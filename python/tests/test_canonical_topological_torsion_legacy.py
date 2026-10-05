"""Pinned source fixed legacy public protocol proposals."""
import ast
from pathlib import Path
import pytest
import cosmolkit as ck

METHODS = ["legacy_topological_torsion_sparse_count_fingerprint", "legacy_topological_torsion_count_fingerprint", "legacy_topological_torsion_fingerprint"]

def test_legacy_source_default_unfolded_and_hashed_fixed_values():
    m = ck.Molecule.from_smiles("CCCCO")
    params = ck.LegacyTopologicalTorsionParams()
    assert (params.torsion_atom_count, params.include_chirality, params.fp_size, params.bits_per_entry) == (4, False, 2048, 4)
    assert params.from_atoms is None and params.ignore_atoms is None and params.custom_atom_invariants is None
    value = m.legacy_topological_torsion_sparse_count_fingerprint()
    assert value.length() == (1 << 36) - 1
    assert value.nonzero_elements() == {4437590048: 1, 12893306913: 1}
    value = m.legacy_topological_torsion_count_fingerprint_with_params(ck.LegacyTopologicalTorsionParams(fp_size=1000))
    assert value.length() == 1000 and value.nonzero_elements() == {24: 1, 288: 1}
    for method in METHODS:
        a, b = getattr(m, method)(), getattr(m, method + "_with_params")(params)
        if hasattr(a, "nonzero_elements"):
            assert a.length() == b.length() and a.nonzero_elements() == b.nonzero_elements()
        else:
            assert a.n_bits() == b.n_bits() and a.on_bits() == b.on_bits()

@pytest.mark.parametrize("method", METHODS)
def test_legacy_empty_roots_and_original_short_custom_error(method):
    m = ck.Molecule.from_smiles("CCCCO"); before = m.to_smiles()
    value = getattr(m, method + "_with_params")(ck.LegacyTopologicalTorsionParams(from_atoms=[]))
    assert value.nonzero_elements() == {} if hasattr(value, "nonzero_elements") else value.on_bits() == []
    with pytest.raises(ck.TopologicalTorsionReadError, match="bad atomInvariants size") as error:
        getattr(m, method + "_with_params")(ck.LegacyTopologicalTorsionParams(custom_atom_invariants=[1]))
    assert error.value.kind == "Generator" and error.value.__cause__ is not None
    assert m.to_smiles() == before
    assert m.legacy_topological_torsion_sparse_count_fingerprint().nonzero_elements() == {4437590048: 1, 12893306913: 1}

@pytest.mark.parametrize("entry", [1, 2, 4, 6])
def test_legacy_four_and_nonfour_source_threshold_projection(entry):
    m = ck.Molecule.from_smiles("CCCCCCCCCCCC"); inv = [7] * 12
    counts = m.legacy_topological_torsion_count_fingerprint_with_params(ck.LegacyTopologicalTorsionParams(fp_size=16, custom_atom_invariants=inv)).nonzero_elements()
    assert len(counts) == 1
    block, count = next(iter(counts.items())); assert count == 9
    value = m.legacy_topological_torsion_fingerprint_with_params(ck.LegacyTopologicalTorsionParams(fp_size=16 * entry, bits_per_entry=entry, custom_atom_invariants=inv))
    assert value.n_bits() == 16 * entry and value.on_bits() == [block * entry + i for i in range(entry)]

def test_legacy_owned_inputs_conversion_and_generated_constructor_protocol():
    roots, ignored, inv = [0], [1], [17, 18, 19, 20, 21]
    p = ck.LegacyTopologicalTorsionParams(from_atoms=roots, ignore_atoms=ignored, custom_atom_invariants=inv)
    roots.append(2); ignored.clear(); inv.clear()
    assert p.from_atoms == [0] and p.ignore_atoms == [1] and p.custom_atom_invariants == [17, 18, 19, 20, 21]
    copy = p.from_atoms; copy.clear(); assert p.from_atoms == [0]
    with pytest.raises(AttributeError): p.fp_size = 10
    for keyword, value, error in [("fp_size", -1, OverflowError), ("torsion_atom_count", 4294967296, OverflowError), ("bits_per_entry", 1.5, TypeError), ("from_atoms", [-1], OverflowError)]:
        with pytest.raises(error): ck.LegacyTopologicalTorsionParams(**{keyword: value})
    tree = ast.parse((Path(__file__).resolve().parents[1] / "cosmolkit.pyi").read_text())
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "LegacyTopologicalTorsionParams")
    methods = {node.name: node for node in cls.body if isinstance(node, ast.FunctionDef)}
    fields = ["torsion_atom_count", "include_chirality", "fp_size", "bits_per_entry", "from_atoms", "ignore_atoms", "custom_atom_invariants"]
    assert set(methods) == {"__new__", *fields}
    for field in fields:
        assert [ast.unparse(d) for d in methods[field].decorator_list] == ["property"]
    constructor = methods["__new__"]
    assert [a.arg for a in constructor.args.kwonlyargs] == fields
    assert [ast.unparse(v) for v in constructor.args.kw_defaults] == ["4", "False", "2048", "4", "None", "None", "None"]
