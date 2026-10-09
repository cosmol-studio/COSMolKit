"""AtomPair canonical transport delivery proposal; source owner supplies chemistry."""
import pytest
import cosmolkit as ck

METHODS = ["fingerprint_atom_pair", "fingerprint_atom_pair_sparse", "fingerprint_atom_pair_count", "fingerprint_atom_pair_sparse_count"]

@pytest.mark.parametrize("method", METHODS)
def test_all_forms_fill_original_output_and_preserve_input(method):
    molecule = ck.Molecule.from_smiles("CC(C)O")
    before = molecule.to_smiles()
    output = ck.FingerprintAdditionalOutput()
    output.allocate_atom_counts(); output.allocate_atom_to_bits()
    output.allocate_bit_info_map(); output.allocate_atoms_per_bit()
    params = ck.AtomPairFingerprintParams(generator=ck.AtomPairParams(fp_size=256))
    result = getattr(molecule, method + "_with_params")(params, output)
    ids = set(result.on_bits() if hasattr(result, "on_bits") else result.nonzero_elements())
    assert ids
    assert output.atom_counts() == [3, 3, 3, 3]
    assert len(output.atom_to_bits()) == molecule.num_atoms()
    assert set(output.bit_info_map()) <= ids
    assert set(output.atoms_per_bit()) <= ids
    assert output.bit_paths() is None
    assert molecule.to_smiles() == before
    counts = output.atom_counts(); counts[0] = 0
    assert output.atom_counts()[0] == 3


def test_empty_roots_optional_values_and_frozen_parameters():
    molecule = ck.Molecule.from_smiles("CCO")
    empty = ck.AtomPairFingerprintParams(from_atoms=[])
    assert empty.from_atoms == [] and empty.ignore_atoms is None
    assert empty.conformer_id == -1
    assert molecule.fingerprint_atom_pair_sparse_count_with_params(empty, None).nonzero_elements() == {}
    assert molecule.fingerprint_atom_pair_sparse_count().nonzero_elements()
    empty.conformer_id = 4
    assert empty.conformer_id == 4
    copy = empty.from_atoms; copy.append(0)
    assert empty.from_atoms == []
    assert ck.AtomPairParams(count_bounds=[]).count_bounds == []


def test_source_generator_precondition_retains_typed_cause_and_input():
    molecule = ck.Molecule.from_smiles("CCO")
    before = molecule.to_smiles()
    bad = ck.AtomPairFingerprintParams(generator=ck.AtomPairParams(min_distance=5, max_distance=4))
    with pytest.raises(ck.AtomPairReadError) as caught:
        molecule.fingerprint_atom_pair_with_params(bad, None)
    assert caught.value.domain == "fingerprints"
    assert caught.value.kind == "Generator"
    assert caught.value.__cause__ is not None
    assert "bad distances provided" in str(caught.value)
    assert molecule.to_smiles() == before


def test_custom_invariants_modulo_and_filters_reach_owner():
    molecule = ck.Molecule.from_smiles("CCO")
    first = ck.AtomPairFingerprintParams(from_atoms=[0], ignore_atoms=[2], custom_atom_invariants=[10,20,30])
    # Original fixed source pair-code regression: 10,20,distance1 => 328001.
    assert molecule.fingerprint_atom_pair_sparse_count_with_params(first,None).nonzero_elements() == {328001:1}
    torsion = ck.AtomPairAtomInvariantsGenerator(topological_torsion_correction=True)
    assert torsion.topological_torsion_correction and not torsion.include_chirality
    custom = ck.AtomPairFingerprintParams(atom_invariants_generator=torsion)
    assert molecule.fingerprint_atom_pair_sparse_count_with_params(custom,None).nonzero_elements() != molecule.fingerprint_atom_pair_sparse_count().nonzero_elements()


@pytest.mark.parametrize("chiral", [False, True])
@pytest.mark.parametrize("correction", [False, True])
def test_invariant_generator_source_metadata_repr_and_frozen_protocol(chiral, correction):
    # Pinned AtomPairGenerator.cpp45-55 includes correction in info_string
    # and both independent fields as quoted bools in Boost JSON.
    import ast
    import json
    from pathlib import Path
    generator = ck.AtomPairAtomInvariantsGenerator(
        include_chirality=chiral, topological_torsion_correction=correction)
    expected = f"AtomPairInvariantGenerator topologicalTorsionCorrection={int(correction)}"
    assert generator.info_string() == expected
    assert json.loads(generator.to_json()) == {
        "type": "AtomPairAtomInvGenerator",
        "includeChirality": str(chiral).lower(),
        "topologicalTorsionCorrection": str(correction).lower(),
    }
    # Original modern Python PyAtomPairAtomInvariantsGenerator repr is a
    # thin formatting projection of the exact source information string.
    assert repr(generator) == f"AtomPairAtomInvariantsGenerator({expected})"
    for name in ("include_chirality", "topological_torsion_correction"):
        original = getattr(generator, name)
        setattr(generator, name, not original)
        assert getattr(generator, name) is not original
        setattr(generator, name, original)
    for name in ("info_string", "to_json", "__repr__"):
        with pytest.raises(TypeError): getattr(generator, name)(0)
    stub = ast.parse((Path(__file__).resolve().parents[1] / "cosmolkit.pyi").read_text())
    node = next(n for n in stub.body if isinstance(n, ast.ClassDef) and n.name == "AtomPairAtomInvariantsGenerator")
    methods = {n.name: n for n in node.body if isinstance(n, ast.FunctionDef)}
    assert set(methods) == {"__new__", "include_chirality", "topological_torsion_correction", "info_string", "to_json", "__repr__"}
    for name in ("info_string", "to_json", "__repr__"):
        assert ast.unparse(methods[name].returns) == "builtins.str"
