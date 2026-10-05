"""Proposed canonical value transform defaults, result and error protocols."""
import ast
from pathlib import Path
import pytest
import cosmolkit as ck

def test_source_carbon_chiral_codes_and_explicit_result_molecule():
    source = ck.Molecule.from_smiles("C[C@H](F)Cl")
    ordinary = source.with_atom_pair_atom_code(1)
    assert ordinary.code == 35
    assert isinstance(ordinary, ck.AtomPairAtomCodeResult)
    assert isinstance(ordinary.molecule, ck.Molecule)
    assert ordinary.molecule.num_atoms() == source.num_atoms() == 4
    for legacy in [True, False]:
        result = source.with_atom_pair_atom_code(1, include_chirality=True, use_legacy_stereo_perception=legacy)
        assert result.code == 547
        assert result.molecule.with_atom_pair_atom_code(1, include_chirality=True, use_legacy_stereo_perception=False).code == 547
        assert result.molecule.with_atom_pair_atom_code(1).code == 35
    assert source.with_atom_pair_atom_code(1).code == 35

def test_source_unsigned_branch_subtract_and_original_pi_bits():
    source = ck.Molecule.from_smiles("C=CC(=O)O")
    assert [source.with_atom_pair_atom_code(i).code for i in range(5)] == [41, 42, 43, 105, 97]
    assert source.with_atom_pair_atom_code(2, 2).code == 41
    assert source.with_atom_pair_atom_code(2, (1 << 32) - 1).code == 40

def test_frozen_result_owns_public_molecule_projection():
    source = ck.Molecule.from_smiles("CC")
    result = source.with_atom_pair_atom_code(0)
    for field in ["code", "molecule", "other"]:
        with pytest.raises(AttributeError): setattr(result, field, None)
    with pytest.raises(TypeError): ck.AtomPairAtomCodeResult()
    first, second = result.molecule, result.molecule
    del source, result
    assert first.with_atom_pair_atom_code(0).code == second.with_atom_pair_atom_code(0).code == 33

def test_boundary_error_remains_structured_and_source_molecule_usable():
    source = ck.Molecule.from_smiles("CC")
    with pytest.raises(ck.OperationError) as raised: source.with_atom_pair_atom_code(2, include_chirality=True, use_legacy_stereo_perception=False)
    error = raised.value
    assert error.domain == "operation" and error.kind == "AtomCode"
    assert error.__cause__ is not None
    assert source.with_atom_pair_atom_code(0).code == 33
    for subtract in [-1, 1 << 32]:
        with pytest.raises(OverflowError): source.with_atom_pair_atom_code(0, subtract)
    with pytest.raises(OverflowError): source.with_atom_pair_atom_code(-1)
    with pytest.raises(TypeError): source.with_atom_pair_atom_code("0")

def test_original_generated_stub_exposes_only_explicit_value_transform():
    tree = ast.parse((Path(__file__).resolve().parents[1] / "cosmolkit.pyi").read_text())
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    result = {node.name: node for node in classes["AtomPairAtomCodeResult"].body if isinstance(node, ast.FunctionDef)}
    assert set(result) == {"code", "molecule"}
    assert ast.unparse(result["code"].returns) == "builtins.int"
    assert ast.unparse(result["molecule"].returns) == "Molecule"
    for declaration in result.values(): assert [ast.unparse(d) for d in declaration.decorator_list] == ["property"]
    methods = {node.name: node for node in classes["Molecule"].body if isinstance(node, ast.FunctionDef)}
    method = methods["with_atom_pair_atom_code"]
    assert [arg.arg for arg in method.args.args] == ["self", "atom_id", "branch_subtract", "include_chirality", "use_legacy_stereo_perception"]
    assert [ast.literal_eval(d) for d in method.args.defaults] == [0, False, True]
    assert ast.unparse(method.returns) == "AtomPairAtomCodeResult"
    assert "get_atom_code" not in methods and "atom_pair_atom_code" not in methods
