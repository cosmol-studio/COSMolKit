"""Small API-boundary regressions; no corpus or external reference files."""
import ast
from pathlib import Path

import cosmolkit as ck
import numpy as np
import pytest


@pytest.mark.parametrize("cls,field,text,member", [
    (ck.SdfReadParams, "coordinate_mode", "require_3d", ck.SdfCoordinateMode.Require3D),
    (ck.MolBlockWriteParams, "format", "v3000", ck.SdfFormat.V3000),
    (ck.Mol2ReadParams, "variant", "corina", ck.Mol2Type.Corina),
    (ck.ValenceParams, "model", "rdkit_like", ck.ValenceModel.RdkitLike),
    (ck.Coordinate2DInputParams, "z_policy", "require_zero", ck.CoordinateZPolicy.RequireZero),
])
def test_enum_inputs_share_constructor_assignment_and_validation(cls, field, text, member):
    value = cls(**{field: text})
    assert type(getattr(value, field)) is type(member)
    assert getattr(value, field) == member
    assert getattr(cls(**{field: member}), field) == member
    for spelling, expected in type(member)._enum_string_values:
        setattr(value, field, spelling)
        assert type(getattr(value, field)) is type(member)
        assert getattr(value, field) == expected
    previous = getattr(value, field)
    for bad, error in [("unknown_value", ValueError), (123, TypeError), (object(), TypeError), (None, TypeError)]:
        with pytest.raises(error):
            cls(**{field: bad})
        with pytest.raises(error):
            setattr(value, field, bad)
        assert getattr(value, field) == previous


def test_enum_strings_forward_through_actual_sdf_operations(tmp_path: Path):
    text = ck.Molecule.from_smiles("CO").to_sdf()
    expected = ck.Molecule.from_sdf(text, coordinate_mode=ck.SdfCoordinateMode.Require3D)
    actual = ck.Molecule.from_sdf(text, coordinate_mode="require_3d")
    assert actual.to_smiles() == expected.to_smiles() == "CO"
    np.testing.assert_array_equal(actual.coordinates_3d(), expected.coordinates_3d())
    assert "V3000" in actual.to_sdf(format="v3000")
    path = tmp_path / "one.sdf"
    path.write_text(text)
    assert ck.Molecule.read_sdf(str(path), coordinate_mode="require_3d").to_smiles() == "CO"
    first = ck.MoleculeBatch.from_sdf_records(text, coordinate_mode="require_3d")[0]
    assert first is not None
    assert first.to_smiles() == "CO"
    with pytest.raises(ValueError):
        _ = ck.Molecule.from_sdf(text, coordinate_mode="automatic")


def test_dynamic_enum_inputs_use_registered_members_not_new_conversion_tables():
    builder = ck.MoleculeBuilder.new()
    first = builder.add_atom(ck.AtomSpec.new(ck.Element.C))
    second = builder.add_atom(ck.AtomSpec.new(ck.Element.C))
    bond = builder.add_bond(ck.BondSpec.new(first, second, "single"))
    assert builder.build().to_smiles() == "CC"
    builder.set_bond_order(bond, "double")
    assert builder.build().to_smiles() == "C=C"
    assert repr(ck.BondSpec(first, second, "double")) == repr(ck.BondSpec(first, second, ck.BondOrder.DOUBLE))
    with pytest.raises(ValueError):
        ck.BondSpec(first, second, "not_a_bond")
    with pytest.raises(TypeError):
        ck.BondSpec(first, second, object())  # pyright: ignore[reportArgumentType] -- invalid input regression
    for name, member in ck.ResidueInfoKind.__members__.items():
        assert ck.expand_one_letter("A", name.lower()) == ck.expand_one_letter("A", member)
    params = ck.BioReadParams(format="pdb")
    assert params.format == ck.BioCoordinateFormat.Pdb
    params.format = "mmcif"
    assert params.format == ck.BioCoordinateFormat.Mmcif
    with pytest.raises(ValueError):
        params.format = "invalid_format"
    assert params.format == ck.BioCoordinateFormat.Mmcif


@pytest.mark.parametrize("method", ["substruct_match", "has_substruct_match"])
def test_configured_single_and_boolean_matching_use_canonical_call_forms(method):
    molecule = ck.Molecule.from_smiles("F[C@](Cl)(Br)I")
    query = ck.parse_smarts("F[C@@](Cl)(Br)I")
    call = getattr(molecule, method)
    before = molecule.to_smiles()
    params = ck.SubstructMatchParams(use_chirality=True)
    assert call(query, params) == call(query, use_chirality=True)
    assert not call(query, params)
    assert call(query)
    with pytest.raises(TypeError, match="mutually exclusive"):
        call(query, params, use_chirality=False)
    assert params.use_chirality
    assert molecule.to_smiles() == before


def test_targeted_option_names_and_fingerprint_discovery():
    for cls in (ck.SmilesParseParams, ck.SdfReadParams, ck.Mol2ReadParams, ck.InchiReadParams):
        params = cls(remove_hs=False)
        assert not params.remove_hs
        assert not hasattr(params, "remove_hydrogens")
    writer = ck.SmilesWriteParams(kekule=True, isomeric_smiles=False)
    assert writer.kekule and not writer.isomeric_smiles
    assert not hasattr(writer, "do_kekule")
    assert not hasattr(writer, "do_isomeric_smiles")
    molecule = ck.Molecule.from_smiles("CCO")
    for family in ("layered", "pattern", "morgan", "atom_pair", "topological", "topological_torsion", "maccs"):
        assert callable(getattr(molecule, f"fingerprint_{family}"))
        assert not hasattr(molecule, f"{family}_fingerprint")
    assert callable(ck.inchi_to_key)
    assert not hasattr(ck, "inchi_to_inchi_key")
    assert hasattr(molecule, "remove_hydrogens_")


def test_generated_stub_describes_enum_input_but_typed_output():
    stub = Path(__file__).resolve().parents[1] / "cosmolkit.pyi"
    tree = ast.parse(stub.read_text())
    params = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "SdfReadParams")
    getter = next(node for node in params.body if isinstance(node, ast.FunctionDef) and node.name == "coordinate_mode" and any(isinstance(d, ast.Name) and d.id == "property" for d in node.decorator_list))
    setter = next(node for node in params.body if isinstance(node, ast.FunctionDef) and node.name == "coordinate_mode" and any(isinstance(d, ast.Attribute) and d.attr == "setter" for d in node.decorator_list))
    assert getter.returns is not None
    annotation = setter.args.args[-1].annotation
    assert annotation is not None
    assert ast.unparse(getter.returns) == "SdfCoordinateMode"
    assert ast.unparse(annotation) == "SdfCoordinateMode | builtins.str"
