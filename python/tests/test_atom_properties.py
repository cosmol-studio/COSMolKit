"""Typed atom metadata and opt-in reaction provenance, with no external oracle."""
import math
import pytest
import cosmolkit as ck


@pytest.mark.parametrize("kind,value", [("int",42), ("uint",2**32-1), ("bool",True), ("bool",False), ("double",1.25), ("double",-0.0), ("string","α\0β"), ("int_vector",[1,-2]), ("string_vector",["a","β"]), ("int_vector",[])])
def test_native_atom_properties_and_cow(kind, value):
    original = ck.Molecule.from_smiles("CCO")
    tagged = getattr(original, f"with_atom_property_{kind}")(0, "tracking_id", value)
    actual = getattr(tagged, f"atom_property_{kind}")(0, "tracking_id")
    assert actual == value
    assert type(actual) is type(value)
    if isinstance(value, float):
        assert math.copysign(1.0, actual) == math.copysign(1.0, value)
    assert getattr(original, f"atom_property_{kind}")(0, "tracking_id") is None
    copied = tagged.with_atom_property_string(0, "copy_note", "copy")
    tagged.set_atom_property_int_(0, "tracking_id", 7)
    assert getattr(copied, f"atom_property_{kind}")(0, "tracking_id") == value
    assert tagged.atom_property_int(0, "tracking_id") == 7
    assert tagged.to_smiles() == original.to_smiles()


def test_atom_property_validation_and_failed_write_atomicity():
    mol = ck.Molecule.from_smiles("C").with_atom_property_int(0, "tracking_id", 42)
    for bad in [None, {}, [True], [1, "a"], [1.5]]:
        with pytest.raises(TypeError):
            mol.set_atom_property_int_(0, "tracking_id", bad)
    with pytest.raises(OverflowError):
        mol.set_atom_property_int_(0, "tracking_id", 2**40)
    for atom, key in [(9, "tracking_id"), (0, "_CIPCode"), (0, "")]:
        with pytest.raises(ck.OperationError):
            mol.set_atom_property_int_(atom, key, 7)
    assert mol.atom_property_int(0, "tracking_id") == 42


def test_reaction_copy_atom_properties_uses_input_origins():
    carbon = ck.Molecule.from_smiles("C").with_atom_property_int(0, "tracking_id", 42)
    oxygen = ck.Molecule.from_smiles("O").with_atom_property_int(0, "tracking_id", 84)
    rxn = ck.Reaction.from_smirks("[C:1].[O:2]>>[C:1][O:2]")
    params = ck.ReactionRunParams(copy_atom_properties=True)
    assert params.copy_atom_properties is True
    product = rxn.run([carbon, oxygen], params)[0][0]
    assert [product.atom_property_int(i, "tracking_id") for i in range(2)] == [42, 84]
    params.copy_atom_properties = False
    plain = rxn.run([carbon, oxygen], params)[0][0]
    assert [plain.atom_property_int(i, "tracking_id") for i in range(2)] == [None, None]
    assert product.to_smiles() == plain.to_smiles() == "CO"
    assert carbon.atom_property_int(0, "tracking_id") == 42
    assert oxygen.atom_property_int(0, "tracking_id") == 84
