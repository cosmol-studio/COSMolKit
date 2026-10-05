import cosmolkit
import pytest

def test_editing_commit_boundary_matches_sanitize_behavior():
    invalid_editor = cosmolkit.Molecule.from_smiles("CC").to_builder()
    oxygen_a = invalid_editor.add_atom(cosmolkit.AtomSpec(cosmolkit.Element.O))
    oxygen_b = invalid_editor.add_atom(cosmolkit.AtomSpec(cosmolkit.Element.O))
    invalid_editor.add_bond(cosmolkit.BondSpec(1, oxygen_a, cosmolkit.BondOrder.DOUBLE))
    invalid_editor.add_bond(cosmolkit.BondSpec(1, oxygen_b, cosmolkit.BondOrder.DOUBLE))

    with pytest.raises(cosmolkit.OperationError, match="sanitization failed") as error:
        _ = invalid_editor.build().sanitize()
    assert isinstance(error.value, ValueError)
    assert error.value.domain == "operation"
    assert error.value.kind == "Sanitize"
    assert error.value.__cause__ is not None

    edited = invalid_editor.build()
    assert len(edited) == 4
    assert [bond.order() for bond in edited.bonds()][-2:] == [
        cosmolkit.BondOrder.DOUBLE,
        cosmolkit.BondOrder.DOUBLE,
    ]

    valid_editor = cosmolkit.Molecule.from_smiles("CC").to_builder()
    oxygen = valid_editor.add_atom(cosmolkit.AtomSpec(cosmolkit.Element.O))
    valid_editor.add_bond(cosmolkit.BondSpec(1, oxygen, cosmolkit.BondOrder.SINGLE))
    valid = valid_editor.build().sanitize()
    assert valid.to_smiles_with_params(cosmolkit.SmilesWriteParams(canonical=False)) == "CCO"

    metal_editor = cosmolkit.Molecule.from_smiles("C").to_builder()
    hg = metal_editor.add_atom(cosmolkit.AtomSpec(cosmolkit.Element.HG))
    metal = metal_editor.build()
    assert metal.atoms()[hg].atomic_number() == 80
