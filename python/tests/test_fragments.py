"""Representative component ownership and exception-conversion checks."""
import cosmolkit as ck
import pytest


def test_fragment_bindings_order_tie_and_nonmutation():
    mol = ck.Molecule.from_smiles("CC.OO.[Na+]")
    before = mol.to_binary()
    assert [f.to_smiles() for f in mol.fragments()] == ["CC", "OO", "[Na+]"]
    assert mol.largest_fragment().to_smiles() == "OO"
    assert mol.to_binary() == before
    assert ck.Molecule.from_smiles("").fragments() == []
    with pytest.raises(ck.OperationError) as caught:
        ck.Molecule.from_smiles("").largest_fragment()
    assert caught.value.kind == "EmptyFragments"
