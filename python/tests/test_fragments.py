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
        _ = ck.Molecule.from_smiles("").largest_fragment()
    assert caught.value.variant.name == "EmptyFragments"


def test_cut_attachment_sides_and_in_place_result() -> None:
    source = ck.Molecule.from_smiles("CCO")
    target = ck.Molecule.from_smiles("CCO")
    cuts = [ck.BondCut.at_bond(1, begin_label=101, end_label=202)]

    product = source.with_bonds_cut(cuts)
    assert source.to_smiles() == "CCO"
    assert target.cut_bonds_(cuts) is None
    assert target.to_smiles() == product.to_smiles() == "[101*]CC.[202*]O"
    assert [part.to_smiles() for part in product.fragments()] == ["[101*]CC", "[202*]O"]


def test_invalid_cut_reports_the_source_bond_without_modifying_the_receiver() -> None:
    mol = ck.Molecule.from_smiles("CCO")
    before = mol.to_binary()
    cuts = [ck.BondCut.at_bond(0), ck.BondCut.at_bond(1, 65536, 8)]

    for apply in (mol.with_bonds_cut, mol.cut_bonds_):
        with pytest.raises(ck.OperationError) as caught:
            _ = apply(cuts)
        detail = caught.value.bond_cut
        assert isinstance(detail, ck.BondCutError)
        assert detail.variant.name == "IsotopeOutOfRange"
        assert (detail.bond, detail.attachment, detail.label) == (1, 1, 65536)
        assert isinstance(caught.value.__cause__, ck.BondCutError)
        assert mol.to_binary() == before

    with pytest.raises(ck.BondCutError) as incomplete:
        _ = ck.BondCut.at_bond(0, end_label=0)
    assert incomplete.value.variant.name == "IncompleteLabels"
