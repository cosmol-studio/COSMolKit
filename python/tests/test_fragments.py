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


def test_brics_sites_convert_to_selections_in_rule_order() -> None:
    source = ck.Molecule.from_smiles("CCCOCCC")
    sites = source.find_brics_bonds()
    assert [(s.bond, s.ends) for s in sites] == [
        (2, ((3, 3), (2, 4))), (3, ((3, 3), (4, 4))),
    ]
    cuts = [site.to_bond_cut() for site in sites]
    assert source.with_bonds_cut(cuts[:1]).to_smiles() == "[3*]OCCC.[4*]CCC"
    assert source.cut_bonds_(cuts) is None
    assert [p.to_smiles() for p in source.fragments()] == ["[4*]CCC", "[3*]O[3*]", "[4*]CCC"]


def test_regular_r_brics_on_octane() -> None:
    mol = ck.Molecule.from_smiles("CCCCCCCC")
    (site,) = mol.find_r_brics_bonds()
    assert isinstance(site, ck.RBricsBondMatch)
    assert site.ends == ((3, 20), (4, 21))
    assert sorted(mol.decompose_r_brics()) == ["[20*]CCCC", "[21*]CCCC"]
    assert mol.find_brics_bonds() == []
    assert mol.with_bonds_cut([site.to_bond_cut()]).to_smiles() == "[20*]CCCC.[21*]CCCC"


def test_recursive_brics_returns_unique_terminal_smiles() -> None:
    source = ck.Molecule.from_smiles("CCCOCCC")
    before = source.to_binary()
    assert sorted(source.decompose_brics()) == ["[3*]O[3*]", "[4*]CCC"]
    assert source.to_binary() == before


def test_decomposition_reports_unsanitizable_product() -> None:
    source = ck.Molecule.from_smiles("O=C([O-])C1=CC=C(N2N=C3[N]4=C2C5=CC=CC=C5[O-][Fe+3]4[O-]C6=CC=CC=C63)C=C1")
    with pytest.raises(ck.BricsDecompositionError) as caught:
        _ = source.decompose_brics()
    assert caught.value.variant.name == "UnkeyedProduct"
    assert isinstance(caught.value.__cause__, ck.SanitizeError)

    with pytest.raises(ck.RecapDecompositionError) as recap_error:
        _ = source.decompose_recap()
    assert recap_error.value.variant.name == "UnkeyedProduct"
    assert isinstance(recap_error.value.__cause__, ck.SanitizeError)



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
    cuts = [ck.BondCut.at_bond(0), ck.BondCut.at_bond(2)]

    for apply in (mol.with_bonds_cut, mol.cut_bonds_):
        with pytest.raises(ck.OperationError) as caught:
            _ = apply(cuts)
        detail = caught.value.bond_cut
        assert isinstance(detail, ck.BondCutError)
        assert detail.variant.name == "BondOutOfRange"
        assert (detail.bond, detail.count) == (2, 2)
        assert isinstance(caught.value.__cause__, ck.BondCutError)
        assert mol.to_binary() == before

    with pytest.raises(ck.BondCutError) as incomplete:
        _ = ck.BondCut.at_bond(0, end_label=0)
    assert incomplete.value.variant.name == "IncompleteLabels"


def test_recap_shared_product_and_missing_key() -> None:
    source = ck.Molecule.from_smiles("CC(=O)Nc1ccc(NC(=O)C)cc1")
    before = source.to_binary()
    graph = source.decompose_recap()
    intermediate = "*Nc1ccc(NC(C)=O)cc1"
    assert sorted(graph.parents("*C(C)=O")) == [intermediate, graph.root]
    assert sorted(graph.children(graph.root)) == ["*C(C)=O", intermediate]
    assert sorted(graph.leaves()) == ["*C(C)=O", "*Nc1ccc(N*)cc1"]
    assert graph.children("*C(C)=O") == []
    assert source.to_binary() == before
    with pytest.raises(ck.RecapLookupError) as missing:
        _ = graph.children("unknown")
    assert missing.value.product == "unknown"

    uncut = ck.Molecule.from_smiles("C").decompose_recap()
    assert uncut.root == "C"
    assert uncut.children(uncut.root) == uncut.parents(uncut.root) == uncut.leaves() == []

