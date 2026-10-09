import pickle

import cosmolkit
import pytest


@pytest.mark.parametrize("smiles", ["", "C[C@H](F)Cl", "F/C=C/F", "c1ccccc1.O", pytest.param(None, id="builder-isotope-max-map")])
def test_native_binary_and_all_python_pickle_protocols_preserve_exact_canonical_state(smiles):
    if smiles is None:
        # The pinned SMILES grammar rejects the entire max/10 prefix.
        # Build this typed boundary value through its canonical checked API.
        builder = cosmolkit.MoleculeBuilder.new()
        carbon = builder.add_atom(cosmolkit.AtomSpec(cosmolkit.Element.C).with_isotope(13).with_atom_map(2147483647).with_explicit_hydrogens(3))
        nitrogen = builder.add_atom(cosmolkit.AtomSpec(cosmolkit.Element.N).with_formal_charge(1).with_explicit_hydrogens(3))
        builder.add_bond(cosmolkit.BondSpec(carbon, nitrogen, cosmolkit.BondOrder.SINGLE))
        molecule = builder.build().sanitize()
    else:
        molecule = cosmolkit.Molecule.from_smiles(smiles)
    source = molecule.to_binary()
    assert isinstance(source, bytes)
    assert source[:12] == b"COSMOL\x00\x00\x02\x00\x00\x00"
    restored = cosmolkit.Molecule.from_binary(source)
    assert restored is not molecule
    assert restored.to_binary() == source
    assert molecule.to_binary() == source
    assert restored.to_smiles() == molecule.to_smiles()
    function, args = molecule.__reduce__()
    assert args == (source,)
    assert function(*args).to_binary() == source
    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        restored = pickle.loads(pickle.dumps(molecule, protocol=protocol))
        assert isinstance(restored, cosmolkit.Molecule)
        assert restored is not molecule
        assert restored.to_binary() == source
        assert molecule.to_binary() == source


def test_native_binary_errors_preserve_canonical_variant_and_carrier_fields():
    assert issubclass(cosmolkit.PickleError, ValueError)
    for data, kind, fields in [
        (b"", "UnexpectedEof", {}),
        (b"\xff", "UnsupportedVersion", {"version": 255}),
        (b"CSMOLPKL\x02\x00\x00\x00\x00\x00", "UnsupportedArchiveVersion", {"major": 2, "minor": 0}),
    ]:
        with pytest.raises(cosmolkit.PickleError) as caught:
            cosmolkit.Molecule.from_binary(data)
        error = caught.value
        assert error.domain == "serialization"
        assert error.kind == kind
        for key, value in fields.items():
            assert getattr(error, key) == value


def test_native_binary_preserves_computed_cip_and_failure_atomicity_after_pickle():
    molecule = cosmolkit.Molecule.from_smiles("C[C@H](F)Cl").with_cip_labels()
    source = molecule.to_binary()
    restored = pickle.loads(pickle.dumps(molecule))
    assert restored.cip_computed() == molecule.cip_computed()
    assert restored.atoms()[1].cip_descriptor() == molecule.atoms()[1].cip_descriptor()
    assert restored.atoms()[1].cip_neighbor_order() == molecule.atoms()[1].cip_neighbor_order()
    assert restored.to_binary() == source


def test_original_max_prefix_smiles_rejection_remains_visible():
    with pytest.raises(cosmolkit.SmilesError, match="invalid atom"):
        cosmolkit.Molecule.from_smiles("[13CH3:2147483647][NH3+]")
