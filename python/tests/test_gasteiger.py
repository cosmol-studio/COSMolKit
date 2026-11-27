"""Gasteiger charges and checked options on the packaged Python API."""

import math

import cosmolkit as ck
import pytest


def test_charge_rows_keep_disconnected_values_and_input_properties() -> None:
    mol = ck.Molecule.from_smiles("[Fe]C.CCO")
    mol.set_atom_property_double_(2, "user_charge", 42.0)
    rows = mol.gasteiger_charges()

    assert math.isnan(rows[0].atom_charge) and rows[0].attached_hydrogen_charge == 0
    assert math.isnan(rows[1].atom_charge) and math.isnan(
        rows[1].attached_hydrogen_charge
    )
    assert rows[4].atom_charge == -0.3966637147812598
    assert rows[4].attached_hydrogen_charge == 0.21002230656786322
    assert math.isnan(mol.max_partial_charge())
    assert mol.peoe_vsa_14() == mol.peoe_vsa()[13]

    with pytest.raises(ck.DescriptorReadError) as exc:
        mol.gasteiger_charges(throw_on_param_failure=True)
    assert isinstance(exc.value.__cause__.__cause__, ck.GasteigerChargeError)
    assert mol.atom_property_double(2, "user_charge") == 42.0
    for index in range(len(rows)):
        assert mol.atom_property_double(index, "_GasteigerCharge") is None
        assert mol.atom_property_double(index, "_GasteigerHCharge") is None


def test_invalid_iterations_cannot_replace_a_valid_setting() -> None:
    acetate = ck.Molecule.from_smiles("CC(=O)[O-]")
    params = ck.GasteigerChargeParams(iterations=0)
    failures = (TypeError, ValueError, OverflowError)
    for value in [-1, 1.25, math.inf, math.nan, 2**31]:
        with pytest.raises(failures):
            ck.GasteigerChargeParams(iterations=value)
        with pytest.raises(failures):
            params.iterations = value
        with pytest.raises(failures):
            acetate.gasteiger_charges(iterations=value)
        assert params.iterations == 0

    rows = acetate.gasteiger_charges_with_params(params)
    assert [row.atom_charge for row in rows] == [0, 0, -0.5, -0.5]
    assert all(row.attached_hydrogen_charge == 0 for row in rows)

    for configured in (
        acetate.gasteiger_charges(params),
        acetate.gasteiger_charges(iterations=0),
    ):
        assert [row.atom_charge for row in configured] == [
            row.atom_charge for row in rows
        ]
        assert all(row.attached_hydrogen_charge == 0 for row in configured)

    with pytest.raises(TypeError):
        acetate.gasteiger_charges(params, throw_on_param_failure=False)
    with pytest.raises(TypeError):
        acetate.gasteiger_charges(params, iterations=None)
