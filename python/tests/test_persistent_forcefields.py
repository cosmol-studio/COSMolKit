"""Local binding regressions; corpus parity belongs to the parity test crate."""

import gc

import cosmolkit as ck
import numpy as np
import pytest


def molecule():
    source = ck.Molecule.from_smiles("CC.CC")
    builder = source.to_builder()
    builder.add_3d_conformer(
        [[0.0, 0.0, 0.0], [1.9, 0.2, 0.0], [5.0, 1.0, 0.0], [6.7, 1.1, 0.3]]
    )
    return builder.build().with_assigned_valence()


@pytest.mark.parametrize("kind", ["mmff", "uff"])
def test_owned_factory_matches_existing_evaluator_and_survives_source(kind):
    mol = molecule()
    handle = getattr(mol, f"{kind}_force_field")()
    result = getattr(mol, f"{kind}_energy_gradient")()
    assert handle.energy() == result.energy()
    # Existing one-shot results are flat xyz; the owned handle exposes N x 3.
    np.testing.assert_array_equal(handle.gradient(), np.asarray(result.gradient()).reshape(4, 3))
    original = handle.positions()
    del mol
    gc.collect()
    assert handle.positions().dtype == np.float64
    assert handle.positions().shape == (4, 3)
    np.testing.assert_array_equal(handle.positions(), original)
    assert np.isfinite(handle.energy())


@pytest.mark.parametrize("kind", ["mmff", "uff"])
def test_coordinate_snapshots_cache_refresh_and_joint_results(kind):
    mol = molecule()
    handle = getattr(mol, f"{kind}_force_field")()
    peer = getattr(mol, f"{kind}_force_field")()
    before = peer.positions()
    energy = handle.energy()
    snapshot = handle.positions()
    snapshot[1, 0] += 0.2
    np.testing.assert_array_equal(handle.positions(), before)
    # Noncontiguous arrays and ordinary sequences both project to Rust rows.
    view = np.zeros((4, 6), dtype=np.float64)
    view[:, ::2] = snapshot
    handle.set_positions_(view[:, ::2])
    assert handle.energy() != energy
    result = handle.energy_gradient()
    assert result.energy == handle.energy()
    np.testing.assert_array_equal(result.gradient, handle.gradient())
    gradient = result.gradient
    gradient[:] = 123.0
    assert not np.all(result.gradient == 123.0)
    handle.set_position_(1, [2.3, 0.2, 0.0])
    assert handle.position(1) == (2.3, 0.2, 0.0)
    np.testing.assert_array_equal(peer.positions(), before)
    assert getattr(mol, f"{kind}_force_field")().energy() == energy


@pytest.mark.parametrize("kind", ["mmff", "uff"])
@pytest.mark.parametrize(
    "operation, arguments, error_kind",
    [
        ("set_position_", (99, [0.0, 0.0, 0.0]), ck.MolecularForceFieldErrorKind.InvalidAtomIndex),
        ("set_position_", (-1, [0.0, 0.0, 0.0]), ck.MolecularForceFieldErrorKind.InvalidAtomIndex),
        ("set_position_", (0, [0.0, float("nan"), 0.0]), ck.MolecularForceFieldErrorKind.NonFiniteCoordinate),
        ("set_position_", (0, [0.0, 0.0]), ck.MolecularForceFieldErrorKind.CoordinateShape),
        ("set_positions_", ([[0.0, 0.0, 0.0]],), ck.MolecularForceFieldErrorKind.CoordinateCount),
        ("set_positions_", ([[0.0, 0.0]],), ck.MolecularForceFieldErrorKind.CoordinateShape),
        ("set_fixed_atoms_", ([0, 99],), ck.MolecularForceFieldErrorKind.InvalidFixedAtom),
        ("set_fixed_atoms_", ([-1],), ck.MolecularForceFieldErrorKind.InvalidFixedAtom),
    ],
)
def test_invalid_updates_preserve_positions_energy_and_pins(kind, operation, arguments, error_kind):
    handle = getattr(molecule(), f"{kind}_force_field")()
    handle.set_fixed_atoms_([1])
    before, energy = handle.positions(), handle.energy()
    with pytest.raises(ck.ForceFieldError) as caught:
        getattr(handle, operation)(*arguments)
    assert caught.value.kind == error_kind
    assert isinstance(caught.value.kind, ck.MolecularForceFieldErrorKind)
    np.testing.assert_array_equal(handle.positions(), before)
    assert handle.fixed_atoms() == (1,)
    assert handle.energy() == energy


@pytest.mark.parametrize("kind", ["mmff", "uff"])
def test_fixed_atom_dragging_and_minimization(kind):
    handle = getattr(molecule(), f"{kind}_force_field")()
    handle.set_fixed_atoms_([1, 0, 1])
    assert handle.fixed_atoms() == (0, 1)
    handle.set_fixed_atoms_([0])
    handle.set_position_(0, [0.1, 0.0, 0.0])
    np.testing.assert_array_equal(handle.gradient()[0], [0.0, 0.0, 0.0])
    before = handle.positions()
    zero = handle.minimize_(max_iterations=0)
    assert zero.iterations == 0
    np.testing.assert_array_equal(handle.positions(), before)
    outcome = handle.minimize_with_params_(ck.ForceFieldMinimizeParams(max_iterations=1))
    assert outcome.iterations == 1
    assert outcome.energy == handle.energy()
    assert isinstance(outcome.converged, bool)
    np.testing.assert_array_equal(handle.positions()[0], before[0])
    handle.set_fixed_atoms_([])
    assert handle.fixed_atoms() == ()


def test_params_explicit_factories_and_missing_conformer_errors():
    mmff = ck.MmffForceFieldParams()
    uff = ck.UffForceFieldParams()
    minimize = ck.ForceFieldMinimizeParams()
    assert (mmff.conformer_id, mmff.mmff_variant, mmff.non_bonded_threshold) == (None, "MMFF94", 100.0)
    assert (uff.conformer_id, uff.vdw_threshold) == (None, 10.0)
    assert (minimize.max_iterations, minimize.force_tolerance, minimize.energy_tolerance) == (200, 1e-4, 1e-6)
    with pytest.raises(AttributeError):
        setattr(mmff, "mmff_variant", "MMFF94s")
    mol = molecule()
    assert mol.mmff_force_field_with_params(mmff).energy() == mol.mmff_force_field().energy()
    assert mol.uff_force_field_with_params(uff).energy() == mol.uff_force_field().energy()
    assert mol.mmff_force_field(mmff_variant="MMFF94s").energy() == mol.mmff_force_field_with_params(
        ck.MmffForceFieldParams(mmff_variant="MMFF94s")
    ).energy()
    empty = ck.Molecule.from_smiles("CC")
    for kind, error in [("mmff", ck.MmffForceFieldError), ("uff", ck.UffForceFieldError)]:
        with pytest.raises(error) as caught:
            getattr(empty, f"{kind}_force_field")()
        assert caught.value.kind == ck.MolecularForceFieldErrorKind.MissingConformer
        with pytest.raises(error) as caught:
            getattr(mol, f"{kind}_force_field")(conformer_id=99)
        assert caught.value.kind == ck.MolecularForceFieldErrorKind.MissingConformer
        assert caught.value.requested == 99


@pytest.mark.parametrize("kind", ["mmff", "uff"])
def test_energy_tolerance_keyword_and_params_share_source_optimizer(kind):
    mol = molecule()
    baseline = getattr(mol, f"{kind}_force_field")()
    expected = baseline.minimize_(max_iterations=2)
    for energy_tolerance in (0.0, 1e-12, 1e6):
        keyword = getattr(mol, f"{kind}_force_field")()
        explicit = getattr(mol, f"{kind}_force_field")()
        params = ck.ForceFieldMinimizeParams(
            max_iterations=2, force_tolerance=1e-4, energy_tolerance=energy_tolerance
        )
        assert params.energy_tolerance == energy_tolerance
        a = keyword.minimize_(max_iterations=2, energy_tolerance=energy_tolerance)
        b = explicit.minimize_with_params_(params)
        assert (a.converged, a.iterations, a.energy) == (
            expected.converged, expected.iterations, expected.energy
        ) == (b.converged, b.iterations, b.energy)
        np.testing.assert_array_equal(keyword.positions(), baseline.positions())
        np.testing.assert_array_equal(explicit.positions(), baseline.positions())


@pytest.mark.parametrize("kind", ["mmff", "uff"])
@pytest.mark.parametrize("force_tolerance", [0.0, -1.0, float("nan")])
def test_typed_tolerance_errors_retain_cause_and_state(kind, force_tolerance):
    handle = getattr(molecule(), f"{kind}_force_field")()
    before, energy = handle.positions(), handle.energy()
    with pytest.raises(ck.ForceFieldError) as caught:
        handle.minimize_(force_tolerance=force_tolerance)
    assert caught.value.kind == ck.MolecularForceFieldErrorKind.InvalidTolerance
    assert caught.value.component == 0
    assert caught.value.__cause__ is not None
    np.testing.assert_array_equal(handle.positions(), before)
    assert handle.energy() == energy
