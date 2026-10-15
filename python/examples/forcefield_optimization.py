"""UFF/MMFF optimization and persistent distance/angle/torsion/position restraints.

Run from the repository root: .venv/bin/python python/examples/forcefield_optimization.py
The installed extension must include the conformer capability (included in full).
"""

from __future__ import annotations

import numpy as np

import cosmolkit as ck

ethanol_3d = """ethanol_3d
  COSMolKit      3D

  9  8  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6000    0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6000   -0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9000   -0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7000    0.0000    1.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9000    1.2000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
"""

mol = ck.mol_from_mol(ethanol_3d, coordinate_mode="require_3d", remove_hs=False)
start = mol.coordinates_3d().copy()

if mol.uff_has_all_molecule_params():
    result = mol.with_uff_optimized(max_iterations=200)
    optimized = result.molecule()

    print("UFF converged:", not result.needs_more())
    print("UFF status code:", result.status_code())
    print("UFF energy:", result.energy())
    print("coordinates changed:", not np.allclose(start, optimized.coordinates_3d()))

if mol.mmff_has_all_molecule_params():
    result = mol.with_mmff_optimized(max_iterations=200)
    print("MMFF94 converged:", not result.needs_more())
    print("MMFF94 status code:", result.status_code())

# These values are explicit example choices, not implicit numerical defaults.
for model in ("UFF", "MMFF94", "MMFF94s"):
    field = mol.uff_force_field() if model == "UFF" else mol.mmff_force_field(mmff_variant=model)

    # Either pass a parameter object or keyword-only configuration to the same method.
    distance = ck.ForceFieldDistanceConstraintParams(
        relative=False, min_distance=1.4, max_distance=1.6, force_constant=50.0,
    )
    distance.force_constant = 100.0  # Configuration fields are writable.
    field.add_distance_constraint_(0, 1, distance)
    # Relative limits capture the current angle/dihedral when the constraint is added.
    field.add_angle_constraint_(
        0, 1, 2, relative=True, min_angle_deg=-5.0, max_angle_deg=5.0, force_constant=1.0,
    )
    field.add_torsion_constraint_(
        3, 0, 1, 2,
        relative=True, min_dihedral_deg=-10.0, max_dihedral_deg=10.0, force_constant=1.0,
    )
    # Distances/displacements use Angstroms; angles/dihedrals use degrees.
    # This soft O restraint keeps its original anchor after coordinate updates.
    field.add_position_constraint_(2, max_displacement=0.25, force_constant=10.0)
    field.set_fixed_atoms_([0])  # Unlike the soft restraint, atom 0 cannot move.
    pinned = field.position(0)
    oxygen = field.position(2)
    field.set_position_(2, [oxygen[0] + 0.3, oxygen[1] + 0.1, oxygen[2] + 0.2])
    initial = field.energy_gradient()
    assert initial.gradient.shape == (mol.num_atoms(), 3)
    outcome = field.minimize_()
    assert np.isfinite(outcome.energy)
    assert outcome.energy <= initial.energy
    assert field.position(0) == pinned
    assert np.array_equal(mol.coordinates_3d(), start)

    # Evaluator edits do not change mol. Install the result explicitly in a new value.
    positions = field.positions()
    relaxed = mol.with_3d_coordinates(positions, conformer_id=0)
    assert np.array_equal(relaxed.coordinates_3d(), positions)
    assert np.array_equal(mol.coordinates_3d(), start)
    print(
        f"{model} restrained: energy={initial.energy:.6f} -> {outcome.energy:.6f},",
        f"converged={outcome.converged}, iterations={outcome.iterations}"
    )
