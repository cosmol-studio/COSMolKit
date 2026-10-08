Force Fields: Current Usage and Persistent API Design
=====================================================

This page separates the current Python API from the proposed persistent
``MolecularForceField`` API. Both MMFF and UFF use coordinates in angstroms,
energies in kcal/mol, and energy gradients in kcal/mol/angstrom.

Current API
-----------

Prepare a molecule with 3D coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Force-field evaluation and optimization require an existing 3D conformer.
They do not generate coordinates or add hydrogens automatically. Supply
coordinates in atom order, with explicit hydrogens when appropriate.

The following small, self-contained water example supplies coordinates
directly; it is not a conformer-generation algorithm:

.. code-block:: python

   import cosmolkit as ck
   import numpy as np

   base = ck.Molecule.from_smiles("O").with_hydrogens()
   builder = base.to_builder()
   builder.add_3d_conformer([
       [0.000, 0.000, 0.000],  # O
       [0.960, 0.000, 0.000],  # H
       [-0.240, 0.930, 0.000], # H
   ])
   mol = builder.build().with_assigned_valence()
   conformer_id = mol.conformers_3d()[0].id()

An application can instead supply a molecule loaded from a 3D molecular file
or produced by a conformer-generation API. Force-field calls below do not
depend on how those coordinates were obtained.

MMFF parameters, energy, and gradient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   if not mol.mmff_has_all_molecule_params():
       raise ValueError("MMFF parameters are unavailable for this molecule")

   properties = mol.mmff_properties_with_params(
       ck.MmffPropertiesParams(mmff_variant="MMFF94")
   )
   print(properties.atom_type(0))
   print(properties.formal_charge(0))
   print(properties.partial_charge(0))

   evaluation = mol.mmff_energy_gradient_with_params(
       ck.MmffEvaluationParams(
           mmff_variant="MMFF94",
           conformer_id=conformer_id,
           non_bonded_threshold=100.0,
           ignore_interfragment_interactions=True,
       )
   )
   if evaluation is None:
       raise ValueError("MMFF typing did not produce a valid force field")

   print(evaluation.energy())
   gradient = np.asarray(evaluation.gradient(), dtype=np.float64).reshape(-1, 3)
   print(gradient.shape)

The short forms ``mol.mmff_properties()`` and ``mol.mmff_energy_gradient()``
use default parameters. The current gradient is a flat Python list in
``[dx0, dy0, dz0, dx1, dy1, dz1, ...]`` order; the reshape above is performed
by the caller. Force is the negative energy gradient.

MMFF single- and multiple-conformer optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   result = mol.with_mmff_optimized_with_params(
       ck.MmffOptimizationParams(
           mmff_variant="MMFF94",
           max_iterations=200,
           conformer_id=conformer_id,
           non_bonded_threshold=100.0,
       )
   )
   optimized = result.molecule()
   print(result.status_code())
   print(result.needs_more())

   # The single-conformer MMFF result has no energy() accessor.
   final_evaluation = optimized.mmff_energy_gradient_with_params(
       ck.MmffEvaluationParams(conformer_id=conformer_id)
   )
   if final_evaluation is not None:
       print(final_evaluation.energy())

   # Optimize every stored 3D conformer independently.
   multiple = mol.with_mmff_optimized_confs_with_params(
       ck.MmffConformerOptimizationParams(
           num_threads=1,
           max_iterations=200,
           mmff_variant="MMFF94",
           non_bonded_threshold=100.0,
       )
   )
   for row in multiple.conformer_results():
       print(row.status_code(), row.needs_more(), row.energy())
   optimized_conformers = multiple.molecule()

These value-style calls leave ``mol`` unchanged. Status ``0`` indicates
convergence, ``1`` indicates that more iterations are needed, and MMFF status
``-1`` indicates invalid parameterization. Therefore ``not needs_more()``
alone is not a success check: inspect ``status_code()`` as well.

Default single-conformer MMFF optimization uses 200 iterations and a
non-bonded threshold of 100.0. Default multiple-conformer optimization uses
1000 iterations and a threshold of 10.0. The example explicitly chooses the
same threshold for both workflows.

UFF evaluation and optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   if not mol.uff_has_all_molecule_params():
       raise ValueError("UFF parameters are unavailable for this molecule")

   evaluation = mol.uff_energy_gradient_with_params(
       ck.UffEvaluationParams(conformer_id=conformer_id, vdw_threshold=10.0)
   )
   print(evaluation.energy())
   gradient = np.asarray(evaluation.gradient(), dtype=np.float64).reshape(-1, 3)

   result = mol.with_uff_optimized_with_params(
       ck.UffOptimizationParams(max_iterations=200, conformer_id=conformer_id)
   )
   optimized = result.molecule()
   print(result.status_code(), result.needs_more(), result.energy())

   multiple = mol.with_uff_optimized_confs_with_params(
       ck.UffConformerOptimizationParams(num_threads=1, max_iterations=200)
   )
   for row in multiple.conformer_results():
       print(row.conformer_id(), row.status_code(), row.energy())

The corresponding default calls are ``uff_energy_gradient()``,
``with_uff_optimized()``, and ``with_uff_optimized_confs()``. Default UFF
optimization uses 1000 iterations; its van der Waals threshold is 10.0.
Malformed inputs and evaluation failures propagate as Python exceptions;
non-convergence is an optimization result, not an exception.

Proposed Persistent API — Not Yet Implemented
---------------------------------------------

.. warning::

   Everything below is a target API design, not a currently callable API.
   ``MolecularForceField`` and its factory, parameter, result, and error types
   are proposed names. This page does not register them, implement them, or
   claim release or parity acceptance.

The goal is to build an owned force field once, then update its coordinates,
evaluate it, fix atoms, and perform short minimization runs without repeating
atom typing and contribution construction on every interaction.

Interactive Python workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # DESIGN EXAMPLE ONLY: these methods do not exist yet.
   ff = mol.mmff_force_field(
       conformer_id=conformer_id,
       mmff_variant="MMFF94",
       non_bonded_threshold=100.0,
       ignore_interfragment_interactions=True,
   )
   # UFF exposes the same handle interface:
   # ff = mol.uff_force_field(conformer_id=conformer_id, vdw_threshold=10.0)

   x, y, z = ff.position(1)
   ff.set_position(1, (x + 0.2, y, z))
   print(ff.energy())
   gradient = ff.gradient()  # independent float64 NumPy array, shape (N, 3)

   ff.set_fixed_atoms([0, 2])
   outcome = ff.minimize(max_iterations=20)
   print(outcome.converged, outcome.iterations, outcome.energy)

   # Another drag, followed by relaxation from the current coordinates.
   ff.set_position(1, (x + 0.4, y, z))
   outcome = ff.minimize(max_iterations=20)

   positions = ff.positions()  # independent writable NumPy snapshot
   positions[:, 2] += 0.1
   ff.set_positions(positions)
   ff.set_fixed_atoms([])      # replace the fixed set with an empty set

   evaluation = ff.energy_gradient()
   print(evaluation.energy, evaluation.gradient.shape)

Factory and configuration surface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Proposed Python factories accept keyword-only configuration:

.. code-block:: python

   mol.mmff_force_field(
       *, conformer_id=None, mmff_variant="MMFF94",
       non_bonded_threshold=100.0, ignore_interfragment_interactions=True,
   ) -> MolecularForceField

   mol.uff_force_field(
       *, conformer_id=None, vdw_threshold=10.0,
       ignore_interfragment_interactions=True,
   ) -> MolecularForceField

Reusable configuration remains available through immutable parameter objects:

.. code-block:: python

   # DESIGN EXAMPLE ONLY.
   params = ck.MmffForceFieldParams(
       conformer_id=conformer_id,
       mmff_variant="MMFF94",
   )
   ff = mol.mmff_force_field_with_params(params)
   outcome = ff.minimize_with_params(
       ck.ForceFieldMinimizeParams(
           max_iterations=20, force_tolerance=1e-4, energy_tolerance=1e-6
       )
   )

Rust retains default factories and explicit ``*_with_params`` methods with
typed parameters. Python keyword calls must construct those same parameter
values, with fields and defaults taken from the same binding contract. They
must not implement a second force-field path or a separate default table.
This keyword convenience projection must be explicitly settled in the API
rules and registry before implementation; it is not an existing overload.

Proposed handle and result API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   class MolecularForceField:
       def position(self, atom_id: int) -> tuple[float, float, float]: ...
       def set_position(self, atom_id: int, position) -> None: ...
       def positions(self) -> np.ndarray: ...
       def set_positions(self, positions) -> None: ...
       def fixed_atoms(self) -> tuple[int, ...]: ...
       def set_fixed_atoms(self, atom_ids) -> None: ...
       def energy(self) -> float: ...
       def gradient(self) -> np.ndarray: ...
       def energy_gradient(self) -> ForceFieldEnergyGradient: ...
       def minimize(self, *, max_iterations: int = 200,
                    force_tolerance: float = 1e-4,
                    energy_tolerance: float = 1e-6) -> ForceFieldMinimizeOutcome: ...
       def minimize_with_params(self, params) -> ForceFieldMinimizeOutcome: ...

   class ForceFieldEnergyGradient:
       energy: float           # read-only
       gradient: np.ndarray    # independent (N, 3) snapshot

   class ForceFieldMinimizeOutcome:
       converged: bool         # read-only
       iterations: int         # actual completed iterations
       energy: float           # final energy, read-only

``set_position`` accepts a three-component sequence. ``set_positions`` accepts
nested numeric sequences or NumPy arrays of shape ``(N, 3)``; input values are
copied into handle-owned storage. Python atom IDs are indices in the fixed
atom order, not IDs transferable between unrelated molecules. Rust uses
``AtomId`` and borrowed coordinate slices; Python never receives a mutable
view into the handle's live storage.

Both keyword minimization and immutable ``ForceFieldMinimizeParams`` expose
``energy_tolerance=1e-6``, alongside ``max_iterations=200`` and
``force_tolerance=1e-4``. Rust constructs these parameters with
``ForceFieldMinimizeParams::new(max_iterations, force_tolerance, energy_tolerance)``.
All three arguments reach the same source-backed optimizer. The pinned RDKit
``BFGSOpt.h`` explicitly ignores ``funcTol`` via ``RDUNUSED_PARAM(funcTol)``:
changing energy tolerance alone therefore does not currently change convergence.
CK preserves that behavior rather than inventing an energy stopping criterion.

Persistent factory and evaluator errors expose a typed
``MolecularForceFieldErrorKind`` through Rust ``kind()`` and Python ``kind``.
For example, compare ``error.kind == ck.MolecularForceFieldErrorKind.MissingConformer``,
not a string. Context remains available through ``requested``, ``atom_index``,
``expected``, ``actual`` and ``component`` where applicable; Rust ``source()``
and Python ``__cause__`` preserve the concrete lower-level cause.
``InvalidTolerance`` identifies the source optimizer's rejected force tolerance
(``component == 0``); energy tolerance is forwarded without a new CK-only constraint.

Ownership, updates, and numerical semantics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The handle owns its coordinates and evaluation state. It can outlive the
  input molecule and never implicitly writes coordinates back to it.
* Factories select an existing 3D conformer. With ``conformer_id=None``, the
  target follows the existing force-field default of selecting the first
  stored 3D conformer. No embedding or hydrogen addition is implicit.
* Topology, atom order, and parameterization are fixed at construction.
  Changing topology requires a new handle.
* Coordinate setters validate shape, atom indices, and finite values before
  committing. Invalid input leaves handle coordinates unchanged.
* Coordinate updates invalidate relevant distance caches automatically.
  Users do not call ``Initialize``. A single-atom update avoids copying all
  positions or retyping atoms; this does not make energy evaluation O(1).
* ``positions()`` and ``gradient()`` return independent NumPy snapshots.
  Editing a snapshot alone does not change the force field.
* ``set_fixed_atoms`` replaces the complete fixed set. Fixed atoms remain in
  energy terms but have zero components in the constrained gradient and do
  not move during minimization. Explicit setters may move them to new anchors.
* Non-bonded contributions are selected at construction using the configured
  threshold and initial geometry. Updating coordinates refreshes evaluation
  caches, not the contribution list. Rebuild the handle if interactions
  excluded at construction must be reconsidered after a large displacement.
* Repeated minimization resumes from current coordinates, but each call starts
  a new optimizer history. Two 20-iteration calls are not promised to equal
  one 40-iteration call. Non-convergence is reported by ``converged=False``.

Proposed construction errors are ``MmffForceFieldError`` and
``UffForceFieldError``; handle update, evaluation, and minimization errors use
``ForceFieldError`` with structured kinds and preserved causes. Invalid MMFF
parameterization must fail construction rather than return a dummy handle.

The existing one-shot APIs remain useful for ordinary molecule evaluation and
optimization. The persistent API adds an interactive evaluator; it does not
replace the chemistry algorithms or introduce a second force-field kernel.
