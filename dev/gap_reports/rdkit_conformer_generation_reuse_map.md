# RDKit Conformer Generation Reuse Map

This is the Step 5 artifact for
`dev/archive/plans/rdkit_conformer_generation_full_port_plan.md`.

This audit compares the Step 3 RDKit conformer-generation source inventory
against current COSMolKit Rust code. It is an audit only; it does not expand
behavior or public API.

## Reuse Verdict

Existing COSMolKit code can be reused for:

- The local `BoundsMatrix` storage API and most RDKit `BoundsMatrix.h`
  behavior in `crates/cosmolkit-core/src/chemistry/distgeom.rs`.
- The DG bounds-builder path from `BoundsMatrixBuilder.cpp`, including
  topological 1-2, 1-3, 1-4, 1-5, VDW lower-bound, and wrapper matrix export
  behavior, subject to marker caveats below.
- Core `ForceField`, point-vector, distance-constraint, position-constraint,
  minimization, scatter/gather, and contribution dispatch machinery in
  `crates/cosmolkit-core/src/chemistry/forcefield/core.rs`.
- CrystalFF torsion preference parsing and torsion contribution machinery in
  `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/`.
- Existing UFF atom typing and parameter helpers reached by DG bounds
  generation.

Existing COSMolKit code cannot yet be reused as a complete conformer generator
because the RDKit `DistGeomUtils` and `Embedder` control-flow symbols are not
ported. In particular, there is no Rust source-backed implementation of
`EmbedMolecule`, `EmbedMultipleConfs`, `EmbedParameters`, ETKDG presets,
failure tracking, random distance matrix generation, metric embedding,
fourth-dimensional minimization, chiral checks, RMS pruning, or conformer
retry scheduling.

## DG Bounds Matrix And Triangle Smoothing

| RDKit symbol | Rust symbol | Status |
|---|---|---|
| `DistGeom::BoundsMatrix` | private `distgeom::BoundsMatrix` | Reusable internally. Storage is `Vec<Vec<f64>>`, not RDKit contiguous matrix storage. |
| `BoundsMatrix::BoundsMatrix(unsigned int)` | `BoundsMatrix::new` | Reusable internally; constructor initializes defaults. |
| `BoundsMatrix::getData()` | `BoundsMatrix::to_vec_vec` plus wrapper copy in `dg_bounds_matrix_with_options` | Partial reuse. Public surface returns nested vectors, not raw contiguous backing. |
| `BoundsMatrix::numRows()` | `BoundsMatrix::num_rows` | Reusable internally. |
| `getUpperBound`, `setUpperBound`, `setUpperBoundIfBetter` | `get_upper`, `set_upper`, `set_upper_if_better` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `getLowerBound`, `setLowerBound`, `setLowerBoundIfBetter` | `get_lower`, `set_lower`, `set_lower_if_better` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `BoundsMatrix::checkValid()` | `BoundsMatrix::check_valid` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `triangleSmoothBounds(BoundsMatPtr, double)` | `triangle_smooth_bounds_shared` | Reusable internally with `RDKit✔️✔️` forwarding behavior. |
| `triangleSmoothBounds(BoundsMatrix*, double)` | `triangle_smooth_bounds_ptr` | Behavior reusable, but marker is `RDKit✔️❌` because current storage is less cache-friendly than RDKit contiguous access. |

Remaining gap: if future conformer-generation code needs RDKit-equivalent
raw square-matrix performance, `BoundsMatrix` storage must be revisited or
the `RDKit✔️❌` performance marker retained.

## BoundsMatrixBuilder

| RDKit symbol/group | Rust symbol(s) | Status |
|---|---|---|
| `Path14Configuration` | private `Path14Configuration` | Reusable internally. |
| `ComputedData` and `visitedBound` | private `ComputedData`, `ComputedData::visited_bound` | Reusable internally. |
| `initBoundsMat` overloads | `init_bounds_mat_ptr`, `init_bounds_mat_shared` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `_checkAndSetBounds` | `BoundsMatrix::check_and_set_bounds_with_mode` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `set12Bounds` | `set_12_bounds` | Partial reuse. Source markers include `RDKit✔️❗`, so performance or modeled-state coverage is not fully settled. |
| `set13Bounds` and `_set13BoundsHelper` | `set_13_bounds`, `set_13_bounds_helper` | Partial reuse. Entry point has `RDKit✔️❗`; helper has stronger local anchors. |
| `set14Bounds` dispatch | `set_14_bounds` | Partial reuse. Entry point is `RDKit❗✔️`, so behavior is not yet fully proven equivalent. |
| `_record14Path` | `record_14_path` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_setInRing14Bounds` | `set_in_ring_14_bounds` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_setTwoInSameRing14Bounds` | `set_two_in_same_ring_14_bounds` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_setTwoInDiffRing14Bounds` | `set_two_in_diff_ring_14_bounds` | Reusable internally with `RDKit✔️✔️`. |
| `_setShareRingBond14Bounds` | `set_share_ring_bond_14_bounds` | Reusable internally with `RDKit✔️✔️`. |
| `_setMacrocycleTwoInSameRing14Bounds` | `set_macrocycle_two_in_same_ring_14_bounds` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_setMacrocycleAllInSameRing14Bounds` | `set_macrocycle_all_in_same_ring_14_bounds` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_setChain14Bounds` | `set_chain_14_bounds` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_checkH2NX3H1OX2` | `check_h2_nx3_h1_ox2` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_checkNhChChNh` | `check_nh_ch_ch_nh` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_checkAmideEster14` | `check_amide_ester_14` | Reusable internally with `RDKit✔️✔️`. |
| `_checkMacrocycleAllInSameRingAmideEster14` | `check_macrocycle_all_in_same_ring_amide_ester_14` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_isCarbonyl` | `is_carbonyl` | Reusable internally with `RDKit✔️✔️`. |
| `_checkAmideEster15` | `check_amide_ester_15` | Behavior reusable, performance marker currently `RDKit✔️❌`. |
| `_compute15Dists*` | `compute_15_dist_cis_cis`, `compute_15_dist_cis_trans`, `compute_15_dist_trans_trans`, `compute_15_dist_trans_cis` | Reusable internally with source-backed tests. |
| `set15Bounds` and `_set15BoundsHelper` | `set_15_bounds`, `set_15_bounds_helper` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `setLowerBoundVDW` overloads | `set_lower_bound_vdw` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `setTopolBounds` overloads | `set_topol_bounds`, `set_topol_bounds_with_outputs` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `collectBondsAndAngles` | `collect_bonds_and_angles` | Reusable internally with `RDKit✔️✔️` source anchors. |

Remaining gap: conformer generation can call this existing bounds-builder
code, but any step claiming full RDKit parity must either preserve existing
non-`✔️✔️` markers or upgrade them only after source-level review and tests.

## Public DG Bounds Surface

| RDKit wrapper symbol | Rust/Python surface | Status |
|---|---|---|
| `rdDistGeom::GetMoleculeBoundsMatrix` / `getMolBoundsMatrix` | `distgeom::dg_bounds_matrix_with_options`, `distgeom::dg_bounds_matrix`, Python `Molecule.dg_bounds_matrix` | Reusable for bounds matrix examples and tests. Marker is `RDKit✔️❌` because return storage differs from NumPy/raw memory behavior. |

Remaining gap: this is not conformer generation. It only computes bounds.

## DistGeom Chiral, Violation, Fourth-Dimension, And Utility Symbols

| RDKit symbol/group | Rust symbol | Status |
|---|---|---|
| `ChiralSetStructureFlags`, `ChiralSet`, `ChiralSetPtr`, `VECT_CHIRALSET` | None found | Missing. |
| `calcChiralVolume` overloads | None found | Missing. |
| `ChiralViolationContribsParams`, `ChiralViolationContribs` | None found | Missing. |
| `DistViolationContribsParams`, `DistViolationContribs` | None with RDKit DistGeom semantics | Missing. `DistanceConstraintContribs` is a ForceField helper but not the same source symbol. |
| `FourthDimContribsParams`, `FourthDimContribs` | None found | Missing. |
| `EIGVAL_TOL`, `KNOWN_DIST_TOL`, `KNOWN_DIST_FORCE_CONSTANT` | None found | Missing. |
| `pickRandomDistMat` overloads | None found | Missing. |
| `computeInitialCoords` overloads | None found | Missing. |
| `computeRandomCoords` overloads | None found | Missing. |
| `constructForceField` | None found | Missing. |
| `addImproperTorsionTerms` | None found in DistGeom context | Missing. |
| `addExperimentalTorsionTerms` | None found in DistGeom context | Missing. |
| `add12Terms`, `add13Terms`, `addLongRangeDistanceConstraints` | None found | Missing. |
| `construct3DForceField`, `constructPlain3DForceField`, `construct3DImproperForceField` | None found | Missing. |

`crates/cosmolkit-core/src/chemistry/forcefield/core.rs` contains reusable
generic force-field pieces:

- `ForceFieldVec3`
- `ForceField`
- `ForceFieldContrib`
- `DistanceConstraintContrib`
- `DistanceConstraintContribs`
- `PositionConstraintContrib`

These are reusable building blocks for future `DistGeomUtils` ports, but they
do not satisfy any missing DistGeom source symbol by themselves.

## CrystalFF Dependencies

| RDKit symbol/group | Rust symbol | Status |
|---|---|---|
| `CrystalFFDetails` | `CrystalFFDetails` | Reusable. |
| `ExpTorsionAngle` | private `ExpTorsionAngle` | Reusable internally. |
| `ExpTorsionAngleCollection` constructor | `ExpTorsionAngleCollection::new` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `ExpTorsionAngleCollection::getParams` | `ExpTorsionAngleCollection::get_params` | Reusable internally with `RDKit✔️✔️` source anchors. |
| `getExperimentalTorsions` overloads | `get_experimental_torsions`, `get_experimental_torsions_without_bonds` | Partial reuse. Main overload has `RDKit✔️❗`, so performance/coverage is unresolved. |
| `torsionPreferences_v1.in` | vendored source read by `torsion_preferences_v1` | Reusable. |
| `torsionPreferences_v2.in` | vendored source read by `torsion_preferences_v2` | Reusable. |
| `torsionPreferences_smallrings.in` | vendored source read by `torsion_preferences_small_rings` | Reusable. |
| `torsionPreferences_macrocycles.in` | vendored source read by `torsion_preferences_macrocycles` | Reusable. |
| `calcTorsionEnergyM6` | `calc_torsion_energy` and `TorsionAngleContribM6::get_energy` path | Needs Step 75/77 audit before declaring complete for conformer generation. |
| `TorsionAngleContribs` | `TorsionAngleContribs` | Reusable candidate, but Step 75 must audit all conformer-generation call requirements. |
| `TorsionAngleContribM6` | `TorsionAngleContribM6` | Reusable candidate, but Step 75 must audit overload parity. |

Remaining gap: CrystalFF code is substantially ported, but Step 75 remains
necessary because `DistGeomUtils` calls CrystalFF through conformer-generation
specific paths that have not been wired or audited end-to-end.

## UFF And MMFF Force-Field Code

Existing reusable UFF symbols include:

- `uff::atom_typer::get_atom_types_for_uff`
- `uff::atom_typer::get_atom_label_for_uff`
- `uff::params::AtomicParams`
- UFF bond, angle, torsion, inversion, nonbonded, and builder modules.

Existing reusable MMFF symbols include the MMFF parameter and contribution
modules under `crates/cosmolkit-core/src/chemistry/forcefield/mmff/`.

These modules support force-field optimization and DG bounds helper behavior,
but RDKit conformer generation's cleanup force field primarily depends on
`DistGeomUtils` plus CrystalFF torsion terms. UFF/MMFF optimization wrappers
are not a substitute for `EmbedMolecule` or `EmbedMultipleConfs`.

## Embedder Symbols

No existing Rust source-backed implementation was found for the RDKit embedder
source symbols:

- `EmbedFailureCauses`
- `EmbedParameters` and all fields
- `KDG`, `ETDG`, `ETDGv2`, `ETKDG`, `ETKDGv2`, `ETKDGv3`, `srETKDGv3`
- `DGeomHelpers::detail::EmbedArgs`
- `haveOppositeSign`
- failure-count mutex helpers
- `_volumeTest`
- `_sameSide`
- `_centerInVolume`
- `_boundsFulfilled`
- `generateInitialCoords`
- `firstMinimization`
- `checkTetrahedralCenters`
- `checkChiralCenters`
- `minimizeFourthDimension`
- `minimizeWithExpTorsions`
- `doubleBondGeometryChecks`
- `doubleBondStereoChecks`
- `finalChiralChecks`
- `embedPoints`
- `findDoubleBonds`
- `findChiralSets`
- `adjustBoundsMatFromCoordMap`
- `initETKDG`
- `setupInitialBoundsMatrix`
- `_fillAtomPositions`
- `_isConfFarFromRest`
- `multiplication_overflows_`
- `embedHelper_`
- `getMolSelfMatches`
- `EmbedMolecule`
- `EmbedMultipleConfs`

Existing `Conformer2D`, `Conformer3D`, and `ConformerStore` in the molecule
model can store generated coordinates, but they are storage primitives, not
embedder behavior.

## Python Wrapper Reuse

Current Python surface exposes:

- `Molecule.dg_bounds_matrix`
- UFF and MMFF parameter checks and optimization APIs for molecules that
  already have 3D conformers.
- `Molecule.num_conformers` and 3D coordinate access for stored conformers.

No Python source-backed conformer-generation wrapper exists yet. Future public
wrappers must expose the eventual Rust conformer-generation API directly and
must not rely on RDKit interop.

## Remaining Source Gaps By Plan Area

The next implementation steps must port, in order:

- `ChiralSet` and aliases.
- `calcChiralVolume`.
- `ChiralViolationContribs`.
- `DistViolationContribs`.
- `FourthDimContribs`.
- `pickRandomDistMat`.
- `computeInitialCoords`.
- `computeRandomCoords`.
- `constructForceField` and all 3D cleanup force-field builders.
- CrystalFF overload gaps identified by Step 75.
- Full `Embedder` parameter, preset, helper, retry, pruning, failure, and
  public wrapper behavior.

No existing COSMolKit symbol satisfies those missing source symbols completely
at this point.
