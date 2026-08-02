# RDKit 2D Coordinate Remaining Source Scan

## Baseline

This audit compares the active Rust 2D coordinate generation surface against the
RDKit depiction baseline required by
`dev/plans/coordinate_2d_rdkit_full_port_checklist.md`:

- `third_party/rdkit/Code/GraphMol/Depictor/RDDepictor.h`
- `third_party/rdkit/Code/GraphMol/Depictor/RDDepictor.cpp`
- `third_party/rdkit/Code/GraphMol/Depictor/DepictUtils.h`
- `third_party/rdkit/Code/GraphMol/Depictor/DepictUtils.cpp`
- `third_party/rdkit/Code/GraphMol/Depictor/EmbeddedFrag.h`
- `third_party/rdkit/Code/GraphMol/Depictor/EmbeddedFrag.cpp`
- `third_party/rdkit/Code/GraphMol/Depictor/Templates.h`
- `third_party/rdkit/Code/GraphMol/Depictor/Templates.cpp`
- `third_party/rdkit/Code/GraphMol/Depictor/Basement/Depictor.cpp`

Primary active Rust surface inspected:

- `crates/cosmolkit-core/src/chemistry/coordinates.rs`
- `crates/cosmolkit-core/src/model/molecule.rs`
- `crates/cosmolkit-core/src/operations/ops.rs`
- `crates/cosmolkit-core/src/properties/batch.rs`
- `crates/cosmolkit-core/src/io/molblock.rs`
- `crates/cosmolkit-core/src/properties/draw.rs`
- `crates/cosmolkit-core/src/support.rs`

## Current Status Summary

The previous early-plan gap report is no longer accurate. The active Rust
surface now includes all major RDKit 2D depiction entrypoints that were missing
when this file was first written:

- explicit `Compute2DCoordParameters`-equivalent parameter handling
- explicit `preferCoordGen` runtime flag handling
- `Add2DCoordsToMol`-equivalent molecule wrapper semantics
- ring-template parsing, runtime model, registry, default loading, set/add APIs
- `compute2DCoordsMimicDistMat(...)`
- constrained depiction overloads for 2D and 3D matching
- `straightenDepiction(...)`
- `normalizeDepiction(...)`
- value-style public exposure through registered `with_2d_coordinates`
  operation, batch helpers, MolBlock writer, and drawing fallbacks

The remaining work is no longer “missing API surface”. The remaining gap state
is concentrated in exact-source audit/marker closure, support-status wording,
and final whole-feature validation.

## Direct Function Inventory

### RDDepictor.h / RDDepictor.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| `preferCoordGen` | Present | Rust models the flag and explicit unsupported CoordGen runtime branch. Final audit should classify this as a non-Rust-runtime wrapper exclusion, not a missing chemistry helper. |
| `Compute2DCoordParameters` | Present | Active Rust parameter surface exists in `coordinates.rs`; only final marker/support audit remains. |
| `compute2DCoords(mol, const Compute2DCoordParameters&)` | Present | Implemented in active Rust call chain. Remaining work is final whole-feature validation and marker review. |
| `compute2DCoords(mol, coordMap, canonOrient, clearConfs, nFlipsPerSample, nSamples, sampleSeed, permuteDeg4Nodes, forceRDKit, useRingTemplates)` | Present | Implemented via explicit Rust options surface. |
| `compute2DCoordsMimicDistMat(...)` | Present | Implemented. Remaining gap is that unsupported/unmodeled branches must stay explicitly documented where Rust still rejects them. |
| `ConstrainedDepictionParams` | Present | Implemented. Some copied-source blocks still intentionally carry non-`✔️✔️` markers and need final audit review instead of new functionality. |
| `generateDepictionMatching2DStructure(...)` overload set | Present | Implemented. Remaining work is marker/support audit, not missing functionality. |
| `generateDepictionMatching3DStructure(...)` | Present | Implemented. Remaining work is final audit coverage statement. |
| `straightenDepiction(...)` | Present | Implemented with targeted tests already added. |
| `normalizeDepiction(...)` | Present | Implemented with targeted tests already added. |
| `DepictorLocal::getRankedAtomNeighbors(...)` | Present | Wired into active non-tetrahedral path. |
| `DepictorLocal::embedSquarePlanar(...)` | Present | Wired into active pipeline. |
| `DepictorLocal::embedTBP(...)` | Present | Wired into active pipeline. |
| `DepictorLocal::embedOctahedral(...)` | Present | Wired into active pipeline. |
| `DepictorLocal::embedNontetrahedralStereo(...)` | Present | Wired into active pipeline. |
| `DepictorLocal::embedFusedSystems(...)` | Present | Final marker/perf review only. |
| `DepictorLocal::embedCisTransSystems(...)` | Present | Final marker/perf review only. |
| `DepictorLocal::getNonEmbeddedAtoms(...)` | Present | Final audit should confirm copied-source closure where helper logic is distributed across Rust helpers. |
| `DepictorLocal::_findLargestFrag(...)` | Present | Final audit should confirm copied-source closure where helper logic is distributed across Rust helpers. |
| `DepictorLocal::_shiftCoords(...)` | Present | Implemented. |
| `computeInitialCoords(...)` | Present | Implemented in active path. Remaining work is final audit classification, not new porting. |
| `copyCoordinate(...)` | Present | Implemented in active path. |
| top-level `compute2DCoords(...)` body | Present | Implemented; public operation/caller exposure now also uses it. |

### DepictUtils.h / DepictUtils.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| globals `BOND_LEN`, `COLLISION_THRES`, `BOND_THRES`, `ANGLE_OPEN`, `MAX_COLL_ITERS`, `HETEROATOM_COLL_SCALE`, `NUM_BONDS_FLIPS` | Present | Final marker/perf review only. |
| `embedRing(...)` | Present | Final audit only. |
| `transformPoints(...)` | Present | Final audit only. |
| `computeBisectPoint(...)` | Present | Final audit only. |
| `reflectPoint(...)` | Present | Final audit only. |
| `reflectPoints(...)` | Present | Final audit only. |
| `setNbrOrder(...)` | Present | Final audit only. |
| `findCoreRings(...)` | Present | Final audit only. |
| `findNextRingToEmbed(...)` | Present | Final audit only. |
| `rankAtomsByRank(...)` | Present for the selected checklist baseline | Keep exact selected-baseline wording in final audit; do not overclaim beyond this checklist scope. |
| `computeSubAngle(...)` | Present | Final audit only. |
| `rotationDir(...)` | Present | Final audit only. |
| `pickFirstRingToEmbed(...)` | Present | Final audit only. |
| `getAllRotatableBonds(...)` | Present | Final audit only. |
| `getRotatableBonds(...)` | Present | Final audit only. |
| `findBondsPairsToPermuteDeg4(...)` | Present | Final audit only. |
| `hasTerminalRGroupOrQueryHydrogen(...)` | Present | Final audit only. |
| `prepareTemplateForRGroups(...)` | Present | Final audit only. |
| `reducedToFullMatches(...)` | Present | Final audit only. |
| `invertWedgingIfMolHasFlipped(...)` | Present | Final audit only. |
| mimic-distance / random-sampling / collision helpers in this file | Present | Final audit only. |

### Templates.h / Templates.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| `CoordinateTemplates` singleton | Present | Rust uses a process-wide registry protected by synchronization primitives; final audit should document this as the active architectural equivalent. |
| `CoordinateTemplates::assertValidTemplate(...)` | Present | Final audit only. |
| `CoordinateTemplates::loadTemplatesFromPath(...)` | Present | Final audit only. |
| `CoordinateTemplates::setRingSystemTemplates(...)` | Present | Final audit only. |
| `CoordinateTemplates::addRingSystemTemplates(...)` | Present | Final audit only. |
| `loadDefaultTemplates()` | Present | Final audit only. |
| `RDDepict::setRingSystemTemplates(...)` wrapper | Present inside active coordinates surface | If a broader public facade wrapper is desired later, that is an API-design choice, not a current chemistry-gap blocker for this checklist. |
| `RDDepict::addRingSystemTemplates(...)` wrapper | Present inside active coordinates surface | Same note as above. |
| `RDDepict::loadDefaultRingSystemTemplates()` | Present inside active coordinates surface | Same note as above. |

### EmbeddedFrag.h / EmbeddedFrag.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| `EmbeddedAtom` / `EmbeddedFrag` helper surface | Present as active Rust equivalents | The final audit must judge source-closure by behavior and copied-source anchors, not by exact C++ type names. |
| ring constructors / dblBond constructors / coord-map constructors | Present as active Rust helper paths | Final audit only. |
| `expandEfrag(...)` | Present | Final audit only. |
| `addNonRingAtom(...)` | Present | Final audit only. |
| `mergeNoCommon(...)` | Present | Final audit only. |
| `mergeWithCommon(...)` | Present | Final audit only. |
| `mergeFragsWithComm(...)` | Present | Final audit only. |
| `setupNewNeighs()` / `updateNewNeighs(...)` / `setupAttachmentPoints()` | Present in active Rust helper surface | Final audit should name the exact Rust helper inventory if needed. |
| collision-density and local fragment transforms | Present | Final audit only. |

### Basement/Depictor.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| `RDKit::Add2DCoordsToMol(ROMol &, bool useDLL)` | Present | The remaining gap is only explicit documentation that CoordGen-backed `useDLL` routing is unsupported in this runtime and fails explicitly instead of silently diverging. |

## Confirmed Remaining Gaps

1. Final audit gap:
   the old “missing function inventory” phase is over, but the selected baseline
   still needs a final source-level audit that states exactly which remaining
   non-`✔️✔️` markers are wrapper/runtime exclusions versus active chemistry gaps.
2. Support-status wording gap:
   `support.rs` and `dev/porting_inventory.md` must describe the new surface
   accurately without overclaiming blanket RDKit parity.
3. Final validation gap:
   full strict `cargo check` / `cargo test` validation and the frozen final
   audit step remain pending at this point in the checklist sequence.

## Execution Guidance For Next Steps

- Do not reopen earlier implementation steps based on this outdated audit model.
- Treat the remaining work as:
  - status synchronization
  - final whole-feature validation
  - final frozen audit wording
- Keep “no heuristic closure” strict: any additional code change after this
  point still requires a direct RDKit source basis.
