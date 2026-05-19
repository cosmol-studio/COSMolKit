# RDKit 2D Coordinate Remaining Source Scan

## Baseline

This audit compares the active Rust 2D coordinate generation surface against the
RDKit depiction baseline required by
`dev/coordinate_2d_rdkit_full_port_checklist.md`:

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
- `crates/cosmolkit-core/src/io/molblock.rs`
- `crates/cosmolkit-core/src/model/molecule.rs`
- `crates/cosmolkit-core/src/operations/ops.rs`
- `crates/cosmolkit-core/src/support.rs`

## Current Status Summary

- Rust already contains a substantial 2D coordinate implementation in
  `coordinates.rs`, including:
  - `embedRing`, `transformPoints`, `computeBisectPoint`, `reflectPoint`
  - `rankAtomsByRank`-named surface, but not full RDKit CIP-rank behavior
  - `computeInitialCoords`-like logic
  - non-tetrahedral helper functions
  - some collision-removal helpers
- The current public entrypoint is still a narrow
  `compute_2d_coords(atoms, bonds)` function, not the full RDKit parameterized
  API surface.
- The file still contains explicit TODO/FIXME notes showing missing source
  closure:
  - CIP-rank tiebreaking is not complete.
  - non-tetrahedral stereo helpers are present but not wired into the active
    placement pipeline.
  - full top-level `compute2DCoords()` control flow is not yet ported.
- No current Rust artifact was found for:
  - `Compute2DCoordParameters`
  - `preferCoordGen`
  - `CoordinateTemplates`
  - `compute2DCoordsMimicDistMat(...)`
  - constrained depiction APIs
  - `straightenDepiction(...)`
  - `normalizeDepiction(...)`
  - `Add2DCoordsToMol(...)`-equivalent wrapper

## Direct Function Inventory

### RDDepictor.h / RDDepictor.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| `preferCoordGen` | Missing | Add explicit runtime routing state and source-backed non-CoordGen behavior. |
| `Compute2DCoordParameters` | Missing | Add explicit Rust parameter type with RDKit defaults and forwarding semantics. |
| `compute2DCoords(mol, const Compute2DCoordParameters&)` | Missing | Add entrypoint equivalent instead of narrow default-only helper. |
| `compute2DCoords(mol, coordMap, canonOrient, clearConfs, nFlipsPerSample, nSamples, sampleSeed, permuteDeg4Nodes, forceRDKit, useRingTemplates)` | Missing | Add explicit overload-equivalent entrypoint and exact default control flow. |
| `compute2DCoordsMimicDistMat(...)` | Missing | No active Rust API or source-backed implementation found. |
| `ConstrainedDepictionParams` | Missing | Add explicit Rust parameter type and semantics. |
| `generateDepictionMatching2DStructure(...)` overload set | Missing | No active Rust constrained-depiction API found. |
| `generateDepictionMatching3DStructure(...)` | Missing | No active Rust API found. |
| `straightenDepiction(...)` | Missing | No active Rust API found. |
| `normalizeDepiction(...)` | Missing | No active Rust API found. |
| `DepictorLocal::getRankedAtomNeighbors(...)` | Present | Needs pipeline integration audit and targeted tests. |
| `DepictorLocal::embedSquarePlanar(...)` | Present | Defined but not wired into active top-level pipeline. |
| `DepictorLocal::embedTBP(...)` | Present | Defined but not wired into active top-level pipeline. |
| `DepictorLocal::embedOctahedral(...)` | Present | Defined but not wired into active top-level pipeline. |
| `DepictorLocal::embedNontetrahedralStereo(...)` | Present | Defined but not wired into active top-level pipeline. |
| `DepictorLocal::embedFusedSystems(...)` | Partial | RDKit-equivalent helper exists only in part; full pipeline integration still open. |
| `DepictorLocal::embedCisTransSystems(...)` | Partial | Seed detection exists in part, but full RDKit entry sequencing remains open. |
| `DepictorLocal::getNonEmbeddedAtoms(...)` | Missing/implicit | No exact source-backed helper artifact found. |
| `DepictorLocal::_findLargestFrag(...)` | Missing/implicit | No exact helper artifact found. |
| `DepictorLocal::_shiftCoords(...)` | Partial | Component shifting exists, but exact RDKit fragment-shift helper not closed. |
| `computeInitialCoords(...)` | Partial | Large body exists, but still lacks full RDKit pre-seeding and control flow closure. |
| `copyCoordinate(...)` | Missing | No exact Rust helper found. |
| top-level `compute2DCoords(...)` body | Partial | Current public helper does not reproduce full RDKit flow. |

### DepictUtils.h / DepictUtils.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| globals `BOND_LEN`, `COLLISION_THRES`, `BOND_THRES`, `ANGLE_OPEN`, `MAX_COLL_ITERS`, `HETEROATOM_COLL_SCALE`, `NUM_BONDS_FLIPS` | Partial | Some constants exist; full source-backed global surface and usage audit still needed. |
| `embedRing(...)` | Present | Needs exact source-anchor audit and targeted tests per checklist. |
| `transformPoints(...)` | Present | Needs exact source-anchor audit and targeted tests per checklist. |
| `computeBisectPoint(...)` | Present | Needs targeted tests per checklist. |
| `reflectPoint(...)` | Present | Needs targeted tests per checklist. |
| `reflectPoints(...)` | Present | Needs targeted tests and call-site audit. |
| `setNbrOrder(...)` | Partial | Current code still has explicit CIP-rank TODOs. |
| `findCoreRings(...)` | Present/partial | Requires exact helper-by-helper source audit and targeted tests. |
| `findNextRingToEmbed(...)` | Present/partial | Requires exact source audit and targeted tests. |
| `rankAtomsByRank(...)` | Partial | Not source-closed; CIP ranking fallback is incomplete. |
| `computeSubAngle(...)` | Present | Needs source-marker/perf review when touched. |
| `rotationDir(...)` | Present | Needs source-marker/perf review when touched. |
| `pickFirstRingToEmbed(...)` | Present/partial | Needs targeted tests and source-marker confirmation. |
| `getAllRotatableBonds(...)` | Missing | No exact active helper found. |
| `getRotatableBonds(...)` | Missing | No exact active helper found. |
| `findBondsPairsToPermuteDeg4(...)` | Missing/partial | No exact source-closed helper found. |
| `hasTerminalRGroupOrQueryHydrogen(...)` | Missing | No active helper found. |
| `prepareTemplateForRGroups(...)` | Missing | No active helper found. |
| `reducedToFullMatches(...)` | Missing | No active helper found. |
| `invertWedgingIfMolHasFlipped(...)` | Missing | No active helper found. |
| mimic-distance / random-sampling / collision helpers in this file | Partial | Some collision helpers exist, but full helper inventory and API closure remain open. |

### Templates.h / Templates.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| `CoordinateTemplates` singleton | Missing | No Rust registry found. |
| `CoordinateTemplates::assertValidTemplate(...)` | Missing | No Rust parser/validator found. |
| `CoordinateTemplates::loadTemplatesFromPath(...)` | Missing | No Rust loader found. |
| `CoordinateTemplates::setRingSystemTemplates(...)` | Missing | No Rust registry API found. |
| `CoordinateTemplates::addRingSystemTemplates(...)` | Missing | No Rust registry API found. |
| `loadDefaultTemplates()` | Missing | No Rust default-template loader found. |
| `RDDepict::setRingSystemTemplates(...)` wrapper | Missing | No Rust public wrapper found. |
| `RDDepict::addRingSystemTemplates(...)` wrapper | Missing | No Rust public wrapper found. |
| `RDDepict::loadDefaultRingSystemTemplates()` | Missing | No Rust public wrapper found. |

### EmbeddedFrag.h / EmbeddedFrag.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| `EmbeddedAtom` class state surface | Missing as exact type | Rust has `TreeEmbeddedAtom`, but not a source-closed `EmbeddedAtom`/`EmbeddedFrag` model. |
| `EmbeddedFrag(unsigned int aid, ...)` | Missing as exact artifact | Current Rust architecture does not expose the RDKit container/helper type. |
| `EmbeddedFrag(coordMap)` | Missing as exact artifact | Prespecified-coordinate fragment path not source-closed. |
| `EmbeddedFrag(fusedRings, useRingTemplates)` | Missing as exact artifact | Fused-ring embedding is not represented as a source-closed fragment type. |
| `EmbeddedFrag(dblBond)` | Missing as exact artifact | Cis/trans fragment constructor not source-closed. |
| `expandEfrag(...)` | Missing | No exact helper artifact found. |
| `addNonRingAtom(...)` | Missing | No exact helper artifact found. |
| `mergeNoCommon(...)` | Missing | No exact helper artifact found. |
| `mergeWithCommon(...)` | Missing | No exact helper artifact found. |
| `mergeFragsWithComm(...)` | Missing | No exact helper artifact found. |
| `setupNewNeighs()` / `updateNewNeighs(...)` / `setupAttachmentPoints()` | Missing/implicit | Logic is partially distributed in Rust helpers, not source-closed as RDKit helper bodies. |
| collision-density and local fragment transforms | Missing/partial | Some equivalent math exists, but fragment-level helper surface is not ported as direct artifacts. |

### Basement/Depictor.cpp

| RDKit function or type | Current Rust state | Remaining gap |
|---|---|---|
| `RDKit::Add2DCoordsToMol(ROMol &, bool useDLL)` | Missing | Need explicit Rust wrapper semantics; on non-CoordGen/non-Windows paths this must still be handled explicitly, not silently omitted. |

## Confirmed High-Priority Gaps Blocking Full Port

1. Parameter/API gap:
   `Compute2DCoordParameters` and both `compute2DCoords` entrypoints are not
   present as explicit Rust APIs.
2. Ranking gap:
   `rankAtomsByRank` is not source-closed because CIP-rank tie-breaking is
   still a TODO in active code.
3. Template gap:
   `CoordinateTemplates` and ring-template load/set/add/default APIs are absent.
4. Pipeline gap:
   non-tetrahedral stereo helpers exist but are not wired into the actual
   placement pipeline.
5. Top-level flow gap:
   current `compute_2d_coords(atoms, bonds)` does not yet implement the full
   RDKit `compute2DCoords(...)` control flow and parameter semantics.
6. Public-surface gap:
   no active Rust equivalents found for mimic-distance embedding, constrained
   depiction, normalization, straightening, or `Add2DCoordsToMol`.

## Execution Guidance For Next Steps

- The checklist start point is valid:
  - Step 5 should begin by introducing `Compute2DCoordParameters` and the two
    `compute2DCoords`-equivalent Rust entrypoints.
  - Step 7 should add targeted tests for default mapping and non-default
    forwarding.
  - Step 9 should run the specified focused test command.
- Later steps must not treat current helper presence as completion; several
  helpers exist only as partial or disconnected ports.
- No heuristic closure is acceptable for the remaining gaps above; each item
  must be completed from the corresponding RDKit source body.
