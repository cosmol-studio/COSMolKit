# Porting Inventory

Per-feature status of checked README items against the redesigned Rust core.
This file is a ledger, not primary evidence. Current Rust source is
authoritative when this file disagrees with code.

## Current Planning Mode

- Scope: checked README claims for `crates/cosmolkit-core/`.
- Out of scope for this pass: Python bindings, Python docs/stubs/tests, strict
  RDKit parity hardening, bit-identical hash/distance/drawing parity.
- Working estimate from current source: about 65-75% complete, use 70% as the
  planning baseline.
- Remaining work: functional closure, unsupported-branch cleanup, support-spec
  consistency, and operation-contract correctness.

## Status Legend

| Status | Meaning |
|--------|---------|
| `implemented` | Public Rust code path exists and is not a placeholder |
| `substantial` | Broad implementation exists, but README-level closure still has known gaps |
| `partial` | Some implementation exists, but important README-visible behavior remains |
| `deferred` | Explicitly outside this pass |

## Cross-Cutting Gaps

| Area | Current Evidence | Required Resolution |
|------|------------------|---------------------|
| Code/ledger drift | Current source is ahead of older plan notes for DG bounds, drawing, batch, Morgan, SDF, and SMILES writer. | Treat Rust code as primary evidence; keep this ledger synchronized after source inspection. |
| Python claims | README contains Python-facing claims and examples. | Deferred for this pass. Do not change Python bindings unless explicitly requested. |
| Strict parity | README often implies RDKit-compatible behavior. | Defer strict parity hardening. Close functional Rust behavior first; record parity gaps separately later. |
| Unsupported branches | SMILES writer `unsupported_stage` guards are gone, and the frozen writer/canon/rings files are now marker-closed. The frozen-scope audit in `dev/gap_reports/smiles_rdkit_remaining_source_scan.md` still shows unresolved parser, kekulize, valence, and sanitize-orchestration copied-source blocks. | Keep marker audits current and close the remaining frozen-scope parity branches with source-backed tests. |
| SMILES dependencies | Reader-side postprocessing now runs through explicit sanitize/valence/kekulize/ring operations, and writer-side unsupported-stage guards are gone. | Keep `support.rs`, checklist evidence, and chemistry-core feature rows synchronized when those dependent operations change status. |
| Operation policy | Topology and derived-state changes must go through registered operations and `OpParts`. | Keep code audits focused on mutation paths, not just public API existence. |

## Feature Status

### Chemistry Core

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| Atom / Bond / Molecule data model | `atom.rs`, `bond.rs`, `molecule.rs` | — | implemented | Continue guarding public mutable storage boundaries. |
| Adjacency representation | `adjacency.rs` | — | implemented | Keep cache invalidation tied to operation contracts. |
| Bond order + formal charge | `bond.rs`, `atom.rs`, `valence.rs` | `VALENCE_FEATURE` | substantial | Query/valence edge cases, radicals, dative branches. |
| SMILES parser | `smiles.rs`, `Molecule::from_smiles` | `SMILES_PARSE_FEATURE` (Experimental) | substantial | Reader port is broad (conformer selection, non-tetrahedral assignment branches, atropisomer mutation, removeHs/query-completion subsets), and the parser's chemistry-core handoff now depends on explicit experimental sanitize/valence/kekulize/rings operations instead of placeholder branches. The frozen-scope audit is still open: `notation/smiles.rs` currently contains 1 `RDKit❌❌`, 2 `RDKit❗❗`, 14 `RDKit✔️❌`, and 713 `RDKit❗✔️` copied-source lines across the remaining parser/helper blockers listed in `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`. |
| Ring perception | `rings.rs` | `RINGS_FEATURE` | substantial | SSSR active-bond filtering, D2 duplicate-candidate handling, D3/extra-ring discovery, symmetrized K4 storage, fastFindRings DFS traversal, and the URF-enabled ring-family/relevant-cycle path now have focused regression coverage. The frozen ring-perception file is marker-closed for the current checklist scope, but the feature remains experimental and should not be described as full RDKit closure outside that audited scope. |
| Valence handling | `valence.rs`, `ops.rs` | `VALENCE_FEATURE` | substantial | RDKit-like assignment is the active dependency for sanitize/kekulize/SMILES postprocessing. The frozen-scope audit is still open: `chemistry/valence.rs` now only retains 4 `RDKit✔️❌` lines in `ValenceContext::new`, and related sanitize/property orchestration in `operations/ops.rs` still carries unresolved copied-source blocks. Remaining work is concentrated in property-cache maintenance, radicals, dative/query edge cases, and entrypoint/orchestration closure. |
| Kekulization | `kekulize.rs`, `ops.rs` | `KEKULIZE_FEATURE` | substantial | Fragment filtering, fused aromatic candidate selection, worker ordering/backtracking, dummy-question permutation, and value-style `KekulizeIfPossible` restoration now have focused parity tests. The frozen-scope audit is still open because `chemistry/kekulize.rs` currently contains 397 `RDKit✔️❌` copied-source lines across the remaining kekulize helper and entrypoint blocks. |
| SMILES writer | `smiles_write.rs`, `Molecule::to_smiles` | `SMILES_WRITE_FEATURE` (Experimental) | substantial | The checklist-closed double-bond-direction, non-tetrahedral writer, and helper/CX targeted tests are covered, and writer-side unsupported-stage guards are removed. The frozen writer file is now marker-closed for the current checklist scope, but the feature still depends on experimental aromatic/valence/kekulize state conventions and should not be described as full end-to-end RDKit parity closure. |
| Atom/bond feature extraction | `Atom`, `Bond`, `query.rs` | — | partial | Public query-inspection surface is not complete. |
| Explicit hydrogen expansion/removal | `hydrogens.rs`, `ops.rs` | `HYDROGENS_FEATURE` | substantial | Coordinates, residue info, and less common AddHs/RemoveHs branches. |
| Tetrahedral stereo representation/perception | `stereo.rs`, `atom.rs`, `bond.rs`, `smiles.rs` | `STEREO_FEATURE` | substantial | CIP ranking (assignAtomCIPRanks), R/S label assignment, pseudo-3D wedge detection, non-tetrahedral stereo infrastructure, bond stereo codes, ring special cases all ported. Full assignAtomChiralTagsFromStructure (3D geometry → ChiralTag) and assignBondStereoCodes (E/Z from bond direction + StereoInfo) remain blocked on Conformer infrastructure. |
| DG bounds matrix | `distgeom.rs`, `Molecule::dg_bounds_matrix` | `DG_BOUNDS_FEATURE` | implemented | The checklist-scoped RDKit DG bounds baseline is behaviorally source-backed in Rust: raw `BoundsMatrix` upper/lower triangle storage, `triangleSmoothBounds`, 1-2/1-3/1-4/1-5 bound setting, `collectBondsAndAngles`, both `setTopolBounds` overloads, and `getMolBoundsMatrix(...)` wrapper defaults are implemented and covered by focused strict tests. Final audit: no known first-axis `RDKit❌*` gap remains in the audited call chain, but residual `RDKit✔️❌`, `RDKit✔️❗`, and `RDKit❗✔️` markers remain intentionally visible for performance and helper-abstraction caveats. This is not a blanket RDKit parity claim outside that audited baseline. |
| MolBlock V2000/V3000 handling | `io/sdf.rs`, `Molecule::from_mol_block` | `MOLBLOCK_IO_FEATURE` (Experimental) | substantial | V3000 reader, complex query predicates, supplier edge cases. |
| Sanitization pipeline | `ops.rs`, `sanitize.rs` | `SANITIZE_FEATURE` | partial | Registered weak-topology sanitize flow exists and is now the explicit reader-side chemistry handoff. The frozen-scope audit is still open in sanitize/property/cleanup orchestration: `operations/ops.rs` currently contains 216 `RDKit✔️❌` copied-source lines across the remaining orchestration and cleanup helpers listed in `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`. |
| Morgan fingerprint | `fingerprint.rs` | `FINGERPRINT_FEATURE` (Experimental) | substantial | Functional options mostly present; strict RDKit hash parity deferred; audit unsupported options. |
| Force fields | `forcefield/`, `uff_*`, `mmff_*` | — | substantial | The source-backed force-field core, UFF, MMFF, and CrystalFF modules have local branch coverage, and `rdkit_forcefield_params_parity.rs` currently passes with `0 ignored` tests for UFF/MMFF parameter coverage, initial energy/gradient, and single-/multi-conformer final-coordinate comparisons on the shared corpus. This is not a blanket claim that every remaining copied RDKit marker in force-field files is closed; marker cleanup and any unsupported source branches must remain explicit. |

### Chemical File I/O

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| MOL reader | `Molecule::from_mol_block`, `io/sdf.rs` | `MOLBLOCK_IO_FEATURE` | partial | V2000 path exists; V3000 and edge parsing need closure. |
| SDF reader | `SdfDataset`, `SdfReader` | `MOLBLOCK_IO_FEATURE` | substantial | Robust multi-record behavior and supplier edge cases need audit. |
| SDF writer | `io/sdf.rs` | `MOLBLOCK_IO_FEATURE` | substantial | Complex query/stereo/SGroup branches need classification and closure. |
| SMILES output | `smiles_write.rs` | `SMILES_WRITE_FEATURE` | substantial | Unsupported-stage guards are removed and the targeted isomeric/directional/non-tetra/helper-CX output regressions are closed. The frozen writer file is marker-closed for the current checklist scope, but the writer still depends on experimental aromatic/valence/kekulize state conventions and is not yet a blanket RDKit-parity claim. |

### Batch Processing

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| Ordered construction from SMILES | `batch.rs` | `BATCH_FEATURE` | implemented | Depends on parser behavior. |
| Batch transforms | `batch.rs` | `BATCH_FEATURE` | substantial | Existing wrappers need audit across all operation errors. |
| Parallel scheduling / `n_jobs` | `batch.rs` | `BATCH_FEATURE` | partial | Full scheduling semantics not closed. |
| Structured errors | `batch.rs` | `BATCH_FEATURE` | substantial | Error taxonomy/reporting needs completion. |
| Validity masks and filtering | `batch.rs` | `BATCH_FEATURE` | substantial | Reports and richer summaries need closure. |
| Selection / indexing | `batch.rs` | `BATCH_FEATURE` | partial | README-level Rust selection API breadth needs decision and implementation. |
| SDF/image export | `batch.rs`, `io/sdf.rs`, `draw.rs` | `BATCH_FEATURE` | partial | Image export exists; SDF export and parallel export need closure. |

### Query, Descriptors, Computation

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| DG bounds matrix | `distgeom.rs` | `DG_BOUNDS_FEATURE` | implemented | The selected RDKit DG bounds call graph is behaviorally source-backed in Rust, including raw matrix semantics, smoothing, helper dispatch, and wrapper defaults. Local strict `distgeom` coverage passes for the audited scope. Final audit found no remaining first-axis `RDKit❌*` gap in the active DG bounds call chain, while performance and helper-abstraction markers remain documented. Keep parity claims separate: this ledger row records audited baseline completion, not blanket parity for every future RDKit-comparison fixture. |
| Morgan fingerprint + Tanimoto | `fingerprint.rs` | `FINGERPRINT_FEATURE` | substantial | Hash parity deferred; unsupported option audit remains. |
| Query atom/bond storage | `query.rs`, `atom.rs`, `bond.rs`, `io/sdf.rs` | — | partial | Internal storage exists; public inspection/matching APIs deferred unless README-core requires them. |

### 2D Coordinates and Drawing

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| 2D coordinate generation | `coordinates.rs`, `ops.rs`, `Molecule::with_2d_coordinates` | `COORDINATE_2D_FEATURE` | substantial | The selected RDKit 2D depiction baseline is now broadly ported in active Rust: parameterized `compute2DCoords` entrypoints, `preferCoordGen`/`forceRDKit` routing, ring-template registry/loaders, mimic-distance embedding, constrained 2D/3D depiction matching, normalization/straightening, and value-style exposure through registered operations, batch helpers, MolBlock writing, and drawing fallbacks are all present. Remaining work is final source-audit closure, support-status wording, and whole-feature strict validation. |
| SVG drawing | `draw.rs`, `Molecule::to_svg` | `DRAWING_FEATURE` | substantial | Annotation overlays (CIP codes, notes, SGroups, brackets, variable bonds, links, highlights, close-contacts), SVG metadata/data-tag/class attributes all ported. Representative output tests and support-status consistency. |
| PNG rendering | `draw.rs`, `Molecule::to_png` | `DRAWING_FEATURE` | substantial | PNG byte-generation tests and batch export integration. |
| Embedded Noto Sans font | `draw.rs`, `assets/fonts/NotoSans-Regular.ttf` | `DRAWING_FEATURE` | implemented | Keep asset path and packaging assumptions audited. |

## Deferred For Later

- Python binding parity with README examples.
- Strict RDKit behavior parity, bit-identical Morgan hashes,
  and visual regression parity for drawings.
- Public query matching/inspection APIs unless the Rust README claim is
  explicitly narrowed to require them in this pass.
