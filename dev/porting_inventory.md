# Porting Inventory

Per-feature status of checked README items against the redesigned Rust core.
This file is a ledger, not primary evidence. Current Rust source is
authoritative when this file disagrees with code.

## Current Planning Mode

- Scope: checked README claims for `crates/cosmolkit-core/`.
- This ledger does not redefine parity evidence. Current exact comparison
  boundaries and corpus results are authoritative in
  [`parity_scope.md`](./parity_scope.md).
- No repository-wide completion percentage is maintained. The feature rows and
  source-backed gap reports are the status evidence.
- Remaining ledger work: functional closure outside the documented parity
  surfaces, capability-boundary documentation, support-spec consistency, and
  operation-contract correctness.
- Executable remaining queues are listed in [`plans/`](./plans/).

## Status Legend

| Status | Meaning |
|--------|---------|
| `implemented` | Public Rust code path exists and is not a placeholder |
| `substantial` | Broad implementation exists, but README-level closure still has known gaps |
| `partial` | Some implementation exists, but important README-visible behavior remains |
| `deferred` | Explicitly outside this pass |
| `unsupported` | A separately named public capability is outside implemented scope and returns no result; this is not a status for failing rows inside a supported boundary |

## Cross-Cutting Gaps

| Area | Current Evidence | Required Resolution |
|------|------------------|---------------------|
| Code/ledger drift | Current source is ahead of older plan notes for DG bounds, drawing, batch, Morgan, SDF, and SMILES writer. | Treat Rust code as primary evidence; keep this ledger synchronized after source inspection. |
| Python claims | README contains Python-facing claims and examples. | Deferred for this pass. Do not change Python bindings unless explicitly requested. |
| Strict parity | Exact covered boundaries and their ChEMBL 37 or maintained-corpus evidence are recorded in `parity_scope.md`. | Keep this implementation ledger synchronized without weakening or broadening those claims. |
| Source-audit scope | SMILES writer `unsupported_stage` guards are gone, and the frozen writer/canon/rings files are marker-closed. The source audit in `dev/gap_reports/smiles_rdkit_remaining_source_scan.md` still tracks parser, kekulize, valence, and sanitize-orchestration copied-source blocks outside its completed marker set; these are audit-closure items, not accepted mismatches in the passing profiles documented by `parity_scope.md`. | Keep marker audits current and close the remaining source-review items without weakening the documented passing boundaries. |
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
| InChI scalar generation, keying, and parsing | `inchi.rs`, `mol_to_inchi`, `mol_to_inchi_key`, `inchi_to_inchi_key`, `mol_from_inchi`, Python `Chem` | `INCHI_FEATURE` (SupportedWithRdkitParity for source-defined behavior) | implemented | The configured five-root official-engine closure and RDKit adapter boundary are source-ported for exactly four scalar APIs. Official InChI v1.07.5 and RDKit 2026.03.1 tests compare source-defined behavior exactly. The official-C undefined `NormalizeAndCompare` initial-allocation path intentionally returns structured `allocation_failed` in Rust. MolBlock, SDF/V3000, IXA, AuxInfo, INCHIGEN, version query, and extended-polymer APIs are separately scoped upstream surfaces outside this four-API entry. |
| Atom/bond feature extraction | `Atom`, `Bond`, `query.rs` | — | partial | Public query-inspection surface is not complete. |
| SMARTS parser, writer, and ordinary-molecule substructure matching | `search/smarts_parse.rs`, `search/smarts_write.rs`, `search/substruct.rs` | `substructure.match` | implemented | The canonical query-bearing `Molecule` path is checked by 162 pinned RDKit 2026.03.1 source rows: 91 accepted parser inputs, nine rejects, and 62 ordered matcher cases. Exact parser acceptance, atom/bond counts, molecule SMARTS output, and ordered atom mappings are golden-compared. Separate strict suites cover fragment/CXSMARTS composition, parameters, callbacks, Python conversion, architecture, and existing consumers. Reaction SMARTS, database/container SMARTS, and the Bison-specific `debug_parse=true` diagnostic stream remain explicitly excluded; debug mode returns a structured unsupported error. |
| Explicit hydrogen expansion/removal | `hydrogens.rs`, `ops.rs` | `HYDROGENS_FEATURE` | substantial | Coordinates, residue info, and less common AddHs/RemoveHs branches. |
| Stereo representation and atom-tag perception | `stereo.rs`, `atom.rs`, `bond.rs`, `smiles.rs`, `Molecule::with_chiral_tags_from_structure` | `STEREO_FEATURE` | substantial | The standalone `assignAtomChiralTagsFromStructure` / `assignChiralTypesFrom3D` scope is source-ported and exact across 77 pinned-RDKit full-state oracle records, including conformer selection, replacement, tetrahedral and enabled non-tetrahedral branches, properties, no-ops, and errors. Typed stereo inspection, CIP ranking/codes, pseudo-3D wedge detection, non-tetrahedral tables, and ring cases are present. The broader `assignStereochemistryFrom3D` workflow and 3D double-bond direction/E-Z assignment remain separate, unclaimed capabilities rather than gaps hidden inside the completed atom-tag operation. |
| Modern CIPLabeler | `ciplabeler.rs`, `Molecule::with_cip_labels*`, Python assignment/query methods | `CIP_LABELER_FEATURE` | implemented | The single source-backed modern core covers full and selected assignment, recursion limits, typed R/S/r/s, E/Z, and M/P/m/p state, lifecycle, Rust/Python APIs, focused fixtures, the maintained 5,000-row gate, and 11,417,448 exact ChEMBL 37 full/selected state comparisons with zero mismatch. Non-tetrahedral configurations outside pinned `findConfigs` and source-equivalent process cancellation are distinct unsupported capabilities, not gaps inside the claimed dispatcher boundary. |
| DG bounds matrix | `distgeom.rs`, `Molecule::dg_bounds_matrix` | `DG_BOUNDS_FEATURE` | implemented | The checklist-scoped RDKit DG bounds baseline is behaviorally source-backed in Rust: raw `BoundsMatrix` upper/lower triangle storage, `triangleSmoothBounds`, 1-2/1-3/1-4/1-5 bound setting, `collectBondsAndAngles`, both `setTopolBounds` overloads, and `getMolBoundsMatrix(...)` wrapper defaults are implemented and covered by focused strict tests. Final audit: no known first-axis `RDKit❌*` gap remains in the audited call chain, but residual `RDKit✔️❌`, `RDKit✔️❗`, and `RDKit❗✔️` markers remain intentionally visible for performance and helper-abstraction caveats. This is not a blanket RDKit parity claim outside that audited baseline. |
| 3D conformer generation | `distgeom.rs`, `Molecule::with_3d_conformer*`, `Molecule::with_3d_conformers*`, Python wrapper | `CONFORMER_GENERATION_FEATURE` | substantial | The active conformer-generation surface is source-backed across the exposed DG/KDG/ETDG/ETKDG presets, Rust and Python value-style entry points, failure tracking, deterministic explicit-seed single-conformer path, deterministic batch seed policy, source-backed unseeded `clock()` RNG path, pruning, coordMap, CPCI, custom bounds-matrix validation, stereo/chiral checks, macrocycle and small-ring torsion paths, examples, and strict parity coverage. Final marker audit is not closed: `dev/gap_reports/rdkit_conformer_generation_final_marker_audit.md` records one remaining first-axis `RDKit❌❌` block for `MolAlign::details::symmetrizeTerminalAtoms()` during symmetry pruning, plus unresolved `RDKit❗✔️`, `RDKit✔️❗`, and `RDKit✔️❌` helper-surface markers in bounds building and CrystalFF torsion preference collection. |
| MolBlock V2000/V3000 handling | `io/sdf.rs`, `Molecule::from_mol_block` | `MOLBLOCK_IO_FEATURE` (Experimental) | substantial | V3000 reader, complex query predicates, supplier edge cases. |
| Sanitization pipeline | `ops.rs`, `sanitize.rs` | `SANITIZE_FEATURE` | partial | Registered weak-topology sanitize flow exists and is now the explicit reader-side chemistry handoff. The frozen-scope audit is still open in sanitize/property/cleanup orchestration: `operations/ops.rs` currently contains 216 `RDKit✔️❌` copied-source lines across the remaining orchestration and cleanup helpers listed in `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`. |
| Morgan fingerprint | `fingerprint.rs` | `FINGERPRINT_FEATURE` (Experimental) | substantial | The documented Morgan profiles and `AdditionalOutput` branches have exact ChEMBL 37 evidence; unrelated upstream fingerprint APIs remain outside this row. |
| Topological Torsion fingerprint | `fingerprint.rs`, `MoleculeBatch::topological_torsion_*` | `TOPOLOGICAL_TORSION_FINGERPRINT_FEATURE` | implemented | The modern four-vector generator, options, selections, chirality preparation, shared `AdditionalOutput`, JSON, scalar/bulk APIs, Rust ordered batch conveniences, three legacy adapters, and Python atom-code/torsion helpers share one source-backed chemistry/vector core. All 152 focused rows pass 12 profiles, all 5,000 maintained-corpus rows pass nine profiles, and all 2,897,804 mutually parseable ChEMBL 37 molecules match across 36 vector and eight complete provenance outputs (127,503,376 exact comparisons) against pinned RDKit 2026.03.1. Atom Pair fingerprints, RDKFingerprint, and unrelated fingerprint families remain outside this distinct entry. |
| RDKFingerprint/topological fingerprint | `fingerprint.rs`, `Molecule::topological_fingerprint` | `FINGERPRINT_FEATURE` | substantial | Source-exact RDKitFP generator and focused exact-bit tests are implemented. All 2,897,804 mutually parseable ChEMBL 37 molecules match across 14 vector and two complete provenance profiles (46,364,864 exact comparisons); broader unmodeled RDKit fingerprint families remain outside this entry. |
| Avalon fingerprint | `avalon_fingerprint.rs`, `Molecule::avalon_fingerprint` | `AVALON_FINGERPRINT_FEATURE` | implemented | Source-exact Avalon/REACCS explicit-bit path and byte-rounded vector semantics match across all 2,897,804 mutually parseable ChEMBL 37 molecules and 23 profiles (66,649,492 exact comparisons). Count/string overloads and unrelated Avalon APIs remain out of scope. |
| Force fields | `forcefield/`, `uff_*`, `mmff_*` | — | substantial | The source-backed force-field core, UFF, MMFF, and CrystalFF modules have local branch coverage, and `rdkit_forcefield_params_parity.rs` currently passes with `0 ignored` tests for UFF/MMFF parameter coverage, initial energy/gradient, and single-/multi-conformer final-coordinate comparisons on the shared corpus. This is not a blanket claim for unrelated upstream force-field APIs; separately scoped capabilities and remaining source markers stay visible without weakening the passing documented boundary. |

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
| 3D conformer generation | `distgeom.rs`, `support.rs`, Rust facade, Python wrapper | `CONFORMER_GENERATION_FEATURE` | substantial | The strict validation slice now passes for the exposed conformer-generation surface, including Rust facade exports, Python wrapper exports, examples, Sphinx docs, whole-workspace strict tests, and source-ported terminal-group symmetrization during symmetry pruning. The final marker audit no longer has any first-axis `RDKit❌❌` block in the audited conformer-generation path, but multiple second-axis/perf-review markers remain in bounds-builder helpers and CrystalFF torsion preference collection, so this row remains substantial instead of implemented. |
| Morgan fingerprint + Tanimoto | `fingerprint.rs` | `FINGERPRINT_FEATURE` | substantial | The enumerated Morgan vector and `AdditionalOutput` profiles have exact ChEMBL 37 evidence. Tanimoto behavior and unrelated upstream overloads retain their separately documented boundaries. |
| Topological Torsion fingerprint | `fingerprint.rs`, `batch.rs` | `TOPOLOGICAL_TORSION_FINGERPRINT_FEATURE` | implemented | Exact focused, 5,000-row maintained-corpus, and complete ChEMBL 37 tests cover scalar and shared-generator outputs; the full-corpus audit performs 127,503,376 exact vector/provenance comparisons with zero mismatch. Focused strict batch tests require ordered scalar equality, indexed errors, 1/2/4-thread determinism, provenance equality, and source immutability. |
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
- Upstream operation families and overloads not enumerated by the current
  feature-specific parity boundaries.
- Public query matching/inspection APIs unless the Rust README claim is
  explicitly narrowed to require them in this pass.
