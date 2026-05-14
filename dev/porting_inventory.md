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
| Unsupported branches | SMILES writer `unsupported_stage` guards have been removed, but the frozen-scope audit in `dev/gap_reports/smiles_rdkit_remaining_source_scan.md` still shows unresolved reader/writer/rings/kekulize/valence/sanitize copied-source blocks. | Keep marker audits current and close remaining frozen-scope parity branches with source-backed tests. |
| SMILES dependencies | Reader-side postprocessing now runs through explicit sanitize/valence/kekulize/ring operations, and writer-side unsupported-stage guards are gone. | Keep `support.rs`, checklist evidence, and chemistry-core feature rows synchronized when those dependent operations change status. |
| Operation policy | Topology and derived-state changes must go through registered operations and `OpParts`. | Keep code audits focused on mutation paths, not just public API existence. |

## Feature Status

### Chemistry Core

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| Atom / Bond / Molecule data model | `atom.rs`, `bond.rs`, `molecule.rs` | — | implemented | Continue guarding public mutable storage boundaries. |
| Adjacency representation | `adjacency.rs` | — | implemented | Keep cache invalidation tied to operation contracts. |
| Bond order + formal charge | `bond.rs`, `atom.rs`, `valence.rs` | `VALENCE_FEATURE` | substantial | Query/valence edge cases, radicals, dative branches. |
| SMILES parser | `smiles.rs`, `Molecule::from_smiles` | `SMILES_PARSE_FEATURE` (Experimental) | substantial | Reader port is broad (conformer selection, non-tetrahedral assignment branches, atropisomer mutation, removeHs/query-completion subsets), and the parser's chemistry-core handoff now depends on explicit experimental sanitize/valence/kekulize/rings operations instead of placeholder branches. The frozen-scope audit is still open: `notation/smiles.rs` currently contains 1 `RDKit❌❌`, 2 `RDKit❗❗`, 14 `RDKit✔️❌`, and 2037 `RDKit❗✔️` copied-source lines across broader parse/CX/helper orchestration areas. |
| Ring perception | `rings.rs` | `RINGS_FEATURE` | substantial | SSSR active-bond filtering, D2 duplicate-candidate handling, D3/extra-ring discovery, symmetrized K4 storage, fastFindRings DFS traversal, and the URF-enabled ring-family/relevant-cycle path now have focused regression coverage. The frozen-scope audit is still open: `chemistry/rings.rs` currently contains 7 `RDKit❌❌` lines in the explicit non-URF branch plus 78 `RDKit❗✔️` lines in broader aggregate `findSSSR` / `symmetrizeSSSR` surfaces. |
| Valence handling | `valence.rs`, `ops.rs` | `VALENCE_FEATURE` | substantial | RDKit-like assignment is the active dependency for sanitize/kekulize/SMILES postprocessing. The frozen-scope audit is still open: `chemistry/valence.rs` currently contains 13 `RDKit❗✔️` lines, and related sanitize/property orchestration in `operations/ops.rs` still carries unresolved copied-source blocks. Remaining work is concentrated in property-cache maintenance, radicals, dative/query edge cases, and entrypoint/orchestration closure. |
| Kekulization | `kekulize.rs`, `ops.rs` | `KEKULIZE_FEATURE` | substantial | Fragment filtering, fused aromatic candidate selection, worker ordering/backtracking, dummy-question permutation, and value-style `KekulizeIfPossible` restoration now have focused parity tests. The frozen-scope audit is still open: `chemistry/kekulize.rs` currently contains 365 `RDKit❗✔️` lines in broader helper and entrypoint blocks. |
| SMILES writer | `smiles_write.rs`, `Molecule::to_smiles` | `SMILES_WRITE_FEATURE` (Experimental) | substantial | The checklist-closed double-bond-direction, non-tetrahedral writer, and helper/CX targeted tests are covered, and writer-side unsupported-stage guards are removed. The feature still depends on experimental aromatic/valence/kekulize state conventions, and the frozen-scope audit remains open because `notation/smiles_write.rs` currently contains 2 `RDKit❗❗` and 582 `RDKit❗✔️` copied-source lines across top-level writer orchestration, double-bond canonicalization, fragment handling, and broader CX/helper sections. |
| Atom/bond feature extraction | `Atom`, `Bond`, `query.rs` | — | partial | Public query-inspection surface is not complete. |
| Explicit hydrogen expansion/removal | `hydrogens.rs`, `ops.rs` | `HYDROGENS_FEATURE` | substantial | Coordinates, residue info, and less common AddHs/RemoveHs branches. |
| Tetrahedral stereo representation/perception | `stereo.rs`, `atom.rs`, `bond.rs`, `smiles.rs` | `STEREO_FEATURE` | substantial | CIP ranking (assignAtomCIPRanks), R/S label assignment, pseudo-3D wedge detection, non-tetrahedral stereo infrastructure, bond stereo codes, ring special cases all ported. Full assignAtomChiralTagsFromStructure (3D geometry → ChiralTag) and assignBondStereoCodes (E/Z from bond direction + StereoInfo) remain blocked on Conformer infrastructure. |
| DG bounds matrix | `distgeom.rs`, `Molecule::dg_bounds_matrix` | `DG_BOUNDS_FEATURE` | substantial | set15Bounds with stereochemistry detection ported. Triangle inequality smoother present. Implementation exists; strict distance parity deferred. |
| MolBlock V2000/V3000 handling | `io/sdf.rs`, `Molecule::from_mol_block` | `MOLBLOCK_IO_FEATURE` (Experimental) | substantial | V3000 reader, complex query predicates, supplier edge cases. |
| Sanitization pipeline | `ops.rs`, `sanitize.rs` | `SANITIZE_FEATURE` | partial | Registered weak-topology sanitize flow exists and is now the explicit reader-side chemistry handoff. The frozen-scope audit is still open in sanitize/property/cleanup orchestration: `operations/ops.rs` currently contains 3 `RDKit❗❗`, 180 `RDKit✔️❌`, and 33 `RDKit❗✔️` copied-source lines in the relevant orchestration surface. |
| Morgan fingerprint | `fingerprint.rs` | `FINGERPRINT_FEATURE` (Experimental) | substantial | Functional options mostly present; strict RDKit hash parity deferred; audit unsupported options. |

### Chemical File I/O

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| MOL reader | `Molecule::from_mol_block`, `io/sdf.rs` | `MOLBLOCK_IO_FEATURE` | partial | V2000 path exists; V3000 and edge parsing need closure. |
| SDF reader | `SdfDataset`, `SdfReader` | `MOLBLOCK_IO_FEATURE` | substantial | Robust multi-record behavior and supplier edge cases need audit. |
| SDF writer | `io/sdf.rs` | `MOLBLOCK_IO_FEATURE` | substantial | Complex query/stereo/SGroup branches need classification and closure. |
| SMILES output | `smiles_write.rs` | `SMILES_WRITE_FEATURE` | substantial | Unsupported-stage guards are removed and the targeted isomeric/directional/non-tetra/helper-CX output regressions are closed. The writer still depends on experimental aromatic/valence/kekulize state conventions, and the frozen-scope audit remains open because `notation/smiles_write.rs` still contains 2 `RDKit❗❗` and 582 `RDKit❗✔️` copied-source lines. |

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
| DG bounds matrix | `distgeom.rs` | `DG_BOUNDS_FEATURE` | substantial | Functional implementation exists; strict parity deferred. |
| Morgan fingerprint + Tanimoto | `fingerprint.rs` | `FINGERPRINT_FEATURE` | substantial | Hash parity deferred; unsupported option audit remains. |
| Query atom/bond storage | `query.rs`, `atom.rs`, `bond.rs`, `io/sdf.rs` | — | partial | Internal storage exists; public inspection/matching APIs deferred unless README-core requires them. |

### 2D Coordinates and Drawing

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| 2D coordinate generation | `coordinates.rs`, `ops.rs`, `Molecule::with_2d_coordinates` | `COORDINATE_2D_FEATURE` | substantial | Unsupported branch audit and functional test breadth. |
| SVG drawing | `draw.rs`, `Molecule::to_svg` | `DRAWING_FEATURE` | substantial | Annotation overlays (CIP codes, notes, SGroups, brackets, variable bonds, links, highlights, close-contacts), SVG metadata/data-tag/class attributes all ported. Representative output tests and support-status consistency. |
| PNG rendering | `draw.rs`, `Molecule::to_png` | `DRAWING_FEATURE` | substantial | PNG byte-generation tests and batch export integration. |
| Embedded Noto Sans font | `draw.rs`, `assets/fonts/NotoSans-Regular.ttf` | `DRAWING_FEATURE` | implemented | Keep asset path and packaging assumptions audited. |

## Deferred For Later

- Python binding parity with README examples.
- Strict RDKit behavior parity, bit-identical Morgan hashes, exact DG bounds,
  and visual regression parity for drawings.
- Public query matching/inspection APIs unless the Rust README claim is
  explicitly narrowed to require them in this pass.
