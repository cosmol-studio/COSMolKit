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
| Unsupported branches | `unsupported_stage` and `UnsupportedFeature` still appear in README-relevant modules. | Classify each as README-relevant or deferred; implement README-relevant branches or fail explicitly with support docs. |
| Operation policy | Topology and derived-state changes must go through registered operations and `OpParts`. | Keep code audits focused on mutation paths, not just public API existence. |

## Feature Status

### Chemistry Core

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| Atom / Bond / Molecule data model | `atom.rs`, `bond.rs`, `molecule.rs` | — | implemented | Continue guarding public mutable storage boundaries. |
| Adjacency representation | `adjacency.rs` | — | implemented | Keep cache invalidation tied to operation contracts. |
| Bond order + formal charge | `bond.rs`, `atom.rs`, `valence.rs` | `VALENCE_FEATURE` | substantial | Query/valence edge cases, radicals, dative branches. |
| SMILES parser | `smiles.rs`, `Molecule::from_smiles` | `SMILES_PARSE_FEATURE` (Experimental) | substantial | User-visible sanitize and unsupported stereo/query edge cases need audit; strict parity deferred. |
| Ring perception | `rings.rs` | `RINGS_FEATURE` | substantial | Utility surface and edge classes need functional closure. |
| Valence handling | `valence.rs`, `ops.rs` | `VALENCE_FEATURE` | substantial | Property cache, radicals, dative/query edge cases. |
| Kekulization | `kekulize.rs`, `ops.rs` | `KEKULIZE_FEATURE` | substantial | Edge failures and operation-state interactions need audit. |
| SMILES writer | `smiles_write.rs`, `Molecule::to_smiles` | `SMILES_WRITE_FEATURE` (Experimental) | substantial | Several `unsupported_stage` branches remain; classify and close README-relevant writer branches. |
| Atom/bond feature extraction | `Atom`, `Bond`, `query.rs` | — | partial | Public query-inspection surface is not complete. |
| Explicit hydrogen expansion/removal | `hydrogens.rs`, `ops.rs` | `HYDROGENS_FEATURE` | substantial | Coordinates, residue info, and less common AddHs/RemoveHs branches. |
| Tetrahedral stereo representation/perception | `stereo.rs`, `atom.rs`, `bond.rs`, `smiles.rs` | `STEREO_FEATURE` | partial | Typed state exists; CIP/geometric perception and some E/Z behavior remain. |
| DG bounds matrix | `distgeom.rs`, `Molecule::dg_bounds_matrix` | `DG_BOUNDS_FEATURE` | substantial | Implementation exists; support docs/tests need audit; strict distance parity deferred. |
| MolBlock V2000/V3000 handling | `io/sdf.rs`, `Molecule::from_mol_block` | `MOLBLOCK_IO_FEATURE` (Experimental) | substantial | V3000 reader, complex query predicates, supplier edge cases. |
| Sanitization pipeline | `ops.rs`, `sanitize.rs` | `SANITIZE_FEATURE` | partial | Full flag/error/cleanup surface not closed. |
| Morgan fingerprint | `fingerprint.rs` | `FINGERPRINT_FEATURE` (Experimental) | substantial | Functional options mostly present; strict RDKit hash parity deferred; audit unsupported options. |

### Chemical File I/O

| Checked Item | Entry Point | Feature | Status | Current Gaps |
|---|---|---|---|---|
| MOL reader | `Molecule::from_mol_block`, `io/sdf.rs` | `MOLBLOCK_IO_FEATURE` | partial | V2000 path exists; V3000 and edge parsing need closure. |
| SDF reader | `SdfDataset`, `SdfReader` | `MOLBLOCK_IO_FEATURE` | substantial | Robust multi-record behavior and supplier edge cases need audit. |
| SDF writer | `io/sdf.rs` | `MOLBLOCK_IO_FEATURE` | substantial | Complex query/stereo/SGroup branches need classification and closure. |
| SMILES output | `smiles_write.rs` | `SMILES_WRITE_FEATURE` | substantial | Remaining unsupported writer branches. |

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
| SVG drawing | `draw.rs`, `Molecule::to_svg` | `DRAWING_FEATURE` | substantial | Representative output tests and support-status consistency. |
| PNG rendering | `draw.rs`, `Molecule::to_png` | `DRAWING_FEATURE` | substantial | PNG byte-generation tests and batch export integration. |
| Embedded Noto Sans font | `draw.rs`, `assets/fonts/NotoSans-Regular.ttf` | `DRAWING_FEATURE` | implemented | Keep asset path and packaging assumptions audited. |

## Deferred For Later

- Python binding parity with README examples.
- Strict RDKit behavior parity, bit-identical Morgan hashes, exact DG bounds,
  and visual regression parity for drawings.
- Public query matching/inspection APIs unless the Rust README claim is
  explicitly narrowed to require them in this pass.
