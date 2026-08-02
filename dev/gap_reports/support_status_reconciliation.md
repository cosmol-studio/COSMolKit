# Public Support Status Reconciliation

This is a working audit for reconciling the root README, public `FeatureSpec`
values, operation-registry metadata, and executable parity evidence. Delete this
report after every row has been resolved into the normative status model and
documentation.

## Decision Rule

A public feature may be stable within a documented scope when:

- the public entry point is implemented and its unsupported boundary is
  explicit;
- focused tests cover its source-defined active branches;
- user-visible behavior has direct or end-to-end parity evidence where an
  external reference defines that behavior.

Complete first- and second-axis source-marker closure is not required for a
stable public scope. Source markers remain port-review evidence and must not be
silently upgraded from a product support decision.

## Evidence Classes

- **Direct parity**: the test invokes the public target operation and compares
  its observable output with pinned external expected data or an independent
  oracle.
- **End-to-end parity**: the target behavior is exercised as part of a public
  workflow and its resulting graph, coordinates, or serialization is compared.
- **Source-focused**: focused tests exercise source-defined branches against
  source-derived expectations, but do not independently invoke an external
  implementation.

The operation registry currently names `parity_profile` values without a
machine-enforced mapping from each profile to its test target. This is a status
infrastructure gap, not evidence that the underlying focused tests do not
exist.

## Status Decisions

| Feature | Current evidence | Decision |
|---|---|---|
| SMILES parse | Graph-field golden parity plus parser error/branch focused tests | Stable with pinned RDKit parity for the documented parse scope |
| SMILES write | Exact common branch matrix and explicit exhaustive matrix | Stable with pinned RDKit parity for supported branches |
| MolBlock/SDF IO | Reader field/state parity and writer output parity | Stable with pinned RDKit parity for supported V2000/V3000 branches |
| MOL2 read | Direct topology, atom, bond, coordinate, chirality, and SMILES parity | Stable with pinned RDKit parity for the exposed Tripos profile |
| Hydrogen operations | End-to-end graph parity plus extensive AddHs/RemoveHs source-focused tests | Stable with pinned RDKit parity for modeled parameters |
| 2D coordinates | Prepared-coordinate/final-depiction parity plus direct source-focused tests | Stable with pinned RDKit parity for the local RDKit coordinate path |
| Coordinate editing | Rust-native operation-contract and invariant tests | Stable COSMolKit API |
| 3D conformer generation | Seeded direct coordinate and parameter parity | Stable with pinned RDKit parity for exposed presets and controls |
| Sanitize | End-to-end graph parity plus extensive per-stage source-focused tests | Stable with pinned RDKit parity for modeled `SanitizeOps` branches |
| Kekulize | Direct flag golden parity plus focused source tests | Stable with pinned RDKit parity for modeled branches |
| Valence/radicals | Graph-field parity plus focused operation/source tests | Stable with pinned RDKit parity for modeled branches |
| Rings/aromaticity | Graph-field parity, focused ring/aromaticity tests, and RDL-backed ring-family tests | Stable with pinned RDKit parity for modeled graph states |
| Morgan/MACCS | Exact bit/count/additional-output parity for enumerated branches | Stable with pinned RDKit parity; split the combined feature descriptor later |
| SVG/PNG drawing | Exact normalized SVG parity; PNG is local deterministic rasterization | Stable; split SVG parity from PNG support later |
| Batch | Order, value semantics, errors, and parallel determinism | Stable COSMolKit API |
| DG bounds | Full matrix parity over the active profile | Stable with pinned RDKit parity |
| BioStructure and structural readers | Extensive source-focused records and public integration fixtures | Stable for documented PDB/mmCIF/mmJSON structural scope |
| InChI scalar APIs | Official C and pinned RDKit exact field/branch parity | Already stable with parity |
| Molecular descriptors | Exact value/bit parity for all exposed options plus matching parse-failure outcomes | Stable with pinned RDKit parity for the exposed descriptor surface |

## Rows Requiring More Work

### Stereo

Split stable typed stereo inspection and covered perception from broader
coordinate-derived assignment branches. Do not use one broad status to hide an
unsupported sub-capability.

### Substructure and SMARTS

Separate molecule-query matching, SMARTS parsing, and SMARTS-query matching.
Current molecule-query parity is small and does not establish complete SMARTS
query semantics.

### PDB/mmCIF Molecule Conversion

PDB conversion has an independent RDKit comparison over the 7SH6 fixture and
focused branch tests. mmCIF conversion uses the Gemmi-primary `BioStructure`
reader followed by the RDKit-compatible molecule conversion profile and has
focused tests. Add a dedicated structured feature descriptor so this public
surface is not incorrectly represented by the BioStructure-reader features.

## Infrastructure Follow-Up

- Add a machine-readable mapping from every parity-sensitive public feature or
  operation profile to its concrete test target and expected-data manifest.
- Generate or check public support tables from `PUBLIC_FEATURES` instead of
  maintaining independent status claims.
- Move historical `frozen` wording to development reports; public feature docs
  should describe only supported and unsupported behavior boundaries.

## Applied Internal Status Changes

The evidence-backed feature descriptors for SMILES parse/write, MolBlock/SDF
IO, MOL2 read, hydrogen operations, 2D coordinates, conformer generation,
sanitize, kekulize, valence/radicals, rings, aromaticity, fingerprints, DG
bounds, drawing, batch, and BioStructure/readers now describe their stable
documented scopes. Registered parity-sensitive molecule operations now use
`RequiredNow` for the covered hydrogen, sanitize, kekulize, valence, radical,
ring, aromaticity, 2D-coordinate, and 3D-conformer surfaces.

`with_2d_coordinates` was incorrectly attached to the generic coordinate-edit
feature. It now references the dedicated `coordinates.2d` feature descriptor.

Stereo, substructure/SMARTS, and bio selection remain experimental for the
specific unresolved boundaries listed above. Descriptor parity now covers
`only_heavy=true`, all four exposed Crippen option combinations, all four
exposed TPSA option combinations, and the COSMolKit parse-failure outcome for
all 12 RDKit-rejected rows in the 150-row profile. Python top-level function
stubs now come from metadata adjacent to each Rust `#[pyfunction]`; a Python
test compares runtime public functions with generated declarations and
`__all__` without maintaining another function list. Hand-written declarations
remain only for runtime-created protocol/module and enum surfaces that are not
top-level Rust `#[pyfunction]` values. No source-reproduction marker was
upgraded as part of this product-status audit.

## Validation

- `cargo fmt --all -- --check`: passed.
- `cargo check -p cosmolkit-core --features op-contracts-strict`: passed with
  existing warnings.
- `cargo check -p cosmolkit-py --no-default-features --features dev-stub`:
  passed with existing warnings.
- `cargo test -p cosmolkit-core --release --features op-contracts-strict
  --test rdkit_molecular_descriptor_parity`: passed with 3 tests executed, 0
  failed, 0 ignored, and 0 filtered. The 150-row profile contains 138
  RDKit-success rows and 12 RDKit parse-failure rows; COSMolKit matched every
  corresponding success value/bit field and failure outcome. The compared
  matrix includes `only_heavy=true`, all four `include_hs`/`force` Crippen
  combinations, and all four `force`/`include_sandp` TPSA combinations.
- `python/tests/test_descriptors.py` and `python/tests/test_stub_surface.py`:
  passed with 3 tests executed. The complete Python suite passed with 531 tests
  and 37 explicit skips. `basedpyright python/tests python/examples` reported
  0 errors and 416 existing warnings.
- Stub generation is reproducible: regenerating `python/cosmolkit.pyi` produced
  the same SHA-256
  `8a87ef008192754b4154639f49e3258406289562f283c7f413bed87e0bdbdc2a`.
- Each of the three focused registry tests executed one test and passed:
  `registered_ops_have_support_and_invariant_entries`,
  `parity_registered_ops_have_parity_entries`, and
  `supported_kekulize_runs_through_operation_pipeline_without_changing_source`.
- `cargo test --workspace --release --features
  cosmolkit-core/op-contracts-strict`: passed. The core unit target executed
  2,529 tests with 45 existing ignored tests; the InChI unit target executed
  1,271 tests with no ignored or filtered tests; the RingDecomposer target
  executed 11 tests. All default integration targets and doctests passed.
- The public RDKit expected-data preparation entrypoint validated all 26
  `smiles_small` outputs across 17 families as exact identity/checksum reuse.
  No `.smiles_small.prepare-*` directory remained after publication. The
  declared development environment now includes `epam.indigo`, which the
  ConfSeq reference generator requires; generation no longer depends on an
  undeclared local package.
- The expensive exhaustive SDF API matrix and exhaustive SMILES writer branch
  matrix remained ignored in this default command. Their representative and
  common-branch parity tests did execute and pass; the exhaustive ignored rows
  are not claimed as executed by this validation run.
