# RDKit 2D Coordinate Frozen Coverage Final Audit

## Baseline

This frozen audit covers the checklist-selected RDKit 2D depiction baseline:

- `third_party/rdkit/Code/GraphMol/Depictor/RDDepictor.h`
- `third_party/rdkit/Code/GraphMol/Depictor/RDDepictor.cpp`
- `third_party/rdkit/Code/GraphMol/Depictor/DepictUtils.h`
- `third_party/rdkit/Code/GraphMol/Depictor/DepictUtils.cpp`
- `third_party/rdkit/Code/GraphMol/Depictor/EmbeddedFrag.h`
- `third_party/rdkit/Code/GraphMol/Depictor/EmbeddedFrag.cpp`
- `third_party/rdkit/Code/GraphMol/Depictor/Templates.h`
- `third_party/rdkit/Code/GraphMol/Depictor/Templates.cpp`
- `third_party/rdkit/Code/GraphMol/Depictor/Basement/Depictor.cpp`

Active Rust surface audited:

- `crates/cosmolkit-core/src/chemistry/coordinates.rs`
- `crates/cosmolkit-core/src/model/molecule.rs`
- `crates/cosmolkit-core/src/operations/ops.rs`
- `crates/cosmolkit-core/src/properties/batch.rs`
- `crates/cosmolkit-core/src/io/molblock.rs`
- `crates/cosmolkit-core/src/properties/draw.rs`
- `crates/cosmolkit-core/src/support.rs`

Validation evidence:

- `cargo fmt --all`
- `cargo check -p cosmolkit-core --features op-contracts-strict`
- `cargo test -p cosmolkit-core --features op-contracts-strict`

All three commands passed in this execution round.

## Verdict

For the selected checklist baseline, the active Rust 2D depiction surface is at
or near `100%` source-port coverage in the sense defined by
`dev/plans/coordinate_2d_rdkit_full_port_checklist.md`:

- the major chemistry entrypoints are implemented
- the direct helper inventory required by the active call chain is present
- public value-style exposure is wired through registered operations and the
  known caller surfaces
- targeted tests and full strict validation pass

No new active chemistry-behavior hole was found in the selected baseline during
this final audit.

## Closed Surface Inventory

The following checklist-critical surfaces are present in the active Rust path:

- `Compute2DCoordParameters`-equivalent parameter surface
- both `compute2DCoords(...)` entrypoint forms
- `preferCoordGen` / `forceRDKit` routing semantics
- `Add2DCoordsToMol`-equivalent wrapper semantics
- ring-template parser, runtime model, registry, set/add/default loading
- `compute2DCoordsMimicDistMat(...)`
- constrained depiction matching against 2D and 3D references
- `straightenDepiction(...)`
- `normalizeDepiction(...)`
- public `with_2d_coordinates` value-style operation with parameterized route
- batch helper exposure
- MolBlock writer fallback generation
- drawing fallback generation

## Residual Exclusions

The remaining exclusions are narrow runtime-wrapper or status-wording items, not
missing chemistry helpers in the selected baseline:

1. CoordGen-backed runtime routing is not available in this Rust build.
   The active Rust behavior is explicit failure when `preferCoordGen=true` or a
   CoordGen-only route would be required. This is a non-Rust-runtime wrapper
   exclusion, not a silent chemistry divergence.
2. Some copied-source blocks still intentionally carry non-`✔️✔️` second-axis
   markers for performance or helper-abstraction reasons. Those markers are not
   evidence of an unimplemented first-axis chemistry gap by themselves.
3. Support status remains `experimental` in `support.rs`. This is an intentional
   public-policy choice to avoid overclaiming blanket RDKit parity outside the
   audited baseline and surrounding dependent feature surfaces.

## Final Statement

Within the selected RDKit 2D depiction baseline used by this checklist, this
port is functionally closed for the active Rust runtime. The remaining caveats
are explicit runtime-wrapper exclusions and support-status conservatism, not
newly discovered missing depiction chemistry behavior.
