# RDKit ForceField Unit Branch Coverage

Date: 2026-05-27

## Scope

This report covers the current force-field port surface under:

- `crates/cosmolkit-core/src/chemistry/forcefield/core.rs`
- `crates/cosmolkit-core/src/chemistry/forcefield/mod.rs`
- `crates/cosmolkit-core/src/chemistry/forcefield/uff/*`
- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/*`
- `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/*`

The verdict here is limited to unit-test branch coverage evidence for the
currently ported Rust force-field code. It does not claim whole-crate closure.

## Coverage Verdict

For the currently ported force-field Rust modules, every module with executable
ported branch logic has at least one in-tree unit-test block that exercises its
ported branches, preconditions, unsupported/error paths, or source-specific
formula variants.

Modules that are export-only and contain no branch-bearing port logic are
explicitly marked as `export-only`.

## Module Inventory

### Force-field core

- `crates/cosmolkit-core/src/chemistry/forcefield/core.rs`
  - Verdict: covered
  - Evidence:
    - local `#[cfg(test)] mod tests`
    - force-field constructor/copy/distance/energy/gradient/constraint tests
    - inventory assertions:
      - `forcefield_core_symbol_inventory_names_required_symbols`
      - `mmff_symbol_inventory_names_required_contrib_and_builder_files`

- `crates/cosmolkit-core/src/chemistry/forcefield/mod.rs`
  - Verdict: export-only
  - Evidence:
    - explicit re-export surface only
    - no local executable branch logic

### UFF

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/mod.rs`
  - Verdict: export-only
  - Evidence: explicit re-export surface only

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/params.rs`
  - Verdict: covered
  - Evidence: local `#[cfg(test)] mod tests` covering source parameter lookup and parse branches

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/utils.rs`
  - Verdict: covered
  - Evidence:
    - local `#[cfg(test)] mod tests`
    - formula and boundary tests for bond, angle, nonbonded, and inversion helpers

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/bond_stretch.rs`
  - Verdict: covered
  - Evidence: local tests cover constructor, term accumulation, energy, gradient, and panic branches

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/angle_bend.rs`
  - Verdict: covered
  - Evidence: local tests cover hybridization/formula branches, gradient paths, and preconditions

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/nonbonded.rs`
  - Verdict: covered
  - Evidence: local tests cover vdW energy/gradient branches and owner/input guards

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/inversion.rs`
  - Verdict: covered
  - Evidence: local tests cover trigonal/group-specific coefficient branches, energy, gradient, and guards

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/torsion_angle.rs`
  - Verdict: covered
  - Evidence:
    - local tests cover all source torsion cases that were ported
    - explicit inventory guard for the audited absence of a fabricated helper:
      - `uff_torsionanglecontrib_source_inventory_rejects_missing_calc_torsion_energy`

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/atom_typer.rs`
  - Verdict: covered
  - Evidence:
    - local tests cover charge flags, atom labels, typing, public parameter APIs, and value-style optimization wrappers

- `crates/cosmolkit-core/src/chemistry/forcefield/uff/builder.rs`
  - Verdict: covered
  - Evidence:
    - local tests cover neighbor-matrix relation branches, force-field builder toggles, torsion/ring filters, unsupported SMARTS, and conformer/precondition paths
    - explicit branch fixtures:
      - `branch_coverage_molecule`
      - `uff_builder_build_neighbor_matrix_covers_all_shared_endpoint_branches`

### MMFF

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/mod.rs`
  - Verdict: export-only
  - Evidence: explicit re-export surface only

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/params.rs`
  - Verdict: covered
  - Evidence: local tests cover default-data parsing, malformed-line failures, lookup/fallback branches, and variant-specific tables

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/mol_properties.rs`
  - Verdict: covered
  - Evidence: local tests cover sanitize, aromaticity/property setup, validity, unsupported-feature boundaries, and public API wrappers

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/bond_stretch.rs`
  - Verdict: covered
  - Evidence: local tests cover constructor/add-term/energy/gradient and owner/input precondition branches

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/angle_bend.rs`
  - Verdict: covered
  - Evidence:
    - local tests cover nonlinear vs linear branches, cosine clipping, multi-term accumulation, gradient singularity clamp, and preconditions

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/stretch_bend.rs`
  - Verdict: covered
  - Evidence: local tests cover force-constant selection, multi-term accumulation, zero-delta branches, gradient accumulation, and preconditions

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/oop_bend.rs`
  - Verdict: covered
  - Evidence: local tests cover out-of-plane geometry branches, multi-term accumulation, gradient branches, and preconditions

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/torsion_angle.rs`
  - Verdict: covered
  - Evidence: local tests cover cosine clipping, energy formula branches, degenerate torsion handling, finite-difference gradient checks, and preconditions

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/nonbonded.rs`
  - Verdict: covered
  - Evidence: local tests cover vdW/electrostatic branches, dielectric-model branches, 1-4 handling, and guards

- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/builder.rs`
  - Verdict: covered
  - Evidence:
    - local tests cover conformer selection, term toggles, threshold/fragment filters, default torsion SMARTS, unsupported custom SMARTS, invalid property rejection, and per-term builder branches

### CrystalFF

- `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/mod.rs`
  - Verdict: export-only
  - Evidence: explicit re-export surface only

- `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/torsion_angle_contribs.rs`
  - Verdict: covered
  - Evidence: local tests cover closed-form energy helper, accumulation, empty/degenerate branches, gradient branches, and owner/input guards

- `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/torsion_angle_m6.rs`
  - Verdict: covered
  - Evidence: local tests cover constructor guards, energy formula, zero-force and degenerate gradient branches, and owner/input guards

## Open Coverage Gaps

None at the module-inventory level for the currently ported force-field scope.

This statement is intentionally narrow:

1. it does not claim zero remaining RDKit parity markers across every
   force-field file,
2. it does not claim whole-crate strict validation is already complete,
3. it does not claim export-only modules need standalone branch tests where no
   executable branch logic exists.
