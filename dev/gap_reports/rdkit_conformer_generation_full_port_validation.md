# RDKit Conformer Generation Full Port Validation

Date: 2026-05-30

## Scope

This validation record covers the conformer-generation execution plan tail:

- final marker audit
- actionable marker-gap closure
- strict core and workspace validation
- wrapper/examples/docs surfaces already landed in earlier conformer-generation
  plan steps

Primary evidence files:

- `dev/gap_reports/rdkit_conformer_generation_final_marker_audit.md`
- `crates/cosmolkit-core/src/chemistry/distgeom.rs`
- `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs`
- `crates/cosmolkit-core/tests/rdkit_conformer_generation_parity.rs`
- `crates/cosmolkit-core/src/support.rs`

## Commands Executed

```bash
cargo test -p cosmolkit-core --features op-contracts-strict clock_seeded_path -- --nocapture
cargo test -p cosmolkit-core --features op-contracts-strict conformer_generation -- --nocapture
cargo test --workspace --features cosmolkit-core/op-contracts-strict conformer_generation -- --nocapture
cargo test --workspace --features cosmolkit-core/op-contracts-strict
cargo fmt --all
```

## Results

### 1. Clock-seeded RNG regression slice

Command:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict clock_seeded_path -- --nocapture
```

Result:

- passed
- `2` targeted tests passed
- verifies that the previously unported `clock()`-seeded paths now execute
  instead of panicking

Covered tests:

- `vector_set_to_random_clock_seeded_path_returns_normalized_vector`
- `power_eigen_solver_accepts_clock_seeded_path`

### 2. Core conformer-generation strict slice

Command:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict conformer_generation -- --nocapture
```

Result:

- passed
- conformer-generation parity/inventory/registry tests passed under strict mode

Observed conformer-generation tests:

- `conformer_generation_fixture_inventory_has_required_rdkit_sources`
- `conformer_generation_fixture_inventory_matches_rdkit_sources`
- `conformer_generation_parameter_parity_covers_public_parameter_controls`
- `conformer_generation_single_parity_embeds_rdkit_simple_torsion_fixtures`
- `conformer_generation_multi_parity_covers_rdkit_multi_conformer_controls`
- `conformer_generation_stereo_parity_covers_source_backed_stereo_and_torsion_paths`

### 3. Workspace conformer-generation slice

Command:

```bash
cargo test --workspace --features cosmolkit-core/op-contracts-strict conformer_generation -- --nocapture
```

Result:

- passed
- Rust facade conformer-generation re-export coverage passed
- core conformer-generation parity tests passed again in the workspace build

### 4. Whole-workspace strict validation

Command:

```bash
cargo test --workspace --features cosmolkit-core/op-contracts-strict
```

Result:

- passed
- `cosmolkit-core` unit tests: `2220 passed`, `0 failed`, `11 ignored`
- conformer-generation parity integration tests: `6 passed`, `0 failed`
- full workspace doctests passed

Notes:

- The emitted warnings are pre-existing broad workspace warnings outside the
  conformer-generation delta.
- No conformer-generation validation failure was reported in the whole-workspace
  strict run.

### 5. Formatting

Command:

```bash
cargo fmt --all
```

Result:

- passed

## Validation Verdict

The current conformer-generation plan tail is internally consistent under the
required strict validation commands.

The surface that is now validated as present:

- source-backed seeded and unseeded RNG setup
- deterministic explicit-seed single-conformer path
- deterministic batch seed policy
- Rust facade and Python wrapper exports
- parity/inventory/registry tests

The surface that is still not marker-closed:

- the final marker audit still records remaining `RDKit❗✔️`, `RDKit✔️❗`, and
  `RDKit✔️❌` helper-surface markers in bounds building and CrystalFF torsion
  preference collection

Accordingly, the feature should continue to be described as experimental and
substantial rather than fully marker-closed.
