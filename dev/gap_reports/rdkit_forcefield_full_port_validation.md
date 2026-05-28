# RDKit ForceField Full Port Validation

Date: 2026-05-27

## Scope

This validation record covers the strict local validation slice used for the
ported force-field code before the whole-crate strict run:

- `crates/cosmolkit-core/src/chemistry/forcefield/core.rs`
- `crates/cosmolkit-core/src/chemistry/forcefield/uff/*`
- `crates/cosmolkit-core/src/chemistry/forcefield/mmff/*`
- `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/*`

## Commands Executed

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --features op-contracts-strict forcefield -- --nocapture
```

## Results

### 1. `cargo check -p cosmolkit-core --features op-contracts-strict`

- Result: passed
- Observed warnings: `230`
- Failure count: `0`

Notes:

- The warning set remains dominated by pre-existing unused-item and visibility
  warnings outside the active force-field export/coverage steps.
- No force-field validation failure was reported by `cargo check`.

### 2. `cargo test -p cosmolkit-core --features op-contracts-strict forcefield -- --nocapture`

- Result: passed
- Unit test summary:
  - `700 passed`
  - `0 failed`
  - `0 ignored`
  - `1237 filtered out`

Notes:

- The emitted panic messages in the log come from expected `#[should_panic]`
  force-field tests and therefore are not failures.
- The command completed successfully with exit code `0`.

## Validation Verdict

The local strict validation slice for the current force-field ported code is
passing.

This verdict is intentionally narrower than the later whole-crate strict run in
Step 1153. It confirms that the current force-field implementation, export
surface, and branch-coverage-report test slice are internally consistent under
strict mode before the full crate test pass is attempted.
