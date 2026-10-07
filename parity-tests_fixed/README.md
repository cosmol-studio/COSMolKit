# Parity tests

Run from the repository root. Each lane has two commands: prepare references,
then compare with Cargo. Ordinary regressions stay in their owning crates.

## Corpus

```bash
cargo run -p cosmolkit-parity-tests-fixed --release -- prepare --corpus smiles_5000 --threads 112
cargo test -p cosmolkit-parity-tests-fixed --release --features cosmolkit/op-contracts-strict --test corpus
```

Available corpora: `smiles_smoke`, `smiles_small`, `smiles_5000`, `bio_small`.
The prepared corpus becomes the test selection; `PARITY_CORPUS` can override it.
Without a task filter, all tasks for that corpus run.
Preparation shows task N/total and a live completed-cases progress bar per task;
parameter combinations for each case remain together. Valid references show `reused`.

To prepare and run just one task:

```bash
cargo run -p cosmolkit-parity-tests-fixed --release -- prepare --corpus smiles_smoke --threads 112 --task smiles_write_smiles
cargo test -p cosmolkit-parity-tests-fixed --release --features cosmolkit/op-contracts-strict --test corpus smiles_write_smiles -- --exact
```

Cargo owns name filtering and scheduling. `smiles_write_smiles` tests 768
profiles per molecule: eight booleans × default/first/last root, matching the
existing writer generator. Random traversal is off. Batch tests are unchanged.

UFF and MMFF single-/multi-conformer optimization use at most **two iterations**.
MMFF covers both MMFF94 and MMFF94s; multi-conformer inputs contain two conformers.
Both libraries receive the same prepared coordinates. MMFF energy/coordinate
comparison uses the existing public MMFF suite's absolute tolerance of `1e-6`;
UFF retains its exact-bit comparison. Parameter-availability queries are included.
Tautomer enumeration and canonicalization are ordinary registered corpus tasks;
the long-conjugated tautomer case is a special regression below.

MACCS, Topological, Layered, Pattern, fuzzy AND and fuzzy OR are registered
SMILES corpus tasks. Each molecule gets **one** reproducible parameter combination
per task (seed `0x434b465020261007`, keyed by task, case ID and SMILES).
Parameters are saved in `expected/` inputs; RDKit and CK consume the same values.
MACCS compares raw 167-bit and public 166-bit results; Layered also compares
seeded atom counts. Fuzzy operations pair each molecule with the next (wrapping
at the end), use Morgan counts, and sample signed counts and 32-/64-bit indices.
No separate fingerprint-pairs corpus is needed. Results compare exactly.

## Special regressions

```bash
cargo run -p cosmolkit-parity-tests-fixed --release -- prepare --special all --threads 112
cargo test -p cosmolkit-parity-tests-fixed --release --features cosmolkit/op-contracts-strict --test special_regression
```

For one special regression, replace `all` with `structure_tags` or
`tautomer_long_conjugated`, then use the same name as Cargo's test filter.

Inputs: `testdata/`. Generated references and checksums: `expected/`.
Results: `reports/`. Preparation uses `.venv/bin/python` with pinned RDKit
and Gemmi packages. Tests validate all selected references before
CK calls; missing/stale references fail, never generate or overwrite expectations.
