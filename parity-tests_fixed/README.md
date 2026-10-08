# Parity tests

This is the current COSMolKit 0.5.0 corpus and special-regression workflow.
Complete 0.5.0 validation is pending; historical 0.3.0 results are recorded
separately in [VALIDATION.md](../VALIDATION.md).

Run from the repository root. Each lane has two commands: prepare references,
then compare with Cargo. Ordinary regressions stay in their owning crates.

## Corpus

```bash
cargo run -p cosmolkit-parity-tests-fixed --profile dev-test -- prepare --corpus smiles_5000 --threads 112
cargo test -p cosmolkit-parity-tests-fixed --profile dev-test --features cosmolkit/op-contracts-strict --test corpus
```

Available corpora: `smiles_smoke`, `smiles_small`, `smiles_5000`, `bio_small`.
The prepared corpus becomes the test selection; `PARITY_CORPUS` can override it.
Without a task filter, all tasks for that corpus run.
Preparation shows task N/total and a live completed-cases progress bar per task;
parameter combinations for each case remain together. Valid references show `reused`.

To prepare and run just one task:

```bash
cargo run -p cosmolkit-parity-tests-fixed --profile dev-test -- prepare --corpus smiles_smoke --threads 112 --task smiles_write_smiles
cargo test -p cosmolkit-parity-tests-fixed --profile dev-test --features cosmolkit/op-contracts-strict --test corpus smiles_write_smiles -- --exact
```

Cargo owns name filtering and scheduling. `smiles_write_smiles` tests 768
profiles per molecule: eight booleans × default/first/last root, matching the
existing writer generator. Random traversal is off. Batch tests are unchanged.

UFF and MMFF single-/multi-conformer optimization use at most **two iterations**.
MMFF covers both MMFF94 and MMFF94s; multi-conformer inputs contain two conformers.
Both libraries receive the same prepared coordinates. MMFF energy/coordinate
comparison uses the existing public MMFF suite's absolute tolerance of `1e-6`;
UFF retains its exact-bit comparison. Parameter-availability queries are included.
`mmff_force_field_smiles` and `uff_force_field_smiles` test the owned persistent
evaluators separately: fixed-seed arbitrary initial coordinates, at most two
minimization iterations, and bitexact initial/final energy, gradient and positions.
Their common coordinate bits are saved in prepared inputs; no embedding or
tolerance comparison is used. Filter preparation and Cargo by either task name.
MMFF parameter unavailability and UFF's source-defined missing TBP center
parameter error are compared separately, not counted as successful minimizations.
Tautomer enumeration and canonicalization are ordinary registered corpus tasks;
the long-conjugated tautomer case is a special regression below.

`molalign_smiles` preserves the six-operation alignment/RMSD comparison on
`smiles_small` (152 molecules) or `smiles_5000`. It compares RMSD, transforms,
atom maps, conformer order/IDs, every coordinate and source preservation.
Floating fields retain the original absolute `1e-8` tolerance; maps and IDs
are exact. The same prepare command supplies native conformers and parameters.
Reference transport preserves native float bits in both preparation and tests
(`serde_json/float_roundtrip`); it must not round the inputs independently of
the native results.

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
cargo run -p cosmolkit-parity-tests-fixed --profile dev-test -- prepare --special all --threads 112
cargo test -p cosmolkit-parity-tests-fixed --profile dev-test --features cosmolkit/op-contracts-strict --test special_regression
```

For one special regression, replace `all` with `structure_tags`,
`tautomer_long_conjugated`, `tautomer_focused`, `molalign_focused`, or
`bio_mmcif_switches`, then use the same name as
Cargo's test filter. The focused tautomer matrix keeps 18 inputs, eight profiles
and all 136 valid enumeration branches; it is not part of ordinary crate tests.
MolAlign retains its 14 fixed boundary calls, including typed errors, through
the same preparation and comparison workflow; no external oracle directory
environment variable is needed.

`bio_mmcif_switches` reuses the existing `bio/sample.pdb` and `bio/sample.cif`
inputs. It compares exact mmCIF output bytes with pinned Gemmi for its 32
Python-exposed single-switch profiles under both global defaults (128 calls).
The `modres` single-switch profile is not covered: Gemmi 0.7.5 does not expose
it in Python. Preparation does not compile or invoke a C++ oracle. This is separate
from `bio_pdb_output_pdb`/`bio_pdb_output_cif`, which both produce PDB output.

Inputs: `testdata/`. Generated references and checksums: `expected/`.
Results: `reports/`. Preparation uses `.venv/bin/python` with pinned RDKit
and Gemmi packages. Tests validate all selected references before
CK calls; missing/stale references fail, never generate or overwrite expectations.
