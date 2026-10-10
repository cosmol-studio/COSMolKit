# Parity tests

This is the current COSMolKit 0.5.0 corpus and special-regression workflow.
No differences have been observed on the known corpus of several hundred
thousand molecules; million-scale 0.5.0 parity validation is pending.
Historical 0.3.0 results are recorded
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

For the audited recovery cache, add `--reuse-from /path/to/recovery-checkout`
to prepare. Only the approved recipe pair and identical complete inputs can
be imported; three renamed input keys are converted without changing native
outputs. Original manifests remain in an import receipt. Different parameters
use normal reference generation; comparison validation is unchanged.

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
Both libraries receive the same prepared coordinates.
UFF/MMFF optimization geometry preparation has a user-selected 60-second
deadline: RDKit's native timeout is backed by a supervised process deadline.
Timed-out inputs remain in the references and reports, with their case IDs and
timeout mechanism, and are excluded from comparison. Reports distinguish
total rows, actual comparisons, mismatches and timeouts; timeouts never count
as matches, and a task with zero actual comparisons fails. Other preparation
errors retain their existing comparisons. Geometry, seeds and optimization
parameters are unchanged for completed inputs.

MMFF energy/coordinate comparison uses the existing public MMFF suite's absolute tolerance of `1e-6`;
UFF retains its exact-bit comparison. Parameter-availability queries are included.
`mmff_force_field_smiles` and `uff_force_field_smiles` test the owned persistent
evaluators separately: fixed-seed arbitrary initial coordinates, at most two
minimization iterations, and bitexact initial/final energy, gradient and positions.
Their common coordinate bits are saved in prepared inputs; no embedding or
tolerance comparison is used. Both sides parse the original SMILES and add
explicit hydrogens in source order; MOL serialization is not part of these
evaluator tests. Filter preparation and Cargo by either task name.
MMFF parameter unavailability and UFF's source-defined missing TBP center
parameter error are compared separately, not counted as successful minimizations.
Rejected SMILES compare at the parsing boundary without coordinates; missing
geometry for an accepted molecule remains an error.
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

MACCS, Topological, Layered, Pattern, Avalon, fuzzy AND and fuzzy OR are registered
SMILES corpus tasks. Each molecule gets **one** reproducible parameter combination
per task (seed `0x434b465020261007`, keyed by task, case ID and SMILES).
Parameters are saved in `expected/` inputs; RDKit and CK consume the same values.
MACCS compares raw 167-bit and public 166-bit results; Layered also compares
seeded atom counts. Fuzzy operations pair each molecule with the next (wrapping
at the end), use Morgan counts, and sample signed counts and 32-/64-bit indices.
No separate fingerprint-pairs corpus is needed. Results compare exactly.
Invalid SMILES compare as parsing rejections, including which operand failed
for fuzzy operations; they are not empty fingerprints or preparation failures.
Avalon samples bit-vector sizes (including non-byte-aligned sizes), query mode
and native feature masks; both sides receive the recorded explicit flags.
`fragments_smiles` compares source-ordered sanitized fragment SMILES;
`largest_fragment_smiles` selects the most atoms, retaining the last tie
(CK's existing convenience rule, not MolStandardize's chooser policy).
The upstream unrooted-linear Layered branch is executed in an isolated
process because the pinned source can crash while treating atom-path indices
as bond indices. Every original case and seeded profile is retained. Native
segmentation faults retain their complete inputs, exit code, process ID,
diagnostics and CK outcome in reports under `UpstreamReferenceCrash`. This
explicitly approved exception applies only to that unrooted-linear Layered
branch with SIGSEGV/Windows access-violation exits. Reports count these rows
separately as `upstream_crashed`, never as matches or preparation timeouts.
Other native failures remain failures; zero actual comparisons still fails.
All successful native results retain the existing exact comparison.

## Special regressions

```bash
cargo run -p cosmolkit-parity-tests-fixed --profile dev-test -- prepare --special all --threads 112
cargo test -p cosmolkit-parity-tests-fixed --profile dev-test --features cosmolkit/op-contracts-strict --test special_regression
```

For one special regression, replace `all` with `structure_tags`,
`tautomer_long_conjugated`, `tautomer_focused`, `molalign_focused`, or
`bio_mmcif_switches`, `forcefield_optimizers`, `mmff_builtin`, `mcs_upstream`, or
`mcs_jnk1`, then use the same name as
Cargo's test filter. The focused tautomer matrix keeps 18 inputs, eight profiles
and all 136 valid enumeration branches; it is not part of ordinary crate tests.
MolAlign retains its 14 fixed boundary calls, including typed errors, through
the same preparation and comparison workflow; no external oracle directory
environment variable is needed.

`conformer_fixed19`, `conformer_library` and `forcefield_properties` retain
the 19 fixed embedding cases, 152 seeded library rows and 152 UFF/MMFF
parameter rows previously dependent on owner-crate oracle caches. Inputs,
parameters and comparison assertions are unchanged: embedding coordinates
use 1e-6, MMFF formal/partial charges use 1e-12, and discrete results are exact.
They use the same prepare/Cargo entrypoints, with no reference-directory
environment variable. Ordinary crate tests need no generated reference files.

`forcefield_optimizers` retains 20 fixed MMFF and four fixed UFF counterexamples.
It compares parameter availability, initial energy/gradient and zero-step status,
then single- and two-conformer optimization with **at most two iterations**.
All energy, gradient and coordinate fields compare **bitexact**; status and
counts compare exactly. The original seed and CXSMILES input geometry are retained.
`mmff_builtin` retains all 2,052 rows of the four upstream dative/hypervalent
matrices, comparing availability and every MMFF94/MMFF94s atom type.
These references are regenerated with the current pin, not copied from old goldens.

MCS uses the same two-command workflow:

```bash
cargo run -p cosmolkit-parity-tests-fixed --profile dev-test -- prepare --special mcs_upstream --threads 112
cargo test -p cosmolkit-parity-tests-fixed --profile dev-test --features cosmolkit/op-contracts-strict --test special_regression mcs_upstream
```

Replace `mcs_upstream` in both commands with `mcs_jnk1` for the complete
**210 pairs** from 21 upstream JNK1 ligands. `mcs_upstream` retains **44 calls**:
42 calls from 14 selected upstream Python test methods, plus the C++
Github9034/StoreAll and JNK1/MaxDistance cases. It is not the entire upstream
callback/type-validation suite. Original SMILES/MOL text, parsing options,
source paths/lines/checksums and comparison parameters are frozen under
`testdata/special/mcs_*.json`; preparation and tests need no external fixture
checkout. Each call requests a 30-second timeout; original upstream timeout
values remain recorded separately.

Each case is an independent Cargo test. Counts, completion, SMARTS text,
degenerate SMARTS keys, serialized query and query-to-input match results
compare exactly; input molecule binaries must remain unchanged. SMARTS are
never normalized to hide differences. Native/CK errors and incomplete results
remain failing observations, even when both sides fail or counts agree.
Per-case results are saved under `reports/mcs_upstream/` and `reports/mcs_jnk1/`.

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
