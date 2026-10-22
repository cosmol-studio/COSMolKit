# Local official InChI comparisons

These two ignored integration tests compare exact InChI and InChIKey strings:

- `official_sdf_outputs`: the official `inchi` and `mcule` CI SDF datasets.
- `official_atom_permutation_invariance`: original inputs plus ten atom-order
  permutations, using the official RDKit helper and seed 42 per molecule. Every
  variant is compared with official InChI and the pinned RDKit pipeline.
  Invariance compares the ten rewritten variants with one another, not with the
  original SDF; upstream noninvariance is reported separately.

CK reads every input through `Molecule::from_sdf_with_params`, preserving explicit
hydrogens, without sanitization, with non-strict parsing. There is no second SDF
parser in this crate. The facade dependency is test-only.

Prepare locally (requires CMake, a C compiler, and Python with RDKit):

```bash
.venv/bin/python crates/cosmolkit-inchi/tests/support/prepare_official.py --jobs 8
```

The command compiles the unmodified pinned `third_party/InChI` shared library and
prints the reference path. Use that path:

```bash
CK_INCHI_OFFICIAL_REFERENCE=/absolute/path/to/reference.jsonl \
  cargo test -p cosmolkit-inchi --release --test official_suite -- --ignored --nocapture
```

To reuse a previously prepared official reference without rebuilding C or
regenerating official outputs:

```bash
.venv/bin/python crates/cosmolkit-inchi/tests/support/prepare_official.py \
  --reuse-official /absolute/path/to/reference.jsonl
```

This writes a separate companion reference containing RDKit 2026.09.1 outcomes
for the exact same MolBlocks; the original reference remains unchanged. Use the
printed companion path for Cargo comparison.

For a bounded smoke run, prepare with `--limit 2` (two original records per
dataset). This is not a full dataset run. References and mismatch reports remain
under `target/inchi-official-special/`; valid completed references are reused.
Cache identity includes upstream contents, preparation code, RDKit version and
record limit. Interrupted preparation leaves a `.partial.jsonl`, never a valid
cache. Preparation process failure, including a native crash, fails the command;
do not use its incomplete reference.

Log/message text and AuxInfo are not parity assertions. Source rejection is
reported separately from successful output comparisons, retaining CK's outcome.
RDKit is the SDF interpretation contract. Exact InChI/InChIKey matches with RDKit
that differ from official direct-SDF parsing are reported separately, not as
errors. Matching generation errors preserve the valence detail or the complete
unkekulized atom list; arbitrary errors are never accepted as equivalent.
Missing permutation inputs, CK panics and CK/RDKit mismatches fail the comparison.
Reports retain the exact input and expected/actual output. Upstream atom-order noninvariance is not
counted as a CK mismatch when CK matches the upstream output for each variant.

`Molecule::to_inchi_result()` and `to_inchi_result_with_params()` expose the InChI,
native return code, AuxInfo, message, log and structured adapter diagnostics.
Python exposes these as readonly properties; JS uses camelCase. The existing
`to_inchi()` string API remains unchanged. Diagnostics are available to callers
without making their wording a chemical parity requirement.
