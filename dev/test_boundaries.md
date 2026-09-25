# Regression Tests and Corpus Parity

COSMolKit uses two kinds of tests:

- **Regression tests** verify the correctness of a local function or behavior
  using focused cases covering its results, boundaries and errors, including
  previously fixed bugs. They live in the crate that owns the behavior.
- **Corpus parity tests** run registered operations on a selected corpus and
  compare results against a pinned reference implementation. They live in the
  top-level `parity-tests/` crate.

Choose by the purpose and scope of the test: verifying a local function's
specific behavior is regression testing; checking agreement with a reference
across a selected corpus is parity testing. Committing data or storing fixed
expectations does not turn a corpus parity test into a regression test.
A feature may need both. The sections below explain how to write each kind.

## Regression tests: local correctness

Keep these in the owning crate's unit tests or `tests/` directory.

- Check a specific behavior, boundary, error, or previously fixed bug. Specify
  fixed inputs, options and expected values or errors, based on the public
  contract, pinned source or a recorded bug.
- Put small inputs and expected values directly in the test. Store larger fixed
  fixtures under repository-level `testdata/`, with source/version, license and
  checksum information when derived from upstream.
- Read or embed fixtures directly. No shared test-support crate is required.
- Do not generate reference data, invoke upstream implementations, parse
  third-party source code, or implement a corpus runner inside these tests.
- Change fixtures explicitly and review their diffs; tests never refresh them.
- For tables, compare required fields, order and row count, not just sample
  names. Compare floating-point bits when exact bits are required.
- For I/O fixtures, assert the fields that must survive a roundtrip and identify
  allowed losses; parsing successfully is not enough.
- For batch or binding edge cases, assert results, error categories, indices
  and ordering. Keep a bug's input and correct expected result after fixing it.

Examples:

- `fuzzy_and({1:5, 3:-2, 8:4}, {1:3, 3:-4, 9:7})` must return `{1:3, 3:-4}`.
- A residue lookup function must return the expected fields for a known residue
  and its documented result for an unknown name.
- A fixed invalid input must return the documented error category.

A fixture is input or expected data used by a test, not a test category.
Include only data needed to verify the local behavior being tested.

## `parity-tests`: selected corpora compared against a reference

The top-level `parity-tests/` crate owns this workflow. Its chemistry dependency
is public `cosmolkit` with `full`, not individual domain crates.

- A Rust registry defines tasks, operations, widths/options and typed comparison
  fields. Apply those declarations to the selected corpus.
- No task selection means all registered tasks. The corpus is selected separately.
- `run` reuses valid reference data or prepares missing/invalid data using a
  pinned reference implementation, then preflights **all selected tasks**.
- Rust operations start only after the complete selection is ready. Preparation
  failure stops the run before any selected operation executes.
- Compare every selected case, including reference errors where error parity
  is required. Declare normalization or numeric tolerance explicitly; never
  omit mismatching fields or accept extra answers to make a test pass.
- Report case identity, options, expected and actual results on mismatch.
  Known failures still execute; report unexpected passes and changed failures.
  An expected failure does not establish parity.
- Never regenerate expectations to fit a mismatch. Keep generated corpora,
  reference caches and reports outside Git.
- Domain crates contain no corpus preparation or parity-runner logic.

Example: run `fuzzy_and` and `fuzzy_or` on 5,000 fingerprint pairs across both
registered index widths, comparing all 20,000 results against pinned RDKit.

Rust is the primary corpus execution path. Python/JS binding checks may use the
full set or a declared subset, with focused coverage for binding, FFI and WASM
differences. This does not claim those adapters or million-row execution are
already implemented.

Do not add a third ad-hoc corpus or timing harness. Performance evidence needs
an explicit baseline; repetition counts and elapsed time alone prove nothing.
