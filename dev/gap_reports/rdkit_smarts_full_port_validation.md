# RDKit SMARTS Full-Port Validation

## Scope and Revision

This report closes the SMARTS plan's validation boundary for the current
workspace source. The Rust package version is `0.2.13`; the pinned RDKit
source/runtime is `2026.03.1` (the vendored tree identifies this release in
`third_party/rdkit/ReleaseNotes.md`). The source revision is
`351f8f378f8ad6bbd517980c38896e66bf907af8`, with archive SHA-256
`cc2e3f5da0f4510754914f63b3f6fa75936ca25be576f58af5c5024fd56f147a`, as
recorded in `dev/source_provenance/rdkit_smarts_2026_03_1.md`.

The claimed boundary is ordinary-molecule SMARTS: parser/compiler, typed query
graph, SMARTS and CXSMARTS writers, recursive query evaluation, direct
substructure matching, Python bindings, and the existing descriptor,
fingerprint, force-field, coordinate, MolBlock, and SDF consumers. Reaction
SMARTS and database/container SMARTS remain explicitly outside the boundary.

## Single-Core Audit

The architecture baseline `crates/cosmolkit-core/tests/smarts_architecture.rs`
passes all four checks: duplicate inventory, existing-consumer regression,
single-core architecture, and canonical-module baseline.

A direct production-tree search also found zero occurrences of the removed
parser/converter names and zero legacy `search::query::parse_smarts` callers;
the only remaining mentions are assertions in this architecture test itself.

The production tree has exactly these ownership points:

| Role | Canonical implementation |
|---|---|
| SMARTS compiler | `crates/cosmolkit-core/src/search/smarts_parse.rs` (`mol_from_smarts`, private parser helpers) |
| Typed query model | `crates/cosmolkit-core/src/search/query.rs` |
| SMARTS/CXSMARTS writer | `crates/cosmolkit-core/src/search/smarts_write.rs` |
| Matcher | `crates/cosmolkit-core/src/search/substruct.rs` |
| Python parser/writer/matcher facade | `python/src/lib.rs`, delegating to the Rust core |

The duplicate search found zero production occurrences of `SmartsMolecule`,
`build_query_molecule`, consumer-local `ParserState`, or imports/calls to an
alternate parser in `search::query`. No consumer-local parser or recursive
reparse path exists. All consumers compile and match through the canonical
query-bearing `Molecule` path.

## Provenance and Fixtures

The committed `testdata/smarts/corpus/rdkit_source_cases.json` contains 162
pinned RDKit 2026.03.1 source rows: 91 distinct accepted parser inputs, nine
rejected inputs, and 62 matcher inputs. Its ignored expected-data profile is
reproduced through `tools/testdata/rdkit/generate_all.py --profile
smarts_source --suite smarts`; coverage CI performs that preparation before
Rust tests. The manifest records generator, input, schema, output checksum,
and pinned reference identity. Rust test execution reads the validated
manifest and never falls back to a runtime RDKit invocation.

The golden comparisons cover parser acceptance, atom/bond counts, exact
molecule SMARTS output, and exact ordered atom mappings. Fragment and
CXSMARTS write-parse checks are composition tests, while parameters,
callbacks, Python conversion, architecture, and existing consumers are
covered by separate strict tests. Those boundaries are not attributed to the
162-row golden corpus.

## Validation Commands

All required strict validation completed successfully:

```text
cargo fmt --all
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --release --features op-contracts-strict
cargo test --workspace --release --features cosmolkit-core/op-contracts-strict
uv sync --group dev
env -u CONDA_PREFIX VIRTUAL_ENV="$PWD/.venv" .venv/bin/maturin develop --manifest-path python/Cargo.toml
cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen
.venv/bin/pytest
.venv/bin/python -m sphinx -W -b html python/docs/source python/docs/build/html
```

The strict core and workspace release runs pass, including the SMARTS corpus,
composition, existing-consumer, and architecture tests. Python completed with
`622 passed, 37 skipped`. The generated stub surface tests pass.
Basedpyright reports `0 errors` when run with the project source path and
`--level error`:

```text
PYTHONPATH=python .venv/bin/basedpyright --level error python/tests python/examples
0 errors, 0 warnings, 0 notes
```

The targeted closure checks also pass: the 162-row SMARTS corpus suite, the
four-test SMARTS architecture suite, the validated support metadata unit
test, and `python/tests/test_stub_surface.py` (`2 passed`). `git diff --check`
reports no whitespace errors.

The default basedpyright invocation emits existing warning-level diagnostics;
those do not represent SMARTS type errors. Three Optional match assertions in
`python/tests/test_python_api_sanity.py` were made explicit so the proven
non-None branch is visible to the type checker.

## Support Boundary

`SUBSTRUCTURE_FEATURE` is now `SupportedWithRdkitParity("2026.03.1")` for the
ordinary-molecule SMARTS/query boundary. `dev/parity_scope.md`,
`dev/porting_inventory.md`, the Rust support metadata, Python docstrings, and
the generated stub describe the same scope. The support claim does not cover
reaction SMARTS, database/container SMARTS, or unrelated RDKit query APIs.
The `debug_parse=true` Bison diagnostic stream is also outside the behavioral
boundary and returns a structured unsupported error; it is never silently
ignored by the hand-written parser.

## Marker Review

The source-reproduction protocol remains in force. SMARTS production functions
retain copied RDKit anchors with the two-axis markers; architecture and corpus
tests are the behavioral evidence for the `✔️` first-axis claims. Performance
markers remain visible where local inspection has not proven equivalence or
where an intentional allocation/algorithmic difference is documented. No
marker was upgraded merely because a test happened to pass.

## Result

The ordinary-molecule SMARTS port has one core implementation, no historical
duplicate parser/matcher/writer blocks, explicit consumer routing checks, and
strict Rust/Python/docs validation. Remaining excluded SMARTS families are
named explicitly rather than silently approximated.
