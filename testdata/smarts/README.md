# SMARTS Test Data

`corpus/rdkit_source_cases.json` is a committed, source-derived input corpus
for ordinary-molecule SMARTS parsing, writing, and direct substructure
matching. The cases come from the pinned RDKit 2026.03.1 source revision
recorded in the corpus and preserve the source test function for every row.

The parser portion contains all 91 distinct `testPass` inputs and all nine
rejecting `testFail` inputs from `Code/GraphMol/SmilesParse/smatest.cpp`.
Repeated success sentinels used only to prove parser recovery are not copied.
The matcher portion contains 62 rows from the foundational `testMatches`,
`testMatches2`, and `testMatches3` branches plus the complete `k`/ring-range
extension table in `smarts_catch_tests.cpp`; explicit-H expansion is tested as
a separate molecule operation.

Generate the ignored RDKit expected data with:

```bash
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python --profile smarts_source --suite smarts
```

The generated `smarts.jsonl` records parse acceptance, atom/bond counts,
exact RDKit SMARTS output, and exact ordered target-atom mappings. The Rust
suite also exercises SMARTS/fragment/CXSMARTS write-parse composition, without
presenting those non-golden composition checks as independent RDKit parity.
Normal test execution only reads the validated manifest and does not invoke
RDKit. Coverage CI prepares this dedicated profile before running Rust tests.
