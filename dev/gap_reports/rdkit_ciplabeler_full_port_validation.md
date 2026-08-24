# RDKit CIPLabeler Full-Port Validation

This report records the complete maintained ChEMBL 37 CIPLabeler phase for
the pinned RDKit `2026.03.1` source (`351f8f378f8ad6bbd517980c38896e66bf907af8c`).
The phase is reproducible with:

```text
.venv/bin/python dev/tools/chembl_parity/run.py \
  --corpus-dir <chembl37-128-v1> \
  --run-dir tmp/chembl37-ciplabeler-full-v8 \
  --workers 96 --phase ciplabeler
```

The committed phase definition is `dev/tools/chembl_parity/profiles/complete.json`.
Workers use deterministic 128-way partitioning and compare complete atom and
bond state for full assignment, selected-atom assignment, selected-bond
assignment, and empty-selection dispatch. Each comparison includes descriptor,
neighbor-order, rank, and molecule `_CIPComputed` state, plus exact success or
error status and message.

## Result

| Branch | Records compared | Matching checks | Blocking mismatches |
|---|---:|---:|---:|
| Full assignment | 2,854,362 | 2,854,362 | 0 |
| Selected atom | 2,854,362 | 2,854,362 | 0 |
| Selected bond | 2,854,362 | 2,854,362 | 0 |
| Empty selection | 2,854,362 | 2,854,362 | 0 |

The source corpus contains 2,897,819 records. Fifteen records are rejected by
both parsers and 43,442 accepted records exceed the configured 80-atom CIP
phase boundary; neither category is counted as a parity comparison. All 128
shards completed successfully, with no retained finding examples and no
informational mismatch counters.

The final aggregate is stored in the untracked run artifact
`tmp/chembl37-ciplabeler-full-v8/ciplabeler/aggregate.json`; generated run
outputs remain outside the repository's committed testdata policy.

The same source-aligned sorter-prefix fix is covered by the focused and
maintained 5,000-row Rust gates. It removes the prior three complex-molecule
`Max Iterations Exceeded` divergences without changing the pinned recursion
budget or adding molecule-specific behavior.
