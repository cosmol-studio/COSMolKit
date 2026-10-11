# Inputs

`corpora/smiles_small.smi` retains the existing 152-record repository smoke
corpus, byte-for-byte from `parity-tests_fixed/testdata/corpora/smiles_small.smi`.
This reuses input data only, not the old runner's architecture or references.
The original corpus and its historical provenance remain unchanged.

SMILES inputs contain one record per nonempty, non-comment line; the first
whitespace-separated field is SMILES. Additional fields are source annotations.
Order and invalid SMILES are retained. Runtime `--offset` and `--limit` select
records, not parameter-expanded rows. Large external selections remain ignored.
