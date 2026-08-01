# SMILES Corpora

## `smiles_small.smi`

This is the 150-row repository-curated regression corpus formerly stored as
`tests/smiles.smi`. Migration preserved every row and its order without
filtering. The historical rationale for selecting individual rows was not
recorded, so this corpus must not be attributed to an external dataset.

## `smiles_5000.smi`

This is the complete 5,000-row `SMILES.csv` attachment from
[`kent-tokyo/chematic#11`](https://github.com/kent-tokyo/chematic/issues/11),
in source order. The reproducible normalization is:

1. Remove the single `SMILES` CSV header row.
2. Remove CR bytes from the attachment's CRLF line endings.
3. Preserve all 5,000 SMILES rows without filtering or reordering.

The source attachment SHA-256 is
`98a272f88ccb4411fd5cad88a77df8ae9964295439dabc0e9e0e82eeaf320e4b`.
The issue describes the attachment as drug-like, natural-product, and
small-molecule edge cases supplied by the issue author.

`source_manifest.jsonl` records record counts, byte lengths, corpus
SHA-256 values, sources, and normalization details.
