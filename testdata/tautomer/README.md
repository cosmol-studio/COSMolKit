# Tautomer Test Data

`fixtures/rdkit/` contains focused source-branch inputs. Generated pinned-RDKit
observations are prepared under the ignored `expected/rdkit/<profile>/`
directories by `tools/testdata/rdkit/generate_all.py`.

`corpus/rdkit/1kPCS_tautomer.csv.gz` and
`corpus/rdkit/100kPCS_tautomer.csv.gz` are exact copies of the corresponding
RDKit PCS tautomer corpora. Each CSV row contains an input SMILES in column one
and its expected canonical tautomer SMILES in column two.

| Field | Value |
|---|---|
| Upstream project | RDKit |
| Source revision | `351f8f378f8ad6bbd517980c38896e66bf907af8c` |
| Source paths | `rdkit/Chem/MolStandardize/test_data/1kPCS_tautomer.csv.gz`; `rdkit/Chem/MolStandardize/test_data/100kPCS_tautomer.csv.gz` |
| Selection | Complete files, byte-for-byte; no row filtering or reordering |
| License context | RDKit BSD 3-Clause license |
| 1k record count | 1,000 |
| 1k SHA-256 | `e1cb1c99502d1644b53af17327cfe2594db22ebfe242751644bec3dcc76e675c` |
| 100k record count | 99,991 |
| 100k SHA-256 | `873790d2e3ab7d555c879d3dad2a71a11d62fbda9f25f1d66c03fa4a95166425` |

Prepare the complete reference observations with:

```bash
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python \
  --profile tautomer_pcs_1k \
  --suite tautomer \
  --jobs 4

.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python \
  --profile tautomer_pcs_100k \
  --suite tautomer \
  --jobs 4
```

The manifest binds each generated family to the corpus checksum, pinned RDKit
identity, generator sources, options, record count, and output checksum.
