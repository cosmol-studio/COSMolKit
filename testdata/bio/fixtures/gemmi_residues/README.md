# Fixed Gemmi residue-table regression fixture

`residues.tsv` is a reviewed static snapshot, not a generated parity cache.
Ordinary Bio tests embed this committed fixture with `include_str!`. They never read
upstream source, invoke Gemmi, or regenerate expected values.

## Provenance

- Project: Gemmi 0.7.5.
- Pinned revision: `5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e`.
- Source: `src/resinfo.cpp`, `static ResidueInfo array[368]`.
- Source URL: https://github.com/project-gemmi/gemmi/blob/5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e/src/resinfo.cpp
- Source SHA-256: `5b97fb945b416f2128ca3f11d43be0b17ac058067fa3d91e2478252623837951`.
- Fixture SHA-256: `ec75f969c2dc9fc0e7d36ec234ce76709add936314b9dfd7b4958dfacc48bb82`.
- Copyright: 2018–2022 Global Phasing Ltd.
- License: Mozilla Public License 2.0; the accompanying `LICENSE.txt` is
  copied from Gemmi. This derived table retains that license context.

## Selection and representation

All 368 initializer rows are retained in source order, including the terminal
unknown entry. No filtering, deduplication, chemistry normalization or values
derived from COSMolKit are used. The header defines seven tab-separated fields:
zero-based index, name, kind, linking type, one-letter code, hydrogen count,
and the decimal f32 weight literal.

The one-letter code may be a literal space; it must not be trimmed. Weight
decimal spellings are copied verbatim with only the C++ `f` suffix removed.
Rust parses them directly as f32, as the previous test did, and compares
`to_bits()` against the implementation. This preserves the existing bitwise
assertion without an intermediate f64 conversion.

## Updating

Updates are explicit reviewed fixture changes, never a test preparation step.
For a deliberate upstream update, identify and hash the new pinned source,
extract the complete initializer in order, remove C++ row punctuation and
`RI::` prefixes, retain the literal field values, and assign consecutive row
indices. Review the table diff and update both identities above. A new source
checkout by itself must not change this fixture or the test expectations.

Validation from the repository root:

```sh
sha256sum testdata/bio/fixtures/gemmi_residues/residues.tsv
cargo test -p cosmolkit-bio --release --test migration_bio_residue
```
