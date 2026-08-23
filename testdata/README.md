# Test Data

Repository-owned reusable test inputs live under `testdata/`, grouped by
stable format or behavior domain. Tests must use these files instead of
`third_party/` trees or crate-local fixture copies.

## Layout

- `smiles/corpus/`: committed `smiles_small` and `smiles_5000` input corpora.
- `<domain>/fixtures/`: committed file inputs with nearby provenance metadata.
- `<domain>/expected/<reference>/<profile>/`: generated, uncommitted expected
  data plus `manifest.json`.
- `inchi/reference/`: committed official-InChI oracle schema and provenance.
- `known_failures/`: explicit executable known-failure declarations.
- `topology/corpus/`: operation-invariant corpora.

Generated expected data is ignored by Git. It may be prepared locally or in
CI and cached, but a restored cache is reusable only when its manifest,
reference version, generator/input/options/platform identity, output checksum,
and record count match exactly.

Coverage CI uses cache schema `v4`. Cache saves have a unique workflow
run/attempt suffix and restores use the stable identity prefix. If preparation
rejects and regenerates a restored candidate, the replacement can therefore be
saved without attempting to overwrite GitHub Actions' immutable cache entry.

Fixture directories and corpus directories carry a `README.md` and the
provenance manifest named by that README. Most use `source_manifest.jsonl`;
the conformer fixture family retains the descriptive
`rdkit_inventory.jsonl` name used by its tests and generator. These manifests
record source classification, selection notes, byte lengths, and SHA-256
values. A manifest entry marked
`historical_provenance_incomplete` is an explicit unresolved provenance gap;
it must not be upgraded by inference from a filename alone.

## RDKit Expected Data

RDKit `2026.03.1` is the pinned generator reference. Set up the project Python
environment and use the single public preparation entrypoint:

```bash
uv sync --group dev
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python \
  --profile smiles_small \
  --suite default \
  --jobs 4
```

The strict 5000-row InChI preparation is:

```bash
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python \
  --profile smiles_5000 \
  --suite inchi \
  --jobs 1
```

The lightweight strict force-field audit deliberately separates broad
parameter coverage from expensive optimizer parity. Prepare its all-row
coverage oracle with:

```bash
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python \
  --profile smiles_5000 \
  --suite forcefield_coverage \
  --jobs 1
```

Run all 5,000 parameter/type/charge comparisons, then the complete curated
energy/gradient/optimizer suite:

```bash
COSMOLKIT_PARITY_PROFILE=smiles_5000 \
  cargo test -p cosmolkit-core --test rdkit_forcefield_coverage_parity

COSMOLKIT_PARITY_PROFILE=smiles_small \
  cargo test -p cosmolkit-core --test rdkit_forcefield_params_parity
```

The first command does not claim optimizer coverage. The second retains every
eligible row, comparison field, and the declared `1e-6` energy, gradient, and
coordinate tolerances.

Omit `--clean` to validate and reuse an exact cache. Use `--clean` only when
an intentional full regeneration is required. Individual `_generate_*.py`
files are internal helpers and are not supported public workflows.

## Gemmi Expected Data

Gemmi `0.7.5` at commit `5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e`
is the pinned structural-writer reference. Prepare the exact mmCIF writer
profile through its single public entrypoint:

```bash
.venv/bin/python tools/testdata/gemmi/generate_all.py \
  --profile bio_mmcif_writer \
  --suite mmcif_writer \
  --jobs 4
```

The generator verifies the pinned source hashes before compiling its narrow
oracle under `target/`. Coverage CI initializes only the Gemmi submodule and
runs this preparation before Rust tests; generated output remains ignored in
accordance with the repository-wide expected-data policy.

`cargo test` never generates expected data. Missing, stale, corrupt, or
incomplete data fails with the corresponding prepare command; tests do not
skip or silently fall back.

## InChI Oracle

The independent official InChI C oracle runner is
`tools/oracles/official_inchi/run.sh`. It verifies the pinned InChI v1.07.5
source identity before building under `target/` and is test-only. Production
Rust never links or invokes it.

The public Python InChI parity boundary is `Molecule.to_inchi()`,
`Molecule.to_inchi_key()`, top-level `inchi_to_key()`, and
`Molecule.from_inchi()`.
Source-defined behavior is compared exactly. Official-C undefined behavior is
mapped to a deterministic structured Rust error and is not presented as an
exact C result.

## Strict Test Commands

```bash
cargo test -p cosmolkit-core --release --features op-contracts-strict

cargo test -p cosmolkit-core --release --features op-contracts-strict \
  --test rdkit_inchi_parity \
  inchi_matches_pinned_rdkit_for_every_active_profile_row -- --exact

COSMOLKIT_PARITY_PROFILE=smiles_5000 \
  cargo test -p cosmolkit-core --release --features op-contracts-strict \
  --test rdkit_inchi_parity \
  inchi_matches_pinned_rdkit_for_every_active_profile_row -- --exact
```

The focused commands above are successful only when the log reports exactly
one executed test and `1 passed`; a filtered run with zero executed tests is
not validation.

Feature-specific comparison boundaries are documented in
[`dev/parity_scope.md`](../dev/parity_scope.md). The normative layout and
cache rules are in
[`dev/repository_organization_policy.md`](../dev/repository_organization_policy.md).
