# Shared Regression Fixtures

In COSMolKit 0.5.0, this directory holds shared ordinary-regression fixtures,
source manifests and historical provenance. It is not a corpus preparation
entrypoint. Existing corpus files remain available to their current consumers;
this documentation change does not move or delete inputs.

Current corpus and designated special-regression inputs belong in
`parity-tests_fixed/testdata/`; generated references belong in that package's
`expected/`, and comparison reports in `reports/`. Use the single
[prepare / Cargo test workflow](../parity-tests_fixed/README.md).
Do not duplicate preparation instructions in individual fixture directories.

The current reference is RDKit `2026.03.6`, recorded in
[`reference/rdkit.json`](reference/rdkit.json). Older versions in fixture
READMEs and source manifests are the original provenance of those fixed
inputs, not a competing current reference pin. Do not rewrite their source
versions or checksums merely to match the current environment.

Ordinary regressions live in their owning crates. They use small inline cases
or fixed shared fixtures, never invoke reference generators, and do not depend
on generated corpus expectations. Pinned `third_party/` fixtures may be read
at test runtime under the [repository policy](../dev/repository_organization_policy.md);
production and package builds must not require them.

Keep source versions, selection notes, licenses and checksums with existing
fixtures. An unresolved provenance gap remains unresolved; do not infer a
source from the filename.

[VALIDATION.md](../VALIDATION.md) records completed COSMolKit 0.5.0 validation:
18,000 ChEMBL 37 source records and 20,505,775 parameterized comparisons against
RDKit 2026.03.6, with zero mismatches under the documented exact, bit-exact and
tolerance rules. Reference timeouts and upstream crashes are reported separately.
Million-scale 0.5.0 validation remains pending.
