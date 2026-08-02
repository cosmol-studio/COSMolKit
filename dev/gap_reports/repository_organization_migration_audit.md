# Repository Organization Migration Audit

This audit records the completion boundary for
`dev/repository_organization_policy.md`. It distinguishes public integration
suites from private source-port tests that must remain inline.

## Generated expected data

- RDKit is the only generated expected-data family currently defined.
- `tools/testdata/rdkit/generate_all.py` is its single public prepare entrypoint.
- The obsolete standalone XYZ built-in helper was removed after its fixture
  coverage moved into `builtin_fixture_migration_golden`; no unregistered
  expected-data generator remains.
- Gemmi and official InChI currently have no generated expected-data family.
  Official InChI behavior checks are live test oracles under
  `tools/oracles/official_inchi/`; no empty generator entrypoint is required.
- Expected outputs are uncommitted, manifest-identified, checksum-validated,
  and cacheable only for exact generation identity.

### Identity and schema validation

- `testdata/reference/rdkit.json` pins the reference implementation identity.
- `tools/testdata/rdkit/_expected_schema.py` defines and validates each JSONL
  output schema independently of the generator that writes the output.
- Every generated profile/domain directory containing outputs has a
  `manifest.json`. Each output entry records generator and dependency
  checksums, input identity, options, schema, record count, byte length, and
  SHA-256.
- `cosmolkit-test-support` rejects missing, stale, incomplete, or corrupt
  requested output entries. Its ordinary lookup validates only the requested
  output and does not scan unrelated outputs.
- The 26 `_generate_*.py` helper files exactly match the 26 entries in
  `GENERATOR_SPECS`; there are no unregistered generation helpers.

### Checksum execution model

- Preparation validates the complete selected generated family once before it
  reuses or publishes that family. This is the only phase that scans every
  selected output.
- An ordinary test binary validates only the manifest entry, checksum, and
  record count for the output it consumes. It does not hash unrelated outputs.
- Therefore a multi-gigabyte expected-data family is not rescanned once per
  test binary. Cache correctness is established by preparation, while each
  consumer still detects corruption of its own requested file.

### Immutable CI cache repair

- Coverage CI uses cache schema `v2` and restores by stable complete-identity
  prefix.
- Saves use a key ending in the unique workflow run id and attempt. A corrupt
  immutable cache can therefore be rejected, regenerated, and saved under a
  new key instead of attempting an impossible overwrite.
- A replacement cache is saved only when preparation reports regeneration;
  exact validated restores are reused without duplicate cache writes.

## Orphan expected-output quarantine

Fifteen files found under active expected-data paths had no valid identity
manifest and were not admitted as current expected data. They were moved,
without byte changes, to:

```text
tmp/unvalidated_rdkit_expected_backup/
  repository_organization_migration_20260731/
```

`original_paths.sha256` and `original_paths.sizes` preserve their former
paths, checksums, and sizes. All 15 backup files still match both records.
The backup is evidence, not an active expected-data family, and must not be
silently reused by tests.

The backup root also contains 22 older top-level JSONL files, separate from
the dated 15-file quarantine. They have no complete generation-identity
manifest or original-path baseline and are therefore retained only as
unvalidated local evidence. They are not counted among the 15 checksum-backed
moves, are not referenced by active tests or tools, and must not be promoted
to expected data based on filename or current checksum alone. At audit time,
the complete backup directory occupied approximately 4.3 GiB and the dated
quarantine approximately 907 MiB. No backup file was deleted.

## Fixture and corpus provenance

- Every fixture/corpus directory has a README and a machine-readable
  provenance manifest. Most use `source_manifest.jsonl`; conformer fixtures
  use the existing `rdkit_inventory.jsonl` name declared by their README.
- Manifest coverage, byte length, and SHA-256 were checked for all 1,107
  committed fixture/corpus files.
- Three historical ligand fixtures retain
  `historical_provenance_incomplete` because their exact upstream source or
  export command cannot be recovered from repository evidence:
  `testdata/mol2/fixtures/3rj7_ligand.mol2`,
  `testdata/stereo/fixtures/1aid_ligand.mol2`, and
  `testdata/stereo/fixtures/1aid_ligand.sdf`.
- Their current bytes and checksums are recorded. This is not a claim of
  complete historical provenance.

## Integration placement

- Public MolBlock V2000 corpus parity is in
  `crates/cosmolkit-core/tests/rdkit_molblock_v2000_write_parity.rs`.
- Other public generated-output parity suites are under
  `crates/cosmolkit-core/tests/`.

## Retained inline parity

The following categories remain inline because moving them would require
widening production visibility solely for tests:

- `properties/descriptors.rs::qed_tests`: the DeleteSubstructs matrix verifies
  private `rdkit_qed_delete_substructs` control flow and its private SMARTS
  construction helper.
- `properties/rdkit_prepared_draw_parity.rs`: verifies the private
  `PreparedDrawMolecule` intermediate snapshot and private batch preparation
  path before SVG serialization.
- `chemistry/forcefield/crystalff/torsion_preferences.rs::tests`: focused tests
  of private CrystalFF source-port functions use committed shared fixtures.
  Repository-root discovery is delegated to `cosmolkit-test-support`.

These tests resolve shared files and generated outputs through
`cosmolkit-test-support` under `cfg(test)` and do not add a production runtime
dependency.

## Python coverage review

- Python parity tests use representative cases rather than the complete Rust
  parity corpora.
- Fingerprint tests cover Python-visible bit vectors and additional-output
  container conversion.
- Substructure tests cover Python-visible atom-mapping tuple conversion and
  ordering.
- `from_rdkit` tests cover the Python object bridge, enum conversion,
  conformer conversion, and Python exception behavior. This boundary can
  independently alter results and therefore retains its representative parity
  cases.
- Complete chemistry corpora remain owned by Rust integration tests.

## Layout result

- Reusable inputs are under `testdata/<domain>/fixtures` or `corpus`.
- Generated outputs are under ignored
  `testdata/<domain>/expected/<reference>/<profile>` directories.
- `tests/` contains no data or executable tests.
- Crate-local and Python-local fixture copies are absent.
- Runtime tests do not read `third_party/rdkit`; upstream source paths remain
  only as provenance and source-reproduction anchors.
- Active instructions, CI configuration, and generator defaults use
  `testdata/`; old paths remain only in historical completed plans.

## Final validation

The following checks passed after the migration:

```text
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --profile smiles_small --suite all
  -> 26 outputs validated; all 17 selected domains reused on the second run

cargo fmt --all -- --check
  -> passed

cargo check -p cosmolkit-core --features op-contracts-strict
  -> passed (existing warnings only)

cargo test -p cosmolkit-core --release --features op-contracts-strict
  -> 2529 unit tests passed, 45 ignored, 0 filtered out;
     all integration tests and 6 doctests passed

cargo test --workspace --release \
  --features cosmolkit-core/op-contracts-strict
  -> passed with exit code 0 across every workspace test binary and doctest
  -> cosmolkit-inchi: 1271 passed, 0 failed, 0 filtered out
  -> official InChI provenance integration test: 1 passed
  -> cosmolkit-ringdecomposer: 11 passed
  -> cosmolkit-test-support: 4 passed

.venv/bin/python -m unittest discover \
  -s tools/testdata/rdkit -p 'test_*.py' -v
  -> 4 passed

.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --profile smiles_5000 --suite inchi --jobs 1
  -> 1 output reused after exact identity and output-checksum validation

coverage.yml YAML parse
  -> valid; rust-coverage job parsed
```

The active-path scan below produced no matches:

```text
rg -n \
  'tests/(fixtures|corpus|golden|known_failures)|crates/cosmolkit-core/tests/fixtures|python/tests/fixtures' \
  AGENTS.md README.md dev/README.md testdata tools .github \
  crates/*/src crates/*/tests python/tests
```

The set of `_generate_*.py` helpers under `tools/testdata/rdkit/` exactly
matches the set declared by `GENERATOR_SPECS`; no helper is an undocumented
standalone generation workflow. The temporary migration section has been
removed from the policy. Layout, generated-data identity, cache behavior, and
test-consumer migration are closed. Historical provenance is not fully closed:
the three ligand fixtures listed above remain explicit unresolved provenance
records.
