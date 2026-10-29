# Development Tools

This directory contains repository-maintenance tools that are not production
runtime dependencies.

## Rust release versions

From the repository root, using Python 3.11 or newer:

```sh
python3 dev/tools/bump_rust_version.py 0.6.0-dev.1 --dry-run
python3 dev/tools/bump_rust_version.py 0.6.0-dev.1
```

Updates the shared Rust crate version, explicitly listed internal dependencies,
and the marked installation example in `crates/cosmolkit/README.md`.
Also updates the Python and WASM binding crate package versions and
project.version in `python/pyproject.toml`.
Accepted forms are `X.Y.Z`, `X.Y.Z-rc.N` and `X.Y.Z-dev.N`; the numeric
suffix has no primality restriction. Python uses `X.Y.ZrcN` and `X.Y.Z.devN`.
Dev and RC packages can be published manually or from matching tags;
stable publication requires a tag push. npm/Conda use `dev` and `rc` labels
for prereleases, and `latest`/`main` for stable releases. GitHub dev/RC
releases are marked prereleases and do not replace Latest.
Then `cargo update --workspace` refreshes both workspace lockfiles
without upgrading already-locked third-party dependencies. Historical records and
unmarked version examples remain unchanged.

The script lists each manifest, section and dependency name in `DEPENDENCIES`.
It does not scan Markdown for version numbers. New replacement locations must
be added explicitly; missing fields or README markers cause an error.
`--dry-run` only prints changes. Cargo failure leaves version edits in place and
returns an error. No build, publish, commit or push is performed.

Validate the updater without modifying release files:

```sh
python3 -m unittest discover -s dev/tools -p test_bump_rust_version.py -v
```

## Other tools

- [`inchi/`](./inchi/): official InChI source inventories, call-graph audit,
  and owned Rust source-type generation.
- [`chembl_parity/`](./chembl_parity/): checksummed corpus preparation and
  resumable, profile-driven ChEMBL 37 differential parity audits.
- [`benchmark_pattern_fingerprint.py`](./benchmark_pattern_fingerprint.py):
  deterministic fresh-process Pattern fingerprint timing, cache-reuse, memory,
  and exact-output comparison against the pinned RDKit build. Machine results
  are written under the ignored `tmp/` tree.
- [`debug/`](./debug/): narrow diagnostic probes that are not test or
  production entrypoints.

Historical plan-generation tools live under [`../archive/tools/`](../archive/tools/).
