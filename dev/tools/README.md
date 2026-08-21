# Development Tools

This directory contains repository-maintenance tools that are not production
runtime dependencies.

- [`inchi/`](./inchi/): official InChI source inventories, call-graph audit,
  and owned Rust source-type generation.
- [`chembl_parity/`](./chembl_parity/): checksummed corpus preparation and
  resumable, profile-driven ChEMBL 37 differential parity audits.
- [`debug/`](./debug/): narrow diagnostic probes that are not test or
  production entrypoints.

Historical plan-generation tools live under [`../archive/tools/`](../archive/tools/).
