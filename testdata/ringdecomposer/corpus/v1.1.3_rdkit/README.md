# RingDecomposerLib v1.1.3_rdkit regression corpus

These files are byte-preserved test inputs from RareyLab RingDecomposerLib.

- Upstream repository: https://github.com/rareylab/RingDecomposerLib
- Pinned tag: `v1.1.3_rdkit`
- Pinned commit: `7b1629781cfb7fda29716d1af14a6110bb553892`
- Fixture source path: `test/molecule_00_0.dimacs` through `test/molecule_25_0.dimacs`
- License source path: `LICENSE`
- Format source path: `test/README`
- Retrieved from the corresponding `raw.githubusercontent.com` paths.

`LICENSE` is the upstream BSD New license. `UPSTREAM_TEST_README` is the
upstream format and corpus description. The DIMACS fixtures are consumed only
by `crates/cosmolkit-ringdecomposer/tests/migration_ring_families.rs`; no
production parser is introduced.
