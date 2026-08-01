# RDKit Built-In Fixtures

This directory is a committed, byte-for-byte test copy of files from the
pinned RDKit `2026.03.1` source tree at revision
`351f8f378f8ad6bbd517980c38896e66bf907af8c`.

The relative `Code/...` and `Data/...` layout is retained so the fixture
migration runner can exercise the same file families without requiring the
RDKit submodule at test time. `source_manifest.jsonl` records the exact
upstream path, byte length, and SHA-256 of every copied file. Selection is the
complete set consumed by COSMolKit's RDKit built-in fixture migration suite;
no case is filtered at test time.

The copied upstream files retain RDKit's BSD-3-Clause license context. The
license text is available from the pinned RDKit release as `license.txt`.
