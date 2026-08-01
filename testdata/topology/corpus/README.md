# Topology Invariant Corpora

These CSV files are repository-authored operation-contract inputs. `core.csv`
contains the shared topology, chemistry-state, coordinate, property, and
stereo cases. `cow_small.csv` is the smaller copy-on-write value-semantics
matrix.

Tests consume every named case in committed order. No external reference
dataset is claimed, and no row is filtered at test time. The files use the
repository's license context. `source_manifest.jsonl` records record counts,
byte lengths, and SHA-256 values.
