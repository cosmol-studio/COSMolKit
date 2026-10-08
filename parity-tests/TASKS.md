# Parity task definitions

COSMolKit 0.5.0 uses the executable registry in
[`parity-tests_fixed/src/registry.rs`](../parity-tests_fixed/src/registry.rs).
The [runner README](../parity-tests_fixed/README.md) defines corpus selection,
task filters, reference preparation and Cargo comparison.

Tasks compare one registered behavior on one declared input type. Parameters,
reference versions and comparison rules belong to executable definitions,
not a separately maintained Markdown task-count table. Registration is not
passing evidence; complete 0.5.0 validation is pending.

The registry in this legacy package records its original implementation only.
It must not be used to infer current missing functions or create a second queue.
