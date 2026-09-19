# musl source provenance (no executable oracle)

The C oracle, ABI adapters and musl differential test have been retired.
This directory retains only source provenance and checksums, not a test runner.
Historical commands in audit reports describe past runs and are not executable
instructions for the current tree.

The implementation reference remains the unmodified `third_party/musl`
submodule: musl `v1.2.5`, commit
`0784374d561435f7c787a555aeab8ede699ed298`, from
https://git.musl-libc.org/git/musl. `source_manifest.json` and `SHA256SUMS`
identify the selected source files; licensing is retained in the submodule's
`COPYRIGHT` and `crates/cosmolkit-io/THIRD_PARTY_NOTICES.md`.

Behavior acceptance uses the pinned RDKit environment, not musl equivalence:

```sh
cargo test -p cosmolkit-io --release --features cosmolkit-core/op-contracts-strict --lib v3k_atom_numbers_coordinate_contract -- --include-ignored --nocapture
```

The ordinary-coordinate contract and exclusions remain documented on
`parse_rdkit_atof_prefix`. Existing Rust regressions and discovered
counterexamples remain; removing the musl oracle does not establish all-input
RDKit parity or resolve the excluded extreme-input differences.
