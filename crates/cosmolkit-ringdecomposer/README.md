# cosmolkit-ringdecomposer

`cosmolkit-ringdecomposer` is a pure Rust port of the RingDecomposerLib graph
algorithms used by RDKit for Unique Ring Families and relevant-cycle
decomposition.

The source-alignment reference is:

- Upstream repository: <https://github.com/rareylab/RingDecomposerLib>
- Tag: `v1.1.3_rdkit`
- Commit: `7b1629781cfb7fda29716d1af14a6110bb553892`
- Local source path: `third_party/RingDecomposerLib/src/RingDecomposerLib`

The crate is independent from `cosmolkit-core`: users provide plain graph nodes
and undirected edges. COSMolKit adapts molecular topology to this graph layer in
its own core crate.

## Porting Markers

RDL source-reproduction comments use the same two-axis marker policy as the
COSMolKit RDKit porting protocol:

- first marker: behavioral reproduction status;
- second marker: performance and algorithmic-complexity status.

RDL markers use `RDL` instead of `RDKit`, for example:

```text
RDL❌❌: not implemented
RDL❗✔️: behavior is structurally reproduced but not fully proven; local
         inspection supports comparable performance and complexity
RDL✔️✔️: behavior, performance, and complexity are considered equivalent after
         local inspection
```

Do not expand behavior without first copying the corresponding C function body
inside a `BEGIN RDL C FUNCTION ...` / `END RDL C FUNCTION ...` frame.
