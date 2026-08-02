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

## Current Port Status

The Rust implementation contains no `RDL❌❌` source lines: every inventoried
source frame has an implementation. Behavioral proof is not yet closed for the
entire crate. At the current revision, 447 copied source lines are marked
`RDL✔️✔️` and 286 remain `RDL❗✔️` because their behavior is implemented but not
fully proven across all source-defined edge states.

The crate backs COSMolKit's stable ring-perception surface for its documented
supported graph states, with focused regression and RDKit parity coverage for
those states. The term "port" here does not claim blanket exact parity for the
remaining `RDL❗✔️` source frames.

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
