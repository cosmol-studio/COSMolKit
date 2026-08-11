# RDKit MMFF Minimizer Fourth-Divergence Probe

## Scope

This probe starts after the default-torsion SMARTS match order was aligned
with pinned RDKit 2026.03.1. It investigates the remaining row-123
multi-conformer coordinate mismatch using conformer 0, MMFF94,
`nonBondedThresh=100.0`, ignored inter-fragment interactions, and the normal
`1e-4` force and `1e-6` energy tolerances.

No tolerance, expected value, molecule filter, force-field term, or optimizer
behavior was changed by the probe. Temporary Rust and Python probes were
removed after collecting the results.

## Reproduced Failure

The public row-123 regression first reports this final-coordinate mismatch:

| Component | COSMolKit | RDKit | Absolute difference |
|---|---:|---:|---:|
| atom 0, axis 0 | `3.3744424137236` | `3.3744173734462777` | `2.504027732216362e-5` |

The mismatch is specific to the multi-conformer fixture path. The preceding
single-conformer initial energy, initial gradient, final energy, and final
coordinates pass their declared parity assertions.

## Probe Log

### Probe 1

- Checked: row-123 conformer-0 initial energy, complete initial gradient,
  `maxIts=0/1/2/3/200` energy bits, and coordinate bits.
- Last matching boundary: molecule topology, conformer count, and decimal
  coordinate values at the fixture boundary.
- First differing boundary: the initial coordinate bit patterns; initial
  energy already differs by 2 ULP (`0x40722552c214d628` in Rust versus
  `0x40722552c214d626` in RDKit), before BFGS runs.
- Reduced range: multi-conformer golden coordinate transport through molecule
  construction, upstream of MMFF contribution construction and minimization.
- Next probe: compare the raw JSON coordinate numbers with each runtime's
  deserialized values and conformer storage.

### Probe 2

- Checked: raw JSON number text, deserialized `f64` bits, and conformer-stored
  coordinate bits on both sides.
- Last matching boundary: the raw decimal literals in
  `forcefield_params.jsonl`.
- First differing boundary: Python `json` and Rust `serde_json` deserialize
  several literals to adjacent `f64` values. Each conformer API preserves the
  value it receives without further modification.
- Reduced range: JSON floating-point deserialization, not RDKit or COSMolKit
  MMFF code.
- Next probe: compare Python parsing, Rust standard-library parsing, and the
  active `serde_json` parsing implementation for representative literals.

### Probe 3

- Checked: Python `json.loads`, Rust `str::parse::<f64>()`, and
  `serde_json::from_str::<f64>()` with the workspace's active serde_json
  features.
- Last matching boundary: Python parsing and Rust standard-library parsing are
  bit-identical for every inspected literal.
- First differing boundary: serde_json 1.0.151 without the
  `float_roundtrip` feature.
- Reduced range: serde_json's non-`float_roundtrip` `f64_from_parts()` branch
  in `src/de.rs`.
- Next probe: none; the first divergent source block is identified.

## Source-Level First Divergence

The Python generator serializes RDKit coordinates as JSON numbers. Python
parses those shortest round-tripping decimal strings back to the original
binary64 values. Rust's standard parser produces the same values:

| JSON literal | Python / Rust standard bits | serde_json default bits |
|---|---:|---:|
| `-1.3006687716310001` | `0xbff4cf8a0ed15695` | `0xbff4cf8a0ed15694` |
| `-0.39018908733716867` | `0xbfd8f8dba657a16b` | `0xbfd8f8dba657a16a` |
| `1.3682138422726557` | `0x3ff5e43432a7edcf` | `0x3ff5e43432a7edce` |
| `1.8158026507898606` | `0x3ffd0d8714921ef7` | `0x3ffd0d8714921ef8` |

In serde_json 1.0.151, `float_roundtrip` is not a default feature. The
non-feature `f64_from_parts()` implementation converts the decimal
significand to `f64` and then multiplies or divides by a power of ten. The
feature-enabled implementation instead calls lexical's concise float parser,
which preserves round trips.

The row-123 Rust molecule therefore starts from coordinates adjacent to the
coordinates used for the pinned RDKit result. MMFF minimization amplifies
those input ULP differences. This is why changing MMFF expressions or BFGS
control flow would be an incorrect fix.

## Required Fix

Enable serde_json's `float_roundtrip` feature in the repository test-support
dependency so parity-test binaries deserialize generated binary64 values
exactly. Add a smallest-boundary regression for representative literals and
rerun row 123, the charge and 1-4 set, the angle set, and the complete MMFF
optimizer tests.

The fix must not change MMFF production behavior, regenerate or hand-edit the
golden coordinates, relax `1e-6`, or special-case row 123.
