# RDKit MMFF Optimizer First-Divergence Probe

## Scope

This probe compares pinned RDKit 2026.03.1 and COSMolKit at the MMFF atom
property, force-field contribution, initial-gradient, minimizer-result, and
termination boundaries. Both sides use the checked-in `smiles_small` CXSMILES
coordinates, MMFF94, a nonbonded threshold of 100.0, inter-fragment
interactions disabled, and the normal 200-iteration minimization limit.

The audited failing rows are 19, 21, 103, and 123 of
`testdata/forcefield/expected/rdkit/smiles_small/forcefield_params.jsonl`.

## First Divergence

The first divergence is before force-field construction and before the
minimizer. `MMFFMolProperties::new` assigns the same audited atom types as
RDKit, but leaves every formal and partial charge at zero because the source
call to `MMFFMolProperties::computeMMFFCharges` is still marked `RDKit❌❌`.

For row 19, RDKit assigns `(atom type, formal charge, partial charge)` as:

```text
[(1, 0.0, 0.0), (1, 0.0, 0.28), (6, 0.0, -0.28),
 (1, 0.0, 0.28), (1, 0.0, 0.0), (6, 0.0, -0.28)]
```

COSMolKit assigns:

```text
[(1, 0.0, 0.0), (1, 0.0, 0.0), (6, 0.0, 0.0),
 (1, 0.0, 0.0), (1, 0.0, 0.0), (6, 0.0, 0.0)]
```

`builder.rs` only constructs an electrostatic pair contribution when both
partial charges are nonzero. The all-zero property state therefore removes
the complete electrostatic contribution rather than exposing an optimizer
algorithm difference.

After restoring charge calculation, a second source divergence became
observable: the MMFF builder reused UFF's neighbor matrix. The pinned UFF
implementation only distinguishes 1-2 and 1-3 pairs, while pinned MMFF builds
a topological distance matrix and additionally marks distance-three pairs as
1-4. Treating those pairs as ordinary 1-X pairs omitted MMFF's `0.75`
electrostatic scaling. The source fix therefore also requires MMFF's own
`Tools::buildNeighborMatrix`, not cross-model reuse.

## Contribution Isolation

Pinned RDKit term toggles produced the following initial energies. In every
row, the observed COSMolKit energy is exactly the sum of the six
non-electrostatic RDKit terms at displayed precision.

| Row | Bond | Angle | Stretch-bend | OOP | Torsion | vdW | Electrostatic | RDKit all | COSMolKit |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 19 | 0.798553858813 | 0.931637719005 | -0.232433413990 | 0 | 1.528499438693 | 2.041670920298 | 5.373040665038 | 10.440969187858 | 5.067928522820 |
| 21 | 2.179936264571 | 5.705018214769 | -1.008141433602 | 0 | 1.283461564315 | 3.740484329678 | -28.531711953750 | -16.630953014019 | 11.900758939731 |
| 103 | 7.859338994953 | 15.744229051948 | -1.301270788942 | 0.000007808377 | 26.937641471461 | 43.796340945488 | -30.955277858265 | 62.081009625021 | 93.036287483286 |
| 123 | 31.786492812572 | 201.071739401839 | -8.235828066150 | 11.171571356136 | 20.071224873936 | 59.480560766427 | -28.071150078494 | 287.274611066266 | 315.345761144760 |

This also explains why failures previously attributed to angle or torsion rows
had differences numerically equal to the absent electrostatic energy.

## Gradient Boundary

For row 19, flat gradient index 6 isolates the same issue:

```text
RDKit all terms:                  9.399621650533480
RDKit electrostatic term:        0.683871485286421
RDKit with electrostatic off:    8.715750165247059
COSMolKit:                       8.715750165247059
```

The non-electrostatic gradient is therefore already identical at this audited
boundary. Divergence begins when charge assignment suppresses electrostatic
term construction.

## Minimizer And Termination Boundary

For row 19, a pinned-RDKit force field with only the electrostatic term disabled
matches the current COSMolKit result through minimization:

| Boundary | RDKit all terms | RDKit electrostatic off | COSMolKit |
|---|---:|---:|---:|
| Initial energy | 10.440969187858046 | 5.067928522820338 | 5.067928522820338 |
| Initial gradient index 6 | 9.399621650533480 | 8.715750165247059 | 8.715750165247059 |
| Termination status after 200 iterations | 0 | 0 | 0 |
| Final energy | 8.586072606774610 | 3.296244900268137 | 3.296244900268137 |

The optimized coordinates from COSMolKit and electrostatic-disabled RDKit
agree to approximately `1.4e-15` maximum absolute difference. Thus neither the
iteration path nor the termination condition is the first divergence for this
failure class.

## Required Source Fix

The source-backed implementation must port the complete pinned-RDKit
`MMFFMolProperties::computeMMFFCharges` function, including:

- every formal-charge atom-type switch branch;
- delocalized ring and conjugated-nitrogen charge sharing;
- terminal O/S neighbor counting and precedence;
- negative-neighbor and atom-type-62 formal-charge adjustments;
- MMFF bond-charge-increment lookup and PBCI fallback; and
- MMFF.V equation 15 partial-charge calculation.

It must also port MMFF's byte-per-cell topological relationship matrix and
retain the distinct `RELATION_1_4` state used by electrostatic contributions.

The existing `MmffPropCollection`, `MmffPbciCollection`, and
`MmffChgCollection` supply the required parameter data. No optimizer tolerance,
minimizer rule, molecule-specific condition, or heuristic correction is
justified by this probe.
