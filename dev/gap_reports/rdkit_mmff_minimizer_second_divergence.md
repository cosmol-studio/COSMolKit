# RDKit MMFF Minimizer Second-Divergence Probe

## Scope

This probe starts after the MMFF charge and 1-4 neighbor-matrix fixes. It uses
the pinned RDKit 2026.03.1 wheel and the same `smiles_small` CXSMILES,
MMFF94 parameters, nonbonded threshold `100.0`, inter-fragment policy, and
`200`-iteration limit as the force-field parity suite.

The audited rows are:

| Row | Reason |
|---:|---|
| 76 | first remaining angle-regression final-energy mismatch |
| 103 | first remaining focused single-conformer coordinate mismatch |
| 123 | first remaining focused multi-conformer coordinate mismatch |

No tolerance, field selection, molecule filtering, or optimizer behavior was
changed by the probe.

## Initial Boundary

Initial energy and the declared `1e-6` initial-gradient comparisons pass on all
three rows. A bit/ULP comparison nevertheless shows that the unscaled gradient
is not bit-identical:

| Row | Dimension | Largest initial-gradient difference | ULPs |
|---:|---:|---:|---:|
| 76 | 105 | `4.440892098500626e-15` at component 15 | 40 |
| 103 | 78 | `6.661338147750939e-15` at component 34 | 15 |
| 123 | 60 | `7.105427357601002e-15` at component 30 | 64 |

The force field contains one aggregate contribution for each enabled MMFF
term family. Comparing each family independently against RDKit gives:

| Contribution | Row 76 | Row 103 | Row 123 |
|---|---|---|---|
| bond stretch | bit-identical | bit-identical | bit-identical |
| angle bend | differs | differs | differs |
| stretch-bend | differs | differs | differs |
| out-of-plane bend | bit-identical | bit-identical | bit-identical |
| torsion | differs | differs | differs |
| nonbonded | bit-identical | bit-identical | bit-identical |

This excludes charge construction, the neighbor matrix, bond gradients, OOP
gradients, and nonbonded gradients as the first remaining divergence.

## Iteration Boundary

The pinned RDKit force field was restarted from the same CXSMILES coordinates
with `maxIts=1`, `2`, and `3`. Those results were compared to snapshots from
COSMolKit's source-matched BFGS implementation.

| Row | First non-bit-identical accepted state | Evidence |
|---:|---:|---|
| 76 | iteration 1 | accepted coordinates differ by single ULPs; energy differs by 2 ULPs |
| 103 | iteration 3 | iterations 1 and 2 have identical coordinate and energy bit hashes |
| 123 | iteration 1 | one near-zero coordinate already differs by 40 ULPs; energy is initially identical |

The line-search acceptance and result codes agree at these checkpoints. The
first divergence is therefore upstream of BFGS control flow: the contribution
gradients feed slightly different directions into the otherwise matching
minimizer. Repeated Hessian updates and line searches amplify those differences
until the final `1e-6` coordinate/energy assertions fail.

## Source-Level First Divergence

The three nonmatching contribution families contain algebraically equivalent
but source-inconsistent floating-point rewrites.

### Angle bend

Pinned RDKit evaluates every component as:

```cpp
g[0][0] += dE_dTheta * dCos_dS[0] / (-sinTheta);
```

COSMolKit currently precomputes `dE_dTheta / (-sinTheta)` and then multiplies
by `dCos_dS`. This changes the operation and rounding order.

### Stretch-bend

Pinned RDKit evaluates the six cosine derivatives as:

```cpp
double dCos_dS1 = 1.0 / dist1 * (p32.x - cosTheta * p12.x);
```

COSMolKit currently evaluates `(p32.x - cosTheta * p12.x) / dist1`. The same
rewrite exists in all six components.

### Torsion

Pinned RDKit evaluates the six torsion cosine derivatives as:

```cpp
double dCos_dT[6] = {1.0 / d[0] * (t[1].x - cosPhi * t[0].x), ...};
```

COSMolKit currently evaluates each numerator divided by `d[0]` or `d[1]`.
Again, the algebra is equivalent but the IEEE-754 rounding path is not.

## Required Fix

The source-backed fix is limited to restoring the pinned RDKit evaluation
order in the complete affected expression sets:

1. all nine angle-bend gradient writes, without a precombined scale;
2. all six stretch-bend cosine derivatives; and
3. all six torsion cosine derivatives.

The fix must retain the copied C++ lines and two-axis markers, add
smallest-boundary evaluation-order regressions, and rerun the row 76/103/123
optimizer regressions plus the complete MMFF optimizer parity tests. It must
not introduce molecule-specific branches, altered tolerances, or result
filtering.
