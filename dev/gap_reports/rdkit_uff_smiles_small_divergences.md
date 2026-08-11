# RDKit UFF smiles_small Divergence Probe

## Scope

This probe covers the two UFF failures remaining after exact golden-float
deserialization and complete MMFF parity on `smiles_small`:

| Row | Surface | First reported mismatch |
|---:|---|---|
| 34 | multi-conformer, conformer 1 | atom 0 axis 0: `0.36058882743683085` versus `0.36050643340176663` |
| 47 | single-conformer | atom 0 axis 0: `1.3370117772556906` versus `1.3606108238464545` |

The probe uses pinned RDKit 2026.03.1, UFF, `vdwThresh=10.0`, ignored
inter-fragment interactions, the generated initial coordinates, and the
declared `1e-6` parity tolerances. No tolerance, fixture row, expected value,
or force-field behavior was changed.

## Probe Log

### Probe 1

- Checked: topology and UFF parameter availability, initial energy, initial
  gradient, single-conformer optimization, and multi-conformer optimization.
- Last matching boundary: both rows have matching topology/parameter
  availability; their initial energy and gradient pass the declared `1e-6`
  checks.
- First differing boundary: optimized coordinates for row 34 conformer 1 and
  row 47 conformer 0.
- Reduced range: ordered contribution construction through minimization.
- Next probe: compare every default torsion SMARTS match, including order and
  orientation, against the current UFF bond-table scan.

### Probe 2

- Checked: complete eligible central-bond membership, order, orientation, and
  the pinned default query's `maxMatches=1000` behavior.
- Last matching boundary: the eligible central-bond set is identical.
- First differing boundary: the ordered and oriented central-bond list passed
  into UFF torsion contribution construction.
- Reduced range: pinned UFF `SubstructMatch` call versus COSMolKit's direct
  bond-table scan in `UFF::Tools::addTorsions`.
- Next probe: none; this contribution boundary precedes all energy, gradient,
  line-search, and minimizer states, so source bisection stops here.

## Source-Level First Divergence

Pinned RDKit UFF `Builder.cpp:499-585` uses the same fixed symmetric query as
MMFF:

```cpp
std::vector<MatchVectType> matchVect;
const ROMol *defaultQuery = DefaultTorsionBondSmarts::query();
unsigned int nHits = SubstructMatch(mol, *query, matchVect);
for (unsigned int i = 0; i < nHits; i++) {
  MatchVectType match = matchVect[i];
  int idx1 = match[0].second;
  int idx2 = match[1].second;
```

COSMolKit UFF instead filters `mol.bonds()` and reuses each bond's stored
begin/end orientation. The resulting rows are:

```text
row 47 RDKit:
(0,1), (0,4), (1,2), (2,3), (3,4)

row 47 current UFF:
(0,1), (1,2), (2,3), (3,4), (4,0)
```

```text
row 34 RDKit:
(1,2), (2,3), (3,4), (3,17), (4,5), (5,6), (6,7), (6,16),
(7,8), (7,15), (8,9), (9,10), (10,11), (10,14), (11,12),
(14,15), (16,17)

row 34 current UFF:
(1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8), (8,9),
(9,10), (10,11), (11,12), (10,14), (14,15), (6,16), (16,17),
(17,3), (15,7)
```

Thus both failures cross a ring-closure ordering/orientation boundary before
UFF torsion terms are appended. The same fixed-query source model was already
audited against pinned RDKit on all 5,000 force-field rows: 120,528 matches
with zero membership, order, orientation, or count mismatches. A global pair
sort is not source-equivalent on high-degree atoms.

## Required Fix

Port the complete fixed default-query match behavior into UFF:

1. traverse first atoms by atom index and second atoms by adjacency order;
2. apply the complete `!$(*#*)&!D1` predicate to both endpoints;
3. preserve accepted query orientation and reject the reverse symmetric atom
   set;
4. stop at the source default of 1000 matches; and
5. pass the ordered endpoint pair and bond index directly into
   `addTorsions`.

Custom SMARTS must remain explicitly unsupported. The fix must add smallest
ring-closure, high-degree adjacency-order, filtering, and match-limit tests,
plus row-34 and row-47 public optimizer regressions. It must not sort pairs,
special-case either molecule, alter UFF equations, or weaken `1e-6`.

## Post-Match-Order Probe

The first public row-34/row-47 regression after the central-bond match-order
port showed that the ordered central matches were necessary but not sufficient:

```text
row 34 conformer 1 atom 0 axis 0
after central-match fix: 0.36073603016092776
RDKit:                  0.36050643340176663
```

The next probe compared every complete torsion tuple, including all four atom
indices and all three bond indices. Row 34 matched pinned RDKit on all 17
central bonds and all 30 constructed torsions. Row 47 matched on the three
torsions whose central atoms have UFF parameters. The two remaining Se-centered
query hits are skipped by both implementations: pinned RDKit 2026.03.1 reports
`UFFTYPER: Unrecognized atom type: Se2+2 (3)` and `has_all=false`, so this is
expected partial-parameter behavior rather than an unported Se torsion.

A pinned-wheel C++ probe then compared the complete row-34 force field against
COSMolKit at the contribution boundary. Both sides constructed exactly 189
contributions in the same type ranges:

```text
0..18    BondStretchContrib
19..42   AngleBendContrib
43..152  vdWContrib
153..182 TorsionAngleContrib
183..188 InversionContrib
```

All 189 individual initial energies matched bit-for-bit. Cumulative initial
gradients matched through contribution 20. The first differing value appeared
after contribution 21, the third `AngleBendContrib`, at flat coordinate 7:

```text
RDKit bits:     c02a6ee2c7313cc5
COSMolKit bits: c02a6ee2c7313cc4
```

The source-level cause is in `UFF::Utils::calcAngleBendGrad`. Pinned RDKit
evaluates each component as:

```cpp
dE_dTheta * dCos_dS[i] / (-sinTheta)
```

COSMolKit precomputed `dE_dTheta / (-sinTheta)` once and then multiplied by
`dCos_dS[i]`. These expressions are algebraically equivalent but have a
different binary64 evaluation order. The complete source-backed fix is to
preserve the C++ multiply-then-divide order independently for all nine gradient
updates. No equation, parameter, tolerance, row, or expected value should
change.

## Post-Angle-Gradient Result

After reproducing the nine `calcAngleBendGrad` multiply-then-divide expressions,
the complete row-34 regression passes parameter coverage, initial energy,
complete initial gradient, single-conformer optimization, and both
multi-conformer optimizations at the declared `1e-6` tolerances.

Row 47 also passes its RDKit `has_all=false` parameter-coverage result, initial
energy, and complete initial gradient. Its single-conformer optimization still
diverges at the first reported coordinate:

```text
row 47 atom 0 axis 0
COSMolKit: 1.3370117772556993
RDKit:     1.3606108238464545
```

The earlier complete contribution probe covered row 34 only. Row 47 therefore
requires its own complete ordered-contribution and minimizer-state probe; the
unsupported Se atom makes the force field underconstrained, so passing a
`1e-6` initial-gradient comparison is not sufficient evidence of source-level
evaluation parity.

## Row-47 Complete Contribution And Minimizer Probe

A pinned-wheel C++ probe and a temporary Rust unit probe compared the complete
fresh aromatic row-47 force fields. Both constructed exactly eight
contributions in the same ranges:

```text
0..2 BondStretchContrib
3..4 AngleBendContrib
5..7 TorsionAngleContrib
```

All eight individual energies and every component of the cumulative gradient
after each contribution matched bit-for-bit. With `snapshotFreq=1`, all 13
accepted BFGS states also matched bit-for-bit in energy and all 15 coordinate
components. Both fresh-input implementations finish with atom-0 x coordinate
`1.3370117772556993` and energy `8.341918182756241e-13`.

The golden generator's `1.3606108238464545` result was reproduced only after
calling `MMFFGetMoleculeProperties` on the shared RDKit `parity_mol` before the
UFF single-conformer helper. For unsupported Se, that call returns `None` but
mutates the molecule from aromatic `c1cc[se]c1` to kekulized
`C1=C[Se]C=C1`. The recorded CXSMILES remains aromatic while the later UFF
golden is therefore computed from a different chemical state.

This is a golden-generation state-isolation defect, not a COSMolKit UFF or
minimizer defect. Each declared force-field surface must receive its own
`Chem.Mol` copy before invoking RDKit APIs that can mutate aromaticity or other
chemical state. The `smiles_small` force-field golden must then be regenerated;
the row-47 public regression must continue to compare every field at `1e-6`.

## Remaining Post-Isolation Failures

After regenerating the complete `smiles_small` force-field family through
`generate_all.py`, the row-47 regression passes. Fourteen of the sixteen
complete force-field parity surfaces pass, while two independently executable
UFF regressions retain coordinate differences above the unchanged `1e-6`
tolerance:

```text
row 81, multi-conformer conformer 0, atom 0 axis 0:
COSMolKit: 1.0514981452102172
RDKit:     1.0515026395495046
absolute difference: 4.494339287397286e-6

row 113, single-conformer, atom 15 axis 0:
COSMolKit: -2.488080635736201
RDKit:     -2.488081666444538
absolute difference: 1.0307083369029966e-6
```

Both focused regressions first compare per-row UFF parameter coverage, initial
energy, and every initial-gradient component, then execute single- and
multi-conformer optimization. Neither row, comparison surface, field, or
tolerance is filtered. Complete contribution and accepted-minimizer-state
instrumentation is required to locate the first source-level divergence.

## Rows 81 And 113 Complete Contribution Probe

A source-built C++ probe loaded the pinned RDKit 2026.03.1 wheel libraries and
compared every ordered contribution energy, every cumulative gradient
component, and every accepted BFGS snapshot against a temporary Rust probe.
Both sides used the exact generated initial-coordinate binary64 values,
`vdwThresh=100.0`, ignored inter-fragment interactions, `maxIters=200`,
`forceTol=1e-4`, and `energyTol=1e-6`.

The force-field shapes and minimizer termination states match:

```text
row 81:  26 points, 394 contributions, result 0, 186 snapshots
row 113: 34 points, 644 contributions, result 1, 200 snapshots
```

Every earlier contribution matches bit-for-bit. The first differences are both
ordinary order-zero `UFF::AngleBendContrib` values:

```text
row 81, global contribution 29, atoms (1, 2, 3):
RDKit energy bits:     3ff435cc85077680
COSMolKit energy bits: 3ff435cc8507765a
first cumulative-gradient difference at flat index 3:
RDKit:                 c03072748829c752
COSMolKit:             c03072748829c748

row 113, global contribution 84, atoms (24, 33, 20):
RDKit energy bits:     3f874af898fbcede
COSMolKit energy bits: 3f874af898fbbe38
first cumulative-gradient difference at flat index 60:
RDKit:                 403c5792b78c06dd
COSMolKit:             403c5792b78c06d9
```

The first accepted row-81 snapshot already differs by one energy bit. The
first accepted row-113 snapshot has equal energy but differs by one coordinate
bit at flat index 62. This proves that the final coordinate failures propagate
from force-field construction/evaluation and do not begin in BFGS.

### Source-Level First Divergence

Pinned RDKit `ForceField/UFF/Params.cpp:65-67` parses the degree value and then
evaluates the conversion left-to-right:

```cpp
paramObj.theta0 = boost::lexical_cast<double>(*token);
paramObj.theta0 = paramObj.theta0 * M_PI / 180.;
```

COSMolKit's runtime parser and build-time default-table generator instead emit
the algebraically equivalent but binary64-distinct expression:

```rust
theta0_degrees * DEG2RAD
```

where `DEG2RAD` was precomputed as `PI / 180.0`. For the shared O_3 value
`104.51`, the complete constructor probes show:

```text
RDKit (104.51 * PI) / 180.0:     3ffd2f4857de2fb3
COSMolKit 104.51 * (PI / 180.0): 3ffd2f4857de2fb4
```

That one-bit parameter divergence changes `C0`, `C1`, and `C2` for both
contributions. Row 81 also changes its angle force constant by one bit; row 113
has the same force-constant bits but still differs through the coefficients.
The complete source-backed fix is to preserve RDKit's multiply-then-divide
evaluation order in both default-table code generation sites and in the custom
parameter parser. It must apply to every UFF row, not to either molecule or
only to `104.51`.

## ChEMBL Row 151 Distance-Cache Probe

After the four deterministic ChEMBL MMFF coverage shards passed, the expanded
`smiles_small` force-field family exposed one additional UFF case containing
phosphorus and two directly attached chlorines. The single-conformer final
energy differed while the multi-conformer conformer-0 final energy differed by
a smaller amount. Both public tests retained the declared `1e-6` comparisons.

A source-built pinned-RDKit probe and a temporary Rust probe compared the full
force-field state for this row. Both implementations constructed exactly 621
ordered contributions. Every initial contribution energy matched bit-for-bit,
every cumulative initial-gradient component matched bit-for-bit, and all 96
optimized coordinate components matched bit-for-bit. The first divergence was
the final `calcEnergy()` evaluation after optimization.

### Source-Level First Divergence

RDKit contributions retain a mutable `ForceField *dp_forceField` and call the
non-const overload of `ForceField::distance()`. That overload writes the
computed distance into the force field's distance cache and reuses the cached
value on subsequent contribution evaluation. COSMolKit's corresponding UFF
and MMFF contributions called `distance_const()`, which recomputed the distance
from the supplied coordinate slice instead of reproducing the source object's
cache semantics. This difference is independent of phosphorus and applies to
every affected bond, angle, stretch-bend, and nonbonded contribution.

The source-level fix replaces `distance_const()` with the mutable cached
`distance()` call in the corresponding UFF bond-stretch, angle-bend, and
nonbonded contributions and MMFF bond-stretch, angle-bend, stretch-bend, and
nonbonded contributions. No molecule-specific branch or tolerance change was
introduced. A smallest-boundary UFF van der Waals regression now proves that a
second energy evaluation through the same source force field reuses the first
cached distance, matching RDKit object semantics.

After the fix, both the complete single-conformer and complete multi-conformer
UFF parity surfaces pass for every eligible `smiles_small` row, including the
curated ChEMBL row 151 case.
