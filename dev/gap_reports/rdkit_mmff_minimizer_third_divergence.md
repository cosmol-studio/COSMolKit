# RDKit MMFF Minimizer Third-Divergence Probe

## Scope

This probe starts after restoring the pinned RDKit evaluation order in the
angle-bend, stretch-bend, and torsion-gradient expressions. It uses pinned
RDKit 2026.03.1 and the row-103 `smiles_small` CXSMILES with MMFF94,
`nonBondedThresh=100.0`, ignored inter-fragment interactions, and the normal
`1e-4` force and `1e-6` energy tolerances.

No tolerance, expected coordinate, molecule selection, or production
force-field behavior was changed by the probe.

## Reproduced Failure

The focused public optimizer regression still fails at row 103. The first
reported final-coordinate difference is atom 1, axis 0:

| COSMolKit | RDKit | Absolute difference |
|---:|---:|---:|
| `-0.6358963962626648` | `-0.6361468191287529` | `2.504228660881402e-4` |

The initial energy is bit-identical (`0x404f0a5e85fd0fee`). The unscaled
initial gradient now differs in only three of 78 components:

| Component | COSMolKit bits | RDKit bits | ULP distance |
|---:|---:|---:|---:|
| 32 | `13847405499541408552` | `13847405499541408554` | 2 |
| 64 | `13822838106945573246` | `13822838106945573238` | 8 |
| 65 | `4623066711331896392` | `4623066711331896394` | 2 |

## Contribution Boundary

Each MMFF family was enabled independently at the same initial coordinates.
Bond stretch, angle bend, stretch-bend, out-of-plane bend, and combined
nonbonded gradients are bit-identical. Only the torsion gradient differs.

The first accepted BFGS state remains bit-identical after gradient scaling and
line search:

| Accepted state | COSMolKit energy bits | RDKit energy bits | Coordinate result |
|---:|---:|---:|---|
| 1 | `0x404b7e22a6974bf4` | `0x404b7e22a6974bf4` | bit-identical |
| 2 | `0x4048fd550eecee13` | `0x4048fd550eecee12` | different bit hash |

The later BFGS mismatch is therefore an amplification of a torsion-gradient
input difference, not a new line-search branch or tolerance difference.

## Source-Level First Divergence

Pinned `Builder.cpp:589-704` obtains central torsion bonds through:

```cpp
std::vector<MatchVectType> matchVect;
const ROMol *defaultQuery = DefaultTorsionBondSmarts::query();
unsigned int nHits = SubstructMatch(mol, *query, matchVect);
for (unsigned int i = 0; i < nHits; ++i) {
  MatchVectType match = matchVect[i];
  int idx2 = match[0].second;
  int idx3 = match[1].second;
```

COSMolKit instead scans the molecule bond table and returns bond-table indices.
That preserves the eligible bond set but not the source match orientation or
iteration order. Row 103 demonstrates the difference:

```text
RDKit prefix:
(2,4), (4,5), (5,7), (7,8), (7,21), (8,9), ...

COSMolKit prefix:
(2,4), (4,5), (5,7), (7,8), (8,9), (9,10), ...
```

RDKit also normalizes later ring-closure matches through the symmetric query;
COSMolKit retains bond-table orientations such as `(21,7)`.

The ordering is defined by the pinned substructure source, not by a global
lexicographic sort:

- `SubstructMatchParameters.uniquify` defaults to `true` and
  `maxMatches` defaults to `1000` in `SubstructMatch.h`.
- `SubstructMatch.cpp` runs `boost::vf2_all`, writes pairs into query-atom
  index order, and rejects a repeated target-atom set in the final checker.
- For this fixed symmetric two-atom query, query atom 0 advances in target
  atom-index order, query atom 1 follows that target atom's adjacency order,
  and the reverse mapping is removed by atom-set uniquification.

A source-specialized model of exactly those rules was compared against pinned
RDKit on all 5,000 rows of `smiles_5000`: 120,528 returned matches, zero order,
orientation, membership, or count mismatches. A plain globally sorted pair
model was rejected because it mismatched high-degree adjacency ordering.

## Required Fix

Port the complete modeled behavior of the fixed default torsion SMARTS match
boundary:

1. iterate eligible first atoms in molecule atom-index order;
2. iterate their bonds in molecule adjacency order;
3. apply the complete `!$(*#*)&!D1` predicate to both endpoints;
4. preserve query-atom orientation and remove the symmetric reverse match by
   the same target-atom-set uniqueness rule;
5. enforce the source default `maxMatches=1000` limit;
6. pass the ordered endpoint pairs directly into `addTorsions`.

Custom torsion SMARTS must remain explicitly unsupported. The fix must retain
the pinned C++ source blocks and two-axis markers, add smallest-boundary tests
for ring-closure orientation, high-degree adjacency ordering, filtering, and
the 1000-match limit, and rerun the row 76/103/123 plus complete MMFF optimizer
parity tests.
