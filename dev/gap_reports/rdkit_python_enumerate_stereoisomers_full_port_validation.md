# RDKit Python EnumerateStereoisomers Full-Port Validation

## Validated Boundary

This report closes the source-reproduction audit for the pinned RDKit
`2026.03.1` Python stereoisomer-enumeration path. The public parity boundary is:

```text
Chem.FindPotentialStereo
GetStereoisomerCount
EnumerateStereoisomers
_AtomFlipper / _BondFlipper / _StereoGroupFlipper
_RangeBitsGenerator / _UniqueRandomBitsGenerator
tryEmbedding
```

The newer C++ `StereoisomerEnumerator` is a different API and is not part of
this claim. The Python path does not enumerate atropisomers; already represented
atropisomer state is preserved where the Python source preserves it. Distance
geometry is reused through COSMolKit's existing source-backed implementation
instead of being duplicated inside the enumerator.

## Source-Line And Ownership Closure

The complete source inventory is fixed in
[`rdkit_python_enumerate_stereoisomers_source_audit.md`](./rdkit_python_enumerate_stereoisomers_source_audit.md).
The two production owners have no incomplete behavioral or complexity markers:

| Owner | Complete source markers | Incomplete markers |
|---|---:|---:|
| `chemistry/potential_stereo.rs` | 1,144 `RDKit✔️✔️` | 0 |
| `chemistry/stereo_enumerate.rs` | 196 `RDKit✔️✔️`; 88 `CPython✔️✔️` | 0 |

The shared SMILES writer still contains incomplete markers for writer behavior
outside this plan. This report therefore does not claim completion of the whole
writer. The cleanup, assignment, direction, and canonical-isomeric-SMILES
subset actually called by enumeration is covered through exact focused,
small-corpus, 5,000-row, and ChEMBL full-state observations.

Architecture guards prove that production has one potential-stereo workspace
owner, one typed flipper model, one configuration source, and one lazy iterator.
The following retired paths are absent from production: independent atom and
bond candidate detectors, local double-bond feasibility checks, eager Boolean
combination builders, `EnumerationStrategy`, `EnumerationParams`, xorshift
sampling, retry caps, and the former 20-center limit. Rust and Python public
APIs are value-style delegates over this single implementation.

## Algorithm And Allocation Closure

Potential-stereo perception keeps source-record order and applies canonical
reranking until the source fixed point is stable. Its mutable work is held in
private typed vectors and molecule state; public source molecules are not
mutated.

Enumeration preserves the source complexity and allocation shape:

```text
exhaustive configurations       lazy arbitrary-width BigUint counter
random configurations           finite unique HashSet exhaustion
theoretical count               arbitrary-width BigUint
output uniqueness               canonical isomeric SMILES HashSet
one emitted configuration       one private workspace/output clone
configuration universe          never materialized as 2^n vectors
custom random-source failure    structured error at the lazy iterator boundary
```

There is no center-count heuristic, unbounded duplicate retry loop, or
molecule-specific feasibility filter. Optional embedding composes the existing
hydrogen and distance-geometry implementations using the exact Python source
seed, one-conformer attempt, failure filtering, heavy-atom coordinate copy, and
successful-output counting semantics.

## Operation And Invariant Compliance

Stereoisomer enumeration is a multi-result value derivation. It clones the
source into a private workspace and never exposes an in-place public topology
mutation, so a new `MoleculeOpSpec` would misclassify the API. Public hydrogen
and coordinate operations used by composition continue through their existing
registered operation boundaries. Every emitted molecule passes the project's
output invariant enforcement, and source preservation is exercised by focused,
small, 5,000-row, ChEMBL, repeated-call, and parallel tests.

## Exact Parity Evidence

All committed expected data is generated from the pinned RDKit and CPython
identities and is protected by repository manifests. Comparisons include
ordered potential-stereo records, cleaned molecule state, theoretical counts,
bounded-state flags, emitted order, canonical isomeric SMILES, atom chiral
tags, bond directions/stereo/controllers, conformer counts, source
preservation, and structured parse outcomes.

| Evidence set | Rows | Profiles and branches | Exact result |
|---|---:|---|---:|
| Focused source fixtures | 23 | all four potential-stereo cleanup/possible combinations plus enumeration, random, group, uniqueness, embedding, and error fixtures | zero mismatch |
| `smiles_small` | 152 | 2 potential-stereo profiles and 4 bounded enumeration profiles | 89,047 comparisons; zero mismatch |
| `smiles_5000` | 5,000 | 2 potential-stereo profiles and 4 bounded enumeration profiles | 4,680,490 comparisons; zero mismatch |
| ChEMBL 37 | 2,897,819 | 4 potential-stereo profiles and 4 bounded enumeration profiles | 89,831,939 comparisons; zero mismatch |

The ChEMBL phase completed all `128/128` deterministic shards, compared every
one of 2,897,819 input rows, fully compared 2,897,804 rows accepted by both
parsers. The remaining 15 rows were rejected by both parsers at parse entry. It recorded no blocking,
informational, or other mismatch and retained no mismatch example. Its exact
comparison total is one parse-acceptance observation per row plus 30
state/output/source observations per row accepted by both parsers.

Focused fixture coverage is reported as cases and exercised branches rather
than added to the corpus leaf totals: its tests use branch-specific structural
assertions, not the corpus JSON leaf-count accounting rule. The two committed
corpora contribute 4,769,537 leaf comparisons; the complete ChEMBL phase adds
89,831,939 phase observations. No row filter, field omission, tolerance,
alternative accepted output, or known-failure allowance is used.

## Benchmark Evidence

The deterministic benchmark runs each engine in a fresh process, verifies
exact output before timing, and records seven rounds at intact public stage
boundaries. Median COSMolKit/RDKit ratios are:

| Profile | Ratio |
|---|---:|
| Candidate discovery | 0.283 |
| Configuration count and setup | 0.362 |
| One-configuration finalization | 0.299 |
| Lazy prefix, 5 outputs | 0.854 |
| Full exhaustive, 16 outputs | 0.909 |
| Bounded random, 16 outputs | 1.015 |
| Embedding, 8 outputs | 1.283 |
| Parallel, 32 jobs | 0.614 |

No enumerator-specific complexity or allocation regression is present. The
measured embedding difference belongs to the already shared distance-geometry
composition path rather than a duplicate or asymptotically weaker enumeration
implementation. The machine-local benchmark and ChEMBL aggregates remain
temporary validation artifacts and are not repository fixtures.

## Validation Conclusion

The pinned Python `EnumerateStereoisomers` behavioral boundary is closed for
the modeled public state space. Source lines, option defaults, potential-stereo
fixed point, typed flippers, arbitrary-width counts, exhaustive and random
configuration order, finalization, uniqueness, embedding composition, lazy
errors, source preservation, scalar/parallel use, and complete ChEMBL 37
coverage have zero unexplained mismatch. The validation does not claim the
separate C++ enumerator, atropisomer enumeration, or unrelated incomplete
SMILES-writer behavior.
