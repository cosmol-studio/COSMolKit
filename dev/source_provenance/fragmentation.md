# Fragmentation source review

RDKit: `fece8caa860bdf6c9c82bdb9253f22a0cf87150c` (2026.09.1).
Regular r-BRICS: `03938c2fd0f16abfcb7c197d4544a23c86a4b530`.
Scope: selected bonds with dummy attachments and the default BRICS, regular
r-BRICS and RECAP decomposition options. This does not implement reBRICS or
the other configurable decomposition modes.

| Boundary | Source | Review |
|---|---|---|
| BRICS vocabulary | `Chem/BRICS.py::environs`, `reactionDefs` | All 19 environments and 46 ordered rules match |
| Regular r-BRICS vocabulary | `rBRICS_public.py::environs`, `reactionDefs` | All 34 environments and 135 ordered rules match |
| RECAP templates | `Chem/Recap.py::reactionDefs` | All 12 ordered templates match |
| Site matching | `FindBRICSBonds`, `FindrBRICSBonds` | Environment prefilter, default match limits, rule order and first-match bond deduplication retained |
| Recursive BRICS | `BRICSDecompose` | FIFO work pool and reaction-group passes retained; stable product-size ordering and accepted-sequence errors retained |
| Recursive regular r-BRICS | `rBRICSDecompose` | Same traversal; a sanitation failure rejects its whole sequence as upstream requires |
| RECAP graph | `RecapDecompose`, `GetLeaves` | Shared canonical products keep all direct parents; forbidden sibling rejects the sequence; root excluded from leaves |
| Product sanitation | `SanitizeMol` | Reuses core sanitation; private reactant provenance is not part of returned string/graph results |
| Selected cuts | `MolFragmenter.cpp::fragmentOnBonds` | Preserves bond type, query, direction, stereo references, incident-order parity and opposite-atom coordinates |
| Isotope width | `Atom.h::setIsotope` | Corrected to unsigned uint16 narrowing, including default atom-index labels; removed the PR-only overflow error |
| Batch deletion | `RWMol.cpp::batchRemoveBonds` | Indexed incident rows and one compaction replace repeated whole-table scans and shifts; source-ordered effects and partial errors retained |
| OPS integration | Existing strong operation lifecycle | Keep required mappings, effects, COW and runtime finalization unchanged |

The RECAP leaf scan avoids repeatedly traversing shared paths. Its result is
the same set of childless reachable products. Interning strings avoids copying
each complete molecule into graph edges. Static compiled rule/reaction tables
are reused. Neither improvement changes chemistry or specified ordering.

Existing reaction, matcher and sanitizer dependencies retain their own source
status; this review does not promote those implementations or claim that every
end-to-end call is as fast as RDKit. Corpus comparison follows source review,
never supplies implementation rules.

## Cost evidence

The initial review overstated the asymptotic difference from RDKit. Its
`ROMol.cpp::getBondWithIdx` explicitly walks the edge iterator to the requested
index; `fragmentOnBonds` and `batchRemoveBonds` call it repeatedly. Therefore
the complete upstream operation is not generally O(V+E+k). CK's former
per-deletion scans and Vec shifts were nevertheless avoidable.

An isolated release benchmark constructs a carbon chain before timing, cuts
every other bond, and checks output atom/bond counts. Median of seven runs,
milliseconds, same machine; no coordinates or SGroups:

| Atoms | Cuts | Original PR | Updated CK | RDKit 2026.09.1 |
|---:|---:|---:|---:|---:|
| 1,000 | 500 | 6.255 | 1.356 | 6.263 |
| 4,000 | 2,000 | 80.789 | 5.780 | 120.367 |
| 16,000 | 8,000 | 2,059.348 | 26.232 | 2,969.061 |

These are selected-bond operation timings, not general BRICS/RECAP throughput
claims. CK uses the same optimized build settings on both sides; the RDKit
column uses the pinned installed wheel and includes Python argument conversion.
Temporary benchmark sources and raw measurements remain outside the checkout.

The final index uses flat offsets and physical row positions, not a per-atom
HashMap. If all stereo-reference lists are empty, their deletion effects are
provably empty: no index or neighbor traversal is needed. Malformed detached
endpoints retain the original scan and error order, without new rejection.

The initial HashMap attempt regressed small molecules and was replaced.
Additional median timings in microseconds per call (seven samples, repeated
calls within each sample):

| Atoms | Cuts | Original PR | Final CK |
|---:|---:|---:|---:|
| 20 | 1 | 4.707 | 4.138 |
| 20 | 10 | 12.540 | 11.606 |
| 100 | 1 | 15.799 | 15.545 |
| 100 | 50 | 70.685 | 49.954 |
| 1,000 | 1 | 148.242 | 142.760 |

Small differences are subject to timing noise; these observations do not prove
a speedup on every input. The removal of repeated whole-edge shifts/scans is
the source-level improvement, with no changed chemistry or OPS permissions.
