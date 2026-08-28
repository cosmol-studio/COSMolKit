# RDKit High-Feasibility Descriptor Source Inventory

## Completed Boundary

- Reference: RDKit `2026.03.1`, revision
  `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- Source root: `third_party/rdkit/`.
- COSMolKit implementation root: `crates/cosmolkit-core/`.
- Included families: connectivity/Chi, Hall-Kier/Kappa/Phi, selected Lipinski
  and ring/stereo counts, MQN, LabuteASA, SlogP-VSA, and SMR-VSA.
- Excluded families: EState, fragment descriptors, graph-complexity
  descriptors, Gasteiger/PEOE-VSA, custom-property VSA, and 3D descriptors.

Every included source function is implemented, source-marked, exposed through
one Rust calculation path, and covered by focused and corpus parity evidence.
This file records ownership and reuse boundaries; the numeric evidence is in
`rdkit_high_feasibility_descriptors_full_port_validation.md`.

## Shared Core Ownership

| Source dependency | Sole COSMolKit owner | Completed disposition |
|---|---|---|
| `findAllPathsOfLengthN` | `chemistry/subgraph.rs` | One neutral source-order atom-path kernel shared by fingerprints and connectivity descriptors. |
| `PeriodicTable::getNouterElecs` and `getRb0` | `chemistry/valence.rs` | One authoritative periodic-table implementation, also used by hydrogen placement. |
| `Atom::getTotalNumHs` | topology/valence cache accessors | Descriptors reuse initialized source-backed total-H state; no family-local hydrogen rules. |
| SSSR and ring membership | `chemistry/rings.rs` plus `lipinski::descriptor_ring_info` | One selected-family acquisition path borrows initialized SSSR-or-better state and performs the source-required cold computation once. |
| SMARTS parsing and matching | `search/smarts_parse.rs` and `search/substruct.rs` | Fixed source SMARTS use retained parsed queries and the shared matcher. |
| Stereo assignment | `chemistry/stereo.rs` | Stereo-center counts reuse typed `stereo_valid` state and the source conditional-clone assignment path. |
| Crippen atom contributions | `properties/descriptors.rs::rdkit_crippen_atom_contribs` | One ordered first-match atom-typing core with retained parsed parameters and a typed contribution cache. |
| Connectivity descriptors | `properties/descriptors/connectivity.rs` | One core owns HK deltas, nVals, alpha, Chi, Kappa, and Phi. |
| Lipinski/ring/stereo descriptors | `properties/descriptors/lipinski.rs` | One core owns the added direct counts and ring/stereo projections while reusing established HBA/HBD, aromatic-ring, fraction-Csp3, and rotor implementations. |
| MQN | `properties/descriptors/mqn.rs` | One fixed `[u32; 42]` core preserves source index and accumulation order. |
| Labute and VSA | `properties/descriptors/mol_surface.rs` | One Labute contribution engine and one VSA `upper_bound` binning loop serve scalar, vector, contribution, and custom-bin projections. |
| Computed descriptor state | `model/molecule.rs::ComputedPropertyCache` | Typed topology-invalidated HK-delta, nVal, Crippen, and Labute entries preserve cold, warm, forced, clone, and parallel-read behavior. |

Architecture tests reject duplicate owners, descriptor-local chemistry,
alternate facade paths, and Python extraction that would discard observational
cache writes.

## Connectivity Ledger

| RDKit source function | Source anchor | COSMolKit implementation |
|---|---|---|
| `detail::hkDeltas` | `ConnectivityDescriptors.cpp:21` | Complete in `connectivity::hk_deltas`, including atomic-number branches, total-H use, reciprocal square roots, typed cache, and `force`. |
| `detail::nVals` | `ConnectivityDescriptors.cpp:50` | Complete in `connectivity::n_vals`, including outer-electron arithmetic, total H, typed cache, and `force`. |
| `detail::getAlpha` | `ConnectivityDescriptors.cpp:71` | Complete in the connectivity core with every atomic-number/hybridization branch and Rb0 fallback path. |
| `calcChiNv`, `calcChiNn` | `ConnectivityDescriptors.cpp:169,188` | Complete through the shared source-order atom-path kernel, including closed-path final-factor omission. |
| `calcChi0v` through `calcChi4v` | `ConnectivityDescriptors.cpp:209-234` | Complete as ordered fixed projections over the valence core. |
| `calcChi0n` through `calcChi4n` | `ConnectivityDescriptors.cpp:238-263` | Complete as ordered fixed projections over the nVal core. |
| Python `Chi0`, `Chi1` | `GraphDescriptors.py:220,237` | Complete graph-degree and bond-degree projections with source filtering and summation order. |
| `calcHallKierAlpha` | `ConnectivityDescriptors.cpp:267` | Complete scalar and atom-contribution forms, including dummy handling and authoritative Rb0 lookup. |
| `kappa1Helper`, `kappa2Helper`, `kappa3Helper` | `ConnectivityDescriptors.cpp:297-324` | Complete private helpers with source polynomial and odd/even behavior. |
| `calcKappa1`, `calcKappa2`, `calcKappa3`, `calcPhi` | `ConnectivityDescriptors.cpp:327-355` | Complete public projections through the shared alpha and path cores, including empty-heavy-atom behavior. |

## Lipinski, Ring, And Stereo Ledger

| RDKit source function | Completed disposition |
|---|---|
| `calcLipinskiHBA`, `calcLipinskiHBD` | Complete as the distinct direct N/O count and N/O total-H sum; they do not replace standard SMARTS HBA/HBD. |
| `calcNumHBA`, `calcNumHBD` | Existing source-backed SMARTS implementations remain the sole standard acceptor/donor cores. |
| `calcNumHeteroatoms`, `calcNumAmideBonds` | Complete through retained source SMARTS and the shared matcher. |
| `calcNumHeavyAtoms`, `calcNumAtoms` | Complete graph projections, including source implicit-H semantics for total atoms. |
| `calcNumRings`, heterocycle/carbocycle, aromatic/aliphatic, and saturated ring families | Complete source-order projections through `descriptor_ring_info`. |
| `calcNumAromaticRings` | Existing source-backed implementation retained and routed through the same ring-state provider. |
| `calcNumSpiroAtoms`, `calcNumBridgeheadAtoms` | Complete with source SSSR intersections, count criteria, and first-observation order. |
| `hasStereoAssigned`, `numAtomStereoCenters`, `numUnspecifiedAtomStereoCenters` | Complete through typed stereo state and conditional source-compatible assignment. |
| `calcFractionCSP3`, `calcNumRotatableBonds` | Existing source-backed implementations retained as the sole cores; MQN reuses the source-default rotor route. |

## MQN Ledger

`calcMQNs` from `MQN.cpp:20` is complete as one `[u32; 42]` implementation.
It preserves the full source index order and accumulation order for element,
bond, polarity, charge, degree, ring-size, multi-ring, and rotatable-bond
entries. Its direct N/O donor/acceptor rules intentionally remain distinct
from the standard SMARTS HBA/HBD descriptors. The Python wrapper converts this
same fixed result and does not introduce scalar MQN implementations.

## Labute And VSA Ledger

| RDKit source function | Completed disposition |
|---|---|
| `getLabuteAtomContribs` | Complete with authoritative Rb0 values, source bond scales, aromatic branch, clamp, implicit-H aggregate, threshold, expression order, and shared-key cache lifecycle. |
| `calcLabuteASA` | Complete cache-first scalar projection through the sole contribution engine. |
| Python `_LabuteHelper` | Complete hydrogen-first projection exposed as a typed Rust result and a Python tuple. |
| `assignContribsToBins` | Complete as the sole binning loop with C++ `upper_bound` equality semantics and source accumulation order. |
| `calcSlogP_VSA` | Complete with the 11 source boundaries, Labute contributions, and cached Crippen logP atom contributions. |
| `calcSMR_VSA` | Complete with the 9 source boundaries, Labute contributions, and cached Crippen MR atom contributions. |
| Python `SlogP_VSA1..12`, `SMR_VSA1..10` | Thin indexes into the corresponding vector core; no duplicate calculations. |

RDKit's Labute computed-property key is shared across `includeHs` choices.
COSMolKit preserves that first-call cache behavior rather than introducing a
second, more intuitive cache key.

## API And Validation Disposition

- Public Rust functions use project-native `snake_case` names and return
  `DescriptorResult` where a source dependency can report an error.
- Fixed families return `[u32; 42]`, `[f64; 12]`, and `[f64; 10]`; custom VSA
  boundaries return `Vec<f64>` through the same binning core.
- Contribution APIs use typed Rust results while preserving the source wrapper
  projection order.
- Python functions borrow the retained `Molecule` with `PyRef`, so computed
  cache writes survive mixed and repeated descriptor calls without mutating
  molecular topology or properties.
- Read-only descriptor calls update only typed observational computed caches
  and therefore do not enter the molecule operation registry.
- The focused behavioral matrix, maintained 5,000-row corpus, and complete
  ChEMBL 37 run compare scalar values, vectors, atom contributions, custom
  bins, and cache sequences exactly against the pinned source with zero
  mismatch.

No selected function remains missing, approximate, duplicated, or blocked.
