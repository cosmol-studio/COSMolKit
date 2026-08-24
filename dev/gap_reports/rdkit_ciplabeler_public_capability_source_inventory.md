# RDKit CIPLabeler Public Capability Source Inventory

## Purpose

This inventory records the source boundary and final closure audit for the
maintained modern RDKit CIP labeler capability. The function ledger preserves
the initial audit route; the closure sections record how every identified gap
was resolved and which adjacent source capabilities remain deliberately
outside the supported dispatcher boundary.

The intended architecture has one modern implementation core:

```text
public Rust/Python assignment operations
                    |
                    v
chemistry::ciplabeler (modern CIPLabeler port)
                    |
          +---------+---------+
          |                   |
          v                   v
shared stereo helpers   molecule property/state APIs
```

Legacy `Chirality.cpp` rank/code helpers remain separate only where a pinned
legacy RDKit call path requires them. They must never implement, wrap, or
silently substitute for the modern public assignment operation.

## Normative Pin

The reference manifest is `testdata/reference/rdkit.json`:

| Field | Value |
|---|---|
| RDKit release | `2026.03.1` |
| Source revision | `351f8f378f8ad6bbd517980c38896e66bf907af8c` |
| Python distribution | `2026.3.1` |
| Python runtime | CPython `3.13.12` |

The source files listed below are relative to
`Code/GraphMol/CIPLabeler/` at that revision. A local extracted source tree is
useful for inspection, but the committed manifest remains the source pin.

## Current COSMolKit State

`crates/cosmolkit-core/src/chemistry/ciplabeler.rs` is the sole crate-private,
source-anchored modern implementation. It owns the digraph, sequence rules,
configurations, assignment dispatcher, selection masks, recursion budget, and
observable source-state transitions. The registered public Rust operation,
Python binding, and fingerprint consumers all route to this core; no second
modern CIP engine exists.

The current public `chemistry::stereo` module also contains legacy behavior:

- `assign_atom_cip_ranks` and `assign_atom_cip_ranks_in_place` reproduce the
  legacy `_CIPRank` path used by older RDKit stereochemistry consumers.
- `assign_legacy_cip_labels` reproduces a legacy `_CIPCode` path.

Those functions are not equivalent to modern
`CIPLabeler::assignCIPLabels`. The reachability audit confirms that the modern
public operation does not call them. They remain only for pinned legacy
consumers that require their distinct behavior.

## Public Source Boundary

Pinned `CIPLabeler.cpp::findConfigs()` constructs configurations only for:

- atoms tagged `CHI_TETRAHEDRAL_CW` or `CHI_TETRAHEDRAL_CCW`;
- bonds tagged `STEREOE`, `STEREOZ`, `STEREOTRANS`, or `STEREOCIS`;
- bonds tagged `STEREOATROPCCW` or `STEREOATROPCW`.

The modern descriptor enum can represent `SP_4`, `TBPY_5`, and `OC_6`, but
the pinned public assignment dispatcher does not construct corresponding
configuration objects. Enum representability is not assignment support.
Focused oracle evidence must lock this boundary instead of claiming broad
non-tetrahedral assignment.

`rules/Rule5.cpp` is present in the source tree but is not part of
`all_rules`; the pinned assignment path uses `Rule5New`. The old Rule5 is
source history unless a separate reachable pinned caller is demonstrated.

The following adjacent stereochemistry capabilities are outside this plan
unless an exact dependency is proved during source audit:

- `findPotentialStereo`;
- full `assignStereochemistry` orchestration;
- `assignStereochemistryFrom3D`;
- 3D stereo perception and coordinate-to-stereo assignment.

## Source-To-Rust Function Ledger

The third column preserves the pre-port audit finding that routed execution.
`present`, `boundary`, and `consolidate` therefore describe the starting point,
not an unresolved current status. Every finding is closed in the final audit
near the end of this document.

### Descriptor, Priority, Sort, And Edge

| RDKit source function or member | Current Rust counterpart | Planning state |
|---|---|---|
| `Descriptor::to_string` | `descriptor_to_string` | present; make public typed conversion derive from this single mapping |
| `Priority::Priority` | `CipPriority::new` | present |
| `CIPLabeler_detail::decrementRemainingCallCountAndCheck` | `CipLabelerContext::decrement_remaining_call_count_and_check` | present; audit zero/underflow and per-call replacement of source thread-local state |
| `Priority::isUnique` | `CipPriority::is_unique` | present |
| `Priority::isPseudoAsymetric` | `CipPriority::is_pseudo_asymetric` | present; preserve source spelling only internally |
| `Sort::Sort` overloads | `CipSort::new`, `from_rules` | present |
| `Sort::getRules` | `CipSort::get_rules` | present |
| `Sort::prioritize` | `CipSort::prioritize` | present |
| `Sort::getGroups` | `CipSort::get_groups` | present |
| `Sort::compare` | `CipSort::compare_substituents` | present |
| `Edge::Edge` | `CipEdge::new` | present |
| `Edge::getOther` | `CipEdge::get_other` | present; exact failure behavior required |
| `Edge::getBeg` / `getEnd` | `get_beg` / `get_end` | present |
| `Edge::getBond` | `get_bond_idx` | present; index-based representation audit required |
| `Edge::getAux` / `setAux` | `get_aux` / `set_aux` | present |
| `Edge::isBeg` / `isEnd` | `is_beg` / `is_end` | present |
| `Edge::flip` | `flip` | present |

### Sequence Rules And Rule Aggregation

| RDKit source function or member | Current Rust counterpart | Planning state |
|---|---|---|
| `SequenceRule::getBondLabel` | `CipSequenceRule::get_bond_label` | present |
| `SequenceRule::getComparision` overloads | `get_comparison`, `get_comparison_with_sort_rules` | present; retain source comparison semantics despite upstream typo |
| `SequenceRule::compare` overloads | `compare`, `compare_with_sort_rules` | present |
| `SequenceRule::getSorter` / `setSorter` | explicit `sort_rules` propagation and `CipSort` construction | present as a structural Rust adaptation; every replacement-rule call path requires equivalence proof |
| `SequenceRule::recursiveCompare` overloads | `recursive_compare`, `recursive_compare_with_sort_rules`, `recursive_compare_sequence_rule` | present; recursion-budget boundary is public-observable |
| `SequenceRule::sort` | `sort`, `sort_edges_for_sequence_rule` | present |
| `SequenceRule::areUpEdges` | `are_up_edges` | present |
| `three_way_comparison` | `three_way_comparison_i32` | present |
| `Rules::Rules` / `add` | `CipRules::new` / `add` | present; null-rule error mapping required |
| `Rules::getNumSubRules` | `get_num_sub_rules` | present |
| `Rules::getSorter` | `get_sorter` | present |
| `Rules::compare` and comparison delegation | `compare`, `get_comparison`, `sort` | present |
| `constitutional_rules` | `cip_constitutional_rules` | present |
| `all_rules` | `cip_all_rules` | present; must contain `Rule5New`, not old `Rule5` |

### Individual Sequence Rules

| RDKit source function or member | Current Rust counterpart | Planning state |
|---|---|---|
| `Rule1a::compare` | `CipRule1a::compare` | present |
| `Rule1b::compare` | `CipRule1b::compare` | present |
| `Rule2::compare` | `CipRule2::compare` | present |
| `Rule3::ord` / `compare` | `CipRule3::ord` / `compare` | present |
| `Rule4a::ord` / `compare` | `CipRule4a::ord` / `compare` | present |
| `PairList::ref` | `CipPairList::ref_descriptor` | present |
| `PairList` constructors | `new`, `with_ref`, `from_head_tail` | present |
| `PairList::getRefDescriptor` | `get_ref_descriptor` | present |
| `PairList::add` / `addAll` | `add` / `add_all` | present |
| `PairList::getPairing` | `get_pairing` | present; fixed-width bit behavior required |
| `PairList::compareTo` | `compare_to` | present |
| `PairList::toString` | `to_rdkit_string` | present |
| `PairList::addAndPair` | `add_and_pair` | present |
| `Rule4b::compare` overloads | `compare`, `compare_with_sort_rules` | present |
| `Rule4b::getReferenceDescriptors` | `get_reference_descriptors` | present |
| `Rule4b::hasDescriptors` | `has_descriptors` | present |
| `Rule4b::getReference` | `get_reference` | present |
| `Rule4b::initialLevel` / `getNextLevel` | `initial_level` / `get_next_level` | present |
| `Rule4b::toNodeList` | `to_node_list` | present |
| `Rule4b::newPairLists` | `new_pair_lists` | present |
| `Rule4b::fillPairs` / `comparePairs` | `fill_pairs` / `compare_pairs` | present |
| `Rule4b::getRefSorter` | `get_ref_sorter` | present |
| `Rule4c::ord` / `compare` | `CipRule4c::ord` / `compare` | present |
| `Rule5New::compare` overloads | `compare`, `compare_with_sort_rules` | present |
| `Rule5New::fillPairs` | `fill_pairs` | present |
| `Rule5New::getRefSorter` | `get_ref_sorter` | present |
| `Rule6::compare` | `CipRule6::compare` | present |
| `Rule5::ord` / `compare` | no assignment-path counterpart required | boundary; prove unreachable from pinned modern dispatcher |

### Node And Digraph

| RDKit source function or member | Current Rust counterpart | Planning state |
|---|---|---|
| `Node::Node` constructors | `CipNode::new`, `new_terminal_child` | present |
| `Node::getDigraph` | `get_digraph` | present |
| `Node::getAtom` / `getAtomIdx` | `atom_idx` / `get_atom_idx` | present; sentinel width audit required |
| `Node::getDistance` | `get_distance` | present |
| `Node::getAtomicNumFraction` | `get_atomic_num_fraction` | present |
| `Node::getAtomicNum` | `get_atomic_num` | present |
| `Node::getMassNum` / `getAtomicMass` | `get_mass_num` / `get_atomic_mass` | present |
| `Node::getAux` / `setAux` | `get_aux` / `set_aux` | present |
| `Node::isSet` / `isDuplicate` / `isDuplicateOrH` | corresponding `is_*` functions | present |
| `Node::isTerminal` / `isExpanded` / `isVisited` | corresponding `is_*` functions | present |
| child and duplicate constructors | `new_child`, `new_bond_duplicate_child`, `new_ring_duplicate_child`, `new_implicit_hydrogen_child` | present |
| `Node::add` | `add` | present |
| `Node::getEdges` overloads | `get_edges`, `get_edges_for_atom` | present |
| `Node::getNonTerminalOutEdges` | `get_non_terminal_out_edges` | present |
| `Digraph::Digraph` | `CipDigraph::new` | present |
| `Digraph::getMol` | `mol` | present |
| `Digraph::getOriginalRoot` / `getCurrentRoot` | corresponding getters | present |
| `Digraph::getNumNodes` / `getNodes` | `get_num_nodes` / `get_nodes` | present |
| `Digraph::addNode` / `addEdge` | `add_node`, `add_existing_node`, `add_edge` | present |
| `Digraph::getEdges` overloads | `node_edges`, `node_edges_for_atom` | present |
| `Digraph::getNonTerminalOutEdges` | `non_terminal_out_edges` | present |
| `Digraph::getRule6Ref` / `setRule6Ref` | corresponding functions | present |
| `Digraph::changeRoot` | `change_root` | present; edge-direction restoration required |
| `Digraph::expand` | `expand` | present; duplicate/ring behavior requires source audit |

### CIPMol And Mancude

| RDKit source function or member | Current Rust counterpart | Planning state |
|---|---|---|
| `CIPMol::CIPMol` | `CipMol::new` | present |
| `CIPMol::getFractionalAtomicNum` | `get_fractional_atomic_num` | present |
| `CIPMol::getNumAtoms` / `getNumBonds` | corresponding getters | present |
| `CIPMol::getAtom` / `atoms` | `atom` / `atoms` | present |
| `CIPMolSpan::CIPMolIter` and span `begin/end` | adjacency-backed atom/bond index vectors | representation adaptation; ordering and invalidation behavior require evidence rather than a duplicate iterator abstraction |
| `CIPMol::getBond` / `getBonds` | `bond` / `bond_indices_for_atom` | present |
| `CIPMol::getNeighbors` | `neighbor_indices` | present |
| `CIPMol::isInRing` | `is_in_ring` | present |
| `CIPMol::getBondOrder` | `get_bond_order` | present; aromatic copy semantics required |
| implicit hydrogen and other-end helpers | `total_num_hs`, `other_atom_idx` | present; source dependency audit required |
| `Mancude::SeedTypes` | `seed_types` | present |
| `Mancude::RelaxTypes` | `relax_types` | present |
| Mancude component traversal | `visit_part`, `visit_parts` | present |
| `Mancude::calcFracAtomNums` | `calc_frac_atom_nums` | present |

### Configuration Base And Tetrahedral

| RDKit source function or member | Current Rust counterpart | Planning state |
|---|---|---|
| `Configuration` constructors | `CipConfiguration::new`, `with_foci` | present |
| `Configuration::getFocus` / `getFoci` | corresponding getters | present |
| `Configuration::setCarriers` / `getCarriers` | corresponding functions | present |
| `Configuration::getDigraph` | `get_digraph` / `get_digraph_mut` | present |
| `Configuration::isInternalEdge` / internal-edge removal | `is_internal_edge`, `find_internal_edge`, `remove_internal_edges` | present |
| duplicate/hydrogen edge filtering | `is_duplicate_or_hydrogen_edge`, `remove_duplicates_and_hs` | present |
| `Configuration::parity4` | `parity4` | present |
| base label/primary-label virtual surface | `CipConfigurationLike` methods | present; dispatch ownership audit required |
| `Tetrahedral::Tetrahedral` | `CipTetrahedral::new` | present |
| Tetrahedral foci/carriers/digraph access | corresponding getters | present |
| `Tetrahedral::setPrimaryLabel` | `set_primary_label` | present; atom property behavior required |
| `Tetrahedral::hasPrimaryLabel` / `resetPrimaryLabel` | corresponding functions | present |
| `Tetrahedral::label` overloads | `label`, `label_with_external_digraph`, `label_node`, `label_node_in_digraph`, `label_node_impl` | present |

### Sp2 And Atropisomer Configurations

| RDKit source function or member | Current Rust counterpart | Planning state |
|---|---|---|
| `Sp2Bond::Sp2Bond` | `CipSp2Bond::new` | present |
| `Chirality::findStereoAtoms` dependency | local constructor block and `find_highest_cip_neighbor_like_rdkit` | present; exact adjacent-source audit required |
| Sp2 foci/carriers/ranked anchors | corresponding getters | present |
| `Sp2Bond::setPrimaryLabel` | `set_primary_label` | present; bond property and stereo-atom behavior required |
| `Sp2Bond::hasPrimaryLabel` / `resetPrimaryLabel` | corresponding functions | present |
| `Sp2Bond::label` overloads | `label`, `label_with_external_digraph`, `label_node`, `label_node_impl` | present |
| `AtropisomerBond::AtropisomerBond` | `CipAtropisomerBond::new` | present |
| `Atropisomers::getAtropisomerAtomsAndBonds` dependency | CIP-local `atropisomer_carriers_like_rdkit` and `chemistry::atropisomer` helper | consolidate into one exact internal helper |
| atropisomer foci/carriers/ranked anchors | corresponding getters | present |
| `AtropisomerBond::setPrimaryLabel` | `set_primary_label` | present |
| `AtropisomerBond::hasPrimaryLabel` / `resetPrimaryLabel` | corresponding functions | present |
| `AtropisomerBond::label` overloads | `label`, `label_with_external_digraph`, `label_with_root`, `label_node_impl` | present |

### Top-Level Assignment

| RDKit source function or member | Current Rust counterpart | Planning state |
|---|---|---|
| `findConfigs` | `cip_find_configs` | present; public supported configuration boundary |
| `labelAux` | `cip_label_aux`, `cip_set_center_node_aux` | present |
| `label` | `cip_label`, `cip_label_with_center_digraph` | present; two-pass budget behavior required |
| selected `assignCIPLabels` overload | `assign_cip_labels_for_indices` | present; crate-private and mask-shaped |
| full `assignCIPLabels` overload | `assign_cip_labels` | present; crate-private |
| selection label clearing/application | `cip_clear_selected_labels`, `cip_apply_primary_labels` | present; source-observable lifecycle audit required |
| `_CIPNeighborOrder` serialization | `cip_neighbor_order_value` | present; currently at risk of using `usize::MAX` instead of `Atom::NOATOM == UINT_MAX` |
| `ControlCHandler::reset/getGotSignal` and `ControlCCaught` | no equivalent project cancellation mechanism | boundary; implement real cancellation or return structured unsupported status without a parity claim |

## Initial Function Audit

The first execution audit compared the complete pinned production source tree
under `Code/GraphMol/CIPLabeler/` (excluding tests and the Python wrapper) with
the Rust source-anchor inventory. This section preserves the initial findings;
all gaps listed here were subsequently revalidated and closed as recorded in
the final closure audit.

| Source unit | Function coverage found in the modern Rust core | Gaps or required proof |
|---|---|---|
| `Descriptor.h` | enum ordering and `to_string` are present | public typed conversion must distinguish descriptors emitted by assignment from enum-only values |
| `Priority.h` | constructor and both accessors are present | inventory mapping was corrected so recursion context is not confused with `Priority` |
| `Edge.{h,cpp}` | constructor and all eight operational methods are present | Rust structured endpoint error must be checked against the source exception boundary |
| `Sort.{h,cpp}` | both constructors, `getRules`, `prioritize`, `getGroups`, and `compareSubstituents` are present | performance marker requires explicit allocation and repeated-comparison review |
| `SequenceRule.{h,cpp}` | constructor/destructor-equivalent ownership, bond labels, comparison overloads, recursive traversal, sorting, and upward-edge handling are present | `getSorter/setSorter` are represented by explicit rule-slice propagation; every Rule4b/Rule5New replacement path must prove the same rule ordering and ownership semantics; cancellation polling is absent |
| `Rules.h` | construction, add, destruction-by-ownership, sorter construction, comparison overloads, and subrule count are present | null rule and dynamic sorter replacement errors require exact tests |
| `Rule1a/1b/2/3` | every constructor-equivalent, ordinal helper, and comparison body is present | IUPAC 2013 constant behavior and mass-table source dependency require oracle coverage |
| `Pairlist.h` | descriptor mapping, all constructors, add/addAll, pairing, comparison, ordering, formatting, and addAndPair are present | pairing storage width and overflow behavior require exact-width review |
| `Rule4a/4b/4c` | all ordinal, traversal, reference, pair construction, pair comparison, replacement sorter, and compare functions are present | replacement-sorter structural adaptation remains unproven outside inline unit examples |
| `Rule5New/Rule6` | all reachable modern assignment functions are present | old `Rule5` is not in `all_rules` and remains deliberately outside the modern dispatcher; Rule6 reference identity requires repeated-root evidence |
| `Node.{h,cpp}` | all constructors, child factories, getters, flags, edge access, lazy expansion, and auxiliary-state methods are present | `Atom::NOATOM` is incorrectly represented by platform-width `usize::MAX`; null-node failure paths use local unsupported errors and need source-boundary review |
| `Digraph.{h,cpp}` | construction, add node/edge, getters, node lookup, Rule6 reference, root change, and expansion are present | vector-plus-stable-ID is a valid representation candidate but complexity markers need review; max-node exception text and root restoration need boundary tests |
| `CIPMol.{h,cpp}` | molecule access, atom/bond iteration substitutes, ring lookup, and cached kekulized bond order are present | the Rust kekulization call catches its full local error result while RDKit catches `MolSanitizeException`; exact failure-class equivalence is unproven; span iteration ordering must be locked |
| `Mancude.{h,cpp}` | type enum, seed, relax, component traversal, and fractional-number calculation are present line-for-line | rational integer overflow and valence/implicit-H source dependencies require corpus and boundary evidence |
| `Configuration.{h,cpp}` | constructors, foci/carriers/digraph access, internal-edge helpers, duplicate filtering, base label, and `parity4` are present | `nullptr` carriers are represented by `usize` sentinels in some paths and need one exact-width model |
| `Tetrahedral.{h,cpp}` | constructor, all label overloads, primary-state functions, parity, carrier construction, and ranked anchors are present | deferred property writes do not currently reproduce repeated-call `resetPrimaryLabel/hasPrimaryLabel` semantics |
| `Sp2Bond.{h,cpp}` | constructor, all label overloads, primary-state functions, stereo writes, and ranked anchors are present | reachable `Chirality::findStereoAtoms/findHighestCIPNeighbor` source dependency must be rechecked; preexisting `_CIPRank` read behavior needs oracle rows |
| `AtropisomerBond.{h,cpp}` | constructor, all label overloads, primary-state functions, parity, and ranked anchors are present | carrier extraction is duplicated locally and in `chemistry::atropisomer`; it must have one source-backed owner |
| `CIPLabeler.cpp` | `findConfigs`, `labelAux`, two-pass `label`, both assignment overloads, selection masks, and recursion budget are present | source interruption behavior is absent; completion-state and repeated-call behavior differ through the deferred-write adapter |
| `TooManyNodesException.h` | mapped to `CipLabelerError::TooManyNodes` | exact message and catch/propagation surface require public error tests |

### Initial Behavioral Gaps

1. `CipTetrahedral::has_primary_label`, `CipSp2Bond::has_primary_label`, and
   `CipAtropisomerBond::has_primary_label` consult preexisting input
   `_CIPCode` in addition to the deferred primary label. Their corresponding
   Rust `reset_primary_label` methods clear only deferred state. On repeated
   assignment, the second pass can therefore skip a selected center because
   the immutable input still has an old code, even though the output adapter
   already cleared that code. The single modern core must model source mutation
   during the algorithm or provide an equivalent working-state layer; merely
   patching the final output is not sufficient.
2. Pinned `resetPrimaryLabel()` clears `_CIPCode` only. The initial adapter
   clears both `_CIPCode` and `_CIPNeighborOrder` before labeling. An unresolved
   selected center therefore has different residual neighbor-order state.
   Full, selected, repeated, and failure cases must use the pinned lifecycle.
3. `_CIPNeighborOrder` is a typed vector of source-width unsigned atom indices
   in RDKit, while COSMolKit currently serializes it into a string property.
   `Atom::NOATOM` is `UINT_MAX` (`4294967295`), not 64-bit `usize::MAX`.
4. `_CIPComputed` is set as a computed property upstream. COSMolKit currently
   stores only the string value `"1"`; it has no demonstrated equivalent
   computed-property lifecycle.
5. `ControlCHandler` polling and interrupted-return behavior remain marked
   `RDKit❌❌`. No public parity claim may include this branch until it has a real
   implementation or an explicit structured unsupported boundary.
6. The embedded `Chirality::findStereoAtoms` closing-brace marker at the
   current Rust line near 4481 is incorrectly marked unreproduced even though
   the enclosing branch is present. The Sp2 audit must determine whether this
   is marker bookkeeping or a behavioral omission before changing it.
7. `Atropisomers::getAtropisomerAtomsAndBonds` has two local reproductions.
   The existing public helper returns only endpoints despite its source name,
   while the CIP-local helper separately determines ordered carrier bonds.
8. The initial replacement of source-owned mutable molecule writes with deferred
   `CipPrimaryLabel` records was not yet proven equivalent for label ordering,
   auxiliary labeling, errors, or partial completion. This adapter must be
   evaluated as part of the core, not treated as harmless facade plumbing.

### Source Items Intentionally Not Reimplemented As Separate Functions

- `CIPMolSpan` iterator mechanics are replaced by the molecule adjacency/index
  representation; source ordering remains part of parity.
- default constructors and destructors with no behavior map to Rust defaults
  and ownership rather than duplicate functions.
- old `Rule5` is not reachable from pinned `all_rules`; `Rule5New` is the sole
  Rule 5 implementation in the modern assignment path.
- pointer identity is represented by stable node/edge IDs. This does not relax
  identity-dependent Rule6, root, or queue semantics.

## Observable State Contract

The pinned modern assignment path can observably affect:

- atom `_CIPCode`;
- bond `_CIPCode`;
- atom `_CIPNeighborOrder`;
- bond `_CIPNeighborOrder`;
- molecule `_CIPComputed` as a computed property;
- Sp2 bond stereo atoms in the source-dependent constructor/setter path;
- recursion-limit failure and interruption completion state.

Modern CIPLabeler does not assign `_CIPRank`. Preexisting `_CIPRank` may still
be read by the Sp2 fallback used when explicit stereo atoms are absent, and its
preservation must be compared. The public capability therefore needs separate
typed descriptor accessors without pretending `_CIPRank` is modern output.

## Modern And Legacy Ownership Audit

The repository has two source-distinct CIP families. They share property names
because RDKit does, but they do not share an assignment algorithm and must not
be merged:

```text
pinned GraphMol/CIPLabeler
  -> chemistry/ciplabeler.rs (sole modern implementation owner)
     -> registered `with_cip_labels_with_options` operation
     -> Morgan bond chirality invariants
     -> Morgan environment chirality
     -> AtomPair chirality atom codes
     -> Topological Torsion chirality through the shared fingerprint path

pinned GraphMol/Chirality.cpp and its direct consumers
  -> chemistry/stereo.rs (sole legacy rank/code implementation owner)
     -> SMILES stereo parse/write and canonical ordering
     -> stereo enumeration
     -> coordinate/stereo orchestration
     -> InChI bridge state
     -> molecular hashes
     -> hydrogen-operation compatibility state
```

The modern fingerprint consumers call the same modern core only when their
pinned source path requires CIP codes and `_CIPComputed` is absent. They consume
the returned value without mutating the caller's molecule. No fingerprint-local
labeling implementation exists.

The `_CIPRank` consumers are retained legacy boundaries. Modern CIPLabeler does
not produce `_CIPRank`; redirecting those call sites to modern assignment would
change pinned SMILES, InChI, depiction, enumeration, and fingerprint behavior.
The one modern exception is the pinned Sp2 fallback that may read a preexisting
`_CIPRank` when stereo atoms are absent. That is an input-state dependency, not
permission for modern assignment to generate legacy ranks.

| Local entry point or state consumer | Boundary decision | Ownership reason |
|---|---|---|
| `chemistry::ciplabeler::assign_cip_labels` | retained behind one registered public operation | sole modern full-molecule core |
| `chemistry::ciplabeler::assign_cip_labels_for_indices` | retained behind the same registered operation | sole modern selected atom/bond core |
| `chemistry::stereo::assign_atom_cip_ranks` | retain as legacy source-backed API | direct legacy rank consumers require `_CIPRank` semantics |
| `chemistry::stereo::assign_atom_chiral_codes` | retain but document as legacy helper | reproduces legacy atom-code assignment, not modern CIPLabeler |
| `chemistry::stereo::assign_bond_stereo_codes` | retain but document as legacy helper | reproduces legacy bond-code assignment, not modern CIPLabeler |
| `chemistry::stereo::assign_legacy_cip_labels` | retained as an explicitly legacy source boundary and is not reachable from the modern operation | legacy behavior is not a substitute modern implementation |
| `chemistry::atropisomer` carrier extraction | shared source-backed helper used by the CIP core | one exact carrier-order implementation owns the behavior |
| generic `_CIPCode` readers in drawing and notation | retain as state observers | they render or serialize assigned state and do not own assignment |
| `_CIPCode` clearing in notation, PDB, and hydrogen operations | lifecycle audit required before changing | these are source-backed state transitions whose interaction with modern computed state must be explicit |

The root facade exposes the project-native typed modern operation and
descriptor queries without reexporting the internal labeler. Legacy functions
remain source-distinct and are not silently reinterpreted.

The implemented lifecycle specifies:

- full versus selected label clearing and preservation;
- when `_CIPComputed` is cleared and set;
- state after recursion-limit errors, interruption, invariant failure, and
  operation failure;
- clone and value-style operation behavior;
- binary and supported text serialization behavior;
- topology/stereo mutation invalidation or remapping;
- interaction with existing computed-property flags and generic property maps.

## Pinned Observable-State Lifecycle Matrix

This matrix records the behavior of the pinned `assignCIPLabels()` path and the
implemented operation/state contract.

| Invocation or transition | `_CIPCode` | `_CIPNeighborOrder` | `_CIPComputed` | `_CIPRank` | bond stereo / stereo atoms |
|---|---|---|---|---|---|
| Full assignment entry | each configuration clears only its own code when reached by the first pass | preserved until a center is successfully relabeled | cleared before `CIPMol` construction | preserved | preserved until a successful Sp2 label write |
| Full assignment success | resolved atom/bond configurations receive the emitted descriptor; unresolved selected configurations have no code | overwritten for each successfully labeled configuration; an unresolved configuration retains any preexisting value because reset does not clear it | set to boolean `true` and marked computed, even when no configurations exist | preserved | a successfully labeled Sp2 bond receives source carrier stereo atoms and normalized CIS/TRANS configuration; atropisomer labeling does not rewrite stereo atoms |
| Selected assignment entry | only selected atoms/bonds that become configurations are reset; unselected codes remain | selected and unselected values remain until successful selected writes | cleared globally | preserved globally | unselected state remains unchanged |
| Selected assignment success | only selected configurations may receive new codes; unselected codes are preserved | only successfully labeled selected configurations are overwritten; all other values are preserved | set globally to computed boolean `true`, matching the pinned selected overload | preserved globally | only successfully labeled selected Sp2 bonds may rewrite stereo atoms/stereo |
| Constitutional-pass recursion exhaustion | the exception is caught per configuration; that configuration remains without a code after reset and later receives a full-rule retry | old neighbor order remains unless the retry resolves and overwrites it | remains absent until the whole call succeeds | preserved | unchanged unless a later retry resolves |
| Full-rule recursion error or other propagated labeling error | source mutation is not transactional: earlier configurations may already have new codes, the failing/current one may have had its old code cleared, and later configurations may be untouched | earlier successful writes remain; reset never clears old neighbor order | remains absent | preserved | earlier successful Sp2 writes remain |
| Interrupt caught by the pinned control handler | partial source writes remain | partial source writes remain | remains absent and the function returns after warning | preserved | partial successful Sp2 writes remain |
| Clone/copy | copied | copied together with its computed-property registration | copied together with its computed-property registration | copied | copied |
| `clearComputedProps()` | preserved because `setPrimaryLabel()` writes it as non-computed | cleared because primary labeling marks it computed | cleared because assignment marks it computed | legacy `_CIPRank` is computed upstream and is cleared | stereo and stereo-atom structural state is preserved by property clearing alone |
| COSMolKit binary roundtrip, required project behavior | must preserve descriptor state | must preserve typed source-width vector state and its computed classification | must preserve boolean value and computed classification | preserve according to the existing legacy archive contract | preserve exactly |
| Text formats without an explicit property channel | no general promise to serialize internal CIP properties; reconstruct only where that format's pinned reader does so | not serialized implicitly | not serialized implicitly | not serialized implicitly | format-defined stereo fields remain authoritative |
| Topology/stereo mutation | operation-specific pinned behavior; never silently treat a stale code as freshly computed | invalidate or remap through the registered operation's source behavior | clear whenever the mutation invalidates completed CIP assignment | retain/clear through the separate legacy source contract | remap, preserve, or clear through the operation's declared stereo behavior |

Three source details control the implementation:

1. `resetPrimaryLabel()` clears `_CIPCode` only. Clearing neighbor order at
   assignment entry is a divergence.
2. `setPrimaryLabel()` writes `_CIPCode` without the computed flag, but writes
   `_CIPNeighborOrder` with the computed flag. The top-level function writes
   `_CIPComputed` as a computed boolean.
3. Failure is source-observably non-transactional. COSMolKit's value-style
   operation may return an error without exposing its discarded working value,
   while the in-place operation must remain structurally valid and may retain
   the same partial mutation allowed by the project operation policy. It must
   not claim rollback parity.

COSMolKit now preserves typed descriptor access, source-width neighbor-order
state, computed-property membership, binary archive classification, and
source-order partial mutation on the in-place error path. The value-style path
discards its failed private working value while leaving the source unchanged.

## Error, Limit, Selection, And Cancellation Audit

| Pinned boundary | Pinned behavior | Required Rust boundary | Required Python boundary |
|---|---|---|---|
| `maxRecursiveIterations == 0` | uses `UINT_MAX`; it means effectively unlimited, not zero recursive comparisons | `u32` option value `0` with source-width wrapping decrement | non-negative integer in the `u32` range; `0` retains pinned meaning |
| constitutional fast pass | resets the budget to 2,000 for every configuration and catches `MaxIterationsExceeded` for that pass only | internal retry behavior, never exposed as an error by itself | same |
| full-rule budget | one source-width counter is shared across the ordered unresolved configuration loop; it is not reset per center | `CipLabelerError::MaxIterationsExceeded` with exact source message `Max Iterations Exceeded in CIP label calculation` | `RuntimeError` with the same message, matching the pinned wrapper translator |
| digraph node cap | expansion at 100,000 nodes throws `TooManyNodesException` with `Digraph generation failed: more than 100000 nodes found.` | a distinct structured error retaining the exact limit and display text | `RuntimeError` with the same display text unless the project defines a narrower public exception class |
| selected atom or bond index at/above count | pinned Python wrapper rejects it before assignment with `list element larger than allowed value`; the C++ core accepts bitsets rather than raw indices | validate typed atom and bond index collections before constructing exact-size masks; return the corresponding structured out-of-range error | `ValueError("list element larger than allowed value")` for parity with the pinned wrapper |
| duplicate selected indexes | dynamic bitset idempotently selects the index once | deduplicate through the exact-size mask without changing order of configuration discovery | accepted and treated as one selection |
| false-valued selections | the pinned wrapper truth-tests both Python objects: any combination of `None` and empty collections selects all when both are false-valued; once either category is non-empty, the other false-valued category is empty | both absent or empty selects all; once either collection is non-empty, an absent or empty category selects no entries | reproduces the same wrapper truth-value dispatch |
| negative or non-integer Python selection | conversion to the unsigned bitset index type fails in wrapper conversion before assignment | not applicable to typed Rust indexes | propagate a normal conversion/type error; do not wrap, clamp, or ignore |
| malformed focus/configuration/stereo carrier state | pinned `CHECK_INVARIANT`, `PRECONDITION`, or runtime error propagates; no fallback configuration is fabricated | retain the existing specific structured variants (`Bad*`, `IncorrectNumberOfStereoAtoms`, `CarrierMismatch`, endpoint/index errors) and expose them through the operation error without panic | map malformed molecule/state variants to `ValueError`; internal contract failures remain `RuntimeError` |
| non-kekulizable Mancude input | pinned `CIPMol` catches `MolSanitizeException` around its kekulized copy and proceeds with the source fallback state | catch only the audited equivalent sanitize class; unrelated internal errors must propagate | follows the resulting Rust structured error or successful fallback |
| cancellation signal | recursive comparison throws `ControlCCaught`; top level catches it, leaves partial writes and `_CIPComputed` absent, logs, and returns; pinned Python wrapper then raises `KeyboardInterrupt("Assign CIP labels cancelled")` | no equivalent project cancellation source exists; cancellation requests return `CipLabelerError::CancellationUnsupported` before mutation | no interrupted-call parity claim; unsupported cancellation is explicit and structured |

Source-width behavior is part of these boundaries. The decrement operation is
unsigned 32-bit pre-decrement: a zero counter wraps to `UINT_MAX`, while a
counter of one becomes zero and throws. Selection masks must always have the
molecule's exact atom/bond lengths; the current crate-private slice API's
short-mask omission and long-mask behavior are not suitable as public
semantics.

`TooManyNodesException`, `MaxIterationsExceeded`, and malformed-state errors
leave the same non-transactional working-state effects described in the
lifecycle matrix. The value-style Rust API can discard that failed working
value, but the in-place operation must preserve structural validity while
documenting that source-order partial CIP state may remain after an error.

## Oracle Schema And Evidence Ownership

The focused oracle is JSON Lines with one independently reproducible case per
row. A row has this logical schema (the generator and schema validator added
later must encode these fields directly rather than flattening away state):

```text
schema_version: 1
case_id: stable unique string
source:
  fixture: repository-relative path or corpus profile plus source row id
  input_kind: smiles | molblock
  input: exact source text
  sanitize: boolean
initial_state: CipStateSnapshot
calls: [
  {
    atoms_to_label: omitted | [u32, ...]
    bonds_to_label: omitted | [u32, ...]
    max_recursive_iterations: u32
    result:
      status: ok | error
      error_type: null or stable pinned exception name
      error_message: null or exact message
    state: CipStateSnapshot
  }, ...
]
```

`null` and an empty list remain distinct in the fixture even though both are
false-valued at the pinned wrapper boundary; this keeps the wrapper dispatch
observable and proves their common all-molecule result when both categories
are false-valued. `CipStateSnapshot` contains:

```text
molecule:
  cip_computed: { present, value, value_type, computed }
atoms: [
  {
    index,
    chiral_tag,
    cip_code: { present, value, computed },
    cip_neighbor_order: { present, value_u32, computed },
    cip_rank: { present, value_u32, computed }
  }, ...
]
bonds: [
  {
    index,
    begin,
    end,
    stereo,
    stereo_atoms_u32,
    cip_code: { present, value, computed },
    cip_neighbor_order: { present, value_u32, computed }
  }, ...
]
```

Every atom and bond row is emitted, including absent values. Lists use source
`unsigned int` values and therefore preserve `Atom::NOATOM == 4294967295`.
Computed membership is obtained from the difference between pinned property
name queries with and without computed properties; it is not inferred from the
property name. The state following an error is recorded because the pinned
operation is non-transactional. Repeated-call cases use multiple entries in
`calls` and never reconstruct the molecule between entries.

Artifact ownership is fixed as follows:

| Artifact | Owner and location | Review rule |
|---|---|---|
| human-authored focused inputs and source citations | `testdata/stereo/fixtures/` and its source manifest | committed; identifiers and pinned upstream provenance are stable |
| generated focused expected state | `testdata/stereo/expected/rdkit/<profile>/ciplabeler.jsonl` | committed; produced only by the pinned generator |
| generator | `tools/testdata/rdkit/_generate_ciplabeler_golden.py` | registered in `generate_all.py`; refuses a non-pinned RDKit identity |
| JSONL schema validation | `tools/testdata/rdkit/_expected_schema.py` | validates every field, exact integer widths, ordering, and no unknown schema version |
| Rust focused parity gate | `crates/cosmolkit-core/tests/rdkit_ciplabeler_parity.rs` | consumes committed expected data; never calls RDKit at test time |
| large-corpus runner and aggregates | `dev/tools/chembl_parity/` | deterministic partitions; report and mismatch artifacts identify source rows |

Maintained evidence tiers are cumulative:

| Tier | Input scope | Required branches and retained output |
|---|---|---|
| Focused | source-derived cases, pinned upstream regressions, and local boundary regressions | full/atom-only/bond-only/empty selection, repeated calls, all emitted descriptor families, unresolved centers, exact-width sentinels, recursion errors, malformed state, and unsupported dispatcher/cancellation boundaries; complete row output committed |
| Small corpus | complete `testdata/smiles/corpus/smiles_small.smi` | full assignment plus deterministic selected profiles and recursion profiles; complete expected JSONL committed |
| 5,000 corpus | every active complete 5,000-row profile, with no convenient-row filtering | same successful observable state plus deterministic selected and bounded-recursion branches; sharded complete expected output committed or represented by the repository's reviewed generated-data manifest policy |
| ChEMBL 37 | complete maintained ChEMBL 37 corpus through deterministic partitions | full and selected calls, configured recursion profiles, complete state comparison in workers; aggregate branch/comparison/error counts and every mismatch row retained in the maintained parity report |

Corpus parse failures and molecules without supported configurations remain
counted rows. They are not dropped: parse outcome and the successful no-config
`_CIPComputed` transition are themselves observable branches. Large-corpus
summaries may aggregate passing rows, but any mismatch artifact must contain
the complete schema row needed to reproduce it locally.

## Public API Shape

The completed project-native surface provides:

- a typed `CipDescriptor` covering the exact descriptor strings emitted by the
  pinned modern path;
- atom and bond descriptor query methods backed by the single stored state;
- registered value-style and trailing-underscore in-place molecule operations;
- options for full assignment, selected atoms, selected bonds, and
  `max_recursive_iterations`;
- Python methods following the same semantics without exposing mutable
  internals or a second Python-only implementation.

Selection semantics must match the pinned wrapper: both false-valued categories
(`None` or empty) mean all; a non-empty atom collection labels only those atoms
when the bond category is false-valued, and conversely for bonds. Atom-local
assignment must still execute the molecular-context CIPLabeler.

## Oracle And Corpus Evidence Boundary

Focused oracle rows must compare more than a final notation string. Each row
must capture:

- success or exact error category/message boundary;
- requested atom and bond selections;
- every atom and bond `_CIPCode`;
- every atom and bond `_CIPNeighborOrder`;
- molecule `_CIPComputed` presence/value/computed status;
- `_CIPRank` before and after;
- bond stereo and stereo-atom state before and after;
- preservation of unselected labels;
- recursion-budget and repeated-call behavior.

The maintained evidence tiers are:

1. focused source-derived and upstream-regression fixtures;
2. the repository `smiles_small` corpus;
3. the complete maintained 5,000-row corpus;
4. a reproducible ChEMBL 37 phase where the measured cost is reasonable.

Fixtures belong under `testdata/stereo/fixtures/` and expected oracle output
under `testdata/stereo/expected/rdkit/<profile>/`. Reference generation belongs
under `tools/testdata/rdkit/` and remains reachable through
`tools/testdata/rdkit/generate_all.py`. Large-scale execution belongs in
`dev/tools/chembl_parity/`, not an ad hoc temporary script.

## Final Closure Audit

| Audit boundary | Closure result |
|---|---|
| Source functions and markers | Every reachable production function in the pinned `CIPLabeler` dispatcher has an inline source anchor and completed two-axis marker. Old `Rule5` remains unreachable because pinned `all_rules` installs `Rule5New`. |
| Exact-width state | Selection indexes and recursion counters use source-width `u32`; `Atom::NOATOM` is `4294967295`; PairList pairing retains the pinned 64-bit storage behavior. Boundary tests cover conversions, masks, sentinels, and decrement wrap behavior. |
| Sorter ownership | Each sequence rule receives the sorter prefix present when `Rules::add()` installed it. Rule4b and Rule5New replacement sorters retain the corresponding prefix. Focused regressions and the complete ChEMBL phase cover the recursive long tail. |
| Single modern implementation | `chemistry::ciplabeler` is the sole modern engine. `with_cip_labels_with_options`, its trailing-underscore form, Python methods, and fingerprint consumers delegate to it. Search found no second public or fingerprint-local modern algorithm. |
| Legacy reachability | `chemistry::stereo` retains source-required legacy rank/code behavior, but the registered modern operation and its public facades do not call it. Modern assignment preserves preexisting `_CIPRank` except where pinned computed-property clearing applies. |
| Atropisomer carriers | The CIP constructor consumes the shared source-backed carrier-order helper; duplicate CIP-local carrier extraction was removed. |
| Observable state | `_CIPCode`, `_CIPNeighborOrder`, `_CIPComputed`, `_CIPRank`, computed membership, chiral tags, bond stereo, and stereo atoms are compared. Clone and binary roundtrip preserve all supported state and computed classification. |
| Operation lifecycle | `CipStatePolicy` is generated with every operation spec. Only `with_cip_labels_with_options` may declare `Assign`; every topology operation is exercised by the generated lifecycle matrix as preserve or clear-computed. In-place failures keep structurally complete storage and source-order partial state. |
| Computed-property semantics | Molecule, atom, and bond computed membership is explicit. `clearComputedProps()` clears registered members, not hard-coded names, so deliberately non-computed properties with the same spelling are preserved like RDKit. The strict operation validator snapshots both value and membership and proves that computed members disappear while same-named non-computed members survive; the lifecycle and MMFF integration regressions exercise both sides. |
| Public Rust and Python APIs | Rust exposes `CipLabelOptions`, typed `CipDescriptor`, atom/bond queries, value-style assignment, and trailing-underscore in-place assignment. Python exposes full and selected assignment plus atom/bond descriptor queries through the same Rust operation. Generated stubs, tests, docs, and examples cover the surface. |
| Errors and unsupported boundaries | Recursion, index, configuration, node-cap, and malformed-state failures are structured with pinned messages where source-observable. COSMolKit has no equivalent process control handler, so requested cancellation returns an explicit unsupported error before mutation. Non-tetrahedral enum values outside pinned `findConfigs` are not claimed as assignment support. |
| Focused and maintained corpus gates | The focused oracle covers all ten emitted descriptor spellings, selection truth-value dispatch, repeated calls, recursion, malformed state, and boundary cases. Complete `smiles_small` and 5,000-row gates are maintained by the expected-data toolchain and CI. |
| Complete ChEMBL 37 phase | 2,854,362 accepted records were compared in each of full, selected-atom, selected-bond, and empty-selection branches: 11,417,448 exact complete-state comparisons, zero mismatch, 128/128 shards complete. Fifteen rows were rejected by both parsers and 43,442 accepted rows exceeded the configured 80-atom phase boundary. |

The public modern CIPLabeler capability is complete for the pinned
`findConfigs` dispatcher boundary. Adjacent full stereochemistry orchestration,
3D perception, source-equivalent process cancellation, and non-tetrahedral
configuration construction remain separate capabilities; their absence does
not represent a mismatch inside the supported modern assignment surface.
