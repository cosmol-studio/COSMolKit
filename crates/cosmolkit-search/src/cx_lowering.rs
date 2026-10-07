//! QueryGraph-native lowering for representation-independent CX records.

use cosmolkit_cx::{
    CxAtomConstraint, CxCoordinateBondKind, CxCountConstraint, CxDataSGroup,
    CxDoubleBondStereoKind, CxEnhancedStereo, CxLinkNode, CxPolymerSGroup, CxRecord,
    CxSGroupHierarchy, CxStereoGroupKind, CxVariableAttachment, CxWedgeBond, CxWedgeDirection,
    ParsedCxExtensions,
};
use cosmolkit_model::{
    AtomId, AtomQueryPredicate, BondDirection, BondId, BondOrder, BondStereo, Conformer3D,
    QueryAtom, QueryAtomIdentity, QueryGraph, QueryGraphError, QueryNode, SGroupData, StereoGroup,
    StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, query_substance_groups,
    replace_query_stereo_groups, replace_query_substance_groups,
};
use cosmolkit_types::{ChiralTag, Hybridization};

const QUERY_SCAN_MAGIC_VALUE: u32 = 0xDEAD_BEEF;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CxQueryLoweringError {
    #[error("CX unsigned property read failed: {0}")]
    Numeric(#[from] cosmolkit_core::PropertyUIntReadError),
    #[error("CX bond property operation failed: {0}")]
    BondProperty(#[from] cosmolkit_model::BondValueError),
    #[error("CX atom property operation failed: {0}")]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error("CX property string conversion failed: {0}")]
    Property(#[from] cosmolkit_core::PropertyStringError),
    #[error("CX molecule property assignment failed: {0}")]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
    #[error("CX atom index {index} is outside the query graph")]
    AtomIndex { index: usize },
    #[error("CX bond index {index} is outside the query graph")]
    BondIndex { index: usize },
    #[error("CX coordinate count {actual} does not match query atom count {expected}")]
    CoordinateCount { actual: usize, expected: usize },
    #[error("CX coordinate bond atom {atom} is not an endpoint of bond {bond}")]
    BondAtomMismatch { atom: usize, bond: usize },
    #[error("CX wedge atom {atom} is not an endpoint of bond {bond}")]
    WedgeAtomMismatch { atom: usize, bond: usize },
    #[error("CX record has no QueryGraph representation: {record}")]
    UnsupportedRecord { record: &'static str },
    #[error("query graph is invalid after CX lowering: {0}")]
    InvalidGraph(String),
}

/// Resolve a validated CX grammar index to its final QueryGraph bond row.
pub(crate) fn query_bond_row_from_source_index(
    graph: &QueryGraph,
    index: usize,
) -> Result<usize, CxQueryLoweringError> {
    // BEGIN COMPLETE PINNED SF184
    // RDKit✔️✔️: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
    // RDKit✔️✔️:   for (auto bnd : mol.bonds()) {
    // RDKit✔️✔️:     unsigned int smilesIdx;
    // RDKit✔️✔️:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
    // RDKit✔️✔️:         smilesIdx == idx) {
    // RDKit✔️✔️:       return bnd;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // END COMPLETE PINNED SF184
    // Complete reached RDProps.h
    // RDKit✔️✔️: bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit✔️✔️:     return d_props.getValIfPresent(key, res);
    // RDKit✔️✔️:   }
    // Complete reached Dict.h
    // RDKit✔️✔️: bool getValIfPresent(const std::string_view what, T &res) const {
    // RDKit✔️✔️:     for (const auto &data : _data) {
    // RDKit✔️✔️:       if (data.key == what) {
    // RDKit✔️✔️:         res = from_rdvalue<T>(data.val);
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // Complete reached RDValue-taggedunion.h
    // RDKit✔️✔️: inline unsigned int rdvalue_cast<unsigned int>(RDValue_cast_t v) {
    // RDKit✔️✔️:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit✔️✔️:     return v.value.u;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (rdvalue_is<int>(v)) {
    // RDKit✔️✔️:     return boost::numeric_cast<unsigned int>(v.value.i);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   throw std::bad_any_cast();
    // RDKit✔️✔️: }
    // Complete reached RDValue.h
    // RDKit✔️✔️: typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    // RDKit✔️✔️:     RDValue_cast_t arg) {
    // RDKit✔️✔️:   T res;
    // RDKit✔️✔️:   if (arg.getTag() == RDTypeTag::StringTag) {
    // RDKit✔️✔️:     Utils::LocaleSwitcher ls;
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       res = rdvalue_cast<T>(arg);
    // RDKit✔️✔️:     } catch (const std::bad_any_cast &exc) {
    // RDKit✔️✔️:       try {
    // RDKit✔️✔️: 	std::string val = rdvalue_cast<std::string>(arg);
    // RDKit✔️✔️: 	// trim only the right characters, this mimics how SD values
    // RDKit✔️✔️: 	//  work on read, they will be trimmed by the MolFile parser
    // RDKit✔️✔️: 	boost::trim_right(val);
    // RDKit✔️✔️:         res = boost::lexical_cast<T>(val);
    // RDKit✔️✔️:       } catch (...) {
    // RDKit✔️✔️:         throw exc;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = rdvalue_cast<T>(arg);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complete reached ROMol.h
    // RDKit✔️✔️: CXXBondIterator<const MolGraph, Bond *const> bonds() const {
    // RDKit✔️✔️:     return {&d_graph};
    // RDKit✔️✔️:   }
    // Complete reached ordered CXXBondIterator
    // RDKit✔️✔️: struct CXXBondIterator {
    // RDKit✔️✔️:   Graph *graph;
    // RDKit✔️✔️:   Iterator vstart, vend;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   struct CXXBondIter {
    // RDKit✔️✔️:     using iterator_category = std::forward_iterator_tag;
    // RDKit✔️✔️:     using difference_type = std::ptrdiff_t;
    // RDKit✔️✔️:     using value_type = Edge;
    // RDKit✔️✔️:     using pointer = Edge *;
    // RDKit✔️✔️:     using reference = Edge &;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     Graph *graph;
    // RDKit✔️✔️:     Iterator pos;
    // RDKit✔️✔️:     Bond *current;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     CXXBondIter(Graph *graph, Iterator pos)
    // RDKit✔️✔️:         : graph(graph), pos(pos), current(nullptr) {}
    // RDKit✔️✔️:
    // RDKit✔️✔️:     reference operator*() {
    // RDKit✔️✔️:       current = (*graph)[*pos];
    // RDKit✔️✔️:       return current;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     CXXBondIter &operator++() {
    // RDKit✔️✔️:       ++pos;
    // RDKit✔️✔️:       return *this;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bool operator==(const CXXBondIter &it) const { return pos == it.pos; }
    // RDKit✔️✔️:     bool operator!=(const CXXBondIter &it) const { return pos != it.pos; }
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CXXBondIterator(Graph *graph) : graph(graph) {
    // RDKit✔️✔️:     auto vs = boost::edges(*graph);
    // RDKit✔️✔️:     vstart = vs.first;
    // RDKit✔️✔️:     vend = vs.second;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   CXXBondIterator(Graph *graph, Iterator start, Iterator end)
    // RDKit✔️✔️:       : graph(graph), vstart(start), vend(end) {};
    // RDKit✔️✔️:   CXXBondIter begin() { return {graph, vstart}; }
    // RDKit✔️✔️:   CXXBondIter end() { return {graph, vend}; }
    // RDKit✔️✔️: }
    // Behavior: absent metadata continues, wrong present types/conversions
    // propagate the canonical structured failure, and first matching row stops
    // the scan before later properties are touched. Null maps to BondIndex at
    // this canonical required-row boundary; it is never replaced by a row ID.
    // Complexity: one ordered O(E) scan with O(1) scratch and no graph clone,
    // secondary index or repeated read. Numeric conversion has its sole CORE
    // owner; comparison is width-independent and cannot hide a conversion error.
    for (row, bond) in graph.bonds().iter().enumerate() {
        let Some(value) = bond.bond().prop("_cxsmilesBondIdx") else {
            continue;
        };
        let source_index = cosmolkit_core::property_value_to_uint(value)?;
        if u128::from(source_index) == index as u128 {
            return Ok(row);
        }
    }
    Err(CxQueryLoweringError::BondIndex { index })
}

pub(crate) struct CxStereoGroupTracker {
    hashes: Vec<u32>,
    first_group_index: usize,
}

impl CxStereoGroupTracker {
    pub(crate) fn new(graph: &QueryGraph) -> Self {
        Self {
            hashes: Vec::new(),
            first_group_index: graph.stereo_groups().len(),
        }
    }
}

pub(crate) fn merge_cx_enhanced_stereo(
    graph: &mut QueryGraph,
    tracker: &mut CxStereoGroupTracker,
    stereo: &CxEnhancedStereo,
) -> Result<(), QueryGraphError> {
    // RDKit source (verbatim; see CXSmilesOps.cpp::VALID_ATIDX and
    // parse_enhanced_stereo):
    /*
    #define VALID_ATIDX(_atidx_) \
      ((_atidx_) >= startAtomIdx && (_atidx_) < startAtomIdx + mol.getNumAtoms())
    if (VALID_ATIDX(aidx)) {
      Atom *atom = mol.getAtomWithIdx(aidx - startAtomIdx);
      if (!atom) {
        BOOST_LOG(rdWarningLog)
            << "Atom " << aidx << " not found!" << std::endl;
        return false;
      }
      atoms.push_back(atom);
    }
    */
    // RDKit✔️✔️: invalid source indexes are skipped; valid indexes append in order.
    // RDKit source (verbatim; see CXSmilesOps.cpp::parse_enhanced_stereo):
    /*
    if (!atoms.empty()) {
      const auto group_hash =
          10 * group_id + static_cast<unsigned int>(group_type);
      std::vector<unsigned int> sgTracker;
      mol.getPropIfPresent(cxsgTracker, sgTracker);
      std::vector<StereoGroup> mol_stereo_groups(mol.getStereoGroups());
      TEST_ASSERT(mol_stereo_groups.size() == sgTracker.size());

      auto iter = std::find(sgTracker.begin(), sgTracker.end(), group_hash);
      if (iter != sgTracker.end()) {
        auto index = iter - sgTracker.begin();
        auto gAtoms = mol_stereo_groups[index].getAtoms();
        gAtoms.insert(gAtoms.end(), atoms.begin(), atoms.end());
        mol_stereo_groups[index] =
            StereoGroup(mol_stereo_groups[index].getGroupType(),
                        std::move(gAtoms), std::move(bonds), group_id);
      } else {
        // not seen this before, create a new stereogroup
        mol_stereo_groups.emplace_back(group_type, std::move(atoms),
                                       std::move(bonds), group_id);
        sgTracker.push_back(group_hash);
        mol.setProp(cxsgTracker, sgTracker);
      }

      mol.setStereoGroups(std::move(mol_stereo_groups));
    }
    */
    // RDKit✔️✔️: the tracker stores first-seen group hashes in parallel order;
    // RDKit✔️✔️: repeated groups append atoms and reconstruct with the incoming ID.
    let (kind, kind_code) = match stereo.kind {
        CxStereoGroupKind::Absolute => (StereoGroupKind::Absolute, 0_u32),
        CxStereoGroupKind::Or => (StereoGroupKind::Or, 1_u32),
        CxStereoGroupKind::And => (StereoGroupKind::And, 2_u32),
    };
    let atoms = stereo
        .atoms
        .iter()
        .filter(|&&index| index < graph.num_atoms())
        .map(|&index| AtomId::new(index))
        .collect::<Vec<_>>();
    if atoms.is_empty() {
        return Ok(());
    }

    let group_hash = stereo.group_id.wrapping_mul(10).wrapping_add(kind_code);
    let mut groups = graph.stereo_groups().to_vec();
    let matched_index = tracker.hashes.iter().position(|hash| *hash == group_hash);
    debug_assert_eq!(
        groups.len(),
        tracker.first_group_index + tracker.hashes.len()
    );
    if let Some(position) = matched_index {
        let group_index = tracker.first_group_index + position;
        let previous_kind = groups[group_index].kind();
        let mut merged_atoms = groups[group_index].atoms().to_vec();
        merged_atoms.extend(atoms);
        groups[group_index] =
            StereoGroup::new(previous_kind, merged_atoms, Vec::new()).with_id(stereo.group_id);
    } else {
        groups.push(StereoGroup::new(kind, atoms, Vec::new()).with_id(stereo.group_id));
    }

    // Complexity review: one tracker scan and one group-vector clone per
    // nonempty record match the source's linear hash lookup and group copy.
    replace_query_stereo_groups(graph, groups)?;
    if matched_index.is_none() {
        tracker.hashes.push(group_hash);
    }
    Ok(())
}

pub(crate) fn cx_link_node_outer_atoms(
    graph: &QueryGraph,
    node: &CxLinkNode,
) -> Result<Option<[usize; 2]>, CxQueryLoweringError> {
    // BEGIN COMPLETE PINNED SF193
    // RDKit✔️❌: bool parse_linknodes(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️❌:                      unsigned int startAtomIdx) {
    // RDKit✔️❌:   // these look like: |LN:1:1.3.2.6,4:1.4.3.6|
    // RDKit✔️❌:   // that's two records:
    // RDKit✔️❌:   //   1:1.3.2.6: 1-3 repeats, atom 1-2, 1-6
    // RDKit✔️❌:   //   4:1.4.3.6: 1-4 repeats, atom 4-3, 4-6
    // RDKit✔️❌:   // which maps to the property value "1 3 2 2 3 2 7|1 4 2 5 4 5 7"
    // RDKit✔️❌:   // If the linking atom only has two neighbors then the outer atom
    // RDKit✔️❌:   // specification (the last two digits) can be left out. So for a molecule
    // RDKit✔️❌:   // where atom 1 has bonds only to atoms 2 and 6 we could have
    // RDKit✔️❌:   // |LN:1:1.3|
    // RDKit✔️❌:   // instead of
    // RDKit✔️❌:   // |LN:1:1.3.2.6|
    // RDKit✔️❌:   if (first >= last || *first != 'L' || first + 1 >= last ||
    // RDKit✔️❌:       *(first + 1) != 'N' || first + 2 >= last || *(first + 2) != ':') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   first += 3;
    // RDKit✔️❌:   std::string accum = "";
    // RDKit✔️❌:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️❌:     unsigned int atidx;
    // RDKit✔️❌:     if (!read_int(first, last, atidx)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // check that we can read at least two more characters:
    // RDKit✔️❌:     if (first + 1 >= last || *first != ':') {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++first;
    // RDKit✔️❌:     unsigned int startReps;
    // RDKit✔️❌:     if (!read_int(first, last, startReps)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first + 1 >= last || *first != '.') {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++first;
    // RDKit✔️❌:     unsigned int endReps;
    // RDKit✔️❌:     if (!read_int(first, last, endReps)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     unsigned int idx1;
    // RDKit✔️❌:     unsigned int idx2;
    // RDKit✔️❌:     if (first < last && *first == '.') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:       if (!read_int(first, last, idx1)) {
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:       if (!read_int(first, last, idx2)) {
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else if (VALID_ATIDX(atidx) &&
    // RDKit✔️❌:                mol.getAtomWithIdx(atidx - startAtomIdx)->getDegree() == 2) {
    // RDKit✔️❌:       auto nbrs =
    // RDKit✔️❌:           mol.getAtomNeighbors(mol.getAtomWithIdx(atidx - startAtomIdx));
    // RDKit✔️❌:       idx1 = *nbrs.first;
    // RDKit✔️❌:       nbrs.first++;
    // RDKit✔️❌:       idx2 = *nbrs.first;
    // RDKit✔️❌:     } else if (VALID_ATIDX(atidx)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first < last && *first == ',') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (VALID_ATIDX(atidx)) {
    // RDKit✔️❌:       if (!accum.empty()) {
    // RDKit✔️❌:         accum += "|";
    // RDKit✔️❌:       }
    // RDKit✔️❌:       accum += (boost::format("%d %d 2 %d %d %d %d") % startReps % endReps %
    // RDKit✔️❌:                 (atidx - startAtomIdx + 1) % (idx1 - startAtomIdx + 1) %
    // RDKit✔️❌:                 (atidx - startAtomIdx + 1) % (idx2 - startAtomIdx + 1))
    // RDKit✔️❌:                    .str();
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!accum.empty()) {
    // RDKit✔️❌:     mol.setProp(common_properties::molFileLinkNodes, accum);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF193
    // BEGIN COMPLETE Atom::getDegree
    // RDKit✔️✔️: unsigned int Atom::getDegree() const {
    // RDKit✔️✔️:   return dp_mol ? getOwningMol().getAtomDegree(this) : 0;
    // RDKit✔️✔️: }
    // END COMPLETE Atom::getDegree
    // BEGIN COMPLETE ROMol::getAtomDegree
    // RDKit✔️✔️: unsigned int ROMol::getAtomDegree(const Atom *at) const {
    // RDKit✔️✔️:   PRECONDITION(at, "no atom");
    // RDKit✔️✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit✔️✔️:                "atom not associated with this molecule");
    // RDKit✔️✔️:   return rdcast<unsigned int>(boost::out_degree(at->getIdx(), d_graph));
    // RDKit✔️✔️: }
    // END COMPLETE ROMol::getAtomDegree
    // BEGIN COMPLETE ROMol::getAtomNeighbors
    // RDKit✔️✔️: ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
    // RDKit✔️✔️:   PRECONDITION(at, "no atom");
    // RDKit✔️✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit✔️✔️:                "atom not associated with this molecule");
    // RDKit✔️✔️:   return boost::adjacent_vertices(at->getIdx(), d_graph);
    // RDKit✔️✔️: }
    // END COMPLETE ROMol::getAtomNeighbors
    // BEGIN COMPLETE Boost1.81::adjacent_vertices
    // Boost✔️✔️: adjacent_vertices(typename Config::vertex_descriptor u,
    // Boost✔️✔️:     const adj_list_helper< Config, Base >& g_)
    // Boost✔️✔️: {
    // Boost✔️✔️:     typedef typename Config::graph_type AdjList;
    // Boost✔️✔️:     const AdjList& cg = static_cast< const AdjList& >(g_);
    // Boost✔️✔️:     AdjList& g = const_cast< AdjList& >(cg);
    // Boost✔️✔️:     typedef typename Config::adjacency_iterator adjacency_iterator;
    // Boost✔️✔️:     typename Config::out_edge_iterator first, last;
    // Boost✔️✔️:     boost::tie(first, last) = out_edges(u, g);
    // Boost✔️✔️:     return std::make_pair(
    // Boost✔️✔️:         adjacency_iterator(first, &g), adjacency_iterator(last, &g));
    // Boost✔️✔️: }
    // END COMPLETE Boost1.81::adjacent_vertices
    // BEGIN COMPLETE Boost1.81::out_edges
    // Boost✔️✔️: out_edges(typename Config::vertex_descriptor u,
    // Boost✔️✔️:     bidirectional_graph_helper_with_property< Config >& g_)
    // Boost✔️✔️: {
    // Boost✔️✔️:     typedef typename Config::global_edgelist_selector EdgeListS;
    // Boost✔️✔️:     BOOST_STATIC_ASSERT((!is_same< EdgeListS, vecS >::value));
    // Boost✔️✔️:
    // Boost✔️✔️:     typedef typename Config::graph_type graph_type;
    // Boost✔️✔️:     typedef typename Config::edge_parallel_category Cat;
    // Boost✔️✔️:     graph_type& g = static_cast< graph_type& >(g_);
    // Boost✔️✔️:     typename Config::OutEdgeList& el = g.out_edge_list(u);
    // Boost✔️✔️:     typename Config::OutEdgeList::iterator ei = el.begin(), ei_end = el.end();
    // Boost✔️✔️:     for (; ei != ei_end; ++ei)
    // Boost✔️✔️:     {
    // Boost✔️✔️:         detail::erase_from_incidence_list(
    // Boost✔️✔️:             in_edge_list(g, (*ei).get_target()), u, Cat());
    // Boost✔️✔️:         g.m_edges.erase((*ei).get_iter());
    // Boost✔️✔️:     }
    // Boost✔️✔️:     g.out_edge_list(u).clear();
    // Boost✔️✔️: }
    // END COMPLETE Boost1.81::out_edges
    // BEGIN COMPLETE Boost1.81::out_degree
    // Boost✔️✔️: inline typename Config::degree_size_type out_degree(
    // Boost✔️✔️:     typename Config::vertex_descriptor u,
    // Boost✔️✔️:     const adj_list_helper< Config, Base >& g_)
    // Boost✔️✔️: {
    // Boost✔️✔️:     typedef typename Config::graph_type AdjList;
    // Boost✔️✔️:     const AdjList& g = static_cast< const AdjList& >(g_);
    // Boost✔️✔️:     return g.out_edge_list(u).size();
    // Boost✔️✔️: }
    // END COMPLETE Boost1.81::out_degree
    // Source global window is zero-based in this canonical local consumer.
    // Skip before neighbor access; explicit indices are not range-validated.
    // O(1) adjacency lookup/degree and the first two insertion-order neighbors.
    if node.atom >= graph.num_atoms() {
        return Ok(None);
    }
    if let Some(outer) = node.outer_atoms {
        return Ok(Some(outer));
    }
    let neighbors = graph.adjacency().get(node.atom).ok_or_else(|| {
        CxQueryLoweringError::InvalidGraph("CX link-node center has no adjacency row".to_owned())
    })?;
    if neighbors.len() != 2 {
        return Err(CxQueryLoweringError::InvalidGraph(format!(
            "CX link-node atom {} has degree {}, expected two when outer atoms are omitted",
            node.atom,
            neighbors.len()
        )));
    }
    Ok(Some([neighbors[0].0, neighbors[1].0]))
}

pub(crate) fn apply_cx_link_nodes_to_query(
    graph: &mut QueryGraph,
    nodes: &[CxLinkNode],
) -> Result<(), CxQueryLoweringError> {
    // BEGIN COMPLETE PINNED SF193
    // RDKit✔️❌: bool parse_linknodes(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️❌:                      unsigned int startAtomIdx) {
    // RDKit✔️❌:   // these look like: |LN:1:1.3.2.6,4:1.4.3.6|
    // RDKit✔️❌:   // that's two records:
    // RDKit✔️❌:   //   1:1.3.2.6: 1-3 repeats, atom 1-2, 1-6
    // RDKit✔️❌:   //   4:1.4.3.6: 1-4 repeats, atom 4-3, 4-6
    // RDKit✔️❌:   // which maps to the property value "1 3 2 2 3 2 7|1 4 2 5 4 5 7"
    // RDKit✔️❌:   // If the linking atom only has two neighbors then the outer atom
    // RDKit✔️❌:   // specification (the last two digits) can be left out. So for a molecule
    // RDKit✔️❌:   // where atom 1 has bonds only to atoms 2 and 6 we could have
    // RDKit✔️❌:   // |LN:1:1.3|
    // RDKit✔️❌:   // instead of
    // RDKit✔️❌:   // |LN:1:1.3.2.6|
    // RDKit✔️❌:   if (first >= last || *first != 'L' || first + 1 >= last ||
    // RDKit✔️❌:       *(first + 1) != 'N' || first + 2 >= last || *(first + 2) != ':') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   first += 3;
    // RDKit✔️❌:   std::string accum = "";
    // RDKit✔️❌:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️❌:     unsigned int atidx;
    // RDKit✔️❌:     if (!read_int(first, last, atidx)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // check that we can read at least two more characters:
    // RDKit✔️❌:     if (first + 1 >= last || *first != ':') {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++first;
    // RDKit✔️❌:     unsigned int startReps;
    // RDKit✔️❌:     if (!read_int(first, last, startReps)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first + 1 >= last || *first != '.') {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++first;
    // RDKit✔️❌:     unsigned int endReps;
    // RDKit✔️❌:     if (!read_int(first, last, endReps)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     unsigned int idx1;
    // RDKit✔️❌:     unsigned int idx2;
    // RDKit✔️❌:     if (first < last && *first == '.') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:       if (!read_int(first, last, idx1)) {
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:       if (!read_int(first, last, idx2)) {
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else if (VALID_ATIDX(atidx) &&
    // RDKit✔️❌:                mol.getAtomWithIdx(atidx - startAtomIdx)->getDegree() == 2) {
    // RDKit✔️❌:       auto nbrs =
    // RDKit✔️❌:           mol.getAtomNeighbors(mol.getAtomWithIdx(atidx - startAtomIdx));
    // RDKit✔️❌:       idx1 = *nbrs.first;
    // RDKit✔️❌:       nbrs.first++;
    // RDKit✔️❌:       idx2 = *nbrs.first;
    // RDKit✔️❌:     } else if (VALID_ATIDX(atidx)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first < last && *first == ',') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (VALID_ATIDX(atidx)) {
    // RDKit✔️❌:       if (!accum.empty()) {
    // RDKit✔️❌:         accum += "|";
    // RDKit✔️❌:       }
    // RDKit✔️❌:       accum += (boost::format("%d %d 2 %d %d %d %d") % startReps % endReps %
    // RDKit✔️❌:                 (atidx - startAtomIdx + 1) % (idx1 - startAtomIdx + 1) %
    // RDKit✔️❌:                 (atidx - startAtomIdx + 1) % (idx2 - startAtomIdx + 1))
    // RDKit✔️❌:                    .str();
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!accum.empty()) {
    // RDKit✔️❌:     mol.setProp(common_properties::molFileLinkNodes, accum);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF193
    // BEGIN COMPLETE source common_properties::molFileLinkNodes
    // RDKit✔️✔️: inline constexpr std::string_view molFileLinkNodes = "_molLinkNodes";
    // END COMPLETE source common_properties::molFileLinkNodes
    // BEGIN COMPLETE Boost1.81::parse_printf_directive
    // Boost✔️✔️:     bool parse_printf_directive(Iter & start, const Iter& last,
    // Boost✔️✔️:                                 detail::format_item<Ch, Tr, Alloc> * fpar,
    // Boost✔️✔️:                                 const Facet& fac,
    // Boost✔️✔️:                                 std::size_t offset, unsigned char exceptions)
    // Boost✔️✔️:     {
    // Boost✔️✔️:         typedef typename basic_format<Ch, Tr, Alloc>::format_item_t format_item_t;
    // Boost✔️✔️:
    // Boost✔️✔️:         fpar->argN_ = format_item_t::argN_no_posit;  // if no positional-directive
    // Boost✔️✔️:         bool precision_set = false;
    // Boost✔️✔️:         bool in_brackets=false;
    // Boost✔️✔️:         Iter start0 = start;
    // Boost✔️✔️:         std::size_t fstring_size = last-start0+offset;
    // Boost✔️✔️:         char mssiz = 0;
    // Boost✔️✔️:
    // Boost✔️✔️:         if(start>= last) { // empty directive : this is a trailing %
    // Boost✔️✔️:                 maybe_throw_exception(exceptions, start-start0 + offset, fstring_size);
    // Boost✔️✔️:                 return false;
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:         if(*start== const_or_not(fac).widen( '|')) {
    // Boost✔️✔️:             in_brackets=true;
    // Boost✔️✔️:             if( ++start >= last ) {
    // Boost✔️✔️:                 maybe_throw_exception(exceptions, start-start0 + offset, fstring_size);
    // Boost✔️✔️:                 return false;
    // Boost✔️✔️:             }
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:         // the flag '0' would be picked as a digit for argument order, but here it's a flag :
    // Boost✔️✔️:         if(*start== const_or_not(fac).widen( '0'))
    // Boost✔️✔️:             goto parse_flags;
    // Boost✔️✔️:
    // Boost✔️✔️:         // handle argument order (%2$d)  or possibly width specification: %2d
    // Boost✔️✔️:         if(wrap_isdigit(fac, *start)) {
    // Boost✔️✔️:             int n;
    // Boost✔️✔️:             start = str2int(start, last, n, fac);
    // Boost✔️✔️:             if( start >= last ) {
    // Boost✔️✔️:                 maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
    // Boost✔️✔️:                 return false;
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:             // %N% case : this is already the end of the directive
    // Boost✔️✔️:             if( *start ==  const_or_not(fac).widen( '%') ) {
    // Boost✔️✔️:                 fpar->argN_ = n-1;
    // Boost✔️✔️:                 ++start;
    // Boost✔️✔️:                 if( in_brackets)
    // Boost✔️✔️:                     maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
    // Boost✔️✔️:                 return true;
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:             if ( *start== const_or_not(fac).widen( '$') ) {
    // Boost✔️✔️:                 fpar->argN_ = n-1;
    // Boost✔️✔️:                 ++start;
    // Boost✔️✔️:             }
    // Boost✔️✔️:             else {
    // Boost✔️✔️:                 // non-positional directive
    // Boost✔️✔️:                 fpar->fmtstate_.width_ = n;
    // Boost✔️✔️:                 fpar->argN_  = format_item_t::argN_no_posit;
    // Boost✔️✔️:                 goto parse_precision;
    // Boost✔️✔️:             }
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:       parse_flags:
    // Boost✔️✔️:         // handle flags
    // Boost✔️✔️:         while (start != last) { // as long as char is one of + - = _ # 0 or ' '
    // Boost✔️✔️:             switch ( wrap_narrow(fac, *start, 0)) {
    // Boost✔️✔️:                 case '\'':
    // Boost✔️✔️:                     break; // no effect yet. (painful to implement)
    // Boost✔️✔️:                 case '-':
    // Boost✔️✔️:                     fpar->fmtstate_.flags_ |= std::ios_base::left;
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case '=':
    // Boost✔️✔️:                     fpar->pad_scheme_ |= format_item_t::centered;
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case '_':
    // Boost✔️✔️:                     fpar->fmtstate_.flags_ |= std::ios_base::internal;
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case ' ':
    // Boost✔️✔️:                     fpar->pad_scheme_ |= format_item_t::spacepad;
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case '+':
    // Boost✔️✔️:                     fpar->fmtstate_.flags_ |= std::ios_base::showpos;
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case '0':
    // Boost✔️✔️:                     fpar->pad_scheme_ |= format_item_t::zeropad;
    // Boost✔️✔️:                     // need to know alignment before really setting flags,
    // Boost✔️✔️:                     // so just add 'zeropad' flag for now, it will be processed later.
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case '#':
    // Boost✔️✔️:                     fpar->fmtstate_.flags_ |= std::ios_base::showpoint | std::ios_base::showbase;
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 default:
    // Boost✔️✔️:                     goto parse_width;
    // Boost✔️✔️:             }
    // Boost✔️✔️:             ++start;
    // Boost✔️✔️:         } // loop on flag.
    // Boost✔️✔️:
    // Boost✔️✔️:         if( start>=last) {
    // Boost✔️✔️:             maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
    // Boost✔️✔️:             return true;
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:       // first skip 'asterisk fields' : * or num (length)
    // Boost✔️✔️:       parse_width:
    // Boost✔️✔️:         if(*start == const_or_not(fac).widen( '*') )
    // Boost✔️✔️:             ++start;
    // Boost✔️✔️:         else if(start!=last && wrap_isdigit(fac, *start))
    // Boost✔️✔️:             start = str2int(start, last, fpar->fmtstate_.width_, fac);
    // Boost✔️✔️:
    // Boost✔️✔️:       parse_precision:
    // Boost✔️✔️:         if( start>= last) {
    // Boost✔️✔️:             maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
    // Boost✔️✔️:             return true;
    // Boost✔️✔️:         }
    // Boost✔️✔️:         // handle precision spec
    // Boost✔️✔️:         if (*start== const_or_not(fac).widen( '.')) {
    // Boost✔️✔️:             ++start;
    // Boost✔️✔️:             if(start != last && *start == const_or_not(fac).widen( '*') )
    // Boost✔️✔️:                 ++start;
    // Boost✔️✔️:             else if(start != last && wrap_isdigit(fac, *start)) {
    // Boost✔️✔️:                 start = str2int(start, last, fpar->fmtstate_.precision_, fac);
    // Boost✔️✔️:                 precision_set = true;
    // Boost✔️✔️:             }
    // Boost✔️✔️:             else
    // Boost✔️✔️:                 fpar->fmtstate_.precision_ =0;
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:       // argument type modifiers
    // Boost✔️✔️:         while (start != last) {
    // Boost✔️✔️:             switch (wrap_narrow(fac, *start, 0)) {
    // Boost✔️✔️:                 case 'h':
    // Boost✔️✔️:                 case 'l':
    // Boost✔️✔️:                 case 'j':
    // Boost✔️✔️:                 case 'z':
    // Boost✔️✔️:                 case 'L':
    // Boost✔️✔️:                     // boost::format ignores argument type modifiers as it relies on
    // Boost✔️✔️:                     // the type of the argument fed into it by operator %
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Note that the ptrdiff_t argument type 't' from C++11 is not honored
    // Boost✔️✔️:                 // because it was already in use as the tabulation specifier in boost::format
    // Boost✔️✔️:                 // case 't':
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Microsoft extensions:
    // Boost✔️✔️:                 // https://msdn.microsoft.com/en-us/library/tcxf1dw6.aspx
    // Boost✔️✔️:
    // Boost✔️✔️:                 case 'w':
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case 'I':
    // Boost✔️✔️:                     mssiz = 'I';
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case '3':
    // Boost✔️✔️:                     if (mssiz != 'I') {
    // Boost✔️✔️:                         maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
    // Boost✔️✔️:                         return true;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                     mssiz = '3';
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case '2':
    // Boost✔️✔️:                     if (mssiz != '3') {
    // Boost✔️✔️:                         maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
    // Boost✔️✔️:                         return true;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                     mssiz = 0x00;
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case '6':
    // Boost✔️✔️:                     if (mssiz != 'I') {
    // Boost✔️✔️:                         maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
    // Boost✔️✔️:                         return true;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                     mssiz = '6';
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 case '4':
    // Boost✔️✔️:                     if (mssiz != '6') {
    // Boost✔️✔️:                         maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
    // Boost✔️✔️:                         return true;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                     mssiz = 0x00;
    // Boost✔️✔️:                     break;
    // Boost✔️✔️:                 default:
    // Boost✔️✔️:                     if (mssiz && mssiz == 'I') {
    // Boost✔️✔️:                         mssiz = 0;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                     goto parse_conversion_specification;
    // Boost✔️✔️:             }
    // Boost✔️✔️:             ++start;
    // Boost✔️✔️:         } // loop on argument type modifiers to pick up 'hh', 'll', and the more complex microsoft ones
    // Boost✔️✔️:
    // Boost✔️✔️:       parse_conversion_specification:
    // Boost✔️✔️:         if (start >= last || mssiz) {
    // Boost✔️✔️:             maybe_throw_exception(exceptions, start - start0 + offset, fstring_size);
    // Boost✔️✔️:             return true;
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:         if( in_brackets && *start== const_or_not(fac).widen( '|') ) {
    // Boost✔️✔️:             ++start;
    // Boost✔️✔️:             return true;
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:         // The default flags are "dec" and "skipws"
    // Boost✔️✔️:         // so if changing the base, need to unset basefield first
    // Boost✔️✔️:
    // Boost✔️✔️:         switch (wrap_narrow(fac, *start, 0))
    // Boost✔️✔️:         {
    // Boost✔️✔️:             // Boolean
    // Boost✔️✔️:             case 'b':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::boolalpha;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:
    // Boost✔️✔️:             // Decimal
    // Boost✔️✔️:             case 'u':
    // Boost✔️✔️:             case 'd':
    // Boost✔️✔️:             case 'i':
    // Boost✔️✔️:                 // Defaults are sufficient
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:
    // Boost✔️✔️:             // Hex
    // Boost✔️✔️:             case 'X':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
    // Boost✔️✔️:                 BOOST_FALLTHROUGH;
    // Boost✔️✔️:             case 'x':
    // Boost✔️✔️:             case 'p': // pointer => set hex.
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::hex;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:
    // Boost✔️✔️:             // Octal
    // Boost✔️✔️:             case 'o':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::oct;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:
    // Boost✔️✔️:             // Floating
    // Boost✔️✔️:             case 'A':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
    // Boost✔️✔️:                 BOOST_FALLTHROUGH;
    // Boost✔️✔️:             case 'a':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ &= ~std::ios_base::basefield;
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::fixed;
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::scientific;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:             case 'E':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
    // Boost✔️✔️:                 BOOST_FALLTHROUGH;
    // Boost✔️✔️:             case 'e':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::scientific;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:             case 'F':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
    // Boost✔️✔️:                 BOOST_FALLTHROUGH;
    // Boost✔️✔️:             case 'f':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::fixed;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:             case 'G':
    // Boost✔️✔️:                 fpar->fmtstate_.flags_ |= std::ios_base::uppercase;
    // Boost✔️✔️:                 BOOST_FALLTHROUGH;
    // Boost✔️✔️:             case 'g':
    // Boost✔️✔️:                 // default flags are correct here
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:
    // Boost✔️✔️:             // Tabulation (a boost::format extension)
    // Boost✔️✔️:             case 'T':
    // Boost✔️✔️:                 ++start;
    // Boost✔️✔️:                 if( start >= last) {
    // Boost✔️✔️:                     maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
    // Boost✔️✔️:                     return false;
    // Boost✔️✔️:                 } else {
    // Boost✔️✔️:                     fpar->fmtstate_.fill_ = *start;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:                 fpar->pad_scheme_ |= format_item_t::tabulation;
    // Boost✔️✔️:                 fpar->argN_ = format_item_t::argN_tabulation;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:             case 't':
    // Boost✔️✔️:                 fpar->fmtstate_.fill_ = const_or_not(fac).widen( ' ');
    // Boost✔️✔️:                 fpar->pad_scheme_ |= format_item_t::tabulation;
    // Boost✔️✔️:                 fpar->argN_ = format_item_t::argN_tabulation;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:
    // Boost✔️✔️:             // Character
    // Boost✔️✔️:             case 'C':
    // Boost✔️✔️:             case 'c':
    // Boost✔️✔️:                 fpar->truncate_ = 1;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:
    // Boost✔️✔️:             // String
    // Boost✔️✔️:             case 'S':
    // Boost✔️✔️:             case 's':
    // Boost✔️✔️:                 if(precision_set) // handle truncation manually, with own parameter.
    // Boost✔️✔️:                     fpar->truncate_ = fpar->fmtstate_.precision_;
    // Boost✔️✔️:                 fpar->fmtstate_.precision_ = 6; // default stream precision.
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:
    // Boost✔️✔️:             // %n is insecure and ignored by boost::format
    // Boost✔️✔️:             case 'n' :
    // Boost✔️✔️:                 fpar->argN_ = format_item_t::argN_ignored;
    // Boost✔️✔️:                 break;
    // Boost✔️✔️:
    // Boost✔️✔️:             default:
    // Boost✔️✔️:                 maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
    // Boost✔️✔️:         }
    // Boost✔️✔️:         ++start;
    // Boost✔️✔️:
    // Boost✔️✔️:         if( in_brackets ) {
    // Boost✔️✔️:             if( start != last && *start== const_or_not(fac).widen( '|') ) {
    // Boost✔️✔️:                 ++start;
    // Boost✔️✔️:                 return true;
    // Boost✔️✔️:             }
    // Boost✔️✔️:             else  maybe_throw_exception(exceptions, start-start0+offset, fstring_size);
    // Boost✔️✔️:         }
    // Boost✔️✔️:         return true;
    // Boost✔️✔️:     }
    // END COMPLETE Boost1.81::parse_printf_directive
    // BEGIN COMPLETE Boost1.81::put
    // Boost✔️✔️:     void put( T x,
    // Boost✔️✔️:               const format_item<Ch, Tr, Alloc>& specs,
    // Boost✔️✔️:               typename basic_format<Ch, Tr, Alloc>::string_type& res,
    // Boost✔️✔️:               typename basic_format<Ch, Tr, Alloc>::internal_streambuf_t & buf,
    // Boost✔️✔️:               io::detail::locale_t *loc_p = NULL)
    // Boost✔️✔️:     {
    // Boost✔️✔️: #ifdef BOOST_MSVC
    // Boost✔️✔️:        // If std::min<unsigned> or std::max<unsigned> are already instantiated
    // Boost✔️✔️:        // at this point then we get a blizzard of warning messages when we call
    // Boost✔️✔️:        // those templates with std::size_t as arguments.  Weird and very annoyning...
    // Boost✔️✔️: #pragma warning(push)
    // Boost✔️✔️: #pragma warning(disable:4267)
    // Boost✔️✔️: #endif
    // Boost✔️✔️:         // does the actual conversion of x, with given params, into a string
    // Boost✔️✔️:         // using the supplied stringbuf.
    // Boost✔️✔️:
    // Boost✔️✔️:         typedef typename basic_format<Ch, Tr, Alloc>::string_type   string_type;
    // Boost✔️✔️:         typedef typename basic_format<Ch, Tr, Alloc>::format_item_t format_item_t;
    // Boost✔️✔️:         typedef typename string_type::size_type size_type;
    // Boost✔️✔️:
    // Boost✔️✔️:         basic_oaltstringstream<Ch, Tr, Alloc>  oss( &buf);
    // Boost✔️✔️:
    // Boost✔️✔️: #if !defined(BOOST_NO_STD_LOCALE)
    // Boost✔️✔️:         if(loc_p != NULL)
    // Boost✔️✔️:             oss.imbue(*loc_p);
    // Boost✔️✔️: #endif
    // Boost✔️✔️:
    // Boost✔️✔️:         specs.fmtstate_.apply_on(oss, loc_p);
    // Boost✔️✔️:
    // Boost✔️✔️:         // the stream format state can be modified by manipulators in the argument :
    // Boost✔️✔️:         put_head( oss, x );
    // Boost✔️✔️:         // in case x is a group, apply the manip part of it,
    // Boost✔️✔️:         // in order to find width
    // Boost✔️✔️:
    // Boost✔️✔️:         const std::ios_base::fmtflags fl=oss.flags();
    // Boost✔️✔️:         const bool internal = (fl & std::ios_base::internal) != 0;
    // Boost✔️✔️:         const std::streamsize w = oss.width();
    // Boost✔️✔️:         const bool two_stepped_padding= internal && (w!=0);
    // Boost✔️✔️:
    // Boost✔️✔️:         res.resize(0);
    // Boost✔️✔️:         if(! two_stepped_padding) {
    // Boost✔️✔️:             if(w>0) // handle padding via mk_str, not natively in stream
    // Boost✔️✔️:                 oss.width(0);
    // Boost✔️✔️:             put_last( oss, x);
    // Boost✔️✔️:             const Ch * res_beg = buf.pbase();
    // Boost✔️✔️:             Ch prefix_space = 0;
    // Boost✔️✔️:             if(specs.pad_scheme_ & format_item_t::spacepad)
    // Boost✔️✔️:                 if(buf.pcount()== 0 ||
    // Boost✔️✔️:                    (res_beg[0] !=oss.widen('+') && res_beg[0] !=oss.widen('-')  ))
    // Boost✔️✔️:                     prefix_space = oss.widen(' ');
    // Boost✔️✔️:             size_type res_size = (std::min)(
    // Boost✔️✔️:                 (static_cast<size_type>((specs.truncate_ & (std::numeric_limits<size_type>::max)())) - !!prefix_space),
    // Boost✔️✔️:                 buf.pcount() );
    // Boost✔️✔️:             mk_str(res, res_beg, res_size, w, oss.fill(), fl,
    // Boost✔️✔️:                    prefix_space, (specs.pad_scheme_ & format_item_t::centered) !=0 );
    // Boost✔️✔️:         }
    // Boost✔️✔️:         else  { // 2-stepped padding
    // Boost✔️✔️:             // internal can be implied by zeropad, or user-set.
    // Boost✔️✔️:             // left, right, and centered alignment overrule internal,
    // Boost✔️✔️:             // but spacepad or truncate might be mixed with internal (using manipulator)
    // Boost✔️✔️:             put_last( oss, x); // may pad
    // Boost✔️✔️:             const Ch * res_beg = buf.pbase();
    // Boost✔️✔️:             size_type res_size = buf.pcount();
    // Boost✔️✔️:             bool prefix_space=false;
    // Boost✔️✔️:             if(specs.pad_scheme_ & format_item_t::spacepad)
    // Boost✔️✔️:                 if(buf.pcount()== 0 ||
    // Boost✔️✔️:                    (res_beg[0] !=oss.widen('+') && res_beg[0] !=oss.widen('-')  ))
    // Boost✔️✔️:                     prefix_space = true;
    // Boost✔️✔️:             if(res_size == static_cast<size_type>(w) && w<=specs.truncate_ && !prefix_space) {
    // Boost✔️✔️:                 // okay, only one thing was printed and padded, so res is fine
    // Boost✔️✔️:                 res.assign(res_beg, res_size);
    // Boost✔️✔️:             }
    // Boost✔️✔️:             else { //   length w exceeded
    // Boost✔️✔️:                 // either it was multi-output with first output padding up all width..
    // Boost✔️✔️:                 // either it was one big arg and we are fine.
    // Boost✔️✔️:                 // Note that res_size<w is possible  (in case of bad user-defined formatting)
    // Boost✔️✔️:                 res.assign(res_beg, res_size);
    // Boost✔️✔️:                 res_beg=NULL;  // invalidate pointers.
    // Boost✔️✔️:
    // Boost✔️✔️:                 // make a new stream, to start re-formatting from scratch :
    // Boost✔️✔️:                 buf.clear_buffer();
    // Boost✔️✔️:                 basic_oaltstringstream<Ch, Tr, Alloc>  oss2( &buf);
    // Boost✔️✔️:                 specs.fmtstate_.apply_on(oss2, loc_p);
    // Boost✔️✔️:                 put_head( oss2, x );
    // Boost✔️✔️:
    // Boost✔️✔️:                 oss2.width(0);
    // Boost✔️✔️:                 if(prefix_space)
    // Boost✔️✔️:                     oss2 << ' ';
    // Boost✔️✔️:                 put_last(oss2, x );
    // Boost✔️✔️:                 if(buf.pcount()==0 && specs.pad_scheme_ & format_item_t::spacepad) {
    // Boost✔️✔️:                     prefix_space =true;
    // Boost✔️✔️:                     oss2 << ' ';
    // Boost✔️✔️:                 }
    // Boost✔️✔️:                 // we now have the minimal-length output
    // Boost✔️✔️:                 const Ch * tmp_beg = buf.pbase();
    // Boost✔️✔️:                 size_type tmp_size = (std::min)(
    // Boost✔️✔️:                     (static_cast<size_type>(specs.truncate_ & (std::numeric_limits<size_type>::max)())),
    // Boost✔️✔️:                     buf.pcount());
    // Boost✔️✔️:
    // Boost✔️✔️:                 if(static_cast<size_type>(w) <= tmp_size) {
    // Boost✔️✔️:                     // minimal length is already >= w, so no padding (cool!)
    // Boost✔️✔️:                         res.assign(tmp_beg, tmp_size);
    // Boost✔️✔️:                 }
    // Boost✔️✔️:                 else { // hum..  we need to pad (multi_output, or spacepad present)
    // Boost✔️✔️:                     //find where we should pad
    // Boost✔️✔️:                     size_type sz = (std::min)(res_size + (prefix_space ? 1 : 0), tmp_size);
    // Boost✔️✔️:                     size_type i = prefix_space;
    // Boost✔️✔️:                     for(; i<sz && tmp_beg[i] == res[i - (prefix_space ? 1 : 0)]; ++i) {}
    // Boost✔️✔️:                     if(i>=tmp_size) i=prefix_space;
    // Boost✔️✔️:                     res.assign(tmp_beg, i);
    // Boost✔️✔️:                                         std::streamsize d = w - static_cast<std::streamsize>(tmp_size);
    // Boost✔️✔️:                                         BOOST_ASSERT(d>0);
    // Boost✔️✔️:                     res.append(static_cast<size_type>( d ), oss2.fill());
    // Boost✔️✔️:                     res.append(tmp_beg+i, tmp_size-i);
    // Boost✔️✔️:                     BOOST_ASSERT(i+(tmp_size-i)+(std::max)(d,(std::streamsize)0)
    // Boost✔️✔️:                                  == static_cast<size_type>(w));
    // Boost✔️✔️:                     BOOST_ASSERT(res.size() == static_cast<size_type>(w));
    // Boost✔️✔️:                 }
    // Boost✔️✔️:             }
    // Boost✔️✔️:         }
    // Boost✔️✔️:         buf.clear_buffer();
    // Boost✔️✔️: #ifdef BOOST_MSVC
    // Boost✔️✔️: #pragma warning(pop)
    // Boost✔️✔️: #endif
    // Boost✔️✔️:     }
    // END COMPLETE Boost1.81::put
    // BEGIN COMPLETE Boost1.81::put_last<unsigned>
    // Boost✔️✔️:     void put_last( BOOST_IO_STD basic_ostream<Ch, Tr> & os, const T& x ) {
    // Boost✔️✔️:         os << x ;
    // Boost✔️✔️:     }
    // END COMPLETE Boost1.81::put_last<unsigned>
    // Accumulate locally; any later item failure keeps the prior property.
    // Native u32 +1 wraps before unsigned decimal formatting. Empty output
    // never clears prior state. Same linear accumulation; source property
    // tree storage and retained parse records add allocation overhead.
    let mut accum = String::new();
    for node in nodes {
        let Some([outer_one, outer_two]) = cx_link_node_outer_atoms(graph, node)? else {
            continue;
        };

        let source_uint = |value: usize| {
            u32::try_from(value).map_err(|_| {
                CxQueryLoweringError::InvalidGraph(
                    "CX link-node integer exceeds the source unsigned range".to_owned(),
                )
            })
        };
        let center_one = source_uint(node.atom)?.wrapping_add(1);
        let outer_one = source_uint(outer_one)?.wrapping_add(1);
        let outer_two = source_uint(outer_two)?.wrapping_add(1);
        let start_repetitions = source_uint(node.start_repetitions)?;
        let end_repetitions = source_uint(node.end_repetitions)?;

        if !accum.is_empty() {
            accum.push('|');
        }
        use std::fmt::Write as _;
        write!(
            &mut accum,
            "{start_repetitions} {end_repetitions} 2 {center_one} {outer_one} {center_one} {outer_two}"
        )
        .expect("writing a link-node property to String cannot fail");
    }

    if !accum.is_empty() {
        graph.set_prop("_molLinkNodes", accum)?;
    }
    Ok(())
}

fn append_atom_predicate(
    graph: &mut QueryGraph,
    atom: usize,
    predicate: AtomQueryPredicate,
) -> Result<(), CxQueryLoweringError> {
    // BEGIN COMPLETE QueryOps::replaceAtomWithQueryAtom
    // RDKit✔️❌: Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom) {
    // RDKit✔️❌:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❌:   PRECONDITION(atom, "bad atom");
    // RDKit✔️❌:   if (atom->hasQuery()) {
    // RDKit✔️❌:     return atom;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   QueryAtom qa(*atom);
    // RDKit✔️❌:   unsigned int idx = atom->getIdx();
    // RDKit✔️❌:
    // RDKit✔️❌:   if (atom->hasProp(common_properties::_hasMassQuery)) {
    // RDKit✔️❌:     qa.expandQuery(makeAtomMassQuery(static_cast<int>(atom->getMass())));
    // RDKit✔️❌:   }
    // RDKit✔️❌:   mol->replaceAtom(idx, &qa);
    // RDKit✔️❌:   return mol->getAtomWithIdx(idx);
    // RDKit✔️❌: }
    // END COMPLETE QueryOps::replaceAtomWithQueryAtom
    // BEGIN COMPLETE makeAtomUnsaturatedQuery
    // RDKit✔️✔️: ATOM_EQUALS_QUERY *makeAtomUnsaturatedQuery() {
    // RDKit✔️✔️:   auto *res =
    // RDKit✔️✔️:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomUnsaturated);
    // RDKit✔️✔️:   res->setDescription("AtomUnsaturated");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END COMPLETE makeAtomUnsaturatedQuery
    // BEGIN COMPLETE RWMol::replaceAtom
    // RDKit✔️❌: void RWMol::replaceAtom(unsigned int idx, Atom *atom_pin, bool,
    // RDKit✔️❌:                         bool preserveProps) {
    // RDKit✔️❌:   PRECONDITION(atom_pin, "bad atom passed to replaceAtom");
    // RDKit✔️❌:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit✔️❌:   auto atom_p = atom_pin->copy();
    // RDKit✔️❌:   atom_p->setOwningMol(this);
    // RDKit✔️❌:   atom_p->setIdx(idx);
    // RDKit✔️❌:   auto vd = boost::vertex(idx, d_graph);
    // RDKit✔️❌:   if (preserveProps) {
    // RDKit✔️❌:     const bool replaceExistingData = false;
    // RDKit✔️❌:     atom_p->updateProps(*d_graph[vd], replaceExistingData);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   const auto orig_p = d_graph[vd];
    // RDKit✔️❌:   delete orig_p;
    // RDKit✔️❌:   d_graph[vd] = atom_p;
    // RDKit✔️❌:
    // RDKit✔️❌:   // handle bookmarks
    // RDKit✔️❌:   for (auto &ab : d_atomBookmarks) {
    // RDKit✔️❌:     for (auto &elem : ab.second) {
    // RDKit✔️❌:       if (elem == orig_p) {
    // RDKit✔️❌:         elem = atom_p;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // handle stereo group
    // RDKit✔️❌:   for (auto &group : d_stereo_groups) {
    // RDKit✔️❌:     auto groupId = group.getReadId();
    // RDKit✔️❌:     auto atoms = group.getAtoms();
    // RDKit✔️❌:     auto bonds = group.getBonds();
    // RDKit✔️❌:     auto aiter = std::find(atoms.begin(), atoms.end(), orig_p);
    // RDKit✔️❌:     while (aiter != atoms.end()) {
    // RDKit✔️❌:       *aiter = atom_p;
    // RDKit✔️❌:       ++aiter;
    // RDKit✔️❌:       aiter = std::find(aiter, atoms.end(), orig_p);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     group = StereoGroup(group.getGroupType(), std::move(atoms),
    // RDKit✔️❌:                         std::move(bonds), groupId);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END COMPLETE RWMol::replaceAtom
    // BEGIN COMPLETE QueryAtom::expandQuery
    // RDKit✔️✔️: void QueryAtom::expandQuery(QUERYATOM_QUERY *what,
    // RDKit✔️✔️:                             Queries::CompositeQueryType how,
    // RDKit✔️✔️:                             bool maintainOrder) {
    // RDKit✔️✔️:   PRECONDITION(dp_query, "Can't expand empty query");
    // RDKit✔️✔️:   bool thisIsNullQuery = dp_query->getDescription() == "AtomNull";
    // RDKit✔️✔️:   bool otherIsNullQuery = what->getDescription() == "AtomNull";
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (thisIsNullQuery || otherIsNullQuery) {
    // RDKit✔️✔️:     mergeNullQueries(dp_query, thisIsNullQuery, what, otherIsNullQuery, how);
    // RDKit✔️✔️:     delete what;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   QUERYATOM_QUERY *origQ = dp_query;
    // RDKit✔️✔️:   std::string descrip;
    // RDKit✔️✔️:   switch (how) {
    // RDKit✔️✔️:     case Queries::COMPOSITE_AND:
    // RDKit✔️✔️:       dp_query = new ATOM_AND_QUERY;
    // RDKit✔️✔️:       descrip = "AtomAnd";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Queries::COMPOSITE_OR:
    // RDKit✔️✔️:       dp_query = new ATOM_OR_QUERY;
    // RDKit✔️✔️:       descrip = "AtomOr";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Queries::COMPOSITE_XOR:
    // RDKit✔️✔️:       dp_query = new ATOM_XOR_QUERY;
    // RDKit✔️✔️:       descrip = "AtomXor";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       UNDER_CONSTRUCTION("unrecognized combination query");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   dp_query->setDescription(descrip);
    // RDKit✔️✔️:   if (maintainOrder) {
    // RDKit✔️✔️:     dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(origQ));
    // RDKit✔️✔️:     dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(what));
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(what));
    // RDKit✔️✔️:     dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(origQ));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END COMPLETE QueryAtom::expandQuery
    // Source hasQuery is the stored explicit origin, never an inferred shape.
    // Ordinary carriers first acquire the native constructor query and optional
    // mass leaf; existing explicit query roots are expanded directly. Native
    // row replacement (only ordinary branch) also reconstructs every stereo
    // group, preserving read ID and resetting write ID through the canonical
    // replacement helper. Carrier members/cache/source facts remain unchanged.
    // Cost: scalar constructor has bounded leaves; ordinary replacement copies
    // typed carrier/property payload and stereo memberships. Cloning a uniform
    // carrier also copies its old non-source derived tree before discarding it,
    // an explicit overhead vs the source ordinary Atom. No graph clone/reparse.
    let source = graph
        .atom(atom)
        .ok_or(CxQueryLoweringError::AtomIndex { index: atom })?;
    if source.predicate_is_carrier_derived() {
        let source_query = crate::query_behavior::query_from_plain_atom_carrier(
            source.atomic_number(),
            source.isotope(),
            source.formal_charge(),
            source.radical_electrons(),
            || -> Result<Option<u16>, CxQueryLoweringError> {
                if source.prop("_hasMassQuery").is_none() {
                    return Ok(None);
                }
                // This conversion is confined to the actual source getMass
                // point. Without a mass query, every raw query atomic number
                // remains intact and never crosses the Element-only boundary.
                // Non-Element numbers are outside the native mass table too;
                // its explicit Atomic number not found precondition is an
                // error, never an unsupported label or guessed mass.
                let carrier = source.try_to_atom().map_err(|_| {
                    CxQueryLoweringError::InvalidGraph(format!(
                        "source mass lookup: Atomic number not found ({})",
                        source.atomic_number(),
                    ))
                })?;
                crate::query_behavior::rdkit_atom_mass(&carrier)
                    .map(|mass| Some(mass as u16))
                    .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))
            },
        )?;
        let mut incoming = source.clone();
        incoming.set_predicate(source_query);
        replace_query_atom_like_rdkit(graph, atom, &incoming, false)?;
    }
    let query_atom = graph
        .atom_mut(atom)
        .ok_or(CxQueryLoweringError::AtomIndex { index: atom })?;
    crate::query_behavior::query_atom_expand_query(
        query_atom.predicate_mut(),
        QueryNode::predicate(predicate),
        crate::query_behavior::CompositeQueryType::And,
        true,
    );
    Ok(())
}

pub(crate) fn apply_cx_query_constraint_item(
    graph: &mut QueryGraph,
    record: &CxRecord,
    item_index: usize,
) -> Result<(), CxQueryLoweringError> {
    // BEGIN COMPLETE PINNED SF192
    // RDKit✔️❌: bool parse_ring_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️❌:                       unsigned int startAtomIdx) {
    // RDKit✔️❌:   if (first >= last || *first != 'r' || first + 1 >= last ||
    // RDKit✔️❌:       *(first + 1) != 'b' || first + 2 >= last || *(first + 2) != ':') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   first += 3;
    // RDKit✔️❌:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️❌:     unsigned int n1;
    // RDKit✔️❌:     if (!read_int(first, last, n1)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // check that we can read at least two more characters:
    // RDKit✔️❌:     if (first + 1 >= last || *first != ':') {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++first;
    // RDKit✔️❌:     unsigned int n2;
    // RDKit✔️❌:     bool gt = false;
    // RDKit✔️❌:     if (*first == '*') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:       n2 = 0xDEADBEEF;
    // RDKit✔️❌:       if (VALID_ATIDX(n1)) {
    // RDKit✔️❌:         mol.setProp(common_properties::_NeedsQueryScan, 1);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       if (!read_int(first, last, n2)) {
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       switch (n2) {
    // RDKit✔️❌:         case 0:
    // RDKit✔️❌:         case 2:
    // RDKit✔️❌:         case 3:
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 4:
    // RDKit✔️❌:           gt = true;
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         default:
    // RDKit✔️❌:           BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:               << "unrecognized rb value: " << n2 << std::endl;
    // RDKit✔️❌:           return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (VALID_ATIDX(n1)) {
    // RDKit✔️❌:       auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
    // RDKit✔️❌:       if (!atom->hasQuery()) {
    // RDKit✔️❌:         atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (!gt) {
    // RDKit✔️❌:         atom->expandQuery(makeAtomRingBondCountQuery(n2),
    // RDKit✔️❌:                           Queries::COMPOSITE_AND);
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         auto q = static_cast<ATOM_EQUALS_QUERY *>(new ATOM_LESSEQUAL_QUERY);
    // RDKit✔️❌:         q->setVal(n2);
    // RDKit✔️❌:         q->setDescription("AtomRingBondCount");
    // RDKit✔️❌:         q->setDataFunc(queryAtomRingBondCount);
    // RDKit✔️❌:         atom->expandQuery(q, Queries::COMPOSITE_AND);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first < last && *first == ',') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF192
    // BEGIN COMPLETE makeAtomRingBondCountQuery
    // RDKit✔️✔️: ATOM_EQUALS_QUERY *makeAtomRingBondCountQuery(int what) {
    // RDKit✔️✔️:   ATOM_EQUALS_QUERY *res = new AtomRingQuery(what);
    // RDKit✔️✔️:   res->setDescription("AtomRingBondCount");
    // RDKit✔️✔️:   res->setDataFunc(queryAtomRingBondCount);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END COMPLETE makeAtomRingBondCountQuery
    // BEGIN COMPLETE PINNED SF191
    // RDKit✔️❌: bool parse_unsaturation(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️❌:                         unsigned int startAtomIdx) {
    // RDKit✔️❌:   if (first + 1 >= last || *first != 'u') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   if (first >= last || *first != ':') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️❌:     unsigned int idx;
    // RDKit✔️❌:     if (!read_int(first, last, idx)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (VALID_ATIDX(idx)) {
    // RDKit✔️❌:       auto atom = mol.getAtomWithIdx(idx - startAtomIdx);
    // RDKit✔️❌:       if (!atom->hasQuery()) {
    // RDKit✔️❌:         atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️❌:       }
    // RDKit✔️❌:       atom->expandQuery(makeAtomUnsaturatedQuery(), Queries::COMPOSITE_AND);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first < last && *first == ',') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF191
    // RDKit source (verbatim; one complete source record item mutates one atom):
    /*
    // RDKit✔️✔️: if (VALID_ATIDX(idx)) {
    // RDKit✔️✔️:   auto atom = mol.getAtomWithIdx(idx - startAtomIdx);
    // RDKit✔️✔️:   if (!atom->hasQuery()) {
    // RDKit✔️✔️:     atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atom->expandQuery(makeAtomUnsaturatedQuery(), Queries::COMPOSITE_AND);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:   auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
    // RDKit✔️✔️:   if (!atom->hasQuery()) {
    // RDKit✔️✔️:     atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!gt) {
    // RDKit✔️✔️:     atom->expandQuery(makeAtomRingBondCountQuery(n2),
    // RDKit✔️✔️:                       Queries::COMPOSITE_AND);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     auto q = static_cast<ATOM_EQUALS_QUERY *>(new ATOM_LESSEQUAL_QUERY);
    // RDKit✔️✔️:     q->setVal(n2);
    // RDKit✔️✔️:     q->setDescription("AtomRingBondCount");
    // RDKit✔️✔️:     q->setDataFunc(queryAtomRingBondCount);
    // RDKit✔️✔️:     atom->expandQuery(q, Queries::COMPOSITE_AND);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:   auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
    // RDKit✔️✔️:   if (!atom->hasQuery()) {
    // RDKit✔️✔️:     atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atom->expandQuery(makeAtomNonHydrogenDegreeQuery(n2),
    // RDKit✔️✔️:                     Queries::COMPOSITE_AND);
    // RDKit✔️✔️: }
     */
    // The parser's source-valid atom-index window is the detached graph range.
    // Skipping before direct indexing retains VALID_ATIDX's source behavior.
    let (atom, predicate) = match record {
        CxRecord::Unsaturation(indices) => (
            *indices.get(item_index).ok_or_else(|| {
                CxQueryLoweringError::InvalidGraph(
                    "CX progress item references a missing unsaturation index".to_owned(),
                )
            })?,
            AtomQueryPredicate::IsUnsaturated,
        ),
        CxRecord::RingBonds(constraints) => {
            let constraint = constraints.get(item_index).ok_or_else(|| {
                CxQueryLoweringError::InvalidGraph(
                    "CX progress item references a missing ring-bond constraint".to_owned(),
                )
            })?;
            let predicate = match constraint.constraint {
                CxCountConstraint::Exact(value) => AtomQueryPredicate::RingBondCount(
                    i32::try_from(value)
                        .expect("CX ring-bond equality is parser-bounded to 0, 2, or 3"),
                ),
                CxCountConstraint::LessEqual(value) => {
                    AtomQueryPredicate::RingBondCountLessEqual(value as u8)
                }
                CxCountConstraint::QueryScan => {
                    AtomQueryPredicate::RingBondCount(QUERY_SCAN_MAGIC_VALUE as i32)
                }
            };
            (constraint.atom, predicate)
        }
        CxRecord::Substitution(constraints) => {
            let CxAtomConstraint { atom, constraint } =
                constraints.get(item_index).ok_or_else(|| {
                    CxQueryLoweringError::InvalidGraph(
                        "CX progress item references a missing substitution constraint".to_owned(),
                    )
                })?;
            let predicate = match constraint {
                CxCountConstraint::Exact(value) => AtomQueryPredicate::NonHydrogenDegree(*value),
                CxCountConstraint::LessEqual(value) => {
                    AtomQueryPredicate::NonHydrogenDegreeLessEqual(*value)
                }
                CxCountConstraint::QueryScan => {
                    AtomQueryPredicate::NonHydrogenDegree(QUERY_SCAN_MAGIC_VALUE)
                }
            };
            (*atom, predicate)
        }
        _ => {
            return Err(CxQueryLoweringError::InvalidGraph(
                "CX progress item references a non-query-constraint record".to_owned(),
            ));
        }
    };
    if atom >= graph.num_atoms() {
        return Ok(());
    }
    if let CxRecord::RingBonds(constraints) = record
        && matches!(
            constraints[item_index].constraint,
            CxCountConstraint::QueryScan
        )
    {
        // RDKit✔️✔️: if (VALID_ATIDX(n1)) {
        // RDKit✔️✔️:   mol.setProp(common_properties::_NeedsQueryScan, 1);
        // RDKit✔️✔️: }
        // Source ordinary int1 is assigned before query conversion/expansion;
        // preserve it if the optional source mass lookup or replacement fails.
        graph
            .set_prop("_NeedsQueryScan", 1_i32)
            .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
    }
    append_atom_predicate(graph, atom, predicate)
}

const CX_LABELS_PROCESSED_PROP: &str = "_cxsmilesLabelsProcessed";

/// Apply RDKit's deferred CX label replacement to the detached query graph.
/// The temporary guard remains set until the enclosing CX application calls
/// `finish_cx_smiles_labels`, so SGroup helpers can invoke this before attach.
pub(crate) fn process_cx_smiles_labels(graph: &mut QueryGraph) -> Result<(), CxQueryLoweringError> {
    // RDKit✔️❌: void processCXSmilesLabels(RWMol &mol) {
    // RDKit✔️❌:   if (mol.hasProp("_cxsmilesLabelsProcessed")) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     std::string symb = "";
    // RDKit✔️❌:     if (atom->getPropIfPresent(common_properties::atomLabel, symb)) {
    // RDKit✔️❌:       atom->clearProp(common_properties::dummyLabel);
    // RDKit✔️❌:       if (symb == "star_e") {
    // RDKit✔️❌:         /* according to the MDL spec, these match anything, but in MARVIN they
    // RDKit✔️❌:         are "unspecified end groups" for polymers */
    // RDKit✔️❌:         addquery(makeAtomNullQuery(), symb, mol, atom->getIdx());
    // RDKit✔️❌:       } else if (symb == "Q_e") {
    // RDKit✔️❌:         addquery(makeQAtomQuery(), symb, mol, atom->getIdx());
    // RDKit✔️❌:       } else if (symb == "QH_p") {
    // RDKit✔️❌:         addquery(makeQHAtomQuery(), symb, mol, atom->getIdx());
    // RDKit✔️❌:       } else if (symb == "AH_p") {  // this seems wrong...
    // RDKit✔️❌:         /* According to the MARVIN Sketch, AH is "any atom, including H" -
    // RDKit✔️❌:         this would be "*" in SMILES - and "A" is "any atom except H".
    // RDKit✔️❌:         The CXSMILES docs say that "A" can be represented normally in SMILES
    // RDKit✔️❌:         and that "AH" needs to be written out as AH_p. I'm going to assume that
    // RDKit✔️❌:         this is a Marvin internal thing and just parse it as they describe it.
    // RDKit✔️❌:         This means that "*" in the SMILES itself needs to be treated
    // RDKit✔️❌:         differently, which we do below. */
    // RDKit✔️❌:         addquery(makeAHAtomQuery(), symb, mol, atom->getIdx());
    // RDKit✔️❌:       } else if (symb == "X_p") {
    // RDKit✔️❌:         addquery(makeXAtomQuery(), symb, mol, atom->getIdx());
    // RDKit✔️❌:       } else if (symb == "XH_p") {
    // RDKit✔️❌:         addquery(makeXHAtomQuery(), symb, mol, atom->getIdx());
    // RDKit✔️❌:       } else if (symb == "M_p") {
    // RDKit✔️❌:         addquery(makeMAtomQuery(), symb, mol, atom->getIdx());
    // RDKit✔️❌:       } else if (symb == "MH_p") {
    // RDKit✔️❌:         addquery(makeMHAtomQuery(), symb, mol, atom->getIdx());
    // RDKit✔️❌:       } else if (std::find(pseudoatoms_p.begin(), pseudoatoms_p.end(), symb) !=
    // RDKit✔️❌:                  pseudoatoms_p.end()) {
    // RDKit✔️❌:         // strip off the "_p":
    // RDKit✔️❌:         atom->setProp(common_properties::dummyLabel,
    // RDKit✔️❌:                       symb.substr(0, symb.size() - 2));
    // RDKit✔️❌:         atom->clearProp(common_properties::atomLabel);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else if (atom->getAtomicNum() == 0 && !atom->hasQuery() &&
    // RDKit✔️❌:                !atom->getIsotope() && atom->getSymbol() == "*") {
    // RDKit✔️❌:       addquery(makeAAtomQuery(), "", mol, atom->getIdx());
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   mol.setProp("_cxsmilesLabelsProcessed", 1, true);
    // RDKit✔️❌: }
    // RDKit✔️❌: std::string Atom::getSymbol() const {
    // RDKit✔️❌:   std::string res;
    // RDKit✔️❌:   // handle dummies differently:
    // RDKit✔️❌:   if (d_atomicNum != 0 ||
    // RDKit✔️❌:       !getPropIfPresent<std::string>(common_properties::dummyLabel, res)) {
    // RDKit✔️❌:     res = PeriodicTable::getTable()->getElementSymbol(d_atomicNum);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // Behavior: preserve guard presence and ordered atom visitation. Label
    // conversion occurs before clearing dummyLabel; all factories and pseudo
    // label stripping retain the source dispatch. The unlabeled branch short
    // circuits atomic number, hasQuery and isotope before getSymbol, including
    // a typed dummyLabel property. Successful completion alone sets the computed
    // guard; earlier atom mutations remain visible if a later conversion fails.
    // Complexity: O(V + copied query/properties + stereo memberships), with
    // source-required allocations; canonical property-tree allocation costs
    // differ from Dict's contiguous rows, hence the separate loss marker.
    // These predicates reuse the existing source-anchored QueryOps factories
    // in query_behavior; this dispatcher does not build parallel query trees.
    if graph.prop(CX_LABELS_PROCESSED_PROP).is_some() {
        return Ok(());
    }

    for atom_index in 0..graph.num_atoms() {
        let label = graph
            .atom(atom_index)
            .and_then(|atom| atom.prop("atomLabel"))
            .map(cosmolkit_core::property_value_to_string)
            .transpose()?;

        if let Some(label) = label {
            graph
                .atom_mut(atom_index)
                .ok_or(CxQueryLoweringError::AtomIndex { index: atom_index })?
                .clear_prop("dummyLabel")?;

            let predicate = match label.as_bytes() {
                b"star_e" => Some(crate::query_behavior::make_atom_null_query()),
                b"Q_e" => Some(crate::query_behavior::make_q_atom_query()),
                b"QH_p" => Some(crate::query_behavior::make_q_h_atom_query()),
                b"AH_p" => Some(crate::query_behavior::make_a_h_atom_query()),
                b"X_p" => Some(crate::query_behavior::make_x_atom_query()),
                b"XH_p" => Some(crate::query_behavior::make_x_h_atom_query()),
                b"M_p" => Some(crate::query_behavior::make_m_atom_query()),
                b"MH_p" => Some(crate::query_behavior::make_m_h_atom_query()),
                _ => None,
            };
            if let Some(predicate) = predicate {
                add_query_like_rdkit(graph, atom_index, predicate, Some(&label))?;
            } else if let Some(dummy_label) = match label.as_bytes() {
                b"Pol_p" => Some("Pol"),
                b"Mod_p" => Some("Mod"),
                _ => None,
            } {
                let atom = graph
                    .atom_mut(atom_index)
                    .ok_or(CxQueryLoweringError::AtomIndex { index: atom_index })?;
                atom.set_prop("dummyLabel", dummy_label)?;
                atom.clear_prop("atomLabel")?;
            }
        } else {
            let atom = &graph.atoms()[atom_index];
            // Atom::getSymbol is reached only after these source guards.
            if atom.atomic_number() == 0
                && atom.predicate_is_carrier_derived()
                && atom.isotope().unwrap_or(0) == 0
            {
                let symbol = atom
                    .prop("dummyLabel")
                    .map(cosmolkit_core::property_value_to_string)
                    .transpose()?
                    .unwrap_or_else(|| cosmolkit_model::PropertyText::from("*"));
                if symbol.as_bytes() == b"*" {
                    add_query_like_rdkit(
                        graph,
                        atom_index,
                        crate::query_behavior::make_a_atom_query(),
                        None,
                    )?;
                }
            }
        }
    }

    // Local complexity review: this loop visits each atom once. A query-label
    // replacement copies that atom's property map and reconstructs the same
    // stereo-group members that RWMol::replaceAtom scans; both costs are
    // linear in the copied properties and affected group membership.
    graph.set_computed_prop(CX_LABELS_PROCESSED_PROP, 1_i32)?;
    Ok(())
}

/// Finish one complete CX application, matching parseCXExtensions' guard
/// cleanup after its final label pass.
pub(crate) fn finish_cx_smiles_labels(graph: &mut QueryGraph) -> Result<(), CxQueryLoweringError> {
    process_cx_smiles_labels(graph)?;
    graph.clear_prop(CX_LABELS_PROCESSED_PROP)?;
    Ok(())
}

pub(crate) fn apply_cx_data_sgroup_to_query(
    graph: &mut QueryGraph,
    data: &CxDataSGroup,
    cx_sequence_id: u32,
) -> Result<(), CxQueryLoweringError> {
    // RDKit source (verbatim; source atom filtering and local group mutation):
    /*
      SubstanceGroup sgroup(&mol, std::string("DAT"));
      sgroup.setProp(cxsmilesindex, nSGroups);
      bool keepSGroup = false;
      for (auto idx : atoms) {
        if (VALID_ATIDX(idx)) {
          keepSGroup = true;
          sgroup.addAtomWithIdx(idx - startAtomIdx);
        }
      }
      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDNAME");
      if (keepSGroup) {
        sgroup.setProp("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
      }
      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "DATAFIELDS", true);
      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "QUERYOP");
      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDINFO");
      parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDTAG");
      if (first < last && *first == '(') {
        std::string coords = read_text_to(first, last, ")");
        ++first;
        if (keepSGroup) {
          sgroup.setProp("COORDS", coords);
        }
      }
      if (keepSGroup) {
        processCXSmilesLabels(mol);
        sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
        addSubstanceGroup(mol, sgroup);
      }
    */
    // RDKit❗✔️: preserve valid input atom order and duplicates, process labels
    // before attachment, and keep CX sequence ID separate from dense storage ID.
    let atoms = data
        .atoms
        .iter()
        .filter(|&&index| index < graph.num_atoms())
        .map(|&index| AtomId::new(index))
        .collect::<Vec<_>>();
    if atoms.is_empty() {
        return Ok(());
    }

    process_cx_smiles_labels(graph)?;

    let mut groups = query_substance_groups(graph).to_vec();
    let dense_id = SubstanceGroupId::new(groups.len());
    let source_index = groups.len() + 1;
    let typed_data = SGroupData {
        field_name: (!data.field_name.is_empty()).then(|| data.field_name.clone()),
        field_info: (!data.field_info.is_empty()).then(|| data.field_info.clone()),
        field_display: Some("    0.0000    0.0000    DR    ALL  0       0".into()),
        query_op: (!data.query_op.is_empty()).then(|| data.query_op.clone()),
        values: (!data.data.is_empty())
            .then(|| vec![data.data.clone()])
            .unwrap_or_default(),
        ..SGroupData::default()
    };
    let mut group = SubstanceGroup::new(dense_id, SubstanceGroupKind::Data)
        .with_rdkit_sequence_id(cx_sequence_id)
        .with_atoms(atoms)
        .with_data(typed_data);
    // BEGIN COMPLETE PINNED SF194 property-field writes
    // RDKit✔️❌: void parse_data_sgroup_attr(Iterator &first, Iterator last,
    // RDKit✔️❌:                             SubstanceGroup &sgroup, bool keepSGroup,
    // RDKit✔️❌:                             std::string fieldName, bool fieldIsArray = false) {
    // RDKit✔️❌:   PRECONDITION(first < last, "parse_data_sgroup_attr: first >= last");
    // RDKit✔️❌:   if (first != last && *first != '|') {
    // RDKit✔️❌:     std::string data = read_text_to(first, last, ":");
    // RDKit✔️❌:     ++first;
    // RDKit✔️❌:     if (!data.empty() && keepSGroup) {
    // RDKit✔️❌:       if (fieldIsArray) {
    // RDKit✔️❌:         std::vector<std::string> dataFields = {data};
    // RDKit✔️❌:         sgroup.setProp(fieldName, dataFields);
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         sgroup.setProp(fieldName, data);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF194 property-field writes
    // BEGIN COMPLETE SubstanceGroup::SubstanceGroup(TYPE)
    // RDKit✔️❌: SubstanceGroup::SubstanceGroup(ROMol *owning_mol, const std::string &type)
    // RDKit✔️❌:     : RDProps(), dp_mol(owning_mol) {
    // RDKit✔️❌:   PRECONDITION(owning_mol, "supplied owning molecule is bad");
    // RDKit✔️❌:
    // RDKit✔️❌:   // TYPE is required to be set , as other properties will depend on it.
    // RDKit✔️❌:   setProp<std::string>("TYPE", type);
    // RDKit✔️❌: }
    // END COMPLETE SubstanceGroup::SubstanceGroup(TYPE)
    // Source property insertion order is retained by the sole store: TYPE,
    // sequence index, FIELDNAME when nonempty, FIELDDISP, remaining fields,
    // optional COORDS, then dense source index at helper completion.
    group.set_prop("TYPE", "DAT")?;
    group.set_prop("_cxsmilesindex", cx_sequence_id)?;
    if !data.field_name.is_empty() {
        group.set_prop("FIELDNAME", data.field_name.clone())?;
    }
    group.set_prop("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0")?;
    if !data.data.is_empty() {
        group.set_prop("DATAFIELDS", vec![data.data.clone()])?;
        group.push_data_field(data.data.clone());
    }
    if !data.query_op.is_empty() {
        group.set_prop("QUERYOP", data.query_op.clone())?;
    }
    if !data.field_info.is_empty() {
        group.set_prop("FIELDINFO", data.field_info.clone())?;
    }
    if !data.field_tag.is_empty() {
        group.set_prop("FIELDTAG", data.field_tag.clone())?;
    }
    if let Some(coordinates) = &data.coordinates {
        group.set_prop("COORDS", coordinates.clone())?;
    }
    group.set_prop("index", source_index as u32)?;
    groups.push(group);
    replace_query_substance_groups(graph, groups)
        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))
}

pub(crate) fn validate_cx_variable_attachment_atom_to_query(
    graph: &QueryGraph,
    attachment: &CxVariableAttachment,
) -> Result<(), CxQueryLoweringError> {
    // RDKit source (verbatim; CXSmilesOps.cpp::parse_variable_attachments):
    // RDKit❗✔️:     if (VALID_ATIDX(at1idx) &&
    // RDKit❗✔️:         mol.getAtomWithIdx(at1idx - startAtomIdx)->getDegree() != 1) {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:           << "position variation bond to atom with more than one bond"
    // RDKit❗✔️:           << std::endl;
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    if attachment.atom >= graph.num_atoms() {
        return Ok(());
    }
    let degree = graph.adjacency().get(attachment.atom).map_or(0, Vec::len);
    if degree != 1 {
        return Err(CxQueryLoweringError::InvalidGraph(
            "position variation bond to atom with more than one bond".to_owned(),
        ));
    }
    Ok(())
}

pub(crate) fn apply_cx_variable_attachment_effect_to_query(
    graph: &mut QueryGraph,
    attachment: &CxVariableAttachment,
) -> Result<(), CxQueryLoweringError> {
    // RDKit source (verbatim; CXSmilesOps.cpp::parse_variable_attachments):
    // RDKit❗✔️:       if (VALID_ATIDX(aidx)) {
    // RDKit❗✔️:         others.push_back(std::to_string(aidx - startAtomIdx + 1));
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (VALID_ATIDX(at1idx)) {
    // RDKit❗✔️:         std::string endPts = "(" + std::to_string(others.size());
    // RDKit❗✔️:         for (auto idx : others) {
    // RDKit❗✔️:           endPts += " " + idx;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         endPts += ")";
    // RDKit❗✔️:         for (auto nbri : boost::make_iterator_range(
    // RDKit❗✔️:                  mol.getAtomBonds(mol.getAtomWithIdx(at1idx - startAtomIdx)))) {
    // RDKit❗✔️:           auto bnd = mol[nbri];
    // RDKit❗✔️:           bnd->setProp(common_properties::_MolFileBondEndPts, endPts);
    // RDKit❗✔️:           bnd->setProp(common_properties::_MolFileBondAttach, std::string("ANY"));
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    if attachment.atom >= graph.num_atoms() {
        return Ok(());
    }
    let endpoints = attachment
        .endpoints
        .iter()
        .filter(|&&index| index < graph.num_atoms())
        .map(|&index| {
            u32::try_from(index)
                .map(|index| index.wrapping_add(1).to_string())
                .map_err(|_| {
                    CxQueryLoweringError::InvalidGraph(
                        "CX variable-attachment endpoint exceeds the source unsigned-int domain"
                            .to_owned(),
                    )
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut end_points = format!("({}", endpoints.len());
    for endpoint in endpoints {
        end_points.push(' ');
        end_points.push_str(&endpoint);
    }
    end_points.push(')');

    let degree = graph.adjacency().get(attachment.atom).map_or(0, Vec::len);
    for position in 0..degree {
        let bond_index = graph
            .adjacency()
            .get(attachment.atom)
            .and_then(|neighbors| neighbors.get(position))
            .map(|neighbor| neighbor.1)
            .ok_or_else(|| {
                CxQueryLoweringError::InvalidGraph(
                    "query adjacency changed during CX variable attachment".to_owned(),
                )
            })?;
        let bond = graph
            .bonds_mut()
            .get_mut(bond_index)
            .ok_or(CxQueryLoweringError::BondIndex { index: bond_index })?;
        bond.bond_mut()
            .set_prop("_MolFileBondEndPts", end_points.clone())?;
        bond.bond_mut().set_prop("_MolFileBondAttach", "ANY")?;
    }
    Ok(())
}

pub(crate) fn apply_cx_variable_attachment_to_query(
    graph: &mut QueryGraph,
    attachment: &CxVariableAttachment,
) -> Result<(), CxQueryLoweringError> {
    validate_cx_variable_attachment_atom_to_query(graph, attachment)?;
    apply_cx_variable_attachment_effect_to_query(graph, attachment)
}

pub(crate) fn apply_cx_wedge_bond_to_query(
    graph: &mut QueryGraph,
    wedge: &CxWedgeBond,
) -> Result<(), CxQueryLoweringError> {
    // RDKit source (verbatim; CXSmilesOps.cpp::parse_wedged_bonds):
    // RDKit❗✔️: template <typename Iterator>
    // RDKit❗✔️: bool parse_wedged_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit❗✔️:                         unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit❗✔️:   // these look like: CC(O)Cl |w:1.0|
    // RDKit❗✔️:   // also wD and wU for down and up wedges.
    // RDKit❗✔️:   //
    // RDKit❗✔️:   // We do not end up using this to set stereochemistry, but the relevant bond
    // RDKit❗✔️:   // properties are set in case client code wants to do something with the
    // RDKit❗✔️:   // information.
    // RDKit❗✔️:   if (first >= last || *first != 'w' || first + 1 >= last) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ++first;
    // RDKit❗✔️:   Bond::BondDir state = Bond::BondDir::NONE;
    // RDKit❗✔️:   unsigned int cfg = 0;
    // RDKit❗✔️:   switch (*first) {
    // RDKit❗✔️:     case ':':
    // RDKit❗✔️:       state = Bond::BondDir::UNKNOWN;
    // RDKit❗✔️:       cfg = 2;
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 'U':
    // RDKit❗✔️:       state = Bond::BondDir::BEGINWEDGE;
    // RDKit❗✔️:       cfg = 1;
    // RDKit❗✔️:       ++first;
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 'D':
    // RDKit❗✔️:       state = Bond::BondDir::BEGINDASH;
    // RDKit❗✔️:       cfg = 3;
    // RDKit❗✔️:       ++first;
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       break;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (state == Bond::BondDir::NONE || first >= last || first + 1 >= last ||
    // RDKit❗✔️:       *first != ':') {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ++first;
    // RDKit❗✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit❗✔️:     unsigned int atomIdx;
    // RDKit❗✔️:     if (!read_int(first, last, atomIdx)) {
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (first < last && *first == '.') {
    // RDKit❗✔️:       ++first;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "improperly formatted w block" << std::endl;
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     unsigned int bondIdx;
    // RDKit❗✔️:     if (!read_int(first, last, bondIdx)) {
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (VALID_ATIDX(atomIdx) && VALID_BNDIDX(bondIdx)) {
    // RDKit❗✔️:       auto atom = mol.getAtomWithIdx(atomIdx - startAtomIdx);
    // RDKit❗✔️:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit❗✔️:       if (!bond) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "bond " << bondIdx << " not found, wedge from atom " << atomIdx
    // RDKit❗✔️:             << " cannot be applied." << std::endl;
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (bond->hasProp(common_properties::_MolFileBondCfg)) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "w block attempts to set wedging on bond " << bond->getIdx()
    // RDKit❗✔️:             << " more than once." << std::endl;
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (atom->getIdx() != bond->getBeginAtomIdx()) {
    // RDKit❗✔️:         if (atom->getIdx() != bond->getEndAtomIdx()) {
    // RDKit❗✔️:           BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:               << "atom " << atomIdx << " is not associated with bond "
    // RDKit❗✔️:               << bondIdx << "(" << bond->getBeginAtomIdx() + startAtomIdx << "-"
    // RDKit❗✔️:               << bond->getEndAtomIdx() + startAtomIdx << ")"
    // RDKit❗✔️:               << " in w block" << std::endl;
    // RDKit❗✔️:           return false;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         auto eidx = bond->getBeginAtomIdx();
    // RDKit❗✔️:         bond->setBeginAtomIdx(atom->getIdx());
    // RDKit❗✔️:         bond->setEndAtomIdx(eidx);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       bond->setProp(common_properties::_MolFileBondCfg, cfg);
    // RDKit❗✔️:       bond->setBondDir(state);
    // RDKit❗✔️:       if (cfg == 2 && canHaveDirection(*bond)) {
    // RDKit❗✔️:         bond->getBeginAtom()->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit❗✔️:         mol.setProp(detail::_needsDetectBondStereo, 1);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if ((cfg == 1 || cfg == 3) && canHaveDirection(*bond)) {
    // RDKit❗✔️:         mol.setProp(detail::_needsDetectAtomStereo, 1);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (first < last && *first == ',') {
    // RDKit❗✔️:       ++first;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // RDKit source helper (verbatim; CXSmilesOps.cpp::get_bond_with_smiles_idx):
    // RDKit❗✔️: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
    // RDKit❗✔️:   for (auto bnd : mol.bonds()) {
    // RDKit❗✔️:     unsigned int smilesIdx;
    // RDKit❗✔️:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
    // RDKit❗✔️:         smilesIdx == idx) {
    // RDKit❗✔️:       return bnd;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return nullptr;
    // RDKit❗✔️: }
    // RDKit source helper (verbatim; Bond.h::canHaveDirection):
    // RDKit❗✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit❗✔️:   auto bondType = bond.getBondType();
    // RDKit❗✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit❗✔️: }
    // Source grammar indices and final bond rows are distinct. Lookup scans
    // the preserved scalar properties in final row order, O(E) as in RDKit.
    if wedge.atom >= graph.num_atoms() || wedge.bond >= graph.num_bonds() {
        return Ok(());
    }
    let row = query_bond_row_from_source_index(graph, wedge.bond)?;
    let bond = graph
        .bonds_mut()
        .get_mut(row)
        .ok_or(CxQueryLoweringError::BondIndex { index: wedge.bond })?;
    if bond.bond().prop("_MolFileBondCfg").is_some() {
        return Err(CxQueryLoweringError::InvalidGraph(format!(
            "w block attempts to set wedging on bond {} more than once.",
            bond.id().index()
        )));
    }
    let atom = AtomId::new(wedge.atom);
    if bond.begin() != atom && bond.end() != atom {
        return Err(CxQueryLoweringError::WedgeAtomMismatch {
            atom: wedge.atom,
            bond: wedge.bond,
        });
    }
    let can_have_direction = matches!(bond.bond().order(), BondOrder::Single | BondOrder::Aromatic);
    if bond.begin() != atom {
        let previous_begin = bond.begin();
        bond.bond_mut().set_endpoints(atom, previous_begin);
    }
    let (configuration, direction) = match wedge.direction {
        CxWedgeDirection::Unknown => (2_u32, BondDirection::Unknown),
        CxWedgeDirection::BeginWedge => (1_u32, BondDirection::BeginWedge),
        CxWedgeDirection::BeginDash => (3_u32, BondDirection::BeginDash),
    };
    bond.bond_mut().set_prop("_MolFileBondCfg", configuration)?;
    bond.bond_mut().set_direction(direction);
    if wedge.direction == CxWedgeDirection::Unknown && can_have_direction {
        graph
            .atom_mut(wedge.atom)
            .ok_or(CxQueryLoweringError::AtomIndex { index: wedge.atom })?
            .set_chiral_tag(ChiralTag::Unspecified);
        graph.set_prop("_needsDetectBondStereo", 1_i32)?;
    }
    if matches!(
        wedge.direction,
        CxWedgeDirection::BeginWedge | CxWedgeDirection::BeginDash
    ) && can_have_direction
    {
        graph.set_prop("_needsDetectAtomStereo", 1_i32)?;
    }
    Ok(())
}

pub(crate) fn apply_cx_double_bond_stereo_to_query(
    graph: &mut QueryGraph,
    bond_index: usize,
    stereo: CxDoubleBondStereoKind,
) -> Result<(), CxQueryLoweringError> {
    // RDKit source (verbatim; CXSmilesOps.cpp::parse_doublebond_stereo):
    // RDKit❗✔️: template <typename Iterator>
    // RDKit❗✔️: bool parse_doublebond_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit❗✔️:                              unsigned int, unsigned int startBondIdx,
    // RDKit❗✔️:                              Bond::BondStereo stereo) {
    // RDKit❗✔️:   // these look like: C1CCCC/C=C/CCC1 |ctu:5|
    // RDKit❗✔️:   // also c and t for cis or trans
    // RDKit❗✔️:   while (first < last && *first != ':') {
    // RDKit❗✔️:     ++first;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (first >= last || *first != ':') {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ++first;
    // RDKit❗✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit❗✔️:     unsigned int bondIdx;
    // RDKit❗✔️:     if (!read_int(first, last, bondIdx)) {
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (VALID_BNDIDX(bondIdx)) {
    // RDKit❗✔️:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit❗✔️:       if (!bond) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "bond " << bondIdx
    // RDKit❗✔️:             << " not found, cannot mark as stereo double bond." << std::endl;
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       bool useCXOrdering = true;
    // RDKit❗✔️:       Chirality::detail::setStereoForBond(mol, bond, stereo, useCXOrdering);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (first < last && *first == ',') {
    // RDKit❗✔️:       ++first;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // RDKit source helper (verbatim; Chirality.cpp::setStereoForBond):
    // RDKit❗✔️: void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo,
    // RDKit❗✔️:                       bool useCXSmilesOrdering) {
    // RDKit❗✔️:   // NOTE:  moved from parse_doublebond_stereo CXSmilesOps
    // RDKit❗✔️:   // IF useCXSmilesOrdering is true, the cis/trans/unknown marker will be
    // RDKit❗✔️:   // assigned relative to the lowest-numbered neighbor of each double bond atom.
    // RDKit❗✔️:   // Otherwise it uses the lowest-numbered neighbor on the lower-numbered atom
    // RDKit❗✔️:   // of the double bond and the highest-numbered neighbor on the higher-numbered
    // RDKit❗✔️:   // atom
    // RDKit❗✔️:   auto begAtom = bond->getBeginAtom();
    // RDKit❗✔️:   auto endAtom = bond->getEndAtom();
    // RDKit❗✔️:   if (begAtom->getIdx() > endAtom->getIdx()) {
    // RDKit❗✔️:     std::swap(begAtom, endAtom);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    // RDKit❗✔️:     unsigned int begControl = mol.getNumAtoms();
    // RDKit❗✔️:     for (auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit❗✔️:       if (nbr == endAtom) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       begControl = std::min(nbr->getIdx(), begControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
    // RDKit❗✔️:     for (auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit❗✔️:       if (nbr == begAtom) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
    // RDKit❗✔️:                                        : std::max(nbr->getIdx(), endControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (begAtom != bond->getBeginAtom()) {
    // RDKit❗✔️:       std::swap(begControl, endControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     bond->setStereoAtoms(begControl, endControl);
    // RDKit❗✔️:     bond->setStereo(stereo);
    // RDKit❗✔️:     mol.setProp("_needsDetectBondStereo", 1);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Source property lookup scans final rows, O(E) as in RDKit;
    // adjacency iteration retains the same source control-atom minima.
    if bond_index >= graph.num_bonds() {
        return Ok(());
    }
    let bond_index = query_bond_row_from_source_index(graph, bond_index)?;
    let (begin, end, begin_degree, end_degree) = {
        let bond = graph
            .bonds()
            .get(bond_index)
            .ok_or(CxQueryLoweringError::BondIndex { index: bond_index })?;
        let begin = bond.begin();
        let end = bond.end();
        let begin_degree = graph.adjacency().get(begin.index()).map_or(0, Vec::len);
        let end_degree = graph.adjacency().get(end.index()).map_or(0, Vec::len);
        (begin, end, begin_degree, end_degree)
    };
    if begin_degree <= 1 || end_degree <= 1 {
        return Ok(());
    }
    let (low, high) = if begin.index() <= end.index() {
        (begin, end)
    } else {
        (end, begin)
    };
    let find_control = |atom: AtomId, other: AtomId| {
        graph
            .adjacency()
            .get(atom.index())
            .into_iter()
            .flatten()
            .filter_map(|&(neighbor, _)| (neighbor != other.index()).then_some(neighbor))
            .min()
            .unwrap_or(graph.num_atoms())
    };
    let mut begin_control = find_control(low, high);
    let mut end_control = find_control(high, low);
    if low != begin {
        std::mem::swap(&mut begin_control, &mut end_control);
    }
    let value = match stereo {
        CxDoubleBondStereoKind::Any => BondStereo::Any,
        CxDoubleBondStereoKind::Cis => BondStereo::Cis,
        CxDoubleBondStereoKind::Trans => BondStereo::Trans,
    };
    let bond = graph
        .bonds_mut()
        .get_mut(bond_index)
        .ok_or(CxQueryLoweringError::BondIndex { index: bond_index })?;
    bond.bond_mut()
        .set_stereo_atoms(Some([AtomId::new(begin_control), AtomId::new(end_control)]));
    bond.bond_mut()
        .set_stereo(value)
        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
    graph.set_prop("_needsDetectBondStereo", 1_i32)?;
    Ok(())
}

pub(crate) fn apply_cx_polymer_sgroup_to_query(
    graph: &mut QueryGraph,
    polymer: &CxPolymerSGroup,
    cx_sequence_id: u32,
) -> Result<(), CxQueryLoweringError> {
    // RDKit source (verbatim; CXSmilesOps.cpp::sgroupTypemap):
    // RDKit✔️✔️: const std::map<std::string, std::string> sgroupTypemap = {
    // RDKit✔️✔️:     {"n", "SRU"},   {"mon", "MON"}, {"mer", "MER"}, {"co", "COP"},
    // RDKit✔️✔️:     {"xl", "CRO"},  {"mod", "MOD"}, {"mix", "MIX"}, {"f", "FOR"},
    // RDKit✔️✔️:     {"any", "ANY"}, {"gen", "GEN"}, {"c", "COM"},   {"grf", "GRA"},
    // RDKit✔️✔️:     {"alt", "COP"}, {"ran", "COP"}, {"blk", "COP"}};
    // RDKit source (verbatim; CXSmilesOps.cpp::parse_polymer_sgroup):
    /*
    bool keepSGroup = false;
    SubstanceGroup sgroup(&mol, type->second);
    sgroup.setProp(cxsmilesindex, nSGroups);
    if (type_code == "alt") {
      sgroup.setProp("SUBTYPE", std::string("ALT"));
    } else if (type_code == "ran") {
      sgroup.setProp("SUBTYPE", std::string("RAN"));
    } else if (type_code == "blk") {
      sgroup.setProp("SUBTYPE", std::string("BLO"));
    }
    for (auto idx : atoms) {
      if (VALID_ATIDX(idx)) {
        sgroup.addAtomWithIdx(idx - startAtomIdx);
        keepSGroup = true;
      }
    }
    if (keepSGroup) {
      processCXSmilesLabels(mol);
      finalizePolymerSGroup(mol, sgroup);
      sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
      addSubstanceGroup(mol, sgroup);
    }
    */
    // RDKit✔️❌: valid atom occurrences retain source order and duplicates; a
    // group with no valid atoms is skipped before labels or crossings run.
    let (kind, source_type) = match polymer.type_code.as_bytes() {
        b"n" => (SubstanceGroupKind::StructuralRepeatUnit, "SRU"),
        b"mon" => (SubstanceGroupKind::Monomer, "MON"),
        b"mer" => (SubstanceGroupKind::Mer, "MER"),
        b"co" => (SubstanceGroupKind::Copolymer, "COP"),
        b"xl" => (SubstanceGroupKind::Crosslink, "CRO"),
        b"mod" => (SubstanceGroupKind::Modification, "MOD"),
        b"mix" => (SubstanceGroupKind::MixtureComponent, "MIX"),
        b"f" => (SubstanceGroupKind::Formulation, "FOR"),
        b"any" => (SubstanceGroupKind::AnyPolymer, "ANY"),
        b"gen" => (SubstanceGroupKind::Generic("GEN".into()), "GEN"),
        b"c" => (SubstanceGroupKind::Generic("COM".into()), "COM"),
        b"grf" => (SubstanceGroupKind::Graft, "GRA"),
        b"alt" | b"ran" | b"blk" => (SubstanceGroupKind::Copolymer, "COP"),
        _ => {
            return Err(CxQueryLoweringError::InvalidGraph(
                "unknown CX polymer SGroup type".to_owned(),
            ));
        }
    };
    let atoms = polymer
        .atoms
        .iter()
        .filter(|&&index| index < graph.num_atoms())
        .map(|&index| AtomId::new(index))
        .collect::<Vec<_>>();
    if atoms.is_empty() {
        return Ok(());
    }

    // RDKit❗✔️: an explicit crossing outside VALID_ATIDX skips the complete
    // local SGroup; a source-valid value is then checked as a bond by core.
    if polymer
        .head_crossings
        .iter()
        .chain(&polymer.tail_crossings)
        .any(|&index| index >= graph.num_atoms())
    {
        return Ok(());
    }
    let head = polymer
        .head_crossings
        .iter()
        .copied()
        .map(BondId::new)
        .collect::<Vec<_>>();
    let tail = polymer
        .tail_crossings
        .iter()
        .copied()
        .map(BondId::new)
        .collect::<Vec<_>>();

    // Source calls processCXSmilesLabels before finalizing or attaching this
    // local group, so a later finalizer error retains label mutations only.
    process_cx_smiles_labels(graph)?;

    let mut groups = query_substance_groups(graph).to_vec();
    let dense_id = SubstanceGroupId::new(groups.len());
    let mut group = SubstanceGroup::new(dense_id, kind)
        .with_rdkit_sequence_id(cx_sequence_id)
        .with_atoms(atoms);
    // RDKit✔️❌: SubstanceGroup::SubstanceGroup(ROMol *owning_mol, const std::string &type)
    // RDKit✔️❌:     : RDProps(), dp_mol(owning_mol) {
    // RDKit✔️❌:   PRECONDITION(owning_mol, "supplied owning molecule is bad");
    // RDKit✔️❌:
    // RDKit✔️❌:   // TYPE is required to be set , as other properties will depend on it.
    // RDKit✔️❌:   setProp<std::string>("TYPE", type);
    // RDKit✔️❌: }
    // Ordinary TYPE assignment uses the source String tag; the sole
    // store preserves order, with its known extra key allocation cost.
    group.set_prop("TYPE", source_type)?;
    group.set_prop("_cxsmilesindex", cx_sequence_id)?;
    match polymer.type_code.as_bytes() {
        b"alt" => {
            group.set_subtype("ALT");
            group.set_prop("SUBTYPE", "ALT")?;
        }
        b"ran" => {
            group.set_subtype("RAN");
            group.set_prop("SUBTYPE", "RAN")?;
        }
        b"blk" => {
            group.set_subtype("BLO");
            group.set_prop("SUBTYPE", "BLO")?;
        }
        _ => {}
    }
    if !polymer.label.is_empty() {
        group.set_label(polymer.label.clone());
        group.set_prop("LABEL", polymer.label.clone())?;
    }
    if !polymer.connect.is_empty() {
        group.set_prop("CONNECT", polymer.connect.clone())?;
    }
    cosmolkit_core::finalize_polymer_sgroup(
        &mut group,
        (!polymer.connect.is_empty()).then_some(polymer.connect.as_bytes()),
        &head,
        &tail,
        graph.num_atoms(),
        graph.num_bonds(),
        |atom| {
            graph
                .adjacency()
                .get(atom.index())
                .into_iter()
                .flatten()
                .map(|&(neighbor, bond)| (AtomId::new(neighbor), BondId::new(bond)))
        },
    )
    .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;

    // RDKit sets the one-based source index only after finalization succeeds.
    let source_index = groups.len() + 1;
    group.set_prop("index", source_index as u32)?;
    groups.push(group);
    replace_query_substance_groups(graph, groups)
        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))
}

pub(crate) fn resolve_cx_sgroup_hierarchy_parent(
    groups: &[SubstanceGroup],
    cx_parent_id: usize,
) -> Result<Option<(SubstanceGroupId, u32)>, CxQueryLoweringError> {
    // RDKit source (verbatim; find_matching_sgroup from CXSmilesOps.cpp):
    // RDKit✔️✔️: std::vector<RDKit::SubstanceGroup>::iterator find_matching_sgroup(
    // RDKit✔️✔️:     std::vector<RDKit::SubstanceGroup> &sgs, unsigned int targetId) {
    // RDKit✔️✔️:   return std::find_if(sgs.begin(), sgs.end(), [targetId](const auto &sg) {
    // RDKit✔️✔️:     unsigned int pval;
    // RDKit✔️✔️:     if (sg.getPropIfPresent(cxsmilesindex, pval)) {
    // RDKit✔️✔️:       if (pval == targetId) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   });
    // RDKit✔️✔️: }
    // RDKit source (verbatim; parent resolution in parse_sgroup_hierarchy):
    // RDKit✔️✔️:     auto psg = find_matching_sgroup(sgs, parentId);
    // RDKit✔️✔️:     if (psg == sgs.end()) {
    // RDKit✔️✔️:       validParent = false;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       psg->getPropIfPresent("index", parentId);
    // RDKit✔️✔️:     }
    let cx_parent_id = u32::try_from(cx_parent_id).map_err(|_| {
        CxQueryLoweringError::InvalidGraph(
            "CX SGroup hierarchy parent id exceeds the source unsigned range".to_owned(),
        )
    })?;
    let Some(parent_index) = find_query_sgroup_by_cx_sequence_id(groups, cx_parent_id)? else {
        return Ok(None);
    };
    let parent = &groups[parent_index];
    let parent_property_id = match parent.props().get(b"index".as_slice()) {
        Some(value) => cosmolkit_core::property_value_to_uint(value)?,
        None => cx_parent_id,
    };
    Ok(Some((parent.id(), parent_property_id)))
}

pub(crate) fn apply_cx_sgroup_hierarchy_child(
    groups: &mut [SubstanceGroup],
    resolved_parent: Option<(SubstanceGroupId, u32)>,
    cx_child_id: usize,
) -> Result<bool, CxQueryLoweringError> {
    // RDKit source (verbatim; child loop in CXSmilesOps.cpp::parse_sgroup_hierarchy):
    // RDKit✔️✔️: std::vector<RDKit::SubstanceGroup>::iterator find_matching_sgroup(
    // RDKit✔️✔️:     std::vector<RDKit::SubstanceGroup> &sgs, unsigned int targetId) {
    // RDKit✔️✔️:   return std::find_if(sgs.begin(), sgs.end(), [targetId](const auto &sg) {
    // RDKit✔️✔️:     unsigned int pval;
    // RDKit✔️✔️:     if (sg.getPropIfPresent(cxsmilesindex, pval)) {
    // RDKit✔️✔️:       if (pval == targetId) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   });
    // RDKit✔️✔️: }
    // RDKit✔️✔️:           for (auto childId : children) {
    // RDKit✔️✔️:             if (childId >= sgs.size()) {
    // RDKit✔️✔️:               throw SmilesParseException(
    // RDKit✔️✔️:                   "child id references non-existent SGroup");
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             auto csg = find_matching_sgroup(sgs, childId);
    // RDKit✔️✔️:             if (csg != sgs.end()) {
    // RDKit✔️✔️:               unsigned int cid;
    // RDKit✔️✔️:               csg->getProp("index", cid);
    // RDKit✔️✔️:               csg->setProp("PARENT", parentId);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    let Some((parent_id, parent_property_id)) = resolved_parent else {
        return Ok(false);
    };
    if cx_child_id >= groups.len() {
        return Err(CxQueryLoweringError::InvalidGraph(
            "child id references non-existent SGroup".to_owned(),
        ));
    }
    let cx_child_id = u32::try_from(cx_child_id).map_err(|_| {
        CxQueryLoweringError::InvalidGraph(
            "CX SGroup hierarchy child id exceeds the source unsigned range".to_owned(),
        )
    })?;
    let Some(child_index) = find_query_sgroup_by_cx_sequence_id(groups, cx_child_id)? else {
        return Ok(false);
    };
    let child = &groups[child_index];
    let child_index_value = child.props().get(b"index".as_slice()).ok_or_else(|| {
        CxQueryLoweringError::InvalidGraph(
            "query SGroup child is missing its source index property".to_owned(),
        )
    })?;
    cosmolkit_core::property_value_to_uint(child_index_value)?;
    let child = &mut groups[child_index];
    // Native parentId is unsigned; preserve its UInt tag and propagate the
    // sole store error before installing the detached semantic projection.
    child.set_prop("PARENT", parent_property_id)?;
    child.set_parent(parent_id);
    Ok(true)
}

fn find_query_sgroup_by_cx_sequence_id(
    groups: &[SubstanceGroup],
    target_id: u32,
) -> Result<Option<usize>, CxQueryLoweringError> {
    // RDKit source (verbatim; find_matching_sgroup from CXSmilesOps.cpp):
    // RDKit✔️✔️: std::vector<RDKit::SubstanceGroup>::iterator find_matching_sgroup(
    // RDKit✔️✔️:     std::vector<RDKit::SubstanceGroup> &sgs, unsigned int targetId) {
    // RDKit✔️✔️:   return std::find_if(sgs.begin(), sgs.end(), [targetId](const auto &sg) {
    // RDKit✔️✔️:     unsigned int pval;
    // RDKit✔️✔️:     if (sg.getPropIfPresent(cxsmilesindex, pval)) {
    // RDKit✔️✔️:       if (pval == targetId) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   });
    // RDKit✔️✔️: }
    for (index, group) in groups.iter().enumerate() {
        let Some(value) = group.props().get(b"_cxsmilesindex".as_slice()) else {
            continue;
        };
        let sequence_id = cosmolkit_core::property_value_to_uint(value)?;
        if sequence_id == target_id {
            return Ok(Some(index));
        }
    }
    Ok(None)
}

pub(crate) fn apply_cx_sgroup_hierarchy_to_query(
    graph: &mut QueryGraph,
    hierarchies: &[CxSGroupHierarchy],
) -> Result<(), CxQueryLoweringError> {
    // RDKit✔️✔️: parent relationships and child references are visited in
    // source order; each successful child sets typed parent state and PARENT.
    let mut groups = query_substance_groups(graph).to_vec();
    let mut dirty = false;
    for hierarchy in hierarchies {
        let resolved_parent = resolve_cx_sgroup_hierarchy_parent(&groups, hierarchy.parent)?;
        for &child_id in &hierarchy.children {
            match apply_cx_sgroup_hierarchy_child(&mut groups, resolved_parent, child_id) {
                Ok(changed) => dirty |= changed,
                Err(error) => {
                    if dirty {
                        replace_query_substance_groups(graph, groups).map_err(|validation| {
                            CxQueryLoweringError::InvalidGraph(validation.to_string())
                        })?;
                    }
                    return Err(error);
                }
            }
        }
    }
    if dirty {
        replace_query_substance_groups(graph, groups)
            .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
    }
    Ok(())
}

fn replace_query_atom_like_rdkit(
    graph: &mut QueryGraph,
    atom_index: usize,
    incoming: &QueryAtom,
    preserve_props: bool,
) -> Result<(), CxQueryLoweringError> {
    // RDKit✔️❌: void RWMol::replaceAtom(unsigned int idx, Atom *atom_pin, bool,
    // RDKit✔️❌:                         bool preserveProps) {
    // RDKit✔️❌:   PRECONDITION(atom_pin, "bad atom passed to replaceAtom");
    // RDKit✔️❌:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit✔️❌:   auto atom_p = atom_pin->copy();
    // RDKit✔️❌:   atom_p->setOwningMol(this);
    // RDKit✔️❌:   atom_p->setIdx(idx);
    // RDKit✔️❌:   auto vd = boost::vertex(idx, d_graph);
    // RDKit✔️❌:   if (preserveProps) {
    // RDKit✔️❌:     const bool replaceExistingData = false;
    // RDKit✔️❌:     atom_p->updateProps(*d_graph[vd], replaceExistingData);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   const auto orig_p = d_graph[vd];
    // RDKit✔️❌:   delete orig_p;
    // RDKit✔️❌:   d_graph[vd] = atom_p;
    // RDKit✔️❌:
    // RDKit✔️❌:   // handle bookmarks
    // RDKit✔️❌:   for (auto &ab : d_atomBookmarks) {
    // RDKit✔️❌:     for (auto &elem : ab.second) {
    // RDKit✔️❌:       if (elem == orig_p) {
    // RDKit✔️❌:         elem = atom_p;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // handle stereo group
    // RDKit✔️❌:   for (auto &group : d_stereo_groups) {
    // RDKit✔️❌:     auto groupId = group.getReadId();
    // RDKit✔️❌:     auto atoms = group.getAtoms();
    // RDKit✔️❌:     auto bonds = group.getBonds();
    // RDKit✔️❌:     auto aiter = std::find(atoms.begin(), atoms.end(), orig_p);
    // RDKit✔️❌:     while (aiter != atoms.end()) {
    // RDKit✔️❌:       *aiter = atom_p;
    // RDKit✔️❌:       ++aiter;
    // RDKit✔️❌:       aiter = std::find(aiter, atoms.end(), orig_p);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     group = StereoGroup(group.getGroupType(), std::move(atoms),
    // RDKit✔️❌:                         std::move(bonds), groupId);
    // RDKit✔️❌:   }
    // RDKit✔️❌: };
    // Source copy helper: Code/GraphMol/QueryAtom.cpp
    // RDKit✔️❌: Atom *QueryAtom::copy() const {
    // RDKit✔️❌:   auto *res = new QueryAtom(*this);
    // RDKit✔️❌:   return static_cast<Atom *>(res);
    // RDKit✔️❌: }
    // Source copy helper: Code/GraphMol/QueryAtom.h
    // RDKit✔️❌:   QueryAtom(const QueryAtom &other) : Atom(other) {
    // RDKit✔️❌:     if (other.dp_query) {
    // RDKit✔️❌:       dp_query = other.dp_query->copy();
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       dp_query = nullptr;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // Source copy helper: Code/GraphMol/Atom.cpp
    // RDKit✔️❌:
    // RDKit✔️❌: void Atom::initFromOther(const Atom &other) {
    // RDKit✔️❌:   RDProps::operator=(other);
    // RDKit✔️❌:   // NOTE: we do *not* copy ownership!
    // RDKit✔️❌:   dp_mol = nullptr;
    // RDKit✔️❌:   d_atomicNum = other.d_atomicNum;
    // RDKit✔️❌:   d_index = 0;
    // RDKit✔️❌:   d_formalCharge = other.d_formalCharge;
    // RDKit✔️❌:   df_noImplicit = other.df_noImplicit;
    // RDKit✔️❌:   df_isAromatic = other.df_isAromatic;
    // RDKit✔️❌:   d_numExplicitHs = other.d_numExplicitHs;
    // RDKit✔️❌:   d_numRadicalElectrons = other.d_numRadicalElectrons;
    // RDKit✔️❌:   d_isotope = other.d_isotope;
    // RDKit✔️❌:   // d_pos = other.d_pos;
    // RDKit✔️❌:   d_chiralTag = other.d_chiralTag;
    // RDKit✔️❌:   d_hybrid = other.d_hybrid;
    // RDKit✔️❌:   d_implicitValence = other.d_implicitValence;
    // RDKit✔️❌:   d_explicitValence = other.d_explicitValence;
    // RDKit✔️❌:   if (other.dp_monomerInfo) {
    // RDKit✔️❌:     dp_monomerInfo = other.dp_monomerInfo->copy();
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     dp_monomerInfo = nullptr;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   d_flags = other.d_flags;
    // RDKit✔️❌: }
    // Source copy helper: Code/Query/Query.h
    // RDKit✔️❌:   virtual Query<MatchFuncArgType, DataFuncArgType, needsConversion> *copy()
    // RDKit✔️❌:       const {
    // RDKit✔️❌:     Query<MatchFuncArgType, DataFuncArgType, needsConversion> *res =
    // RDKit✔️❌:         new Query<MatchFuncArgType, DataFuncArgType, needsConversion>();
    // RDKit✔️❌:     for (auto iter = this->beginChildren(); iter != this->endChildren();
    // RDKit✔️❌:          ++iter) {
    // RDKit✔️❌:       res->addChild(CHILD_TYPE(iter->get()->copy()));
    // RDKit✔️❌:     }
    // RDKit✔️❌:     res->d_val = this->d_val;
    // RDKit✔️❌:     res->d_tol = this->d_tol;
    // RDKit✔️❌:     res->df_negate = this->df_negate;
    // RDKit✔️❌:     res->d_matchFunc = this->d_matchFunc;
    // RDKit✔️❌:     res->d_dataFunc = this->d_dataFunc;
    // RDKit✔️❌:     res->d_description = this->d_description;
    // RDKit✔️❌:     res->d_queryType = this->d_queryType;
    // RDKit✔️❌:     return res;
    // RDKit✔️❌:   }
    // Source copy helper: Code/Query/AndQuery.h
    // RDKit✔️❌:   Query<MatchFuncArgType, DataFuncArgType, needsConversion> *copy()
    // RDKit✔️❌:       const override {
    // RDKit✔️❌:     AndQuery<MatchFuncArgType, DataFuncArgType, needsConversion> *res =
    // RDKit✔️❌:         new AndQuery<MatchFuncArgType, DataFuncArgType, needsConversion>();
    // RDKit✔️❌:     typename BASE::CHILD_VECT_CI i;
    // RDKit✔️❌:     for (i = this->beginChildren(); i != this->endChildren(); ++i) {
    // RDKit✔️❌:       res->addChild(typename BASE::CHILD_TYPE(i->get()->copy()));
    // RDKit✔️❌:     }
    // RDKit✔️❌:     res->setNegation(this->getNegation());
    // RDKit✔️❌:     res->d_description = this->d_description;
    // RDKit✔️❌:     res->d_queryType = this->d_queryType;
    // RDKit✔️❌:     return res;
    // RDKit✔️❌:   }
    // Source copy helper: Code/Query/OrQuery.h
    // RDKit✔️❌:   Query<MatchFuncArgType, DataFuncArgType, needsConversion> *copy()
    // RDKit✔️❌:       const override {
    // RDKit✔️❌:     OrQuery<MatchFuncArgType, DataFuncArgType, needsConversion> *res =
    // RDKit✔️❌:         new OrQuery<MatchFuncArgType, DataFuncArgType, needsConversion>();
    // RDKit✔️❌:
    // RDKit✔️❌:     typename BASE::CHILD_VECT_CI i;
    // RDKit✔️❌:     for (i = this->beginChildren(); i != this->endChildren(); ++i) {
    // RDKit✔️❌:       res->addChild(typename BASE::CHILD_TYPE(i->get()->copy()));
    // RDKit✔️❌:     }
    // RDKit✔️❌:     res->setNegation(this->getNegation());
    // RDKit✔️❌:     res->d_description = this->d_description;
    // RDKit✔️❌:     res->d_queryType = this->d_queryType;
    // RDKit✔️❌:     return res;
    // RDKit✔️❌:   }
    // Source copy helper: Code/Query/XOrQuery.h
    // RDKit✔️❌:   Query<MatchFuncArgType, DataFuncArgType, needsConversion> *copy()
    // RDKit✔️❌:       const override {
    // RDKit✔️❌:     XOrQuery<MatchFuncArgType, DataFuncArgType, needsConversion> *res =
    // RDKit✔️❌:         new XOrQuery<MatchFuncArgType, DataFuncArgType, needsConversion>();
    // RDKit✔️❌:
    // RDKit✔️❌:     typename BASE::CHILD_VECT_CI i;
    // RDKit✔️❌:     for (i = this->beginChildren(); i != this->endChildren(); ++i) {
    // RDKit✔️❌:       res->addChild(typename BASE::CHILD_TYPE(i->get()->copy()));
    // RDKit✔️❌:     }
    // RDKit✔️❌:     res->setNegation(this->getNegation());
    // RDKit✔️❌:     res->d_description = this->d_description;
    // RDKit✔️❌:     res->d_queryType = this->d_queryType;
    // RDKit✔️❌:     return res;
    // RDKit✔️❌:   }
    // Source copy helper: Code/Query/EqualityQuery.h
    // RDKit✔️❌:   Query<MatchFuncArgType, DataFuncArgType, needsConversion> *copy()
    // RDKit✔️❌:       const override {
    // RDKit✔️❌:     EqualityQuery<MatchFuncArgType, DataFuncArgType, needsConversion> *res =
    // RDKit✔️❌:         new EqualityQuery<MatchFuncArgType, DataFuncArgType, needsConversion>();
    // RDKit✔️❌:     res->setNegation(this->getNegation());
    // RDKit✔️❌:     res->setVal(this->d_val);
    // RDKit✔️❌:     res->setTol(this->d_tol);
    // RDKit✔️❌:     res->setDataFunc(this->d_dataFunc);
    // RDKit✔️❌:     res->d_description = this->d_description;
    // RDKit✔️❌:     res->d_queryType = this->d_queryType;
    // RDKit✔️❌:     return res;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // Behavior: Rust's non-null borrow supplies the source atom precondition.
    // Check the target row before copying. Stable row IDs carry all references
    // that the source repairs by replacing pointers (including duplicates).
    // QueryGraph has no pointer-owning bookmarks or owning-molecule handles.
    // Copy incoming atom members/query/cache facts, then replace only its
    // dictionary when requested. Both branches retain incoming monomer info,
    // temporary flags and valence storage. Install before rebuilding groups.
    // Reconstruct every group in its original slot: no merging, sorting,
    // validation, or new error path. The read ID survives; the write ID resets.
    // Complexity: O(incoming query/property payload + preserved properties +
    // stereo memberships). PropertyStore's tree uses per-key allocations versus
    // the source contiguous Dict vector, hence a conservative cost-loss marker.
    let source = graph
        .atom(atom_index)
        .ok_or(CxQueryLoweringError::AtomIndex { index: atom_index })?;
    let mut replacement = incoming.clone().with_id(AtomId::new(atom_index));
    if preserve_props {
        replacement.replace_source_properties_from(source);
    }
    graph.atoms_mut()[atom_index] = replacement;
    for group in graph.stereo_groups_mut() {
        let mut replacement =
            StereoGroup::new(group.kind(), group.atoms().to_vec(), group.bonds().to_vec());
        if let Some(read_id) = group.id() {
            replacement = replacement.with_id(read_id);
        }
        *group = replacement;
    }
    Ok(())
}

fn add_query_like_rdkit(
    graph: &mut QueryGraph,
    atom_index: usize,
    predicate: QueryNode<AtomQueryPredicate>,
    label: Option<&cosmolkit_model::PropertyText>,
) -> Result<(), CxQueryLoweringError> {
    // Source: Code/GraphMol/SmilesParse/CXSmilesOps.cpp
    // RDKit✔️❌: void addquery(Q *qry, std::string symbol, RDKit::RWMol &mol, unsigned int idx) {
    // RDKit✔️❌:   PRECONDITION(qry, "bad query");
    // RDKit✔️❌:   auto *qa = new QueryAtom(0);
    // RDKit✔️❌:   qa->setQuery(qry);
    // RDKit✔️❌:   qa->setNoImplicit(true);
    // RDKit✔️❌:   bool updateLabel = false;
    // RDKit✔️❌:   bool preserveProps = true;
    // RDKit✔️❌:   mol.replaceAtom(idx, qa, updateLabel, preserveProps);
    // RDKit✔️❌:   if (symbol != "") {
    // RDKit✔️❌:     mol.getAtomWithIdx(idx)->setProp(RDKit::common_properties::atomLabel,
    // RDKit✔️❌:                                      symbol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   delete qa;
    // RDKit✔️❌: }
    // Source: Code/GraphMol/QueryAtom.h
    // RDKit✔️❌:   explicit QueryAtom(int num) : Atom(num), dp_query(makeAtomNumQuery(num)) {}
    // Source: Code/GraphMol/QueryAtom.h
    // RDKit✔️❌:   //! returns the label associated to this query
    // RDKit✔️❌:   std::string getQueryType() const override { return dp_query->getTypeLabel(); }
    // RDKit✔️❌:
    // RDKit✔️❌:   //! replaces our current query with the value passed in
    // RDKit✔️❌:   void setQuery(QUERYATOM_QUERY *what) override {
    // Source: Code/GraphMol/Atom.cpp
    // RDKit✔️❌: Atom::Atom(unsigned int num) : RDProps() {
    // RDKit✔️❌:   d_atomicNum = num;
    // RDKit✔️❌:   initAtom();
    // RDKit✔️❌: };
    // Source: Code/GraphMol/Atom.cpp
    // RDKit✔️❌: void Atom::initAtom() {
    // RDKit✔️❌:   df_isAromatic = false;
    // RDKit✔️❌:   df_noImplicit = false;
    // RDKit✔️❌:   d_numExplicitHs = 0;
    // RDKit✔️❌:   d_numRadicalElectrons = 0;
    // RDKit✔️❌:   d_formalCharge = 0;
    // RDKit✔️❌:   d_index = 0;
    // RDKit✔️❌:   d_isotope = 0;
    // RDKit✔️❌:   d_chiralTag = CHI_UNSPECIFIED;
    // RDKit✔️❌:   d_hybrid = UNSPECIFIED;
    // RDKit✔️❌:   dp_mol = nullptr;
    // RDKit✔️❌:   dp_monomerInfo = nullptr;
    // RDKit✔️❌:
    // RDKit✔️❌:   d_implicitValence = -1;
    // RDKit✔️❌:   d_explicitValence = -1;
    // RDKit✔️❌: }
    // Behavior: a typed owned predicate satisfies the non-null query
    // precondition. A fresh atomic-number-zero QueryAtom carries source Atom
    // defaults (including -1 valence facts and zero flags); setting the supplied
    // query replaces the constructor's temporary atomic-number predicate.
    // Delegate replacement with preserveProps=true, then overwrite atomLabel
    // only for a nonempty supplied symbol. RAII destroys the temporary query.
    // Complexity: one temporary query value and the source-required replacement
    // copy. The constructor's immediately discarded predicate need not allocate
    // in Rust. Overall cost includes PropertyStore's tree allocation difference
    // from the source contiguous dictionary, marked separately as a cost loss.
    let mut incoming = QueryAtom::from_identity_parts(
        AtomId::new(0),
        QueryAtomIdentity::from_atomic_number(0),
        predicate,
    );
    incoming.set_no_implicit(true);
    replace_query_atom_like_rdkit(graph, atom_index, &incoming, true)?;
    if let Some(label) = label.filter(|label| !label.is_empty()) {
        graph.atoms_mut()[atom_index].set_prop("atomLabel", label)?;
    }
    Ok(())
}

pub(crate) fn apply_cx_coordinate_bond(
    graph: &mut QueryGraph,
    reference: cosmolkit_cx::CxBondReference,
    kind: CxCoordinateBondKind,
) -> Result<(), CxQueryLoweringError> {
    // BEGIN COMPLETE PINNED SF189
    // RDKit✔️❌: bool parse_coordinate_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️❌:                             Bond::BondType typ, unsigned int startAtomIdx,
    // RDKit✔️❌:                             unsigned int startBondIdx) {
    // RDKit✔️❌:   if (first >= last || (*first != 'C' && *first != 'H')) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   if (first >= last || *first != ':') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️❌:     unsigned int aidx;
    // RDKit✔️❌:     unsigned int bidx;
    // RDKit✔️❌:     if (read_int_pair(first, last, aidx, bidx)) {
    // RDKit✔️❌:       if (VALID_ATIDX(aidx) && VALID_BNDIDX(bidx)) {
    // RDKit✔️❌:         auto bnd = get_bond_with_smiles_idx(mol, bidx - startBondIdx);
    // RDKit✔️❌:         if (!bnd || (bnd->getBeginAtomIdx() != aidx - startAtomIdx &&
    // RDKit✔️❌:                      bnd->getEndAtomIdx() != aidx - startAtomIdx)) {
    // RDKit✔️❌:           BOOST_LOG(rdWarningLog) << "BOND NOT FOUND! " << bidx
    // RDKit✔️❌:                                   << " involving atom " << aidx << std::endl;
    // RDKit✔️❌:           return false;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         bnd->setBondType(typ);
    // RDKit✔️❌:         if (bnd->getBeginAtomIdx() != aidx - startAtomIdx) {
    // RDKit✔️❌:           unsigned int tmp = bnd->getBeginAtomIdx();
    // RDKit✔️❌:           bnd->setBeginAtomIdx(aidx - startAtomIdx);
    // RDKit✔️❌:           bnd->setEndAtomIdx(tmp);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first < last && *first == ',') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF189
    // BEGIN COMPLETE PINNED SF184
    // RDKit✔️✔️: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
    // RDKit✔️✔️:   for (auto bnd : mol.bonds()) {
    // RDKit✔️✔️:     unsigned int smilesIdx;
    // RDKit✔️✔️:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
    // RDKit✔️✔️:         smilesIdx == idx) {
    // RDKit✔️✔️:       return bnd;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // END COMPLETE PINNED SF184
    // BEGIN COMPLETE REACHED Code/GraphMol/Bond.h void setBondType(BondType bT)
    // RDKit✔️✔️:   void setBondType(BondType bT) { d_bondType = bT; }
    // END COMPLETE REACHED Code/GraphMol/Bond.h void setBondType(BondType bT)
    // BEGIN COMPLETE REACHED Code/GraphMol/Bond.h unsigned int getBeginAtomIdx() const
    // RDKit✔️✔️:   unsigned int getBeginAtomIdx() const { return d_beginAtomIdx; }
    // END COMPLETE REACHED Code/GraphMol/Bond.h unsigned int getBeginAtomIdx() const
    // BEGIN COMPLETE REACHED Code/GraphMol/Bond.h unsigned int getEndAtomIdx() const
    // RDKit✔️✔️:   unsigned int getEndAtomIdx() const { return d_endAtomIdx; }
    // END COMPLETE REACHED Code/GraphMol/Bond.h unsigned int getEndAtomIdx() const
    // BEGIN COMPLETE REACHED Code/GraphMol/Bond.cpp void Bond::setBeginAtomIdx(unsigned int what)
    // RDKit✔️✔️: void Bond::setBeginAtomIdx(unsigned int what) {
    // RDKit✔️✔️:   if (dp_mol) {
    // RDKit✔️✔️:     URANGE_CHECK(what, getOwningMol().getNumAtoms());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   d_beginAtomIdx = what;
    // RDKit✔️✔️: }
    // END COMPLETE REACHED Code/GraphMol/Bond.cpp void Bond::setBeginAtomIdx(unsigned int what)
    // BEGIN COMPLETE REACHED Code/GraphMol/Bond.cpp void Bond::setEndAtomIdx(unsigned int what)
    // RDKit✔️✔️: void Bond::setEndAtomIdx(unsigned int what) {
    // RDKit✔️✔️:   if (dp_mol) {
    // RDKit✔️✔️:     URANGE_CHECK(what, getOwningMol().getNumAtoms());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   d_endAtomIdx = what;
    // RDKit✔️✔️: }
    // END COMPLETE REACHED Code/GraphMol/Bond.cpp void Bond::setEndAtomIdx(unsigned int what)
    // Behavior: graph-local startAtomIdx/startBondIdx are both zero. Native
    // u32 pair values were validated by the syntax owner, with no row guessing.
    // Compare both ranges before reading any source-index property. Missing
    // properties are skipped, wrong present property conversion propagates,
    // and first matching source index short circuits through canonical SF184.
    // Missing/mismatched bonds log the source payload before failure; only a
    // matching bond is mutated, type first and then begin/end reversal. Query
    // predicates and adjacency remain the same undirected edge's source state.
    // Cost: one O(E) borrowed source-index scan, no graph clone or reparse.
    if reference.atom >= graph.num_atoms() || reference.bond >= graph.num_bonds() {
        return Ok(());
    }
    let row = match query_bond_row_from_source_index(graph, reference.bond) {
        Ok(row) => row,
        Err(error @ CxQueryLoweringError::BondIndex { .. }) => {
            warn_cx_coordinate_bond_not_found(reference);
            return Err(error);
        }
        Err(error) => return Err(error),
    };
    let bond = graph
        .bonds_mut()
        .get_mut(row)
        .ok_or(CxQueryLoweringError::BondIndex {
            index: reference.bond,
        })?;
    let begin = bond.begin();
    if begin.index() != reference.atom && bond.end().index() != reference.atom {
        warn_cx_coordinate_bond_not_found(reference);
        return Err(CxQueryLoweringError::BondAtomMismatch {
            atom: reference.atom,
            bond: reference.bond,
        });
    }
    bond.bond_mut().set_order(match kind {
        CxCoordinateBondKind::Dative => BondOrder::Dative,
        CxCoordinateBondKind::Hydrogen => BondOrder::Hydrogen,
    });
    if begin.index() != reference.atom {
        bond.bond_mut()
            .set_endpoints(AtomId::new(reference.atom), begin);
    }
    Ok(())
}

fn warn_cx_coordinate_bond_not_found(reference: cosmolkit_cx::CxBondReference) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog) << "BOND NOT FOUND! " << bidx
    // RDKit✔️✔️:                                   << " involving atom " << aidx << std::endl;
    // Reproduce this source payload/default output plus newline/flush, with
    // non-throwing I/O like the source ostream default. Independent RDLog
    // enable flags, alternate sinks and prefix formatting are unmodeled.
    use std::io::Write;
    let mut output = std::io::stderr().lock();
    let _ = writeln!(
        output,
        "BOND NOT FOUND! {} involving atom {}",
        reference.bond, reference.atom
    )
    .and_then(|_| output.flush());
}

pub(crate) fn apply_cx_zero_bond(
    graph: &mut QueryGraph,
    index: usize,
) -> Result<(), CxQueryLoweringError> {
    // BEGIN COMPLETE PINNED SF190
    // RDKit✔️❌: bool parse_zero_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️❌:                       unsigned int, unsigned int startBondIdx) {
    // RDKit✔️❌:   // these look like: C1CCCCC~CCCC1 |Z:5|
    // RDKit✔️❌:   if (first >= last || *first != 'Z') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   if (first >= last || *first != ':') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:
    // RDKit✔️❌:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️❌:     unsigned int bondIdx;
    // RDKit✔️❌:     if (!read_int(first, last, bondIdx)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (VALID_BNDIDX(bondIdx)) {
    // RDKit✔️❌:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit✔️❌:
    // RDKit✔️❌:       if (!bond) {
    // RDKit✔️❌:         BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:             << "bond " << bondIdx
    // RDKit✔️❌:             << " not found, cannot mark as zero order bond." << std::endl;
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       bond->setBondType(Bond::ZERO);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first < last && *first == ',') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF190
    // BEGIN COMPLETE PINNED SF184
    // RDKit✔️✔️: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
    // RDKit✔️✔️:   for (auto bnd : mol.bonds()) {
    // RDKit✔️✔️:     unsigned int smilesIdx;
    // RDKit✔️✔️:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
    // RDKit✔️✔️:         smilesIdx == idx) {
    // RDKit✔️✔️:       return bnd;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // END COMPLETE PINNED SF184
    // RDKit✔️✔️: void setBondType(BondType bT) { d_bondType = bT; }
    // Source-local startBondIdx=0: skip an index outside the final bond count
    // before consulting properties. Use the first source index match; present
    // wrong-kind numeric properties propagate, with no warning or fallback.
    // A source-index hole logs the native payload before failure. Valid writes
    // change only d_bondType; query predicates, endpoints and properties stay.
    // One borrowed O(E) canonical scan and one assignment, without graph copy.
    if index >= graph.num_bonds() {
        return Ok(());
    }
    let row = match query_bond_row_from_source_index(graph, index) {
        Ok(row) => row,
        Err(error @ CxQueryLoweringError::BondIndex { .. }) => {
            warn_cx_zero_bond_not_found(index);
            return Err(error);
        }
        Err(error) => return Err(error),
    };
    let bond = graph
        .bonds_mut()
        .get_mut(row)
        .ok_or(CxQueryLoweringError::BondIndex { index })?;
    bond.bond_mut().set_order(BondOrder::Zero);
    Ok(())
}

fn warn_cx_zero_bond_not_found(index: usize) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "bond " << bondIdx
    // RDKit✔️✔️:             << " not found, cannot mark as zero order bond." << std::endl;
    // Source payload/default stderr, decimal native index, newline and flush;
    // default ostream I/O does not throw. Independent RDLog enable flags,
    // alternate sinks and prefixes remain unmodeled, not inferred here.
    use std::io::Write;
    let mut output = std::io::stderr().lock();
    let _ = writeln!(
        output,
        "bond {} not found, cannot mark as zero order bond.",
        index
    )
    .and_then(|_| output.flush());
}

/// Apply parsed CX records directly to the canonical query value.
///
/// Parsing remains owned by `cosmolkit-cx`; this function owns only the
/// destination semantics. It never projects query data through a concrete
/// molecule.
pub fn apply_cx_to_query_graph(
    graph: &mut QueryGraph,
    parsed: &ParsedCxExtensions,
) -> Result<(), CxQueryLoweringError> {
    let mut stereo_tracker = CxStereoGroupTracker::new(graph);
    let mut cx_sequence_id = 0_u32;
    for record in parsed.records() {
        match record {
            CxRecord::Coordinates(coordinates) => {
                // BEGIN COMPLETE PINNED SF188
                // RDKit✔️❌: bool parse_coords(Iterator &first, Iterator last, RDKit::RWMol &mol,
                // RDKit✔️❌:                   unsigned int startAtomIdx, unsigned int confIdx) {
                // RDKit✔️❌:   if (first >= last || *first != '(') {
                // RDKit✔️❌:     return false;
                // RDKit✔️❌:   }
                // RDKit✔️❌:
                // RDKit✔️❌:   auto *conf = new Conformer(mol.getNumAtoms());
                // RDKit✔️❌:   mol.addConformer(conf);
                // RDKit✔️❌:   conf->setId(confIdx);
                // RDKit✔️❌:   ++first;
                // RDKit✔️❌:   unsigned int atIdx = 0;
                // RDKit✔️❌:   bool is3D = false;
                // RDKit✔️❌:   while (first <= last && *first != ')') {
                // RDKit✔️❌:     RDGeom::Point3D pt;
                // RDKit✔️❌:     std::string tkn = read_text_to(first, last, ";)");
                // RDKit✔️❌:     if (VALID_ATIDX(atIdx)) {
                // RDKit✔️❌:       if (!tkn.empty()) {
                // RDKit✔️❌:         std::vector<std::string> tokens;
                // RDKit✔️❌:         boost::split(tokens, tkn, boost::is_any_of(std::string(",")));
                // RDKit✔️❌:         if (tokens.size() >= 1 && tokens[0].size()) {
                // RDKit✔️❌:           pt.x = boost::lexical_cast<double>(tokens[0]);
                // RDKit✔️❌:         }
                // RDKit✔️❌:         if (tokens.size() >= 2 && tokens[1].size()) {
                // RDKit✔️❌:           pt.y = boost::lexical_cast<double>(tokens[1]);
                // RDKit✔️❌:         }
                // RDKit✔️❌:         if (tokens.size() >= 3 && tokens[2].size()) {
                // RDKit✔️❌:           pt.z = boost::lexical_cast<double>(tokens[2]);
                // RDKit✔️❌:           is3D = true;
                // RDKit✔️❌:         }
                // RDKit✔️❌:       }
                // RDKit✔️❌:
                // RDKit✔️❌:       conf->setAtomPos(atIdx - startAtomIdx, pt);
                // RDKit✔️❌:     }
                // RDKit✔️❌:     ++atIdx;
                // RDKit✔️❌:     if (first <= last && *first != ')') {
                // RDKit✔️❌:       ++first;
                // RDKit✔️❌:     }
                // RDKit✔️❌:   }
                // RDKit✔️❌:   // make sure that the conformer really is 3D!
                // RDKit✔️❌:   if (is3D && hasNonZeroZCoords(*conf)) {
                // RDKit✔️❌:     conf->set3D(true);
                // RDKit✔️❌:   } else {
                // RDKit✔️❌:     conf->set3D(false);
                // RDKit✔️❌:   }
                // RDKit✔️❌:   if (first >= last || *first != ')') {
                // RDKit✔️❌:     return false;
                // RDKit✔️❌:   }
                // RDKit✔️❌:   ++first;
                // RDKit✔️❌:   return true;
                // RDKit✔️❌: }
                // END COMPLETE PINNED SF188
                let atom_count = graph.num_atoms();
                let source_atom_count = u32::try_from(atom_count).map_err(|_| {
                    CxQueryLoweringError::InvalidGraph(
                        "CX atom count exceeds source unsigned32 domain".to_owned(),
                    )
                })?;
                let mut values = vec![[0.0; 3]; atom_count];
                for (slot, value) in coordinates.values.iter().enumerate() {
                    let index = slot as u32;
                    if index < source_atom_count {
                        values[index as usize] = value.unwrap_or([0.0; 3]);
                    }
                }
                // This complete-record lowerer also accepts syntax-only records;
                // derive the final flag from the actual mapped destination rows.
                let is_3d = coordinates.is_3d && values.iter().any(|point| point[2].abs() > 1e-3);
                graph
                    .add_conformer_3d(Conformer3D::new(coordinates.conformer, values, is_3d))
                    .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
            }
            CxRecord::AtomLabels(values) => {
                // BEGIN COMPLETE PINNED SF187 GRAPH WRITE
                // RDKit✔️✔️: bool parse_atom_labels(Iterator &first, Iterator last, RDKit::RWMol &mol,
                // RDKit✔️✔️:                        unsigned int startAtomIdx) {
                // RDKit✔️✔️:   if (first >= last || *first != '$') {
                // RDKit✔️✔️:     return false;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   ++first;
                // RDKit✔️✔️:   unsigned int atIdx = 0;
                // RDKit✔️✔️:   while (first <= last && *first != '$') {
                // RDKit✔️✔️:     std::string tkn = read_text_to(first, last, ";$");
                // RDKit✔️✔️:     if (!tkn.empty() && VALID_ATIDX(atIdx)) {
                // RDKit✔️✔️:       mol.getAtomWithIdx(atIdx - startAtomIdx)
                // RDKit✔️✔️:           ->setProp(RDKit::common_properties::atomLabel, tkn);
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:     ++atIdx;
                // RDKit✔️✔️:     if (first <= last && *first != '$') {
                // RDKit✔️✔️:       ++first;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   if (first >= last || *first != '$') {
                // RDKit✔️✔️:     return false;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   ++first;
                // RDKit✔️✔️:   return true;
                // RDKit✔️✔️: }
                // END COMPLETE PINNED SF187 GRAPH WRITE
                // BEGIN RDKIT CPP FUNCTION parse_atom_labels
                // RDKit✔️✔️: std::string tkn = read_text_to(first, last, ";$");
                // RDKit✔️✔️: if (!tkn.empty() && VALID_ATIDX(atIdx)) {
                // RDKit✔️✔️:   mol.getAtomWithIdx(atIdx - startAtomIdx)
                // RDKit✔️✔️:       ->setProp(RDKit::common_properties::atomLabel, tkn);
                // RDKit✔️✔️: }
                // RDKit✔️✔️: ++atIdx;
                // END RDKIT CPP FUNCTION parse_atom_labels
                // `cosmolkit-cx` has decoded numeric character entities before
                // creating these ordered slots; keep each nonempty source value
                // unchanged and skip slots outside the detached query graph.
                // RDKit✔️✔️: unsigned int atIdx = 0;
                // RDKit✔️✔️: ++atIdx;
                // Transport slot remains usize; source atom ordinal is u32.
                for (slot, value) in values.iter().enumerate() {
                    let index = (slot as u32) as usize;
                    if index >= graph.num_atoms() {
                        continue;
                    }
                    let Some(value) = value.as_ref().filter(|value| !value.is_empty()) else {
                        continue;
                    };
                    let atom = graph
                        .atom_mut(index)
                        .ok_or(CxQueryLoweringError::AtomIndex { index })?;
                    atom.set_prop("atomLabel", value)?;
                }
            }
            CxRecord::AtomValues(values) => {
                // BEGIN COMPLETE PINNED SF185 GRAPH WRITE
                // RDKit✔️✔️: bool parse_atom_values(Iterator &first, Iterator last, RDKit::RWMol &mol,
                // RDKit✔️✔️:                        unsigned int startAtomIdx) {
                // RDKit✔️✔️:   if (first >= last || *first != ':') {
                // RDKit✔️✔️:     return false;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   ++first;
                // RDKit✔️✔️:   unsigned int atIdx = 0;
                // RDKit✔️✔️:   while (first <= last && *first != '$') {
                // RDKit✔️✔️:     std::string tkn = read_text_to(first, last, ";$");
                // RDKit✔️✔️:     if (tkn != "" && VALID_ATIDX(atIdx)) {
                // RDKit✔️✔️:       mol.getAtomWithIdx(atIdx)->setProp(RDKit::common_properties::molFileValue,
                // RDKit✔️✔️:                                          tkn);
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:     ++atIdx;
                // RDKit✔️✔️:     if (first <= last && *first != '$') {
                // RDKit✔️✔️:       ++first;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   if (first >= last || *first != '$') {
                // RDKit✔️✔️:     return false;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   ++first;
                // RDKit✔️✔️:   return true;
                // RDKit✔️✔️: }
                // END COMPLETE PINNED SF185 GRAPH WRITE
                // BEGIN RDKIT CPP FUNCTION parse_atom_values
                // RDKit✔️✔️: std::string tkn = read_text_to(first, last, ";$");
                // RDKit✔️✔️: if (tkn != "" && VALID_ATIDX(atIdx)) {
                // RDKit✔️✔️:   mol.getAtomWithIdx(atIdx)->setProp(
                // RDKit✔️✔️:       RDKit::common_properties::molFileValue, tkn);
                // RDKit✔️✔️: }
                // RDKit✔️✔️: ++atIdx;
                // END RDKIT CPP FUNCTION parse_atom_values
                // Values use the same already-decoded ordered slots as labels.
                // RDKit✔️✔️: unsigned int atIdx = 0;
                // RDKit✔️✔️: ++atIdx;
                // Record-slot ordinals are transport positions; chemistry uses
                // the source 32-bit unsigned atom index, including wraparound.
                for (slot, value) in values.iter().enumerate() {
                    let index = (slot as u32) as usize;
                    if index >= graph.num_atoms() {
                        continue;
                    }
                    let Some(value) = value.as_ref().filter(|value| !value.is_empty()) else {
                        continue;
                    };
                    let atom = graph
                        .atom_mut(index)
                        .ok_or(CxQueryLoweringError::AtomIndex { index })?;
                    atom.set_prop("molFileValue", value)?;
                }
            }
            CxRecord::AtomProperties(properties) => {
                // BEGIN RDKIT CPP FUNCTION parse_atom_props
                // RDKit✔️✔️: std::string pname = read_text_to(first, last, ".");
                // RDKit✔️✔️: if (!pname.empty()) {
                // RDKit✔️✔️:   std::string pval = read_text_to(first, last, ":|,");
                // RDKit✔️✔️:   if (VALID_ATIDX(atIdx) && !pval.empty()) {
                // RDKit✔️✔️:     mol.getAtomWithIdx(atIdx - startAtomIdx)
                // RDKit✔️✔️:         ->setProp(pname, pval);
                // RDKit✔️✔️:   }
                // RDKit✔️✔️: }
                // END RDKIT CPP FUNCTION parse_atom_props
                // Iterate the typed items in source order: repeated names
                // overwrite in that order, and invalid source atom indices are
                // skipped as they are by VALID_ATIDX.
                for property in properties {
                    if property.atom >= graph.num_atoms()
                        || property.name.is_empty()
                        || property.value.is_empty()
                    {
                        continue;
                    }
                    let atom =
                        graph
                            .atom_mut(property.atom)
                            .ok_or(CxQueryLoweringError::AtomIndex {
                                index: property.atom,
                            })?;
                    atom.set_prop(property.name.clone(), property.value.clone())?;
                }
            }
            CxRecord::CoordinateBonds(annotation) => {
                for reference in &annotation.bonds {
                    apply_cx_coordinate_bond(graph, *reference, annotation.kind)?;
                }
            }
            CxRecord::ZeroBonds(indices) => {
                for &index in indices {
                    apply_cx_zero_bond(graph, index)?;
                }
            }
            CxRecord::Unsaturation(indices) => {
                for item_index in 0..indices.len() {
                    apply_cx_query_constraint_item(graph, record, item_index)?;
                }
            }
            CxRecord::RingBonds(constraints) => {
                for item_index in 0..constraints.len() {
                    apply_cx_query_constraint_item(graph, record, item_index)?;
                }
            }
            CxRecord::Substitution(constraints) => {
                for item_index in 0..constraints.len() {
                    apply_cx_query_constraint_item(graph, record, item_index)?;
                }
            }
            CxRecord::EnhancedStereo(stereo) => {
                merge_cx_enhanced_stereo(graph, &mut stereo_tracker, stereo)
                    .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
            }
            CxRecord::WedgedBonds(wedges) => {
                for wedge in wedges {
                    apply_cx_wedge_bond_to_query(graph, wedge)?;
                }
            }
            CxRecord::DoubleBondStereo(stereo) => {
                for &index in &stereo.bonds {
                    apply_cx_double_bond_stereo_to_query(graph, index, stereo.stereo)?;
                }
            }
            CxRecord::Radicals(radicals) => {
                for radical in radicals {
                    // RDKit❗✔️: if (VALID_ATIDX(atIdx)) {
                    // RDKit❗✔️:   mol.getAtomWithIdx(atIdx - startAtomIdx)
                    // RDKit❗✔️:       ->setNumRadicalElectrons(numRadicalElectrons);
                    // RDKit❗✔️: }
                    if let Some(atom) = graph.atom_mut(radical.atom) {
                        atom.set_radical_electrons(radical.electrons);
                    }
                }
            }
            CxRecord::LinkNodes(nodes) => apply_cx_link_nodes_to_query(graph, nodes)?,
            CxRecord::DataSGroup(data) => {
                apply_cx_data_sgroup_to_query(graph, data, cx_sequence_id)?;
                cx_sequence_id = cx_sequence_id.wrapping_add(1);
            }
            CxRecord::SGroupHierarchy(hierarchies) => {
                apply_cx_sgroup_hierarchy_to_query(graph, hierarchies)?;
            }
            CxRecord::PolymerSGroup(polymer) => {
                apply_cx_polymer_sgroup_to_query(graph, polymer, cx_sequence_id)?;
                cx_sequence_id = cx_sequence_id.wrapping_add(1);
            }
            CxRecord::VariableAttachments(attachments) => {
                for attachment in attachments {
                    apply_cx_variable_attachment_to_query(graph, attachment)?;
                }
            }
            CxRecord::Unknown(_) => {}
        }
    }
    finish_cx_smiles_labels(graph)?;
    graph
        .validate()
        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))
}

#[cfg(test)]
mod tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::*;
    use crate::query_behavior::{
        make_atom_in_ring_of_size_query, make_atom_min_ring_size_query,
        make_atom_ring_bond_count_query, make_atom_ring_query,
    };
    use cosmolkit_cx::{
        CxAtomConstraint, CxBondReference, CxCoordinateBondKind, CxCoordinateBonds,
        CxCountConstraint, CxLinkNode, CxRingBond, CxSGroupHierarchy, ParsedCxExtensions,
    };
    use cosmolkit_model::{
        Atom, AtomSpec, BondId, BondSpec, PropertyValue, QueryAtom, SGroupBondRole, SGroupCState,
        SGroupData, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, query_substance_groups,
        replace_query_substance_groups,
    };
    use cosmolkit_types::Element;

    fn graph() -> QueryGraph {
        let atoms = vec![
            cosmolkit_model::QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
            cosmolkit_model::QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![cosmolkit_model::QueryBond::new(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        QueryGraph::from_parts(
            atoms,
            bonds,
            Vec::<(
                cosmolkit_model::PropertyText,
                cosmolkit_model::PropertyValue,
            )>::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    #[test]
    fn query_sgroups_cx_labels_replace_special_query_atoms_and_preserve_typed_state() {
        let labels = [
            ("star_e", crate::query_behavior::make_atom_null_query()),
            ("Q_e", crate::query_behavior::make_q_atom_query()),
            ("QH_p", crate::query_behavior::make_q_h_atom_query()),
            ("AH_p", crate::query_behavior::make_a_h_atom_query()),
            ("X_p", crate::query_behavior::make_x_atom_query()),
            ("XH_p", crate::query_behavior::make_x_h_atom_query()),
            ("M_p", crate::query_behavior::make_m_atom_query()),
            ("MH_p", crate::query_behavior::make_m_h_atom_query()),
        ];

        for (label, expected_predicate) in labels {
            let mut query = graph();
            let atom = query.atom_mut(0).expect("first query atom");
            atom.set_formal_charge(1);
            atom.set_isotope(Some(13));
            atom.set_atom_map(Some(9));
            atom.set_prop("dummyLabel", "stale")
                .expect("dummy label property");
            atom.set_prop("sourceProperty", "retained")
                .expect("source property");

            let substance_group =
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_rdkit_sequence_id(23)
                    .with_external_id(41)
                    .with_atoms(vec![AtomId::new(0), AtomId::new(1), AtomId::new(0)])
                    .with_bonds(vec![BondId::new(0), BondId::new(0)])
                    .with_bond_role(BondId::new(0), SGroupBondRole::Contained)
                    .with_head_crossing_bonds(vec![BondId::new(0), BondId::new(0)])
                    .with_crossing_bond_correspondence(vec![BondId::new(0)])
                    .with_parent_atoms(vec![AtomId::new(1)])
                    .with_label("typed polymer label")
                    .with_data(SGroupData {
                        field_name: Some("FIELD".into()),
                        field_type: Some("S".into()),
                        values: vec!["one".into(), "two".into()],
                        ..SGroupData::default()
                    })
                    .with_cstates(vec![SGroupCState::new(BondId::new(0), [0.25, 0.5, 0.75])])
                    .with_prop("origin", "existing")
                    .unwrap()
                    .with_data_field("first source row")
                    .with_data_field("second source row");
            replace_query_substance_groups(&mut query, vec![substance_group.clone()])
                .expect("typed SGroup references are valid");

            let stereo_group = StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(0), AtomId::new(1), AtomId::new(0)],
                vec![BondId::new(0), BondId::new(0)],
            )
            .with_id(71);
            replace_query_stereo_groups(&mut query, vec![stereo_group.clone()])
                .expect("stereo group references are valid");

            let parsed = ParsedCxExtensions::new(
                vec![CxRecord::AtomLabels(vec![Some(label.into()), None])],
                0,
            );
            apply_cx_to_query_graph(&mut query, &parsed).expect("source CX label lowering");

            let atom = query.atom(0).expect("replacement query atom");
            assert_eq!(
                atom.identity(),
                QueryAtomIdentity::Element(Element::DUMMY),
                "{label}"
            );
            assert_eq!(atom.predicate(), &expected_predicate, "{label}");
            assert!(atom.no_implicit(), "{label}");
            assert_eq!(atom.formal_charge(), 0, "{label}");
            assert_eq!(atom.isotope(), None, "{label}");
            assert_eq!(atom.atom_map(), Some(9), "{label}");
            assert_eq!(
                atom.prop("atomLabel"),
                Some(&cosmolkit_model::PropertyValue::String(label.into())),
                "{label}"
            );
            assert_eq!(atom.prop("dummyLabel"), None, "{label}");
            assert_eq!(
                atom.prop("sourceProperty"),
                Some(&cosmolkit_model::PropertyValue::String("retained".into())),
                "{label}"
            );
            assert_eq!(query_substance_groups(&query), &[substance_group]);
            assert_eq!(query.stereo_groups(), &[stereo_group]);
            assert_eq!(query.prop(CX_LABELS_PROCESSED_PROP), None);
        }
    }

    #[test]
    fn query_sgroups_cx_labels_keep_ordinary_names_and_strip_only_pseudo_suffixes() {
        for (label, expected_label, expected_dummy) in [
            ("ordinary", Some("ordinary"), None),
            ("Pol", Some("Pol"), None),
            ("Pol_p", None, Some("Pol")),
            ("Mod_p", None, Some("Mod")),
        ] {
            let mut query = graph();
            let atom = query.atom_mut(0).expect("first query atom");
            atom.set_formal_charge(1);
            atom.set_prop("dummyLabel", "stale")
                .expect("dummy label property");
            let source_predicate = atom.predicate().clone();
            let source_identity = atom.identity();
            let parsed = ParsedCxExtensions::new(
                vec![CxRecord::AtomLabels(vec![Some(label.into()), None])],
                0,
            );

            apply_cx_to_query_graph(&mut query, &parsed).expect("source CX label lowering");

            let atom = query.atom(0).expect("unreplaced query atom");
            assert_eq!(atom.identity(), source_identity, "{label}");
            assert_eq!(atom.predicate(), &source_predicate, "{label}");
            assert_eq!(atom.formal_charge(), 1, "{label}");
            assert_eq!(
                atom.prop("atomLabel"),
                expected_label
                    .map(cosmolkit_model::PropertyValue::from)
                    .as_ref(),
                "{label}"
            );
            assert_eq!(
                atom.prop("dummyLabel"),
                expected_dummy
                    .map(cosmolkit_model::PropertyValue::from)
                    .as_ref(),
                "{label}"
            );
        }
    }

    #[test]
    fn query_sgroups_cx_labels_apply_atomprop_precedence_after_label_records() {
        let mut query = graph();
        let original_predicate = query.atom(0).expect("first atom").predicate().clone();
        query
            .atom_mut(0)
            .expect("first atom")
            .set_prop("dummyLabel", "stale")
            .expect("dummy label property");
        let parsed = ParsedCxExtensions::new(
            vec![
                CxRecord::AtomLabels(vec![Some("Q_e".into()), None]),
                CxRecord::AtomProperties(vec![cosmolkit_cx::CxAtomProperty {
                    atom: 0,
                    name: "atomLabel".into(),
                    value: "ordinary after atomProp".into(),
                }]),
            ],
            0,
        );

        apply_cx_to_query_graph(&mut query, &parsed).expect("source property ordering");

        let atom = query.atom(0).expect("first atom");
        assert_eq!(
            atom.prop("atomLabel"),
            Some(&cosmolkit_model::PropertyValue::String(
                "ordinary after atomProp".into()
            ))
        );
        assert_eq!(atom.prop("dummyLabel"), None);
        assert_eq!(atom.predicate(), &original_predicate);
        assert_eq!(atom.identity(), QueryAtomIdentity::Element(Element::C));
    }

    #[test]
    fn query_sgroups_cx_labels_guard_repeated_processing_until_finish() {
        let mut query = graph();
        query
            .atom_mut(0)
            .expect("first atom")
            .set_prop("atomLabel", "Q_e")
            .expect("special label");
        process_cx_smiles_labels(&mut query).expect("first source label pass");
        let first_predicate = query.atom(0).expect("replacement atom").predicate().clone();
        assert_eq!(first_predicate, crate::query_behavior::make_q_atom_query());

        let atom = query.atom_mut(0).expect("replacement atom");
        atom.set_prop("atomLabel", "Pol_p")
            .expect("later source label");
        atom.set_prop("dummyLabel", "later property")
            .expect("later dummy label");
        process_cx_smiles_labels(&mut query).expect("guarded repeated source pass");

        let atom = query.atom(0).expect("guarded atom");
        assert_eq!(
            atom.prop("atomLabel"),
            Some(&cosmolkit_model::PropertyValue::String("Pol_p".into()))
        );
        assert_eq!(
            atom.prop("dummyLabel"),
            Some(&cosmolkit_model::PropertyValue::String(
                "later property".into()
            ))
        );
        assert_eq!(atom.predicate(), &first_predicate);
        assert_eq!(
            query.prop(CX_LABELS_PROCESSED_PROP),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );

        finish_cx_smiles_labels(&mut query).expect("outer parse clears the guard");
        assert_eq!(query.prop(CX_LABELS_PROCESSED_PROP), None);
    }

    #[test]
    fn query_sgroups_cx_labels_run_before_dat_and_polymer_group_attachment() {
        let mut query = graph();
        query
            .atom_mut(0)
            .expect("first atom")
            .set_prop("atomLabel", "star_e")
            .expect("special label");
        process_cx_smiles_labels(&mut query).expect("source processes labels before attachment");

        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(5)
                .with_atoms(vec![AtomId::new(0)])
                .with_label("DAT after label pass")
                .with_data_field("typed data"),
            SubstanceGroup::new(
                SubstanceGroupId::new(1),
                SubstanceGroupKind::StructuralRepeatUnit,
            )
            .with_parent(SubstanceGroupId::new(0))
            .with_rdkit_sequence_id(9)
            .with_atoms(vec![AtomId::new(0), AtomId::new(1), AtomId::new(0)])
            .with_bonds(vec![BondId::new(0), BondId::new(0)])
            .with_head_crossing_bonds(vec![BondId::new(0), BondId::new(0)])
            .with_crossing_bond_correspondence(vec![BondId::new(0)])
            .with_label("polymer after label pass"),
        ];
        replace_query_substance_groups(&mut query, groups.clone())
            .expect("DAT and polymer references use unchanged graph IDs");

        finish_cx_smiles_labels(&mut query).expect("outer parse clears the label guard");

        assert_eq!(query_substance_groups(&query), groups);
        assert_eq!(
            query.atom(0).unwrap().predicate(),
            &crate::query_behavior::make_atom_null_query()
        );
    }

    #[test]
    fn query_sgroups_cx_labels_convert_unqueried_dummy_atoms_to_a_query() {
        let atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::DUMMY));
        let query_atom = QueryAtom::from_carrier_parts(
            atom,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(0)),
        );
        let mut query = QueryGraph::from_parts(
            vec![query_atom],
            Vec::new(),
            Vec::<(
                cosmolkit_model::PropertyText,
                cosmolkit_model::PropertyValue,
            )>::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("ordinary dummy carrier graph");
        let parsed = ParsedCxExtensions::new(Vec::new(), 0);

        apply_cx_to_query_graph(&mut query, &parsed).expect("source fallback label pass");

        let atom = query.atom(0).expect("fallback query atom");
        assert_eq!(atom.identity(), QueryAtomIdentity::Element(Element::DUMMY));
        assert_eq!(
            atom.predicate(),
            &crate::query_behavior::make_a_atom_query()
        );
        assert!(atom.no_implicit());
    }

    #[test]
    fn lowers_query_constraints_without_concrete_projection() {
        let parsed = ParsedCxExtensions::new(
            vec![
                CxRecord::Unsaturation(vec![0]),
                CxRecord::Substitution(vec![CxAtomConstraint {
                    atom: 1,
                    constraint: CxCountConstraint::Exact(1),
                }]),
                CxRecord::AtomLabels(vec![Some("left".into()), None]),
            ],
            8,
        );
        let mut query = graph();
        apply_cx_to_query_graph(&mut query, &parsed).unwrap();
        assert_eq!(
            query.atom(0).unwrap().prop("atomLabel"),
            Some(&PropertyValue::String("left".into()))
        );
        assert!(matches!(
            query.atom(0).unwrap().predicate(),
            QueryNode::And(children) if children.len() == 2
        ));
        assert!(matches!(
            query.atom(1).unwrap().predicate(),
            QueryNode::And(children) if children.len() == 2
        ));
    }

    #[test]
    fn cx_progress_bonds_record_lowering_orients_pairs_and_skips_invalid_indices() {
        let mut query = graph();
        // CXSmilesOps::get_bond_with_smiles_idx resolves parser properties,
        // never final row IDs. Source smarts.yy assigns this UInt on insertion.
        // Preserve every original orientation/order/invalid-index assertion.
        query.bonds_mut()[0]
            .bond_mut()
            .set_prop("_cxsmilesBondIdx", PropertyValue::UInt(0))
            .expect("source parser bond-index property");
        let coordinate = ParsedCxExtensions::new(
            vec![CxRecord::CoordinateBonds(CxCoordinateBonds {
                kind: CxCoordinateBondKind::Dative,
                bonds: vec![
                    CxBondReference { atom: 1, bond: 0 },
                    CxBondReference { atom: 9, bond: 9 },
                    CxBondReference { atom: 9, bond: 0 },
                    CxBondReference { atom: 0, bond: 9 },
                ],
            })],
            0,
        );
        apply_cx_to_query_graph(&mut query, &coordinate).unwrap();
        assert_eq!(query.bonds_mut()[0].bond().order(), BondOrder::Dative);
        assert_eq!(query.bonds_mut()[0].endpoints(), (1, 0));

        let zero = ParsedCxExtensions::new(vec![CxRecord::ZeroBonds(vec![0, 9])], 0);
        apply_cx_to_query_graph(&mut query, &zero).unwrap();
        assert_eq!(query.bonds_mut()[0].bond().order(), BondOrder::Zero);
        assert_eq!(query.bonds_mut()[0].endpoints(), (1, 0));
    }

    #[test]
    fn q07e_ring_factories_and_cx_lowering_preserve_i32_targets_and_sentinels() {
        let maximum = 2_147_483_639;
        for target in [0, 255, 256, maximum] {
            assert_eq!(
                make_atom_ring_query(target),
                QueryNode::predicate(AtomQueryPredicate::NumAtomRings(target))
            );
            assert_eq!(
                make_atom_in_ring_of_size_query(target),
                QueryNode::predicate(AtomQueryPredicate::InRingOfSize(target))
            );
            assert_eq!(
                make_atom_min_ring_size_query(target),
                QueryNode::predicate(AtomQueryPredicate::SmallestRingSize(target))
            );
            assert_eq!(
                make_atom_ring_bond_count_query(target),
                QueryNode::predicate(AtomQueryPredicate::RingBondCount(target))
            );
        }
        assert_eq!(
            make_atom_ring_query(-1),
            QueryNode::predicate(AtomQueryPredicate::NumAtomRings(-1))
        );
        assert_eq!(
            make_atom_ring_bond_count_query(i32::MIN),
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(i32::MIN))
        );

        fn contains_predicate(
            node: &QueryNode<AtomQueryPredicate>,
            expected: &AtomQueryPredicate,
        ) -> bool {
            match node {
                QueryNode::Predicate(predicate) => predicate == expected,
                QueryNode::And(children) | QueryNode::Or(children) => children
                    .iter()
                    .any(|child| contains_predicate(child, expected)),
                _ => false,
            }
        }

        let mut query = graph();
        let carrier_before = query.atom(0).unwrap().try_to_atom().unwrap();
        let parsed = ParsedCxExtensions::new(
            vec![CxRecord::RingBonds(vec![
                CxRingBond {
                    atom: 0,
                    constraint: CxCountConstraint::Exact(3),
                },
                CxRingBond {
                    atom: 0,
                    constraint: CxCountConstraint::QueryScan,
                },
                CxRingBond {
                    atom: 1,
                    constraint: CxCountConstraint::LessEqual(4),
                },
            ])],
            8,
        );
        apply_cx_to_query_graph(&mut query, &parsed).unwrap();

        let atom_zero = query.atom(0).unwrap();
        assert!(contains_predicate(
            atom_zero.predicate(),
            &AtomQueryPredicate::RingBondCount(3)
        ));
        assert!(contains_predicate(
            atom_zero.predicate(),
            &AtomQueryPredicate::RingBondCount(0xDEAD_BEEF_u32 as i32)
        ));
        assert_eq!(atom_zero.try_to_atom().unwrap(), carrier_before);
        assert!(!atom_zero.predicate_is_carrier_derived());
        assert!(contains_predicate(
            query.atom(1).unwrap().predicate(),
            &AtomQueryPredicate::RingBondCountLessEqual(4)
        ));
    }

    #[test]
    fn cx_progress_linknodes_empty_record_preserves_existing_property() {
        let parsed = ParsedCxExtensions::new(vec![CxRecord::LinkNodes(Vec::new())], 4);
        let mut query = graph().with_prop("molFileLinkNodes", "prior").unwrap();
        apply_cx_to_query_graph(&mut query, &parsed).expect("empty source accumulator");
        assert_eq!(
            query.prop("molFileLinkNodes").map(fixture_value),
            Some("prior")
        );
    }

    #[test]
    fn cx_progress_linknodes_lowering_preserves_order_and_neighbor_order() {
        let atoms = (0..3)
            .map(|index| {
                cosmolkit_model::QueryAtom::new(AtomId::new(index), AtomSpec::new(Element::C))
            })
            .collect();
        let bonds = vec![
            cosmolkit_model::QueryBond::new(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            cosmolkit_model::QueryBond::new(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ),
        ];
        let mut query = QueryGraph::from_parts(
            atoms,
            bonds,
            Vec::<(
                cosmolkit_model::PropertyText,
                cosmolkit_model::PropertyValue,
            )>::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("degree-two center graph");
        let parsed = ParsedCxExtensions::new(
            vec![CxRecord::LinkNodes(vec![
                CxLinkNode {
                    atom: 0,
                    start_repetitions: 1,
                    end_repetitions: 3,
                    outer_atoms: Some([1, 2]),
                },
                CxLinkNode {
                    atom: 0,
                    start_repetitions: 2,
                    end_repetitions: 4,
                    outer_atoms: None,
                },
                CxLinkNode {
                    atom: 9,
                    start_repetitions: 7,
                    end_repetitions: 8,
                    outer_atoms: None,
                },
            ])],
            0,
        );

        apply_cx_to_query_graph(&mut query, &parsed).expect("source link-node projection");

        assert_eq!(
            query.prop("_molLinkNodes").map(fixture_value),
            Some("1 3 2 1 2 1 3|2 4 2 1 2 1 3")
        );
    }

    #[test]
    fn cx_progress_linknodes_later_degree_error_keeps_prior_property() {
        let parsed = ParsedCxExtensions::new(
            vec![CxRecord::LinkNodes(vec![
                CxLinkNode {
                    atom: 0,
                    start_repetitions: 1,
                    end_repetitions: 2,
                    outer_atoms: Some([1, 0]),
                },
                CxLinkNode {
                    atom: 1,
                    start_repetitions: 3,
                    end_repetitions: 4,
                    outer_atoms: None,
                },
            ])],
            0,
        );
        let mut query = graph().with_prop("molFileLinkNodes", "prior").unwrap();

        assert!(apply_cx_to_query_graph(&mut query, &parsed).is_err());
        assert_eq!(
            query.prop("molFileLinkNodes").map(fixture_value),
            Some("prior")
        );
    }

    #[test]
    fn cx_progress_radicals_source_skips_out_of_range_atom_indices() {
        let parsed = ParsedCxExtensions::new(
            vec![CxRecord::Radicals(vec![
                cosmolkit_cx::CxRadical {
                    atom: 0,
                    electrons: 2,
                },
                cosmolkit_cx::CxRadical {
                    atom: 3,
                    electrons: 1,
                },
            ])],
            4,
        );
        let mut query = graph();
        apply_cx_to_query_graph(&mut query, &parsed).expect("source-skipped radical index");
        assert_eq!(query.atom(0).unwrap().radical_electrons(), 2);
        assert_eq!(query.atom(1).unwrap().radical_electrons(), 0);
    }

    #[test]
    fn cx_progress_stereo_merge_reconstruction_clears_previous_bonds() {
        let mut query = graph();
        query.add_stereo_group(
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(0)],
                vec![BondId::new(0)],
            )
            .with_id(4),
        );
        let mut tracker = CxStereoGroupTracker {
            hashes: vec![41],
            first_group_index: 0,
        };
        let incoming = CxEnhancedStereo {
            kind: CxStereoGroupKind::Or,
            group_id: 4,
            atoms: vec![1],
        };

        merge_cx_enhanced_stereo(&mut query, &mut tracker, &incoming)
            .expect("reconstruct tracked source group");

        assert_eq!(query.stereo_groups().len(), 1);
        assert_eq!(query.stereo_groups()[0].kind(), StereoGroupKind::Or);
        assert_eq!(query.stereo_groups()[0].id(), Some(4));
        assert_eq!(
            query.stereo_groups()[0].atoms(),
            &[AtomId::new(0), AtomId::new(1)]
        );
        assert!(query.stereo_groups()[0].bonds().is_empty());
    }

    #[test]
    fn cx_progress_hierarchy_maps_source_sequences_to_canonical_groups() {
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(5)
                .with_prop("_cxsmilesindex", "5")
                .unwrap()
                .with_prop("index", "71")
                .unwrap(),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(0)
                .with_prop("_cxsmilesindex", "0")
                .unwrap()
                .with_prop("index", "92")
                .unwrap(),
            SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(1)
                .with_prop("_cxsmilesindex", "1")
                .unwrap()
                .with_prop("index", "93")
                .unwrap(),
        ];
        let parsed =
            cosmolkit_cx::parse_cx_extensions("|SgH:5:0.1.0|").expect("source hierarchy syntax");
        let mut query = graph();
        replace_query_substance_groups(&mut query, groups).expect("initial query SGroups");

        apply_cx_to_query_graph(&mut query, &parsed).expect("source hierarchy lowering");

        let groups = query_substance_groups(&query);
        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0].rdkit_sequence_id(), Some(5));
        assert_eq!(groups[1].rdkit_sequence_id(), Some(0));
        assert_eq!(groups[2].rdkit_sequence_id(), Some(1));
        assert_eq!(groups[0].parent(), None);
        assert_eq!(groups[1].parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            groups[1]
                .props()
                .get("PARENT".as_bytes())
                .map(|value| value.as_uint().expect("source unsigned property kind")),
            Some(71_u32)
        );
        assert_eq!(groups[2].parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            groups[2]
                .props()
                .get("PARENT".as_bytes())
                .map(|value| value.as_uint().expect("source unsigned property kind")),
            Some(71_u32)
        );
    }

    #[test]
    fn cx_progress_hierarchy_skips_unmatched_parent_before_child_range_check() {
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(5)
                .with_prop("_cxsmilesindex", "5")
                .unwrap()
                .with_prop("index", "71")
                .unwrap(),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(0)
                .with_prop("_cxsmilesindex", "0")
                .unwrap()
                .with_prop("index", "92")
                .unwrap(),
        ];
        let parsed = cosmolkit_cx::parse_cx_extensions("|SgH:99:4294967295|")
            .expect("source hierarchy syntax");
        let mut query = graph();
        replace_query_substance_groups(&mut query, groups.clone()).expect("initial query SGroups");

        apply_cx_to_query_graph(&mut query, &parsed)
            .expect("source skips every child when its parent is missing");

        assert_eq!(query_substance_groups(&query), groups);
    }

    #[test]
    fn cx_progress_hierarchy_keeps_prior_child_when_later_child_fails() {
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(5)
                .with_prop("_cxsmilesindex", "5")
                .unwrap()
                .with_prop("index", "71")
                .unwrap(),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(0)
                .with_prop("_cxsmilesindex", "0")
                .unwrap()
                .with_prop("index", "92")
                .unwrap(),
            SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(1)
                .with_prop("_cxsmilesindex", "1")
                .unwrap(),
        ];
        let parsed =
            cosmolkit_cx::parse_cx_extensions("|SgH:5:0.1|").expect("source hierarchy syntax");
        let mut query = graph();
        replace_query_substance_groups(&mut query, groups).expect("initial query SGroups");

        let error = apply_cx_to_query_graph(&mut query, &parsed)
            .expect_err("matched child without source index property fails");

        assert!(error.to_string().contains("source index property"));
        let groups = query_substance_groups(&query);
        assert_eq!(groups[1].parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            groups[1]
                .props()
                .get("PARENT".as_bytes())
                .map(|value| value.as_uint().expect("source unsigned property kind")),
            Some(71_u32)
        );
        assert_eq!(groups[2].parent(), None);
        assert_eq!(groups[2].props().get("PARENT".as_bytes()), None);
    }

    #[test]
    fn cx_progress_hierarchy_without_parent_index_uses_cx_parent_id() {
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(17)
                .with_prop("_cxsmilesindex", "17")
                .unwrap(),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(0)
                .with_prop("_cxsmilesindex", "0")
                .unwrap()
                .with_prop("index", "4")
                .unwrap(),
        ];
        let parsed = ParsedCxExtensions::new(
            vec![CxRecord::SGroupHierarchy(vec![CxSGroupHierarchy {
                parent: 17,
                children: vec![0],
            }])],
            0,
        );
        let mut query = graph();
        replace_query_substance_groups(&mut query, groups).expect("initial query SGroups");

        apply_cx_to_query_graph(&mut query, &parsed).expect("optional parent index uses CX id");

        let child = &query_substance_groups(&query)[1];
        assert_eq!(child.parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            child
                .props()
                .get("PARENT".as_bytes())
                .map(|value| value.as_uint().expect("source unsigned property kind")),
            Some(17_u32)
        );
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    fn query_graph(props: Vec<cosmolkit_model::PropertyValue>) -> cosmolkit_model::QueryGraph {
        let atoms = (0..props.len() + 1)
            .map(|i| {
                cosmolkit_model::QueryAtom::new(
                    cosmolkit_model::AtomId::new(i),
                    cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::C),
                )
            })
            .collect();
        let bonds = props
            .into_iter()
            .enumerate()
            .map(|(i, v)| {
                cosmolkit_model::QueryBond::new(
                    cosmolkit_model::BondId::new(i),
                    cosmolkit_model::BondSpec::new(
                        cosmolkit_model::AtomId::new(i),
                        cosmolkit_model::AtomId::new(i + 1),
                        cosmolkit_types::BondOrder::Single,
                    )
                    .with_prop("_cxsmilesBondIdx", v)
                    .unwrap(),
                )
            })
            .collect();
        cosmolkit_model::QueryGraph::from_parts(
            atoms,
            bonds,
            Vec::<(
                cosmolkit_model::PropertyText,
                cosmolkit_model::PropertyValue,
            )>::new(),
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }

    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_search/CXlower_0
    #[test]
    fn uint_cell_unsigned_consumer_search_cxlower_0_cx_lowering() {
        let g = query_graph(vec![cosmolkit_model::PropertyValue::UInt(0_u32)]);
        let before = g.clone();
        assert_eq!(query_bond_row_from_source_index(&g, 0_usize), Ok(0));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_search/CXlower_1
    #[test]
    fn uint_cell_unsigned_consumer_search_cxlower_1_cx_lowering() {
        let g = query_graph(vec![cosmolkit_model::PropertyValue::UInt(1_u32)]);
        let before = g.clone();
        assert_eq!(query_bond_row_from_source_index(&g, 1_usize), Ok(0));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_search/CXlower_2147483646
    #[test]
    fn uint_cell_unsigned_consumer_search_cxlower_2147483646_cx_lowering() {
        let g = query_graph(vec![cosmolkit_model::PropertyValue::UInt(2147483646_u32)]);
        let before = g.clone();
        assert_eq!(
            query_bond_row_from_source_index(&g, 2147483646_usize),
            Ok(0)
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_search/CXlower_2147483647
    #[test]
    fn uint_cell_unsigned_consumer_search_cxlower_2147483647_cx_lowering() {
        let g = query_graph(vec![cosmolkit_model::PropertyValue::UInt(2147483647_u32)]);
        let before = g.clone();
        assert_eq!(
            query_bond_row_from_source_index(&g, 2147483647_usize),
            Ok(0)
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_search/CXlower_2147483648
    #[test]
    fn uint_cell_unsigned_consumer_search_cxlower_2147483648_cx_lowering() {
        let g = query_graph(vec![cosmolkit_model::PropertyValue::UInt(2147483648_u32)]);
        let before = g.clone();
        assert_eq!(
            query_bond_row_from_source_index(&g, 2147483648_usize),
            Ok(0)
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UNSIGNED_CONSUMER_search/CXlower_4294967295
    #[test]
    fn uint_cell_unsigned_consumer_search_cxlower_4294967295_cx_lowering() {
        let g = query_graph(vec![cosmolkit_model::PropertyValue::UInt(4294967295_u32)]);
        let before = g.clone();
        assert_eq!(
            query_bond_row_from_source_index(&g, 4294967295_usize),
            Ok(0)
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: CX_FIRST_MATCH
    #[test]
    fn uint_cell_cx_first_match_cx_lowering() {
        let g = query_graph(vec![
            cosmolkit_model::PropertyValue::UInt(0),
            cosmolkit_model::PropertyValue::IntVector(vec![]),
        ]);
        let before = g.clone();
        assert_eq!(query_bond_row_from_source_index(&g, 0), Ok(0));
        assert_eq!(g, before);
    }
}
