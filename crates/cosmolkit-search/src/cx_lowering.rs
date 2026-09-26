//! QueryGraph-native lowering for representation-independent CX records.

use cosmolkit_cx::{
    CxAtomConstraint, CxCoordinateBondKind, CxCountConstraint, CxDataSGroup,
    CxDoubleBondStereoKind, CxEnhancedStereo, CxLinkNode, CxPolymerSGroup, CxRecord,
    CxSGroupHierarchy, CxStereoGroupKind, CxVariableAttachment, CxWedgeBond, CxWedgeDirection,
    ParsedCxExtensions,
};
use cosmolkit_model::{
    AtomId, AtomQueryPredicate, BondDirection, BondId, BondOrder, BondStereo, Conformer3D,
    QueryAtomIdentity, QueryGraph, QueryGraphError, QueryNode, SGroupData, StereoGroup,
    StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, query_substance_groups,
    replace_query_stereo_groups, replace_query_substance_groups,
};
use cosmolkit_types::{ChiralTag, Hybridization};

const QUERY_SCAN_MAGIC_VALUE: u32 = 0xDEAD_BEEF;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CxQueryLoweringError {
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

pub(crate) fn apply_cx_link_nodes_to_query(
    graph: &mut QueryGraph,
    nodes: &[CxLinkNode],
) -> Result<(), CxQueryLoweringError> {
    // RDKit source (verbatim; CXSmilesOps.cpp::parse_linknodes):
    /*
    if (first < last && *first == '.') {
      ++first;
      if (!read_int(first, last, idx1)) {
        return false;
      }
      ++first;
      if (!read_int(first, last, idx2)) {
        return false;
      }
    } else if (VALID_ATIDX(atidx) &&
               mol.getAtomWithIdx(atidx - startAtomIdx)->getDegree() == 2) {
      auto nbrs =
          mol.getAtomNeighbors(mol.getAtomWithIdx(atidx - startAtomIdx));
      idx1 = *nbrs.first;
      nbrs.first++;
      idx2 = *nbrs.first;
    } else if (VALID_ATIDX(atidx)) {
      return false;
    }
    if (first < last && *first == ',') {
      ++first;
    }
    if (VALID_ATIDX(atidx)) {
      if (!accum.empty()) {
        accum += "|";
      }
      accum += (boost::format("%d %d 2 %d %d %d %d") % startReps % endReps %
                (atidx - startAtomIdx + 1) % (idx1 - startAtomIdx + 1) %
                (atidx - startAtomIdx + 1) % (idx2 - startAtomIdx + 1))
                   .str();
    }
    if (!accum.empty()) {
      mol.setProp(common_properties::molFileLinkNodes, accum);
    }
    */
    // RDKit❗✔️: preserve source neighbor order and one-based unsigned values;
    // build locally so a later degree error leaves the prior property intact.
    let mut accum = String::new();
    for node in nodes {
        if node.atom >= graph.num_atoms() {
            continue;
        }

        let [outer_one, outer_two] = match node.outer_atoms {
            Some(outer_atoms) => outer_atoms,
            None => {
                let neighbors = graph.adjacency().get(node.atom).ok_or_else(|| {
                    CxQueryLoweringError::InvalidGraph(
                        "CX link-node center has no adjacency row".to_owned(),
                    )
                })?;
                if neighbors.len() != 2 {
                    return Err(CxQueryLoweringError::InvalidGraph(format!(
                        "CX link-node atom {} has degree {}, expected two when outer atoms are omitted",
                        node.atom,
                        neighbors.len()
                    )));
                }
                [neighbors[0].0, neighbors[1].0]
            }
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
        graph.set_prop("molFileLinkNodes", accum);
    }
    Ok(())
}

fn append_atom_predicate(
    graph: &mut QueryGraph,
    atom: usize,
    predicate: AtomQueryPredicate,
) -> Result<(), CxQueryLoweringError> {
    // BEGIN RDKIT CPP FUNCTION QueryOps::replaceAtomWithQueryAtom
    // RDKit✔️✔️: Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   if (atom->hasQuery()) {
    // RDKit✔️✔️:     return atom;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   QueryAtom qa(*atom);
    // RDKit✔️✔️:   unsigned int idx = atom->getIdx();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (atom->hasProp(common_properties::_hasMassQuery)) {
    // RDKit✔️✔️:     qa.expandQuery(makeAtomMassQuery(static_cast<int>(atom->getMass())));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol->replaceAtom(idx, &qa);
    // RDKit✔️✔️:   return mol->getAtomWithIdx(idx);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryOps::replaceAtomWithQueryAtom
    // The detached QueryGraph already stores the carrier, its predicate, and
    // whether that predicate came from an ordinary atom in the corresponding
    // typed value. Mutating the predicate marks that origin explicit without
    // replacing the raw query identity or carrier state.
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
    append_atom_predicate(graph, atom, predicate)
}

const CX_LABELS_PROCESSED_PROP: &str = "_cxsmilesLabelsProcessed";

/// Apply RDKit's deferred CX label replacement to the detached query graph.
/// The temporary guard remains set until the enclosing CX application calls
/// `finish_cx_smiles_labels`, so SGroup helpers can invoke this before attach.
pub(crate) fn process_cx_smiles_labels(graph: &mut QueryGraph) -> Result<(), CxQueryLoweringError> {
    // BEGIN RDKIT CPP FUNCTION processCXSmilesLabels
    // RDKit❗✔️: void processCXSmilesLabels(RWMol &mol) {
    // RDKit❗✔️:   if (mol.hasProp("_cxsmilesLabelsProcessed")) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (auto atom : mol.atoms()) {
    // RDKit❗✔️:     std::string symb = "";
    // RDKit❗✔️:     if (atom->getPropIfPresent(common_properties::atomLabel, symb)) {
    // RDKit❗✔️:       atom->clearProp(common_properties::dummyLabel);
    // RDKit❗✔️:       if (symb == "star_e") {
    // RDKit❗✔️:         addquery(makeAtomNullQuery(), symb, mol, atom->getIdx());
    // RDKit❗✔️:       } else if (symb == "Q_e") {
    // RDKit❗✔️:         addquery(makeQAtomQuery(), symb, mol, atom->getIdx());
    // RDKit❗✔️:       } else if (symb == "QH_p") {
    // RDKit❗✔️:         addquery(makeQHAtomQuery(), symb, mol, atom->getIdx());
    // RDKit❗✔️:       } else if (symb == "AH_p") {
    // RDKit❗✔️:         addquery(makeAHAtomQuery(), symb, mol, atom->getIdx());
    // RDKit❗✔️:       } else if (symb == "X_p") {
    // RDKit❗✔️:         addquery(makeXAtomQuery(), symb, mol, atom->getIdx());
    // RDKit❗✔️:       } else if (symb == "XH_p") {
    // RDKit❗✔️:         addquery(makeXHAtomQuery(), symb, mol, atom->getIdx());
    // RDKit❗✔️:       } else if (symb == "M_p") {
    // RDKit❗✔️:         addquery(makeMAtomQuery(), symb, mol, atom->getIdx());
    // RDKit❗✔️:       } else if (symb == "MH_p") {
    // RDKit❗✔️:         addquery(makeMHAtomQuery(), symb, mol, atom->getIdx());
    // RDKit❗✔️:       } else if (std::find(pseudoatoms_p.begin(), pseudoatoms_p.end(), symb) !=
    // RDKit❗✔️:                  pseudoatoms_p.end()) {
    // RDKit❗✔️:         atom->setProp(common_properties::dummyLabel,
    // RDKit❗✔️:                       symb.substr(0, symb.size() - 2));
    // RDKit❗✔️:         atom->clearProp(common_properties::atomLabel);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (atom->getAtomicNum() == 0 && !atom->hasQuery() &&
    // RDKit❗✔️:                !atom->getIsotope() && atom->getSymbol() == "*") {
    // RDKit❗✔️:       addquery(makeAAtomQuery(), "", mol, atom->getIdx());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   mol.setProp("_cxsmilesLabelsProcessed", 1, true);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION processCXSmilesLabels
    // These predicates reuse the existing source-anchored QueryOps factories
    // in query_behavior; this dispatcher does not build parallel query trees.
    if graph.prop(CX_LABELS_PROCESSED_PROP).is_some() {
        return Ok(());
    }

    for atom_index in 0..graph.num_atoms() {
        let label = graph
            .atom(atom_index)
            .and_then(|atom| atom.prop("atomLabel"))
            .map(str::to_owned);

        if let Some(label) = label {
            graph
                .atom_mut(atom_index)
                .ok_or(CxQueryLoweringError::AtomIndex { index: atom_index })?
                .clear_prop("dummyLabel");

            let predicate = match label.as_str() {
                "star_e" => Some(crate::query_behavior::make_atom_null_query()),
                "Q_e" => Some(crate::query_behavior::make_q_atom_query()),
                "QH_p" => Some(crate::query_behavior::make_q_h_atom_query()),
                "AH_p" => Some(crate::query_behavior::make_a_h_atom_query()),
                "X_p" => Some(crate::query_behavior::make_x_atom_query()),
                "XH_p" => Some(crate::query_behavior::make_x_h_atom_query()),
                "M_p" => Some(crate::query_behavior::make_m_atom_query()),
                "MH_p" => Some(crate::query_behavior::make_m_h_atom_query()),
                _ => None,
            };
            if let Some(predicate) = predicate {
                replace_query_atom_like_rdkit(graph, atom_index, predicate, Some(&label))?;
            } else if let Some(dummy_label) = match label.as_str() {
                "Pol_p" => Some("Pol"),
                "Mod_p" => Some("Mod"),
                _ => None,
            } {
                let atom = graph
                    .atom_mut(atom_index)
                    .ok_or(CxQueryLoweringError::AtomIndex { index: atom_index })?;
                atom.set_prop("dummyLabel", dummy_label)
                    .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
                atom.clear_prop("atomLabel");
            }
        } else if graph.atom(atom_index).is_some_and(|atom| {
            atom.atomic_number() == 0
                && atom.isotope().is_none()
                && atom.predicate_is_carrier_derived()
        }) {
            replace_query_atom_like_rdkit(
                graph,
                atom_index,
                crate::query_behavior::make_a_atom_query(),
                None,
            )?;
        }
    }

    // Local complexity review: this loop visits each atom once. A query-label
    // replacement copies that atom's property map and reconstructs the same
    // stereo-group members that RWMol::replaceAtom scans; both costs are
    // linear in the copied properties and affected group membership.
    graph.set_prop(CX_LABELS_PROCESSED_PROP, "1");
    Ok(())
}

/// Finish one complete CX application, matching parseCXExtensions' guard
/// cleanup after its final label pass.
pub(crate) fn finish_cx_smiles_labels(graph: &mut QueryGraph) -> Result<(), CxQueryLoweringError> {
    process_cx_smiles_labels(graph)?;
    graph.clear_prop(CX_LABELS_PROCESSED_PROP);
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
        field_display: Some("    0.0000    0.0000    DR    ALL  0       0".to_owned()),
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
    group.set_prop("_cxsmilesindex", cx_sequence_id.to_string());
    group.set_prop("index", source_index.to_string());
    group.set_prop("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
    if !data.field_name.is_empty() {
        group.set_prop("FIELDNAME", data.field_name.clone());
    }
    if !data.data.is_empty() {
        group.set_prop("DATAFIELDS", data.data.clone());
        group.push_data_field(data.data.clone());
    }
    if !data.query_op.is_empty() {
        group.set_prop("QUERYOP", data.query_op.clone());
    }
    if !data.field_info.is_empty() {
        group.set_prop("FIELDINFO", data.field_info.clone());
    }
    if !data.field_tag.is_empty() {
        group.set_prop("FIELDTAG", data.field_tag.clone());
    }
    if let Some(coordinates) = &data.coordinates {
        group.set_prop("COORDS", coordinates.clone());
    }
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
            .set_prop("_MolFileBondEndPts", end_points.clone())
            .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
        bond.bond_mut()
            .set_prop("_MolFileBondAttach", "ANY")
            .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
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
    // RDKit❗🔝: template <typename Iterator>
    // RDKit❗🔝: bool parse_wedged_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit❗🔝:                         unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit❗🔝:   // these look like: CC(O)Cl |w:1.0|
    // RDKit❗🔝:   // also wD and wU for down and up wedges.
    // RDKit❗🔝:   //
    // RDKit❗🔝:   // We do not end up using this to set stereochemistry, but the relevant bond
    // RDKit❗🔝:   // properties are set in case client code wants to do something with the
    // RDKit❗🔝:   // information.
    // RDKit❗🔝:   if (first >= last || *first != 'w' || first + 1 >= last) {
    // RDKit❗🔝:     return false;
    // RDKit❗🔝:   }
    // RDKit❗🔝:   ++first;
    // RDKit❗🔝:   Bond::BondDir state = Bond::BondDir::NONE;
    // RDKit❗🔝:   unsigned int cfg = 0;
    // RDKit❗🔝:   switch (*first) {
    // RDKit❗🔝:     case ':':
    // RDKit❗🔝:       state = Bond::BondDir::UNKNOWN;
    // RDKit❗🔝:       cfg = 2;
    // RDKit❗🔝:       break;
    // RDKit❗🔝:     case 'U':
    // RDKit❗🔝:       state = Bond::BondDir::BEGINWEDGE;
    // RDKit❗🔝:       cfg = 1;
    // RDKit❗🔝:       ++first;
    // RDKit❗🔝:       break;
    // RDKit❗🔝:     case 'D':
    // RDKit❗🔝:       state = Bond::BondDir::BEGINDASH;
    // RDKit❗🔝:       cfg = 3;
    // RDKit❗🔝:       ++first;
    // RDKit❗🔝:       break;
    // RDKit❗🔝:     default:
    // RDKit❗🔝:       break;
    // RDKit❗🔝:   }
    // RDKit❗🔝:   if (state == Bond::BondDir::NONE || first >= last || first + 1 >= last ||
    // RDKit❗🔝:       *first != ':') {
    // RDKit❗🔝:     return false;
    // RDKit❗🔝:   }
    // RDKit❗🔝:   ++first;
    // RDKit❗🔝:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit❗🔝:     unsigned int atomIdx;
    // RDKit❗🔝:     if (!read_int(first, last, atomIdx)) {
    // RDKit❗🔝:       return false;
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (first < last && *first == '.') {
    // RDKit❗🔝:       ++first;
    // RDKit❗🔝:     } else {
    // RDKit❗🔝:       BOOST_LOG(rdWarningLog) << "improperly formatted w block" << std::endl;
    // RDKit❗🔝:       return false;
    // RDKit❗🔝:     }
    // RDKit❗🔝:     unsigned int bondIdx;
    // RDKit❗🔝:     if (!read_int(first, last, bondIdx)) {
    // RDKit❗🔝:       return false;
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (VALID_ATIDX(atomIdx) && VALID_BNDIDX(bondIdx)) {
    // RDKit❗🔝:       auto atom = mol.getAtomWithIdx(atomIdx - startAtomIdx);
    // RDKit❗🔝:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit❗🔝:       if (!bond) {
    // RDKit❗🔝:         BOOST_LOG(rdWarningLog)
    // RDKit❗🔝:             << "bond " << bondIdx << " not found, wedge from atom " << atomIdx
    // RDKit❗🔝:             << " cannot be applied." << std::endl;
    // RDKit❗🔝:         return false;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       if (bond->hasProp(common_properties::_MolFileBondCfg)) {
    // RDKit❗🔝:         BOOST_LOG(rdWarningLog)
    // RDKit❗🔝:             << "w block attempts to set wedging on bond " << bond->getIdx()
    // RDKit❗🔝:             << " more than once." << std::endl;
    // RDKit❗🔝:         return false;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       if (atom->getIdx() != bond->getBeginAtomIdx()) {
    // RDKit❗🔝:         if (atom->getIdx() != bond->getEndAtomIdx()) {
    // RDKit❗🔝:           BOOST_LOG(rdWarningLog)
    // RDKit❗🔝:               << "atom " << atomIdx << " is not associated with bond "
    // RDKit❗🔝:               << bondIdx << "(" << bond->getBeginAtomIdx() + startAtomIdx << "-"
    // RDKit❗🔝:               << bond->getEndAtomIdx() + startAtomIdx << ")"
    // RDKit❗🔝:               << " in w block" << std::endl;
    // RDKit❗🔝:           return false;
    // RDKit❗🔝:         }
    // RDKit❗🔝:         auto eidx = bond->getBeginAtomIdx();
    // RDKit❗🔝:         bond->setBeginAtomIdx(atom->getIdx());
    // RDKit❗🔝:         bond->setEndAtomIdx(eidx);
    // RDKit❗🔝:       }
    // RDKit❗🔝:       bond->setProp(common_properties::_MolFileBondCfg, cfg);
    // RDKit❗🔝:       bond->setBondDir(state);
    // RDKit❗🔝:       if (cfg == 2 && canHaveDirection(*bond)) {
    // RDKit❗🔝:         bond->getBeginAtom()->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit❗🔝:         mol.setProp(detail::_needsDetectBondStereo, 1);
    // RDKit❗🔝:       }
    // RDKit❗🔝:       if ((cfg == 1 || cfg == 3) && canHaveDirection(*bond)) {
    // RDKit❗🔝:         mol.setProp(detail::_needsDetectAtomStereo, 1);
    // RDKit❗🔝:       }
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (first < last && *first == ',') {
    // RDKit❗🔝:       ++first;
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return true;
    // RDKit❗🔝: }
    // RDKit source helper (verbatim; CXSmilesOps.cpp::get_bond_with_smiles_idx):
    // RDKit❗🔝: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
    // RDKit❗🔝:   for (auto bnd : mol.bonds()) {
    // RDKit❗🔝:     unsigned int smilesIdx;
    // RDKit❗🔝:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
    // RDKit❗🔝:         smilesIdx == idx) {
    // RDKit❗🔝:       return bnd;
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return nullptr;
    // RDKit❗🔝: }
    // RDKit source helper (verbatim; Bond.h::canHaveDirection):
    // RDKit❗🔝: inline bool canHaveDirection(const Bond &bond) {
    // RDKit❗🔝:   auto bondType = bond.getBondType();
    // RDKit❗🔝:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit❗🔝: }
    // CX source IDs equal QueryGraph bond insertion order; direct indexing
    // replaces RDKit's property scan with one vector lookup.
    if wedge.atom >= graph.num_atoms() || wedge.bond >= graph.num_bonds() {
        return Ok(());
    }
    let bond = graph
        .bonds_mut()
        .get_mut(wedge.bond)
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
        CxWedgeDirection::Unknown => ("2", BondDirection::Unknown),
        CxWedgeDirection::BeginWedge => ("1", BondDirection::BeginWedge),
        CxWedgeDirection::BeginDash => ("3", BondDirection::BeginDash),
    };
    bond.bond_mut()
        .set_prop("_MolFileBondCfg", configuration)
        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
    bond.bond_mut().set_direction(direction);
    if wedge.direction == CxWedgeDirection::Unknown && can_have_direction {
        graph
            .atom_mut(wedge.atom)
            .ok_or(CxQueryLoweringError::AtomIndex { index: wedge.atom })?
            .set_chiral_tag(ChiralTag::Unspecified);
        graph.set_prop("_needsDetectBondStereo", "1");
    }
    if matches!(
        wedge.direction,
        CxWedgeDirection::BeginWedge | CxWedgeDirection::BeginDash
    ) && can_have_direction
    {
        graph.set_prop("_needsDetectAtomStereo", "1");
    }
    Ok(())
}

pub(crate) fn apply_cx_double_bond_stereo_to_query(
    graph: &mut QueryGraph,
    bond_index: usize,
    stereo: CxDoubleBondStereoKind,
) -> Result<(), CxQueryLoweringError> {
    // RDKit source (verbatim; CXSmilesOps.cpp::parse_doublebond_stereo):
    // RDKit❗🔝: template <typename Iterator>
    // RDKit❗🔝: bool parse_doublebond_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit❗🔝:                              unsigned int, unsigned int startBondIdx,
    // RDKit❗🔝:                              Bond::BondStereo stereo) {
    // RDKit❗🔝:   // these look like: C1CCCC/C=C/CCC1 |ctu:5|
    // RDKit❗🔝:   // also c and t for cis or trans
    // RDKit❗🔝:   while (first < last && *first != ':') {
    // RDKit❗🔝:     ++first;
    // RDKit❗🔝:   }
    // RDKit❗🔝:   if (first >= last || *first != ':') {
    // RDKit❗🔝:     return false;
    // RDKit❗🔝:   }
    // RDKit❗🔝:   ++first;
    // RDKit❗🔝:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit❗🔝:     unsigned int bondIdx;
    // RDKit❗🔝:     if (!read_int(first, last, bondIdx)) {
    // RDKit❗🔝:       return false;
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (VALID_BNDIDX(bondIdx)) {
    // RDKit❗🔝:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit❗🔝:       if (!bond) {
    // RDKit❗🔝:         BOOST_LOG(rdWarningLog)
    // RDKit❗🔝:             << "bond " << bondIdx
    // RDKit❗🔝:             << " not found, cannot mark as stereo double bond." << std::endl;
    // RDKit❗🔝:         return false;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       bool useCXOrdering = true;
    // RDKit❗🔝:       Chirality::detail::setStereoForBond(mol, bond, stereo, useCXOrdering);
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (first < last && *first == ',') {
    // RDKit❗🔝:       ++first;
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return true;
    // RDKit❗🔝: }
    // RDKit source helper (verbatim; Chirality.cpp::setStereoForBond):
    // RDKit❗🔝: void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo,
    // RDKit❗🔝:                       bool useCXSmilesOrdering) {
    // RDKit❗🔝:   // NOTE:  moved from parse_doublebond_stereo CXSmilesOps
    // RDKit❗🔝:   // IF useCXSmilesOrdering is true, the cis/trans/unknown marker will be
    // RDKit❗🔝:   // assigned relative to the lowest-numbered neighbor of each double bond atom.
    // RDKit❗🔝:   // Otherwise it uses the lowest-numbered neighbor on the lower-numbered atom
    // RDKit❗🔝:   // of the double bond and the highest-numbered neighbor on the higher-numbered
    // RDKit❗🔝:   // atom
    // RDKit❗🔝:   auto begAtom = bond->getBeginAtom();
    // RDKit❗🔝:   auto endAtom = bond->getEndAtom();
    // RDKit❗🔝:   if (begAtom->getIdx() > endAtom->getIdx()) {
    // RDKit❗🔝:     std::swap(begAtom, endAtom);
    // RDKit❗🔝:   }
    // RDKit❗🔝:   if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    // RDKit❗🔝:     unsigned int begControl = mol.getNumAtoms();
    // RDKit❗🔝:     for (auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit❗🔝:       if (nbr == endAtom) {
    // RDKit❗🔝:         continue;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       begControl = std::min(nbr->getIdx(), begControl);
    // RDKit❗🔝:     }
    // RDKit❗🔝:     unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
    // RDKit❗🔝:     for (auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit❗🔝:       if (nbr == begAtom) {
    // RDKit❗🔝:         continue;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
    // RDKit❗🔝:                                        : std::max(nbr->getIdx(), endControl);
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (begAtom != bond->getBeginAtom()) {
    // RDKit❗🔝:       std::swap(begControl, endControl);
    // RDKit❗🔝:     }
    // RDKit❗🔝:     bond->setStereoAtoms(begControl, endControl);
    // RDKit❗🔝:     bond->setStereo(stereo);
    // RDKit❗🔝:     mol.setProp("_needsDetectBondStereo", 1);
    // RDKit❗🔝:   }
    // RDKit❗🔝: }
    // CX source IDs equal QueryGraph insertion order, so direct indexing
    // removes the source scan; adjacency iteration retains the same minima.
    if bond_index >= graph.num_bonds() {
        return Ok(());
    }
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
    graph.set_prop("_needsDetectBondStereo", "1");
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
    let kind = match polymer.type_code.as_str() {
        "n" => SubstanceGroupKind::StructuralRepeatUnit,
        "mon" => SubstanceGroupKind::Monomer,
        "mer" => SubstanceGroupKind::Mer,
        "co" => SubstanceGroupKind::Copolymer,
        "xl" => SubstanceGroupKind::Crosslink,
        "mod" => SubstanceGroupKind::Modification,
        "mix" => SubstanceGroupKind::MixtureComponent,
        "f" => SubstanceGroupKind::Formulation,
        "any" => SubstanceGroupKind::AnyPolymer,
        "gen" => SubstanceGroupKind::Generic("GEN".to_owned()),
        "c" => SubstanceGroupKind::Generic("COM".to_owned()),
        "grf" => SubstanceGroupKind::Graft,
        "alt" | "ran" | "blk" => SubstanceGroupKind::Copolymer,
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
    group.set_prop("_cxsmilesindex", cx_sequence_id.to_string());
    if !polymer.label.is_empty() {
        group.set_label(polymer.label.clone());
        group.set_prop("LABEL", polymer.label.clone());
    }
    if !polymer.connect.is_empty() {
        group.set_prop("CONNECT", polymer.connect.clone());
    }
    match polymer.type_code.as_str() {
        "alt" => {
            group.set_subtype("ALT");
            group.set_prop("SUBTYPE", "ALT");
        }
        "ran" => {
            group.set_subtype("RAN");
            group.set_prop("SUBTYPE", "RAN");
        }
        "blk" => {
            group.set_subtype("BLO");
            group.set_prop("SUBTYPE", "BLO");
        }
        _ => {}
    }
    cosmolkit_core::finalize_polymer_sgroup(
        &mut group,
        (!polymer.connect.is_empty()).then_some(polymer.connect.as_str()),
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
    group.set_prop("index", source_index.to_string());
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
    let parent_property_id = match parent.props().get("index") {
        Some(value) => value.parse::<u32>().map_err(|_| {
            CxQueryLoweringError::InvalidGraph(
                "query SGroup index property is not a source unsigned integer".to_owned(),
            )
        })?,
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
    child
        .props()
        .get("index")
        .ok_or_else(|| {
            CxQueryLoweringError::InvalidGraph(
                "query SGroup child is missing its source index property".to_owned(),
            )
        })?
        .parse::<u32>()
        .map_err(|_| {
            CxQueryLoweringError::InvalidGraph(
                "query SGroup child index property is not a source unsigned integer".to_owned(),
            )
        })?;
    let child = &mut groups[child_index];
    child.set_parent(parent_id);
    child.set_prop("PARENT", parent_property_id.to_string());
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
        let Some(value) = group.props().get("_cxsmilesindex") else {
            continue;
        };
        let sequence_id = value.parse::<u32>().map_err(|_| {
            CxQueryLoweringError::InvalidGraph(
                "query SGroup CX sequence property is not a source unsigned integer".to_owned(),
            )
        })?;
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
    predicate: QueryNode<AtomQueryPredicate>,
    label: Option<&str>,
) -> Result<(), CxQueryLoweringError> {
    // BEGIN RDKIT CPP FUNCTION CXSmilesOps::addquery
    // RDKit❗✔️: void addquery(Q *qry, std::string symbol, RDKit::RWMol &mol, unsigned int idx) {
    // RDKit❗✔️:   PRECONDITION(qry, "bad query");
    // RDKit❗✔️:   auto *qa = new QueryAtom(0);
    // RDKit❗✔️:   qa->setQuery(qry);
    // RDKit❗✔️:   qa->setNoImplicit(true);
    // RDKit❗✔️:   bool updateLabel = false;
    // RDKit❗✔️:   bool preserveProps = true;
    // RDKit❗✔️:   mol.replaceAtom(idx, qa, updateLabel, preserveProps);
    // RDKit❗✔️:   if (symbol != "") {
    // RDKit❗✔️:     mol.getAtomWithIdx(idx)->setProp(RDKit::common_properties::atomLabel,
    // RDKit❗✔️:                                      symbol);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   delete qa;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION CXSmilesOps::addquery
    // BEGIN RDKIT CPP FUNCTION RWMol::replaceAtom
    // RDKit❗✔️: if (preserveProps) atom_p->updateProps(*d_graph[vd], false);
    // RDKit❗✔️: for (auto &group : d_stereo_groups) {
    // RDKit❗✔️:   auto groupId = group.getReadId();
    // RDKit❗✔️:   auto atoms = group.getAtoms();
    // RDKit❗✔️:   auto bonds = group.getBonds();
    // RDKit❗✔️:     auto aiter = std::find(atoms.begin(), atoms.end(), orig_p);
    // RDKit❗✔️:     while (aiter != atoms.end()) {
    // RDKit❗✔️:       *aiter = atom_p;
    // RDKit❗✔️:       ++aiter;
    // RDKit❗✔️:       aiter = std::find(aiter, atoms.end(), orig_p);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     group = StereoGroup(group.getGroupType(), std::move(atoms),
    // RDKit❗✔️:                         std::move(bonds), groupId);
    // END RDKIT CPP FUNCTION RWMol::replaceAtom
    let source = graph
        .atom(atom_index)
        .ok_or(CxQueryLoweringError::AtomIndex { index: atom_index })?;
    let id = source.id();
    let mut replacement = source
        .clone()
        .with_identity(QueryAtomIdentity::from_atomic_number(0));
    replacement.set_predicate(predicate);
    replacement.set_formal_charge(0);
    replacement.set_explicit_hydrogens(0);
    replacement.set_chiral_tag(ChiralTag::Unspecified);
    replacement.set_chiral_permutation(None);
    replacement.set_unknown_stereo(false);
    replacement.set_mol_parity(None);
    replacement.set_mol_inversion_flag(None);
    replacement.set_implicit_hydrogen(false);
    replacement.set_tracked_isotopic_hydrogens(Vec::new());
    replacement.set_aromatic(false);
    replacement.set_isotope(None);
    replacement.set_atom_map(None);
    replacement.set_no_implicit(true);
    replacement.set_radical_electrons(0);
    replacement.set_hybridization(Hybridization::Unspecified);
    replacement.set_pdb_residue_info(None);
    replacement.set_template_attachment_order(None);
    if let Some(label) = label.filter(|label| !label.is_empty()) {
        replacement
            .set_prop("atomLabel", label)
            .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
    }
    debug_assert_eq!(replacement.id(), id);
    graph.atoms_mut()[atom_index] = replacement;

    let stereo_groups = graph
        .stereo_groups()
        .iter()
        .map(|group| {
            let replacement =
                StereoGroup::new(group.kind(), group.atoms().to_vec(), group.bonds().to_vec());
            if let Some(read_id) = group.id() {
                replacement.with_id(read_id)
            } else {
                replacement
            }
        })
        .collect();
    replace_query_stereo_groups(graph, stereo_groups)
        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;

    // The model currently exposes the stereo read ID only. RDKit reconstructs
    // each group with getReadId(), which resets its separate write ID to zero;
    // keep this reconstruction so tab4's upcoming write-ID field gets the
    // same zero-default behavior when integrated.
    // Local complexity review: the source also copies each stereo group's
    // ordered atom/bond vectors per replacement. This reconstruction has the
    // same O(groups + memberships) work and preserves duplicate order.
    Ok(())
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
                // BEGIN RDKIT CPP FUNCTION parse_coords
                // RDKit❗✔️: auto *conf = new Conformer(mol.getNumAtoms());
                // RDKit❗✔️: mol.addConformer(conf);
                // RDKit❗✔️: conf->setId(confIdx);
                // RDKit❗✔️: unsigned int atIdx = 0;
                // RDKit❗✔️: while (first <= last && *first != ')') {
                // RDKit❗✔️:   if (VALID_ATIDX(atIdx)) {
                // RDKit❗✔️:     conf->setAtomPos(atIdx - startAtomIdx, pt);
                // RDKit❗✔️:   }
                // RDKit❗✔️:   ++atIdx;
                // RDKit❗✔️: }
                // RDKit❗✔️: if (is3D && hasNonZeroZCoords(*conf)) {
                // RDKit❗✔️:   conf->set3D(true);
                // RDKit❗✔️: } else {
                // RDKit❗✔️:   conf->set3D(false);
                // RDKit❗✔️: }
                // END RDKIT CPP FUNCTION parse_coords
                // BEGIN RDKIT CPP MACRO VALID_ATIDX
                // RDKit❗✔️: #define VALID_ATIDX(_atidx_) \
                // RDKit❗✔️:   ((_atidx_) >= startAtomIdx && (_atidx_) < startAtomIdx + mol.getNumAtoms())
                // END RDKIT CPP MACRO VALID_ATIDX
                // BEGIN RDKIT CPP FUNCTION Conformer::Conformer(unsigned int)
                // RDKit❗✔️: Conformer(unsigned int numAtoms)
                // RDKit❗✔️:     : d_positions(numAtoms, RDGeom::Point3D(0.0, 0.0, 0.0)) {}
                // END RDKIT CPP FUNCTION Conformer::Conformer(unsigned int)
                // BEGIN RDKIT CPP FUNCTION Conformer::setAtomPos
                // RDKit❗✔️: inline void setAtomPos(unsigned int atomId, const RDGeom::Point3D &position) {
                // RDKit❗✔️:   if (atomId == std::numeric_limits<unsigned int>::max()) {
                // RDKit❗✔️:     throw ValueErrorException("atom index overflow");
                // RDKit❗✔️:   }
                // RDKit❗✔️:   if (atomId >= d_positions.size()) {
                // RDKit❗✔️:     d_positions.resize(atomId + 1, RDGeom::Point3D(0.0, 0.0, 0.0));
                // RDKit❗✔️:   }
                // RDKit❗✔️:   d_positions[atomId] = position;
                // RDKit❗✔️: }
                // END RDKIT CPP FUNCTION Conformer::setAtomPos
                // BEGIN RDKIT CPP FUNCTION hasNonZeroZCoords
                // RDKit❗✔️: inline bool hasNonZeroZCoords(const Conformer &conf) {
                // RDKit❗✔️:   constexpr double zeroTol = 1e-3;
                // RDKit❗✔️:   for (auto p : conf.getPositions()) {
                // RDKit❗✔️:     if (std::abs(p.z) > zeroTol) {
                // RDKit❗✔️:       return true;
                // RDKit❗✔️:     }
                // RDKit❗✔️:   }
                // RDKit❗✔️:   return false;
                // RDKit❗✔️: }
                // END RDKIT CPP FUNCTION hasNonZeroZCoords
                // RDKit uses a zero-filled atom-sized conformer and increments
                // the source slot for every row. Here row ordinals are already
                // decoded, and query graphs use startAtomIdx == 0.
                let atom_count = graph.num_atoms();
                let mut values = vec![[0.0; 3]; atom_count];
                // This allocates one atom-sized destination and scans each
                // parsed row once; rows rejected by VALID_ATIDX do not write.
                for (index, value) in coordinates.values.iter().enumerate() {
                    if index >= atom_count {
                        continue;
                    }
                    if let Some(value) = value {
                        values[index] = *value;
                    }
                }
                // Ignore z values belonging to source rows the conformer
                // rejected as invalid atom indices before setting its 3D flag.
                let is_3d = coordinates.is_3d && values.iter().any(|point| point[2].abs() > 1e-3);
                graph
                    .add_conformer_3d(Conformer3D::new(coordinates.conformer, values, is_3d))
                    .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
            }
            CxRecord::AtomLabels(values) => {
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
                for (index, value) in values.iter().enumerate() {
                    if index >= graph.num_atoms() {
                        continue;
                    }
                    let Some(value) = value.as_ref().filter(|value| !value.is_empty()) else {
                        continue;
                    };
                    let atom = graph
                        .atom_mut(index)
                        .ok_or(CxQueryLoweringError::AtomIndex { index })?;
                    atom.set_prop("atomLabel", value)
                        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
                }
            }
            CxRecord::AtomValues(values) => {
                // BEGIN RDKIT CPP FUNCTION parse_atom_values
                // RDKit✔️✔️: std::string tkn = read_text_to(first, last, ";$");
                // RDKit✔️✔️: if (tkn != "" && VALID_ATIDX(atIdx)) {
                // RDKit✔️✔️:   mol.getAtomWithIdx(atIdx)->setProp(
                // RDKit✔️✔️:       RDKit::common_properties::molFileValue, tkn);
                // RDKit✔️✔️: }
                // RDKit✔️✔️: ++atIdx;
                // END RDKIT CPP FUNCTION parse_atom_values
                // Values use the same already-decoded ordered slots as labels.
                for (index, value) in values.iter().enumerate() {
                    if index >= graph.num_atoms() {
                        continue;
                    }
                    let Some(value) = value.as_ref().filter(|value| !value.is_empty()) else {
                        continue;
                    };
                    let atom = graph
                        .atom_mut(index)
                        .ok_or(CxQueryLoweringError::AtomIndex { index })?;
                    atom.set_prop("molFileValue", value)
                        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
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
                    atom.set_prop(property.name.clone(), property.value.clone())
                        .map_err(|error| CxQueryLoweringError::InvalidGraph(error.to_string()))?;
                }
            }
            CxRecord::CoordinateBonds(annotation) => {
                let order = match annotation.kind {
                    CxCoordinateBondKind::Dative => BondOrder::Dative,
                    CxCoordinateBondKind::Hydrogen => BondOrder::Hydrogen,
                };
                for reference in &annotation.bonds {
                    // RDKit❗🔝: if (VALID_ATIDX(aidx) && VALID_BNDIDX(bidx)) {
                    // RDKit❗🔝:   auto bnd = get_bond_with_smiles_idx(mol, bidx - startBondIdx);
                    // RDKit❗🔝:   if (!bnd || (bnd->getBeginAtomIdx() != aidx - startAtomIdx &&
                    // RDKit❗🔝:                bnd->getEndAtomIdx() != aidx - startAtomIdx)) {
                    // RDKit❗🔝:     return false;
                    // RDKit❗🔝:   }
                    // RDKit❗🔝:   bnd->setBondType(typ);
                    // RDKit❗🔝:   if (bnd->getBeginAtomIdx() != aidx - startAtomIdx) {
                    // RDKit❗🔝:     unsigned int tmp = bnd->getBeginAtomIdx();
                    // RDKit❗🔝:     bnd->setBeginAtomIdx(aidx - startAtomIdx);
                    // RDKit❗🔝:     bnd->setEndAtomIdx(tmp);
                    // RDKit❗🔝:   }
                    // RDKit❗🔝: }
                    // RDKit❗🔝: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
                    // RDKit❗🔝:   for (auto bnd : mol.bonds()) {
                    // RDKit❗🔝:     unsigned int smilesIdx;
                    // RDKit❗🔝:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
                    // RDKit❗🔝:         smilesIdx == idx) {
                    // RDKit❗🔝:       return bnd;
                    // RDKit❗🔝:     }
                    // RDKit❗🔝:   }
                    // RDKit❗🔝:   return nullptr;
                    // RDKit❗🔝: }
                    // The detached parser assigns source indices in bond
                    // insertion order; direct indexing removes the source scan.
                    if reference.atom >= graph.num_atoms() || reference.bond >= graph.num_bonds() {
                        continue;
                    }
                    let bond = graph.bonds_mut().get_mut(reference.bond).ok_or(
                        CxQueryLoweringError::BondIndex {
                            index: reference.bond,
                        },
                    )?;
                    if bond.begin().index() != reference.atom
                        && bond.end().index() != reference.atom
                    {
                        return Err(CxQueryLoweringError::BondAtomMismatch {
                            atom: reference.atom,
                            bond: reference.bond,
                        });
                    }
                    let begin = bond.begin();
                    bond.bond_mut().set_order(order);
                    if begin.index() != reference.atom {
                        bond.bond_mut()
                            .set_endpoints(AtomId::new(reference.atom), begin);
                    }
                }
            }
            CxRecord::ZeroBonds(indices) => {
                for &index in indices {
                    // RDKit❗🔝: if (VALID_BNDIDX(bondIdx)) {
                    // RDKit❗🔝:   auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
                    // RDKit❗🔝:   if (!bond) {
                    // RDKit❗🔝:     return false;
                    // RDKit❗🔝:   }
                    // RDKit❗🔝:   bond->setBondType(Bond::ZERO);
                    // RDKit❗🔝: }
                    // RDKit❗🔝: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
                    // RDKit❗🔝:   for (auto bnd : mol.bonds()) {
                    // RDKit❗🔝:     unsigned int smilesIdx;
                    // RDKit❗🔝:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
                    // RDKit❗🔝:         smilesIdx == idx) {
                    // RDKit❗🔝:       return bnd;
                    // RDKit❗🔝:     }
                    // RDKit❗🔝:   }
                    // RDKit❗🔝:   return nullptr;
                    // RDKit❗🔝: }
                    // Parser-order bond IDs make direct vector lookup remove
                    // the pinned source's linear scan without changing order.
                    // The detached parser assigns this source index in bond
                    // insertion order; direct indexing avoids this linear scan.
                    if index >= graph.num_bonds() {
                        continue;
                    }
                    let bond = graph
                        .bonds_mut()
                        .get_mut(index)
                        .ok_or(CxQueryLoweringError::BondIndex { index })?;
                    bond.bond_mut().set_order(BondOrder::Zero);
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
        Atom, AtomSpec, BondId, BondSpec, QueryAtom, SGroupBondRole, SGroupCState, SGroupData,
        SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, query_substance_groups,
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
            Default::default(),
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
                        field_name: Some("FIELD".to_owned()),
                        field_type: Some("S".to_owned()),
                        values: vec!["one".to_owned(), "two".to_owned()],
                        ..SGroupData::default()
                    })
                    .with_cstates(vec![SGroupCState::new(BondId::new(0), [0.25, 0.5, 0.75])])
                    .with_prop("origin", "existing")
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
                vec![CxRecord::AtomLabels(vec![Some(label.to_owned()), None])],
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
            assert_eq!(atom.atom_map(), None, "{label}");
            assert_eq!(atom.prop("atomLabel"), Some(label), "{label}");
            assert_eq!(atom.prop("dummyLabel"), None, "{label}");
            assert_eq!(atom.prop("sourceProperty"), Some("retained"), "{label}");
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
                vec![CxRecord::AtomLabels(vec![Some(label.to_owned()), None])],
                0,
            );

            apply_cx_to_query_graph(&mut query, &parsed).expect("source CX label lowering");

            let atom = query.atom(0).expect("unreplaced query atom");
            assert_eq!(atom.identity(), source_identity, "{label}");
            assert_eq!(atom.predicate(), &source_predicate, "{label}");
            assert_eq!(atom.formal_charge(), 1, "{label}");
            assert_eq!(atom.prop("atomLabel"), expected_label, "{label}");
            assert_eq!(atom.prop("dummyLabel"), expected_dummy, "{label}");
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
                CxRecord::AtomLabels(vec![Some("Q_e".to_owned()), None]),
                CxRecord::AtomProperties(vec![cosmolkit_cx::CxAtomProperty {
                    atom: 0,
                    name: "atomLabel".to_owned(),
                    value: "ordinary after atomProp".to_owned(),
                }]),
            ],
            0,
        );

        apply_cx_to_query_graph(&mut query, &parsed).expect("source property ordering");

        let atom = query.atom(0).expect("first atom");
        assert_eq!(atom.prop("atomLabel"), Some("ordinary after atomProp"));
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
        assert_eq!(atom.prop("atomLabel"), Some("Pol_p"));
        assert_eq!(atom.prop("dummyLabel"), Some("later property"));
        assert_eq!(atom.predicate(), &first_predicate);
        assert_eq!(query.prop(CX_LABELS_PROCESSED_PROP), Some("1"));

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
            Default::default(),
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
                CxRecord::AtomLabels(vec![Some("left".to_owned()), None]),
            ],
            8,
        );
        let mut query = graph();
        apply_cx_to_query_graph(&mut query, &parsed).unwrap();
        assert_eq!(query.atom(0).unwrap().prop("atomLabel"), Some("left"));
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
        let mut query = graph().with_prop("molFileLinkNodes", "prior");
        apply_cx_to_query_graph(&mut query, &parsed).expect("empty source accumulator");
        assert_eq!(query.prop("molFileLinkNodes"), Some("prior"));
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
            Default::default(),
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
            query.prop("molFileLinkNodes"),
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
        let mut query = graph().with_prop("molFileLinkNodes", "prior");

        assert!(apply_cx_to_query_graph(&mut query, &parsed).is_err());
        assert_eq!(query.prop("molFileLinkNodes"), Some("prior"));
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
                .with_prop("index", "71"),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(0)
                .with_prop("_cxsmilesindex", "0")
                .with_prop("index", "92"),
            SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(1)
                .with_prop("_cxsmilesindex", "1")
                .with_prop("index", "93"),
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
            groups[1].props().get("PARENT").map(String::as_str),
            Some("71")
        );
        assert_eq!(groups[2].parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            groups[2].props().get("PARENT").map(String::as_str),
            Some("71")
        );
    }

    #[test]
    fn cx_progress_hierarchy_skips_unmatched_parent_before_child_range_check() {
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(5)
                .with_prop("_cxsmilesindex", "5")
                .with_prop("index", "71"),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(0)
                .with_prop("_cxsmilesindex", "0")
                .with_prop("index", "92"),
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
                .with_prop("index", "71"),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(0)
                .with_prop("_cxsmilesindex", "0")
                .with_prop("index", "92"),
            SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(1)
                .with_prop("_cxsmilesindex", "1"),
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
            groups[1].props().get("PARENT").map(String::as_str),
            Some("71")
        );
        assert_eq!(groups[2].parent(), None);
        assert_eq!(groups[2].props().get("PARENT"), None);
    }

    #[test]
    fn cx_progress_hierarchy_without_parent_index_uses_cx_parent_id() {
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(17)
                .with_prop("_cxsmilesindex", "17"),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(0)
                .with_prop("_cxsmilesindex", "0")
                .with_prop("index", "4"),
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
        assert_eq!(child.props().get("PARENT").map(String::as_str), Some("17"));
    }
}
