//! Query-graph post-processing owned by the search implementation.

use cosmolkit_model::{AtomId, BondId, PropertyValue, QueryGraph};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo};

pub(crate) const SMILES_START_PROP: &str = "_SmilesStart";
pub(crate) const CXSMILES_BOND_IDX_PROP: &str = "_cxsmilesBondIdx";
pub(crate) const UNSPECIFIED_ORDER_PROP: &str = "_unspecifiedOrder";

/// Complete delayed parser cleanup after reaction-wide CX lowering.
#[doc(hidden)]
pub fn cleanup_query_graph_parser_state(
    graph: &mut QueryGraph,
) -> Result<(), crate::SmartsParseError> {
    cleanup_query_parser_state(graph)
}

fn neighboring_directed_bond(
    graph: &QueryGraph,
    atom: AtomId,
) -> Result<Option<BondId>, crate::SmartsParseError> {
    // RDKit✔️✔️: const Bond *getNeighboringDirectedBond(const ROMol &mol, const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "no atom");
    // RDKit✔️✔️:   for (const auto &bondIdx :
    // RDKit✔️✔️:        boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit✔️✔️:     const Bond *bond = mol[bondIdx];
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit✔️✔️:         hasStereoBondDir(bond)) {
    // RDKit✔️✔️:       return bond;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // Canonical QueryGraph adjacency is physical bond-row order, validated
    // against that order by its constructor. No neighbor sorting or copying.
    // A missing row is an invariant failure, never a source nullptr fallback.
    let neighbors = graph
        .adjacency()
        .get(atom.index())
        .ok_or(cosmolkit_model::QueryGraphError::AdjacencyMismatch)?;
    for &(_, bond_index) in neighbors {
        let bond = graph
            .bonds()
            .get(bond_index)
            .ok_or(cosmolkit_model::QueryGraphError::AdjacencyMismatch)?;
        if bond.bond().order() != BondOrder::Double
            && matches!(
                bond.bond().direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            )
        {
            return Ok(Some(bond.id()));
        }
    }
    Ok(None)
}

fn opposite_direction(direction: BondDirection) -> BondDirection {
    // RDKit✔️✔️: Bond::BondDir getOppositeBondDir(Bond::BondDir dir) {
    // RDKit✔️✔️:   PRECONDITION(dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT,
    // RDKit✔️✔️:                "bad bond direction");
    // RDKit✔️✔️:   switch (dir) {
    // RDKit✔️✔️:     case Bond::ENDDOWNRIGHT:
    // RDKit✔️✔️:       return Bond::ENDUPRIGHT;
    // RDKit✔️✔️:     case Bond::ENDUPRIGHT:
    // RDKit✔️✔️:       return Bond::ENDDOWNRIGHT;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       return Bond::NONE;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Source callers pass only the two directed-bond states. Enum dispatch
    // preserves those inversions in O(1), without allocation or extra lookup.
    match direction {
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        other => other,
    }
}

pub(crate) fn set_bond_stereo_from_directions(
    graph: &mut QueryGraph,
) -> Result<(), crate::SmartsParseError> {
    // RDKit✔️✔️: void setBondStereoFromDirections(ROMol &mol) {
    // RDKit✔️✔️:   mol.clearProp("_needsDetectBondStereo");
    // RDKit✔️✔️:   for (Bond *bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:         bond->getStereo() != Bond::STEREOANY) {
    // RDKit✔️✔️:       const Atom *stereoBondBeginAtom = bond->getBeginAtom();
    // RDKit✔️✔️:       const Atom *stereoBondEndAtom = bond->getEndAtom();
    // RDKit✔️✔️:
    // RDKit✔️✔️:       const Bond *directedBondAtBegin =
    // RDKit✔️✔️:           Chirality::getNeighboringDirectedBond(mol, stereoBondBeginAtom);
    // RDKit✔️✔️:       const Bond *directedBondAtEnd =
    // RDKit✔️✔️:           Chirality::getNeighboringDirectedBond(mol, stereoBondEndAtom);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (directedBondAtBegin != nullptr && directedBondAtEnd != nullptr) {
    // RDKit✔️✔️:         unsigned beginSideStereoAtom =
    // RDKit✔️✔️:             directedBondAtBegin->getOtherAtomIdx(stereoBondBeginAtom->getIdx());
    // RDKit✔️✔️:         unsigned endSideStereoAtom =
    // RDKit✔️✔️:             directedBondAtEnd->getOtherAtomIdx(stereoBondEndAtom->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:         bond->setStereoAtoms(beginSideStereoAtom, endSideStereoAtom);
    // RDKit✔️✔️:
    // RDKit✔️✔️:         auto beginSideBondDirection = directedBondAtBegin->getBondDir();
    // RDKit✔️✔️:         if (directedBondAtBegin->getBeginAtom() == stereoBondBeginAtom) {
    // RDKit✔️✔️:           beginSideBondDirection = getOppositeBondDir(beginSideBondDirection);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         auto endSideBondDirection = directedBondAtEnd->getBondDir();
    // RDKit✔️✔️:         if (directedBondAtEnd->getEndAtom() == stereoBondEndAtom) {
    // RDKit✔️✔️:           endSideBondDirection = getOppositeBondDir(endSideBondDirection);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if (beginSideBondDirection == endSideBondDirection) {
    // RDKit✔️✔️:           bond->setStereo(Bond::STEREOTRANS);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           bond->setStereo(Bond::STEREOCIS);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // Same source clear-first and physical row mutation order. Each row's
    // source-selected neighboring directions are borrowed before mutation;
    // only a two-ID scalar result is retained, not a graph/update buffer.
    // Source absence of a directed neighbor is the only skip branch; every
    // reached property/invariant/setStereo failure propagates structurally.
    graph.clear_prop("_needsDetectBondStereo")?;
    for index in 0..graph.num_bonds() {
        let bond = &graph.bonds()[index];
        if bond.bond().order() != BondOrder::Double || bond.bond().stereo() == BondStereo::Any {
            continue;
        }
        let begin = bond.begin();
        let end = bond.end();
        let begin_id = neighboring_directed_bond(graph, begin)?;
        let end_id = neighboring_directed_bond(graph, end)?;
        let (Some(begin_id), Some(end_id)) = (begin_id, end_id) else {
            continue;
        };
        let begin_bond = graph
            .bond(begin_id.index())
            .ok_or(cosmolkit_model::QueryGraphError::AdjacencyMismatch)?;
        let end_bond = graph
            .bond(end_id.index())
            .ok_or(cosmolkit_model::QueryGraphError::AdjacencyMismatch)?;
        let begin_atom = if begin_bond.begin() == begin {
            begin_bond.end()
        } else if begin_bond.end() == begin {
            begin_bond.begin()
        } else {
            return Err(cosmolkit_model::QueryGraphError::AdjacencyMismatch.into());
        };
        let end_atom = if end_bond.begin() == end {
            end_bond.end()
        } else if end_bond.end() == end {
            end_bond.begin()
        } else {
            return Err(cosmolkit_model::QueryGraphError::AdjacencyMismatch.into());
        };
        let begin_direction = if begin_bond.begin() == begin {
            opposite_direction(begin_bond.bond().direction())
        } else {
            begin_bond.bond().direction()
        };
        let end_direction = if end_bond.end() == end {
            opposite_direction(end_bond.bond().direction())
        } else {
            end_bond.bond().direction()
        };
        let stereo = if begin_direction == end_direction {
            BondStereo::Trans
        } else {
            BondStereo::Cis
        };
        let carrier = graph.bonds_mut()[index].bond_mut();
        carrier.set_stereo_atoms(Some([begin_atom, end_atom]));
        carrier.set_stereo(stereo)?;
    }
    Ok(())
}

pub(crate) fn check_chiral_permutation(tag: cosmolkit_types::ChiralTag, permutation: i32) -> bool {
    // RDKit✔️✔️: if (chiralTag > RDKit::Atom::ChiralType::CHI_OTHER &&
    // RDKit✔️✔️:     permutationLimits.find(chiralTag) != permutationLimits.end() &&
    // RDKit✔️✔️:     (permutation < 0 || permutation > permutationLimits.at(chiralTag))) {
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return true;
    let limit = match tag {
        cosmolkit_types::ChiralTag::Tetrahedral | cosmolkit_types::ChiralTag::Allene => Some(2),
        cosmolkit_types::ChiralTag::SquarePlanar => Some(3),
        cosmolkit_types::ChiralTag::TrigonalBipyramidal => Some(20),
        cosmolkit_types::ChiralTag::Octahedral => Some(30),
        _ => None,
    };
    limit.is_none_or(|limit| permutation >= 0 && permutation <= limit)
}

/// Apply shared source parser chirality finalization to detached query carriers.
/// Syntax flags and source-order ring bond IDs must already be materialized.
#[doc(hidden)]
pub(crate) fn finalize_query_parser_chirality(
    graph: &mut QueryGraph,
) -> Result<(), cosmolkit_core::parser_helpers::ParserCarrierError> {
    let rings = graph
        .atoms()
        .iter()
        .map(|atom| match atom.prop("_RingClosures") {
            None => Ok(Vec::new()),
            Some(PropertyValue::IntVector(ids)) => {
                ids.iter()
                    .map(|&id| {
                        usize::try_from(id).map(BondId::new).map_err(|_| {
                            cosmolkit_core::parser_helpers::ParserCarrierError::Model(
                "query ring closure remained unresolved during chirality adjustment".into())
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
            }
            Some(_) => Err(cosmolkit_core::parser_helpers::ParserCarrierError::Model(
                "query ring closure property has the wrong type".into(),
            )),
        })
        .collect::<Result<Vec<_>, _>>()?;
    let starts = graph
        .atoms()
        .iter()
        .map(|atom| atom.prop("_SmilesStart").is_some())
        .collect::<Vec<_>>();
    let assignments = cosmolkit_core::parser_helpers::parser_chirality_assignments(
        graph.atoms(),
        graph.bonds(),
        |index| {
            graph.adjacency()[index]
                .iter()
                .map(|&(atom, bond)| (atom, BondId::new(bond)))
        },
        &rings,
        &starts,
    )?;
    for (atom, (tag, permutation)) in graph.atoms_mut().iter_mut().zip(assignments) {
        atom.set_chiral_tag(tag);
        atom.set_chiral_permutation(permutation);
    }
    Ok(())
}

/// Cleanup source parser state on the canonical detached query value.
#[doc(hidden)]
fn cleanup_query_parser_state(graph: &mut QueryGraph) -> Result<(), crate::SmartsParseError> {
    // RDKit✔️❌: void CleanupAfterParsing(RWMol *mol) {
    // RDKit✔️❌:   PRECONDITION(mol, "no molecule");
    // RDKit✔️❌:   for (auto atom : mol->atoms()) {
    // RDKit✔️❌:     atom->clearProp(common_properties::_RingClosures);
    // RDKit✔️❌:     atom->clearProp(common_properties::_SmilesStart);
    // RDKit✔️❌:     std::string label;
    // RDKit✔️❌:     if (atom->getAtomicNum() == 0 &&
    // RDKit✔️❌:         atom->getPropIfPresent(common_properties::atomLabel, label)) {
    // RDKit✔️❌:       // marvinsketch can output higher labels than _AP1 and _AP2, but they
    // RDKit✔️❌:       // aren't part of the MOL file spec so we don't treat them as attachment
    // RDKit✔️❌:       // points
    // RDKit✔️❌:       if (label == "_AP1") {
    // RDKit✔️❌:         atom->setProp(common_properties::_fromAttachPoint, 1);
    // RDKit✔️❌:       } else if (label == "_AP2") {
    // RDKit✔️❌:         atom->setProp(common_properties::_fromAttachPoint, 2);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto bond : mol->bonds()) {
    // RDKit✔️❌:     bond->clearProp(common_properties::_unspecifiedOrder);
    // RDKit✔️❌:     bond->clearProp("_cxsmilesBondIdx");
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto sg : RDKit::getSubstanceGroups(*mol)) {
    // RDKit✔️❌:     sg.clearProp("_cxsmilesindex");
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!Chirality::getAllowNontetrahedralChirality()) {
    // RDKit✔️❌:     bool needWarn = false;
    // RDKit✔️❌:     for (auto atom : mol->atoms()) {
    // RDKit✔️❌:       if (atom->hasProp(common_properties::_chiralPermutation)) {
    // RDKit✔️❌:         needWarn = true;
    // RDKit✔️❌:         atom->clearProp(common_properties::_chiralPermutation);
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER) {
    // RDKit✔️❌:         needWarn = true;
    // RDKit✔️❌:         atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (needWarn) {
    // RDKit✔️❌:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:           << "ignoring non-tetrahedral stereo specification since setAllowNontetrahedralChirality() is false."
    // RDKit✔️❌:           << std::endl;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // Propagate errors in source atom, bond, copied-group, then final
    // non-tetrahedral pass order. Shared CORE owns both atom algorithms.
    // Known cost: each local SGroup clone includes tree/order key storage.

    cosmolkit_core::parser_helpers::cleanup_parser_atoms(graph.atoms_mut())?;
    for bond in graph.bonds_mut() {
        bond.bond_mut().clear_prop("_unspecifiedOrder")?;
        bond.bond_mut().clear_prop(CXSMILES_BOND_IDX_PROP)?;
    }
    cosmolkit_core::parser_helpers::cleanup_parser_substance_groups(
        cosmolkit_model::query_substance_groups(graph),
    )?;
    cosmolkit_core::parser_helpers::cleanup_parser_nontetrahedral_atoms(graph.atoms_mut())?;
    Ok(())
}
