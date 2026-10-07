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

fn neighboring_directed_bond(graph: &QueryGraph, atom: AtomId) -> Option<BondId> {
    // RDKit✔️✔️: for (const auto &bondIdx :
    // RDKit✔️✔️:      boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit✔️✔️:   const Bond *bond = mol[bondIdx];
    // RDKit✔️✔️:   if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit✔️✔️:       hasStereoBondDir(bond)) { return bond; }
    // RDKit✔️✔️: }
    for (_, bond_index) in graph.adjacency().get(atom.index())?.iter().copied() {
        let bond = graph.bonds().get(bond_index)?;
        if bond.bond().order() != BondOrder::Double
            && matches!(
                bond.bond().direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            )
        {
            return Some(bond.id());
        }
    }
    None
}

fn opposite_direction(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        other => other,
    }
}

pub(crate) fn set_bond_stereo_from_directions(
    graph: &mut QueryGraph,
) -> Result<(), crate::SmartsParseError> {
    // RDKit✔️✔️: mol.clearProp("_needsDetectBondStereo");
    // RDKit✔️✔️: if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:     bond->getStereo() != Bond::STEREOANY) {
    // RDKit✔️✔️:   const Bond *directedBondAtBegin =
    // RDKit✔️✔️:       Chirality::getNeighboringDirectedBond(mol, stereoBondBeginAtom);
    // RDKit✔️✔️:   const Bond *directedBondAtEnd =
    // RDKit✔️✔️:       Chirality::getNeighboringDirectedBond(mol, stereoBondEndAtom);
    // RDKit✔️✔️:   if (beginSideBondDirection == endSideBondDirection) {
    // RDKit✔️✔️:     bond->setStereo(Bond::STEREOTRANS);
    // RDKit✔️✔️:   } else { bond->setStereo(Bond::STEREOCIS); }
    // RDKit✔️✔️: }
    graph.clear_prop("_needsDetectBondStereo")?;
    let mut updates = Vec::new();
    for bond in graph.bonds() {
        if bond.bond().order() != BondOrder::Double || bond.bond().stereo() == BondStereo::Any {
            continue;
        }
        let begin = bond.begin();
        let end = bond.end();
        let (Some(begin_id), Some(end_id)) = (
            neighboring_directed_bond(graph, begin),
            neighboring_directed_bond(graph, end),
        ) else {
            continue;
        };
        let (Some(begin_bond), Some(end_bond)) =
            (graph.bond(begin_id.index()), graph.bond(end_id.index()))
        else {
            continue;
        };
        let begin_atom = if begin_bond.begin() == begin {
            begin_bond.end()
        } else {
            begin_bond.begin()
        };
        let end_atom = if end_bond.begin() == end {
            end_bond.end()
        } else {
            end_bond.begin()
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
        updates.push((bond.id(), [begin_atom, end_atom], stereo));
    }
    for (bond_id, stereo_atoms, stereo) in updates {
        if let Some(bond) = graph.bonds_mut().get_mut(bond_id.index()) {
            bond.bond_mut().set_stereo_atoms(Some(stereo_atoms));
            bond.bond_mut().set_stereo(stereo)?;
        }
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
    // RDKit❗❌: void CleanupAfterParsing(RWMol *mol) {
    // RDKit❗❌:   PRECONDITION(mol, "no molecule");
    // RDKit❗❌:   for (auto atom : mol->atoms()) {
    // RDKit❗❌:     atom->clearProp(common_properties::_RingClosures);
    // RDKit❗❌:     atom->clearProp(common_properties::_SmilesStart);
    // RDKit❗❌:     std::string label;
    // RDKit❗❌:     if (atom->getAtomicNum() == 0 &&
    // RDKit❗❌:         atom->getPropIfPresent(common_properties::atomLabel, label)) {
    // RDKit❗❌:       // marvinsketch can output higher labels than _AP1 and _AP2, but they
    // RDKit❗❌:       // aren't part of the MOL file spec so we don't treat them as attachment
    // RDKit❗❌:       // points
    // RDKit❗❌:       if (label == "_AP1") {
    // RDKit❗❌:         atom->setProp(common_properties::_fromAttachPoint, 1);
    // RDKit❗❌:       } else if (label == "_AP2") {
    // RDKit❗❌:         atom->setProp(common_properties::_fromAttachPoint, 2);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto bond : mol->bonds()) {
    // RDKit❗❌:     bond->clearProp(common_properties::_unspecifiedOrder);
    // RDKit❗❌:     bond->clearProp("_cxsmilesBondIdx");
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto sg : RDKit::getSubstanceGroups(*mol)) {
    // RDKit❗❌:     sg.clearProp("_cxsmilesindex");
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!Chirality::getAllowNontetrahedralChirality()) {
    // RDKit❗❌:     bool needWarn = false;
    // RDKit❗❌:     for (auto atom : mol->atoms()) {
    // RDKit❗❌:       if (atom->hasProp(common_properties::_chiralPermutation)) {
    // RDKit❗❌:         needWarn = true;
    // RDKit❗❌:         atom->clearProp(common_properties::_chiralPermutation);
    // RDKit❗❌:       }
    // RDKit❗❌:       if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER) {
    // RDKit❗❌:         needWarn = true;
    // RDKit❗❌:         atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (needWarn) {
    // RDKit❗❌:       BOOST_LOG(rdWarningLog)
    // RDKit❗❌:           << "ignoring non-tetrahedral stereo specification since setAllowNontetrahedralChirality() is false."
    // RDKit❗❌:           << std::endl;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Propagate reached property errors in atom, bond, then copied-group
    // order. Full non-tetrahedral cleanup remains its existing CORE boundary.
    // Known cost: each local SGroup clone includes tree/order key storage.

    cosmolkit_core::parser_helpers::cleanup_parser_atoms(graph.atoms_mut())?;
    for bond in graph.bonds_mut() {
        bond.bond_mut().clear_prop("_unspecifiedOrder")?;
        bond.bond_mut().clear_prop(CXSMILES_BOND_IDX_PROP)?;
    }
    cosmolkit_core::parser_helpers::cleanup_parser_substance_groups(
        cosmolkit_model::query_substance_groups(graph),
    )?;
    Ok(())
}

#[cfg(test)]
mod source_property_failure_tests {
    use super::*;
    use crate::{SmartsParseError, SmartsParseParams, parse_smarts};

    #[test]
    fn directional_stereo_retains_reserved_property_failure_before_bond_updates() {
        for (input, stereo) in [
            ("C/C=C/C", BondStereo::Trans),
            (r"C/C=C\C", BondStereo::Cis),
        ] {
            for marker_present in [false, true] {
                let mut graph = parse_smarts(input, &SmartsParseParams::default()).unwrap();
                graph.bonds_mut()[1]
                    .bond_mut()
                    .set_stereo(BondStereo::None)
                    .unwrap();
                graph.bonds_mut()[1].bond_mut().set_stereo_atoms(None);
                if marker_present {
                    graph.set_prop("_needsDetectBondStereo", true).unwrap();
                }
                graph
                    .set_prop("__computedProps", PropertyValue::Int(7))
                    .unwrap();
                let before = graph.clone();
                let error = set_bond_stereo_from_directions(&mut graph).unwrap_err();
                assert!(matches!(
                    &error,
                    SmartsParseError::MoleculeProperty(
                        cosmolkit_model::MoleculePropertyError::ComputedListKind(_)
                    )
                ));
                assert!(
                    std::error::Error::source(&error)
                        .unwrap()
                        .downcast_ref::<cosmolkit_model::MoleculePropertyError>()
                        .is_some()
                );
                assert_eq!(
                    graph, before,
                    "source property failure precedes every stereo update"
                );

                let names = if marker_present {
                    vec![cosmolkit_model::PropertyText::from(
                        "_needsDetectBondStereo",
                    )]
                } else {
                    Vec::new()
                };
                graph
                    .set_prop("__computedProps", PropertyValue::StringVector(names))
                    .unwrap();
                set_bond_stereo_from_directions(&mut graph).unwrap();
                assert_eq!(graph.prop("_needsDetectBondStereo"), None);
                assert_eq!(
                    graph.prop("__computedProps"),
                    Some(&PropertyValue::StringVector(Vec::new()))
                );
                assert_eq!(graph.bonds()[1].bond().stereo(), stereo);
                assert_eq!(
                    graph.bonds()[1].bond().stereo_atoms(),
                    Some([AtomId::new(0), AtomId::new(3)])
                );
            }
        }
    }
}
