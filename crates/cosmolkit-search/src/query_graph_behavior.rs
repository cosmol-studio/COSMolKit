//! Query-graph post-processing owned by the search implementation.

use cosmolkit_model::{AtomId, BondId, PropertyValue, QueryGraph};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo};

pub(crate) const SMILES_START_PROP: &str = "_SmilesStart";
pub(crate) const CXSMILES_BOND_IDX_PROP: &str = "_cxsmilesBondIdx";
pub(crate) const UNSPECIFIED_ORDER_PROP: &str = "_unspecifiedOrder";

pub(crate) fn cleanup_query_graph_parser_state(graph: &mut QueryGraph) {
    cleanup_query_parser_state(graph);
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

pub(crate) fn set_bond_stereo_from_directions(graph: &mut QueryGraph) {
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
    graph.clear_prop("_needsDetectBondStereo");
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
            bond.bond_mut().set_stereo(stereo);
        }
    }
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
fn cleanup_query_parser_state(graph: &mut QueryGraph) {
    cosmolkit_core::parser_helpers::cleanup_parser_atoms(graph.atoms_mut());
    for bond in graph.bonds_mut() {
        bond.bond_mut().clear_prop("_unspecifiedOrder");
        bond.bond_mut().clear_prop(CXSMILES_BOND_IDX_PROP);
    }
    cosmolkit_core::parser_helpers::cleanup_parser_substance_groups(
        cosmolkit_model::query_substance_groups(graph),
    );
}
