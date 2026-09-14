use std::collections::BTreeMap;

use cosmolkit_core::{ValenceAssignment, ValenceModel, fast_find_rings_from_parts};
use cosmolkit_model::{Atom, AtomId, Bond, BondId, TopologyBlock};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag};

use crate::{SmilesParseError, SmilesRecord, canonical_rank, stereo};

mod direction;

const MAX_NATOMS: i64 = 5000;
const MAX_BONDTYPE: i64 = 32;
const MAX_CYCLES: usize = 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AtomColor {
    White,
    Grey,
    Black,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MolStackElem {
    Atom(usize),
    Bond { bond: BondId, atom_to_left: usize },
    Ring(usize),
    BranchOpen,
    BranchClose,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Possible {
    rank: i64,
    atom: usize,
    bond: BondId,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct ChiralAdjustment {
    chiral_tag_override: Option<ChiralTag>,
    invert_tetrahedral: bool,
    nontetrahedral_permutation: Option<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SmilesWriteOutput {
    pub(crate) text: String,
    pub(crate) atom_order: Vec<AtomId>,
    pub(crate) bond_order: Vec<BondId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FragmentWriteOutput {
    text: String,
    atom_order: Vec<AtomId>,
    bond_order: Vec<BondId>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmilesWriteParams {
    pub canonical: bool,
}

impl Default for SmilesWriteParams {
    fn default() -> Self {
        Self { canonical: true }
    }
}

/// Writes canonical SMILES from detached values using RDKit-compatible atom
/// ranking and traversal.
pub fn write_smiles(record: &SmilesRecord) -> Result<String, SmilesParseError> {
    write_smiles_with_params(record, &SmilesWriteParams::default())
}

/// Writes SMILES from detached values with explicit canonicalization policy.
pub fn write_smiles_with_params(
    record: &SmilesRecord,
    params: &SmilesWriteParams,
) -> Result<String, SmilesParseError> {
    write_smiles_output(record, params, false).map(|output| output.text)
}

pub(crate) fn write_smiles_for_cx(
    record: &SmilesRecord,
    params: &SmilesWriteParams,
) -> Result<SmilesWriteOutput, SmilesParseError> {
    write_smiles_output(record, params, true)
}

fn write_smiles_output(
    record: &SmilesRecord,
    params: &SmilesWriteParams,
    doing_cx_smiles: bool,
) -> Result<SmilesWriteOutput, SmilesParseError> {
    record
        .topology
        .validate()
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;
    reject_unmodeled_stereochemical_writing(&record.topology, doing_cx_smiles)?;
    if record.topology.atoms.is_empty() {
        return Ok(SmilesWriteOutput {
            text: String::new(),
            atom_order: Vec::new(),
            bond_order: Vec::new(),
        });
    }

    // RDKit's public writer performs stereo preparation and direction
    // canonicalization on a private molecule copy. Keep the detached record
    // immutable and make those same temporary changes on this working block.
    let mut topology = record.topology.clone();
    let source_valence = if doing_cx_smiles {
        Some(
            cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                ValenceModel::RdkitLike,
                false,
            )
            .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?,
        )
    } else {
        None
    };
    let mut dative_donors = Vec::new();
    let mut hydrogen_bond_atoms = Vec::new();
    // RDKit✔️✔️: if (!doingCXSmiles || !includeStereoGroups) {
    // RDKit✔️✔️:   std::vector<StereoGroup> noStereoGroups;
    // RDKit✔️✔️:   tmol->setStereoGroups(noStereoGroups);
    // RDKit✔️✔️: }
    if !doing_cx_smiles {
        topology.stereo_groups.clear();
    } else {
        // RDKit✔️✔️: if (doingCXSmiles || !params.includeDativeBonds) {
        // RDKit✔️✔️:   for (auto bond : tmol->bonds()) {
        // RDKit✔️✔️:     if (bond->getBondType() == Bond::DATIVE) {
        // RDKit✔️✔️:       bond->setBondType(Bond::SINGLE);
        // RDKit✔️✔️:       bond->getBeginAtom()->calcExplicitValence(false);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: if (doingCXSmiles) {
        // RDKit✔️✔️:   for (auto bond : tmol->bonds()) {
        // RDKit✔️✔️:     if (bond->getBondType() == Bond::HYDROGEN) {
        // RDKit✔️✔️:       bond->setBondType(Bond::SINGLE);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        for bond in &mut topology.bonds {
            if bond.order() == BondOrder::Dative {
                dative_donors.push(bond.begin());
                bond.set_order(BondOrder::Single);
            } else if bond.order() == BondOrder::Hydrogen {
                hydrogen_bond_atoms.push(bond.begin());
                hydrogen_bond_atoms.push(bond.end());
                bond.set_order(BondOrder::Single);
            }
            if !matches!(
                bond.direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            ) {
                bond.set_direction(BondDirection::None);
            }
            if bond.stereo() == BondStereo::Any {
                bond.set_stereo(BondStereo::None);
            }
        }
        topology.adjacency =
            cosmolkit_model::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
    }
    // RDKit✔️✔️:     for (auto atom : tmol->atoms()) {
    // RDKit✔️✔️:       atom->updatePropertyCache(false);
    // RDKit✔️✔️:     }
    let mut valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &topology,
        ValenceModel::RdkitLike,
        false,
    )
    .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?;
    if let Some(source_valence) = &source_valence {
        // RDKit✔️✔️: bond->setBondType(Bond::SINGLE);
        // RDKit✔️✔️: // update the explicit valence of the begin atom since the implicit
        // RDKit✔️✔️: // valence will no longer be properly perceived
        // RDKit✔️✔️: bond->getBeginAtom()->calcExplicitValence(false);
        // The source deliberately retains the donor's already perceived
        // implicit valence while recalculating only its explicit valence.
        for atom in &dative_donors {
            valence.implicit_hydrogens[atom.index()] =
                source_valence.implicit_hydrogens[atom.index()];
        }
        // RDKit converts hydrogen bonds to single bonds without recalculating
        // either endpoint's cached valence.
        for atom in &hydrogen_bond_atoms {
            valence.explicit_valence[atom.index()] = source_valence.explicit_valence[atom.index()];
            valence.implicit_hydrogens[atom.index()] =
                source_valence.implicit_hydrogens[atom.index()];
        }
    }
    let rings =
        fast_find_rings_from_parts(topology.atoms.len(), &topology.bonds, &topology.adjacency)
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
    if topology.bonds.iter().any(|bond| {
        matches!(
            bond.direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        )
    }) {
        let ranks = cosmolkit_core::assign_atom_cip_ranks(&topology, &valence)
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        topology = cosmolkit_core::assign_directional_double_bond_stereo(topology, &ranks, &rings)
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?
            .topology;
    }
    let mut ranking_topology = topology.clone();
    if doing_cx_smiles {
        for atom in dative_donors
            .iter()
            .chain(hydrogen_bond_atoms.iter())
            .copied()
        {
            let total_hydrogens =
                usize::from(ranking_topology.atoms[atom.index()].explicit_hydrogens())
                    + usize::try_from(valence.implicit_hydrogens[atom.index()].max(0))
                        .unwrap_or(usize::MAX);
            ranking_topology.atoms[atom.index()]
                .set_explicit_hydrogens(u8::try_from(total_hydrogens).unwrap_or(u8::MAX));
            ranking_topology.atoms[atom.index()].set_no_implicit(true);
        }
    }
    let ring_bonds = find_ring_bonds(&topology);
    let mut colors = vec![AtomColor::White; topology.atoms.len()];
    let mut fragments = Vec::new();
    let components = connected_components(&topology);

    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles non-canonical fragments
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       std::iota(ranks.begin(), ranks.end(), 0);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     subSmi = SmilesWrite::FragmentSmilesConstruct(
    // RDKit✔️✔️:         *tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (params.canonical) {
    // RDKit✔️✔️:     std::sort(tmp.begin(), tmp.end());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (unsigned i = 0; i < vfragsmi.size(); ++i) {
    // RDKit✔️✔️:       result += vfragsmi[i];
    // RDKit✔️✔️:       if (i < vfragsmi.size() - 1) {
    // RDKit✔️✔️:         result += ".";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles non-canonical fragments
    for component in components {
        let ranks = if params.canonical {
            canonical_rank::rank_component_atoms(&ranking_topology, &component)
                .map_err(|error| SmilesParseError::CanonicalRank(error.to_string()))?
                .into_iter()
                .map(|rank| rank as i64)
                .collect::<Vec<_>>()
        } else {
            (0..topology.atoms.len())
                .map(|index| index as i64)
                .collect::<Vec<_>>()
        };
        let start = if params.canonical {
            component
                .iter()
                .copied()
                .min_by_key(|atom| ranks[*atom])
                .expect("connected component is nonempty")
        } else {
            component[0]
        };
        let mut cycle_colors = colors.clone();
        let mut ring_closures = vec![Vec::new(); topology.atoms.len()];
        dfs_find_cycles(
            &topology,
            start,
            None,
            &mut cycle_colors,
            &ranks,
            &ring_bonds,
            &mut ring_closures,
        );
        let mut stack = Vec::with_capacity(topology.atoms.len() + topology.bonds.len());
        let mut ring_ids = vec![None; topology.bonds.len()];
        let mut traversal_ring_closure_bonds = vec![false; topology.bonds.len()];
        let mut available_ring_ids = vec![true; MAX_CYCLES];
        let mut atom_traversal_bond_order = vec![Vec::new(); topology.atoms.len()];
        dfs_build_stack(
            &topology,
            start,
            None,
            &mut colors,
            &ranks,
            &ring_bonds,
            &ring_closures,
            &mut ring_ids,
            &mut available_ring_ids,
            &mut stack,
            &mut atom_traversal_bond_order,
            &mut traversal_ring_closure_bonds,
        )?;
        let chiral_adjustments = compute_chiral_adjustments(
            &topology,
            &valence,
            start,
            &ring_closures,
            &atom_traversal_bond_order,
            &stack,
            doing_cx_smiles && params.canonical,
        )?;
        direction::canonicalize_double_bond_directions_for_writer(
            &mut topology,
            &stack,
            &traversal_ring_closure_bonds,
        )?;
        let text = write_mol_stack(&topology, &valence, &stack, &chiral_adjustments)?;
        let atom_order = stack
            .iter()
            .filter_map(|element| match element {
                MolStackElem::Atom(atom) => Some(AtomId::new(*atom)),
                _ => None,
            })
            .collect();
        let bond_order = stack
            .iter()
            .filter_map(|element| match element {
                MolStackElem::Bond { bond, .. } => Some(*bond),
                _ => None,
            })
            .collect();
        fragments.push(FragmentWriteOutput {
            text,
            atom_order,
            bond_order,
        });
    }
    if params.canonical {
        fragments.sort_by(|left, right| left.text.cmp(&right.text));
    }
    let text = fragments
        .iter()
        .map(|fragment| fragment.text.as_str())
        .collect::<Vec<_>>()
        .join(".");
    let atom_order = fragments
        .iter()
        .flat_map(|fragment| fragment.atom_order.iter().copied())
        .collect();
    let bond_order = fragments
        .iter()
        .flat_map(|fragment| fragment.bond_order.iter().copied())
        .collect();
    Ok(SmilesWriteOutput {
        text,
        atom_order,
        bond_order,
    })
}

fn connected_components(topology: &TopologyBlock) -> Vec<Vec<usize>> {
    let mut seen = vec![false; topology.atoms.len()];
    let mut components = Vec::new();
    for start in 0..topology.atoms.len() {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut component = Vec::new();
        let mut pending = vec![start];
        while let Some(atom) = pending.pop() {
            component.push(atom);
            for neighbor in topology.adjacency.neighbors_of(atom) {
                if !seen[neighbor.atom_index] {
                    seen[neighbor.atom_index] = true;
                    pending.push(neighbor.atom_index);
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }
    components
}

fn reject_unmodeled_stereochemical_writing(
    topology: &TopologyBlock,
    doing_cx_smiles: bool,
) -> Result<(), SmilesParseError> {
    if topology.atoms.iter().any(|atom| {
        !matches!(
            atom.chiral_tag(),
            ChiralTag::Unspecified
                | ChiralTag::TetrahedralCw
                | ChiralTag::TetrahedralCcw
                | ChiralTag::SquarePlanar
                | ChiralTag::TrigonalBipyramidal
                | ChiralTag::Octahedral
        ) || atom.unknown_stereo()
    }) {
        return Err(SmilesParseError::UnsupportedWriter(
            "allene, generic, or unknown atom stereochemistry is not modeled by the detached writer",
        ));
    }
    if topology.bonds.iter().any(|bond| {
        (!doing_cx_smiles
            && (!matches!(
                bond.direction(),
                BondDirection::None | BondDirection::EndDownRight | BondDirection::EndUpRight
            ) || !matches!(
                bond.stereo(),
                BondStereo::None
                    | BondStereo::E
                    | BondStereo::Z
                    | BondStereo::Cis
                    | BondStereo::Trans
            ) || bond.unknown_stereo()))
            || (doing_cx_smiles
                && (!matches!(
                    bond.direction(),
                    BondDirection::None
                        | BondDirection::BeginWedge
                        | BondDirection::BeginDash
                        | BondDirection::EndDownRight
                        | BondDirection::EndUpRight
                        | BondDirection::EitherDouble
                        | BondDirection::Unknown
                ) || !matches!(
                    bond.stereo(),
                    BondStereo::None
                        | BondStereo::Any
                        | BondStereo::E
                        | BondStereo::Z
                        | BondStereo::Cis
                        | BondStereo::Trans
                )))
    }) {
        return Err(SmilesParseError::UnsupportedWriter(
            "unknown, atropisomeric, or non-directional bond stereochemistry is not modeled by the detached writer",
        ));
    }
    Ok(())
}

fn dfs_find_cycles(
    topology: &TopologyBlock,
    atom: usize,
    incoming_bond: Option<BondId>,
    colors: &mut [AtomColor],
    ranks: &[i64],
    ring_bonds: &[bool],
    atom_ring_closures: &mut [Vec<BondId>],
) {
    // BEGIN RDKIT CPP FUNCTION Canon::dfsFindCycles
    // RDKit✔️✔️:   colors[atomIdx] = GREY_NODE;
    // RDKit✔️✔️:   std::vector<PossibleType> possibles;
    // RDKit✔️✔️:   for (auto &possible : possibles) {
    // RDKit✔️✔️:     int possibleIdx = std::get<1>(possible);
    // RDKit✔️✔️:     Bond *bond = std::get<2>(possible);
    // RDKit✔️✔️:     switch (colors[possibleIdx]) {
    // RDKit✔️✔️:       case WHITE_NODE:
    // RDKit✔️✔️:         dfsFindCycles(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit✔️✔️:                       atomRingClosures, bondsInPlay, bondSymbols, doRandom);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case GREY_NODE:
    // RDKit✔️✔️:         atomRingClosures[possibleIdx].push_back(bond->getIdx());
    // RDKit✔️✔️:         atomRingClosures[atomIdx].push_back(bond->getIdx());
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   colors[atomIdx] = BLACK_NODE;
    // END RDKIT CPP FUNCTION Canon::dfsFindCycles
    colors[atom] = AtomColor::Grey;
    let mut possibles = topology
        .adjacency
        .neighbors_of(atom)
        .iter()
        .copied()
        .filter(|neighbor| Some(neighbor.bond) != incoming_bond)
        .map(|neighbor| {
            let bond = &topology.bonds[neighbor.bond.index()];
            let mut rank = ranks[neighbor.atom_index];
            if colors[neighbor.atom_index] == AtomColor::Grey {
                rank -= (MAX_BONDTYPE + 1) * MAX_NATOMS * MAX_NATOMS;
                rank += (MAX_BONDTYPE - rdkit_bond_type_code(bond.order())) * MAX_NATOMS;
            } else if ring_bonds[neighbor.bond.index()] {
                rank +=
                    (MAX_BONDTYPE - rdkit_bond_type_code(bond.order())) * MAX_NATOMS * MAX_NATOMS;
            }
            Possible {
                rank,
                atom: neighbor.atom_index,
                bond: neighbor.bond,
            }
        })
        .collect::<Vec<_>>();
    possibles.sort_by_key(|possible| possible.rank);

    for possible in possibles {
        match colors[possible.atom] {
            AtomColor::White => dfs_find_cycles(
                topology,
                possible.atom,
                Some(possible.bond),
                colors,
                ranks,
                ring_bonds,
                atom_ring_closures,
            ),
            AtomColor::Grey => {
                atom_ring_closures[possible.atom].push(possible.bond);
                atom_ring_closures[atom].push(possible.bond);
            }
            AtomColor::Black => {}
        }
    }
    colors[atom] = AtomColor::Black;
}

#[allow(clippy::too_many_arguments)]
fn dfs_build_stack(
    topology: &TopologyBlock,
    atom: usize,
    incoming_bond: Option<BondId>,
    colors: &mut [AtomColor],
    ranks: &[i64],
    ring_bonds: &[bool],
    atom_ring_closures: &[Vec<BondId>],
    ring_ids: &mut [Option<usize>],
    available_ring_ids: &mut [bool],
    stack: &mut Vec<MolStackElem>,
    atom_traversal_bond_order: &mut [Vec<BondId>],
    traversal_ring_closure_bonds: &mut [bool],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Canon::dfsBuildStack ring closures and branches
    // RDKit✔️✔️:   molStack.push_back(MolStackElem(atom));
    // RDKit✔️✔️:   colors[atomIdx] = GREY_NODE;
    // RDKit✔️✔️:   if (!atomRingClosures[atomIdx].empty()) {
    // RDKit✔️✔️:     std::vector<unsigned int> ringsClosed;
    // RDKit✔️✔️:     for (auto bIdx : atomRingClosures[atomIdx]) {
    // RDKit✔️✔️:       Bond *bond = mol.getBondWithIdx(bIdx);
    // RDKit✔️✔️:       if (bond->getPropIfPresent(common_properties::_TraversalRingClosureBond,
    // RDKit✔️✔️:                                  ringIdx)) {
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(ringIdx));
    // RDKit✔️✔️:         ringsClosed.push_back(ringIdx - 1);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         auto lowestRingIdx = cyclesAvailable.find_first();
    // RDKit✔️✔️:         cyclesAvailable.set(lowestRingIdx, false);
    // RDKit✔️✔️:         ++lowestRingIdx;
    // RDKit✔️✔️:         bond->setProp(common_properties::_TraversalRingClosureBond,
    // RDKit✔️✔️:                       static_cast<unsigned int>(lowestRingIdx));
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(lowestRingIdx));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (auto ringIdx : ringsClosed) {
    // RDKit✔️✔️:       cyclesAvailable.set(ringIdx);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto possiblesIt = possibles.begin(); possiblesIt != possibles.end();
    // RDKit✔️✔️:        ++possiblesIt) {
    // RDKit✔️✔️:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit✔️✔️:       molStack.push_back(
    // RDKit✔️✔️:           MolStackElem("(", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit✔️✔️:     dfsBuildStack(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit✔️✔️:                   cyclesAvailable, molStack, atomRingClosures,
    // RDKit✔️✔️:                   atomTraversalBondOrder, bondsInPlay, bondSymbols, doRandom);
    // RDKit✔️✔️:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit✔️✔️:       molStack.push_back(
    // RDKit✔️✔️:           MolStackElem(")", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION Canon::dfsBuildStack ring closures and branches
    stack.push(MolStackElem::Atom(atom));
    colors[atom] = AtomColor::Grey;
    let mut traversal_order = Vec::new();
    if let Some(incoming_bond) = incoming_bond {
        traversal_order.push(incoming_bond);
    }
    let mut seen_from_here = vec![false; topology.atoms.len()];
    seen_from_here[atom] = true;
    let mut closed = Vec::new();
    for &bond_id in &atom_ring_closures[atom] {
        traversal_order.push(bond_id);
        let bond = &topology.bonds[bond_id.index()];
        seen_from_here[other_atom(bond, atom)?] = true;
        if let Some(ring_id) = ring_ids[bond_id.index()] {
            stack.push(MolStackElem::Bond {
                bond: bond_id,
                atom_to_left: atom,
            });
            stack.push(MolStackElem::Ring(ring_id));
            closed.push(ring_id - 1);
        } else {
            let Some(slot) = available_ring_ids.iter().position(|available| *available) else {
                return Err(SmilesParseError::UnsupportedWriter(
                    "more than 1024 traversal ring closures are open",
                ));
            };
            available_ring_ids[slot] = false;
            let ring_id = slot + 1;
            ring_ids[bond_id.index()] = Some(ring_id);
            traversal_ring_closure_bonds[bond_id.index()] = true;
            stack.push(MolStackElem::Ring(ring_id));
        }
    }
    for slot in closed {
        available_ring_ids[slot] = true;
    }

    let mut possibles = topology
        .adjacency
        .neighbors_of(atom)
        .iter()
        .copied()
        .filter(|neighbor| Some(neighbor.bond) != incoming_bond)
        .filter(|neighbor| {
            colors[neighbor.atom_index] == AtomColor::White && !seen_from_here[neighbor.atom_index]
        })
        .map(|neighbor| {
            let bond = &topology.bonds[neighbor.bond.index()];
            let mut rank = ranks[neighbor.atom_index];
            if ring_bonds[neighbor.bond.index()] {
                rank +=
                    (MAX_BONDTYPE - rdkit_bond_type_code(bond.order())) * MAX_NATOMS * MAX_NATOMS;
            }
            Possible {
                rank,
                atom: neighbor.atom_index,
                bond: neighbor.bond,
            }
        })
        .collect::<Vec<_>>();
    possibles.sort_by_key(|possible| possible.rank);

    for (position, possible) in possibles.iter().copied().enumerate() {
        if colors[possible.atom] != AtomColor::White {
            continue;
        }
        let is_branch = position + 1 != possibles.len();
        if is_branch {
            stack.push(MolStackElem::BranchOpen);
        }
        stack.push(MolStackElem::Bond {
            bond: possible.bond,
            atom_to_left: atom,
        });
        traversal_order.push(possible.bond);
        dfs_build_stack(
            topology,
            possible.atom,
            Some(possible.bond),
            colors,
            ranks,
            ring_bonds,
            atom_ring_closures,
            ring_ids,
            available_ring_ids,
            stack,
            atom_traversal_bond_order,
            traversal_ring_closure_bonds,
        )?;
        if is_branch {
            stack.push(MolStackElem::BranchClose);
        }
    }
    // RDKit✔️✔️:   atomTraversalBondOrder[atom->getIdx()] = travList;
    atom_traversal_bond_order[atom] = traversal_order;
    colors[atom] = AtomColor::Black;
    Ok(())
}

fn compute_chiral_adjustments(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    start_atom: usize,
    atom_ring_closures: &[Vec<BondId>],
    atom_traversal_bond_order: &[Vec<BondId>],
    stack: &[MolStackElem],
    include_stereo_groups: bool,
) -> Result<Vec<ChiralAdjustment>, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral traversal section
    // RDKit✔️✔️: const INT_LIST &trueOrder = atomTraversalBondOrder[atom->getIdx()];
    // RDKit✔️✔️: int nSwaps = 0;
    // RDKit✔️✔️: if (trueOrder.size() < atom->getDegree()) {
    // RDKit✔️✔️:   INT_LIST tOrder = trueOrder;
    // RDKit✔️✔️:   for (const auto bnd : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (std::find(trueOrder.begin(), trueOrder.end(), bnd->getIdx()) ==
    // RDKit✔️✔️:         trueOrder.end()) {
    // RDKit✔️✔️:       tOrder.push_back(bnd->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   nSwaps = atom->getPerturbationOrder(tOrder);
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   nSwaps = atom->getPerturbationOrder(trueOrder);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (!perm) {
    // RDKit✔️✔️:   nSwaps = atom->getPerturbationOrder(tOrder);
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
    // RDKit✔️✔️:   perm = Chirality::getChiralPermutation(atom, tOrder);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (doChiralInversions &&
    // RDKit✔️✔️:     chiralAtomNeedsTagInversion(mol, atom, firstInPart,
    // RDKit✔️✔️:                                 atomRingClosures[atom->getIdx()].size())) {
    // RDKit✔️✔️:   ++nSwaps;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (nSwaps % 2) {
    // RDKit✔️✔️:   numSwapsChiralAtoms.set(atom->getIdx());
    // RDKit✔️✔️: }
    // RDKit✔️✔️: atomPermutationIndices[atom->getIdx()] = perm;
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral traversal section
    let mut adjustments = vec![ChiralAdjustment::default(); topology.atoms.len()];
    for atom_index in 0..topology.atoms.len() {
        let atom = &topology.atoms[atom_index];
        if atom.chiral_tag() == ChiralTag::Unspecified {
            continue;
        }
        let incident = topology
            .adjacency
            .neighbors_of(atom_index)
            .iter()
            .map(|neighbor| neighbor.bond)
            .collect::<Vec<_>>();
        let mut traversal = atom_traversal_bond_order[atom_index].clone();
        if traversal.len() < incident.len() {
            for bond in &incident {
                if !traversal.contains(bond) {
                    traversal.push(*bond);
                }
            }
        }
        match atom.chiral_tag() {
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                if traversal.is_empty() {
                    continue;
                }
                let mut swaps = stereo::count_swaps_to_interconvert(&traversal, incident)
                    .ok_or_else(|| {
                        SmilesParseError::Model(
                            "writer traversal and storage bond orderings are not permutations"
                                .into(),
                        )
                    })?;
                let unsaturated =
                    topology
                        .adjacency
                        .neighbors_of(atom_index)
                        .iter()
                        .any(|neighbor| {
                            stereo::bond_order_as_double(
                                topology.bonds[neighbor.bond.index()].order(),
                            ) > 1.0
                        });
                if stereo::chiral_atom_needs_tag_inversion(
                    topology.adjacency.neighbors_of(atom_index).len(),
                    atom.explicit_hydrogens(),
                    atom_index == start_atom,
                    stereo::atom_has_fourth_valence(
                        atom.explicit_hydrogens(),
                        valence.implicit_hydrogens[atom_index] == 1,
                    ),
                    atom_ring_closures[atom_index].len(),
                    unsaturated,
                ) {
                    swaps += 1;
                }
                adjustments[atom_index].invert_tetrahedral = swaps % 2 == 1;
            }
            ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral => {
                let mut probe = traversal.into_iter().map(Some).collect::<Vec<_>>();
                stereo::insert_implicit_nontetrahedral_neighbors(
                    &mut probe,
                    atom.chiral_tag(),
                    atom_index == start_atom,
                );
                let permutation = stereo::nontetrahedral_chiral_permutation(
                    atom.chiral_permutation().unwrap_or(0),
                    atom.chiral_tag(),
                    topology.bonds.len(),
                    &incident,
                    &probe,
                    false,
                )
                .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
                if permutation != 0 {
                    adjustments[atom_index].nontetrahedral_permutation = Some(permutation);
                }
            }
            _ => {}
        }
    }
    apply_relative_chiral_adjustments(topology, stack, include_stereo_groups, &mut adjustments)?;
    Ok(adjustments)
}

fn apply_relative_chiral_adjustments(
    topology: &TopologyBlock,
    stack: &[MolStackElem],
    include_stereo_groups: bool,
    adjustments: &mut [ChiralAdjustment],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral post-processing section
    // RDKit✔️✔️: boost::dynamic_bitset<> ringStereoChemAdjusted(nAtoms);
    // RDKit✔️✔️: for (auto &msI : molStack) {
    // RDKit✔️✔️:   if (msI.type == MOL_STACK_ATOM &&
    // RDKit✔️✔️:       msI.obj.atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:       !msI.obj.atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit✔️✔️:     if (msI.obj.atom->hasProp(common_properties::_ringStereoAtoms)) {
    // RDKit✔️✔️:       if (!ringStereoChemAdjusted[msI.obj.atom->getIdx()]) {
    // RDKit✔️✔️:         msI.obj.atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:         ringStereoChemAdjusted.set(msI.obj.atom->getIdx());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       const INT_VECT &ringStereoAtoms = msI.obj.atom->getProp<INT_VECT>(
    // RDKit✔️✔️:           common_properties::_ringStereoAtoms);
    // RDKit✔️✔️:       for (auto nbrV : ringStereoAtoms) {
    // RDKit✔️✔️:         int nbrIdx = abs(nbrV) - 1;
    // RDKit✔️✔️:         if (!ringStereoChemAdjusted[nbrIdx] &&
    // RDKit✔️✔️:             atomVisitOrders[nbrIdx] >
    // RDKit✔️✔️:                 atomVisitOrders[msI.obj.atom->getIdx()]) {
    // RDKit✔️✔️:           mol.getAtomWithIdx(nbrIdx)->setChiralTag(
    // RDKit✔️✔️:               msI.obj.atom->getChiralTag());
    // RDKit✔️✔️:           if (nbrV < 0) {
    // RDKit✔️✔️:             mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           if (numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
    // RDKit✔️✔️:             if (!numSwapsChiralAtoms[nbrIdx]) {
    // RDKit✔️✔️:               mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             if (numSwapsChiralAtoms[nbrIdx]) {
    // RDKit✔️✔️:               mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           ringStereoChemAdjusted.set(nbrIdx);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (size_t sgidx;
    // RDKit✔️✔️:                msI.obj.atom->getPropIfPresent("_stereoGroup", sgidx) &&
    // RDKit✔️✔️:                mol.getStereoGroups().size() > sgidx) {
    // RDKit✔️✔️:       auto &sg = mol.getStereoGroups()[sgidx];
    // RDKit✔️✔️:       bool swapIt =
    // RDKit✔️✔️:           msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW;
    // RDKit✔️✔️:       if (swapIt) {
    // RDKit✔️✔️:         msI.obj.atom->invertChirality();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (swapIt || numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
    // RDKit✔️✔️:         for (auto at : sg.getAtoms()) {
    // RDKit✔️✔️:           if (at == msI.obj.atom) {
    // RDKit✔️✔️:             continue;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           at->invertChirality();
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if (msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:           msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:         if ((numSwapsChiralAtoms[msI.obj.atom->getIdx()])) {
    // RDKit✔️✔️:           msI.obj.atom->invertChirality();
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else if (atomPermutationIndices[msI.obj.atom->getIdx()]) {
    // RDKit✔️✔️:         msI.obj.atom->setProp(
    // RDKit✔️✔️:             common_properties::_chiralPermutation,
    // RDKit✔️✔️:             atomPermutationIndices[msI.obj.atom->getIdx()]);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral post-processing section
    let mut atom_visit_order = vec![usize::MAX; topology.atoms.len()];
    let mut visited_atoms = Vec::new();
    for (position, element) in stack.iter().enumerate() {
        if let MolStackElem::Atom(atom) = element {
            atom_visit_order[*atom] = position;
            visited_atoms.push(*atom);
        }
    }

    let mut stereo_group_reference = BTreeMap::new();
    if include_stereo_groups {
        for (group_index, group) in topology.stereo_groups.iter().enumerate() {
            if let Some(reference) = group
                .atoms()
                .iter()
                .map(|atom| atom.index())
                .filter(|atom| atom_visit_order[*atom] != usize::MAX)
                .min_by_key(|atom| atom_visit_order[*atom])
            {
                stereo_group_reference.insert(reference, group_index);
            }
        }
    }

    let mut ring_stereo_adjusted = vec![false; topology.atoms.len()];
    for atom_index in visited_atoms {
        let atom = &topology.atoms[atom_index];
        if atom.chiral_tag() == ChiralTag::Unspecified {
            continue;
        }
        if let Some(encoded) = atom.prop("_ringStereoAtoms") {
            let relations = parse_ring_stereo_atoms(encoded, topology.atoms.len())?;
            let source_inverted = adjustments[atom_index].invert_tetrahedral;
            if !ring_stereo_adjusted[atom_index] {
                adjustments[atom_index].chiral_tag_override = Some(ChiralTag::TetrahedralCcw);
                adjustments[atom_index].invert_tetrahedral = false;
                ring_stereo_adjusted[atom_index] = true;
            }
            let source_tag = adjustments[atom_index]
                .chiral_tag_override
                .unwrap_or(atom.chiral_tag());
            for (same_orientation, neighbor_index) in relations {
                if ring_stereo_adjusted[neighbor_index]
                    || atom_visit_order[neighbor_index] <= atom_visit_order[atom_index]
                {
                    continue;
                }
                let mut neighbor_tag = if same_orientation {
                    source_tag
                } else {
                    stereo::invert_tetrahedral_tag(source_tag)
                };
                if source_inverted != adjustments[neighbor_index].invert_tetrahedral {
                    neighbor_tag = stereo::invert_tetrahedral_tag(neighbor_tag);
                }
                adjustments[neighbor_index].chiral_tag_override = Some(neighbor_tag);
                adjustments[neighbor_index].invert_tetrahedral = false;
                ring_stereo_adjusted[neighbor_index] = true;
            }
        } else if let Some(group_index) = stereo_group_reference.get(&atom_index).copied() {
            let current_tag = adjustments[atom_index]
                .chiral_tag_override
                .unwrap_or(atom.chiral_tag());
            let swap_group = current_tag == ChiralTag::TetrahedralCw;
            if swap_group {
                adjustments[atom_index].chiral_tag_override =
                    Some(stereo::invert_tetrahedral_tag(current_tag));
            }
            if swap_group || adjustments[atom_index].invert_tetrahedral {
                for member in topology.stereo_groups[group_index].atoms() {
                    if member.index() == atom_index
                        || atom_visit_order[member.index()] == usize::MAX
                    {
                        continue;
                    }
                    let member_tag = adjustments[member.index()]
                        .chiral_tag_override
                        .unwrap_or(topology.atoms[member.index()].chiral_tag());
                    adjustments[member.index()].chiral_tag_override =
                        Some(stereo::invert_tetrahedral_tag(member_tag));
                }
            }
            adjustments[atom_index].invert_tetrahedral = false;
        }
    }
    Ok(())
}

fn parse_ring_stereo_atoms(
    encoded: &str,
    atom_count: usize,
) -> Result<Vec<(bool, usize)>, SmilesParseError> {
    let mut result = Vec::new();
    for token in encoded.split(',').filter(|token| !token.is_empty()) {
        let value = token.parse::<i64>().map_err(|_| {
            SmilesParseError::WriterStereo("`_ringStereoAtoms` contains a non-integer entry".into())
        })?;
        if value == 0 {
            return Err(SmilesParseError::WriterStereo(
                "`_ringStereoAtoms` cannot contain zero".into(),
            ));
        }
        let index = usize::try_from(value.unsigned_abs() - 1).map_err(|_| {
            SmilesParseError::WriterStereo("`_ringStereoAtoms` index is out of range".into())
        })?;
        if index >= atom_count {
            return Err(SmilesParseError::WriterStereo(
                "`_ringStereoAtoms` index is out of range".into(),
            ));
        }
        result.push((value > 0, index));
    }
    Ok(result)
}

fn write_mol_stack(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    stack: &[MolStackElem],
    chiral_adjustments: &[ChiralAdjustment],
) -> Result<String, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    // RDKit✔️✔️:   for (auto &mSE : molStack) {
    // RDKit✔️✔️:     switch (mSE.type) {
    // RDKit✔️✔️:       case Canon::MOL_STACK_ATOM:
    // RDKit✔️✔️:         for (auto rclosure : ringClosuresToErase) {
    // RDKit✔️✔️:           ringClosureMap.erase(rclosure);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         ringClosuresToErase.clear();
    // RDKit✔️✔️:         res << GetAtomSmiles(mSE.obj.atom, params);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Canon::MOL_STACK_BOND:
    // RDKit✔️✔️:         res << GetBondSmiles(bond, params, mSE.number);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Canon::MOL_STACK_RING:
    // RDKit✔️✔️:         if (ringClosureMap.count(ringIdx)) {
    // RDKit✔️✔️:           closureVal = ringClosureMap[ringIdx];
    // RDKit✔️✔️:           ringClosuresToErase.push_back(ringIdx);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           closureVal = 1;
    // RDKit✔️✔️:           ringClosureMap[ringIdx] = closureVal;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Canon::MOL_STACK_BRANCH_OPEN:
    // RDKit✔️✔️:         res << "(";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Canon::MOL_STACK_BRANCH_CLOSE:
    // RDKit✔️✔️:         res << ")";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    let mut output = String::new();
    let mut display_digits = BTreeMap::<usize, usize>::new();
    let mut closures_to_erase = Vec::new();
    for element in stack {
        match *element {
            MolStackElem::Atom(atom) => {
                for ring_id in closures_to_erase.drain(..) {
                    display_digits.remove(&ring_id);
                }
                output.push_str(&atom_text(
                    topology,
                    valence,
                    atom,
                    chiral_adjustments[atom],
                )?);
            }
            MolStackElem::Bond { bond, atom_to_left } => output.push_str(&bond_text(
                topology,
                &topology.bonds[bond.index()],
                atom_to_left,
            )?),
            MolStackElem::Ring(ring_id) => {
                let display_digit = if let Some(&digit) = display_digits.get(&ring_id) {
                    closures_to_erase.push(ring_id);
                    digit
                } else {
                    let digit = (1..)
                        .find(|candidate| !display_digits.values().any(|used| used == candidate))
                        .expect("positive ring-label space is not exhaustible");
                    display_digits.insert(ring_id, digit);
                    digit
                };
                write_ring_label(&mut output, display_digit);
            }
            MolStackElem::BranchOpen => output.push('('),
            MolStackElem::BranchClose => output.push(')'),
        }
    }
    Ok(output)
}

fn atom_text(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_index: usize,
    chiral_adjustment: ChiralAdjustment,
) -> Result<String, SmilesParseError> {
    let atom = &topology.atoms[atom_index];
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles modeled non-stereo atom emission
    // RDKit❗✔️:   bool needsBracket = true;
    // RDKit❗✔️:   if (!hasCustomSymbol && !params.allHsExplicit) {
    // RDKit❗✔️:     needsBracket = atomNeedsBracket(atom, atString, params);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (needsBracket) {
    // RDKit❗✔️:     res += "[";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (isotope && params.doIsomericSmiles) {
    // RDKit❗✔️:     res += std::to_string(isotope);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   res += symb;
    // RDKit❗✔️:   if (needsBracket) {
    // RDKit❗✔️:     unsigned int totNumHs = atom->getTotalNumHs();
    // RDKit❗✔️:     if (totNumHs > 0) {
    // RDKit❗✔️:       res += "H";
    // RDKit❗✔️:       if (totNumHs > 1) {
    // RDKit❗✔️:         res += std::to_string(totNumHs);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetAtomSmiles modeled non-stereo atom emission
    // Detached records do not yet carry the source writer's calculated total
    // valence cache, so bracket selection is exact for syntax-owned fields and
    // deliberately remains partial for cache-derived nonstandard valence.
    let custom_symbol = atom.prop("smilesSymbol");
    let raw_symbol = custom_symbol.unwrap_or_else(|| atom.element().symbol());
    let symbol = if atom.is_aromatic()
        && matches!(
            atom.atomic_number(),
            5 | 6 | 7 | 8 | 14 | 15 | 16 | 33 | 34 | 52
        ) {
        let mut lowered = raw_symbol.to_owned();
        lowered.get_mut(0..1).map(str::make_ascii_lowercase);
        lowered
    } else {
        raw_symbol.to_owned()
    };
    // BEGIN RDKIT CPP FUNCTION getAtomChiralityInfo
    // RDKit✔️✔️: std::string getAtomChiralityInfo(const Atom *atom) {
    // RDKit✔️✔️:   auto allowNontet = Chirality::getAllowNontetrahedralChirality();
    // RDKit✔️✔️:   std::string atString;
    // RDKit✔️✔️:   switch (atom->getChiralTag()) {
    // RDKit✔️✔️:     case Atom::CHI_TETRAHEDRAL_CW: atString = "@@"; break;
    // RDKit✔️✔️:     case Atom::CHI_TETRAHEDRAL_CCW: atString = "@"; break;
    // RDKit✔️✔️:     default: break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atString.empty() && allowNontet) {
    // RDKit✔️✔️:     switch (atom->getChiralTag()) {
    // RDKit✔️✔️:       case Atom::CHI_SQUAREPLANAR: atString = "@SP"; break;
    // RDKit✔️✔️:       case Atom::CHI_TRIGONALBIPYRAMIDAL: atString = "@TB"; break;
    // RDKit✔️✔️:       case Atom::CHI_OCTAHEDRAL: atString = "@OH"; break;
    // RDKit✔️✔️:       default: break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!atString.empty()) {
    // RDKit✔️✔️:       int permutation = 0;
    // RDKit✔️✔️:       if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER &&
    // RDKit✔️✔️:           atom->getPropIfPresent(common_properties::_chiralPermutation,
    // RDKit✔️✔️:                                  permutation) &&
    // RDKit✔️✔️:           !SmilesParseOps::checkChiralPermutation(atom->getChiralTag(),
    // RDKit✔️✔️:                                                   permutation)) {
    // RDKit✔️✔️:         throw ValueErrorException("bad chirality spec");
    // RDKit✔️✔️:       } else if (permutation) {
    // RDKit✔️✔️:         atString += std::to_string(permutation);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return atString;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getAtomChiralityInfo
    let base_chiral_tag = chiral_adjustment
        .chiral_tag_override
        .unwrap_or(atom.chiral_tag());
    let chiral_tag = if chiral_adjustment.invert_tetrahedral {
        stereo::invert_tetrahedral_tag(base_chiral_tag)
    } else {
        base_chiral_tag
    };
    let mut chirality = match chiral_tag {
        ChiralTag::TetrahedralCw => "@@".to_owned(),
        ChiralTag::TetrahedralCcw => "@".to_owned(),
        ChiralTag::SquarePlanar => "@SP".to_owned(),
        ChiralTag::TrigonalBipyramidal => "@TB".to_owned(),
        ChiralTag::Octahedral => "@OH".to_owned(),
        _ => String::new(),
    };
    if matches!(
        chiral_tag,
        ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral
    ) {
        let permutation = chiral_adjustment
            .nontetrahedral_permutation
            .or_else(|| atom.chiral_permutation())
            .unwrap_or(0);
        let limit = match chiral_tag {
            ChiralTag::SquarePlanar => 3,
            ChiralTag::TrigonalBipyramidal => 20,
            ChiralTag::Octahedral => 30,
            _ => unreachable!(),
        };
        if permutation > limit {
            return Err(SmilesParseError::WriterStereo(format!(
                "invalid {} permutation {permutation}; maximum is {limit}",
                chiral_tag.rdkit_name()
            )));
        }
        if permutation != 0 {
            chirality.push_str(&permutation.to_string());
        }
    }
    let total_num_hydrogens = usize::from(atom.explicit_hydrogens())
        + usize::try_from(valence.implicit_hydrogens[atom_index].max(0)).unwrap_or(usize::MAX);
    let total_valence =
        valence.explicit_valence[atom_index] + valence.implicit_hydrogens[atom_index].max(0);
    let nonstandard_valence = if atom.radical_electrons() != 0
        || (matches!(atom.atomic_number(), 7 | 15)
            && atom.is_aromatic()
            && atom.explicit_hydrogens() > 0)
    {
        true
    } else {
        cosmolkit_core::rdkit_valence_list(atom.atomic_number())
            .ok()
            .flatten()
            .and_then(|values| values.first())
            .is_some_and(|default| total_valence != *default && total_num_hydrogens != 0)
    };
    let bonded_to_metal = topology
        .adjacency
        .neighbors_of(atom_index)
        .iter()
        .any(|neighbor| {
            rdkit_query_ops_is_metal(topology.atoms[neighbor.atom_index].atomic_number())
        });
    let needs_bracket = custom_symbol.is_some()
        || !in_organic_subset(atom.atomic_number())
        || atom.formal_charge() != 0
        || atom.isotope().is_some()
        || atom.atom_map().is_some()
        || !chirality.is_empty()
        || nonstandard_valence
        || bonded_to_metal;
    if !needs_bracket {
        return Ok(append_supplemental_label(atom, symbol));
    }

    let mut output = String::from("[");
    if let Some(isotope) = atom.isotope() {
        output.push_str(&isotope.to_string());
    }
    output.push_str(&symbol);
    output.push_str(&chirality);
    if total_num_hydrogens > 0 {
        output.push('H');
        if total_num_hydrogens > 1 {
            output.push_str(&total_num_hydrogens.to_string());
        }
    }
    match atom.formal_charge() {
        0 => {}
        1 => output.push('+'),
        -1 => output.push('-'),
        charge if charge > 1 => {
            output.push('+');
            output.push_str(&charge.to_string());
        }
        charge => output.push_str(&charge.to_string()),
    }
    if let Some(atom_map) = atom.atom_map() {
        output.push(':');
        output.push_str(&atom_map.to_string());
    }
    output.push(']');
    Ok(append_supplemental_label(atom, output))
}

fn append_supplemental_label(atom: &Atom, mut text: String) -> String {
    if let Some(label) = atom.prop("_supplementalSmilesLabel") {
        text.push_str(label);
    }
    text
}

fn bond_text(
    topology: &TopologyBlock,
    bond: &Bond,
    atom_to_left: usize,
) -> Result<&'static str, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION GetBondSmiles modeled non-stereo bond emission
    // RDKit✔️✔️:     case Bond::SINGLE:
    // RDKit✔️✔️:       if (params.allBondsExplicit) {
    // RDKit✔️✔️:         res = "-";
    // RDKit✔️✔️:       } else if (aromatic && !bond->getIsAromatic()) {
    // RDKit✔️✔️:         res = "-";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::DOUBLE:
    // RDKit✔️✔️:       if (!aromatic || !bond->getIsAromatic() || params.allBondsExplicit) {
    // RDKit✔️✔️:         res = "=";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::TRIPLE:
    // RDKit✔️✔️:       res = "#";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::QUADRUPLE:
    // RDKit✔️✔️:       res = "$";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::AROMATIC:
    // RDKit✔️✔️:       if (params.allBondsExplicit || !aromatic) {
    // RDKit✔️✔️:         res = ":";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::DATIVE:
    // RDKit✔️✔️:       if (atomToLeftIdx >= 0 &&
    // RDKit✔️✔️:           bond->getBeginAtomIdx() == static_cast<unsigned int>(atomToLeftIdx)) {
    // RDKit✔️✔️:         res = "->";
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         res = "<-";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       res = "~";
    // END RDKIT CPP FUNCTION GetBondSmiles modeled non-stereo bond emission
    let other = other_atom(bond, atom_to_left)?;
    let aromatic_context = matches!(
        bond.order(),
        BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
    ) && topology.atoms[atom_to_left].is_aromatic()
        && topology.atoms[other].is_aromatic()
        && (topology.atoms[atom_to_left].atomic_number() != 0
            || topology.atoms[other].atomic_number() != 0);
    if matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic) {
        // RDKit✔️✔️:       if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
        // RDKit✔️✔️:         switch (dir) {
        // RDKit✔️✔️:           case Bond::ENDDOWNRIGHT:
        // RDKit✔️✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
        // RDKit✔️✔️:               res = "\\\\";
        // RDKit✔️✔️:             }
        // RDKit✔️✔️:             break;
        // RDKit✔️✔️:           case Bond::ENDUPRIGHT:
        // RDKit✔️✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
        // RDKit✔️✔️:               res = "/";
        // RDKit✔️✔️:             }
        // RDKit✔️✔️:             break;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        match bond.direction() {
            BondDirection::EndDownRight => return Ok("\\"),
            BondDirection::EndUpRight => return Ok("/"),
            _ => {}
        }
    }
    Ok(match bond.order() {
        BondOrder::Single if aromatic_context && !bond.is_aromatic() => "-",
        BondOrder::Single => "",
        BondOrder::Double if aromatic_context && bond.is_aromatic() => "",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Quadruple => "$",
        BondOrder::Aromatic if aromatic_context => "",
        BondOrder::Aromatic => ":",
        BondOrder::Dative if bond.begin().index() == atom_to_left => "->",
        BondOrder::Dative => "<-",
        _ => "~",
    })
}

fn other_atom(bond: &Bond, atom: usize) -> Result<usize, SmilesParseError> {
    if bond.begin().index() == atom {
        Ok(bond.end().index())
    } else if bond.end().index() == atom {
        Ok(bond.begin().index())
    } else {
        Err(SmilesParseError::Model(format!(
            "bond {} is not incident to atom {atom}",
            bond.id().index()
        )))
    }
}

fn in_organic_subset(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION inOrganicSubset
    // RDKit✔️✔️: const int atomicSmiles[] = {0, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53, -1};
    // RDKit✔️✔️: bool inOrganicSubset(int atomicNumber) {
    // RDKit✔️✔️:   unsigned int idx = 0;
    // RDKit✔️✔️:   while (atomicSmiles[idx] < atomicNumber && atomicSmiles[idx] != -1) {
    // RDKit✔️✔️:     ++idx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return atomicSmiles[idx] == atomicNumber;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION inOrganicSubset
    matches!(
        atomic_number,
        0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53
    )
}

fn rdkit_query_ops_is_metal(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::makeMAtomQuery / QueryOps::isMetal
    // RDKit✔️✔️: // !#0!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️: bool isMetal(const Atom &atom) {
    // RDKit✔️✔️:   static const std::unique_ptr<ATOM_OR_QUERY> q(makeMAtomQuery());
    // RDKit✔️✔️:   return q->Match(&atom);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryOps::makeMAtomQuery / QueryOps::isMetal
    !matches!(
        atomic_number,
        0 | 1
            | 2
            | 5
            | 6
            | 7
            | 8
            | 9
            | 10
            | 14
            | 15
            | 16
            | 17
            | 18
            | 33
            | 34
            | 35
            | 36
            | 52
            | 53
            | 54
            | 85
            | 86
    )
}

fn write_ring_label(output: &mut String, label: usize) {
    // RDKit✔️✔️:         if (closureVal < 10) {
    // RDKit✔️✔️:           res << (char)(closureVal + '0');
    // RDKit✔️✔️:         } else if (closureVal < 100) {
    // RDKit✔️✔️:           res << '%' << closureVal;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           res << "%(" << closureVal << ')';
    // RDKit✔️✔️:         }
    if label < 10 {
        output.push(char::from(b'0' + label as u8));
    } else if label < 100 {
        output.push('%');
        output.push_str(&label.to_string());
    } else {
        output.push_str("%(");
        output.push_str(&label.to_string());
        output.push(')');
    }
}

fn rdkit_bond_type_code(order: BondOrder) -> i64 {
    match order {
        BondOrder::Unspecified => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Quintuple => 5,
        BondOrder::Hextuple => 6,
        BondOrder::OneAndHalf => 7,
        BondOrder::TwoAndHalf => 8,
        BondOrder::ThreeAndHalf => 9,
        BondOrder::FourAndHalf => 10,
        BondOrder::FiveAndHalf => 11,
        BondOrder::Aromatic => 12,
        BondOrder::Ionic => 13,
        BondOrder::Hydrogen => 14,
        BondOrder::ThreeCenter => 15,
        BondOrder::DativeOne => 16,
        BondOrder::Dative => 17,
        BondOrder::DativeLeft => 18,
        BondOrder::DativeRight => 19,
        BondOrder::Other => 20,
        BondOrder::Zero => 21,
    }
}

fn find_ring_bonds(topology: &TopologyBlock) -> Vec<bool> {
    // BEGIN RDKIT CPP FUNCTION findSSSR active bond selection
    // RDKit✔️🔝:   // Zero-order bonds are not candidates for rings, and dative bonds and
    // RDKit✔️🔝:   // hydrogen bonds may also be out
    // RDKit✔️🔝:   boost::dynamic_bitset<> activeBonds(nbnds);
    // RDKit✔️🔝:   activeBonds.set();
    // RDKit✔️🔝:   for (auto bond : mol.bonds()) {
    // RDKit✔️🔝:     if (auto bt = bond->getBondType();
    // RDKit✔️🔝:         bt == Bond::ZERO || (!includeDativeBonds && isDative(bt)) ||
    // RDKit✔️🔝:         (!includeHydrogenBonds && bt == Bond::HYDROGEN)) {
    // RDKit✔️🔝:       activeBonds[bond->getIdx()] = 0;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // END RDKIT CPP FUNCTION findSSSR active bond selection
    // A bond belongs to some cycle exactly when it is not a bridge in the
    // source-selected active-bond graph. The source later materializes an
    // SSSR only to query `numBondRings() != 0`; computing that boolean directly
    // is O(V+E), avoids ring-vector allocation, and preserves its semantics.
    let mut discovery = vec![usize::MAX; topology.atoms.len()];
    let mut low = vec![usize::MAX; topology.atoms.len()];
    let mut bridges = vec![false; topology.bonds.len()];
    let mut next_time = 0;

    fn visit(
        topology: &TopologyBlock,
        atom: usize,
        parent_bond: Option<BondId>,
        discovery: &mut [usize],
        low: &mut [usize],
        bridges: &mut [bool],
        next_time: &mut usize,
    ) {
        discovery[atom] = *next_time;
        low[atom] = *next_time;
        *next_time += 1;
        for neighbor in topology.adjacency.neighbors_of(atom) {
            if Some(neighbor.bond) == parent_bond {
                continue;
            }
            if !ring_perception_eligible(topology.bonds[neighbor.bond.index()].order()) {
                continue;
            }
            if discovery[neighbor.atom_index] == usize::MAX {
                visit(
                    topology,
                    neighbor.atom_index,
                    Some(neighbor.bond),
                    discovery,
                    low,
                    bridges,
                    next_time,
                );
                low[atom] = low[atom].min(low[neighbor.atom_index]);
                if low[neighbor.atom_index] > discovery[atom] {
                    bridges[neighbor.bond.index()] = true;
                }
            } else {
                low[atom] = low[atom].min(discovery[neighbor.atom_index]);
            }
        }
    }

    for atom in 0..topology.atoms.len() {
        if discovery[atom] == usize::MAX {
            visit(
                topology,
                atom,
                None,
                &mut discovery,
                &mut low,
                &mut bridges,
                &mut next_time,
            );
        }
    }
    bridges
        .into_iter()
        .enumerate()
        .map(|(index, is_bridge)| {
            ring_perception_eligible(topology.bonds[index].order()) && !is_bridge
        })
        .collect()
}

fn ring_perception_eligible(order: BondOrder) -> bool {
    !matches!(
        order,
        BondOrder::Zero
            | BondOrder::Dative
            | BondOrder::DativeOne
            | BondOrder::DativeLeft
            | BondOrder::DativeRight
            | BondOrder::Hydrogen
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SmilesParseParams, parse_smiles};

    fn roundtrip(input: &str) -> String {
        let record = parse_smiles(input, &SmilesParseParams::default()).expect("parse");
        write_smiles_with_params(&record, &SmilesWriteParams { canonical: false }).expect("write")
    }

    #[test]
    fn source_order_writer_places_all_but_last_child_in_branches() {
        assert_eq!(roundtrip("C(O)N"), "C(O)N");
        assert_eq!(roundtrip("OC.C"), "OC.C");
    }

    #[test]
    fn writer_serializes_aliphatic_and_aromatic_cycles() {
        assert_eq!(roundtrip("C1CCCCC1"), "C1CCCCC1");
        assert_eq!(roundtrip("c1ccccc1"), "c1ccccc1");
        assert_eq!(roundtrip("c1ccccc-1"), "c1ccccc-1");
    }

    #[test]
    fn writer_matches_source_order_for_fused_bridged_and_spiro_cycles() {
        for smiles in [
            "C1CCC2CCCCC2C1",
            "C12(CCCCC1)CCCCC2",
            "C1C2C3C1C2C3",
            "C1CC2CCC1C2",
            "C1(C2CC2)CC1",
            "C1CC2(CC1)CCC2",
        ] {
            assert_eq!(roundtrip(smiles), smiles);
        }
    }

    #[test]
    fn writer_preserves_dative_ring_orientation() {
        assert_eq!(roundtrip("N1CC->1"), "N1CC->1");
        assert_eq!(roundtrip("N1CC<-1"), "N1CC<-1");
    }

    #[test]
    fn writer_emits_bracket_fields_without_loss() {
        assert_eq!(roundtrip("[13CH3+:7]C"), "[13CH3+:7]C");
        assert_eq!(roundtrip("[Na+]Cl"), "[Na+][Cl]");
    }

    #[test]
    fn writer_matches_rdkit_tetrahedral_traversal_inversions() {
        for (input, expected) in [
            ("[C@H](F)(Cl)Br", "F[C@@H](Cl)Br"),
            ("[C@@H](F)(Cl)Br", "F[C@H](Cl)Br"),
            ("F[C@H](Cl)Br", "F[C@H](Cl)Br"),
            ("Br[C@@H](Cl)F", "F[C@H](Cl)Br"),
            ("N[C@](F)(Cl)Br", "N[C@](F)(Cl)Br"),
            ("C[C@H](F)Cl", "C[C@H](F)Cl"),
            ("F[C@]1(Br)CCO1", "F[C@]1(Br)CCO1"),
            ("F[C@]1(CCO1)Br", "F[C@@]1(Br)CCO1"),
            ("[C@H](F)(Cl)Br.[C@@H](I)(N)O", "F[C@@H](Cl)Br.N[C@H](O)I"),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            assert_eq!(write_smiles(&record).unwrap(), expected, "{input}");
        }

        let bond_stereo = parse_smiles("F/C=C/F", &Default::default()).unwrap();
        assert_eq!(write_smiles(&bond_stereo).unwrap(), "F/C=C/F");
    }

    #[test]
    fn writer_matches_rdkit_nontetrahedral_permutation_reordering() {
        for (input, expected_noncanonical, expected_canonical) in [
            (
                "[Pt@SP1](F)(Cl)(Br)I",
                "[Pt@SP1]([F])([Cl])([Br])[I]",
                "[F][Pt@SP1]([Cl])([Br])[I]",
            ),
            (
                "[Pt@SP2](F)(Cl)(Br)I",
                "[Pt@SP2]([F])([Cl])([Br])[I]",
                "[F][Pt@SP2]([Cl])([Br])[I]",
            ),
            (
                "I[Pt@SP1](Br)(Cl)F",
                "[I][Pt@SP1]([Br])([Cl])[F]",
                "[F][Pt@SP1]([Cl])([Br])[I]",
            ),
            (
                "[Pt@SP1](F)(Cl)Br",
                "[Pt@SP1]([F])([Cl])[Br]",
                "[F][Pt@SP3]([Cl])[Br]",
            ),
            (
                "[P@TB1](F)(Cl)(Br)(I)N",
                "[P@TB1](F)(Cl)(Br)(I)N",
                "N[P@TB8](F)(Cl)(Br)I",
            ),
            (
                "[P@TB20](F)(Cl)(Br)(I)N",
                "[P@TB20](F)(Cl)(Br)(I)N",
                "N[P@TB3](F)(Cl)(Br)I",
            ),
            (
                "N[P@TB1](I)(Br)(Cl)F",
                "N[P@TB1](I)(Br)(Cl)F",
                "N[P@TB8](F)(Cl)(Br)I",
            ),
            (
                "[P@TB1](F)(Cl)(Br)I",
                "[P@TB1](F)(Cl)(Br)I",
                "F[P@TB9](Cl)(Br)I",
            ),
            (
                "[Co@OH1](F)(Cl)(Br)(I)(N)O",
                "[Co@OH1]([F])([Cl])([Br])([I])([NH2])[OH]",
                "[NH2][Co@OH9]([OH])([F])([Cl])([Br])[I]",
            ),
            (
                "[Co@OH30](F)(Cl)(Br)(I)(N)O",
                "[Co@OH30]([F])([Cl])([Br])([I])([NH2])[OH]",
                "[NH2][Co@OH15]([OH])([F])([Cl])([Br])[I]",
            ),
            (
                "O[Co@OH1](N)(I)(Br)(Cl)F",
                "[OH][Co@OH1]([NH2])([I])([Br])([Cl])[F]",
                "[NH2][Co@OH9]([OH])([F])([Cl])([Br])[I]",
            ),
            (
                "[Co@OH1](F)(Cl)(Br)(I)N",
                "[Co@OH1]([F])([Cl])([Br])([I])[NH2]",
                "[NH2][Co@OH30]([F])([Cl])([Br])[I]",
            ),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            assert_eq!(
                write_smiles_with_params(&record, &SmilesWriteParams { canonical: false }).unwrap(),
                expected_noncanonical,
                "non-canonical {input}"
            );
            assert_eq!(
                write_smiles(&record).unwrap(),
                expected_canonical,
                "canonical {input}"
            );
        }
    }

    #[test]
    fn writer_canonicalizes_ring_relative_stereo_like_rdkit() {
        for (input, relation, expected) in [
            (
                "C1[C@H](F)CC[C@H](Cl)C1",
                ("6", "2"),
                "F[C@H]1CC[C@@H](Cl)CC1",
            ),
            (
                "C1[C@H](F)CC[C@@H](Cl)C1",
                ("-6", "-2"),
                "F[C@H]1CC[C@H](Cl)CC1",
            ),
        ] {
            let mut record = parse_smiles(input, &Default::default()).unwrap();
            record.topology.atoms[1].set_prop("_ringStereoAtoms", relation.0);
            record.topology.atoms[5].set_prop("_ringStereoAtoms", relation.1);
            assert_eq!(write_smiles(&record).unwrap(), expected, "{input}");
        }
    }

    #[test]
    fn writer_rejects_malformed_ring_relative_stereo_state() {
        let mut record = parse_smiles("C1[C@H](F)CC[C@H](Cl)C1", &Default::default()).unwrap();
        record.topology.atoms[1].set_prop("_ringStereoAtoms", "0");
        assert!(matches!(
            write_smiles(&record),
            Err(SmilesParseError::WriterStereo(_))
        ));
    }

    #[test]
    fn canonical_writer_matches_rdkit_double_bond_directions() {
        let mut mismatches = Vec::new();
        for (input, expected) in [
            ("F/C=C/F", "F/C=C/F"),
            ("F/C=C\\F", "F/C=C\\F"),
            ("Br\\C=C/F", "F/C=C\\Br"),
            ("F\\C=C\\Br", "F/C=C/Br"),
            ("F\\C=C(/Cl)\\Br", "F/C=C(\\Cl)Br"),
            ("F/C=C(\\Cl)/Br", "F/C=C(\\Cl)Br"),
            ("C/C=C/C", "C/C=C/C"),
            ("C/C=C\\C", "C/C=C\\C"),
            ("C/C=C/C=C\\C", "C/C=C\\C=C\\C"),
            ("Cl/C=C(/C=C/C)\\C=C\\Br", "C/C=C/C(=C/Cl)/C=C/Br"),
            ("C(\\C/C=C/Cl)=C/O", "O/C=C/C/C=C/Cl"),
            ("O=C\\C=C/F", "O=C/C=C\\F"),
            ("C(=O)\\C=C/Br", "O=C/C=C\\Br"),
            ("CC(=O)\\C=C/Br", "CC(=O)/C=C\\Br"),
            ("C/C=C(/C)C", "CC=C(C)C"),
            ("C/C=C(/C(F))C(Cl)", "C/C=C(/CF)CCl"),
            ("C/C=C(/C(F))C(F)", "CC=C(CF)CF"),
            ("C/C=C(/CO)CN", "C/C=C(\\CN)CO"),
            ("F/C=C(/F)Cl", "F/C=C(/F)Cl"),
            ("F/C(Cl)=C(/Br)I", "F/C(Cl)=C(/Br)I"),
            ("C/C=C(/C)\\C", "CC=C(C)C"),
            ("C1COC/C=C\\CCC1", "C1=C\\COCCCCC/1"),
            ("C1COC/C=C/CCC1", "C1=C/COCCCCC/1"),
            ("C1CC/C=C/C=C/CCC1", "C1=C/CCCCCC/C=C/1"),
            ("C/1=C/C=C/CCCCCC1", "C1=C\\CCCCCC/C=C/1"),
            ("C1COC/C=C/C=C/C1", "C1=C/CCCOC/C=C/1"),
            ("C1=C/OCC/C=C\\CC\\1", "C1=C\\CCO/C=C\\CC/1"),
            ("CO/C1=C/C=C\\C=C/C=N\\1", "COC1=C/C=C\\C=C/C=N\\1"),
            ("C1C/C=C/CCCCCCCC1", "C1=C/CCCCCCCCCC/1"),
            ("C1/C=C/C=C/CCCCCCCCC1", "C1=C/CCCCCCCCCC/C=C/1"),
            ("C1/C=C\\CCCCC1", "C1=C\\CCCCCC/1"),
            // This crate's parse boundary is deliberately unsanitized. RDKit
            // with sanitize=false/removeHs=false preserves these directions;
            // its later full sanitize pipeline removes them as non-potential
            // fused-ring stereo.
            ("C1/C=C\\C2CCCCC2C1", "C1=C\\C2CCCCC2CC/1"),
            ("C1C/C=C/CC2CCCCC12", "C1=C/CC2CCCCC2CC/1"),
            ("F/C=C/F.C/C=C\\C", "C/C=C\\C.F/C=C/F"),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            let actual = write_smiles(&record).unwrap();
            if actual != expected {
                mismatches.push((input, expected, actual));
            }
        }
        assert!(mismatches.is_empty(), "{mismatches:#?}");
    }

    #[test]
    fn detached_legacy_cip_ranks_match_rdkit_compute_atom_cip_ranks() {
        for (input, expected) in [
            ("CCO", vec![0, 1, 2]),
            ("Cl/C=C(/C=C/C)\\C=C\\Br", vec![7, 5, 4, 2, 1, 0, 3, 6, 8]),
            ("C1CCCCC1", vec![0, 0, 0, 0, 0, 0]),
            ("c1ccccc1", vec![0, 0, 0, 0, 0, 0]),
            ("CC(C)O", vec![0, 1, 0, 2]),
            ("N#CC(=O)O", vec![2, 0, 1, 4, 3]),
            ("F[C@H](Cl)Br", vec![1, 0, 2, 3]),
            ("OP(=O)(O)O", vec![1, 2, 0, 1, 1]),
            ("[Na+].[O-]C=O", vec![3, 1, 0, 2]),
            ("[CH3:7]CO", vec![1, 0, 2]),
            ("[13CH3]C", vec![1, 0]),
            ("[12CH3][13CH3]", vec![0, 1]),
            ("[35Cl]C[37Cl]", vec![1, 0, 2]),
            ("[127I]C[129I]", vec![1, 0, 2]),
            ("[128Te]C[130Te]", vec![1, 0, 2]),
            ("C1/C=C\\C2CCCCC2C1", vec![5, 8, 9, 7, 4, 1, 0, 2, 6, 3]),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            let valence = cosmolkit_core::assign_valence_with_options_for_topology(
                &record.topology,
                ValenceModel::RdkitLike,
                false,
            )
            .unwrap();
            assert_eq!(
                cosmolkit_core::assign_atom_cip_ranks(&record.topology, &valence).unwrap(),
                expected,
                "{input}"
            );
        }
    }

    #[test]
    fn ring_label_format_matches_smiles_writer_extensions() {
        let mut output = String::new();
        write_ring_label(&mut output, 9);
        write_ring_label(&mut output, 10);
        write_ring_label(&mut output, 100);
        assert_eq!(output, "9%10%(100)");
    }

    #[test]
    fn canonical_writer_matches_rdkit_for_non_stereo_components_and_cycles() {
        for (input, expected) in [
            ("OCC", "CCO"),
            ("C(C)O", "CCO"),
            ("OC.C", "C.CO"),
            ("C.O.CC", "C.CC.O"),
            ("C1CCCCC1", "C1CCCCC1"),
            ("C1CCC2CCCCC2C1", "C1CCC2CCCCC2C1"),
            ("C12(CCCCC1)CCCCC2", "C1CCC2(CC1)CCCCC2"),
            ("C1C2C3C1C2C3", "C1C2C3CC2C13"),
            ("OC(C)C", "CC(C)O"),
            ("n1ccccc1", "c1ccncc1"),
            ("N#C", "C#N"),
            ("[CH3:7]CO", "OC[CH3:7]"),
            ("[Na+]Cl", "[Na+][Cl]"),
            ("CCN(CC)CC", "CCN(CC)CC"),
            ("CC(C)(C)O", "CC(C)(C)O"),
            ("O=C(O)C", "CC(=O)O"),
            ("C#CC=C", "C#CC=C"),
            ("N=C=O", "N=C=O"),
            ("[NH4+]", "[NH4+]"),
            ("[O-]C=O", "O=C[O-]"),
            ("[13CH3]C", "C[13CH3]"),
            ("[*:1]CC", "CC[*:1]"),
            ("B(O)O", "OBO"),
            ("OP(=O)(O)O", "O=P(O)(O)O"),
            ("c1cc[nH]c1", "c1cc[nH]c1"),
            ("[nH]1cccc1", "c1cc[nH]c1"),
            ("c1ccccc1O", "Oc1ccccc1"),
            ("C1=CC2=CC=CC=C2C=C1", "C1=CC=C2C=CC=CC2=C1"),
            ("C1CC1C2CC2", "C1CC1C1CC1"),
            ("N->B", "B<-N"),
            ("N1CC->1", "C1C->N1"),
            ("BrCCl", "ClCBr"),
            ("[SiH4]", "[SiH4]"),
            ("[AsH]", "[AsH]"),
            ("[se]1cccc1", "c1cc[se]c1"),
            ("[C]", "C"),
            ("[CH3]", "[CH3]"),
            ("[NH2-]", "[NH2-]"),
            ("[O:2]=[C:1]O", "O[C:1]=[O:2]"),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            assert_eq!(write_smiles(&record).unwrap(), expected, "{input}");
        }
    }
}
