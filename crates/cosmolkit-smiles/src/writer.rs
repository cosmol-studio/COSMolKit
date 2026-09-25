use std::collections::BTreeMap;

use cosmolkit_core::{KekulizeParams, ValenceAssignment, ValenceModel, fast_find_rings_from_parts};
use cosmolkit_model::{Atom, AtomId, Bond, BondId, StereoGroup, TopologyBlock};
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

#[derive(Debug)]
struct WriterStereoFragment {
    topology: TopologyBlock,
    source_atoms: Vec<AtomId>,
    source_bonds: Vec<BondId>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmilesWriteParams {
    pub do_isomeric_smiles: bool,
    pub do_kekule: bool,
    pub canonical: bool,
    pub clean_stereo: bool,
    pub rooted_at_atom: Option<AtomId>,
    pub all_bonds_explicit: bool,
    pub all_hydrogens_explicit: bool,
    pub include_dative_bonds: bool,
    pub ignore_atom_map_numbers: bool,
}

impl Default for SmilesWriteParams {
    fn default() -> Self {
        Self {
            do_isomeric_smiles: true,
            do_kekule: false,
            canonical: true,
            clean_stereo: true,
            rooted_at_atom: None,
            all_bonds_explicit: false,
            all_hydrogens_explicit: false,
            include_dative_bonds: true,
            ignore_atom_map_numbers: false,
        }
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
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles root validation
    // RDKit✔️❌:   if (!mol.getNumAtoms()) {
    // RDKit✔️❌:     return "";
    // RDKit✔️❌:   }
    // RDKit✔️❌:   PRECONDITION(
    // RDKit✔️❌:       params.rootedAtAtom < 0 ||
    // RDKit✔️❌:           static_cast<unsigned int>(params.rootedAtAtom) < mol.getNumAtoms(),
    // RDKit✔️❌:       "rootedAtAtom must be less than the number of atoms");
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles root validation
    // Detached records need structural validation before indexed access; this
    // adds an O(V + E) pass before the source's constant-time root precondition.
    record
        .topology
        .validate()
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;
    if record.topology.atoms.is_empty() {
        return Ok(SmilesWriteOutput {
            text: String::new(),
            atom_order: Vec::new(),
            bond_order: Vec::new(),
        });
    }
    if let Some(root) = params.rooted_at_atom {
        if root.index() >= record.topology.atoms.len() {
            return Err(SmilesParseError::WriterRootAtomOutOfRange {
                atom_index: root.index(),
                atom_count: record.topology.atoms.len(),
            });
        }
    }
    reject_unmodeled_stereochemical_writing(&record.topology, doing_cx_smiles)?;

    // RDKit's public writer performs stereo preparation and direction
    // canonicalization on a private molecule copy. Keep the detached record
    // immutable and make those same temporary changes on this working block.
    let mut topology = record.topology.clone();
    let original_atom_maps = if params.ignore_atom_map_numbers {
        let original = params.canonical.then(|| {
            topology
                .atoms
                .iter()
                .map(|atom| atom.atom_map().filter(|map| *map != 0))
                .collect::<Vec<_>>()
        });
        // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles map setup
        // RDKit✔️✔️:       if (params.ignoreAtomMapNumbers) {
        // RDKit✔️✔️:         atomMapNums[atom->getIdx()] = atom->getAtomMapNum();
        // RDKit✔️✔️:         atom->setAtomMapNum(0);
        // RDKit✔️✔️:       }
        // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles map setup
        for atom in &mut topology.atoms {
            atom.set_atom_map(None);
        }
        original
    } else {
        None
    };
    let source_valence = if doing_cx_smiles || !params.include_dative_bonds {
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
    let components = connected_components(&topology);
    let stereochem_done_marker_is_computed = record
        .properties
        .prop("_StereochemDone")
        .map(|_| record.properties.is_prop_computed("_StereochemDone"));
    prepare_writer_stereochemistry(
        &mut topology,
        stereochem_done_marker_is_computed,
        params,
        &components,
    )?;
    let mut dative_donors = Vec::new();
    let mut hydrogen_bond_atoms = Vec::new();
    // RDKit✔️✔️: if (!doingCXSmiles || !includeStereoGroups) {
    // RDKit✔️✔️:   std::vector<StereoGroup> noStereoGroups;
    // RDKit✔️✔️:   tmol->setStereoGroups(noStereoGroups);
    // RDKit✔️✔️: }
    if doing_cx_smiles {
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
    if !doing_cx_smiles && !params.include_dative_bonds {
        // BEGIN RDKIT CPP FUNCTION detail::MolToSmiles includeDativeBonds conversion
        // RDKit❗✔️:     if (doingCXSmiles || !params.includeDativeBonds) {
        // RDKit❗✔️:       for (auto bond : tmol->bonds()) {
        // RDKit❗✔️:         if (bond->getBondType() == Bond::DATIVE) {
        // RDKit❗✔️:           bond->setBondType(Bond::SINGLE);
        // RDKit❗✔️:           bond->getBeginAtom()->calcExplicitValence(false);
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // END RDKIT CPP FUNCTION detail::MolToSmiles includeDativeBonds conversion
        for bond in &mut topology.bonds {
            if bond.order() == BondOrder::Dative {
                dative_donors.push(bond.begin());
                bond.set_order(BondOrder::Single);
            }
        }
    }
    if !doing_cx_smiles {
        topology.stereo_groups.clear();
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
    let mut component_ranks = Vec::with_capacity(components.len());
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles canonical rank policy
    // RDKit✔️✔️:       const bool breakTies = true;
    // RDKit✔️✔️:       const bool includeChiralPresence = false;
    // RDKit✔️✔️:       const bool includeIsotopes = params.doIsomericSmiles;
    // RDKit✔️✔️:       const bool includeChirality = params.doIsomericSmiles;
    // RDKit✔️✔️:       const bool includeStereoGroups = params.doIsomericSmiles;
    // RDKit✔️✔️:       const bool useNonStereoRanks = false;
    // RDKit✔️✔️:       const bool includeAtomMaps = true;
    // RDKit✔️✔️:       Canon::rankMolAtoms(*tmol, ranks, breakTies, includeChirality,
    // RDKit✔️✔️:                           includeIsotopes, includeAtomMaps,
    // RDKit✔️✔️:                           includeChiralPresence, includeStereoGroups,
    // RDKit✔️✔️:                           useNonStereoRanks);
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles canonical rank policy
    // Complexity review: policy propagation adds constant-time field assignments per
    // component. Fragment construction, rank refinement, allocations, loop nesting,
    // and component mapping are unchanged; no scan, clone, or lookup was added.
    let mut kekulize_fragments = if params.do_kekule {
        // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags fragment copies before ranking
        // RDKit❗❌:   if (nFrags == 1) {
        // RDKit❗❌:     res.emplace_back(new RWMol(mol));
        // RDKit❗❌:   } else {
        // RDKit❗❌:     if (comp.size() == 1 ||
        // RDKit❗❌:         (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
        // RDKit❗❌:       auto submol = copyMolSubset(mol, atoms, info, opts);
        // RDKit❗❌:       res.push_back(std::move(submol));
        // RDKit❗❌:     } else {
        // RDKit❗❌:       res.emplace_back(new RWMol(mol));
        // RDKit❗❌:       auto &frag = res.back();
        // RDKit❗❌:       frag->beginBatchEdit();
        // RDKit❗❌:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
        // RDKit❗❌:         if (!atomsInFrag[idx]) { frag->removeAtom(idx); }
        // RDKit❗❌:       }
        // RDKit❗❌:       frag->commitBatchEdit();
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // END RDKIT CPP FUNCTION MolOps::getTheFrags fragment copies before ranking
        // Local adapters keep source maps while allocating a membership mask and
        // copying each component topology before any source fragment is ranked.
        components
            .iter()
            .map(|component| {
                if components.len() == 1 {
                    Ok(WriterStereoFragment {
                        topology: topology.clone(),
                        source_atoms: (0..topology.atoms.len()).map(AtomId::new).collect(),
                        source_bonds: (0..topology.bonds.len()).map(BondId::new).collect(),
                    })
                } else {
                    let mut inside = vec![false; topology.atoms.len()];
                    for source_atom in component {
                        inside[*source_atom] = true;
                    }
                    let clears_computed_props = component.len() == 1
                        || (components.len() > 3
                            && !fragment_has_challenging_features(&topology, component, &inside));
                    if clears_computed_props {
                        extract_writer_subset_fragment(&topology, component, &inside)
                    } else {
                        extract_writer_preserving_fragment(&topology, &inside)
                    }
                }
            })
            .collect::<Result<Vec<_>, _>>()?
    } else {
        Vec::new()
    };
    for (component_index, component) in components.iter().enumerate() {
        component_ranks.push(if params.canonical {
            let mut policy = canonical_rank::CanonicalRankPolicy::default();
            policy.include_chirality = params.do_isomeric_smiles;
            policy.include_isotopes = params.do_isomeric_smiles;
            policy.include_stereo_groups = params.do_isomeric_smiles;
            canonical_rank::rank_component_atoms_with_policy(&ranking_topology, component, policy)
                .map_err(|error| SmilesParseError::CanonicalRank(error.to_string()))?
                .into_iter()
                .map(|rank| rank as i64)
                .collect::<Vec<_>>()
        } else {
            (0..topology.atoms.len())
                .map(|index| index as i64)
                .collect::<Vec<_>>()
        });
        // RDKit restores ignored map values after this fragment's rank and
        // before constructing that fragment's SMILES. Keep the source rows
        // mapped to their saved values while later component ranks continue to
        // use the unchanged, map-cleared `ranking_topology`.
        if let Some(original_atom_maps) = original_atom_maps.as_ref() {
            // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles map restoration
            // RDKit✔️✔️:       if (params.ignoreAtomMapNumbers) {
            // RDKit✔️✔️:         for (auto atom : tmol->atoms()) {
            // RDKit✔️✔️:           atom->setAtomMapNum(atomMapNums[atom->getIdx()]);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles map restoration
            // Source map number zero clears the map property on restoration.
            for source_atom in component {
                topology.atoms[*source_atom].set_atom_map(original_atom_maps[*source_atom]);
            }
        }
        if params.do_kekule {
            // BEGIN RDKIT CPP FUNCTION FragmentSmilesConstruct Kekulize branch
            // RDKit❗❌:   if (params.doKekule) {
            // RDKit❗❌:     if (atomsInPlay && bondsInPlay) {
            // RDKit❗❌:       MolOps::details::KekulizeFragment(static_cast<RWMol &>(mol), *atomsInPlay, *bondsInPlay);
            // RDKit❗❌:     } else {
            // RDKit❗❌:       MolOps::Kekulize(static_cast<RWMol &>(mol));
            // RDKit❗❌:     }
            // RDKit❗❌:   }
            // END RDKIT CPP FUNCTION FragmentSmilesConstruct Kekulize branch
            // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment atom error
            // RDKit✔️✔️:           throw AtomKekulizeException(msg, atom->getIdx());
            // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment atom error
            // FragmentSmilesConstruct receives an extracted component from
            // getTheFrags, so its full Kekulize call reports component-local
            // atom rows. Reuse the writer's source-mapped extraction branches
            // and project only core Kekulize's atom/bond fields back to the
            // private full-topology copy.
            let fragment = &mut kekulize_fragments[component_index];
            if let Some(original_atom_maps) = original_atom_maps.as_ref() {
                for (local_index, source_atom) in fragment.source_atoms.iter().copied().enumerate()
                {
                    fragment.topology.atoms[local_index]
                        .set_atom_map(original_atom_maps[source_atom.index()]);
                }
            }
            let kekulized =
                cosmolkit_core::kekulize(&fragment.topology, &KekulizeParams::default())
                    .map_err(SmilesParseError::WriterKekulize)?;
            for (local_index, source_atom) in fragment.source_atoms.iter().copied().enumerate() {
                let kekulized_atom = &kekulized.topology.atoms[local_index];
                let target = &mut topology.atoms[source_atom.index()];
                target.set_aromatic(kekulized_atom.is_aromatic());
                target.set_explicit_hydrogens(kekulized_atom.explicit_hydrogens());
                target.set_no_implicit(kekulized_atom.no_implicit());
            }
            for (local_index, source_bond) in fragment.source_bonds.iter().copied().enumerate() {
                let kekulized_bond = &kekulized.topology.bonds[local_index];
                let target = &mut topology.bonds[source_bond.index()];
                target.set_order(kekulized_bond.order());
                target.set_aromatic(kekulized_bond.is_aromatic());
            }
        }
    }
    if params.canonical && !params.do_isomeric_smiles {
        prepare_canonical_nonisomeric_stereo_fallback(
            &mut topology,
            stereochem_done_marker_is_computed,
            &components,
        )?;
    }
    let ring_bonds = find_ring_bonds(&topology);
    let mut colors = vec![AtomColor::White; topology.atoms.len()];
    let mut fragments = Vec::new();

    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles canonical and non-canonical fragments
    // RDKit✔️✔️:     // Sort the vfragsmi, but also sort the atom and bond order vectors into
    // RDKit✔️✔️:     // the same order
    // RDKit✔️✔️:     typedef std::tuple<std::string, std::vector<unsigned int>,
    // RDKit✔️✔️:                        std::vector<unsigned int>> tplType;
    // RDKit✔️✔️:     std::vector<tplType> tmp(vfragsmi.size());
    // RDKit✔️✔️:     for (unsigned int ti = 0; ti < vfragsmi.size(); ++ti) {
    // RDKit✔️✔️:       tmp[ti] = std::make_tuple(vfragsmi[ti], allAtomOrdering[ti],
    // RDKit✔️✔️:                                 allBondOrdering[ti]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::sort(tmp.begin(), tmp.end());
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
    for (component, ranks) in components.into_iter().zip(component_ranks) {
        // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles rooted component
        // RDKit✔️🔝:     rootedAtAtom = -1;
        // RDKit✔️🔝:     if (params.rootedAtAtom >= 0 && atsPresent[params.rootedAtAtom]) {
        // RDKit✔️🔝:       rootedAtAtom = params.rootedAtAtom - atsPresent.find_first();
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     if (rootedAtAtom >= 0) {
        // RDKit✔️🔝:       nextAtomIdx = rootedAtAtom;
        // RDKit✔️🔝:       rootedAtAtom = -1;
        // RDKit✔️🔝:     } else {
        // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles rooted component
        // Sorted component rows make root membership/remapping O(V) across the
        // writer; the source rebuilds an O(V)-bit set once per component.
        let start = if let Some(root) = params
            .rooted_at_atom
            .filter(|root| component.contains(&root.index()))
        {
            // RDKit subtracts the fragment's first source atom index, then
            // uses that value in the compact fragment's atom-index domain.
            let rooted_fragment_index = root.index() - component[0];
            *component.get(rooted_fragment_index).ok_or(
                SmilesParseError::WriterRootAtomOutOfRange {
                    atom_index: rooted_fragment_index,
                    atom_count: component.len(),
                },
            )?
        } else if params.canonical {
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
        let text = write_mol_stack(&topology, &valence, &stack, &chiral_adjustments, params)?;
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
        fragments.sort_by(|left, right| {
            left.text
                .cmp(&right.text)
                .then_with(|| left.atom_order.cmp(&right.atom_order))
                .then_with(|| left.bond_order.cmp(&right.bond_order))
        });
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

fn prepare_writer_stereochemistry(
    topology: &mut TopologyBlock,
    stereochem_done_marker_is_computed: Option<bool>,
    params: &SmilesWriteParams,
    components: &[Vec<usize>],
) -> Result<(), SmilesParseError> {
    // Keep malformed caller-supplied ring-relative state observable even when
    // clean stereo would clear computed ring relations before serialization.
    for atom in &topology.atoms {
        if let Some(encoded) = atom.prop("_ringStereoAtoms") {
            parse_ring_stereo_atoms(encoded, topology.atoms.len())?;
        }
    }
    if !params.do_isomeric_smiles {
        return Ok(());
    }

    // BEGIN RDKIT CPP FUNCTION MolOps::assignStereochemistry defaults
    // RDKit❗✔️: RDKIT_GRAPHMOL_EXPORT void assignStereochemistry(
    // RDKit❗✔️:     ROMol &mol, bool cleanIt = false, bool force = false,
    // RDKit❗✔️:     bool flagPossibleStereoCenters = false);
    // END RDKIT CPP FUNCTION MolOps::assignStereochemistry defaults
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles stereo preparation
    // RDKit❗✔️:     if (params.doIsomericSmiles) {
    // RDKit❗✔️:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles stereo preparation
    // Behavior review: `clean_it=false` uses the existing false/false wrapper;
    // `clean_it=true` currently uses true/true and fails the preserved pinned
    // modern large-ring case. No wrapper is selected by molecule contents.
    prepare_writer_stereo_components(
        topology,
        stereochem_done_marker_is_computed,
        params.clean_stereo,
        components,
    )
}

fn prepare_canonical_nonisomeric_stereo_fallback(
    topology: &mut TopologyBlock,
    stereochem_done_marker_is_computed: Option<bool>,
    components: &[Vec<usize>],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment stereo fallback
    // RDKit❗✔️:   if (!mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:     MolOps::assignStereochemistry(mol, false);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment stereo fallback
    // The caller reaches this only after the writer has ranked components and
    // performed FragmentSmilesConstruct's optional kekulization.
    // Complexity review: this reuses the writer's existing component mapping
    // and core assignment; the detached adapter still materializes fragment
    // topology and mapping vectors that RDKit already owns as mutable fragments.
    prepare_writer_stereo_components(
        topology,
        stereochem_done_marker_is_computed,
        false,
        components,
    )
}

fn prepare_writer_stereo_components(
    topology: &mut TopologyBlock,
    stereochem_done_marker_is_computed: Option<bool>,
    clean_it: bool,
    components: &[Vec<usize>],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags computed-property copy selection
    // RDKit❗❌:       if (comp.size() == 1 ||
    // RDKit❗❌:           (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
    // RDKit❗❌:         SubsetOptions opts{.sanitize = sanitizeFrags,
    // RDKit❗❌:                            .clearComputedProps = true,
    // RDKit❗❌:                            .copyCoordinates = copyConformers,
    // RDKit❗❌:                            .method = SubsetMethod::BONDS_BETWEEN_ATOMS};
    // RDKit❗❌:         auto submol = copyMolSubset(mol, atoms, info, opts);
    // RDKit❗❌:         res.push_back(std::move(submol));
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res.emplace_back(new RWMol(mol));
    // RDKit❗❌:         auto &frag = res.back();
    // RDKit❗❌:         frag->beginBatchEdit();
    // RDKit❗❌:         for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:           if (!atomsInFrag[idx]) {
    // RDKit❗❌:             frag->removeAtom(idx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         frag->commitBatchEdit();
    // RDKit❗❌:       }
    // END RDKIT CPP FUNCTION MolOps::getTheFrags computed-property copy selection
    // BEGIN RDKIT CPP FUNCTION Subset::copyMolSubset molecule-property transport
    // RDKit❗❌:   auto extracted_mol = std::make_unique<RWMol>();
    // RDKit❗❌:   if (options.clearComputedProps) {
    // RDKit❗❌:     extracted_mol->clearComputedProps();
    // RDKit❗❌:   } else {
    // RDKit❗❌:     copyComputedProps(mol, *extracted_mol);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION Subset::copyMolSubset molecule-property transport
    // BEGIN RDKIT CPP FUNCTION RWMol::commitBatchEdit computed-property clearing
    // RDKit❗❌:   batchRemoveBonds();
    // RDKit❗❌:   batchRemoveAtoms();
    // RDKit❗❌:   dp_ringInfo->reset();
    // RDKit❗❌:   clearComputedProps(true);
    // RDKit❗❌:   dp_delBonds.reset();
    // RDKit❗❌:   dp_delAtoms.reset();
    // END RDKIT CPP FUNCTION RWMol::commitBatchEdit computed-property clearing
    // BEGIN RDKIT CPP FUNCTION ROMol::clearComputedProps
    // RDKit❗❌:   RDProps::clearComputedProps();
    // RDKit❗❌:   for (auto atom : atoms()) {
    // RDKit❗❌:     atom->clearComputedProps();
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto bond : bonds()) {
    // RDKit❗❌:     bond->clearComputedProps();
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION ROMol::clearComputedProps
    // BEGIN RDKIT CPP FUNCTION RDProps::clearComputedProps
    // RDKit❗❌:   STR_VECT compLst;
    // RDKit❗❌:   if (getPropIfPresent(RDKit::detail::computedPropName, compLst) &&
    // RDKit❗❌:       !compLst.empty()) {
    // RDKit❗❌:     for (const auto &sv : compLst) {
    // RDKit❗❌:       d_props.clearVal(sv);
    // RDKit❗❌:     }
    // RDKit❗❌:     compLst.clear();
    // RDKit❗❌:     d_props.setVal(RDKit::detail::computedPropName, compLst);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION RDProps::clearComputedProps
    // Behavior review: the single-component clone preserves the marker;
    // subset extraction starts from an empty molecule, while clone/prune
    // clears computed but retains ordinary molecule properties at batch
    // commit. The source guard survives only when the extracted fragment still
    // has `_StereochemDone`.
    // Complexity review: both implementations scan the source graph per
    // fragment. The Rust preserving branch clones source and working topology,
    // then materializes mapped output rows, while RDKit clones once and prunes
    // in place; that branch has a clear extra full-topology copy.
    for (component_index, component) in components.iter().enumerate() {
        let mut inside = vec![false; topology.atoms.len()];
        for atom in component {
            inside[*atom] = true;
        }
        let clears_computed_props = components.len() > 1
            && (component.len() == 1
                || (components.len() > 3
                    && !fragment_has_challenging_features(topology, component, &inside)));
        let marker_survives_extraction = components.len() == 1
            || (!clears_computed_props && stereochem_done_marker_is_computed == Some(false));
        if stereochem_done_marker_is_computed.is_some() && marker_survives_extraction {
            continue;
        }

        let fragment = if components.len() == 1 {
            // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags one-component path
            // RDKit❗❌:   if (nFrags == 1) {
            // RDKit❗❌:     res.emplace_back(new RWMol(mol));
            // RDKit❗❌:     if (fragsMolAtomMapping) {
            // RDKit❗❌:       INT_VECT comp;
            // RDKit❗❌:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
            // RDKit❗❌:         comp.push_back(idx);
            // RDKit❗❌:       }
            // RDKit❗❌:       (*fragsMolAtomMapping).push_back(comp);
            // RDKit❗❌:     }
            // RDKit❗❌:   }
            // END RDKIT CPP FUNCTION MolOps::getTheFrags one-component path
            // The source deep-clones the sole component and maps each atom to
            // its unchanged source index. This topology clone also retains
            // ordinary/computed atom and bond properties; the writer's
            // detached adapter adds an atom/bond transport map and merge pass.
            WriterStereoFragment {
                topology: topology.clone(),
                source_atoms: (0..topology.atoms.len()).map(AtomId::new).collect(),
                source_bonds: (0..topology.bonds.len()).map(BondId::new).collect(),
            }
        } else if clears_computed_props {
            extract_writer_subset_fragment(topology, component, &inside)?
        } else {
            extract_writer_preserving_fragment(topology, &inside)?
        };

        let valence = cosmolkit_core::assign_valence_with_options_for_topology(
            &fragment.topology,
            ValenceModel::RdkitLike,
            false,
        )
        .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?;
        let rings = fast_find_rings_from_parts(
            fragment.topology.atoms.len(),
            &fragment.topology.bonds,
            &fragment.topology.adjacency,
        )
        .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        let assigned = if clean_it {
            // RDKit❗✔️: the writer passes cleanStereo while the omitted
            // flagPossibleStereoCenters argument defaults to false; the current
            // cleanup wrapper selects true for that flag. The serialized-output
            // regression below checks this wrapper independently.
            cosmolkit_core::assign_legacy_stereochemistry(fragment.topology, &valence, &rings)
        } else {
            cosmolkit_core::assign_legacy_stereochemistry_for_depiction(
                fragment.topology,
                &valence,
                &rings,
            )
        }
        .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;

        merge_writer_stereo_fragment(
            topology,
            &assigned,
            &fragment.source_atoms,
            &fragment.source_bonds,
        )?;

        // Keep the source's component walk visible in this local anchor. The
        // index is intentionally consumed by debug assertions rather than
        // changing the serialized traversal order.
        debug_assert!(component_index < components.len());
    }
    Ok(())
}

fn fragment_has_challenging_features(
    topology: &TopologyBlock,
    component: &[usize],
    inside: &[bool],
) -> bool {
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags fragmentHasChallengingFeatures
    // RDKit✔️✔️:       auto fragmentHasChallengingFeatures =
    // RDKit✔️✔️:           [&](const INT_VECT &comp,
    // RDKit✔️✔️:               const boost::dynamic_bitset<> &atomsInFrag) -> bool {
    // RDKit✔️✔️:         for (auto idx : comp) {
    // RDKit✔️✔️:           const auto atom = mol.getAtomWithIdx(idx);
    // RDKit✔️✔️:           if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:               atom->getChiralTag() != Atom::ChiralType::CHI_OTHER) {
    // RDKit✔️✔️:             return true;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           for (auto bnd : mol.atomBonds(atom)) {
    // RDKit✔️✔️:             if (atomsInFrag[bnd->getOtherAtomIdx(idx)]) {
    // RDKit✔️✔️:               if (bnd->getStereo() != Bond::BondStereo::STEREONONE &&
    // RDKit✔️✔️:                   bnd->getStereo() != Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:                 return true;
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         for (auto sgroup : getSubstanceGroups(mol)) {
    // RDKit✔️✔️:           for (auto aid : sgroup.getAtoms()) {
    // RDKit✔️✔️:             if (atomsInFrag[aid]) {
    // RDKit✔️✔️:               return true;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           for (auto aid : sgroup.getParentAtoms()) {
    // RDKit✔️✔️:             if (atomsInFrag[aid]) {
    // RDKit✔️✔️:               return true;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         for (auto stereoGroup : mol.getStereoGroups()) {
    // RDKit✔️✔️:           for (auto atom : stereoGroup.getAtoms()) {
    // RDKit✔️✔️:             if (atomsInFrag[atom->getIdx()]) {
    // RDKit✔️✔️:               return true;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           for (auto bond : stereoGroup.getBonds()) {
    // RDKit✔️✔️:             if (atomsInFrag[bond->getBeginAtomIdx()] &&
    // RDKit✔️✔️:                 atomsInFrag[bond->getEndAtomIdx()]) {
    // RDKit✔️✔️:               return true;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       };
    // END RDKIT CPP FUNCTION MolOps::getTheFrags fragmentHasChallengingFeatures
    // Behavior review: each modeled challenge category is checked on the
    // selected component. The source and Rust loops inspect atom tags,
    // internal bond stereo, SGroup atom/parent membership, and stereo-group
    // atom plus fully internal bond membership.
    // Complexity review: one component traversal and source-like scans of
    // incident bonds and group memberships; no nested whole-graph rescan is
    // added beyond the source lambda.
    for atom_index in component {
        let atom = &topology.atoms[*atom_index];
        if !matches!(atom.chiral_tag(), ChiralTag::Unspecified | ChiralTag::Other) {
            return true;
        }
        for neighbor in topology.adjacency.neighbors_of(*atom_index) {
            if inside[neighbor.atom_index] {
                let stereo = topology.bonds[neighbor.bond.index()].stereo();
                if !matches!(stereo, BondStereo::None | BondStereo::Any) {
                    return true;
                }
            }
        }
    }
    if topology.substance_groups.iter().any(|group| {
        group
            .atoms()
            .iter()
            .chain(group.parent_atoms())
            .any(|atom| inside[atom.index()])
    }) {
        return true;
    }
    topology.stereo_groups.iter().any(|group| {
        group.atoms().iter().any(|atom| inside[atom.index()])
            || group.bonds().iter().any(|bond| {
                let bond = &topology.bonds[bond.index()];
                inside[bond.begin().index()] && inside[bond.end().index()]
            })
    })
}

fn extract_writer_subset_fragment(
    topology: &TopologyBlock,
    component: &[usize],
    inside: &[bool],
) -> Result<WriterStereoFragment, SmilesParseError> {
    // Behavior review: for a multi-atom connected component, `component_bonds`
    // below contains every bond with both endpoints in `inside`; this is the
    // source `BONDS_BETWEEN_ATOMS` selection. `subtopology_from_path` creates
    // source-ordered rows, clears atom/bond computed properties, clears cis or
    // trans stereo when either stereo reference is unmappable, and returns
    // validated `new_to_old` maps. The singleton branch below supplies the
    // no-bond case directly and copies source-selected groups. The model has
    // additional typed SGroup references beyond RDKit's three selected index
    // lists, so those edge cases remain explicitly incomplete.
    // Complexity review: each component allocates a V-sized mask, scans all
    // source bonds to collect its internal edges, then the core subset path
    // scans atom and bond rows and allocates old/new maps. RDKit's
    // `copyMolSubset` consumes selected bitsets and scans its source rows once;
    // the local edge-inventory pass and mapping vectors are additional work.
    if component.len() == 1 {
        // BEGIN RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds singleton
        // RDKit❗❗:   for (const auto &ref_atom : reference_mol.atoms()) {
        // RDKit❗❗:     if (!selectedAtoms[ref_atom->getIdx()]) {
        // RDKit❗❗:       continue;
        // RDKit❗❗:     }
        // RDKit❗❗:     std::unique_ptr<Atom> extracted_atom{
        // RDKit❗❗:         options.copyAsQuery ? new QueryAtom(*ref_atom) : ref_atom->copy()};
        // RDKit❗❗:     extracted_atom->clearComputedProps();
        // RDKit❗❗:   }
        // END RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds singleton
        // BEGIN RDKIT CPP FUNCTION Subset::isSelectedSGroup/copySelectedSubstanceGroups
        // RDKit❗❗:   auto is_selected_component = [](auto &indices, auto &selection_test) {
        // RDKit❗❗:     return indices.empty() ||
        // RDKit❗❗:            std::all_of(indices.begin(), indices.end(), selection_test);
        // RDKit❗❗:   };
        // RDKit❗❗:   return is_selected_component(sgroup.getAtoms(), atom_test) &&
        // RDKit❗❗:          is_selected_component(sgroup.getBonds(), bond_test) &&
        // RDKit❗❗:          is_selected_component(sgroup.getParentAtoms(), atom_test);
        // RDKit❗❗:     if (!isSelectedSGroup(sgroup, selection_info)) {
        // RDKit❗❗:       continue;
        // RDKit❗❗:     }
        // RDKit❗❗:     SubstanceGroup extracted_sgroup(sgroup);
        // RDKit❗❗:     extracted_sgroup.setOwningMol(&extracted_mol);
        // RDKit❗❗:     update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getAtoms),
        // RDKit❗❗:                    std::mem_fn(&SubstanceGroup::setAtoms), atomMapping);
        // RDKit❗❗:     update_indices(extracted_sgroup,
        // RDKit❗❗:                    std::mem_fn(&SubstanceGroup::getParentAtoms),
        // RDKit❗❗:                    std::mem_fn(&SubstanceGroup::setParentAtoms), atomMapping);
        // RDKit❗❗:     update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getBonds),
        // RDKit❗❗:                    std::mem_fn(&SubstanceGroup::setBonds), bondMapping);
        // RDKit❗❗:     addSubstanceGroup(extracted_mol, std::move(extracted_sgroup));
        // END RDKIT CPP FUNCTION Subset::isSelectedSGroup/copySelectedSubstanceGroups
        // BEGIN RDKIT CPP FUNCTION Subset::copySelectedStereoGroups
        // RDKit❗❗:   auto is_selected_component = [](auto &objects, auto &selected_indices) {
        // RDKit❗❗:     return objects.empty() ||
        // RDKit❗❗:            std::any_of(objects.begin(), objects.end(), [&](auto &object) {
        // RDKit❗❗:              return selected_indices[object->getIdx()];
        // RDKit❗❗:            });
        // RDKit❗❗:   };
        // RDKit❗❗:   auto is_selected_stereo_group = [&](const auto &stereo_group) {
        // RDKit❗❗:     return is_selected_component(stereo_group.getAtoms(),
        // RDKit❗❗:                                  selection_info.selectedAtoms) &&
        // RDKit❗❗:            is_selected_component(stereo_group.getBonds(),
        // RDKit❗❗:                                  selection_info.selectedBonds);
        // RDKit❗❗:   };
        // RDKit❗❗:   if (!is_selected_stereo_group(stereo_group)) { continue; }
        // END RDKIT CPP FUNCTION Subset::copySelectedStereoGroups
        // BEGIN RDKIT CPP FUNCTION Subset::copyMolSubset computed properties
        // RDKit❗❗:   if (options.clearComputedProps) {
        // RDKit❗❗:     extracted_mol->clearComputedProps();
        // RDKit❗❗:   } else {
        // RDKit❗❗:     copyComputedProps(mol, *extracted_mol);
        // RDKit❗❗:   }
        // END RDKIT CPP FUNCTION Subset::copyMolSubset computed properties
        let source_atom = AtomId::new(component[0]);
        let mut atom_old_to_new = vec![None; topology.atoms.len()];
        atom_old_to_new[source_atom.index()] = Some(AtomId::new(0));
        let mut atom = topology.atoms[source_atom.index()]
            .clone()
            .with_id(AtomId::new(0));
        atom.clear_computed_props();
        atom.remap_template_attachment_order(&atom_old_to_new)
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        let stereo_groups = topology
            .stereo_groups
            .iter()
            .filter_map(|group| {
                let atom_side_selected = group.atoms().is_empty()
                    || group.atoms().iter().any(|member| inside[member.index()]);
                let bond_side_selected = group.bonds().is_empty();
                if !atom_side_selected || !bond_side_selected {
                    return None;
                }
                let atoms = group
                    .atoms()
                    .iter()
                    .filter(|member| inside[member.index()])
                    .map(|_| AtomId::new(0))
                    .collect();
                let local = StereoGroup::new(group.kind(), atoms, Vec::new());
                Some(match group.id() {
                    Some(id) => local.with_id(id),
                    None => local,
                })
            })
            .collect();
        let bond_old_to_new = vec![None; topology.bonds.len()];
        let mut selected_sgroups = topology
            .substance_groups
            .iter()
            .map(|group| {
                let source_members_selected = group
                    .atoms()
                    .iter()
                    .all(|member| atom_old_to_new[member.index()].is_some())
                    && group
                        .bonds()
                        .iter()
                        .all(|member| bond_old_to_new[member.index()].is_some())
                    && group
                        .parent_atoms()
                        .iter()
                        .all(|member| atom_old_to_new[member.index()].is_some());
                source_members_selected
                    && group.can_remap_without_parent(&atom_old_to_new, &bond_old_to_new)
            })
            .collect::<Vec<_>>();
        loop {
            let mut changed = false;
            for (index, group) in topology.substance_groups.iter().enumerate() {
                if selected_sgroups[index]
                    && group.parent().is_some_and(|parent| {
                        !selected_sgroups
                            .get(parent.index())
                            .copied()
                            .unwrap_or(false)
                    })
                {
                    selected_sgroups[index] = false;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        let mut sgroup_old_to_new = vec![None; topology.substance_groups.len()];
        let mut next_sgroup_index = 0;
        for (old_index, keep) in selected_sgroups.iter().copied().enumerate() {
            if keep {
                sgroup_old_to_new[old_index] =
                    Some(cosmolkit_model::SubstanceGroupId::new(next_sgroup_index));
                next_sgroup_index += 1;
            }
        }
        let substance_groups = topology
            .substance_groups
            .iter()
            .enumerate()
            .filter(|(index, _)| selected_sgroups[*index])
            .map(|(index, group)| {
                group
                    .remapped(
                        sgroup_old_to_new[index].expect("selected SGroup has a new id"),
                        &atom_old_to_new,
                        &bond_old_to_new,
                        &sgroup_old_to_new,
                    )
                    .ok_or_else(|| {
                        SmilesParseError::WriterStereo(
                            "selected singleton SGroup could not be remapped".into(),
                        )
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        // The model's typed SGroup state can reference attach/crossing atoms
        // and CState bonds outside the three source membership lists. Such a
        // group is excluded above because its references cannot survive this
        // subset. Source-selected memberships and all representable references
        // are remapped into the singleton row before topology validation.
        let fragment =
            TopologyBlock::try_from_parts(vec![atom], Vec::new(), substance_groups, stereo_groups)
                .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        return Ok(WriterStereoFragment {
            topology: fragment,
            source_atoms: vec![source_atom],
            source_bonds: Vec::new(),
        });
    }

    // BEGIN RDKIT CPP FUNCTION Subset::getSubsetInfo BONDS_BETWEEN_ATOMS
    // RDKit✔️❌:   if (options.method == SubsetMethod::BONDS_BETWEEN_ATOMS) {
    // RDKit✔️❌:     for (const auto &atom_idx : path) {
    // RDKit✔️❌:       if (atom_idx < num_atoms) {
    // RDKit✔️❌:         selectedAtoms.set(atom_idx);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     for (const auto &bond : mol.bonds()) {
    // RDKit✔️❌:       if (selectedAtoms[bond->getBeginAtomIdx()] &&
    // RDKit✔️❌:           selectedAtoms[bond->getEndAtomIdx()]) {
    // RDKit✔️❌:         selectedBonds.set(bond->getIdx());
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION Subset::getSubsetInfo BONDS_BETWEEN_ATOMS
    // BEGIN RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds atom old/new rows
    // RDKit✔️❌:   for (const auto &ref_atom : reference_mol.atoms()) {
    // RDKit✔️❌:     if (!selectedAtoms[ref_atom->getIdx()]) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     std::unique_ptr<Atom> extracted_atom{
    // RDKit✔️❌:         options.copyAsQuery ? new QueryAtom(*ref_atom) : ref_atom->copy()};
    // RDKit✔️❌:     extracted_atom->clearComputedProps();
    // RDKit✔️❌:     atomMapping[ref_atom->getIdx()] = extracted_mol.addAtom(
    // RDKit✔️❌:         extracted_atom.release(), updateLabel, takeOwnership);
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds atom old/new rows
    // BEGIN RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds bond old/new rows
    // RDKit✔️❌:     auto num_bonds =
    // RDKit✔️❌:         extracted_mol.addBond(extracted_bond.release(), takeOwnership);
    // RDKit✔️❌:     bondMapping[ref_bond->getIdx()] = num_bonds - 1;
    // END RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds bond old/new rows
    // BEGIN RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds mapped bonds
    // RDKit❗❗:   for (const auto &ref_bond : reference_mol.bonds()) {
    // RDKit❗❗:     if (!selectedBonds[ref_bond->getIdx()]) {
    // RDKit❗❗:       continue;
    // RDKit❗❗:     }
    // RDKit❗❗:     if (atomMapping.find(ref_bond->getBeginAtomIdx()) == atomMapping.end() ||
    // RDKit❗❗:         atomMapping.find(ref_bond->getEndAtomIdx()) == atomMapping.end()) {
    // RDKit❗❗:       throw ValueErrorException("copyMolSubset: subset bonds contain atoms not contained in subset atoms");
    // RDKit❗❗:     }
    // RDKit❗❗:     std::unique_ptr<Bond> extracted_bond{
    // RDKit❗❗:         options.copyAsQuery ? new QueryBond(*ref_bond) : ref_bond->copy()};
    // RDKit❗❗:     auto &atoms = extracted_bond->getStereoAtoms();
    // RDKit❗❗:     if (atoms.size() == 2) {
    // RDKit❗❗:       auto map1 = atomMapping.find(atoms[0]);
    // RDKit❗❗:       auto map2 = atomMapping.find(atoms[1]);
    // RDKit❗❗:       if (map1 != atomMapping.end() && map2 != atomMapping.end()) {
    // RDKit❗❗:         atoms[0] = map1->second;
    // RDKit❗❗:         atoms[1] = map2->second;
    // RDKit❗❗:       } else {
    // RDKit❗❗:         atoms.clear();  // We couldn't map the stereo atoms
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:     for (auto &atomidx : atoms) {
    // RDKit❗❗:       auto map = atomMapping.find(atomidx);
    // RDKit❗❗:       if (map != atomMapping.end()) {
    // RDKit❗❗:         atomidx = map->second;
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // END RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds mapped bonds
    let component_bonds = topology
        .bonds
        .iter()
        .filter(|bond| inside[bond.begin().index()] && inside[bond.end().index()])
        .map(Bond::id)
        .collect::<Vec<_>>();
    let result = cosmolkit_core::subtopology_from_path(
        topology,
        &component_bonds,
        &cosmolkit_core::SubtopologyParams::default(),
    )
    .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
    let cosmolkit_core::DetachedPathSubgraph::Concrete(fragment) = result.subgraph else {
        return Err(SmilesParseError::WriterStereo(
            "writer component extraction returned a query graph".into(),
        ));
    };
    Ok(WriterStereoFragment {
        topology: fragment,
        source_atoms: result
            .mapping
            .atoms
            .new_to_old
            .into_iter()
            .flatten()
            .collect(),
        source_bonds: result
            .mapping
            .bonds
            .new_to_old
            .into_iter()
            .flatten()
            .collect(),
    })
}

fn extract_writer_preserving_fragment(
    topology: &TopologyBlock,
    inside: &[bool],
) -> Result<WriterStereoFragment, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags clone-and-prune branch
    // RDKit❗❌:         res.emplace_back(new RWMol(mol));
    // RDKit❗❌:         auto &frag = res.back();
    // RDKit❗❌:         frag->beginBatchEdit();
    // RDKit❗❌:         for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:           if (!atomsInFrag[idx]) {
    // RDKit❗❌:             frag->removeAtom(idx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         frag->commitBatchEdit();
    // END RDKIT CPP FUNCTION MolOps::getTheFrags clone-and-prune branch
    // Behavior review: row order and atom/bond mappings follow source order;
    // ordinary atom/bond properties and typed groups survive the validated
    // mapping. RDKit's `commitBatchEdit` clears computed properties on the
    // molecule, atoms, and bonds, so this topology-only adapter clears the
    // computed atom/bond properties after `finish`. The source removes
    // nonmembers in ascending source atom order and remaps surviving rows.
    // Complexity review: clearing computed properties scans retained atom
    // and bond rows, matching the source's linear property-clear pass;
    // `begin_batch_edit` still clones both source and working blocks, then
    // `finish` creates remapped rows. These two full topology copies are a
    // confirmed local cost increase over RDKit's one-copy in-place branch.
    let mut edit = topology
        .begin_batch_edit()
        .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
    for (index, selected) in inside.iter().copied().enumerate() {
        if !selected {
            edit.remove_atom(AtomId::new(index))
                .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        }
    }
    let (mut fragment, mapping) = edit
        .finish()
        .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
    // BEGIN RDKIT CPP FUNCTION RWMol::commitBatchEdit computed-property clearing
    // RDKit❗❌:   // fix properties
    // RDKit❗❌:   clearComputedProps(true);
    // END RDKIT CPP FUNCTION RWMol::commitBatchEdit computed-property clearing
    // BEGIN RDKIT CPP FUNCTION ROMol::clearComputedProps
    // RDKit❗❌:   RDProps::clearComputedProps();
    // RDKit❗❌:   for (auto atom : atoms()) {
    // RDKit❗❌:     atom->clearComputedProps();
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto bond : bonds()) {
    // RDKit❗❌:     bond->clearComputedProps();
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION ROMol::clearComputedProps
    // TopologyBlock has no molecule-level property block or persisted ring
    // cache; those source-owned fields remain handled by their caller. Its
    // retained atom and bond rows carry the source property-clearing effect.
    for atom in &mut fragment.atoms {
        atom.clear_computed_props();
    }
    for bond in &mut fragment.bonds {
        bond.clear_computed_props();
    }
    Ok(WriterStereoFragment {
        topology: fragment,
        source_atoms: mapping.atoms.new_to_old.into_iter().flatten().collect(),
        source_bonds: mapping.bonds.new_to_old.into_iter().flatten().collect(),
    })
}

fn merge_writer_stereo_fragment(
    topology: &mut TopologyBlock,
    assigned: &TopologyBlock,
    source_atoms: &[AtomId],
    source_bonds: &[BondId],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles assignment producer
    // RDKit❗❌:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit❗❌:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit❗❌:       }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles assignment producer
    // Rust-only detached transport: RDKit mutates each temporary fragment in
    // place and has no merge function. The core assignment returns an owned
    // topology, so this adapter maps each assigned row through the extraction
    // `new_to_old` vectors: atoms carry chiral tag, explicit-H/no-implicit
    // state, and computed properties; bonds carry direction, stereo enum,
    // stereo-atom references, and computed properties. Ordinary target props
    // and the original typed groups stay on the full topology. `_ringStereoAtoms`
    // is the one computed property whose atom indices are fragment-relative;
    // signed one-based local IDs are converted to signed one-based source IDs.
    // D10's regression covers both signs and ordinary/computed property
    // preservation. The extraction mapping is validated before these rows are
    // returned by the subset/batch-edit owners.
    // Complexity review: source mutation has no copy-back pass. This adapter
    // scans each mapped atom/bond and clones computed property strings; that
    // is an extra O(V+E) transport and allocation cost.
    for (local_index, source_id) in source_atoms.iter().copied().enumerate() {
        let prepared = &assigned.atoms[local_index];
        let target = &mut topology.atoms[source_id.index()];
        target.set_chiral_tag(prepared.chiral_tag());
        target.set_explicit_hydrogens(prepared.explicit_hydrogens());
        target.set_no_implicit(prepared.no_implicit());
        target.clear_computed_props();
        for (key, value) in prepared.props() {
            if !prepared.is_prop_computed(key) {
                continue;
            }
            let value = if key == "_ringStereoAtoms" {
                parse_ring_stereo_atoms(value, assigned.atoms.len())?
                    .into_iter()
                    .map(|(same_orientation, local_atom)| {
                        let source_atom = source_atoms[local_atom].index() + 1;
                        let signed = if same_orientation {
                            i64::try_from(source_atom).unwrap_or(i64::MAX)
                        } else {
                            -i64::try_from(source_atom).unwrap_or(i64::MAX)
                        };
                        signed.to_string()
                    })
                    .collect::<Vec<_>>()
                    .join(",")
            } else {
                value.clone()
            };
            target
                .set_computed_prop(key.clone(), value)
                .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        }
    }
    for (local_index, source_id) in source_bonds.iter().copied().enumerate() {
        let prepared = &assigned.bonds[local_index];
        let target = &mut topology.bonds[source_id.index()];
        target.set_direction(prepared.direction());
        let stereo_atoms = prepared.stereo_atoms().map(|atoms| {
            [
                source_atoms[atoms[0].index()],
                source_atoms[atoms[1].index()],
            ]
        });
        target.set_stereo_atoms(stereo_atoms);
        target
            .set_stereo(prepared.stereo())
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        target.clear_computed_props();
        for (key, value) in prepared.props() {
            if prepared.is_prop_computed(key) {
                target
                    .set_computed_prop(key.clone(), value.clone())
                    .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
            }
        }
    }
    Ok(())
}

fn connected_components(topology: &TopologyBlock) -> Vec<Vec<usize>> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getMolFrags component mapping
    // RDKit✔️❌: unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
    // RDKit✔️❌:   unsigned int natms = mol.getNumAtoms();
    // RDKit✔️❌:   mapping.resize(natms);
    // RDKit✔️❌:   return natms ? boost::connected_components(mol.getTopology(), &mapping[0])
    // RDKit✔️❌:                    : 0;
    // RDKit✔️❌: };
    // RDKit✔️❌: unsigned int getMolFrags(const ROMol &mol, VECT_INT_VECT &frags) {
    // RDKit✔️❌:   frags.clear();
    // RDKit✔️❌:   INT_VECT mapping;
    // RDKit✔️❌:   getMolFrags(mol, mapping);
    // RDKit✔️❌:
    // RDKit✔️❌:   INT_INT_VECT_MAP comMap;
    // RDKit✔️❌:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit✔️❌:     int mi = mapping[i];
    // RDKit✔️❌:     if (comMap.find(mi) == comMap.end()) {
    // RDKit✔️❌:       INT_VECT comp;
    // RDKit✔️❌:       comMap[mi] = comp;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     comMap[mi].push_back(i);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (INT_INT_VECT_MAP_CI mci = comMap.begin(); mci != comMap.end(); mci++) {
    // RDKit✔️❌:     frags.push_back((*mci).second);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return rdcast<unsigned int>(frags.size());
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::getMolFrags component mapping
    // Complexity review: local DFS is O(V + E), but sorting each component adds
    // O(sum(|component| log |component|)); RDKit emits each component's atom rows
    // in source order while grouping the connected-component labels.
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
    // RDKit❗❌: const INT_LIST &trueOrder = atomTraversalBondOrder[atom->getIdx()];
    // RDKit❗❌: int nSwaps = 0;
    // RDKit❗❌: if (trueOrder.size() < atom->getDegree()) {
    // RDKit❗❌:   INT_LIST tOrder = trueOrder;
    // RDKit❗❌:   for (const auto bnd : mol.atomBonds(atom)) {
    // RDKit❗❌:     if (std::find(trueOrder.begin(), trueOrder.end(), bnd->getIdx()) ==
    // RDKit❗❌:         trueOrder.end()) {
    // RDKit❗❌:       tOrder.push_back(bnd->getIdx());
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   nSwaps = atom->getPerturbationOrder(tOrder);
    // RDKit❗❌: } else {
    // RDKit❗❌:   nSwaps = atom->getPerturbationOrder(trueOrder);
    // RDKit❗❌: }
    // RDKit❗❌: if (!perm) {
    // RDKit❗❌:   nSwaps = atom->getPerturbationOrder(tOrder);
    // RDKit❗❌: } else {
    // RDKit❗❌:   insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
    // RDKit❗❌:   perm = Chirality::getChiralPermutation(atom, tOrder);
    // RDKit❗❌: }
    // RDKit❗❌: if (doChiralInversions &&
    // RDKit❗❌:     chiralAtomNeedsTagInversion(mol, atom, firstInPart,
    // RDKit❗❌:                                 atomRingClosures[atom->getIdx()].size())) {
    // RDKit❗❌:   ++nSwaps;
    // RDKit❗❌: }
    // RDKit❗❌: if (nSwaps % 2) {
    // RDKit❗❌:   numSwapsChiralAtoms.set(atom->getIdx());
    // RDKit❗❌: }
    // RDKit❗❌: atomPermutationIndices[atom->getIdx()] = perm;
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
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment stack visit order
    // RDKit❗❌:   std::vector<unsigned int> atomVisitOrders(mol.getNumAtoms());
    // RDKit❗❌:   std::vector<unsigned int> bondVisitOrders(mol.getNumBonds());
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int pos = 0;
    // RDKit❗❌:   for (const auto &msI : molStack) {
    // RDKit❗❌:     if (msI.type == MOL_STACK_ATOM) {
    // RDKit❗❌:       atomVisitOrders[msI.obj.atom->getIdx()] = pos;
    // RDKit❗❌:     } else if (msI.type == MOL_STACK_BOND) {
    // RDKit❗❌:       bondVisitOrders[msI.obj.bond->getIdx()] = pos;
    // RDKit❗❌:       auto dir = msI.obj.bond->getBondDir();
    // RDKit❗❌:       if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
    // RDKit❗❌:         msI.obj.bond->setBondDir(Bond::NONE);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     ++pos;
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment stack visit order
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral post-processing section
    // RDKit❗❌: boost::dynamic_bitset<> ringStereoChemAdjusted(nAtoms);
    // RDKit❗❌: for (auto &msI : molStack) {
    // RDKit❗❌:   if (msI.type == MOL_STACK_ATOM &&
    // RDKit❗❌:       msI.obj.atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❌:       !msI.obj.atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗❌:     if (msI.obj.atom->hasProp(common_properties::_ringStereoAtoms)) {
    // RDKit❗❌:       if (!ringStereoChemAdjusted[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:         msI.obj.atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit❗❌:         ringStereoChemAdjusted.set(msI.obj.atom->getIdx());
    // RDKit❗❌:       }
    // RDKit❗❌:       const INT_VECT &ringStereoAtoms = msI.obj.atom->getProp<INT_VECT>(
    // RDKit❗❌:           common_properties::_ringStereoAtoms);
    // RDKit❗❌:       for (auto nbrV : ringStereoAtoms) {
    // RDKit❗❌:         int nbrIdx = abs(nbrV) - 1;
    // RDKit❗❌:         if (!ringStereoChemAdjusted[nbrIdx] &&
    // RDKit❗❌:             atomVisitOrders[nbrIdx] >
    // RDKit❗❌:                 atomVisitOrders[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:           mol.getAtomWithIdx(nbrIdx)->setChiralTag(
    // RDKit❗❌:               msI.obj.atom->getChiralTag());
    // RDKit❗❌:           if (nbrV < 0) {
    // RDKit❗❌:             mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit❗❌:           }
    // RDKit❗❌:           if (numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:             if (!numSwapsChiralAtoms[nbrIdx]) {
    // RDKit❗❌:               mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit❗❌:             }
    // RDKit❗❌:           } else {
    // RDKit❗❌:             if (numSwapsChiralAtoms[nbrIdx]) {
    // RDKit❗❌:               mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:           ringStereoChemAdjusted.set(nbrIdx);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     } else if (size_t sgidx;
    // RDKit❗❌:                msI.obj.atom->getPropIfPresent("_stereoGroup", sgidx) &&
    // RDKit❗❌:                mol.getStereoGroups().size() > sgidx) {
    // RDKit❗❌:       auto &sg = mol.getStereoGroups()[sgidx];
    // RDKit❗❌:       bool swapIt =
    // RDKit❗❌:           msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW;
    // RDKit❗❌:       if (swapIt) {
    // RDKit❗❌:         msI.obj.atom->invertChirality();
    // RDKit❗❌:       }
    // RDKit❗❌:       if (swapIt || numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:         for (auto at : sg.getAtoms()) {
    // RDKit❗❌:           if (at == msI.obj.atom) {
    // RDKit❗❌:             continue;
    // RDKit❗❌:           }
    // RDKit❗❌:           at->invertChirality();
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       if (msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗❌:           msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗❌:         if ((numSwapsChiralAtoms[msI.obj.atom->getIdx()])) {
    // RDKit❗❌:           msI.obj.atom->invertChirality();
    // RDKit❗❌:         }
    // RDKit❗❌:       } else if (atomPermutationIndices[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:         msI.obj.atom->setProp(
    // RDKit❗❌:             common_properties::_chiralPermutation,
    // RDKit❗❌:             atomPermutationIndices[msI.obj.atom->getIdx()]);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral post-processing section
    // This detached path stores stack atom IDs again, uses byte-per-entry Vec<bool>,
    // and builds a BTreeMap for stereo-group references; RDKit iterates molStack
    // directly and uses a compact dynamic bitset. Its per-component storage is larger.

    // Source value-initializes absent atomVisitOrders to zero. Keep an explicit
    // sentinel and reject absent component members before the strict later test.
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
        if atom.chiral_tag() == ChiralTag::Unspecified || atom.prop("_brokenChirality").is_some() {
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
                    || atom_visit_order[neighbor_index] == usize::MAX
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
    // BEGIN RDKIT CPP TYPE RDGeneral/types.h INT_VECT
    // RDKit❗❌: typedef std::vector<int> INT_VECT;
    // END RDKIT CPP TYPE RDGeneral/types.h INT_VECT
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment ring-relative references
    // RDKit❗❌:           const INT_VECT &ringStereoAtoms = atom->getProp<INT_VECT>(
    // RDKit❗❌:               common_properties::_ringStereoAtoms);
    // RDKit❗❌:           for (auto nbrV : ringStereoAtoms) {
    // RDKit❗❌:             int nbrIdx = abs(nbrV) - 1;
    // RDKit❗❌:               if (nbrV < 0) {
    // RDKit❗❌:                 mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment ring-relative references
    // The detached model stores this computed property as text; decode only
    // the source's signed, one-based int entries and preserve their order.
    let mut result = Vec::new();
    for token in encoded.split(',') {
        let value = token.parse::<i32>().map_err(|_| {
            SmilesParseError::WriterStereo(
                "`_ringStereoAtoms` is not a signed source INT_VECT value (bad_any_cast)".into(),
            )
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
    params: &SmilesWriteParams,
) -> Result<String, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    // RDKit✔️❌:   for (auto &mSE : molStack) {
    // RDKit✔️❌:     switch (mSE.type) {
    // RDKit✔️❌:       case Canon::MOL_STACK_ATOM:
    // RDKit✔️❌:         for (auto rclosure : ringClosuresToErase) {
    // RDKit✔️❌:           ringClosureMap.erase(rclosure);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         ringClosuresToErase.clear();
    // RDKit✔️❌:         res << GetAtomSmiles(mSE.obj.atom, params);
    // RDKit✔️❌:         atomOrdering.push_back(mSE.obj.atom->getIdx());
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Canon::MOL_STACK_BOND:
    // RDKit✔️❌:         bond = mSE.obj.bond;
    // RDKit✔️❌:         res << GetBondSmiles(bond, params, mSE.number);
    // RDKit✔️❌:         bondOrdering.push_back(bond->getIdx());
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Canon::MOL_STACK_RING:
    // RDKit✔️❌:         ringIdx = mSE.number;
    // RDKit✔️❌:         if (ringClosureMap.count(ringIdx)) {
    // RDKit✔️❌:           closureVal = ringClosureMap[ringIdx];
    // RDKit✔️❌:           ringClosuresToErase.push_back(ringIdx);
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           closureVal = 1;
    // RDKit✔️❌:           bool done = false;
    // RDKit✔️❌:           while (!done) {
    // RDKit✔️❌:             std::map<int, int>::iterator mapIt;
    // RDKit✔️❌:             for (mapIt = ringClosureMap.begin();
    // RDKit✔️❌:                  mapIt != ringClosureMap.end(); ++mapIt) {
    // RDKit✔️❌:               if (mapIt->second == closureVal) {
    // RDKit✔️❌:                 break;
    // RDKit✔️❌:               }
    // RDKit✔️❌:             }
    // RDKit✔️❌:             if (mapIt == ringClosureMap.end()) {
    // RDKit✔️❌:               done = true;
    // RDKit✔️❌:             } else {
    // RDKit✔️❌:               closureVal += 1;
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    // RDKit✔️❌:           ringClosureMap[ringIdx] = closureVal;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (closureVal < 10) {
    // RDKit✔️❌:           res << (char)(closureVal + '0');
    // RDKit✔️❌:         } else if (closureVal < 100) {
    // RDKit✔️❌:           res << '%' << closureVal;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res << "%(" << closureVal << ')';
    // RDKit✔️❌:         }
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Canon::MOL_STACK_BRANCH_OPEN:
    // RDKit✔️❌:         res << "(";
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Canon::MOL_STACK_BRANCH_CLOSE:
    // RDKit✔️❌:         res << ")";
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       default:
    // RDKit✔️❌:         break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    // `write_smiles_output` collects atom and bond output-order rows in two
    // separate scans after this emission pass; RDKit records both in this loop.
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
                    params,
                )?);
            }
            MolStackElem::Bond { bond, atom_to_left } => output.push_str(&bond_text(
                topology,
                &topology.bonds[bond.index()],
                atom_to_left,
                params,
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
    params: &SmilesWriteParams,
) -> Result<String, SmilesParseError> {
    let atom = &topology.atoms[atom_index];
    // BEGIN RDKIT CPP FUNCTION atomNeedsBracket
    // RDKit❗✔️: bool atomNeedsBracket(const Atom *atom, const std::string &atString,
    // RDKit❗✔️:                       const SmilesWriteParams &params) {
    // RDKit❗✔️:   PRECONDITION(atom, "null atom");
    // RDKit❗✔️:   auto num = atom->getAtomicNum();
    // RDKit❗✔️:   if (!inOrganicSubset(num)) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (atom->getFormalCharge()) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (params.doIsomericSmiles && (atom->getIsotope() || !atString.empty())) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (atom->hasProp(common_properties::molAtomMapNumber)) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const INT_VECT &defaultVs = PeriodicTable::getTable()->getValenceList(num);
    // RDKit❗✔️:   int totalValence = atom->getTotalValence();
    // RDKit❗✔️:   bool nonStandard = false;
    // RDKit❗✔️:   if (atom->getNumRadicalElectrons()) {
    // RDKit❗✔️:     nonStandard = true;
    // RDKit❗✔️:   } else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
    // RDKit❗✔️:              atom->getNumExplicitHs()) {
    // RDKit❗✔️:     nonStandard = true;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     nonStandard = (totalValence != defaultVs.front() && atom->getTotalNumHs());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (nonStandard) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (atom->hasOwningMol()) {
    // RDKit❗✔️:     for (const auto bond : atom->getOwningMol().atomBonds(atom)) {
    // RDKit❗✔️:       auto oatom = bond->getOtherAtom(atom);
    // RDKit❗✔️:       if (QueryOps::isMetal(*oatom)) {
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION atomNeedsBracket
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
    // RDKit❗✔️:     if (fc > 0) {
    // RDKit❗✔️:       res += "+";
    // RDKit❗✔️:       if (fc > 1) {
    // RDKit❗✔️:         res += std::to_string(fc);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (fc < 0) {
    // RDKit❗✔️:       if (fc < -1) {
    // RDKit❗✔️:         res += std::to_string(fc);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         res += "-";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetAtomSmiles modeled non-stereo atom emission
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles custom symbol selection and aromatic spelling
    // RDKit❗✔️:   std::string symb;
    // RDKit❗✔️:   bool hasCustomSymbol =
    // RDKit❗✔️:       atom->getPropIfPresent(common_properties::smilesSymbol, symb);
    // RDKit❗✔️:   if (!hasCustomSymbol) {
    // RDKit❗✔️:     symb = PeriodicTable::getTable()->getElementSymbol(num);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // this was originally only done for the organic subset,
    // RDKit❗✔️:   // applying it to other atom-types is a fix for Issue 3152751:
    // RDKit❗✔️:   // Only accept for atom->getAtomicNum() in [5, 6, 7, 8, 14, 15, 16, 33, 34,
    // RDKit❗✔️:   // 52]
    // RDKit❗✔️:   if (!params.doKekule && atom->getIsAromatic() && symb[0] >= 'A' &&
    // RDKit❗✔️:       symb[0] <= 'Z') {
    // RDKit❗✔️:     switch (atom->getAtomicNum()) {
    // RDKit❗✔️:       case 5:
    // RDKit❗✔️:       case 6:
    // RDKit❗✔️:       case 7:
    // RDKit❗✔️:       case 8:
    // RDKit❗✔️:       case 14:
    // RDKit❗✔️:       case 15:
    // RDKit❗✔️:       case 16:
    // RDKit❗✔️:       case 33:
    // RDKit❗✔️:       case 34:
    // RDKit❗✔️:       case 52:
    // RDKit❗✔️:         symb[0] -= ('A' - 'a');
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetAtomSmiles custom symbol selection and aromatic spelling
    let custom_symbol = atom.prop("smilesSymbol");
    let raw_symbol = custom_symbol.unwrap_or_else(|| atom.element().symbol());
    let symbol = if !params.do_kekule
        && atom.is_aromatic()
        && raw_symbol
            .as_bytes()
            .first()
            .is_some_and(u8::is_ascii_uppercase)
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
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles chirality selection
    // RDKit❗✔️:   if (params.doIsomericSmiles) {
    // RDKit❗✔️:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:         !atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗✔️:       atString = getAtomChiralityInfo(atom);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetAtomSmiles chirality selection
    let mut chirality = String::new();
    if params.do_isomeric_smiles && atom.prop("_brokenChirality").is_none() {
        let base_chiral_tag = chiral_adjustment
            .chiral_tag_override
            .unwrap_or(atom.chiral_tag());
        let chiral_tag = if chiral_adjustment.invert_tetrahedral {
            stereo::invert_tetrahedral_tag(base_chiral_tag)
        } else {
            base_chiral_tag
        };
        chirality = match chiral_tag {
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
    }
    let needs_bracket = if custom_symbol.is_some() || params.all_hydrogens_explicit {
        true
    } else {
        if !in_organic_subset(atom.atomic_number()) {
            true
        } else if atom.formal_charge() != 0 {
            true
        } else if params.do_isomeric_smiles && (atom.isotope().is_some() || !chirality.is_empty()) {
            true
        } else if atom.atom_map().is_some() {
            true
        } else {
            let default_valences = cosmolkit_core::rdkit_valence_list(atom.atomic_number())
                .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?
                .ok_or_else(|| {
                    SmilesParseError::WriterValence(format!(
                        "RDKit valence list is unavailable for atomic number {}",
                        atom.atomic_number()
                    ))
                })?;
            let default_valence = default_valences.first().copied().ok_or_else(|| {
                SmilesParseError::WriterValence(format!(
                    "RDKit valence list is empty for atomic number {}",
                    atom.atomic_number()
                ))
            })?;
            let total_valence = valence.explicit_valence[atom_index]
                + valence.implicit_hydrogens[atom_index].max(0);
            let total_num_hydrogens = usize::from(atom.explicit_hydrogens())
                + usize::try_from(valence.implicit_hydrogens[atom_index].max(0))
                    .unwrap_or(usize::MAX);
            let nonstandard_valence = if atom.radical_electrons() != 0
                || (matches!(atom.atomic_number(), 7 | 15)
                    && atom.is_aromatic()
                    && atom.explicit_hydrogens() > 0)
            {
                true
            } else {
                total_valence != default_valence && total_num_hydrogens > 0
            };
            nonstandard_valence
                || topology
                    .adjacency
                    .neighbors_of(atom_index)
                    .iter()
                    .any(|neighbor| {
                        rdkit_query_ops_is_metal(
                            topology.atoms[neighbor.atom_index].atomic_number(),
                        )
                    })
        }
    };
    if !needs_bracket {
        return Ok(append_supplemental_label(atom, symbol));
    }

    let mut output = String::from("[");
    if params.do_isomeric_smiles {
        if let Some(isotope) = atom.isotope() {
            output.push_str(&isotope.to_string());
        }
    }
    output.push_str(&symbol);
    output.push_str(&chirality);
    let total_num_hydrogens = usize::from(atom.explicit_hydrogens())
        + usize::try_from(valence.implicit_hydrogens[atom_index].max(0)).unwrap_or(usize::MAX);
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
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles supplemental label
    // RDKit❗✔️:   // If the atom has this property, the contained string will
    // RDKit❗✔️:   // be inserted directly in the SMILES:
    // RDKit❗✔️:   std::string label;
    // RDKit❗✔️:   if (atom->getPropIfPresent(common_properties::_supplementalSmilesLabel,
    // RDKit❗✔️:                              label)) {
    // RDKit❗✔️:     res += label;
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetAtomSmiles supplemental label
    if let Some(label) = atom.prop("_supplementalSmilesLabel") {
        text.push_str(label);
    }
    text
}

fn bond_text(
    topology: &TopologyBlock,
    bond: &Bond,
    atom_to_left: usize,
    params: &SmilesWriteParams,
) -> Result<&'static str, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION GetBondSmiles aromatic-context selection
    // RDKit❗✔️:   bool aromatic = false;
    // RDKit❗✔️:   if (!params.doKekule && (bond->getBondType() == Bond::SINGLE ||
    // RDKit❗✔️:                            bond->getBondType() == Bond::DOUBLE ||
    // RDKit❗✔️:                            bond->getBondType() == Bond::AROMATIC)) {
    // RDKit❗✔️:     if (bond->hasOwningMol()) {
    // RDKit❗✔️:       auto a1 = bond->getOwningMol().getAtomWithIdx(atomToLeftIdx);
    // RDKit❗✔️:       auto a2 = bond->getOwningMol().getAtomWithIdx(
    // RDKit❗✔️:           bond->getOtherAtomIdx(atomToLeftIdx));
    // RDKit❗✔️:       if ((a1->getIsAromatic() && a2->getIsAromatic()) &&
    // RDKit❗✔️:           (a1->getAtomicNum() || a2->getAtomicNum())) {
    // RDKit❗✔️:         aromatic = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       aromatic = false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetBondSmiles aromatic-context selection
    // BEGIN RDKIT CPP FUNCTION GetBondSmiles direction selection
    // RDKit❗✔️:   Bond::BondDir dir = bond->getBondDir();
    // RDKit❗✔️:   switch (bond->getBondType()) {
    // RDKit❗✔️:     case Bond::SINGLE:
    // RDKit❗✔️:       if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
    // RDKit❗✔️:         switch (dir) {
    // RDKit❗✔️:           case Bond::ENDDOWNRIGHT:
    // RDKit❗✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
    // RDKit❗✔️:               res = "\\";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           case Bond::ENDUPRIGHT:
    // RDKit❗✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
    // RDKit❗✔️:               res = "/";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           default:
    // RDKit❗✔️:             if (params.allBondsExplicit) {
    // RDKit❗✔️:               res = "-";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         if (params.allBondsExplicit) {
    // RDKit❗✔️:           res = "-";
    // RDKit❗✔️:         } else if (aromatic && !bond->getIsAromatic()) {
    // RDKit❗✔️:           res = "-";
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::DOUBLE:
    // RDKit❗✔️:       if (!aromatic || !bond->getIsAromatic() || params.allBondsExplicit) {
    // RDKit❗✔️:         res = "=";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::TRIPLE:
    // RDKit❗✔️:       res = "#";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::QUADRUPLE:
    // RDKit❗✔️:       res = "$";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::AROMATIC:
    // RDKit❗✔️:       if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
    // RDKit❗✔️:         switch (dir) {
    // RDKit❗✔️:           case Bond::ENDDOWNRIGHT:
    // RDKit❗✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
    // RDKit❗✔️:               res = "\\";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           case Bond::ENDUPRIGHT:
    // RDKit❗✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
    // RDKit❗✔️:               res = "/";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           default:
    // RDKit❗✔️:             if (params.allBondsExplicit || !aromatic) {
    // RDKit❗✔️:               res = ":";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else if (params.allBondsExplicit || !aromatic) {
    // RDKit❗✔️:         res = ":";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::DATIVE:
    // RDKit❗✔️:       if (atomToLeftIdx >= 0 &&
    // RDKit❗✔️:           bond->getBeginAtomIdx() == static_cast<unsigned int>(atomToLeftIdx)) {
    // RDKit❗✔️:         res = "->";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         res = "<-";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       res = "~";
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetBondSmiles direction selection
    let other = other_atom(bond, atom_to_left)?;
    let aromatic_context = !params.do_kekule
        && matches!(
            bond.order(),
            BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
        )
        && topology.atoms[atom_to_left].is_aromatic()
        && topology.atoms[other].is_aromatic()
        && (topology.atoms[atom_to_left].atomic_number() != 0
            || topology.atoms[other].atomic_number() != 0);
    let direction = bond.direction();
    let direction_is_specified = !matches!(direction, BondDirection::None | BondDirection::Unknown);
    match bond.order() {
        BondOrder::Single => {
            if direction_is_specified {
                Ok(match direction {
                    BondDirection::EndDownRight
                        if params.all_bonds_explicit || params.do_isomeric_smiles =>
                    {
                        "\\"
                    }
                    BondDirection::EndUpRight
                        if params.all_bonds_explicit || params.do_isomeric_smiles =>
                    {
                        "/"
                    }
                    _ if params.all_bonds_explicit => "-",
                    _ => "",
                })
            } else if params.all_bonds_explicit || (aromatic_context && !bond.is_aromatic()) {
                Ok("-")
            } else {
                Ok("")
            }
        }
        BondOrder::Double => Ok(
            if !aromatic_context || !bond.is_aromatic() || params.all_bonds_explicit {
                "="
            } else {
                ""
            },
        ),
        BondOrder::Triple => Ok("#"),
        BondOrder::Quadruple => Ok("$"),
        BondOrder::Aromatic => {
            if direction_is_specified {
                Ok(match direction {
                    BondDirection::EndDownRight
                        if params.all_bonds_explicit || params.do_isomeric_smiles =>
                    {
                        "\\"
                    }
                    BondDirection::EndUpRight
                        if params.all_bonds_explicit || params.do_isomeric_smiles =>
                    {
                        "/"
                    }
                    _ if params.all_bonds_explicit || !aromatic_context => ":",
                    _ => "",
                })
            } else {
                Ok(if params.all_bonds_explicit || !aromatic_context {
                    ":"
                } else {
                    ""
                })
            }
        }
        BondOrder::Dative if bond.begin().index() == atom_to_left => Ok("->"),
        BondOrder::Dative => Ok("<-"),
        _ => Ok("~"),
    }
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
    // RDKit✔️❌:         if (closureVal < 10) {
    // RDKit✔️❌:           res << (char)(closureVal + '0');
    // RDKit✔️❌:         } else if (closureVal < 100) {
    // RDKit✔️❌:           res << '%' << closureVal;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res << "%(" << closureVal << ')';
    // RDKit✔️❌:         }
    // The multi-digit branches allocate a temporary with `to_string`; RDKit
    // streams the integer into its existing output buffer.
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
    // BEGIN RDKIT CPP FUNCTION RingInfo::numBondRings
    // RDKit✔️🔝: unsigned int RingInfo::numBondRings(unsigned int idx) const {
    // RDKit✔️🔝:   PRECONDITION(df_init, "RingInfo not initialized");
    // RDKit✔️🔝:
    // RDKit✔️🔝:   if (idx < d_bondMembers.size()) {
    // RDKit✔️🔝:     return rdcast<unsigned int>(d_bondMembers[idx].size());
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return 0;
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION RingInfo::numBondRings
    // END RDKIT CPP FUNCTION findSSSR active bond selection
    // A cycle basis spans the source-selected active graph's cycle space, so
    // every bond in a cycle has a stored ring membership. Tarjan's bridge test
    // computes that same `numBondRings() != 0` predicate directly in O(V+E),
    // avoiding SSSR ring-vector construction without changing the result.
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
        write_smiles_with_params(
            &record,
            &SmilesWriteParams {
                canonical: false,
                ..Default::default()
            },
        )
        .expect("write")
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
                write_smiles_with_params(
                    &record,
                    &SmilesWriteParams {
                        canonical: false,
                        ..Default::default()
                    }
                )
                .unwrap(),
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

    #[test]
    fn writer_subset_and_clone_prune_preserve_source_property_and_group_rules() {
        let mut record =
            parse_smiles("C[C@H](F)Cl.C[C@H](O)Br", &SmilesParseParams::default()).unwrap();
        record.topology.atoms[1]
            .set_prop("ordinary_atom", "kept")
            .unwrap();
        record.topology.atoms[1]
            .set_computed_prop("_computed_atom", "cleared")
            .unwrap();
        record.topology.bonds[0]
            .set_prop("ordinary_bond", "kept")
            .unwrap();
        record.topology.bonds[0]
            .set_computed_prop("_computed_bond", "cleared")
            .unwrap();
        record.topology.stereo_groups.push(
            StereoGroup::new(
                cosmolkit_model::StereoGroupKind::Or,
                vec![AtomId::new(1), AtomId::new(5)],
                vec![BondId::new(0), BondId::new(3)],
            )
            .with_id(17),
        );
        let before = record.clone();
        let inside = [true, true, true, true, false, false, false, false];

        let subset =
            extract_writer_subset_fragment(&record.topology, &[0, 1, 2, 3], &inside).unwrap();
        assert_eq!(
            subset
                .source_atoms
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>(),
            [0, 1, 2, 3]
        );
        assert_eq!(
            subset
                .source_bonds
                .iter()
                .map(|bond| bond.index())
                .collect::<Vec<_>>(),
            [0, 1, 2]
        );
        assert_eq!(subset.topology.atoms[1].prop("ordinary_atom"), Some("kept"));
        assert_eq!(subset.topology.atoms[1].prop("_computed_atom"), None);
        assert_eq!(subset.topology.bonds[0].prop("ordinary_bond"), Some("kept"));
        assert_eq!(subset.topology.bonds[0].prop("_computed_bond"), None);
        assert_eq!(subset.topology.stereo_groups.len(), 1);
        assert_eq!(subset.topology.stereo_groups[0].id(), Some(17));
        assert_eq!(
            subset.topology.stereo_groups[0].kind(),
            cosmolkit_model::StereoGroupKind::Or
        );
        assert_eq!(subset.topology.stereo_groups[0].atoms(), &[AtomId::new(1)]);
        assert_eq!(subset.topology.stereo_groups[0].bonds(), &[BondId::new(0)]);

        let clone_prune = extract_writer_preserving_fragment(&record.topology, &inside).unwrap();
        assert_eq!(
            clone_prune
                .source_atoms
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>(),
            [0, 1, 2, 3]
        );
        assert_eq!(
            clone_prune
                .source_bonds
                .iter()
                .map(|bond| bond.index())
                .collect::<Vec<_>>(),
            [0, 1, 2]
        );
        assert_eq!(
            clone_prune.topology.atoms[1].prop("ordinary_atom"),
            Some("kept")
        );
        assert_eq!(clone_prune.topology.atoms[1].prop("_computed_atom"), None);
        assert!(!clone_prune.topology.atoms[1].is_prop_computed("_computed_atom"));
        assert_eq!(
            clone_prune.topology.bonds[0].prop("ordinary_bond"),
            Some("kept")
        );
        assert_eq!(clone_prune.topology.bonds[0].prop("_computed_bond"), None);
        assert!(!clone_prune.topology.bonds[0].is_prop_computed("_computed_bond"));
        assert_eq!(clone_prune.topology.stereo_groups.len(), 1);
        assert_eq!(clone_prune.topology.stereo_groups[0].id(), Some(17));
        assert_eq!(
            clone_prune.topology.stereo_groups[0].kind(),
            cosmolkit_model::StereoGroupKind::Or
        );
        assert_eq!(
            clone_prune.topology.stereo_groups[0].atoms(),
            &[AtomId::new(1)]
        );
        assert_eq!(
            clone_prune.topology.stereo_groups[0].bonds(),
            &[BondId::new(0)]
        );
        assert_eq!(
            record, before,
            "fragment extraction leaves its source unchanged"
        );
    }

    #[test]
    fn writer_singleton_subset_clears_computed_props_and_keeps_partial_group_membership() {
        let mut record = parse_smiles("C.CC", &SmilesParseParams::default()).unwrap();
        record.topology.atoms[0]
            .set_prop("ordinary_atom", "kept")
            .unwrap();
        record.topology.atoms[0]
            .set_computed_prop("_computed_atom", "cleared")
            .unwrap();
        record.topology.substance_groups.push(
            cosmolkit_model::SubstanceGroup::new(
                cosmolkit_model::SubstanceGroupId::new(0),
                cosmolkit_model::SubstanceGroupKind::Data,
            )
            .with_atoms(vec![AtomId::new(0)])
            .with_label("singleton"),
        );
        record.topology.stereo_groups.push(
            StereoGroup::new(
                cosmolkit_model::StereoGroupKind::Or,
                vec![AtomId::new(0), AtomId::new(1)],
                Vec::new(),
            )
            .with_id(23),
        );
        let singleton =
            extract_writer_subset_fragment(&record.topology, &[0], &[true, false, false]).unwrap();

        assert_eq!(singleton.source_atoms, [AtomId::new(0)]);
        assert!(singleton.source_bonds.is_empty());
        assert_eq!(
            singleton.topology.atoms[0].prop("ordinary_atom"),
            Some("kept")
        );
        assert_eq!(singleton.topology.atoms[0].prop("_computed_atom"), None);
        assert_eq!(singleton.topology.substance_groups.len(), 1);
        assert_eq!(
            singleton.topology.substance_groups[0].id(),
            cosmolkit_model::SubstanceGroupId::new(0)
        );
        assert_eq!(
            singleton.topology.substance_groups[0].atoms(),
            &[AtomId::new(0)]
        );
        assert_eq!(
            singleton.topology.substance_groups[0].label(),
            Some("singleton")
        );
        assert_eq!(singleton.topology.stereo_groups.len(), 1);
        assert_eq!(singleton.topology.stereo_groups[0].id(), Some(23));
        assert_eq!(
            singleton.topology.stereo_groups[0].kind(),
            cosmolkit_model::StereoGroupKind::Or
        );
        assert_eq!(
            singleton.topology.stereo_groups[0].atoms(),
            &[AtomId::new(0)]
        );
        assert!(singleton.topology.stereo_groups[0].bonds().is_empty());
    }

    #[test]
    fn writer_stereo_merge_remaps_signed_ring_relative_ids_and_only_computed_props() {
        let mut target = parse_smiles("CCCCCC", &SmilesParseParams::default()).unwrap();
        let mut assigned = parse_smiles("CCCC", &SmilesParseParams::default()).unwrap();
        target.topology.atoms[2]
            .set_prop("target_atom", "source")
            .unwrap();
        target.topology.atoms[2]
            .set_computed_prop("_stale_atom", "remove")
            .unwrap();
        target.topology.bonds[2]
            .set_prop("target_bond", "source")
            .unwrap();
        target.topology.bonds[2]
            .set_computed_prop("_stale_bond", "remove")
            .unwrap();
        assigned.topology.atoms[0]
            .set_prop("incoming_atom_ordinary", "do-not-copy")
            .unwrap();
        assigned.topology.atoms[0]
            .set_computed_prop("_ringStereoAtoms", "-2,4")
            .unwrap();
        assigned.topology.atoms[0]
            .set_computed_prop("_incoming_atom_computed", "copy")
            .unwrap();
        assigned.topology.bonds[0]
            .set_prop("incoming_bond_ordinary", "do-not-copy")
            .unwrap();
        assigned.topology.bonds[0]
            .set_computed_prop("_incoming_bond_computed", "copy")
            .unwrap();
        assigned.topology.bonds[0].set_direction(BondDirection::EndUpRight);

        merge_writer_stereo_fragment(
            &mut target.topology,
            &assigned.topology,
            &[
                AtomId::new(2),
                AtomId::new(3),
                AtomId::new(4),
                AtomId::new(5),
            ],
            &[BondId::new(2), BondId::new(3), BondId::new(4)],
        )
        .unwrap();

        assert_eq!(target.topology.atoms[2].prop("target_atom"), Some("source"));
        assert_eq!(
            target.topology.atoms[2].prop("incoming_atom_ordinary"),
            None
        );
        assert_eq!(target.topology.atoms[2].prop("_stale_atom"), None);
        assert_eq!(
            target.topology.atoms[2].prop("_ringStereoAtoms"),
            Some("-4,6")
        );
        assert!(target.topology.atoms[2].is_prop_computed("_ringStereoAtoms"));
        assert_eq!(
            target.topology.atoms[2].prop("_incoming_atom_computed"),
            Some("copy")
        );
        assert!(target.topology.atoms[2].is_prop_computed("_incoming_atom_computed"));
        assert_eq!(target.topology.bonds[2].prop("target_bond"), Some("source"));
        assert_eq!(
            target.topology.bonds[2].prop("incoming_bond_ordinary"),
            None
        );
        assert_eq!(target.topology.bonds[2].prop("_stale_bond"), None);
        assert_eq!(
            target.topology.bonds[2].prop("_incoming_bond_computed"),
            Some("copy")
        );
        assert!(target.topology.bonds[2].is_prop_computed("_incoming_bond_computed"));
        assert_eq!(
            target.topology.bonds[2].direction(),
            BondDirection::EndUpRight
        );
    }

    #[test]
    fn canonical_nonisomeric_stereo_fallback_is_deferred_until_after_ranking() {
        let mut record = parse_smiles("F[C@H](Cl)Br", &SmilesParseParams::default()).unwrap();
        for key in [
            "_CIPCode",
            "_CIPRank",
            "_ChiralityPossible",
            "_ringStereoAtoms",
        ] {
            record.topology.atoms[1].clear_prop(key);
        }
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: true,
            ..Default::default()
        };
        let components = connected_components(&record.topology);

        prepare_writer_stereochemistry(&mut record.topology, None, &params, &components).unwrap();

        assert!(
            record
                .topology
                .atoms
                .iter()
                .all(|atom| atom.prop("_CIPCode").is_none()),
            "canonical non-isomeric stereo assignment belongs after canonical ranking"
        );
    }

    fn assert_s33_standard_component_order_row(
        record: &SmilesRecord,
        params: &SmilesWriteParams,
        expected_text: &str,
        expected_atom_order: &[usize],
        expected_bond_order: &[usize],
        case: &str,
    ) {
        let before = record.clone();
        let output = write_smiles_output(record, params, false).expect(case);
        assert_eq!(output.text, expected_text, "{case}: text");
        assert_eq!(
            output.atom_order,
            expected_atom_order
                .iter()
                .copied()
                .map(AtomId::new)
                .collect::<Vec<_>>(),
            "{case}: atom_order"
        );
        assert_eq!(
            output.bond_order,
            expected_bond_order
                .iter()
                .copied()
                .map(BondId::new)
                .collect::<Vec<_>>(),
            "{case}: bond_order"
        );
        assert_eq!(record, &before, "{case}: input mutation");
    }

    fn s33_standard_params(isomeric: bool) -> SmilesWriteParams {
        SmilesWriteParams {
            do_isomeric_smiles: isomeric,
            do_kekule: false,
            canonical: true,
            clean_stereo: true,
            rooted_at_atom: None,
            all_bonds_explicit: false,
            all_hydrogens_explicit: false,
            include_dative_bonds: true,
            ignore_atom_map_numbers: false,
        }
    }

    fn parse_s33_standard_input(input: &str) -> SmilesRecord {
        let parse_params = SmilesParseParams::default();
        let parsed = parse_smiles(input, &parse_params).unwrap();
        let remove_params = cosmolkit_core::RemoveHsParams {
            update_explicit_count: true,
            sanitize: parse_params.sanitize,
            ..cosmolkit_core::RemoveHsParams::default()
        };
        let prepared = cosmolkit_core::remove_hydrogens_with_params(
            parsed.topology,
            parsed.coordinates,
            parsed.properties,
            &remove_params,
        )
        .unwrap();
        crate::finalize_smiles_stereo(
            SmilesRecord {
                topology: prepared.topology,
                coordinates: prepared.coordinates,
                properties: prepared.properties,
            },
            &parse_params,
        )
        .unwrap()
    }

    fn s33_interleave_oc_co_components(record: &SmilesRecord) -> SmilesRecord {
        assert_eq!(record.topology.atoms.len(), 4);
        assert_eq!(record.topology.bonds.len(), 2);
        let new_to_old = [0, 2, 1, 3];
        let mut old_to_new = [0; 4];
        for (new_index, old_index) in new_to_old.into_iter().enumerate() {
            old_to_new[old_index] = new_index;
        }

        let mut interleaved = record.clone();
        interleaved.topology.atoms = new_to_old
            .into_iter()
            .enumerate()
            .map(|(new_index, old_index)| {
                record.topology.atoms[old_index]
                    .clone()
                    .with_id(AtomId::new(new_index))
            })
            .collect();
        for bond in &mut interleaved.topology.bonds {
            bond.set_endpoints(
                AtomId::new(old_to_new[bond.begin().index()]),
                AtomId::new(old_to_new[bond.end().index()]),
            );
        }
        interleaved.topology.adjacency = cosmolkit_model::AdjacencyList::from_topology(
            interleaved.topology.atoms.len(),
            &interleaved.topology.bonds,
        );
        interleaved.topology.validate().unwrap();
        interleaved
    }

    #[test]
    fn s33_standard_component_order_matches_all_pinned_text_and_map_rows() {
        // The ten rows are the pinned standard-writer oracle in
        // S33_component_order_oracle.md. The CX profile is deliberately not used.
        let unmapped = [
            (
                "C[C@H](F)Cl.[13CH3]O",
                true,
                "C[C@H](F)Cl.[13CH3]O",
                &[0, 1, 2, 3, 4, 5][..],
                &[0, 1, 2, 3][..],
                "unmapped isomeric source order",
            ),
            (
                "[13CH3]O.C[C@H](F)Cl",
                true,
                "C[C@H](F)Cl.[13CH3]O",
                &[2, 3, 4, 5, 0, 1][..],
                &[1, 2, 3, 0][..],
                "unmapped isomeric reversed components",
            ),
            (
                "C[C@H](F)Cl.[13CH3]O",
                false,
                "CC(F)Cl.CO",
                &[0, 1, 2, 3, 4, 5][..],
                &[0, 1, 2, 3][..],
                "unmapped nonisomeric source order",
            ),
            (
                "[13CH3]O.C[C@H](F)Cl",
                false,
                "CC(F)Cl.CO",
                &[2, 3, 4, 5, 0, 1][..],
                &[1, 2, 3, 0][..],
                "unmapped nonisomeric reversed components",
            ),
        ];
        for (input, isomeric, text, atoms, bonds, case) in unmapped {
            let record = parse_s33_standard_input(input);
            assert_s33_standard_component_order_row(
                &record,
                &s33_standard_params(isomeric),
                text,
                atoms,
                bonds,
                case,
            );
        }

        let mapped = [
            (
                "[CH3:10][C@H:11]([F:12])[Cl:13].[13CH3:20][O:21]",
                true,
                "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
                &[4, 5, 0, 1, 2, 3][..],
                &[3, 0, 1, 2][..],
                "mapped isomeric source order",
            ),
            (
                "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
                true,
                "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
                &[0, 1, 2, 3, 4, 5][..],
                &[0, 1, 2, 3][..],
                "mapped isomeric reversed components",
            ),
            (
                "[CH3:10][C@H:11]([F:12])[Cl:13].[13CH3:20][O:21]",
                false,
                "[CH3:10][CH:11]([F:12])[Cl:13].[CH3:20][O:21]",
                &[0, 1, 2, 3, 4, 5][..],
                &[0, 1, 2, 3][..],
                "mapped nonisomeric source order",
            ),
            (
                "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
                false,
                "[CH3:10][CH:11]([F:12])[Cl:13].[CH3:20][O:21]",
                &[2, 3, 4, 5, 0, 1][..],
                &[1, 2, 3, 0][..],
                "mapped nonisomeric reversed components",
            ),
        ];
        for (input, isomeric, text, atoms, bonds, case) in mapped {
            let record = parse_s33_standard_input(input);
            assert_s33_standard_component_order_row(
                &record,
                &s33_standard_params(isomeric),
                text,
                atoms,
                bonds,
                case,
            );
        }

        let source_order = parse_s33_standard_input("OC.CO");
        assert_s33_standard_component_order_row(
            &source_order,
            &s33_standard_params(true),
            "CO.CO",
            &[1, 0, 2, 3],
            &[0, 1],
            "equal-text components in source order",
        );
        let interleaved = s33_interleave_oc_co_components(&source_order);
        assert_eq!(
            interleaved
                .topology
                .bonds
                .iter()
                .map(|bond| (bond.begin().index(), bond.end().index()))
                .collect::<Vec<_>>(),
            [(0, 2), (1, 3)],
            "the [0,2,1,3] input permutation interleaves component rows"
        );
        assert_s33_standard_component_order_row(
            &interleaved,
            &s33_standard_params(true),
            "CO.CO",
            &[1, 3, 2, 0],
            &[1, 0],
            "equal-text interleaved components use atom/bond map tie-break",
        );
    }
}
