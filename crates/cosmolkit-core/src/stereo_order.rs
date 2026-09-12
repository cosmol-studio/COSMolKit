//! Detached tetrahedral ligand ordering and mapping primitives.
//!
//! This module owns parity-preserving order conversion. It does not perceive
//! stereo, mutate topology, or accept a live molecule/runtime capability.

use cosmolkit_model::{
    AtomId, Bond, BondId, MappingValidationError, TopologyBlock, TopologyMapping,
    TopologyValidationError,
};
use cosmolkit_types::{BondOrder, ChiralTag};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TetrahedralLigand {
    Bond(BondId),
    Implicit,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralRemap {
    pub center: AtomId,
    pub ligands: Vec<TetrahedralLigand>,
    pub tag: ChiralTag,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum StereoOrderError {
    #[error("reference order has {reference} entries but probe order has {probe}")]
    PermutationLength { reference: usize, probe: usize },
    #[error(
        "probe order does not contain the value required at reference position {reference_position}"
    )]
    MissingProbeValue { reference_position: usize },
    #[error("invalid topology: {0}")]
    InvalidTopology(TopologyValidationError),
    #[error("invalid topology mapping: {0}")]
    InvalidMapping(MappingValidationError),
    #[error("stereo center {center} is out of range for {atom_count} atoms")]
    CenterOutOfRange { center: AtomId, atom_count: usize },
    #[error("bond {bond} is out of range for {bond_count} bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
    #[error("bond {bond} is not incident to stereo center {center}")]
    BondNotIncident { bond: BondId, center: AtomId },
    #[error("tetrahedral order has {actual} ligands; maximum is {maximum}")]
    TooManyLigands { actual: usize, maximum: usize },
    #[error("bond {bond} is duplicated at tetrahedral ligand position {position}")]
    DuplicateBondLigand { position: usize, bond: BondId },
    #[error("implicit ligand is duplicated at positions {first} and {second}")]
    MultipleImplicitLigands { first: usize, second: usize },
    #[error("chiral tag {tag:?} is not a supported explicit tetrahedral CW/CCW tag")]
    UnsupportedTetrahedralTag { tag: ChiralTag },
    #[error("stereo center {center} was removed by the topology mapping")]
    RemovedCenter { center: AtomId },
    #[error("tetrahedral ligand bond {bond} was removed by the topology mapping")]
    RemovedLigand { bond: BondId },
    #[error("implicit replacement bond {bond} is out of range for {old_bond_count} source bonds")]
    ImplicitReplacementOutOfRange { bond: BondId, old_bond_count: usize },
    #[error("implicit replacement bond {bond} was retained as {mapped}")]
    ImplicitReplacementRetained { bond: BondId, mapped: BondId },
    #[error("implicit replacement bond {bond} is not present in the source ligand order")]
    ImplicitReplacementNotInSource { bond: BondId },
}

/// Count the first-match swaps used by RDKit to convert `probe` to `reference`.
pub fn count_swaps_to_interconvert<T: Copy + Eq>(
    reference: &[T],
    probe: &[T],
) -> Result<usize, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION RDGeneral::countSwapsToInterconvert
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
    // RDKit✔️✔️:   PRECONDITION(ref.size() == probe.size(), "size mismatch");
    // RDKit✔️✔️:   typename T::const_iterator refIt = ref.begin();
    // RDKit✔️✔️:   typename T::iterator probeIt = probe.begin();
    // RDKit✔️✔️:   typename T::iterator probeIt2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int nSwaps = 0;
    // RDKit✔️✔️:   while (refIt != ref.end()) {
    // RDKit✔️✔️:     if ((*probeIt) != (*refIt)) {
    // RDKit✔️✔️:       bool foundIt = false;
    // RDKit✔️✔️:       probeIt2 = probeIt;
    // RDKit✔️✔️:       while ((*probeIt2) != (*refIt) && probeIt2 != probe.end()) {
    // RDKit✔️✔️:         ++probeIt2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (probeIt2 != probe.end()) {
    // RDKit✔️✔️:         foundIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       CHECK_INVARIANT(foundIt, "could not find probe element");
    // RDKit✔️✔️:
    // RDKit✔️✔️:       std::swap(*probeIt, *probeIt2);
    // RDKit✔️✔️:       nSwaps++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++probeIt;
    // RDKit✔️✔️:     ++refIt;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nSwaps;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDGeneral::countSwapsToInterconvert
    if reference.len() != probe.len() {
        return Err(StereoOrderError::PermutationLength {
            reference: reference.len(),
            probe: probe.len(),
        });
    }

    let mut probe = probe.to_vec();
    let mut swaps = 0;
    for (reference_position, expected) in reference.iter().enumerate() {
        if probe[reference_position] == *expected {
            continue;
        }
        let Some(offset) = probe[reference_position..]
            .iter()
            .position(|candidate| candidate == expected)
        else {
            return Err(StereoOrderError::MissingProbeValue { reference_position });
        };
        probe.swap(reference_position, reference_position + offset);
        swaps += 1;
    }
    Ok(swaps)
}

pub fn invert_tetrahedral_tag(tag: ChiralTag) -> Result<ChiralTag, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Atom::invertChirality tetrahedral CW/CCW branches
    // RDKit✔️✔️:     case CHI_TETRAHEDRAL_CW:
    // RDKit✔️✔️:       setChiralTag(CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     case CHI_TETRAHEDRAL_CCW:
    // RDKit✔️✔️:       setChiralTag(CHI_TETRAHEDRAL_CW);
    // RDKit✔️✔️:       return true;
    // END RDKIT CPP FUNCTION Atom::invertChirality tetrahedral CW/CCW branches
    match tag {
        ChiralTag::TetrahedralCw => Ok(ChiralTag::TetrahedralCcw),
        ChiralTag::TetrahedralCcw => Ok(ChiralTag::TetrahedralCw),
        tag => Err(StereoOrderError::UnsupportedTetrahedralTag { tag }),
    }
}

pub fn bond_affects_atom_chirality(bond: &Bond, center: AtomId) -> Result<bool, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::detail::bondAffectsAtomChirality
    // RDKit✔️✔️: bool bondAffectsAtomChirality(const Bond *bond, const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond pointer");
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom pointer");
    // RDKit✔️✔️:   if (bond->getBondType() == Bond::BondType::UNSPECIFIED ||
    // RDKit✔️✔️:       bond->getBondType() == Bond::BondType::ZERO ||
    // RDKit✔️✔️:       (bond->getBondType() == Bond::BondType::DATIVE &&
    // RDKit✔️✔️:        bond->getBeginAtomIdx() == atom->getIdx())) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::detail::bondAffectsAtomChirality
    if bond.begin() != center && bond.end() != center {
        return Err(StereoOrderError::BondNotIncident {
            bond: bond.id(),
            center,
        });
    }
    Ok(
        !matches!(bond.order(), BondOrder::Unspecified | BondOrder::Zero)
            && !(bond.order() == BondOrder::Dative && bond.begin() == center),
    )
}

pub fn atom_nonzero_degree(
    topology: &TopologyBlock,
    center: AtomId,
) -> Result<usize, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::detail::getAtomNonzeroDegree
    // RDKit✔️❌: unsigned int getAtomNonzeroDegree(const Atom *atom) {
    // RDKit✔️❌:   PRECONDITION(atom, "bad pointer");
    // RDKit✔️❌:   PRECONDITION(atom->hasOwningMol(), "no owning molecule");
    // RDKit✔️❌:   unsigned int res = 0;
    // RDKit✔️❌:   for (auto bond : atom->getOwningMol().atomBonds(atom)) {
    // RDKit✔️❌:     if (!bondAffectsAtomChirality(bond, atom)) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++res;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Chirality::detail::getAtomNonzeroDegree
    // The detached boundary validates the full topology before the source
    // degree scan, changing O(degree) to O(V+E); behavior is preserved but the
    // second marker is intentionally negative.
    topology
        .validate()
        .map_err(StereoOrderError::InvalidTopology)?;
    if center.index() >= topology.atoms.len() {
        return Err(StereoOrderError::CenterOutOfRange {
            center,
            atom_count: topology.atoms.len(),
        });
    }
    let mut degree = 0;
    for neighbor in topology.adjacency.neighbors_of(center.index()) {
        let bond =
            topology
                .bonds
                .get(neighbor.bond.index())
                .ok_or(StereoOrderError::BondOutOfRange {
                    bond: neighbor.bond,
                    bond_count: topology.bonds.len(),
                })?;
        if bond_affects_atom_chirality(bond, center)? {
            degree += 1;
        }
    }
    Ok(degree)
}

pub fn incident_tetrahedral_bond_order(
    topology: &TopologyBlock,
    center: AtomId,
) -> Result<Vec<TetrahedralLigand>, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION Atom::getPerturbationOrder incident storage order
    // RDKit✔️❌:   INT_LIST ref;
    // RDKit✔️❌:   for (const auto bnd : getOwningMol().atomBonds(this)) {
    // RDKit✔️❌:     ref.push_back(bnd->getIdx());
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return static_cast<int>(countSwapsToInterconvert(probe, ref));
    // END RDKIT CPP FUNCTION Atom::getPerturbationOrder incident storage order
    // The detached topology validation adds O(V+E) work before retaining the
    // source incident-bond order, so the second marker is intentionally
    // negative. Non-chirality-affecting bonds are then removed by the separate
    // source-backed predicate above.
    topology
        .validate()
        .map_err(StereoOrderError::InvalidTopology)?;
    if center.index() >= topology.atoms.len() {
        return Err(StereoOrderError::CenterOutOfRange {
            center,
            atom_count: topology.atoms.len(),
        });
    }

    let mut ligands = Vec::new();
    for neighbor in topology.adjacency.neighbors_of(center.index()) {
        let bond =
            topology
                .bonds
                .get(neighbor.bond.index())
                .ok_or(StereoOrderError::BondOutOfRange {
                    bond: neighbor.bond,
                    bond_count: topology.bonds.len(),
                })?;
        if bond_affects_atom_chirality(bond, center)? {
            ligands.push(TetrahedralLigand::Bond(bond.id()));
        }
    }
    validate_ligand_order(&ligands)?;
    Ok(ligands)
}

pub fn tetrahedral_tag_after_order_change(
    tag: ChiralTag,
    reference: &[TetrahedralLigand],
    probe: &[TetrahedralLigand],
) -> Result<ChiralTag, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION FindStereo tetrahedral controlling-order normalization
    // RDKit✔️✔️:       unsigned nSwaps =
    // RDKit✔️✔️:           countSwapsToInterconvert(origNbrOrder, sinfo.controllingAtoms);
    // RDKit✔️✔️:       if (nSwaps % 2) {
    // RDKit✔️✔️:         stereo = (stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CCW
    // RDKit✔️✔️:                       ? Atom::ChiralType::CHI_TETRAHEDRAL_CW
    // RDKit✔️✔️:                       : Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:       }
    // END RDKIT CPP FUNCTION FindStereo tetrahedral controlling-order normalization
    validate_tetrahedral_tag(tag)?;
    validate_ligand_order(reference)?;
    validate_ligand_order(probe)?;
    let swaps = count_swaps_to_interconvert(reference, probe)?;
    if swaps % 2 == 1 {
        invert_tetrahedral_tag(tag)
    } else {
        Ok(tag)
    }
}

#[allow(clippy::too_many_arguments)]
pub fn remap_tetrahedral_center(
    source_center: AtomId,
    source_ligands: &[TetrahedralLigand],
    source_tag: ChiralTag,
    mapping: &TopologyMapping,
    old_atom_count: usize,
    new_atom_count: usize,
    old_bond_count: usize,
    new_bond_count: usize,
    removed_bond_as_implicit: Option<BondId>,
    target_ligands: &[TetrahedralLigand],
) -> Result<TetrahedralRemap, StereoOrderError> {
    // BEGIN RDKIT CPP FUNCTION AddHs remove-H tetrahedral order adjustment
    // RDKit✔️❌:       INT_LIST neighborIndices;
    // RDKit✔️❌:       for (const auto &nbnd : mol.atomBonds(heavyAtom)) {
    // RDKit✔️❌:         if (nbnd->getIdx() != bond->getIdx()) {
    // RDKit✔️❌:           neighborIndices.push_back(nbnd->getIdx());
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       neighborIndices.push_back(bond->getIdx());
    // RDKit✔️❌:
    // RDKit✔️❌:       int nSwaps = heavyAtom->getPerturbationOrder(neighborIndices);
    // RDKit✔️❌:       if (nSwaps % 2) {
    // RDKit✔️❌:         heavyAtom->invertChirality();
    // RDKit✔️❌:       }
    // END RDKIT CPP FUNCTION AddHs remove-H tetrahedral order adjustment
    // The complete detached mapping validation is O(V+E), while the source
    // adjustment only scans one center's incident bonds. This intentional
    // fail-closed boundary accounts for the negative performance marker.
    mapping
        .validate_for_counts(
            old_atom_count,
            new_atom_count,
            old_bond_count,
            new_bond_count,
        )
        .map_err(StereoOrderError::InvalidMapping)?;
    validate_tetrahedral_tag(source_tag)?;
    validate_ligand_order(source_ligands)?;
    validate_ligand_order(target_ligands)?;
    validate_bond_ranges(source_ligands, old_bond_count)?;
    validate_bond_ranges(target_ligands, new_bond_count)?;

    if source_center.index() >= old_atom_count {
        return Err(StereoOrderError::CenterOutOfRange {
            center: source_center,
            atom_count: old_atom_count,
        });
    }
    let target_center = mapping.atoms().old_to_new()[source_center.index()].ok_or(
        StereoOrderError::RemovedCenter {
            center: source_center,
        },
    )?;

    if let Some(replacement) = removed_bond_as_implicit {
        if replacement.index() >= old_bond_count {
            return Err(StereoOrderError::ImplicitReplacementOutOfRange {
                bond: replacement,
                old_bond_count,
            });
        }
        if !source_ligands.contains(&TetrahedralLigand::Bond(replacement)) {
            return Err(StereoOrderError::ImplicitReplacementNotInSource { bond: replacement });
        }
        if let Some(mapped) = mapping.bonds().old_to_new()[replacement.index()] {
            return Err(StereoOrderError::ImplicitReplacementRetained {
                bond: replacement,
                mapped,
            });
        }
    }

    let mut mapped_ligands = Vec::with_capacity(source_ligands.len());
    for ligand in source_ligands {
        match *ligand {
            TetrahedralLigand::Implicit => mapped_ligands.push(TetrahedralLigand::Implicit),
            TetrahedralLigand::Bond(bond) => match mapping.bonds().old_to_new()[bond.index()] {
                Some(mapped) => mapped_ligands.push(TetrahedralLigand::Bond(mapped)),
                None if removed_bond_as_implicit == Some(bond) => {
                    mapped_ligands.push(TetrahedralLigand::Implicit);
                }
                None => return Err(StereoOrderError::RemovedLigand { bond }),
            },
        }
    }
    validate_ligand_order(&mapped_ligands)?;

    let tag = tetrahedral_tag_after_order_change(source_tag, &mapped_ligands, target_ligands)?;
    Ok(TetrahedralRemap {
        center: target_center,
        ligands: target_ligands.to_vec(),
        tag,
    })
}

fn validate_tetrahedral_tag(tag: ChiralTag) -> Result<(), StereoOrderError> {
    if matches!(tag, ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw) {
        Ok(())
    } else {
        Err(StereoOrderError::UnsupportedTetrahedralTag { tag })
    }
}

fn validate_ligand_order(ligands: &[TetrahedralLigand]) -> Result<(), StereoOrderError> {
    const MAXIMUM: usize = 4;
    if ligands.len() > MAXIMUM {
        return Err(StereoOrderError::TooManyLigands {
            actual: ligands.len(),
            maximum: MAXIMUM,
        });
    }
    let mut implicit_position = None;
    for (position, ligand) in ligands.iter().copied().enumerate() {
        match ligand {
            TetrahedralLigand::Implicit => {
                if let Some(first) = implicit_position {
                    return Err(StereoOrderError::MultipleImplicitLigands {
                        first,
                        second: position,
                    });
                }
                implicit_position = Some(position);
            }
            TetrahedralLigand::Bond(bond) => {
                if ligands[..position].contains(&TetrahedralLigand::Bond(bond)) {
                    return Err(StereoOrderError::DuplicateBondLigand { position, bond });
                }
            }
        }
    }
    Ok(())
}

fn validate_bond_ranges(
    ligands: &[TetrahedralLigand],
    bond_count: usize,
) -> Result<(), StereoOrderError> {
    for ligand in ligands {
        if let TetrahedralLigand::Bond(bond) = ligand
            && bond.index() >= bond_count
        {
            return Err(StereoOrderError::BondOutOfRange {
                bond: *bond,
                bond_count,
            });
        }
    }
    Ok(())
}
