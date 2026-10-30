//! Project-native stereo reads adapted from preserved d892 source.
use cosmolkit_core::{
    PotentialStereoCenter, PotentialStereoDescriptor, PotentialStereoParams,
    PotentialStereoSpecified, PotentialStereoType, potential_stereo,
};
use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{
    AtomId, BondId, LigandRef, TetrahedralStereo, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::BondOrder;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum StereoReadError {
    #[error(transparent)]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("valence preparation failed: {0}")]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error("ring preparation failed: {0}")]
    Rings(#[from] cosmolkit_core::RingFindingError),
    #[error("potential stereo perception failed: {0}")]
    PotentialStereo(#[from] cosmolkit_core::PotentialStereoError),
    #[error("CIP label assignment failed: {0}")]
    CipLabeler(#[from] crate::CipLabelerError),
    #[error("CIP label property conversion failed: {0}")]
    PropertyString(#[from] cosmolkit_core::PropertyStringError),
    #[error("CIP label at {atom:?} is not UTF-8: {source}")]
    CipLabelEncoding {
        atom: AtomId,
        #[source]
        source: std::string::FromUtf8Error,
    },
}

pub fn tetrahedral_stereo(
    topology: &TopologyBlock,
    source_valence: Option<&ValenceAssignment>,
    source_rings: Option<&RingInfo>,
    include_unassigned: bool,
) -> Result<Vec<TetrahedralStereo>, StereoReadError> {
    let (valence, rings) =
        crate::chiral_centers::prepare_chiral_state(topology, source_valence, source_rings)?;
    let perception = potential_stereo(
        topology,
        &valence,
        &rings,
        &PotentialStereoParams::default(),
    )?;
    let mut result = Vec::new();
    // RDKit✔️✔️:       for si in itms:
    // RDKit✔️✔️:         if si.type == StereoType.Atom_Tetrahedral and (includeUnassigned or si.specified
    // RDKit✔️✔️:                                                        == StereoSpecified.Specified):
    // FindMolChiralCenters modern filtering; reuse its perception owner without
    // computing CIP labels. Ordered ligands below are COSMolKit's value convention.
    for info in perception.stereo {
        if info.stereo_type != PotentialStereoType::AtomTetrahedral
            || (!include_unassigned && info.specified != PotentialStereoSpecified::Specified)
        {
            continue;
        }
        let PotentialStereoCenter::Atom(center) = info.centered_on else {
            unreachable!("tetrahedral perception uses atom identifiers");
        };
        if info.specified != PotentialStereoSpecified::Specified {
            result.push(TetrahedralStereo {
                center,
                ligands: None,
            });
            continue;
        }
        // FindStereo.cpp::getStereoInfo already folds neighbor-order parity
        // into Tet_CW/Tet_CCW relative to the sorted controlling atoms.
        // RDKit✔️✔️:   std::sort(sinfo.controllingAtoms.begin(), sinfo.controllingAtoms.end());
        // The source permits P/As with two atom neighbors plus H and a lone
        // pair, and pyramidal N/P/As/S/Se with three occupied positions.
        // Eligibility is decided only by the perception owner, never here.
        let hydrogens = cosmolkit_core::total_hydrogen_count_from_validated(
            topology, &*valence, center, false,
        )? as usize;
        let occupied = info.controlling_atoms.len() + hydrogens;
        if !(3..=4).contains(&occupied) || info.controlling_atoms.iter().any(Option::is_none) {
            return Err(cosmolkit_core::PotentialStereoError::InvalidAtomDegree {
                atom: center,
                degree: occupied,
            }
            .into());
        }
        let mut ligands = [LigandRef::LonePair; 4];
        for (slot, neighbor) in ligands.iter_mut().zip(&info.controlling_atoms) {
            *slot = LigandRef::Atom(neighbor.expect("checked controlling atom"));
        }
        for slot in &mut ligands[info.controlling_atoms.len()..occupied] {
            *slot = LigandRef::ImplicitHydrogen;
        }
        match info.descriptor {
            PotentialStereoDescriptor::TetrahedralClockwise => ligands.swap(0, 1),
            PotentialStereoDescriptor::TetrahedralCounterclockwise => {}
            _ => {
                return Err(
                    cosmolkit_core::PotentialStereoError::InvalidChiralPermutation {
                        atom: center,
                        permutation: info.permutation,
                    }
                    .into(),
                );
            }
        }
        result.push(TetrahedralStereo {
            center,
            ligands: Some(canonicalize_tetrahedral_ligands(ligands)),
        });
    }
    Ok(result)
}

fn canonicalize_tetrahedral_ligands(ligands: [LigandRef; 4]) -> [LigandRef; 4] {
    // Keep virtual positions last; even permutations preserve handedness.
    let last = if ligands.contains(&LigandRef::LonePair) {
        Some(LigandRef::LonePair)
    } else if ligands.contains(&LigandRef::ImplicitHydrogen) {
        Some(LigandRef::ImplicitHydrogen)
    } else {
        None
    };
    let both_virtual =
        ligands.contains(&LigandRef::LonePair) && ligands.contains(&LigandRef::ImplicitHydrogen);
    let mut best = None;
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let perm = [a, b, c, d];
                    if has_duplicate_indices(perm) || !is_even_permutation(perm) {
                        continue;
                    }
                    let candidate = [ligands[a], ligands[b], ligands[c], ligands[d]];
                    if last.is_some_and(|last| candidate[3] != last)
                        || (both_virtual && candidate[2] != LigandRef::ImplicitHydrogen)
                    {
                        continue;
                    }
                    if best.is_none_or(|best| candidate < best) {
                        best = Some(candidate);
                    }
                }
            }
        }
    }
    best.expect("an even permutation places virtual ligands last")
}

fn has_duplicate_indices(indices: [usize; 4]) -> bool {
    indices[0] == indices[1]
        || indices[0] == indices[2]
        || indices[0] == indices[3]
        || indices[1] == indices[2]
        || indices[1] == indices[3]
        || indices[2] == indices[3]
}

fn is_even_permutation(indices: [usize; 4]) -> bool {
    let mut inversions = 0;
    for i in 0..4 {
        for j in (i + 1)..4 {
            if indices[i] > indices[j] {
                inversions += 1;
            }
        }
    }
    inversions % 2 == 0
}

fn should_detect_double_bond_stereo(
    topology: &TopologyBlock,
    rings: Option<&RingInfo>,
    bond: BondId,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION shouldDetectDoubleBondStereo (Chirality.cpp)
    // RDKit✔️✔️: bool shouldDetectDoubleBondStereo(const Bond *bond) {
    // RDKit✔️✔️:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
    // RDKit✔️✔️:   return (!ri->numBondRings(bond->getIdx()) ||
    // RDKit✔️✔️:           ri->minBondRingSize(bond->getIdx()) >=
    // RDKit✔️✔️:               Chirality::minRingSizeForDoubleBondStereo);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    // Preserve the original project caller's Double/Aromatic gate and absent
    // ring assignment branch. Queries borrow the existing indexed ring owner.
    let bond = &topology.bonds[bond.index()];
    if bond.order() != BondOrder::Double && bond.order() != BondOrder::Aromatic {
        return false;
    }
    let Some(rings) = rings else {
        return true;
    };
    rings.num_bond_rings(bond.id()) == 0 || rings.min_bond_ring_size(bond.id()) >= 8
}

pub fn perceive_stereochemistry(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
    rings: Option<&RingInfo>,
) -> Result<(), StereoReadError> {
    // COSMolKit-native source d892 stereo.rs:392. This is a read-only validation
    // query: no chiral tag, CIP descriptor or geometry is assigned.
    let _ = tetrahedral_stereo(topology, valence, rings, false)?;
    for bond in &topology.bonds {
        let _ = should_detect_double_bond_stereo(topology, rings, bond.id());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn tetrahedral_stereo_canonicalizes_even_ligand_permutations() {
        let base = [
            LigandRef::Atom(AtomId::new(0)),
            LigandRef::Atom(AtomId::new(2)),
            LigandRef::Atom(AtomId::new(3)),
            LigandRef::Atom(AtomId::new(4)),
        ];

        for a in 0..4 {
            for b in 0..4 {
                for c in 0..4 {
                    for d in 0..4 {
                        let perm = [a, b, c, d];
                        if has_duplicate_indices(perm) {
                            continue;
                        }
                        let candidate = [base[a], base[b], base[c], base[d]];
                        if is_even_permutation(perm) {
                            assert_eq!(canonicalize_tetrahedral_ligands(candidate), base);
                        } else {
                            assert_ne!(canonicalize_tetrahedral_ligands(candidate), base);
                        }
                    }
                }
            }
        }
    }
}
