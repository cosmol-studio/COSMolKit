//! Project-native stereo reads adapted from preserved d892 source.
use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{
    AtomId, BondId, LigandRef, TetrahedralStereo, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, ChiralTag};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum StereoReadError {
    #[error(transparent)]
    InvalidTopology(#[from] TopologyValidationError),
}

pub fn tetrahedral_stereo(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
) -> Result<Vec<TetrahedralStereo>, StereoReadError> {
    // COSMolKit-native source d892 chemistry/stereo.rs:209-347. This is the
    // project ordered-ligand contract, not RDKit assignStereochemistry.
    // Cost review: borrow existing adjacency instead of cloning it; retain one
    // neighbor vector per tagged center and the same fixed 4^4 permutation walk.
    topology.validate()?;
    // Detect tetrahedral stereo centers from the typed atom state.
    // Atoms with ChiralTag::TetrahedralCw or TetrahedralCcw are stereo
    // centers. The four ligands are the explicit neighbors plus an
    // implicit hydrogen if degree < 4.
    let adjacency = &topology.adjacency;

    let mut result = Vec::new();
    for atom in &topology.atoms {
        let tag = atom.chiral_tag();
        if tag != ChiralTag::TetrahedralCw && tag != ChiralTag::TetrahedralCcw {
            continue;
        }
        let center = atom.id();
        let nbrs: Vec<AtomId> = adjacency
            .neighbors_of(center.index())
            .iter()
            .map(|n| n.atom_index)
            .map(AtomId::new)
            .collect();
        let degree = nbrs.len();
        let implicit_hs = valence
            .and_then(|v| v.implicit_hydrogens.get(center.index()).copied())
            .unwrap_or(0) as usize;
        let hydrogen_ligands = atom.explicit_hydrogens() as usize + implicit_hs;

        // Must have 4 ligands total (explicit + implicit)
        if degree + hydrogen_ligands != 4 {
            continue;
        }

        // Build ligand list: explicit neighbors first, then implicit H
        let mut ligands: [LigandRef; 4] = [
            LigandRef::ImplicitHydrogen,
            LigandRef::ImplicitHydrogen,
            LigandRef::ImplicitHydrogen,
            LigandRef::ImplicitHydrogen,
        ];
        for (i, &nbr) in nbrs.iter().enumerate().take(4) {
            ligands[i] = LigandRef::Atom(nbr);
        }

        // COSMolKit defines tetrahedral stereo as center + ordered ligands.
        // RDKit-style CW/CCW tags are compatibility input state, so fold their
        // parity into an odd ligand permutation instead of carrying a second
        // orientation flag.
        let perm = atom.chiral_permutation().unwrap_or(0);
        if matches!(
            (tag, perm % 2),
            (ChiralTag::TetrahedralCw, 0) | (ChiralTag::TetrahedralCcw, 1)
        ) {
            ligands.swap(0, 1);
        }

        let ligands = canonicalize_tetrahedral_ligands(ligands);

        result.push(TetrahedralStereo { center, ligands });
    }
    Ok(result)
}

fn canonicalize_tetrahedral_ligands(ligands: [LigandRef; 4]) -> [LigandRef; 4] {
    if ligands.contains(&LigandRef::ImplicitHydrogen) {
        return canonicalize_tetrahedral_ligands_with_implicit_hydrogen(ligands);
    }

    let mut best = ligands;

    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let perm = [a, b, c, d];
                    if has_duplicate_indices(perm) || !is_even_permutation(perm) {
                        continue;
                    }
                    let candidate = [
                        ligands[perm[0]],
                        ligands[perm[1]],
                        ligands[perm[2]],
                        ligands[perm[3]],
                    ];
                    if candidate < best {
                        best = candidate;
                    }
                }
            }
        }
    }

    best
}

fn canonicalize_tetrahedral_ligands_with_implicit_hydrogen(
    ligands: [LigandRef; 4],
) -> [LigandRef; 4] {
    let mut best = ligands;

    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let perm = [a, b, c, d];
                    if has_duplicate_indices(perm) || !is_even_permutation(perm) {
                        continue;
                    }
                    let candidate = [
                        ligands[perm[0]],
                        ligands[perm[1]],
                        ligands[perm[2]],
                        ligands[perm[3]],
                    ];
                    if candidate[3] == LigandRef::ImplicitHydrogen && candidate < best {
                        best = candidate;
                    }
                }
            }
        }
    }

    best
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
    let _ = tetrahedral_stereo(topology, valence)?;
    for bond in &topology.bonds {
        let _ = should_detect_double_bond_stereo(topology, rings, bond.id());
    }
    Ok(())
}

pub fn find_chiral_centers(
    topology: &TopologyBlock,
    include_unassigned: bool,
) -> Vec<(usize, String)> {
    // COSMolKit-native source: preserved Python find_chiral_centers filter.
    // Retain all unspecified atoms when requested and the original tag labels;
    // this query does not claim to calculate CIP labels or potential centers.
    topology
        .atoms
        .iter()
        .filter_map(|atom| match atom.chiral_tag() {
            ChiralTag::Unspecified if include_unassigned => {
                Some((atom.id().index(), "?".to_owned()))
            }
            ChiralTag::TetrahedralCw => Some((atom.id().index(), "CHI_TETRAHEDRAL_CW".to_owned())),
            ChiralTag::TetrahedralCcw => {
                Some((atom.id().index(), "CHI_TETRAHEDRAL_CCW".to_owned()))
            }
            ChiralTag::TrigonalBipyramidal => {
                Some((atom.id().index(), "CHI_TRIGONALBIPYRAMIDAL".to_owned()))
            }
            _ => None,
        })
        .collect()
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
