//! Original COSMolKit native 64-bit molecular hash over detached state.
//! This is the retained project algorithm, not RDKit MolHash or a replacement digest.
use crate::hash::hash_combine;
use cosmolkit_core::{CipRankError, RingInfo, ValenceAssignment, assign_atom_cip_ranks};
use cosmolkit_model::{
    AtomId, BondOrder, ChiralTag, NeighborRef, TopologyBlock, TopologyValidationError,
};
use std::fmt;

/// Structured original hashing and detached-input failures.
#[derive(Debug)]
pub enum MoleculeHashError {
    EmptyMolecule,
    MissingPreparedValence,
    CipRanks(CipRankError),
    InvalidTopology(TopologyValidationError),
    RankCount { actual: usize, atom_count: usize },
}
impl fmt::Display for MoleculeHashError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMolecule => f.write_str("cannot hash an empty molecule"),
            Self::MissingPreparedValence => {
                f.write_str("CIP ranking requires a prepared valence assignment")
            }
            Self::CipRanks(e) => write!(f, "molecular hash CIP ranking failed: {e}"),
            Self::InvalidTopology(e) => write!(f, "invalid molecular hash topology: {e}"),
            Self::RankCount { actual, atom_count } => write!(
                f,
                "molecular hash rank row has length {actual}, expected {atom_count}"
            ),
        }
    }
}
impl std::error::Error for MoleculeHashError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::CipRanks(e) => Some(e),
            Self::InvalidTopology(e) => Some(e),
            _ => None,
        }
    }
}

/// Native hash using the original legacy CIP ranks and already prepared values.
/// No valence or ring perception runs when its explicit input is absent.
pub fn molecule_hash(
    topology: &TopologyBlock,
    valence: Option<&ValenceAssignment>,
    rings: Option<&RingInfo>,
) -> Result<u64, MoleculeHashError> {
    // COSMolKit❗✔️: pub fn mol_hash(mol: &Molecule) -> Result<u64, HashError> {
    // COSMolKit❗✔️:     if mol.num_atoms() == 0 {
    // COSMolKit❗✔️:         return Err(HashError::EmptyMolecule);
    // COSMolKit❗✔️:     }
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:     let ranks = assign_atom_cip_ranks(mol)?;
    // COSMolKit❗✔️:     mol_hash_with_ranks(mol, &ranks)
    // COSMolKit❗✔️: }
    if topology.atoms.is_empty() {
        return Err(MoleculeHashError::EmptyMolecule);
    }
    let valence = valence.ok_or(MoleculeHashError::MissingPreparedValence)?;
    let ranks = assign_atom_cip_ranks(topology, valence).map_err(MoleculeHashError::CipRanks)?;
    molecule_hash_with_ranks(topology, rings, &ranks)
}

/// Native hash using one caller-supplied rank per atom.
/// Original rank/index and stable neighbor-rank sorting, charge wrapping, ring
/// absence and both 32-bit folds are retained. Borrow adjacency instead of
/// cloning it: the read-only algorithm cannot change ordering or graph state.
pub fn molecule_hash_with_ranks(
    topology: &TopologyBlock,
    rings: Option<&RingInfo>,
    ranks: &[u32],
) -> Result<u64, MoleculeHashError> {
    // COSMolKit❗✔️: pub fn mol_hash_with_ranks(mol: &Molecule, ranks: &[u32]) -> Result<u64, HashError> {
    // COSMolKit❗✔️:     if mol.num_atoms() == 0 {
    // COSMolKit❗✔️:         return Err(HashError::EmptyMolecule);
    // COSMolKit❗✔️:     }
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:     let n = mol.num_atoms();
    // COSMolKit❗✔️:     let adjacency = mol.topology_block().adjacency.clone();
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:     let ri = mol.derived_cache().rings.as_ref();
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:     // Build (rank, original_index) pairs and sort by rank
    // COSMolKit❗✔️:     let mut order: Vec<(u32, usize)> = ranks
    // COSMolKit❗✔️:         .iter()
    // COSMolKit❗✔️:         .copied()
    // COSMolKit❗✔️:         .enumerate()
    // COSMolKit❗✔️:         .map(|(i, r)| (r, i))
    // COSMolKit❗✔️:         .collect();
    // COSMolKit❗✔️:     order.sort();
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:     let mut hash_hi: u32 = 0;
    // COSMolKit❗✔️:     let mut hash_lo: u32 = 0;
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:     for &(_rank, atom_idx) in &order {
    // COSMolKit❗✔️:         let atom = &mol.atoms()[atom_idx];
    // COSMolKit❗✔️:         let mut atom_hash: u32 = 0;
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:         // Atomic number
    // COSMolKit❗✔️:         hash_combine(&mut atom_hash, atom.atomic_number() as u32);
    // COSMolKit❗✔️:         // Formal charge (offset by 8 to make it unsigned)
    // COSMolKit❗✔️:         hash_combine(&mut atom_hash, atom.formal_charge().wrapping_add(8) as u32);
    // COSMolKit❗✔️:         // Isotope
    // COSMolKit❗✔️:         hash_combine(&mut atom_hash, atom.isotope().unwrap_or(0) as u32);
    // COSMolKit❗✔️:         // Chirality
    // COSMolKit❗✔️:         let chirality_code = chiral_tag_to_hash_code(atom.chiral_tag());
    // COSMolKit❗✔️:         hash_combine(&mut atom_hash, chirality_code);
    // COSMolKit❗✔️:         // Aromaticity
    // COSMolKit❗✔️:         hash_combine(&mut atom_hash, if atom.is_aromatic() { 1u32 } else { 0u32 });
    // COSMolKit❗✔️:         // Ring membership
    // COSMolKit❗✔️:         let in_ring = ri.map_or(false, |rinfo| {
    // COSMolKit❗✔️:             rinfo.num_atom_rings(AtomId::new(atom_idx)) > 0
    // COSMolKit❗✔️:         });
    // COSMolKit❗✔️:         hash_combine(&mut atom_hash, if in_ring { 1u32 } else { 0u32 });
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:         // Neighbors: sorted by neighbor rank
    // COSMolKit❗✔️:         let nbrs = adjacency.neighbors_of(atom_idx);
    // COSMolKit❗✔️:         let mut nbr_entries: Vec<(u32, &NeighborRef)> = nbrs
    // COSMolKit❗✔️:             .iter()
    // COSMolKit❗✔️:             .map(|nbr| (ranks[nbr.atom_index], nbr))
    // COSMolKit❗✔️:             .collect();
    // COSMolKit❗✔️:         nbr_entries.sort_by_key(|&(r, _)| r);
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:         for (_nbr_rank, nbr) in &nbr_entries {
    // COSMolKit❗✔️:             let bond = &mol.bonds()[nbr.bond.index()];
    // COSMolKit❗✔️:             let bond_val = bond_order_to_hash_val(bond.order());
    // COSMolKit❗✔️:             hash_combine(&mut atom_hash, bond_val);
    // COSMolKit❗✔️:             let nbr_atom = &mol.atoms()[nbr.atom_index];
    // COSMolKit❗✔️:             hash_combine(&mut atom_hash, nbr_atom.atomic_number() as u32);
    // COSMolKit❗✔️:             hash_combine(&mut atom_hash, if bond.is_aromatic() { 1u32 } else { 0u32 });
    // COSMolKit❗✔️:         }
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:         // Fold atom_hash into the overall 64-bit result
    // COSMolKit❗✔️:         hash_combine(&mut hash_hi, atom_hash);
    // COSMolKit❗✔️:         hash_combine(&mut hash_lo, atom_hash.wrapping_mul(0x9e3779b9));
    // COSMolKit❗✔️:     }
    // COSMolKit❗✔️:
    // COSMolKit❗✔️:     Ok(((hash_hi as u64) << 32) | (hash_lo as u64))
    // COSMolKit❗✔️: }
    if topology.atoms.len() == 0 {
        return Err(MoleculeHashError::EmptyMolecule);
    }

    topology
        .validate()
        .map_err(MoleculeHashError::InvalidTopology)?;
    if ranks.len() != topology.atoms.len() {
        return Err(MoleculeHashError::RankCount {
            actual: ranks.len(),
            atom_count: topology.atoms.len(),
        });
    }
    let adjacency = &topology.adjacency;

    let ri = rings;

    // Build (rank, original_index) pairs and sort by rank
    let mut order: Vec<(u32, usize)> = ranks
        .iter()
        .copied()
        .enumerate()
        .map(|(i, r)| (r, i))
        .collect();
    order.sort();

    let mut hash_hi: u32 = 0;
    let mut hash_lo: u32 = 0;

    for &(_rank, atom_idx) in &order {
        let atom = &topology.atoms[atom_idx];
        let mut atom_hash: u32 = 0;

        // Atomic number
        hash_combine(&mut atom_hash, atom.atomic_number() as u32);
        // Formal charge (offset by 8 to make it unsigned)
        hash_combine(&mut atom_hash, atom.formal_charge().wrapping_add(8) as u32);
        // Isotope
        hash_combine(&mut atom_hash, atom.isotope().unwrap_or(0) as u32);
        // Chirality
        let chirality_code = chiral_tag_to_hash_code(atom.chiral_tag());
        hash_combine(&mut atom_hash, chirality_code);
        // Aromaticity
        hash_combine(&mut atom_hash, if atom.is_aromatic() { 1u32 } else { 0u32 });
        // Ring membership
        let in_ring = ri.map_or(false, |rinfo| {
            rinfo.num_atom_rings(AtomId::new(atom_idx)) > 0
        });
        hash_combine(&mut atom_hash, if in_ring { 1u32 } else { 0u32 });

        // Neighbors: sorted by neighbor rank
        let nbrs = adjacency.neighbors_of(atom_idx);
        let mut nbr_entries: Vec<(u32, &NeighborRef)> = nbrs
            .iter()
            .map(|nbr| (ranks[nbr.atom_index], nbr))
            .collect();
        nbr_entries.sort_by_key(|&(r, _)| r);

        for (_nbr_rank, nbr) in &nbr_entries {
            let bond = &topology.bonds[nbr.bond.index()];
            let bond_val = bond_order_to_hash_val(bond.order());
            hash_combine(&mut atom_hash, bond_val);
            let nbr_atom = &topology.atoms[nbr.atom_index];
            hash_combine(&mut atom_hash, nbr_atom.atomic_number() as u32);
            hash_combine(&mut atom_hash, if bond.is_aromatic() { 1u32 } else { 0u32 });
        }

        // Fold atom_hash into the overall 64-bit result
        hash_combine(&mut hash_hi, atom_hash);
        hash_combine(&mut hash_lo, atom_hash.wrapping_mul(0x9e3779b9));
    }

    Ok(((hash_hi as u64) << 32) | (hash_lo as u64))
}

fn bond_order_to_hash_val(order: BondOrder) -> u32 {
    // COSMolKit❗✔️: fn bond_order_to_hash_val(order: BondOrder) -> u32 {
    // COSMolKit❗✔️:     match order {
    // COSMolKit❗✔️:         BondOrder::Single => 1,
    // COSMolKit❗✔️:         BondOrder::Double => 2,
    // COSMolKit❗✔️:         BondOrder::Triple => 3,
    // COSMolKit❗✔️:         BondOrder::Quadruple => 4,
    // COSMolKit❗✔️:         BondOrder::Aromatic => 5,
    // COSMolKit❗✔️:         BondOrder::OneAndHalf => 6,
    // COSMolKit❗✔️:         BondOrder::TwoAndHalf => 7,
    // COSMolKit❗✔️:         BondOrder::ThreeAndHalf => 8,
    // COSMolKit❗✔️:         BondOrder::Dative | BondOrder::DativeOne => 9,
    // COSMolKit❗✔️:         BondOrder::Hydrogen => 10,
    // COSMolKit❗✔️:         BondOrder::Ionic => 11,
    // COSMolKit❗✔️:         BondOrder::Other | BondOrder::Unspecified | BondOrder::Zero | BondOrder::Null => 0,
    // COSMolKit❗✔️:         _ => 0,
    // COSMolKit❗✔️:     }
    // COSMolKit❗✔️: }
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 5,
        BondOrder::OneAndHalf => 6,
        BondOrder::TwoAndHalf => 7,
        BondOrder::ThreeAndHalf => 8,
        BondOrder::Dative | BondOrder::DativeOne => 9,
        BondOrder::Hydrogen => 10,
        BondOrder::Ionic => 11,
        BondOrder::Other | BondOrder::Unspecified | BondOrder::Zero => 0,
        _ => 0,
    }
}

fn chiral_tag_to_hash_code(tag: ChiralTag) -> u32 {
    // COSMolKit❗✔️: fn chiral_tag_to_hash_code(tag: ChiralTag) -> u32 {
    // COSMolKit❗✔️:     match tag {
    // COSMolKit❗✔️:         ChiralTag::Unspecified => 0,
    // COSMolKit❗✔️:         ChiralTag::TetrahedralCw => 3,  // R
    // COSMolKit❗✔️:         ChiralTag::TetrahedralCcw => 2, // S
    // COSMolKit❗✔️:         ChiralTag::Other => 1,
    // COSMolKit❗✔️:         ChiralTag::Tetrahedral => 4,
    // COSMolKit❗✔️:         ChiralTag::Allene => 5,
    // COSMolKit❗✔️:         ChiralTag::SquarePlanar => 6,
    // COSMolKit❗✔️:         ChiralTag::TrigonalBipyramidal => 7,
    // COSMolKit❗✔️:         ChiralTag::Octahedral => 8,
    // COSMolKit❗✔️:     }
    // COSMolKit❗✔️: }
    match tag {
        ChiralTag::Unspecified => 0,
        ChiralTag::TetrahedralCw => 3,  // R
        ChiralTag::TetrahedralCcw => 2, // S
        ChiralTag::Other => 1,
        ChiralTag::Tetrahedral => 4,
        ChiralTag::Allene => 5,
        ChiralTag::SquarePlanar => 6,
        ChiralTag::TrigonalBipyramidal => 7,
        ChiralTag::Octahedral => 8,
    }
}
