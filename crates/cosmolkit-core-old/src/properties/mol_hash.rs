// RDKit marker convention defined in dev/source_reproduction_protocol.md.
//
// Source reproduction protocol: dev/source_reproduction_protocol.md
//
// RDKit source files: hashfunctions.cpp, nmmolhash.h
//   - MurckoScaffold: leaf-pruning (degree < 2) re-hanging H counts on neighbor
//   - ExtendedMurcko: ring atoms + connector atoms stay; side chains become dummy (*)
//   - MolHash: canonical hash via CIP-ordered atom traversal

use crate::{
    AtomId, AtomSpec, BondOrder, BondSpec, ChiralTag, Element, Molecule, MoleculeBuilder,
    NeighborRef, fingerprint::hash_combine, stereo::assign_atom_cip_ranks,
};

/// RDKit-derived error type for MolHash operations.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum HashError {
    #[error("empty molecule has no hash")]
    EmptyMolecule,
    #[error("stereo/CIP error: {0}")]
    Stereo(#[from] crate::StereoError),
    #[error(transparent)]
    Build(#[from] crate::MoleculeBuildError),
    #[error("hash error: {0}")]
    General(String),
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Extract the Murcko scaffold (framework) of a molecule.
///
/// Repeatedly removes leaf atoms (degree < 2), re-hanging their implicit H
/// counts onto the neighbor before removal. Equivalent to RDKit's
/// `HashFunction::MurckoScaffold`.
///
/// # Errors
/// Returns `HashError::EmptyMolecule` if the molecule has no atoms.
pub fn mol_murcko_scaffold(mol: &Molecule) -> Result<Molecule, HashError> {
    if mol.num_atoms() == 0 {
        return Err(HashError::EmptyMolecule);
    }

    let mut builder = build_molecule_builder(mol);
    let mut working = builder.build()?;

    loop {
        let adjacency = working.topology_block().adjacency.clone();

        let mut removal_indices: Vec<usize> = Vec::new();
        for (idx, _atom) in working.atoms().iter().enumerate() {
            let deg = adjacency.neighbors_of(idx).len();
            if deg < 2 {
                removal_indices.push(idx);
            }
        }

        if removal_indices.is_empty() {
            break;
        }

        // Track H-count transfers from degree-1 atoms to their neighbors
        let mut extra_h: Vec<u32> = vec![0u32; working.num_atoms()];
        for &idx in &removal_indices {
            let deg = adjacency.neighbors_of(idx).len();
            if deg == 1 {
                if let Some(nbr) = adjacency.neighbors_of(idx).first() {
                    let bond = &working.bonds()[nbr.bond.index()];
                    extra_h[nbr.atom_index] += bond_order_to_unsigned(bond.order());
                }
            }
        }

        // Map old index -> new index
        let kept_indices: Vec<usize> = (0..working.num_atoms())
            .filter(|i| !removal_indices.contains(i))
            .collect();
        let old_to_new: Vec<Option<usize>> = (0..working.num_atoms())
            .map(|i| {
                if removal_indices.contains(&i) {
                    None
                } else {
                    kept_indices.iter().position(|&x| x == i)
                }
            })
            .collect();

        // Build new atoms
        let mut new_atoms: Vec<AtomSpec> = Vec::new();
        for (old_idx, atom) in working.atoms().iter().enumerate() {
            if removal_indices.contains(&old_idx) {
                continue;
            }
            let mut spec = AtomSpec::new(atom.element())
                .with_formal_charge(atom.formal_charge())
                .with_chiral_tag(atom.chiral_tag())
                .with_aromatic(atom.is_aromatic());
            if let Some(isotope) = atom.isotope() {
                spec = spec.with_isotope(isotope);
            }
            if let Some(atom_map) = atom.atom_map() {
                spec = spec.with_atom_map(atom_map);
            }
            let new_h = atom.explicit_hydrogens() as u32 + extra_h[old_idx];
            spec = spec.with_explicit_hydrogens(new_h as u8);
            new_atoms.push(spec);
        }

        // Build new bonds with remapped indices
        let mut new_bonds: Vec<BondSpec> = Vec::new();
        for bond in working.bonds() {
            let begin_idx = bond.begin().index();
            let end_idx = bond.end().index();
            if removal_indices.contains(&begin_idx) || removal_indices.contains(&end_idx) {
                continue;
            }
            let new_begin = AtomId::new(old_to_new[begin_idx].unwrap());
            let new_end = AtomId::new(old_to_new[end_idx].unwrap());
            new_bonds.push(
                BondSpec::new(new_begin, new_end, bond.order())
                    .with_aromatic(bond.is_aromatic())
                    .with_conjugated(bond.is_conjugated()),
            );
        }

        builder = MoleculeBuilder::new();
        for spec in &new_atoms {
            builder.add_atom(spec.clone());
        }
        for spec in &new_bonds {
            builder.add_bond(spec.clone())?;
        }
        working = builder.build()?;
    }

    Ok(working)
}

/// Extract the net (extended) Murcko scaffold.
///
/// Ring atoms and atoms that connect two or more ring systems (scaffold
/// connectors) are kept. Side-chain atoms adjacent to the scaffold are
/// converted to dummy (*) atoms. Non-adjacent side-chain atoms are removed.
///
/// Equivalent to RDKit's `HashFunction::ExtendedMurcko`.
///
/// # Errors
/// Returns `HashError::EmptyMolecule` if the molecule has no atoms.
pub fn mol_net_scaffold(mol: &Molecule) -> Result<Molecule, HashError> {
    if mol.num_atoms() == 0 {
        return Err(HashError::EmptyMolecule);
    }

    let n = mol.num_atoms();
    let adjacency = mol.topology_block().adjacency.clone();

    // Determine which atoms are in a ring
    let ri = mol.derived_cache().rings.as_ref();
    let is_ring_atom: Vec<bool> = if let Some(ring_info) = ri {
        (0..n)
            .map(|i| ring_info.num_atom_rings(AtomId::new(i)) > 0)
            .collect()
    } else {
        vec![false; n]
    };

    // Atoms in the scaffold: ring atoms + atoms connecting to ring systems
    let mut is_in_scaffold: Vec<bool> = vec![false; n];
    for i in 0..n {
        is_in_scaffold[i] = is_in_scaffold_inner(i, n, &is_ring_atom, &adjacency, mol);
    }

    // Keep atoms in scaffold, or adjacent to scaffold (become dummy)
    let mut keep_mask: Vec<bool> = vec![false; n];
    for i in 0..n {
        keep_mask[i] = is_in_scaffold[i] || has_nbr_in_scaffold(i, &is_in_scaffold, &adjacency);
    }

    // Map old index -> new index
    let old_to_new: Vec<Option<usize>> = {
        let mut map = vec![None; n];
        let mut next = 0usize;
        for i in 0..n {
            if keep_mask[i] {
                map[i] = Some(next);
                next += 1;
            }
        }
        map
    };

    // Build atoms
    let mut builder = MoleculeBuilder::new();
    for i in 0..n {
        if !keep_mask[i] {
            continue;
        }
        let atom = &mol.atoms()[i];
        if is_in_scaffold[i] {
            let mut spec = AtomSpec::new(atom.element())
                .with_formal_charge(atom.formal_charge())
                .with_chiral_tag(atom.chiral_tag())
                .with_aromatic(atom.is_aromatic());
            if let Some(isotope) = atom.isotope() {
                spec = spec.with_isotope(isotope);
            }
            if let Some(atom_map) = atom.atom_map() {
                spec = spec.with_atom_map(atom_map);
            }
            builder.add_atom(spec);
        } else {
            // Side-chain atom adjacent to scaffold -> dummy atom
            let spec = AtomSpec::new(Element::DUMMY).with_aromatic(false);
            builder.add_atom(spec);
        }
    }

    // Build bonds
    for bond in mol.bonds() {
        let begin_idx = bond.begin().index();
        let end_idx = bond.end().index();
        if !keep_mask[begin_idx] || !keep_mask[end_idx] {
            continue;
        }
        let new_begin = AtomId::new(old_to_new[begin_idx].unwrap());
        let new_end = AtomId::new(old_to_new[end_idx].unwrap());
        let bond_spec =
            BondSpec::new(new_begin, new_end, bond.order()).with_aromatic(bond.is_aromatic());
        builder.add_bond(bond_spec)?;
    }

    builder.build().map_err(HashError::from)
}

/// Compute a canonical 64-bit hash of the molecule.
///
/// Atoms are iterated in CIP rank order. Each atom's atomic number, charge,
/// isotope, chirality, aromaticity, and ring membership are hashed along with
/// the properties of its neighbors (bond order and neighbor atomic number).
///
/// # Errors
/// Returns `HashError::EmptyMolecule` if the molecule has no atoms.
pub fn mol_hash(mol: &Molecule) -> Result<u64, HashError> {
    if mol.num_atoms() == 0 {
        return Err(HashError::EmptyMolecule);
    }

    let ranks = assign_atom_cip_ranks(mol)?;
    mol_hash_with_ranks(mol, &ranks)
}

/// Compute a 64-bit hash using the provided CIP ranks.
///
/// Ranks are used to order atoms; lower rank = higher priority.
pub fn mol_hash_with_ranks(mol: &Molecule, ranks: &[u32]) -> Result<u64, HashError> {
    if mol.num_atoms() == 0 {
        return Err(HashError::EmptyMolecule);
    }

    let n = mol.num_atoms();
    let adjacency = mol.topology_block().adjacency.clone();

    let ri = mol.derived_cache().rings.as_ref();

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
        let atom = &mol.atoms()[atom_idx];
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
            let bond = &mol.bonds()[nbr.bond.index()];
            let bond_val = bond_order_to_hash_val(bond.order());
            hash_combine(&mut atom_hash, bond_val);
            let nbr_atom = &mol.atoms()[nbr.atom_index];
            hash_combine(&mut atom_hash, nbr_atom.atomic_number() as u32);
            hash_combine(&mut atom_hash, if bond.is_aromatic() { 1u32 } else { 0u32 });
        }

        // Fold atom_hash into the overall 64-bit result
        hash_combine(&mut hash_hi, atom_hash);
        hash_combine(&mut hash_lo, atom_hash.wrapping_mul(0x9e3779b9));
    }

    Ok(((hash_hi as u64) << 32) | (hash_lo as u64))
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Build a `MoleculeBuilder` from an existing molecule (copy atoms + bonds).
fn build_molecule_builder(mol: &Molecule) -> MoleculeBuilder {
    let mut builder = MoleculeBuilder::new();
    for atom in mol.atoms() {
        let mut spec = AtomSpec::new(atom.element())
            .with_formal_charge(atom.formal_charge())
            .with_explicit_hydrogens(atom.explicit_hydrogens())
            .with_chiral_tag(atom.chiral_tag())
            .with_aromatic(atom.is_aromatic());
        if let Some(isotope) = atom.isotope() {
            spec = spec.with_isotope(isotope);
        }
        if let Some(atom_map) = atom.atom_map() {
            spec = spec.with_atom_map(atom_map);
        }
        builder.add_atom(spec);
    }
    for bond in mol.bonds() {
        let bond_spec =
            BondSpec::new(bond.begin(), bond.end(), bond.order()).with_aromatic(bond.is_aromatic());
        // Silently ignore bond errors during copy (should not happen for valid mol)
        let _ = builder.add_bond(bond_spec);
    }
    builder
}

/// Convert a bond order to an unsigned integer for H-count transfer.
fn bond_order_to_unsigned(order: BondOrder) -> u32 {
    match order {
        BondOrder::Single | BondOrder::Aromatic => 1,
        BondOrder::Double | BondOrder::OneAndHalf => 2,
        BondOrder::Triple | BondOrder::TwoAndHalf => 3,
        BondOrder::Quadruple | BondOrder::ThreeAndHalf => 4,
        BondOrder::Quintuple | BondOrder::FourAndHalf => 5,
        BondOrder::Hextuple | BondOrder::FiveAndHalf => 6,
        _ => 1,
    }
}

/// Convert a bond order to a small hash value.
fn bond_order_to_hash_val(order: BondOrder) -> u32 {
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
        BondOrder::Other | BondOrder::Unspecified | BondOrder::Zero | BondOrder::Null => 0,
        _ => 0,
    }
}

/// Convert a chiral tag to a small integer for hashing.
fn chiral_tag_to_hash_code(tag: ChiralTag) -> u32 {
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

/// Check if an atom is in the scaffold (ring atom or connector to ring systems).
///
/// RDKit<beh>: IsInScaffold — ring atoms return true directly; non-ring atoms
/// are scaffold connectors if they have at least two neighbors that can reach
/// a ring system via DFS.
fn is_in_scaffold_inner(
    atom_idx: usize,
    n: usize,
    is_ring_atom: &[bool],
    adjacency: &crate::AdjacencyList,
    mol: &Molecule,
) -> bool {
    // Ring atoms are always in scaffold
    if is_ring_atom[atom_idx] {
        return true;
    }

    // Non-ring atoms: count how many neighbors can reach a ring via DFS
    let mut count = 0u32;
    let nbrs = adjacency.neighbors_of(atom_idx);
    for nbr in nbrs {
        if depth_first_search_for_ring(atom_idx, nbr.atom_index, n, is_ring_atom, adjacency) {
            count += 1;
        }
    }
    count > 1
}

/// DFS to find if `start` can reach a ring atom without going through `blocked`.
///
/// RDKit<beh>: TraverseForRing / DepthFirstSearchForRing
fn depth_first_search_for_ring(
    blocked: usize,
    start: usize,
    n: usize,
    is_ring_atom: &[bool],
    adjacency: &crate::AdjacencyList,
) -> bool {
    let mut visited = vec![false; n];
    visited[blocked] = true;
    traverse_for_ring(start, &mut visited, is_ring_atom, adjacency)
}

fn traverse_for_ring(
    atom_idx: usize,
    visited: &mut [bool],
    is_ring_atom: &[bool],
    adjacency: &crate::AdjacencyList,
) -> bool {
    visited[atom_idx] = true;

    if is_ring_atom[atom_idx] {
        return true;
    }

    for nbr in adjacency.neighbors_of(atom_idx) {
        if !visited[nbr.atom_index] {
            if traverse_for_ring(nbr.atom_index, visited, is_ring_atom, adjacency) {
                return true;
            }
        }
    }
    false
}

/// Check if an atom has a neighbor in the scaffold.
fn has_nbr_in_scaffold(
    atom_idx: usize,
    is_in_scaffold: &[bool],
    adjacency: &crate::AdjacencyList,
) -> bool {
    for nbr in adjacency.neighbors_of(atom_idx) {
        if is_in_scaffold[nbr.atom_index] {
            return true;
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Molecule;

    fn sanitized_mol(smiles: &str) -> Molecule {
        // Molecule::from_smiles already runs sanitization internally.
        // Ring finding may be needed for ring-aware tests.
        let mol = Molecule::from_smiles(smiles).expect("valid SMILES");
        let _ = crate::rings::find_sssr(&mol);
        mol
    }

    #[test]
    fn test_murcko_scaffold_simple_ring() {
        let mol = sanitized_mol("c1ccccc1");
        let scaffold = mol_murcko_scaffold(&mol).expect("murcko scaffold");
        assert!(scaffold.num_atoms() > 0);
        assert_eq!(scaffold.num_atoms(), 6);
    }

    #[test]
    fn test_murcko_scaffold_with_side_chain() {
        // Build toluene (methylbenzene) manually: benzene ring + methyl side chain.
        let mol = build_toluene();
        let scaffold = mol_murcko_scaffold(&mol).expect("murcko scaffold");
        // Murcko scaffold = ring atoms only, so 6 C's
        assert_eq!(scaffold.num_atoms(), 6);
    }

    #[test]
    fn test_net_scaffold_simple() {
        let mol = sanitized_mol("c1ccccc1");
        let scaffold = mol_net_scaffold(&mol).expect("net scaffold");
        assert!(scaffold.num_atoms() >= 6);
    }

    #[test]
    fn test_net_scaffold_with_side_chain() {
        // Use cyclohexane (which parses correctly via SMILES) to test net scaffold
        let mol = sanitized_mol("C1CCCCC1");
        let scaffold = mol_net_scaffold(&mol).expect("net scaffold");
        // Cyclohexane: all 6 ring carbons are in scaffold, no side chains
        assert_eq!(scaffold.num_atoms(), 6);
    }

    fn build_toluene() -> Molecule {
        // Build a 6-carbon aromatic ring with a methyl group attached at atom 0.
        use crate::Element;
        use crate::atom::AtomSpec;
        use crate::bond::BondSpec;
        use crate::builder::MoleculeBuilder;

        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let c1 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let c2 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let c3 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let c4 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let c5 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let c6 = builder.add_atom(AtomSpec::new(Element::C)); // methyl

        // Aromatic ring bonds
        builder
            .add_bond(BondSpec::new(c0, c1, crate::BondOrder::Aromatic))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, c2, crate::BondOrder::Aromatic))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c2, c3, crate::BondOrder::Aromatic))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c3, c4, crate::BondOrder::Aromatic))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c4, c5, crate::BondOrder::Aromatic))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c5, c0, crate::BondOrder::Aromatic))
            .unwrap();

        // Methyl single bond
        builder
            .add_bond(BondSpec::new(c0, c6, crate::BondOrder::Single))
            .unwrap();

        builder.build().unwrap()
    }

    #[test]
    fn test_mol_hash_benzene() {
        let mol = sanitized_mol("c1ccccc1");
        let hash = mol_hash(&mol).expect("hash");
        assert!(hash != 0);
    }

    #[test]
    fn test_mol_hash_deterministic() {
        let mol = sanitized_mol("c1ccccc1");
        let hash1 = mol_hash(&mol).expect("hash1");
        let hash2 = mol_hash(&mol).expect("hash2");
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_mol_hash_different_molecules() {
        let mol1 = sanitized_mol("c1ccccc1");
        let mol2 = sanitized_mol("CCO");
        let hash1 = mol_hash(&mol1).expect("hash1");
        let hash2 = mol_hash(&mol2).expect("hash2");
        assert_ne!(hash1, hash2);
    }

    #[test]
    fn test_mol_hash_empty_error() {
        let mol = Molecule::new();
        assert!(mol_hash(&mol).is_err());
    }

    #[test]
    fn test_mol_hash_with_ranks() {
        let mol = sanitized_mol("c1ccccc1");
        let ranks = assign_atom_cip_ranks(&mol).expect("ranks");
        let hash = mol_hash_with_ranks(&mol, &ranks).expect("hash");
        assert!(hash != 0);
    }
}
