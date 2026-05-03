use std::collections::{HashMap, VecDeque};

use crate::{
    AdjacencyList, Bond, BondDirection, BondOrder, ChiralTag, Molecule, ValenceModel,
    assign_valence,
};
use glam::DVec3;

/// A ligand participating in a tetrahedral stereocenter.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum LigandRef {
    /// A graph atom by molecule atom index.
    Atom(usize),
    /// A hydrogen carried as an atom hydrogen count rather than a graph atom.
    ImplicitH,
}

/// Ordered tetrahedral stereochemistry for one atom center.
///
/// Specification: `tetrahedral_stereo_representation.md` at repository root.
///
/// For RDKit-style tetrahedral tags, the ordering is derived from the atom's
/// incoming bond order and chiral tag.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub center: usize,
    pub ligands: [LigandRef; 4],
}

/// Hidden helper struct exposing RDKit legacy stereochemistry perception
/// fields for parity tests.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegacyStereoAtomProps {
    pub cip_code: Option<String>,
    pub cip_rank: Option<i64>,
    pub chirality_possible: bool,
}

fn bond_affects_atom_chirality(bond: &Bond, atom_index: usize) -> bool {
    !(matches!(bond.order, BondOrder::Null)
        || matches!(bond.order, BondOrder::Dative) && bond.begin_atom == atom_index)
}

fn atom_nonzero_degree(mol: &Molecule, atom_index: usize) -> usize {
    mol.bonds()
        .iter()
        .filter(|bond| {
            (bond.begin_atom == atom_index || bond.end_atom == atom_index)
                && bond_affects_atom_chirality(bond, atom_index)
        })
        .count()
}

fn atom_nonzero_degree_with_adjacency(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    atom_index: usize,
) -> usize {
    adjacency
        .neighbors_of(atom_index)
        .iter()
        .filter(|neighbor| {
            bond_affects_atom_chirality(&mol.bonds()[neighbor.bond_index], atom_index)
        })
        .count()
}

fn has_protium_neighbor_with_adjacency(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    atom_index: usize,
) -> bool {
    adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
        let atom = &mol.atoms()[neighbor.atom_index];
        atom.atomic_num == 1 && atom.isotope.is_none()
    })
}

fn bond_is_conjugated_for_stereo_with_adjacency(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    bond_index: usize,
) -> bool {
    let bond = &mol.bonds()[bond_index];
    if bond.is_aromatic || matches!(bond.order, BondOrder::Double | BondOrder::Triple) {
        return true;
    }
    if !matches!(bond.order, BondOrder::Single) {
        return false;
    }
    [bond.begin_atom, bond.end_atom].into_iter().any(|center| {
        adjacency.neighbors_of(center).iter().any(|neighbor| {
            neighbor.bond_index != bond_index && {
                let other = &mol.bonds()[neighbor.bond_index];
                other.is_aromatic || matches!(other.order, BondOrder::Double | BondOrder::Triple)
            }
        })
    })
}

fn atom_has_conjugated_bond_with_adjacency(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    atom_index: usize,
) -> bool {
    adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
        bond_is_conjugated_for_stereo_with_adjacency(mol, adjacency, neighbor.bond_index)
    })
}

fn sub3(lhs: DVec3, rhs: DVec3) -> DVec3 {
    DVec3::new(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z)
}

fn unit_or_zero(v: DVec3) -> DVec3 {
    let len = v.length();
    if len == 0.0 { DVec3::ZERO } else { v / len }
}

fn graph_h_neighbor_count_with_adjacency(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    atom_index: usize,
) -> usize {
    adjacency
        .neighbors_of(atom_index)
        .iter()
        .filter(|neighbor| mol.atoms()[neighbor.atom_index].atomic_num == 1)
        .count()
}

fn bond_valence_for_total_hs(bond: &Bond, atom_index: usize) -> usize {
    if !bond_affects_atom_chirality(bond, atom_index) {
        return 0;
    }
    match bond.order {
        BondOrder::Single | BondOrder::Aromatic | BondOrder::Dative => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Null => 0,
    }
}

fn rdkit_total_num_hs_for_3d_with_adjacency(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    atom_index: usize,
    implicit_hydrogens: &[u8],
) -> usize {
    let graph_hs = graph_h_neighbor_count_with_adjacency(mol, adjacency, atom_index);
    let assigned = mol.atoms()[atom_index].explicit_hydrogens as usize
        + implicit_hydrogens[atom_index] as usize
        + graph_hs;
    if assigned != 0 {
        return assigned;
    }
    let atom = &mol.atoms()[atom_index];
    if atom.no_implicit || atom.formal_charge != 0 || atom.is_aromatic {
        return graph_hs;
    }
    let explicit_valence = adjacency
        .neighbors_of(atom_index)
        .iter()
        .map(|neighbor| bond_valence_for_total_hs(&mol.bonds()[neighbor.bond_index], atom_index))
        .sum::<usize>();
    let default_valence: usize = match atom.atomic_num {
        6 => 4,
        7 | 15 | 33 => 3,
        8 | 16 | 34 => 2,
        _ => return graph_hs,
    };
    graph_hs + default_valence.saturating_sub(explicit_valence)
}

pub(crate) fn assign_chiral_types_from_3d_rdkit_subset(mol: &mut Molecule) {
    // Source mapping: RDKit Chirality.cpp::assignChiralTypesFrom3D().
    const ZERO_VOLUME_TOL: f64 = 0.1;
    let Some(coords) = mol.coords_3d() else {
        return;
    };
    if coords.len() != mol.atoms().len() {
        return;
    }
    let coords = coords.to_vec();
    let Ok(assignment) = assign_valence(mol, ValenceModel::RdkitLike) else {
        return;
    };
    let adjacency = AdjacencyList::from_topology(mol.atoms().len(), mol.bonds());
    for atom_index in 0..mol.atoms().len() {
        mol.atoms_mut()[atom_index].chiral_tag = ChiralTag::Unspecified;

        let nz_degree = atom_nonzero_degree_with_adjacency(mol, &adjacency, atom_index);
        let total_num_hs = rdkit_total_num_hs_for_3d_with_adjacency(
            mol,
            &adjacency,
            atom_index,
            &assignment.implicit_hydrogens,
        );
        let tnz_degree = nz_degree + total_num_hs;

        if nz_degree < 3 || tnz_degree > 6 {
            continue;
        }
        if assign_nontetrahedral_chiral_type_from_3d_with_adjacency(
            mol, &adjacency, atom_index, &coords,
        ) {
            continue;
        }
        if tnz_degree > 4 {
            continue;
        }
        let atomic_num = mol.atoms()[atom_index].atomic_num;
        if matches!(atomic_num, 7 | 15 | 33) {
            continue;
        }
        if !matches!(atomic_num, 16 | 34) && tnz_degree != 4 {
            continue;
        }

        let center = coords[atom_index];
        let mut nbrs = Vec::with_capacity(4);
        for neighbor in adjacency.neighbors_of(atom_index) {
            let bond = &mol.bonds()[neighbor.bond_index];
            if !bond_affects_atom_chirality(bond, atom_index) {
                continue;
            }
            nbrs.push(coords[neighbor.atom_index]);
        }
        if nbrs.len() < 3 {
            continue;
        }

        let v1 = sub3(nbrs[0], center);
        let v2 = sub3(nbrs[1], center);
        let v3 = sub3(nbrs[2], center);
        let mut chiral_vol = v1.dot(v2.cross(v3));
        if chiral_vol < -ZERO_VOLUME_TOL {
            mol.atoms_mut()[atom_index].chiral_tag = ChiralTag::TetrahedralCw;
        } else if chiral_vol > ZERO_VOLUME_TOL {
            mol.atoms_mut()[atom_index].chiral_tag = ChiralTag::TetrahedralCcw;
        } else if nbrs.len() == 4 {
            let v4 = sub3(nbrs[3], center);
            chiral_vol = -v1.dot(v2.cross(v4));
            if chiral_vol < -ZERO_VOLUME_TOL {
                mol.atoms_mut()[atom_index].chiral_tag = ChiralTag::TetrahedralCw;
            } else if chiral_vol > ZERO_VOLUME_TOL {
                mol.atoms_mut()[atom_index].chiral_tag = ChiralTag::TetrahedralCcw;
            }
        }
    }

    let legacy_props = mol.rdkit_legacy_stereo_atom_props(false);
    let cleanup_ranks = legacy_props
        .iter()
        .map(|props| props.cip_rank.unwrap_or(0))
        .collect::<Vec<_>>();
    let has_cip_code = legacy_props
        .iter()
        .map(|props| props.cip_code.is_some())
        .collect::<Vec<_>>();
    let mut cycle_cache = StereoBondCycleCache::new(mol.bonds().len());
    for (atom_index, legacy_prop) in legacy_props.iter().enumerate() {
        if !matches!(
            mol.atoms()[atom_index].chiral_tag,
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            continue;
        }
        if legacy_prop.cip_code.is_some() {
            continue;
        }
        let has_ring_stereo_partner = has_ring_stereo_partner_with_cache(
            mol,
            &adjacency,
            &mut cycle_cache,
            atom_index,
            &cleanup_ranks,
            &has_cip_code,
        );
        if has_ring_stereo_partner {
            continue;
        }
        mol.atoms_mut()[atom_index].chiral_tag = ChiralTag::Unspecified;
        if mol.atoms()[atom_index].explicit_hydrogens == 1
            && mol.atoms()[atom_index].formal_charge == 0
            && !mol.atoms()[atom_index].is_aromatic
        {
            mol.atoms_mut()[atom_index].explicit_hydrogens = 0;
            mol.atoms_mut()[atom_index].no_implicit = false;
        }
    }
}

fn assign_nontetrahedral_chiral_type_from_3d_with_adjacency(
    mol: &mut Molecule,
    adjacency: &AdjacencyList,
    atom_index: usize,
    coords: &[DVec3],
) -> bool {
    // Source mapping: RDKit Chirality.cpp::assignNontetrahedralChiralTypeFrom3D(),
    // count=5 trigonal-bipyramidal branch.
    if mol.atoms()[atom_index].atomic_num < 15 {
        return false;
    }
    if adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
        matches!(
            mol.bonds()[neighbor.bond_index].stereo,
            crate::BondStereo::Any
        )
    }) {
        return false;
    }

    let center = coords[atom_index];
    let mut vectors = Vec::new();
    for neighbor in adjacency.neighbors_of(atom_index) {
        if vectors.len() == 6 {
            return false;
        }
        let direction = coords[neighbor.atom_index] - center;
        let length = direction.length();
        if length == 0.0 {
            return false;
        }
        vectors.push(direction / length);
    }
    if vectors.len() != 4 && vectors.len() != 5 {
        return false;
    }

    assign_nontetrahedral_chiral_type_from_vectors(mol, atom_index, &vectors)
}

fn assign_nontetrahedral_chiral_type_from_vectors(
    mol: &mut Molecule,
    atom_index: usize,
    vectors: &[DVec3],
) -> bool {
    let mut pair = [0usize; 6];
    let mut pairs = 0usize;
    for i in 0..vectors.len() {
        for j in (i + 1)..vectors.len() {
            if vectors[i].dot(vectors[j]) < -(1.0 - 0.1) {
                if pair[i] != 0 || pair[j] != 0 {
                    return false;
                }
                pair[i] = j + 1;
                pair[j] = i + 1;
                pairs += 1;
            }
        }
    }
    if pairs != 1 {
        return false;
    }

    let voltest =
        |x: usize, y: usize, z: usize| vectors[x].dot(vectors[y].cross(vectors[z])) >= 0.0;
    let perm = if vectors.len() == 4 {
        if pair[0] == 2 {
            if angle_between(vectors[2], vectors[3]) < 100.0_f64.to_radians() {
                return false;
            }
            if voltest(0, 2, 3) { 7 } else { 8 }
        } else if pair[0] == 3 {
            if angle_between(vectors[1], vectors[3]) < 100.0_f64.to_radians() {
                return false;
            }
            if voltest(0, 1, 3) { 5 } else { 6 }
        } else if pair[0] == 4 {
            if angle_between(vectors[1], vectors[2]) < 100.0_f64.to_radians() {
                return false;
            }
            if voltest(0, 1, 2) { 3 } else { 4 }
        } else if pair[1] == 3 {
            if angle_between(vectors[0], vectors[3]) < 100.0_f64.to_radians() {
                return false;
            }
            if voltest(1, 0, 3) { 13 } else { 14 }
        } else if pair[1] == 4 {
            if angle_between(vectors[0], vectors[2]) < 100.0_f64.to_radians() {
                return false;
            }
            if voltest(1, 0, 2) { 10 } else { 12 }
        } else {
            if angle_between(vectors[0], vectors[1]) < 100.0_f64.to_radians() {
                return false;
            }
            if voltest(3, 0, 1) { 16 } else { 19 }
        }
    } else if pair[0] == 2 {
        if voltest(0, 2, 3) { 7 } else { 8 }
    } else if pair[0] == 3 {
        if voltest(0, 1, 3) { 5 } else { 6 }
    } else if pair[0] == 4 {
        if voltest(0, 1, 2) { 3 } else { 4 }
    } else if pair[0] == 5 {
        if voltest(0, 1, 2) { 1 } else { 2 }
    } else if pair[1] == 3 {
        if voltest(1, 0, 3) { 13 } else { 14 }
    } else if pair[1] == 4 {
        if voltest(1, 0, 2) { 10 } else { 12 }
    } else if pair[1] == 5 {
        if voltest(1, 0, 2) { 9 } else { 11 }
    } else if pair[2] == 4 {
        if voltest(2, 0, 1) { 16 } else { 19 }
    } else if pair[2] == 5 {
        if voltest(2, 0, 1) { 15 } else { 20 }
    } else if voltest(3, 0, 1) {
        17
    } else {
        18
    };
    mol.atoms_mut()[atom_index].chiral_tag = ChiralTag::TrigonalBipyramidal;
    mol.atoms_mut()[atom_index]
        .props
        .insert("_chiralPermutation".to_owned(), perm.to_string());
    true
}

fn angle_between(a: DVec3, b: DVec3) -> f64 {
    a.dot(b).clamp(-1.0, 1.0).acos()
}

struct StereoBondCycleCache {
    min_cycle_sizes: Vec<Option<Option<usize>>>,
}

impl StereoBondCycleCache {
    fn new(bond_count: usize) -> Self {
        Self {
            min_cycle_sizes: vec![None; bond_count],
        }
    }

    fn min_cycle_size_for_bond(
        &mut self,
        mol: &Molecule,
        adjacency: &AdjacencyList,
        bond_index: usize,
    ) -> Option<usize> {
        if let Some(cached) = self.min_cycle_sizes[bond_index] {
            return cached;
        }
        let result = min_cycle_size_for_bond_with_adjacency(mol, adjacency, bond_index);
        self.min_cycle_sizes[bond_index] = Some(result);
        result
    }
}

fn min_cycle_size_for_bond_with_adjacency(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    bond_index: usize,
) -> Option<usize> {
    let bond = &mol.bonds()[bond_index];
    let start = bond.begin_atom;
    let goal = bond.end_atom;
    let mut queue = VecDeque::from([(start, 0usize)]);
    let mut seen = vec![false; mol.atoms().len()];
    seen[start] = true;
    while let Some((atom, dist)) = queue.pop_front() {
        for neighbor in adjacency.neighbors_of(atom) {
            if neighbor.bond_index == bond_index {
                continue;
            }
            let next_atom = neighbor.atom_index;
            if next_atom == goal {
                return Some(dist + 2);
            }
            if !seen[next_atom] {
                seen[next_atom] = true;
                queue.push_back((next_atom, dist + 1));
            }
        }
    }
    None
}

fn has_ring_stereo_partner_with_cache(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    cycle_cache: &mut StereoBondCycleCache,
    atom_index: usize,
    ranks: &[i64],
    has_cip_code: &[bool],
) -> bool {
    if !atom_is_candidate_for_ring_stereochem_with_cache(
        mol,
        adjacency,
        cycle_cache,
        atom_index,
        ranks,
    ) {
        return false;
    }
    let atom_ring_neighbor_count = adjacency
        .neighbors_of(atom_index)
        .iter()
        .filter(|neighbor| {
            cycle_cache
                .min_cycle_size_for_bond(mol, adjacency, neighbor.bond_index)
                .is_some()
        })
        .count();
    let mut seen_atoms = vec![false; mol.atoms().len()];
    let mut seen_bonds = vec![false; mol.bonds().len()];
    let mut queue = VecDeque::new();
    seen_atoms[atom_index] = true;
    for neighbor in adjacency.neighbors_of(atom_index) {
        seen_bonds[neighbor.bond_index] = true;
        if cycle_cache
            .min_cycle_size_for_bond(mol, adjacency, neighbor.bond_index)
            .is_some()
        {
            let other = neighbor.atom_index;
            if !seen_atoms[other] {
                seen_atoms[other] = true;
                queue.push_back(other);
            }
        }
    }
    while let Some(current) = queue.pop_front() {
        if !matches!(mol.atoms()[current].chiral_tag, ChiralTag::Unspecified)
            && !has_cip_code.get(current).copied().unwrap_or(false)
            && (atom_ring_neighbor_count == 4
                || atom_is_candidate_for_ring_stereochem_with_cache(
                    mol,
                    adjacency,
                    cycle_cache,
                    current,
                    ranks,
                ))
        {
            return true;
        }
        for neighbor in adjacency.neighbors_of(current) {
            if seen_bonds[neighbor.bond_index] {
                continue;
            }
            seen_bonds[neighbor.bond_index] = true;
            if cycle_cache
                .min_cycle_size_for_bond(mol, adjacency, neighbor.bond_index)
                .is_some()
            {
                let other = neighbor.atom_index;
                if !seen_atoms[other] {
                    seen_atoms[other] = true;
                    queue.push_back(other);
                }
            }
        }
    }
    false
}

fn atom_is_candidate_for_ring_stereochem_with_cache(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    cycle_cache: &mut StereoBondCycleCache,
    atom_index: usize,
    ranks: &[i64],
) -> bool {
    if !adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
        cycle_cache
            .min_cycle_size_for_bond(mol, adjacency, neighbor.bond_index)
            .is_some()
    }) {
        return false;
    }
    let atom = &mol.atoms()[atom_index];
    if atom.atomic_num == 7
        && atom_nonzero_degree_with_adjacency(mol, adjacency, atom_index) == 3
        && !adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
            cycle_cache
                .min_cycle_size_for_bond(mol, adjacency, neighbor.bond_index)
                .is_some_and(|size| size == 3)
        })
        && !is_atom_bridgehead_with_cache(mol, adjacency, cycle_cache, atom_index)
    {
        return false;
    }

    let mut non_ring_neighbors = Vec::new();
    let mut ring_neighbor_ranks = std::collections::BTreeSet::new();
    let mut ring_neighbor_count = 0usize;
    for neighbor in adjacency.neighbors_of(atom_index) {
        if cycle_cache
            .min_cycle_size_for_bond(mol, adjacency, neighbor.bond_index)
            .is_some()
        {
            ring_neighbor_count += 1;
            ring_neighbor_ranks.insert(ranks[neighbor.atom_index]);
        } else {
            non_ring_neighbors.push(neighbor.atom_index);
        }
    }
    match non_ring_neighbors.len() {
        2 => {
            ranks[non_ring_neighbors[0]] != ranks[non_ring_neighbors[1]]
                && ring_neighbor_count != ring_neighbor_ranks.len()
        }
        1 => ring_neighbor_count > ring_neighbor_ranks.len(),
        0 => {
            (ring_neighbor_count == 4 && ring_neighbor_ranks.len() == 3)
                || (ring_neighbor_count == 3 && ring_neighbor_ranks.len() == 2)
        }
        _ => false,
    }
}

fn is_atom_bridgehead_with_cache(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    cycle_cache: &mut StereoBondCycleCache,
    atom_index: usize,
) -> bool {
    if atom_nonzero_degree_with_adjacency(mol, adjacency, atom_index) < 3 {
        return false;
    }
    let atom_ring_bonds = adjacency
        .neighbors_of(atom_index)
        .iter()
        .filter(|neighbor| {
            cycle_cache
                .min_cycle_size_for_bond(mol, adjacency, neighbor.bond_index)
                .is_some()
        })
        .count();
    atom_ring_bonds >= 3
}

pub(crate) fn assign_chiral_types_from_bond_dirs_rdkit_subset(mol: &mut Molecule) {
    // Source mapping: RDKit Chirality.cpp::assignChiralTypesFromBondDirs()
    // tetrahedral pseudo-3D branch. This covers molfile wedge/dash input and
    // deliberately leaves unsupported/ambiguous layouts unspecified.
    const ZERO_TOL: f64 = 1e-3;
    const COORD_ZERO_TOL: f64 = 1e-4;
    const T_SHAPE_TOL: f64 = 0.00031;
    const PSEUDO_3D_OFFSET: f64 = 0.1;
    const VOLUME_TOLERANCE: f64 = 0.00174;

    let Some(coords_2d) = mol.coords_2d() else {
        return;
    };
    if coords_2d.len() != mol.atoms().len() {
        return;
    }
    let coords = coords_2d
        .iter()
        .map(|coord| DVec3::new(coord.x, coord.y, 0.0))
        .collect::<Vec<_>>();
    let mut atoms_set = vec![false; mol.atoms().len()];

    for bond_index in 0..mol.bonds().len() {
        let direction = mol.bonds()[bond_index].direction;
        if !matches!(
            direction,
            BondDirection::EndUpRight | BondDirection::EndDownRight
        ) {
            continue;
        }
        let atom_index = mol.bonds()[bond_index].begin_atom;
        if atoms_set[atom_index] || atom_nonzero_degree(mol, atom_index) > 4 {
            continue;
        }

        let center = coords[atom_index];
        let wedged_atom = mol.bonds()[bond_index].end_atom;
        let ref_length = (center - coords[wedged_atom]).length();
        let z_offset = if matches!(direction, BondDirection::EndUpRight) {
            PSEUDO_3D_OFFSET
        } else {
            -PSEUDO_3D_OFFSET
        } * if ref_length == 0.0 { 1.0 } else { ref_length };

        let mut ref_idx = usize::MAX;
        let mut bond_vects = Vec::<DVec3>::new();
        let mut all_single = true;
        for (neighbor_idx, neighbor_bond) in mol
            .bonds()
            .iter()
            .filter(|bond| bond.begin_atom == atom_index || bond.end_atom == atom_index)
            .enumerate()
        {
            let other = if neighbor_bond.begin_atom == atom_index {
                neighbor_bond.end_atom
            } else {
                neighbor_bond.begin_atom
            };
            let mut point = coords[other];
            if neighbor_bond.index == bond_index {
                ref_idx = neighbor_idx;
                point.z = z_offset;
            } else if neighbor_bond.begin_atom == atom_index
                && matches!(
                    neighbor_bond.direction,
                    BondDirection::EndUpRight | BondDirection::EndDownRight
                )
            {
                point.z = if matches!(neighbor_bond.direction, BondDirection::EndUpRight) {
                    PSEUDO_3D_OFFSET
                } else {
                    -PSEUDO_3D_OFFSET
                } * if ref_length == 0.0 { 1.0 } else { ref_length };
            }
            if neighbor_bond.index != bond_index && (center - point).length_squared() < ZERO_TOL {
                ref_idx = usize::MAX;
                break;
            }
            if !matches!(neighbor_bond.order, BondOrder::Single) {
                all_single = false;
            }
            bond_vects.push(unit_or_zero(point - center));
        }
        if ref_idx == usize::MAX {
            continue;
        }
        let n_nbrs = bond_vects.len();
        if !(3..=4).contains(&n_nbrs) {
            continue;
        }
        if !all_single && !matches!(mol.atoms()[atom_index].atomic_num, 15 | 16) {
            continue;
        }

        let mut order = [0usize, 1, 2, 3];
        let mut prefactor = 1.0;
        if ref_idx != 0 {
            order.swap(0, ref_idx);
            prefactor *= -1.0;
        }
        if n_nbrs > 3
            && bond_vects[order[1]]
                .cross(bond_vects[order[2]])
                .length_squared()
                < 10.0 * ZERO_TOL
            && bond_vects[order[1]]
                .cross(bond_vects[order[0]])
                .length_squared()
                > 10.0 * ZERO_TOL
        {
            bond_vects[order[1]].z = -bond_vects[order[0]].z;
        }
        if n_nbrs == 3 {
            let cp01 = bond_vects[order[0]].cross(bond_vects[order[1]]);
            let cp02 = bond_vects[order[0]].cross(bond_vects[order[2]]);
            let dp01 = bond_vects[order[0]].dot(bond_vects[order[1]]);
            let dp02 = bond_vects[order[0]].dot(bond_vects[order[2]]);
            if needs_swap_pseudo_3d_rdkit_subset(cp01, cp02, dp01, dp02, ZERO_TOL) {
                order.swap(1, 2);
                prefactor *= -1.0;
            }
        } else {
            let mut ordered = Vec::with_capacity(3);
            for i in 1..4 {
                let cp0i = bond_vects[order[0]].cross(bond_vects[order[i]]);
                let sgn = if cp0i.z < -ZERO_TOL { -1.0 } else { 1.0 };
                let dp0i = bond_vects[order[0]].dot(bond_vects[order[i]]);
                ordered.push((sgn, sgn * dp0i, order[i]));
            }
            ordered.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
            let mut n_changed = 0;
            for i in 1..4 {
                let next = ordered[i - 1].2;
                if order[i] != next {
                    order[i] = next;
                    n_changed += 1;
                }
            }
            if n_changed == 2 {
                prefactor *= -1.0;
            }
        }

        let mut opposing_conflict = false;
        'opposing: for i in 0..n_nbrs {
            for j in (i + 1)..n_nbrs {
                if bond_vects[order[i]].z * bond_vects[order[j]].z < -ZERO_TOL {
                    let cp = bond_vects[order[i]]
                        .cross(bond_vects[order[j]])
                        .length_squared();
                    if cp < 0.01 {
                        if n_nbrs == 4
                            && (bond_vects[order[i]].dot(bond_vects[order[j]]) + 1.0).abs()
                                < ZERO_TOL
                            && (j - i == 1 || (i == 0 && j == 3))
                        {
                            bond_vects[order[j]].z = 0.0;
                            continue;
                        }
                        opposing_conflict = true;
                        break 'opposing;
                    }
                }
            }
        }
        if opposing_conflict {
            continue;
        }

        if n_nbrs == 3 {
            let conflict = if bond_vects[order[1]].z * bond_vects[order[0]].z < -COORD_ZERO_TOL
                && bond_vects[order[2]].z.abs() < COORD_ZERO_TOL
            {
                bond_vects[order[2]].cross(bond_vects[order[0]]).z
                    * bond_vects[order[2]].cross(bond_vects[order[1]]).z
                    < -1e-4
            } else if bond_vects[order[2]].z * bond_vects[order[0]].z < -COORD_ZERO_TOL
                && bond_vects[order[1]].z.abs() < COORD_ZERO_TOL
            {
                bond_vects[order[1]].cross(bond_vects[order[0]]).z
                    * bond_vects[order[1]].cross(bond_vects[order[2]]).z
                    < -COORD_ZERO_TOL
            } else {
                false
            };
            if conflict {
                continue;
            }
        }

        let mut bv1 = bond_vects[order[1]];
        bv1.z = 0.0;
        let mut bv2 = bond_vects[order[2]];
        bv2.z = 0.0;
        let mut crossp1 = bv1.cross(bv2);
        if n_nbrs == 3 {
            if crossp1.length_squared() < T_SHAPE_TOL {
                bv1.z = -bond_vects[order[0]].z;
                bv2.z = -bond_vects[order[0]].z;
                crossp1 = bv1.cross(bv2);
            }
        } else if crossp1.length_squared() < 10.0 * ZERO_TOL
            && bond_vects[order[3]].z.abs() < COORD_ZERO_TOL
        {
            bond_vects[order[3]].z = -bond_vects[order[0]].z;
        }
        let mut vol = crossp1.dot(bond_vects[order[0]]);
        if n_nbrs == 4 {
            let dotp1 = bond_vects[order[1]].dot(bond_vects[order[2]]);
            let mut bv3 = bond_vects[order[3]];
            bv3.z = 0.0;
            let crossp2 = bv1.cross(bv3);
            let dotp2 = bond_vects[order[1]].dot(bond_vects[order[3]]);
            let vol2 = crossp2.dot(bond_vects[order[0]]);
            if vol.abs() < ZERO_TOL {
                if vol2.abs() < ZERO_TOL {
                    continue;
                }
                vol = vol2;
                prefactor *= -1.0;
            } else if vol * vol2 > 0.0 && vol2.abs() > VOLUME_TOLERANCE && dotp1 < dotp2 {
                vol = vol2;
                prefactor *= -1.0;
            } else if vol.abs() < VOLUME_TOLERANCE && vol2.abs() > VOLUME_TOLERANCE {
                if vol * vol2 < 0.0 {
                    prefactor *= -1.0;
                }
                vol = vol2;
            }
        }
        vol *= prefactor;
        if vol > VOLUME_TOLERANCE {
            mol.atoms_mut()[atom_index].chiral_tag = ChiralTag::TetrahedralCcw;
            atoms_set[atom_index] = true;
        } else if vol < -VOLUME_TOLERANCE {
            mol.atoms_mut()[atom_index].chiral_tag = ChiralTag::TetrahedralCw;
            atoms_set[atom_index] = true;
        }
    }
}

fn needs_swap_pseudo_3d_rdkit_subset(
    cp01: DVec3,
    cp02: DVec3,
    dp01: f64,
    dp02: f64,
    zero_tol: f64,
) -> bool {
    if dp01.abs() - 1.0 > -zero_tol {
        return cp02.z < 0.0;
    }
    if dp02.abs() - 1.0 > -zero_tol && cp01.z < 0.0 {
        return true;
    }
    if cp01.z * cp02.z < -zero_tol {
        return cp01.z < cp02.z;
    }
    if dp01 * dp02 < -zero_tol {
        return dp01 < dp02;
    }
    dp01.abs() > dp02.abs()
}

fn count_swaps_to_interconvert(probe: &[usize], reference: &[usize]) -> usize {
    let mut probe = probe.to_vec();
    let mut positions: HashMap<usize, usize> =
        probe.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let mut swaps = 0usize;
    for (i, &target) in reference.iter().enumerate() {
        if probe[i] == target {
            continue;
        }
        let j = *positions
            .get(&target)
            .expect("reference element missing from probe permutation");
        let current = probe[i];
        probe.swap(i, j);
        positions.insert(current, j);
        positions.insert(target, i);
        swaps += 1;
    }
    swaps
}

fn perturbation_order_for_atom_with_adjacency(
    adjacency: &AdjacencyList,
    atom_index: usize,
    probe: &[usize],
) -> usize {
    let reference = adjacency
        .neighbors_of(atom_index)
        .iter()
        .map(|neighbor| neighbor.bond_index)
        .collect::<Vec<_>>();
    count_swaps_to_interconvert(probe, &reference)
}

pub(crate) fn min_cycle_size_for_bond(mol: &Molecule, bond_index: usize) -> Option<usize> {
    let bond = &mol.bonds()[bond_index];
    let start = bond.begin_atom;
    let goal = bond.end_atom;
    let mut queue = VecDeque::from([(start, 0usize)]);
    let mut seen = vec![false; mol.atoms().len()];
    seen[start] = true;
    while let Some((atom, dist)) = queue.pop_front() {
        for (idx, edge) in mol.bonds().iter().enumerate() {
            if idx == bond_index {
                continue;
            }
            let next = if edge.begin_atom == atom {
                Some(edge.end_atom)
            } else if edge.end_atom == atom {
                Some(edge.begin_atom)
            } else {
                None
            };
            let Some(next_atom) = next else {
                continue;
            };
            if next_atom == goal {
                return Some(dist + 2);
            }
            if !seen[next_atom] {
                seen[next_atom] = true;
                queue.push_back((next_atom, dist + 1));
            }
        }
    }
    None
}

pub(crate) fn should_detect_double_bond_stereo(mol: &Molecule, bond_index: usize) -> bool {
    min_cycle_size_for_bond(mol, bond_index).is_none_or(|size| size >= 8)
}

fn atom_has_directional_bond_to_other_than_with_adjacency(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    atom_index: usize,
    excluded_neighbor: usize,
) -> bool {
    adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
        neighbor.atom_index != excluded_neighbor
            && matches!(
                mol.bonds()[neighbor.bond_index].direction,
                crate::BondDirection::EndDownRight | crate::BondDirection::EndUpRight
            )
    })
}

fn legacy_assign_bond_stereo_would_leave_unassigned_with_cache(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    cycle_cache: &mut StereoBondCycleCache,
) -> bool {
    // Source mapping: RDKit Chirality.cpp::assignBondStereoCodes(). During
    // clean legacy stereo perception RDKit clears double-bond stereo first,
    // then counts eligible double bonds that cannot be assigned from adjacent
    // directional single bonds.
    for (bond_index, bond) in mol.bonds().iter().enumerate() {
        if !matches!(bond.order, BondOrder::Double)
            || cycle_cache
                .min_cycle_size_for_bond(mol, adjacency, bond_index)
                .is_some_and(|size| size < 8)
        {
            continue;
        }
        let begin_degree = atom_nonzero_degree_with_adjacency(mol, adjacency, bond.begin_atom);
        let end_degree = atom_nonzero_degree_with_adjacency(mol, adjacency, bond.end_atom);
        if !matches!(begin_degree, 2 | 3) || !matches!(end_degree, 2 | 3) {
            continue;
        }
        let begin_has_dir = atom_has_directional_bond_to_other_than_with_adjacency(
            mol,
            adjacency,
            bond.begin_atom,
            bond.end_atom,
        );
        let end_has_dir = atom_has_directional_bond_to_other_than_with_adjacency(
            mol,
            adjacency,
            bond.end_atom,
            bond.begin_atom,
        );
        if !(begin_has_dir && end_has_dir) {
            return true;
        }
    }
    false
}

pub(crate) fn is_atom_bridgehead(mol: &Molecule, atom_index: usize) -> bool {
    if atom_nonzero_degree(mol, atom_index) < 3 {
        return false;
    }
    let atom_ring_bonds = mol
        .bonds()
        .iter()
        .filter(|bond| bond.begin_atom == atom_index || bond.end_atom == atom_index)
        .filter(|bond| min_cycle_size_for_bond(mol, bond.index).is_some())
        .count();
    atom_ring_bonds >= 3
}

fn is_atom_potential_chiral_center_with_cache(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    cycle_cache: &mut StereoBondCycleCache,
    atom_index: usize,
    ranks: &[i64],
    explicit_valence: &[u8],
    implicit_hydrogens: &[u8],
) -> (bool, bool, Vec<(i64, usize)>) {
    let atom = &mol.atoms()[atom_index];
    let nz_degree = atom_nonzero_degree_with_adjacency(mol, adjacency, atom_index);
    let total_num_hs = atom.explicit_hydrogens as usize
        + implicit_hydrogens[atom_index] as usize
        + graph_h_neighbor_count_with_adjacency(mol, adjacency, atom_index);
    let tnz_degree = nz_degree + total_num_hs;
    let mut legal_center = true;
    let mut has_dupes = false;
    let mut nbrs = Vec::new();

    if tnz_degree > 4 || tnz_degree < 3 || (nz_degree < 3 && !matches!(atom.atomic_num, 15 | 33)) {
        legal_center = false;
    } else if nz_degree == 3 {
        if total_num_hs == 1 {
            if has_protium_neighbor_with_adjacency(mol, adjacency, atom_index) {
                legal_center = false;
            }
        } else {
            legal_center = false;
            match atom.atomic_num {
                7 => {
                    let in_three_ring = adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
                        cycle_cache.min_cycle_size_for_bond(mol, adjacency, neighbor.bond_index)
                            == Some(3)
                    });
                    if !atom_has_conjugated_bond_with_adjacency(mol, adjacency, atom_index)
                        && (in_three_ring
                            || is_atom_bridgehead_with_cache(
                                mol,
                                adjacency,
                                cycle_cache,
                                atom_index,
                            ))
                    {
                        legal_center = true;
                    }
                }
                15 | 33 => {
                    legal_center = true;
                }
                16 | 34 => {
                    let ev = explicit_valence[atom_index] as i32;
                    if ev == 4 || (ev == 3 && atom.formal_charge == 1) {
                        legal_center = true;
                    }
                }
                _ => {}
            }
        }
    }

    if legal_center && !ranks.is_empty() {
        let mut seen_ranks = HashMap::<i64, ()>::new();
        for neighbor in adjacency.neighbors_of(atom_index) {
            let bond = &mol.bonds()[neighbor.bond_index];
            nbrs.push((ranks[neighbor.atom_index], neighbor.bond_index));
            if !bond_affects_atom_chirality(bond, atom_index) {
                continue;
            }
            if seen_ranks.insert(ranks[neighbor.atom_index], ()).is_some() {
                has_dupes = true;
                break;
            }
        }
    }

    (legal_center, has_dupes, nbrs)
}

impl Molecule {
    /// Return ordered tetrahedral stereocenters derived from RDKit-style atom
    /// chiral tags and molecule bond order.
    #[must_use]
    pub fn tetrahedral_stereo(&self) -> Vec<TetrahedralStereo> {
        let mut out = Vec::new();
        for atom in self.atoms() {
            let tag = atom.chiral_tag;
            if !matches!(tag, ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw) {
                continue;
            }

            let mut ligands = Vec::with_capacity(4);
            for bond in self.bonds() {
                if bond.begin_atom == atom.index {
                    ligands.push(LigandRef::Atom(bond.end_atom));
                } else if bond.end_atom == atom.index {
                    ligands.push(LigandRef::Atom(bond.begin_atom));
                }
            }

            if ligands.len() == 3 {
                ligands.push(LigandRef::ImplicitH);
            }
            let Ok(mut ligands) = <[LigandRef; 4]>::try_from(ligands) else {
                continue;
            };

            // RDKit DistGeom treats CHI_TETRAHEDRAL_CCW as positive volume
            // for atom bond order, and CHI_TETRAHEDRAL_CW as negative. COSMolKit
            // stores ordered ligands with positive orientation, so CW is an odd
            // permutation of the RDKit atom-bond order.
            if matches!(tag, ChiralTag::TetrahedralCw) {
                ligands.swap(0, 1);
            }

            out.push(TetrahedralStereo {
                center: atom.index,
                ligands,
            });
        }
        out
    }

    /// Hidden helper that mirrors the atom-property side of RDKit legacy
    /// stereochemistry perception closely enough for parity testing.
    #[doc(hidden)]
    #[must_use]
    pub fn rdkit_legacy_stereo_atom_props(
        &self,
        flag_possible_stereo_centers: bool,
    ) -> Vec<LegacyStereoAtomProps> {
        let Ok(assignment) = assign_valence(self, ValenceModel::RdkitLike) else {
            return vec![
                LegacyStereoAtomProps {
                    cip_code: None,
                    cip_rank: None,
                    chirality_possible: false,
                };
                self.atoms().len()
            ];
        };
        let adjacency = AdjacencyList::from_topology(self.atoms().len(), self.bonds());
        let mut cycle_cache = StereoBondCycleCache::new(self.bonds().len());

        let mut has_stereo_atoms = false;
        let mut has_potential_stereo_atoms = false;
        for atom in self.atoms() {
            if !has_stereo_atoms && !matches!(atom.chiral_tag, ChiralTag::Unspecified) {
                has_stereo_atoms = true;
            } else if flag_possible_stereo_centers && !has_potential_stereo_atoms {
                has_potential_stereo_atoms = is_atom_potential_chiral_center_with_cache(
                    self,
                    &adjacency,
                    &mut cycle_cache,
                    atom.index,
                    &[],
                    &assignment.explicit_valence,
                    &assignment.implicit_hydrogens,
                )
                .0;
            }
        }

        let mut has_stereo_bonds = false;
        let mut has_potential_stereo_bonds = false;
        for (bond_index, bond) in self.bonds().iter().enumerate() {
            if !matches!(bond.order, BondOrder::Double) {
                continue;
            }
            let begin = bond.begin_atom;
            let end = bond.end_atom;
            let is_specified = self.bonds().iter().any(|other| {
                (other.begin_atom == begin
                    || other.end_atom == begin
                    || other.begin_atom == end
                    || other.end_atom == end)
                    && matches!(
                        other.direction,
                        crate::BondDirection::EndUpRight | crate::BondDirection::EndDownRight
                    )
            });
            if is_specified {
                has_stereo_bonds = true;
            } else if flag_possible_stereo_centers
                && !has_potential_stereo_bonds
                && cycle_cache
                    .min_cycle_size_for_bond(self, &adjacency, bond_index)
                    .is_none_or(|size| size >= 8)
            {
                has_potential_stereo_bonds = true;
            }
            if has_stereo_bonds && has_potential_stereo_bonds {
                break;
            }
        }

        let mut keep_going = has_stereo_atoms || has_stereo_bonds;
        if !keep_going {
            keep_going = flag_possible_stereo_centers
                && (has_potential_stereo_atoms || has_potential_stereo_bonds);
        }
        if !keep_going {
            return vec![
                LegacyStereoAtomProps {
                    cip_code: None,
                    cip_rank: None,
                    chirality_possible: false,
                };
                self.atoms().len()
            ];
        }

        let ranks = crate::io::molblock::rdkit_cip_ranks_for_depict(self);
        let mut out: Vec<LegacyStereoAtomProps> = ranks
            .iter()
            .map(|&rank| LegacyStereoAtomProps {
                cip_code: None,
                cip_rank: Some(rank),
                chirality_possible: false,
            })
            .collect();

        for atom in self.atoms() {
            let tag = atom.chiral_tag;
            if !flag_possible_stereo_centers && matches!(tag, ChiralTag::Unspecified) {
                continue;
            }

            let (legal_center, has_dupes, mut nbrs) = is_atom_potential_chiral_center_with_cache(
                self,
                &adjacency,
                &mut cycle_cache,
                atom.index,
                &ranks,
                &assignment.explicit_valence,
                &assignment.implicit_hydrogens,
            );
            if legal_center && !has_dupes && flag_possible_stereo_centers {
                out[atom.index].chirality_possible = true;
            }
            if !legal_center || has_dupes || matches!(tag, ChiralTag::Unspecified) {
                continue;
            }

            nbrs.sort();
            let probe: Vec<usize> = nbrs.into_iter().map(|(_, bond_index)| bond_index).collect();
            let mut swaps =
                perturbation_order_for_atom_with_adjacency(&adjacency, atom.index, &probe);
            let total_num_hs = atom.explicit_hydrogens as usize
                + assignment.implicit_hydrogens[atom.index] as usize
                + graph_h_neighbor_count_with_adjacency(self, &adjacency, atom.index);
            if probe.len() == 3 && total_num_hs == 1 {
                swaps += 1;
            }
            let mut final_tag = tag;
            if swaps % 2 == 1 {
                final_tag = match final_tag {
                    ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                    ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                    ChiralTag::TrigonalBipyramidal => ChiralTag::TrigonalBipyramidal,
                    ChiralTag::Unspecified => ChiralTag::Unspecified,
                };
            }
            out[atom.index].cip_code = Some(
                if matches!(final_tag, ChiralTag::TetrahedralCcw) {
                    "S"
                } else {
                    "R"
                }
                .to_owned(),
            );
        }
        let cip_codes: Vec<Option<String>> =
            out.iter().map(|props| props.cip_code.clone()).collect();
        if cip_codes.iter().any(Option::is_some)
            && legacy_assign_bond_stereo_would_leave_unassigned_with_cache(
                self,
                &adjacency,
                &mut cycle_cache,
            )
        {
            let reranked =
                crate::io::molblock::rdkit_cip_reranks_with_legacy_stereo(self, &ranks, &cip_codes);
            for (props, rank) in out.iter_mut().zip(reranked) {
                props.cip_rank = Some(rank);
            }
        }
        out
    }

    /// Hidden helper mirroring RDKit's legacy small-ring gate for double-bond
    /// stereochemistry.
    #[doc(hidden)]
    #[must_use]
    pub fn rdkit_should_detect_double_bond_stereo(&self, bond_index: usize) -> bool {
        should_detect_double_bond_stereo(self, bond_index)
    }
}

pub(crate) fn cache_rdkit_legacy_cip_ranks(mol: &mut Molecule) {
    let ranks = mol.rdkit_legacy_stereo_atom_props(true);
    for (atom, props) in mol.atoms_mut().iter_mut().zip(ranks) {
        atom.rdkit_cip_rank = props.cip_rank;
    }
}

#[cfg(test)]
mod tests {
    use super::{LigandRef, TetrahedralStereo};
    use crate::Molecule;

    #[test]
    fn tetrahedral_stereo_is_derived_from_chiral_tags() {
        let ccw = Molecule::from_smiles("F[C@](Cl)(Br)I").expect("parse chiral smiles");
        let cw = Molecule::from_smiles("F[C@@](Cl)(Br)I").expect("parse chiral smiles");

        let ccw_stereo = ccw.tetrahedral_stereo();
        let cw_stereo = cw.tetrahedral_stereo();

        assert_eq!(ccw_stereo.len(), 1);
        assert_eq!(cw_stereo.len(), 1);
        assert_eq!(ccw_stereo[0].center, 1);
        assert_eq!(cw_stereo[0].center, 1);
        assert_eq!(
            ccw_stereo[0].ligands,
            [
                LigandRef::Atom(0),
                LigandRef::Atom(2),
                LigandRef::Atom(3),
                LigandRef::Atom(4)
            ]
        );
        assert_eq!(
            cw_stereo[0].ligands,
            [
                LigandRef::Atom(2),
                LigandRef::Atom(0),
                LigandRef::Atom(3),
                LigandRef::Atom(4)
            ]
        );
    }

    #[test]
    fn tetrahedral_stereo_places_implicit_hydrogen_as_fourth_ligand() {
        let mol = Molecule::from_smiles("[13CH3:7][C@H](F)Cl").expect("parse chiral smiles");
        let stereo = mol.tetrahedral_stereo();

        assert_eq!(
            stereo,
            vec![TetrahedralStereo {
                center: 1,
                ligands: [
                    LigandRef::Atom(0),
                    LigandRef::Atom(2),
                    LigandRef::Atom(3),
                    LigandRef::ImplicitH
                ],
            }]
        );
    }
}
