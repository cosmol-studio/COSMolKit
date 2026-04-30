use crate::smiles_write::SmilesWriteError;
use crate::{BondOrder, Molecule};
use std::cmp::Ordering;

const MAX_NATOMS: i64 = 5000;
const MAX_CYCLES: usize = 1024;
const MAX_BONDTYPE: i64 = 32;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum AtomColor {
    White,
    Grey,
    Black,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum MolStackType {
    Atom,
    Bond,
    Ring,
    BranchOpen,
    BranchClose,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MolStackElem {
    Atom {
        atom_idx: usize,
    },
    Bond {
        bond_idx: usize,
        atom_to_left_idx: usize,
        traversal_ring_closure: bool,
    },
    Ring {
        ring_idx: usize,
    },
    BranchOpen {
        branch_idx: usize,
    },
    BranchClose {
        branch_idx: usize,
    },
}

impl MolStackElem {
    pub fn kind(&self) -> MolStackType {
        match self {
            Self::Atom { .. } => MolStackType::Atom,
            Self::Bond { .. } => MolStackType::Bond,
            Self::Ring { .. } => MolStackType::Ring,
            Self::BranchOpen { .. } => MolStackType::BranchOpen,
            Self::BranchClose { .. } => MolStackType::BranchClose,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FragmentTraversal {
    pub atom_colors: Vec<AtomColor>,
    pub atom_ranks: Vec<u32>,
    pub mol_stack: Vec<MolStackElem>,
    pub atom_ring_closure_counts: Vec<usize>,
    pub atom_traversal_bond_order: Vec<Vec<usize>>,
    pub chiral_tag_overrides: Vec<Option<crate::ChiralTag>>,
    pub bond_direction_overrides: Vec<Option<crate::BondDirection>>,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
struct Possible {
    rank: i64,
    other_idx: usize,
    bond_idx: usize,
}

fn possible_compare(lhs: &Possible, rhs: &Possible) -> std::cmp::Ordering {
    lhs.rank.cmp(&rhs.rank)
}

fn rdkit_bond_type_value(order: BondOrder) -> i64 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 12,
        BondOrder::Dative => 17,
        BondOrder::Null => 0,
    }
}

fn adjacency(mol: &Molecule) -> crate::AdjacencyList {
    mol.adjacency
        .clone()
        .unwrap_or_else(|| crate::AdjacencyList::from_topology(mol.atoms.len(), &mol.bonds))
}

fn bond_in_any_cycle(mol: &Molecule, bond_idx: usize) -> bool {
    let bond = &mol.bonds[bond_idx];
    let adj = adjacency(mol);
    let mut seen = vec![false; mol.atoms.len()];
    let mut queue = std::collections::VecDeque::new();
    seen[bond.begin_atom] = true;
    queue.push_back(bond.begin_atom);
    while let Some(atom_idx) = queue.pop_front() {
        for nbr in adj.neighbors_of(atom_idx) {
            if nbr.bond_index == bond_idx {
                continue;
            }
            if nbr.atom_index == bond.end_atom {
                return true;
            }
            if !seen[nbr.atom_index] {
                seen[nbr.atom_index] = true;
                queue.push_back(nbr.atom_index);
            }
        }
    }
    false
}

fn build_possibles(
    mol: &Molecule,
    atom_idx: usize,
    in_bond_idx: Option<usize>,
    colors: &[AtomColor],
    ranks: &[u32],
    do_random: bool,
) -> Vec<Possible> {
    let adj = adjacency(mol);
    let mut possibles = Vec::new();
    for nbr in adj.neighbors_of(atom_idx) {
        if in_bond_idx.is_some_and(|idx| idx == nbr.bond_index) {
            continue;
        }
        let mut rank = i64::from(ranks[nbr.atom_index]);
        if !do_random {
            if colors[nbr.atom_index] == AtomColor::Grey {
                rank -= (MAX_BONDTYPE + 1) * MAX_NATOMS * MAX_NATOMS;
                rank += (MAX_BONDTYPE - rdkit_bond_type_value(mol.bonds[nbr.bond_index].order))
                    * MAX_NATOMS;
            } else if bond_in_any_cycle(mol, nbr.bond_index) {
                rank += (MAX_BONDTYPE - rdkit_bond_type_value(mol.bonds[nbr.bond_index].order))
                    * MAX_NATOMS
                    * MAX_NATOMS;
            }
        } else {
            return Vec::new();
        }
        possibles.push(Possible {
            rank,
            other_idx: nbr.atom_index,
            bond_idx: nbr.bond_index,
        });
    }
    possibles.sort_by(possible_compare);
    possibles
}

fn dfs_find_cycles(
    mol: &Molecule,
    atom_idx: usize,
    in_bond_idx: Option<usize>,
    colors: &mut [AtomColor],
    ranks: &[u32],
    atom_ring_closures: &mut [Vec<usize>],
    do_random: bool,
) {
    colors[atom_idx] = AtomColor::Grey;
    for possible in build_possibles(mol, atom_idx, in_bond_idx, colors, ranks, do_random) {
        match colors[possible.other_idx] {
            AtomColor::White => {
                dfs_find_cycles(
                    mol,
                    possible.other_idx,
                    Some(possible.bond_idx),
                    colors,
                    ranks,
                    atom_ring_closures,
                    do_random,
                );
            }
            AtomColor::Grey => {
                atom_ring_closures[possible.other_idx].push(possible.bond_idx);
                atom_ring_closures[atom_idx].push(possible.bond_idx);
            }
            AtomColor::Black => {}
        }
    }
    colors[atom_idx] = AtomColor::Black;
}

fn dfs_build_stack(
    mol: &Molecule,
    atom_idx: usize,
    in_bond_idx: Option<usize>,
    colors: &mut [AtomColor],
    ranks: &[u32],
    cycles_available: &mut [bool],
    ring_openings: &mut std::collections::HashMap<usize, usize>,
    mol_stack: &mut Vec<MolStackElem>,
    atom_ring_closures: &mut [Vec<usize>],
    atom_traversal_bond_order: &mut [Vec<usize>],
    do_random: bool,
) -> Result<(), SmilesWriteError> {
    let adj = adjacency(mol);
    let mut seen_from_here = vec![false; mol.atoms.len()];
    seen_from_here[atom_idx] = true;
    mol_stack.push(MolStackElem::Atom { atom_idx });
    colors[atom_idx] = AtomColor::Grey;
    let mut trav_list = Vec::<usize>::new();
    if let Some(in_bond_idx) = in_bond_idx {
        trav_list.push(in_bond_idx);
    }

    if !atom_ring_closures[atom_idx].is_empty() {
        let mut rings_closed = Vec::new();
        let closures = atom_ring_closures[atom_idx].clone();
        for bond_idx in closures {
            trav_list.push(bond_idx);
            let bond = &mol.bonds[bond_idx];
            let other = if bond.begin_atom == atom_idx {
                bond.end_atom
            } else {
                bond.begin_atom
            };
            seen_from_here[other] = true;
            if let Some(existing_ring_idx) = ring_openings.remove(&bond_idx) {
                mol_stack.push(MolStackElem::Bond {
                    bond_idx,
                    atom_to_left_idx: atom_idx,
                    traversal_ring_closure: true,
                });
                mol_stack.push(MolStackElem::Ring {
                    ring_idx: existing_ring_idx,
                });
                rings_closed.push(existing_ring_idx - 1);
            } else if let Some(lowest) = cycles_available.iter().position(|v| *v) {
                cycles_available[lowest] = false;
                ring_openings.insert(bond_idx, lowest + 1);
                mol_stack.push(MolStackElem::Ring {
                    ring_idx: lowest + 1,
                });
            } else {
                return Err(SmilesWriteError::UnsupportedPath(
                    "Too many rings open at once. SMILES cannot be generated.",
                ));
            }
        }
        for ring_idx in rings_closed {
            cycles_available[ring_idx] = true;
        }
    }

    let mut possibles = Vec::new();
    for nbr in adj.neighbors_of(atom_idx) {
        if in_bond_idx.is_some_and(|idx| idx == nbr.bond_index) {
            continue;
        }
        if colors[nbr.atom_index] != AtomColor::White || seen_from_here[nbr.atom_index] {
            continue;
        }
        let mut rank = i64::from(ranks[nbr.atom_index]);
        if !do_random {
            if bond_in_any_cycle(mol, nbr.bond_index) {
                rank += (MAX_BONDTYPE - rdkit_bond_type_value(mol.bonds[nbr.bond_index].order))
                    * MAX_NATOMS
                    * MAX_NATOMS;
            }
        } else {
            return Err(SmilesWriteError::UnsupportedPath(
                "RDKit doRandom traversal path is not ported yet",
            ));
        }
        possibles.push(Possible {
            rank,
            other_idx: nbr.atom_index,
            bond_idx: nbr.bond_index,
        });
    }
    possibles.sort_by(possible_compare);

    for (branch_idx, possible) in possibles.iter().enumerate() {
        if colors[possible.other_idx] != AtomColor::White {
            continue;
        }
        trav_list.push(possible.bond_idx);
        if branch_idx + 1 != possibles.len() {
            mol_stack.push(MolStackElem::BranchOpen { branch_idx });
        }
        mol_stack.push(MolStackElem::Bond {
            bond_idx: possible.bond_idx,
            atom_to_left_idx: atom_idx,
            traversal_ring_closure: false,
        });
        dfs_build_stack(
            mol,
            possible.other_idx,
            Some(possible.bond_idx),
            colors,
            ranks,
            cycles_available,
            ring_openings,
            mol_stack,
            atom_ring_closures,
            atom_traversal_bond_order,
            do_random,
        )?;
        if branch_idx + 1 != possibles.len() {
            mol_stack.push(MolStackElem::BranchClose { branch_idx });
        }
    }

    colors[atom_idx] = AtomColor::Black;
    atom_traversal_bond_order[atom_idx] = trav_list;
    Ok(())
}

#[derive(Debug, Clone)]
struct CanonBondHolder {
    bond_type: u8,
    bond_stereo: u8,
    nbr_sym_class: usize,
    nbr_idx: usize,
    stype: crate::BondStereo,
    controlling_atoms: [Option<usize>; 4],
}

#[derive(Debug, Clone)]
struct CanonRankAtom {
    index: usize,
    degree: usize,
    total_num_hs: usize,
    atomic_num: u8,
    isotope: u16,
    atom_map_num: u32,
    formal_charge: i8,
    chiral_presence: bool,
    chiral_tag: crate::ChiralTag,
    has_ring_nbr: bool,
    is_ring_stereo_atom: bool,
    nbr_ids: Vec<usize>,
    bonds: Vec<CanonBondHolder>,
}

fn rdkit_bond_type_rank(order: BondOrder) -> u8 {
    match order {
        BondOrder::Null => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 12,
        BondOrder::Dative => 17,
    }
}

fn rdkit_bond_stereo_rank(stereo: crate::BondStereo) -> u8 {
    match stereo {
        crate::BondStereo::None => 0,
        crate::BondStereo::Any => 1,
        crate::BondStereo::Cis => 4,
        crate::BondStereo::Trans => 5,
    }
}

fn flip_if_needed(
    stereo: crate::BondStereo,
    controlling_atoms: &[Option<usize>; 4],
    atoms: &[CanonRankAtom],
) -> crate::BondStereo {
    let mut flip = false;
    let first0 = controlling_atoms[0].expect("missing controlling atom 0");
    let first2 = controlling_atoms[2].expect("missing controlling atom 2");
    if let Some(ctrl1) = controlling_atoms[1]
        && atoms[ctrl1].index > atoms[first0].index
    {
        flip = !flip;
    }
    if let Some(ctrl3) = controlling_atoms[3]
        && atoms[ctrl3].index > atoms[first2].index
    {
        flip = !flip;
    }
    if flip {
        match stereo {
            crate::BondStereo::Cis => crate::BondStereo::Trans,
            crate::BondStereo::Trans => crate::BondStereo::Cis,
            other => other,
        }
    } else {
        stereo
    }
}

fn init_canon_rank_atoms(
    mol: &Molecule,
    include_chirality: bool,
) -> Result<Vec<CanonRankAtom>, SmilesWriteError> {
    let assignment = crate::assign_valence(mol, crate::ValenceModel::RdkitLike).map_err(|_| {
        SmilesWriteError::UnsupportedPath(
            "RDKit canonical ranking requires RDKit-like valence assignment",
        )
    })?;
    let ring_stereo_atoms = if include_chirality {
        find_chiral_atom_special_cases(mol)?
    } else {
        vec![Vec::new(); mol.atoms.len()]
    };
    let adjacency = adjacency(mol);
    let mut atoms = Vec::with_capacity(mol.atoms.len());
    for atom in &mol.atoms {
        let has_ring_nbr = adjacency.neighbors_of(atom.index).iter().any(|nbr| {
            !ring_stereo_atoms[nbr.atom_index].is_empty()
                && matches!(
                    mol.atoms[nbr.atom_index].chiral_tag,
                    crate::ChiralTag::TetrahedralCw | crate::ChiralTag::TetrahedralCcw
                )
        });
        atoms.push(CanonRankAtom {
            index: atom.index,
            degree: adjacency.neighbors_of(atom.index).len(),
            total_num_hs: atom.explicit_hydrogens as usize
                + assignment.implicit_hydrogens[atom.index] as usize,
            atomic_num: atom.atomic_num,
            isotope: atom.isotope.unwrap_or(0),
            atom_map_num: atom.atom_map_num.unwrap_or(0),
            formal_charge: atom.formal_charge,
            chiral_presence: !matches!(atom.chiral_tag, crate::ChiralTag::Unspecified),
            chiral_tag: atom.chiral_tag,
            has_ring_nbr,
            is_ring_stereo_atom: !ring_stereo_atoms[atom.index].is_empty()
                && matches!(
                    atom.chiral_tag,
                    crate::ChiralTag::TetrahedralCw | crate::ChiralTag::TetrahedralCcw
                ),
            nbr_ids: adjacency
                .neighbors_of(atom.index)
                .iter()
                .map(|nbr| nbr.atom_index)
                .collect(),
            bonds: Vec::new(),
        });
    }
    for atom_idx in 0..mol.atoms.len() {
        let mut bonds = Vec::new();
        for nbr in adjacency.neighbors_of(atom_idx) {
            let bond = &mol.bonds[nbr.bond_index];
            bonds.push(CanonBondHolder {
                bond_type: rdkit_bond_type_rank(bond.order),
                bond_stereo: if include_chirality {
                    rdkit_bond_stereo_rank(bond.stereo)
                } else {
                    0
                },
                nbr_sym_class: 0,
                nbr_idx: nbr.atom_index,
                stype: if include_chirality {
                    bond.stereo
                } else {
                    crate::BondStereo::None
                },
                controlling_atoms: canonical_bond_controlling_atoms(mol, bond),
            });
        }
        bonds.sort_by(canon_bond_greater);
        atoms[atom_idx].bonds = bonds;
    }
    Ok(atoms)
}

fn canonical_bond_controlling_atoms(mol: &Molecule, bond: &crate::Bond) -> [Option<usize>; 4] {
    let mut out = [None, None, None, None];
    if !matches!(
        bond.stereo,
        crate::BondStereo::Cis | crate::BondStereo::Trans
    ) {
        return out;
    }
    if bond.stereo_atoms.len() >= 2 {
        out[0] = Some(bond.stereo_atoms[0]);
        out[2] = Some(bond.stereo_atoms[1]);
    } else {
        return out;
    }
    let begin_degree = adjacency(mol).neighbors_of(bond.begin_atom).len();
    if begin_degree > 2 {
        for nbr in adjacency(mol).neighbors_of(bond.begin_atom) {
            if nbr.atom_index != bond.end_atom && Some(nbr.atom_index) != out[0] {
                out[1] = Some(nbr.atom_index);
            }
        }
    }
    let end_degree = adjacency(mol).neighbors_of(bond.end_atom).len();
    if end_degree > 2 {
        for nbr in adjacency(mol).neighbors_of(bond.end_atom) {
            if nbr.atom_index != bond.begin_atom && Some(nbr.atom_index) != out[2] {
                out[3] = Some(nbr.atom_index);
            }
        }
    }
    out
}

fn canon_bond_compare(lhs: &CanonBondHolder, rhs: &CanonBondHolder) -> Ordering {
    lhs.bond_type
        .cmp(&rhs.bond_type)
        .then_with(|| lhs.bond_stereo.cmp(&rhs.bond_stereo))
        .then_with(|| lhs.nbr_sym_class.cmp(&rhs.nbr_sym_class))
}

fn canon_bond_greater(lhs: &CanonBondHolder, rhs: &CanonBondHolder) -> Ordering {
    canon_bond_compare(rhs, lhs)
}

fn canon_bond_compare_with_atoms(
    atoms: &[CanonRankAtom],
    lhs: &CanonBondHolder,
    rhs: &CanonBondHolder,
) -> Ordering {
    let mut ord = lhs.bond_type.cmp(&rhs.bond_type);
    if ord != Ordering::Equal {
        return ord;
    }
    ord = lhs.bond_stereo.cmp(&rhs.bond_stereo);
    if ord != Ordering::Equal {
        return ord;
    }
    ord = lhs.nbr_sym_class.cmp(&rhs.nbr_sym_class);
    if ord != Ordering::Equal {
        return ord;
    }
    if lhs.bond_stereo != 0 && rhs.bond_stereo != 0 {
        let st1 = lhs.stype;
        let st2 = rhs.stype;
        let cmp = compare_bond_stereo(
            atoms,
            st1,
            &lhs.controlling_atoms,
            st2,
            &rhs.controlling_atoms,
        );
        if cmp != 0 {
            return if cmp < 0 {
                Ordering::Less
            } else {
                Ordering::Greater
            };
        }
    }
    Ordering::Equal
}

fn compare_bond_stereo(
    atoms: &[CanonRankAtom],
    st1: crate::BondStereo,
    c1: &[Option<usize>; 4],
    st2: crate::BondStereo,
    c2: &[Option<usize>; 4],
) -> i32 {
    use crate::BondStereo::{Any, Cis, None, Trans};
    match (st1, st2) {
        (None, None) => return 0,
        (None, _) => return -1,
        (_, None) => return 1,
        (Any, Any) => return 0,
        (Any, _) => return -1,
        (_, Any) => return 1,
        _ => {}
    }
    let fst1 = flip_if_needed(st1, c1, atoms);
    let fst2 = flip_if_needed(st2, c2, atoms);
    let rank = |s| match s {
        Cis => 4u8,
        Trans => 5u8,
        None => 0u8,
        Any => 1u8,
    };
    rank(fst1).cmp(&rank(fst2)) as i32
}

fn update_atom_neighbor_index(atoms: &mut [CanonRankAtom], atom_idx: usize) {
    let classes: Vec<(usize, usize)> = atoms[atom_idx]
        .bonds
        .iter()
        .map(|bond| (bond.nbr_idx, atoms[bond.nbr_idx].index))
        .collect();
    for (bond, (_, cls)) in atoms[atom_idx].bonds.iter_mut().zip(classes.into_iter()) {
        bond.nbr_sym_class = cls;
    }
    atoms[atom_idx]
        .bonds
        .sort_by(|lhs, rhs| canon_bond_greater(lhs, rhs));
}

#[derive(Debug, Clone, Copy)]
struct AtomCompareFunctor {
    use_nbrs: bool,
    use_isotopes: bool,
    use_chirality: bool,
    use_atom_maps: bool,
    use_non_stereo_ranks: bool,
    use_chiral_presence: bool,
}

impl AtomCompareFunctor {
    fn get_atom_ring_nbr_code(&self, atoms: &[CanonRankAtom], i: usize) -> u32 {
        if !atoms[i].has_ring_nbr {
            return 0;
        }
        let mut code = 0u32;
        for &nbr in &atoms[i].nbr_ids {
            if atoms[nbr].is_ring_stereo_atom {
                code = code.saturating_add((atoms[nbr].index as u32) * 10000 + 1);
            }
        }
        code
    }

    fn get_chiral_rank(&self, atoms: &[CanonRankAtom], i: usize) -> u32 {
        let mut perm = Vec::with_capacity(atoms[i].degree);
        for &nbr in &atoms[i].nbr_ids {
            let rank = atoms[nbr].index;
            if perm.contains(&rank) {
                return 0;
            }
            perm.push(rank);
        }
        if perm.len() != atoms[i].degree {
            return 0;
        }
        match atoms[i].chiral_tag {
            crate::ChiralTag::TetrahedralCw | crate::ChiralTag::TetrahedralCcw => {
                let mut sorted = perm.clone();
                sorted.sort_unstable();
                let swaps = count_swaps(&perm, &sorted);
                let mut res = if matches!(atoms[i].chiral_tag, crate::ChiralTag::TetrahedralCw) {
                    2
                } else {
                    1
                };
                if swaps % 2 == 1 {
                    res = if res == 2 { 1 } else { 2 };
                }
                res
            }
            crate::ChiralTag::Unspecified => 0,
        }
    }

    fn base_compare(&self, atoms: &[CanonRankAtom], i: usize, j: usize) -> i32 {
        let ai = &atoms[i];
        let aj = &atoms[j];
        macro_rules! cmp_field {
            ($lhs:expr, $rhs:expr) => {
                if $lhs < $rhs {
                    return -1;
                } else if $lhs > $rhs {
                    return 1;
                }
            };
        }
        cmp_field!(ai.index, aj.index);
        if self.use_non_stereo_ranks {
            return 0;
        }
        if self.use_atom_maps || ai.atomic_num == 0 || aj.atomic_num == 0 {
            cmp_field!(ai.atom_map_num, aj.atom_map_num);
        }
        cmp_field!(ai.degree, aj.degree);
        cmp_field!(ai.atomic_num, aj.atomic_num);
        if self.use_isotopes {
            cmp_field!(ai.isotope, aj.isotope);
        }
        cmp_field!(ai.total_num_hs, aj.total_num_hs);
        cmp_field!(ai.formal_charge as u32, aj.formal_charge as u32);
        if self.use_chiral_presence {
            cmp_field!(ai.chiral_presence, aj.chiral_presence);
        }
        if self.use_chirality {
            let i_has = !matches!(ai.chiral_tag, crate::ChiralTag::Unspecified);
            let j_has = !matches!(aj.chiral_tag, crate::ChiralTag::Unspecified);
            cmp_field!(i_has, j_has);
            if i_has && j_has {
                cmp_field!(
                    self.get_chiral_rank(atoms, i),
                    self.get_chiral_rank(atoms, j)
                );
            }
            cmp_field!(
                self.get_atom_ring_nbr_code(atoms, i),
                self.get_atom_ring_nbr_code(atoms, j)
            );
        }
        0
    }

    fn compare(&self, atoms: &mut [CanonRankAtom], i: usize, j: usize) -> i32 {
        let v = self.base_compare(atoms, i, j);
        if v != 0 {
            return v;
        }
        if self.use_nbrs {
            update_atom_neighbor_index(atoms, i);
            update_atom_neighbor_index(atoms, j);
            let nb = atoms[i].bonds.len().min(atoms[j].bonds.len());
            for idx in 0..nb {
                match canon_bond_compare_with_atoms(
                    atoms,
                    &atoms[i].bonds[idx],
                    &atoms[j].bonds[idx],
                ) {
                    Ordering::Less => return -1,
                    Ordering::Greater => return 1,
                    Ordering::Equal => {}
                }
            }
            if atoms[i].bonds.len() < atoms[j].bonds.len() {
                return -1;
            } else if atoms[i].bonds.len() > atoms[j].bonds.len() {
                return 1;
            }
        }
        0
    }
}

fn count_swaps(lhs: &[usize], rhs: &[usize]) -> usize {
    let mut probe = lhs.to_vec();
    let mut positions: std::collections::BTreeMap<usize, usize> =
        probe.iter().enumerate().map(|(idx, &v)| (v, idx)).collect();
    let mut swaps = 0usize;
    for (idx, &target) in rhs.iter().enumerate() {
        if probe[idx] == target {
            continue;
        }
        let j = *positions
            .get(&target)
            .expect("target missing in count_swaps permutation");
        let current = probe[idx];
        probe.swap(idx, j);
        positions.insert(current, j);
        positions.insert(target, idx);
        swaps += 1;
    }
    swaps
}

fn hanoi_sort_partition(
    order: &mut [usize],
    count: &mut [usize],
    changed: &mut [bool],
    atoms: &mut [CanonRankAtom],
    cmp: AtomCompareFunctor,
) {
    fn hanoi(
        base: &mut [usize],
        temp: &mut [usize],
        count: &mut [usize],
        changed: &mut [bool],
        atoms: &mut [CanonRankAtom],
        cmp: AtomCompareFunctor,
    ) -> bool {
        let nel = base.len();
        if nel == 1 {
            count[base[0]] = 1;
            return false;
        } else if nel == 2 {
            let n1 = base[0];
            let n2 = base[1];
            let stat = if changed[n1] || changed[n2] {
                cmp.compare(atoms, n1, n2)
            } else {
                0
            };
            if stat == 0 {
                count[n1] = 2;
                count[n2] = 0;
                return false;
            } else if stat < 0 {
                count[n1] = 1;
                count[n2] = 1;
                return false;
            } else {
                count[n1] = 1;
                count[n2] = 1;
                base[0] = n2;
                base[1] = n1;
                return false;
            }
        }
        let n1 = nel / 2;
        let n2 = nel - n1;
        let (b1, b2) = base.split_at_mut(n1);
        let (t1, t2) = temp.split_at_mut(n1);
        let (s1, result) = if hanoi(b1, t1, count, changed, atoms, cmp) {
            (t1.to_vec(), false)
        } else {
            (b1.to_vec(), true)
        };
        let s2 = if hanoi(b2, t2, count, changed, atoms, cmp) {
            t2.to_vec()
        } else {
            b2.to_vec()
        };
        let mut ptr = Vec::with_capacity(nel);
        let mut i1 = 0usize;
        let mut i2 = 0usize;
        let mut rem1 = n1;
        let mut rem2 = n2;
        loop {
            let stat = if changed[s1[i1]] || changed[s2[i2]] {
                cmp.compare(atoms, s1[i1], s2[i2])
            } else {
                0
            };
            let len1 = count[s1[i1]];
            let len2 = count[s2[i2]];
            if stat == 0 {
                count[s1[i1]] = len1 + len2;
                count[s2[i2]] = 0;
                ptr.extend_from_slice(&s1[i1..i1 + len1]);
                i1 += len1;
                rem1 -= len1;
                if rem1 == 0 {
                    ptr.extend_from_slice(&s2[i2..]);
                    break;
                }
                ptr.extend_from_slice(&s2[i2..i2 + len2]);
                i2 += len2;
                rem2 -= len2;
                if rem2 == 0 {
                    ptr.extend_from_slice(&s1[i1..]);
                    break;
                }
            } else if stat < 0 {
                ptr.extend_from_slice(&s1[i1..i1 + len1]);
                i1 += len1;
                rem1 -= len1;
                if rem1 == 0 {
                    ptr.extend_from_slice(&s2[i2..]);
                    break;
                }
            } else {
                ptr.extend_from_slice(&s2[i2..i2 + len2]);
                i2 += len2;
                rem2 -= len2;
                if rem2 == 0 {
                    ptr.extend_from_slice(&s1[i1..]);
                    break;
                }
            }
        }
        if result {
            temp.copy_from_slice(&ptr);
        } else {
            base.copy_from_slice(&ptr);
        }
        result
    }
    let mut temp = vec![0usize; order.len()];
    if hanoi(order, &mut temp, count, changed, atoms, cmp) {
        order.copy_from_slice(&temp);
    }
}

fn activate_partitions(
    n_atoms: usize,
    order: &[usize],
    count: &[usize],
    next: &mut [isize],
    changed: &mut [bool],
) -> isize {
    let _ = n_atoms;
    let mut activeset = -1isize;
    next.fill(-2);
    let mut i = 0usize;
    while i < order.len() {
        let j = order[i];
        if count[j] > 1 {
            next[j] = activeset;
            activeset = j as isize;
            i += count[j];
        } else {
            i += 1;
        }
    }
    for &j in order {
        changed[j] = true;
    }
    activeset
}

fn refine_partitions(
    atoms: &mut [CanonRankAtom],
    order: &mut [usize],
    count: &mut [usize],
    activeset: &mut isize,
    next: &mut [isize],
    changed: &mut [bool],
    touched: &mut [bool],
    cmp: AtomCompareFunctor,
) {
    while *activeset != -1 {
        let partition = *activeset as usize;
        *activeset = next[partition];
        next[partition] = -2;
        let len = count[partition];
        let offset = atoms[partition].index;
        hanoi_sort_partition(&mut order[offset..offset + len], count, changed, atoms, cmp);
        let start = &order[offset..offset + len];
        for &idx in start {
            changed[idx] = false;
        }
        let mut index = start[0];
        let mut i = count[index];
        let mut symclass = 0usize;
        while i < len {
            index = start[i];
            if count[index] != 0 {
                symclass = offset + i;
            }
            atoms[index].index = symclass;
            for &nbr in &atoms[index].nbr_ids {
                changed[nbr] = true;
            }
            i += 1;
        }
        index = start[0];
        i = count[index];
        while i < len {
            index = start[i];
            for &nbr in &atoms[index].nbr_ids {
                touched[atoms[nbr].index] = true;
            }
            i += 1;
        }
        for ii in 0..atoms.len() {
            if touched[ii] {
                let part = order[ii];
                if count[part] > 1 && next[part] == -2 {
                    next[part] = *activeset;
                    *activeset = part as isize;
                }
                touched[ii] = false;
            }
        }
    }
}

fn break_ties_partitions(
    atoms: &mut [CanonRankAtom],
    order: &mut [usize],
    count: &mut [usize],
    activeset: &mut isize,
    next: &mut [isize],
    changed: &mut [bool],
    touched: &mut [bool],
    cmp: AtomCompareFunctor,
) {
    let n_atoms = atoms.len();
    let mut i = 0usize;
    while i < n_atoms {
        let partition = order[i];
        let old_part = atoms[partition].index;
        while count[partition] > 1 {
            let len = count[partition];
            let offset = atoms[partition].index + len - 1;
            let index = order[offset];
            atoms[index].index = offset;
            count[partition] = len - 1;
            count[index] = 1;
            if atoms[index].degree < 1 {
                continue;
            }
            for &nbr in &atoms[index].nbr_ids {
                touched[atoms[nbr].index] = true;
                changed[nbr] = true;
            }
            for ii in 0..n_atoms {
                if touched[ii] {
                    let npart = order[ii];
                    if count[npart] > 1 && next[npart] == -2 {
                        next[npart] = *activeset;
                        *activeset = npart as isize;
                    }
                    touched[ii] = false;
                }
            }
            refine_partitions(atoms, order, count, activeset, next, changed, touched, cmp);
        }
        if atoms[partition].index == old_part {
            i += 1;
        }
    }
}

pub fn rank_mol_atoms(
    mol: &Molecule,
    break_ties: bool,
    include_chirality: bool,
    include_isotopes: bool,
    include_atom_maps: bool,
    include_chiral_presence: bool,
    _include_stereo_groups: bool,
    use_non_stereo_ranks: bool,
    _include_ring_stereo: bool,
) -> Result<Vec<u32>, SmilesWriteError> {
    let _has_bond_stereo = mol.bonds.iter().any(|bond| {
        !matches!(bond.stereo, crate::BondStereo::None) || !bond.stereo_atoms.is_empty()
    });
    let _has_ring_bond = (0..mol.bonds.len()).any(|bond_idx| bond_in_any_cycle(mol, bond_idx));
    let mut atoms = init_canon_rank_atoms(mol, include_chirality)?;
    let n = atoms.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    for atom in &mut atoms {
        atom.index = 0;
    }
    let cmp = AtomCompareFunctor {
        use_nbrs: true,
        use_isotopes: include_isotopes,
        use_chirality: include_chirality,
        use_atom_maps: include_atom_maps,
        use_non_stereo_ranks,
        use_chiral_presence: include_chiral_presence,
    };
    let mut order: Vec<usize> = (0..n).collect();
    let mut count = vec![0usize; n];
    let mut next = vec![-2isize; n];
    let mut changed = vec![true; n];
    let mut touched = vec![false; n];
    count[0] = n;
    let mut activeset = activate_partitions(n, &order, &count, &mut next, &mut changed);
    refine_partitions(
        &mut atoms,
        &mut order,
        &mut count,
        &mut activeset,
        &mut next,
        &mut changed,
        &mut touched,
        cmp,
    );
    if break_ties {
        break_ties_partitions(
            &mut atoms,
            &mut order,
            &mut count,
            &mut activeset,
            &mut next,
            &mut changed,
            &mut touched,
            cmp,
        );
    }
    let mut res = vec![0u32; n];
    for &idx in &order {
        res[idx] = atoms[idx].index as u32;
    }
    Ok(res)
}

pub fn build_noncanonical_fragment(
    mol: &Molecule,
    atom_idx: usize,
    atom_ranks: &[u32],
    do_random: bool,
) -> Result<FragmentTraversal, SmilesWriteError> {
    let mut cycle_colors = vec![AtomColor::White; mol.atoms.len()];
    let mut atom_ring_closures = vec![Vec::<usize>::new(); mol.atoms.len()];
    dfs_find_cycles(
        mol,
        atom_idx,
        None,
        &mut cycle_colors,
        atom_ranks,
        &mut atom_ring_closures,
        do_random,
    );
    let mut colors = vec![AtomColor::White; mol.atoms.len()];
    let mut cycles_available = vec![true; MAX_CYCLES];
    let mut ring_openings = std::collections::HashMap::new();
    let mut atom_traversal_bond_order = vec![Vec::<usize>::new(); mol.atoms.len()];
    let mut mol_stack = Vec::with_capacity(mol.atoms.len() + mol.bonds.len());
    dfs_build_stack(
        mol,
        atom_idx,
        None,
        &mut colors,
        atom_ranks,
        &mut cycles_available,
        &mut ring_openings,
        &mut mol_stack,
        &mut atom_ring_closures,
        &mut atom_traversal_bond_order,
        do_random,
    )?;
    Ok(FragmentTraversal {
        atom_colors: colors,
        atom_ranks: atom_ranks.to_vec(),
        mol_stack,
        atom_ring_closure_counts: atom_ring_closures.iter().map(Vec::len).collect(),
        atom_traversal_bond_order,
        chiral_tag_overrides: vec![None; mol.atoms.len()],
        bond_direction_overrides: vec![None; mol.bonds.len()],
    })
}

pub fn canonicalize_fragment(
    mol: &Molecule,
    atom_idx: usize,
    atom_ranks: &[u32],
    do_isomeric_smiles: bool,
    do_random: bool,
    do_chiral_inversions: bool,
) -> Result<FragmentTraversal, SmilesWriteError> {
    let mut traversal = build_noncanonical_fragment(mol, atom_idx, atom_ranks, do_random)?;
    if !do_isomeric_smiles {
        return Ok(traversal);
    }
    let ring_stereo_atoms = find_chiral_atom_special_cases(mol)?;
    let assignment = crate::assign_valence(mol, crate::ValenceModel::RdkitLike).map_err(|_| {
        SmilesWriteError::UnsupportedPath(
            "canonical fragment chirality handling requires RDKit-like valence assignment",
        )
    })?;
    let first_idx = traversal
        .mol_stack
        .iter()
        .find_map(|elem| match elem {
            MolStackElem::Atom { atom_idx } => Some(*atom_idx),
            _ => None,
        })
        .ok_or(SmilesWriteError::UnsupportedPath(
            "empty MolStack in canonical fragment traversal",
        ))?;
    for atom in &mol.atoms {
        if matches!(atom.chiral_tag, crate::ChiralTag::Unspecified) {
            continue;
        }
        let mut true_order = traversal.atom_traversal_bond_order[atom.index].clone();
        let degree = mol
            .bonds
            .iter()
            .filter(|bond| bond.begin_atom == atom.index || bond.end_atom == atom.index)
            .count();
        if true_order.len() < degree {
            for bond in &mol.bonds {
                if (bond.begin_atom == atom.index || bond.end_atom == atom.index)
                    && !true_order.contains(&bond.index)
                {
                    true_order.push(bond.index);
                }
            }
        }
        let reference: Vec<usize> = mol
            .bonds
            .iter()
            .filter(|bond| bond.begin_atom == atom.index || bond.end_atom == atom.index)
            .map(|bond| bond.index)
            .collect();
        let mut n_swaps = count_swaps(&true_order, &reference);
        if do_chiral_inversions
            && chiral_atom_needs_tag_inversion(
                mol,
                atom.index,
                atom.index == first_idx,
                traversal.atom_ring_closure_counts[atom.index],
                &assignment.implicit_hydrogens,
            )
        {
            n_swaps += 1;
        }
        if n_swaps % 2 == 1 {
            traversal.chiral_tag_overrides[atom.index] = Some(match atom.chiral_tag {
                crate::ChiralTag::TetrahedralCw => crate::ChiralTag::TetrahedralCcw,
                crate::ChiralTag::TetrahedralCcw => crate::ChiralTag::TetrahedralCw,
                crate::ChiralTag::Unspecified => crate::ChiralTag::Unspecified,
            });
        }
    }
    let mut atom_visit_orders = vec![usize::MAX; mol.atoms.len()];
    let mut bond_visit_orders = vec![usize::MAX; mol.bonds.len()];
    let mut traversal_ring_closure_bonds = vec![false; mol.bonds.len()];
    let mut pos = 0usize;
    for elem in &traversal.mol_stack {
        match elem {
            MolStackElem::Atom { atom_idx } => {
                atom_visit_orders[*atom_idx] = pos;
            }
            MolStackElem::Bond {
                bond_idx,
                traversal_ring_closure,
                ..
            } => {
                bond_visit_orders[*bond_idx] = pos;
                traversal_ring_closure_bonds[*bond_idx] = *traversal_ring_closure;
            }
            _ => {}
        }
        pos += 1;
    }
    let mut ring_stereo_adjusted = vec![false; mol.atoms.len()];
    for elem in &traversal.mol_stack {
        let MolStackElem::Atom { atom_idx } = elem else {
            continue;
        };
        let atom = &mol.atoms[*atom_idx];
        if matches!(atom.chiral_tag, crate::ChiralTag::Unspecified)
            || ring_stereo_atoms[*atom_idx].is_empty()
        {
            continue;
        }
        if !ring_stereo_adjusted[*atom_idx] {
            traversal.chiral_tag_overrides[*atom_idx] = Some(crate::ChiralTag::TetrahedralCcw);
            ring_stereo_adjusted[*atom_idx] = true;
        }
        for &nbr_v in &ring_stereo_atoms[*atom_idx] {
            let nbr_idx = nbr_v.unsigned_abs() as usize - 1;
            if ring_stereo_adjusted[nbr_idx]
                || atom_visit_orders[nbr_idx] <= atom_visit_orders[*atom_idx]
            {
                continue;
            }
            let mut tag = traversal.chiral_tag_overrides[*atom_idx].unwrap_or(atom.chiral_tag);
            if nbr_v < 0 {
                tag = match tag {
                    crate::ChiralTag::TetrahedralCw => crate::ChiralTag::TetrahedralCcw,
                    crate::ChiralTag::TetrahedralCcw => crate::ChiralTag::TetrahedralCw,
                    crate::ChiralTag::Unspecified => crate::ChiralTag::Unspecified,
                };
            }
            let first_swapped =
                traversal.chiral_tag_overrides[*atom_idx].is_some_and(|t| t != atom.chiral_tag);
            let nbr_swapped = traversal.chiral_tag_overrides[nbr_idx]
                .is_some_and(|t| t != mol.atoms[nbr_idx].chiral_tag);
            if first_swapped ^ nbr_swapped {
                tag = match tag {
                    crate::ChiralTag::TetrahedralCw => crate::ChiralTag::TetrahedralCcw,
                    crate::ChiralTag::TetrahedralCcw => crate::ChiralTag::TetrahedralCw,
                    crate::ChiralTag::Unspecified => crate::ChiralTag::Unspecified,
                };
            }
            traversal.chiral_tag_overrides[nbr_idx] = Some(tag);
            ring_stereo_adjusted[nbr_idx] = true;
        }
    }

    let mut tmp_mol = mol.clone();
    let mut bond_dir_counts = vec![0i8; tmp_mol.bonds.len()];
    let mut atom_dir_counts = vec![0i8; tmp_mol.atoms.len()];
    canonicalize_double_bonds(
        &mut tmp_mol,
        &bond_visit_orders,
        &atom_visit_orders,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        &traversal.mol_stack,
        &traversal_ring_closure_bonds,
    )?;
    remove_unwanted_bond_dir_specs(
        &mut tmp_mol,
        &traversal.mol_stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        &bond_visit_orders,
    );
    remove_redundant_bond_dir_specs(
        &mut tmp_mol,
        &traversal.mol_stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
    );
    for bond in &tmp_mol.bonds {
        if !matches!(bond.direction, crate::BondDirection::None) {
            traversal.bond_direction_overrides[bond.index] = Some(bond.direction);
        }
    }
    Ok(traversal)
}

fn chiral_atom_needs_tag_inversion(
    mol: &Molecule,
    atom_idx: usize,
    is_atom_first: bool,
    num_closures: usize,
    implicit_hydrogens: &[u8],
) -> bool {
    let degree = mol
        .bonds
        .iter()
        .filter(|bond| bond.begin_atom == atom_idx || bond.end_atom == atom_idx)
        .count();
    if degree != 3 {
        return false;
    }
    let atom = &mol.atoms[atom_idx];
    (is_atom_first && atom.explicit_hydrogens == 1)
        || (!atom_has_fourth_valence(atom, implicit_hydrogens[atom_idx])
            && num_closures == 1
            && !is_unsaturated_atom(mol, atom_idx))
}

fn atom_has_fourth_valence(atom: &crate::Atom, implicit_hydrogens: u8) -> bool {
    atom.explicit_hydrogens == 1 || implicit_hydrogens == 1
}

fn is_unsaturated_atom(mol: &Molecule, atom_idx: usize) -> bool {
    mol.bonds.iter().any(|bond| {
        (bond.begin_atom == atom_idx || bond.end_atom == atom_idx)
            && matches!(
                bond.order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
            )
    })
}

fn flip_stereo_bond_dir(direction: crate::BondDirection) -> crate::BondDirection {
    match direction {
        crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
        crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
        crate::BondDirection::None => crate::BondDirection::None,
    }
}

fn can_set_double_bond_stereo(order: BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Single | BondOrder::Aromatic | BondOrder::Dative
    )
}

fn can_have_direction(order: BondOrder) -> bool {
    matches!(order, BondOrder::Single | BondOrder::Aromatic)
}

fn bond_other_atom(bond: &crate::Bond, atom_idx: usize) -> usize {
    if bond.begin_atom == atom_idx {
        bond.end_atom
    } else {
        bond.begin_atom
    }
}

fn set_direction_from_neighboring_bond(
    mol: &mut Molecule,
    source_bond_idx: usize,
    source_flipped: bool,
    target_bond_idx: usize,
    target_flipped: bool,
) {
    let mut dir = mol.bonds[source_bond_idx].direction;
    if source_flipped == target_flipped {
        dir = flip_stereo_bond_dir(dir);
    }
    mol.bonds[target_bond_idx].direction = dir;
}

fn get_reference_direction(
    mol: &Molecule,
    dbl_bond_idx: usize,
    ref_atom_idx: usize,
    target_atom_idx: usize,
    ref_controlling_bond_idx: usize,
    ref_is_flipped: bool,
    target_bond_idx: usize,
    target_is_flipped: bool,
) -> Result<crate::BondDirection, SmilesWriteError> {
    let dbl_bond = &mol.bonds[dbl_bond_idx];
    let mut dir = match dbl_bond.stereo {
        crate::BondStereo::Trans => mol.bonds[ref_controlling_bond_idx].direction,
        crate::BondStereo::Cis => {
            flip_stereo_bond_dir(mol.bonds[ref_controlling_bond_idx].direction)
        }
        _ => crate::BondDirection::None,
    };
    if matches!(dir, crate::BondDirection::None) {
        return Err(SmilesWriteError::UnsupportedPath(
            "RDKit Canon::getReferenceDirection reached bond with unset stereo direction",
        ));
    }

    let adj = adjacency(mol);
    let ref_ctrl_other = bond_other_atom(&mol.bonds[ref_controlling_bond_idx], ref_atom_idx);
    if adj.neighbors_of(ref_atom_idx).len() == 3 && !dbl_bond.stereo_atoms.contains(&ref_ctrl_other)
    {
        dir = flip_stereo_bond_dir(dir);
    }
    let target_other = bond_other_atom(&mol.bonds[target_bond_idx], target_atom_idx);
    if adj.neighbors_of(target_atom_idx).len() == 3
        && !dbl_bond.stereo_atoms.contains(&target_other)
    {
        dir = flip_stereo_bond_dir(dir);
    }
    if ref_is_flipped != target_is_flipped {
        dir = flip_stereo_bond_dir(dir);
    }
    Ok(dir)
}

fn same_side_dirs_are_compatible(
    mol: &Molecule,
    first_bond_idx: usize,
    second_bond_idx: usize,
    first_flipped: bool,
    second_flipped: bool,
) -> bool {
    let dirs_should_match = first_flipped != second_flipped;
    let dirs_match = mol.bonds[first_bond_idx].direction == mol.bonds[second_bond_idx].direction;
    dirs_match == dirs_should_match
}

fn get_neighboring_stereo_bond(
    mol: &Molecule,
    dbl_bond_atom_idx: usize,
    nbr_bond_idx: usize,
) -> Option<usize> {
    let other_atom_idx = bond_other_atom(&mol.bonds[nbr_bond_idx], dbl_bond_atom_idx);
    for nbr in adjacency(mol).neighbors_of(other_atom_idx) {
        if nbr.bond_index != nbr_bond_idx
            && matches!(mol.bonds[nbr.bond_index].order, BondOrder::Double)
            && matches!(
                mol.bonds[nbr.bond_index].stereo,
                crate::BondStereo::Cis | crate::BondStereo::Trans
            )
        {
            return Some(nbr.bond_index);
        }
    }
    None
}

fn find_neighbor_bonds(
    mol: &Molecule,
    dbl_bond_idx: usize,
    atom_idx: usize,
    bond_dir_counts: &[i8],
    bond_visit_orders: &[usize],
) -> (Option<usize>, Option<usize>, bool) {
    let mut first_neighbor_bond = None;
    let mut second_neighbor_bond = None;
    let mut dir_set = false;
    let mut first_visit_order = mol.bonds.len() + 1;
    for nbr in adjacency(mol).neighbors_of(atom_idx) {
        if nbr.bond_index == dbl_bond_idx
            || !can_set_double_bond_stereo(mol.bonds[nbr.bond_index].order)
        {
            continue;
        }
        if bond_dir_counts[nbr.bond_index] > 0 {
            dir_set = true;
        }
        if first_neighbor_bond.is_none() || bond_visit_orders[nbr.bond_index] < first_visit_order {
            if let Some(prev) = first_neighbor_bond {
                second_neighbor_bond = Some(prev);
            }
            first_neighbor_bond = Some(nbr.bond_index);
            first_visit_order = bond_visit_orders[nbr.bond_index];
        } else {
            second_neighbor_bond = Some(nbr.bond_index);
        }
    }
    (first_neighbor_bond, second_neighbor_bond, dir_set)
}

#[allow(clippy::too_many_arguments)]
fn fix_conflict_across_double_bond(
    mol: &mut Molecule,
    dbl_bond_idx: usize,
    atom_idx: usize,
    first_bond_idx: usize,
    first_flipped: bool,
    second_bond_idx: usize,
    second_flipped: bool,
    ref_atom_idx: usize,
    ref_bond_idx: usize,
    ref_flipped: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> Result<bool, SmilesWriteError> {
    for (bond_idx, is_flipped, other_bond_idx) in [
        (first_bond_idx, first_flipped, second_bond_idx),
        (second_bond_idx, second_flipped, first_bond_idx),
    ] {
        let other_idx = bond_other_atom(&mol.bonds[other_bond_idx], atom_idx);
        if atom_dir_counts[other_idx] != 2 {
            continue;
        }
        let expected_atom_dir = get_reference_direction(
            mol,
            dbl_bond_idx,
            ref_atom_idx,
            atom_idx,
            ref_bond_idx,
            ref_flipped,
            bond_idx,
            is_flipped,
        )?;
        if expected_atom_dir == mol.bonds[bond_idx].direction {
            bond_dir_counts[other_bond_idx] = 0;
            atom_dir_counts[atom_idx] -= 1;
            atom_dir_counts[other_idx] -= 1;
            return Ok(true);
        }
    }
    Ok(false)
}

#[allow(clippy::too_many_arguments)]
fn handle_dir_conflicts_across_double_bond(
    mol: &mut Molecule,
    dbl_bond_idx: usize,
    atom1_idx: usize,
    atom1_dirs_are_consistent: bool,
    first_from_atom1_idx: usize,
    is_first_from_atom1_flipped: bool,
    second_from_atom1_idx: usize,
    is_second_from_atom1_flipped: bool,
    atom2_idx: usize,
    atom2_dirs_are_consistent: bool,
    first_from_atom2_idx: usize,
    is_first_from_atom2_flipped: bool,
    second_from_atom2_idx: usize,
    is_second_from_atom2_flipped: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> Result<bool, SmilesWriteError> {
    if atom1_dirs_are_consistent && atom2_dirs_are_consistent {
        let expected = get_reference_direction(
            mol,
            dbl_bond_idx,
            atom1_idx,
            atom2_idx,
            first_from_atom1_idx,
            is_first_from_atom1_flipped,
            first_from_atom2_idx,
            is_first_from_atom2_flipped,
        )?;
        Ok(expected == mol.bonds[first_from_atom2_idx].direction)
    } else if !atom2_dirs_are_consistent && atom1_dirs_are_consistent {
        fix_conflict_across_double_bond(
            mol,
            dbl_bond_idx,
            atom2_idx,
            first_from_atom2_idx,
            is_first_from_atom2_flipped,
            second_from_atom2_idx,
            is_second_from_atom2_flipped,
            atom1_idx,
            first_from_atom1_idx,
            is_first_from_atom1_flipped,
            bond_dir_counts,
            atom_dir_counts,
        )
    } else if !atom1_dirs_are_consistent && atom2_dirs_are_consistent {
        fix_conflict_across_double_bond(
            mol,
            dbl_bond_idx,
            atom1_idx,
            first_from_atom1_idx,
            is_first_from_atom1_flipped,
            second_from_atom1_idx,
            is_second_from_atom1_flipped,
            atom2_idx,
            first_from_atom2_idx,
            is_first_from_atom2_flipped,
            bond_dir_counts,
            atom_dir_counts,
        )
    } else {
        for (atom1_bond_idx, atom1_flipped, atom1_other_idx) in [
            (
                first_from_atom1_idx,
                is_first_from_atom1_flipped,
                second_from_atom1_idx,
            ),
            (
                second_from_atom1_idx,
                is_second_from_atom1_flipped,
                first_from_atom1_idx,
            ),
        ] {
            for (atom2_bond_idx, atom2_flipped, atom2_other_idx) in [
                (
                    first_from_atom2_idx,
                    is_first_from_atom2_flipped,
                    second_from_atom2_idx,
                ),
                (
                    second_from_atom2_idx,
                    is_second_from_atom2_flipped,
                    first_from_atom2_idx,
                ),
            ] {
                let expected = get_reference_direction(
                    mol,
                    dbl_bond_idx,
                    atom1_idx,
                    atom2_idx,
                    atom1_bond_idx,
                    atom1_flipped,
                    atom2_bond_idx,
                    atom2_flipped,
                )?;
                if expected != mol.bonds[atom2_bond_idx].direction {
                    continue;
                }
                let atom1_other_atom_idx = bond_other_atom(&mol.bonds[atom1_other_idx], atom1_idx);
                if atom_dir_counts[atom1_other_atom_idx] != 2 {
                    continue;
                }
                let atom2_other_atom_idx = bond_other_atom(&mol.bonds[atom2_other_idx], atom2_idx);
                if atom1_other_atom_idx == atom2_other_atom_idx
                    || atom_dir_counts[atom2_other_atom_idx] != 2
                {
                    continue;
                }
                bond_dir_counts[atom1_other_idx] = 0;
                atom_dir_counts[atom1_idx] -= 1;
                atom_dir_counts[atom1_other_atom_idx] -= 1;
                bond_dir_counts[atom2_other_idx] = 0;
                atom_dir_counts[atom2_idx] -= 1;
                atom_dir_counts[atom2_other_atom_idx] -= 1;
                return Ok(true);
            }
        }
        Ok(false)
    }
}

fn clear_bond_dirs(
    mol: &mut Molecule,
    ref_bond_idx: usize,
    from_atom_idx: usize,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    let clear_direction = |mol: &mut Molecule,
                           bond_idx: usize,
                           from_atom_idx: usize,
                           bond_dir_counts: &mut [i8],
                           atom_dir_counts: &mut [i8]| {
        bond_dir_counts[bond_idx] -= 1;
        if bond_dir_counts[bond_idx] == 0 {
            mol.bonds[bond_idx].direction = crate::BondDirection::None;
            atom_dir_counts[from_atom_idx] -= 1;
            let other_atom_idx = bond_other_atom(&mol.bonds[bond_idx], from_atom_idx);
            if atom_dir_counts[other_atom_idx] > 0 {
                atom_dir_counts[other_atom_idx] -= 1;
            }
        }
    };

    for nbr in adjacency(mol).neighbors_of(from_atom_idx) {
        if nbr.bond_index != ref_bond_idx && can_have_direction(mol.bonds[nbr.bond_index].order) {
            if bond_dir_counts[nbr.bond_index] >= bond_dir_counts[ref_bond_idx]
                && atom_dir_counts[mol.bonds[nbr.bond_index].begin_atom] != 1
                && atom_dir_counts[mol.bonds[nbr.bond_index].end_atom] != 1
            {
                clear_direction(
                    mol,
                    nbr.bond_index,
                    from_atom_idx,
                    bond_dir_counts,
                    atom_dir_counts,
                );
            } else if atom_dir_counts[mol.bonds[ref_bond_idx].begin_atom] != 1
                && atom_dir_counts[mol.bonds[ref_bond_idx].end_atom] != 1
            {
                clear_direction(
                    mol,
                    ref_bond_idx,
                    from_atom_idx,
                    bond_dir_counts,
                    atom_dir_counts,
                );
            }
            break;
        }
    }
}

fn canonicalize_double_bond(
    mol: &mut Molecule,
    dbl_bond_idx: usize,
    bond_visit_orders: &[usize],
    atom_visit_orders: &[usize],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
    traversal_ring_closure_bonds: &[bool],
) -> Result<(), SmilesWriteError> {
    let dbl_bond = &mol.bonds[dbl_bond_idx];
    if !matches!(dbl_bond.order, BondOrder::Double)
        || !matches!(
            dbl_bond.stereo,
            crate::BondStereo::Cis | crate::BondStereo::Trans
        )
        || dbl_bond.stereo_atoms.len() < 2
    {
        return Ok(());
    }
    let adj = adjacency(mol);
    let begin_degree = adj.neighbors_of(dbl_bond.begin_atom).len();
    let end_degree = adj.neighbors_of(dbl_bond.end_atom).len();
    if !matches!(begin_degree, 2 | 3) || !matches!(end_degree, 2 | 3) {
        return Ok(());
    }

    let (mut atom1_idx, mut atom2_idx) = (dbl_bond.begin_atom, dbl_bond.end_atom);
    if atom_visit_orders[dbl_bond.begin_atom] >= atom_visit_orders[dbl_bond.end_atom] {
        std::mem::swap(&mut atom1_idx, &mut atom2_idx);
    }

    let (first_from_atom1, second_from_atom1, dir1_set) = find_neighbor_bonds(
        mol,
        dbl_bond_idx,
        atom1_idx,
        bond_dir_counts,
        bond_visit_orders,
    );
    let (first_from_atom2, second_from_atom2, dir2_set) = find_neighbor_bonds(
        mol,
        dbl_bond_idx,
        atom2_idx,
        bond_dir_counts,
        bond_visit_orders,
    );
    let (Some(first_from_atom1_idx), Some(first_from_atom2_idx)) =
        (first_from_atom1, first_from_atom2)
    else {
        return Ok(());
    };

    let is_first_from_atom1_flipped = {
        let anchor_idx = bond_other_atom(&mol.bonds[first_from_atom1_idx], atom1_idx);
        (atom_visit_orders[atom1_idx] < atom_visit_orders[anchor_idx])
            != traversal_ring_closure_bonds[first_from_atom1_idx]
    };
    let is_second_from_atom1_flipped = if let Some(second_from_atom1_idx) = second_from_atom1 {
        let anchor_idx = bond_other_atom(&mol.bonds[second_from_atom1_idx], atom1_idx);
        (atom_visit_orders[atom1_idx] < atom_visit_orders[anchor_idx])
            != traversal_ring_closure_bonds[second_from_atom1_idx]
    } else {
        false
    };
    let is_first_from_atom2_flipped = {
        let anchor_idx = bond_other_atom(&mol.bonds[first_from_atom2_idx], atom2_idx);
        (atom_visit_orders[anchor_idx] < atom_visit_orders[atom2_idx])
            != traversal_ring_closure_bonds[first_from_atom2_idx]
    };
    let is_second_from_atom2_flipped = if let Some(second_from_atom2_idx) = second_from_atom2 {
        let anchor_idx = bond_other_atom(&mol.bonds[second_from_atom2_idx], atom2_idx);
        (atom_visit_orders[anchor_idx] < atom_visit_orders[atom2_idx])
            != traversal_ring_closure_bonds[second_from_atom2_idx]
    } else {
        false
    };

    if dir1_set && dir2_set {
        let mut atom1_dirs_are_consistent = true;
        if let Some(second_from_atom1_idx) = second_from_atom1 {
            if bond_dir_counts[first_from_atom1_idx] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    second_from_atom1_idx,
                    is_second_from_atom1_flipped,
                    first_from_atom1_idx,
                    is_first_from_atom1_flipped,
                );
            } else if bond_dir_counts[second_from_atom1_idx] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    first_from_atom1_idx,
                    is_first_from_atom1_flipped,
                    second_from_atom1_idx,
                    is_second_from_atom1_flipped,
                );
            } else {
                atom1_dirs_are_consistent = same_side_dirs_are_compatible(
                    mol,
                    first_from_atom1_idx,
                    second_from_atom1_idx,
                    is_first_from_atom1_flipped,
                    is_second_from_atom1_flipped,
                );
            }
            bond_dir_counts[second_from_atom1_idx] += 1;
            atom_dir_counts[atom1_idx] += 1;
        }
        bond_dir_counts[first_from_atom1_idx] += 1;
        atom_dir_counts[atom1_idx] += 1;

        let mut atom2_dirs_are_consistent = true;
        if let Some(second_from_atom2_idx) = second_from_atom2 {
            if bond_dir_counts[first_from_atom2_idx] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    second_from_atom2_idx,
                    is_second_from_atom2_flipped,
                    first_from_atom2_idx,
                    is_first_from_atom2_flipped,
                );
            } else if bond_dir_counts[second_from_atom2_idx] == 0 {
                set_direction_from_neighboring_bond(
                    mol,
                    first_from_atom2_idx,
                    is_first_from_atom2_flipped,
                    second_from_atom2_idx,
                    is_second_from_atom2_flipped,
                );
            } else {
                atom2_dirs_are_consistent = same_side_dirs_are_compatible(
                    mol,
                    first_from_atom2_idx,
                    second_from_atom2_idx,
                    is_first_from_atom2_flipped,
                    is_second_from_atom2_flipped,
                );
            }
            bond_dir_counts[second_from_atom2_idx] += 1;
            atom_dir_counts[atom2_idx] += 1;
        }
        bond_dir_counts[first_from_atom2_idx] += 1;
        atom_dir_counts[atom2_idx] += 1;

        if let (Some(second_from_atom1_idx), Some(second_from_atom2_idx)) =
            (second_from_atom1, second_from_atom2)
        {
            let _ = handle_dir_conflicts_across_double_bond(
                mol,
                dbl_bond_idx,
                atom1_idx,
                atom1_dirs_are_consistent,
                first_from_atom1_idx,
                is_first_from_atom1_flipped,
                second_from_atom1_idx,
                is_second_from_atom1_flipped,
                atom2_idx,
                atom2_dirs_are_consistent,
                first_from_atom2_idx,
                is_first_from_atom2_flipped,
                second_from_atom2_idx,
                is_second_from_atom2_flipped,
                bond_dir_counts,
                atom_dir_counts,
            )?;
        }
        return Ok(());
    }

    let mut set_from_bond1 = true;
    let mut atom1_controlling_bond_idx = first_from_atom1_idx;
    let mut atom2_controlling_bond_idx = first_from_atom2_idx;
    if !dir1_set && !dir2_set {
        mol.bonds[first_from_atom1_idx].direction = crate::BondDirection::EndUpRight;
        bond_dir_counts[first_from_atom1_idx] += 1;
        atom_dir_counts[atom1_idx] += 1;
    } else if !dir2_set {
        if bond_dir_counts[first_from_atom1_idx] > 0 {
            bond_dir_counts[first_from_atom1_idx] += 1;
            atom_dir_counts[atom1_idx] += 1;
            if let Some(second_from_atom1_idx) = second_from_atom1
                && bond_dir_counts[second_from_atom1_idx] > 0
            {
                let _ = same_side_dirs_are_compatible(
                    mol,
                    first_from_atom1_idx,
                    second_from_atom1_idx,
                    is_first_from_atom1_flipped,
                    is_second_from_atom1_flipped,
                );
                bond_dir_counts[second_from_atom1_idx] += 1;
                atom_dir_counts[atom1_idx] += 1;
            }
        } else {
            let Some(second_from_atom1_idx) = second_from_atom1 else {
                return Err(SmilesWriteError::UnsupportedPath(
                    "RDKit Canon::canonicalizeDoubleBond hit inconsistent atom1 stereo state",
                ));
            };
            if bond_dir_counts[second_from_atom1_idx] <= 0 {
                return Err(SmilesWriteError::UnsupportedPath(
                    "RDKit Canon::canonicalizeDoubleBond expected atom1 second bond direction",
                ));
            }
            set_direction_from_neighboring_bond(
                mol,
                second_from_atom1_idx,
                is_second_from_atom1_flipped,
                first_from_atom1_idx,
                is_first_from_atom1_flipped,
            );
            bond_dir_counts[second_from_atom1_idx] += 1;
            bond_dir_counts[first_from_atom1_idx] += 1;
            atom_dir_counts[atom1_idx] += 2;
            atom1_controlling_bond_idx = second_from_atom1_idx;
        }
    } else {
        set_from_bond1 = false;
        if bond_dir_counts[first_from_atom2_idx] > 0 {
            bond_dir_counts[first_from_atom2_idx] += 1;
            atom_dir_counts[atom2_idx] += 1;
            if let Some(second_from_atom2_idx) = second_from_atom2
                && bond_dir_counts[second_from_atom2_idx] > 0
            {
                let _ = same_side_dirs_are_compatible(
                    mol,
                    first_from_atom2_idx,
                    second_from_atom2_idx,
                    is_first_from_atom2_flipped,
                    is_second_from_atom2_flipped,
                );
                bond_dir_counts[second_from_atom2_idx] += 1;
                atom_dir_counts[atom2_idx] += 1;
            }
        } else {
            let Some(second_from_atom2_idx) = second_from_atom2 else {
                return Err(SmilesWriteError::UnsupportedPath(
                    "RDKit Canon::canonicalizeDoubleBond hit inconsistent atom2 stereo state",
                ));
            };
            if bond_dir_counts[second_from_atom2_idx] <= 0 {
                return Err(SmilesWriteError::UnsupportedPath(
                    "RDKit Canon::canonicalizeDoubleBond expected atom2 second bond direction",
                ));
            }
            set_direction_from_neighboring_bond(
                mol,
                second_from_atom2_idx,
                is_second_from_atom2_flipped,
                first_from_atom2_idx,
                is_first_from_atom2_flipped,
            );
            bond_dir_counts[second_from_atom2_idx] += 1;
            bond_dir_counts[first_from_atom2_idx] += 1;
            atom_dir_counts[atom2_idx] += 2;
            atom2_controlling_bond_idx = second_from_atom2_idx;
        }
    }

    if set_from_bond1 {
        let is_ctrl_flipped = if atom1_controlling_bond_idx == first_from_atom1_idx {
            is_first_from_atom1_flipped
        } else {
            is_second_from_atom1_flipped
        };
        let atom2_dir = get_reference_direction(
            mol,
            dbl_bond_idx,
            atom1_idx,
            atom2_idx,
            atom1_controlling_bond_idx,
            is_ctrl_flipped,
            first_from_atom2_idx,
            is_first_from_atom2_flipped,
        )?;
        mol.bonds[first_from_atom2_idx].direction = atom2_dir;
        bond_dir_counts[first_from_atom2_idx] += 1;
        atom_dir_counts[atom2_idx] += 1;
    } else {
        let is_ctrl_flipped = if atom2_controlling_bond_idx == first_from_atom2_idx {
            is_first_from_atom2_flipped
        } else {
            is_second_from_atom2_flipped
        };
        let atom1_dir = get_reference_direction(
            mol,
            dbl_bond_idx,
            atom2_idx,
            atom1_idx,
            atom2_controlling_bond_idx,
            is_ctrl_flipped,
            first_from_atom1_idx,
            is_first_from_atom1_flipped,
        )?;
        mol.bonds[first_from_atom1_idx].direction = atom1_dir;
        bond_dir_counts[first_from_atom1_idx] += 1;
        atom_dir_counts[atom1_idx] += 1;
    }

    if begin_degree == 3
        && let Some(second_from_atom1_idx) = second_from_atom1
        && bond_dir_counts[second_from_atom1_idx] == 0
    {
        set_direction_from_neighboring_bond(
            mol,
            first_from_atom1_idx,
            is_first_from_atom1_flipped,
            second_from_atom1_idx,
            is_second_from_atom1_flipped,
        );
        bond_dir_counts[second_from_atom1_idx] += 1;
        atom_dir_counts[atom1_idx] += 1;
    }

    if end_degree == 3
        && let Some(second_from_atom2_idx) = second_from_atom2
        && bond_dir_counts[second_from_atom2_idx] == 0
    {
        set_direction_from_neighboring_bond(
            mol,
            first_from_atom2_idx,
            is_first_from_atom2_flipped,
            second_from_atom2_idx,
            is_second_from_atom2_flipped,
        );
        bond_dir_counts[second_from_atom2_idx] += 1;
        atom_dir_counts[atom2_idx] += 1;
    }

    Ok(())
}

fn canonicalize_double_bonds(
    mol: &mut Molecule,
    bond_visit_orders: &[usize],
    atom_visit_orders: &[usize],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
    mol_stack: &[MolStackElem],
    traversal_ring_closure_bonds: &[bool],
) -> Result<(), SmilesWriteError> {
    use std::cmp::Reverse;
    use std::collections::{BinaryHeap, VecDeque};

    let mut stereo_bond_nbrs = vec![Vec::<usize>::new(); mol.bonds.len()];
    let mut heap = BinaryHeap::<(usize, Reverse<usize>, usize)>::new();
    for ms in mol_stack {
        let MolStackElem::Bond { bond_idx, .. } = ms else {
            continue;
        };

        if matches!(
            mol.bonds[*bond_idx].direction,
            crate::BondDirection::EndDownRight | crate::BondDirection::EndUpRight
        ) {
            mol.bonds[*bond_idx].direction = crate::BondDirection::None;
        }

        if !matches!(mol.bonds[*bond_idx].order, BondOrder::Double)
            || !matches!(
                mol.bonds[*bond_idx].stereo,
                crate::BondStereo::Cis | crate::BondStereo::Trans
            )
            || mol.bonds[*bond_idx].stereo_atoms.len() < 2
        {
            mol.bonds[*bond_idx].stereo = crate::BondStereo::None;
            continue;
        }

        for &dbl_bond_atom_idx in &[
            mol.bonds[*bond_idx].begin_atom,
            mol.bonds[*bond_idx].end_atom,
        ] {
            for nbr in adjacency(mol).neighbors_of(dbl_bond_atom_idx) {
                if !can_have_direction(mol.bonds[nbr.bond_index].order) {
                    continue;
                }
                if let Some(nbr_stereo_bond_idx) =
                    get_neighboring_stereo_bond(mol, dbl_bond_atom_idx, nbr.bond_index)
                {
                    stereo_bond_nbrs[*bond_idx].push(nbr_stereo_bond_idx);
                }
            }
        }
        stereo_bond_nbrs[*bond_idx].sort_by_key(|&idx| bond_visit_orders[idx]);
        heap.push((
            stereo_bond_nbrs[*bond_idx].len(),
            Reverse(bond_visit_orders[*bond_idx]),
            *bond_idx,
        ));
    }

    let mut seen_bonds = vec![false; mol.bonds.len()];
    while let Some((_nbr_count, _visit_order, bond_idx)) = heap.pop() {
        if seen_bonds[bond_idx] {
            continue;
        }

        let mut connected_bonds_q = VecDeque::new();
        connected_bonds_q.push_back(bond_idx);
        while let Some(current_bond_idx) = connected_bonds_q.pop_front() {
            if seen_bonds[current_bond_idx] {
                continue;
            }
            canonicalize_double_bond(
                mol,
                current_bond_idx,
                bond_visit_orders,
                atom_visit_orders,
                bond_dir_counts,
                atom_dir_counts,
                traversal_ring_closure_bonds,
            )?;
            seen_bonds[current_bond_idx] = true;
            for &nbr_stereo_bond_idx in &stereo_bond_nbrs[current_bond_idx] {
                if !seen_bonds[nbr_stereo_bond_idx] {
                    connected_bonds_q.push_back(nbr_stereo_bond_idx);
                }
            }
        }
    }

    Ok(())
}

fn remove_unwanted_bond_dir_specs(
    mol: &mut Molecule,
    mol_stack: &[MolStackElem],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
    bond_visit_orders: &[usize],
) {
    for ms in mol_stack {
        let MolStackElem::Bond { bond_idx, .. } = ms else {
            continue;
        };

        if !matches!(mol.bonds[*bond_idx].order, BondOrder::Double)
            || matches!(
                mol.bonds[*bond_idx].stereo,
                crate::BondStereo::Cis | crate::BondStereo::Trans
            )
        {
            continue;
        }

        let first_atom = mol.bonds[*bond_idx].begin_atom;
        let second_atom = mol.bonds[*bond_idx].end_atom;
        if adjacency(mol).neighbors_of(first_atom).len() == 1
            || adjacency(mol).neighbors_of(second_atom).len() == 1
        {
            continue;
        }

        let mut removal_candidates = Vec::new();
        for nbr in adjacency(mol).neighbors_of(first_atom) {
            if bond_dir_counts[nbr.bond_index] > 0 {
                removal_candidates.push(nbr.bond_index);
            }
        }
        if removal_candidates.is_empty() {
            continue;
        }
        if atom_dir_counts[first_atom] > 0 {
            removal_candidates.clear();
        }

        let mut candidates_on_second_end = 0u8;
        for nbr in adjacency(mol).neighbors_of(second_atom) {
            if bond_dir_counts[nbr.bond_index] > 0 {
                removal_candidates.push(nbr.bond_index);
                candidates_on_second_end += 1;
            }
        }
        if candidates_on_second_end == 0 {
            continue;
        }
        if atom_dir_counts[second_atom] > 0 {
            continue;
        }

        removal_candidates.sort_by_key(|&idx| bond_visit_orders[idx]);
        for candidate_bond_idx in removal_candidates {
            let other_atom = if mol.bonds[candidate_bond_idx].begin_atom == first_atom
                || mol.bonds[candidate_bond_idx].end_atom == first_atom
            {
                bond_other_atom(&mol.bonds[candidate_bond_idx], first_atom)
            } else {
                bond_other_atom(&mol.bonds[candidate_bond_idx], second_atom)
            };
            if atom_dir_counts[other_atom] == 2 {
                bond_dir_counts[candidate_bond_idx] = 0;
                mol.bonds[candidate_bond_idx].direction = crate::BondDirection::None;
                atom_dir_counts[other_atom] -= 1;
                break;
            }
        }
    }
}

fn remove_redundant_bond_dir_specs(
    mol: &mut Molecule,
    mol_stack: &[MolStackElem],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    for ms in mol_stack {
        let MolStackElem::Bond {
            bond_idx,
            atom_to_left_idx,
            ..
        } = ms
        else {
            continue;
        };

        if can_have_direction(mol.bonds[*bond_idx].order) && bond_dir_counts[*bond_idx] > 0 {
            let canon_begin_atom = *atom_to_left_idx;
            let canon_end_atom = bond_other_atom(&mol.bonds[*bond_idx], *atom_to_left_idx);
            if atom_dir_counts[canon_begin_atom] >= 2
                && adjacency(mol)
                    .neighbors_of(canon_begin_atom)
                    .iter()
                    .any(|nbr| {
                        nbr.bond_index != *bond_idx
                            && matches!(mol.bonds[nbr.bond_index].order, BondOrder::Double)
                            && matches!(
                                mol.bonds[nbr.bond_index].stereo,
                                crate::BondStereo::Cis | crate::BondStereo::Trans
                            )
                    })
            {
                clear_bond_dirs(
                    mol,
                    *bond_idx,
                    canon_begin_atom,
                    bond_dir_counts,
                    atom_dir_counts,
                );
            }
            if atom_dir_counts[canon_end_atom] >= 2
                && adjacency(mol)
                    .neighbors_of(canon_end_atom)
                    .iter()
                    .any(|nbr| {
                        nbr.bond_index != *bond_idx
                            && matches!(mol.bonds[nbr.bond_index].order, BondOrder::Double)
                            && matches!(
                                mol.bonds[nbr.bond_index].stereo,
                                crate::BondStereo::Cis | crate::BondStereo::Trans
                            )
                    })
            {
                clear_bond_dirs(
                    mol,
                    *bond_idx,
                    canon_end_atom,
                    bond_dir_counts,
                    atom_dir_counts,
                );
            }
        } else if !matches!(mol.bonds[*bond_idx].direction, crate::BondDirection::None) {
            mol.bonds[*bond_idx].direction = crate::BondDirection::None;
        }
    }
}

fn shortest_path_ignoring_bond(
    mol: &Molecule,
    begin: usize,
    end: usize,
    skip_bond: usize,
) -> Option<Vec<usize>> {
    let adj = adjacency(mol);
    let mut prev = vec![None::<usize>; mol.atoms.len()];
    let mut seen = vec![false; mol.atoms.len()];
    let mut q = std::collections::VecDeque::new();
    seen[begin] = true;
    q.push_back(begin);
    while let Some(curr) = q.pop_front() {
        for nbr in adj.neighbors_of(curr) {
            if nbr.bond_index == skip_bond || seen[nbr.atom_index] {
                continue;
            }
            seen[nbr.atom_index] = true;
            prev[nbr.atom_index] = Some(curr);
            if nbr.atom_index == end {
                let mut path = vec![end];
                let mut p = end;
                while let Some(parent) = prev[p] {
                    path.push(parent);
                    if parent == begin {
                        break;
                    }
                    p = parent;
                }
                path.reverse();
                return Some(path);
            }
            q.push_back(nbr.atom_index);
        }
    }
    None
}

fn cycle_key(cycle: &[usize]) -> Vec<usize> {
    let mut key = cycle.to_vec();
    key.sort_unstable();
    key
}

fn all_cycle_candidates(mol: &Molecule) -> Vec<Vec<usize>> {
    let mut out = Vec::<Vec<usize>>::new();
    let mut seen = std::collections::BTreeSet::<Vec<usize>>::new();
    for bond in &mol.bonds {
        if let Some(path) =
            shortest_path_ignoring_bond(mol, bond.begin_atom, bond.end_atom, bond.index)
        {
            if path.len() < 3 {
                continue;
            }
            let key = cycle_key(&path);
            if seen.insert(key) {
                out.push(path);
            }
        }
    }
    out
}

fn atom_is_in_ring_size(mol: &Molecule, atom_idx: usize, ring_size: usize) -> bool {
    all_cycle_candidates(mol)
        .into_iter()
        .any(|cycle| cycle.len() == ring_size && cycle.contains(&atom_idx))
}

fn atom_ring_count(mol: &Molecule, atom_idx: usize) -> usize {
    all_cycle_candidates(mol)
        .into_iter()
        .filter(|cycle| cycle.contains(&atom_idx))
        .count()
}

fn bond_ring_count(mol: &Molecule, bond_idx: usize) -> usize {
    let bond = &mol.bonds[bond_idx];
    all_cycle_candidates(mol)
        .into_iter()
        .filter(|cycle| {
            cycle.windows(2).any(|w| {
                (w[0] == bond.begin_atom && w[1] == bond.end_atom)
                    || (w[0] == bond.end_atom && w[1] == bond.begin_atom)
            }) || {
                let first = cycle[0];
                let last = *cycle.last().unwrap_or(&first);
                (first == bond.begin_atom && last == bond.end_atom)
                    || (first == bond.end_atom && last == bond.begin_atom)
            }
        })
        .count()
}

fn atom_is_candidate_for_ring_stereochem(
    mol: &Molecule,
    atom_idx: usize,
    atom_ranks: &[u32],
) -> bool {
    if atom_ring_count(mol, atom_idx) == 0 {
        return false;
    }
    let atom = &mol.atoms[atom_idx];
    if atom.atomic_num == 7 {
        let total_degree =
            adjacency(mol).neighbors_of(atom_idx).len() + atom.explicit_hydrogens as usize;
        if total_degree == 3 && !atom_is_in_ring_size(mol, atom_idx, 3) {
            return false;
        }
    }
    let mut non_ring_nbrs = Vec::<usize>::new();
    let mut ring_nbrs = Vec::<usize>::new();
    let mut ring_nbr_ranks = std::collections::BTreeSet::<u32>::new();
    let adj = adjacency(mol);
    for nbr in adj.neighbors_of(atom_idx) {
        if bond_ring_count(mol, nbr.bond_index) == 0 {
            non_ring_nbrs.push(nbr.atom_index);
        } else {
            ring_nbrs.push(nbr.atom_index);
            ring_nbr_ranks.insert(atom_ranks[nbr.atom_index]);
        }
    }
    match non_ring_nbrs.len() {
        2 => {
            atom_ranks[non_ring_nbrs[0]] != atom_ranks[non_ring_nbrs[1]]
                && ring_nbrs.len() != ring_nbr_ranks.len()
        }
        1 => ring_nbrs.len() > ring_nbr_ranks.len(),
        0 => {
            (ring_nbrs.len() == 4 && ring_nbr_ranks.len() == 3)
                || (ring_nbrs.len() == 3 && ring_nbr_ranks.len() == 2)
        }
        _ => false,
    }
}

fn find_chiral_atom_special_cases(mol: &Molecule) -> Result<Vec<Vec<i32>>, SmilesWriteError> {
    let legacy = mol.rdkit_legacy_stereo_atom_props(false);
    let atom_ranks: Vec<u32> = legacy
        .iter()
        .map(|p| p.cip_rank.unwrap_or(0) as u32)
        .collect();
    let mut possible_special_cases = vec![Vec::<i32>::new(); mol.atoms.len()];
    let mut atoms_seen = vec![false; mol.atoms.len()];
    let mut atoms_used = vec![false; mol.atoms.len()];
    let mut bonds_seen = vec![false; mol.bonds.len()];
    let adj = adjacency(mol);
    for atom in &mol.atoms {
        if atoms_seen[atom.index] {
            continue;
        }
        if matches!(atom.chiral_tag, crate::ChiralTag::Unspecified)
            || legacy[atom.index].cip_code.is_some()
            || atom_ring_count(mol, atom.index) == 0
            || !atom_is_candidate_for_ring_stereochem(mol, atom.index, &atom_ranks)
        {
            continue;
        }
        let mut next_atoms = std::collections::VecDeque::<usize>::new();
        for nbr in adj.neighbors_of(atom.index) {
            if !bonds_seen[nbr.bond_index] {
                bonds_seen[nbr.bond_index] = true;
                if bond_ring_count(mol, nbr.bond_index) > 0 && !atoms_seen[nbr.atom_index] {
                    next_atoms.push_back(nbr.atom_index);
                    atoms_used[nbr.atom_index] = true;
                }
            }
        }
        let mut ring_stereo_atoms = possible_special_cases[atom.index].clone();
        while let Some(ratom_idx) = next_atoms.pop_front() {
            atoms_seen[ratom_idx] = true;
            let ratom = &mol.atoms[ratom_idx];
            if !matches!(ratom.chiral_tag, crate::ChiralTag::Unspecified)
                && legacy[ratom_idx].cip_code.is_none()
                && atom_is_candidate_for_ring_stereochem(mol, ratom_idx, &atom_ranks)
            {
                let same = if ratom.chiral_tag == atom.chiral_tag {
                    1
                } else {
                    -1
                };
                ring_stereo_atoms.push(same * (ratom_idx as i32 + 1));
                possible_special_cases[ratom_idx].push(same * (atom.index as i32 + 1));
            }
            for nbr in adj.neighbors_of(ratom_idx) {
                if !bonds_seen[nbr.bond_index] {
                    bonds_seen[nbr.bond_index] = true;
                    if bond_ring_count(mol, nbr.bond_index) > 0
                        && !atoms_seen[nbr.atom_index]
                        && !atoms_used[nbr.atom_index]
                    {
                        next_atoms.push_back(nbr.atom_index);
                        atoms_used[nbr.atom_index] = true;
                    }
                }
            }
        }
        if !ring_stereo_atoms.is_empty() {
            possible_special_cases[atom.index] = ring_stereo_atoms.clone();
            for i in 0..ring_stereo_atoms.len() {
                let ring_atom_entry = ring_stereo_atoms[i];
                let ring_atom_idx = ring_atom_entry.unsigned_abs() as usize - 1;
                let mut lringatoms = possible_special_cases[ring_atom_idx].clone();
                for oring_atom_entry in ring_stereo_atoms.iter().skip(i + 1).copied() {
                    let oring_atom_idx = oring_atom_entry.unsigned_abs() as usize - 1;
                    let these_different = (ring_atom_entry < 0) ^ (oring_atom_entry < 0);
                    lringatoms.push(if these_different {
                        -((oring_atom_idx as i32) + 1)
                    } else {
                        oring_atom_idx as i32 + 1
                    });
                    let mut olringatoms = possible_special_cases[oring_atom_idx].clone();
                    olringatoms.push(if these_different {
                        -((ring_atom_idx as i32) + 1)
                    } else {
                        ring_atom_idx as i32 + 1
                    });
                    possible_special_cases[oring_atom_idx] = olringatoms;
                }
                possible_special_cases[ring_atom_idx] = lringatoms;
            }
        }
        atoms_seen[atom.index] = true;
    }
    Ok(possible_special_cases)
}
