use std::collections::{BTreeSet, VecDeque};

use crate::{BondOrder, Molecule, ValenceModel, assign_valence, rdkit_valence_list};

const DIST12_DELTA: f64 = 0.01;
const DIST13_TOL: f64 = 0.04;
const GEN_DIST_TOL: f64 = 0.06;
const MAX_UPPER: f64 = 1000.0;
const H_BOND_LENGTH: f64 = 1.8;
const VDW_SCALE_15: f64 = 0.7;
const UFF_LAMBDA: f64 = 0.1332;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DgBoundsError {
    #[error("molecule has no atoms")]
    EmptyMolecule,
    #[error("Too many bonds in the molecule, cannot compute 1-4 bounds")]
    TooManyBonds,
    #[error("DG bounds generation unsupported: {0}")]
    Unsupported(String),
}

#[derive(Debug, Clone, PartialEq)]
struct BoundsMatrix {
    data: Vec<f64>,
    n: usize,
}

impl BoundsMatrix {
    fn new(n: usize) -> Self {
        Self {
            data: vec![0.0; n * n],
            n,
        }
    }

    fn num_rows(&self) -> usize {
        self.n
    }

    fn idx(&self, row: usize, col: usize) -> usize {
        row * self.n + col
    }

    fn get_val_unchecked(&self, row: usize, col: usize) -> f64 {
        self.data[self.idx(row, col)]
    }

    fn set_val_unchecked(&mut self, row: usize, col: usize, val: f64) {
        let idx = self.idx(row, col);
        self.data[idx] = val;
    }

    fn get_upper_bound(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.get_val_unchecked(i, j)
        } else {
            self.get_val_unchecked(j, i)
        }
    }

    fn set_upper_bound(&mut self, i: usize, j: usize, val: f64) {
        assert!(val >= 0.0, "Negative upper bound");
        if i < j {
            self.set_val_unchecked(i, j, val);
        } else {
            self.set_val_unchecked(j, i, val);
        }
    }

    fn get_lower_bound(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.get_val_unchecked(j, i)
        } else {
            self.get_val_unchecked(i, j)
        }
    }

    fn set_lower_bound(&mut self, i: usize, j: usize, val: f64) {
        assert!(val >= 0.0, "Negative lower bound");
        if i < j {
            self.set_val_unchecked(j, i, val);
        } else {
            self.set_val_unchecked(i, j, val);
        }
    }

    fn into_rows(self) -> Vec<Vec<f64>> {
        self.data.chunks(self.n).map(|row| row.to_vec()).collect()
    }
}

#[derive(Debug, Clone)]
struct ComputedData {
    bond_lengths: Vec<f64>,
    visited12_bounds: Vec<bool>,
    visited13_bounds: Vec<bool>,
    visited14_bounds: Vec<bool>,
    set15_atoms: Vec<bool>,
    bond_adj: Vec<i32>,
    bond_angles: Vec<f64>,
    paths14: Vec<Path14Configuration>,
    cis_paths: Vec<u64>,
    trans_paths: Vec<u64>,
}

impl ComputedData {
    fn new(n_atoms: usize, n_bonds: usize) -> Self {
        Self {
            bond_lengths: vec![0.0; n_bonds],
            visited12_bounds: vec![false; n_atoms * n_atoms],
            visited13_bounds: vec![false; n_atoms * n_atoms],
            visited14_bounds: vec![false; n_atoms * n_atoms],
            set15_atoms: vec![false; n_atoms * n_atoms],
            bond_adj: vec![-1; n_bonds * n_bonds],
            bond_angles: vec![-1.0; n_bonds * n_bonds],
            paths14: Vec::new(),
            cis_paths: Vec::new(),
            trans_paths: Vec::new(),
        }
    }

    fn bond_mat_idx(&self, n_bonds: usize, i: usize, j: usize) -> usize {
        i * n_bonds + j
    }

    fn set_bond_adj(&mut self, n_bonds: usize, i: usize, j: usize, value: i32) {
        let idx = self.bond_mat_idx(n_bonds, i, j);
        let rev = self.bond_mat_idx(n_bonds, j, i);
        self.bond_adj[idx] = value;
        self.bond_adj[rev] = value;
    }

    fn get_bond_adj(&self, n_bonds: usize, i: usize, j: usize) -> i32 {
        self.bond_adj[self.bond_mat_idx(n_bonds, i, j)]
    }

    fn set_bond_angle(&mut self, n_bonds: usize, i: usize, j: usize, value: f64) {
        let idx = self.bond_mat_idx(n_bonds, i, j);
        let rev = self.bond_mat_idx(n_bonds, j, i);
        self.bond_angles[idx] = value;
        self.bond_angles[rev] = value;
    }

    fn get_bond_angle(&self, n_bonds: usize, i: usize, j: usize) -> f64 {
        self.bond_angles[self.bond_mat_idx(n_bonds, i, j)]
    }

    fn visited_bound(&self, pid: usize, up_to_14: bool) -> bool {
        self.visited12_bounds[pid]
            || self.visited13_bounds[pid]
            || (up_to_14 && self.visited14_bounds[pid])
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum Path14Type {
    Cis,
    Trans,
    Other,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
struct Path14Configuration {
    bid1: usize,
    bid2: usize,
    bid3: usize,
    kind: Path14Type,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum RdkitHybridization {
    S,
    Sp,
    Sp2,
    Sp3,
    Sp3d,
    Sp3d2,
    Unspecified,
}

#[derive(Debug, Copy, Clone)]
struct UffAtomicParams {
    r1: f64,
    gmp_xi: f64,
}

const UFF_PARAMS: &[(&str, f64, f64)] =
    cosmolkit_macros::rdkit_uff_params!("src/data/rdkit_uff_params.tsv");

const PERIODIC_RVDW: &[(u8, f64)] =
    cosmolkit_macros::rdkit_periodic_rvdw!("src/data/rdkit_periodic_rvdw.tsv");

fn uff_params_for_label(label: &str) -> Option<UffAtomicParams> {
    UFF_PARAMS.iter().find_map(|(key, r1, gmp_xi)| {
        (*key == label).then_some(UffAtomicParams {
            r1: *r1,
            gmp_xi: *gmp_xi,
        })
    })
}

fn init_bounds_mat(mmat: &mut BoundsMatrix, default_min: f64, default_max: f64) {
    let n = mmat.num_rows();
    for i in 1..n {
        for j in 0..i {
            mmat.set_upper_bound(i, j, default_max);
            mmat.set_lower_bound(i, j, default_min);
        }
    }
}

fn triangle_smooth_bounds(bounds_mat: &mut BoundsMatrix, tol: f64) -> bool {
    let npt = bounds_mat.num_rows();
    for k in 0..npt {
        for i in 0..(npt.saturating_sub(1)) {
            if i == k {
                continue;
            }
            let mut ii = i;
            let mut ik = k;
            if ii > ik {
                std::mem::swap(&mut ii, &mut ik);
            }
            let uik = bounds_mat.get_val_unchecked(ii, ik);
            let lik = bounds_mat.get_val_unchecked(ik, ii);
            for j in (i + 1)..npt {
                if j == k {
                    continue;
                }
                let mut jj = j;
                let mut jk = k;
                if jj > jk {
                    std::mem::swap(&mut jj, &mut jk);
                }
                let ukj = bounds_mat.get_val_unchecked(jj, jk);
                let sum_uik_ukj = uik + ukj;
                if bounds_mat.get_val_unchecked(i, j) > sum_uik_ukj {
                    bounds_mat.set_val_unchecked(i, j, sum_uik_ukj);
                }

                let diff_lik_ujk = lik - ukj;
                let diff_ljk_uik = bounds_mat.get_val_unchecked(jk, jj) - uik;
                if bounds_mat.get_val_unchecked(j, i) < diff_lik_ujk {
                    bounds_mat.set_val_unchecked(j, i, diff_lik_ujk);
                } else if bounds_mat.get_val_unchecked(j, i) < diff_ljk_uik {
                    bounds_mat.set_val_unchecked(j, i, diff_ljk_uik);
                }
                let l_bound = bounds_mat.get_val_unchecked(j, i);
                let u_bound = bounds_mat.get_val_unchecked(i, j);
                if tol > 0.0
                    && (l_bound - u_bound) / l_bound > 0.0
                    && (l_bound - u_bound) / l_bound < tol
                {
                    bounds_mat.set_val_unchecked(i, j, l_bound);
                } else if l_bound - u_bound > 0.0 {
                    return false;
                }
            }
        }
    }
    true
}

fn check_and_set_bounds(i: usize, j: usize, lb: f64, ub: f64, mmat: &mut BoundsMatrix) {
    let clb = mmat.get_lower_bound(i, j);
    let cub = mmat.get_upper_bound(i, j);
    assert!(ub > lb, "upper bound not greater than lower bound");
    assert!(lb > DIST12_DELTA || clb > DIST12_DELTA, "bad lower bound");
    if clb <= DIST12_DELTA {
        mmat.set_lower_bound(i, j, lb);
    } else if lb < clb && lb > DIST12_DELTA {
        mmat.set_lower_bound(i, j, lb);
    }
    if cub >= MAX_UPPER {
        mmat.set_upper_bound(i, j, ub);
    } else if ub > cub && ub < MAX_UPPER {
        mmat.set_upper_bound(i, j, ub);
    }
}

fn graph_distance_matrix(mol: &Molecule) -> Vec<f64> {
    let n = mol.atoms.len();
    let adjacency = mol
        .adjacency
        .clone()
        .unwrap_or_else(|| crate::AdjacencyList::from_topology(n, &mol.bonds));
    let mut out = vec![f64::INFINITY; n * n];
    for start in 0..n {
        let mut q = VecDeque::new();
        let mut seen = vec![false; n];
        q.push_back((start, 0usize));
        seen[start] = true;
        out[start * n + start] = 0.0;
        while let Some((u, d)) = q.pop_front() {
            for nb in adjacency.neighbors_of(u) {
                if seen[nb.atom_index] {
                    continue;
                }
                seen[nb.atom_index] = true;
                out[start * n + nb.atom_index] = (d + 1) as f64;
                q.push_back((nb.atom_index, d + 1));
            }
        }
    }
    out
}

fn rdkit_default_valence(atomic_num: u8) -> Option<i32> {
    let vals = rdkit_valence_list(atomic_num)?;
    vals.iter().copied().find(|v| *v >= 0)
}

fn rdkit_n_outer_electrons(atomic_num: u8) -> Option<i32> {
    match atomic_num {
        0 => Some(0),
        1 => Some(1),
        2 => Some(2),
        3 => Some(1),
        4 => Some(2),
        5 => Some(3),
        6 => Some(4),
        7 => Some(5),
        8 => Some(6),
        9 => Some(7),
        10 => Some(8),
        11 => Some(1),
        12 => Some(2),
        13 => Some(3),
        14 => Some(4),
        15 => Some(5),
        16 => Some(6),
        17 => Some(7),
        18 => Some(8),
        19 => Some(1),
        20 => Some(2),
        21 => Some(3),
        22 => Some(4),
        23 => Some(5),
        24 => Some(6),
        25 => Some(7),
        26 => Some(8),
        27 => Some(9),
        28 => Some(10),
        29 => Some(11),
        30 => Some(2),
        31 => Some(3),
        32 => Some(4),
        33 => Some(5),
        34 => Some(6),
        35 => Some(7),
        36 => Some(8),
        37 => Some(1),
        38 => Some(2),
        39 => Some(3),
        40 => Some(4),
        41 => Some(5),
        42 => Some(6),
        43 => Some(7),
        44 => Some(8),
        45 => Some(9),
        46 => Some(10),
        47 => Some(11),
        48 => Some(2),
        49 => Some(3),
        50 => Some(4),
        51 => Some(5),
        52 => Some(6),
        53 => Some(7),
        54 => Some(8),
        55 => Some(1),
        56 => Some(2),
        57 => Some(3),
        58 => Some(4),
        59 => Some(3),
        60 => Some(4),
        61 => Some(5),
        62 => Some(6),
        63 => Some(7),
        64 => Some(8),
        65 => Some(9),
        66 => Some(10),
        67 => Some(11),
        68 => Some(12),
        69 => Some(13),
        70 => Some(14),
        71 => Some(15),
        72 => Some(4),
        73 => Some(5),
        74 => Some(6),
        75 => Some(7),
        76 => Some(8),
        77 => Some(9),
        78 => Some(10),
        79 => Some(11),
        80 => Some(2),
        81 => Some(3),
        82 => Some(4),
        83 => Some(5),
        84 => Some(6),
        85 => Some(7),
        86 => Some(8),
        87 => Some(1),
        88 => Some(2),
        89 => Some(3),
        90 => Some(4),
        91 => Some(3),
        92 => Some(4),
        93 => Some(5),
        94 => Some(6),
        95 => Some(7),
        96 => Some(8),
        97 => Some(9),
        98 => Some(10),
        99 => Some(11),
        100 => Some(12),
        101 => Some(13),
        102 => Some(14),
        103 => Some(15),
        104..=118 => Some(2),
        _ => None,
    }
}

fn bond_valence_contrib_for_atom(bond: &crate::Bond, atom_index: usize) -> f64 {
    if bond.begin_atom != atom_index && bond.end_atom != atom_index {
        return 0.0;
    }
    match bond.order {
        BondOrder::Null => 0.0,
        BondOrder::Single => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Aromatic => 1.5,
        BondOrder::Dative => {
            if bond.end_atom == atom_index {
                1.0
            } else {
                0.0
            }
        }
    }
}

fn count_atom_electrons_rdkit(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
    atom_index: usize,
) -> i32 {
    let atom = &mol.atoms[atom_index];
    let Some(dv) = rdkit_default_valence(atom.atomic_num) else {
        return -1;
    };
    if dv <= 1 {
        return -1;
    }
    let mut degree = atom_degree[atom_index] as i32
        + atom.explicit_hydrogens as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    for b in &mol.bonds {
        if (b.begin_atom == atom_index || b.end_atom == atom_index)
            && bond_valence_contrib_for_atom(b, atom_index) == 0.0
        {
            degree -= 1;
        }
    }
    if degree > 3 {
        return -1;
    }
    let Some(nouter) = rdkit_n_outer_electrons(atom.atomic_num) else {
        return -1;
    };
    let nlp = (nouter - dv - atom.formal_charge as i32).max(0);
    let radicals = atom.num_radical_electrons as i32;
    let mut res = (dv - degree) + nlp - radicals;
    if res > 1 {
        let n_unsaturations =
            assignment.explicit_valence[atom_index] as i32 - atom_degree[atom_index] as i32;
        if n_unsaturations > 1 {
            res = 1;
        }
    }
    res
}

fn is_atom_conjug_cand(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
    atom_index: usize,
) -> bool {
    let at = &mol.atoms[atom_index];
    if let Some(vals) = rdkit_valence_list(at.atomic_num)
        && at.formal_charge == 0
        && !vals.is_empty()
        && vals[0] >= 0
    {
        let total_valence = assignment.explicit_valence[atom_index] as i32
            + assignment.implicit_hydrogens[atom_index] as i32;
        if total_valence > vals[0] {
            return false;
        }
    }
    let nouter = rdkit_n_outer_electrons(at.atomic_num).unwrap_or(0);
    let row_ok = at.atomic_num <= 10
        || (nouter != 5 && nouter != 6)
        || (nouter == 6 && atom_degree[atom_index] < 2);
    row_ok && count_atom_electrons_rdkit(mol, assignment, atom_degree, atom_index) > 0
}

fn compute_conjugated_bonds(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
) -> Vec<bool> {
    let mut conjugated = vec![false; mol.bonds.len()];
    for (bi, b) in mol.bonds.iter().enumerate() {
        conjugated[bi] = b.is_aromatic;
    }
    for at in 0..mol.atoms.len() {
        if !is_atom_conjug_cand(mol, assignment, atom_degree, at) {
            continue;
        }
        let sbo = atom_degree[at]
            + mol.atoms[at].explicit_hydrogens as usize
            + assignment.implicit_hydrogens[at] as usize;
        if !(2..=3).contains(&sbo) {
            continue;
        }
        let bnds: Vec<usize> = mol
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(bi, b)| {
                if b.begin_atom == at || b.end_atom == at {
                    Some(bi)
                } else {
                    None
                }
            })
            .collect();
        for &b1 in &bnds {
            let bond1 = &mol.bonds[b1];
            if bond_valence_contrib_for_atom(bond1, at) < 1.5 {
                continue;
            }
            let o1 = if bond1.begin_atom == at {
                bond1.end_atom
            } else {
                bond1.begin_atom
            };
            if !is_atom_conjug_cand(mol, assignment, atom_degree, o1) {
                continue;
            }
            for &b2 in &bnds {
                if b1 == b2 {
                    continue;
                }
                let bond2 = &mol.bonds[b2];
                let o2 = if bond2.begin_atom == at {
                    bond2.end_atom
                } else {
                    bond2.begin_atom
                };
                let sbo2 = atom_degree[o2]
                    + mol.atoms[o2].explicit_hydrogens as usize
                    + assignment.implicit_hydrogens[o2] as usize;
                if sbo2 > 3 {
                    continue;
                }
                if is_atom_conjug_cand(mol, assignment, atom_degree, o2) {
                    conjugated[b1] = true;
                    conjugated[b2] = true;
                }
            }
        }
    }
    conjugated
}

fn compute_hybridizations(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
    atom_has_conjugated_bond: &[bool],
) -> Vec<RdkitHybridization> {
    let mut out = Vec::with_capacity(mol.atoms.len());
    for (atom_index, atom) in mol.atoms.iter().enumerate() {
        if atom.atomic_num == 1 {
            if !atom.no_implicit && atom.isotope.is_none() {
                out.push(RdkitHybridization::Unspecified);
                continue;
            }
        }
        if atom.atomic_num == 0 {
            out.push(RdkitHybridization::Unspecified);
            continue;
        }
        let mut deg = atom_degree[atom_index] as i32
            + atom.explicit_hydrogens as i32
            + assignment.implicit_hydrogens[atom_index] as i32;
        for b in &mol.bonds {
            if (b.begin_atom == atom_index || b.end_atom == atom_index)
                && matches!(b.order, BondOrder::Dative)
                && b.end_atom != atom_index
            {
                deg -= 1;
            }
        }
        let hyb = if atom.atomic_num <= 1 {
            match deg {
                0 | 1 => RdkitHybridization::S,
                2 => RdkitHybridization::Sp,
                3 => RdkitHybridization::Sp2,
                4 => RdkitHybridization::Sp3,
                5 => RdkitHybridization::Sp3d,
                6 => RdkitHybridization::Sp3d2,
                _ => RdkitHybridization::Unspecified,
            }
        } else {
            let nouter = rdkit_n_outer_electrons(atom.atomic_num).unwrap_or(0);
            let total_valence = assignment.explicit_valence[atom_index] as i32
                + assignment.implicit_hydrogens[atom_index] as i32;
            let num_free = nouter - (total_valence + atom.formal_charge as i32);
            let norbs = if total_valence + nouter - (atom.formal_charge as i32) < 8 {
                let radicals = atom.num_radical_electrons as i32;
                let lone_pairs = (num_free - radicals) / 2;
                deg + lone_pairs + radicals
            } else {
                let lone_pairs = num_free / 2;
                deg + lone_pairs
            };
            match norbs {
                0 | 1 => RdkitHybridization::S,
                2 => RdkitHybridization::Sp,
                3 => RdkitHybridization::Sp2,
                4 => {
                    let total_degree = atom_degree[atom_index]
                        + atom.explicit_hydrogens as usize
                        + assignment.implicit_hydrogens[atom_index] as usize;
                    if total_degree > 3 || !atom_has_conjugated_bond[atom_index] {
                        RdkitHybridization::Sp3
                    } else {
                        RdkitHybridization::Sp2
                    }
                }
                5 => RdkitHybridization::Sp3d,
                6 => RdkitHybridization::Sp3d2,
                _ => RdkitHybridization::Unspecified,
            }
        };
        out.push(hyb);
    }
    out
}

fn atom_total_valence(
    _mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_index: usize,
) -> i32 {
    assignment.explicit_valence[atom_index] as i32
        + assignment.implicit_hydrogens[atom_index] as i32
}

fn add_atom_charge_flags(
    atom_index: usize,
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_key: &mut String,
) {
    let atom = &mol.atoms[atom_index];
    let total_valence = atom_total_valence(mol, assignment, atom_index);
    match atom.atomic_num {
        29 | 47 => {
            if total_valence == 1 || atom.formal_charge == 1 {
                atom_key.push_str("+1");
            }
        }
        4 | 20 | 25 | 26 | 28 | 46 | 78 => {
            if total_valence == 2 || atom.formal_charge == 2 {
                atom_key.push_str("+2");
            }
        }
        21 | 24 | 27 | 79 | 89 | 96..=103 => {
            if total_valence == 3 || atom.formal_charge == 3 {
                atom_key.push_str("+3");
            }
        }
        2 | 18 | 22 | 36 | 54 | 90 | 92 | 93 | 94 | 95 => {
            if total_valence == 4 || atom.formal_charge == 4 {
                atom_key.push_str("+4");
            }
        }
        23 | 41 | 43 | 73 => {
            if total_valence == 5 || atom.formal_charge == 5 {
                atom_key.push_str("+5");
            }
        }
        42 => {
            if total_valence == 6 || atom.formal_charge == 6 {
                atom_key.push_str("+6");
            }
        }
        12 | 30 | 48 => {
            if total_valence == 2 {
                atom_key.push_str("+2");
            }
        }
        15 => match total_valence {
            3 => atom_key.push_str("+3"),
            5 => atom_key.push_str("+5"),
            _ => {}
        },
        16 => {
            if compute_single_hybridization(mol, assignment, atom_index) != RdkitHybridization::Sp2
            {
                match total_valence {
                    2 => atom_key.push_str("+2"),
                    4 => atom_key.push_str("+4"),
                    6 => atom_key.push_str("+6"),
                    _ => {}
                }
            }
        }
        31 | 33 | 49 | 51 | 81 | 82 | 83 => {
            if total_valence == 3 {
                atom_key.push_str("+3");
            }
        }
        34 | 52 | 84 => {
            if total_valence == 2 {
                atom_key.push_str("+2");
            }
        }
        _ => {
            if (57..=71).contains(&atom.atomic_num) && total_valence == 6 {
                atom_key.push_str("+3");
            }
        }
    }
}

fn compute_single_hybridization(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_index: usize,
) -> RdkitHybridization {
    let mut atom_degree = vec![0usize; mol.atoms.len()];
    for b in &mol.bonds {
        atom_degree[b.begin_atom] += 1;
        atom_degree[b.end_atom] += 1;
    }
    let conjugated = compute_conjugated_bonds(mol, assignment, &atom_degree);
    let mut atom_has_conjugated_bond = vec![false; mol.atoms.len()];
    for (bi, bond) in mol.bonds.iter().enumerate() {
        if conjugated[bi] {
            atom_has_conjugated_bond[bond.begin_atom] = true;
            atom_has_conjugated_bond[bond.end_atom] = true;
        }
    }
    compute_hybridizations(mol, assignment, &atom_degree, &atom_has_conjugated_bond)[atom_index]
}

fn atom_symbol(atomic_num: u8) -> Option<&'static str> {
    match atomic_num {
        0 => Some("*"),
        1 => Some("H"),
        5 => Some("B"),
        6 => Some("C"),
        7 => Some("N"),
        8 => Some("O"),
        9 => Some("F"),
        11 => Some("Na"),
        14 => Some("Si"),
        15 => Some("P"),
        16 => Some("S"),
        17 => Some("Cl"),
        29 => Some("Cu"),
        34 => Some("Se"),
        35 => Some("Br"),
        45 => Some("Rh"),
        53 => Some("I"),
        _ => None,
    }
}

fn get_atom_label(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    hybridizations: &[RdkitHybridization],
    atom_has_conjugated_bond: &[bool],
    atom_index: usize,
) -> Result<String, DgBoundsError> {
    let atom = &mol.atoms[atom_index];
    let mut atom_key = atom_symbol(atom.atomic_num)
        .ok_or_else(|| {
            DgBoundsError::Unsupported(format!(
                "UFF atom symbol lookup missing for atomic number {}",
                atom.atomic_num
            ))
        })?
        .to_string();
    if atom_key.len() == 1 {
        atom_key.push('_');
    }
    if atom.atomic_num != 0 {
        let nouter = rdkit_n_outer_electrons(atom.atomic_num).unwrap_or(0);
        if rdkit_default_valence(atom.atomic_num) == Some(-1) || (nouter != 1 && nouter != 7) {
            match atom.atomic_num {
                12 | 13 | 14 | 15 | 50 | 51 | 52 | 81 | 82 | 83 | 84 => atom_key.push('3'),
                80 => atom_key.push('1'),
                _ => match hybridizations[atom_index] {
                    RdkitHybridization::S => {}
                    RdkitHybridization::Sp => atom_key.push('1'),
                    RdkitHybridization::Sp2 => {
                        if (atom.is_aromatic || atom_has_conjugated_bond[atom_index])
                            && matches!(atom.atomic_num, 6 | 7 | 8 | 16)
                        {
                            atom_key.push('R');
                        } else {
                            atom_key.push('2');
                        }
                    }
                    RdkitHybridization::Sp3 => atom_key.push('3'),
                    RdkitHybridization::Sp3d => atom_key.push('5'),
                    RdkitHybridization::Sp3d2 => atom_key.push('6'),
                    RdkitHybridization::Unspecified => {}
                },
            }
        }
    }
    add_atom_charge_flags(atom_index, mol, assignment, &mut atom_key);
    Ok(atom_key)
}

fn calc_bond_rest_length(bond_order: f64, end1: &UffAtomicParams, end2: &UffAtomicParams) -> f64 {
    let ri = end1.r1;
    let rj = end2.r1;
    let r_bo = -UFF_LAMBDA * (ri + rj) * bond_order.ln();
    let xi = end1.gmp_xi;
    let xj = end2.gmp_xi;
    let root_delta = xi.sqrt() - xj.sqrt();
    let r_en = ri * rj * root_delta * root_delta / (xi * ri + xj * rj);
    ri + rj + r_bo - r_en
}

fn shortest_path_ignoring_edge(
    adj: &[Vec<(usize, usize)>],
    src: usize,
    dst: usize,
    forbidden_edge: usize,
) -> Option<Vec<usize>> {
    let mut q = VecDeque::new();
    let mut seen = vec![false; adj.len()];
    let mut prev = vec![usize::MAX; adj.len()];
    seen[src] = true;
    q.push_back(src);
    while let Some(v) = q.pop_front() {
        if v == dst {
            break;
        }
        for &(nb, eidx) in &adj[v] {
            if eidx == forbidden_edge || seen[nb] {
                continue;
            }
            seen[nb] = true;
            prev[nb] = v;
            q.push_back(nb);
        }
    }
    if !seen[dst] {
        return None;
    }
    let mut path = Vec::new();
    let mut cur = dst;
    path.push(cur);
    while cur != src {
        cur = prev[cur];
        path.push(cur);
    }
    path.reverse();
    Some(path)
}

fn canonical_cycle(mut c: Vec<usize>) -> Vec<usize> {
    let min_pos = c
        .iter()
        .enumerate()
        .min_by_key(|(_, v)| **v)
        .map(|(i, _)| i)
        .unwrap_or(0);
    c.rotate_left(min_pos);
    let mut rev = c.clone();
    rev.reverse();
    let min_pos_rev = rev
        .iter()
        .enumerate()
        .min_by_key(|(_, v)| **v)
        .map(|(i, _)| i)
        .unwrap_or(0);
    rev.rotate_left(min_pos_rev);
    if rev < c { rev } else { c }
}

fn adjacency_vec(mol: &Molecule) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); mol.atoms.len()];
    for bond in &mol.bonds {
        adjacency[bond.begin_atom].push(bond.end_atom);
        adjacency[bond.end_atom].push(bond.begin_atom);
    }
    adjacency
}

fn connected_components(mol: &Molecule, adjacency: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; mol.atoms.len()];
    let mut comps = Vec::new();
    for start in 0..mol.atoms.len() {
        if seen[start] {
            continue;
        }
        let mut comp = Vec::new();
        let mut q = VecDeque::new();
        q.push_back(start);
        seen[start] = true;
        while let Some(cur) = q.pop_front() {
            comp.push(cur);
            for &nb in &adjacency[cur] {
                if !seen[nb] {
                    seen[nb] = true;
                    q.push_back(nb);
                }
            }
        }
        comp.sort_unstable();
        comps.push(comp);
    }
    comps.sort_by_key(|c| c[0]);
    comps
}

fn rdkit_single_ring_order(
    mol: &Molecule,
    ring_set: &BTreeSet<usize>,
    degrees_after_trim: &[usize],
) -> Option<Vec<usize>> {
    let root = ring_set
        .iter()
        .copied()
        .filter(|&aid| degrees_after_trim[aid] == 2)
        .min()
        .or_else(|| ring_set.iter().copied().min())?;
    const WHITE: u8 = 0;
    const GRAY: u8 = 1;
    const BLACK: u8 = 2;

    let mut done = vec![WHITE; mol.atoms.len()];
    let mut parents = vec![usize::MAX; mol.atoms.len()];
    let mut depths = vec![0usize; mol.atoms.len()];
    let mut bfsq = VecDeque::<usize>::new();
    bfsq.push_back(root);

    let mut rings = Vec::<Vec<usize>>::new();
    let mut cur_size = usize::MAX;
    while let Some(curr) = bfsq.pop_front() {
        done[curr] = BLACK;
        let depth = depths[curr] + 1;
        if depth > cur_size {
            break;
        }
        for bond in &mol.bonds {
            let nbr = if bond.begin_atom == curr {
                bond.end_atom
            } else if bond.end_atom == curr {
                bond.begin_atom
            } else {
                continue;
            };
            if !ring_set.contains(&nbr) {
                continue;
            }
            if done[nbr] == BLACK || parents[curr] == nbr {
                continue;
            }
            if done[nbr] == WHITE {
                parents[nbr] = curr;
                done[nbr] = GRAY;
                depths[nbr] = depth;
                bfsq.push_back(nbr);
            } else {
                let mut ring = vec![nbr];
                let mut parent = parents[nbr];
                while parent != usize::MAX && parent != root {
                    ring.push(parent);
                    parent = parents[parent];
                }
                ring.insert(0, curr);
                parent = parents[curr];
                while parent != usize::MAX {
                    if ring.contains(&parent) {
                        ring.clear();
                        break;
                    }
                    ring.insert(0, parent);
                    parent = parents[parent];
                }
                if ring.len() > 1 {
                    if ring.len() <= cur_size {
                        cur_size = ring.len();
                        rings.push(ring);
                    } else {
                        break;
                    }
                }
            }
        }
    }
    rings.into_iter().next()
}

fn rdkit_smallest_rings_bfs(
    mol: &Molecule,
    root: usize,
    active_bonds: &[bool],
    forbidden: &[bool],
) -> Vec<Vec<usize>> {
    const WHITE: u8 = 0;
    const GRAY: u8 = 1;
    const BLACK: u8 = 2;

    let mut done = vec![WHITE; mol.atoms.len()];
    for (idx, &is_forbidden) in forbidden.iter().enumerate() {
        if is_forbidden {
            done[idx] = BLACK;
        }
    }
    let mut parents = vec![usize::MAX; mol.atoms.len()];
    let mut depths = vec![0usize; mol.atoms.len()];
    let mut bfsq = VecDeque::<usize>::new();
    bfsq.push_back(root);

    let mut rings = Vec::<Vec<usize>>::new();
    let mut cur_size = usize::MAX;
    while let Some(curr) = bfsq.pop_front() {
        done[curr] = BLACK;
        let depth = depths[curr] + 1;
        if depth > cur_size {
            break;
        }

        for (bond_idx, bond) in mol.bonds.iter().enumerate() {
            if !active_bonds[bond_idx] {
                continue;
            }
            let nbr = if bond.begin_atom == curr {
                bond.end_atom
            } else if bond.end_atom == curr {
                bond.begin_atom
            } else {
                continue;
            };
            if done[nbr] == BLACK || parents[curr] == nbr {
                continue;
            }
            if done[nbr] == WHITE {
                parents[nbr] = curr;
                done[nbr] = GRAY;
                depths[nbr] = depth;
                bfsq.push_back(nbr);
            } else {
                let mut ring = vec![nbr];
                let mut parent = parents[nbr];
                while parent != usize::MAX && parent != root {
                    ring.push(parent);
                    parent = parents[parent];
                }
                ring.insert(0, curr);
                parent = parents[curr];
                while parent != usize::MAX {
                    if ring.contains(&parent) {
                        ring.clear();
                        break;
                    }
                    ring.insert(0, parent);
                    parent = parents[parent];
                }
                if ring.len() > 1 {
                    if ring.len() <= cur_size {
                        cur_size = ring.len();
                        rings.push(ring);
                    } else {
                        return rings;
                    }
                }
            }
        }
    }
    rings
}

fn rdkit_trim_bonds(
    mol: &Molecule,
    cand: usize,
    changed: &mut VecDeque<usize>,
    atom_degrees: &mut [usize],
    active_bonds: &mut [bool],
) {
    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        if !active_bonds[bond_idx] {
            continue;
        }
        let other = if bond.begin_atom == cand {
            bond.end_atom
        } else if bond.end_atom == cand {
            bond.begin_atom
        } else {
            continue;
        };
        if atom_degrees[other] <= 2 {
            changed.push_back(other);
        }
        active_bonds[bond_idx] = false;
        atom_degrees[other] = atom_degrees[other].saturating_sub(1);
        atom_degrees[cand] = atom_degrees[cand].saturating_sub(1);
    }
}

fn rdkit_mark_useless_d2s(
    mol: &Molecule,
    root: usize,
    forb: &mut [bool],
    atom_degrees: &[usize],
    active_bonds: &[bool],
) {
    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        if !active_bonds[bond_idx] {
            continue;
        }
        let other = if bond.begin_atom == root {
            bond.end_atom
        } else if bond.end_atom == root {
            bond.begin_atom
        } else {
            continue;
        };
        if !forb[other] && atom_degrees[other] == 2 {
            forb[other] = true;
            rdkit_mark_useless_d2s(mol, other, forb, atom_degrees, active_bonds);
        }
    }
}

fn rdkit_pick_d2_nodes(
    mol: &Molecule,
    comp: &[usize],
    atom_degrees: &[usize],
    active_bonds: &[bool],
) -> Vec<usize> {
    let mut d2nodes = Vec::new();
    let mut forb = vec![false; mol.atoms.len()];
    loop {
        let mut root = None;
        for &idx in comp {
            if atom_degrees[idx] == 2 && !forb[idx] {
                root = Some(idx);
                d2nodes.push(idx);
                forb[idx] = true;
                break;
            }
        }
        let Some(root) = root else {
            break;
        };
        rdkit_mark_useless_d2s(mol, root, &mut forb, atom_degrees, active_bonds);
    }
    d2nodes
}

fn rdkit_find_sssr_orders(mol: &Molecule, comp: &[usize]) -> Vec<Vec<usize>> {
    let mut atom_degrees = vec![0usize; mol.atoms.len()];
    for bond in &mol.bonds {
        atom_degrees[bond.begin_atom] += 1;
        atom_degrees[bond.end_atom] += 1;
    }
    let mut active_bonds = vec![true; mol.bonds.len()];
    let mut changed = VecDeque::<usize>::new();
    for &idx in comp {
        if atom_degrees[idx] < 2 {
            changed.push_back(idx);
        }
    }

    let mut done_atoms = vec![false; mol.atoms.len()];
    let mut n_atoms_done = 0usize;
    let mut rings = Vec::<Vec<usize>>::new();
    let mut invariants = BTreeSet::<Vec<usize>>::new();

    while n_atoms_done <= comp.len().saturating_sub(3) {
        while let Some(cand) = changed.pop_front() {
            if !done_atoms[cand] {
                done_atoms[cand] = true;
                n_atoms_done += 1;
                rdkit_trim_bonds(
                    mol,
                    cand,
                    &mut changed,
                    &mut atom_degrees,
                    &mut active_bonds,
                );
            }
        }

        let d2nodes = rdkit_pick_d2_nodes(mol, comp, &atom_degrees, &active_bonds);
        if d2nodes.is_empty() {
            break;
        }

        for &cand in &d2nodes {
            let forbidden = vec![false; mol.atoms.len()];
            for ring in rdkit_smallest_rings_bfs(mol, cand, &active_bonds, &forbidden) {
                let mut inv = ring.clone();
                inv.sort_unstable();
                if invariants.insert(inv) {
                    rings.push(ring);
                }
            }
        }

        for cand in d2nodes {
            done_atoms[cand] = true;
            n_atoms_done += 1;
            rdkit_trim_bonds(
                mol,
                cand,
                &mut changed,
                &mut atom_degrees,
                &mut active_bonds,
            );
        }
    }
    rings
}

fn rdkit_atom_rings(mol: &Molecule) -> Vec<Vec<usize>> {
    let adjacency = adjacency_vec(mol);
    let comps = connected_components(mol, &adjacency);
    let mut out = Vec::new();
    for comp in comps {
        let comp_set: BTreeSet<usize> = comp.iter().copied().collect();
        let mut deg_in_comp = vec![0usize; mol.atoms.len()];
        for &a in &comp {
            deg_in_comp[a] = adjacency[a].iter().filter(|n| comp_set.contains(n)).count();
        }
        let mut removed = BTreeSet::<usize>::new();
        let mut queue: Vec<usize> = comp
            .iter()
            .copied()
            .filter(|&a| deg_in_comp[a] <= 1)
            .collect();
        while let Some(a) = queue.pop() {
            if removed.contains(&a) {
                continue;
            }
            removed.insert(a);
            for &nb in &adjacency[a] {
                if !comp_set.contains(&nb) || removed.contains(&nb) {
                    continue;
                }
                if deg_in_comp[nb] > 0 {
                    deg_in_comp[nb] -= 1;
                    if deg_in_comp[nb] == 1 {
                        queue.push(nb);
                    }
                }
            }
        }
        let ring_atoms: Vec<usize> = comp
            .iter()
            .copied()
            .filter(|a| !removed.contains(a))
            .collect();
        if ring_atoms.is_empty() {
            continue;
        }
        let ring_set: BTreeSet<usize> = ring_atoms.iter().copied().collect();
        let ring_bond_count = mol
            .bonds
            .iter()
            .filter(|b| ring_set.contains(&b.begin_atom) && ring_set.contains(&b.end_atom))
            .count();
        let target_cycle_count = ring_bond_count + 1 - ring_atoms.len();
        if target_cycle_count == 0 {
            continue;
        }

        let mut cycles = rdkit_find_sssr_orders(mol, &comp);
        if target_cycle_count == 1 && cycles.is_empty() {
            if let Some(cycle) = rdkit_single_ring_order(mol, &ring_set, &deg_in_comp) {
                cycles.push(cycle);
            }
        } else if cycles.len() < target_cycle_count {
            let mut seen = BTreeSet::<Vec<usize>>::new();
            let mut all_cycles = Vec::<Vec<usize>>::new();
            let mut ring_atoms_sorted = ring_atoms.clone();
            ring_atoms_sorted.sort_unstable();
            const MAX_RING_ENUM: usize = 8;
            for &start in &ring_atoms_sorted {
                let mut path = vec![start];
                let mut used = BTreeSet::<usize>::new();
                used.insert(start);
                fn dfs_cycles(
                    start: usize,
                    cur: usize,
                    path: &mut Vec<usize>,
                    used: &mut BTreeSet<usize>,
                    adjacency: &[Vec<usize>],
                    ring_set: &BTreeSet<usize>,
                    seen: &mut BTreeSet<Vec<usize>>,
                    cycles: &mut Vec<Vec<usize>>,
                ) {
                    if path.len() >= 3 && adjacency[cur].contains(&start) {
                        let c = canonical_cycle(path.clone());
                        if seen.insert(c.clone()) {
                            cycles.push(c);
                        }
                    }
                    if path.len() == MAX_RING_ENUM {
                        return;
                    }
                    let mut nbs: Vec<usize> = adjacency[cur]
                        .iter()
                        .copied()
                        .filter(|n| ring_set.contains(n) && !used.contains(n))
                        .collect();
                    nbs.sort_unstable();
                    for nb in nbs {
                        if nb < start {
                            continue;
                        }
                        used.insert(nb);
                        path.push(nb);
                        dfs_cycles(start, nb, path, used, adjacency, ring_set, seen, cycles);
                        path.pop();
                        used.remove(&nb);
                    }
                }
                dfs_cycles(
                    start,
                    start,
                    &mut path,
                    &mut used,
                    &adjacency,
                    &ring_set,
                    &mut seen,
                    &mut all_cycles,
                );
            }
            all_cycles.retain(|cyc| {
                (0..cyc.len()).all(|i| {
                    let a = cyc[i];
                    let b = cyc[(i + 1) % cyc.len()];
                    adjacency[a].contains(&b) && ring_set.contains(&a) && ring_set.contains(&b)
                })
            });
            all_cycles.sort_by_key(|cyc| (cyc.len(), cyc.clone()));
            let mut covered_edges = BTreeSet::<(usize, usize)>::new();
            for cyc in all_cycles {
                let edges: Vec<(usize, usize)> = (0..cyc.len())
                    .map(|i| {
                        let a = cyc[i];
                        let b = cyc[(i + 1) % cyc.len()];
                        if a <= b { (a, b) } else { (b, a) }
                    })
                    .collect();
                if cycles.is_empty() || edges.iter().any(|edge| !covered_edges.contains(edge)) {
                    for edge in edges {
                        covered_edges.insert(edge);
                    }
                    cycles.push(cyc);
                }
                if cycles.len() == target_cycle_count {
                    break;
                }
            }
        }
        out.extend(cycles);
    }
    out.sort_by_key(|r| r.len());
    out
}

fn rdkit_bond_rings(mol: &Molecule) -> Vec<Vec<usize>> {
    rdkit_atom_rings(mol)
        .into_iter()
        .filter_map(|ring| cycle_bond_indices(mol, &ring))
        .collect()
}

fn is_atom_in_ring_of_size(rings: &[Vec<usize>], atom_index: usize, size: usize) -> bool {
    rings
        .iter()
        .any(|ring| ring.len() == size && ring.contains(&atom_index))
}

fn cycle_bond_indices(mol: &Molecule, ring: &[usize]) -> Option<Vec<usize>> {
    let mut out = Vec::with_capacity(ring.len());
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        let bond = mol.bonds.iter().find(|bond| {
            (bond.begin_atom == a && bond.end_atom == b)
                || (bond.begin_atom == b && bond.end_atom == a)
        })?;
        out.push(bond.index);
    }
    Some(out)
}

fn set_ring_angle(hyb: RdkitHybridization, ring_size: usize) -> f64 {
    if (hyb == RdkitHybridization::Sp2 && ring_size <= 8) || ring_size == 3 || ring_size == 4 {
        std::f64::consts::PI * (1.0 - 2.0 / ring_size as f64)
    } else if hyb == RdkitHybridization::Sp3 {
        if ring_size == 5 {
            104.0_f64.to_radians()
        } else {
            109.5_f64.to_radians()
        }
    } else if hyb == RdkitHybridization::Sp3d {
        105.0_f64.to_radians()
    } else if hyb == RdkitHybridization::Sp3d2 {
        90.0_f64.to_radians()
    } else {
        120.0_f64.to_radians()
    }
}

fn bond_in_ring_of_size(mol: &Molecule, bond_index: usize, size: usize) -> bool {
    let mut adj = vec![Vec::<(usize, usize)>::new(); mol.atoms.len()];
    for (bi, b) in mol.bonds.iter().enumerate() {
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }
    let bond = &mol.bonds[bond_index];
    shortest_path_ignoring_edge(&adj, bond.begin_atom, bond.end_atom, bond_index)
        .is_some_and(|path| path.len() == size)
}

fn bond_order_as_double(order: BondOrder) -> f64 {
    match order {
        BondOrder::Null => 0.0,
        BondOrder::Single | BondOrder::Dative => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Aromatic => 1.5,
    }
}

fn rvdw(atomic_num: u8) -> f64 {
    PERIODIC_RVDW
        .iter()
        .find_map(|(num, value)| (*num == atomic_num).then_some(*value))
        .unwrap_or(1.7)
}

fn set12_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
) -> Result<(), DgBoundsError> {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).map_err(|err| {
        DgBoundsError::Unsupported(format!(
            "RDKit UFF atom typing valence assignment failed: {err}"
        ))
    })?;
    let mut atom_degree = vec![0usize; mol.atoms.len()];
    for bond in &mol.bonds {
        atom_degree[bond.begin_atom] += 1;
        atom_degree[bond.end_atom] += 1;
    }
    let conjugated = compute_conjugated_bonds(mol, &assignment, &atom_degree);
    let mut atom_has_conjugated_bond = vec![false; mol.atoms.len()];
    for (bi, bond) in mol.bonds.iter().enumerate() {
        if conjugated[bi] {
            atom_has_conjugated_bond[bond.begin_atom] = true;
            atom_has_conjugated_bond[bond.end_atom] = true;
        }
    }
    let hybridizations =
        compute_hybridizations(mol, &assignment, &atom_degree, &atom_has_conjugated_bond);
    let mut atom_params = Vec::with_capacity(mol.atoms.len());
    for atom_index in 0..mol.atoms.len() {
        let label = get_atom_label(
            mol,
            &assignment,
            &hybridizations,
            &atom_has_conjugated_bond,
            atom_index,
        )?;
        atom_params.push(uff_params_for_label(&label));
    }

    let mut squish_atoms = vec![false; mol.atoms.len()];
    for (bond_index, bond) in mol.bonds.iter().enumerate() {
        if conjugated[bond_index]
            && (mol.atoms[bond.begin_atom].atomic_num > 10
                || mol.atoms[bond.end_atom].atomic_num > 10)
            && bond_in_ring_of_size(mol, bond_index, 5)
        {
            squish_atoms[bond.begin_atom] = true;
            squish_atoms[bond.end_atom] = true;
        }
    }

    for bond in &mol.bonds {
        let beg = bond.begin_atom;
        let end = bond.end_atom;
        let bond_order = bond_order_as_double(bond.order);
        if let (Some(p1), Some(p2)) = (atom_params[beg], atom_params[end]) {
            if bond_order <= 0.0 {
                return Err(DgBoundsError::Unsupported(format!(
                    "RDKit set12Bounds encountered non-positive bond order at bond {}",
                    bond.index
                )));
            }
            let bl = calc_bond_rest_length(bond_order, &p1, &p2);
            let extra_squish = if squish_atoms[beg] || squish_atoms[end] {
                0.2
            } else {
                0.0
            };
            accum_data.bond_lengths[bond.index] = bl;
            mmat.set_upper_bound(beg, end, bl + extra_squish + DIST12_DELTA);
            mmat.set_lower_bound(beg, end, bl - extra_squish - DIST12_DELTA);
        } else {
            let bl = (rvdw(mol.atoms[beg].atomic_num) + rvdw(mol.atoms[end].atomic_num)) / 2.0;
            accum_data.bond_lengths[bond.index] = bl;
            mmat.set_upper_bound(beg, end, 1.5 * bl);
            mmat.set_lower_bound(beg, end, 0.5 * bl);
        }
        let pid = beg.min(end) * mol.atoms.len() + beg.max(end);
        accum_data.visited12_bounds[pid] = true;
    }
    Ok(())
}

fn compute_13_dist(d1: f64, d2: f64, angle: f64) -> f64 {
    (d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * angle.cos()).sqrt()
}

fn compute_14_dist_3d(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64, tor_ang: f64) -> f64 {
    let p1x = d1 * ang12.cos();
    let p1y = d1 * ang12.sin();
    let p4x = d2 - d3 * ang23.cos();
    let p4y = d3 * ang23.sin() * tor_ang.cos();
    let p4z = d3 * ang23.sin() * tor_ang.sin();
    let dx = p4x - p1x;
    let dy = p4y - p1y;
    let dz = p4z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn compute_14_dist_cis(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() - d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

fn compute_14_dist_trans(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() + d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

fn is_larger_sp2_atom(
    mol: &Molecule,
    hybridizations: &[RdkitHybridization],
    atom_rings: &[Vec<usize>],
    atom_index: usize,
) -> bool {
    mol.atoms[atom_index].atomic_num > 13
        && hybridizations[atom_index] == RdkitHybridization::Sp2
        && atom_rings.iter().any(|ring| ring.contains(&atom_index))
}

fn set13_bounds_helper(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &ComputedData,
    hybridizations: &[RdkitHybridization],
    atom_rings: &[Vec<usize>],
    aid1: usize,
    aid2: usize,
    aid3: usize,
    bid1: usize,
    bid2: usize,
    angle: f64,
) {
    let dl = compute_13_dist(
        accum_data.bond_lengths[bid1],
        accum_data.bond_lengths[bid2],
        angle,
    );
    let mut dist_tol = DIST13_TOL;
    if is_larger_sp2_atom(mol, hybridizations, atom_rings, aid1) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom(mol, hybridizations, atom_rings, aid2) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom(mol, hybridizations, atom_rings, aid3) {
        dist_tol *= 2.0;
    }
    check_and_set_bounds(aid1, aid3, dl - dist_tol, dl + dist_tol, mmat);
}

fn set13_bounds_nonring(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
) -> Result<(), DgBoundsError> {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).map_err(|err| {
        DgBoundsError::Unsupported(format!("set13Bounds valence assignment failed: {err}"))
    })?;
    let mut atom_degree = vec![0usize; mol.atoms.len()];
    for bond in &mol.bonds {
        atom_degree[bond.begin_atom] += 1;
        atom_degree[bond.end_atom] += 1;
    }
    let conjugated = compute_conjugated_bonds(mol, &assignment, &atom_degree);
    let mut atom_has_conjugated_bond = vec![false; mol.atoms.len()];
    for (bi, bond) in mol.bonds.iter().enumerate() {
        if conjugated[bi] {
            atom_has_conjugated_bond[bond.begin_atom] = true;
            atom_has_conjugated_bond[bond.end_atom] = true;
        }
    }
    let hybridizations =
        compute_hybridizations(mol, &assignment, &atom_degree, &atom_has_conjugated_bond);
    let atom_rings = rdkit_atom_rings(mol);

    let nb = mol.bonds.len();
    let mut visited = vec![0usize; mol.atoms.len()];
    let mut angle_taken = vec![0.0f64; mol.atoms.len()];
    let mut done_paths = Vec::<(usize, usize)>::new();
    for ring in &atom_rings {
        let r_size = ring.len();
        let mut aid1 = ring[r_size - 1];
        for i in 0..r_size {
            let aid2 = ring[i];
            let aid3 = if i == r_size - 1 {
                ring[0]
            } else {
                ring[i + 1]
            };
            let bid1 = mol
                .bonds
                .iter()
                .find(|b| {
                    (b.begin_atom == aid1 && b.end_atom == aid2)
                        || (b.begin_atom == aid2 && b.end_atom == aid1)
                })
                .map(|b| b.index)
                .ok_or_else(|| {
                    DgBoundsError::Unsupported("set13Bounds ring edge lookup failed".to_string())
                })?;
            let bid2 = mol
                .bonds
                .iter()
                .find(|b| {
                    (b.begin_atom == aid2 && b.end_atom == aid3)
                        || (b.begin_atom == aid3 && b.end_atom == aid2)
                })
                .map(|b| b.index)
                .ok_or_else(|| {
                    DgBoundsError::Unsupported("set13Bounds ring edge lookup failed".to_string())
                })?;
            if done_paths.contains(&(bid1, bid2)) || done_paths.contains(&(bid2, bid1)) {
                aid1 = aid2;
                continue;
            }
            let angle = set_ring_angle(hybridizations[aid2], r_size);
            let pid = aid1.min(aid3) * mol.atoms.len() + aid1.max(aid3);
            if !accum_data.visited12_bounds[pid] {
                set13_bounds_helper(
                    mol,
                    mmat,
                    accum_data,
                    &hybridizations,
                    &atom_rings,
                    aid1,
                    aid2,
                    aid3,
                    bid1,
                    bid2,
                    angle,
                );
                accum_data.visited13_bounds[pid] = true;
            }
            accum_data.set_bond_angle(nb, bid1, bid2, angle);
            accum_data.set_bond_adj(nb, bid1, bid2, aid2 as i32);
            visited[aid2] += 1;
            angle_taken[aid2] += angle;
            done_paths.push((bid1, bid2));
            aid1 = aid2;
        }
    }

    for aid2 in 0..mol.atoms.len() {
        let neighbors: Vec<(usize, usize)> = mol
            .bonds
            .iter()
            .filter_map(|bond| {
                if bond.begin_atom == aid2 {
                    Some((bond.index, bond.end_atom))
                } else if bond.end_atom == aid2 {
                    Some((bond.index, bond.begin_atom))
                } else {
                    None
                }
            })
            .collect();
        for i in 0..neighbors.len() {
            for j in 0..i {
                let (bid1, aid1) = neighbors[i];
                let (bid2, aid3) = neighbors[j];
                let pid = aid1.min(aid3) * mol.atoms.len() + aid1.max(aid3);
                if accum_data.visited12_bounds[pid]
                    || accum_data.get_bond_angle(nb, bid1, bid2) > 0.0
                {
                    continue;
                }
                let deg = neighbors.len();
                let angle = if visited[aid2] >= 1 {
                    match hybridizations[aid2] {
                        RdkitHybridization::Sp2 => {
                            (2.0 * std::f64::consts::PI - angle_taken[aid2])
                                / ((deg * (deg - 1) / 2 - visited[aid2]) as f64)
                        }
                        RdkitHybridization::Sp3 => {
                            if is_atom_in_ring_of_size(&atom_rings, aid2, 3) {
                                116.0_f64.to_radians()
                            } else if is_atom_in_ring_of_size(&atom_rings, aid2, 4) {
                                112.0_f64.to_radians()
                            } else {
                                109.5_f64.to_radians()
                            }
                        }
                        _ => {
                            if deg == 5 {
                                105.0_f64.to_radians()
                            } else if deg == 6 {
                                135.0_f64.to_radians()
                            } else {
                                120.0_f64.to_radians()
                            }
                        }
                    }
                } else {
                    match hybridizations[aid2] {
                        RdkitHybridization::Sp => std::f64::consts::PI,
                        RdkitHybridization::Sp2 => 2.0 * std::f64::consts::PI / 3.0,
                        RdkitHybridization::Sp3 => 109.5_f64.to_radians(),
                        RdkitHybridization::Sp3d => 105.0_f64.to_radians(),
                        RdkitHybridization::Sp3d2 => 135.0_f64.to_radians(),
                        RdkitHybridization::S | RdkitHybridization::Unspecified => {
                            120.0_f64.to_radians()
                        }
                    }
                };
                set13_bounds_helper(
                    mol,
                    mmat,
                    accum_data,
                    &hybridizations,
                    &atom_rings,
                    aid1,
                    aid2,
                    aid3,
                    bid1,
                    bid2,
                    angle,
                );
                accum_data.visited13_bounds[pid] = true;
                accum_data.set_bond_angle(nb, bid1, bid2, angle);
                accum_data.set_bond_adj(nb, bid1, bid2, aid2 as i32);
                visited[aid2] += 1;
                angle_taken[aid2] += angle;
            }
        }
    }
    Ok(())
}

fn get_atom_stereo(bond: &crate::Bond, aid1: usize, aid4: usize) -> crate::BondStereo {
    let mut stype = bond.stereo;
    if bond.stereo_atoms.len() >= 2
        && ((bond.stereo_atoms[0] != aid1) ^ (bond.stereo_atoms[1] != aid4))
    {
        stype = match stype {
            crate::BondStereo::Cis => crate::BondStereo::Trans,
            crate::BondStereo::Trans => crate::BondStereo::Cis,
            other => other,
        };
    }
    stype
}

fn record_path_flag(path_store: &mut Vec<u64>, path_id: u64) {
    if !path_store.contains(&path_id) {
        path_store.push(path_id);
    }
}

fn has_path_flag(path_store: &[u64], path_id: u64) -> bool {
    path_store.contains(&path_id)
}

fn record_14_path(
    nb: usize,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    hybridizations: &[RdkitHybridization],
    accum_data: &mut ComputedData,
) {
    let aid2 = accum_data.get_bond_adj(nb, bid1, bid2) as usize;
    let aid3 = accum_data.get_bond_adj(nb, bid2, bid3) as usize;
    let mut path14 = Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind: Path14Type::Other,
    };
    if hybridizations[aid2] == RdkitHybridization::Sp2
        && hybridizations[aid3] == RdkitHybridization::Sp2
    {
        path14.kind = Path14Type::Cis;
        let path_id = bid1 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid3 as u64;
        record_path_flag(&mut accum_data.cis_paths, path_id);
        record_path_flag(
            &mut accum_data.cis_paths,
            bid3 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid1 as u64,
        );
    }
    accum_data.paths14.push(path14);
}

fn set_in_ring_14_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    hybridizations: &[RdkitHybridization],
    bond_ring_counts: &[usize],
    bond_rings: &[Vec<usize>],
    dmat: &[f64],
    bid1: usize,
    bid2: usize,
    bid3: usize,
    ring_size: usize,
) {
    let nb = mol.bonds.len();
    let aid2 = accum_data.get_bond_adj(nb, bid1, bid2) as usize;
    let aid3 = accum_data.get_bond_adj(nb, bid2, bid3) as usize;
    let aid1 = if mol.bonds[bid1].begin_atom == aid2 {
        mol.bonds[bid1].end_atom
    } else {
        mol.bonds[bid1].begin_atom
    };
    let aid4 = if mol.bonds[bid3].begin_atom == aid3 {
        mol.bonds[bid3].end_atom
    } else {
        mol.bonds[bid3].begin_atom
    };
    let pid = aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4);
    if accum_data.visited_bound(pid, false) {
        return;
    }
    if dmat[aid1.max(aid4) * mol.atoms.len() + aid1.min(aid4)] < 2.9 {
        return;
    }

    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(nb, bid1, bid2);
    let ba23 = accum_data.get_bond_angle(nb, bid2, bid3);

    let mut path14 = Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind: Path14Type::Other,
    };
    let stype = get_atom_stereo(&mol.bonds[bid2], aid1, aid4);
    let mut prefer_cis = false;
    let mut prefer_trans = false;
    if ring_size <= 8
        && hybridizations[aid2] == RdkitHybridization::Sp2
        && hybridizations[aid3] == RdkitHybridization::Sp2
        && !matches!(stype, crate::BondStereo::Trans)
    {
        if bond_ring_counts[bid2] > 1 {
            if bond_ring_counts[bid1] == 1 && bond_ring_counts[bid3] == 1 {
                for br in bond_rings {
                    if br.contains(&bid1) {
                        if br.contains(&bid3) {
                            prefer_cis = true;
                        }
                        break;
                    }
                }
            }
        } else {
            prefer_cis = true;
        }
    } else if matches!(stype, crate::BondStereo::Cis) {
        prefer_cis = true;
    } else if matches!(stype, crate::BondStereo::Trans) {
        prefer_trans = true;
    }

    let (dl, du) = if prefer_cis {
        path14.kind = Path14Type::Cis;
        let path_id = bid1 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid3 as u64;
        record_path_flag(&mut accum_data.cis_paths, path_id);
        record_path_flag(
            &mut accum_data.cis_paths,
            bid3 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid1 as u64,
        );
        let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
        (dl, dl + 2.0 * GEN_DIST_TOL)
    } else if prefer_trans {
        path14.kind = Path14Type::Trans;
        let path_id = bid1 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid3 as u64;
        record_path_flag(&mut accum_data.trans_paths, path_id);
        record_path_flag(
            &mut accum_data.trans_paths,
            bid3 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid1 as u64,
        );
        let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
        (dl, dl + 2.0 * GEN_DIST_TOL)
    } else {
        let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if du < dl {
            std::mem::swap(&mut dl, &mut du);
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        (dl, du)
    };
    accum_data.paths14.push(path14);
    accum_data.visited14_bounds[pid] = true;
    check_and_set_bounds(aid1, aid4, dl, du, mmat);
}

fn set_two_in_same_ring_14_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    hybridizations: &[RdkitHybridization],
    dmat: &[f64],
    bid1: usize,
    bid2: usize,
    bid3: usize,
) {
    let nb = mol.bonds.len();
    let aid2 = accum_data.get_bond_adj(nb, bid1, bid2) as usize;
    let aid3 = accum_data.get_bond_adj(nb, bid2, bid3) as usize;
    let aid1 = if mol.bonds[bid1].begin_atom == aid2 {
        mol.bonds[bid1].end_atom
    } else {
        mol.bonds[bid1].begin_atom
    };
    let aid4 = if mol.bonds[bid3].begin_atom == aid3 {
        mol.bonds[bid3].end_atom
    } else {
        mol.bonds[bid3].begin_atom
    };
    let pid = aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4);
    if accum_data.visited_bound(pid, false) {
        return;
    }
    if dmat[aid1.max(aid4) * mol.atoms.len() + aid1.min(aid4)] < 2.9 {
        return;
    }
    if mol.bonds.iter().any(|b| {
        (b.begin_atom == aid1 && b.end_atom == aid3)
            || (b.end_atom == aid1 && b.begin_atom == aid3)
            || (b.begin_atom == aid4 && b.end_atom == aid2)
            || (b.end_atom == aid4 && b.begin_atom == aid2)
    }) {
        return;
    }

    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(nb, bid1, bid2);
    let ba23 = accum_data.get_bond_angle(nb, bid2, bid3);

    let mut path14 = Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind: Path14Type::Other,
    };
    let (dl, du) = if hybridizations[aid2] == RdkitHybridization::Sp2
        && hybridizations[aid3] == RdkitHybridization::Sp2
    {
        let du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        path14.kind = Path14Type::Trans;
        let path_id = bid1 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid3 as u64;
        record_path_flag(&mut accum_data.trans_paths, path_id);
        record_path_flag(
            &mut accum_data.trans_paths,
            bid3 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid1 as u64,
        );
        (du - GEN_DIST_TOL, du + GEN_DIST_TOL)
    } else {
        let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if du < dl {
            std::mem::swap(&mut dl, &mut du);
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        (dl, du)
    };
    check_and_set_bounds(aid1, aid4, dl, du, mmat);
    accum_data.paths14.push(path14);
    accum_data.visited14_bounds[pid] = true;
}

fn is_carbonyl(mol: &Molecule, atom_index: usize) -> bool {
    let atom = &mol.atoms[atom_index];
    if atom.atomic_num != 6 {
        return false;
    }
    let degree = mol
        .bonds
        .iter()
        .filter(|b| b.begin_atom == atom_index || b.end_atom == atom_index)
        .count();
    if degree <= 2 {
        return false;
    }
    mol.bonds.iter().any(|bond| {
        let nbr = if bond.begin_atom == atom_index {
            Some(bond.end_atom)
        } else if bond.end_atom == atom_index {
            Some(bond.begin_atom)
        } else {
            None
        };
        nbr.is_some_and(|nbr_idx| {
            let z = mol.atoms[nbr_idx].atomic_num;
            (z == 8 || z == 7) && matches!(bond.order, BondOrder::Double)
        })
    })
}

fn check_amide_ester14(
    mol: &Molecule,
    bnd1: &crate::Bond,
    bnd3: &crate::Bond,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> bool {
    let a2_num = mol.atoms[atm2].atomic_num;
    let a3_num = mol.atoms[atm3].atomic_num;
    let a4_num = mol.atoms[atm4].atomic_num;
    a3_num == 6
        && matches!(bnd3.order, BondOrder::Double)
        && (a4_num == 8 || a4_num == 7)
        && matches!(bnd1.order, BondOrder::Single)
        && (a2_num == 8 || (a2_num == 7 && atom_total_num_hs(mol, atm2) == 1))
}

fn check_amide_ester15(
    mol: &Molecule,
    bnd1: &crate::Bond,
    bnd3: &crate::Bond,
    atm2: usize,
    atm3: usize,
) -> bool {
    let a2_num = mol.atoms[atm2].atomic_num;
    (a2_num == 8 || (a2_num == 7 && atom_total_num_hs(mol, atm2) == 1))
        && matches!(bnd1.order, BondOrder::Single)
        && mol.atoms[atm3].atomic_num == 6
        && matches!(bnd3.order, BondOrder::Single)
        && is_carbonyl(mol, atm3)
}

fn compute_15_dists_cis_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    compute_13_dist(d14, d4, ang34 - ang143)
}

fn compute_15_dists_cis_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    compute_13_dist(d14, d4, ang34 + ang143)
}

fn compute_15_dists_trans_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    compute_13_dist(d14, d4, ang34 + ang143)
}

fn compute_15_dists_trans_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    compute_13_dist(d14, d4, ang34 - ang143)
}

fn set_chain14_bounds_nonring(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
) -> Result<(), DgBoundsError> {
    let nb = mol.bonds.len();
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).map_err(|err| {
        DgBoundsError::Unsupported(format!("set14Bounds valence assignment failed: {err}"))
    })?;
    let mut atom_degree = vec![0usize; mol.atoms.len()];
    for bond in &mol.bonds {
        atom_degree[bond.begin_atom] += 1;
        atom_degree[bond.end_atom] += 1;
    }
    let conjugated = compute_conjugated_bonds(mol, &assignment, &atom_degree);
    let mut atom_has_conjugated_bond = vec![false; mol.atoms.len()];
    for (bi, bond) in mol.bonds.iter().enumerate() {
        if conjugated[bi] {
            atom_has_conjugated_bond[bond.begin_atom] = true;
            atom_has_conjugated_bond[bond.end_atom] = true;
        }
    }
    let hybridizations =
        compute_hybridizations(mol, &assignment, &atom_degree, &atom_has_conjugated_bond);
    let bond_rings = rdkit_bond_rings(mol);
    let mut ring_bond_pairs = BTreeSet::<u64>::new();
    let mut done_paths = BTreeSet::<u64>::new();
    let mut bond_ring_counts = vec![0usize; mol.bonds.len()];
    for ring in &bond_rings {
        for &bid in ring {
            bond_ring_counts[bid] += 1;
        }
    }
    for bring in &bond_rings {
        let r_size = bring.len();
        if r_size < 3 {
            continue;
        }
        let mut bid1 = bring[r_size - 1];
        for i in 0..r_size {
            let bid2 = bring[i];
            let bid3 = bring[(i + 1) % r_size];
            ring_bond_pairs.insert(bid1 as u64 * nb as u64 + bid2 as u64);
            ring_bond_pairs.insert(bid2 as u64 * nb as u64 + bid1 as u64);
            done_paths.insert(
                bid1 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid3 as u64,
            );
            done_paths.insert(
                bid3 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid1 as u64,
            );
            if r_size > 5 {
                set_in_ring_14_bounds(
                    mol,
                    mmat,
                    accum_data,
                    &hybridizations,
                    &bond_ring_counts,
                    &bond_rings,
                    &graph_distance_matrix(mol),
                    bid1,
                    bid2,
                    bid3,
                    r_size,
                );
            } else {
                record_14_path(nb, bid1, bid2, bid3, &hybridizations, accum_data);
            }
            bid1 = bid2;
        }
    }
    let dmat = graph_distance_matrix(mol);
    for bnd2 in &mol.bonds {
        let bid2 = bnd2.index;
        let aid2 = bnd2.begin_atom;
        let aid3 = bnd2.end_atom;
        let left_bonds: Vec<&crate::Bond> = mol
            .bonds
            .iter()
            .filter(|bond| bond.index != bid2 && (bond.begin_atom == aid2 || bond.end_atom == aid2))
            .collect();
        let right_bonds: Vec<&crate::Bond> = mol
            .bonds
            .iter()
            .filter(|bond| bond.index != bid2 && (bond.begin_atom == aid3 || bond.end_atom == aid3))
            .collect();
        for bnd1 in &left_bonds {
            for bnd3 in &right_bonds {
                let bid1 = bnd1.index;
                let bid3 = bnd3.index;
                let id1 =
                    bid1 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid3 as u64;
                let id2 =
                    bid3 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid1 as u64;
                if done_paths.contains(&id1) || done_paths.contains(&id2) {
                    continue;
                }
                let aid1 = if bnd1.begin_atom == aid2 {
                    bnd1.end_atom
                } else {
                    bnd1.begin_atom
                };
                let aid4 = if bnd3.begin_atom == aid3 {
                    bnd3.end_atom
                } else {
                    bnd3.begin_atom
                };
                if aid1 == aid4 {
                    continue;
                }
                let pid = aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4);
                if accum_data.visited_bound(pid, false) {
                    continue;
                }
                let pid1 = bid1 as u64 * nb as u64 + bid2 as u64;
                let pid2 = bid2 as u64 * nb as u64 + bid1 as u64;
                let pid3 = bid2 as u64 * nb as u64 + bid3 as u64;
                let pid4 = bid3 as u64 * nb as u64 + bid2 as u64;
                if ring_bond_pairs.contains(&pid1)
                    || ring_bond_pairs.contains(&pid2)
                    || ring_bond_pairs.contains(&pid3)
                    || ring_bond_pairs.contains(&pid4)
                {
                    set_two_in_same_ring_14_bounds(
                        mol,
                        mmat,
                        accum_data,
                        &hybridizations,
                        &dmat,
                        bid1,
                        bid2,
                        bid3,
                    );
                    continue;
                }
                if (bond_ring_counts[bid1] > 0 && bond_ring_counts[bid2] > 0)
                    || (bond_ring_counts[bid2] > 0 && bond_ring_counts[bid3] > 0)
                    || bond_ring_counts[bid2] > 0
                {
                    set_in_ring_14_bounds(
                        mol,
                        mmat,
                        accum_data,
                        &hybridizations,
                        &bond_ring_counts,
                        &bond_rings,
                        &dmat,
                        bid1,
                        bid2,
                        bid3,
                        0,
                    );
                    continue;
                }
                let bl1 = accum_data.bond_lengths[bid1];
                let bl2 = accum_data.bond_lengths[bid2];
                let bl3 = accum_data.bond_lengths[bid3];
                let ba12 = accum_data.get_bond_angle(nb, bid1, bid2);
                let ba23 = accum_data.get_bond_angle(nb, bid2, bid3);
                if !(ba12 > 0.0 && ba23 > 0.0) {
                    return Err(DgBoundsError::Unsupported(
                        "set14Bounds missing prerequisite bond angles from set13Bounds".to_string(),
                    ));
                }
                let mut path14 = Path14Configuration {
                    bid1,
                    bid2,
                    bid3,
                    kind: Path14Type::Other,
                };
                let (mut dl, mut du) = match bnd2.order {
                    BondOrder::Double => {
                        if matches!(bnd1.order, BondOrder::Double)
                            || matches!(bnd3.order, BondOrder::Double)
                        {
                            let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                            let du = dl + 2.0 * GEN_DIST_TOL;
                            path14.kind = Path14Type::Cis;
                            record_path_flag(&mut accum_data.cis_paths, id1);
                            record_path_flag(&mut accum_data.cis_paths, id2);
                            (dl, du)
                        } else if matches!(
                            bnd2.stereo,
                            crate::BondStereo::Cis | crate::BondStereo::Trans
                        ) {
                            match get_atom_stereo(bnd2, aid1, aid4) {
                                crate::BondStereo::Cis => {
                                    let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23)
                                        - GEN_DIST_TOL;
                                    let du = dl + 2.0 * GEN_DIST_TOL;
                                    path14.kind = Path14Type::Cis;
                                    record_path_flag(&mut accum_data.cis_paths, id1);
                                    record_path_flag(&mut accum_data.cis_paths, id2);
                                    (dl, du)
                                }
                                crate::BondStereo::Trans => {
                                    let du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                                    let dl = du - GEN_DIST_TOL;
                                    let du = du + GEN_DIST_TOL;
                                    path14.kind = Path14Type::Trans;
                                    record_path_flag(&mut accum_data.trans_paths, id1);
                                    record_path_flag(&mut accum_data.trans_paths, id2);
                                    (dl, du)
                                }
                                crate::BondStereo::Any | crate::BondStereo::None => unreachable!(),
                            }
                        } else {
                            (
                                compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                                compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                            )
                        }
                    }
                    BondOrder::Single => {
                        let atom2 = &mol.atoms[aid2];
                        let atom3 = &mol.atoms[aid3];
                        let atom1 = &mol.atoms[aid1];
                        let atom4 = &mol.atoms[aid4];
                        let deg2 = mol
                            .bonds
                            .iter()
                            .filter(|b| b.begin_atom == aid2 || b.end_atom == aid2)
                            .count();
                        let deg3 = mol
                            .bonds
                            .iter()
                            .filter(|b| b.begin_atom == aid3 || b.end_atom == aid3)
                            .count();
                        if atom2.atomic_num == 16
                            && atom3.atomic_num == 16
                            && deg2 == 2
                            && deg3 == 2
                        {
                            let dl = compute_14_dist_3d(
                                bl1,
                                bl2,
                                bl3,
                                ba12,
                                ba23,
                                std::f64::consts::PI / 2.0,
                            ) - GEN_DIST_TOL;
                            (dl, dl + 2.0 * GEN_DIST_TOL)
                        } else if check_amide_ester14(mol, bnd1, bnd3, aid2, aid3, aid4)
                            || check_amide_ester14(mol, bnd3, bnd1, aid3, aid2, aid1)
                        {
                            let secondary_amide_h = (atom1.atomic_num == 1
                                && atom2.atomic_num == 7
                                && deg2 == 3
                                && atom_total_num_hs(mol, aid2) == 1)
                                || (atom4.atomic_num == 1
                                    && atom3.atomic_num == 7
                                    && deg3 == 3
                                    && atom_total_num_hs(mol, aid3) == 1);
                            if secondary_amide_h {
                                let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                                path14.kind = Path14Type::Trans;
                                record_path_flag(&mut accum_data.trans_paths, id1);
                                record_path_flag(&mut accum_data.trans_paths, id2);
                                (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL)
                            } else {
                                let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                                path14.kind = Path14Type::Cis;
                                record_path_flag(&mut accum_data.cis_paths, id1);
                                record_path_flag(&mut accum_data.cis_paths, id2);
                                (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL)
                            }
                        } else if check_amide_ester15(mol, bnd1, bnd3, aid2, aid3)
                            || check_amide_ester15(mol, bnd3, bnd1, aid3, aid2)
                        {
                            if atom2.atomic_num == 7
                                && deg2 == 3
                                && atom1.atomic_num == 1
                                && atom_total_num_hs(mol, aid2) == 1
                            {
                                continue;
                            }
                            let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                            path14.kind = Path14Type::Trans;
                            record_path_flag(&mut accum_data.trans_paths, id1);
                            record_path_flag(&mut accum_data.trans_paths, id2);
                            (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL)
                        } else {
                            (
                                compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                                compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                            )
                        }
                    }
                    _ => (
                        compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                        compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                    ),
                };
                if (du - dl).abs() < DIST12_DELTA {
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                }
                check_and_set_bounds(aid1, aid4, dl, du, mmat);
                accum_data.paths14.push(path14);
                accum_data.visited14_bounds[pid] = true;
            }
        }
    }
    Ok(())
}

fn set15_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
) {
    let nb = mol.bonds.len();
    let na = mol.atoms.len();
    let paths14 = accum_data.paths14.clone();
    for path in paths14 {
        for &(bid1, bid2, bid3, kind) in &[
            (path.bid1, path.bid2, path.bid3, path.kind),
            (path.bid3, path.bid2, path.bid1, path.kind),
        ] {
            let aid2 = accum_data.get_bond_adj(nb, bid1, bid2) as usize;
            let aid1 = if mol.bonds[bid1].begin_atom == aid2 {
                mol.bonds[bid1].end_atom
            } else {
                mol.bonds[bid1].begin_atom
            };
            let aid3 = accum_data.get_bond_adj(nb, bid2, bid3) as usize;
            let aid4 = if mol.bonds[bid3].begin_atom == aid3 {
                mol.bonds[bid3].end_atom
            } else {
                mol.bonds[bid3].begin_atom
            };
            let d1 = accum_data.bond_lengths[bid1];
            let d2 = accum_data.bond_lengths[bid2];
            let d3 = accum_data.bond_lengths[bid3];
            let ang12 = accum_data.get_bond_angle(nb, bid1, bid2);
            let ang23 = accum_data.get_bond_angle(nb, bid2, bid3);
            for i in 0..nb {
                if accum_data.get_bond_adj(nb, bid3, i) != aid4 as i32 {
                    continue;
                }
                let aid5 = if mol.bonds[i].begin_atom == aid4 {
                    mol.bonds[i].end_atom
                } else {
                    mol.bonds[i].begin_atom
                };
                let pid = aid1.min(aid5) * na + aid1.max(aid5);
                if accum_data.visited_bound(pid, true) {
                    continue;
                }
                if dmat[aid1.max(aid5) * na + aid1.min(aid5)] < 3.9 || aid1 == aid5 {
                    continue;
                }
                let pid1 = aid1 * na + aid5;
                let pid2 = aid5 * na + aid1;
                if !(mmat.get_lower_bound(aid1, aid5) < DIST12_DELTA
                    || accum_data.set15_atoms[pid1]
                    || accum_data.set15_atoms[pid2])
                {
                    continue;
                }
                let d4 = accum_data.bond_lengths[i];
                let ang34 = accum_data.get_bond_angle(nb, bid3, i);
                let path_id =
                    bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + i as u64;
                let (mut dl, mut du) = match kind {
                    Path14Type::Cis => {
                        if has_path_flag(&accum_data.cis_paths, path_id) {
                            let base =
                                compute_15_dists_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                            (base - 0.08, base + 0.08)
                        } else if has_path_flag(&accum_data.trans_paths, path_id) {
                            let base =
                                compute_15_dists_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                            (base - 0.08, base + 0.08)
                        } else {
                            (
                                compute_15_dists_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34)
                                    - 0.08,
                                compute_15_dists_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                                    + 0.08,
                            )
                        }
                    }
                    Path14Type::Trans => {
                        if has_path_flag(&accum_data.cis_paths, path_id) {
                            let base =
                                compute_15_dists_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                            (base - 0.08, base + 0.08)
                        } else if has_path_flag(&accum_data.trans_paths, path_id) {
                            let base =
                                compute_15_dists_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                            (base - 0.08, base + 0.08)
                        } else {
                            (
                                compute_15_dists_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34)
                                    - 0.08,
                                compute_15_dists_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                                    + 0.08,
                            )
                        }
                    }
                    Path14Type::Other => {
                        if has_path_flag(&accum_data.cis_paths, path_id) {
                            (
                                compute_15_dists_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12)
                                    - 0.08,
                                compute_15_dists_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                                    + 0.08,
                            )
                        } else if has_path_flag(&accum_data.trans_paths, path_id) {
                            (
                                compute_15_dists_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12)
                                    - 0.08,
                                compute_15_dists_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                                    + 0.08,
                            )
                        } else {
                            let vw1 = rvdw(mol.atoms[aid1].atomic_num);
                            let vw5 = rvdw(mol.atoms[aid5].atomic_num);
                            (VDW_SCALE_15 * (vw1 + vw5), MAX_UPPER)
                        }
                    }
                };
                if du < 0.0 {
                    du = MAX_UPPER;
                }
                if dl < 0.0 {
                    dl = 0.0;
                }
                let debug_pair = std::env::var("COSMOLKIT_DEBUG_DG_PAIR").ok().and_then(|s| {
                    let (a, b) = s.split_once(',')?;
                    Some((a.parse::<usize>().ok()?, b.parse::<usize>().ok()?))
                });
                if std::env::var_os("COSMOLKIT_DEBUG_DG").is_some()
                    && debug_pair
                        .map(|(a, b)| (aid1 == a && aid5 == b) || (aid1 == b && aid5 == a))
                        .unwrap_or(false)
                {
                    let dbg_ciscis = compute_15_dists_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                    let dbg_cistrans =
                        compute_15_dists_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                    let dbg_transcis =
                        compute_15_dists_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                    let dbg_transtrans =
                        compute_15_dists_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                    eprintln!(
                        "DEBUG set15 0-4 kind={:?} path=({}, {}, {}) next={} path_id={} cis?={} trans?={} cc={} ct={} tc={} tt={} dl={} du={}",
                        kind,
                        bid1,
                        bid2,
                        bid3,
                        i,
                        path_id,
                        has_path_flag(&accum_data.cis_paths, path_id),
                        has_path_flag(&accum_data.trans_paths, path_id),
                        dbg_ciscis,
                        dbg_cistrans,
                        dbg_transcis,
                        dbg_transtrans,
                        dl,
                        du
                    );
                }
                check_and_set_bounds(aid1, aid5, dl, du, mmat);
                accum_data.set15_atoms[pid1] = true;
                accum_data.set15_atoms[pid2] = true;
            }
        }
    }
}

fn is_hbond_acceptor(atom: &crate::Atom) -> bool {
    atom.atomic_num == 7 || atom.atomic_num == 8
}

fn atom_total_num_hs(mol: &Molecule, atom_index: usize) -> usize {
    mol.atoms[atom_index].explicit_hydrogens as usize
        + match assign_valence(mol, ValenceModel::RdkitLike) {
            Ok(v) => v.implicit_hydrogens[atom_index] as usize,
            Err(_) => 0,
        }
        + mol
            .bonds
            .iter()
            .filter(|b| {
                (b.begin_atom == atom_index && mol.atoms[b.end_atom].atomic_num == 1)
                    || (b.end_atom == atom_index && mol.atoms[b.begin_atom].atomic_num == 1)
            })
            .count()
}

fn is_h_in_hbond_donor(atom_index: usize, mol: &Molecule) -> bool {
    if mol.atoms[atom_index].atomic_num != 1 {
        return false;
    }
    mol.bonds.iter().any(|b| {
        let nbr = if b.begin_atom == atom_index {
            Some(b.end_atom)
        } else if b.end_atom == atom_index {
            Some(b.begin_atom)
        } else {
            None
        };
        nbr.is_some_and(|n| {
            let z = mol.atoms[n].atomic_num;
            z == 7 || z == 8
        })
    })
}

fn set_lower_bound_vdw(mol: &Molecule, mmat: &mut BoundsMatrix, dmat: &[f64]) {
    let n = mol.atoms.len();
    let mut hin_hbond_donors = vec![false; n];
    let mut hbond_acceptors = vec![false; n];
    for i in 1..n {
        let atom_i = &mol.atoms[i];
        if is_h_in_hbond_donor(i, mol) {
            hin_hbond_donors[i] = true;
        }
        if is_hbond_acceptor(atom_i) {
            hbond_acceptors[i] = true;
        }
        let vw1 = rvdw(atom_i.atomic_num);
        for j in 0..i {
            let atom_j = &mol.atoms[j];
            let vw2 = rvdw(atom_j.atomic_num);
            if mmat.get_lower_bound(i, j) < DIST12_DELTA {
                if (hin_hbond_donors[i] && hbond_acceptors[j])
                    || (hbond_acceptors[i] && hin_hbond_donors[j])
                {
                    mmat.set_lower_bound(i, j, H_BOND_LENGTH);
                } else {
                    let dij = dmat[i * n + j];
                    if dij == 4.0 {
                        mmat.set_lower_bound(i, j, VDW_SCALE_15 * (vw1 + vw2));
                    } else if dij == 5.0 {
                        mmat.set_lower_bound(
                            i,
                            j,
                            (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2),
                        );
                    } else {
                        mmat.set_lower_bound(i, j, vw1 + vw2);
                    }
                }
            }
        }
    }
}

pub fn dg_bounds_matrix(mol: &Molecule) -> Result<Vec<Vec<f64>>, DgBoundsError> {
    let na = mol.atoms.len();
    let nb = mol.bonds.len();
    if na == 0 {
        return Err(DgBoundsError::EmptyMolecule);
    }
    let max_num_bonds = (u64::MAX as f64).powf(1.0 / 3.0) as usize;
    if nb >= max_num_bonds {
        return Err(DgBoundsError::TooManyBonds);
    }

    let mut mmat = BoundsMatrix::new(na);
    init_bounds_mat(&mut mmat, 0.0, MAX_UPPER);
    let dmat = graph_distance_matrix(mol);
    let mut accum_data = ComputedData::new(na, nb);
    set12_bounds(mol, &mut mmat, &mut accum_data)?;

    if na > 2 {
        set13_bounds_nonring(mol, &mut mmat, &mut accum_data)?;
    }

    let needs_14 = dmat.iter().any(|&d| d > 2.0 && d.is_finite());
    if needs_14 {
        set_chain14_bounds_nonring(mol, &mut mmat, &mut accum_data)?;
    }

    let needs_15 = dmat.iter().any(|&d| d > 3.0 && d.is_finite());
    if needs_15 {
        set15_bounds(mol, &mut mmat, &mut accum_data, &dmat);
    }

    set_lower_bound_vdw(mol, &mut mmat, &dmat);
    if !triangle_smooth_bounds(&mut mmat, 0.0) {
        return Err(DgBoundsError::Unsupported(
            "triangle smoothing failed for current DG bounds subset".to_string(),
        ));
    }
    Ok(mmat.into_rows())
}
