use core::fmt;
use std::f64::consts::PI;

use crate::{BondOrder, BondStereo, Molecule};

#[derive(Debug)]
pub enum MolWriteError {
    UnsupportedSubset(&'static str),
}

impl fmt::Display for MolWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSubset(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for MolWriteError {}

pub fn mol_to_v2000_block_minimal(mol: &Molecule) -> Result<String, MolWriteError> {
    // Mirror RDKit MolToMolBlock default behavior for V2000 output:
    // write through a kekulized copy when aromatic bonds are present.
    // (Compute2DCoords is done on the original molecule state.)
    let mut write_mol = mol.clone();
    if write_mol
        .bonds
        .iter()
        .any(|b| matches!(b.order, BondOrder::Aromatic))
    {
        crate::kekulize::kekulize_in_place(&mut write_mol)
            .map_err(|_| MolWriteError::UnsupportedSubset("kekulize before v2000 write failed"))?;
    }

    let coords = compute_2d_coords_minimal(mol)?;

    let mut out = String::new();
    out.push('\n');
    out.push_str("     COSMolKit      2D\n");
    out.push('\n');
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        write_mol.atoms.len(),
        write_mol.bonds.len()
    ));

    for (idx, atom) in write_mol.atoms.iter().enumerate() {
        let (x, y) = coords[idx];
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            x,
            y,
            0.0,
            atom_symbol(atom.atomic_num)
        ));
    }

    for bond in &write_mol.bonds {
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}  0\n",
            bond.begin_atom + 1,
            bond.end_atom + 1,
            bond_type_code(bond.order)
        ));
    }
    out.push_str("M  END\n");
    Ok(out)
}

pub fn mol_to_sdf_record_minimal(mol: &Molecule) -> Result<String, MolWriteError> {
    let mut block = mol_to_v2000_block_minimal(mol)?;
    block.push_str("$$$$\n");
    Ok(block)
}

fn compute_2d_coords_minimal(mol: &Molecule) -> Result<Vec<(f64, f64)>, MolWriteError> {
    let n = mol.atoms.len();
    if n == 0 {
        return Err(MolWriteError::UnsupportedSubset(
            "empty molecule is unsupported",
        ));
    }
    if let Some(coords) = strict_rdkit_subset_coords(mol) {
        return Ok(coords);
    }
    compute_2d_coords_rdkit_like(mol)
}

// Source-aligned entrypoint for 2D depiction.
// Target mapping (RDKit):
//   compute2DCoords()
//     -> computeInitialCoords()
//     -> removeCollisionsBondFlip()
//     -> removeCollisionsOpenAngles()
//     -> removeCollisionsShortenBonds()
//     -> canonicalizeOrientation()
//     -> _shiftCoords()
//
// Stages are explicit so we can replace internals with one-to-one source
// replication without changing call sites.
fn compute_2d_coords_rdkit_like(mol: &Molecule) -> Result<Vec<(f64, f64)>, MolWriteError> {
    let staged = rdkit_stage_compute_initial_coords_subset(mol)?;
    let staged = rdkit_stage_cleanup_subset(mol, staged)?;
    let staged = rdkit_stage_orient_and_shift_subset(staged)?;
    Ok(staged)
}

// Stage 1 (interim): delegated to the existing implementation.
// This is temporary while we port computeInitialCoords()/EmbeddedFrag flow.
fn rdkit_stage_compute_initial_coords_subset(
    mol: &Molecule,
) -> Result<Vec<(f64, f64)>, MolWriteError> {
    rdkit_compute_initial_coords_strict(mol)
}

// Stage 2 (interim no-op): collision cleanup currently happens inside the
// legacy function and will be split out during source-level porting.
fn rdkit_stage_cleanup_subset(
    _mol: &Molecule,
    coords: Vec<(f64, f64)>,
) -> Result<Vec<(f64, f64)>, MolWriteError> {
    Ok(coords)
}

// Stage 3 (interim no-op): canonical orientation and component shifting also
// currently happen inside the legacy function.
fn rdkit_stage_orient_and_shift_subset(
    coords: Vec<(f64, f64)>,
) -> Result<Vec<(f64, f64)>, MolWriteError> {
    Ok(coords)
}

fn rdkit_compute_initial_coords_strict(mol: &Molecule) -> Result<Vec<(f64, f64)>, MolWriteError> {
    let n = mol.atoms.len();
    let mut adjacency = vec![Vec::<usize>::new(); n];
    let mut degree = vec![0usize; n];
    for b in &mol.bonds {
        adjacency[b.begin_atom].push(b.end_atom);
        adjacency[b.end_atom].push(b.begin_atom);
        degree[b.begin_atom] += 1;
        degree[b.end_atom] += 1;
    }
    for nbs in &mut adjacency {
        nbs.sort_unstable();
        nbs.dedup();
    }
    let hybridizations = rdkit_hybridizations_for_depict(mol, &degree)?;
    let cip_ranks = rdkit_cip_ranks_for_depict(mol);

    let mut visited = vec![false; n];
    let mut components: Vec<Vec<usize>> = Vec::new();
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut stack = vec![start];
        visited[start] = true;
        let mut comp = Vec::new();
        while let Some(v) = stack.pop() {
            comp.push(v);
            for &nb in &adjacency[v] {
                if !visited[nb] {
                    visited[nb] = true;
                    stack.push(nb);
                }
            }
        }
        comp.sort_unstable();
        components.push(comp);
    }
    components.sort_by_key(|c| c[0]);

    let mut local_components: Vec<Vec<(usize, (f64, f64))>> = Vec::new();

    for comp in components {
        let k = comp.len();
        let local: Option<Vec<(usize, (f64, f64))>> = if k == 1 {
            Some(vec![(comp[0], (0.0, 0.0))])
        } else if k == 2 {
            let bond_order = mol
                .bonds
                .iter()
                .find(|b| {
                    (b.begin_atom == comp[0] && b.end_atom == comp[1])
                        || (b.begin_atom == comp[1] && b.end_atom == comp[0])
                })
                .map(|b| b.order);
            if matches!(bond_order, Some(BondOrder::Null)) {
                Some(vec![(comp[0], (0.7500, 0.0)), (comp[1], (-0.7500, 0.0))])
            } else {
                Some(vec![(comp[0], (-0.7500, 0.0)), (comp[1], (0.7500, 0.0))])
            }
        } else if k == 3 {
            // Keep explicit support for the common linear dative triad:
            // [donor]->[metal]<-[donor]
            // This pattern is present in the regression corpus and RDKit
            // places it linearly in 2D.
            let dative_bonds: Vec<_> = mol
                .bonds
                .iter()
                .filter(|b| {
                    matches!(b.order, BondOrder::Dative)
                        && comp.contains(&b.begin_atom)
                        && comp.contains(&b.end_atom)
                })
                .collect();
            if dative_bonds.len() == 2 {
                let center = comp.iter().copied().find(|&a| degree[a] == 2);
                if let Some(c) = center {
                    let mut ends: Vec<usize> = comp.iter().copied().filter(|&x| x != c).collect();
                    ends.sort_unstable();
                    Some(vec![
                        (ends[0], (-1.5000, 0.0)),
                        (c, (0.0, 0.0)),
                        (ends[1], (1.5000, 0.0)),
                    ])
                } else {
                    None
                }
            } else {
                place_acyclic_tree_rdkit_like(
                    mol,
                    &comp,
                    &adjacency,
                    &degree,
                    &hybridizations,
                    &cip_ranks,
                    None,
                )
            }
        } else {
            let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();
            let bond_count_in_comp = mol
                .bonds
                .iter()
                .filter(|b| comp_set.contains(&b.begin_atom) && comp_set.contains(&b.end_atom))
                .count();
            let cyclomatic = bond_count_in_comp as isize - k as isize + 1;

            let is_tree = bond_count_in_comp + 1 == k;
            let branch_try = if is_tree {
                let centers: Vec<_> = comp.iter().copied().filter(|&a| degree[a] > 1).collect();
                if centers.len() == 1 {
                    let center = centers[0];
                    let mut leaves: Vec<_> =
                        comp.iter().copied().filter(|&a| a != center).collect();
                    if degree[center] == 4
                        && leaves.len() == 4
                        && matches!(mol.atoms[center].chiral_tag, crate::ChiralTag::Unspecified)
                        && leaves.iter().all(|&a| degree[a] == 1)
                    {
                        leaves.sort_by_key(|&idx| {
                            let depict = atom_depict_rank(mol.atoms[idx].atomic_num, degree[idx]);
                            depict * mol.atoms.len() + idx
                        });
                        Some(vec![
                            (center, (0.0, 0.0)),
                            (leaves[0], (-1.2990, -0.7500)),
                            (leaves[1], (1.2990, 0.7500)),
                            (leaves[2], (0.7500, -1.2990)),
                            (leaves[3], (-0.7500, 1.2990)),
                        ])
                    } else {
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            };
            let chosen = branch_try.or_else(|| {
                if cyclomatic <= 0 {
                    place_acyclic_tree_rdkit_like(
                        mol,
                        &comp,
                        &adjacency,
                        &degree,
                        &hybridizations,
                        &cip_ranks,
                        None,
                    )
                } else if cyclomatic == 1 {
                    place_multiring_nonfused_component_rdkit_like(mol, &comp, &adjacency, &degree)
                } else {
                    place_multiring_nonfused_component_rdkit_like(mol, &comp, &adjacency, &degree)
                }
            });
            chosen
        };

        let mut local = local.ok_or(MolWriteError::UnsupportedSubset(
            "strict RDKit computeInitialCoords branch missing for this component; no heuristic fallback",
        ))?;
        if std::env::var_os("COSMOLKIT_DEBUG_PRECANON").is_some() {
            eprintln!("DEBUG precanon local: {:?}", local);
        }
        canonicalize_component_rdkit_like(&mut local);
        local_components.push(local);
    }

    let mut out = vec![(0.0f64, 0.0f64); n];
    if local_components.is_empty() {
        return Ok(out);
    }
    let box_of = |component: &Vec<(usize, (f64, f64))>| -> (f64, f64, f64, f64) {
        let mut max_x = f64::NEG_INFINITY;
        let mut min_x = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        for &(_, (x, y)) in component {
            max_x = max_x.max(x);
            min_x = min_x.min(x);
            max_y = max_y.max(y);
            min_y = min_y.min(y);
        }
        (
            max_x.max(0.0),
            (-min_x).max(0.0),
            max_y.max(0.0),
            (-min_y).max(0.0),
        )
    };
    let (mut xmax, xmin, mut ymax, ymin) = box_of(&local_components[0]);
    for &(idx, (x, y)) in &local_components[0] {
        out[idx] = (x, y);
    }
    for component in local_components.iter().skip(1) {
        let (xp, xn, yp, yn) = box_of(component);
        let mut shift_x = 0.0;
        let mut shift_y = 0.0;
        let mut xshift = true;
        if xmax + xmin > ymax + ymin {
            xshift = false;
        }
        if xshift {
            shift_x = xmax + xn + 1.0;
            xmax += xp + xn + 1.0;
        } else {
            shift_y = ymax + yn + 1.0;
            ymax += yp + yn + 1.0;
        }
        for &(idx, (x, y)) in component {
            out[idx] = (x + shift_x, y + shift_y);
        }
    }

    Ok(out)
}

fn strict_rdkit_subset_coords(mol: &Molecule) -> Option<Vec<(f64, f64)>> {
    let n = mol.atoms.len();
    if n == 1 {
        return Some(vec![(0.0, 0.0)]);
    }
    if n > 3 {
        return None;
    }
    if mol.bonds.len() != n - 1 {
        return None;
    }
    if mol.bonds.iter().any(|b| {
        !matches!(
            b.order,
            BondOrder::Single | BondOrder::Double | BondOrder::Triple
        )
    }) {
        return None;
    }
    let mut degree = vec![0usize; n];
    for b in &mol.bonds {
        degree[b.begin_atom] += 1;
        degree[b.end_atom] += 1;
    }
    if degree.iter().any(|d| *d > 2) {
        return None;
    }
    if n == 2 {
        return Some(vec![(-0.7500, 0.0), (0.7500, -0.0)]);
    }
    None
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum RdkitHybridization {
    Unspecified,
    S,
    Sp,
    Sp2,
    Sp3,
    Sp3d,
    Sp3d2,
}

fn canonicalize_component_rdkit_like(component: &mut Vec<(usize, (f64, f64))>) {
    if component.len() <= 1 {
        return;
    }
    let n = component.len() as f64;
    let (mut cx, mut cy) = (0.0f64, 0.0f64);
    for &(_, (x, y)) in component.iter() {
        cx += x;
        cy += y;
    }
    cx /= n;
    cy /= n;

    let (mut xx, mut xy, mut yy) = (0.0f64, 0.0f64, 0.0f64);
    let mut centered: Vec<(f64, f64)> = Vec::with_capacity(component.len());
    for &(_, (x, y)) in component.iter() {
        let px = x - cx;
        let py = y - cy;
        centered.push((px, py));
        xx += px * px;
        xy += px * py;
        yy += py * py;
    }

    let d = ((xx - yy) * (xx - yy) + 4.0 * xy * xy).sqrt();
    let mut eig1 = (2.0 * xy, (yy - xx) + d);
    let eig1_len = norm(eig1);
    if eig1_len <= 1e-4 {
        for (i, (_, pos)) in component.iter_mut().enumerate() {
            *pos = centered[i];
        }
        return;
    }
    let e_val1 = (xx + yy + d) / 2.0;
    eig1 = (eig1.0 / eig1_len, eig1.1 / eig1_len);

    let mut eig2 = (2.0 * xy, (yy - xx) - d);
    let e_val2 = (xx + yy - d) / 2.0;
    let eig2_len = norm(eig2);
    if eig2_len > 1e-4 {
        eig2 = (eig2.0 / eig2_len, eig2.1 / eig2_len);
        if e_val2 > e_val1 {
            std::mem::swap(&mut eig1, &mut eig2);
        }
    }

    for (i, (_, pos)) in component.iter_mut().enumerate() {
        let (px, py) = centered[i];
        // trans = [[eig1.x, eig1.y],[-eig1.y,eig1.x]]
        let rx = px * eig1.0 + py * eig1.1;
        let ry = -px * eig1.1 + py * eig1.0;
        *pos = (rx, ry);
    }
}

#[derive(Clone)]
struct TreeEmbeddedAtom {
    loc: (f64, f64),
    normal: (f64, f64),
    ccw: bool,
    cis_trans_nbr: Option<usize>,
    angle: f64,
    nbr1: Option<usize>,
    nbr2: Option<usize>,
    rot_dir: i32,
    pending: Vec<usize>,
}

fn atom_depict_rank(atomic_num: u8, degree: usize) -> usize {
    let anum = if atomic_num == 1 {
        1000usize
    } else {
        atomic_num as usize
    };
    100 * anum + degree
}

fn n_outer_electrons_rdkit(atomic_num: u8) -> Option<i32> {
    // RDKit PeriodicTable::getNouterElecs(), values mirrored from chem-core's
    // RDKit 2025.03.5 valence support table.
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

fn rdkit_num_bonds_plus_lone_pairs(
    mol: &Molecule,
    atom_index: usize,
    degree: usize,
    explicit_valence: u8,
    implicit_hydrogens: u8,
) -> Option<i32> {
    // RDKit ConjugHybrid.cpp::numBondsPlusLonePairs().
    let atom = &mol.atoms[atom_index];
    let mut deg = degree as i32 + atom.explicit_hydrogens as i32 + implicit_hydrogens as i32;
    for bond in &mol.bonds {
        if bond.begin_atom != atom_index && bond.end_atom != atom_index {
            continue;
        }
        if matches!(bond.order, BondOrder::Null)
            || (matches!(bond.order, BondOrder::Dative) && bond.end_atom != atom_index)
        {
            deg -= 1;
        }
    }

    if atom.atomic_num <= 1 {
        return Some(deg);
    }
    let nouter = n_outer_electrons_rdkit(atom.atomic_num)?;
    let total_valence = explicit_valence as i32 + implicit_hydrogens as i32;
    let charge = atom.formal_charge as i32;
    let free_electrons = nouter - (total_valence + charge);
    if total_valence + nouter - charge < 8 {
        let radicals = atom.num_radical_electrons as i32;
        Some(deg + (free_electrons - radicals) / 2 + radicals)
    } else {
        Some(deg + free_electrons / 2)
    }
}

fn rdkit_hybridizations_for_depict(
    mol: &Molecule,
    degree: &[usize],
) -> Result<Vec<RdkitHybridization>, MolWriteError> {
    let assignment = crate::assign_valence(mol, crate::ValenceModel::RdkitLike).map_err(|_| {
        MolWriteError::UnsupportedSubset("RDKit hybridization valence assignment failed")
    })?;

    let mut out = Vec::with_capacity(mol.atoms.len());
    for (idx, atom) in mol.atoms.iter().enumerate() {
        if atom.atomic_num == 0 {
            out.push(RdkitHybridization::Unspecified);
            continue;
        }
        let norbs = if atom.atomic_num < 89 {
            rdkit_num_bonds_plus_lone_pairs(
                mol,
                idx,
                degree[idx],
                assignment.explicit_valence[idx],
                assignment.implicit_hydrogens[idx],
            )
            .ok_or(MolWriteError::UnsupportedSubset(
                "RDKit hybridization outer-electron lookup failed",
            ))?
        } else {
            degree[idx] as i32
                + atom.explicit_hydrogens as i32
                + assignment.implicit_hydrogens[idx] as i32
        };

        out.push(match norbs {
            0 | 1 => RdkitHybridization::S,
            2 => RdkitHybridization::Sp,
            3 => RdkitHybridization::Sp2,
            4 => RdkitHybridization::Sp3,
            5 => RdkitHybridization::Sp3d,
            6 => RdkitHybridization::Sp3d2,
            _ => RdkitHybridization::Unspecified,
        });
    }
    Ok(out)
}

fn rdkit_twice_bond_type(order: BondOrder) -> i64 {
    // RDKit Bond.cpp::getTwiceBondType().
    match order {
        BondOrder::Null => 0,
        BondOrder::Single => 2,
        BondOrder::Double => 4,
        BondOrder::Triple => 6,
        BondOrder::Quadruple => 8,
        BondOrder::Aromatic => 3,
        BondOrder::Dative => 2,
    }
}

fn rdkit_cip_invariants(mol: &Molecule) -> Vec<i64> {
    // RDKit Chirality.cpp::buildCIPInvariants() without isotope/map-number
    // variation; both are currently represented as their default values here.
    const MASS_BITS: i64 = 10;
    const MAX_MASS: i64 = 1 << MASS_BITS;
    mol.atoms
        .iter()
        .map(|atom| {
            let mut invariant = i64::from(atom.atomic_num % 128);
            let mass = MAX_MASS / 2;
            invariant = (invariant << MASS_BITS) | mass;
            let mapnum = 0;
            (invariant << 10) | mapnum
        })
        .collect()
}

fn find_cip_segments(sorted: &mut [(Vec<i64>, usize, i64)]) -> (Vec<(usize, usize)>, usize) {
    let mut segments = Vec::new();
    if sorted.is_empty() {
        return (segments, 0);
    }

    let mut num_independent = sorted.len();
    let mut running_rank = 0i64;
    sorted[0].2 = running_rank;
    let mut current = 0usize;
    let mut in_equal_section = false;

    for i in 1..sorted.len() {
        if sorted[current].0 == sorted[i].0 {
            sorted[i].2 = running_rank;
            num_independent -= 1;
            if !in_equal_section {
                in_equal_section = true;
                segments.push((i - 1, 0));
            }
        } else {
            running_rank += 1;
            sorted[i].2 = running_rank;
            current = i;

            if in_equal_section {
                if let Some((_, last)) = segments.last_mut() {
                    *last = i;
                }
                in_equal_section = false;
            }
        }
    }

    if in_equal_section {
        if let Some((_, last)) = segments.last_mut() {
            *last = sorted.len() - 1;
        }
    }

    (segments, num_independent)
}

fn recompute_cip_ranks(sorted: &[(Vec<i64>, usize, i64)], ranks: &mut [i64]) {
    for entry in sorted {
        ranks[entry.1] = entry.2;
    }
}

pub(crate) fn rdkit_cip_ranks_for_depict(mol: &Molecule) -> Vec<i64> {
    // RDKit legacy Chirality.cpp::assignAtomCIPRanks()/iterateCIPRanks() enough
    // for depictor ordering. This preserves duplicated neighbor entries for
    // multiple bonds through getTwiceBondType().
    let n = mol.atoms.len();
    let invars = rdkit_cip_invariants(mol);
    let mut cip_entries: Vec<Vec<i64>> = invars.iter().copied().map(|v| vec![v]).collect();
    let mut sortable: Vec<(Vec<i64>, usize, i64)> = cip_entries
        .iter()
        .cloned()
        .enumerate()
        .map(|(idx, entry)| (entry, idx, -1))
        .collect();

    sortable.sort_by(|a, b| a.0.cmp(&b.0));
    let (mut needs_sorting, mut num_ranks) = find_cip_segments(&mut sortable);
    let mut ranks = vec![0i64; n];
    recompute_cip_ranks(&sortable, &mut ranks);

    for i in 0..n {
        cip_entries[i][0] = i64::from(mol.atoms[i].atomic_num);
        cip_entries[i].push(ranks[i]);
    }
    for item in &mut sortable {
        item.0 = cip_entries[item.1].clone();
    }

    let mut bond_features = vec![Vec::<(i64, usize)>::new(); n];
    for atom_idx in 0..n {
        for bond in &mol.bonds {
            let nbr_idx = if bond.begin_atom == atom_idx {
                bond.end_atom
            } else if bond.end_atom == atom_idx {
                bond.begin_atom
            } else {
                continue;
            };
            bond_features[atom_idx].push((rdkit_twice_bond_type(bond.order), nbr_idx));
        }
    }

    let max_its = n / 2 + 1;
    let mut num_its = 0usize;
    let mut last_num_ranks: Option<usize> = None;

    while !needs_sorting.is_empty()
        && num_its < max_its
        && last_num_ranks.map_or(true, |last| last < num_ranks)
    {
        for index in 0..n {
            let mut features = bond_features[index].clone();
            if features.len() > 1 {
                features.sort_by(|a, b| ranks[b.1].cmp(&ranks[a.1]));
            }
            for (count, nbr_idx) in features {
                for _ in 0..count {
                    cip_entries[index].push(ranks[nbr_idx] + 1);
                }
            }
        }
        for item in &mut sortable {
            item.0 = cip_entries[item.1].clone();
        }

        last_num_ranks = Some(num_ranks);
        for &(first, last) in &needs_sorting {
            sortable[first..=last].sort_by(|a, b| a.0.cmp(&b.0));
        }
        let found = find_cip_segments(&mut sortable);
        needs_sorting = found.0;
        num_ranks = found.1;
        recompute_cip_ranks(&sortable, &mut ranks);

        if Some(num_ranks) != last_num_ranks {
            for i in 0..n {
                cip_entries[i].resize(3, 0);
                cip_entries[i][2] = ranks[i];
            }
            for item in &mut sortable {
                item.0 = cip_entries[item.1].clone();
            }
        }

        num_its += 1;
    }

    ranks
}

fn rdkit_rank_atoms_by_rank(
    mol: &Molecule,
    atoms: &mut [usize],
    degree: &[usize],
    cip_ranks: &[i64],
) {
    atoms.sort_by_key(|&idx| {
        let rank = cip_ranks.get(idx).copied().unwrap_or_else(|| idx as i64);
        (rank, idx)
    });

    // Mirrors DepictUtils.cpp::rankAtomsByRank() fallback when no CIP rank is
    // present. In this code path CIP ranks are always generated above, so this
    // is intentionally only used for length mismatches.
    if cip_ranks.len() != mol.atoms.len() {
        atoms.sort_by_key(|&idx| {
            let depict = atom_depict_rank(mol.atoms[idx].atomic_num, degree[idx]);
            depict * mol.atoms.len() + idx
        });
    }
}

fn rdkit_set_nbr_order(
    aid: usize,
    nbrs: &[usize],
    mol: &Molecule,
    adjacency: &[Vec<usize>],
    degree: &[usize],
    cip_ranks: &[i64],
) -> Vec<usize> {
    let mut ref_atom: Option<usize> = None;
    for &nb in &adjacency[aid] {
        if !nbrs.contains(&nb) {
            ref_atom = Some(nb);
        }
    }

    let mut thold = nbrs.to_vec();
    if let Some(r) = ref_atom {
        thold.push(r);
    }
    if thold.len() <= 3 {
        let mut out = nbrs.to_vec();
        rdkit_rank_atoms_by_rank(mol, &mut out, degree, cip_ranks);
        return out;
    }

    rdkit_rank_atoms_by_rank(mol, &mut thold, degree, cip_ranks);

    let ln = thold.len();
    thold.swap(ln - 3, ln - 2);

    if let Some(r) = ref_atom {
        if let Some(pos) = thold.iter().position(|&a| a == r) {
            let mut out = Vec::with_capacity(nbrs.len());
            out.extend_from_slice(&thold[pos + 1..]);
            out.extend_from_slice(&thold[..pos]);
            return out;
        }
    }
    thold
}

fn compute_sub_angle(degree: usize, hybridization: RdkitHybridization) -> f64 {
    // RDKit DepictUtils.h::computeSubAngle().
    match hybridization {
        RdkitHybridization::Unspecified | RdkitHybridization::Sp3 => {
            if degree == 4 {
                PI / 2.0
            } else {
                2.0 * PI / 3.0
            }
        }
        RdkitHybridization::Sp2 => 2.0 * PI / 3.0,
        RdkitHybridization::S
        | RdkitHybridization::Sp
        | RdkitHybridization::Sp3d
        | RdkitHybridization::Sp3d2 => 2.0 * PI / degree as f64,
    }
}

fn rotate(v: (f64, f64), angle: f64) -> (f64, f64) {
    let (c, s) = (angle.cos(), angle.sin());
    (v.0 * c - v.1 * s, v.0 * s + v.1 * c)
}

fn rotate_around(p: (f64, f64), center: (f64, f64), angle: f64) -> (f64, f64) {
    let rel = (p.0 - center.0, p.1 - center.1);
    let rot = rotate(rel, angle);
    (center.0 + rot.0, center.1 + rot.1)
}

fn norm(v: (f64, f64)) -> f64 {
    (v.0 * v.0 + v.1 * v.1).sqrt()
}

fn normalize(v: (f64, f64)) -> (f64, f64) {
    let l = norm(v);
    if l <= 1e-12 {
        (1.0, 0.0)
    } else {
        (v.0 / l, v.1 / l)
    }
}

fn rotation_dir(center: (f64, f64), loc1: (f64, f64), loc2: (f64, f64), rem_angle: f64) -> i32 {
    let pt1 = (loc1.0 - center.0, loc1.1 - center.1);
    let pt2 = (loc2.0 - center.0, loc2.1 - center.1);
    let cross = (pt1.0 * pt2.1 - pt1.1 * pt2.0) * (PI - rem_angle);
    if cross >= 0.0 { -1 } else { 1 }
}

fn place_acyclic_tree_rdkit_like(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
    degree: &[usize],
    hybridizations: &[RdkitHybridization],
    cip_ranks: &[i64],
    seeded_coords: Option<&std::collections::HashMap<usize, (f64, f64)>>,
) -> Option<Vec<(usize, (f64, f64))>> {
    if comp.is_empty() {
        return None;
    }
    let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();
    let bond_count_in_comp = mol
        .bonds
        .iter()
        .filter(|b| comp_set.contains(&b.begin_atom) && comp_set.contains(&b.end_atom))
        .count();
    let is_tree_comp = bond_count_in_comp + 1 == comp.len();
    let n = mol.atoms.len();
    let mut states: Vec<Option<TreeEmbeddedAtom>> = vec![None; n];
    let mut unembedded: std::collections::HashSet<usize> = comp.iter().copied().collect();
    let mut attach_pts = std::collections::VecDeque::<usize>::new();
    let mut attach_in_queue = std::collections::HashSet::<usize>::new();

    let mut root: Option<usize> = None;
    let mut stereo_seeded = false;
    if let Some(seeded) = seeded_coords {
        if !seeded.is_empty() {
            let mut cx = 0.0f64;
            let mut cy = 0.0f64;
            let mut cnt = 0usize;
            for (&aid, &(x, y)) in seeded {
                if comp_set.contains(&aid) {
                    cx += x;
                    cy += y;
                    cnt += 1;
                }
            }
            if cnt > 0 {
                cx /= cnt as f64;
                cy /= cnt as f64;
            }
            for (&aid, &(x, y)) in seeded {
                if !comp_set.contains(&aid) {
                    continue;
                }
                let mut nrm = normalize((x - cx, y - cy));
                if norm(nrm) <= 1e-12 {
                    nrm = (1.0, 0.0);
                }
                states[aid] = Some(TreeEmbeddedAtom {
                    loc: (x, y),
                    normal: nrm,
                    ccw: true,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: None,
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                });
                unembedded.remove(&aid);
            }
            root = comp
                .iter()
                .copied()
                .filter(|a| states[*a].is_some())
                .min_by_key(|&idx| {
                    let depict = atom_depict_rank(mol.atoms[idx].atomic_num, degree[idx]);
                    depict * mol.atoms.len() + idx
                });
        }
    }
    if root.is_none() && seeded_coords.is_none() {
        if let Some(bond) = mol.bonds.iter().find(|bond| {
            comp_set.contains(&bond.begin_atom)
                && comp_set.contains(&bond.end_atom)
                && matches!(bond.order, BondOrder::Double)
                && matches!(bond.stereo, BondStereo::Cis | BondStereo::Trans)
                && bond.stereo_atoms.len() == 2
        }) {
            let begin = bond.begin_atom;
            let end = bond.end_atom;
            let begin_control = bond.stereo_atoms[0];
            let end_control = bond.stereo_atoms[1];

            states[begin] = Some(TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (0.0, -1.0),
                ccw: false,
                cis_trans_nbr: Some(begin_control),
                angle: -1.0,
                nbr1: Some(end),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            });
            let (end_normal, end_ccw) = if matches!(bond.stereo, BondStereo::Cis) {
                ((0.0, -1.0), true)
            } else {
                ((0.0, 1.0), false)
            };
            states[end] = Some(TreeEmbeddedAtom {
                loc: (1.5, 0.0),
                normal: end_normal,
                ccw: end_ccw,
                cis_trans_nbr: Some(end_control),
                angle: -1.0,
                nbr1: Some(begin),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            });
            unembedded.remove(&begin);
            unembedded.remove(&end);
            root = Some(begin);
            stereo_seeded = true;
        }
    }
    if root.is_none() {
        let r = *comp.iter().min_by_key(|&&idx| {
            let depict = atom_depict_rank(mol.atoms[idx].atomic_num, degree[idx]);
            depict * mol.atoms.len() + idx
        })?;
        states[r] = Some(TreeEmbeddedAtom {
            loc: (0.0, 0.0),
            normal: (1.0, 0.0),
            ccw: true,
            cis_trans_nbr: None,
            angle: -1.0,
            nbr1: None,
            nbr2: None,
            rot_dir: 0,
            pending: Vec::new(),
        });
        unembedded.remove(&r);
        root = Some(r);
    }

    let find_num_neigh =
        |states: &Vec<Option<TreeEmbeddedAtom>>, pt: (f64, f64), radius: f64| -> i32 {
            states
                .iter()
                .flatten()
                .filter(|st| {
                    let dx = st.loc.0 - pt.0;
                    let dy = st.loc.1 - pt.1;
                    (dx * dx + dy * dy).sqrt() < radius
                })
                .count() as i32
        };

    let compute_angle = |center: (f64, f64), p1: (f64, f64), p2: (f64, f64)| -> f64 {
        let v1 = (p1.0 - center.0, p1.1 - center.1);
        let v2 = (p2.0 - center.0, p2.1 - center.1);
        let d1 = (v1.0 * v1.0 + v1.1 * v1.1).sqrt();
        let d2 = (v2.0 * v2.0 + v2.1 * v2.1).sqrt();
        if d1 <= 1e-12 || d2 <= 1e-12 {
            return PI;
        }
        let mut c = (v1.0 * v2.0 + v1.1 * v2.1) / (d1 * d2);
        c = c.clamp(-1.0, 1.0);
        c.acos()
    };

    let update_new_neighs =
        |aid: usize,
         states: &mut Vec<Option<TreeEmbeddedAtom>>,
         unembedded: &std::collections::HashSet<usize>,
         attach_pts: &mut std::collections::VecDeque<usize>,
         attach_in_queue: &mut std::collections::HashSet<usize>| {
            if !is_tree_comp && states[aid].is_some() {
                let st_loc = states[aid].as_ref().map(|s| s.loc).unwrap_or((0.0, 0.0));
                let done_nbrs: Vec<usize> = adjacency[aid]
                    .iter()
                    .copied()
                    .filter(|nb| {
                        comp_set.contains(nb) && !unembedded.contains(nb) && states[*nb].is_some()
                    })
                    .collect();
                if done_nbrs.len() >= 3 {
                    let mut angle_pairs: Vec<(f64, (usize, usize))> = Vec::new();
                    for i in 0..done_nbrs.len() {
                        for j in i + 1..done_nbrs.len() {
                            let a = done_nbrs[i];
                            let b = done_nbrs[j];
                            let pa = states[a]
                                .as_ref()
                                .map(|n| n.loc)
                                .unwrap_or((st_loc.0 - 1.0, st_loc.1));
                            let pb = states[b]
                                .as_ref()
                                .map(|n| n.loc)
                                .unwrap_or((st_loc.0 + 1.0, st_loc.1));
                            angle_pairs.push((compute_angle(st_loc, pa, pb), (a, b)));
                        }
                    }
                    angle_pairs
                        .sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
                    if let Some((w_ang, (wnb1, wnb2))) = angle_pairs.last().copied() {
                        let mut nb1 = wnb1;
                        let mut nb2 = wnb2;
                        for &(_, (a, b)) in angle_pairs.iter().rev() {
                            if a == wnb1 && b != wnb2 {
                                nb1 = b;
                                nb2 = wnb1;
                                break;
                            } else if b == wnb1 && a != wnb2 {
                                nb1 = a;
                                nb2 = wnb1;
                                break;
                            } else if a == wnb2 && b != wnb1 {
                                nb1 = b;
                                nb2 = wnb2;
                                break;
                            } else if b == wnb2 && a != wnb1 {
                                nb1 = a;
                                nb2 = wnb2;
                                break;
                            }
                        }
                        let p1 = states[nb1]
                            .as_ref()
                            .map(|n| n.loc)
                            .unwrap_or((st_loc.0 - 1.0, st_loc.1));
                        let p2 = states[nb2]
                            .as_ref()
                            .map(|n| n.loc)
                            .unwrap_or((st_loc.0 + 1.0, st_loc.1));
                        if let Some(st) = states[aid].as_mut() {
                            st.nbr1 = Some(nb1);
                            st.nbr2 = Some(nb2);
                            st.angle = 2.0 * PI - w_ang;
                            st.rot_dir = rotation_dir(st_loc, p1, p2, w_ang);
                        }
                    }
                }
            }

            let mut neighs: Vec<usize> = adjacency[aid]
                .iter()
                .copied()
                .filter(|nb| comp_set.contains(nb) && unembedded.contains(nb))
                .collect();
            if !neighs.is_empty() && (degree[aid] < 4 || neighs.len() < 3) {
                rdkit_rank_atoms_by_rank(mol, &mut neighs, degree, cip_ranks);
            } else if degree[aid] >= 4 && neighs.len() >= 3 {
                neighs = rdkit_set_nbr_order(aid, &neighs, mol, adjacency, degree, cip_ranks);
            }
            if let Some(st) = states[aid].as_mut() {
                st.pending = neighs;
                if !st.pending.is_empty() && !attach_in_queue.contains(&aid) {
                    attach_pts.push_back(aid);
                    attach_in_queue.insert(aid);
                }
            }
        };

    if stereo_seeded {
        for &aid in comp {
            if states[aid].is_some() {
                update_new_neighs(
                    aid,
                    &mut states,
                    &unembedded,
                    &mut attach_pts,
                    &mut attach_in_queue,
                );
            }
        }
    } else if let Some(seeded) = seeded_coords {
        if !seeded.is_empty() {
            for &aid in comp {
                if states[aid].is_some() {
                    update_new_neighs(
                        aid,
                        &mut states,
                        &unembedded,
                        &mut attach_pts,
                        &mut attach_in_queue,
                    );
                }
            }
        } else {
            update_new_neighs(
                root.expect("root initialized"),
                &mut states,
                &unembedded,
                &mut attach_pts,
                &mut attach_in_queue,
            );
        }
    } else {
        update_new_neighs(
            root.expect("root initialized"),
            &mut states,
            &unembedded,
            &mut attach_pts,
            &mut attach_in_queue,
        );
    }

    while let Some(aid) = attach_pts.pop_front() {
        attach_in_queue.remove(&aid);
        let neighs = states[aid]
            .as_ref()
            .map(|st| st.pending.clone())
            .unwrap_or_default();
        if let Some(st) = states[aid].as_mut() {
            st.pending.clear();
        }
        let total = neighs.len();
        for (i, nb) in neighs.into_iter().enumerate() {
            if !unembedded.contains(&nb) {
                continue;
            }

            let ref_state = states[aid].as_ref()?.clone();
            let nnbr = (total - i) as f64;
            if ref_state.angle > 0.0 {
                let rem_angle = 2.0 * PI - ref_state.angle;
                let mut curr_angle = rem_angle / (1.0 + nnbr);
                if let Some(st) = states[aid].as_mut() {
                    st.angle += curr_angle;
                }
                let nb1 = states[ref_state.nbr1?].as_ref()?.loc;
                let nb2 = states[ref_state.nbr2?].as_ref()?.loc;
                if states[aid].as_ref()?.rot_dir == 0 {
                    let rd = rotation_dir(ref_state.loc, nb1, nb2, rem_angle);
                    if let Some(st) = states[aid].as_mut() {
                        st.rot_dir = rd;
                    }
                }
                curr_angle *= states[aid].as_ref()?.rot_dir as f64;
                let mut curr_loc = rotate_around(nb2, ref_state.loc, curr_angle);
                if rem_angle.abs() - PI < 1e-3 {
                    let curr_loc2 = rotate_around(nb2, ref_state.loc, -curr_angle);
                    if find_num_neigh(&states, curr_loc, 0.5)
                        > find_num_neigh(&states, curr_loc2, 0.5)
                    {
                        curr_loc = curr_loc2;
                    }
                }
                if let Some(st) = states[aid].as_mut() {
                    st.nbr2 = Some(nb);
                }

                let tpt = (curr_loc.0 - ref_state.loc.0, curr_loc.1 - ref_state.loc.1);
                let mut normal = normalize((-tpt.1, tpt.0));
                let tp1 = (curr_loc.0 + normal.0, curr_loc.1 + normal.1);
                let tp2 = (curr_loc.0 - normal.0, curr_loc.1 - normal.1);
                let nccw = find_num_neigh(&states, tp1, 2.5);
                let ncw = find_num_neigh(&states, tp2, 2.5);
                let (ccw, out_normal) = if nccw < ncw {
                    (false, normal)
                } else {
                    normal = (-normal.0, -normal.1);
                    (true, normal)
                };
                states[nb] = Some(TreeEmbeddedAtom {
                    loc: curr_loc,
                    normal: out_normal,
                    ccw,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: Some(aid),
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                });
            } else {
                let mut ref_ccw = ref_state.ccw;
                let mut curr_vec = ref_state.normal;
                if ref_state.cis_trans_nbr.is_some_and(|ct| ct != nb) {
                    ref_ccw = !ref_ccw;
                    curr_vec = (-curr_vec.0, -curr_vec.1);
                }
                let degree_here = degree[aid];
                let mut angle = compute_sub_angle(degree_here, hybridizations[aid]);
                let mut flip_norm = false;

                if states[aid].as_ref()?.nbr1.is_some() {
                    if let Some(st) = states[aid].as_mut() {
                        st.angle = angle;
                        st.nbr2 = Some(nb);
                    }
                } else {
                    let norm0 = states[aid].as_ref()?.normal;
                    let rot = rotate(norm0, angle);
                    if let Some(st) = states[aid].as_mut() {
                        st.normal = rot;
                        st.nbr1 = Some(nb);
                    }
                    flip_norm = true;
                }

                angle -= PI / 2.0;
                if !ref_ccw {
                    angle *= -1.0;
                }
                curr_vec = rotate(curr_vec, angle);
                let curr_loc = (
                    ref_state.loc.0 + 1.5 * curr_vec.0,
                    ref_state.loc.1 + 1.5 * curr_vec.1,
                );

                let tpt = (ref_state.loc.0 - curr_loc.0, ref_state.loc.1 - curr_loc.1);
                let mut normal = (-tpt.1, tpt.0);
                if ref_ccw ^ flip_norm {
                    normal = (-normal.0, -normal.1);
                }
                normal = normalize(normal);
                states[nb] = Some(TreeEmbeddedAtom {
                    loc: curr_loc,
                    normal,
                    ccw: (!ref_ccw) ^ flip_norm,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: Some(aid),
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                });
            }

            unembedded.remove(&nb);
            update_new_neighs(
                nb,
                &mut states,
                &unembedded,
                &mut attach_pts,
                &mut attach_in_queue,
            );
        }
        update_new_neighs(
            aid,
            &mut states,
            &unembedded,
            &mut attach_pts,
            &mut attach_in_queue,
        );
    }

    if !unembedded.is_empty() {
        return None;
    }

    let mut out: Vec<(usize, (f64, f64))> = comp
        .iter()
        .copied()
        .filter_map(|a| states[a].as_ref().map(|st| (a, st.loc)))
        .collect();
    if out.len() != comp.len() {
        return None;
    }
    out.sort_by_key(|(a, _)| *a);
    Some(out)
}

fn place_unicyclic_component(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
) -> Option<Vec<(usize, (f64, f64))>> {
    let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();

    let mut deg_in_comp = vec![0usize; mol.atoms.len()];
    for &a in comp {
        deg_in_comp[a] = adjacency[a].iter().filter(|n| comp_set.contains(n)).count();
    }

    // Peel leaves to get ring-core atoms.
    let mut removed = std::collections::HashSet::<usize>::new();
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

    let mut ring_atoms: Vec<usize> = comp
        .iter()
        .copied()
        .filter(|a| !removed.contains(a))
        .collect();
    if ring_atoms.len() < 3 {
        return None;
    }

    // Ensure ring core is a simple cycle.
    let ring_set: std::collections::HashSet<usize> = ring_atoms.iter().copied().collect();
    for &a in &ring_atoms {
        let d = adjacency[a].iter().filter(|n| ring_set.contains(n)).count();
        if d != 2 {
            return None;
        }
    }

    // Traverse cycle order deterministically.
    ring_atoms.sort_unstable();
    let start = ring_atoms[0];
    let mut order = Vec::with_capacity(ring_atoms.len());
    order.push(start);
    let mut prev = usize::MAX;
    let mut cur = start;
    loop {
        let mut next_candidates: Vec<usize> = adjacency[cur]
            .iter()
            .copied()
            .filter(|n| ring_set.contains(n) && *n != prev)
            .collect();
        if next_candidates.is_empty() {
            return None;
        }
        next_candidates.sort_unstable();
        let nxt = next_candidates[0];
        if nxt == start {
            break;
        }
        order.push(nxt);
        prev = cur;
        cur = nxt;
        if order.len() > ring_atoms.len() + 1 {
            return None;
        }
    }
    if order.len() != ring_atoms.len() {
        return None;
    }

    let n = order.len();
    let bond_len = 1.5f64;
    let radius = bond_len / (2.0 * (PI / n as f64).sin());
    let mut out: Vec<(usize, (f64, f64))> = Vec::new();
    let mut placed = std::collections::HashSet::<usize>::new();

    // TODO: replace this traversal with RDKit MolOps::symmetrizeSSSR ring order.
    let aromatic_ring = mol.bonds.iter().any(|b| {
        ring_set.contains(&b.begin_atom)
            && ring_set.contains(&b.end_atom)
            && matches!(b.order, BondOrder::Aromatic)
    });
    let direction = if aromatic_ring && n == 5 { 1.0 } else { -1.0 };
    for (i, &a) in order.iter().enumerate() {
        let theta = direction * 2.0 * PI * i as f64 / n as f64;
        out.push((a, (radius * theta.cos(), radius * theta.sin())));
        placed.insert(a);
    }

    let pos_map = |idx: usize, data: &Vec<(usize, (f64, f64))>| -> (f64, f64) {
        data.iter()
            .find(|(a, _)| *a == idx)
            .map(|(_, p)| *p)
            .expect("placed atom")
    };
    let normalize = |x: f64, y: f64| -> (f64, f64) {
        let l = (x * x + y * y).sqrt();
        if l < 1e-12 {
            (1.0, 0.0)
        } else {
            (x / l, y / l)
        }
    };
    let rotate = |dx: f64, dy: f64, deg: f64| -> (f64, f64) {
        let rad = deg * PI / 180.0;
        let (c, s) = (rad.cos(), rad.sin());
        (dx * c - dy * s, dx * s + dy * c)
    };

    fn place_branch(
        mol: &Molecule,
        adjacency: &[Vec<usize>],
        comp_set: &std::collections::HashSet<usize>,
        placed: &mut std::collections::HashSet<usize>,
        out: &mut Vec<(usize, (f64, f64))>,
        parent: usize,
        current: usize,
        in_dir: (f64, f64),
        depth: usize,
    ) {
        let current_pos = out
            .iter()
            .find(|(a, _)| *a == current)
            .map(|(_, p)| *p)
            .expect("current placed");
        let mut children: Vec<usize> = adjacency[current]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n) && *n != parent && !placed.contains(n))
            .collect();
        children.sort_unstable();
        if children.is_empty() {
            return;
        }

        let angle_set: Vec<f64> = match children.len() {
            1 => vec![if depth % 2 == 0 { -60.0 } else { 60.0 }],
            2 => vec![-60.0, 60.0],
            3 => vec![-90.0, 0.0, 90.0],
            _ => {
                let mut v = Vec::with_capacity(children.len());
                let step = 180.0 / (children.len() as f64 - 1.0);
                for i in 0..children.len() {
                    v.push(-90.0 + step * i as f64);
                }
                v
            }
        };

        for (i, child) in children.into_iter().enumerate() {
            let angle = angle_set[i];
            let rad = angle * PI / 180.0;
            let (c, s) = (rad.cos(), rad.sin());
            let new_dir = (in_dir.0 * c - in_dir.1 * s, in_dir.0 * s + in_dir.1 * c);
            let child_pos = (
                current_pos.0 + 1.5 * new_dir.0,
                current_pos.1 + 1.5 * new_dir.1,
            );
            out.push((child, child_pos));
            placed.insert(child);
            place_branch(
                mol,
                adjacency,
                comp_set,
                placed,
                out,
                current,
                child,
                new_dir,
                depth + 1,
            );
        }
    }

    for &ra in &order {
        let (x, y) = pos_map(ra, &out);
        let (dx, dy) = normalize(x, y);
        let mut nbrs: Vec<usize> = adjacency[ra]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n) && !ring_set.contains(n) && !placed.contains(n))
            .collect();
        nbrs.sort_unstable();
        for (i, nb) in nbrs.into_iter().enumerate() {
            let angle = if i == 0 { 0.0 } else { (i as f64) * 30.0 };
            let (ndx, ndy) = rotate(dx, dy, angle);
            let pos = (x + 1.5 * ndx, y + 1.5 * ndy);
            out.push((nb, pos));
            placed.insert(nb);
            place_branch(
                mol,
                adjacency,
                &comp_set,
                &mut placed,
                &mut out,
                ra,
                nb,
                (ndx, ndy),
                0,
            );
        }
    }

    if comp.iter().any(|a| !placed.contains(a)) {
        return None;
    }
    out.sort_by_key(|(a, _)| *a);
    Some(out)
}

fn regularize_isolated_aromatic_six_cycles(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
    local: &mut Vec<(usize, (f64, f64))>,
) {
    let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();
    let mut coord_map = std::collections::HashMap::<usize, (f64, f64)>::new();
    for &(a, p) in local.iter() {
        coord_map.insert(a, p);
    }

    let mut aromatic_adj = std::collections::HashMap::<usize, Vec<usize>>::new();
    for &a in comp {
        aromatic_adj.insert(a, Vec::new());
    }
    for b in &mol.bonds {
        if !comp_set.contains(&b.begin_atom) || !comp_set.contains(&b.end_atom) {
            continue;
        }
        if !matches!(b.order, BondOrder::Aromatic) {
            continue;
        }
        aromatic_adj
            .entry(b.begin_atom)
            .or_default()
            .push(b.end_atom);
        aromatic_adj
            .entry(b.end_atom)
            .or_default()
            .push(b.begin_atom);
    }
    for nbs in aromatic_adj.values_mut() {
        nbs.sort_unstable();
        nbs.dedup();
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

    let mut seen = std::collections::HashSet::<Vec<usize>>::new();
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    let mut comp_sorted = comp.to_vec();
    comp_sorted.sort_unstable();
    for &start in &comp_sorted {
        let mut path = vec![start];
        let mut used = std::collections::HashSet::<usize>::new();
        used.insert(start);
        fn dfs6(
            start: usize,
            cur: usize,
            depth: usize,
            path: &mut Vec<usize>,
            used: &mut std::collections::HashSet<usize>,
            aromatic_adj: &std::collections::HashMap<usize, Vec<usize>>,
            seen: &mut std::collections::HashSet<Vec<usize>>,
            cycles: &mut Vec<Vec<usize>>,
        ) {
            if depth == 6 {
                if aromatic_adj
                    .get(&cur)
                    .is_some_and(|nbs| nbs.contains(&start))
                {
                    let c = canonical_cycle(path.clone());
                    if seen.insert(c.clone()) {
                        cycles.push(c);
                    }
                }
                return;
            }
            let mut nbs = aromatic_adj.get(&cur).cloned().unwrap_or_default();
            nbs.sort_unstable();
            for nb in nbs {
                if used.contains(&nb) {
                    continue;
                }
                if nb < start {
                    continue;
                }
                used.insert(nb);
                path.push(nb);
                dfs6(start, nb, depth + 1, path, used, aromatic_adj, seen, cycles);
                path.pop();
                used.remove(&nb);
            }
        }
        dfs6(
            start,
            start,
            1,
            &mut path,
            &mut used,
            &aromatic_adj,
            &mut seen,
            &mut cycles,
        );
    }
    if cycles.is_empty() {
        return;
    }

    let mut atom_cycle_count = std::collections::HashMap::<usize, usize>::new();
    for cyc in &cycles {
        for &a in cyc {
            *atom_cycle_count.entry(a).or_insert(0) += 1;
        }
    }

    for cyc in cycles {
        if cyc.len() != 6 {
            continue;
        }
        if cyc
            .iter()
            .any(|a| atom_cycle_count.get(a).copied().unwrap_or(0) != 1)
        {
            continue;
        }
        if cyc.iter().any(|a| !coord_map.contains_key(a)) {
            continue;
        }

        let ext_deg = |a: usize| -> usize {
            adjacency[a]
                .iter()
                .filter(|n| comp_set.contains(n) && !cyc.contains(n))
                .count()
        };

        let mut opposite_pairs: Vec<(usize, usize)> = Vec::new();
        for i in 0..6 {
            let j = (i + 3) % 6;
            if ext_deg(cyc[i]) > 0 && ext_deg(cyc[j]) > 0 {
                opposite_pairs.push((i, j));
            }
        }
        if opposite_pairs.is_empty() {
            continue;
        }

        let score_for =
            |candidate: &std::collections::HashMap<usize, (f64, f64)>, order: &[usize]| -> f64 {
                let mut score = 0.0;
                for &a in order {
                    let pa = candidate[&a];
                    for &nb in &adjacency[a] {
                        if !comp_set.contains(&nb) || order.contains(&nb) {
                            continue;
                        }
                        let pb = *coord_map.get(&nb).expect("neighbor coords");
                        let dx = pa.0 - pb.0;
                        let dy = pa.1 - pb.1;
                        let d = (dx * dx + dy * dy).sqrt();
                        score += (d - 1.5).abs() * 10.0;
                    }
                }
                // Favor orientations where ring-entry/exit bond directions are
                // aligned with external substituent directions at opposite anchors.
                let a0 = order[0];
                let a1 = order[1];
                let a5 = order[5];
                let a3 = order[3];
                let a2 = order[2];
                let a4 = order[4];
                let p0 = candidate[&a0];
                let p1 = candidate[&a1];
                let p5 = candidate[&a5];
                let p3 = candidate[&a3];
                let p2 = candidate[&a2];
                let p4 = candidate[&a4];

                let d01 = normalize((p1.0 - p0.0, p1.1 - p0.1));
                let d05 = normalize((p5.0 - p0.0, p5.1 - p0.1));
                let d32 = normalize((p2.0 - p3.0, p2.1 - p3.1));
                let d34 = normalize((p4.0 - p3.0, p4.1 - p3.1));

                let mut anchor_penalty = 0.0;
                for &nb in &adjacency[a0] {
                    if !comp_set.contains(&nb) || order.contains(&nb) {
                        continue;
                    }
                    let pb = *coord_map.get(&nb).expect("neighbor coords");
                    let v = normalize((pb.0 - p0.0, pb.1 - p0.1));
                    let c = (v.0 * d01.0 + v.1 * d01.1).max(v.0 * d05.0 + v.1 * d05.1);
                    anchor_penalty += (1.0 - c).abs();
                }
                for &nb in &adjacency[a3] {
                    if !comp_set.contains(&nb) || order.contains(&nb) {
                        continue;
                    }
                    let pb = *coord_map.get(&nb).expect("neighbor coords");
                    let v = normalize((pb.0 - p3.0, pb.1 - p3.1));
                    let c = (v.0 * d32.0 + v.1 * d32.1).max(v.0 * d34.0 + v.1 * d34.1);
                    anchor_penalty += (1.0 - c).abs();
                }
                score += anchor_penalty * 20.0;
                for &a in order {
                    let pa = candidate[&a];
                    for &b in comp {
                        if order.contains(&b) {
                            continue;
                        }
                        let pb = coord_map[&b];
                        let dx = pa.0 - pb.0;
                        let dy = pa.1 - pb.1;
                        let d = (dx * dx + dy * dy).sqrt();
                        if d < 1.0 {
                            score += (1.0 - d) * (1.0 - d) * 50.0;
                        }
                    }
                }
                for &a in order {
                    let po = coord_map[&a];
                    let pc = candidate[&a];
                    let dx = pc.0 - po.0;
                    let dy = pc.1 - po.1;
                    score += (dx * dx + dy * dy) * 0.5;
                }
                score
            };

        let mut best_coords: Option<std::collections::HashMap<usize, (f64, f64)>> = None;
        let mut best_score = f64::INFINITY;
        for (i0, i3) in opposite_pairs {
            let mut order = cyc.clone();
            order.rotate_left(i0);
            if order[3] != cyc[i3] {
                order.reverse();
                let pos = order.iter().position(|a| *a == cyc[i0]).unwrap_or(0);
                order.rotate_left(pos);
                if order[3] != cyc[i3] {
                    continue;
                }
            }
            let a0 = order[0];
            let a3 = order[3];
            let p0 = *coord_map.get(&a0).expect("coord exists");
            let p3 = *coord_map.get(&a3).expect("coord exists");
            let center = ((p0.0 + p3.0) * 0.5, (p0.1 + p3.1) * 0.5);
            let axis = normalize((p3.0 - p0.0, p3.1 - p0.1));
            let p0_t = (center.0 - 1.5 * axis.0, center.1 - 1.5 * axis.1);
            let p3_t = (center.0 + 1.5 * axis.0, center.1 + 1.5 * axis.1);

            for sign in [1.0f64, -1.0f64] {
                let mut cand = std::collections::HashMap::<usize, (f64, f64)>::new();
                cand.insert(a0, p0_t);
                for k in 1..6 {
                    let rot = rotate_around(p0_t, center, sign * (k as f64) * PI / 3.0);
                    cand.insert(order[k], rot);
                }
                if let Some(&p3c) = cand.get(&a3) {
                    let dx = p3c.0 - p3_t.0;
                    let dy = p3c.1 - p3_t.1;
                    if (dx * dx + dy * dy).sqrt() > 1e-3 {
                        continue;
                    }
                }
                let s = score_for(&cand, &order);
                if s < best_score {
                    best_score = s;
                    best_coords = Some(cand);
                }
            }
        }
        let Some(chosen) = best_coords else { continue };
        for (&a, &p) in &chosen {
            coord_map.insert(a, p);
        }
    }

    for item in local.iter_mut() {
        if let Some(&p) = coord_map.get(&item.0) {
            item.1 = p;
        }
    }
}

fn remove_collisions_bond_flip_like(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
    local: &mut Vec<(usize, (f64, f64))>,
) {
    const COLLISION_THRES: f64 = 0.70;
    const HETEROATOM_COLL_SCALE: f64 = 1.3;
    const MAX_COLL_ITERS: usize = 15;
    const NUM_BONDS_FLIPS: usize = 3;

    let comp_set: std::collections::BTreeSet<usize> = comp.iter().copied().collect();
    let mut pos = std::collections::BTreeMap::<usize, (f64, f64)>::new();
    for &(a, p) in local.iter() {
        pos.insert(a, p);
    }

    let is_ring_bond = |bond_idx: usize| -> bool {
        let b = &mol.bonds[bond_idx];
        let u = b.begin_atom;
        let v = b.end_atom;
        let mut stack = vec![u];
        let mut seen = std::collections::BTreeSet::<usize>::new();
        seen.insert(u);
        while let Some(cur) = stack.pop() {
            for (idx, bond) in mol.bonds.iter().enumerate() {
                if idx == bond_idx {
                    continue;
                }
                let nb = if bond.begin_atom == cur {
                    bond.end_atom
                } else if bond.end_atom == cur {
                    bond.begin_atom
                } else {
                    continue;
                };
                if !comp_set.contains(&nb) {
                    continue;
                }
                if nb == v {
                    return true;
                }
                if seen.insert(nb) {
                    stack.push(nb);
                }
            }
        }
        false
    };

    let shortest_path = |src: usize, dst: usize| -> Option<Vec<usize>> {
        let mut q = std::collections::VecDeque::<usize>::new();
        let mut prev = std::collections::BTreeMap::<usize, usize>::new();
        q.push_back(src);
        prev.insert(src, src);
        while let Some(u) = q.pop_front() {
            if u == dst {
                break;
            }
            for &v in &adjacency[u] {
                if !comp_set.contains(&v) {
                    continue;
                }
                if let std::collections::btree_map::Entry::Vacant(e) = prev.entry(v) {
                    e.insert(u);
                    q.push_back(v);
                }
            }
        }
        if !prev.contains_key(&dst) {
            return None;
        }
        let mut path = vec![dst];
        let mut cur = dst;
        while cur != src {
            cur = prev[&cur];
            path.push(cur);
        }
        path.reverse();
        Some(path)
    };

    let get_rotatable_bonds = |aid1: usize, aid2: usize| -> Vec<usize> {
        let Some(mut path) = shortest_path(aid1, aid2) else {
            return Vec::new();
        };
        if path.len() < 4 {
            return Vec::new();
        }
        path.remove(0);
        path.pop();
        let mut res = Vec::new();
        let mut pid = path[0];
        for aid in path {
            if let Some((bond_idx, _)) = mol.bonds.iter().enumerate().find(|(_, bond)| {
                (bond.begin_atom == pid && bond.end_atom == aid)
                    || (bond.begin_atom == aid && bond.end_atom == pid)
            }) {
                if !is_ring_bond(bond_idx) {
                    res.push(bond_idx);
                }
            }
            pid = aid;
        }
        res
    };

    let reflect_point = |p: (f64, f64), a: (f64, f64), b: (f64, f64)| -> (f64, f64) {
        let vx = b.0 - a.0;
        let vy = b.1 - a.1;
        let v2 = vx * vx + vy * vy;
        if v2 <= 1e-12 {
            return p;
        }
        let t = ((p.0 - a.0) * vx + (p.1 - a.1) * vy) / v2;
        let proj = (a.0 + t * vx, a.1 + t * vy);
        (2.0 * proj.0 - p.0, 2.0 * proj.1 - p.1)
    };

    let find_collisions =
        |pos: &std::collections::BTreeMap<usize, (f64, f64)>| -> (Vec<(usize, usize)>, f64) {
            let mut collisions = Vec::new();
            let mut density = std::collections::BTreeMap::<usize, f64>::new();
            let mut atoms = comp.to_vec();
            atoms.sort_unstable();
            let col_thres2 = COLLISION_THRES * COLLISION_THRES;
            for i in 0..atoms.len() {
                let a = atoms[i];
                let factor_a = if mol.atoms[a].atomic_num != 6 {
                    HETEROATOM_COLL_SCALE
                } else {
                    1.0
                };
                for &b in atoms.iter().take(i) {
                    let factor_b = if mol.atoms[b].atomic_num != 6 {
                        HETEROATOM_COLL_SCALE
                    } else {
                        1.0
                    };
                    let pa = pos[&a];
                    let pb = pos[&b];
                    let dx = pb.0 - pa.0;
                    let dy = pb.1 - pa.1;
                    let d2_raw = dx * dx + dy * dy;
                    let add = if d2_raw > 1.0e-3 {
                        1.0 / d2_raw
                    } else {
                        1000.0
                    };
                    *density.entry(a).or_insert(0.0) += add;
                    *density.entry(b).or_insert(0.0) += add;
                    if d2_raw / (factor_a * factor_b) < col_thres2 {
                        collisions.push((a, b));
                    }
                }
            }
            (collisions, density.values().sum())
        };

    let flip_about_bond = |pos: &mut std::collections::BTreeMap<usize, (f64, f64)>,
                           bond_idx: usize,
                           flip_end: bool| {
        let bond = &mol.bonds[bond_idx];
        let mut beg = bond.begin_atom;
        let mut end = bond.end_atom;
        if !flip_end {
            std::mem::swap(&mut beg, &mut end);
        }
        let mut end_side = Vec::new();
        let mut stack = vec![end];
        let mut seen = std::collections::BTreeSet::<usize>::new();
        seen.insert(end);
        while let Some(cur) = stack.pop() {
            end_side.push(cur);
            for &nb in &adjacency[cur] {
                if !comp_set.contains(&nb) || nb == beg {
                    continue;
                }
                if seen.insert(nb) {
                    stack.push(nb);
                }
            }
        }
        let end_side_set: std::collections::BTreeSet<usize> = end_side.iter().copied().collect();
        let mut end_side_flip = true;
        if comp.len().saturating_sub(end_side.len()) < end_side.len() {
            end_side_flip = false;
        }
        let beg_loc = pos[&beg];
        let end_loc = pos[&end];
        for &aid in comp {
            let in_end_side = end_side_set.contains(&aid);
            if end_side_flip ^ !in_end_side {
                let p = pos[&aid];
                pos.insert(aid, reflect_point(p, beg_loc, end_loc));
            }
        }
    };

    let (mut colls, _) = find_collisions(&pos);
    let mut done_bonds = std::collections::BTreeMap::<usize, usize>::new();
    let mut iter = 0usize;
    while iter < MAX_COLL_ITERS && !colls.is_empty() {
        let ncols = colls.len();
        let (aid1, aid2) = colls[0];
        let rot_bonds = get_rotatable_bonds(aid1, aid2);
        let (_, prev_density) = find_collisions(&pos);
        for ri in rot_bonds {
            let done = *done_bonds.get(&ri).unwrap_or(&0);
            if done >= NUM_BONDS_FLIPS {
                continue;
            }
            done_bonds.insert(ri, done + 1);

            let before = pos.clone();
            flip_about_bond(&mut pos, ri, true);
            let (new_colls, new_density) = find_collisions(&pos);
            if new_colls.len() < ncols {
                done_bonds.insert(ri, NUM_BONDS_FLIPS);
                colls = new_colls;
                break;
            } else if new_colls.len() == ncols && new_density < prev_density {
                colls = new_colls;
                break;
            }

            pos = before;
            flip_about_bond(&mut pos, ri, false);
            let (new_colls, new_density) = find_collisions(&pos);
            if new_colls.len() < ncols {
                done_bonds.insert(ri, NUM_BONDS_FLIPS);
                colls = new_colls;
                break;
            } else if new_colls.len() == ncols && new_density < prev_density {
                colls = new_colls;
                break;
            } else {
                flip_about_bond(&mut pos, ri, false);
                let (reset_colls, _) = find_collisions(&pos);
                colls = reset_colls;
            }
        }
        iter += 1;
    }

    for item in local.iter_mut() {
        if let Some(&p) = pos.get(&item.0) {
            item.1 = p;
        }
    }
}

#[allow(dead_code)]
fn remove_collisions_bond_flip_non_source_dead(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
    local: &mut Vec<(usize, (f64, f64))>,
) {
    let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();
    let mut pos = std::collections::HashMap::<usize, (f64, f64)>::new();
    for &(a, p) in local.iter() {
        pos.insert(a, p);
    }

    let mut bond_pairs = std::collections::HashSet::<(usize, usize)>::new();
    for b in &mol.bonds {
        if !comp_set.contains(&b.begin_atom) || !comp_set.contains(&b.end_atom) {
            continue;
        }
        let (a, c) = if b.begin_atom <= b.end_atom {
            (b.begin_atom, b.end_atom)
        } else {
            (b.end_atom, b.begin_atom)
        };
        bond_pairs.insert((a, c));
    }

    let ring_edge = |u: usize, v: usize| -> bool {
        let mut stack = vec![u];
        let mut seen = std::collections::HashSet::<usize>::new();
        seen.insert(u);
        while let Some(cur) = stack.pop() {
            for &nb in &adjacency[cur] {
                if !comp_set.contains(&nb) {
                    continue;
                }
                if (cur == u && nb == v) || (cur == v && nb == u) {
                    continue;
                }
                if nb == v {
                    return true;
                }
                if seen.insert(nb) {
                    stack.push(nb);
                }
            }
        }
        false
    };

    let score = |pos: &std::collections::HashMap<usize, (f64, f64)>| -> f64 {
        let mut atoms = comp.to_vec();
        atoms.sort_unstable();
        let mut s = 0.0f64;
        for i in 0..atoms.len() {
            for j in (i + 1)..atoms.len() {
                let a = atoms[i];
                let b = atoms[j];
                let pair = if a <= b { (a, b) } else { (b, a) };
                if bond_pairs.contains(&pair) {
                    continue;
                }
                let pa = pos[&a];
                let pb = pos[&b];
                let dx = pa.0 - pb.0;
                let dy = pa.1 - pb.1;
                let d = (dx * dx + dy * dy).sqrt();
                let dd = d.max(1e-3);
                s += 1.0 / (dd * dd);
            }
        }
        for b in &mol.bonds {
            if !comp_set.contains(&b.begin_atom) || !comp_set.contains(&b.end_atom) {
                continue;
            }
            let pa = pos[&b.begin_atom];
            let pb = pos[&b.end_atom];
            let dx = pa.0 - pb.0;
            let dy = pa.1 - pb.1;
            let d = (dx * dx + dy * dy).sqrt();
            s += (d - 1.5).abs() * 0.2;
        }
        s
    };

    let reflect_point = |p: (f64, f64), a: (f64, f64), b: (f64, f64)| -> (f64, f64) {
        let vx = b.0 - a.0;
        let vy = b.1 - a.1;
        let v2 = vx * vx + vy * vy;
        if v2 <= 1e-12 {
            return p;
        }
        let t = ((p.0 - a.0) * vx + (p.1 - a.1) * vy) / v2;
        let proj = (a.0 + t * vx, a.1 + t * vy);
        (2.0 * proj.0 - p.0, 2.0 * proj.1 - p.1)
    };

    let find_collisions =
        |pos: &std::collections::HashMap<usize, (f64, f64)>| -> Vec<(usize, usize)> {
            let mut out = Vec::new();
            let mut atoms = comp.to_vec();
            atoms.sort_unstable();
            for i in 0..atoms.len() {
                let a = atoms[i];
                let fa = if mol.atoms[a].atomic_num != 6 {
                    1.3
                } else {
                    1.0
                };
                for &b in atoms.iter().take(i) {
                    let fb = if mol.atoms[b].atomic_num != 6 {
                        1.3
                    } else {
                        1.0
                    };
                    let pa = pos[&a];
                    let pb = pos[&b];
                    let dx = pa.0 - pb.0;
                    let dy = pa.1 - pb.1;
                    let mut d2 = dx * dx + dy * dy;
                    d2 /= fa * fb;
                    if d2 < 0.70f64 * 0.70f64 {
                        out.push((a, b));
                    }
                }
            }
            out
        };

    let shortest_path = |src: usize, dst: usize| -> Option<Vec<usize>> {
        let mut q = std::collections::VecDeque::<usize>::new();
        let mut prev = std::collections::HashMap::<usize, usize>::new();
        q.push_back(src);
        prev.insert(src, src);
        while let Some(u) = q.pop_front() {
            if u == dst {
                break;
            }
            for &v in &adjacency[u] {
                if !comp_set.contains(&v) {
                    continue;
                }
                if let std::collections::hash_map::Entry::Vacant(e) = prev.entry(v) {
                    e.insert(u);
                    q.push_back(v);
                }
            }
        }
        if !prev.contains_key(&dst) {
            return None;
        }
        let mut path = vec![dst];
        let mut cur = dst;
        while cur != src {
            cur = prev[&cur];
            path.push(cur);
        }
        path.reverse();
        Some(path)
    };

    let mut done_bonds = std::collections::HashMap::<(usize, usize), usize>::new();
    let mut iter = 0usize;
    while iter < 15 {
        let colls = find_collisions(&pos);
        if colls.is_empty() {
            break;
        }
        let (a, b) = colls[0];
        let Some(path) = shortest_path(a, b) else {
            break;
        };
        if path.len() < 4 {
            iter += 1;
            continue;
        }
        let internal = &path[1..path.len() - 1];
        let mut rot_bonds: Vec<(usize, usize)> = Vec::new();
        let mut pid = internal[0];
        for &aid in internal.iter() {
            if aid == pid {
                continue;
            }
            let mut u = pid;
            let mut v = aid;
            if u > v {
                std::mem::swap(&mut u, &mut v);
            }
            let bond = mol.bonds.iter().find(|bd| {
                (bd.begin_atom == u && bd.end_atom == v) || (bd.begin_atom == v && bd.end_atom == u)
            });
            let Some(bd) = bond else {
                pid = aid;
                continue;
            };
            if !matches!(bd.order, BondOrder::Single) {
                pid = aid;
                continue;
            }
            if ring_edge(u, v) {
                pid = aid;
                continue;
            }
            rot_bonds.push((u, v));
            pid = aid;
        }

        let mut accepted = false;
        let ncols = colls.len();
        let prev_density = score(&pos);
        for &(u, v) in &rot_bonds {
            let cnt = *done_bonds.get(&(u, v)).unwrap_or(&0);
            if cnt >= 3 {
                continue;
            }
            done_bonds.insert((u, v), cnt + 1);

            let try_flip = |from: usize,
                            to: usize,
                            pos_in: &std::collections::HashMap<usize, (f64, f64)>|
             -> Option<std::collections::HashMap<usize, (f64, f64)>> {
                let mut side = std::collections::HashSet::<usize>::new();
                let mut stack = vec![to];
                side.insert(to);
                while let Some(cur) = stack.pop() {
                    for &nb in &adjacency[cur] {
                        if !comp_set.contains(&nb) {
                            continue;
                        }
                        if cur == to && nb == from {
                            continue;
                        }
                        if side.insert(nb) {
                            stack.push(nb);
                        }
                    }
                }
                if side.is_empty() || side.len() == comp.len() {
                    return None;
                }
                let pu = pos_in[&from];
                let pv = pos_in[&to];
                let mut trial = pos_in.clone();
                for &x in &side {
                    if x == from {
                        continue;
                    }
                    let p = trial[&x];
                    trial.insert(x, reflect_point(p, pu, pv));
                }
                Some(trial)
            };

            for (from, to) in [(u, v), (v, u)] {
                let Some(trial) = try_flip(from, to, &pos) else {
                    continue;
                };
                let c2 = find_collisions(&trial).len();
                let den2 = score(&trial);
                if c2 < ncols || (c2 == ncols && den2 < prev_density) {
                    pos = trial;
                    if c2 < ncols {
                        done_bonds.insert((u, v), 3);
                    }
                    accepted = true;
                    break;
                }
            }
            if accepted {
                break;
            }
        }
        if !accepted {
            break;
        }
        iter += 1;
    }

    for item in local.iter_mut() {
        if let Some(&p) = pos.get(&item.0) {
            item.1 = p;
        }
    }
}

fn remove_collisions_open_angles_like(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
    local: &mut Vec<(usize, (f64, f64))>,
) {
    const COLLISION_THRES: f64 = 0.70;
    const HETEROATOM_COLL_SCALE: f64 = 1.3;
    const ANGLE_OPEN: f64 = 0.1222;

    if comp.len() < 3 {
        return;
    }
    let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();
    let mut pos = std::collections::HashMap::<usize, (f64, f64)>::new();
    for &(a, p) in local.iter() {
        pos.insert(a, p);
    }

    let n = mol.atoms.len();
    let mut dmat = vec![vec![usize::MAX; n]; n];
    for &s in comp {
        let mut q = std::collections::VecDeque::<usize>::new();
        dmat[s][s] = 0;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            let du = dmat[s][u];
            for &v in &adjacency[u] {
                if !comp_set.contains(&v) {
                    continue;
                }
                if dmat[s][v] == usize::MAX {
                    dmat[s][v] = du + 1;
                    q.push_back(v);
                }
            }
        }
    }

    let find_collisions =
        |pos: &std::collections::HashMap<usize, (f64, f64)>| -> Vec<(usize, usize)> {
            let mut out = Vec::new();
            let mut atoms = comp.to_vec();
            atoms.sort_unstable();
            for i in 0..atoms.len() {
                let a = atoms[i];
                let fa = if mol.atoms[a].atomic_num != 6 {
                    HETEROATOM_COLL_SCALE
                } else {
                    1.0
                };
                for &b in atoms.iter().take(i) {
                    let fb = if mol.atoms[b].atomic_num != 6 {
                        HETEROATOM_COLL_SCALE
                    } else {
                        1.0
                    };
                    let pa = pos[&a];
                    let pb = pos[&b];
                    let dx = pa.0 - pb.0;
                    let dy = pa.1 - pb.1;
                    let mut d2 = dx * dx + dy * dy;
                    d2 /= fa * fb;
                    if d2 < COLLISION_THRES * COLLISION_THRES {
                        out.push((a, b));
                    }
                }
            }
            out
        };

    let find_deg1_neighbor = |aid: usize| -> Option<usize> {
        let nbs: Vec<usize> = adjacency[aid]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .collect();
        if nbs.len() == 1 { Some(nbs[0]) } else { None }
    };

    let find_closest_neighbor = |aid1: usize, aid2: usize| -> Option<usize> {
        adjacency[aid2]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .min_by_key(|n| dmat[aid1][*n])
    };

    let rotate_point = |p: (f64, f64), c: (f64, f64), ang: f64| -> (f64, f64) {
        let v = (p.0 - c.0, p.1 - c.1);
        let r = rotate(v, ang);
        (c.0 + r.0, c.1 + r.1)
    };

    for (aid1, aid2) in find_collisions(&pos) {
        let deg1 = adjacency[aid1]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .count();
        let deg2 = adjacency[aid2]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .count();
        if deg1 > 1 && deg2 > 1 {
            continue;
        }

        let (aid_a, aid_b, kind) = if deg1 == 1 && deg2 == 1 {
            (
                find_deg1_neighbor(aid1),
                find_deg1_neighbor(aid2),
                1, // both endpoints
            )
        } else if deg1 == 1 {
            (
                find_deg1_neighbor(aid1),
                find_closest_neighbor(find_deg1_neighbor(aid1).unwrap_or(aid1), aid2),
                2,
            )
        } else {
            (
                find_closest_neighbor(find_deg1_neighbor(aid2).unwrap_or(aid2), aid1),
                find_deg1_neighbor(aid2),
                3,
            )
        };
        let (Some(aid_a), Some(aid_b)) = (aid_a, aid_b) else {
            continue;
        };

        let v2 = (pos[&aid1].0 - pos[&aid_a].0, pos[&aid1].1 - pos[&aid_a].1);
        let v1 = (pos[&aid_b].0 - pos[&aid_a].0, pos[&aid_b].1 - pos[&aid_a].1);
        let cross = v1.0 * v2.1 - v1.1 * v2.0;
        match kind {
            1 => {
                let mut angle = ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let p1 = rotate_point(pos[&aid1], pos[&aid_a], angle);
                let p2 = rotate_point(pos[&aid2], pos[&aid_b], -angle);
                pos.insert(aid1, p1);
                pos.insert(aid2, p2);
            }
            2 => {
                let mut angle = 2.0 * ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let p1 = rotate_point(pos[&aid1], pos[&aid_a], angle);
                pos.insert(aid1, p1);
            }
            3 => {
                let mut angle = -2.0 * ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let p2 = rotate_point(pos[&aid2], pos[&aid_b], angle);
                pos.insert(aid2, p2);
            }
            _ => {}
        }
    }

    for item in local.iter_mut() {
        if let Some(&p) = pos.get(&item.0) {
            item.1 = p;
        }
    }
}

fn optimize_aromatic_ring_reflections(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
    local: &mut Vec<(usize, (f64, f64))>,
) {
    let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();
    let mut pos = std::collections::HashMap::<usize, (f64, f64)>::new();
    for &(a, p) in local.iter() {
        pos.insert(a, p);
    }

    let mut aromatic_adj = std::collections::HashMap::<usize, Vec<usize>>::new();
    for &a in comp {
        aromatic_adj.insert(a, Vec::new());
    }
    for b in &mol.bonds {
        if !comp_set.contains(&b.begin_atom) || !comp_set.contains(&b.end_atom) {
            continue;
        }
        if !matches!(b.order, BondOrder::Aromatic) {
            continue;
        }
        aromatic_adj
            .entry(b.begin_atom)
            .or_default()
            .push(b.end_atom);
        aromatic_adj
            .entry(b.end_atom)
            .or_default()
            .push(b.begin_atom);
    }
    for nbs in aromatic_adj.values_mut() {
        nbs.sort_unstable();
        nbs.dedup();
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

    let mut seen = std::collections::HashSet::<Vec<usize>>::new();
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    let mut comp_sorted = comp.to_vec();
    comp_sorted.sort_unstable();
    for &start in &comp_sorted {
        let mut path = vec![start];
        let mut used = std::collections::HashSet::<usize>::new();
        used.insert(start);
        fn dfs6(
            start: usize,
            cur: usize,
            depth: usize,
            path: &mut Vec<usize>,
            used: &mut std::collections::HashSet<usize>,
            aromatic_adj: &std::collections::HashMap<usize, Vec<usize>>,
            seen: &mut std::collections::HashSet<Vec<usize>>,
            cycles: &mut Vec<Vec<usize>>,
        ) {
            if depth == 6 {
                if aromatic_adj
                    .get(&cur)
                    .is_some_and(|nbs| nbs.contains(&start))
                {
                    let c = canonical_cycle(path.clone());
                    if seen.insert(c.clone()) {
                        cycles.push(c);
                    }
                }
                return;
            }
            let mut nbs = aromatic_adj.get(&cur).cloned().unwrap_or_default();
            nbs.sort_unstable();
            for nb in nbs {
                if used.contains(&nb) || nb < start {
                    continue;
                }
                used.insert(nb);
                path.push(nb);
                dfs6(start, nb, depth + 1, path, used, aromatic_adj, seen, cycles);
                path.pop();
                used.remove(&nb);
            }
        }
        dfs6(
            start,
            start,
            1,
            &mut path,
            &mut used,
            &aromatic_adj,
            &mut seen,
            &mut cycles,
        );
    }
    if cycles.is_empty() {
        return;
    }

    let mut atom_cycle_count = std::collections::HashMap::<usize, usize>::new();
    for cyc in &cycles {
        for &a in cyc {
            *atom_cycle_count.entry(a).or_insert(0) += 1;
        }
    }

    let score = |pos: &std::collections::HashMap<usize, (f64, f64)>| -> f64 {
        let mut atoms = comp.to_vec();
        atoms.sort_unstable();
        let mut s = 0.0f64;
        for i in 0..atoms.len() {
            for j in (i + 1)..atoms.len() {
                let a = atoms[i];
                let b = atoms[j];
                let pa = pos[&a];
                let pb = pos[&b];
                let dx = pa.0 - pb.0;
                let dy = pa.1 - pb.1;
                let d = (dx * dx + dy * dy).sqrt().max(1e-3);
                s += 1.0 / (d * d);
            }
        }
        s
    };

    let reflect_point = |p: (f64, f64), a: (f64, f64), b: (f64, f64)| -> (f64, f64) {
        let vx = b.0 - a.0;
        let vy = b.1 - a.1;
        let v2 = vx * vx + vy * vy;
        if v2 <= 1e-12 {
            return p;
        }
        let t = ((p.0 - a.0) * vx + (p.1 - a.1) * vy) / v2;
        let proj = (a.0 + t * vx, a.1 + t * vy);
        (2.0 * proj.0 - p.0, 2.0 * proj.1 - p.1)
    };

    let mut best = score(&pos);
    for cyc in cycles {
        if cyc.len() != 6 {
            continue;
        }
        if cyc
            .iter()
            .any(|a| atom_cycle_count.get(a).copied().unwrap_or(0) != 1)
        {
            continue;
        }
        let ext_deg = |a: usize| -> usize {
            adjacency[a]
                .iter()
                .filter(|n| comp_set.contains(n) && !cyc.contains(n))
                .count()
        };
        let mut opposite_pairs: Vec<(usize, usize)> = Vec::new();
        for i in 0..6 {
            let j = (i + 3) % 6;
            if ext_deg(cyc[i]) > 0 && ext_deg(cyc[j]) > 0 {
                opposite_pairs.push((i, j));
            }
        }
        for (i0, i3) in opposite_pairs {
            let a0 = cyc[i0];
            let a3 = cyc[i3];
            let Some(&p0) = pos.get(&a0) else { continue };
            let Some(&p3) = pos.get(&a3) else { continue };

            let mut trial = pos.clone();
            for &a in &cyc {
                if a == a0 || a == a3 {
                    continue;
                }
                let pp = trial[&a];
                trial.insert(a, reflect_point(pp, p0, p3));
            }
            let s = score(&trial);
            if s + 1e-8 < best {
                pos = trial;
                best = s;
            }
        }
    }

    for item in local.iter_mut() {
        if let Some(&p) = pos.get(&item.0) {
            item.1 = p;
        }
    }
}

fn place_fused_hex_component(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
) -> Option<Vec<(usize, (f64, f64))>> {
    let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();
    let mut deg_in_comp = vec![0usize; mol.atoms.len()];
    for &a in comp {
        deg_in_comp[a] = adjacency[a].iter().filter(|n| comp_set.contains(n)).count();
    }

    // Ring core by leaf peeling (same approach as unicyclic).
    let mut removed = std::collections::HashSet::<usize>::new();
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
    if ring_atoms.len() < 6 {
        return None;
    }
    let ring_set: std::collections::HashSet<usize> = ring_atoms.iter().copied().collect();

    // RDKit feeds EmbeddedFrag::embedFusedRings() from symmetrizeSSSR(); for
    // this strict subset we use the smallest independent cycles from the ring
    // core as the corresponding fused ring set.
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

    let mut seen = std::collections::HashSet::<Vec<usize>>::new();
    let mut all_cycles: Vec<Vec<usize>> = Vec::new();
    let mut ring_atoms_sorted = ring_atoms.clone();
    ring_atoms_sorted.sort_unstable();
    for &start in &ring_atoms_sorted {
        let mut path = vec![start];
        let mut used = std::collections::HashSet::<usize>::new();
        used.insert(start);
        fn dfs_cycles(
            start: usize,
            cur: usize,
            path: &mut Vec<usize>,
            used: &mut std::collections::HashSet<usize>,
            adjacency: &[Vec<usize>],
            ring_set: &std::collections::HashSet<usize>,
            seen: &mut std::collections::HashSet<Vec<usize>>,
            cycles: &mut Vec<Vec<usize>>,
        ) {
            if path.len() >= 3 && adjacency[cur].contains(&start) {
                let c = canonical_cycle(path.clone());
                if seen.insert(c.clone()) {
                    cycles.push(c);
                }
            }
            if path.len() == 8 {
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
            adjacency,
            &ring_set,
            &mut seen,
            &mut all_cycles,
        );
    }
    if all_cycles.is_empty() {
        return None;
    }

    all_cycles.retain(|cyc| {
        (0..cyc.len()).all(|i| {
            let a = cyc[i];
            let b = cyc[(i + 1) % cyc.len()];
            adjacency[a].contains(&b) && ring_set.contains(&a) && ring_set.contains(&b)
        })
    });
    if all_cycles.is_empty() {
        return None;
    }

    all_cycles.sort_by_key(|cyc| (cyc.len(), cyc.clone()));
    let ring_bond_count = mol
        .bonds
        .iter()
        .filter(|b| ring_set.contains(&b.begin_atom) && ring_set.contains(&b.end_atom))
        .count();
    let target_cycle_count = ring_bond_count + 1 - ring_atoms.len();
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    let mut covered_edges = std::collections::HashSet::<(usize, usize)>::new();
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
    if cycles.is_empty() {
        return None;
    }
    // Detect whether ring cycles are fused (share an edge).
    let mut has_shared_edge = false;
    for i in 0..cycles.len() {
        for j in (i + 1)..cycles.len() {
            let ci = &cycles[i];
            let cj = &cycles[j];
            'edge_scan: for k in 0..ci.len() {
                let a = ci[k];
                let b = ci[(k + 1) % ci.len()];
                for l in 0..cj.len() {
                    let x = cj[l];
                    let y = cj[(l + 1) % cj.len()];
                    if (a == x && b == y) || (a == y && b == x) {
                        has_shared_edge = true;
                        break 'edge_scan;
                    }
                }
            }
            if has_shared_edge {
                break;
            }
        }
        if has_shared_edge {
            break;
        }
    }
    let mut ring_centers = vec![(0.0f64, 0.0f64); cycles.len()];
    let mut ring_done = vec![false; cycles.len()];
    let mut coords: std::collections::HashMap<usize, (f64, f64)> = std::collections::HashMap::new();

    if !has_shared_edge {
        // route non-fused multi-ring systems to the dedicated seeded-growth path
        return place_multiring_nonfused_component_rdkit_like(mol, comp, adjacency, &deg_in_comp);
    }
    if has_shared_edge {
        // Pick first ring using a simplified RDKit heuristic: fewer substituents, then larger size.
        let mut first_id = 0usize;
        let mut best = (usize::MAX, 0usize, usize::MAX);
        for (i, cyc) in cycles.iter().enumerate() {
            let subs = cyc
                .iter()
                .copied()
                .filter(|&a| {
                    adjacency[a]
                        .iter()
                        .any(|n| comp_set.contains(n) && !ring_set.contains(n))
                })
                .count();
            let key = (subs, usize::MAX - cyc.len(), i);
            if key < best {
                best = key;
                first_id = i;
            }
        }

        // embedRing from RDKit DepictUtils: regular polygon around origin.
        let n0 = cycles[first_id].len();
        let ang = 2.0 * PI / n0 as f64;
        let radius = 1.5f64 / (2.0 * (PI / n0 as f64).sin());
        for (k, &a) in cycles[first_id].iter().enumerate() {
            coords.insert(
                a,
                (
                    radius * (k as f64 * ang).cos(),
                    radius * (k as f64 * ang).sin(),
                ),
            );
        }
        ring_done[first_id] = true;
        let c0 = cycles[first_id]
            .iter()
            .map(|a| coords[a])
            .fold((0.0, 0.0), |acc, p| (acc.0 + p.0, acc.1 + p.1));
        ring_centers[first_id] = (c0.0 / n0 as f64, c0.1 / n0 as f64);

        // BFS-like fused placement by shared edge.
        loop {
            let mut progressed = false;
            for rid in 0..cycles.len() {
                if ring_done[rid] {
                    continue;
                }
                let ring = &cycles[rid];
                let common: Vec<usize> = ring
                    .iter()
                    .copied()
                    .filter(|a| coords.contains_key(a))
                    .collect();
                if common.len() >= 3 {
                    let n = ring.len();
                    let radius = 1.5f64 / (2.0 * (PI / n as f64).sin());
                    let local: std::collections::HashMap<usize, (f64, f64)> = ring
                        .iter()
                        .enumerate()
                        .map(|(k, &a)| {
                            let theta = 2.0 * PI * k as f64 / n as f64;
                            (a, (radius * theta.cos(), radius * theta.sin()))
                        })
                        .collect();
                    let common_set: std::collections::HashSet<usize> =
                        common.iter().copied().collect();
                    let mut endpoints: Vec<usize> = common
                        .iter()
                        .copied()
                        .filter(|&a| {
                            adjacency[a]
                                .iter()
                                .filter(|n| common_set.contains(n))
                                .count()
                                == 1
                        })
                        .collect();
                    endpoints.sort_unstable();
                    let (aid1, aid2) = if endpoints.len() >= 2 {
                        (endpoints[0], endpoints[1])
                    } else {
                        (common[0], *common.last().expect("common len checked"))
                    };
                    let mut best_candidate: Option<std::collections::HashMap<usize, (f64, f64)>> =
                        None;
                    let mut best_score = f64::INFINITY;
                    for (pin1, pin2) in [(aid1, aid2), (aid2, aid1)] {
                        let src1 = local[&pin1];
                        let src2 = local[&pin2];
                        let dst1 = coords[&pin1];
                        let dst2 = coords[&pin2];
                        for reflect in [false, true] {
                            let map_point = |p: (f64, f64)| -> Option<(f64, f64)> {
                                let p = if reflect { (p.0, -p.1) } else { p };
                                let s1 = if reflect { (src1.0, -src1.1) } else { src1 };
                                let s2 = if reflect { (src2.0, -src2.1) } else { src2 };
                                let sv = (s2.0 - s1.0, s2.1 - s1.1);
                                let tv = (dst2.0 - dst1.0, dst2.1 - dst1.1);
                                let sl = norm(sv);
                                let tl = norm(tv);
                                if sl <= 1e-12 || tl <= 1e-12 {
                                    return None;
                                }
                                let scale = tl / sl;
                                let cos_t = (sv.0 * tv.0 + sv.1 * tv.1) / (sl * tl);
                                let sin_t = (sv.0 * tv.1 - sv.1 * tv.0) / (sl * tl);
                                let rel = (p.0 - s1.0, p.1 - s1.1);
                                Some((
                                    dst1.0 + scale * (cos_t * rel.0 - sin_t * rel.1),
                                    dst1.1 + scale * (sin_t * rel.0 + cos_t * rel.1),
                                ))
                            };

                            let mut candidate =
                                std::collections::HashMap::<usize, (f64, f64)>::new();
                            let mut score = 0.0;
                            let mut ok = true;
                            for &a in ring {
                                let Some(p) = map_point(local[&a]) else {
                                    ok = false;
                                    break;
                                };
                                if let Some(existing) = coords.get(&a) {
                                    let dx = existing.0 - p.0;
                                    let dy = existing.1 - p.1;
                                    let d2 = dx * dx + dy * dy;
                                    if d2 > 1e-4 {
                                        ok = false;
                                        break;
                                    }
                                    score += d2;
                                } else {
                                    for (&existing_atom, &existing) in &coords {
                                        if ring.contains(&existing_atom) {
                                            continue;
                                        }
                                        let dx = existing.0 - p.0;
                                        let dy = existing.1 - p.1;
                                        let d2 = dx * dx + dy * dy;
                                        if d2 < 1e-8 {
                                            score += 1.0e9;
                                        } else {
                                            score += 1.0 / d2;
                                        }
                                    }
                                }
                                candidate.insert(a, p);
                            }
                            if ok && score < best_score {
                                best_score = score;
                                best_candidate = Some(candidate);
                            }
                        }
                    }
                    if let Some(chosen) = best_candidate {
                        for (&k, &v) in &chosen {
                            coords.entry(k).or_insert(v);
                        }
                        let c = ring
                            .iter()
                            .map(|a| coords[a])
                            .fold((0.0, 0.0), |acc, p| (acc.0 + p.0, acc.1 + p.1));
                        ring_centers[rid] = (c.0 / ring.len() as f64, c.1 / ring.len() as f64);
                        ring_done[rid] = true;
                        progressed = true;
                        continue;
                    }
                }
                // find a done ring sharing a bond (2 adjacent atoms)
                let mut anchor: Option<(usize, usize, usize)> = None; // (done_ring, a, b)
                for did in 0..cycles.len() {
                    if !ring_done[did] {
                        continue;
                    }
                    let done = &cycles[did];
                    for i in 0..ring.len() {
                        let a = ring[i];
                        let b = ring[(i + 1) % ring.len()];
                        if done
                            .iter()
                            .enumerate()
                            .any(|(j, &x)| x == a && done[(j + 1) % done.len()] == b)
                            || done
                                .iter()
                                .enumerate()
                                .any(|(j, &x)| x == b && done[(j + 1) % done.len()] == a)
                        {
                            if coords.contains_key(&a) && coords.contains_key(&b) {
                                anchor = Some((did, a, b));
                                break;
                            }
                        }
                    }
                    if anchor.is_some() {
                        break;
                    }
                }
                let Some((did, a, b)) = anchor else { continue };
                let pa = *coords.get(&a)?;
                let pb = *coords.get(&b)?;
                let shared_mid = ((pa.0 + pb.0) * 0.5, (pa.1 + pb.1) * 0.5);
                let src_center = ring_centers[did];
                let side_src = (pb.0 - pa.0) * (src_center.1 - shared_mid.1)
                    - (pb.1 - pa.1) * (src_center.0 - shared_mid.0);

                let try_build =
                    |reverse: bool,
                     sign: f64|
                     -> Option<std::collections::HashMap<usize, (f64, f64)>> {
                        let mut order = ring.clone();
                        if reverse {
                            order.reverse();
                        }
                        let n = order.len();
                        let mut i0 = None;
                        for i in 0..n {
                            if order[i] == a && order[(i + 1) % n] == b {
                                i0 = Some(i);
                                break;
                            }
                        }
                        let i0 = i0?;
                        let mut local = std::collections::HashMap::<usize, (f64, f64)>::new();
                        local.insert(a, pa);
                        local.insert(b, pb);
                        let turn = sign * (2.0 * PI / n as f64);
                        let mut prev = pb;
                        let mut dir = (pb.0 - pa.0, pb.1 - pa.1);
                        for step in 2..n {
                            dir = rotate(dir, turn);
                            let idx = order[(i0 + step) % n];
                            let p = (prev.0 + dir.0, prev.1 + dir.1);
                            local.insert(idx, p);
                            prev = p;
                        }
                        Some(local)
                    };

                let mut best_candidate: Option<std::collections::HashMap<usize, (f64, f64)>> = None;
                let mut best_score = f64::INFINITY;
                for reverse in [false, true] {
                    for sign in [1.0, -1.0] {
                        let Some(cand) = try_build(reverse, sign) else {
                            continue;
                        };
                        let c = ring
                            .iter()
                            .map(|a| cand[a])
                            .fold((0.0, 0.0), |acc, p| (acc.0 + p.0, acc.1 + p.1));
                        let center = (c.0 / ring.len() as f64, c.1 / ring.len() as f64);
                        let side = (pb.0 - pa.0) * (center.1 - shared_mid.1)
                            - (pb.1 - pa.1) * (center.0 - shared_mid.0);
                        if side * side_src >= 0.0 {
                            continue;
                        }
                        let mut score = 0.0;
                        let mut ok = true;
                        for &ra in ring {
                            if let Some(&pp) = coords.get(&ra) {
                                let cp = cand[&ra];
                                let dx = pp.0 - cp.0;
                                let dy = pp.1 - cp.1;
                                let d2 = dx * dx + dy * dy;
                                if d2 > 1e-4 {
                                    ok = false;
                                    break;
                                }
                                score += d2;
                            }
                        }
                        if ok && score < best_score {
                            best_score = score;
                            best_candidate = Some(cand);
                        }
                    }
                }
                let Some(chosen) = best_candidate else {
                    continue;
                };
                for (&k, &v) in &chosen {
                    coords.entry(k).or_insert(v);
                }
                let c = ring
                    .iter()
                    .map(|a| coords[a])
                    .fold((0.0, 0.0), |acc, p| (acc.0 + p.0, acc.1 + p.1));
                ring_centers[rid] = (c.0 / ring.len() as f64, c.1 / ring.len() as f64);
                ring_done[rid] = true;
                progressed = true;
            }
            if !progressed {
                break;
            }
        }
        if ring_done.iter().any(|d| !*d) {
            return None;
        }
    } else {
        // Non-fused multi-ring system: place each ring independently as a
        // template, then continue with substituent growth.
        let mut x_shift = 0.0f64;
        for (rid, ring) in cycles.iter().enumerate() {
            let n0 = ring.len();
            let ang = 2.0 * PI / n0 as f64;
            let radius = 1.5f64 / (2.0 * (PI / n0 as f64).sin());
            for (k, &a) in ring.iter().enumerate() {
                coords.entry(a).or_insert((
                    x_shift + radius * (k as f64 * ang).cos(),
                    radius * (k as f64 * ang).sin(),
                ));
            }
            let c0 = ring
                .iter()
                .map(|a| coords[a])
                .fold((0.0, 0.0), |acc, p| (acc.0 + p.0, acc.1 + p.1));
            ring_centers[rid] = (c0.0 / n0 as f64, c0.1 / n0 as f64);
            ring_done[rid] = true;
            x_shift += 2.0 * radius + 2.5;
        }
    }
    if ring_atoms.iter().any(|a| !coords.contains_key(a)) {
        return None;
    }

    // Place non-ring substituents outward from nearest ring atom.
    let core_center = ring_atoms
        .iter()
        .map(|a| coords[a])
        .fold((0.0, 0.0), |acc, p| (acc.0 + p.0, acc.1 + p.1));
    let core_center = (
        core_center.0 / ring_atoms.len() as f64,
        core_center.1 / ring_atoms.len() as f64,
    );

    let mut stack: Vec<(usize, usize)> = Vec::new(); // (parent, child)
    for &ra in &ring_atoms {
        let mut outs: Vec<usize> = adjacency[ra]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n) && !ring_set.contains(n))
            .collect();
        outs.sort_unstable();
        let pa = coords[&ra];
        let ring_nbs: Vec<usize> = adjacency[ra]
            .iter()
            .copied()
            .filter(|n| ring_set.contains(n))
            .collect();
        let mut base = if ring_nbs.len() >= 2 {
            let p1 = coords[&ring_nbs[0]];
            let p2 = coords[&ring_nbs[1]];
            let v1 = normalize((p1.0 - pa.0, p1.1 - pa.1));
            let v2 = normalize((p2.0 - pa.0, p2.1 - pa.1));
            let bis = normalize((v1.0 + v2.0, v1.1 + v2.1));
            (-bis.0, -bis.1)
        } else {
            normalize((pa.0 - core_center.0, pa.1 - core_center.1))
        };
        if norm(base) < 1e-8 {
            base = (1.0, 0.0);
        }
        for (i, ch) in outs.into_iter().enumerate() {
            if coords.contains_key(&ch) {
                continue;
            }
            let delta = rotate(base, (i as f64) * (PI / 6.0));
            let p = (pa.0 + 1.5 * delta.0, pa.1 + 1.5 * delta.1);
            coords.insert(ch, p);
            stack.push((ra, ch));
        }
    }
    while let Some((parent, cur)) = stack.pop() {
        let pc = coords[&cur];
        let pp = coords[&parent];
        let in_dir = normalize((pc.0 - pp.0, pc.1 - pp.1));
        let mut children: Vec<usize> = adjacency[cur]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n) && *n != parent && !coords.contains_key(n))
            .collect();
        children.sort_unstable();
        let angles: Vec<f64> = match children.len() {
            0 => vec![],
            1 => vec![PI / 3.0],
            2 => vec![-PI / 3.0, PI / 3.0],
            3 => vec![-PI / 2.0, 0.0, PI / 2.0],
            m => {
                let step = PI / (m as f64 - 1.0);
                (0..m).map(|i| -PI / 2.0 + step * i as f64).collect()
            }
        };
        for (i, ch) in children.into_iter().enumerate() {
            let d = rotate(in_dir, angles[i]);
            let p = (pc.0 + 1.5 * d.0, pc.1 + 1.5 * d.1);
            coords.insert(ch, p);
            stack.push((cur, ch));
        }
    }

    if comp.iter().any(|a| !coords.contains_key(a)) {
        return None;
    }
    let mut out: Vec<(usize, (f64, f64))> = comp.iter().map(|a| (*a, coords[a])).collect();
    out.sort_by_key(|(a, _)| *a);
    Some(out)
}

fn place_multiring_nonfused_component_rdkit_like(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
    degree: &[usize],
) -> Option<Vec<(usize, (f64, f64))>> {
    #[derive(Clone)]
    struct DepictFrag {
        atoms: std::collections::BTreeMap<usize, TreeEmbeddedAtom>,
        attach_pts: std::collections::VecDeque<usize>,
        done: bool,
    }

    fn update_new_neighs_for_frag(
        frag: &mut DepictFrag,
        mol: &Molecule,
        comp_set: &std::collections::BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        aid: usize,
    ) {
        let mut heavy = Vec::new();
        let mut hydrogens = Vec::new();
        for &nb in &adjacency[aid] {
            if !comp_set.contains(&nb) || frag.atoms.contains_key(&nb) {
                continue;
            }
            if mol.atoms[nb].atomic_num == 1 {
                hydrogens.push(nb);
            } else {
                heavy.push(nb);
            }
        }
        heavy.extend(hydrogens);

        if !heavy.is_empty() && (degree[aid] < 4 || heavy.len() < 3) {
            let cip_ranks = rdkit_cip_ranks_for_depict(mol);
            rdkit_rank_atoms_by_rank(mol, &mut heavy, degree, &cip_ranks);
        } else if degree[aid] >= 4 && heavy.len() >= 3 {
            let cip_ranks = rdkit_cip_ranks_for_depict(mol);
            heavy = rdkit_set_nbr_order(aid, &heavy, mol, adjacency, degree, &cip_ranks);
        }

        if let Some(st) = frag.atoms.get_mut(&aid) {
            st.pending = heavy;
            if !st.pending.is_empty() && !frag.attach_pts.iter().any(|&x| x == aid) {
                frag.attach_pts.push_back(aid);
            }
        }
    }

    fn setup_new_neighs_for_frag(
        frag: &mut DepictFrag,
        mol: &Molecule,
        comp_set: &std::collections::BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
    ) {
        frag.attach_pts.clear();
        let atoms: Vec<usize> = frag.atoms.keys().copied().collect();
        for aid in atoms {
            update_new_neighs_for_frag(frag, mol, comp_set, adjacency, degree, aid);
        }
        let cip_ranks = rdkit_cip_ranks_for_depict(mol);
        let attach = frag.attach_pts.make_contiguous();
        rdkit_rank_atoms_by_rank(mol, attach, degree, &cip_ranks);
    }

    fn find_num_neigh_frag(frag: &DepictFrag, pt: (f64, f64), radius: f64) -> i32 {
        frag.atoms
            .values()
            .filter(|st| {
                let dx = st.loc.0 - pt.0;
                let dy = st.loc.1 - pt.1;
                (dx * dx + dy * dy).sqrt() < radius
            })
            .count() as i32
    }

    fn add_atom_to_atom_with_ang_frag(
        frag: &mut DepictFrag,
        aid: usize,
        to_aid: usize,
    ) -> Option<()> {
        let ref_state = frag.atoms.get(&to_aid)?.clone();
        let nnbr = ref_state.pending.len() as f64;
        let rem_angle = 2.0 * PI - ref_state.angle;
        let mut curr_angle = rem_angle / (1.0 + nnbr);
        if let Some(st) = frag.atoms.get_mut(&to_aid) {
            st.angle += curr_angle;
        }

        let nb1 = frag.atoms.get(&ref_state.nbr1?)?.loc;
        let nb2 = frag.atoms.get(&ref_state.nbr2?)?.loc;
        if frag.atoms.get(&to_aid)?.rot_dir == 0 {
            let rd = rotation_dir(ref_state.loc, nb1, nb2, rem_angle);
            if let Some(st) = frag.atoms.get_mut(&to_aid) {
                st.rot_dir = rd;
            }
        }
        curr_angle *= frag.atoms.get(&to_aid)?.rot_dir as f64;
        let mut curr_loc = rotate_around(nb2, ref_state.loc, curr_angle);
        if rem_angle.abs() - PI < 1e-3 {
            let curr_loc2 = rotate_around(nb2, ref_state.loc, -curr_angle);
            if find_num_neigh_frag(frag, curr_loc, 0.5) > find_num_neigh_frag(frag, curr_loc2, 0.5)
            {
                curr_loc = curr_loc2;
            }
        }
        if let Some(st) = frag.atoms.get_mut(&to_aid) {
            st.nbr2 = Some(aid);
        }

        let tpt = (curr_loc.0 - ref_state.loc.0, curr_loc.1 - ref_state.loc.1);
        let probe_normal = (-tpt.1, tpt.0);
        let tp1 = (curr_loc.0 + probe_normal.0, curr_loc.1 + probe_normal.1);
        let tp2 = (curr_loc.0 - probe_normal.0, curr_loc.1 - probe_normal.1);
        let nccw = find_num_neigh_frag(frag, tp1, 2.5);
        let ncw = find_num_neigh_frag(frag, tp2, 2.5);
        let mut normal = normalize(probe_normal);
        let (ccw, out_normal) = if nccw < ncw {
            (false, normal)
        } else {
            normal = (-normal.0, -normal.1);
            (true, normal)
        };
        frag.atoms.insert(
            aid,
            TreeEmbeddedAtom {
                loc: curr_loc,
                normal: out_normal,
                ccw,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: Some(to_aid),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            },
        );
        Some(())
    }

    fn add_atom_to_atom_with_no_ang_frag(
        frag: &mut DepictFrag,
        mol: &Molecule,
        degree: &[usize],
        aid: usize,
        to_aid: usize,
    ) -> Option<()> {
        let ref_state = frag.atoms.get(&to_aid)?.clone();
        let mut ref_ccw = ref_state.ccw;
        let mut curr_vec = ref_state.normal;
        if ref_state.cis_trans_nbr.is_some_and(|ct| ct != aid) {
            ref_ccw = !ref_ccw;
            curr_vec = (-curr_vec.0, -curr_vec.1);
        }
        let hybridizations = rdkit_hybridizations_for_depict(mol, degree).ok()?;
        let mut angle = compute_sub_angle(degree[to_aid], hybridizations[to_aid]);
        let mut flip_norm = false;

        if frag.atoms.get(&to_aid)?.nbr1.is_some() {
            if let Some(st) = frag.atoms.get_mut(&to_aid) {
                st.angle = angle;
                st.nbr2 = Some(aid);
            }
        } else {
            if mol.atoms[to_aid].atomic_num == 0 {
                unimplemented!("query/dummy atom depict degree handling");
            }
            let rot = rotate(frag.atoms.get(&to_aid)?.normal, angle);
            if let Some(st) = frag.atoms.get_mut(&to_aid) {
                st.normal = rot;
                st.nbr1 = Some(aid);
            }
            flip_norm = true;
        }

        angle -= PI / 2.0;
        if !ref_ccw {
            angle *= -1.0;
        }
        curr_vec = rotate(curr_vec, angle);
        let curr_loc = (
            ref_state.loc.0 + 1.5 * curr_vec.0,
            ref_state.loc.1 + 1.5 * curr_vec.1,
        );

        let tpt = (ref_state.loc.0 - curr_loc.0, ref_state.loc.1 - curr_loc.1);
        let mut normal = (-tpt.1, tpt.0);
        if ref_ccw ^ flip_norm {
            normal = (-normal.0, -normal.1);
        }
        normal = normalize(normal);
        frag.atoms.insert(
            aid,
            TreeEmbeddedAtom {
                loc: curr_loc,
                normal,
                ccw: (!ref_ccw) ^ flip_norm,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: Some(to_aid),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            },
        );
        Some(())
    }

    fn add_non_ring_atom_frag(
        frag: &mut DepictFrag,
        mol: &Molecule,
        comp_set: &std::collections::BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        aid: usize,
        to_aid: usize,
    ) -> Option<()> {
        if frag.atoms.contains_key(&aid) || !frag.atoms.contains_key(&to_aid) {
            return None;
        }
        if frag.atoms.get(&to_aid)?.angle > 0.0 {
            add_atom_to_atom_with_ang_frag(frag, aid, to_aid)?;
        } else {
            add_atom_to_atom_with_no_ang_frag(frag, mol, degree, aid, to_aid)?;
        }
        if let Some(st) = frag.atoms.get_mut(&to_aid) {
            st.pending.retain(|&x| x != aid);
        }
        update_new_neighs_for_frag(frag, mol, comp_set, adjacency, degree, aid);
        Some(())
    }

    fn transform_frag_two_point(
        frag: &mut DepictFrag,
        ref1: (f64, f64),
        ref2: (f64, f64),
        oth1: (f64, f64),
        oth2: (f64, f64),
    ) -> Option<()> {
        let rv = (ref2.0 - ref1.0, ref2.1 - ref1.1);
        let ov = (oth2.0 - oth1.0, oth2.1 - oth1.1);
        let rlen = norm(rv);
        let olen = norm(ov);
        if rlen <= 1e-12 || olen <= 1e-12 {
            return None;
        }
        let cos_t = (ov.0 * rv.0 + ov.1 * rv.1) / (olen * rlen);
        let sin_t = (ov.0 * rv.1 - ov.1 * rv.0) / (olen * rlen);
        for st in frag.atoms.values_mut() {
            let x = st.loc.0 - oth1.0;
            let y = st.loc.1 - oth1.1;
            st.loc = (
                ref1.0 + cos_t * x - sin_t * y,
                ref1.1 + sin_t * x + cos_t * y,
            );
            let n = st.normal;
            st.normal = normalize((cos_t * n.0 - sin_t * n.1, sin_t * n.0 + cos_t * n.1));
        }
        Some(())
    }

    fn reflect_point_line(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
        let vx = b.0 - a.0;
        let vy = b.1 - a.1;
        let len2 = vx * vx + vy * vy;
        if len2 <= 1e-12 {
            return p;
        }
        let wx = p.0 - a.0;
        let wy = p.1 - a.1;
        let proj = (wx * vx + wy * vy) / len2;
        let q = (a.0 + proj * vx, a.1 + proj * vy);
        (2.0 * q.0 - p.0, 2.0 * q.1 - p.1)
    }

    fn reflect_frag(frag: &mut DepictFrag, a: (f64, f64), b: (f64, f64)) {
        let vx = b.0 - a.0;
        let vy = b.1 - a.1;
        let len2 = vx * vx + vy * vy;
        if len2 <= 1e-12 {
            return;
        }
        for st in frag.atoms.values_mut() {
            let old = st.loc;
            st.loc = reflect_point_line(old, a, b);
            let n_end = (old.0 + st.normal.0, old.1 + st.normal.1);
            let n_ref = reflect_point_line(n_end, a, b);
            st.normal = normalize((n_ref.0 - st.loc.0, n_ref.1 - st.loc.1));
            st.ccw = !st.ccw;
        }
    }

    fn reflect_if_necessary_density(
        master: &DepictFrag,
        other: &mut DepictFrag,
        aid1: usize,
        aid2: usize,
    ) {
        let pin1 = master.atoms[&aid1].loc;
        let pin2 = master.atoms[&aid2].loc;
        let mut density_normal = 0.0f64;
        let mut density_reflect = 0.0f64;
        for (&oa, ost) in &other.atoms {
            if master.atoms.contains_key(&oa) {
                continue;
            }
            let loc = ost.loc;
            let rloc = reflect_point_line(loc, pin1, pin2);
            for tst in master.atoms.values() {
                let d = norm((tst.loc.0 - loc.0, tst.loc.1 - loc.1));
                let rd = norm((tst.loc.0 - rloc.0, tst.loc.1 - rloc.1));
                density_normal += if d > 1e-3 { 1.0 / d } else { 1000.0 };
                density_reflect += if rd > 1e-3 { 1.0 / rd } else { 1000.0 };
            }
        }
        if std::env::var_os("COSMOLKIT_DEBUG_FUSED_DENSITY").is_some() {
            eprintln!(
                "DEBUG density pins=({aid1},{aid2}) normal={density_normal:.8} reflect={density_reflect:.8} choose_reflect={}",
                density_normal - density_reflect > 1e-4
            );
        }
        if density_normal - density_reflect > 1e-4 {
            reflect_frag(other, pin1, pin2);
        }
    }

    fn merge_no_common_frag(
        master: &mut DepictFrag,
        other: &mut DepictFrag,
        mol: &Molecule,
        comp_set: &std::collections::BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        to_aid: usize,
        nbr_aid: usize,
    ) -> Option<()> {
        add_non_ring_atom_frag(master, mol, comp_set, adjacency, degree, nbr_aid, to_aid)?;
        add_non_ring_atom_frag(other, mol, comp_set, adjacency, degree, to_aid, nbr_aid)?;
        let ref1 = master.atoms.get(&to_aid)?.loc;
        let ref2 = master.atoms.get(&nbr_aid)?.loc;
        let oth1 = other.atoms.get(&to_aid)?.loc;
        let oth2 = other.atoms.get(&nbr_aid)?.loc;
        transform_frag_two_point(other, ref1, ref2, oth1, oth2)?;
        reflect_if_necessary_density(master, other, to_aid, nbr_aid);

        let common = [to_aid, nbr_aid];
        let mut to_update = Vec::new();
        for (&aid, ost) in &other.atoms {
            if !common.contains(&aid) {
                master.atoms.insert(aid, ost.clone());
                if !ost.pending.is_empty() && !master.attach_pts.iter().any(|&x| x == aid) {
                    master.attach_pts.push_back(aid);
                }
            } else if let Some(mst) = master.atoms.get_mut(&aid) {
                if ost.cis_trans_nbr.is_some() {
                    mst.cis_trans_nbr = ost.cis_trans_nbr;
                    mst.normal = ost.normal;
                    mst.ccw = ost.ccw;
                }
                if ost.angle > 0.0 {
                    mst.angle = ost.angle;
                    mst.nbr1 = ost.nbr1;
                    mst.nbr2 = ost.nbr2;
                }
            }
            to_update.push(aid);
        }
        for aid in common {
            to_update.push(aid);
        }
        to_update.sort_unstable();
        to_update.dedup();
        for aid in to_update {
            if master.atoms.contains_key(&aid) {
                update_new_neighs_for_frag(master, mol, comp_set, adjacency, degree, aid);
            }
        }
        other.done = true;
        Some(())
    }

    fn init_ring_frag_from_order(ring: &[usize]) -> DepictFrag {
        let n = ring.len();
        let radius = 1.5f64 / (2.0 * (PI / n as f64).sin());
        let largest_angle = PI * (1.0 - (2.0 / n as f64));
        let mut frag = DepictFrag {
            atoms: std::collections::BTreeMap::new(),
            attach_pts: std::collections::VecDeque::new(),
            done: false,
        };
        for (k, &a) in ring.iter().enumerate() {
            let theta = 2.0 * PI * k as f64 / n as f64;
            let prev = if k == 0 { ring[n - 1] } else { ring[k - 1] };
            let next = ring[(k + 1) % n];
            frag.atoms.insert(
                a,
                TreeEmbeddedAtom {
                    loc: (radius * theta.cos(), radius * theta.sin()),
                    normal: (1.0, 0.0),
                    ccw: true,
                    cis_trans_nbr: None,
                    angle: largest_angle,
                    nbr1: Some(prev),
                    nbr2: Some(next),
                    rot_dir: 0,
                    pending: Vec::new(),
                },
            );
        }
        frag
    }

    fn merge_ring_frag(
        master: &mut DepictFrag,
        emb_ring: &DepictFrag,
        n_common: usize,
        pin_atoms: &[usize],
    ) {
        for (&aid, ost) in &emb_ring.atoms {
            if !master.atoms.contains_key(&aid) {
                master.atoms.insert(aid, ost.clone());
            } else if n_common <= 2 && pin_atoms.contains(&aid) {
                if let Some(mst) = master.atoms.get_mut(&aid) {
                    mst.angle += ost.angle;
                    if mst.nbr1 == ost.nbr1 {
                        mst.nbr1 = ost.nbr2;
                    } else if mst.nbr1 == ost.nbr2 {
                        mst.nbr1 = ost.nbr1;
                    } else if mst.nbr2 == ost.nbr1 {
                        mst.nbr2 = ost.nbr2;
                    } else if mst.nbr2 == ost.nbr2 {
                        mst.nbr2 = ost.nbr1;
                    }
                }
            }
        }
    }

    fn share_ring_edge(a: &[usize], b: &[usize]) -> bool {
        for i in 0..a.len() {
            let a1 = a[i];
            let a2 = a[(i + 1) % a.len()];
            for j in 0..b.len() {
                let b1 = b[j];
                let b2 = b[(j + 1) % b.len()];
                if (a1 == b1 && a2 == b2) || (a1 == b2 && a2 == b1) {
                    return true;
                }
            }
        }
        false
    }

    fn share_ring_atom(a: &[usize], b: &[usize]) -> bool {
        a.iter().any(|x| b.contains(x))
    }

    fn pick_first_ring_to_embed(degree: &[usize], cycles: &[Vec<usize>], group: &[usize]) -> usize {
        let mut res = group[0];
        let mut min_subs = usize::MAX;
        let mut max_size = 0usize;
        for &rid in group {
            let ring = &cycles[rid];
            let subs = ring.iter().filter(|&&a| degree[a] > 2).count();
            if subs < min_subs || (subs == min_subs && ring.len() > max_size) {
                res = rid;
                min_subs = subs;
                max_size = ring.len();
            }
        }
        res
    }

    fn build_fused_group_frag(
        mol: &Molecule,
        comp_set: &std::collections::BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cycles: &[Vec<usize>],
        group: &[usize],
    ) -> Option<DepictFrag> {
        let first = pick_first_ring_to_embed(degree, cycles, group);
        let mut master = init_ring_frag_from_order(&cycles[first]);
        let union: std::collections::BTreeSet<usize> = group
            .iter()
            .flat_map(|&rid| cycles[rid].iter().copied())
            .collect();
        let mut done = std::collections::BTreeSet::<usize>::new();
        done.insert(first);

        while master.atoms.len() < union.len() {
            let done_atoms: std::collections::BTreeSet<usize> = done
                .iter()
                .flat_map(|&rid| cycles[rid].iter().copied())
                .collect();
            let mut best: Option<(usize, Vec<usize>)> = None;
            let mut max_common = 0usize;
            for &rid in group {
                if done.contains(&rid) {
                    continue;
                }
                let common: Vec<usize> = cycles[rid]
                    .iter()
                    .copied()
                    .filter(|a| done_atoms.contains(a))
                    .collect();
                if common.is_empty() {
                    continue;
                }
                if common.len() == 2 {
                    best = Some((rid, common));
                    break;
                }
                if common.len() > max_common {
                    max_common = common.len();
                    best = Some((rid, common));
                }
            }
            let (rid, mut common) = best?;
            let cmn_lst = common
                .iter()
                .zip(cycles[rid].iter())
                .take_while(|(a, b)| *a == *b)
                .count();
            if cmn_lst > 0 && cmn_lst < common.len() {
                common.rotate_left(cmn_lst);
            }
            if std::env::var_os("COSMOLKIT_DEBUG_FUSED_DENSITY").is_some() {
                eprintln!("DEBUG next ring rid={rid} common={common:?}");
            }
            let mut emb_ring = init_ring_frag_from_order(&cycles[rid]);
            let mut pin_atoms = Vec::new();
            if common.len() == 1 {
                let aid = common[0];
                pin_atoms.push(aid);
                let rcr = master.atoms.get(&aid)?.loc;
                let oeatm = emb_ring.atoms.get(&aid)?.clone();
                let ccr = oeatm.loc;
                let onb1 = oeatm.nbr1?;
                let onb2 = oeatm.nbr2?;
                let onb1_loc = emb_ring.atoms.get(&onb1)?.loc;
                let onb2_loc = emb_ring.atoms.get(&onb2)?.loc;
                let mid_pt = (
                    (onb1_loc.0 + onb2_loc.0) * 0.5,
                    (onb1_loc.1 + onb2_loc.1) * 0.5,
                );

                let nb1 = master.atoms.get(&aid)?.nbr1?;
                let nb2 = master.atoms.get(&aid)?.nbr2?;
                let nbp1 = master.atoms.get(&nb1)?.loc;
                let nbp2 = master.atoms.get(&nb2)?.loc;
                let ang = master.atoms.get(&aid)?.angle;
                let largest_angle = 2.0 * PI - ang;
                let mut bpt = ((nbp1.0 + nbp2.0) * 0.5, (nbp1.1 + nbp2.1) * 0.5);
                if largest_angle > PI {
                    bpt = (2.0 * rcr.0 - bpt.0, 2.0 * rcr.1 - bpt.1);
                }
                transform_frag_two_point(&mut emb_ring, rcr, bpt, ccr, mid_pt)?;
            } else {
                let aid1 = *common.first()?;
                let aid2 = *common.last()?;
                pin_atoms.push(aid1);
                pin_atoms.push(aid2);
                let ref1 = master.atoms.get(&aid1)?.loc;
                let ref2 = master.atoms.get(&aid2)?.loc;
                let oth1 = emb_ring.atoms.get(&aid1)?.loc;
                let oth2 = emb_ring.atoms.get(&aid2)?.loc;
                transform_frag_two_point(&mut emb_ring, ref1, ref2, oth1, oth2)?;
                reflect_if_necessary_density(&master, &mut emb_ring, aid1, aid2);
            }
            merge_ring_frag(&mut master, &emb_ring, common.len(), &pin_atoms);
            done.insert(rid);
        }

        setup_new_neighs_for_frag(&mut master, mol, comp_set, adjacency, degree);
        Some(master)
    }

    let comp_set: std::collections::BTreeSet<usize> = comp.iter().copied().collect();
    let mut deg_in_comp = vec![0usize; mol.atoms.len()];
    for &a in comp {
        deg_in_comp[a] = adjacency[a].iter().filter(|n| comp_set.contains(n)).count();
    }
    let mut removed = std::collections::BTreeSet::<usize>::new();
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
    let ring_set: std::collections::BTreeSet<usize> = ring_atoms.iter().copied().collect();

    fn rdkit_single_ring_order(
        mol: &Molecule,
        ring_set: &std::collections::BTreeSet<usize>,
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
        let mut bfsq = std::collections::VecDeque::<usize>::new();
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
        let mut bfsq = std::collections::VecDeque::<usize>::new();
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
        changed: &mut std::collections::BTreeSet<usize>,
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
                changed.insert(other);
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

    fn rdkit_find_sssr_orders(
        mol: &Molecule,
        comp: &[usize],
        target_cycle_count: usize,
    ) -> Vec<Vec<usize>> {
        let mut atom_degrees = vec![0usize; mol.atoms.len()];
        for bond in &mol.bonds {
            atom_degrees[bond.begin_atom] += 1;
            atom_degrees[bond.end_atom] += 1;
        }
        let mut active_bonds = vec![true; mol.bonds.len()];
        let mut changed = std::collections::BTreeSet::<usize>::new();
        for &idx in comp {
            if atom_degrees[idx] < 2 {
                changed.insert(idx);
            }
        }

        let mut done_atoms = vec![false; mol.atoms.len()];
        let mut n_atoms_done = 0usize;
        let mut rings = Vec::<Vec<usize>>::new();
        let mut invariants = std::collections::BTreeSet::<Vec<usize>>::new();

        while n_atoms_done < comp.len() {
            while let Some(cand) = changed.pop_first() {
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

            if rings.len() >= target_cycle_count {
                break;
            }
        }
        rings
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

    let ring_bond_count = mol
        .bonds
        .iter()
        .filter(|b| ring_set.contains(&b.begin_atom) && ring_set.contains(&b.end_atom))
        .count();
    let target_cycle_count = ring_bond_count + 1 - ring_atoms.len();

    let mut seen = std::collections::BTreeSet::<Vec<usize>>::new();
    let mut all_cycles: Vec<Vec<usize>> = Vec::new();
    let mut ring_atoms_sorted = ring_atoms.clone();
    ring_atoms_sorted.sort_unstable();
    const MAX_RING_ENUM: usize = 8;
    for &start in &ring_atoms_sorted {
        let mut path = vec![start];
        let mut used = std::collections::BTreeSet::<usize>::new();
        used.insert(start);
        fn dfs_cycles(
            start: usize,
            cur: usize,
            path: &mut Vec<usize>,
            used: &mut std::collections::BTreeSet<usize>,
            adjacency: &[Vec<usize>],
            ring_set: &std::collections::BTreeSet<usize>,
            seen: &mut std::collections::BTreeSet<Vec<usize>>,
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
            adjacency,
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
    if all_cycles.is_empty() {
        return None;
    }

    all_cycles.sort_by_key(|cyc| (cyc.len(), cyc.clone()));
    let mut cycles: Vec<Vec<usize>> = rdkit_find_sssr_orders(mol, comp, target_cycle_count);
    if target_cycle_count == 1 && cycles.is_empty() {
        cycles.push(rdkit_single_ring_order(mol, &ring_set, &deg_in_comp)?);
    } else if cycles.len() < target_cycle_count {
        let mut covered_edges = std::collections::BTreeSet::<(usize, usize)>::new();
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

    if cycles.is_empty() {
        return None;
    }
    if std::env::var_os("COSMOLKIT_DEBUG_RING_CYCLES").is_some() {
        eprintln!("DEBUG cycles: {:?}", cycles);
    }
    let embedded_ring_set: std::collections::BTreeSet<usize> =
        cycles.iter().flatten().copied().collect();
    let mut parent: Vec<usize> = (0..cycles.len()).collect();
    fn find_parent(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            let root = find_parent(parent, parent[x]);
            parent[x] = root;
        }
        parent[x]
    }
    fn union_parent(parent: &mut [usize], a: usize, b: usize) {
        let ra = find_parent(parent, a);
        let rb = find_parent(parent, b);
        if ra != rb {
            parent[rb] = ra;
        }
    }
    for i in 0..cycles.len() {
        for j in i + 1..cycles.len() {
            if share_ring_edge(&cycles[i], &cycles[j]) || share_ring_atom(&cycles[i], &cycles[j]) {
                union_parent(&mut parent, i, j);
            }
        }
    }
    let mut grouped = std::collections::BTreeMap::<usize, Vec<usize>>::new();
    for i in 0..cycles.len() {
        let root = find_parent(&mut parent, i);
        grouped.entry(root).or_default().push(i);
    }

    let mut frags = Vec::<DepictFrag>::new();
    for group in grouped.values() {
        let frag = if group.len() == 1 {
            let ring = &cycles[group[0]];
            let ring_set: std::collections::BTreeSet<usize> = ring.iter().copied().collect();
            let ring_order = rdkit_single_ring_order(mol, &ring_set, &deg_in_comp)
                .unwrap_or_else(|| ring.clone());
            let mut frag = init_ring_frag_from_order(&ring_order);
            setup_new_neighs_for_frag(&mut frag, mol, &comp_set, adjacency, degree);
            frag
        } else {
            build_fused_group_frag(mol, &comp_set, adjacency, degree, &cycles, group)?
        };
        frags.push(frag);
    }

    let mut non_ring_atoms: std::collections::BTreeSet<usize> = comp
        .iter()
        .copied()
        .filter(|a| !embedded_ring_set.contains(a))
        .collect();

    let mut master_idx = frags
        .iter()
        .enumerate()
        .max_by_key(|(_, f)| f.atoms.len())
        .map(|(idx, _)| idx)?;
    // RDKit _findLargestFrag keeps the first fragment on size ties.
    for (idx, frag) in frags.iter().enumerate() {
        if !frag.done && frag.atoms.len() == frags[master_idx].atoms.len() {
            master_idx = idx;
            break;
        }
    }
    let mut master = frags.remove(master_idx);
    master.done = true;

    while !master.attach_pts.is_empty() || !non_ring_atoms.is_empty() {
        if master.attach_pts.is_empty() {
            unimplemented!("RDKit fallback to start a new non-ring fragment after ring expansion");
        }
        let aid = master.attach_pts.pop_front()?;
        let nbrs = master
            .atoms
            .get(&aid)
            .map(|st| st.pending.clone())
            .unwrap_or_default();
        for nb in nbrs {
            if non_ring_atoms.contains(&nb) {
                add_non_ring_atom_frag(&mut master, mol, &comp_set, adjacency, degree, nb, aid)?;
                non_ring_atoms.remove(&nb);
            } else {
                let mut found = None;
                for (idx, frag) in frags.iter().enumerate() {
                    if !frag.done && frag.atoms.contains_key(&nb) {
                        found = Some(idx);
                        break;
                    }
                }
                if let Some(idx) = found {
                    let mut other = frags.remove(idx);
                    merge_no_common_frag(
                        &mut master,
                        &mut other,
                        mol,
                        &comp_set,
                        adjacency,
                        degree,
                        aid,
                        nb,
                    )?;
                }
            }
        }
        if let Some(st) = master.atoms.get_mut(&aid) {
            st.pending.clear();
        }
    }

    if comp.iter().any(|a| !master.atoms.contains_key(a)) {
        return None;
    }

    let mut local: Vec<(usize, (f64, f64))> = comp
        .iter()
        .copied()
        .map(|a| (a, master.atoms[&a].loc))
        .collect();
    local.sort_by_key(|(idx, _)| *idx);
    if std::env::var_os("COSMOLKIT_DEBUG_CLEANUP").is_some() {
        eprintln!("DEBUG before cleanup: {:?}", local);
    }

    remove_collisions_bond_flip_like(mol, comp, adjacency, &mut local);
    if std::env::var_os("COSMOLKIT_DEBUG_CLEANUP").is_some() {
        eprintln!("DEBUG after bond_flip: {:?}", local);
    }
    remove_collisions_open_angles_like(mol, comp, adjacency, &mut local);
    if std::env::var_os("COSMOLKIT_DEBUG_CLEANUP").is_some() {
        eprintln!("DEBUG after open_angles: {:?}", local);
    }
    remove_collisions_shorten_bonds_like(mol, comp, adjacency, &mut local);
    if std::env::var_os("COSMOLKIT_DEBUG_CLEANUP").is_some() {
        eprintln!("DEBUG after cleanup: {:?}", local);
    }
    Some(local)
}

fn remove_collisions_shorten_bonds_like(
    mol: &Molecule,
    comp: &[usize],
    adjacency: &[Vec<usize>],
    local: &mut Vec<(usize, (f64, f64))>,
) {
    const COLLISION_THRES: f64 = 0.70;
    const HETEROATOM_COLL_SCALE: f64 = 1.3;
    const MAX_COLL_ITERS: usize = 15;

    let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();
    let mut pos = std::collections::HashMap::<usize, (f64, f64)>::new();
    for &(a, p) in local.iter() {
        pos.insert(a, p);
    }

    let find_collisions =
        |pos: &std::collections::HashMap<usize, (f64, f64)>| -> Vec<(usize, usize)> {
            let mut out = Vec::new();
            let mut atoms = comp.to_vec();
            atoms.sort_unstable();
            for i in 0..atoms.len() {
                let a = atoms[i];
                let fa = if mol.atoms[a].atomic_num != 6 {
                    HETEROATOM_COLL_SCALE
                } else {
                    1.0
                };
                for &b in atoms.iter().take(i) {
                    let fb = if mol.atoms[b].atomic_num != 6 {
                        HETEROATOM_COLL_SCALE
                    } else {
                        1.0
                    };
                    let pa = pos[&a];
                    let pb = pos[&b];
                    let dx = pa.0 - pb.0;
                    let dy = pa.1 - pb.1;
                    let mut d2 = dx * dx + dy * dy;
                    d2 /= fa * fb;
                    if d2 < COLLISION_THRES * COLLISION_THRES {
                        out.push((a, b));
                    }
                }
            }
            out
        };

    let degree_in_comp =
        |a: usize| -> usize { adjacency[a].iter().filter(|n| comp_set.contains(n)).count() };
    let find_deg1_neighbor = |aid: usize| -> Option<usize> {
        let nbs: Vec<usize> = adjacency[aid]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .collect();
        if nbs.len() == 1 { Some(nbs[0]) } else { None }
    };
    let edge_is_ring_like = |u: usize, v: usize| -> bool {
        let mut stack = vec![u];
        let mut seen = std::collections::HashSet::<usize>::new();
        seen.insert(u);
        while let Some(cur) = stack.pop() {
            for &nb in &adjacency[cur] {
                if !comp_set.contains(&nb) {
                    continue;
                }
                if (cur == u && nb == v) || (cur == v && nb == u) {
                    continue;
                }
                if nb == v {
                    return true;
                }
                if seen.insert(nb) {
                    stack.push(nb);
                }
            }
        }
        false
    };
    let shortest_path = |src: usize, dst: usize| -> Option<Vec<usize>> {
        let mut q = std::collections::VecDeque::<usize>::new();
        let mut prev = std::collections::HashMap::<usize, usize>::new();
        q.push_back(src);
        prev.insert(src, src);
        while let Some(u) = q.pop_front() {
            if u == dst {
                break;
            }
            for &v in &adjacency[u] {
                if !comp_set.contains(&v) {
                    continue;
                }
                if let std::collections::hash_map::Entry::Vacant(e) = prev.entry(v) {
                    e.insert(u);
                    q.push_back(v);
                }
            }
        }
        if !prev.contains_key(&dst) {
            return None;
        }
        let mut path = vec![dst];
        let mut cur = dst;
        while cur != src {
            cur = prev[&cur];
            path.push(cur);
        }
        path.reverse();
        Some(path)
    };

    let mut colls = find_collisions(&pos);
    let mut iter = 0usize;
    while !colls.is_empty() && iter < MAX_COLL_ITERS {
        let (mut aid1, mut aid2) = colls[0];
        let mut deg1 = degree_in_comp(aid1);
        let mut deg2 = degree_in_comp(aid2);
        if deg2 > deg1 {
            std::mem::swap(&mut aid1, &mut aid2);
            std::mem::swap(&mut deg1, &mut deg2);
        }
        if let Some(mut path) = shortest_path(aid1, aid2) {
            if !path.is_empty() {
                // pop head (aid1), mirror RDKit logic
                path.remove(0);
            }
            let mut has_non_ring = false;
            let mut prev = aid1;
            for &next in &path {
                if !edge_is_ring_like(prev, next) {
                    has_non_ring = true;
                    break;
                }
                prev = next;
            }
            if has_non_ring {
                if deg1 == 1 {
                    if let Some(a) = find_deg1_neighbor(aid1) {
                        let pa = pos[&a];
                        let mut v = (pos[&aid1].0 - pa.0, pos[&aid1].1 - pa.1);
                        v.0 *= 0.9;
                        v.1 *= 0.9;
                        let len = (v.0 * v.0 + v.1 * v.1).sqrt();
                        if len > 0.75 {
                            pos.insert(aid1, (pa.0 + v.0, pa.1 + v.1));
                        }
                    }
                }
                if deg2 == 1 {
                    if let Some(a) = find_deg1_neighbor(aid2) {
                        let pa = pos[&a];
                        let mut v = (pos[&aid2].0 - pa.0, pos[&aid2].1 - pa.1);
                        v.0 *= 0.9;
                        v.1 *= 0.9;
                        let len = (v.0 * v.0 + v.1 * v.1).sqrt();
                        if len > 0.75 {
                            pos.insert(aid2, (pa.0 + v.0, pa.1 + v.1));
                        }
                    }
                }
            }
        }
        colls = find_collisions(&pos);
        iter += 1;
    }

    for item in local.iter_mut() {
        if let Some(&p) = pos.get(&item.0) {
            item.1 = p;
        }
    }
}

fn atom_symbol(atomic_num: u8) -> &'static str {
    match atomic_num {
        1 => "H",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        11 => "Na",
        14 => "Si",
        15 => "P",
        16 => "S",
        17 => "Cl",
        29 => "Cu",
        34 => "Se",
        35 => "Br",
        45 => "Rh",
        53 => "I",
        _ => "*",
    }
}

fn bond_type_code(order: BondOrder) -> usize {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Dative => 9,
        BondOrder::Null => 8,
        BondOrder::Quadruple => 0,
    }
}

#[cfg(test)]
mod tests {
    use core::fmt;
    use std::collections::HashSet;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::{Path, PathBuf};
    use std::process::Command;

    use super::mol_to_v2000_block_minimal;
    use crate::Molecule;
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    struct GoldenRecord {
        smiles: String,
        parse_ok: bool,
        parse_error: Option<String>,
        v2000_ok: bool,
        v2000_body: Option<String>,
        v2000_error: Option<String>,
        v3000_ok: bool,
        v3000_body: Option<String>,
        v3000_error: Option<String>,
    }

    #[derive(Debug, Deserialize)]
    struct KekulizeGoldenRecord {
        smiles: String,
        parse_ok: bool,
        parse_error: Option<String>,
        kekulize_ok: bool,
        kekulize_error: Option<String>,
        v2000_ok: bool,
        v2000_body: Option<String>,
        v2000_error: Option<String>,
        v3000_ok: bool,
        v3000_body: Option<String>,
        v3000_error: Option<String>,
    }

    #[derive(Debug)]
    enum TestDataError {
        Io(std::io::Error),
        Json {
            line_no: usize,
            source: serde_json::Error,
        },
    }

    impl fmt::Display for TestDataError {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                Self::Io(err) => write!(f, "{err}"),
                Self::Json { line_no, source } => {
                    write!(f, "invalid jsonl at line {line_no}: {source}")
                }
            }
        }
    }

    impl std::error::Error for TestDataError {}

    impl From<std::io::Error> for TestDataError {
        fn from(value: std::io::Error) -> Self {
            Self::Io(value)
        }
    }

    fn repo_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("crates/")
            .parent()
            .expect("repo root")
            .to_path_buf()
    }

    fn body(block: &str) -> String {
        let lines: Vec<_> = block.lines().collect();
        lines[3..].join("\n")
    }

    fn normalize_signed_zero(text: &str) -> String {
        text.replace("-0.0000", " 0.0000")
    }

    fn atom_symbol_from_v2000_line(line: &str) -> String {
        if line.len() >= 34 {
            return line[31..34].trim().to_owned();
        }
        String::new()
    }

    fn atom_symbols_equivalent(ours: &str, expected: &str) -> bool {
        if ours == expected {
            return true;
        }
        // RDKit writes mapped dummies as "R" atoms in V2000 while our minimal
        // subset keeps "*" placeholders.
        (ours == "*" && expected == "R") || (ours == "R" && expected == "*")
    }

    fn parse_bond_line(line: &str) -> (usize, usize, usize) {
        let a = line[0..3].trim().parse::<usize>().expect("bond a index");
        let b = line[3..6].trim().parse::<usize>().expect("bond b index");
        let order = line[6..9].trim().parse::<usize>().expect("bond order");
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        (lo, hi, order)
    }

    #[derive(Debug)]
    struct ParsedBody {
        atoms: Vec<String>,
        coords: Vec<(f64, f64)>,
        bonds: Vec<(usize, usize, usize)>,
    }

    fn parse_v2000_body(block_body: &str) -> Option<ParsedBody> {
        let lines: Vec<_> = block_body.lines().collect();
        let counts = *lines.first()?;
        let toks: Vec<_> = counts.split_whitespace().collect();
        let atom_count = toks.first()?.parse::<usize>().ok()?;
        let bond_count = toks.get(1)?.parse::<usize>().ok()?;
        if lines.len() < 1 + atom_count + bond_count {
            return None;
        }

        let mut atoms = Vec::with_capacity(atom_count);
        let mut coords = Vec::with_capacity(atom_count);
        for i in 0..atom_count {
            let line = lines[1 + i];
            atoms.push(atom_symbol_from_v2000_line(line));
            let x = line[0..10].trim().parse::<f64>().ok()?;
            let y = line[10..20].trim().parse::<f64>().ok()?;
            coords.push((x, y));
        }

        let mut bonds = Vec::with_capacity(bond_count);
        for i in 0..bond_count {
            bonds.push(parse_bond_line(lines[1 + atom_count + i]));
        }
        Some(ParsedBody {
            atoms,
            coords,
            bonds,
        })
    }

    fn parse_v3000_body(block_body: &str) -> Option<ParsedBody> {
        let mut in_atom = false;
        let mut in_bond = false;
        let mut atoms_raw: Vec<(usize, String, f64, f64)> = Vec::new();
        let mut bonds: Vec<(usize, usize, usize)> = Vec::new();

        for line in block_body.lines() {
            let t = line.trim();
            if t == "M  V30 BEGIN ATOM" {
                in_atom = true;
                continue;
            }
            if t == "M  V30 END ATOM" {
                in_atom = false;
                continue;
            }
            if t == "M  V30 BEGIN BOND" {
                in_bond = true;
                continue;
            }
            if t == "M  V30 END BOND" {
                in_bond = false;
                continue;
            }

            if in_atom && t.starts_with("M  V30 ") {
                let toks: Vec<_> = t.split_whitespace().collect();
                if toks.len() < 5 {
                    return None;
                }
                let idx = toks[2].parse::<usize>().ok()?;
                let symbol = toks[3].to_owned();
                let x = toks[4].parse::<f64>().ok()?;
                let y = toks[5].parse::<f64>().ok()?;
                atoms_raw.push((idx, symbol, x, y));
            } else if in_bond && t.starts_with("M  V30 ") {
                let toks: Vec<_> = t.split_whitespace().collect();
                if toks.len() < 6 {
                    return None;
                }
                let order = toks[3].parse::<usize>().ok()?;
                let a = toks[4].parse::<usize>().ok()?;
                let b = toks[5].parse::<usize>().ok()?;
                let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
                bonds.push((lo, hi, order));
            }
        }

        if atoms_raw.is_empty() {
            return None;
        }
        atoms_raw.sort_by_key(|(idx, _, _, _)| *idx);
        let mut atoms = Vec::with_capacity(atoms_raw.len());
        let mut coords = Vec::with_capacity(atoms_raw.len());
        for (_, s, x, y) in atoms_raw {
            atoms.push(s);
            coords.push((x, y));
        }
        Some(ParsedBody {
            atoms,
            coords,
            bonds,
        })
    }

    fn parse_body_for_compare(block_body: &str) -> Option<ParsedBody> {
        let first = block_body.lines().next().unwrap_or_default();
        if first.contains("V3000") {
            parse_v3000_body(block_body)
        } else {
            parse_v2000_body(block_body)
        }
    }

    fn canonical_bonds_for_compare(
        bonds: &[(usize, usize, usize)],
        flexible_pairs: &HashSet<(usize, usize)>,
    ) -> Vec<(usize, usize, usize)> {
        let mut out = Vec::with_capacity(bonds.len());
        for &(a, b, mut order) in bonds {
            if flexible_pairs.contains(&(a, b)) {
                order = 0;
            }
            out.push((a, b, order));
        }
        out.sort_unstable();
        out
    }

    fn coords_match_strict(ours: &[(f64, f64)], expected: &[(f64, f64)]) -> bool {
        if ours.len() != expected.len() {
            return false;
        }
        let coord_tol = 1e-3f64;
        for i in 0..ours.len() {
            let (ox, oy) = ours[i];
            let (ex, ey) = expected[i];
            if (ox - ex).abs() > coord_tol || (oy - ey).abs() > coord_tol {
                return false;
            }
        }
        true
    }

    fn align_rigid_2d(
        source: &[(f64, f64)],
        target: &[(f64, f64)],
    ) -> Option<(f64, f64, f64, f64)> {
        if source.len() != target.len() || source.is_empty() {
            return None;
        }

        let n = source.len() as f64;
        let src_cx = source.iter().map(|(x, _)| x).sum::<f64>() / n;
        let src_cy = source.iter().map(|(_, y)| y).sum::<f64>() / n;
        let tgt_cx = target.iter().map(|(x, _)| x).sum::<f64>() / n;
        let tgt_cy = target.iter().map(|(_, y)| y).sum::<f64>() / n;

        let mut a = 0.0f64;
        let mut b = 0.0f64;
        for i in 0..source.len() {
            let sx = source[i].0 - src_cx;
            let sy = source[i].1 - src_cy;
            let tx = target[i].0 - tgt_cx;
            let ty = target[i].1 - tgt_cy;
            a += sx * tx + sy * ty;
            b += sx * ty - sy * tx;
        }

        let norm = (a * a + b * b).sqrt();
        if norm <= 1e-12 {
            return None;
        }
        let cos_t = a / norm;
        let sin_t = b / norm;
        let tx = tgt_cx - (cos_t * src_cx - sin_t * src_cy);
        let ty = tgt_cy - (sin_t * src_cx + cos_t * src_cy);
        Some((cos_t, sin_t, tx, ty))
    }

    fn apply_rigid_2d(p: (f64, f64), tf: (f64, f64, f64, f64)) -> (f64, f64) {
        let (cos_t, sin_t, tx, ty) = tf;
        let x = cos_t * p.0 - sin_t * p.1 + tx;
        let y = sin_t * p.0 + cos_t * p.1 + ty;
        (x, y)
    }

    fn rmsd_for_indices(
        ours: &[(f64, f64)],
        expected: &[(f64, f64)],
        indices_1based: &[usize],
        tf: (f64, f64, f64, f64),
    ) -> f64 {
        let mut sse = 0.0f64;
        for &idx in indices_1based {
            let i = idx - 1;
            let po = apply_rigid_2d(ours[i], tf);
            let pe = expected[i];
            let dx = po.0 - pe.0;
            let dy = po.1 - pe.1;
            sse += dx * dx + dy * dy;
        }
        (sse / (indices_1based.len() as f64)).sqrt()
    }

    fn non_coordinate_sections_match(ours_body: &str, expected_body: &str, mol: &Molecule) -> bool {
        let Some(ours) = parse_body_for_compare(ours_body) else {
            return false;
        };
        let Some(expected) = parse_body_for_compare(expected_body) else {
            return false;
        };

        if ours.atoms.len() != expected.atoms.len()
            || ours.coords.len() != expected.coords.len()
            || ours.bonds.len() != expected.bonds.len()
        {
            return false;
        }
        let coords_ok = coords_match_strict(&ours.coords, &expected.coords);

        for i in 0..ours.atoms.len() {
            let ours_symbol = &ours.atoms[i];
            let expected_symbol = &expected.atoms[i];
            if !atom_symbols_equivalent(&ours_symbol, &expected_symbol) {
                return false;
            }
        }

        let flexible_pairs: HashSet<(usize, usize)> = mol
            .bonds
            .iter()
            .filter(|b| {
                matches!(
                    b.order,
                    crate::BondOrder::Aromatic | crate::BondOrder::Dative | crate::BondOrder::Null
                )
            })
            .map(|b| {
                let a = b.begin_atom + 1;
                let c = b.end_atom + 1;
                if a <= c { (a, c) } else { (c, a) }
            })
            .collect();

        canonical_bonds_for_compare(&ours.bonds, &flexible_pairs)
            == canonical_bonds_for_compare(&expected.bonds, &flexible_pairs)
            && coords_ok
    }

    fn compare_against_expected(
        ours_body: &str,
        expected_body: &str,
        mol: &Molecule,
        smiles: &str,
        row_idx_1based: usize,
        variant: &str,
    ) {
        let ours_norm = normalize_signed_zero(ours_body);
        let expected_norm = normalize_signed_zero(expected_body);
        let mismatch_detail = {
            let ours = parse_body_for_compare(&ours_norm);
            let expected = parse_body_for_compare(&expected_norm);
            match (ours, expected) {
                (Some(ours), Some(expected)) => {
                    let mut detail = String::new();
                    if ours.atoms.len() == expected.atoms.len() {
                        if let Some((idx, (o, e))) = ours
                            .coords
                            .iter()
                            .zip(expected.coords.iter())
                            .enumerate()
                            .find(|(_, (o, e))| {
                                (o.0 - e.0).abs() > 1e-3f64 || (o.1 - e.1).abs() > 1e-3f64
                            })
                        {
                            detail.push_str(&format!(
                                "first coordinate mismatch at atom {}: ours=({:.4},{:.4}) expected=({:.4},{:.4})\n",
                                idx + 1,
                                o.0,
                                o.1,
                                e.0,
                                e.1
                            ));
                        }
                    }
                    detail.push_str("ours body:\n");
                    detail.push_str(&ours_norm);
                    detail.push_str("\nexpected body:\n");
                    detail.push_str(&expected_norm);
                    detail
                }
                _ => format!(
                    "failed to parse ours/expected body\nours:\n{ours_norm}\nexpected:\n{expected_norm}"
                ),
            }
        };
        assert!(
            non_coordinate_sections_match(&ours_norm, &expected_norm, mol),
            "molblock mismatch (including coordinates) at row {} ({}) against {}\n{}",
            row_idx_1based,
            smiles,
            variant,
            mismatch_detail
        );
    }

    fn bonds_and_atoms_match_strict(ours_body: &str, expected_body: &str) -> bool {
        let Some(ours) = parse_body_for_compare(ours_body) else {
            return false;
        };
        let Some(expected) = parse_body_for_compare(expected_body) else {
            return false;
        };

        if ours.atoms.len() != expected.atoms.len() || ours.bonds.len() != expected.bonds.len() {
            return false;
        }
        for i in 0..ours.atoms.len() {
            if !atom_symbols_equivalent(&ours.atoms[i], &expected.atoms[i]) {
                return false;
            }
        }
        let mut ours_bonds = ours.bonds;
        let mut expected_bonds = expected.bonds;
        ours_bonds.sort_unstable();
        expected_bonds.sort_unstable();
        ours_bonds == expected_bonds
    }

    fn load_smiles() -> Result<Vec<String>, TestDataError> {
        let path = repo_root().join("tests/smiles.smi");
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut rows = Vec::new();
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            rows.push(trimmed.to_owned());
        }
        Ok(rows)
    }

    fn load_golden() -> Result<Vec<GoldenRecord>, TestDataError> {
        let path = repo_root().join("tests/golden/molblock_v2000_minimal.jsonl");
        ensure_golden_exists(&path);
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut rows = Vec::new();
        for (idx, line) in reader.lines().enumerate() {
            let line_no = idx + 1;
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let record = serde_json::from_str::<GoldenRecord>(&line)
                .map_err(|source| TestDataError::Json { line_no, source })?;
            rows.push(record);
        }
        Ok(rows)
    }

    fn load_kekulize_golden() -> Result<Vec<KekulizeGoldenRecord>, TestDataError> {
        let path = repo_root().join("tests/golden/molblock_v2000_kekulized.jsonl");
        ensure_kekulize_golden_exists(&path);
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut rows = Vec::new();
        for (idx, line) in reader.lines().enumerate() {
            let line_no = idx + 1;
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let record = serde_json::from_str::<KekulizeGoldenRecord>(&line)
                .map_err(|source| TestDataError::Json { line_no, source })?;
            rows.push(record);
        }
        Ok(rows)
    }

    fn ensure_golden_exists(golden_path: &Path) {
        if golden_path.exists() {
            return;
        }
        let repo = repo_root();
        let script = repo.join("tests/scripts/gen_rdkit_v2000_minimal_golden.py");
        let input = repo.join("tests/smiles.smi");

        let candidates = [
            std::env::var("COSMOLKIT_PYTHON").ok(),
            Some(repo.join(".venv/bin/python").display().to_string()),
            Some(String::from("python3")),
        ];

        let mut last_error = String::new();
        for candidate in candidates.iter().flatten() {
            let output = Command::new(candidate)
                .arg(&script)
                .arg("--input")
                .arg(&input)
                .arg("--output")
                .arg(golden_path)
                .output();
            match output {
                Ok(out) if out.status.success() => return,
                Ok(out) => {
                    last_error = format!(
                        "python={} exit={} stderr={}",
                        candidate,
                        out.status,
                        String::from_utf8_lossy(&out.stderr)
                    );
                }
                Err(err) => {
                    last_error = format!("python={} spawn error={}", candidate, err);
                }
            }
        }

        panic!(
            "golden file missing and auto-generation failed.\n\
             expected: {}\n\
             tried COSMOLKIT_PYTHON, .venv/bin/python, python3.\n\
             last error: {}\n\
             please run:\n\
             uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_v2000_minimal_golden.py --input tests/smiles.smi --output tests/golden/molblock_v2000_minimal.jsonl",
            golden_path.display(),
            last_error
        );
    }

    fn ensure_kekulize_golden_exists(golden_path: &Path) {
        if golden_path.exists() {
            return;
        }
        let repo = repo_root();
        let script = repo.join("tests/scripts/gen_rdkit_kekulize_molblock_golden.py");
        let input = repo.join("tests/smiles.smi");

        let candidates = [
            std::env::var("COSMOLKIT_PYTHON").ok(),
            Some(repo.join(".venv/bin/python").display().to_string()),
            Some(String::from("python3")),
        ];

        let mut last_error = String::new();
        for candidate in candidates.iter().flatten() {
            let output = Command::new(candidate)
                .arg(&script)
                .arg("--input")
                .arg(&input)
                .arg("--output")
                .arg(golden_path)
                .output();
            match output {
                Ok(out) if out.status.success() => return,
                Ok(out) => {
                    last_error = format!(
                        "python={} exit={} stderr={}",
                        candidate,
                        out.status,
                        String::from_utf8_lossy(&out.stderr)
                    );
                }
                Err(err) => {
                    last_error = format!("python={} spawn error={}", candidate, err);
                }
            }
        }

        panic!(
            "kekulize golden file missing and auto-generation failed.\n\
             expected: {}\n\
             tried COSMOLKIT_PYTHON, .venv/bin/python, python3.\n\
             last error: {}\n\
             please run:\n\
             uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_kekulize_molblock_golden.py --input tests/smiles.smi --output tests/golden/molblock_v2000_kekulized.jsonl",
            golden_path.display(),
            last_error
        );
    }

    #[test]
    fn molblock_v2000_golden_has_one_record_per_smiles() {
        let smiles = load_smiles().expect("read tests/smiles.smi");
        let golden = load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
        assert_eq!(
            golden.len(),
            smiles.len(),
            "golden rows must match input smiles rows"
        );

        for (idx, (record, input_smiles)) in golden.iter().zip(smiles.iter()).enumerate() {
            assert_eq!(
                record.smiles,
                *input_smiles,
                "smiles mismatch at row {}",
                idx + 1
            );

            if record.parse_ok {
                assert!(
                    record.parse_error.is_none(),
                    "parse_ok=true should not carry parse_error at row {}",
                    idx + 1
                );
            } else {
                assert!(
                    record.parse_error.is_some(),
                    "parse_ok=false should carry parse_error at row {}",
                    idx + 1
                );
            }

            if record.v2000_ok {
                assert!(
                    record.v2000_body.is_some(),
                    "v2000_ok=true requires v2000_body at row {}",
                    idx + 1
                );
                assert!(
                    record.v2000_error.is_none(),
                    "v2000_ok=true should not carry v2000_error at row {}",
                    idx + 1
                );
            } else {
                assert!(
                    record.v2000_body.is_none(),
                    "v2000_ok=false should not carry v2000_body at row {}",
                    idx + 1
                );
                assert!(
                    record.v2000_error.is_some(),
                    "v2000_ok=false should carry v2000_error at row {}",
                    idx + 1
                );
            }

            if record.v3000_ok {
                assert!(
                    record.v3000_body.is_some(),
                    "v3000_ok=true requires v3000_body at row {}",
                    idx + 1
                );
                assert!(
                    record.v3000_error.is_none(),
                    "v3000_ok=true should not carry v3000_error at row {}",
                    idx + 1
                );
            } else {
                assert!(
                    record.v3000_body.is_none(),
                    "v3000_ok=false should not carry v3000_body at row {}",
                    idx + 1
                );
                assert!(
                    record.v3000_error.is_some(),
                    "v3000_ok=false should carry v3000_error at row {}",
                    idx + 1
                );
            }
        }
    }

    #[test]
    fn molblock_v2000_body_matches_rdkit_coordinates_and_topology() {
        let golden = load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
        for (idx, record) in golden.iter().enumerate() {
            let parsed = Molecule::from_smiles(&record.smiles);
            if record.parse_ok {
                assert!(
                    parsed.is_ok(),
                    "parse should succeed at row {} ({})",
                    idx + 1,
                    record.smiles
                );
            } else {
                assert!(
                    parsed.is_err(),
                    "parse should fail at row {} ({})",
                    idx + 1,
                    record.smiles
                );
                continue;
            }
            let mol = parsed.expect("parse checked above");
            let ours = mol_to_v2000_block_minimal(&mol)
                .unwrap_or_else(|e| panic!("write failed at row {}: {}", idx + 1, e));
            let ours_body = body(&ours);

            if record.v2000_ok {
                let expected_v2000 = record
                    .v2000_body
                    .as_ref()
                    .expect("v2000_ok=true requires v2000_body");
                compare_against_expected(
                    &ours_body,
                    expected_v2000,
                    &mol,
                    &record.smiles,
                    idx + 1,
                    "v2000",
                );
            }
            if record.v3000_ok {
                let expected_v3000 = record
                    .v3000_body
                    .as_ref()
                    .expect("v3000_ok=true requires v3000_body");
                compare_against_expected(
                    &ours_body,
                    expected_v3000,
                    &mol,
                    &record.smiles,
                    idx + 1,
                    "v3000",
                );
            }
        }
    }

    #[test]
    fn molblock_kekulized_golden_has_one_record_per_smiles() {
        let smiles = load_smiles().expect("read tests/smiles.smi");
        let golden =
            load_kekulize_golden().expect("read tests/golden/molblock_v2000_kekulized.jsonl");
        assert_eq!(
            golden.len(),
            smiles.len(),
            "kekulize golden rows must match input smiles rows"
        );
        for (idx, (record, input_smiles)) in golden.iter().zip(smiles.iter()).enumerate() {
            assert_eq!(
                record.smiles,
                *input_smiles,
                "kekulize smiles mismatch at row {}",
                idx + 1
            );
            if record.parse_ok {
                assert!(
                    record.parse_error.is_none(),
                    "kekulize parse_ok=true should not carry parse_error at row {}",
                    idx + 1
                );
            } else {
                assert!(
                    record.parse_error.is_some(),
                    "kekulize parse_ok=false should carry parse_error at row {}",
                    idx + 1
                );
            }
            if record.kekulize_ok {
                assert!(
                    record.kekulize_error.is_none(),
                    "kekulize_ok=true should not carry kekulize_error at row {}",
                    idx + 1
                );
                if record.v2000_ok {
                    assert!(
                        record.v2000_body.is_some() && record.v2000_error.is_none(),
                        "kekulize v2000_ok=true requires body and no error at row {}",
                        idx + 1
                    );
                } else {
                    assert!(
                        record.v2000_body.is_none() && record.v2000_error.is_some(),
                        "kekulize v2000_ok=false requires error and no body at row {}",
                        idx + 1
                    );
                }
                if record.v3000_ok {
                    assert!(
                        record.v3000_body.is_some() && record.v3000_error.is_none(),
                        "kekulize v3000_ok=true requires body and no error at row {}",
                        idx + 1
                    );
                } else {
                    assert!(
                        record.v3000_body.is_none() && record.v3000_error.is_some(),
                        "kekulize v3000_ok=false requires error and no body at row {}",
                        idx + 1
                    );
                }
            } else {
                assert!(
                    record.kekulize_error.is_some(),
                    "kekulize_ok=false should carry kekulize_error at row {}",
                    idx + 1
                );
                assert!(
                    !record.v2000_ok && !record.v3000_ok,
                    "kekulize_ok=false should not have successful molblock variants at row {}",
                    idx + 1
                );
            }
        }
    }

    #[test]
    fn molblock_kekulized_topology_matches_rdkit_golden() {
        let golden =
            load_kekulize_golden().expect("read tests/golden/molblock_v2000_kekulized.jsonl");
        for (idx, record) in golden.iter().enumerate() {
            let parsed = Molecule::from_smiles(&record.smiles);
            if record.parse_ok {
                assert!(
                    parsed.is_ok(),
                    "parse should succeed at row {} ({})",
                    idx + 1,
                    record.smiles
                );
            } else {
                assert!(
                    parsed.is_err(),
                    "parse should fail at row {} ({})",
                    idx + 1,
                    record.smiles
                );
                continue;
            }
            let mut mol = parsed.expect("parse checked above");

            let ours_kek = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                crate::kekulize::kekulize_in_place(&mut mol)
            }));
            let kek_ok = matches!(ours_kek, Ok(Ok(())));

            if record.kekulize_ok {
                assert!(
                    kek_ok,
                    "kekulize should succeed at row {} ({})",
                    idx + 1,
                    record.smiles
                );
            } else {
                assert!(
                    !kek_ok,
                    "kekulize should fail at row {} ({})",
                    idx + 1,
                    record.smiles
                );
                continue;
            }

            let ours = mol_to_v2000_block_minimal(&mol)
                .unwrap_or_else(|e| panic!("write failed at row {}: {}", idx + 1, e));
            let ours_body = body(&ours);

            if record.v2000_ok {
                let expected_v2000 = record
                    .v2000_body
                    .as_ref()
                    .expect("v2000_ok=true requires v2000_body");
                assert!(
                    bonds_and_atoms_match_strict(
                        &normalize_signed_zero(&ours_body),
                        &normalize_signed_zero(expected_v2000)
                    ),
                    "kekulized bond block mismatch at row {} ({}) against v2000",
                    idx + 1,
                    record.smiles
                );
            }
            if record.v3000_ok {
                let expected_v3000 = record
                    .v3000_body
                    .as_ref()
                    .expect("v3000_ok=true requires v3000_body");
                assert!(
                    bonds_and_atoms_match_strict(
                        &normalize_signed_zero(&ours_body),
                        &normalize_signed_zero(expected_v3000)
                    ),
                    "kekulized bond block mismatch at row {} ({}) against v3000",
                    idx + 1,
                    record.smiles
                );
            }
        }
    }

    #[test]
    fn long_biaryl_ester_anchor_aligned_ring_rmsd_diagnostic() {
        let golden = load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
        let diagnostic_smiles = "CCCCCCCCOc1ccc(C(=O)COC(=O)c2ccc(CCCCCCC)cc2)cc1";
        let record = golden
            .iter()
            .find(|record| record.smiles == diagnostic_smiles)
            .expect("diagnostic molecule should exist in golden");
        assert!(
            record.v2000_ok,
            "diagnostic molecule should have v2000 golden body"
        );

        let mol = Molecule::from_smiles(&record.smiles).expect("parse diagnostic molecule");
        let ours_block = mol_to_v2000_block_minimal(&mol).expect("write diagnostic molecule");
        let ours = parse_v2000_body(&body(&ours_block)).expect("parse ours diagnostic body");
        let expected = parse_v2000_body(
            record
                .v2000_body
                .as_ref()
                .expect("diagnostic v2000_ok=true requires body"),
        )
        .expect("parse rdkit diagnostic body");

        // O=C-CO-C=O six-atom anchor: [O, C, C, O, C, O]
        let anchor = [17usize, 16, 18, 19, 20, 21];
        let ring_a = [10usize, 11, 12, 13, 33, 34];
        let ring_b = [20usize, 21, 22, 23, 31, 32];

        let src_anchor: Vec<(f64, f64)> = anchor.iter().map(|&i| ours.coords[i - 1]).collect();
        let tgt_anchor: Vec<(f64, f64)> = anchor.iter().map(|&i| expected.coords[i - 1]).collect();
        let tf = align_rigid_2d(&src_anchor, &tgt_anchor).expect("anchor rigid align");
        let anchor_rmsd = rmsd_for_indices(&ours.coords, &expected.coords, &anchor, tf);
        let ring_a_rmsd = rmsd_for_indices(&ours.coords, &expected.coords, &ring_a, tf);
        let ring_b_rmsd = rmsd_for_indices(&ours.coords, &expected.coords, &ring_b, tf);

        assert!(
            anchor_rmsd <= 1e-4 && ring_a_rmsd <= 1e-4 && ring_b_rmsd <= 1e-4,
            "long biaryl ester anchor-fit RMSD: anchor={:.6}, ring1={:.6}, ring2={:.6}",
            anchor_rmsd,
            ring_a_rmsd,
            ring_b_rmsd
        );
    }

    #[test]
    fn long_biaryl_ester_coordinate_diagnostic() {
        let golden = load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
        let diagnostic_smiles = "CCCCCCCCOc1ccc(C(=O)COC(=O)c2ccc(CCCCCCC)cc2)cc1";
        let record = golden
            .iter()
            .find(|record| record.smiles == diagnostic_smiles)
            .expect("diagnostic molecule should exist in golden");
        let mol = Molecule::from_smiles(&record.smiles).expect("parse diagnostic molecule");
        let ours_block = mol_to_v2000_block_minimal(&mol).expect("write diagnostic molecule");
        let ours = parse_v2000_body(&body(&ours_block)).expect("parse ours diagnostic body");
        let expected = parse_v2000_body(
            record
                .v2000_body
                .as_ref()
                .expect("diagnostic v2000_ok=true requires body"),
        )
        .expect("parse rdkit diagnostic body");

        assert_eq!(ours.atoms, expected.atoms, "diagnostic atom symbols differ");
        let mut ours_bonds = ours.bonds.clone();
        let mut expected_bonds = expected.bonds.clone();
        ours_bonds.sort_unstable();
        expected_bonds.sort_unstable();
        assert_eq!(ours_bonds, expected_bonds, "diagnostic bonds differ");

        let mut sse = 0.0f64;
        let mut sse_rot = 0.0f64;
        let mut max = (0usize, 0.0f64, (0.0f64, 0.0f64), (0.0f64, 0.0f64));
        let mut max_rot = (0usize, 0.0f64, (0.0f64, 0.0f64), (0.0f64, 0.0f64));
        for i in 0..ours.coords.len() {
            let dx = ours.coords[i].0 - expected.coords[i].0;
            let dy = ours.coords[i].1 - expected.coords[i].1;
            let d = (dx * dx + dy * dy).sqrt();
            sse += d * d;
            if d > max.1 {
                max = (i + 1, d, ours.coords[i], expected.coords[i]);
            }
            let rot = (-ours.coords[i].0, -ours.coords[i].1);
            let dxr = rot.0 - expected.coords[i].0;
            let dyr = rot.1 - expected.coords[i].1;
            let dr = (dxr * dxr + dyr * dyr).sqrt();
            sse_rot += dr * dr;
            if dr > max_rot.1 {
                max_rot = (i + 1, dr, rot, expected.coords[i]);
            }
        }
        let rmsd = (sse / ours.coords.len() as f64).sqrt();
        assert!(
            rmsd <= 1e-4,
            "long biaryl ester direct RMSD={:.6}; rot180 RMSD={:.6}; direct max atom {} d={:.6} ours=({:.4},{:.4}) rdkit=({:.4},{:.4}); rot max atom {} d={:.6} ours=({:.4},{:.4}) rdkit=({:.4},{:.4})",
            rmsd,
            (sse_rot / ours.coords.len() as f64).sqrt(),
            max.0,
            max.1,
            max.2.0,
            max.2.1,
            max.3.0,
            max.3.1,
            max_rot.0,
            max_rot.1,
            max_rot.2.0,
            max_rot.2.1,
            max_rot.3.0,
            max_rot.3.1
        );
    }

    #[test]
    #[ignore = "temporary row 71 layout diagnostic"]
    fn row71_precanonical_coordinate_diagnostic() {
        let smiles = "c1sccc1[C@@](C)(O)[C@H]2CCCC(C23C)=Cc4c(C3)cnn4-c5ccc(F)cc5";
        let mol = Molecule::from_smiles(smiles).expect("parse row 71 diagnostic molecule");
        let coords = super::rdkit_compute_initial_coords_strict(&mol)
            .expect("compute row 71 diagnostic coordinates");
        for (idx, (x, y)) in coords.iter().enumerate() {
            println!("{idx}: {x:.4} {y:.4}");
        }
    }

    #[test]
    #[ignore = "temporary row 71 ring-order diagnostic"]
    fn row71_ring_order_diagnostic() {
        let smiles = "c1sccc1[C@@](C)(O)[C@H]2CCCC(C23C)=Cc4c(C3)cnn4-c5ccc(F)cc5";
        let mol = Molecule::from_smiles(smiles).expect("parse row 71 diagnostic molecule");
        let n = mol.atoms.len();
        let mut adjacency = vec![Vec::<usize>::new(); n];
        for b in &mol.bonds {
            adjacency[b.begin_atom].push(b.end_atom);
            adjacency[b.end_atom].push(b.begin_atom);
        }
        for nbs in &mut adjacency {
            nbs.sort_unstable();
            nbs.dedup();
        }
        let comp: Vec<usize> = (0..n).collect();
        let mut deg_in_comp = vec![0usize; n];
        for &a in &comp {
            deg_in_comp[a] = adjacency[a].len();
        }
        let mut removed = std::collections::HashSet::<usize>::new();
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
                if removed.contains(&nb) {
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
        let ring_set: std::collections::HashSet<usize> = ring_atoms.iter().copied().collect();
        let ring_bond_count = mol
            .bonds
            .iter()
            .filter(|b| ring_set.contains(&b.begin_atom) && ring_set.contains(&b.end_atom))
            .count();
        let target_cycle_count = ring_bond_count + 1 - ring_atoms.len();
        let cycles = super::place_multiring_nonfused_component_rdkit_like;
        let _ = cycles;
        println!("target_cycle_count={target_cycle_count}");
        // keep this test as a compile-time anchor; use the existing failure output
        // for detailed coordinate diagnostics.
    }
}
