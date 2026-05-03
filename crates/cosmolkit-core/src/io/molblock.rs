use std::collections::HashMap;
use std::f64::consts::PI;

#[cfg(any())]
use std::cmp::Ordering;

use crate::{BondDirection, BondOrder, BondStereo, ChiralTag, Molecule};

unsafe extern "C" {
    #[link_name = "cos"]
    fn c_libm_cos(x: f64) -> f64;
    #[link_name = "sin"]
    fn c_libm_sin(x: f64) -> f64;
    #[link_name = "acos"]
    fn c_libm_acos(x: f64) -> f64;
    #[link_name = "sqrt"]
    fn c_libm_sqrt(x: f64) -> f64;
}

fn rdkit_cos(x: f64) -> f64 {
    unsafe { c_libm_cos(x) }
}

fn rdkit_sin(x: f64) -> f64 {
    unsafe { c_libm_sin(x) }
}

fn rdkit_acos(x: f64) -> f64 {
    unsafe { c_libm_acos(x) }
}

fn rdkit_sqrt(x: f64) -> f64 {
    unsafe { c_libm_sqrt(x) }
}

#[derive(Debug, thiserror::Error)]
pub enum MolWriteError {
    #[error("{0}")]
    UnsupportedSubset(&'static str),
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SdfFormat {
    Auto,
    V2000,
    V3000,
}

pub fn mol_to_v2000_2d_block(mol: &Molecule) -> Result<String, MolWriteError> {
    let coords = mol.coords_2d().ok_or(MolWriteError::UnsupportedSubset(
        "2D coordinates are required; call Molecule::compute_2d_coords() before writing 2D SDF",
    ))?;
    let coords = coords
        .iter()
        .map(|coord| (coord.x, coord.y, 0.0))
        .collect::<Vec<_>>();
    mol_to_v2000_block_with_coords(mol, &coords, "2D")
}

pub fn mol_to_v2000_3d_block(mol: &Molecule) -> Result<String, MolWriteError> {
    let coords = mol.coords_3d().ok_or(MolWriteError::UnsupportedSubset(
        "3D coordinates are required; parse a 3D molfile or add a 3D conformer before writing 3D SDF",
    ))?;
    let coords = coords
        .iter()
        .map(|coord| (coord.x, coord.y, coord.z))
        .collect::<Vec<_>>();
    mol_to_v2000_block_with_coords(mol, &coords, "3D")
}

pub fn mol_to_v2000_block(mol: &Molecule) -> Result<String, MolWriteError> {
    if mol.coords_2d().is_some() {
        mol_to_v2000_2d_block(mol)
    } else {
        mol_to_v2000_3d_block(mol)
    }
}

fn mol_to_v2000_block_with_coords(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    coord_label: &str,
) -> Result<String, MolWriteError> {
    // Mirror the golden generation path: tests/scripts first calls
    // RDDepict::Compute2DCoords() on the parsed molecule, then MolToMolBlock()
    // writes through a kekulized copy while preserving the existing conformer.
    let mut write_mol = mol.clone();
    if write_mol.bonds.iter().any(|b| b.is_aromatic) {
        crate::kekulize::kekulize_in_place(&mut write_mol, true)
            .map_err(|_| MolWriteError::UnsupportedSubset("kekulize before v2000 write failed"))?;
    }

    if coords.len() != write_mol.atoms.len() {
        return Err(MolWriteError::UnsupportedSubset(
            "coordinate count does not match atom count",
        ));
    }

    let mut out = String::new();
    out.push('\n');
    out.push_str(&format!("     COSMolKit      {coord_label}\n"));
    out.push('\n');
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        write_mol.atoms.len(),
        write_mol.bonds.len()
    ));

    for (idx, atom) in write_mol.atoms.iter().enumerate() {
        let (x, y, z) = coords[idx];
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0{:>3}{:>3}  0  0  0  0  0  0  0  0  0\n",
            x,
            y,
            z,
            atom_symbol(atom.atomic_num),
            0,
            0
        ));
    }

    let coords_2d = coords
        .iter()
        .map(|(x, y, _)| glam::DVec2::new(*x, *y))
        .collect::<Vec<_>>();
    let wedge = pick_bonds_to_wedge_rdkit_subset(&write_mol);
    for bond in &write_mol.bonds {
        let (dir_code, reverse) = if matches!(bond.order, BondOrder::Single) {
            let dir = determine_bond_wedge_state_rdkit_subset(bond, &wedge, &coords_2d, &write_mol);
            let reverse = wedge
                .get(&bond.index)
                .is_some_and(|from_atom| *from_atom != bond.begin_atom)
                && !matches!(dir, BondDirection::None);
            (bond_stereo_code_with_direction(bond, dir), reverse)
        } else {
            (bond_stereo_code(bond), false)
        };
        let (begin_atom, end_atom) = if reverse {
            (bond.end_atom + 1, bond.begin_atom + 1)
        } else {
            (bond.begin_atom + 1, bond.end_atom + 1)
        };
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}{:>3}\n",
            begin_atom,
            end_atom,
            bond_type_code(bond.order),
            dir_code
        ));
    }
    append_v2000_property_lines(&mut out, &write_mol);
    out.push_str("M  END\n");
    Ok(out)
}

pub fn mol_to_2d_sdf_record(mol: &Molecule, format: SdfFormat) -> Result<String, MolWriteError> {
    let mut block = match format {
        SdfFormat::Auto => {
            if mol.atoms.len() > 999
                || mol.bonds.len() > 999
                || mol
                    .bonds
                    .iter()
                    .any(|bond| matches!(bond.order, BondOrder::Dative))
            {
                mol_to_v3000_2d_block(mol)?
            } else {
                mol_to_v2000_2d_block(mol)?
            }
        }
        SdfFormat::V2000 => mol_to_v2000_2d_block(mol)?,
        SdfFormat::V3000 => mol_to_v3000_2d_block(mol)?,
    };
    block.push_str("$$$$\n");
    Ok(block)
}

pub fn mol_to_3d_sdf_record(mol: &Molecule, format: SdfFormat) -> Result<String, MolWriteError> {
    let mut block = match format {
        SdfFormat::Auto => {
            if mol.atoms.len() > 999
                || mol.bonds.len() > 999
                || mol
                    .bonds
                    .iter()
                    .any(|bond| matches!(bond.order, BondOrder::Dative))
            {
                mol_to_v3000_3d_block(mol)?
            } else {
                mol_to_v2000_3d_block(mol)?
            }
        }
        SdfFormat::V2000 => mol_to_v2000_3d_block(mol)?,
        SdfFormat::V3000 => mol_to_v3000_3d_block(mol)?,
    };
    block.push_str("$$$$\n");
    Ok(block)
}

pub fn mol_to_v3000_2d_block(mol: &Molecule) -> Result<String, MolWriteError> {
    let coords = mol.coords_2d().ok_or(MolWriteError::UnsupportedSubset(
        "2D coordinates are required; call Molecule::compute_2d_coords() before writing 2D SDF",
    ))?;
    let coords = coords
        .iter()
        .map(|coord| (coord.x, coord.y, 0.0))
        .collect::<Vec<_>>();
    mol_to_v3000_block_with_coords(mol, &coords, "2D")
}

pub fn mol_to_v3000_3d_block(mol: &Molecule) -> Result<String, MolWriteError> {
    let coords = mol.coords_3d().ok_or(MolWriteError::UnsupportedSubset(
        "3D coordinates are required; parse a 3D molfile or add a 3D conformer before writing 3D SDF",
    ))?;
    let coords = coords
        .iter()
        .map(|coord| (coord.x, coord.y, coord.z))
        .collect::<Vec<_>>();
    mol_to_v3000_block_with_coords(mol, &coords, "3D")
}

fn mol_to_v3000_block_with_coords(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    coord_label: &str,
) -> Result<String, MolWriteError> {
    if coords.len() != mol.atoms.len() {
        return Err(MolWriteError::UnsupportedSubset(
            "coordinate count does not match atom count",
        ));
    }

    let mut out = String::new();
    out.push('\n');
    out.push_str(&format!("     COSMolKit      {coord_label}\n"));
    out.push('\n');
    out.push_str("  0  0  0  0  0  0  0  0  0  0999 V3000\n");
    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {} 0 0 0\n",
        mol.atoms.len(),
        mol.bonds.len()
    ));
    out.push_str("M  V30 BEGIN ATOM\n");
    for (idx, atom) in mol.atoms.iter().enumerate() {
        let (x, y, z) = coords[idx];
        out.push_str(&format!(
            "M  V30 {} {} {:.6} {:.6} {:.6} 0",
            idx + 1,
            atom_symbol(atom.atomic_num),
            x,
            y,
            z
        ));
        if atom.formal_charge != 0 {
            out.push_str(&format!(" CHG={}", atom.formal_charge));
        }
        if let Some(isotope) = atom.isotope {
            out.push_str(&format!(" MASS={isotope}"));
        }
        for (key, value) in &atom.props {
            out.push_str(&format!(" {key}={value}"));
        }
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
    out.push_str("M  V30 BEGIN BOND\n");
    for (idx, bond) in mol.bonds.iter().enumerate() {
        out.push_str(&format!(
            "M  V30 {} {} {} {}",
            idx + 1,
            bond_type_code(bond.order),
            bond.begin_atom + 1,
            bond.end_atom + 1
        ));
        if let Some(cfg) = bond_v3000_cfg_code(bond) {
            out.push_str(&format!(" CFG={cfg}"));
        }
        out.push('\n');
    }
    out.push_str("M  V30 END BOND\n");
    out.push_str("M  V30 END CTAB\n");
    out.push_str("M  END\n");
    Ok(out)
}

pub(crate) fn compute_2d_coords(mol: &Molecule) -> Result<Vec<(f64, f64)>, MolWriteError> {
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
    let hybridizations = rdkit_hybridizations_for_depict(mol, &degree)?;
    let cip_ranks = mol
        .rdkit_legacy_stereo_atom_props(true)
        .into_iter()
        .enumerate()
        .map(|(idx, props)| props.cip_rank.unwrap_or(idx as i64))
        .collect::<Vec<_>>();

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

    let mut local_components: Vec<(usize, usize, Vec<(usize, (f64, f64))>)> = Vec::new();

    for comp in components {
        let k = comp.len();
        let comp_order = if k > 1 { 0usize } else { 1usize };
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

            let chosen = if cyclomatic <= 0 {
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
            };
            chosen
        };

        let local = local.ok_or(MolWriteError::UnsupportedSubset(
            "strict RDKit computeInitialCoords branch missing for this component; no heuristic fallback",
        ))?;
        let has_ring_or_stereo_seed = if k > 1 {
            let comp_set: std::collections::HashSet<usize> = comp.iter().copied().collect();
            let bond_count_in_comp = mol
                .bonds
                .iter()
                .filter(|b| comp_set.contains(&b.begin_atom) && comp_set.contains(&b.end_atom))
                .count();
            bond_count_in_comp + 1 > k
                || mol.bonds.iter().any(|bond| {
                    comp_set.contains(&bond.begin_atom)
                        && comp_set.contains(&bond.end_atom)
                        && matches!(bond.order, BondOrder::Double)
                        && matches!(bond.stereo, BondStereo::Cis | BondStereo::Trans)
                        && bond.stereo_atoms.len() == 2
                })
        } else {
            false
        };
        let rdkit_order = if has_ring_or_stereo_seed {
            0usize
        } else {
            1usize
        };
        let rdkit_rank = comp
            .iter()
            .copied()
            .map(|idx| atom_depict_rank(mol.atoms[idx].atomic_num, degree[idx]) * n + idx)
            .min()
            .unwrap_or(usize::MAX);
        local_components.push((rdkit_order, rdkit_rank.max(comp_order), local));
    }
    local_components.sort_by_key(|(rdkit_order, rank, component)| {
        (
            *rdkit_order,
            *rank,
            std::cmp::Reverse(component.len()),
            component[0].0,
        )
    });
    for (_, _, component) in &mut local_components {
        canonicalize_component_rdkit_like(component);
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
    let (mut xmax, xmin, mut ymax, ymin) = box_of(&local_components[0].2);
    for &(idx, (x, y)) in &local_components[0].2 {
        out[idx] = (x, y);
    }
    for (_, _, component) in local_components.iter().skip(1) {
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
    let inv_n = 1.0 / n;
    cx *= inv_n;
    cy *= inv_n;

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

#[derive(Clone, Debug)]
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

fn rdkit_ring_radius(ring_size: usize, bond_len: f64) -> f64 {
    let ang = 2.0 * PI / ring_size as f64;
    bond_len / rdkit_sqrt(2.0 * (1.0 - rdkit_cos(ang)))
}

fn n_outer_electrons_rdkit(atomic_num: u8) -> Option<i32> {
    // RDKit PeriodicTable::getNouterElecs(), values mirrored from chem-core's
    // RDKit 2026.03.1 valence support table.
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
    // map-number variation; atom maps are not currently stored in the graph.
    const MASS_BITS: i64 = 10;
    const MAX_MASS: i64 = 1 << MASS_BITS;

    fn most_common_isotope(atomic_num: u8) -> i32 {
        match atomic_num {
            1 => 1,
            5 => 11,
            6 => 12,
            7 => 14,
            8 => 16,
            9 => 19,
            11 => 23,
            12 => 24,
            13 => 27,
            14 => 28,
            15 => 31,
            16 => 32,
            17 => 35,
            29 => 63,
            34 => 80,
            35 => 79,
            53 => 127,
            _ => 0,
        }
    }

    mol.atoms
        .iter()
        .map(|atom| {
            let mut invariant = i64::from(atom.atomic_num % 128);
            let mut mass = 0i64;
            if let Some(isotope) = atom.isotope {
                mass = i64::from(isotope) - i64::from(most_common_isotope(atom.atomic_num));
                if mass >= 0 {
                    mass += 1;
                }
            }
            mass += MAX_MASS / 2;
            if mass < 0 {
                mass = 0;
            } else {
                mass %= MAX_MASS;
            }
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

fn rdkit_cip_ranks_from_invariants(
    mol: &Molecule,
    invars: Vec<i64>,
    seed_with_invars: bool,
) -> Vec<i64> {
    let n = mol.atoms.len();
    let assignment = crate::assign_valence(mol, crate::ValenceModel::RdkitLike).ok();
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
        if seed_with_invars {
            cip_entries[i][0] = invars[i];
        } else {
            cip_entries[i][0] = i64::from(mol.atoms[i].atomic_num);
            cip_entries[i].push(ranks[i]);
        }
    }
    for item in &mut sortable {
        item.0 = cip_entries[item.1].clone();
    }
    let cip_rank_index = if seed_with_invars { 1 } else { 2 };

    const K_MAX_BONDS: usize = 16;
    let mut bond_features = vec![(0i64, 0usize); n * K_MAX_BONDS];
    let mut num_neighbors = vec![0usize; n];
    for atom_idx in 0..n {
        let mut index_offset = atom_idx * K_MAX_BONDS;
        for bond in &mol.bonds {
            let nbr_idx = if bond.begin_atom == atom_idx {
                bond.end_atom
            } else if bond.end_atom == atom_idx {
                bond.begin_atom
            } else {
                continue;
            };
            let count = if matches!(bond.order, BondOrder::Double) {
                let nbr = &mol.atoms[nbr_idx];
                let nbr_degree = mol
                    .bonds
                    .iter()
                    .filter(|other| other.begin_atom == nbr_idx || other.end_atom == nbr_idx)
                    .count();
                if nbr.atomic_num == 15 && matches!(nbr_degree, 3 | 4) {
                    1
                } else {
                    rdkit_twice_bond_type(bond.order)
                }
            } else {
                rdkit_twice_bond_type(bond.order)
            };
            if index_offset < (atom_idx + 1) * K_MAX_BONDS {
                bond_features[index_offset] = (count, nbr_idx);
            }
            num_neighbors[atom_idx] += 1;
            index_offset += 1;
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
            let index_offset = K_MAX_BONDS * index;
            let feature_len = num_neighbors[index].min(K_MAX_BONDS);
            let mut features = bond_features[index_offset..index_offset + feature_len].to_vec();
            if num_neighbors[index] > 1 {
                features.sort_by(|a, b| ranks[b.1].cmp(&ranks[a.1]));
            }
            for (count, nbr_idx) in features {
                for _ in 0..count {
                    cip_entries[index].push(ranks[nbr_idx] + 1);
                }
            }
            let total_hs = if let Some(assignment) = &assignment {
                mol.atoms[index].explicit_hydrogens as usize
                    + assignment.implicit_hydrogens[index] as usize
            } else {
                mol.atoms[index].explicit_hydrogens as usize
            };
            cip_entries[index].extend(std::iter::repeat_n(0, total_hs));
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
                cip_entries[i].resize(cip_rank_index + 1, 0);
                cip_entries[i][cip_rank_index] = ranks[i];
            }
            for item in &mut sortable {
                item.0 = cip_entries[item.1].clone();
            }
        }

        num_its += 1;
    }

    ranks
}

#[cfg(any())]
const ATNUM_CLASS_OFFSET: usize = 10000;

#[cfg(any())]
#[derive(Debug, Clone)]
struct ChiralBondHolder {
    bond_type: u8,
    bond_stereo: u8,
    nbr_sym_class: usize,
    nbr_idx: usize,
}

#[cfg(any())]
#[derive(Debug, Clone)]
struct ChiralCanonAtom {
    index: usize,
    degree: usize,
    atomic_num: u8,
    isotope: u16,
    cip_code_rank: u8,
    nbr_ids: Vec<usize>,
    bonds: Vec<ChiralBondHolder>,
}

#[cfg(any())]
fn rdkit_chiral_bond_stereo_rank(stereo: BondStereo) -> u8 {
    match stereo {
        BondStereo::Cis => 1,
        BondStereo::Trans => 2,
        BondStereo::None | BondStereo::Any => 0,
    }
}

#[cfg(any())]
fn rdkit_chiral_bond_repetitions(
    order: BondOrder,
    neighbor_atomic_num: u8,
    neighbor_degree: usize,
) -> usize {
    if matches!(order, BondOrder::Double)
        && neighbor_atomic_num == 15
        && matches!(neighbor_degree, 3 | 4)
    {
        return 1;
    }
    match order {
        BondOrder::Null => 0,
        BondOrder::Single | BondOrder::Dative => 2,
        BondOrder::Aromatic => 3,
        BondOrder::Double => 4,
        BondOrder::Triple => 6,
        BondOrder::Quadruple => 8,
    }
}

#[cfg(any())]
fn chiral_bond_compare(lhs: &ChiralBondHolder, rhs: &ChiralBondHolder, div: usize) -> Ordering {
    lhs.bond_type
        .cmp(&rhs.bond_type)
        .then_with(|| lhs.bond_stereo.cmp(&rhs.bond_stereo))
        .then_with(|| (lhs.nbr_sym_class / div).cmp(&(rhs.nbr_sym_class / div)))
}

#[cfg(any())]
fn chiral_bond_greater(lhs: &ChiralBondHolder, rhs: &ChiralBondHolder) -> Ordering {
    chiral_bond_compare(rhs, lhs, 1)
}

#[cfg(any())]
fn init_chiral_canon_atoms(mol: &Molecule) -> Vec<ChiralCanonAtom> {
    let assignment = crate::assign_valence(mol, crate::ValenceModel::RdkitLike).ok();
    let mut graph_degree = vec![0usize; mol.atoms.len()];
    for bond in &mol.bonds {
        graph_degree[bond.begin_atom] += 1;
        graph_degree[bond.end_atom] += 1;
    }

    let mut atoms: Vec<ChiralCanonAtom> = mol
        .atoms
        .iter()
        .map(|atom| ChiralCanonAtom {
            index: atom.index,
            degree: graph_degree[atom.index],
            atomic_num: atom.atomic_num,
            isotope: atom.isotope.unwrap_or(0),
            cip_code_rank: 0,
            nbr_ids: Vec::with_capacity(graph_degree[atom.index]),
            bonds: Vec::new(),
        })
        .collect();

    for bond in &mol.bonds {
        atoms[bond.begin_atom].nbr_ids.push(bond.end_atom);
        atoms[bond.end_atom].nbr_ids.push(bond.begin_atom);

        for (atom_idx, nbr_idx) in [
            (bond.begin_atom, bond.end_atom),
            (bond.end_atom, bond.begin_atom),
        ] {
            let nbr = &mol.atoms[nbr_idx];
            let reps =
                rdkit_chiral_bond_repetitions(bond.order, nbr.atomic_num, graph_degree[nbr_idx]);
            let holder = ChiralBondHolder {
                bond_type: 1,
                bond_stereo: rdkit_chiral_bond_stereo_rank(bond.stereo),
                nbr_sym_class: usize::from(nbr.atomic_num) * ATNUM_CLASS_OFFSET + nbr_idx + 1,
                nbr_idx,
            };
            atoms[atom_idx]
                .bonds
                .extend(std::iter::repeat_n(holder, reps));
        }
    }

    for atom in &mut atoms {
        let total_hs = atom.degree.checked_sub(atom.degree).unwrap_or(0)
            + mol.atoms[atom.index].explicit_hydrogens as usize
            + assignment
                .as_ref()
                .map(|valence| valence.implicit_hydrogens[atom.index] as usize)
                .unwrap_or(0);
        for _ in 0..total_hs {
            atom.bonds.push(ChiralBondHolder {
                bond_type: 1,
                bond_stereo: 0,
                nbr_sym_class: ATNUM_CLASS_OFFSET,
                nbr_idx: ATNUM_CLASS_OFFSET,
            });
            atom.bonds.push(ChiralBondHolder {
                bond_type: 1,
                bond_stereo: 0,
                nbr_sym_class: ATNUM_CLASS_OFFSET,
                nbr_idx: ATNUM_CLASS_OFFSET,
            });
        }
        atom.bonds.sort_by(chiral_bond_greater);
    }

    atoms
}

#[cfg(any())]
fn update_chiral_atom_neighborhood(atoms: &mut [ChiralCanonAtom], atom_idx: usize) {
    let updated: Vec<(usize, usize)> = atoms[atom_idx]
        .bonds
        .iter()
        .map(|bond| {
            if bond.nbr_idx == ATNUM_CLASS_OFFSET {
                (ATNUM_CLASS_OFFSET, ATNUM_CLASS_OFFSET)
            } else {
                let nbr = &atoms[bond.nbr_idx];
                (
                    bond.nbr_idx,
                    usize::from(nbr.atomic_num) * ATNUM_CLASS_OFFSET + nbr.index + 1,
                )
            }
        })
        .collect();
    for (bond, (_, sym_class)) in atoms[atom_idx].bonds.iter_mut().zip(updated) {
        bond.nbr_sym_class = sym_class;
    }
    atoms[atom_idx].bonds.sort_by(chiral_bond_greater);
}

#[cfg(any())]
fn chiral_atom_compare(atoms: &mut [ChiralCanonAtom], i: usize, j: usize) -> Ordering {
    let base = atoms[i]
        .index
        .cmp(&atoms[j].index)
        .then_with(|| atoms[i].atomic_num.cmp(&atoms[j].atomic_num))
        .then_with(|| atoms[i].isotope.cmp(&atoms[j].isotope))
        .then_with(|| atoms[i].cip_code_rank.cmp(&atoms[j].cip_code_rank));
    if base != Ordering::Equal {
        return base;
    }

    update_chiral_atom_neighborhood(atoms, i);
    update_chiral_atom_neighborhood(atoms, j);

    let shared = atoms[i].bonds.len().min(atoms[j].bonds.len());
    for k in 0..shared {
        let cmp = chiral_bond_compare(&atoms[i].bonds[k], &atoms[j].bonds[k], ATNUM_CLASS_OFFSET);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    for k in 0..shared {
        let cmp = chiral_bond_compare(&atoms[i].bonds[k], &atoms[j].bonds[k], 1);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    atoms[i].bonds.len().cmp(&atoms[j].bonds.len())
}

#[cfg(any())]
fn chiral_atom_compare_i32(atoms: &mut [ChiralCanonAtom], i: usize, j: usize) -> i32 {
    match chiral_atom_compare(atoms, i, j) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

#[cfg(any())]
fn hanoi_sort_chiral_partition(
    order: &mut [usize],
    count: &mut [usize],
    changed: &mut [bool],
    atoms: &mut [ChiralCanonAtom],
) {
    fn hanoi(
        base: &mut [usize],
        temp: &mut [usize],
        count: &mut [usize],
        changed: &mut [bool],
        atoms: &mut [ChiralCanonAtom],
    ) -> bool {
        let nel = base.len();
        if nel == 1 {
            count[base[0]] = 1;
            return false;
        }
        if nel == 2 {
            let n1 = base[0];
            let n2 = base[1];
            let stat = if changed[n1] || changed[n2] {
                chiral_atom_compare_i32(atoms, n1, n2)
            } else {
                0
            };
            if stat == 0 {
                count[n1] = 2;
                count[n2] = 0;
                return false;
            }
            count[n1] = 1;
            count[n2] = 1;
            if stat > 0 {
                base[0] = n2;
                base[1] = n1;
            }
            return false;
        }

        let n1_len = nel / 2;
        let n2_len = nel - n1_len;
        let (b1, b2) = base.split_at_mut(n1_len);
        let (t1, t2) = temp.split_at_mut(n1_len);

        let left_temp = hanoi(b1, t1, count, changed, atoms);
        let right_temp = hanoi(b2, t2, count, changed, atoms);
        let s1 = if left_temp { t1.to_vec() } else { b1.to_vec() };
        let s2 = if right_temp { t2.to_vec() } else { b2.to_vec() };
        let result = !left_temp;
        let mut merged = Vec::with_capacity(nel);
        let mut i1 = 0usize;
        let mut i2 = 0usize;

        loop {
            let a = s1[i1];
            let b = s2[i2];
            let stat = if changed[a] || changed[b] {
                chiral_atom_compare_i32(atoms, a, b)
            } else {
                0
            };
            let len1 = count[a];
            let len2 = count[b];

            if stat == 0 {
                count[a] = len1 + len2;
                count[b] = 0;
                merged.extend_from_slice(&s1[i1..i1 + len1]);
                i1 += len1;
                if i1 == n1_len {
                    merged.extend_from_slice(&s2[i2..]);
                    break;
                }
                merged.extend_from_slice(&s2[i2..i2 + len2]);
                i2 += len2;
                if i2 == n2_len {
                    merged.extend_from_slice(&s1[i1..]);
                    break;
                }
            } else if stat < 0 {
                merged.extend_from_slice(&s1[i1..i1 + len1]);
                i1 += len1;
                if i1 == n1_len {
                    merged.extend_from_slice(&s2[i2..]);
                    break;
                }
            } else {
                merged.extend_from_slice(&s2[i2..i2 + len2]);
                i2 += len2;
                if i2 == n2_len {
                    merged.extend_from_slice(&s1[i1..]);
                    break;
                }
            }
        }

        if result {
            temp.copy_from_slice(&merged);
        } else {
            base.copy_from_slice(&merged);
        }
        result
    }

    let mut temp = vec![0usize; order.len()];
    if hanoi(order, &mut temp, count, changed, atoms) {
        order.copy_from_slice(&temp);
    }
}

#[cfg(any())]
fn activate_chiral_partitions(
    order: &[usize],
    count: &[usize],
    next: &mut [isize],
    changed: &mut [bool],
) -> isize {
    next.fill(-2);
    let mut activeset = -1isize;
    let mut i = 0usize;
    while i < order.len() {
        let atom_idx = order[i];
        if count[atom_idx] > 1 {
            next[atom_idx] = activeset;
            activeset = atom_idx as isize;
            i += count[atom_idx];
        } else {
            i += 1;
        }
    }
    for &atom_idx in order {
        changed[atom_idx] = true;
    }
    activeset
}

#[cfg(any())]
fn refine_chiral_partitions(
    atoms: &mut [ChiralCanonAtom],
    order: &mut [usize],
    count: &mut [usize],
    next: &mut [isize],
    changed: &mut [bool],
    touched: &mut [bool],
    mut activeset: isize,
) {
    while activeset != -1 {
        let partition = activeset as usize;
        activeset = next[partition];
        next[partition] = -2;

        let len = count[partition];
        let offset = atoms[partition].index;
        let start = &mut order[offset..offset + len];
        hanoi_sort_chiral_partition(start, count, changed, atoms);

        for &atom_idx in start.iter() {
            changed[atom_idx] = false;
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
            for nbr in atoms[index].nbr_ids.clone() {
                changed[nbr] = true;
            }
            i += 1;
        }

        index = start[0];
        let mut i = count[index];
        while i < len {
            index = start[i];
            for nbr in atoms[index].nbr_ids.clone() {
                touched[atoms[nbr].index] = true;
            }
            i += 1;
        }

        for touched_idx in 0..atoms.len() {
            if touched[touched_idx] {
                let partition = order[touched_idx];
                if count[partition] > 1 && next[partition] == -2 {
                    next[partition] = activeset;
                    activeset = partition as isize;
                }
                touched[touched_idx] = false;
            }
        }
    }
}

#[cfg(any())]
fn rdkit_chiral_rank_mol_atoms(mol: &Molecule) -> Vec<i64> {
    // Source mapping: RDKit 2026.03.1
    //   Code/GraphMol/Chirality.cpp::assignAtomCIPRanks()
    //   Code/GraphMol/new_canon.cpp::Canon::chiralRankMolAtoms()
    //   Code/GraphMol/new_canon.cpp::detail::getChiralBonds()
    //   Code/GraphMol/new_canon.h::ChiralAtomCompareFunctor
    let n = mol.atoms.len();
    if n == 0 {
        return Vec::new();
    }

    let mut atoms = init_chiral_canon_atoms(mol);
    let mut order: Vec<usize> = (0..n).collect();
    let mut count = vec![0usize; n];
    let mut next = vec![-2isize; n];
    let mut changed = vec![true; n];
    let mut touched = vec![false; n];

    for atom in &mut atoms {
        atom.index = 0;
    }
    count[0] = n;

    let activeset = activate_chiral_partitions(&order, &count, &mut next, &mut changed);
    refine_chiral_partitions(
        &mut atoms,
        &mut order,
        &mut count,
        &mut next,
        &mut changed,
        &mut touched,
        activeset,
    );

    let mut ranks = vec![0i64; n];
    for atom_idx in order {
        ranks[atom_idx] = atoms[atom_idx].index as i64;
    }
    ranks
}

pub(crate) fn rdkit_cip_ranks_for_depict(mol: &Molecule) -> Vec<i64> {
    rdkit_cip_ranks_from_invariants(mol, rdkit_cip_invariants(mol), false)
}

pub(crate) fn rdkit_cip_reranks_with_legacy_stereo(
    mol: &Molecule,
    ranks: &[i64],
    cip_codes: &[Option<String>],
) -> Vec<i64> {
    let mut factor = 100i64;
    while factor < mol.atoms.len() as i64 {
        factor *= 10;
    }
    let mut invars = vec![0i64; mol.atoms.len()];
    for atom in &mol.atoms {
        let mut invariant = ranks[atom.index] * factor;
        if let Some(code) = &cip_codes[atom.index] {
            if code == "S" {
                invariant += 10;
            } else if code == "R" {
                invariant += 20;
            }
        }
        for bond in &mol.bonds {
            if (bond.begin_atom == atom.index || bond.end_atom == atom.index)
                && matches!(bond.order, BondOrder::Double)
            {
                invariant += match bond.stereo {
                    BondStereo::Trans => 1,
                    BondStereo::Cis => 2,
                    BondStereo::None | BondStereo::Any => 0,
                };
            }
        }
        invars[atom.index] = invariant;
    }
    rdkit_cip_ranks_from_invariants(mol, invars, true)
}

fn rdkit_rank_atoms_by_rank(
    mol: &Molecule,
    atoms: &mut [usize],
    degree: &[usize],
    chiral_atom_ranks: Option<&[i64]>,
) {
    atoms.sort_by_key(|&idx| {
        let rank = mol.atoms[idx]
            .rdkit_cip_rank
            .or_else(|| chiral_atom_ranks.and_then(|ranks| ranks.get(idx).copied()))
            .unwrap_or_else(|| {
                (atom_depict_rank(mol.atoms[idx].atomic_num, degree[idx]) * mol.atoms.len() + idx)
                    as i64
            });
        (rank, idx as i64)
    });
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
        rdkit_rank_atoms_by_rank(mol, &mut out, degree, Some(cip_ranks));
        return out;
    }

    rdkit_rank_atoms_by_rank(mol, &mut thold, degree, Some(cip_ranks));

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

fn transform2d_mul3(lhs: [f64; 9], rhs: [f64; 9]) -> [f64; 9] {
    let mut out = [0.0; 9];
    for i in 0..3 {
        let id_a = i * 3;
        let id_c = id_a;
        for j in 0..3 {
            let id_ct = id_c + j;
            out[id_ct] = 0.0;
            for k in 0..3 {
                let id_at = id_a + k;
                let id_b = k * 3 + j;
                out[id_ct] += lhs[id_at] * rhs[id_b];
            }
        }
    }
    out
}

fn transform2d_to_affine(data: [f64; 9]) -> [f64; 6] {
    [data[0], data[1], data[2], data[3], data[4], data[5]]
}

fn transform2d_set_transform_center_angle(pt: (f64, f64), angle: f64) -> [f64; 6] {
    let trans1 = [1.0, 0.0, -pt.0, 0.0, 1.0, -pt.1, 0.0, 0.0, 1.0];
    let rot = [
        rdkit_cos(angle),
        -rdkit_sin(angle),
        0.0,
        rdkit_sin(angle),
        rdkit_cos(angle),
        0.0,
        0.0,
        0.0,
        1.0,
    ];
    let this = transform2d_mul3(rot, trans1);
    let trans2 = [1.0, 0.0, pt.0, 0.0, 1.0, pt.1, 0.0, 0.0, 1.0];
    transform2d_to_affine(transform2d_mul3(trans2, this))
}

fn rotate(v: (f64, f64), angle: f64) -> (f64, f64) {
    transform2d_point(v, transform2d_set_transform_center_angle((0.0, 0.0), angle))
}

fn rotate_around(p: (f64, f64), center: (f64, f64), angle: f64) -> (f64, f64) {
    transform2d_point(p, transform2d_set_transform_center_angle(center, angle))
}

fn transform2d_point(pt: (f64, f64), data: [f64; 6]) -> (f64, f64) {
    (
        data[0] * pt.0 + data[1] * pt.1 + data[2],
        data[3] * pt.0 + data[4] * pt.1 + data[5],
    )
}

fn transform2d_set_transform_two_point(
    ref1: (f64, f64),
    ref2: (f64, f64),
    pt1: (f64, f64),
    pt2: (f64, f64),
) -> [f64; 6] {
    let rvec = (ref2.0 - ref1.0, ref2.1 - ref1.1);
    let pvec = (pt2.0 - pt1.0, pt2.1 - pt1.1);
    let dp = rvec.0 * pvec.0 + rvec.1 * pvec.1;
    let lp = norm(rvec) * norm(pvec);
    if lp <= 0.0 {
        return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    }
    let mut cval = dp / lp;
    cval = cval.clamp(-1.0, 1.0);
    let mut ang = rdkit_acos(cval);
    let cross = pvec.0 * rvec.1 - pvec.1 * rvec.0;
    if cross < 0.0 {
        ang *= -1.0;
    }
    let mut data = [
        rdkit_cos(ang),
        -rdkit_sin(ang),
        0.0,
        rdkit_sin(ang),
        rdkit_cos(ang),
        0.0,
    ];
    let npt1 = transform2d_point(pt1, data);
    data[2] = ref1.0 - npt1.0;
    data[5] = ref1.1 - npt1.1;
    data
}

fn rdkit_reflect_point(point: (f64, f64), loc1: (f64, f64), loc2: (f64, f64)) -> (f64, f64) {
    let org = (0.0, 0.0);
    let xaxis = (1.0, 0.0);
    let mut cent = (loc1.0 + loc2.0, loc1.1 + loc2.1);
    cent.0 *= 0.5;
    cent.1 *= 0.5;

    let trans = transform2d_set_transform_two_point(org, xaxis, cent, loc1);
    let itrans = transform2d_set_transform_two_point(cent, loc1, org, xaxis);

    let mut res = transform2d_point(point, trans);
    res.1 = -res.1;
    transform2d_point(res, itrans)
}

fn norm(v: (f64, f64)) -> f64 {
    rdkit_sqrt(v.0 * v.0 + v.1 * v.1)
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
        rdkit_acos(c)
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
                rdkit_rank_atoms_by_rank(mol, &mut neighs, degree, Some(cip_ranks));
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
    remove_collisions_bond_flip_like(mol, comp, adjacency, &mut out);
    remove_collisions_open_angles_like(mol, comp, adjacency, &mut out);
    remove_collisions_shorten_bonds_like(mol, comp, adjacency, &mut out);
    Some(out)
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
            for bond in &mol.bonds {
                let v = if bond.begin_atom == u {
                    bond.end_atom
                } else if bond.end_atom == u {
                    bond.begin_atom
                } else {
                    continue;
                };
                if comp_set.contains(&v) {
                    if let std::collections::btree_map::Entry::Vacant(e) = prev.entry(v) {
                        e.insert(u);
                        q.push_back(v);
                    }
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
                if matches!(
                    mol.bonds[bond_idx].stereo,
                    BondStereo::None | BondStereo::Any
                ) && !is_ring_bond(bond_idx)
                {
                    res.push(bond_idx);
                }
            }
            pid = aid;
        }
        res
    };

    let graph_distance_matrix = || -> Vec<Vec<f64>> {
        let n = mol.atoms.len();
        let mut out = vec![vec![f64::INFINITY; n]; n];
        for start in 0..n {
            let mut q = std::collections::VecDeque::<usize>::new();
            out[start][start] = 0.0;
            q.push_back(start);
            while let Some(cur) = q.pop_front() {
                let next_dist = out[start][cur] + 1.0;
                for &nb in &adjacency[cur] {
                    if !comp_set.contains(&nb) {
                        continue;
                    }
                    if out[start][nb].is_infinite() {
                        out[start][nb] = next_dist;
                        q.push_back(nb);
                    }
                }
            }
        }
        out
    };
    let dmat = graph_distance_matrix();

    let reflect_point = rdkit_reflect_point;

    let find_collisions =
        |pos: &std::collections::BTreeMap<usize, (f64, f64)>| -> (Vec<(usize, usize)>, f64) {
            let mut collisions = Vec::new();
            let mut density = std::collections::BTreeMap::<usize, f64>::new();
            let mut atoms = comp.to_vec();
            atoms.sort_unstable();
            let col_thres2 = COLLISION_THRES * COLLISION_THRES;
            let bond_thres2 = 0.50f64 * 0.50f64;
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
            for (bid1, b1) in mol.bonds.iter().enumerate() {
                let beg1 = b1.begin_atom;
                let end1 = b1.end_atom;
                if !pos.contains_key(&beg1) || !pos.contains_key(&end1) {
                    continue;
                }
                let p_beg1 = pos[&beg1];
                let p_end1 = pos[&end1];
                let v1 = (p_end1.0 - p_beg1.0, p_end1.1 - p_beg1.1);
                let avg1 = ((p_end1.0 + p_beg1.0) * 0.5, (p_end1.1 + p_beg1.1) * 0.5);
                for (bid2, b2) in mol.bonds.iter().enumerate() {
                    if bid2 <= bid1 {
                        continue;
                    }
                    let beg2 = b2.begin_atom;
                    let end2 = b2.end_atom;
                    if !pos.contains_key(&beg2) || !pos.contains_key(&end2) {
                        continue;
                    }
                    let p_beg2 = pos[&beg2];
                    let p_end2 = pos[&end2];
                    let avg2 = ((p_end2.0 + p_beg2.0) * 0.5, (p_end2.1 + p_beg2.1) * 0.5);
                    let avg_delta = (avg2.0 - avg1.0, avg2.1 - avg1.1);
                    let avg_d2 = avg_delta.0 * avg_delta.0 + avg_delta.1 * avg_delta.1;
                    if avg_d2 < 0.5 && avg_d2 < bond_thres2 {
                        let v2 = (p_beg2.0 - p_beg1.0, p_beg2.1 - p_beg1.1);
                        let v3 = (p_end2.0 - p_beg1.0, p_end2.1 - p_beg1.1);
                        let cross_v1_v2 = v1.0 * v2.1 - v1.1 * v2.0;
                        let cross_v1_v3 = v1.0 * v3.1 - v1.1 * v3.0;
                        if cross_v1_v2 * cross_v1_v3 < -1.0e-6 {
                            let candidates = [
                                (dmat[beg1][beg2], (beg1, beg2)),
                                (dmat[beg1][end2], (beg1, end2)),
                                (dmat[end1][beg2], (end1, beg2)),
                                (dmat[end1][end2], (end1, end2)),
                            ];
                            let mut best = candidates[0];
                            for candidate in candidates.iter().skip(1) {
                                if candidate.0 < best.0 {
                                    best = *candidate;
                                }
                            }
                            collisions.push(best.1);
                        }
                    }
                }
            }
            let total_density = atoms.iter().fold(0.0, |accum, aid| {
                density.get(aid).copied().unwrap_or(0.0) + accum
            });
            (collisions, total_density)
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
        let mut atoms: Vec<usize> = pos.keys().copied().collect();
        atoms.retain(|aid| comp_set.contains(aid));
        for aid in atoms {
            let in_end_side = end_side_set.contains(&aid);
            if end_side_flip ^ !in_end_side {
                let p = pos[&aid];
                let beg_loc = pos[&beg];
                let end_loc = pos[&end];
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

            flip_about_bond(&mut pos, ri, true);
            let _ = find_collisions(&pos);
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

    let reflect_point = rdkit_reflect_point;

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
            rdkit_rank_atoms_by_rank(mol, &mut heavy, degree, Some(&cip_ranks));
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
        let chiral_atom_ranks;
        let rank_fallback = if mol.atoms.iter().any(|atom| atom.atomic_num == 1) {
            chiral_atom_ranks = rdkit_cip_ranks_for_depict(mol);
            Some(chiral_atom_ranks.as_slice())
        } else {
            None
        };
        let attach = frag.attach_pts.make_contiguous();
        rdkit_rank_atoms_by_rank(mol, attach, degree, rank_fallback);
    }

    fn find_num_neigh_frag(frag: &DepictFrag, pt: (f64, f64), radius: f64) -> i32 {
        frag.atoms
            .values()
            .filter(|st| {
                let dx = st.loc.0 - pt.0;
                let dy = st.loc.1 - pt.1;
                rdkit_sqrt(dx * dx + dy * dy) < radius
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
        curr_vec.0 *= 1.5;
        curr_vec.1 *= 1.5;
        curr_vec.0 += ref_state.loc.0;
        curr_vec.1 += ref_state.loc.1;
        let curr_loc = curr_vec;

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
        let trans = transform2d_set_transform_two_point(ref1, ref2, oth1, oth2);
        let ov = (oth2.0 - oth1.0, oth2.1 - oth1.1);
        if norm(ov) <= 1e-12 {
            return None;
        }
        for st in frag.atoms.values_mut() {
            let loc = st.loc;
            st.loc = transform2d_point(loc, trans);
            let temp = (loc.0 + st.normal.0, loc.1 + st.normal.1);
            let temp = transform2d_point(temp, trans);
            st.normal = (temp.0 - st.loc.0, temp.1 - st.loc.1);
        }
        Some(())
    }

    fn reflect_point_line(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
        rdkit_reflect_point(p, a, b)
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
            st.normal = (n_ref.0 - st.loc.0, n_ref.1 - st.loc.1);
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
        if density_normal - density_reflect > 1e-4 {
            reflect_frag(other, pin1, pin2);
        }
    }

    fn reflect_if_necessary_cis_trans(
        master: &DepictFrag,
        other: &mut DepictFrag,
        ct_case: u8,
        aid1: usize,
        aid2: usize,
    ) -> Option<()> {
        let p1_loc = master.atoms.get(&aid1)?.loc;
        let (p1_norm, ring_atom) = if ct_case == 1 {
            let aid1_state = other.atoms.get(&aid1)?;
            (aid1_state.normal, aid1_state.cis_trans_nbr?)
        } else {
            let aid1_state = master.atoms.get(&aid1)?;
            (aid1_state.normal, aid1_state.cis_trans_nbr?)
        };
        let r_atm_loc = if ct_case == 1 {
            master.atoms.get(&ring_atom)?.loc
        } else {
            other.atoms.get(&ring_atom)?.loc
        };
        let r_rel = (r_atm_loc.0 - p1_loc.0, r_atm_loc.1 - p1_loc.1);
        let dot = r_rel.0 * p1_norm.0 + r_rel.1 * p1_norm.1;
        if dot < 0.0 {
            let p2_loc = master.atoms.get(&aid2)?.loc;
            reflect_frag(other, p1_loc, p2_loc);
        }
        Some(())
    }

    fn merge_with_common_frag(
        master: &mut DepictFrag,
        other: &mut DepictFrag,
        mol: &Molecule,
        comp_set: &std::collections::BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        mut common: Vec<usize>,
    ) -> Option<()> {
        let mut ct_case = 0u8;
        if common.len() == 1 {
            let comm_aid = common[0];
            let other_atom;
            if master
                .atoms
                .get(&comm_aid)
                .and_then(|st| st.cis_trans_nbr)
                .is_some()
            {
                ct_case = 2;
                other_atom = master.atoms.get(&comm_aid)?.nbr1;
                if let Some(aid) = other_atom {
                    add_non_ring_atom_frag(other, mol, comp_set, adjacency, degree, aid, comm_aid)?;
                }
            } else if other
                .atoms
                .get(&comm_aid)
                .and_then(|st| st.cis_trans_nbr)
                .is_some()
            {
                ct_case = 1;
                other_atom = other.atoms.get(&comm_aid)?.nbr1;
                if let Some(aid) = other_atom {
                    add_non_ring_atom_frag(
                        master, mol, comp_set, adjacency, degree, aid, comm_aid,
                    )?;
                }
            } else {
                other_atom = master.atoms.get(&comm_aid)?.nbr1;
                if let Some(aid) = other_atom {
                    add_non_ring_atom_frag(other, mol, comp_set, adjacency, degree, aid, comm_aid)?;
                }
            }
            if let Some(aid) = other_atom {
                common.push(aid);
            }
        }

        if common.len() == 1 {
            unimplemented!("RDKit one-atom transform merge is not yet ported");
        }
        let aid1 = common[0];
        let aid2 = common[1];
        let ref1 = master.atoms.get(&aid1)?.loc;
        let ref2 = master.atoms.get(&aid2)?.loc;
        let oth1 = other.atoms.get(&aid1)?.loc;
        let oth2 = other.atoms.get(&aid2)?.loc;
        transform_frag_two_point(other, ref1, ref2, oth1, oth2)?;
        if common.len() >= 2 {
            if ct_case > 0 {
                reflect_if_necessary_cis_trans(master, other, ct_case, aid1, aid2)?;
            } else if common.len() == 2 {
                reflect_if_necessary_density(master, other, aid1, aid2);
            } else {
                unimplemented!("RDKit third-point reflection merge is not yet ported");
            }
        }

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
        }
        for aid in common {
            if master.atoms.contains_key(&aid) {
                update_new_neighs_for_frag(master, mol, comp_set, adjacency, degree, aid);
            }
        }
        other.done = true;
        Some(())
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
        merge_with_common_frag(
            master,
            other,
            mol,
            comp_set,
            adjacency,
            degree,
            vec![to_aid, nbr_aid],
        )
    }

    fn find_common_atoms(master: &DepictFrag, other: &DepictFrag) -> Vec<usize> {
        let mut common = Vec::new();
        for aid in master.atoms.keys() {
            if other.atoms.contains_key(aid) {
                common.push(*aid);
            }
        }
        common
    }

    fn merge_frags_with_common(
        master: &mut DepictFrag,
        frags: &mut Vec<DepictFrag>,
        mol: &Molecule,
        comp_set: &std::collections::BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
    ) -> Option<()> {
        loop {
            let mut found = None;
            for (idx, frag) in frags.iter().enumerate() {
                if frag.done {
                    continue;
                }
                let common = find_common_atoms(master, frag);
                if !common.is_empty() {
                    found = Some((idx, common));
                    break;
                }
            }
            let Some((idx, common)) = found else { break };
            let mut other = frags.remove(idx);
            let common_for_cleanup = common.clone();
            merge_with_common_frag(master, &mut other, mol, comp_set, adjacency, degree, common)?;
            for cai in common_for_cleanup {
                let remove_attach = master
                    .atoms
                    .get(&cai)
                    .is_some_and(|st| st.pending.is_empty())
                    && master.attach_pts.iter().any(|&x| x == cai);
                if remove_attach {
                    master.attach_pts.retain(|&x| x != cai);
                }
            }
        }
        Some(())
    }

    fn mirror_trans_ring_atoms(
        mol: &Molecule,
        ring: &[usize],
        coords: &mut std::collections::BTreeMap<usize, (f64, f64)>,
    ) {
        for i in 0..ring.len() {
            let atom1 = ring[i];
            let atom2 = ring[(i + 1) % ring.len()];
            let Some(bond) = mol.bonds.iter().find(|bond| {
                (bond.begin_atom == atom1 && bond.end_atom == atom2)
                    || (bond.begin_atom == atom2 && bond.end_atom == atom1)
            }) else {
                continue;
            };
            if !matches!(bond.order, BondOrder::Double) {
                continue;
            }
            if matches!(bond.stereo, BondStereo::None | BondStereo::Any)
                || bond.stereo_atoms.len() != 2
            {
                continue;
            }
            let left_is_in = ring.contains(&bond.stereo_atoms[0]);
            let right_is_in = ring.contains(&bond.stereo_atoms[1]);
            let is_trans = if matches!(bond.stereo, BondStereo::Trans) {
                left_is_in == right_is_in
            } else {
                left_is_in != right_is_in
            };
            if !is_trans {
                continue;
            }

            let left = ring[(i + ring.len() - 1) % ring.len()];
            let right = atom2;
            let last = coords[&left];
            let ref_pt = coords[&right];
            let interest = coords[&atom1];
            let d = (last.0 - ref_pt.0, last.1 - ref_pt.1);
            let dot = d.0 * d.0 + d.1 * d.1;
            let a = (d.0 * d.0 - d.1 * d.1) / dot;
            let b = 2.0 * d.0 * d.1 / dot;
            let x = a * (interest.0 - ref_pt.0) + b * (interest.1 - ref_pt.1) + ref_pt.0;
            let y = b * (interest.0 - ref_pt.0) - a * (interest.1 - ref_pt.1) + ref_pt.1;
            coords.insert(atom1, (x, y));
        }
    }

    fn init_ring_frag_from_order(mol: &Molecule, ring: &[usize]) -> DepictFrag {
        let n = ring.len();
        let radius = rdkit_ring_radius(n, 1.5);
        let largest_angle = PI * (1.0 - (2.0 / n as f64));
        let mut coords = std::collections::BTreeMap::<usize, (f64, f64)>::new();
        let ang = 2.0 * PI / n as f64;
        for (k, &a) in ring.iter().enumerate() {
            let theta = k as f64 * ang;
            coords.insert(a, (radius * rdkit_cos(theta), radius * rdkit_sin(theta)));
        }
        mirror_trans_ring_atoms(mol, ring, &mut coords);
        let mut frag = DepictFrag {
            atoms: std::collections::BTreeMap::new(),
            attach_pts: std::collections::VecDeque::new(),
            done: false,
        };
        for (k, &a) in ring.iter().enumerate() {
            let prev = if k == 0 { ring[n - 1] } else { ring[k - 1] };
            let next = ring[(k + 1) % n];
            frag.atoms.insert(
                a,
                TreeEmbeddedAtom {
                    loc: coords[&a],
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

    fn init_cis_trans_frag_from_bond(bond: &crate::Bond) -> DepictFrag {
        let begin = bond.begin_atom;
        let end = bond.end_atom;
        let begin_control = bond.stereo_atoms[0];
        let end_control = bond.stereo_atoms[1];

        let (end_normal, end_ccw) = if matches!(bond.stereo, BondStereo::Cis) {
            ((0.0, -1.0), true)
        } else {
            ((0.0, 1.0), false)
        };

        let mut frag = DepictFrag {
            atoms: std::collections::BTreeMap::new(),
            attach_pts: std::collections::VecDeque::new(),
            done: false,
        };
        frag.atoms.insert(
            begin,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (0.0, -1.0),
                ccw: false,
                cis_trans_nbr: Some(begin_control),
                angle: -1.0,
                nbr1: Some(end),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            },
        );
        frag.atoms.insert(
            end,
            TreeEmbeddedAtom {
                loc: (1.5, 0.0),
                normal: end_normal,
                ccw: end_ccw,
                cis_trans_nbr: Some(end_control),
                angle: -1.0,
                nbr1: Some(begin),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            },
        );
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
        let mut master = init_ring_frag_from_order(mol, &cycles[first]);
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
            let mut emb_ring = init_ring_frag_from_order(mol, &cycles[rid]);
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
    let mut cycles = rdkit_find_sssr_orders(mol, comp, 0);
    let ring_atoms: Vec<usize> = cycles
        .iter()
        .flatten()
        .copied()
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    if ring_atoms.is_empty() {
        return None;
    }
    let ring_set: std::collections::BTreeSet<usize> = ring_atoms.iter().copied().collect();
    let ring_bond_count = mol
        .bonds
        .iter()
        .filter(|b| ring_set.contains(&b.begin_atom) && ring_set.contains(&b.end_atom))
        .count();

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
        changed: &mut std::collections::VecDeque<usize>,
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

    fn rdkit_find_sssr_orders(
        mol: &Molecule,
        comp: &[usize],
        _target_cycle_count: usize,
    ) -> Vec<Vec<usize>> {
        fn add_ring_edges_and_atoms(
            mol: &Molecule,
            ring: &[usize],
            ring_bonds: &mut [bool],
            ring_atoms: &mut [bool],
        ) {
            for i in 0..ring.len() {
                let a = ring[i];
                let b = ring[(i + 1) % ring.len()];
                if let Some(bond) = mol.bonds.iter().find(|bond| {
                    (bond.begin_atom == a && bond.end_atom == b)
                        || (bond.begin_atom == b && bond.end_atom == a)
                }) {
                    ring_bonds[bond.index] = true;
                    ring_atoms[a] = true;
                }
            }
            if let Some(&last) = ring.last() {
                ring_atoms[last] = true;
            }
        }

        fn rdkit_find_rings_d2_nodes(
            mol: &Molecule,
            rings: &mut Vec<Vec<usize>>,
            invariants: &mut std::collections::BTreeSet<Vec<usize>>,
            d2nodes: &[usize],
            atom_degrees: &mut [usize],
            active_bonds: &mut [bool],
            ring_bonds: &mut [bool],
            ring_atoms: &mut [bool],
        ) {
            let mut dup_d2_cands = std::collections::BTreeMap::<Vec<usize>, Vec<usize>>::new();
            let mut dup_map = std::collections::BTreeMap::<usize, Vec<usize>>::new();

            for &cand in d2nodes {
                let srings = rdkit_smallest_rings_bfs(
                    mol,
                    cand,
                    active_bonds,
                    &vec![false; mol.atoms.len()],
                );
                for ring in &srings {
                    let mut inv = ring.clone();
                    inv.sort_unstable();
                    let duplicate_invars = dup_d2_cands.entry(inv.clone()).or_default();
                    if !invariants.contains(&inv) {
                        rings.push(ring.clone());
                        invariants.insert(inv);
                        add_ring_edges_and_atoms(mol, ring, ring_bonds, ring_atoms);
                    } else {
                        for &other_cand in duplicate_invars.iter() {
                            dup_map.entry(cand).or_default().push(other_cand);
                            dup_map.entry(other_cand).or_default().push(cand);
                        }
                    }
                    duplicate_invars.push(cand);
                }

                if srings.is_empty() {
                    let mut changed = std::collections::VecDeque::<usize>::new();
                    changed.push_back(cand);
                    while let Some(local_cand) = changed.pop_front() {
                        rdkit_trim_bonds(mol, local_cand, &mut changed, atom_degrees, active_bonds);
                    }
                }
            }
        }

        fn rdkit_find_rings_d3_node(
            mol: &Molecule,
            rings: &mut Vec<Vec<usize>>,
            invariants: &mut std::collections::BTreeSet<Vec<usize>>,
            cand: usize,
            active_bonds: &[bool],
        ) {
            let srings =
                rdkit_smallest_rings_bfs(mol, cand, active_bonds, &vec![false; mol.atoms.len()]);
            for ring in &srings {
                let mut inv = ring.clone();
                inv.sort_unstable();
                if invariants.insert(inv) {
                    rings.push(ring.clone());
                }
            }
            if srings.len() >= 3 {
                return;
            }

            let mut active_neighbors = Vec::<usize>::new();
            for bond in &mol.bonds {
                if !active_bonds[bond.index] {
                    continue;
                }
                if bond.begin_atom == cand {
                    active_neighbors.push(bond.end_atom);
                } else if bond.end_atom == cand {
                    active_neighbors.push(bond.begin_atom);
                }
                if active_neighbors.len() == 3 {
                    break;
                }
            }
            if active_neighbors.len() < 3 {
                return;
            }
            let n1 = active_neighbors[0];
            let n2 = active_neighbors[1];
            let n3 = active_neighbors[2];

            if srings.len() == 2 {
                let f = if srings[0].contains(&n1) && srings[1].contains(&n1) {
                    Some(n1)
                } else if srings[0].contains(&n2) && srings[1].contains(&n2) {
                    Some(n2)
                } else if srings[0].contains(&n3) && srings[1].contains(&n3) {
                    Some(n3)
                } else {
                    None
                };
                if let Some(f) = f {
                    let mut forb = vec![false; mol.atoms.len()];
                    forb[f] = true;
                    for ring in rdkit_smallest_rings_bfs(mol, cand, active_bonds, &forb) {
                        let mut inv = ring.clone();
                        inv.sort_unstable();
                        if invariants.insert(inv) {
                            rings.push(ring);
                        }
                    }
                }
            } else if srings.len() == 1 {
                let (f1, f2) = if !srings[0].contains(&n1) {
                    (n2, n3)
                } else if !srings[0].contains(&n2) {
                    (n1, n3)
                } else if !srings[0].contains(&n3) {
                    (n1, n2)
                } else {
                    return;
                };
                let mut forb = vec![false; mol.atoms.len()];
                forb[f2] = true;
                for ring in rdkit_smallest_rings_bfs(mol, cand, active_bonds, &forb) {
                    let mut inv = ring.clone();
                    inv.sort_unstable();
                    if invariants.insert(inv) {
                        rings.push(ring);
                    }
                }
                let mut forb = vec![false; mol.atoms.len()];
                forb[f1] = true;
                for ring in rdkit_smallest_rings_bfs(mol, cand, active_bonds, &forb) {
                    let mut inv = ring.clone();
                    inv.sort_unstable();
                    if invariants.insert(inv) {
                        rings.push(ring);
                    }
                }
            }
        }

        let mut atom_degrees = vec![0usize; mol.atoms.len()];
        for bond in &mol.bonds {
            atom_degrees[bond.begin_atom] += 1;
            atom_degrees[bond.end_atom] += 1;
        }
        let mut active_bonds = vec![true; mol.bonds.len()];
        let mut changed = std::collections::VecDeque::<usize>::new();
        for &idx in comp {
            if atom_degrees[idx] < 2 {
                changed.push_back(idx);
            }
        }

        let mut done_atoms = vec![false; mol.atoms.len()];
        let mut n_atoms_done = 0usize;
        let mut rings = Vec::<Vec<usize>>::new();
        let mut invariants = std::collections::BTreeSet::<Vec<usize>>::new();
        let mut ring_bonds = vec![false; mol.bonds.len()];
        let mut ring_atoms = vec![false; mol.atoms.len()];

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
                let d3node = comp.iter().copied().find(|&cand| {
                    !done_atoms[cand]
                        && atom_degrees[cand] == 3
                        && mol.bonds.iter().any(|bond| {
                            active_bonds[bond.index]
                                && (bond.begin_atom == cand || bond.end_atom == cand)
                        })
                });
                let Some(cand) = d3node else { break };
                rdkit_find_rings_d3_node(mol, &mut rings, &mut invariants, cand, &active_bonds);

                done_atoms[cand] = true;
                n_atoms_done += 1;
                rdkit_trim_bonds(
                    mol,
                    cand,
                    &mut changed,
                    &mut atom_degrees,
                    &mut active_bonds,
                );
                continue;
            }
            rdkit_find_rings_d2_nodes(
                mol,
                &mut rings,
                &mut invariants,
                &d2nodes,
                &mut atom_degrees,
                &mut active_bonds,
                &mut ring_bonds,
                &mut ring_atoms,
            );

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

    let target_cycle_count = ring_bond_count + 1 - ring_atoms.len();
    fn cycle_edges(cycle: &[usize]) -> Vec<(usize, usize)> {
        (0..cycle.len())
            .map(|i| {
                let a = cycle[i];
                let b = cycle[(i + 1) % cycle.len()];
                if a <= b { (a, b) } else { (b, a) }
            })
            .collect()
    }
    fn find_ring_connecting_atoms_rdkit_like(
        mol: &Molecule,
        start: usize,
        end: usize,
        ring_set: &std::collections::BTreeSet<usize>,
        cycles: &mut Vec<Vec<usize>>,
        ring_bond_set: &mut std::collections::BTreeSet<(usize, usize)>,
    ) -> bool {
        let mut invariants = cycles
            .iter()
            .map(|ring| {
                let mut inv = ring.clone();
                inv.sort_unstable();
                inv
            })
            .collect::<std::collections::BTreeSet<_>>();
        let mut queue = std::collections::VecDeque::<Vec<usize>>::new();
        queue.push_back(vec![start]);
        while let Some(path) = queue.pop_front() {
            let Some(&curr) = path.last() else { continue };
            for nbr in atom_neighbors(mol, curr) {
                if nbr == end {
                    if curr == start {
                        continue;
                    }
                    let mut ring = path.clone();
                    ring.push(nbr);
                    if let Some(min_pos) = ring
                        .iter()
                        .enumerate()
                        .min_by_key(|(_, atom)| **atom)
                        .map(|(idx, _)| idx)
                    {
                        ring.rotate_left(min_pos);
                        let mut reversed = ring.clone();
                        reversed.reverse();
                        if let Some(rev_min_pos) = reversed
                            .iter()
                            .enumerate()
                            .min_by_key(|(_, atom)| **atom)
                            .map(|(idx, _)| idx)
                        {
                            reversed.rotate_left(rev_min_pos);
                        }
                        if reversed.get(1) > ring.get(1) {
                            ring = reversed;
                        }
                    }
                    let mut inv = ring.clone();
                    inv.sort_unstable();
                    if invariants.insert(inv) {
                        for edge in cycle_edges(&ring) {
                            ring_bond_set.insert(edge);
                        }
                        cycles.push(ring);
                        return true;
                    }
                } else if ring_set.contains(&nbr) && !path.contains(&nbr) {
                    let mut next = path.clone();
                    next.push(nbr);
                    queue.push_back(next);
                }
            }
        }
        false
    }
    if target_cycle_count == 1 && cycles.is_empty() {
        if let Some(cycle) = rdkit_single_ring_order(mol, &ring_set, &deg_in_comp) {
            cycles.push(cycle);
        }
    }
    let mut ring_bond_set = std::collections::BTreeSet::<(usize, usize)>::new();
    for cycle in &cycles {
        for edge in cycle_edges(cycle) {
            ring_bond_set.insert(edge);
        }
    }
    if cycles.len() < target_cycle_count {
        let mut dead_bonds = std::collections::BTreeSet::<usize>::new();
        loop {
            let possible = mol.bonds.iter().find(|bond| {
                !ring_bond_set.contains(&if bond.begin_atom <= bond.end_atom {
                    (bond.begin_atom, bond.end_atom)
                } else {
                    (bond.end_atom, bond.begin_atom)
                }) && !dead_bonds.contains(&bond.index)
                    && ring_set.contains(&bond.begin_atom)
                    && ring_set.contains(&bond.end_atom)
            });
            let Some(bond) = possible else { break };
            let found = find_ring_connecting_atoms_rdkit_like(
                mol,
                bond.begin_atom,
                bond.end_atom,
                &ring_set,
                &mut cycles,
                &mut ring_bond_set,
            );
            if !found {
                dead_bonds.insert(bond.index);
            }
            if cycles.len() >= target_cycle_count {
                break;
            }
        }
        if cycles.len() < target_cycle_count {
            return None;
        }
    }
    cycles.retain(|ring| {
        ring.iter()
            .all(|a| comp_set.contains(a) && ring_set.contains(a))
    });
    if cycles.is_empty() {
        return None;
    }
    // RDKit Depictor/RDDepictor.cpp::embedFusedSystems() builds a ring
    // neighbor map and then uses RingUtils::pickFusedRings(), which is a DFS
    // over ring indices in insertion order. The resulting traversal order is
    // passed directly to EmbeddedFrag::embedFusedRings() and affects
    // pickFirstRingToEmbed() plus subsequent merge orientation.
    let mut ring_neighbors = vec![Vec::<usize>::new(); cycles.len()];
    for i in 0..cycles.len() {
        for j in i + 1..cycles.len() {
            if share_ring_atom(&cycles[i], &cycles[j]) {
                ring_neighbors[i].push(j);
                ring_neighbors[j].push(i);
            }
        }
    }
    fn pick_fused_rings_rdkit_order(
        curr: usize,
        ring_neighbors: &[Vec<usize>],
        done: &mut [bool],
        out: &mut Vec<usize>,
    ) {
        done[curr] = true;
        out.push(curr);
        for &nb in &ring_neighbors[curr] {
            if !done[nb] {
                pick_fused_rings_rdkit_order(nb, ring_neighbors, done, out);
            }
        }
    }
    let mut grouped = Vec::<Vec<usize>>::new();
    let mut done_rings = vec![false; cycles.len()];
    let mut curr = 0usize;
    while curr < cycles.len() {
        let mut group = Vec::new();
        pick_fused_rings_rdkit_order(curr, &ring_neighbors, &mut done_rings, &mut group);
        grouped.push(group);
        let mut next = None;
        for (idx, done) in done_rings.iter().enumerate() {
            if !*done {
                next = Some(idx);
                break;
            }
        }
        let Some(next_idx) = next else { break };
        curr = next_idx;
    }

    let mut frags = Vec::<DepictFrag>::new();
    for group in &grouped {
        let frag = if group.len() == 1 {
            let ring = &cycles[group[0]];
            let mut frag = init_ring_frag_from_order(mol, ring);
            setup_new_neighs_for_frag(&mut frag, mol, &comp_set, adjacency, degree);
            frag
        } else {
            build_fused_group_frag(mol, &comp_set, adjacency, degree, &cycles, group)?
        };
        frags.push(frag);
    }

    for bond in &mol.bonds {
        if !comp_set.contains(&bond.begin_atom) || !comp_set.contains(&bond.end_atom) {
            continue;
        }
        if !matches!(bond.order, BondOrder::Double)
            || !matches!(bond.stereo, BondStereo::Cis | BondStereo::Trans)
            || bond.stereo_atoms.len() != 2
        {
            continue;
        }
        let edge = if bond.begin_atom <= bond.end_atom {
            (bond.begin_atom, bond.end_atom)
        } else {
            (bond.end_atom, bond.begin_atom)
        };
        if ring_bond_set.contains(&edge) {
            continue;
        }
        let mut frag = init_cis_trans_frag_from_bond(bond);
        setup_new_neighs_for_frag(&mut frag, mol, &comp_set, adjacency, degree);
        frags.push(frag);
    }

    let embedded_atom_set: std::collections::BTreeSet<usize> = frags
        .iter()
        .flat_map(|frag| frag.atoms.keys().copied())
        .collect();

    let mut non_ring_atoms: std::collections::BTreeSet<usize> = comp
        .iter()
        .copied()
        .filter(|a| !embedded_atom_set.contains(a))
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
        merge_frags_with_common(&mut master, &mut frags, mol, &comp_set, adjacency, degree)?;
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
                    if other.atoms.contains_key(&aid) {
                        merge_with_common_frag(
                            &mut master,
                            &mut other,
                            mol,
                            &comp_set,
                            adjacency,
                            degree,
                            vec![aid],
                        )?;
                    } else {
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
        }
        if let Some(st) = master.atoms.get_mut(&aid) {
            st.pending.clear();
        }
        merge_frags_with_common(&mut master, &mut frags, mol, &comp_set, adjacency, degree)?;
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

    remove_collisions_bond_flip_like(mol, comp, adjacency, &mut local);
    remove_collisions_open_angles_like(mol, comp, adjacency, &mut local);
    remove_collisions_shorten_bonds_like(mol, comp, adjacency, &mut local);
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
            } else {
                fn recurse_deg_two_ring_atoms(
                    aid: usize,
                    adjacency: &[Vec<usize>],
                    comp_set: &std::collections::HashSet<usize>,
                    edge_is_ring_like: &dyn Fn(usize, usize) -> bool,
                    r_path: &mut Vec<usize>,
                    nbr_map: &mut std::collections::HashMap<usize, Vec<usize>>,
                ) {
                    let nbrs: Vec<usize> = adjacency[aid]
                        .iter()
                        .copied()
                        .filter(|n| comp_set.contains(n) && edge_is_ring_like(aid, *n))
                        .collect();
                    if nbrs.len() != 2 {
                        return;
                    }
                    r_path.push(aid);
                    nbr_map.insert(aid, nbrs.clone());
                    for nbr in nbrs {
                        if !r_path.contains(&nbr) {
                            recurse_deg_two_ring_atoms(
                                nbr,
                                adjacency,
                                comp_set,
                                edge_is_ring_like,
                                r_path,
                                nbr_map,
                            );
                        }
                    }
                }

                let mut r_path = Vec::<usize>::new();
                let mut nbr_map = std::collections::HashMap::<usize, Vec<usize>>::new();
                recurse_deg_two_ring_atoms(
                    aid1,
                    adjacency,
                    &comp_set,
                    &edge_is_ring_like,
                    &mut r_path,
                    &mut nbr_map,
                );
                if r_path.is_empty() {
                    recurse_deg_two_ring_atoms(
                        aid2,
                        adjacency,
                        &comp_set,
                        &edge_is_ring_like,
                        &mut r_path,
                        &mut nbr_map,
                    );
                }
                let mut move_map = std::collections::HashMap::<usize, (f64, f64)>::new();
                for &rpi in &r_path {
                    let Some(nbrs) = nbr_map.get(&rpi) else {
                        continue;
                    };
                    let p0 = pos[&rpi];
                    let p1 = pos[&nbrs[0]];
                    let p2 = pos[&nbrs[1]];
                    let mut mv = ((p1.0 + p2.0) * 0.5 - p0.0, (p1.1 + p2.1) * 0.5 - p0.1);
                    let len = (mv.0 * mv.0 + mv.1 * mv.1).sqrt();
                    if len > 1e-12 {
                        mv.0 = mv.0 / len * COLLISION_THRES;
                        mv.1 = mv.1 / len * COLLISION_THRES;
                        move_map.insert(rpi, mv);
                    }
                }
                for rpi in r_path {
                    if let Some(mv) = move_map.get(&rpi) {
                        let p = pos[&rpi];
                        pos.insert(rpi, (p.0 + mv.0, p.1 + mv.1));
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

fn bond_stereo_code(bond: &crate::Bond) -> usize {
    match (bond.direction, bond.stereo) {
        (crate::BondDirection::EndUpRight, _) => 1,
        (crate::BondDirection::EndDownRight, _) => 6,
        (_, BondStereo::Any) => 3,
        _ => 0,
    }
}

fn bond_stereo_code_with_direction(bond: &crate::Bond, direction: BondDirection) -> usize {
    match (direction, bond.stereo) {
        (BondDirection::EndUpRight, _) => 1,
        (BondDirection::EndDownRight, _) => 6,
        (_, BondStereo::Any) => 3,
        _ => 0,
    }
}

pub(crate) fn pick_bonds_to_wedge_rdkit_subset(mol: &Molecule) -> HashMap<usize, usize> {
    const NO_NBRS: i32 = 100;
    let atom_rings = rdkit_atom_rings_for_wedging(mol);
    let bond_rings = rdkit_bond_rings_for_wedging(mol, &atom_rings);
    let mut n_chiral_nbrs = vec![NO_NBRS; mol.atoms.len()];

    for atom in &mol.atoms {
        if !matches!(
            atom.chiral_tag,
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            continue;
        }
        n_chiral_nbrs[atom.index] = 0;
        for nb in atom_neighbors(mol, atom.index) {
            if mol.atoms[nb].atomic_num == 1 {
                n_chiral_nbrs[atom.index] -= 10;
                continue;
            }
            if matches!(
                mol.atoms[nb].chiral_tag,
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) {
                n_chiral_nbrs[atom.index] -= 1;
            }
        }
    }

    let mut indices: Vec<usize> = (0..mol.atoms.len()).collect();
    indices.sort_by_key(|i| n_chiral_nbrs[*i]);
    let mut wedge = HashMap::new();

    for idx in indices {
        if n_chiral_nbrs[idx] > NO_NBRS {
            continue;
        }
        let atom = &mol.atoms[idx];
        if !matches!(
            atom.chiral_tag,
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            break;
        }
        if let Some(bid) = pick_bond_to_wedge_rdkit_subset(
            idx,
            mol,
            &n_chiral_nbrs,
            &wedge,
            NO_NBRS,
            &atom_rings,
            &bond_rings,
        ) {
            wedge.insert(bid, idx);
        }
    }
    wedge
}

fn rdkit_atom_rings_for_wedging(mol: &Molecule) -> Vec<Vec<usize>> {
    crate::distgeom::rdkit_atom_rings(mol)
}

fn rdkit_bond_rings_for_wedging(mol: &Molecule, atom_rings: &[Vec<usize>]) -> Vec<Vec<usize>> {
    atom_rings
        .iter()
        .map(|ring| {
            let mut out = Vec::new();
            for i in 0..ring.len() {
                let a1 = ring[i];
                let a2 = ring[(i + 1) % ring.len()];
                if let Some(bond) = mol.bonds.iter().find(|bond| {
                    (bond.begin_atom == a1 && bond.end_atom == a2)
                        || (bond.begin_atom == a2 && bond.end_atom == a1)
                }) {
                    out.push(bond.index);
                }
            }
            out
        })
        .collect()
}

fn pick_bond_to_wedge_rdkit_subset(
    atom_idx: usize,
    mol: &Molecule,
    n_chiral_nbrs: &[i32],
    already: &HashMap<usize, usize>,
    no_nbrs: i32,
    atom_rings: &[Vec<usize>],
    bond_rings: &[Vec<usize>],
) -> Option<usize> {
    let mut best: Option<(i32, usize)> = None;
    for bond in &mol.bonds {
        if !matches!(bond.order, BondOrder::Single) {
            continue;
        }
        let other = if bond.begin_atom == atom_idx {
            Some(bond.end_atom)
        } else if bond.end_atom == atom_idx {
            Some(bond.begin_atom)
        } else {
            None
        };
        let Some(other_idx) = other else { continue };
        if already.contains_key(&bond.index) {
            continue;
        }
        let other_atom = &mol.atoms[other_idx];
        let score = if other_atom.atomic_num == 1 {
            -1_000_000
        } else {
            let mut s = other_atom.atomic_num as i32
                + 100 * atom_degree(mol, other_idx) as i32
                + 1000 * i32::from(!matches!(other_atom.chiral_tag, ChiralTag::Unspecified));
            if n_chiral_nbrs[other_idx] < no_nbrs {
                s -= 100_000 * n_chiral_nbrs[other_idx];
            }
            s += 10_000
                * atom_rings
                    .iter()
                    .filter(|ring| ring.contains(&other_idx))
                    .count() as i32;
            s += 20_000
                * bond_rings
                    .iter()
                    .filter(|ring| ring.contains(&bond.index))
                    .count() as i32;
            let (has_double, has_known_double, has_any_double) =
                get_double_bond_presence(mol, other_idx, usize::MAX);
            s += 11_000 * i32::from(has_double);
            s += 12_000 * i32::from(has_known_double);
            s += 23_000 * i32::from(has_any_double);
            s
        };
        match best {
            None => best = Some((score, bond.index)),
            Some((bscore, _)) if score < bscore => best = Some((score, bond.index)),
            _ => {}
        }
    }
    best.map(|(_, bid)| bid)
}

pub(crate) fn determine_bond_wedge_state_rdkit_subset(
    bond: &crate::Bond,
    wedge: &HashMap<usize, usize>,
    coords: &[glam::DVec2],
    mol: &Molecule,
) -> BondDirection {
    let Some(&from_atom) = wedge.get(&bond.index) else {
        return bond.direction;
    };
    if !matches!(bond.order, BondOrder::Single) {
        return bond.direction;
    }
    let center = from_atom;
    let other = if bond.begin_atom == center {
        bond.end_atom
    } else {
        bond.begin_atom
    };
    let chiral = mol.atoms[center].chiral_tag;
    if !matches!(chiral, ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw) {
        return bond.direction;
    }

    let center_loc = coords[center];
    let ref_vec = coords[other] - center_loc;
    let ref_len2 = ref_vec.length_squared();
    if ref_len2 <= 1e-18 {
        return bond.direction;
    }

    let mut neighbor_bonds = vec![(bond.index, 0.0f64)];
    for nb in atom_bond_indices(mol, center) {
        if nb == bond.index {
            continue;
        }
        let nbond = &mol.bonds[nb];
        let other_atom = if nbond.begin_atom == center {
            nbond.end_atom
        } else {
            nbond.begin_atom
        };
        let tmp_vec = coords[other_atom] - center_loc;
        if tmp_vec.length_squared() <= 1e-18 {
            return bond.direction;
        }
        let mut angle = signed_angle_2d(ref_vec, tmp_vec);
        if angle < 0.0 {
            angle += 2.0 * PI;
        }
        neighbor_bonds.push((nb, angle));
    }
    neighbor_bonds.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let canonical = atom_bond_indices(mol, center);
    let mut perm: Vec<usize> = neighbor_bonds.iter().map(|(bid, _)| *bid).collect();
    let n_swaps = permutation_parity_swaps(&mut perm, &canonical);

    let mut swaps = n_swaps;
    if neighbor_bonds.len() == 3 {
        let angle1 = neighbor_bonds[1].1;
        let angle2 = neighbor_bonds[2].1;
        let angle_tol = PI * 1.9 / 180.0;
        if angle2 - angle1 >= (PI - angle_tol) {
            swaps += 1;
        }
    }

    match chiral {
        ChiralTag::TetrahedralCcw => {
            if swaps % 2 == 1 {
                BondDirection::EndDownRight
            } else {
                BondDirection::EndUpRight
            }
        }
        ChiralTag::TetrahedralCw => {
            if swaps % 2 == 1 {
                BondDirection::EndUpRight
            } else {
                BondDirection::EndDownRight
            }
        }
        ChiralTag::Unspecified => BondDirection::None,
    }
}

fn atom_neighbors(mol: &Molecule, atom_idx: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for b in &mol.bonds {
        if b.begin_atom == atom_idx {
            out.push(b.end_atom);
        } else if b.end_atom == atom_idx {
            out.push(b.begin_atom);
        }
    }
    out
}

fn atom_bond_indices(mol: &Molecule, atom_idx: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for b in &mol.bonds {
        if b.begin_atom == atom_idx || b.end_atom == atom_idx {
            out.push(b.index);
        }
    }
    out
}

fn atom_degree(mol: &Molecule, atom_idx: usize) -> usize {
    atom_bond_indices(mol, atom_idx).len()
}

fn get_double_bond_presence(
    mol: &Molecule,
    atom_idx: usize,
    ignore_bond: usize,
) -> (bool, bool, bool) {
    let mut has_double = false;
    let mut has_known_double = false;
    let mut has_any_double = false;
    for b in &mol.bonds {
        if b.index == ignore_bond {
            continue;
        }
        if b.begin_atom != atom_idx && b.end_atom != atom_idx {
            continue;
        }
        if matches!(b.order, BondOrder::Double) {
            has_double = true;
            if matches!(b.stereo, BondStereo::Any) {
                has_any_double = true;
            } else if matches!(b.stereo, BondStereo::Cis | BondStereo::Trans) {
                has_known_double = true;
            }
        }
    }
    (has_double, has_known_double, has_any_double)
}

fn signed_angle_2d(a: glam::DVec2, b: glam::DVec2) -> f64 {
    let cross = a.x * b.y - a.y * b.x;
    let dot = a.x * b.x + a.y * b.y;
    cross.atan2(dot)
}

fn permutation_parity_swaps(perm: &mut [usize], canonical: &[usize]) -> i32 {
    if perm.len() != canonical.len() {
        return 0;
    }
    let mut pos = HashMap::new();
    for (i, v) in canonical.iter().enumerate() {
        pos.insert(*v, i);
    }
    let mut target = Vec::with_capacity(perm.len());
    for v in perm.iter() {
        let Some(&p) = pos.get(v) else { return 0 };
        target.push(p);
    }
    let mut seen = vec![false; target.len()];
    let mut swaps = 0i32;
    for i in 0..target.len() {
        if seen[i] {
            continue;
        }
        let mut j = i;
        let mut cycle = 0usize;
        while !seen[j] {
            seen[j] = true;
            j = target[j];
            cycle += 1;
        }
        if cycle > 0 {
            swaps += (cycle as i32) - 1;
        }
    }
    swaps
}

fn bond_v3000_cfg_code(bond: &crate::Bond) -> Option<usize> {
    match bond.direction {
        crate::BondDirection::EndUpRight => Some(1),
        crate::BondDirection::EndDownRight => Some(3),
        crate::BondDirection::None => {
            if matches!(bond.stereo, BondStereo::Any) {
                Some(2)
            } else {
                None
            }
        }
    }
}

fn append_v2000_property_lines(out: &mut String, mol: &Molecule) {
    let charges: Vec<(usize, i8)> = mol
        .atoms
        .iter()
        .filter_map(|atom| {
            if atom.formal_charge == 0 {
                None
            } else {
                Some((atom.index + 1, atom.formal_charge))
            }
        })
        .collect();
    for chunk in charges.chunks(8) {
        out.push_str(&format!("M  CHG{:>3}", chunk.len()));
        for (idx, charge) in chunk {
            out.push_str(&format!("{idx:>4}{charge:>4}"));
        }
        out.push('\n');
    }

    let isotopes: Vec<(usize, u16)> = mol
        .atoms
        .iter()
        .filter_map(|atom| atom.isotope.map(|isotope| (atom.index + 1, isotope)))
        .collect();
    for chunk in isotopes.chunks(8) {
        out.push_str(&format!("M  ISO{:>3}", chunk.len()));
        for (idx, isotope) in chunk {
            out.push_str(&format!("{idx:>4}{isotope:>4}"));
        }
        out.push('\n');
    }

    let radicals: Vec<(usize, u8)> = mol
        .atoms
        .iter()
        .filter_map(|atom| {
            if atom.num_radical_electrons == 0 {
                None
            } else {
                let code = if atom.num_radical_electrons % 2 == 1 {
                    2
                } else {
                    3
                };
                Some((atom.index + 1, code))
            }
        })
        .collect();
    for chunk in radicals.chunks(8) {
        out.push_str(&format!("M  RAD{:>3}", chunk.len()));
        for (idx, radical) in chunk {
            out.push_str(&format!("{idx:>4}{radical:>4}"));
        }
        out.push('\n');
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::{Path, PathBuf};
    use std::process::Command;

    use super::mol_to_v2000_block;
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

    #[derive(Debug, thiserror::Error)]
    enum TestDataError {
        #[error(transparent)]
        Io(#[from] std::io::Error),
        #[error("invalid jsonl at line {line_no}: {source}")]
        Json {
            line_no: usize,
            #[source]
            source: serde_json::Error,
        },
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
        let atom_count = counts.get(0..3)?.trim().parse::<usize>().ok()?;
        let bond_count = counts.get(3..6)?.trim().parse::<usize>().ok()?;
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

    fn canonical_bonds_for_compare(bonds: &[(usize, usize, usize)]) -> Vec<(usize, usize, usize)> {
        let mut out = bonds.to_vec();
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

    fn non_coordinate_sections_match(ours_body: &str, expected_body: &str) -> bool {
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
            if !atom_symbols_equivalent(ours_symbol, expected_symbol) {
                return false;
            }
        }

        canonical_bonds_for_compare(&ours.bonds) == canonical_bonds_for_compare(&expected.bonds)
            && coords_ok
    }

    fn compare_against_expected(
        ours_body: &str,
        expected_body: &str,
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
            non_coordinate_sections_match(&ours_norm, &expected_norm),
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
            let mut mol = parsed.expect("parse checked above");
            mol.compute_2d_coords().unwrap_or_else(|e| {
                panic!("2D coordinate generation failed at row {}: {}", idx + 1, e)
            });
            let ours = mol_to_v2000_block(&mol)
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
                crate::kekulize::kekulize_in_place(&mut mol, true)
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

            mol.compute_2d_coords().unwrap_or_else(|e| {
                panic!("2D coordinate generation failed at row {}: {}", idx + 1, e)
            });
            let ours = mol_to_v2000_block(&mol)
                .unwrap_or_else(|e| panic!("write failed at row {}: {}", idx + 1, e));
            let ours_body = body(&ours);

            if record.v2000_ok {
                let expected_v2000 = record
                    .v2000_body
                    .as_ref()
                    .expect("v2000_ok=true requires v2000_body");
                let ours_norm = normalize_signed_zero(&ours_body);
                let expected_norm = normalize_signed_zero(expected_v2000);
                assert!(
                    bonds_and_atoms_match_strict(&ours_norm, &expected_norm),
                    "kekulized bond block mismatch at row {} ({}) against v2000\nours:\n{}\nexpected:\n{}",
                    idx + 1,
                    record.smiles,
                    ours_norm,
                    expected_norm
                );
            }
            if record.v3000_ok {
                let expected_v3000 = record
                    .v3000_body
                    .as_ref()
                    .expect("v3000_ok=true requires v3000_body");
                let ours_norm = normalize_signed_zero(&ours_body);
                let expected_norm = normalize_signed_zero(expected_v3000);
                assert!(
                    bonds_and_atoms_match_strict(&ours_norm, &expected_norm),
                    "kekulized bond block mismatch at row {} ({}) against v3000\nours:\n{}\nexpected:\n{}",
                    idx + 1,
                    record.smiles,
                    ours_norm,
                    expected_norm
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

        let mut mol = Molecule::from_smiles(&record.smiles).expect("parse diagnostic molecule");
        mol.compute_2d_coords()
            .expect("compute diagnostic 2D coordinates");
        let ours_block = mol_to_v2000_block(&mol).expect("write diagnostic molecule");
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
        let mut mol = Molecule::from_smiles(&record.smiles).expect("parse diagnostic molecule");
        mol.compute_2d_coords()
            .expect("compute diagnostic 2D coordinates");
        let ours_block = mol_to_v2000_block(&mol).expect("write diagnostic molecule");
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
}
