use core::fmt;
use std::f64::consts::PI;

use cosmolkit_chem_core::{BondOrder, Molecule};

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
    let coords = compute_2d_coords_minimal(mol)?;

    let mut out = String::new();
    out.push('\n');
    out.push_str("     COSMolKit      2D\n");
    out.push('\n');
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        mol.atoms.len(),
        mol.bonds.len()
    ));

    for (idx, atom) in mol.atoms.iter().enumerate() {
        let (x, y) = coords[idx];
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            x,
            y,
            0.0,
            atom_symbol(atom.atomic_num)
        ));
    }

    for bond in &mol.bonds {
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
    Ok(general_coords_fallback(mol))
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
    Some(vec![(-1.2990, -0.2500), (0.0, 0.5000), (1.2990, -0.2500)])
}

fn general_coords_fallback(mol: &Molecule) -> Vec<(f64, f64)> {
    let n = mol.atoms.len();
    let mut adjacency = vec![Vec::<usize>::new(); n];
    let mut degree = vec![0usize; n];
    for b in &mol.bonds {
        adjacency[b.begin_atom].push(b.end_atom);
        adjacency[b.end_atom].push(b.begin_atom);
        degree[b.begin_atom] += 1;
        degree[b.end_atom] += 1;
    }

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

    let mut coords = vec![(0.0, 0.0); n];
    let mut local_components: Vec<Vec<(usize, (f64, f64))>> = Vec::new();
    for comp in components {
        let k = comp.len();
        if k == 1 {
            local_components.push(vec![(comp[0], (0.0, 0.0))]);
            continue;
        }
        if k == 2 {
            let b = mol
                .bonds
                .iter()
                .find(|b| {
                    (b.begin_atom == comp[0] && b.end_atom == comp[1])
                        || (b.begin_atom == comp[1] && b.end_atom == comp[0])
                })
                .map(|b| b.order);
            let pair = if matches!(b, Some(BondOrder::Null)) {
                vec![(comp[0], (0.7500, 0.0)), (comp[1], (-0.7500, 0.0))]
            } else {
                vec![(comp[0], (-0.7500, 0.0)), (comp[1], (0.7500, 0.0))]
            };
            local_components.push(pair);
            continue;
        }

        // Dative linear triad template (e.g. [NH3]->[Cu+2]<-[NH3]).
        if k == 3 {
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
                let mut center = None;
                for &a in &comp {
                    if degree[a] == 2 {
                        center = Some(a);
                        break;
                    }
                }
                if let Some(c) = center {
                    let mut ends: Vec<usize> = comp.iter().copied().filter(|&x| x != c).collect();
                    ends.sort_unstable();
                    let mut local = Vec::new();
                    local.push((ends[0], (1.5000, 0.0)));
                    local.push((c, (0.0, 0.0)));
                    local.push((ends[1], (-1.5000, 0.0)));
                    local.sort_by_key(|(idx, _)| *idx);
                    local_components.push(local);
                    continue;
                }
            }
        }

        // Pure simple cycle: regular polygon.
        let bond_count_in_comp = mol
            .bonds
            .iter()
            .filter(|b| comp.contains(&b.begin_atom) && comp.contains(&b.end_atom))
            .count();

        // One-ring + tree substituents.
        let cyclomatic = bond_count_in_comp as isize - k as isize + 1;
        if cyclomatic == 1 {
            if let Some(local) = place_unicyclic_component(mol, &comp, &adjacency) {
                local_components.push(local);
                continue;
            }
        } else if cyclomatic > 1 {
            if let Some(local) = place_fused_hex_component(mol, &comp, &adjacency) {
                local_components.push(local);
                continue;
            }
            if let Some(mut local) = place_acyclic_tree_rdkit_like(mol, &comp, &adjacency, &degree)
            {
                regularize_isolated_aromatic_six_cycles(mol, &comp, &adjacency, &mut local);
                remove_collisions_bond_flip_like(mol, &comp, &adjacency, &mut local);
                remove_collisions_open_angles_like(mol, &comp, &adjacency, &mut local);
                optimize_aromatic_ring_reflections(mol, &comp, &adjacency, &mut local);
                local_components.push(local);
                continue;
            }
        }

        // Tree star templates (common in early special-case corpus).
        let is_tree = bond_count_in_comp + 1 == k;
        if is_tree {
            let centers: Vec<_> = comp.iter().copied().filter(|&a| degree[a] > 1).collect();
            if centers.len() == 1 {
                let center = centers[0];
                let mut leaves: Vec<_> = comp.iter().copied().filter(|&a| a != center).collect();
                if degree[center] == 3
                    && leaves.len() == 3
                    && leaves.iter().all(|&a| degree[a] == 1)
                {
                    leaves.sort_by_key(|&a| {
                        let atom = &mol.atoms[a];
                        let anum_rank = if atom.atomic_num == 1 {
                            1000
                        } else {
                            atom.atomic_num as usize
                        };
                        (100 * anum_rank + degree[a], a)
                    });
                    let local = vec![
                        (center, (0.0, 0.0)),
                        (leaves[0], (-1.2990, -0.7500)),
                        (leaves[1], (1.2990, -0.7500)),
                        (leaves[2], (0.0000, 1.5000)),
                    ];
                    let mut local = local;
                    local.sort_by_key(|(idx, _)| *idx);
                    local_components.push(local);
                    continue;
                }
                if degree[center] == 4
                    && leaves.len() == 4
                    && leaves.iter().all(|&a| degree[a] == 1)
                {
                    leaves.sort_unstable();
                    let local = vec![
                        (center, (0.0, 0.0)),
                        (leaves[0], (-1.2990, -0.7500)),
                        (leaves[1], (1.2990, 0.7500)),
                        (leaves[2], (0.7500, -1.2990)),
                        (leaves[3], (-0.7500, 1.2990)),
                    ];
                    let mut local = local;
                    local.sort_by_key(|(idx, _)| *idx);
                    local_components.push(local);
                    continue;
                }
            }
            if let Some(local) = place_acyclic_tree_rdkit_like(mol, &comp, &adjacency, &degree) {
                local_components.push(local);
                continue;
            }
        }

        let is_simple_cycle = bond_count_in_comp == k && comp.iter().all(|&a| degree[a] == 2);
        if is_simple_cycle {
            let bond_len = 1.5f64;
            let radius = bond_len / (2.0 * (PI / k as f64).sin());
            let aromatic_component = mol.bonds.iter().any(|b| {
                comp.contains(&b.begin_atom)
                    && comp.contains(&b.end_atom)
                    && matches!(b.order, BondOrder::Aromatic)
            });
            let direction = if aromatic_component { 1.0 } else { -1.0 };
            let mut local = Vec::new();
            for (rank, &atom_idx) in comp.iter().enumerate() {
                let theta = direction * 2.0 * PI * rank as f64 / k as f64;
                local.push((atom_idx, (radius * theta.cos(), radius * theta.sin())));
            }
            local.sort_by_key(|(idx, _)| *idx);
            local_components.push(local);
            continue;
        }

        let bond_len = 1.5f64;
        let radius = bond_len / (2.0 * (PI / k as f64).sin());
        let mut local = Vec::new();
        for (rank, &atom_idx) in comp.iter().enumerate() {
            let theta = 2.0 * PI * rank as f64 / k as f64;
            local.push((atom_idx, (radius * theta.cos(), radius * theta.sin())));
        }
        local.sort_by_key(|(idx, _)| *idx);
        local_components.push(local);
    }

    for component in &mut local_components {
        canonicalize_component_rdkit_like(component);
    }

    // Mimic RDKit DepictorLocal::_shiftCoords().
    if local_components.is_empty() {
        return coords;
    }
    let box_of = |component: &Vec<(usize, (f64, f64))>| -> (f64, f64, f64, f64) {
        let mut max_x = f64::NEG_INFINITY;
        let mut min_x = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        for &(_, (x, y)) in component {
            if x > max_x {
                max_x = x;
            }
            if x < min_x {
                min_x = x;
            }
            if y > max_y {
                max_y = y;
            }
            if y < min_y {
                min_y = y;
            }
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
        coords[idx] = (x, y);
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
            coords[idx] = (x + shift_x, y + shift_y);
        }
    }

    coords
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
    // RDKit runs this on higher-precision coordinates; with our minimal writer
    // many templates are rounded to 4 decimals, so near-degenerate covariance
    // can spuriously rotate symmetric motifs. Treat tiny-d as degenerate.
    if d <= 1e-3 {
        return;
    }
    let mut eigx = 2.0 * xy;
    let mut eigy = (yy - xx) + d;
    let eig_len = (eigx * eigx + eigy * eigy).sqrt();
    if eig_len <= 1e-4 {
        return;
    }
    eigx /= eig_len;
    eigy /= eig_len;

    for (i, (_, pos)) in component.iter_mut().enumerate() {
        let (px, py) = centered[i];
        // trans = [[eig1.x, eig1.y],[-eig1.y,eig1.x]]
        let rx = px * eigx + py * eigy;
        let ry = -px * eigy + py * eigx;
        *pos = (rx, ry);
    }
}

#[derive(Clone)]
struct TreeEmbeddedAtom {
    loc: (f64, f64),
    normal: (f64, f64),
    ccw: bool,
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

fn compute_sub_angle_from_degree(degree: usize) -> f64 {
    if degree == 4 {
        PI / 2.0
    } else {
        2.0 * PI / 3.0
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

    let root = *comp.iter().min_by_key(|&&idx| {
        let depict = atom_depict_rank(mol.atoms[idx].atomic_num, degree[idx]);
        depict * mol.atoms.len() + idx
    })?;

    states[root] = Some(TreeEmbeddedAtom {
        loc: (0.0, 0.0),
        normal: (1.0, 0.0),
        ccw: true,
        angle: -1.0,
        nbr1: None,
        nbr2: None,
        rot_dir: 0,
        pending: Vec::new(),
    });
    unembedded.remove(&root);

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
            neighs.sort_by_key(|&idx| {
                let depict = atom_depict_rank(mol.atoms[idx].atomic_num, degree[idx]);
                depict * mol.atoms.len() + idx
            });
            if let Some(st) = states[aid].as_mut() {
                st.pending = neighs;
                if !st.pending.is_empty() && !attach_in_queue.contains(&aid) {
                    attach_pts.push_back(aid);
                    attach_in_queue.insert(aid);
                }
            }
        };

    update_new_neighs(
        root,
        &mut states,
        &unembedded,
        &mut attach_pts,
        &mut attach_in_queue,
    );

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
                if (rem_angle.abs() - PI).abs() < 1e-3 {
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
                    angle: -1.0,
                    nbr1: Some(aid),
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                });
            } else {
                let ref_ccw = ref_state.ccw;
                let mut curr_vec = ref_state.normal;
                let degree_here = degree[aid];
                let mut angle = compute_sub_angle_from_degree(degree_here);
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

    // Keep orientation deterministic; aromatic 5-rings tend to align better with +theta.
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

    // In case of any leftover atom (unexpected), signal fallback.
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

        let mut best_pair: Option<(usize, usize)> = None;
        for i in 0..6 {
            let j = (i + 3) % 6;
            if ext_deg(cyc[i]) > 0 && ext_deg(cyc[j]) > 0 {
                let key = (cyc[i].min(cyc[j]), cyc[i].max(cyc[j]), i.min(j), i.max(j));
                let cur = best_pair.map(|(bi, bj)| {
                    let a = cyc[bi];
                    let b = cyc[bj];
                    (a.min(b), a.max(b), bi.min(bj), bi.max(bj))
                });
                if cur.is_none() || Some(key) < cur {
                    best_pair = Some((i, j));
                }
            }
        }
        let Some((i0, i3)) = best_pair else { continue };

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
        let score_for = |candidate: &std::collections::HashMap<usize, (f64, f64)>| -> f64 {
            let mut score = 0.0;
            for &a in &order {
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
            for &a in &order {
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
            for &a in &order {
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
            let s = score_for(&cand);
            if s < best_score {
                best_score = s;
                best_coords = Some(cand);
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

    let mut candidates: Vec<(usize, usize)> = Vec::new();
    for b in &mol.bonds {
        if !matches!(b.order, BondOrder::Single) {
            continue;
        }
        if !comp_set.contains(&b.begin_atom) || !comp_set.contains(&b.end_atom) {
            continue;
        }
        let u = b.begin_atom;
        let v = b.end_atom;
        if ring_edge(u, v) {
            continue;
        }
        candidates.push((u, v));
    }

    let mut best_score = score(&pos);
    for _ in 0..3 {
        let mut improved = false;
        for &(u, v) in &candidates {
            let mut side = std::collections::HashSet::<usize>::new();
            let mut stack = vec![v];
            side.insert(v);
            while let Some(cur) = stack.pop() {
                for &nb in &adjacency[cur] {
                    if !comp_set.contains(&nb) {
                        continue;
                    }
                    if cur == v && nb == u {
                        continue;
                    }
                    if side.insert(nb) {
                        stack.push(nb);
                    }
                }
            }
            if side.is_empty() || side.len() == comp.len() {
                continue;
            }
            let pu = pos[&u];
            let pv = pos[&v];
            let mut trial = pos.clone();
            for &a in &side {
                if a == u {
                    continue;
                }
                let p = trial[&a];
                trial.insert(a, reflect_point(p, pu, pv));
            }
            let sc = score(&trial);
            if sc + 1e-8 < best_score {
                pos = trial;
                best_score = sc;
                improved = true;
            }
        }
        if !improved {
            break;
        }
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
        let mut best_pair: Option<(usize, usize)> = None;
        for i in 0..6 {
            let j = (i + 3) % 6;
            if ext_deg(cyc[i]) > 0 && ext_deg(cyc[j]) > 0 {
                best_pair = Some((i, j));
                break;
            }
        }
        let Some((i0, i3)) = best_pair else { continue };
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

    // Enumerate unique 6-cycles from ring core.
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
    let mut ring_atoms_sorted = ring_atoms.clone();
    ring_atoms_sorted.sort_unstable();
    for &start in &ring_atoms_sorted {
        let mut path = vec![start];
        let mut used = std::collections::HashSet::<usize>::new();
        used.insert(start);
        fn dfs6(
            start: usize,
            cur: usize,
            depth: usize,
            path: &mut Vec<usize>,
            used: &mut std::collections::HashSet<usize>,
            adjacency: &[Vec<usize>],
            ring_set: &std::collections::HashSet<usize>,
            seen: &mut std::collections::HashSet<Vec<usize>>,
            cycles: &mut Vec<Vec<usize>>,
        ) {
            if depth == 6 {
                if adjacency[cur].contains(&start) {
                    let c = canonical_cycle(path.clone());
                    if seen.insert(c.clone()) {
                        cycles.push(c);
                    }
                }
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
                dfs6(
                    start,
                    nb,
                    depth + 1,
                    path,
                    used,
                    adjacency,
                    ring_set,
                    seen,
                    cycles,
                );
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
            adjacency,
            &ring_set,
            &mut seen,
            &mut cycles,
        );
    }
    if cycles.is_empty() {
        return None;
    }

    // Keep only cycles whose all edges exist in ring-core.
    cycles.retain(|cyc| {
        (0..cyc.len()).all(|i| {
            let a = cyc[i];
            let b = cyc[(i + 1) % cyc.len()];
            adjacency[a].contains(&b) && ring_set.contains(&a) && ring_set.contains(&b)
        })
    });
    if cycles.is_empty() {
        return None;
    }

    // This routine is only for fused-ring systems. If no pair of cycles shares
    // an edge, route to the non-fused pipeline instead.
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
    if !has_shared_edge {
        return None;
    }

    let mut ring_centers = vec![(0.0f64, 0.0f64); cycles.len()];
    let mut ring_done = vec![false; cycles.len()];
    let mut coords: std::collections::HashMap<usize, (f64, f64)> = std::collections::HashMap::new();

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

            let try_build = |reverse: bool,
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
    use cosmolkit_chem_core::Molecule;
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
                    cosmolkit_chem_core::BondOrder::Aromatic
                        | cosmolkit_chem_core::BondOrder::Dative
                        | cosmolkit_chem_core::BondOrder::Null
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
        assert!(
            non_coordinate_sections_match(
                &normalize_signed_zero(ours_body),
                &normalize_signed_zero(expected_body),
                mol
            ),
            "molblock mismatch (including coordinates) at row {} ({}) against {}",
            row_idx_1based,
            smiles,
            variant
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
        let path = repo_root().join("tests/molblock_minimal_smiles.txt");
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
        let input = repo.join("tests/molblock_minimal_smiles.txt");

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
             uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_v2000_minimal_golden.py --input tests/molblock_minimal_smiles.txt --output tests/golden/molblock_v2000_minimal.jsonl",
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
        let input = repo.join("tests/molblock_minimal_smiles.txt");

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
             uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_kekulize_molblock_golden.py --input tests/molblock_minimal_smiles.txt --output tests/golden/molblock_v2000_kekulized.jsonl",
            golden_path.display(),
            last_error
        );
    }

    #[test]
    fn molblock_minimal_golden_has_entry_for_each_smiles() {
        let smiles = load_smiles().expect("read tests/molblock_minimal_smiles.txt");
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
    fn molblock_minimal_matches_rdkit_golden() {
        let golden = load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
        for (idx, record) in golden.iter().enumerate() {
            let mol = Molecule::from_smiles(&record.smiles)
                .unwrap_or_else(|e| panic!("parse failed at row {}: {}", idx + 1, e));
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
    fn molblock_kekulize_golden_has_entry_for_each_smiles() {
        let smiles = load_smiles().expect("read tests/molblock_minimal_smiles.txt");
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
    fn molblock_kekulized_bonds_match_rdkit_golden() {
        let golden =
            load_kekulize_golden().expect("read tests/golden/molblock_v2000_kekulized.jsonl");
        for (idx, record) in golden.iter().enumerate() {
            let mut mol = Molecule::from_smiles(&record.smiles)
                .unwrap_or_else(|e| panic!("parse failed at row {}: {}", idx + 1, e));

            let ours_kek = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                cosmolkit_chem_core::kekulize::kekulize_in_place(&mut mol)
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
}
