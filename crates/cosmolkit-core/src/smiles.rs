use std::collections::{HashMap, HashSet, VecDeque};

use crate::{
    Atom, Bond, BondDirection, BondOrder, BondStereo, ChiralTag, Molecule, SmilesParseError,
    ValenceModel, assign_valence, valence::valence_list,
};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum PendingBond {
    Unspecified,
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    DirectionalSingleUp,
    DirectionalSingleDown,
    DativeForward,
    DativeBackward,
    Null,
}

#[derive(Debug, Clone)]
struct RingOpen {
    atom_index: usize,
    bond: PendingBond,
    opening_slot: usize,
}

pub(crate) fn parse_smiles(smiles: &str) -> Result<Molecule, SmilesParseError> {
    let mut parser = Parser::new(smiles);
    let mut mol = parser.parse_molecule()?;
    parser.skip_ascii_whitespace();
    if !parser.is_eof() {
        return Err(parser.error("unexpected trailing characters"));
    }
    parser.reorder_ring_closure_bonds_rdkit_like(&mut mol);
    adjust_atom_chirality_flags(
        &mut mol,
        &parser.smiles_start_atoms,
        &parser.atom_ring_closure_slots,
        &parser.ring_closure_bond_numbers,
    );
    validate_aromatic_atoms_are_in_rings(&mol)?;
    perceive_aromaticity_rdkit_like(&mut mol, &parser.ring_closure_bonds)?;
    prune_noncyclic_aromatic_bonds_rdkit_like(&mut mol);
    assign_double_bond_stereo(&mut mol);
    canonicalize_double_bond_stereo_by_cip_rdkit_like(&mut mol);
    mol.rebuild_adjacency();
    Ok(mol)
}

/// Assign double-bond stereochemistry from stored directional single bonds.
///
/// This mirrors the same RDKit-aligned path used after SMILES parsing, without
/// requiring the molecule to have originated from a SMILES string.
pub fn assign_double_bond_stereo_from_directions(mol: &mut Molecule) {
    assign_double_bond_stereo(mol);
    canonicalize_double_bond_stereo_by_cip_rdkit_like(mol);
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

fn is_unsaturated_atom(mol: &Molecule, atom_index: usize) -> bool {
    mol.bonds.iter().any(|b| {
        (b.begin_atom == atom_index || b.end_atom == atom_index)
            && bond_order_as_double(b.order) > 1.0
    })
}

fn atom_has_fourth_valence(
    atom: &Atom,
    atom_index: usize,
    assignment: Option<&crate::ValenceAssignment>,
) -> bool {
    if atom.explicit_hydrogens == 1 {
        return true;
    }
    if let Some(a) = assignment {
        return a.implicit_hydrogens[atom_index] == 1;
    }
    false
}

fn count_swaps_to_interconvert(probe: &[usize], reference: &[usize]) -> Option<usize> {
    if probe.len() != reference.len() {
        return None;
    }
    let mut work = probe.to_vec();
    let mut swaps = 0usize;
    for i in 0..work.len() {
        if work[i] == reference[i] {
            continue;
        }
        let mut found = None;
        for (j, value) in work.iter().enumerate().skip(i + 1) {
            if *value == reference[i] {
                found = Some(j);
                break;
            }
        }
        let j = found?;
        work.swap(i, j);
        swaps += 1;
    }
    Some(swaps)
}

fn smiles_bond_ordering_for_atom(
    mol: &Molecule,
    atom_index: usize,
    ring_closures: &[usize],
) -> Vec<usize> {
    let mut neighbors = Vec::<(usize, usize)>::new();
    neighbors.push((atom_index, usize::MAX));
    for b in &mol.bonds {
        if b.begin_atom != atom_index && b.end_atom != atom_index {
            continue;
        }
        if ring_closures.contains(&b.index) {
            continue;
        }
        let nbr = if b.begin_atom == atom_index {
            b.end_atom
        } else {
            b.begin_atom
        };
        neighbors.push((nbr, b.index));
    }
    neighbors.sort_by_key(|&(nbr, _)| nbr);

    let mut self_pos = 0usize;
    if neighbors.first().is_some_and(|(idx, _)| *idx != atom_index) {
        self_pos = 1;
    }

    let mut bond_ordering = Vec::<usize>::new();
    for (i, &(_, bond_idx)) in neighbors.iter().enumerate() {
        if i == self_pos {
            bond_ordering.extend_from_slice(ring_closures);
        } else {
            bond_ordering.push(bond_idx);
        }
    }
    bond_ordering
}

fn atom_bond_storage_order_rdkit_like(
    mol: &Molecule,
    atom_index: usize,
    ring_closures: &[usize],
    ring_closure_bond_numbers: &[Option<u32>],
) -> Vec<usize> {
    let mut out = Vec::new();
    for b in &mol.bonds {
        if b.begin_atom == atom_index || b.end_atom == atom_index {
            if ring_closures.contains(&b.index) {
                continue;
            }
            out.push(b.index);
        }
    }
    let mut sorted_ring_closures = ring_closures.to_vec();
    sorted_ring_closures.sort_by_key(|&bond_idx| {
        ring_closure_bond_numbers
            .get(bond_idx)
            .and_then(|n| *n)
            .unwrap_or(u32::MAX)
    });
    for &bond_idx in &sorted_ring_closures {
        if let Some(b) = mol.bonds.get(bond_idx)
            && (b.begin_atom == atom_index || b.end_atom == atom_index)
        {
            out.push(bond_idx);
        }
    }
    out
}

fn invert_atom_chiral_tag(atom: &mut Atom) {
    atom.chiral_tag = match atom.chiral_tag {
        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
        ChiralTag::Unspecified => ChiralTag::Unspecified,
    };
}

fn adjust_atom_chirality_flags(
    mol: &mut Molecule,
    smiles_start_atoms: &[bool],
    atom_ring_closure_slots: &[Vec<Option<usize>>],
    ring_closure_bond_numbers: &[Option<u32>],
) {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).ok();
    for atom_idx in 0..mol.atoms.len() {
        if !matches!(
            mol.atoms[atom_idx].chiral_tag,
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            continue;
        }
        let ring_closures: Vec<usize> = atom_ring_closure_slots
            .get(atom_idx)
            .into_iter()
            .flat_map(|v| v.iter())
            .filter_map(|slot| *slot)
            .collect();
        let bond_ordering = smiles_bond_ordering_for_atom(mol, atom_idx, &ring_closures);
        let storage_order = atom_bond_storage_order_rdkit_like(
            mol,
            atom_idx,
            &ring_closures,
            ring_closure_bond_numbers,
        );
        let mut n_swaps = count_swaps_to_interconvert(&bond_ordering, &storage_order).unwrap_or(0);

        let atom = &mol.atoms[atom_idx];
        let needs_inversion = storage_order.len() == 3
            && ((smiles_start_atoms.get(atom_idx).copied().unwrap_or(false)
                && atom.explicit_hydrogens == 1)
                || (!atom_has_fourth_valence(atom, atom_idx, assignment.as_ref())
                    && ring_closures.len() == 1
                    && !is_unsaturated_atom(mol, atom_idx)));
        if needs_inversion {
            n_swaps += 1;
        }
        if n_swaps % 2 == 1 {
            invert_atom_chiral_tag(&mut mol.atoms[atom_idx]);
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum ElectronDonorType {
    Vacant,
    One,
    Two,
    Any,
    No,
}

fn bond_valence_contrib_for_atom(bond: &Bond, atom_index: usize) -> i32 {
    if bond.begin_atom != atom_index && bond.end_atom != atom_index {
        return 0;
    }
    match bond.order {
        BondOrder::Null => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 1,
        BondOrder::Dative => {
            if bond.end_atom == atom_index {
                1
            } else {
                0
            }
        }
    }
}

fn n_outer_electrons(atomic_num: u8) -> Option<i32> {
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

fn more_electronegative(anum1: u8, anum2: u8) -> bool {
    // Mirrors RDKit PeriodicTable::moreElectroNegative().
    let ne1 = n_outer_electrons(anum1).unwrap_or(0);
    let ne2 = n_outer_electrons(anum2).unwrap_or(0);
    ne1 > ne2 || (ne1 == ne2 && anum1 < anum2)
}

fn default_valence_rdkit(atomic_num: u8) -> Option<i32> {
    let vals = valence_list(atomic_num)?;
    for &v in vals {
        if v >= 0 {
            return Some(v);
        }
    }
    None
}

fn count_atom_electrons_rdkit(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_index: usize,
    atom_degree: &[usize],
) -> i32 {
    let atom = &mol.atoms[atom_index];
    let Some(default_valence) = default_valence_rdkit(atom.atomic_num) else {
        return -1;
    };
    if default_valence <= 1 {
        return -1;
    }

    let mut degree = atom_degree[atom_index] as i32
        + atom.explicit_hydrogens as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    for bond in &mol.bonds {
        if (bond.begin_atom == atom_index || bond.end_atom == atom_index)
            && bond_valence_contrib_for_atom(bond, atom_index) == 0
        {
            degree -= 1;
        }
    }
    if degree > 3 {
        return -1;
    }

    let Some(nouter) = n_outer_electrons(atom.atomic_num) else {
        return -1;
    };
    let nlp = (nouter - default_valence - atom.formal_charge as i32).max(0);
    let radicals = atom.num_radical_electrons as i32;

    let mut res = (default_valence - degree) + nlp - radicals;
    if res > 1 {
        let explicit_valence = assignment.explicit_valence[atom_index] as i32;
        let n_unsaturations = explicit_valence - atom_degree[atom_index] as i32;
        if n_unsaturations > 1 {
            res = 1;
        }
    }
    res
}

fn is_multiple_bond(order: BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
    )
}

fn incident_non_cyclic_multiple_bond(
    mol: &Molecule,
    atom_index: usize,
    ring_bonds: &HashSet<usize>,
) -> Option<usize> {
    for bond in &mol.bonds {
        if bond.begin_atom != atom_index && bond.end_atom != atom_index {
            continue;
        }
        if ring_bonds.contains(&bond.index) || !is_multiple_bond(bond.order) {
            continue;
        }
        return Some(if bond.begin_atom == atom_index {
            bond.end_atom
        } else {
            bond.begin_atom
        });
    }
    None
}

fn incident_cyclic_multiple_bond(
    mol: &Molecule,
    atom_index: usize,
    ring_bonds: &HashSet<usize>,
) -> bool {
    for bond in &mol.bonds {
        if bond.begin_atom != atom_index && bond.end_atom != atom_index {
            continue;
        }
        if ring_bonds.contains(&bond.index) && is_multiple_bond(bond.order) {
            return true;
        }
    }
    false
}

fn incident_multiple_bond(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_index: usize,
    atom_degree: &[usize],
) -> bool {
    let atom = &mol.atoms[atom_index];
    let mut degree = atom_degree[atom_index] as i32 + atom.explicit_hydrogens as i32;
    for bond in &mol.bonds {
        if (bond.begin_atom == atom_index || bond.end_atom == atom_index)
            && bond_valence_contrib_for_atom(bond, atom_index) == 0
        {
            degree -= 1;
        }
    }
    assignment.explicit_valence[atom_index] as i32 != degree
}

fn get_atom_donor_type_rdkit(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_index: usize,
    atom_degree: &[usize],
    ring_bonds: &HashSet<usize>,
) -> ElectronDonorType {
    let atom = &mol.atoms[atom_index];
    if atom.atomic_num == 0 {
        return if incident_cyclic_multiple_bond(mol, atom_index, ring_bonds) {
            ElectronDonorType::One
        } else {
            ElectronDonorType::Any
        };
    }
    let nelec = count_atom_electrons_rdkit(mol, assignment, atom_index, atom_degree);
    if nelec < 0 {
        return ElectronDonorType::No;
    }
    if nelec == 0 {
        return if incident_non_cyclic_multiple_bond(mol, atom_index, ring_bonds).is_some() {
            ElectronDonorType::Vacant
        } else if incident_cyclic_multiple_bond(mol, atom_index, ring_bonds) {
            ElectronDonorType::One
        } else {
            ElectronDonorType::No
        };
    }
    if nelec == 1 {
        if let Some(other) = incident_non_cyclic_multiple_bond(mol, atom_index, ring_bonds) {
            if more_electronegative(mol.atoms[other].atomic_num, atom.atomic_num) {
                return ElectronDonorType::Vacant;
            }
            return ElectronDonorType::One;
        }
        if incident_multiple_bond(mol, assignment, atom_index, atom_degree) {
            return ElectronDonorType::One;
        }
        return if atom.formal_charge == 1 {
            ElectronDonorType::Vacant
        } else {
            ElectronDonorType::No
        };
    }
    let mut adjusted = nelec;
    if let Some(other) = incident_non_cyclic_multiple_bond(mol, atom_index, ring_bonds)
        && more_electronegative(mol.atoms[other].atomic_num, atom.atomic_num)
    {
        adjusted -= 1;
    }
    if adjusted % 2 == 1 {
        ElectronDonorType::One
    } else {
        ElectronDonorType::Two
    }
}

fn ring_electron_bounds(
    ring: &[usize],
    edon: &[ElectronDonorType],
) -> Option<(usize, usize, usize)> {
    let mut low = 0usize;
    let mut high = 0usize;
    let mut any_count = 0usize;
    for &idx in ring {
        match edon[idx] {
            ElectronDonorType::Any => {
                any_count += 1;
                low += 1;
                high += 2;
            }
            ElectronDonorType::One => {
                low += 1;
                high += 1;
            }
            ElectronDonorType::Two => {
                low += 2;
                high += 2;
            }
            ElectronDonorType::Vacant | ElectronDonorType::No => {}
        }
    }
    if any_count > 1 {
        return None;
    }
    Some((low, high, any_count))
}

fn apply_huckel_rule(ring: &[usize], edon: &[ElectronDonorType], min_ring_size: usize) -> bool {
    if ring.len() < min_ring_size {
        return false;
    }
    let Some((low, high, _)) = ring_electron_bounds(ring, edon) else {
        return false;
    };
    if high == 2 {
        return true;
    }
    if high < 6 {
        return false;
    }
    let mut e = low;
    while e <= high {
        if (e + 2) % 4 == 0 {
            return true;
        }
        e += 1;
    }
    false
}

fn is_atom_candidate_for_aromaticity_rdkit(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_index: usize,
    donor_type: ElectronDonorType,
) -> bool {
    let atom = &mol.atoms[atom_index];

    // Mirrors RDKit isAtomCandForArom() default model constraints.
    if atom.atomic_num > 18 && atom.atomic_num != 34 && atom.atomic_num != 52 {
        return false;
    }
    if !matches!(
        donor_type,
        ElectronDonorType::Vacant
            | ElectronDonorType::One
            | ElectronDonorType::Two
            | ElectronDonorType::Any
    ) {
        return false;
    }

    // Atoms above default valence are excluded from aromatic candidates.
    if let Some(def_val) = default_valence_rdkit(atom.atomic_num)
        && def_val > 0
    {
        let adjusted_atomic_num = atom.atomic_num as i32 - atom.formal_charge as i32;
        if adjusted_atomic_num > 0
            && let Some(adjusted_def_val) = default_valence_rdkit(adjusted_atomic_num as u8)
        {
            let total_valence = assignment.explicit_valence[atom_index] as i32
                + assignment.implicit_hydrogens[atom_index] as i32;
            if total_valence > adjusted_def_val {
                return false;
            }
        }
    }

    if atom.num_radical_electrons > 0 && (atom.atomic_num != 6 || atom.formal_charge != 0) {
        return false;
    }

    // Disallow atoms with more than one multiple (double/triple) bond when
    // unsaturation exceeds one, mirroring RDKit sf.net bug 1934360 handling.
    let atom_degree = mol
        .bonds
        .iter()
        .filter(|b| b.begin_atom == atom_index || b.end_atom == atom_index)
        .count() as i32;
    let n_unsaturations = assignment.explicit_valence[atom_index] as i32 - atom_degree;
    if n_unsaturations > 1 {
        let mut n_mult = 0usize;
        for bond in &mol.bonds {
            if bond.begin_atom != atom_index && bond.end_atom != atom_index {
                continue;
            }
            if matches!(bond.order, BondOrder::Double | BondOrder::Triple) {
                n_mult += 1;
                if n_mult > 1 {
                    return false;
                }
            }
        }
    }

    true
}

fn canonical_cycle_key(cycle: &[usize]) -> Vec<usize> {
    let n = cycle.len();
    let mut best = cycle.to_vec();
    for start in 0..n {
        let mut rotated = Vec::with_capacity(n);
        for i in 0..n {
            rotated.push(cycle[(start + i) % n]);
        }
        if rotated < best {
            best = rotated.clone();
        }
        rotated.reverse();
        if rotated < best {
            best = rotated;
        }
    }
    best
}

fn shortest_path_ignoring_edge(
    adj: &[Vec<(usize, usize)>],
    src: usize,
    dst: usize,
    forbidden_edge: usize,
) -> Option<Vec<usize>> {
    let n = adj.len();
    let mut prev = vec![usize::MAX; n];
    let mut seen = vec![false; n];
    let mut q = VecDeque::new();
    q.push_back(src);
    seen[src] = true;

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

fn find_candidate_rings(mol: &Molecule, max_ring_size: usize) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::<(usize, usize)>::new(); mol.atoms.len()];
    for (bi, b) in mol.bonds.iter().enumerate() {
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }

    let mut seen = HashSet::<Vec<usize>>::new();
    let mut rings = Vec::<Vec<usize>>::new();
    for (bi, b) in mol.bonds.iter().enumerate() {
        let Some(path) = shortest_path_ignoring_edge(&adj, b.begin_atom, b.end_atom, bi) else {
            continue;
        };
        if path.len() < 3 || path.len() > max_ring_size {
            continue;
        }
        let key = canonical_cycle_key(&path);
        if seen.insert(key) {
            rings.push(path);
        }
    }
    rings
}

fn graph_component_count(mol: &Molecule) -> usize {
    let mut seen = vec![false; mol.atoms.len()];
    let mut comps = 0usize;
    for start in 0..mol.atoms.len() {
        if seen[start] {
            continue;
        }
        comps += 1;
        let mut q = VecDeque::new();
        q.push_back(start);
        seen[start] = true;
        while let Some(u) = q.pop_front() {
            for b in &mol.bonds {
                let v = if b.begin_atom == u {
                    b.end_atom
                } else if b.end_atom == u {
                    b.begin_atom
                } else {
                    continue;
                };
                if !seen[v] {
                    seen[v] = true;
                    q.push_back(v);
                }
            }
        }
    }
    comps
}

fn highest_set_bit(words: &[u64]) -> Option<usize> {
    for (wi, &w) in words.iter().enumerate().rev() {
        if w != 0 {
            let bit = 63usize - w.leading_zeros() as usize;
            return Some(wi * 64 + bit);
        }
    }
    None
}

fn reduce_to_min_cycle_basis_indices(mol: &Molecule, rings: &[Vec<usize>]) -> Vec<usize> {
    if rings.is_empty() || mol.bonds.is_empty() {
        return Vec::new();
    }
    let cyclomatic = mol.bonds.len() + graph_component_count(mol) - mol.atoms.len();
    if cyclomatic == 0 {
        return Vec::new();
    }
    let words = mol.bonds.len().div_ceil(64);

    let mut entries: Vec<(usize, Vec<usize>, Vec<u64>)> = rings
        .iter()
        .enumerate()
        .map(|(idx, ring)| {
            let mut bits = vec![0u64; words];
            let mut edges = cycle_bond_indices(mol, ring);
            edges.sort_unstable();
            for &e in &edges {
                bits[e / 64] ^= 1u64 << (e % 64);
            }
            (idx, ring.clone(), bits)
        })
        .collect();
    entries.sort_by(|a, b| {
        let ka = canonical_cycle_key(&a.1);
        let kb = canonical_cycle_key(&b.1);
        (a.1.len(), ka).cmp(&(b.1.len(), kb))
    });

    let mut basis: Vec<(usize, Vec<u64>)> = Vec::new();
    let mut selected = Vec::<usize>::new();
    for (ring_idx, _ring, mut row) in entries {
        for (pivot, brow) in &basis {
            if (row[pivot / 64] >> (pivot % 64)) & 1 == 1 {
                for wi in 0..row.len() {
                    row[wi] ^= brow[wi];
                }
            }
        }
        let Some(pivot) = highest_set_bit(&row) else {
            continue;
        };
        basis.push((pivot, row));
        basis.sort_by(|l, r| r.0.cmp(&l.0));
        selected.push(ring_idx);
        if selected.len() >= cyclomatic {
            break;
        }
    }
    selected
}

fn cycle_bond_indices(mol: &Molecule, cycle: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(cycle.len());
    for i in 0..cycle.len() {
        let a = cycle[i];
        let b = cycle[(i + 1) % cycle.len()];
        if let Some(bond) = mol.bonds.iter().find(|bond| {
            (bond.begin_atom == a && bond.end_atom == b)
                || (bond.begin_atom == b && bond.end_atom == a)
        }) {
            out.push(bond.index);
        }
    }
    out
}

fn build_ring_neighbor_map(ring_bonds: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut neighs = vec![Vec::<usize>::new(); ring_bonds.len()];
    for i in 0..ring_bonds.len() {
        let set_i: HashSet<usize> = ring_bonds[i].iter().copied().collect();
        for j in (i + 1)..ring_bonds.len() {
            if ring_bonds[j].iter().any(|b| set_i.contains(b)) {
                neighs[i].push(j);
                neighs[j].push(i);
            }
        }
    }
    neighs
}

fn connected_ring_components(neighs: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; neighs.len()];
    let mut out = Vec::<Vec<usize>>::new();
    for start in 0..neighs.len() {
        if seen[start] {
            continue;
        }
        let mut q = VecDeque::new();
        q.push_back(start);
        seen[start] = true;
        let mut comp = Vec::<usize>::new();
        while let Some(u) = q.pop_front() {
            comp.push(u);
            for &v in &neighs[u] {
                if !seen[v] {
                    seen[v] = true;
                    q.push_back(v);
                }
            }
        }
        out.push(comp);
    }
    out
}

fn is_connected_ring_subset(subset: &[usize], neighs: &[Vec<usize>]) -> bool {
    if subset.len() <= 1 {
        return true;
    }
    let allowed: HashSet<usize> = subset.iter().copied().collect();
    let mut seen = HashSet::<usize>::new();
    let mut q = VecDeque::new();
    q.push_back(subset[0]);
    seen.insert(subset[0]);
    while let Some(u) = q.pop_front() {
        for &v in &neighs[u] {
            if !allowed.contains(&v) || seen.contains(&v) {
                continue;
            }
            seen.insert(v);
            q.push_back(v);
        }
    }
    seen.len() == subset.len()
}

fn for_each_combination<F>(items: &[usize], choose: usize, f: &mut F)
where
    F: FnMut(&[usize]),
{
    fn rec<F>(items: &[usize], choose: usize, start: usize, acc: &mut Vec<usize>, f: &mut F)
    where
        F: FnMut(&[usize]),
    {
        if acc.len() == choose {
            f(acc);
            return;
        }
        let remaining = choose - acc.len();
        if remaining == 0 {
            f(acc);
            return;
        }
        if start >= items.len() || items.len() - start < remaining {
            return;
        }
        for i in start..=items.len() - remaining {
            acc.push(items[i]);
            rec(items, choose, i + 1, acc, f);
            acc.pop();
        }
    }

    if choose == 0 || choose > items.len() {
        return;
    }
    let mut acc = Vec::<usize>::with_capacity(choose);
    rec(items, choose, 0, &mut acc, f);
}

fn mark_aromatic_subset_rdkit_like(
    mol: &mut Molecule,
    rings: &[Vec<usize>],
    ring_bonds: &[Vec<usize>],
    subset: &[usize],
    done_bonds: &mut HashSet<usize>,
) {
    for &ring_idx in subset {
        for &atom_idx in &rings[ring_idx] {
            mol.atoms[atom_idx].is_aromatic = true;
        }
    }

    let mut bond_counts = std::collections::BTreeMap::<usize, usize>::new();
    for &ring_idx in subset {
        for &bond_idx in &ring_bonds[ring_idx] {
            *bond_counts.entry(bond_idx).or_insert(0) += 1;
        }
    }
    for (bond_idx, count) in bond_counts {
        if count != 1 {
            continue;
        }
        if let Some(bond) = mol.bonds.get_mut(bond_idx) {
            if matches!(bond.order, BondOrder::Single | BondOrder::Double) {
                bond.order = BondOrder::Aromatic;
            }
            bond.is_aromatic = true;
            done_bonds.insert(bond_idx);
        }
    }
}

fn perceive_aromaticity_rdkit_like(
    mol: &mut Molecule,
    ring_closure_bonds: &[usize],
) -> Result<(), SmilesParseError> {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike)
        .map_err(|e| SmilesParseError::ParseError(format!("valence assignment failed: {:?}", e)))?;
    let mut atom_degree = vec![0usize; mol.atoms.len()];
    for b in &mol.bonds {
        atom_degree[b.begin_atom] += 1;
        atom_degree[b.end_atom] += 1;
    }
    // Mirrors RDKit maxFusedAromaticRingSize constant used for candidate rings.
    // Candidate rings are derived from explicit SMILES ring closures, matching
    // the graph construction path used by the parser.
    let _ = ring_closure_bonds;
    let all_rings = find_candidate_rings(mol, 24);
    let ring_basis_ids = reduce_to_min_cycle_basis_indices(mol, &all_rings);
    let rings: Vec<Vec<usize>> = if ring_basis_ids.is_empty() {
        all_rings
    } else {
        ring_basis_ids
            .into_iter()
            .map(|idx| all_rings[idx].clone())
            .collect()
    };

    let mut ring_bonds = HashSet::<usize>::new();
    for ring in &rings {
        for bi in cycle_bond_indices(mol, ring) {
            ring_bonds.insert(bi);
        }
    }
    let donor_types: Vec<ElectronDonorType> = (0..mol.atoms.len())
        .map(|idx| get_atom_donor_type_rdkit(mol, &assignment, idx, &atom_degree, &ring_bonds))
        .collect();
    let aromatic_candidates: Vec<bool> = (0..mol.atoms.len())
        .map(|idx| is_atom_candidate_for_aromaticity_rdkit(mol, &assignment, idx, donor_types[idx]))
        .collect();

    // Candidate rings cRings in RDKit aromaticityHelper().
    let mut candidate_ring_ids = Vec::<usize>::new();
    for (ring_idx, ring) in rings.iter().enumerate() {
        let all_dummy = ring.iter().all(|&idx| mol.atoms[idx].atomic_num == 0);
        if all_dummy {
            continue;
        }
        if !ring.iter().all(|&idx| aromatic_candidates[idx]) {
            continue;
        }
        candidate_ring_ids.push(ring_idx);
    }

    if candidate_ring_ids.is_empty() {
        return Ok(());
    }

    let candidate_rings: Vec<Vec<usize>> = candidate_ring_ids
        .iter()
        .map(|&idx| rings[idx].clone())
        .collect();
    let candidate_ring_bonds: Vec<Vec<usize>> = candidate_rings
        .iter()
        .map(|ring| cycle_bond_indices(mol, ring))
        .collect();
    let neighs = build_ring_neighbor_map(&candidate_ring_bonds);
    let components = connected_ring_components(&neighs);

    // Mirrors applyHuckelToFused():
    // check connected ring combinations from size 1..=min(nrings, 6),
    // then mark atoms + outer bonds for aromatic combinations.
    for component in components {
        let max_choose = component.len().min(6);
        let mut done_bonds = HashSet::<usize>::new();
        let mut n_ring_bonds = HashSet::<usize>::new();
        for &ring_idx in &component {
            for &b in &candidate_ring_bonds[ring_idx] {
                n_ring_bonds.insert(b);
            }
        }

        for choose in 1..=max_choose {
            let mut combo = Vec::<usize>::new();
            let mut walker = |picked: &[usize]| {
                combo.clear();
                combo.extend_from_slice(picked);
                if !is_connected_ring_subset(&combo, &neighs) {
                    return;
                }

                let mut ats_in_system = vec![0usize; mol.atoms.len()];
                for &ring_idx in &combo {
                    for &a in &candidate_rings[ring_idx] {
                        ats_in_system[a] += 1;
                    }
                }
                let unon: Vec<usize> = ats_in_system
                    .iter()
                    .enumerate()
                    .filter_map(|(a, &count)| {
                        if count == 1 || count == 2 {
                            Some(a)
                        } else {
                            None
                        }
                    })
                    .collect();
                if apply_huckel_rule(&unon, &donor_types, 0) {
                    mark_aromatic_subset_rdkit_like(
                        mol,
                        &candidate_rings,
                        &candidate_ring_bonds,
                        &combo,
                        &mut done_bonds,
                    );
                }
            };
            for_each_combination(&component, choose, &mut walker);
            if done_bonds.len() >= n_ring_bonds.len() {
                break;
            }
        }
    }
    Ok(())
}

fn opposite_bond_direction(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        BondDirection::None => BondDirection::None,
    }
}

fn has_stereo_bond_dir(direction: BondDirection) -> bool {
    matches!(
        direction,
        BondDirection::EndUpRight | BondDirection::EndDownRight
    )
}

fn assign_double_bond_stereo(mol: &mut Molecule) {
    // Mirrors RDKit Chirality.cpp::setBondStereoFromDirections().
    let bond_count = mol.bonds.len();
    for db_idx in 0..bond_count {
        if !matches!(mol.bonds[db_idx].order, BondOrder::Double) {
            continue;
        }
        let begin = mol.bonds[db_idx].begin_atom;
        let end = mol.bonds[db_idx].end_atom;

        let begin_control = mol.bonds.iter().enumerate().find_map(|(idx, bond)| {
            if idx == db_idx || matches!(bond.order, BondOrder::Double) {
                return None;
            }
            if !has_stereo_bond_dir(bond.direction) {
                return None;
            }
            if bond.begin_atom == begin {
                let mut dir = bond.direction;
                // RDKit begin-side rule: flip if the directed bond begins at
                // the begin atom of the stereo double bond.
                if bond.begin_atom == begin {
                    dir = opposite_bond_direction(dir);
                }
                Some((bond.end_atom, dir))
            } else if bond.end_atom == begin {
                Some((bond.begin_atom, bond.direction))
            } else {
                None
            }
        });

        let end_control = mol.bonds.iter().enumerate().find_map(|(idx, bond)| {
            if idx == db_idx || matches!(bond.order, BondOrder::Double) {
                return None;
            }
            if !has_stereo_bond_dir(bond.direction) {
                return None;
            }
            if bond.begin_atom == end {
                Some((bond.end_atom, bond.direction))
            } else if bond.end_atom == end {
                let mut dir = bond.direction;
                // RDKit end-side rule: flip if the directed bond ends at the
                // end atom of the stereo double bond.
                if bond.end_atom == end {
                    dir = opposite_bond_direction(dir);
                }
                Some((bond.begin_atom, dir))
            } else {
                None
            }
        });

        let (Some((begin_atom, begin_dir)), Some((end_atom, end_dir))) =
            (begin_control, end_control)
        else {
            continue;
        };

        mol.bonds[db_idx].stereo_atoms = vec![begin_atom, end_atom];
        // RDKit: equal directional markers means trans, different means cis.
        mol.bonds[db_idx].stereo = if begin_dir == end_dir {
            BondStereo::Trans
        } else {
            BondStereo::Cis
        };
    }
}

fn highest_cip_neighbor_excluding(
    mol: &Molecule,
    atom_index: usize,
    skip_atom_index: usize,
    cip_ranks: &[i64],
) -> Option<usize> {
    let mut best_atom = None::<usize>;
    let mut best_rank = i64::MIN;
    for bond in &mol.bonds {
        let nbr = if bond.begin_atom == atom_index {
            bond.end_atom
        } else if bond.end_atom == atom_index {
            bond.begin_atom
        } else {
            continue;
        };
        if nbr == skip_atom_index {
            continue;
        }
        let rank = *cip_ranks.get(nbr)?;
        if best_atom.is_none() || rank > best_rank {
            best_atom = Some(nbr);
            best_rank = rank;
        } else if rank == best_rank {
            // Mirrors RDKit findHighestCIPNeighbor(): ties invalidate result.
            best_atom = None;
        }
    }
    best_atom
}

fn canonicalize_double_bond_stereo_by_cip_rdkit_like(mol: &mut Molecule) {
    // Mirrors the post-direction double-bond stereo normalization in RDKit
    // assignStereochemistry: controlling atoms are chosen by highest CIP rank,
    // and cis/trans is inverted when controlling atoms flip parity.
    let cip_ranks = crate::io::molblock::rdkit_cip_ranks_for_depict(mol);
    for bond_idx in 0..mol.bonds.len() {
        if !matches!(mol.bonds[bond_idx].order, BondOrder::Double) {
            continue;
        }
        if !matches!(
            mol.bonds[bond_idx].stereo,
            BondStereo::Cis | BondStereo::Trans
        ) {
            continue;
        }
        let begin = mol.bonds[bond_idx].begin_atom;
        let end = mol.bonds[bond_idx].end_atom;
        let Some(begin_ctrl) = highest_cip_neighbor_excluding(mol, begin, end, &cip_ranks) else {
            continue;
        };
        let Some(end_ctrl) = highest_cip_neighbor_excluding(mol, end, begin, &cip_ranks) else {
            continue;
        };

        let mut flips = 0usize;
        if mol.bonds[bond_idx].stereo_atoms.len() == 2 {
            if mol.bonds[bond_idx].stereo_atoms[0] != begin_ctrl {
                flips += 1;
            }
            if mol.bonds[bond_idx].stereo_atoms[1] != end_ctrl {
                flips += 1;
            }
        }
        mol.bonds[bond_idx].stereo_atoms = vec![begin_ctrl, end_ctrl];
        if flips % 2 == 1 {
            mol.bonds[bond_idx].stereo = match mol.bonds[bond_idx].stereo {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                other => other,
            };
        }
    }
}

fn prune_noncyclic_aromatic_bonds_rdkit_like(mol: &mut Molecule) {
    // Aromatic bonds must participate in an aromatic cycle. Any aromatic edge
    // that is a bridge in the aromatic-only subgraph is demoted.
    let n_atoms = mol.atoms.len();
    if n_atoms == 0 {
        return;
    }
    let mut adj = vec![Vec::<(usize, usize)>::new(); n_atoms];
    for (bi, b) in mol.bonds.iter().enumerate() {
        if !b.is_aromatic {
            continue;
        }
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }
    let mut disc = vec![usize::MAX; n_atoms];
    let mut low = vec![usize::MAX; n_atoms];
    let mut time = 0usize;
    let mut is_bridge = vec![false; mol.bonds.len()];

    fn dfs(
        u: usize,
        parent_edge: Option<usize>,
        adj: &[Vec<(usize, usize)>],
        disc: &mut [usize],
        low: &mut [usize],
        time: &mut usize,
        is_bridge: &mut [bool],
    ) {
        disc[u] = *time;
        low[u] = *time;
        *time += 1;
        for &(v, eidx) in &adj[u] {
            if Some(eidx) == parent_edge {
                continue;
            }
            if disc[v] == usize::MAX {
                dfs(v, Some(eidx), adj, disc, low, time, is_bridge);
                low[u] = low[u].min(low[v]);
                if low[v] > disc[u] {
                    is_bridge[eidx] = true;
                }
            } else {
                low[u] = low[u].min(disc[v]);
            }
        }
    }

    for u in 0..n_atoms {
        if disc[u] == usize::MAX {
            dfs(
                u,
                None,
                &adj,
                &mut disc,
                &mut low,
                &mut time,
                &mut is_bridge,
            );
        }
    }
    for (bi, is_br) in is_bridge.into_iter().enumerate() {
        if !is_br {
            continue;
        }
        if let Some(bond) = mol.bonds.get_mut(bi)
            && bond.is_aromatic
        {
            bond.is_aromatic = false;
            bond.order = BondOrder::Single;
        }
    }
}

fn validate_aromatic_atoms_are_in_rings(mol: &Molecule) -> Result<(), SmilesParseError> {
    // Mirror RDKit sanitize behavior for obvious aromatic-invalid inputs:
    // aromatic atoms must belong to at least one cycle (ring).
    let n = mol.atoms.len();
    if n == 0 || mol.bonds.is_empty() {
        return Ok(());
    }

    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }

    // Tarjan bridges: an edge is in a ring iff it is not a bridge.
    let mut disc = vec![usize::MAX; n];
    let mut low = vec![usize::MAX; n];
    let mut time = 0usize;
    let mut is_bridge = vec![false; mol.bonds.len()];

    fn dfs(
        u: usize,
        parent_edge: Option<usize>,
        adj: &[Vec<(usize, usize)>],
        disc: &mut [usize],
        low: &mut [usize],
        time: &mut usize,
        is_bridge: &mut [bool],
    ) {
        disc[u] = *time;
        low[u] = *time;
        *time += 1;
        for &(v, eidx) in &adj[u] {
            if Some(eidx) == parent_edge {
                continue;
            }
            if disc[v] == usize::MAX {
                dfs(v, Some(eidx), adj, disc, low, time, is_bridge);
                low[u] = low[u].min(low[v]);
                if low[v] > disc[u] {
                    is_bridge[eidx] = true;
                }
            } else {
                low[u] = low[u].min(disc[v]);
            }
        }
    }

    for u in 0..n {
        if disc[u] == usize::MAX {
            dfs(
                u,
                None,
                &adj,
                &mut disc,
                &mut low,
                &mut time,
                &mut is_bridge,
            );
        }
    }

    for a in 0..n {
        if !mol.atoms[a].is_aromatic {
            continue;
        }
        let in_ring = adj[a].iter().any(|&(_, eidx)| !is_bridge[eidx]);
        if !in_ring {
            return Err(SmilesParseError::ParseError(
                "non-ring atom marked aromatic".to_string(),
            ));
        }
    }
    Ok(())
}

struct Parser<'a> {
    input: &'a str,
    pos: usize,
    ring_opens: HashMap<u32, RingOpen>,
    ring_closure_bonds: Vec<usize>,
    ring_closure_bond_numbers: Vec<Option<u32>>,
    smiles_start_atoms: Vec<bool>,
    atom_ring_closure_slots: Vec<Vec<Option<usize>>>,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Self {
        Self {
            input,
            pos: 0,
            ring_opens: HashMap::new(),
            ring_closure_bonds: Vec::new(),
            ring_closure_bond_numbers: Vec::new(),
            smiles_start_atoms: Vec::new(),
            atom_ring_closure_slots: Vec::new(),
        }
    }

    fn add_atom_with_tracking(
        &mut self,
        mol: &mut Molecule,
        atom: Atom,
        smiles_start: bool,
    ) -> usize {
        let idx = mol.add_atom(atom);
        if self.smiles_start_atoms.len() <= idx {
            self.smiles_start_atoms.resize(idx + 1, false);
        }
        if self.atom_ring_closure_slots.len() <= idx {
            self.atom_ring_closure_slots.resize_with(idx + 1, Vec::new);
        }
        self.smiles_start_atoms[idx] = smiles_start;
        idx
    }

    fn parse_molecule(&mut self) -> Result<Molecule, SmilesParseError> {
        self.skip_ascii_whitespace();
        let mut mol = Molecule::new();
        let first_atom = self.parse_atom()?;
        let current = self.add_atom_with_tracking(&mut mol, first_atom, true);
        self.parse_chain(&mut mol, current)?;
        while self.consume_if('.') {
            let fragment_atom = self.parse_atom()?;
            let fragment_idx = self.add_atom_with_tracking(&mut mol, fragment_atom, true);
            self.parse_chain(&mut mol, fragment_idx)?;
        }
        if !self.ring_opens.is_empty() {
            return Err(self.error("unclosed ring"));
        }
        Ok(mol)
    }

    fn reorder_ring_closure_bonds_rdkit_like(&mut self, mol: &mut Molecule) {
        if self.ring_closure_bonds.is_empty() {
            return;
        }

        let ring_set: HashSet<usize> = self.ring_closure_bonds.iter().copied().collect();
        let mut old_order = Vec::with_capacity(mol.bonds.len());
        old_order.extend(
            mol.bonds
                .iter()
                .map(|bond| bond.index)
                .filter(|idx| !ring_set.contains(idx)),
        );

        let mut closures = self.ring_closure_bonds.clone();
        closures.sort_by_key(|&bond_idx| {
            (
                self.ring_closure_bond_numbers
                    .get(bond_idx)
                    .and_then(|n| *n)
                    .unwrap_or(u32::MAX),
                bond_idx,
            )
        });
        old_order.extend(closures);

        let mut old_to_new = vec![usize::MAX; mol.bonds.len()];
        for (new_idx, &old_idx) in old_order.iter().enumerate() {
            old_to_new[old_idx] = new_idx;
        }

        let old_bonds = mol.bonds.clone();
        mol.bonds = old_order
            .iter()
            .enumerate()
            .map(|(new_idx, &old_idx)| {
                let mut bond = old_bonds[old_idx].clone();
                bond.index = new_idx;
                bond
            })
            .collect();

        self.ring_closure_bonds = self
            .ring_closure_bonds
            .iter()
            .map(|&old_idx| old_to_new[old_idx])
            .collect();

        let old_numbers = self.ring_closure_bond_numbers.clone();
        self.ring_closure_bond_numbers = vec![None; mol.bonds.len()];
        for (old_idx, maybe_number) in old_numbers.into_iter().enumerate() {
            if let Some(number) = maybe_number {
                self.ring_closure_bond_numbers[old_to_new[old_idx]] = Some(number);
            }
        }

        for slots in &mut self.atom_ring_closure_slots {
            for slot in slots {
                if let Some(old_idx) = *slot {
                    *slot = Some(old_to_new[old_idx]);
                }
            }
        }
    }

    fn parse_chain(
        &mut self,
        mol: &mut Molecule,
        mut current_atom: usize,
    ) -> Result<(), SmilesParseError> {
        loop {
            self.skip_ascii_whitespace();
            if self.is_eof() || self.peek_char() == Some(')') || self.peek_char() == Some('.') {
                return Ok(());
            }

            if self.peek_char() == Some('(') {
                self.consume_char();
                let branch_bond = self.parse_optional_bond();
                let branch_atom = self.parse_atom()?;
                let branch_atom_idx = self.add_atom_with_tracking(mol, branch_atom, false);
                self.add_resolved_bond(
                    mol,
                    current_atom,
                    branch_atom_idx,
                    branch_bond.unwrap_or(PendingBond::Unspecified),
                );
                self.parse_chain(mol, branch_atom_idx)?;
                self.expect_char(')')?;
                continue;
            }

            if let Some(ring_number) = self.parse_optional_ring_number() {
                self.handle_ring_closure(mol, current_atom, PendingBond::Unspecified, ring_number)?;
                continue;
            }

            if let Some(pending_bond) = self.parse_optional_bond() {
                if let Some(ring_number) = self.parse_optional_ring_number() {
                    self.handle_ring_closure(mol, current_atom, pending_bond, ring_number)?;
                    continue;
                }

                let next_atom = self.parse_atom()?;
                let next_idx = self.add_atom_with_tracking(mol, next_atom, false);
                self.add_resolved_bond(mol, current_atom, next_idx, pending_bond);
                current_atom = next_idx;
                continue;
            }

            let next_atom = self.parse_atom()?;
            let next_idx = self.add_atom_with_tracking(mol, next_atom, false);
            self.add_resolved_bond(mol, current_atom, next_idx, PendingBond::Unspecified);
            current_atom = next_idx;
        }
    }

    fn handle_ring_closure(
        &mut self,
        mol: &mut Molecule,
        current_atom: usize,
        bond: PendingBond,
        ring_number: u32,
    ) -> Result<(), SmilesParseError> {
        if let Some(open) = self.ring_opens.remove(&ring_number) {
            if open.atom_index == current_atom {
                return Err(self.error("duplicated ring closure bonds atom to itself"));
            }
            if mol.bonds.iter().any(|b| {
                (b.begin_atom == open.atom_index && b.end_atom == current_atom)
                    || (b.begin_atom == current_atom && b.end_atom == open.atom_index)
            }) {
                return Err(self.error("ring closure duplicates existing bond"));
            }
            let (begin_atom, end_atom, chosen) = if open.bond != PendingBond::Unspecified {
                (open.atom_index, current_atom, open.bond)
            } else {
                // RDKit CloseMolRings() keeps the first specified partial
                // bond. If the opening side was unspecified, the closing
                // partial bond is retained, so the final bond begins at the
                // current atom and ends at the opening atom.
                (current_atom, open.atom_index, bond)
            };
            let bi = self.add_resolved_bond(mol, begin_atom, end_atom, chosen);
            self.ring_closure_bonds.push(bi);
            if self.ring_closure_bond_numbers.len() <= bi {
                self.ring_closure_bond_numbers.resize(bi + 1, None);
            }
            self.ring_closure_bond_numbers[bi] = Some(ring_number);
            if self.atom_ring_closure_slots.len() <= current_atom {
                self.atom_ring_closure_slots
                    .resize_with(current_atom + 1, Vec::new);
            }
            self.atom_ring_closure_slots[current_atom].push(Some(bi));
            if self.atom_ring_closure_slots.len() <= open.atom_index {
                self.atom_ring_closure_slots
                    .resize_with(open.atom_index + 1, Vec::new);
            }
            if self.atom_ring_closure_slots[open.atom_index].len() <= open.opening_slot {
                self.atom_ring_closure_slots[open.atom_index].resize(open.opening_slot + 1, None);
            }
            self.atom_ring_closure_slots[open.atom_index][open.opening_slot] = Some(bi);
            Ok(())
        } else {
            // Mirrors RDKit SmilesParseOps::CheckRingClosureBranchStatus().
            // This handles chirality inversions required when a ring-opening
            // digit appears after branch constructs.
            let degree = mol
                .bonds
                .iter()
                .filter(|b| b.begin_atom == current_atom || b.end_atom == current_atom)
                .count();
            let is_tetra = matches!(
                mol.atoms[current_atom].chiral_tag,
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            );
            if current_atom != mol.atoms.len().saturating_sub(1)
                && (degree == 1
                    || (degree == 2 && current_atom != 0)
                    || (degree == 3 && current_atom == 0))
                && is_tetra
            {
                mol.atoms[current_atom].chiral_tag = match mol.atoms[current_atom].chiral_tag {
                    ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                    ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                    ChiralTag::Unspecified => ChiralTag::Unspecified,
                };
            }
            if self.atom_ring_closure_slots.len() <= current_atom {
                self.atom_ring_closure_slots
                    .resize_with(current_atom + 1, Vec::new);
            }
            let opening_slot = self.atom_ring_closure_slots[current_atom].len();
            self.atom_ring_closure_slots[current_atom].push(None);
            self.ring_opens.insert(
                ring_number,
                RingOpen {
                    atom_index: current_atom,
                    bond,
                    opening_slot,
                },
            );
            Ok(())
        }
    }

    fn add_resolved_bond(
        &self,
        mol: &mut Molecule,
        begin_atom: usize,
        end_atom: usize,
        pending: PendingBond,
    ) -> usize {
        let mut bgn = begin_atom;
        let mut end = end_atom;
        let (order, is_aromatic) = match pending {
            PendingBond::Single => (BondOrder::Single, false),
            PendingBond::Double => (BondOrder::Double, false),
            PendingBond::Triple => (BondOrder::Triple, false),
            PendingBond::Quadruple => (BondOrder::Quadruple, false),
            PendingBond::Aromatic => (BondOrder::Aromatic, true),
            PendingBond::DirectionalSingleUp | PendingBond::DirectionalSingleDown => {
                (BondOrder::Single, false)
            }
            PendingBond::DativeForward => (BondOrder::Dative, false),
            PendingBond::DativeBackward => {
                // RDKit keeps dative direction in begin/end topology:
                // "<-" means right atom donates to left atom.
                bgn = end_atom;
                end = begin_atom;
                (BondOrder::Dative, false)
            }
            PendingBond::Null => (BondOrder::Null, false),
            PendingBond::Unspecified => {
                let atom1 = &mol.atoms[begin_atom];
                let atom2 = &mol.atoms[end_atom];
                if atom1.is_aromatic && atom2.is_aromatic {
                    (BondOrder::Aromatic, true)
                } else {
                    (BondOrder::Single, false)
                }
            }
        };
        let direction = match pending {
            PendingBond::DirectionalSingleUp => BondDirection::EndUpRight,
            PendingBond::DirectionalSingleDown => BondDirection::EndDownRight,
            _ => BondDirection::None,
        };
        mol.add_bond(Bond {
            index: 0,
            begin_atom: bgn,
            end_atom: end,
            order,
            is_aromatic,
            direction,
            stereo: BondStereo::None,
            stereo_atoms: Vec::new(),
        })
    }

    fn parse_atom(&mut self) -> Result<Atom, SmilesParseError> {
        self.skip_ascii_whitespace();
        match self.peek_char() {
            Some('[') => self.parse_bracket_atom(),
            Some(_) => self.parse_simple_atom(),
            None => Err(self.error("expected atom")),
        }
    }

    fn parse_simple_atom(&mut self) -> Result<Atom, SmilesParseError> {
        if let Some((atomic_num, aromatic, consumed)) = self.match_simple_atom() {
            self.pos += consumed;
            return Ok(Atom {
                index: 0,
                atomic_num,
                is_aromatic: aromatic,
                formal_charge: 0,
                explicit_hydrogens: 0,
                no_implicit: false,
                num_radical_electrons: 0,
                chiral_tag: ChiralTag::Unspecified,
                isotope: None,
                atom_map_num: None,
            });
        }
        Err(self.error("unsupported atom token"))
    }

    fn parse_bracket_atom(&mut self) -> Result<Atom, SmilesParseError> {
        self.expect_char('[')?;
        let isotope = self.parse_optional_number().map(|v| v as u16);
        let mut atom = if self.consume_if('H') {
            Atom {
                index: 0,
                atomic_num: 1,
                is_aromatic: false,
                formal_charge: 0,
                explicit_hydrogens: 0,
                no_implicit: true,
                num_radical_electrons: 0,
                chiral_tag: ChiralTag::Unspecified,
                isotope,
                atom_map_num: None,
            }
        } else if self.consume_if('*') {
            Atom {
                index: 0,
                atomic_num: 0,
                is_aromatic: false,
                formal_charge: 0,
                explicit_hydrogens: 0,
                no_implicit: true,
                num_radical_electrons: 0,
                chiral_tag: ChiralTag::Unspecified,
                isotope,
                atom_map_num: None,
            }
        } else if let Some((atomic_num, aromatic, consumed)) = self.match_bracket_atom_symbol() {
            self.pos += consumed;
            Atom {
                index: 0,
                atomic_num,
                is_aromatic: aromatic,
                formal_charge: 0,
                explicit_hydrogens: 0,
                no_implicit: true,
                num_radical_electrons: 0,
                chiral_tag: ChiralTag::Unspecified,
                isotope,
                atom_map_num: None,
            }
        } else {
            return Err(self.error("unsupported bracket atom"));
        };

        if self.consume_if('@') {
            if self.consume_if('@') {
                atom.chiral_tag = ChiralTag::TetrahedralCw;
            } else {
                atom.chiral_tag = ChiralTag::TetrahedralCcw;
            }
        }
        if self.consume_if('H') {
            atom.explicit_hydrogens = self.parse_optional_number().unwrap_or(1) as u8;
        }
        if self.consume_if('+') {
            atom.formal_charge = if self.consume_if('+') {
                2
            } else {
                self.parse_optional_number().unwrap_or(1) as i8
            };
        } else if self.consume_if('-') {
            atom.formal_charge = if self.consume_if('-') {
                -2
            } else {
                -(self.parse_optional_number().unwrap_or(1) as i8)
            };
        }
        if self.consume_if(':') {
            atom.atom_map_num = Some(self.parse_required_number()? as u32);
        }
        self.expect_char(']')?;
        Ok(atom)
    }

    fn parse_optional_bond(&mut self) -> Option<PendingBond> {
        match self.peek_char()? {
            '<' if self.input[self.pos..].starts_with("<-") => {
                self.pos += 2;
                Some(PendingBond::DativeBackward)
            }
            '-' if self.input[self.pos..].starts_with("->") => {
                self.pos += 2;
                Some(PendingBond::DativeForward)
            }
            '-' => {
                self.consume_char();
                Some(PendingBond::Single)
            }
            '=' => {
                self.consume_char();
                Some(PendingBond::Double)
            }
            '#' => {
                self.consume_char();
                Some(PendingBond::Triple)
            }
            '$' => {
                self.consume_char();
                Some(PendingBond::Quadruple)
            }
            ':' => {
                self.consume_char();
                Some(PendingBond::Aromatic)
            }
            '~' => {
                self.consume_char();
                Some(PendingBond::Null)
            }
            '/' => {
                self.consume_char();
                Some(PendingBond::DirectionalSingleUp)
            }
            '\\' => {
                self.consume_char();
                Some(PendingBond::DirectionalSingleDown)
            }
            _ => None,
        }
    }

    fn parse_optional_ring_number(&mut self) -> Option<u32> {
        if self.consume_if('%') {
            if self.consume_if('(') {
                let number = self.parse_required_number().ok()?;
                if !self.consume_if(')') {
                    return None;
                }
                return Some(number as u32);
            }
            let d1 = self.parse_single_digit()? as u32;
            let d2 = self.parse_single_digit()? as u32;
            return Some(d1 * 10 + d2);
        }
        self.parse_single_digit().map(|v| v as u32)
    }

    fn parse_optional_number(&mut self) -> Option<u32> {
        let start = self.pos;
        let mut value = 0u32;
        let mut consumed = false;
        while let Some(ch) = self.peek_char() {
            if let Some(digit) = ch.to_digit(10) {
                consumed = true;
                value = value.saturating_mul(10).saturating_add(digit);
                self.consume_char();
            } else {
                break;
            }
        }
        if consumed {
            Some(value)
        } else {
            self.pos = start;
            None
        }
    }

    fn parse_required_number(&mut self) -> Result<u32, SmilesParseError> {
        self.parse_optional_number()
            .ok_or_else(|| self.error("expected number"))
    }

    fn parse_single_digit(&mut self) -> Option<u8> {
        let ch = self.peek_char()?;
        if ch.is_ascii_digit() {
            self.consume_char();
            return Some((ch as u8) - b'0');
        }
        None
    }

    fn match_simple_atom(&self) -> Option<(u8, bool, usize)> {
        let rest = &self.input[self.pos..];
        for (token, atomic_num, aromatic) in [
            ("Cl", 17, false),
            ("Br", 35, false),
            ("*", 0, false),
            ("B", 5, false),
            ("C", 6, false),
            ("N", 7, false),
            ("O", 8, false),
            ("P", 15, false),
            ("S", 16, false),
            ("F", 9, false),
            ("I", 53, false),
            ("b", 5, true),
            ("c", 6, true),
            ("n", 7, true),
            ("o", 8, true),
            ("p", 15, true),
            ("s", 16, true),
        ] {
            if rest.starts_with(token) {
                return Some((atomic_num, aromatic, token.len()));
            }
        }
        None
    }

    fn match_bracket_atom_symbol(&self) -> Option<(u8, bool, usize)> {
        let rest = &self.input[self.pos..];
        for (token, atomic_num, aromatic) in [
            ("se", 34, true),
            ("as", 33, true),
            ("te", 52, true),
            ("si", 14, true),
            ("Rh", 45, false),
            ("Cu", 29, false),
            ("Cl", 17, false),
            ("Br", 35, false),
            ("Si", 14, false),
            ("As", 33, false),
            ("Se", 34, false),
            ("Li", 3, false),
            ("Na", 11, false),
            ("Mg", 12, false),
            ("Al", 13, false),
            ("Ca", 20, false),
            ("Fe", 26, false),
            ("Zn", 30, false),
            ("*", 0, false),
            ("b", 5, true),
            ("c", 6, true),
            ("n", 7, true),
            ("o", 8, true),
            ("p", 15, true),
            ("s", 16, true),
            ("B", 5, false),
            ("C", 6, false),
            ("N", 7, false),
            ("O", 8, false),
            ("P", 15, false),
            ("S", 16, false),
            ("F", 9, false),
            ("I", 53, false),
        ] {
            if rest.starts_with(token) {
                return Some((atomic_num, aromatic, token.len()));
            }
        }
        None
    }

    fn expect_char(&mut self, expected: char) -> Result<(), SmilesParseError> {
        if self.consume_if(expected) {
            Ok(())
        } else {
            Err(self.error("unexpected token"))
        }
    }

    fn consume_if(&mut self, expected: char) -> bool {
        if self.peek_char() == Some(expected) {
            self.consume_char();
            true
        } else {
            false
        }
    }

    fn consume_char(&mut self) -> Option<char> {
        let ch = self.peek_char()?;
        self.pos += ch.len_utf8();
        Some(ch)
    }

    fn peek_char(&self) -> Option<char> {
        self.input[self.pos..].chars().next()
    }

    fn skip_ascii_whitespace(&mut self) {
        while matches!(self.peek_char(), Some(ch) if ch.is_ascii_whitespace()) {
            self.consume_char();
        }
    }

    fn is_eof(&self) -> bool {
        self.pos >= self.input.len()
    }

    fn error(&self, message: &str) -> SmilesParseError {
        let _ = self.pos;
        SmilesParseError::ParseError(message.to_owned())
    }
}
