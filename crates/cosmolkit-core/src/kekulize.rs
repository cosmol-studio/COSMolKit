use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};

use crate::valence::valence_list;
use crate::valence::{ValenceModel, assign_valence};
use crate::{BondOrder, Molecule};

/// Kekulization errors.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum KekulizeError {
    #[error("impossible aromatic assignment")]
    ImpossibleAromaticAssignment,
    #[error("kekulization path is not implemented")]
    NotImplemented,
}

/// Convert aromatic bond representation to one concrete alternating form.
pub fn kekulize_in_place(molecule: &mut Molecule) -> Result<(), KekulizeError> {
    if molecule.bonds.is_empty() {
        return Ok(());
    }

    let mut aromatic_bonds: Vec<usize> = Vec::new();
    for (bi, b) in molecule.bonds.iter().enumerate() {
        if matches!(b.order, BondOrder::Aromatic) {
            aromatic_bonds.push(bi);
        }
    }
    if aromatic_bonds.is_empty() {
        return Ok(());
    }

    let mut bond_lookup: HashMap<(usize, usize), usize> = HashMap::new();
    for (bi, b) in molecule.bonds.iter().enumerate() {
        let key = if b.begin_atom <= b.end_atom {
            (b.begin_atom, b.end_atom)
        } else {
            (b.end_atom, b.begin_atom)
        };
        bond_lookup.insert(key, bi);
    }
    let aromatic_set: HashSet<usize> = aromatic_bonds.iter().copied().collect();
    let aromatic_components: Vec<Vec<usize>> = fused_ring_atom_components(molecule)
        .into_iter()
        .filter(|comp| {
            let comp_set: HashSet<usize> = comp.iter().copied().collect();
            aromatic_set.iter().any(|&bi| {
                let bond = &molecule.bonds[bi];
                comp_set.contains(&bond.begin_atom) && comp_set.contains(&bond.end_atom)
            })
        })
        .collect();
    if aromatic_components.is_empty() {
        return Err(KekulizeError::ImpossibleAromaticAssignment);
    }

    let atom_ranks = rank_fragment_atoms_for_kekulize(molecule);

    for all_atms in aromatic_components {
        if !aromatic_component_has_cycle(molecule, &all_atms, &aromatic_set) {
            // RDKit rejects aromatic atoms/bonds outside ring systems in
            // sanitization/kekulization flow.
            return Err(KekulizeError::ImpossibleAromaticAssignment);
        }
        let mut d_bond_cands = vec![false; molecule.atoms.len()];
        let mut done = Vec::<usize>::new();
        let mut questions = Vec::<usize>::new();
        mark_dbond_cands(
            molecule,
            &all_atms,
            &aromatic_set,
            &mut d_bond_cands,
            &mut questions,
            &mut done,
        );

        let mut d_bond_adds = vec![false; molecule.bonds.len()];
        let ok = kekulize_worker(
            molecule,
            &all_atms,
            &aromatic_set,
            &bond_lookup,
            &atom_ranks,
            &mut d_bond_cands,
            &mut d_bond_adds,
            &mut done,
            100,
        );
        let success = if ok {
            true
        } else if !questions.is_empty() {
            permute_dummies_and_kekulize(
                molecule,
                &all_atms,
                &aromatic_set,
                &bond_lookup,
                &atom_ranks,
                &d_bond_cands,
                &questions,
                100,
            )
        } else {
            false
        };
        if !success {
            return Err(KekulizeError::ImpossibleAromaticAssignment);
        }
    }

    // clear aromatic flags in a "clearAromaticFlags=True"-like behavior
    // and ensure no aromatic bond type remains.
    for atom in &mut molecule.atoms {
        atom.is_aromatic = false;
    }
    for b in &mut molecule.bonds {
        if matches!(b.order, BondOrder::Aromatic) {
            b.order = BondOrder::Single;
        }
    }

    Ok(())
}

fn bond_order_contrib(order: BondOrder) -> i32 {
    match order {
        BondOrder::Single | BondOrder::Dative => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 1,
        BondOrder::Null => 0,
    }
}

fn bond_valence_contrib_for_atom(mol: &Molecule, bond_index: usize, atom_idx: usize) -> i32 {
    let b = &mol.bonds[bond_index];
    if b.begin_atom != atom_idx && b.end_atom != atom_idx {
        return 0;
    }
    if matches!(b.order, BondOrder::Dative) {
        // Mirror RDKit valence contribution for dative bonds:
        // donor(begin)=0, acceptor(end)=1.
        if b.end_atom == atom_idx {
            return 1;
        }
        return 0;
    }
    bond_order_contrib(b.order)
}

fn default_valence(atomic_num: u8) -> i32 {
    match atomic_num {
        0 => 0,
        5 => 3,  // B
        6 => 4,  // C
        7 => 3,  // N
        8 => 2,  // O
        14 => 4, // Si
        15 => 3, // P
        16 => 2, // S
        33 => 3, // As
        34 => 2, // Se
        _ => 0,
    }
}

fn is_early_atom(atomic_num: u8) -> bool {
    atomic_num <= 10
}

fn total_valence_for_atom(mol: &Molecule, atom_idx: usize) -> i32 {
    let mut total = mol.atoms[atom_idx].explicit_hydrogens as i32;
    for bi in 0..mol.bonds.len() {
        total += bond_valence_contrib_for_atom(mol, bi, atom_idx);
    }
    total
}

fn degree_for_atom(mol: &Molecule, atom_idx: usize) -> i32 {
    let mut degree = 0;
    for b in &mol.bonds {
        if b.begin_atom == atom_idx || b.end_atom == atom_idx {
            degree += 1;
        }
    }
    degree
}

fn neighbors_of(mol: &Molecule, atom_idx: usize) -> Vec<usize> {
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

fn all_bond_adjacency(mol: &Molecule) -> Vec<Vec<(usize, usize)>> {
    let mut adj = vec![Vec::<(usize, usize)>::new(); mol.atoms.len()];
    for (bi, b) in mol.bonds.iter().enumerate() {
        adj[b.begin_atom].push((b.end_atom, bi));
        adj[b.end_atom].push((b.begin_atom, bi));
    }
    adj
}

fn all_cycle_candidates(mol: &Molecule) -> Vec<Vec<usize>> {
    let adj = all_bond_adjacency(mol);
    let mut seen = HashSet::<Vec<usize>>::new();
    let mut rings = Vec::<Vec<usize>>::new();
    for (bi, bond) in mol.bonds.iter().enumerate() {
        let Some(path) = shortest_path_ignoring_edge(&adj, bond.begin_atom, bond.end_atom, bi)
        else {
            continue;
        };
        if path.len() < 3 {
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
    for (wi, &word) in words.iter().enumerate().rev() {
        if word != 0 {
            return Some(wi * 64 + (63usize - word.leading_zeros() as usize));
        }
    }
    None
}

fn cycle_bond_indices(mol: &Molecule, cycle: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(cycle.len());
    for i in 0..cycle.len() {
        let a = cycle[i];
        let b = cycle[(i + 1) % cycle.len()];
        if let Some((bi, _)) = mol.bonds.iter().enumerate().find(|(_, bond)| {
            (bond.begin_atom == a && bond.end_atom == b)
                || (bond.begin_atom == b && bond.end_atom == a)
        }) {
            out.push(bi);
        }
    }
    out
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
            for edge in edges {
                bits[edge / 64] ^= 1u64 << (edge % 64);
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
        for (pivot, basis_row) in &basis {
            if (row[pivot / 64] >> (pivot % 64)) & 1 == 1 {
                for wi in 0..row.len() {
                    row[wi] ^= basis_row[wi];
                }
            }
        }
        let Some(pivot) = highest_set_bit(&row) else {
            continue;
        };
        basis.push((pivot, row));
        basis.sort_by(|left, right| right.0.cmp(&left.0));
        selected.push(ring_idx);
        if selected.len() >= cyclomatic {
            break;
        }
    }
    selected
}

fn ring_neighbor_components(ring_bonds: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut neighbors = vec![Vec::<usize>::new(); ring_bonds.len()];
    for i in 0..ring_bonds.len() {
        let left: HashSet<usize> = ring_bonds[i].iter().copied().collect();
        for j in (i + 1)..ring_bonds.len() {
            if ring_bonds[j].iter().any(|bond| left.contains(bond)) {
                neighbors[i].push(j);
                neighbors[j].push(i);
            }
        }
    }

    let mut seen = vec![false; ring_bonds.len()];
    let mut components = Vec::<Vec<usize>>::new();
    for start in 0..ring_bonds.len() {
        if seen[start] {
            continue;
        }
        let mut q = VecDeque::new();
        q.push_back(start);
        seen[start] = true;
        let mut comp = Vec::new();
        while let Some(ring_idx) = q.pop_front() {
            comp.push(ring_idx);
            for &nbr in &neighbors[ring_idx] {
                if !seen[nbr] {
                    seen[nbr] = true;
                    q.push_back(nbr);
                }
            }
        }
        components.push(comp);
    }
    components
}

fn fused_ring_atom_components(mol: &Molecule) -> Vec<Vec<usize>> {
    // Source mapping: RDKit 2026.03.1 Kekulize.cpp::KekulizeFragment
    // uses findSSSR(), convertToBonds(), makeRingNeighborMap(), and
    // pickFusedRings() before calling kekulizeFused().
    let all_rings = all_cycle_candidates(mol);
    let basis_ids = reduce_to_min_cycle_basis_indices(mol, &all_rings);
    let rings: Vec<Vec<usize>> = basis_ids
        .into_iter()
        .map(|idx| all_rings[idx].clone())
        .collect();
    let ring_bonds: Vec<Vec<usize>> = rings
        .iter()
        .map(|ring| cycle_bond_indices(mol, ring))
        .collect();
    ring_neighbor_components(&ring_bonds)
        .into_iter()
        .map(|component| {
            let mut atoms = Vec::<usize>::new();
            for ring_idx in component {
                atoms.extend(rings[ring_idx].iter().copied());
            }
            atoms.sort_unstable();
            atoms.dedup();
            atoms
        })
        .collect()
}

#[derive(Debug, Clone)]
struct CanonBond {
    bond_type: u8,
    stereo: u8,
    nbr_sym_class: usize,
    nbr_idx: usize,
}

#[derive(Debug, Clone)]
struct CanonAtom {
    index: usize,
    degree: usize,
    total_num_hs: usize,
    atomic_num: u8,
    isotope: u16,
    formal_charge: i8,
    chiral_presence: bool,
    nbr_ids: Vec<usize>,
    bonds: Vec<CanonBond>,
}

fn rdkit_bond_type_rank(order: BondOrder) -> u8 {
    match order {
        BondOrder::Null => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        // RDKit Bond::AROMATIC has a larger enum value than ordinary bond
        // orders. Only relative ordering matters for canonical ranking.
        BondOrder::Aromatic => 12,
        BondOrder::Dative => 17,
    }
}

fn rdkit_bond_stereo_rank(stereo: crate::BondStereo) -> u8 {
    // RDKit 2026.03.1 Code/GraphMol/Bond.h::BondStereo:
    // STEREONONE=0, STEREOANY=1, STEREOZ=2, STEREOE=3,
    // STEREOCIS=4, STEREOTRANS=5.
    match stereo {
        crate::BondStereo::None => 0,
        crate::BondStereo::Any => 1,
        crate::BondStereo::Cis => 4,
        crate::BondStereo::Trans => 5,
    }
}

fn init_fragment_canon_atoms(mol: &Molecule) -> Vec<CanonAtom> {
    let mut atoms: Vec<CanonAtom> = mol
        .atoms
        .iter()
        .map(|atom| CanonAtom {
            index: atom.index,
            degree: 0,
            total_num_hs: atom.explicit_hydrogens as usize,
            atomic_num: atom.atomic_num,
            isotope: atom.isotope.unwrap_or(0),
            formal_charge: atom.formal_charge,
            chiral_presence: !matches!(atom.chiral_tag, crate::ChiralTag::Unspecified),
            nbr_ids: Vec::new(),
            bonds: Vec::new(),
        })
        .collect();

    for bond in &mol.bonds {
        let begin = bond.begin_atom;
        let end = bond.end_atom;
        atoms[begin].nbr_ids.push(end);
        atoms[end].nbr_ids.push(begin);
        atoms[begin].degree += 1;
        atoms[end].degree += 1;

        atoms[begin].bonds.push(CanonBond {
            bond_type: rdkit_bond_type_rank(bond.order),
            stereo: rdkit_bond_stereo_rank(bond.stereo),
            nbr_sym_class: 0,
            nbr_idx: end,
        });
        atoms[end].bonds.push(CanonBond {
            bond_type: rdkit_bond_type_rank(bond.order),
            stereo: rdkit_bond_stereo_rank(bond.stereo),
            nbr_sym_class: 0,
            nbr_idx: begin,
        });
    }

    for atom in &mut atoms {
        atom.bonds.sort_by(canon_bond_greater);
    }
    atoms
}

fn canon_bond_compare(lhs: &CanonBond, rhs: &CanonBond) -> Ordering {
    lhs.bond_type
        .cmp(&rhs.bond_type)
        .then_with(|| lhs.stereo.cmp(&rhs.stereo))
        .then_with(|| lhs.nbr_sym_class.cmp(&rhs.nbr_sym_class))
}

fn canon_bond_greater(lhs: &CanonBond, rhs: &CanonBond) -> Ordering {
    canon_bond_compare(rhs, lhs)
}

fn update_atom_neighbor_index(atoms: &mut [CanonAtom], atom_idx: usize) {
    let neighbor_classes: Vec<(usize, usize)> = atoms[atom_idx]
        .bonds
        .iter()
        .map(|bond| (bond.nbr_idx, atoms[bond.nbr_idx].index))
        .collect();
    for (bond, (_, nbr_class)) in atoms[atom_idx]
        .bonds
        .iter_mut()
        .zip(neighbor_classes.into_iter())
    {
        bond.nbr_sym_class = nbr_class;
    }
    atoms[atom_idx].bonds.sort_by(canon_bond_greater);
}

fn canon_atom_compare(atoms: &[CanonAtom], i: usize, j: usize) -> Ordering {
    atoms[i]
        .index
        .cmp(&atoms[j].index)
        .then_with(|| atoms[i].degree.cmp(&atoms[j].degree))
        .then_with(|| atoms[i].atomic_num.cmp(&atoms[j].atomic_num))
        .then_with(|| atoms[i].isotope.cmp(&atoms[j].isotope))
        .then_with(|| atoms[i].total_num_hs.cmp(&atoms[j].total_num_hs))
        .then_with(|| atoms[i].formal_charge.cmp(&atoms[j].formal_charge))
        .then_with(|| atoms[i].chiral_presence.cmp(&atoms[j].chiral_presence))
        .then_with(|| {
            let n = atoms[i].bonds.len().min(atoms[j].bonds.len());
            for k in 0..n {
                let cmp = canon_bond_compare(&atoms[i].bonds[k], &atoms[j].bonds[k]);
                if cmp != Ordering::Equal {
                    return cmp;
                }
            }
            atoms[i].bonds.len().cmp(&atoms[j].bonds.len())
        })
}

fn canon_atom_compare_i32(atoms: &mut [CanonAtom], i: usize, j: usize) -> i32 {
    update_atom_neighbor_index(atoms, i);
    update_atom_neighbor_index(atoms, j);
    match canon_atom_compare(atoms, i, j) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

fn hanoi_sort_partition(
    order: &mut [usize],
    count: &mut [usize],
    changed: &mut [bool],
    atoms: &mut [CanonAtom],
) {
    fn hanoi(
        base: &mut [usize],
        temp: &mut [usize],
        count: &mut [usize],
        changed: &mut [bool],
        atoms: &mut [CanonAtom],
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
                canon_atom_compare_i32(atoms, n1, n2)
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

        while i1 < n1_len && i2 < n2_len {
            let a = s1[i1];
            let b = s2[i2];
            let stat = if changed[a] || changed[b] {
                canon_atom_compare_i32(atoms, a, b)
            } else {
                0
            };
            let len1 = count[a];
            let len2 = count[b];
            debug_assert!(len1 > 0);
            debug_assert!(len2 > 0);

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

fn activate_partitions(
    order: &[usize],
    count: &[usize],
    next: &mut [isize],
    changed: &mut [bool],
) -> isize {
    next.fill(-2);
    let mut activeset = -1isize;
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
    atoms: &mut [CanonAtom],
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
        hanoi_sort_partition(start, count, changed, atoms);

        for k in offset..(offset + len) {
            changed[order[k]] = false;
        }

        let mut index = order[offset];
        let mut i = count[index];
        let mut symclass = 0usize;
        while i < len {
            index = order[offset + i];
            if count[index] != 0 {
                symclass = offset + i;
            }
            atoms[index].index = symclass;
            for &nbr in &atoms[index].nbr_ids {
                changed[nbr] = true;
            }
            i += 1;
        }

        index = order[offset];
        let mut i = count[index];
        while i < len {
            index = order[offset + i];
            for &nbr in &atoms[index].nbr_ids {
                touched[atoms[nbr].index] = true;
            }
            i += 1;
        }

        for idx in 0..atoms.len() {
            if touched[idx] {
                let partition = order[idx];
                if count[partition] > 1 && next[partition] == -2 {
                    next[partition] = activeset;
                    activeset = partition as isize;
                }
                touched[idx] = false;
            }
        }
    }
}

fn break_ties(
    atoms: &mut [CanonAtom],
    order: &mut [usize],
    count: &mut [usize],
    next: &mut [isize],
    changed: &mut [bool],
    touched: &mut [bool],
) {
    let n = atoms.len();
    let mut i = 0usize;
    while i < n {
        let partition = order[i];
        let old_part = atoms[partition].index;
        while count[partition] > 1 {
            let len = count[partition];
            let offset = atoms[partition].index + len - 1;
            let idx = order[offset];
            atoms[idx].index = offset;
            count[partition] = len - 1;
            count[idx] = 1;

            if atoms[idx].degree < 1 {
                continue;
            }
            for &nbr in &atoms[idx].nbr_ids {
                touched[atoms[nbr].index] = true;
                changed[nbr] = true;
            }

            let mut activeset = -1isize;
            for part_idx in 0..n {
                if touched[part_idx] {
                    let part = order[part_idx];
                    if count[part] > 1 && next[part] == -2 {
                        next[part] = activeset;
                        activeset = part as isize;
                    }
                    touched[part_idx] = false;
                }
            }
            refine_partitions(atoms, order, count, next, changed, touched, activeset);
        }
        if atoms[partition].index != old_part {
            i = i.saturating_sub(1);
        } else {
            i += 1;
        }
    }
}

fn rank_fragment_atoms_for_kekulize(mol: &Molecule) -> Vec<usize> {
    // Source mapping: RDKit 2026.03.1
    //   Code/GraphMol/Kekulize.cpp::KekulizeFragment(canonical=true)
    //   Code/GraphMol/new_canon.cpp::rankFragmentAtoms()
    //   Code/GraphMol/new_canon.h::{RefinePartitions,BreakTies}
    let n = mol.atoms.len();
    if n == 0 {
        return Vec::new();
    }

    let mut atoms = init_fragment_canon_atoms(mol);
    let mut order: Vec<usize> = (0..n).collect();
    let mut count = vec![0usize; n];
    let mut next = vec![-2isize; n];
    let mut changed = vec![true; n];
    let mut touched = vec![false; n];

    for atom in &mut atoms {
        atom.index = 0;
    }
    count[0] = n;

    let activeset = activate_partitions(&order, &count, &mut next, &mut changed);
    refine_partitions(
        &mut atoms,
        &mut order,
        &mut count,
        &mut next,
        &mut changed,
        &mut touched,
        activeset,
    );
    break_ties(
        &mut atoms,
        &mut order,
        &mut count,
        &mut next,
        &mut changed,
        &mut touched,
    );

    let mut ranks = vec![0usize; n];
    for idx in 0..n {
        ranks[idx] = atoms[idx].index;
    }
    ranks
}

fn bond_index_between(
    bond_lookup: &HashMap<(usize, usize), usize>,
    a: usize,
    b: usize,
) -> Option<usize> {
    let key = if a <= b { (a, b) } else { (b, a) };
    bond_lookup.get(&key).copied()
}

fn aromatic_component_bond_lookup(
    mol: &Molecule,
    all_atms: &[usize],
    aromatic_set: &HashSet<usize>,
) -> (Vec<Vec<(usize, usize)>>, HashMap<(usize, usize), usize>) {
    let mut local_adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); mol.atoms.len()];
    let in_all: HashSet<usize> = all_atms.iter().copied().collect();
    let mut bmap = HashMap::new();
    for &bi in aromatic_set {
        let b = &mol.bonds[bi];
        if !in_all.contains(&b.begin_atom) || !in_all.contains(&b.end_atom) {
            continue;
        }
        local_adj[b.begin_atom].push((b.end_atom, bi));
        local_adj[b.end_atom].push((b.begin_atom, bi));
        let key = if b.begin_atom <= b.end_atom {
            (b.begin_atom, b.end_atom)
        } else {
            (b.end_atom, b.begin_atom)
        };
        bmap.insert(key, bi);
    }
    (local_adj, bmap)
}

fn shortest_path_ignoring_edge(
    local_adj: &[Vec<(usize, usize)>],
    start: usize,
    goal: usize,
    ignored_bond: usize,
) -> Option<Vec<usize>> {
    let mut prev = vec![usize::MAX; local_adj.len()];
    let mut seen = vec![false; local_adj.len()];
    let mut q = VecDeque::new();
    seen[start] = true;
    q.push_back(start);
    while let Some(v) = q.pop_front() {
        if v == goal {
            break;
        }
        for &(nb, bi) in &local_adj[v] {
            if bi == ignored_bond || seen[nb] {
                continue;
            }
            seen[nb] = true;
            prev[nb] = v;
            q.push_back(nb);
        }
    }
    if !seen[goal] {
        return None;
    }
    let mut path = vec![goal];
    let mut cur = goal;
    while cur != start {
        cur = prev[cur];
        if cur == usize::MAX {
            return None;
        }
        path.push(cur);
    }
    path.reverse();
    Some(path)
}

fn aromatic_rings_for_component(
    mol: &Molecule,
    all_atms: &[usize],
    aromatic_set: &HashSet<usize>,
) -> Vec<Vec<usize>> {
    // We do not have RDKit RingInfo/SSSR state in this data model.
    // Build a deterministic cycle set from aromatic edges by:
    // for each aromatic edge (u,v), find shortest path u->v without that edge.
    // path + edge forms one simple cycle candidate.
    let (local_adj, _bmap) = aromatic_component_bond_lookup(mol, all_atms, aromatic_set);
    let mut rings = Vec::<Vec<usize>>::new();
    let mut seen = HashSet::<Vec<usize>>::new();
    for &bi in aromatic_set {
        let b = &mol.bonds[bi];
        if let Some(path) = shortest_path_ignoring_edge(&local_adj, b.begin_atom, b.end_atom, bi) {
            if path.len() < 3 {
                continue;
            }
            // Build canonical ring atom key from path (drop duplicated closure atom).
            let mut ring = path;
            ring.sort_unstable();
            ring.dedup();
            if ring.len() < 3 {
                continue;
            }
            if seen.insert(ring.clone()) {
                rings.push(ring);
            }
        }
    }

    // Empty ring-set is allowed here; caller decides whether component is valid.
    rings
}

fn aromatic_component_has_cycle(
    mol: &Molecule,
    all_atms: &[usize],
    aromatic_set: &HashSet<usize>,
) -> bool {
    let in_all: HashSet<usize> = all_atms.iter().copied().collect();
    let mut aromatic_atoms = HashSet::<usize>::new();
    let e = aromatic_set
        .iter()
        .filter(|&&bi| {
            let b = &mol.bonds[bi];
            let in_component = in_all.contains(&b.begin_atom) && in_all.contains(&b.end_atom);
            if in_component {
                aromatic_atoms.insert(b.begin_atom);
                aromatic_atoms.insert(b.end_atom);
            }
            in_component
        })
        .count();
    // For the aromatic subgraph inside the fused ring component: cycle iff
    // E >= V for one connected undirected component.
    e >= aromatic_atoms.len()
}

fn ring_not_candidates(
    mol: &Molecule,
    rings: &[Vec<usize>],
) -> (Vec<bool>, Vec<usize>, Vec<Vec<usize>>) {
    let mut num_atom_rings = vec![0usize; mol.atoms.len()];
    let mut atom_members = vec![Vec::<usize>::new(); mol.atoms.len()];
    for (ri, ring) in rings.iter().enumerate() {
        for &ai in ring {
            num_atom_rings[ai] += 1;
            atom_members[ai].push(ri);
        }
    }

    // RDKit markDbondCands:
    // ring is "not candidate" unless it has at least one aromatic atom
    // that belongs to exactly one ring.
    let mut is_ring_not_cand = vec![true; rings.len()];
    for (ri, ring) in rings.iter().enumerate() {
        for &ai in ring {
            let at = &mol.atoms[ai];
            if at.is_aromatic && num_atom_rings[ai] == 1 {
                is_ring_not_cand[ri] = false;
                break;
            }
        }
    }
    (is_ring_not_cand, num_atom_rings, atom_members)
}

fn rdkit_wedged_bond_priority_reachable(_mol: &Molecule) -> bool {
    // RDKit branch in kekulizeWorker uses BondDir BEGINWEDGE/BEGINDASH
    // to push wedged neighbors after normal options.
    // Current Bond model does not contain bond direction/stereo-wedge state.
    false
}

fn can_receive_double(
    mol: &Molecule,
    atom_idx: usize,
    aromatic_set: &HashSet<usize>,
    total_hs: &[i32],
) -> bool {
    let atom = &mol.atoms[atom_idx];
    if atom.atomic_num == 0 {
        return true;
    }
    if !atom.is_aromatic {
        return false;
    }

    let mut base_valence = total_hs[atom_idx];
    let mut n_to_ignore = 0i32;
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.begin_atom != atom_idx && b.end_atom != atom_idx {
            continue;
        }
        if aromatic_set.contains(&bi) {
            base_valence += 1;
        } else {
            let contrib = bond_valence_contrib_for_atom(mol, bi, atom_idx);
            base_valence += contrib;
            if contrib == 0 {
                n_to_ignore += 1;
            }
        }
    }

    let mut dv = default_valence(atom.atomic_num);
    if dv <= 0 {
        return false;
    }

    let mut chrg = atom.formal_charge as i32;
    if is_early_atom(atom.atomic_num) {
        chrg = -chrg;
    }
    if atom.atomic_num == 6 && chrg > 0 {
        chrg = -chrg;
    }
    dv += chrg;

    let total_valence = total_valence_for_atom(mol, atom_idx);
    let n_radicals = atom.num_radical_electrons as i32;
    let total_degree = degree_for_atom(mol, atom_idx) + total_hs[atom_idx] - n_to_ignore;

    if let Some(vlist) = valence_list(atom.atomic_num) {
        let mut vi = 1usize;
        while total_valence > dv && vi < vlist.len() && vlist[vi] > 0 {
            dv = vlist[vi] + chrg;
            vi += 1;
        }
    }

    // RDKit special case for aromatic N-oxides:
    // O=n1ccccc1 family can require dv=5 on N/P/As.
    if total_valence == 5
        && base_valence == 4
        && dv == 3
        && total_degree == 3
        && n_radicals == 0
        && chrg == 0
        && total_hs[atom_idx] == 0
        && matches!(atom.atomic_num, 7 | 15 | 33)
    {
        dv = 5;
    }

    if total_degree + n_radicals >= dv {
        return false;
    }
    if dv == base_valence + 1 + n_radicals {
        return true;
    }
    if n_radicals == 0 && atom.no_implicit && dv == base_valence + 2 {
        return true;
    }
    false
}

fn mark_dbond_cands(
    mol: &mut Molecule,
    all_atms: &[usize],
    aromatic_set: &HashSet<usize>,
    d_bond_cands: &mut [bool],
    questions: &mut Vec<usize>,
    done: &mut Vec<usize>,
) {
    let has_aromatic_or_dummy = all_atms.iter().any(|&ai| {
        let at = &mol.atoms[ai];
        at.atomic_num == 0 || at.is_aromatic
    });
    if !has_aromatic_or_dummy {
        return;
    }
    if rdkit_wedged_bond_priority_reachable(mol) {
        unimplemented!(
            "RDKit wedged-bond option ordering in KekulizeWorker is not representable without bond direction field"
        );
    }

    let in_all: HashSet<usize> = all_atms.iter().copied().collect();
    let rings = aromatic_rings_for_component(mol, all_atms, aromatic_set);
    let (is_ring_not_cand, num_atom_rings, atom_members) = ring_not_candidates(mol, &rings);
    let mut make_single: Vec<usize> = Vec::new();
    let total_hs: Vec<i32> = match assign_valence(mol, ValenceModel::RdkitLike) {
        Ok(v) => mol
            .atoms
            .iter()
            .enumerate()
            .map(|(i, at)| at.explicit_hydrogens as i32 + v.implicit_hydrogens[i] as i32)
            .collect(),
        Err(_) => mol
            .atoms
            .iter()
            .map(|at| at.explicit_hydrogens as i32)
            .collect(),
    };

    for &ai in all_atms {
        let at = &mol.atoms[ai];
        if at.atomic_num != 0 && !at.is_aromatic {
            done.push(ai);
            continue;
        }

        let mut non_ar_non_dummy_nbr = 0usize;
        for b in &mol.bonds {
            if b.begin_atom != ai && b.end_atom != ai {
                continue;
            }
            let nbr = if b.begin_atom == ai {
                b.end_atom
            } else {
                b.begin_atom
            };
            if in_all.contains(&nbr) {
                let other = &mol.atoms[nbr];
                if other.atomic_num != 0 && !other.is_aromatic {
                    non_ar_non_dummy_nbr += 1;
                }
            }
            if aromatic_set.contains(&b.index)
                && matches!(
                    b.order,
                    BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
                )
            {
                make_single.push(b.index);
            }
        }

        if at.atomic_num == 0 {
            let n_atom_rings = num_atom_rings[ai];
            let n_non_cand_rings = atom_members[ai]
                .iter()
                .filter(|&&ri| is_ring_not_cand[ri])
                .count();
            if non_ar_non_dummy_nbr < n_atom_rings && n_non_cand_rings < n_atom_rings {
                d_bond_cands[ai] = true;
                questions.push(ai);
            }
        } else if can_receive_double(mol, ai, aromatic_set, &total_hs) {
            d_bond_cands[ai] = true;
        }
    }

    for bi in make_single {
        let b = &mut mol.bonds[bi];
        if in_all.contains(&b.begin_atom) && in_all.contains(&b.end_atom) {
            b.order = BondOrder::Single;
        }
    }
}

fn backtrack(
    mol: &mut Molecule,
    options: &mut HashMap<usize, VecDeque<usize>>,
    last_opt: usize,
    done: &mut Vec<usize>,
    aqueue: &mut VecDeque<usize>,
    d_bond_cands: &mut [bool],
    d_bond_adds: &mut [bool],
) {
    let Some(pos) = done.iter().position(|x| *x == last_opt) else {
        return;
    };
    let tdone: Vec<usize> = done[..pos].to_vec();
    let rollback: Vec<usize> = done[pos..].to_vec();
    for &a in rollback.iter().rev() {
        aqueue.push_front(a);
    }
    let tdone_set: HashSet<usize> = tdone.iter().copied().collect();
    for bi in 0..mol.bonds.len() {
        if !d_bond_adds[bi] {
            continue;
        }
        let b = &mol.bonds[bi];
        if !tdone_set.contains(&b.begin_atom) && !tdone_set.contains(&b.end_atom) {
            d_bond_adds[bi] = false;
            let bb = &mut mol.bonds[bi];
            bb.order = BondOrder::Single;
            d_bond_cands[bb.begin_atom] = true;
            d_bond_cands[bb.end_atom] = true;
        }
    }
    *done = tdone;
    options.retain(|k, _| done.contains(k) || *k == last_opt);
}

fn kekulize_worker(
    mol: &mut Molecule,
    all_atms: &[usize],
    aromatic_set: &HashSet<usize>,
    bond_lookup: &HashMap<(usize, usize), usize>,
    atom_ranks: &[usize],
    d_bond_cands: &mut [bool],
    d_bond_adds: &mut [bool],
    done: &mut Vec<usize>,
    max_backtracks: usize,
) -> bool {
    let all_set: HashSet<usize> = all_atms.iter().copied().collect();
    let mut astack: VecDeque<usize> = VecDeque::new();
    let mut options: HashMap<usize, VecDeque<usize>> = HashMap::new();
    let mut btmoves: Vec<usize> = Vec::new();
    let mut last_opt: Option<usize> = None;
    let mut local_bonds_added = vec![false; mol.bonds.len()];
    let mut num_bt = 0usize;
    let mut sorted_atms = all_atms.to_vec();
    sorted_atms.sort_by_key(|&a| (atom_ranks.get(a).copied().unwrap_or(usize::MAX), a));

    while done.len() < sorted_atms.len() || !astack.is_empty() {
        let curr = if let Some(v) = astack.pop_front() {
            v
        } else {
            let mut next_start = None;
            for &a in &sorted_atms {
                if !done.contains(&a) {
                    next_start = Some(a);
                    break;
                }
            }
            match next_start {
                Some(v) => v,
                None => return true,
            }
        };
        done.push(curr);

        let c_cand = d_bond_cands[curr];
        let mut existing_option_frame = false;
        let mut opts = if let Some(saved) = options.get(&curr) {
            existing_option_frame = true;
            saved.clone()
        } else {
            let mut lstack: VecDeque<usize> = VecDeque::new();
            let mut new_opts: VecDeque<usize> = VecDeque::new();
            let mut neighbors = neighbors_of(mol, curr);
            neighbors.sort_by_key(|&nbr| (atom_ranks.get(nbr).copied().unwrap_or(usize::MAX), nbr));
            for nbr in neighbors {
                if done.contains(&nbr) || !all_set.contains(&nbr) {
                    continue;
                }
                if !astack.contains(&nbr) {
                    lstack.push_back(nbr);
                }
                if !c_cand || !d_bond_cands[nbr] {
                    continue;
                }
                let Some(bi) = bond_index_between(bond_lookup, curr, nbr) else {
                    continue;
                };
                let ok_edge = aromatic_set.contains(&bi)
                    || mol.atoms[curr].atomic_num == 0
                    || mol.atoms[nbr].atomic_num == 0;
                if ok_edge {
                    new_opts.push_back(nbr);
                }
            }
            astack.extend(lstack);
            new_opts
        };

        if c_cand {
            if let Some(ncnd) = opts.pop_front() {
                let Some(bi) = bond_index_between(bond_lookup, curr, ncnd) else {
                    return false;
                };
                mol.bonds[bi].order = BondOrder::Double;
                d_bond_cands[curr] = false;
                d_bond_cands[ncnd] = false;
                d_bond_adds[bi] = true;
                local_bonds_added[bi] = true;

                if existing_option_frame {
                    if opts.is_empty() {
                        options.remove(&curr);
                        btmoves.pop();
                        last_opt = btmoves.last().copied();
                    } else {
                        options.insert(curr, opts);
                    }
                } else if !opts.is_empty() {
                    last_opt = Some(curr);
                    btmoves.push(curr);
                    options.insert(curr, opts);
                }
            } else if mol.atoms[curr].atomic_num != 0 {
                if let Some(lo) = last_opt {
                    if num_bt < max_backtracks {
                        backtrack(
                            mol,
                            &mut options,
                            lo,
                            done,
                            &mut astack,
                            d_bond_cands,
                            d_bond_adds,
                        );
                        num_bt += 1;
                    } else {
                        for (bi, added) in local_bonds_added.iter().copied().enumerate() {
                            if added {
                                mol.bonds[bi].order = BondOrder::Single;
                            }
                        }
                        return false;
                    }
                } else {
                    for (bi, added) in local_bonds_added.iter().copied().enumerate() {
                        if added {
                            mol.bonds[bi].order = BondOrder::Single;
                        }
                    }
                    return false;
                }
            }
        }
    }
    true
}

fn permute_dummies_and_kekulize(
    mol: &mut Molecule,
    all_atms: &[usize],
    aromatic_set: &HashSet<usize>,
    bond_lookup: &HashMap<(usize, usize), usize>,
    atom_ranks: &[usize],
    d_bond_cands: &[bool],
    questions: &[usize],
    max_backtracks: usize,
) -> bool {
    if questions.is_empty() {
        return false;
    }
    let all_set: HashSet<usize> = all_atms.iter().copied().collect();
    let mut mask = 1usize;
    while mask < (1usize << questions.len()) {
        // reset aromatic edges in this fragment
        for (bi, b) in mol.bonds.iter_mut().enumerate() {
            if aromatic_set.contains(&bi)
                && all_set.contains(&b.begin_atom)
                && all_set.contains(&b.end_atom)
            {
                if !matches!(b.order, BondOrder::Single) {
                    b.order = BondOrder::Single;
                }
            }
        }

        let mut t_cands = d_bond_cands.to_vec();
        for (i, &q) in questions.iter().enumerate() {
            if (mask & (1usize << i)) != 0 {
                t_cands[q] = false;
            }
        }
        let mut d_bond_adds = vec![false; mol.bonds.len()];
        let mut done = Vec::new();
        let ok = kekulize_worker(
            mol,
            all_atms,
            aromatic_set,
            bond_lookup,
            atom_ranks,
            &mut t_cands,
            &mut d_bond_adds,
            &mut done,
            max_backtracks,
        );
        if ok {
            return true;
        }
        mask += 1;
    }
    false
}
