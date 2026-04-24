use std::collections::{HashMap, HashSet, VecDeque};

use crate::valence::valence_list;
use crate::valence::{ValenceModel, assign_valence};
use crate::{BondOrder, Molecule};

/// Kekulization errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KekulizeError {
    ImpossibleAromaticAssignment,
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

    let mut aromatic_atom_graph: Vec<Vec<usize>> = vec![Vec::new(); molecule.atoms.len()];
    for &bi in &aromatic_bonds {
        let b = &molecule.bonds[bi];
        aromatic_atom_graph[b.begin_atom].push(b.end_atom);
        aromatic_atom_graph[b.end_atom].push(b.begin_atom);
    }

    let mut visited = vec![false; molecule.atoms.len()];
    let mut aromatic_components: Vec<Vec<usize>> = Vec::new();
    for start in 0..molecule.atoms.len() {
        if visited[start] || aromatic_atom_graph[start].is_empty() {
            continue;
        }
        let mut stack = vec![start];
        visited[start] = true;
        let mut comp = Vec::new();
        while let Some(v) = stack.pop() {
            comp.push(v);
            for &nb in &aromatic_atom_graph[v] {
                if !visited[nb] {
                    visited[nb] = true;
                    stack.push(nb);
                }
            }
        }
        comp.sort_unstable();
        aromatic_components.push(comp);
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
    let e = aromatic_set
        .iter()
        .filter(|&&bi| {
            let b = &mol.bonds[bi];
            in_all.contains(&b.begin_atom) && in_all.contains(&b.end_atom)
        })
        .count();
    // For one connected undirected component: cycle iff E >= V.
    e >= all_atms.len()
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

    while done.len() < all_atms.len() || !astack.is_empty() {
        let curr = if let Some(v) = astack.pop_front() {
            v
        } else {
            let mut fallback = None;
            for &a in all_atms {
                if !done.contains(&a) {
                    fallback = Some(a);
                    break;
                }
            }
            match fallback {
                Some(v) => v,
                None => return true,
            }
        };
        done.push(curr);

        let c_cand = d_bond_cands[curr];
        let mut opts = if let Some(saved) = options.remove(&curr) {
            saved
        } else {
            let mut lstack: VecDeque<usize> = VecDeque::new();
            let mut new_opts: VecDeque<usize> = VecDeque::new();
            let neighbors = neighbors_of(mol, curr);
            for nbr in neighbors {
                if done.contains(&nbr) || !all_set.contains(&nbr) {
                    continue;
                }
                if !astack.contains(&nbr) {
                    lstack.push_front(nbr);
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

                if !opts.is_empty() {
                    last_opt = Some(curr);
                    if btmoves.last().copied() != Some(curr) {
                        btmoves.push(curr);
                    }
                    options.insert(curr, opts);
                } else if btmoves.last().copied() == Some(curr) {
                    let _ = btmoves.pop();
                    last_opt = btmoves.last().copied();
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
