use super::*;

pub(super) fn canonicalize_double_bond_directions_for_writer(
    topology: &mut TopologyBlock,
    stack: &[MolStackElem],
    traversal_ring_closure_bonds: &[bool],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment double-bond direction section
    // RDKit✔️✔️: std::vector<unsigned int> atomVisitOrders(mol.getNumAtoms());
    // RDKit✔️✔️: std::vector<unsigned int> bondVisitOrders(mol.getNumBonds());
    // RDKit✔️✔️: for (const auto &msI : molStack) {
    // RDKit✔️✔️:   if (msI.type == MOL_STACK_ATOM) { atomVisitOrders[msI.obj.atom->getIdx()] = pos; }
    // RDKit✔️✔️:   else if (msI.type == MOL_STACK_BOND) {
    // RDKit✔️✔️:     bondVisitOrders[msI.obj.bond->getIdx()] = pos;
    // RDKit✔️✔️:     auto dir = msI.obj.bond->getBondDir();
    // RDKit✔️✔️:     if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
    // RDKit✔️✔️:       msI.obj.bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: canonicalizeDoubleBonds(mol, bondVisitOrders, atomVisitOrders,
    // RDKit✔️✔️:                         bondDirCounts, atomDirCounts, molStack);
    // RDKit✔️✔️: Canon::removeUnwantedBondDirSpecs(mol, molStack, bondDirCounts,
    // RDKit✔️✔️:                                     atomDirCounts, bondVisitOrders);
    // RDKit✔️✔️: Canon::removeRedundantBondDirSpecs(mol, molStack, bondDirCounts,
    // RDKit✔️✔️:                                    atomDirCounts);
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment double-bond direction section
    let mut atom_visit_orders = vec![usize::MAX; topology.atoms.len()];
    let mut bond_visit_orders = vec![usize::MAX; topology.bonds.len()];
    for (pos, item) in stack.iter().enumerate() {
        match *item {
            MolStackElem::Atom(atom) => atom_visit_orders[atom] = pos,
            MolStackElem::Bond { bond, .. } => {
                bond_visit_orders[bond.index()] = pos;
                if matches!(
                    topology.bonds[bond.index()].direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                ) {
                    topology.bonds[bond.index()].set_direction(BondDirection::None);
                }
            }
            MolStackElem::Ring(_) => {}
            MolStackElem::BranchOpen | MolStackElem::BranchClose => {}
        }
    }

    let mut bond_dir_counts = vec![0i8; topology.bonds.len()];
    let mut atom_dir_counts = vec![0i8; topology.atoms.len()];
    canonicalize_double_bonds_for_writer(
        topology,
        &bond_visit_orders,
        &atom_visit_orders,
        &traversal_ring_closure_bonds,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        stack,
    );
    remove_unwanted_bond_dir_specs_for_writer(
        topology,
        stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        &bond_visit_orders,
    );
    remove_redundant_bond_dir_specs_for_writer(
        topology,
        stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
    );
    Ok(())
}

pub(super) fn canonicalize_double_bonds_for_writer(
    topology: &mut TopologyBlock,
    bond_visit_orders: &[usize],
    atom_visit_orders: &[usize],
    traversal_ring_closure_bonds: &[bool],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
    stack: &[MolStackElem],
) {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeDoubleBonds
    // RDKit✔️✔️: for (auto &msI : molStack) {
    // RDKit✔️✔️:   if (msI.type != MOL_STACK_BOND) { continue; }
    // RDKit✔️✔️:   if (bond->getBondType() != Bond::DOUBLE ||
    // RDKit✔️✔️:       bond->getStereo() <= Bond::STEREOANY ||
    // RDKit✔️✔️:       bond->getStereoAtoms().size() < 2) {
    // RDKit✔️✔️:     bond->setStereo(Bond::STEREONONE);
    // RDKit✔️✔️:     continue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ... prioritize by neighboring stereo bonds and molStack order ...
    // RDKit✔️✔️: }
    // RDKit✔️✔️: while (!q.empty()) { Canon::canonicalizeDoubleBond(...); }
    // END RDKIT CPP FUNCTION Canon::canonicalizeDoubleBonds
    let mut stereo_bond_neighbors = BTreeMap::<BondId, Vec<BondId>>::new();
    let mut candidates = Vec::new();
    for item in stack {
        let MolStackElem::Bond { bond, .. } = *item else {
            continue;
        };
        let bond_ref = &topology.bonds[bond.index()];
        if !is_writer_stereo_double_bond(bond_ref) {
            if bond_ref.order() == BondOrder::Double {
                let bond_mut = &mut topology.bonds[bond.index()];
                bond_mut.set_stereo_atoms(None);
                bond_mut.set_stereo(BondStereo::None);
            }
            continue;
        }
        let mut stereo_nbrs = neighboring_stereo_double_bonds_for_writer(topology, bond);
        stereo_nbrs.sort_by_key(|neighbor| bond_visit_orders[neighbor.index()]);
        stereo_bond_neighbors.insert(bond, stereo_nbrs.clone());
        candidates.push((
            usize::MAX - stereo_nbrs.len(),
            bond_visit_orders[bond.index()],
            bond,
        ));
    }
    candidates.sort_by_key(|(stereo_rank, visit_order, _)| (*stereo_rank, *visit_order));

    let mut seen = vec![false; topology.bonds.len()];
    for (_, _, start_bond) in candidates {
        if seen[start_bond.index()] {
            continue;
        }
        let mut queue = std::collections::VecDeque::from([start_bond]);
        while let Some(bond) = queue.pop_front() {
            if seen[bond.index()] {
                continue;
            }
            canonicalize_double_bond_for_writer(
                topology,
                bond,
                bond_visit_orders,
                atom_visit_orders,
                traversal_ring_closure_bonds,
                bond_dir_counts,
                atom_dir_counts,
            );
            seen[bond.index()] = true;
            for &nbr in stereo_bond_neighbors.get(&bond).into_iter().flatten() {
                if !seen[nbr.index()] {
                    queue.push_back(nbr);
                }
            }
        }
    }
}

pub(super) fn canonicalize_double_bond_for_writer(
    topology: &mut TopologyBlock,
    dbl_bond: BondId,
    bond_visit_orders: &[usize],
    atom_visit_orders: &[usize],
    traversal_ring_closure_bonds: &[bool],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeDoubleBond
    // RDKit✔️✔️: void canonicalizeDoubleBond(Bond *dblBond,
    // RDKit✔️✔️:   const UINT_VECT &bondVisitOrders, const UINT_VECT &atomVisitOrders,
    // RDKit✔️✔️:   std::vector<int8_t> &bondDirCounts,
    // RDKit✔️✔️:   std::vector<int8_t> &atomDirCounts) {
    // RDKit✔️✔️:   Atom *atom1 = dblBond->getBeginAtom();
    // RDKit✔️✔️:   Atom *atom2 = dblBond->getEndAtom();
    // RDKit✔️✔️:   if ((atom1->getDegree() != 2 && atom1->getDegree() != 3) ||
    // RDKit✔️✔️:       (atom2->getDegree() != 2 && atom2->getDegree() != 3)) { return; }
    // RDKit✔️✔️:   if (atomVisitOrders[dblBond->getBeginAtomIdx()] >=
    // RDKit✔️✔️:       atomVisitOrders[dblBond->getEndAtomIdx()]) { std::swap(atom1, atom2); }
    // RDKit✔️✔️:   ... find lowest visit order neighbor bonds from both ends ...
    // RDKit✔️✔️:   bool isFirstFromAtom1Flipped = [&]() {
    // RDKit✔️✔️:     auto anchorIdx = firstFromAtom1->getOtherAtom(atom1)->getIdx();
    // RDKit✔️✔️:     return (atomVisitOrders[atom1->getIdx()] < atomVisitOrders[anchorIdx]) !=
    // RDKit✔️✔️:            firstFromAtom1->hasProp(
    // RDKit✔️✔️:                common_properties::_TraversalRingClosureBond);
    // RDKit✔️✔️:   }();
    // RDKit✔️✔️:   bool isSecondFromAtom1Flipped = [&]() {
    // RDKit✔️✔️:     if (secondFromAtom1 == nullptr) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto anchorIdx = secondFromAtom1->getOtherAtom(atom1)->getIdx();
    // RDKit✔️✔️:     return (atomVisitOrders[atom1->getIdx()] < atomVisitOrders[anchorIdx]) !=
    // RDKit✔️✔️:            secondFromAtom1->hasProp(
    // RDKit✔️✔️:                common_properties::_TraversalRingClosureBond);
    // RDKit✔️✔️:   }();
    // RDKit✔️✔️:   bool isFirstFromAtom2Flipped = [&]() {
    // RDKit✔️✔️:     auto anchorIdx = firstFromAtom2->getOtherAtom(atom2)->getIdx();
    // RDKit✔️✔️:     return (atomVisitOrders[anchorIdx] < atomVisitOrders[atom2->getIdx()]) !=
    // RDKit✔️✔️:            firstFromAtom2->hasProp(
    // RDKit✔️✔️:                common_properties::_TraversalRingClosureBond);
    // RDKit✔️✔️:   }();
    // RDKit✔️✔️:   bool isSecondFromAtom2Flipped = [&]() {
    // RDKit✔️✔️:     if (secondFromAtom2 == nullptr) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto anchorIdx = secondFromAtom2->getOtherAtom(atom2)->getIdx();
    // RDKit✔️✔️:     return (atomVisitOrders[anchorIdx] < atomVisitOrders[atom2->getIdx()]) !=
    // RDKit✔️✔️:            secondFromAtom2->hasProp(
    // RDKit✔️✔️:                common_properties::_TraversalRingClosureBond);
    // RDKit✔️✔️:   }();
    // RDKit✔️✔️:   ... set arbitrary ENDUPRIGHT if neither side is constrained ...
    // RDKit✔️✔️:   ... propagate direction with getReferenceDirection ...
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::canonicalizeDoubleBond
    if !is_writer_stereo_double_bond(&topology.bonds[dbl_bond.index()]) {
        return;
    }
    let mut atom1 = topology.bonds[dbl_bond.index()].begin();
    let mut atom2 = topology.bonds[dbl_bond.index()].end();
    let deg1 = incident_bonds(topology, atom1).len();
    let deg2 = incident_bonds(topology, atom2).len();
    if !(deg1 == 2 || deg1 == 3) || !(deg2 == 2 || deg2 == 3) {
        return;
    }
    if atom_visit_orders[atom1.index()] >= atom_visit_orders[atom2.index()] {
        std::mem::swap(&mut atom1, &mut atom2);
    }

    let (first1, second1, dir1_set) = find_double_bond_neighbor_bonds_for_writer(
        topology,
        dbl_bond,
        atom1,
        bond_visit_orders,
        bond_dir_counts,
    );
    let (first2, second2, dir2_set) = find_double_bond_neighbor_bonds_for_writer(
        topology,
        dbl_bond,
        atom2,
        bond_visit_orders,
        bond_dir_counts,
    );
    let (Some(first1), Some(first2)) = (first1, first2) else {
        return;
    };

    let first1_flipped = writer_bond_is_flipped_from_atom1(
        topology,
        atom1,
        first1,
        atom_visit_orders,
        traversal_ring_closure_bonds,
    );
    let second1_flipped = second1
        .map(|bond| {
            writer_bond_is_flipped_from_atom1(
                topology,
                atom1,
                bond,
                atom_visit_orders,
                traversal_ring_closure_bonds,
            )
        })
        .unwrap_or(false);
    let first2_flipped = writer_bond_is_flipped_from_atom2(
        topology,
        atom2,
        first2,
        atom_visit_orders,
        traversal_ring_closure_bonds,
    );
    let second2_flipped = second2
        .map(|bond| {
            writer_bond_is_flipped_from_atom2(
                topology,
                atom2,
                bond,
                atom_visit_orders,
                traversal_ring_closure_bonds,
            )
        })
        .unwrap_or(false);

    if dir1_set && dir2_set {
        let atom1_consistent = account_same_side_dirs_for_writer(
            topology,
            atom1,
            first1,
            first1_flipped,
            second1,
            second1_flipped,
            bond_dir_counts,
            atom_dir_counts,
        );
        let atom2_consistent = account_same_side_dirs_for_writer(
            topology,
            atom2,
            first2,
            first2_flipped,
            second2,
            second2_flipped,
            bond_dir_counts,
            atom_dir_counts,
        );
        let _ = handle_dir_conflicts_across_double_bond_for_writer(
            topology,
            dbl_bond,
            atom1,
            atom1_consistent,
            first1,
            first1_flipped,
            second1,
            second1_flipped,
            atom2,
            atom2_consistent,
            first2,
            first2_flipped,
            second2,
            second2_flipped,
            bond_dir_counts,
            atom_dir_counts,
        );
        return;
    }

    let mut set_from_atom1 = true;
    let mut atom1_controlling_bond = first1;
    let mut atom2_controlling_bond = first2;
    if !dir1_set && !dir2_set {
        set_bond_direction(topology, first1, BondDirection::EndUpRight);
        bond_dir_counts[first1.index()] += 1;
        atom_dir_counts[atom1.index()] += 1;
    } else if !dir2_set {
        if bond_dir_counts[first1.index()] > 0 {
            bond_dir_counts[first1.index()] += 1;
            atom_dir_counts[atom1.index()] += 1;
            if let Some(second1) = second1
                && bond_dir_counts[second1.index()] > 0
            {
                bond_dir_counts[second1.index()] += 1;
                atom_dir_counts[atom1.index()] += 1;
            }
        } else if let Some(second1) = second1 {
            set_direction_from_neighboring_bond_for_writer(
                topology,
                second1,
                second1_flipped,
                first1,
                first1_flipped,
            );
            bond_dir_counts[second1.index()] += 1;
            bond_dir_counts[first1.index()] += 1;
            atom_dir_counts[atom1.index()] += 2;
            atom1_controlling_bond = second1;
        }
    } else {
        set_from_atom1 = false;
        if bond_dir_counts[first2.index()] > 0 {
            bond_dir_counts[first2.index()] += 1;
            atom_dir_counts[atom2.index()] += 1;
            if let Some(second2) = second2
                && bond_dir_counts[second2.index()] > 0
            {
                bond_dir_counts[second2.index()] += 1;
                atom_dir_counts[atom2.index()] += 1;
            }
        } else if let Some(second2) = second2 {
            set_direction_from_neighboring_bond_for_writer(
                topology,
                second2,
                second2_flipped,
                first2,
                first2_flipped,
            );
            bond_dir_counts[second2.index()] += 1;
            bond_dir_counts[first2.index()] += 1;
            atom_dir_counts[atom2.index()] += 2;
            atom2_controlling_bond = second2;
        }
    }

    if set_from_atom1 {
        let controlling_flipped = if atom1_controlling_bond == first1 {
            first1_flipped
        } else {
            second1_flipped
        };
        let dir = get_reference_direction_for_writer(
            topology,
            dbl_bond,
            atom1,
            atom2,
            atom1_controlling_bond,
            controlling_flipped,
            first2,
            first2_flipped,
        );
        set_bond_direction(topology, first2, dir);
        bond_dir_counts[first2.index()] += 1;
        atom_dir_counts[atom2.index()] += 1;
    } else {
        let controlling_flipped = if atom2_controlling_bond == first2 {
            first2_flipped
        } else {
            second2_flipped
        };
        let dir = get_reference_direction_for_writer(
            topology,
            dbl_bond,
            atom2,
            atom1,
            atom2_controlling_bond,
            controlling_flipped,
            first1,
            first1_flipped,
        );
        set_bond_direction(topology, first1, dir);
        bond_dir_counts[first1.index()] += 1;
        atom_dir_counts[atom1.index()] += 1;
    }

    if incident_bonds(topology, atom1).len() == 3
        && let Some(second1) = second1
        && bond_dir_counts[second1.index()] == 0
    {
        set_direction_from_neighboring_bond_for_writer(
            topology,
            first1,
            first1_flipped,
            second1,
            second1_flipped,
        );
        bond_dir_counts[second1.index()] += 1;
        atom_dir_counts[atom1.index()] += 1;
    }
    if incident_bonds(topology, atom2).len() == 3
        && let Some(second2) = second2
        && bond_dir_counts[second2.index()] == 0
    {
        set_direction_from_neighboring_bond_for_writer(
            topology,
            first2,
            first2_flipped,
            second2,
            second2_flipped,
        );
        bond_dir_counts[second2.index()] += 1;
        atom_dir_counts[atom2.index()] += 1;
    }
}

pub(super) fn find_double_bond_neighbor_bonds_for_writer(
    topology: &TopologyBlock,
    dbl_bond: BondId,
    atom: AtomId,
    bond_visit_orders: &[usize],
    bond_dir_counts: &[i8],
) -> (Option<BondId>, Option<BondId>, bool) {
    let mut first = None;
    let mut second = None;
    let mut first_visit = usize::MAX;
    let mut dir_set = false;
    for bond in incident_bonds(topology, atom) {
        if bond == dbl_bond
            || !can_set_double_bond_stereo_for_writer(topology.bonds[bond.index()].order())
        {
            continue;
        }
        if bond_dir_counts[bond.index()] > 0 {
            dir_set = true;
        }
        if first.is_none() || bond_visit_orders[bond.index()] < first_visit {
            if first.is_some() {
                second = first;
            }
            first = Some(bond);
            first_visit = bond_visit_orders[bond.index()];
        } else {
            second = Some(bond);
        }
    }
    (first, second, dir_set)
}

pub(super) fn account_same_side_dirs_for_writer(
    topology: &mut TopologyBlock,
    atom: AtomId,
    first: BondId,
    first_flipped: bool,
    second: Option<BondId>,
    second_flipped: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> bool {
    let mut consistent = true;
    if let Some(second) = second {
        if bond_dir_counts[first.index()] == 0 {
            set_direction_from_neighboring_bond_for_writer(
                topology,
                second,
                second_flipped,
                first,
                first_flipped,
            );
        } else if bond_dir_counts[second.index()] == 0 {
            set_direction_from_neighboring_bond_for_writer(
                topology,
                first,
                first_flipped,
                second,
                second_flipped,
            );
        } else {
            consistent = same_side_dirs_are_compatible_for_writer(
                topology,
                first,
                first_flipped,
                second,
                second_flipped,
            );
        }
        bond_dir_counts[second.index()] += 1;
        atom_dir_counts[atom.index()] += 1;
    }
    bond_dir_counts[first.index()] += 1;
    atom_dir_counts[atom.index()] += 1;
    consistent
}

#[allow(clippy::too_many_arguments)]
pub(super) fn handle_dir_conflicts_across_double_bond_for_writer(
    topology: &mut TopologyBlock,
    dbl_bond: BondId,
    atom1: AtomId,
    atom1_consistent: bool,
    first1: BondId,
    first1_flipped: bool,
    second1: Option<BondId>,
    second1_flipped: bool,
    atom2: AtomId,
    atom2_consistent: bool,
    first2: BondId,
    first2_flipped: bool,
    second2: Option<BondId>,
    second2_flipped: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> bool {
    // BEGIN RDKIT CPP FUNCTION handleDirConflictsAcrossDoubleBond
    // RDKit✔️✔️: if (atom1DirsAreConsistent && atom2DirsAreConsistent) { ... }
    // RDKit✔️✔️: else if (!atom2DirsAreConsistent && atom1DirsAreConsistent) {
    // RDKit✔️✔️:   return fixConflictAcrossDoubleBond(...);
    // RDKit✔️✔️: } else if (!atom1DirsAreConsistent && atom2DirsAreConsistent) {
    // RDKit✔️✔️:   return fixConflictAcrossDoubleBond(...);
    // RDKit✔️✔️: } else { ... try all four direction combinations ... }
    // END RDKIT CPP FUNCTION handleDirConflictsAcrossDoubleBond
    if atom1_consistent && atom2_consistent {
        return get_reference_direction_for_writer(
            topology,
            dbl_bond,
            atom1,
            atom2,
            first1,
            first1_flipped,
            first2,
            first2_flipped,
        ) == topology.bonds[first2.index()].direction();
    }
    if !atom2_consistent && atom1_consistent {
        if let Some(second2) = second2 {
            return fix_conflict_across_double_bond_for_writer(
                topology,
                dbl_bond,
                atom2,
                first2,
                first2_flipped,
                second2,
                second2_flipped,
                atom1,
                first1,
                first1_flipped,
                bond_dir_counts,
                atom_dir_counts,
            );
        }
        return false;
    }
    if !atom1_consistent && atom2_consistent {
        if let Some(second1) = second1 {
            return fix_conflict_across_double_bond_for_writer(
                topology,
                dbl_bond,
                atom1,
                first1,
                first1_flipped,
                second1,
                second1_flipped,
                atom2,
                first2,
                first2_flipped,
                bond_dir_counts,
                atom_dir_counts,
            );
        }
        return false;
    }
    let (Some(second1), Some(second2)) = (second1, second2) else {
        return false;
    };
    for (atom1_bond, atom1_flipped, atom1_other) in [
        (first1, first1_flipped, second1),
        (second1, second1_flipped, first1),
    ] {
        for (atom2_bond, atom2_flipped, atom2_other) in [
            (first2, first2_flipped, second2),
            (second2, second2_flipped, first2),
        ] {
            let expected = get_reference_direction_for_writer(
                topology,
                dbl_bond,
                atom1,
                atom2,
                atom1_bond,
                atom1_flipped,
                atom2_bond,
                atom2_flipped,
            );
            if expected != topology.bonds[atom2_bond.index()].direction() {
                continue;
            }
            let Some(atom1_other_idx) =
                bond_other_atom(&topology.bonds[atom1_other.index()], atom1)
            else {
                continue;
            };
            let Some(atom2_other_idx) =
                bond_other_atom(&topology.bonds[atom2_other.index()], atom2)
            else {
                continue;
            };
            if atom1_other_idx == atom2_other_idx
                || atom_dir_counts[atom1_other_idx.index()] != 2
                || atom_dir_counts[atom2_other_idx.index()] != 2
            {
                continue;
            }
            bond_dir_counts[atom1_other.index()] = 0;
            atom_dir_counts[atom1.index()] -= 1;
            atom_dir_counts[atom1_other_idx.index()] -= 1;
            bond_dir_counts[atom2_other.index()] = 0;
            atom_dir_counts[atom2.index()] -= 1;
            atom_dir_counts[atom2_other_idx.index()] -= 1;
            return true;
        }
    }
    false
}

#[allow(clippy::too_many_arguments)]
pub(super) fn fix_conflict_across_double_bond_for_writer(
    topology: &TopologyBlock,
    dbl_bond: BondId,
    atom: AtomId,
    first: BondId,
    first_flipped: bool,
    second: BondId,
    second_flipped: bool,
    ref_atom: AtomId,
    ref_bond: BondId,
    ref_flipped: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> bool {
    for (bond, flipped, other_bond) in [
        (first, first_flipped, second),
        (second, second_flipped, first),
    ] {
        let Some(other_idx) = bond_other_atom(&topology.bonds[other_bond.index()], atom) else {
            continue;
        };
        if atom_dir_counts[other_idx.index()] != 2 {
            continue;
        }
        let expected = get_reference_direction_for_writer(
            topology,
            dbl_bond,
            ref_atom,
            atom,
            ref_bond,
            ref_flipped,
            bond,
            flipped,
        );
        if expected == topology.bonds[bond.index()].direction() {
            bond_dir_counts[other_bond.index()] = 0;
            atom_dir_counts[atom.index()] -= 1;
            atom_dir_counts[other_idx.index()] -= 1;
            return true;
        }
    }
    false
}

pub(super) fn remove_unwanted_bond_dir_specs_for_writer(
    topology: &mut TopologyBlock,
    stack: &[MolStackElem],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
    bond_visit_orders: &[usize],
) {
    // BEGIN RDKIT CPP FUNCTION Canon::removeUnwantedBondDirSpecs
    // RDKit✔️✔️: for (auto &msI : molStack) {
    // RDKit✔️✔️:   if (msI.type != MOL_STACK_BOND) { continue; }
    // RDKit✔️✔️:   if (msI.obj.bond->getBondType() != Bond::DOUBLE ||
    // RDKit✔️✔️:       msI.obj.bond->getStereo() > Bond::STEREOANY) { continue; }
    // RDKit✔️✔️:   ... collect removable direction bonds on both ends ...
    // RDKit✔️✔️:   if (atomDirCounts[otherAtom->getIdx()] == 2) {
    // RDKit✔️✔️:     bondDirCounts[candidateBond->getIdx()] = 0;
    // RDKit✔️✔️:     candidateBond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:     atomDirCounts[otherAtom->getIdx()] -= 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::removeUnwantedBondDirSpecs
    for item in stack {
        let MolStackElem::Bond { bond: dbl_bond, .. } = *item else {
            continue;
        };
        let bond = &topology.bonds[dbl_bond.index()];
        if bond.order() != BondOrder::Double || is_writer_stereo_double_bond(bond) {
            continue;
        }
        let first = bond.begin();
        let second = bond.end();
        if incident_bonds(topology, first).len() == 1 || incident_bonds(topology, second).len() == 1
        {
            continue;
        }
        let mut candidates = Vec::new();
        for bond in incident_bonds(topology, first) {
            if bond_dir_counts[bond.index()] > 0 {
                candidates.push(bond);
            }
        }
        if candidates.is_empty() {
            continue;
        }
        if atom_dir_counts[first.index()] > 0 {
            candidates.clear();
        }
        let mut second_end_candidates = 0usize;
        for bond in incident_bonds(topology, second) {
            if bond_dir_counts[bond.index()] > 0 {
                candidates.push(bond);
                second_end_candidates += 1;
            }
        }
        if second_end_candidates == 0 || atom_dir_counts[second.index()] > 0 {
            continue;
        }
        candidates.sort_by_key(|bond| bond_visit_orders[bond.index()]);
        for candidate in candidates {
            let other = if bond_other_atom(&topology.bonds[candidate.index()], first).is_some() {
                bond_other_atom(&topology.bonds[candidate.index()], first)
            } else {
                bond_other_atom(&topology.bonds[candidate.index()], second)
            };
            let Some(other) = other else {
                continue;
            };
            if atom_dir_counts[other.index()] == 2 {
                bond_dir_counts[candidate.index()] = 0;
                topology.bonds[candidate.index()].set_direction(BondDirection::None);
                atom_dir_counts[other.index()] -= 1;
                break;
            }
        }
    }
}

pub(super) fn remove_redundant_bond_dir_specs_for_writer(
    topology: &mut TopologyBlock,
    stack: &[MolStackElem],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    // BEGIN RDKIT CPP FUNCTION Canon::removeRedundantBondDirSpecs / clearBondDirs
    // RDKit✔️✔️: auto clearBondDirsFromAtom = [&mol, &bondDirCounts, &atomDirCounts](
    // RDKit✔️✔️:                                  Bond *tBond, const Atom *atom) {
    // RDKit✔️✔️:   if (atomDirCounts[atom->getIdx()] < 2) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (bond != tBond && bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:         bond->getStereo() > Bond::STEREOANY) {
    // RDKit✔️✔️:       clearBondDirs(mol, tBond, atom, bondDirCounts, atomDirCounts);
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // RDKit✔️✔️: if (canHaveDirection(*tBond) && bondDirCounts[tBond->getIdx()]) {
    // RDKit✔️✔️:   clearBondDirsFromAtom(tBond, canonBeginAtom);
    // RDKit✔️✔️:   clearBondDirsFromAtom(tBond, canonEndAtom);
    // RDKit✔️✔️: } else if (tBond->getBondDir() != Bond::NONE) {
    // RDKit✔️✔️:   tBond->setBondDir(Bond::NONE);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::removeRedundantBondDirSpecs / clearBondDirs
    for item in stack {
        let MolStackElem::Bond {
            bond: t_bond,
            atom_to_left,
        } = *item
        else {
            continue;
        };
        let t_bond_ref = &topology.bonds[t_bond.index()];
        if can_have_direction_for_writer(t_bond_ref.order()) && bond_dir_counts[t_bond.index()] > 0
        {
            let canon_begin_atom = AtomId::new(atom_to_left);
            let Some(canon_end_atom) = bond_other_atom(t_bond_ref, canon_begin_atom) else {
                continue;
            };

            clear_redundant_bond_dirs_from_atom_for_writer(
                topology,
                t_bond,
                canon_begin_atom,
                bond_dir_counts,
                atom_dir_counts,
            );
            clear_redundant_bond_dirs_from_atom_for_writer(
                topology,
                t_bond,
                canon_end_atom,
                bond_dir_counts,
                atom_dir_counts,
            );
        } else if topology.bonds[t_bond.index()].direction() != BondDirection::None {
            topology.bonds[t_bond.index()].set_direction(BondDirection::None);
        }
    }
}

pub(super) fn clear_redundant_bond_dirs_from_atom_for_writer(
    topology: &mut TopologyBlock,
    ref_bond: BondId,
    atom: AtomId,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    if atom_dir_counts[atom.index()] < 2 {
        return;
    }
    for bond in incident_bonds(topology, atom) {
        if bond != ref_bond
            && topology.bonds[bond.index()].order() == BondOrder::Double
            && matches!(
                topology.bonds[bond.index()].stereo(),
                BondStereo::Z
                    | BondStereo::E
                    | BondStereo::Cis
                    | BondStereo::Trans
                    | BondStereo::AtropCw
                    | BondStereo::AtropCcw
            )
        {
            clear_bond_dirs_from_atom_for_writer(
                topology,
                ref_bond,
                atom,
                bond_dir_counts,
                atom_dir_counts,
            );
            return;
        }
    }
}

pub(super) fn clear_bond_dirs_from_atom_for_writer(
    topology: &mut TopologyBlock,
    ref_bond: BondId,
    atom: AtomId,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    // BEGIN RDKIT CPP FUNCTION clearBondDirs
    // RDKit✔️✔️: void clearBondDirs(ROMol &mol, Bond *refBond, const Atom *fromAtom,
    // RDKit✔️✔️:                    std::vector<int8_t> &bondDirCounts,
    // RDKit✔️✔️:                    std::vector<int8_t> &atomDirCounts) {
    // RDKit✔️✔️:   auto clearDirection = [&atomDirCounts, &bondDirCounts](Bond *bond,
    // RDKit✔️✔️:                                                          const Atom *fromAtom) {
    // RDKit✔️✔️:     --bondDirCounts[bond->getIdx()];
    // RDKit✔️✔️:     if (!bondDirCounts[bond->getIdx()]) {
    // RDKit✔️✔️:       bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:       --atomDirCounts[fromAtom->getIdx()];
    // RDKit✔️✔️:       if (auto otherAtom = bond->getOtherAtom(fromAtom);
    // RDKit✔️✔️:           atomDirCounts[otherAtom->getIdx()]) {
    // RDKit✔️✔️:         --atomDirCounts[otherAtom->getIdx()];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   for (auto oBond : mol.atomBonds(fromAtom)) {
    // RDKit✔️✔️:     if (oBond != refBond && canHaveDirection(*oBond)) {
    // RDKit✔️✔️:       if ((bondDirCounts[oBond->getIdx()] >=
    // RDKit✔️✔️:            bondDirCounts[refBond->getIdx()]) &&
    // RDKit✔️✔️:           atomDirCounts[oBond->getBeginAtomIdx()] != 1 &&
    // RDKit✔️✔️:           atomDirCounts[oBond->getEndAtomIdx()] != 1) {
    // RDKit✔️✔️:         clearDirection(oBond, fromAtom);
    // RDKit✔️✔️:       } else if (atomDirCounts[refBond->getBeginAtomIdx()] != 1 &&
    // RDKit✔️✔️:                  atomDirCounts[refBond->getEndAtomIdx()] != 1) {
    // RDKit✔️✔️:         clearDirection(refBond, fromAtom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION clearBondDirs
    if atom_dir_counts[atom.index()] < 2 {
        return;
    }
    for other_bond in incident_bonds(topology, atom) {
        if other_bond == ref_bond
            || !can_have_direction_for_writer(topology.bonds[other_bond.index()].order())
        {
            continue;
        }
        let ref_count = bond_dir_counts[ref_bond.index()];
        let other_count = bond_dir_counts[other_bond.index()];
        let other_begin = topology.bonds[other_bond.index()].begin();
        let other_end = topology.bonds[other_bond.index()].end();
        if other_count >= ref_count
            && atom_dir_counts[other_begin.index()] != 1
            && atom_dir_counts[other_end.index()] != 1
        {
            clear_direction_counted_for_writer(
                topology,
                other_bond,
                atom,
                bond_dir_counts,
                atom_dir_counts,
            );
        } else {
            let ref_begin = topology.bonds[ref_bond.index()].begin();
            let ref_end = topology.bonds[ref_bond.index()].end();
            if atom_dir_counts[ref_begin.index()] != 1 && atom_dir_counts[ref_end.index()] != 1 {
                clear_direction_counted_for_writer(
                    topology,
                    ref_bond,
                    atom,
                    bond_dir_counts,
                    atom_dir_counts,
                );
            }
        }
        break;
    }
}

pub(super) fn clear_direction_counted_for_writer(
    topology: &mut TopologyBlock,
    bond: BondId,
    from_atom: AtomId,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    bond_dir_counts[bond.index()] -= 1;
    if bond_dir_counts[bond.index()] == 0 {
        topology.bonds[bond.index()].set_direction(BondDirection::None);
        atom_dir_counts[from_atom.index()] -= 1;
        if let Some(other) = bond_other_atom(&topology.bonds[bond.index()], from_atom)
            && atom_dir_counts[other.index()] > 0
        {
            atom_dir_counts[other.index()] -= 1;
        }
    }
}

pub(super) fn neighboring_stereo_double_bonds_for_writer(
    topology: &TopologyBlock,
    dbl_bond: BondId,
) -> Vec<BondId> {
    let mut out = Vec::new();
    let bond = &topology.bonds[dbl_bond.index()];
    for atom in [bond.begin(), bond.end()] {
        for nbr_bond in incident_bonds(topology, atom) {
            if !can_have_direction_for_writer(topology.bonds[nbr_bond.index()].order()) {
                continue;
            }
            let Some(other_atom) = bond_other_atom(&topology.bonds[nbr_bond.index()], atom) else {
                continue;
            };
            for other_bond in incident_bonds(topology, other_atom) {
                if other_bond != nbr_bond
                    && is_writer_stereo_double_bond(&topology.bonds[other_bond.index()])
                    && !out.contains(&other_bond)
                {
                    out.push(other_bond);
                }
            }
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
pub(super) fn get_reference_direction_for_writer(
    topology: &TopologyBlock,
    dbl_bond: BondId,
    ref_atom: AtomId,
    target_atom: AtomId,
    ref_controlling_bond: BondId,
    ref_is_flipped: bool,
    target_bond: BondId,
    target_is_flipped: bool,
) -> BondDirection {
    // BEGIN RDKIT CPP FUNCTION getReferenceDirection
    // RDKit❗✔️: if (dblBond.getStereo() == Bond::STEREOE ||
    // RDKit❗✔️:     dblBond.getStereo() == Bond::STEREOTRANS) {
    // RDKit❗✔️:   dir = refControllingBond.getBondDir();
    // RDKit❗✔️: } else if (dblBond.getStereo() == Bond::STEREOZ ||
    // RDKit❗✔️:            dblBond.getStereo() == Bond::STEREOCIS) {
    // RDKit❗✔️:   dir = flipStereoBondDir(refControllingBond.getBondDir());
    // RDKit❗✔️: }
    // RDKit❗✔️: if (refAtom.getDegree() == 3 && ... stereoAtoms ... ) { dir = flipStereoBondDir(dir); }
    // RDKit❗✔️: if (targetAtom.getDegree() == 3 && ... stereoAtoms ... ) { dir = flipStereoBondDir(dir); }
    // RDKit❗✔️: if (refIsFlipped != targetIsFlipped) { dir = flipStereoBondDir(dir); }
    // END RDKIT CPP FUNCTION getReferenceDirection
    let dbl = &topology.bonds[dbl_bond.index()];
    let mut dir = match dbl.stereo() {
        BondStereo::E | BondStereo::Trans => {
            topology.bonds[ref_controlling_bond.index()].direction()
        }
        BondStereo::Z | BondStereo::Cis => flip_stereo_bond_dir_for_writer(
            topology.bonds[ref_controlling_bond.index()].direction(),
        ),
        _ => BondDirection::None,
    };
    if let Some(stereo_atoms) = dbl.stereo_atoms() {
        if incident_bonds(topology, ref_atom).len() == 3
            && !stereo_atoms.contains(
                &bond_other_atom(&topology.bonds[ref_controlling_bond.index()], ref_atom)
                    .unwrap_or(ref_atom),
            )
        {
            dir = flip_stereo_bond_dir_for_writer(dir);
        }
        if incident_bonds(topology, target_atom).len() == 3
            && !stereo_atoms.contains(
                &bond_other_atom(&topology.bonds[target_bond.index()], target_atom)
                    .unwrap_or(target_atom),
            )
        {
            dir = flip_stereo_bond_dir_for_writer(dir);
        }
    }
    if ref_is_flipped != target_is_flipped {
        dir = flip_stereo_bond_dir_for_writer(dir);
    }
    dir
}

pub(super) fn set_direction_from_neighboring_bond_for_writer(
    topology: &mut TopologyBlock,
    source: BondId,
    source_flipped: bool,
    target: BondId,
    target_flipped: bool,
) {
    // BEGIN RDKIT CPP FUNCTION setDirectionFromNeighboringBond
    // RDKit❗✔️: auto dir = sourceBond.getBondDir();
    // RDKit❗✔️: if (isSourceBondFlipped == isTargetBondFlipped) {
    // RDKit❗✔️:   dir = flipStereoBondDir(dir);
    // RDKit❗✔️: }
    // RDKit❗✔️: targetBond.setBondDir(dir);
    // END RDKIT CPP FUNCTION setDirectionFromNeighboringBond
    let mut dir = topology.bonds[source.index()].direction();
    if source_flipped == target_flipped {
        dir = flip_stereo_bond_dir_for_writer(dir);
    }
    set_bond_direction(topology, target, dir);
}

pub(super) fn same_side_dirs_are_compatible_for_writer(
    topology: &TopologyBlock,
    first: BondId,
    first_flipped: bool,
    second: BondId,
    second_flipped: bool,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION sameSideDirsAreCompatible
    // RDKit❗✔️: auto dirsShouldMatch = isFirstBondFlipped != isSecondBondFlipped;
    // RDKit❗✔️: auto dirsMatch = firstBond.getBondDir() == secondBond.getBondDir();
    // RDKit❗✔️: return dirsMatch == dirsShouldMatch;
    // END RDKIT CPP FUNCTION sameSideDirsAreCompatible
    let dirs_should_match = first_flipped != second_flipped;
    let dirs_match =
        topology.bonds[first.index()].direction() == topology.bonds[second.index()].direction();
    dirs_match == dirs_should_match
}

pub(super) fn writer_bond_is_flipped_from_atom1(
    topology: &TopologyBlock,
    atom: AtomId,
    bond: BondId,
    atom_visit_orders: &[usize],
    traversal_ring_closure_bonds: &[bool],
) -> bool {
    let anchor = bond_other_atom(&topology.bonds[bond.index()], atom).unwrap_or(atom);
    (atom_visit_orders[atom.index()] < atom_visit_orders[anchor.index()])
        != traversal_ring_closure_bonds[bond.index()]
}

pub(super) fn writer_bond_is_flipped_from_atom2(
    topology: &TopologyBlock,
    atom: AtomId,
    bond: BondId,
    atom_visit_orders: &[usize],
    traversal_ring_closure_bonds: &[bool],
) -> bool {
    let anchor = bond_other_atom(&topology.bonds[bond.index()], atom).unwrap_or(atom);
    (atom_visit_orders[anchor.index()] < atom_visit_orders[atom.index()])
        != traversal_ring_closure_bonds[bond.index()]
}

pub(super) fn flip_stereo_bond_dir_for_writer(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        other => other,
    }
}

pub(super) fn is_writer_stereo_double_bond(bond: &Bond) -> bool {
    bond.order() == BondOrder::Double
        && matches!(
            bond.stereo(),
            BondStereo::E | BondStereo::Z | BondStereo::Cis | BondStereo::Trans
        )
        && bond.stereo_atoms().is_some()
}

pub(super) fn can_have_direction_for_writer(order: BondOrder) -> bool {
    matches!(order, BondOrder::Single | BondOrder::Aromatic)
}

pub(super) fn can_set_double_bond_stereo_for_writer(order: BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Single
            | BondOrder::Aromatic
            | BondOrder::Dative
            | BondOrder::DativeOne
            | BondOrder::DativeLeft
            | BondOrder::DativeRight
    )
}

pub(super) fn set_bond_direction(
    topology: &mut TopologyBlock,
    bond: BondId,
    direction: BondDirection,
) {
    topology.bonds[bond.index()].set_direction(direction);
}

fn incident_bonds(topology: &TopologyBlock, atom: AtomId) -> Vec<BondId> {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .map(|neighbor| neighbor.bond)
        .collect()
}

fn bond_other_atom(bond: &Bond, atom: AtomId) -> Option<AtomId> {
    if bond.begin() == atom {
        Some(bond.end())
    } else if bond.end() == atom {
        Some(bond.begin())
    } else {
        None
    }
}
