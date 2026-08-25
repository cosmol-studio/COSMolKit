//! RDKit generic-group matching for ordinary substructure matches.

use std::collections::{BTreeSet, VecDeque};

use crate::search::query::{
    CompositeQueryType, make_atom_explicit_degree_query, make_atom_null_query,
    query_atom_expand_query, replace_atom_with_query_atom,
};
use crate::{
    Atom, AtomQueryPredicate, Bond, Molecule, QueryNode, SubstanceGroup, SubstanceGroupId,
    SubstanceGroupKind,
};

type AtomMatcher<'a> = dyn Fn(&Atom) -> bool + 'a;
type BondMatcher<'a> = dyn Fn(&Bond) -> bool + 'a;

fn is_hydrogen(molecule: &Molecule, atom_index: usize, mut ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool IsHydrogen(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                 boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (atom.getAtomicNum() == 1 && mol.getAtomDegree(&atom) == 1) {
    // RDKit✔️✔️:     ignore.set(atom.getIdx());  // just an H atom
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // Complexity review: both implementations perform two O(1) indexed
    // reads and, on success, one O(1) bit write. The by-value ignore bitmap
    // is copied by the caller boundary in both implementations.
    let Some(atom) = molecule.atoms().get(atom_index) else {
        return false;
    };
    if atom.atomic_number() == 1
        && molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom_index)
            .len()
            == 1
    {
        if let Some(bit) = ignore.get_mut(atom_index) {
            *bit = true;
        }
        return true;
    }
    false
}

#[allow(clippy::too_many_arguments)]
fn all_atoms_match(
    molecule: &Molecule,
    atom_index: usize,
    mut ignore: Vec<bool>,
    matcher: Option<&AtomMatcher<'_>>,
    bond_matcher: Option<&BondMatcher<'_>>,
    at_least_one_atom: Option<&AtomMatcher<'_>>,
    at_least_one_bond: Option<&BondMatcher<'_>>,
) -> bool {
    // RDKit✔️✔️: bool AllAtomsMatch(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                    boost::dynamic_bitset<> ignore, AtomMatcherFunc matcher,
    // RDKit✔️✔️:                    BondMatcherFunc bondMatcher = nullptr,
    // RDKit✔️✔️:                    AtomMatcherFunc atLeastOneAtom = nullptr,
    // RDKit✔️✔️:                    BondMatcherFunc atLeastOneBond = nullptr) {
    // RDKit✔️✔️:   PRECONDITION(&atom.getOwningMol() == &mol, "atom not owned by molecule");
    // RDKit✔️✔️:   if (matcher && !matcher(atom)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool atomAtLeast = atLeastOneAtom == nullptr;
    // RDKit✔️✔️:   bool bondAtLeast = atLeastOneBond == nullptr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::deque<const Atom *> nbrs;
    // RDKit✔️✔️:   nbrs.push_back(&atom);
    // RDKit✔️✔️:   while (!nbrs.empty()) {
    // RDKit✔️✔️:     const auto atm = nbrs.front();
    // RDKit✔️✔️:     if (!atomAtLeast && atLeastOneAtom(*atm)) {
    // RDKit✔️✔️:       atomAtLeast = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     nbrs.pop_front();
    // RDKit✔️✔️:     ignore.set(atm->getIdx());
    // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atm)) {
    // RDKit✔️✔️:       if (ignore[nbr->getIdx()]) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (matcher && !matcher(*nbr)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (bondMatcher || !bondAtLeast) {
    // RDKit✔️✔️:         const auto bnd = mol.getBondBetweenAtoms(atm->getIdx(), nbr->getIdx());
    // RDKit✔️✔️:         if (bondMatcher && !(bondMatcher)(*bnd)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (!bondAtLeast && atLeastOneBond(*bnd)) {
    // RDKit✔️✔️:           bondAtLeast = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       nbrs.push_back(nbr);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return atomAtLeast && bondAtLeast;
    // RDKit✔️✔️: }
    // Complexity review: Rust uses the molecule's CSR adjacency and the bond
    // id stored in each neighbor row, so it retains RDKit's O(V + E) traversal,
    // O(V) copied bitmap, and O(V) worst-case queue. It deliberately marks on
    // dequeue, matching RDKit's callback and duplicate-enqueue behavior.
    let Some(atom) = molecule.atoms().get(atom_index) else {
        return false;
    };
    if matcher.is_some_and(|matches| !matches(atom)) {
        return false;
    }

    let mut atom_at_least = at_least_one_atom.is_none();
    let mut bond_at_least = at_least_one_bond.is_none();
    let mut neighbors = VecDeque::from([atom_index]);
    while let Some(current_index) = neighbors.front().copied() {
        let current = &molecule.atoms()[current_index];
        if !atom_at_least && at_least_one_atom.is_some_and(|matches| matches(current)) {
            atom_at_least = true;
        }
        neighbors.pop_front();
        ignore[current_index] = true;

        for neighbor in molecule
            .topology_block()
            .adjacency
            .neighbors_of(current_index)
        {
            if ignore[neighbor.atom_index] {
                continue;
            }
            let neighbor_atom = &molecule.atoms()[neighbor.atom_index];
            if matcher.is_some_and(|matches| !matches(neighbor_atom)) {
                return false;
            }
            if bond_matcher.is_some() || !bond_at_least {
                let bond = &molecule.bonds()[neighbor.bond.index()];
                if bond_matcher.is_some_and(|matches| !matches(bond)) {
                    return false;
                }
                if !bond_at_least && at_least_one_bond.is_some_and(|matches| matches(bond)) {
                    bond_at_least = true;
                }
            }
            neighbors.push_back(neighbor.atom_index);
        }
    }
    atom_at_least && bond_at_least
}

fn group_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool GroupAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                       boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = nullptr;
    // RDKit✔️✔️:   auto bondMatcher = nullptr;
    // RDKit✔️✔️:   auto atLeastMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() != 1;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher,
    // RDKit✔️✔️:                        atLeastMatcher);
    // RDKit✔️✔️: }
    // Complexity review: both paths perform one O(V + E) fast ring pass when
    // no suitable ring state exists, then one O(V + E) component traversal.
    // Rust discards the immutable ring result because this matcher does not
    // inspect it; no second traversal implementation is introduced.
    let _ring_info = crate::fast_find_rings(molecule).ok();
    let heavy_atom = |atom: &Atom| atom.atomic_number() != 1;
    all_atoms_match(
        molecule,
        atom_index,
        ignore,
        None,
        None,
        Some(&heavy_atom),
        None,
    )
}

fn group_h_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool GroupHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                        boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return GroupAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: both versions add one O(1) hydrogen check before the
    // same canonical group matcher and retain its O(V + E) worst case.
    is_hydrogen(molecule, atom_index, ignore.clone())
        || group_atom_matcher(molecule, atom_index, ignore)
}

fn group_star_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool GroupStarAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                           boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = nullptr;
    // RDKit✔️✔️:   auto bondMatcher = nullptr;
    // RDKit✔️✔️:   auto atLeastMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() != 1;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto atLeastBondMatcher = [](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return queryIsBondInRing(&bnd);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher,
    // RDKit✔️✔️:                        atLeastMatcher, atLeastBondMatcher);
    // RDKit✔️✔️: }
    // Complexity review: both implementations build/read fast ring membership
    // in O(V + E), then perform the one O(V + E) canonical component traversal
    // with O(1) indexed bond-ring membership checks.
    let Ok(ring_info) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let heavy_atom = |atom: &Atom| atom.atomic_number() != 1;
    let ring_bond = |bond: &Bond| ring_info.num_bond_rings(bond.id()) != 0;
    all_atoms_match(
        molecule,
        atom_index,
        ignore,
        None,
        None,
        Some(&heavy_atom),
        Some(&ring_bond),
    )
}

fn group_star_h_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool GroupStarHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                            boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return GroupStarAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: both versions add one O(1) hydrogen check before the
    // same canonical ring-containing group matcher.
    is_hydrogen(molecule, atom_index, ignore.clone())
        || group_star_atom_matcher(molecule, atom_index, ignore)
}

fn alkyl_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkylAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                       boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return !at.getIsAromatic() &&
    // RDKit✔️✔️:            (at.getAtomicNum() == 6 || at.getAtomicNum() == 1);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto atLeastMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() == 6;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto bondMatcher = [](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return bnd.getBondType() == Bond::BondType::SINGLE &&
    // RDKit✔️✔️:            !bnd.getIsAromatic() && !queryIsBondInRing(&bnd);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher,
    // RDKit✔️✔️:                        atLeastMatcher);
    // RDKit✔️✔️: }
    // Complexity review: both versions perceive fast rings and traverse the
    // exposed component once in O(V + E), using O(1) predicates.
    let Ok(ring_info) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let atoms = |atom: &Atom| !atom.is_aromatic() && matches!(atom.atomic_number(), 6 | 1);
    let carbon = |atom: &Atom| atom.atomic_number() == 6;
    let bonds = |bond: &Bond| {
        bond.order() == crate::BondOrder::Single
            && !bond.is_aromatic()
            && ring_info.num_bond_rings(bond.id()) == 0
    };
    all_atoms_match(
        molecule,
        atom_index,
        ignore,
        Some(&atoms),
        Some(&bonds),
        Some(&carbon),
        None,
    )
}

fn alkyl_h_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                        boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return !at.getIsAromatic() &&
    // RDKit✔️✔️:            (at.getAtomicNum() == 6 || at.getAtomicNum() == 1);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto bondMatcher = [](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return bnd.getBondType() == Bond::BondType::SINGLE &&
    // RDKit✔️✔️:            !bnd.getIsAromatic() && !queryIsBondInRing(&bnd);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher);
    // RDKit✔️✔️: }
    // Complexity review: both versions use O(V + E) ring perception and the
    // one O(V + E) component traversal, without the carbon witness.
    let Ok(ring_info) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let atoms = |atom: &Atom| !atom.is_aromatic() && matches!(atom.atomic_number(), 6 | 1);
    let bonds = |bond: &Bond| {
        bond.order() == crate::BondOrder::Single
            && !bond.is_aromatic()
            && ring_info.num_bond_rings(bond.id()) == 0
    };
    all_atoms_match(
        molecule,
        atom_index,
        ignore,
        Some(&atoms),
        Some(&bonds),
        None,
        None,
    )
}

fn unsat_alk_x_atom_matcher(
    molecule: &Molecule,
    atom_index: usize,
    ignore: Vec<bool>,
    extra_bond_type: crate::BondOrder,
) -> bool {
    // RDKit✔️✔️: bool UnsatAlkXAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                           boost::dynamic_bitset<> ignore,
    // RDKit✔️✔️:                           Bond::BondType extraBondType) {
    // RDKit✔️✔️:   // nominally requires at least two Cs, but since it can only
    // RDKit✔️✔️:   // contain Cs and Hs and since a multiple bond is required, that condition is
    // RDKit✔️✔️:   // redundant
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return !at.getIsAromatic() &&
    // RDKit✔️✔️:            (at.getAtomicNum() == 6 || at.getAtomicNum() == 1);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto bondMatcher = [extraBondType](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return (bnd.getBondType() == Bond::BondType::SINGLE ||
    // RDKit✔️✔️:             bnd.getBondType() == extraBondType) &&
    // RDKit✔️✔️:            !bnd.getIsAromatic() && !queryIsBondInRing(&bnd);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto atLeastMatcher = [extraBondType](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return bnd.getBondType() == extraBondType;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   AtomMatcherFunc atomAtLeast = nullptr;
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher, atomAtLeast,
    // RDKit✔️✔️:                        atLeastMatcher);
    // RDKit✔️✔️: }
    // Complexity review: both implementations do O(V + E) ring perception
    // and one O(V + E) traversal with constant-time bond predicates.
    let Ok(ring_info) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let atoms = |atom: &Atom| !atom.is_aromatic() && matches!(atom.atomic_number(), 6 | 1);
    let bonds = |bond: &Bond| {
        matches!(bond.order(), crate::BondOrder::Single) || bond.order() == extra_bond_type
    };
    let acyclic_bonds = |bond: &Bond| {
        bonds(bond) && !bond.is_aromatic() && ring_info.num_bond_rings(bond.id()) == 0
    };
    let required_bond = |bond: &Bond| bond.order() == extra_bond_type;
    all_atoms_match(
        molecule,
        atom_index,
        ignore,
        Some(&atoms),
        Some(&acyclic_bonds),
        None,
        Some(&required_bond),
    )
}

fn alkenyl_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkenylAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                         boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   return UnsatAlkXAtomMatcher(mol, atom, ignore, Bond::BondType::DOUBLE);
    // RDKit✔️✔️: }
    // Complexity review: both are O(1) wrappers around the O(V + E) core.
    unsat_alk_x_atom_matcher(molecule, atom_index, ignore, crate::BondOrder::Double)
}

fn alkenyl_h_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkenylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                          boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return UnsatAlkXAtomMatcher(mol, atom, ignore, Bond::BondType::DOUBLE);
    // RDKit✔️✔️: }
    // Complexity review: both add O(1) before the O(V + E) core.
    is_hydrogen(molecule, atom_index, ignore.clone())
        || unsat_alk_x_atom_matcher(molecule, atom_index, ignore, crate::BondOrder::Double)
}

fn alkynyl_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkynylAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                         boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   return UnsatAlkXAtomMatcher(mol, atom, ignore, Bond::BondType::TRIPLE);
    // RDKit✔️✔️: }
    // Complexity review: both are O(1) wrappers around the O(V + E) core.
    unsat_alk_x_atom_matcher(molecule, atom_index, ignore, crate::BondOrder::Triple)
}

fn alkynyl_h_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkynylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                          boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return UnsatAlkXAtomMatcher(mol, atom, ignore, Bond::BondType::TRIPLE);
    // RDKit✔️✔️: }
    // Complexity review: both add O(1) before the O(V + E) core.
    is_hydrogen(molecule, atom_index, ignore.clone())
        || unsat_alk_x_atom_matcher(molecule, atom_index, ignore, crate::BondOrder::Triple)
}

fn acyclic_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AcyclicAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                         boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) == 0;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto atLeastMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() != 1;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher, nullptr, atLeastMatcher);
    // RDKit✔️✔️: }
    // Complexity review: both perceive rings in O(V + E), then perform the
    // canonical O(V + E) traversal with O(1) indexed ring membership.
    let Ok(rings) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let acyclic = |atom: &Atom| rings.num_atom_rings(atom.id()) == 0;
    let heavy = |atom: &Atom| atom.atomic_number() != 1;
    all_atoms_match(
        molecule,
        atom_index,
        ignore,
        Some(&acyclic),
        None,
        Some(&heavy),
        None,
    )
}

fn acyclic_h_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AcyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                          boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) == 0;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher);
    // RDKit✔️✔️: }
    // Complexity review: equivalent O(V + E) ring perception and traversal.
    let Ok(rings) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let acyclic = |atom: &Atom| rings.num_atom_rings(atom.id()) == 0;
    all_atoms_match(
        molecule,
        atom_index,
        ignore,
        Some(&acyclic),
        None,
        None,
        None,
    )
}

fn carboacyclic_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarboacyclicAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                              boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return (at.getAtomicNum() == 6 || at.getAtomicNum() == 1) &&
    // RDKit✔️✔️:            at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) == 0;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto atLeastMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() == 6;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher, nullptr, atLeastMatcher);
    // RDKit✔️✔️: }
    // Complexity review: equivalent O(V + E) ring perception and traversal.
    let Ok(rings) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let atoms =
        |atom: &Atom| matches!(atom.atomic_number(), 6 | 1) && rings.num_atom_rings(atom.id()) == 0;
    let carbon = |atom: &Atom| atom.atomic_number() == 6;
    all_atoms_match(
        molecule,
        atom_index,
        ignore,
        Some(&atoms),
        None,
        Some(&carbon),
        None,
    )
}

fn carboacyclic_h_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarboacyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                               boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return (at.getAtomicNum() == 6 || at.getAtomicNum() == 1) &&
    // RDKit✔️✔️:            at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) == 0;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher);
    // RDKit✔️✔️: }
    // Complexity review: equivalent O(V + E) ring perception and traversal.
    let Ok(rings) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let atoms =
        |atom: &Atom| matches!(atom.atomic_number(), 6 | 1) && rings.num_atom_rings(atom.id()) == 0;
    all_atoms_match(molecule, atom_index, ignore, Some(&atoms), None, None, None)
}

fn heteroacyclic_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool HeteroacyclicAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                               boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) == 0;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto atLeastOne = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() != 6 && at.getAtomicNum() != 1;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   BondMatcherFunc bondMatcher = nullptr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return AllAtomsMatch(mol, atom, ignore, atomMatcher, bondMatcher, atLeastOne);
    // RDKit✔️✔️: }
    // Complexity review: equivalent O(V + E) ring perception and traversal.
    let Ok(rings) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let acyclic = |atom: &Atom| rings.num_atom_rings(atom.id()) == 0;
    let hetero = |atom: &Atom| !matches!(atom.atomic_number(), 6 | 1);
    all_atoms_match(
        molecule,
        atom_index,
        ignore,
        Some(&acyclic),
        None,
        Some(&hetero),
        None,
    )
}

fn heteroacyclic_h_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool HeteroacyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                                boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return HeteroacyclicAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: both add O(1) before the same O(V + E) matcher.
    is_hydrogen(molecule, atom_index, ignore.clone())
        || heteroacyclic_atom_matcher(molecule, atom_index, ignore)
}

fn alkoxyacyclic_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkoxyacyclicAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                               boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atom.getDegree() != 2 || atom.getAtomicNum() != 8) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const Atom *nnbr = nullptr;
    // RDKit✔️✔️:   for (const auto *nbr : mol.atomNeighbors(&atom)) {
    // RDKit✔️✔️:     if (!ignore[nbr->getIdx()]) {
    // RDKit✔️✔️:       nnbr = nbr;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (nnbr == nullptr) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return !at.getIsAromatic() &&
    // RDKit✔️✔️:            (at.getAtomicNum() == 6 || at.getAtomicNum() == 1);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto atLeastMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() == 6;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto bondMatcher = [](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return bnd.getBondType() == Bond::BondType::SINGLE &&
    // RDKit✔️✔️:            !bnd.getIsAromatic() && !queryIsBondInRing(&bnd);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return AllAtomsMatch(mol, *nnbr, ignore, atomMatcher, bondMatcher,
    // RDKit✔️✔️:                        atLeastMatcher);
    // RDKit✔️✔️: }
    // Complexity review: both perform O(V + E) ring perception, O(degree)
    // attachment selection, and one O(V + E) canonical component traversal.
    let Ok(rings) = crate::fast_find_rings(molecule) else {
        return false;
    };
    let Some(atom) = molecule.atoms().get(atom_index) else {
        return false;
    };
    let neighbors = molecule.topology_block().adjacency.neighbors_of(atom_index);
    if neighbors.len() != 2 || atom.atomic_number() != 8 {
        return false;
    }
    let Some(next) = neighbors
        .iter()
        .find(|neighbor| !ignore[neighbor.atom_index])
    else {
        return false;
    };
    let atoms = |atom: &Atom| !atom.is_aromatic() && matches!(atom.atomic_number(), 6 | 1);
    let carbon = |atom: &Atom| atom.atomic_number() == 6;
    let bonds = |bond: &Bond| {
        bond.order() == crate::BondOrder::Single
            && !bond.is_aromatic()
            && rings.num_bond_rings(bond.id()) == 0
    };
    all_atoms_match(
        molecule,
        next.atom_index,
        ignore,
        Some(&atoms),
        Some(&bonds),
        Some(&carbon),
        None,
    )
}

fn alkoxyacyclic_h_atom_matcher(molecule: &Molecule, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkoxyacyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                                boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return AlkoxyacyclicAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: both add O(1) before the same O(V + E) matcher.
    is_hydrogen(molecule, atom_index, ignore.clone())
        || alkoxyacyclic_atom_matcher(molecule, atom_index, ignore)
}

fn check_atom_ring(
    molecule: &Molecule,
    atom_index: usize,
    ignore: &[bool],
    ring: &[crate::AtomId],
    matcher: Option<&AtomMatcher<'_>>,
    at_least_one: Option<&AtomMatcher<'_>>,
) -> bool {
    // RDKit✔️✔️: bool checkAtomRing(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                    const boost::dynamic_bitset<> &ignore,
    // RDKit✔️✔️:                    const std::vector<int> &ring, AtomMatcherFunc matcher,
    // RDKit✔️✔️:                    AtomMatcherFunc atLeastOne) {
    // RDKit✔️✔️:   bool atLeast = atLeastOne == nullptr;
    // RDKit✔️✔️:   for (auto aidx : ring) {
    // RDKit✔️✔️:     if (aidx != static_cast<int>(atom.getIdx()) &&
    // RDKit✔️✔️:         (ignore[aidx] || (matcher && !matcher(*mol.getAtomWithIdx(aidx))))) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!atLeast && atLeastOne(*mol.getAtomWithIdx(aidx))) {
    // RDKit✔️✔️:       atLeast = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return atLeast;
    // RDKit✔️✔️: }
    // Complexity review: both scan one ring in O(r), with O(1) indexed reads.
    let mut found = at_least_one.is_none();
    for id in ring {
        let index = id.index();
        let atom = &molecule.atoms()[index];
        if index != atom_index && (ignore[index] || matcher.is_some_and(|matches| !matches(atom))) {
            return false;
        }
        if !found && at_least_one.is_some_and(|matches| matches(atom)) {
            found = true;
        }
    }
    found
}

fn check_bond_ring(
    molecule: &Molecule,
    ring: &[crate::BondId],
    matcher: Option<&BondMatcher<'_>>,
    at_least_one: Option<&BondMatcher<'_>>,
) -> bool {
    // RDKit✔️✔️: bool checkBondRing(const ROMol &mol, const std::vector<int> &bring,
    // RDKit✔️✔️:                    BondMatcherFunc matcher, BondMatcherFunc atLeastOne) {
    // RDKit✔️✔️:   bool atLeast = atLeastOne == nullptr;
    // RDKit✔️✔️:   for (auto bidx : bring) {
    // RDKit✔️✔️:     if (matcher && !matcher(*mol.getBondWithIdx(bidx))) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!atLeast && atLeastOne(*mol.getBondWithIdx(bidx))) {
    // RDKit✔️✔️:       atLeast = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return atLeast;
    // RDKit✔️✔️: }
    // Complexity review: both scan one ring in O(r), with O(1) indexed reads.
    let mut found = at_least_one.is_none();
    for id in ring {
        let bond = &molecule.bonds()[id.index()];
        if matcher.is_some_and(|matches| !matches(bond)) {
            return false;
        }
        if !found && at_least_one.is_some_and(|matches| matches(bond)) {
            found = true;
        }
    }
    found
}

#[allow(clippy::too_many_arguments)]
fn fused_ring_match(
    molecule: &Molecule,
    atom_index: usize,
    ignore: Vec<bool>,
    atom_matcher: Option<&AtomMatcher<'_>>,
    bond_matcher: Option<&BondMatcher<'_>>,
    atom_per_ring: Option<&AtomMatcher<'_>>,
    bond_per_ring: Option<&BondMatcher<'_>>,
    at_least_one_atom: Option<&AtomMatcher<'_>>,
) -> bool {
    // RDKit✔️✔️: bool FusedRingMatch(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                     boost::dynamic_bitset<> ignore,
    // RDKit✔️✔️:                     AtomMatcherFunc atomMatcher = nullptr,
    // RDKit✔️✔️:                     BondMatcherFunc bondMatcher = nullptr,
    // RDKit✔️✔️:                     AtomMatcherFunc atLeastOneAtomPerRing = nullptr,
    // RDKit✔️✔️:                     BondMatcherFunc atLeastOneBondPerRing = nullptr,
    // RDKit✔️✔️:                     AtomMatcherFunc atLeastOneAtom = nullptr) {
    // RDKit✔️✔️:   PRECONDITION(&atom.getOwningMol() == &mol, "atom not owned by molecule");
    // RDKit✔️✔️:   if (atomMatcher && !atomMatcher(atom)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     MolOps::findSSSR(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!mol.getRingInfo()->numAtomRings(atom.getIdx())) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::set<int> ringAtoms;
    // RDKit✔️✔️:   for (auto i = 0u; i < mol.getRingInfo()->numRings(); ++i) {
    // RDKit✔️✔️:     const auto &ring = mol.getRingInfo()->atomRings()[i];
    // RDKit✔️✔️:     if (std::find(ring.begin(), ring.end(), atom.getIdx()) != ring.end()) {
    // RDKit✔️✔️:       if (!checkAtomRing(mol, atom, ignore, ring, atomMatcher,
    // RDKit✔️✔️:                          atLeastOneAtomPerRing)) { return false; }
    // RDKit✔️✔️:       if (!checkBondRing(mol, mol.getRingInfo()->bondRings()[i], bondMatcher,
    // RDKit✔️✔️:                          atLeastOneBondPerRing)) { return false; }
    // RDKit✔️✔️:       ringAtoms.insert(ring.begin(), ring.end());
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto i = 0u; i < mol.getRingInfo()->numRings(); ++i) {
    // RDKit✔️✔️:     const auto &ring = mol.getRingInfo()->atomRings()[i];
    // RDKit✔️✔️:     std::set<int> sring(ring.begin(), ring.end());
    // RDKit✔️✔️:     std::vector<int> diff(sring.size());
    // RDKit✔️✔️:     auto dit = std::set_difference(sring.begin(), sring.end(),
    // RDKit✔️✔️:                                    ringAtoms.begin(), ringAtoms.end(), diff.begin());
    // RDKit✔️✔️:     auto numNewAtoms = dit - diff.begin();
    // RDKit✔️✔️:     if (!numNewAtoms || sring.size() - numNewAtoms < 2) { continue; }
    // RDKit✔️✔️:     if (!checkAtomRing(mol, atom, ignore, ring, atomMatcher,
    // RDKit✔️✔️:                        atLeastOneAtomPerRing)) { return false; }
    // RDKit✔️✔️:     if (!checkBondRing(mol, mol.getRingInfo()->bondRings()[i], bondMatcher,
    // RDKit✔️✔️:                        atLeastOneBondPerRing)) { return false; }
    // RDKit✔️✔️:     ringAtoms.insert(diff.begin(), dit);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atLeastOneAtom) {
    // RDKit✔️✔️:     return std::find_if(ringAtoms.begin(), ringAtoms.end(),
    // RDKit✔️✔️:       [&mol, atLeastOneAtom](auto idx) -> bool {
    // RDKit✔️✔️:         return atLeastOneAtom(*mol.getAtomWithIdx(idx));
    // RDKit✔️✔️:       }) != ringAtoms.end();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // Complexity review: both use ordered sets, two ring scans, and identical
    // O(R * r log A) set construction/difference behavior after SSSR.
    let Some(atom) = molecule.atoms().get(atom_index) else {
        return false;
    };
    if atom_matcher.is_some_and(|matches| !matches(atom)) {
        return false;
    }
    let Ok(rings) = crate::find_sssr(molecule) else {
        return false;
    };
    if rings.num_atom_rings(atom.id()) == 0 {
        return false;
    }
    let mut seen = BTreeSet::new();
    for (i, ring) in rings.atom_rings().iter().enumerate() {
        if ring.iter().any(|id| id.index() == atom_index) {
            if !check_atom_ring(
                molecule,
                atom_index,
                &ignore,
                ring,
                atom_matcher,
                atom_per_ring,
            ) || !check_bond_ring(
                molecule,
                &rings.bond_rings()[i],
                bond_matcher,
                bond_per_ring,
            ) {
                return false;
            }
            seen.extend(ring.iter().map(|id| id.index()));
            break;
        }
    }
    for (i, ring) in rings.atom_rings().iter().enumerate() {
        let current: BTreeSet<_> = ring.iter().map(|id| id.index()).collect();
        let difference: Vec<_> = current.difference(&seen).copied().collect();
        if difference.is_empty() || current.len() - difference.len() < 2 {
            continue;
        }
        if !check_atom_ring(
            molecule,
            atom_index,
            &ignore,
            ring,
            atom_matcher,
            atom_per_ring,
        ) || !check_bond_ring(
            molecule,
            &rings.bond_rings()[i],
            bond_matcher,
            bond_per_ring,
        ) {
            return false;
        }
        seen.extend(difference);
    }
    at_least_one_atom
        .is_none_or(|matches| seen.iter().any(|index| matches(&molecule.atoms()[*index])))
}

fn carbocycloalkyl_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarbocycloalkylAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                                 boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return !at.getIsAromatic() && at.getAtomicNum() == 6;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto bondMatcher = [](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return !bnd.getIsAromatic() && bnd.getBondType() == Bond::BondType::SINGLE;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher);
    // RDKit✔️✔️: }
    // Complexity review: O(1) predicate setup plus the canonical fused-ring core.
    let atoms = |a: &Atom| !a.is_aromatic() && a.atomic_number() == 6;
    let bonds = |b: &Bond| !b.is_aromatic() && b.order() == crate::BondOrder::Single;
    fused_ring_match(
        mol,
        atom,
        ignore,
        Some(&atoms),
        Some(&bonds),
        None,
        None,
        None,
    )
}

fn carbocycloalkyl_h_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarbocycloalkylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                                  boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CarbocycloalkylAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || carbocycloalkyl_atom_matcher(mol, atom, ignore)
}

fn carbocycloalkenyl_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarbocycloalkenylAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                                   boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool { return at.getAtomicNum() == 6; };
    // RDKit✔️✔️:   auto atLeastOneBond = [](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return bnd.getIsAromatic() || bnd.getBondType() == Bond::BondType::DOUBLE ||
    // RDKit✔️✔️:            bnd.getBondType() == Bond::BondType::AROMATIC;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   AtomMatcherFunc atLeastOne = nullptr;
    // RDKit✔️✔️:   BondMatcherFunc bondMatcher = nullptr;
    // RDKit✔️✔️:   return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher, atLeastOne,
    // RDKit✔️✔️:                         atLeastOneBond);
    // RDKit✔️✔️: }
    // Complexity review: O(1) predicate setup plus the canonical fused-ring core.
    let atoms = |a: &Atom| a.atomic_number() == 6;
    let unsaturated = |b: &Bond| {
        b.is_aromatic()
            || matches!(
                b.order(),
                crate::BondOrder::Double | crate::BondOrder::Aromatic
            )
    };
    fused_ring_match(
        mol,
        atom,
        ignore,
        Some(&atoms),
        None,
        None,
        Some(&unsaturated),
        None,
    )
}

fn carbocycloalkenyl_h_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarbocycloalkenylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                                    boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CarbocycloalkenylAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || carbocycloalkenyl_atom_matcher(mol, atom, ignore)
}

fn carboaryl_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarboarylAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                           boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getIsAromatic() && at.getAtomicNum() == 6;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto bondMatcher = [](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return bnd.getIsAromatic() || bnd.getBondType() == Bond::BondType::AROMATIC;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher);
    // RDKit✔️✔️: }
    // Complexity review: O(1) predicate setup plus the canonical fused-ring core.
    let atoms = |a: &Atom| a.is_aromatic() && a.atomic_number() == 6;
    let bonds = |b: &Bond| b.is_aromatic() || b.order() == crate::BondOrder::Aromatic;
    fused_ring_match(
        mol,
        atom,
        ignore,
        Some(&atoms),
        Some(&bonds),
        None,
        None,
        None,
    )
}

fn carboaryl_h_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarboarylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                            boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CarboarylAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || carboaryl_atom_matcher(mol, atom, ignore)
}

fn carbocyclic_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarbocyclicAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                             boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool { return at.getAtomicNum() == 6; };
    // RDKit✔️✔️:   return FusedRingMatch(mol, atom, ignore, atomMatcher);
    // RDKit✔️✔️: }
    // Complexity review: O(1) predicate setup plus the canonical fused-ring core.
    let atoms = |a: &Atom| a.atomic_number() == 6;
    fused_ring_match(mol, atom, ignore, Some(&atoms), None, None, None, None)
}

fn carbocyclic_h_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarbocyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                              boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CarbocyclicAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || carbocyclic_atom_matcher(mol, atom, ignore)
}

fn no_carbon_ring_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool NoCarbonRingAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                              boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() != 6;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return FusedRingMatch(mol, atom, ignore, atomMatcher);
    // RDKit✔️✔️: }
    // Complexity review: O(1) predicate setup plus the canonical fused-ring core.
    let atoms = |a: &Atom| a.atomic_number() != 6;
    fused_ring_match(mol, atom, ignore, Some(&atoms), None, None, None, None)
}

fn no_carbon_ring_h_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool NoCarbonRingHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                               boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return NoCarbonRingAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || no_carbon_ring_atom_matcher(mol, atom, ignore)
}

fn heterocyclic_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool HeterocyclicAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                              boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   auto atLeastOne = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() != 6 && at.getAtomicNum() != 1;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   AtomMatcherFunc atomMatcher = nullptr;
    // RDKit✔️✔️:   AtomMatcherFunc oneAtomPerRing = nullptr;
    // RDKit✔️✔️:   BondMatcherFunc bondMatcher = nullptr;
    // RDKit✔️✔️:   BondMatcherFunc oneBondPerRing = nullptr;
    // RDKit✔️✔️:   return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher,
    // RDKit✔️✔️:                         oneAtomPerRing, oneBondPerRing, atLeastOne);
    // RDKit✔️✔️: }
    // Complexity review: O(1) predicate setup plus the canonical fused-ring core.
    let hetero = |a: &Atom| !matches!(a.atomic_number(), 6 | 1);
    fused_ring_match(mol, atom, ignore, None, None, None, None, Some(&hetero))
}

fn heterocyclic_h_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool HeterocyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                               boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return HeterocyclicAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || heterocyclic_atom_matcher(mol, atom, ignore)
}

fn heteroaryl_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool HeteroarylAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                            boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool { return at.getIsAromatic(); };
    // RDKit✔️✔️:   auto bondMatcher = [](const Bond &bnd) -> bool {
    // RDKit✔️✔️:     return bnd.getIsAromatic() || bnd.getBondType() == Bond::BondType::AROMATIC;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   auto atLeastOne = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getAtomicNum() != 6 && at.getAtomicNum() != 1;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   AtomMatcherFunc oneAtomPerRing = nullptr;
    // RDKit✔️✔️:   BondMatcherFunc oneBondPerRing = nullptr;
    // RDKit✔️✔️:   return FusedRingMatch(mol, atom, ignore, atomMatcher, bondMatcher,
    // RDKit✔️✔️:                         oneAtomPerRing, oneBondPerRing, atLeastOne);
    // RDKit✔️✔️: }
    // Complexity review: O(1) predicate setup plus the canonical fused-ring core.
    let aromatic_atom = |a: &Atom| a.is_aromatic();
    let aromatic_bond = |b: &Bond| b.is_aromatic() || b.order() == crate::BondOrder::Aromatic;
    let hetero = |a: &Atom| !matches!(a.atomic_number(), 6 | 1);
    fused_ring_match(
        mol,
        atom,
        ignore,
        Some(&aromatic_atom),
        Some(&aromatic_bond),
        None,
        None,
        Some(&hetero),
    )
}

fn heteroaryl_h_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool HeteroarylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                             boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return HeteroarylAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || heteroaryl_atom_matcher(mol, atom, ignore)
}

fn cyclic_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CyclicAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                        boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool {
    // RDKit✔️✔️:     return at.getOwningMol().getRingInfo()->numAtomRings(at.getIdx()) > 0;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return FusedRingMatch(mol, atom, ignore, atomMatcher);
    // RDKit✔️✔️: }
    // Complexity review: one O(V+E) fast-ring pass plus the canonical core.
    let Ok(rings) = crate::fast_find_rings(mol) else {
        return false;
    };
    let cyclic = |a: &Atom| rings.num_atom_rings(a.id()) > 0;
    fused_ring_match(mol, atom, ignore, Some(&cyclic), None, None, None, None)
}

fn cyclic_h_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                         boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CyclicAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || cyclic_atom_matcher(mol, atom, ignore)
}

fn d_atom_matcher(_mol: &Molecule, atom: usize, mut ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool DAtomMatcher(const ROMol &, const Atom &atom,
    // RDKit✔️✔️:                   boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (atom.getAtomicNum() == 1 && atom.getIsotope() == 2) {
    // RDKit✔️✔️:     ignore.set(atom.getIdx());
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // Complexity review: both implementations perform two O(1) atom-field
    // reads and, on success, one O(1) write to the by-value ignore bitmap.
    let Some(atom) = _mol.atoms().get(atom) else {
        return false;
    };
    if atom.atomic_number() == 1 && atom.isotope() == Some(2) {
        if let Some(bit) = ignore.get_mut(atom.id().index()) {
            *bit = true;
        }
        return true;
    }
    false
}

fn t_atom_matcher(_mol: &Molecule, atom: usize, mut ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool TAtomMatcher(const ROMol &, const Atom &atom,
    // RDKit✔️✔️:                   boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (atom.getAtomicNum() == 1 && atom.getIsotope() == 3) {
    // RDKit✔️✔️:     ignore.set(atom.getIdx());
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // Complexity review: both implementations perform two O(1) atom-field
    // reads and, on success, one O(1) write to the by-value ignore bitmap.
    let Some(atom) = _mol.atoms().get(atom) else {
        return false;
    };
    if atom.atomic_number() == 1 && atom.isotope() == Some(3) {
        if let Some(bit) = ignore.get_mut(atom.id().index()) {
            *bit = true;
        }
        return true;
    }
    false
}

fn hplus_atom_matcher(_mol: &Molecule, atom: usize, mut ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool HplusAtomMatcher(const ROMol &, const Atom &atom,
    // RDKit✔️✔️:                       boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (atom.getAtomicNum() == 1 && atom.getFormalCharge() == 1) {
    // RDKit✔️✔️:     ignore.set(atom.getIdx());
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // Complexity review: both implementations perform two O(1) atom-field
    // reads and, on success, one O(1) write to the by-value ignore bitmap.
    let Some(atom) = _mol.atoms().get(atom) else {
        return false;
    };
    if atom.atomic_number() == 1 && atom.formal_charge() == 1 {
        if let Some(bit) = ignore.get_mut(atom.id().index()) {
            *bit = true;
        }
        return true;
    }
    false
}

fn pol_atom_matcher(_mol: &Molecule, atom: usize, mut ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool PolAtomMatcher(const ROMol &, const Atom &atom,
    // RDKit✔️✔️:                     boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   std::string label;
    // RDKit✔️✔️:   if (atom.getPropIfPresent(common_properties::atomLabel, label) &&
    // RDKit✔️✔️:       label == "Pol") {
    // RDKit✔️✔️:     ignore.set(atom.getIdx());
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // Complexity review: COSMolKit's ordered property-map lookup is O(log P),
    // matching RDKit's property-map lookup class; the success path adds one
    // O(1) write to the copied ignore bitmap and allocates no label string.
    let Some(atom) = _mol.atoms().get(atom) else {
        return false;
    };
    if atom.prop("atomLabel") == Some("Pol") {
        if let Some(bit) = ignore.get_mut(atom.id().index()) {
            *bit = true;
        }
        return true;
    }
    false
}

fn r_atom_matcher(mol: &Molecule, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool RAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                   boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   return GroupHAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: both implementations delegate directly to the same
    // canonical GroupH matcher without another traversal or allocation.
    group_h_atom_matcher(mol, atom, ignore)
}

type GenericMatcher = fn(&Molecule, usize, Vec<bool>) -> bool;

fn generic_matcher_for_label(label: &str) -> Option<GenericMatcher> {
    // RDKit✔️✔️: const static std::map<
    // RDKit✔️✔️:     std::string,
    // RDKit✔️✔️:     std::function<bool(const ROMol &, const Atom &, boost::dynamic_bitset<>)>>
    // RDKit✔️✔️:     genericMatchers = {
    // RDKit✔️✔️:         {"Group", Matchers::GroupAtomMatcher},
    // RDKit✔️✔️:         {"G", Matchers::GroupAtomMatcher},
    // RDKit✔️✔️:         {"GroupH", Matchers::GroupHAtomMatcher},
    // RDKit✔️✔️:         {"GH", Matchers::GroupHAtomMatcher},
    // RDKit✔️✔️:         {"Group*", Matchers::GroupStarAtomMatcher},
    // RDKit✔️✔️:         {"G*", Matchers::GroupStarAtomMatcher},
    // RDKit✔️✔️:         {"GroupH*", Matchers::GroupStarHAtomMatcher},
    // RDKit✔️✔️:         {"GH*", Matchers::GroupStarHAtomMatcher},
    // RDKit✔️✔️:         {"Alkyl", Matchers::AlkylAtomMatcher},
    // RDKit✔️✔️:         {"ALK", Matchers::AlkylAtomMatcher},
    // RDKit✔️✔️:         {"AlkylH", Matchers::AlkylHAtomMatcher},
    // RDKit✔️✔️:         {"ALH", Matchers::AlkylHAtomMatcher},
    // RDKit✔️✔️:         {"Alkenyl", Matchers::AlkenylAtomMatcher},
    // RDKit✔️✔️:         {"AEL", Matchers::AlkenylAtomMatcher},
    // RDKit✔️✔️:         {"AlkenylH", Matchers::AlkenylHAtomMatcher},
    // RDKit✔️✔️:         {"AEH", Matchers::AlkenylHAtomMatcher},
    // RDKit✔️✔️:         {"Alkynyl", Matchers::AlkynylAtomMatcher},
    // RDKit✔️✔️:         {"AYL", Matchers::AlkynylAtomMatcher},
    // RDKit✔️✔️:         {"AlkynylH", Matchers::AlkynylHAtomMatcher},
    // RDKit✔️✔️:         {"AYH", Matchers::AlkynylHAtomMatcher},
    // RDKit✔️✔️:         {"Carbocyclic", Matchers::CarbocyclicAtomMatcher},
    // RDKit✔️✔️:         {"CBC", Matchers::CarbocyclicAtomMatcher},
    // RDKit✔️✔️:         {"CarbocyclicH", Matchers::CarbocyclicHAtomMatcher},
    // RDKit✔️✔️:         {"CBH", Matchers::CarbocyclicHAtomMatcher},
    // RDKit✔️✔️:         {"Carbocycloalkyl", Matchers::CarbocycloalkylAtomMatcher},
    // RDKit✔️✔️:         {"CAL", Matchers::CarbocycloalkylAtomMatcher},
    // RDKit✔️✔️:         {"CarbocycloalkylH", Matchers::CarbocycloalkylHAtomMatcher},
    // RDKit✔️✔️:         {"CAH", Matchers::CarbocycloalkylHAtomMatcher},
    // RDKit✔️✔️:         {"Carbocycloalkenyl", Matchers::CarbocycloalkenylAtomMatcher},
    // RDKit✔️✔️:         {"CEL", Matchers::CarbocycloalkenylAtomMatcher},
    // RDKit✔️✔️:         {"CarbocycloalkenylH", Matchers::CarbocycloalkenylHAtomMatcher},
    // RDKit✔️✔️:         {"CEH", Matchers::CarbocycloalkenylHAtomMatcher},
    // RDKit✔️✔️:         {"Carboaryl", Matchers::CarboarylAtomMatcher},
    // RDKit✔️✔️:         {"ARY", Matchers::CarboarylAtomMatcher},
    // RDKit✔️✔️:         {"CarboarylH", Matchers::CarboarylHAtomMatcher},
    // RDKit✔️✔️:         {"ARH", Matchers::CarboarylHAtomMatcher},
    // RDKit✔️✔️:         {"Cyclic", Matchers::CyclicAtomMatcher},
    // RDKit✔️✔️:         {"CYC", Matchers::CyclicAtomMatcher},
    // RDKit✔️✔️:         {"CyclicH", Matchers::CyclicHAtomMatcher},
    // RDKit✔️✔️:         {"CYH", Matchers::CyclicHAtomMatcher},
    // RDKit✔️✔️:         {"Acyclic", Matchers::AcyclicAtomMatcher},
    // RDKit✔️✔️:         {"ACY", Matchers::AcyclicAtomMatcher},
    // RDKit✔️✔️:         {"AcyclicH", Matchers::AcyclicHAtomMatcher},
    // RDKit✔️✔️:         {"ACH", Matchers::AcyclicHAtomMatcher},
    // RDKit✔️✔️:         {"Carboacyclic", Matchers::CarboacyclicAtomMatcher},
    // RDKit✔️✔️:         {"ABC", Matchers::CarboacyclicAtomMatcher},
    // RDKit✔️✔️:         {"CarboacyclicH", Matchers::CarboacyclicHAtomMatcher},
    // RDKit✔️✔️:         {"ABH", Matchers::CarboacyclicHAtomMatcher},
    // RDKit✔️✔️:         {"Heteroacyclic", Matchers::HeteroacyclicAtomMatcher},
    // RDKit✔️✔️:         {"AHC", Matchers::HeteroacyclicAtomMatcher},
    // RDKit✔️✔️:         {"HeteroacyclicH", Matchers::HeteroacyclicHAtomMatcher},
    // RDKit✔️✔️:         {"AHH", Matchers::HeteroacyclicHAtomMatcher},
    // RDKit✔️✔️:         {"Alkoxy", Matchers::AlkoxyacyclicAtomMatcher},
    // RDKit✔️✔️:         {"AOX", Matchers::AlkoxyacyclicAtomMatcher},
    // RDKit✔️✔️:         {"AlkoxyH", Matchers::AlkoxyacyclicHAtomMatcher},
    // RDKit✔️✔️:         {"AOH", Matchers::AlkoxyacyclicHAtomMatcher},
    // RDKit✔️✔️:         {"Heterocyclic", Matchers::HeterocyclicAtomMatcher},
    // RDKit✔️✔️:         {"CHC", Matchers::HeterocyclicAtomMatcher},
    // RDKit✔️✔️:         {"HeterocyclicH", Matchers::HeterocyclicHAtomMatcher},
    // RDKit✔️✔️:         {"CHH", Matchers::HeterocyclicHAtomMatcher},
    // RDKit✔️✔️:         {"Heteroaryl", Matchers::HeteroarylAtomMatcher},
    // RDKit✔️✔️:         {"HAR", Matchers::HeteroarylAtomMatcher},
    // RDKit✔️✔️:         {"HeteroarylH", Matchers::HeteroarylHAtomMatcher},
    // RDKit✔️✔️:         {"HAH", Matchers::HeteroarylHAtomMatcher},
    // RDKit✔️✔️:         {"NoCarbonRing", Matchers::NoCarbonRingAtomMatcher},
    // RDKit✔️✔️:         {"CXX", Matchers::NoCarbonRingAtomMatcher},
    // RDKit✔️✔️:         {"NoCarbonRingH", Matchers::NoCarbonRingHAtomMatcher},
    // RDKit✔️✔️:         {"CXH", Matchers::NoCarbonRingHAtomMatcher}};
    // Complexity review: lookup has bounded constant work over the same fixed
    // label table and returns a plain function pointer without allocation.
    match label {
        "Group" | "G" => Some(group_atom_matcher),
        "GroupH" | "GH" => Some(group_h_atom_matcher),
        "Group*" | "G*" => Some(group_star_atom_matcher),
        "GroupH*" | "GH*" => Some(group_star_h_atom_matcher),
        "Alkyl" | "ALK" => Some(alkyl_atom_matcher),
        "AlkylH" | "ALH" => Some(alkyl_h_atom_matcher),
        "Alkenyl" | "AEL" => Some(alkenyl_atom_matcher),
        "AlkenylH" | "AEH" => Some(alkenyl_h_atom_matcher),
        "Alkynyl" | "AYL" => Some(alkynyl_atom_matcher),
        "AlkynylH" | "AYH" => Some(alkynyl_h_atom_matcher),
        "Carbocyclic" | "CBC" => Some(carbocyclic_atom_matcher),
        "CarbocyclicH" | "CBH" => Some(carbocyclic_h_atom_matcher),
        "Carbocycloalkyl" | "CAL" => Some(carbocycloalkyl_atom_matcher),
        "CarbocycloalkylH" | "CAH" => Some(carbocycloalkyl_h_atom_matcher),
        "Carbocycloalkenyl" | "CEL" => Some(carbocycloalkenyl_atom_matcher),
        "CarbocycloalkenylH" | "CEH" => Some(carbocycloalkenyl_h_atom_matcher),
        "Carboaryl" | "ARY" => Some(carboaryl_atom_matcher),
        "CarboarylH" | "ARH" => Some(carboaryl_h_atom_matcher),
        "Cyclic" | "CYC" => Some(cyclic_atom_matcher),
        "CyclicH" | "CYH" => Some(cyclic_h_atom_matcher),
        "Acyclic" | "ACY" => Some(acyclic_atom_matcher),
        "AcyclicH" | "ACH" => Some(acyclic_h_atom_matcher),
        "Carboacyclic" | "ABC" => Some(carboacyclic_atom_matcher),
        "CarboacyclicH" | "ABH" => Some(carboacyclic_h_atom_matcher),
        "Heteroacyclic" | "AHC" => Some(heteroacyclic_atom_matcher),
        "HeteroacyclicH" | "AHH" => Some(heteroacyclic_h_atom_matcher),
        "Alkoxy" | "AOX" => Some(alkoxyacyclic_atom_matcher),
        "AlkoxyH" | "AOH" => Some(alkoxyacyclic_h_atom_matcher),
        "Heterocyclic" | "CHC" => Some(heterocyclic_atom_matcher),
        "HeterocyclicH" | "CHH" => Some(heterocyclic_h_atom_matcher),
        "Heteroaryl" | "HAR" => Some(heteroaryl_atom_matcher),
        "HeteroarylH" | "HAH" => Some(heteroaryl_h_atom_matcher),
        "NoCarbonRing" | "CXX" => Some(no_carbon_ring_atom_matcher),
        "NoCarbonRingH" | "CXH" => Some(no_carbon_ring_h_atom_matcher),
        _ => None,
    }
}

#[derive(Debug, thiserror::Error)]
enum GenericGroupsError {
    #[error(transparent)]
    Operation(#[from] crate::OperationError),
    #[error(transparent)]
    RingFinding(#[from] crate::RingFindingError),
}

fn set_generic_queries_from_properties(
    mol: &mut Molecule,
    use_atom_labels: bool,
    use_sgroups: bool,
) {
    // RDKit✔️✔️: void setGenericQueriesFromProperties(ROMol &mol, bool useAtomLabels,
    // RDKit✔️✔️:                                      bool useSGroups) {
    // RDKit✔️✔️:   if (useAtomLabels) {
    // RDKit✔️✔️:     for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:       std::string label = "";
    // RDKit✔️✔️:       if (!atom->getPropIfPresent(common_properties::atomLabel, label)) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // special case for generic groups from mol blocks:
    // RDKit✔️✔️:       if (atom->hasQuery() && !atom->getAtomicNum() &&
    // RDKit✔️✔️:           atom->getQuery()->getDescription() == "AtomAtomicNum" &&
    // RDKit✔️✔️:           !atom->getQuery()->getNegation()) {
    // RDKit✔️✔️:         atom->setQuery(makeAtomNullQuery());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // pseudoatom labels from CXSMILES end with "_p"... strip that if
    // RDKit✔️✔️:       // present
    // RDKit✔️✔️:       if (label.size() > 4 && label.compare(label.size() - 2, 2, "_p") == 0) {
    // RDKit✔️✔️:         label = label.substr(0, label.size() - 2);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (genericMatchers.find(label) != genericMatchers.end()) {
    // RDKit✔️✔️:         atom->setProp(common_properties::_QueryAtomGenericLabel, label);
    // RDKit✔️✔️:         atom->clearProp(common_properties::atomLabel);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (useSGroups) {
    // RDKit✔️✔️:     auto &sgs = getSubstanceGroups(mol);
    // RDKit✔️✔️:     auto iter = sgs.begin();
    // RDKit✔️✔️:     while (iter != sgs.end()) {
    // RDKit✔️✔️:       const auto &sgroup = *iter;
    // RDKit✔️✔️:       if (sgroup.getProp<std::string>("TYPE") == "SUP") {
    // RDKit✔️✔️:         std::string label;
    // RDKit✔️✔️:         if (sgroup.getPropIfPresent("LABEL", label) &&
    // RDKit✔️✔️:             genericMatchers.find(label) != genericMatchers.end()) {
    // RDKit✔️✔️:           for (auto aidx : sgroup.getAtoms()) {
    // RDKit✔️✔️:             mol.getAtomWithIdx(aidx)->setProp(
    // RDKit✔️✔️:                 common_properties::_QueryAtomGenericLabel, label);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           iter = sgs.erase(iter);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           ++iter;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         ++iter;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Complexity review: both paths scan atoms and SGroups once. Each label
    // uses the one fixed-table lookup above; Rust retains consumed groups in a
    // single O(S) buffer and resets typed row ids once, avoiding repeated
    // vector erases while preserving source order and overwrite precedence.
    let topology = mol.topology_block_mut();
    if use_atom_labels {
        for atom in &mut topology.atoms {
            let Some(mut label) = atom.prop("atomLabel").map(str::to_owned) else {
                continue;
            };
            if atom.atomic_number() == 0
                && matches!(
                    atom.query(),
                    Some(QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(0)))
                )
            {
                atom.set_query(Some(make_atom_null_query()));
            }
            if label.len() > 4 && label.ends_with("_p") {
                label.truncate(label.len() - 2);
            }
            if generic_matcher_for_label(&label).is_some() {
                atom.set_prop("_QueryAtomGenericLabel", label);
                atom.clear_prop("atomLabel");
            }
        }
    }
    if use_sgroups {
        let source_groups = std::mem::take(&mut topology.substance_groups);
        let mut retained = Vec::with_capacity(source_groups.len());
        for group in source_groups {
            let label = group.label().map(str::to_owned);
            if group.kind() == &SubstanceGroupKind::Superatom
                && label
                    .as_deref()
                    .is_some_and(|value| generic_matcher_for_label(value).is_some())
            {
                let label = label.expect("validated SGroup label");
                for atom_id in group.atoms() {
                    topology.atoms[atom_id.index()]
                        .set_prop("_QueryAtomGenericLabel", label.clone());
                }
            } else {
                retained.push(group);
            }
        }
        for (index, group) in retained.iter_mut().enumerate() {
            group.set_id(SubstanceGroupId::new(index));
        }
        topology.substance_groups = retained;
    }
}

fn convert_generic_queries_to_substance_groups(mol: &mut Molecule) {
    // RDKit✔️✔️: void convertGenericQueriesToSubstanceGroups(ROMol &mol) {
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     std::string label;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::_QueryAtomGenericLabel,
    // RDKit✔️✔️:                                label)) {
    // RDKit✔️✔️:       SubstanceGroup sg(&mol, "SUP");
    // RDKit✔️✔️:       sg.setProp("LABEL", label);
    // RDKit✔️✔️:       sg.addAtomWithIdx(atom->getIdx());
    // RDKit✔️✔️:       addSubstanceGroup(mol, sg);
    // RDKit✔️✔️:       atom->clearProp(common_properties::_QueryAtomGenericLabel);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Complexity review: both perform one O(A) atom scan, append one typed
    // SGroup per labeled atom, and preserve encounter order without graph
    // traversal or repeated lookup.
    let topology = mol.topology_block_mut();
    for atom_index in 0..topology.atoms.len() {
        let Some(label) = topology.atoms[atom_index]
            .prop("_QueryAtomGenericLabel")
            .map(str::to_owned)
        else {
            continue;
        };
        let id = SubstanceGroupId::new(topology.substance_groups.len());
        topology.substance_groups.push(
            SubstanceGroup::new(id, SubstanceGroupKind::Superatom)
                .with_label(label)
                .with_atoms(vec![topology.atoms[atom_index].id()]),
        );
        topology.atoms[atom_index].clear_prop("_QueryAtomGenericLabel");
    }
}

fn adjust_query_properties_with_generic_groups(
    mol: &Molecule,
) -> Result<Molecule, GenericGroupsError> {
    // RDKit✔️✔️: ROMol *adjustQueryPropertiesWithGenericGroups(
    // RDKit✔️✔️:     const ROMol &mol, const MolOps::AdjustQueryParameters *inParams) {
    // RDKit✔️✔️:   auto *res = new RWMol(mol);
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     adjustQueryProperties(*res, inParams);
    // RDKit✔️✔️:     GenericGroups::setGenericQueriesFromProperties(*res);
    // RDKit✔️✔️:   } catch (MolSanitizeException &se) {
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     throw se;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<ROMol *>(res);
    // RDKit✔️✔️: }
    //
    // RDKit's null `inParams` path uses the default AdjustQueryParameters:
    // aromatize with SYMMRINGS|SETAROMATICITY, convert unlabeled dummies to
    // null queries, and add explicit-degree queries only to non-dummy ring
    // atoms. This private Rust entry models that complete default surface.
    // Complexity review: both clone molecule state, run the same ring and
    // aromaticity passes, then make bounded O(A) atom passes. Rust's COW clone
    // delays topology allocation until sanitize mutates it.
    let mut result =
        mol.sanitize_with_ops(crate::SanitizeOps::SYMMRINGS | crate::SanitizeOps::SET_AROMATICITY)?;
    let rings = crate::symmetrize_sssr(&result)?;
    let degrees: Vec<u8> = (0..result.num_atoms())
        .map(|index| {
            u8::try_from(result.topology_block().adjacency.neighbors_of(index).len())
                .expect("atom degree must fit the modeled u8 query domain")
        })
        .collect();
    let topology = result.topology_block_mut();
    for (index, atom) in topology.atoms.iter_mut().enumerate() {
        let atomic_number = atom.atomic_number();
        if atomic_number == 0 && atom.query().is_none() && atom.isotope().unwrap_or(0) == 0 {
            atom.set_query(Some(make_atom_null_query()));
        }
        if atomic_number != 0 && rings.num_atom_rings(atom.id()) != 0 {
            replace_atom_with_query_atom(atom);
            let query = atom.query_mut().expect("query atom was established");
            query_atom_expand_query(
                query,
                make_atom_explicit_degree_query(degrees[index]),
                CompositeQueryType::And,
                false,
            );
        }
    }
    set_generic_queries_from_properties(&mut result, true, true);
    Ok(result)
}

pub(super) fn generic_atom_matcher(mol: &Molecule, query: &Molecule, atom_match: &[usize]) -> bool {
    // RDKit✔️✔️: bool genericAtomMatcher(const ROMol &mol, const ROMol &query,
    // RDKit✔️✔️:                         const std::span<const unsigned int> &match) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> ignore(mol.getNumAtoms());
    // RDKit✔️✔️:   for (const auto idx : match) {
    // RDKit✔️✔️:     ignore.set(idx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto atom : query.atoms()) {
    // RDKit✔️✔️:     if (atom->getDegree() > 1) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::string genericLabel;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::_QueryAtomGenericLabel,
    // RDKit✔️✔️:                                genericLabel)) {
    // RDKit✔️✔️:       auto found = genericMatchers.find(genericLabel);
    // RDKit✔️✔️:       if (found != genericMatchers.end() &&
    // RDKit✔️✔️:           !found->second(mol, *mol.getAtomWithIdx(match[atom->getIdx()]),
    // RDKit✔️✔️:                          ignore)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // Complexity review: both implementations allocate one O(V) ignore
    // bitmap and scan O(Q) query atoms. The exhaustive Rust match has bounded
    // work over RDKit's same fixed label table, introduces no dynamic registry,
    // and delegates every label and alias to exactly one canonical matcher.
    let mut ignore = vec![false; mol.num_atoms()];
    for &index in atom_match {
        ignore[index] = true;
    }

    for atom in query.atoms() {
        let query_index = atom.id().index();
        if query
            .topology_block()
            .adjacency
            .neighbors_of(query_index)
            .len()
            > 1
        {
            continue;
        }
        let Some(label) = atom.prop("_QueryAtomGenericLabel") else {
            continue;
        };
        let matcher = generic_matcher_for_label(label);
        if matcher.is_some_and(|matches| !matches(mol, atom_match[query_index], ignore.clone())) {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use crate::{
        AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder, SubstanceGroup, SubstanceGroupId,
        SubstanceGroupKind,
    };

    use super::{
        acyclic_atom_matcher, acyclic_h_atom_matcher, adjust_query_properties_with_generic_groups,
        alkenyl_atom_matcher, alkenyl_h_atom_matcher, alkoxyacyclic_atom_matcher,
        alkoxyacyclic_h_atom_matcher, alkyl_atom_matcher, alkyl_h_atom_matcher,
        alkynyl_atom_matcher, alkynyl_h_atom_matcher, all_atoms_match, carboacyclic_atom_matcher,
        carboacyclic_h_atom_matcher, carboaryl_atom_matcher, carboaryl_h_atom_matcher,
        carbocyclic_atom_matcher, carbocyclic_h_atom_matcher, carbocycloalkenyl_atom_matcher,
        carbocycloalkenyl_h_atom_matcher, carbocycloalkyl_atom_matcher,
        carbocycloalkyl_h_atom_matcher, check_atom_ring, check_bond_ring,
        convert_generic_queries_to_substance_groups, cyclic_atom_matcher, cyclic_h_atom_matcher,
        d_atom_matcher, fused_ring_match, generic_atom_matcher, group_atom_matcher,
        group_h_atom_matcher, group_star_atom_matcher, group_star_h_atom_matcher,
        heteroacyclic_atom_matcher, heteroacyclic_h_atom_matcher, heteroaryl_atom_matcher,
        heteroaryl_h_atom_matcher, heterocyclic_atom_matcher, heterocyclic_h_atom_matcher,
        hplus_atom_matcher, is_hydrogen, no_carbon_ring_atom_matcher,
        no_carbon_ring_h_atom_matcher, pol_atom_matcher, r_atom_matcher,
        set_generic_queries_from_properties, t_atom_matcher,
    };

    fn explicit_hydrogen_molecule() -> crate::Molecule {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let isolated_hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .expect("add C-H bond");
        let molecule = builder.build().expect("build explicit-H molecule");
        assert_eq!(isolated_hydrogen.index(), 2);
        molecule
    }

    #[test]
    fn smarts_generic__is_hydrogen() {
        let molecule = explicit_hydrogen_molecule();
        let ignore = vec![false; molecule.num_atoms()];
        assert!(!is_hydrogen(&molecule, 0, ignore.clone()));
        assert!(is_hydrogen(&molecule, 1, ignore.clone()));
        assert!(!is_hydrogen(&molecule, 2, ignore));
    }

    #[test]
    fn smarts_generic__all_atoms_match() {
        let molecule = crate::Molecule::from_smiles("CC=C").expect("parse propene");
        let ignore = vec![false; molecule.num_atoms()];
        let carbon = |atom: &crate::Atom| atom.atomic_number() == 6;
        let single_or_double =
            |bond: &crate::Bond| matches!(bond.order(), BondOrder::Single | BondOrder::Double);
        let double = |bond: &crate::Bond| bond.order() == BondOrder::Double;
        let oxygen = |atom: &crate::Atom| atom.atomic_number() == 8;

        assert!(all_atoms_match(
            &molecule,
            0,
            ignore.clone(),
            Some(&carbon),
            Some(&single_or_double),
            Some(&carbon),
            Some(&double),
        ));
        assert!(!all_atoms_match(
            &molecule,
            0,
            ignore.clone(),
            Some(&carbon),
            Some(&single_or_double),
            Some(&oxygen),
            None,
        ));

        let mut cut_at_attachment = ignore;
        cut_at_attachment[1] = true;
        assert!(!all_atoms_match(
            &molecule,
            0,
            cut_at_attachment,
            Some(&carbon),
            Some(&single_or_double),
            None,
            Some(&double),
        ));
    }

    #[test]
    fn smarts_generic__group_atom_matcher() {
        let ethane = crate::Molecule::from_smiles("CC").expect("parse ethane");
        assert!(group_atom_matcher(
            &ethane,
            0,
            vec![false; ethane.num_atoms()]
        ));

        let hydrogen = explicit_hydrogen_molecule();
        let mut ignore = vec![false; hydrogen.num_atoms()];
        ignore[0] = true;
        assert!(!group_atom_matcher(&hydrogen, 1, ignore));
    }

    #[test]
    fn smarts_generic__group_h_atom_matcher() {
        let molecule = explicit_hydrogen_molecule();
        let ignore = vec![false; molecule.num_atoms()];
        assert!(group_h_atom_matcher(&molecule, 1, ignore.clone()));
        assert!(!group_h_atom_matcher(&molecule, 2, ignore));
    }

    #[test]
    fn smarts_generic__group_star_atom_matcher() {
        let ring = crate::Molecule::from_smiles("C1CC1").expect("parse cyclopropane");
        assert!(group_star_atom_matcher(
            &ring,
            0,
            vec![false; ring.num_atoms()]
        ));

        let chain = crate::Molecule::from_smiles("CCC").expect("parse propane");
        assert!(!group_star_atom_matcher(
            &chain,
            0,
            vec![false; chain.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__group_star_h_atom_matcher() {
        let molecule = explicit_hydrogen_molecule();
        assert!(group_star_h_atom_matcher(
            &molecule,
            1,
            vec![false; molecule.num_atoms()]
        ));
        assert!(!group_star_h_atom_matcher(
            &molecule,
            2,
            vec![false; molecule.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__alkyl_atom_matcher() {
        let chain = crate::Molecule::from_smiles("CCC").expect("parse propane");
        assert!(alkyl_atom_matcher(
            &chain,
            0,
            vec![false; chain.num_atoms()]
        ));
        let alkene = crate::Molecule::from_smiles("C=C").expect("parse ethene");
        assert!(!alkyl_atom_matcher(
            &alkene,
            0,
            vec![false; alkene.num_atoms()]
        ));
        let ring = crate::Molecule::from_smiles("C1CC1").expect("parse cyclopropane");
        assert!(!alkyl_atom_matcher(&ring, 0, vec![false; ring.num_atoms()]));
    }

    #[test]
    fn smarts_generic__alkyl_h_atom_matcher() {
        let molecule = explicit_hydrogen_molecule();
        let mut ignore = vec![false; molecule.num_atoms()];
        ignore[0] = true;
        assert!(alkyl_h_atom_matcher(&molecule, 1, ignore));
    }

    #[test]
    fn smarts_generic__alkenyl_atom_matcher() {
        let alkene = crate::Molecule::from_smiles("CC=C").expect("parse propene");
        assert!(alkenyl_atom_matcher(
            &alkene,
            0,
            vec![false; alkene.num_atoms()]
        ));
        let alkane = crate::Molecule::from_smiles("CCC").expect("parse propane");
        assert!(!alkenyl_atom_matcher(
            &alkane,
            0,
            vec![false; alkane.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__alkenyl_h_atom_matcher() {
        let molecule = explicit_hydrogen_molecule();
        assert!(alkenyl_h_atom_matcher(
            &molecule,
            1,
            vec![false; molecule.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__alkynyl_atom_matcher() {
        let alkyne = crate::Molecule::from_smiles("CC#C").expect("parse propyne");
        assert!(alkynyl_atom_matcher(
            &alkyne,
            0,
            vec![false; alkyne.num_atoms()]
        ));
        let alkene = crate::Molecule::from_smiles("CC=C").expect("parse propene");
        assert!(!alkynyl_atom_matcher(
            &alkene,
            0,
            vec![false; alkene.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__alkynyl_h_atom_matcher() {
        let molecule = explicit_hydrogen_molecule();
        assert!(alkynyl_h_atom_matcher(
            &molecule,
            1,
            vec![false; molecule.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__acyclic_atom_matcher() {
        let chain = crate::Molecule::from_smiles("CO").unwrap();
        assert!(acyclic_atom_matcher(
            &chain,
            0,
            vec![false; chain.num_atoms()]
        ));
        let ring = crate::Molecule::from_smiles("C1CC1").unwrap();
        assert!(!acyclic_atom_matcher(
            &ring,
            0,
            vec![false; ring.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__acyclic_h_atom_matcher() {
        let molecule = explicit_hydrogen_molecule();
        let mut ignore = vec![false; molecule.num_atoms()];
        ignore[0] = true;
        assert!(acyclic_h_atom_matcher(&molecule, 1, ignore));
    }

    #[test]
    fn smarts_generic__carboacyclic_atom_matcher() {
        let chain = crate::Molecule::from_smiles("CC").unwrap();
        assert!(carboacyclic_atom_matcher(
            &chain,
            0,
            vec![false; chain.num_atoms()]
        ));
        let hetero = crate::Molecule::from_smiles("CO").unwrap();
        assert!(!carboacyclic_atom_matcher(
            &hetero,
            0,
            vec![false; hetero.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__carboacyclic_h_atom_matcher() {
        let molecule = explicit_hydrogen_molecule();
        let mut ignore = vec![false; molecule.num_atoms()];
        ignore[0] = true;
        assert!(carboacyclic_h_atom_matcher(&molecule, 1, ignore));
    }

    #[test]
    fn smarts_generic__heteroacyclic_atom_matcher() {
        let hetero = crate::Molecule::from_smiles("CO").unwrap();
        assert!(heteroacyclic_atom_matcher(
            &hetero,
            0,
            vec![false; hetero.num_atoms()]
        ));
        let carbon = crate::Molecule::from_smiles("CC").unwrap();
        assert!(!heteroacyclic_atom_matcher(
            &carbon,
            0,
            vec![false; carbon.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__heteroacyclic_h_atom_matcher() {
        let molecule = explicit_hydrogen_molecule();
        assert!(heteroacyclic_h_atom_matcher(
            &molecule,
            1,
            vec![false; molecule.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__alkoxyacyclic_atom_matcher() {
        let ether = crate::Molecule::from_smiles("COC").unwrap();
        let mut ignore = vec![false; ether.num_atoms()];
        ignore[0] = true;
        ignore[1] = true;
        assert!(alkoxyacyclic_atom_matcher(&ether, 1, ignore));

        let alcohol = crate::Molecule::from_smiles("CO").unwrap();
        let mut ignore = vec![false; alcohol.num_atoms()];
        ignore[0] = true;
        assert!(!alkoxyacyclic_atom_matcher(&alcohol, 1, ignore));
    }

    #[test]
    fn smarts_generic__alkoxyacyclic_h_atom_matcher() {
        let molecule = explicit_hydrogen_molecule();
        assert!(alkoxyacyclic_h_atom_matcher(
            &molecule,
            1,
            vec![false; molecule.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic_check_atom_ring() {
        let molecule = crate::Molecule::from_smiles("C1CO1").unwrap();
        let rings = crate::find_sssr(&molecule).unwrap();
        let carbon = |atom: &crate::Atom| atom.atomic_number() == 6;
        let oxygen = |atom: &crate::Atom| atom.atomic_number() == 8;
        assert!(check_atom_ring(
            &molecule,
            2,
            &[false; 3],
            &rings.atom_rings()[0],
            None,
            Some(&oxygen)
        ));
        assert!(!check_atom_ring(
            &molecule,
            0,
            &[false; 3],
            &rings.atom_rings()[0],
            Some(&carbon),
            None
        ));
    }

    #[test]
    fn smarts_generic_check_bond_ring() {
        let molecule = crate::Molecule::from_smiles("C1=CC1").unwrap();
        let rings = crate::find_sssr(&molecule).unwrap();
        let allowed =
            |bond: &crate::Bond| matches!(bond.order(), BondOrder::Single | BondOrder::Double);
        let double = |bond: &crate::Bond| bond.order() == BondOrder::Double;
        assert!(check_bond_ring(
            &molecule,
            &rings.bond_rings()[0],
            Some(&allowed),
            Some(&double)
        ));
    }

    #[test]
    fn smarts_generic__fused_ring_match() {
        let fused = crate::Molecule::from_smiles("C1CCC2CCCCC2C1").unwrap();
        let carbon = |atom: &crate::Atom| atom.atomic_number() == 6;
        assert!(fused_ring_match(
            &fused,
            0,
            vec![false; fused.num_atoms()],
            Some(&carbon),
            None,
            None,
            None,
            None
        ));
        let chain = crate::Molecule::from_smiles("CCC").unwrap();
        assert!(!fused_ring_match(
            &chain,
            0,
            vec![false; chain.num_atoms()],
            None,
            None,
            None,
            None,
            None
        ));
    }

    #[test]
    fn smarts_generic__carbocycloalkyl_atom_matcher() {
        let mol = crate::Molecule::from_smiles("C1CC1").unwrap();
        assert!(carbocycloalkyl_atom_matcher(
            &mol,
            0,
            vec![false; mol.num_atoms()]
        ));
        let unsat = crate::Molecule::from_smiles("C1=CC1").unwrap();
        assert!(!carbocycloalkyl_atom_matcher(
            &unsat,
            0,
            vec![false; unsat.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__carbocycloalkyl_h_atom_matcher() {
        let mol = explicit_hydrogen_molecule();
        assert!(carbocycloalkyl_h_atom_matcher(
            &mol,
            1,
            vec![false; mol.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__carbocycloalkenyl_atom_matcher() {
        let mol = crate::Molecule::from_smiles("C1=CC1").unwrap();
        assert!(carbocycloalkenyl_atom_matcher(
            &mol,
            0,
            vec![false; mol.num_atoms()]
        ));
        let saturated = crate::Molecule::from_smiles("C1CC1").unwrap();
        assert!(!carbocycloalkenyl_atom_matcher(
            &saturated,
            0,
            vec![false; saturated.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__carbocycloalkenyl_h_atom_matcher() {
        let mol = explicit_hydrogen_molecule();
        assert!(carbocycloalkenyl_h_atom_matcher(
            &mol,
            1,
            vec![false; mol.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__carboaryl_atom_matcher() {
        let mol = crate::Molecule::from_smiles("c1ccccc1").unwrap();
        assert!(carboaryl_atom_matcher(
            &mol,
            0,
            vec![false; mol.num_atoms()]
        ));
        let hetero = crate::Molecule::from_smiles("c1ccncc1").unwrap();
        assert!(!carboaryl_atom_matcher(
            &hetero,
            0,
            vec![false; hetero.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__carboaryl_h_atom_matcher() {
        let mol = explicit_hydrogen_molecule();
        assert!(carboaryl_h_atom_matcher(
            &mol,
            1,
            vec![false; mol.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__carbocyclic_atom_matcher() {
        for smiles in ["C1CC1", "c1ccccc1"] {
            let mol = crate::Molecule::from_smiles(smiles).unwrap();
            assert!(carbocyclic_atom_matcher(
                &mol,
                0,
                vec![false; mol.num_atoms()]
            ));
        }
    }

    #[test]
    fn smarts_generic__carbocyclic_h_atom_matcher() {
        let mol = explicit_hydrogen_molecule();
        assert!(carbocyclic_h_atom_matcher(
            &mol,
            1,
            vec![false; mol.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__no_carbon_ring_atom_matcher() {
        let mol = crate::Molecule::from_smiles("N1NN1").unwrap();
        assert!(no_carbon_ring_atom_matcher(
            &mol,
            0,
            vec![false; mol.num_atoms()]
        ));
        let carbon = crate::Molecule::from_smiles("C1NN1").unwrap();
        assert!(!no_carbon_ring_atom_matcher(
            &carbon,
            1,
            vec![false; carbon.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__no_carbon_ring_h_atom_matcher() {
        let mol = explicit_hydrogen_molecule();
        assert!(no_carbon_ring_h_atom_matcher(
            &mol,
            1,
            vec![false; mol.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__heterocyclic_atom_matcher() {
        let mol = crate::Molecule::from_smiles("C1CO1").unwrap();
        assert!(heterocyclic_atom_matcher(
            &mol,
            0,
            vec![false; mol.num_atoms()]
        ));
        let carbon = crate::Molecule::from_smiles("C1CC1").unwrap();
        assert!(!heterocyclic_atom_matcher(
            &carbon,
            0,
            vec![false; carbon.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__heterocyclic_h_atom_matcher() {
        let mol = explicit_hydrogen_molecule();
        assert!(heterocyclic_h_atom_matcher(
            &mol,
            1,
            vec![false; mol.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__heteroaryl_atom_matcher() {
        let mol = crate::Molecule::from_smiles("c1ccncc1").unwrap();
        assert!(heteroaryl_atom_matcher(
            &mol,
            0,
            vec![false; mol.num_atoms()]
        ));
        let carbon = crate::Molecule::from_smiles("c1ccccc1").unwrap();
        assert!(!heteroaryl_atom_matcher(
            &carbon,
            0,
            vec![false; carbon.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__heteroaryl_h_atom_matcher() {
        let mol = explicit_hydrogen_molecule();
        assert!(heteroaryl_h_atom_matcher(
            &mol,
            1,
            vec![false; mol.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__cyclic_atom_matcher() {
        let mol = crate::Molecule::from_smiles("C1CO1").unwrap();
        assert!(cyclic_atom_matcher(&mol, 0, vec![false; mol.num_atoms()]));
        let chain = crate::Molecule::from_smiles("CCO").unwrap();
        assert!(!cyclic_atom_matcher(
            &chain,
            0,
            vec![false; chain.num_atoms()]
        ));
    }

    #[test]
    fn smarts_generic__cyclic_h_atom_matcher() {
        let mol = explicit_hydrogen_molecule();
        assert!(cyclic_h_atom_matcher(&mol, 1, vec![false; mol.num_atoms()]));
    }

    #[test]
    fn smarts_generic__d_atom_matcher() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::H).with_isotope(2));
        builder.add_atom(AtomSpec::new(Element::H).with_isotope(3));
        builder.add_atom(AtomSpec::new(Element::C).with_isotope(2));
        let mol = builder.build().expect("build isotope matcher molecule");
        let ignore = vec![false; mol.num_atoms()];
        assert!(d_atom_matcher(&mol, 0, ignore.clone()));
        assert!(!d_atom_matcher(&mol, 1, ignore.clone()));
        assert!(!d_atom_matcher(&mol, 2, ignore));
    }

    #[test]
    fn smarts_generic__t_atom_matcher() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::H).with_isotope(3));
        builder.add_atom(AtomSpec::new(Element::H).with_isotope(2));
        builder.add_atom(AtomSpec::new(Element::C).with_isotope(3));
        let mol = builder.build().expect("build isotope matcher molecule");
        let ignore = vec![false; mol.num_atoms()];
        assert!(t_atom_matcher(&mol, 0, ignore.clone()));
        assert!(!t_atom_matcher(&mol, 1, ignore.clone()));
        assert!(!t_atom_matcher(&mol, 2, ignore));
    }

    #[test]
    fn smarts_generic__hplus_atom_matcher() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::H).with_formal_charge(1));
        builder.add_atom(AtomSpec::new(Element::H));
        builder.add_atom(AtomSpec::new(Element::C).with_formal_charge(1));
        let mol = builder.build().expect("build charged matcher molecule");
        let ignore = vec![false; mol.num_atoms()];
        assert!(hplus_atom_matcher(&mol, 0, ignore.clone()));
        assert!(!hplus_atom_matcher(&mol, 1, ignore.clone()));
        assert!(!hplus_atom_matcher(&mol, 2, ignore));
    }

    #[test]
    fn smarts_generic__pol_atom_matcher() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C).with_prop("atomLabel", "Pol"));
        builder.add_atom(AtomSpec::new(Element::C).with_prop("atomLabel", "POL"));
        builder.add_atom(AtomSpec::new(Element::C).with_prop("other", "Pol"));
        let mol = builder.build().expect("build labeled matcher molecule");
        let ignore = vec![false; mol.num_atoms()];
        assert!(pol_atom_matcher(&mol, 0, ignore.clone()));
        assert!(!pol_atom_matcher(&mol, 1, ignore.clone()));
        assert!(!pol_atom_matcher(&mol, 2, ignore));
    }

    #[test]
    fn smarts_generic__r_atom_matcher() {
        let mol = explicit_hydrogen_molecule();
        let ignore = vec![false; mol.num_atoms()];
        assert!(r_atom_matcher(&mol, 0, ignore.clone()));
        assert!(r_atom_matcher(&mol, 1, ignore.clone()));
        assert!(!r_atom_matcher(&mol, 2, ignore));
    }

    #[test]
    fn smarts_generic_generic_atom_matcher() {
        fn terminal_query(label: &str) -> crate::Molecule {
            let mut builder = MoleculeBuilder::new();
            let attachment = builder.add_atom(AtomSpec::new(Element::C));
            let generic = builder
                .add_atom(AtomSpec::new(Element::C).with_prop("_QueryAtomGenericLabel", label));
            builder
                .add_bond(BondSpec::new(attachment, generic, BondOrder::Single))
                .expect("add query bond");
            builder.build().expect("build terminal generic query")
        }

        let alkyl_query = terminal_query("Alkyl");
        let alias_query = terminal_query("ALK");
        let saturated = crate::Molecule::from_smiles("CCC").expect("parse propane");
        let unsaturated = crate::Molecule::from_smiles("CC=C").expect("parse propene");
        assert!(generic_atom_matcher(&saturated, &alkyl_query, &[0, 1]));
        assert!(generic_atom_matcher(&saturated, &alias_query, &[0, 1]));
        assert!(!generic_atom_matcher(&unsaturated, &alkyl_query, &[0, 1]));

        let unknown_query = terminal_query("NotARegisteredGenericGroup");
        assert!(generic_atom_matcher(&unsaturated, &unknown_query, &[0, 1]));

        let mut builder = MoleculeBuilder::new();
        let left = builder.add_atom(AtomSpec::new(Element::C));
        let middle = builder
            .add_atom(AtomSpec::new(Element::C).with_prop("_QueryAtomGenericLabel", "Alkynyl"));
        let right = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(left, middle, BondOrder::Single))
            .expect("add left query bond");
        builder
            .add_bond(BondSpec::new(middle, right, BondOrder::Single))
            .expect("add right query bond");
        let internal_label = builder.build().expect("build internal-label query");
        assert!(generic_atom_matcher(
            &saturated,
            &internal_label,
            &[0, 1, 2]
        ));
    }

    #[test]
    fn smarts_generic_adjust_query_properties_with_generic_groups() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C).with_prop("atomLabel", "ALK"));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a2, a0, BondOrder::Single))
            .unwrap();
        let source = builder.build().expect("build labeled ring query");
        let adjusted = adjust_query_properties_with_generic_groups(&source)
            .expect("adjust default generic query properties");

        assert_eq!(source.atoms()[0].prop("atomLabel"), Some("ALK"));
        assert_eq!(
            adjusted.atoms()[0].prop("_QueryAtomGenericLabel"),
            Some("ALK")
        );
        assert_eq!(adjusted.atoms()[0].prop("atomLabel"), None);
        assert!(adjusted.atoms().iter().all(|atom| atom.query().is_some()));
    }

    #[test]
    fn smarts_generic_convert_generic_queries_to_substance_groups() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C).with_prop("_QueryAtomGenericLabel", "Alkyl"));
        builder.add_atom(AtomSpec::new(Element::O));
        let mut mol = builder.build().expect("build generic query molecule");
        convert_generic_queries_to_substance_groups(&mut mol);

        assert_eq!(mol.atoms()[0].prop("_QueryAtomGenericLabel"), None);
        assert_eq!(mol.substance_groups().len(), 1);
        assert_eq!(
            mol.substance_groups()[0].kind(),
            &SubstanceGroupKind::Superatom
        );
        assert_eq!(mol.substance_groups()[0].label(), Some("Alkyl"));
        assert_eq!(mol.substance_groups()[0].atoms(), &[mol.atoms()[0].id()]);
    }

    #[test]
    fn smarts_generic_set_generic_queries_from_properties() {
        let mut builder = MoleculeBuilder::new();
        let first = builder.add_atom(AtomSpec::new(Element::C).with_prop("atomLabel", "Alkyl_p"));
        let second = builder.add_atom(AtomSpec::new(Element::C).with_prop("atomLabel", "Unknown"));
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_label("Alkenyl")
                    .with_atoms(vec![first]),
            )
            .unwrap();
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                    .with_atoms(vec![second]),
            )
            .unwrap();
        let mut mol = builder.build().expect("build property-backed query");
        set_generic_queries_from_properties(&mut mol, true, true);

        assert_eq!(
            mol.atoms()[0].prop("_QueryAtomGenericLabel"),
            Some("Alkenyl")
        );
        assert_eq!(mol.atoms()[0].prop("atomLabel"), None);
        assert_eq!(mol.atoms()[1].prop("atomLabel"), Some("Unknown"));
        assert_eq!(mol.substance_groups().len(), 1);
        assert_eq!(mol.substance_groups()[0].id(), SubstanceGroupId::new(0));
        assert_eq!(mol.substance_groups()[0].kind(), &SubstanceGroupKind::Data);
    }
}
