//! RDKit generic-group matching for ordinary substructure matches.

use std::borrow::Cow;
use std::collections::{BTreeSet, VecDeque};

use crate::{QueryGraph, SearchTarget, SearchTargetAccess};
use cosmolkit_core::{RingFindingError, RingInfo};
use cosmolkit_model::{Atom, AtomId, Bond, BondId};
use cosmolkit_types::BondOrder;

fn fast_find_rings_detached<'a>(
    target: &'a SearchTarget<'_>,
) -> Result<Cow<'a, RingInfo>, RingFindingError> {
    if let Some(rings) = target.ring_info()
        && rings.is_find_fast_or_better()
    {
        return Ok(Cow::Borrowed(rings));
    }
    cosmolkit_core::fast_find_rings_from_parts(
        target.num_atoms(),
        target.bonds(),
        target.adjacency(),
    )
    .map(Cow::Owned)
}

fn find_sssr_detached<'a>(
    target: &'a SearchTarget<'_>,
) -> Result<Cow<'a, RingInfo>, RingFindingError> {
    if let Some(rings) = target.ring_info()
        && rings.is_sssr_or_better()
    {
        return Ok(Cow::Borrowed(rings));
    }
    cosmolkit_core::find_sssr_from_parts(target.num_atoms(), target.bonds(), target.adjacency())
        .map(Cow::Owned)
}

type AtomMatcher<'a> = dyn Fn(&Atom) -> bool + 'a;
type BondMatcher<'a> = dyn Fn(&Bond) -> bool + 'a;

fn is_hydrogen(molecule: &SearchTarget<'_>, atom_index: usize, mut ignore: Vec<bool>) -> bool {
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
    molecule: &SearchTarget<'_>,
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

fn group_atom_matcher(molecule: &SearchTarget<'_>, atom_index: usize, ignore: Vec<bool>) -> bool {
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
    let _ring_info = fast_find_rings_detached(molecule).ok();
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

fn group_h_atom_matcher(molecule: &SearchTarget<'_>, atom_index: usize, ignore: Vec<bool>) -> bool {
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

fn group_star_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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
    let Ok(ring_info) = fast_find_rings_detached(molecule) else {
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

fn group_star_h_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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

fn alkyl_atom_matcher(molecule: &SearchTarget<'_>, atom_index: usize, ignore: Vec<bool>) -> bool {
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
    let Ok(ring_info) = fast_find_rings_detached(molecule) else {
        return false;
    };
    let atoms = |atom: &Atom| !atom.is_aromatic() && matches!(atom.atomic_number(), 6 | 1);
    let carbon = |atom: &Atom| atom.atomic_number() == 6;
    let bonds = |bond: &Bond| {
        bond.order() == BondOrder::Single
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

fn alkyl_h_atom_matcher(molecule: &SearchTarget<'_>, atom_index: usize, ignore: Vec<bool>) -> bool {
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
    let Ok(ring_info) = fast_find_rings_detached(molecule) else {
        return false;
    };
    let atoms = |atom: &Atom| !atom.is_aromatic() && matches!(atom.atomic_number(), 6 | 1);
    let bonds = |bond: &Bond| {
        bond.order() == BondOrder::Single
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
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
    extra_bond_type: BondOrder,
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
    let Ok(ring_info) = fast_find_rings_detached(molecule) else {
        return false;
    };
    let atoms = |atom: &Atom| !atom.is_aromatic() && matches!(atom.atomic_number(), 6 | 1);
    let bonds =
        |bond: &Bond| matches!(bond.order(), BondOrder::Single) || bond.order() == extra_bond_type;
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

fn alkenyl_atom_matcher(molecule: &SearchTarget<'_>, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkenylAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                         boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   return UnsatAlkXAtomMatcher(mol, atom, ignore, Bond::BondType::DOUBLE);
    // RDKit✔️✔️: }
    // Complexity review: both are O(1) wrappers around the O(V + E) core.
    unsat_alk_x_atom_matcher(molecule, atom_index, ignore, BondOrder::Double)
}

fn alkenyl_h_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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
        || unsat_alk_x_atom_matcher(molecule, atom_index, ignore, BondOrder::Double)
}

fn alkynyl_atom_matcher(molecule: &SearchTarget<'_>, atom_index: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool AlkynylAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                         boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   return UnsatAlkXAtomMatcher(mol, atom, ignore, Bond::BondType::TRIPLE);
    // RDKit✔️✔️: }
    // Complexity review: both are O(1) wrappers around the O(V + E) core.
    unsat_alk_x_atom_matcher(molecule, atom_index, ignore, BondOrder::Triple)
}

fn alkynyl_h_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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
        || unsat_alk_x_atom_matcher(molecule, atom_index, ignore, BondOrder::Triple)
}

fn acyclic_atom_matcher(molecule: &SearchTarget<'_>, atom_index: usize, ignore: Vec<bool>) -> bool {
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
    let Ok(rings) = fast_find_rings_detached(molecule) else {
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

fn acyclic_h_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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
    let Ok(rings) = fast_find_rings_detached(molecule) else {
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

fn carboacyclic_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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
    let Ok(rings) = fast_find_rings_detached(molecule) else {
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

fn carboacyclic_h_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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
    let Ok(rings) = fast_find_rings_detached(molecule) else {
        return false;
    };
    let atoms =
        |atom: &Atom| matches!(atom.atomic_number(), 6 | 1) && rings.num_atom_rings(atom.id()) == 0;
    all_atoms_match(molecule, atom_index, ignore, Some(&atoms), None, None, None)
}

fn heteroacyclic_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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
    let Ok(rings) = fast_find_rings_detached(molecule) else {
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

fn heteroacyclic_h_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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

fn alkoxyacyclic_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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
    let Ok(rings) = fast_find_rings_detached(molecule) else {
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
        bond.order() == BondOrder::Single
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

fn alkoxyacyclic_h_atom_matcher(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: Vec<bool>,
) -> bool {
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
    molecule: &SearchTarget<'_>,
    atom_index: usize,
    ignore: &[bool],
    ring: &[AtomId],
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
    molecule: &SearchTarget<'_>,
    ring: &[BondId],
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
    molecule: &SearchTarget<'_>,
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
    let Ok(rings) = find_sssr_detached(molecule) else {
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

fn carbocycloalkyl_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
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
    let bonds = |b: &Bond| !b.is_aromatic() && b.order() == BondOrder::Single;
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

fn carbocycloalkyl_h_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarbocycloalkylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                                  boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CarbocycloalkylAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || carbocycloalkyl_atom_matcher(mol, atom, ignore)
}

fn carbocycloalkenyl_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
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
    let unsaturated =
        |b: &Bond| b.is_aromatic() || matches!(b.order(), BondOrder::Double | BondOrder::Aromatic);
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

fn carbocycloalkenyl_h_atom_matcher(
    mol: &SearchTarget<'_>,
    atom: usize,
    ignore: Vec<bool>,
) -> bool {
    // RDKit✔️✔️: bool CarbocycloalkenylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                                    boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CarbocycloalkenylAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || carbocycloalkenyl_atom_matcher(mol, atom, ignore)
}

fn carboaryl_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
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
    let bonds = |b: &Bond| b.is_aromatic() || b.order() == BondOrder::Aromatic;
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

fn carboaryl_h_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarboarylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                            boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CarboarylAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || carboaryl_atom_matcher(mol, atom, ignore)
}

fn carbocyclic_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarbocyclicAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                             boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   auto atomMatcher = [](const Atom &at) -> bool { return at.getAtomicNum() == 6; };
    // RDKit✔️✔️:   return FusedRingMatch(mol, atom, ignore, atomMatcher);
    // RDKit✔️✔️: }
    // Complexity review: O(1) predicate setup plus the canonical fused-ring core.
    let atoms = |a: &Atom| a.atomic_number() == 6;
    fused_ring_match(mol, atom, ignore, Some(&atoms), None, None, None, None)
}

fn carbocyclic_h_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CarbocyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                              boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CarbocyclicAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || carbocyclic_atom_matcher(mol, atom, ignore)
}

fn no_carbon_ring_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
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

fn no_carbon_ring_h_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool NoCarbonRingHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                               boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return NoCarbonRingAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || no_carbon_ring_atom_matcher(mol, atom, ignore)
}

fn heterocyclic_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
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

fn heterocyclic_h_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool HeterocyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                               boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return HeterocyclicAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || heterocyclic_atom_matcher(mol, atom, ignore)
}

fn heteroaryl_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
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
    let aromatic_bond = |b: &Bond| b.is_aromatic() || b.order() == BondOrder::Aromatic;
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

fn heteroaryl_h_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool HeteroarylHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                             boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return HeteroarylAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || heteroaryl_atom_matcher(mol, atom, ignore)
}

fn cyclic_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
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
    let Ok(rings) = fast_find_rings_detached(mol) else {
        return false;
    };
    let cyclic = |a: &Atom| rings.num_atom_rings(a.id()) > 0;
    fused_ring_match(mol, atom, ignore, Some(&cyclic), None, None, None, None)
}

fn cyclic_h_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool CyclicHAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                         boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   if (IsHydrogen(mol, atom, ignore)) { return true; }
    // RDKit✔️✔️:   return CyclicAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: O(1) hydrogen check plus the canonical matcher.
    is_hydrogen(mol, atom, ignore.clone()) || cyclic_atom_matcher(mol, atom, ignore)
}

fn d_atom_matcher(_mol: &SearchTarget<'_>, atom: usize, mut ignore: Vec<bool>) -> bool {
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

fn t_atom_matcher(_mol: &SearchTarget<'_>, atom: usize, mut ignore: Vec<bool>) -> bool {
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

fn hplus_atom_matcher(_mol: &SearchTarget<'_>, atom: usize, mut ignore: Vec<bool>) -> bool {
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

fn pol_atom_matcher(_mol: &SearchTarget<'_>, atom: usize, mut ignore: Vec<bool>) -> bool {
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

fn r_atom_matcher(mol: &SearchTarget<'_>, atom: usize, ignore: Vec<bool>) -> bool {
    // RDKit✔️✔️: bool RAtomMatcher(const ROMol &mol, const Atom &atom,
    // RDKit✔️✔️:                   boost::dynamic_bitset<> ignore) {
    // RDKit✔️✔️:   return GroupHAtomMatcher(mol, atom, ignore);
    // RDKit✔️✔️: }
    // Complexity review: both implementations delegate directly to the same
    // canonical GroupH matcher without another traversal or allocation.
    group_h_atom_matcher(mol, atom, ignore)
}

type GenericMatcher = for<'a> fn(&SearchTarget<'a>, usize, Vec<bool>) -> bool;

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

pub(super) fn generic_atom_matcher(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    atom_match: &[usize],
) -> bool {
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
        if query.adjacency().get(query_index).map_or(0, Vec::len) > 1 {
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
