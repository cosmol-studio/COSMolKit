//! RDKit-aligned conjugation assignment over detached topology values.

use cosmolkit_model::{Atom, AtomId, BondId, TopologyBlock, TopologyValidationError};

use crate::{
    AromaticityError, ValenceAssignment, ValenceError, bond_valence_contrib,
    periodic_table_outer_electrons, required_valence_list,
};

#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum ConjugationError {
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("valence assignment field {field} has {actual} rows; expected {expected}")]
    ValenceAssignmentLength {
        field: &'static str,
        actual: usize,
        expected: usize,
    },
    #[error("invalid valence row for atom {atom}: {field}={value}")]
    InvalidValenceRow {
        atom: AtomId,
        field: &'static str,
        value: i32,
    },
    #[error("atom {atom} is outside {atom_count} topology atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("bond {bond} is outside {bond_count} topology bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
    #[error("integer overflow while computing {field}")]
    IntegerOverflow { field: &'static str },
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

fn validate_valence(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<(), ConjugationError> {
    let expected = topology.atoms.len();
    for (field, rows) in [
        ("explicit_valence", &valence.explicit_valence),
        ("implicit_hydrogens", &valence.implicit_hydrogens),
    ] {
        if rows.len() != expected {
            return Err(ConjugationError::ValenceAssignmentLength {
                field,
                actual: rows.len(),
                expected,
            });
        }
        if let Some((row, &value)) = rows.iter().enumerate().find(|(_, value)| **value < 0) {
            return Err(ConjugationError::InvalidValenceRow {
                atom: AtomId::new(row),
                field,
                value,
            });
        }
    }
    Ok(())
}

fn atom(topology: &TopologyBlock, atom_id: AtomId) -> Result<&Atom, ConjugationError> {
    topology
        .atoms
        .get(atom_id.index())
        .ok_or(ConjugationError::AtomOutOfRange {
            atom: atom_id,
            atom_count: topology.atoms.len(),
        })
}

fn total_valence(valence: &ValenceAssignment, atom_id: AtomId) -> Result<i32, ConjugationError> {
    let explicit = *valence.explicit_valence.get(atom_id.index()).ok_or(
        ConjugationError::ValenceAssignmentLength {
            field: "explicit_valence",
            actual: valence.explicit_valence.len(),
            expected: atom_id.index() + 1,
        },
    )?;
    let implicit = *valence.implicit_hydrogens.get(atom_id.index()).ok_or(
        ConjugationError::ValenceAssignmentLength {
            field: "implicit_hydrogens",
            actual: valence.implicit_hydrogens.len(),
            expected: atom_id.index() + 1,
        },
    )?;
    explicit
        .checked_add(implicit)
        .ok_or(ConjugationError::IntegerOverflow {
            field: "total valence",
        })
}

fn total_substitutions(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_id: AtomId,
) -> Result<usize, ConjugationError> {
    let hydrogens = usize::try_from(crate::hcount::total_hydrogen_count_from_validated(
        topology, valence, atom_id, false,
    )?)
    .map_err(|_| ConjugationError::IntegerOverflow {
        field: "degree plus total hydrogens",
    })?;
    topology
        .adjacency
        .neighbors_of(atom_id.index())
        .len()
        .checked_add(hydrogens)
        .ok_or(ConjugationError::IntegerOverflow {
            field: "degree plus total hydrogens",
        })
}

fn is_atom_conjugation_candidate(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_id: AtomId,
) -> Result<bool, ConjugationError> {
    // BEGIN RDKIT CPP FUNCTION isAtomConjugCand
    // RDKit✔️✔️: bool isAtomConjugCand(const Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   // return false for neutral atoms where the current valence exceeds the
    // RDKit✔️✔️:   // minimal valence for the atom. logic: if we're hypervalent we aren't
    // RDKit✔️✔️:   // conjugated
    // RDKit✔️✔️:   const auto &vals =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
    // RDKit✔️✔️:   if (!at->getFormalCharge() && vals.front() >= 0 &&
    // RDKit✔️✔️:       at->getTotalValence() > static_cast<unsigned int>(vals.front())) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // the second check here is for Issue211, where the c-P bonds in
    // RDKit✔️✔️:   // Pc1ccccc1 were being marked as conjugated.  This caused the P atom
    // RDKit✔️✔️:   // itself to be SP2 hybridized.  This is wrong.  For now we'll do a quick
    // RDKit✔️✔️:   // hack and forbid this check from adding conjugation to anything out of
    // RDKit✔️✔️:   // the first row of the periodic table.  (Conjugation in aromatic rings
    // RDKit✔️✔️:   // has already been attended to, so this is safe.)
    // RDKit✔️✔️:   int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
    // RDKit✔️✔️:   auto res = ((at->getAtomicNum() <= 10) || (nouter != 5 && nouter != 6) ||
    // RDKit✔️✔️:               (nouter == 6 && at->getTotalDegree() < 2u)) &&
    // RDKit✔️✔️:              MolOps::countAtomElec(at) > 0;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isAtomConjugCand
    let current = atom(topology, atom_id)?;
    let minimum_valence = required_valence_list(current.atomic_number())?[0];
    if current.formal_charge() == 0
        && minimum_valence >= 0
        && total_valence(valence, atom_id)? > minimum_valence
    {
        return Ok(false);
    }
    let outer = periodic_table_outer_electrons(current.atomic_number())?;
    let issue_211_gate = current.atomic_number() <= 10
        || (outer != 5 && outer != 6)
        || (outer == 6 && total_substitutions(topology, valence, atom_id)? < 2);
    Ok(issue_211_gate && crate::aromaticity::count_atom_electrons(topology, valence, atom_id)? > 0)
}

// BEGIN RECOVERY CHEM-15 SOURCE ConjAtomInfo
// RDKit❗❌: struct ConjAtomInfo {
// RDKit❗❌:   unsigned int numSubstituents;
// RDKit❗❌:   bool isCandidate;
// RDKit❗❌: };
// END RECOVERY CHEM-15 SOURCE ConjAtomInfo
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ConjugationAtomInfo {
    num_substituents: usize,
    is_candidate: bool,
}

fn build_conjugation_atom_info(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<Vec<ConjugationAtomInfo>, ConjugationError> {
    // BEGIN RECOVERY CHEM-15 SOURCE setConjugation_cache_transform
    // RDKit❗❌:   std::vector<ConjAtomInfo> atomInfo;
    // RDKit❗❌:   atomInfo.reserve(mol.getNumAtoms());
    // RDKit❗❌:   std::ranges::transform(
    // RDKit❗❌:       mol.atoms(), std::back_inserter(atomInfo), [](const auto atom) {
    // RDKit❗❌:         const auto isCandidate = isAtomConjugCand(atom);
    // RDKit❗❌:         return ConjAtomInfo{
    // RDKit❗❌:             isCandidate ? atom->getDegree() + atom->getTotalNumHs() : 0u,
    // RDKit❗❌:             isCandidate};
    // RDKit❗❌:       });
    // END RECOVERY CHEM-15 SOURCE setConjugation_cache_transform
    // One reserved Vec and one index-order evaluation per atom, matching the
    // source transform. Preserve existing typed valence/H-count boundaries;
    // no candidate re-evaluation occurs inside either nested neighbor loop.
    let mut atom_info = Vec::with_capacity(topology.atoms.len());
    for atom_index in 0..topology.atoms.len() {
        let atom_id = AtomId::new(atom_index);
        let is_candidate = is_atom_conjugation_candidate(topology, valence, atom_id)?;
        let num_substituents = if is_candidate {
            total_substitutions(topology, valence, atom_id)?
        } else {
            0
        };
        atom_info.push(ConjugationAtomInfo {
            num_substituents,
            is_candidate,
        });
        #[cfg(test)]
        recovery_chem15_trace::record(atom_index, num_substituents, is_candidate);
    }
    Ok(atom_info)
}

fn mark_conjugated_atom_bonds(
    source: &TopologyBlock,
    atom_info: &[ConjugationAtomInfo],
    write_flag: &mut impl FnMut(usize, bool),
    atom_id: AtomId,
) -> Result<(), ConjugationError> {
    // BEGIN RECOVERY CHEM-15 SOURCE markConjAtomBonds_cached
    // RDKit❗❌: void markConjAtomBonds(Atom *at,
    // RDKit❗❌:                        const std::span<const ConjAtomInfo> atomInfo) {
    // RDKit❗❌:   PRECONDITION(at, "bad atom");
    // RDKit❗❌:   const auto &info = atomInfo[at->getIdx()];
    // RDKit❗❌:   if (!info.isCandidate) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   auto &mol = at->getOwningMol();
    // RDKit❗❌:
    // RDKit❗❌:   const auto atx = at->getIdx();
    // RDKit❗❌:   // make sure that have either 2 or 3 substitutions on this atom
    // RDKit❗❌:   if ((info.numSubstituents < 2) || (info.numSubstituents > 3)) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto bnd1 : mol.atomBonds(at)) {
    // RDKit❗❌:     if (bnd1->getValenceContrib(at) < 1.5 ||
    // RDKit❗❌:         !atomInfo[bnd1->getOtherAtomIdx(atx)].isCandidate) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     for (const auto bnd2 : mol.atomBonds(at)) {
    // RDKit❗❌:       if (bnd1 == bnd2) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       const auto at2Idx = bnd2->getOtherAtomIdx(atx);
    // RDKit❗❌:       const auto &at2Info = atomInfo[at2Idx];
    // RDKit❗❌:       if (at2Info.numSubstituents > 3) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (at2Info.isCandidate) {
    // RDKit❗❌:         bnd1->setIsConjugated(true);
    // RDKit❗❌:         bnd2->setIsConjugated(true);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RECOVERY CHEM-15 SOURCE markConjAtomBonds_cached
    let info = &atom_info[atom_id.index()];
    if !info.is_candidate {
        return Ok(());
    }
    if !(2..=3).contains(&info.num_substituents) {
        return Ok(());
    }
    let incident = source.adjacency.neighbors_of(atom_id.index());
    for first in incident {
        let first_bond =
            source
                .bonds
                .get(first.bond.index())
                .ok_or(ConjugationError::BondOutOfRange {
                    bond: first.bond,
                    bond_count: source.bonds.len(),
                })?;
        if bond_valence_contrib(first_bond, atom_id)? < 1.5
            || !atom_info[first.atom_index].is_candidate
        {
            continue;
        }
        for second in incident {
            if first.bond == second.bond {
                continue;
            }
            let second_info = &atom_info[second.atom_index];
            if second_info.num_substituents > 3 {
                continue;
            }
            if second_info.is_candidate {
                write_flag(first.bond.index(), true);
                write_flag(second.bond.index(), true);
            }
        }
    }
    Ok(())
}

pub fn atom_has_conjugated_bond(
    topology: &TopologyBlock,
    atom_id: AtomId,
) -> Result<bool, ConjugationError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::atomHasConjugatedBond
    // RDKit✔️❌: bool atomHasConjugatedBond(const Atom *at) {
    // RDKit✔️❌:   PRECONDITION(at, "bad atom");
    // RDKit✔️❌:
    // RDKit✔️❌:   auto &mol = at->getOwningMol();
    // RDKit✔️❌:   for (const auto bnd : mol.atomBonds(at)) {
    // RDKit✔️❌:     if (bnd->getIsConjugated()) {
    // RDKit✔️❌:       return true;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return false;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::atomHasConjugatedBond
    // The detached predicate validates the complete topology before trusting
    // adjacency, adding O(V + E) work to the source's O(degree) scan.
    topology.validate()?;
    atom(topology, atom_id)?;
    Ok(atom_has_conjugated_bond_from_validated(topology, atom_id))
}

pub(crate) fn atom_has_conjugated_bond_from_validated(
    topology: &TopologyBlock,
    atom_id: AtomId,
) -> bool {
    topology
        .adjacency
        .neighbors_of(atom_id.index())
        .iter()
        .any(|neighbor| topology.bonds[neighbor.bond.index()].is_conjugated())
}

// The original detached transform retains its one owned output and validation order.
pub fn assign_conjugation(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<TopologyBlock, ConjugationError> {
    topology.validate()?;
    validate_valence(topology, valence)?;
    let mut working = topology.clone();
    assign_conjugation_with_writer(topology, valence, &mut |index, value| {
        working.bonds[index].set_conjugated(value)
    })?;
    working.validate()?;
    Ok(working)
}
/// Detached bond flags for consumers that do not mutate topology.
/// Uses the same source loop and candidate/valence rules as the topology transform.
pub fn assign_conjugation_flags(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<Vec<bool>, ConjugationError> {
    topology.validate()?;
    validate_valence(topology, valence)?;
    let mut flags = vec![false; topology.bonds.len()];
    assign_conjugation_with_writer(topology, valence, &mut |index, value| flags[index] = value)?;
    Ok(flags)
}
fn assign_conjugation_with_writer(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    writer: &mut impl FnMut(usize, bool),
) -> Result<(), ConjugationError> {
    // BEGIN RECOVERY CHEM-15 SOURCE MolOps_setConjugation
    // RDKit❗❌: void setConjugation(ROMol &mol) {
    // RDKit❗❌:   // start with all bonds being marked unconjugated
    // RDKit❗❌:   // except for aromatic bonds
    // RDKit❗❌:   for (auto bond : mol.bonds()) {
    // RDKit❗❌:     bond->setIsConjugated(bond->getIsAromatic());
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<ConjAtomInfo> atomInfo;
    // RDKit❗❌:   atomInfo.reserve(mol.getNumAtoms());
    // RDKit❗❌:   std::ranges::transform(
    // RDKit❗❌:       mol.atoms(), std::back_inserter(atomInfo), [](const auto atom) {
    // RDKit❗❌:         const auto isCandidate = isAtomConjugCand(atom);
    // RDKit❗❌:         return ConjAtomInfo{
    // RDKit❗❌:             isCandidate ? atom->getDegree() + atom->getTotalNumHs() : 0u,
    // RDKit❗❌:             isCandidate};
    // RDKit❗❌:       });
    // RDKit❗❌:
    // RDKit❗❌:
    // RDKit❗❌:   // loop over each atom and check if the bonds connecting to it can
    // RDKit❗❌:   // be conjugated
    // RDKit❗❌:   for (auto atom : mol.atoms()) {
    // RDKit❗❌:     markConjAtomBonds(atom, atomInfo);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RECOVERY CHEM-15 SOURCE MolOps_setConjugation
    // Both public result forms and all conformer/DG callers share this owner.
    // Source order is reset every bond, build all atom rows, then mark atoms.
    // Extra detached validation/cloning and typed-error behavior remain baseline;
    // cache removes repeated candidate/H-count work only, not all old costs.
    for (index, bond) in topology.bonds.iter().enumerate() {
        writer(index, bond.is_aromatic());
    }
    let atom_info = build_conjugation_atom_info(topology, valence)?;
    for atom_index in 0..topology.atoms.len() {
        mark_conjugated_atom_bonds(topology, &atom_info, writer, AtomId::new(atom_index))?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "tests/conjugation.rs"]
mod source_tests;

#[cfg(test)]
#[path = "tests/conformer_shared_conjugation_assignments.rs"]
mod shared_assignment_tests;

// Passive private observation is absent from release builds and does not
// control production behavior. Only completed cache rows are observed.
#[cfg(test)]
mod recovery_chem15_trace {
    use std::cell::RefCell;
    std::thread_local! {
        static ROWS: RefCell<Option<Vec<(usize, usize, bool)>>> = const { RefCell::new(None) };
    }
    pub(super) fn record(index: usize, count: usize, candidate: bool) {
        ROWS.with(|rows| {
            if let Some(rows) = &mut *rows.borrow_mut() {
                rows.push((index, count, candidate));
            }
        });
    }
    pub(super) fn capture<T>(f: impl FnOnce() -> T) -> (T, Vec<(usize, usize, bool)>) {
        struct Reset;
        impl Drop for Reset {
            fn drop(&mut self) {
                ROWS.with(|rows| *rows.borrow_mut() = None);
            }
        }
        ROWS.with(|rows| {
            assert!(rows.borrow().is_none());
            *rows.borrow_mut() = Some(Vec::new());
        });
        let reset = Reset;
        let output = f();
        let rows = ROWS.with(|rows| rows.borrow_mut().take().unwrap());
        drop(reset);
        (output, rows)
    }
}

#[cfg(test)]
mod recovery_chem15 {
    use super::*;
    use cosmolkit_model::{AtomSpec, Bond, BondSpec};
    use cosmolkit_types::{BondOrder, Element};
    fn graph(specs: Vec<AtomSpec>, edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
        let atoms = specs
            .into_iter()
            .enumerate()
            .map(|(i, s)| Atom::from_spec(AtomId::new(i), s))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b, o))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), o),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn assignment(explicit: &[i32], implicit: &[i32]) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: explicit.to_vec(),
            implicit_hydrogens: implicit.to_vec(),
        }
    }
    fn check(g: TopologyBlock, v: ValenceAssignment, rows: &[(usize, bool)], flags: &[bool]) {
        let before = g.clone();
        let wanted = rows
            .iter()
            .enumerate()
            .map(|(i, &(n, c))| (i, n, c))
            .collect::<Vec<_>>();
        let (result, observed) =
            recovery_chem15_trace::capture(|| assign_conjugation(&g, &v).unwrap());
        assert_eq!(
            observed, wanted,
            "one cache row per actual atom in source index order"
        );
        let (bits, observed_flags) =
            recovery_chem15_trace::capture(|| assign_conjugation_flags(&g, &v).unwrap());
        assert_eq!(
            observed_flags, wanted,
            "shared conformer/DG entry consumes same single cache"
        );
        assert_eq!(bits, flags);
        let mut expected = before.clone();
        for (bond, &flag) in expected.bonds.iter_mut().zip(flags) {
            bond.set_conjugated(flag);
        }
        assert_eq!(result, expected);
        assert_eq!(g, before);
    }
    #[test]
    fn all_atoms_cached_once_including_isolated_and_false_zero_fields() {
        check(
            graph(
                vec![AtomSpec::new(Element::C); 4]
                    .into_iter()
                    .chain([AtomSpec::new(Element::HE)])
                    .collect(),
                &[
                    (0, 1, BondOrder::Double),
                    (1, 2, BondOrder::Single),
                    (2, 3, BondOrder::Double),
                ],
            ),
            assignment(&[2, 3, 3, 2, 0], &[2, 1, 1, 2, 0]),
            &[(3, true), (3, true), (3, true), (3, true), (0, false)],
            &[true, true, true],
        );
        check(
            graph(vec![AtomSpec::new(Element::HE)], &[]),
            assignment(&[0], &[0]),
            &[(0, false)],
            &[],
        );
        check(
            graph(
                vec![AtomSpec::new(Element::C); 2],
                &[(0, 1, BondOrder::Single)],
            ),
            assignment(&[1, 1], &[3, 3]),
            &[(0, false), (0, false)],
            &[false],
        );
    }
    #[test]
    fn lone_pair_candidates_share_cached_substitutions_for_both_outputs() {
        check(
            graph(
                vec![
                    AtomSpec::new(Element::N),
                    AtomSpec::new(Element::C),
                    AtomSpec::new(Element::O),
                ],
                &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Double)],
            ),
            assignment(&[1, 3, 2], &[2, 1, 0]),
            &[(3, true), (3, true), (1, true)],
            &[true, true],
        );
    }
    #[test]
    fn branching_reuses_center_and_neighbor_rows_without_recomputation() {
        check(
            graph(
                vec![AtomSpec::new(Element::C); 6],
                &[
                    (0, 1, BondOrder::Double),
                    (1, 2, BondOrder::Single),
                    (2, 3, BondOrder::Double),
                    (1, 4, BondOrder::Single),
                    (4, 5, BondOrder::Double),
                ],
            ),
            assignment(&[2, 4, 3, 2, 3, 2], &[2, 0, 1, 2, 1, 2]),
            &[(3, true); 6],
            &[true; 5],
        );
    }
    #[test]
    fn second_row_and_hypervalent_guards_store_zero_and_preserve_aromatic_reset() {
        let mut specs = vec![AtomSpec::new(Element::P)];
        specs.extend((0..6).map(|_| AtomSpec::new(Element::C).with_aromatic(true)));
        let mut edges = vec![(0, 1, BondOrder::Single)];
        edges.extend((0..6).map(|i| (1 + i, 1 + (i + 1) % 6, BondOrder::Aromatic)));
        let mut g = graph(specs, &edges);
        for b in &mut g.bonds[1..] {
            b.set_aromatic(true);
        }
        check(
            g,
            assignment(&[1, 4, 3, 3, 3, 3, 3], &[2, 0, 1, 1, 1, 1, 1]),
            &[
                (0, false),
                (3, true),
                (3, true),
                (3, true),
                (3, true),
                (3, true),
                (3, true),
            ],
            &[false, true, true, true, true, true, true],
        );
        check(
            graph(
                vec![
                    AtomSpec::new(Element::P),
                    AtomSpec::new(Element::O),
                    AtomSpec::new(Element::O),
                    AtomSpec::new(Element::O),
                    AtomSpec::new(Element::O),
                ],
                &[
                    (0, 1, BondOrder::Double),
                    (0, 2, BondOrder::Single),
                    (0, 3, BondOrder::Single),
                    (0, 4, BondOrder::Single),
                ],
            ),
            assignment(&[5, 2, 1, 1, 1], &[0, 0, 1, 1, 1]),
            &[(0, false), (1, true), (2, true), (2, true), (2, true)],
            &[false; 4],
        );
    }
    #[test]
    fn complete_cache_failure_precedes_any_marking_but_follows_all_resets() {
        // Controlled invalid carried-valence state, not a naturally sanitized molecule.
        // Source transform evaluates the late isolated row before calling mark.
        let g = graph(
            vec![
                AtomSpec::new(Element::N),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::O),
                AtomSpec::new(Element::C),
            ],
            &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Double)],
        );
        let before = g.clone();
        let mut writes = Vec::new();
        let result = assign_conjugation_with_writer(
            &g,
            &assignment(&[1, 3, 2, i32::MAX], &[2, 1, 0, 1]),
            &mut |i, flag| writes.push((i, flag)),
        );
        assert!(matches!(
            result,
            Err(ConjugationError::IntegerOverflow {
                field: "total valence"
            })
        ));
        assert_eq!(writes, [(0, false), (1, false)]);
        assert_eq!(g, before);
    }
    #[test]
    fn empty_input_and_repeated_assignment_preserve_original_fields() {
        check(TopologyBlock::default(), assignment(&[], &[]), &[], &[]);
        let mut g = graph(
            vec![AtomSpec::new(Element::C); 2],
            &[(0, 1, BondOrder::Single)],
        );
        g.bonds[0].set_conjugated(true);
        let v = assignment(&[1, 1], &[3, 3]);
        let first = assign_conjugation(&g, &v).unwrap();
        let second = assign_conjugation(&first, &v).unwrap();
        assert_eq!(first, second);
        assert!(!first.bonds[0].is_conjugated());
        assert!(g.bonds[0].is_conjugated());
    }
}
