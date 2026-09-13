//! RDKit-aligned conjugation assignment over detached topology values.

use cosmolkit_model::{Atom, AtomId, BondId, TopologyBlock, TopologyValidationError};

use crate::{
    AromaticityError, ValenceAssignment, ValenceError, bond_valence_contrib,
    periodic_table_outer_electrons, required_valence_list,
};

#[derive(Clone, Debug, PartialEq, thiserror::Error)]
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

fn mark_conjugated_atom_bonds(
    source: &TopologyBlock,
    valence: &ValenceAssignment,
    working: &mut TopologyBlock,
    atom_id: AtomId,
) -> Result<(), ConjugationError> {
    // BEGIN RDKIT CPP FUNCTION markConjAtomBonds
    // RDKit✔️✔️: void markConjAtomBonds(Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   if (!isAtomConjugCand(at)) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto &mol = at->getOwningMol();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int atx = at->getIdx();
    // RDKit✔️✔️:   // make sure that have either 2 or 3 substitutions on this atom
    // RDKit✔️✔️:   int sbo = at->getDegree() + at->getTotalNumHs();
    // RDKit✔️✔️:   if ((sbo < 2) || (sbo > 3)) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto bnd1 : mol.atomBonds(at)) {
    // RDKit✔️✔️:     if (bnd1->getValenceContrib(at) < 1.5 ||
    // RDKit✔️✔️:         !isAtomConjugCand(bnd1->getOtherAtom(at))) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto bnd2 : mol.atomBonds(at)) {
    // RDKit✔️✔️:       if (bnd1 == bnd2) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       auto at2 = mol.getAtomWithIdx(bnd2->getOtherAtomIdx(atx));
    // RDKit✔️✔️:       sbo = at2->getDegree() + at2->getTotalNumHs();
    // RDKit✔️✔️:       if (sbo > 3) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (isAtomConjugCand(at2)) {
    // RDKit✔️✔️:         bnd1->setIsConjugated(true);
    // RDKit✔️✔️:         bnd2->setIsConjugated(true);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION markConjAtomBonds
    if !is_atom_conjugation_candidate(source, valence, atom_id)? {
        return Ok(());
    }
    if !(2..=3).contains(&total_substitutions(source, valence, atom_id)?) {
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
            || !is_atom_conjugation_candidate(source, valence, AtomId::new(first.atom_index))?
        {
            continue;
        }
        for second in incident {
            if first.bond == second.bond {
                continue;
            }
            let second_atom = AtomId::new(second.atom_index);
            if total_substitutions(source, valence, second_atom)? > 3 {
                continue;
            }
            if is_atom_conjugation_candidate(source, valence, second_atom)? {
                working.bonds[first.bond.index()].set_conjugated(true);
                working.bonds[second.bond.index()].set_conjugated(true);
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

pub fn assign_conjugation(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<TopologyBlock, ConjugationError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::setConjugation
    // RDKit✔️❌: void setConjugation(ROMol &mol) {
    // RDKit✔️❌:   // start with all bonds being marked unconjugated
    // RDKit✔️❌:   // except for aromatic bonds
    // RDKit✔️❌:   for (auto bond : mol.bonds()) {
    // RDKit✔️❌:     bond->setIsConjugated(bond->getIsAromatic());
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // loop over each atom and check if the bonds connecting to it can
    // RDKit✔️❌:   // be conjugated
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     markConjAtomBonds(atom);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::setConjugation
    // The detached atomic API clones the complete topology and validates both
    // boundaries; source RDKit mutates in place, so allocation cost is known
    // to be higher even though traversal complexity remains linear.
    topology.validate()?;
    validate_valence(topology, valence)?;
    let mut working = topology.clone();
    for bond in &mut working.bonds {
        bond.set_conjugated(bond.is_aromatic());
    }
    for atom_index in 0..topology.atoms.len() {
        mark_conjugated_atom_bonds(topology, valence, &mut working, AtomId::new(atom_index))?;
    }
    working.validate()?;
    Ok(working)
}
