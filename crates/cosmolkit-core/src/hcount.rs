//! RDKit-aligned hydrogen-count composition over detached topology values.

use cosmolkit_model::{Atom, AtomId, TopologyBlock, TopologyValidationError};

use crate::{
    ValenceAssignment, ValenceError, assign_explicit_valence_for_atom_from_parts,
    assign_implicit_valence_for_atom_from_parts_with_explicit_valence,
};

/// Detached result of RDKit's post-aromaticity hydrogen adjustment.
#[derive(Debug, Clone, PartialEq)]
pub struct AdjustHsAssignment {
    pub topology: TopologyBlock,
    pub valence: ValenceAssignment,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AdjustHsError {
    #[error("invalid topology: {source}")]
    InvalidTopology { source: TopologyValidationError },
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
    #[error(
        "explicit hydrogen adjustment overflows at atom {atom}: explicit={original_explicit_hydrogens}, original implicit={original_implicit_valence}, recalculated implicit={recalculated_implicit_valence}"
    )]
    ExplicitHydrogenOverflow {
        atom: AtomId,
        original_explicit_hydrogens: u8,
        original_implicit_valence: i32,
        recalculated_implicit_valence: i32,
    },
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

fn validate_adjust_hs_valence(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<(), AdjustHsError> {
    let expected = topology.atoms.len();
    for (field, rows) in [
        ("explicit_valence", &valence.explicit_valence),
        ("implicit_hydrogens", &valence.implicit_hydrogens),
    ] {
        if rows.len() != expected {
            return Err(AdjustHsError::ValenceAssignmentLength {
                field,
                actual: rows.len(),
                expected,
            });
        }
        if let Some((row, &value)) = rows.iter().enumerate().find(|(_, value)| **value < 0) {
            return Err(AdjustHsError::InvalidValenceRow {
                atom: AtomId::new(row),
                field,
                value,
            });
        }
    }
    Ok(())
}

/// Transfer implicit hydrogens lost after aromaticity perception into the
/// explicit-hydrogen field, preserving the source atom-row order.
pub fn adjust_hs(
    topology: &TopologyBlock,
    original_valence: &ValenceAssignment,
) -> Result<AdjustHsAssignment, AdjustHsError> {
    topology
        .validate()
        .map_err(|source| AdjustHsError::InvalidTopology { source })?;
    validate_adjust_hs_valence(topology, original_valence)?;

    // BEGIN RDKIT CPP FUNCTION MolOps::adjustHs
    // RDKit✔️❌: void adjustHs(RWMol &mol) {
    // RDKit✔️❌:   //
    // RDKit✔️❌:   //  Go through and adjust the number of implicit and explicit Hs
    // RDKit✔️❌:   //  on each atom in the molecule.
    // RDKit✔️❌:   //
    // RDKit✔️❌:   //  Atoms that do not *need* explicit Hs
    // RDKit✔️❌:   //
    // RDKit✔️❌:   //  Assumptions: this is called after the molecule has been
    // RDKit✔️❌:   //  sanitized, aromaticity has been perceived, and the implicit
    // RDKit✔️❌:   //  valence of everything has been calculated.
    // RDKit✔️❌:   //
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     int origImplicitV = atom->getValence(Atom::ValenceType::IMPLICIT);
    // RDKit✔️❌:     atom->calcExplicitValence(false);
    // RDKit✔️❌:     int origExplicitV = atom->getNumExplicitHs();
    // RDKit✔️❌:
    // RDKit✔️❌:     int newImplicitV = atom->calcImplicitValence(false);
    // RDKit✔️❌:     //
    // RDKit✔️❌:     //  Case 1: The disappearing Hydrogen
    // RDKit✔️❌:     //    Smiles:  O=C1NC=CC2=C1C=CC=C2
    // RDKit✔️❌:     //
    // RDKit✔️❌:     //    after perception is done, the N atom has two aromatic
    // RDKit✔️❌:     //    bonds to it and a single implicit H.  When the Smiles is
    // RDKit✔️❌:     //    written, we get: n1ccc2ccccc2c1=O.  Here the nitrogen has
    // RDKit✔️❌:     //    no implicit Hs (because there are two aromatic bonds to
    // RDKit✔️❌:     //    it, giving it a valence of 3).  Also: this SMILES is bogus
    // RDKit✔️❌:     //    (un-kekulizable).  The correct SMILES would be:
    // RDKit✔️❌:     //    [nH]1ccc2ccccc2c1=O.  So we need to loop through the atoms
    // RDKit✔️❌:     //    and find those that have lost implicit H; we'll add those
    // RDKit✔️❌:     //    back as explicit Hs.
    // RDKit✔️❌:     //
    // RDKit✔️❌:     //    <phew> that takes way longer to comment than it does to
    // RDKit✔️❌:     //    write:
    // RDKit✔️❌:     if (newImplicitV < origImplicitV) {
    // RDKit✔️❌:       atom->setNumExplicitHs(origExplicitV + (origImplicitV - newImplicitV));
    // RDKit✔️❌:       atom->calcExplicitValence(false);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::adjustHs
    // The detached boundary must clone the full topology to preserve input
    // immutability, unlike the source's in-place mutation. Per-atom valence
    // calculations retain the source O(V + E) traversal shape, but the full
    // topology clone is a material allocation difference.

    let mut working = topology.clone();
    let mut explicit_valence = Vec::with_capacity(working.atoms.len());
    let mut implicit_hydrogens = Vec::with_capacity(working.atoms.len());

    for atom_row in 0..working.atoms.len() {
        let atom = AtomId::new(atom_row);
        let original_implicit_valence = original_valence.implicit_hydrogens[atom_row];
        let current_explicit_valence = assign_explicit_valence_for_atom_from_parts(
            &working.atoms,
            &working.bonds,
            &working.adjacency,
            atom,
            false,
        )?;
        let original_explicit_hydrogens = working.atoms[atom_row].explicit_hydrogens();
        let recalculated_implicit_valence =
            assign_implicit_valence_for_atom_from_parts_with_explicit_valence(
                &working.atoms,
                &working.bonds,
                &working.adjacency,
                atom,
                current_explicit_valence,
                false,
            )?;

        let final_explicit_valence = if recalculated_implicit_valence < original_implicit_valence {
            let lost = original_implicit_valence
                .checked_sub(recalculated_implicit_valence)
                .ok_or(AdjustHsError::ExplicitHydrogenOverflow {
                    atom,
                    original_explicit_hydrogens,
                    original_implicit_valence,
                    recalculated_implicit_valence,
                })?;
            let adjusted = i32::from(original_explicit_hydrogens)
                .checked_add(lost)
                .and_then(|value| u8::try_from(value).ok())
                .ok_or(AdjustHsError::ExplicitHydrogenOverflow {
                    atom,
                    original_explicit_hydrogens,
                    original_implicit_valence,
                    recalculated_implicit_valence,
                })?;
            working.atoms[atom_row].set_explicit_hydrogens(adjusted);
            assign_explicit_valence_for_atom_from_parts(
                &working.atoms,
                &working.bonds,
                &working.adjacency,
                atom,
                false,
            )?
        } else {
            current_explicit_valence
        };

        explicit_valence.push(final_explicit_valence);
        implicit_hydrogens.push(recalculated_implicit_valence);
    }

    working
        .validate()
        .map_err(|source| AdjustHsError::InvalidTopology { source })?;
    Ok(AdjustHsAssignment {
        topology: working,
        valence: ValenceAssignment {
            explicit_valence,
            implicit_hydrogens,
        },
    })
}

fn explicit_hydrogen_count(atom: &Atom) -> u32 {
    // BEGIN RDKIT CPP FUNCTION Atom::getNumExplicitHs
    // RDKit✔️✔️: unsigned int getNumExplicitHs() const { return d_numExplicitHs; }
    // END RDKIT CPP FUNCTION Atom::getNumExplicitHs
    u32::from(atom.explicit_hydrogens())
}

fn implicit_hydrogen_count(atom: &Atom, valence: &ValenceAssignment) -> Result<u32, ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::getNumImplicitHs
    // RDKit✔️✔️: unsigned int Atom::getNumImplicitHs() const {
    // RDKit✔️✔️:   if (df_noImplicit) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(d_implicitValence > -1,
    // RDKit✔️✔️:                "getNumImplicitHs() called without preceding call to "
    // RDKit✔️✔️:                "calcImplicitValence()");
    // RDKit✔️✔️:   return getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getNumImplicitHs
    if atom.no_implicit() {
        return Ok(0);
    }
    let implicit = valence
        .implicit_hydrogens
        .get(atom.id().index())
        .copied()
        .filter(|value| *value >= 0)
        .ok_or(ValenceError::ImplicitValenceCacheNotInitialized { atom: atom.id() })?;
    Ok(implicit as u32)
}

/// Return an atom's explicit plus implicit hydrogen count and, when requested,
/// its explicit hydrogen-atom neighbors.
pub fn total_hydrogen_count(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_id: AtomId,
    include_neighbors: bool,
) -> Result<u32, ValenceError> {
    topology
        .validate()
        .map_err(|source| ValenceError::InvalidTopology { source })?;
    total_hydrogen_count_from_validated(topology, valence, atom_id, include_neighbors)
}

/// Compose hydrogen counts after the caller has validated topology and
/// assignment dimensions. This is the unique O(degree) implementation used by
/// chemistry phases that would otherwise repeat a whole-topology scan.
pub(crate) fn total_hydrogen_count_from_validated(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_id: AtomId,
    include_neighbors: bool,
) -> Result<u32, ValenceError> {
    let atom = topology
        .atoms
        .get(atom_id.index())
        .ok_or(ValenceError::AtomOutOfRange {
            atom: atom_id,
            atom_count: topology.atoms.len(),
        })?;

    // BEGIN RDKIT CPP FUNCTION Atom::getTotalNumHs
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getTotalNumHs
    // Validation is owned by the public wrapper or the already-validating
    // chemistry phase; this inner composition retains the source O(degree)
    // shape without duplicating the hydrogen algorithm.
    let explicit = explicit_hydrogen_count(atom);
    let implicit = implicit_hydrogen_count(atom, valence)?;
    let neighbor_hydrogens = if include_neighbors {
        topology
            .adjacency
            .neighbors_of(atom_id.index())
            .iter()
            .filter(|neighbor| topology.atoms[neighbor.atom_index].atomic_number() == 1)
            .count()
    } else {
        0
    };
    let overflow = || ValenceError::HydrogenCountOverflow {
        atom: atom_id,
        explicit,
        implicit,
        neighbor_hydrogens,
    };
    let explicit = i32::try_from(explicit).map_err(|_| overflow())?;
    let implicit = i32::try_from(implicit).map_err(|_| overflow())?;
    let neighbor_count = i32::try_from(neighbor_hydrogens).map_err(|_| overflow())?;
    let total = explicit
        .checked_add(implicit)
        .and_then(|count| count.checked_add(neighbor_count))
        .ok_or_else(overflow)?;
    u32::try_from(total).map_err(|_| ValenceError::HydrogenCountOverflow {
        atom: atom_id,
        explicit: explicit as u32,
        implicit: implicit as u32,
        neighbor_hydrogens,
    })
}
