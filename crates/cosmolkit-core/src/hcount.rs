//! RDKit-aligned hydrogen-count composition over detached topology values.

use crate::stereo_graph::{StereoAtomAccess, StereoGraphAccess};
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
        // Source reads only the old implicit field at this point. Use the
        // existing signed-width getter and its noImplicit branch, in atom
        // order; an old explicit sentinel is overwritten without being read.
        let original_implicit_valence =
            implicit_hydrogen_count(&working.atoms[atom_row], original_valence)? as i32;
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
            // RDKit✔️✔️: void setNumExplicitHs(unsigned int what) { d_numExplicitHs = what; }
            // RDKit✔️✔️: std::uint8_t d_numExplicitHs;
            // Both implicit fields have signed int8 source width, so this
            // promoted sum fits i32. The setter stores its low eight bits,
            // including wraparound; source defines no overflow rejection.
            // O(1) arithmetic, no allocation or fallback.
            let adjusted = (i32::from(original_explicit_hydrogens)
                + (original_implicit_valence - recalculated_implicit_valence))
                as u8;
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

fn explicit_hydrogen_count<A: StereoAtomAccess>(atom: &A) -> u32 {
    // BEGIN RDKIT CPP FUNCTION Atom::getNumExplicitHs
    // RDKit✔️✔️: unsigned int getNumExplicitHs() const { return d_numExplicitHs; }
    // END RDKIT CPP FUNCTION Atom::getNumExplicitHs
    u32::from(atom.explicit_hydrogens())
}

pub(crate) fn implicit_hydrogen_count<A: StereoAtomAccess>(
    atom: &A,
    valence: &ValenceAssignment,
) -> Result<u32, ValenceError> {
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
    // RDKit✔️✔️: unsigned int Atom::getValence(ValenceType which) const {
    // RDKit✔️✔️:   if (!dp_mol) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::EXPLICIT || df_noImplicit ||
    // RDKit✔️✔️:        d_implicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
    // RDKit✔️✔️:   if (which == ValenceType::EXPLICIT) {
    // RDKit✔️✔️:     return d_explicitValence;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return df_noImplicit ? 0 : d_implicitValence;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Source getValence reach is IMPLICIT on a topology-owned atom: the
    // explicit branch and no-owning-molecule branch are unreachable here.
    // The noImplicit return precedes cache access; negative/missing actual
    // cache rows propagate the source precondition, never the legacy H flag.
    // Complexity review: one flag, one indexed cache read, constant stack work.
    if atom.no_implicit() {
        return Ok(0);
    }
    // RDKit✔️✔️: std::int8_t d_implicitValence, d_explicitValence;
    // Caller-supplied assignments can carry calculation-width i32 values; the source cached getter
    // observes the signed eight-bit field before its initialization check.
    // One indexed read and cast, with no allocation or repeated graph scan.
    let implicit = valence
        .implicit_hydrogens
        .get(atom.id().index())
        .copied()
        .map(|value| i32::from(value as i8))
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
pub fn total_hydrogen_count_from_validated<G: StereoGraphAccess>(
    topology: &G,
    valence: &ValenceAssignment,
    atom_id: AtomId,
    include_neighbors: bool,
) -> Result<u32, ValenceError> {
    let atom = topology
        .atoms()
        .get(atom_id.index())
        .ok_or(ValenceError::AtomOutOfRange {
            atom: atom_id,
            atom_count: topology.atoms().len(),
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
            .adjacency()
            .neighbors_of(atom_id.index())
            .iter()
            .filter(|neighbor| topology.atoms()[neighbor.atom_index].atomic_number() == 1)
            .count()
    } else {
        0
    };
    // Both source getters return unsigned int: their initial sum is uint32.
    // The native int initializer preserves the low bits on the pinned two's
    // complement target. count_if returns ptrdiff_t; its addition is performed
    // in that wider signed type before conversion back to int. The returned
    // unsigned int therefore preserves the final low 32 bits. This is not a
    // signed-int addition-overflow check or a chemical-count fallback.
    // Actual native cache is int8 and explicit H count uint8; the detached
    // count carrier can also exercise the widened arithmetic boundary.
    Ok(explicit
        .wrapping_add(implicit)
        .wrapping_add(neighbor_hydrogens as u32))
}

#[cfg(test)]
mod source646_implicit_hydrogen_tests {
    use super::*;
    use cosmolkit_model::AtomSpec;
    use cosmolkit_types::Element;
    #[test]
    fn source646_no_implicit_short_circuits_missing_or_negative_actual_cache() {
        let atom = Atom::from_spec(
            AtomId::new(2),
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_implicit_hydrogen(true),
        );
        for rows in [vec![], vec![-1, -1, -1], vec![1, 2, 127]] {
            assert_eq!(
                implicit_hydrogen_count(
                    &atom,
                    &ValenceAssignment {
                        explicit_valence: vec![],
                        implicit_hydrogens: rows
                    }
                ),
                Ok(0)
            );
        }
    }
    #[test]
    fn source646_actual_cache_count_preconditions_and_flag_independence() {
        for flag in [false, true] {
            let atom = Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C).with_implicit_hydrogen(flag),
            );
            for count in [0, 1, 7, 127] {
                assert_eq!(
                    implicit_hydrogen_count(
                        &atom,
                        &ValenceAssignment {
                            explicit_valence: vec![],
                            implicit_hydrogens: vec![count]
                        }
                    ),
                    Ok(count as u32)
                );
            }
            for rows in [vec![], vec![-1], vec![-128]] {
                assert_eq!(
                    implicit_hydrogen_count(
                        &atom,
                        &ValenceAssignment {
                            explicit_valence: vec![],
                            implicit_hydrogens: rows
                        }
                    ),
                    Err(ValenceError::ImplicitValenceCacheNotInitialized {
                        atom: AtomId::new(0)
                    })
                );
            }
        }
    }
}
