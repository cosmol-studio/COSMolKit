//! RDKit-aligned hydrogen-count composition over detached topology values.

use cosmolkit_model::{Atom, AtomId, TopologyBlock};

use crate::{ValenceAssignment, ValenceError};

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
    let atom = topology
        .atoms
        .get(atom_id.index())
        .ok_or(ValenceError::AtomOutOfRange {
            atom: atom_id,
            atom_count: topology.atoms.len(),
        })?;

    // BEGIN RDKIT CPP FUNCTION Atom::getTotalNumHs
    // RDKit✔️❌: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️❌:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️❌:   if (includeNeighbors && dp_mol) {
    // RDKit✔️❌:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️❌:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️❌:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️❌:     });
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Atom::getTotalNumHs
    // The detached public boundary validates the complete topology before
    // trusting adjacency. This preserves behavior on valid input but adds an
    // O(V + E) scan to the source's O(degree) lookup, so the second marker is
    // deliberately `❌`.
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
