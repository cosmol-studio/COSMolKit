//! Atom-count descriptor kernels (`RDKit::Descriptors::calcNumHeavyAtoms`,
//! `RDKit::Descriptors::calcNumAtoms`, `RDKit::ROMol::getNumHeavyAtoms`,
//! `RDKit::ROMol::getNumAtoms`).

use crate::{
    DescriptorError, DescriptorInput, DescriptorResult, prepared_valence, validate_topology,
};
use cosmolkit_core::{ValenceAssignment, total_hydrogen_count_from_validated};
use cosmolkit_model::TopologyBlock;

/// Single-pass heavy-atom count over explicit atom rows.
///
/// Behavior owner for `RDKit::Descriptors::calcNumHeavyAtoms`
/// (MolDescriptors.cpp:28-31) and `RDKit::ROMol::getNumHeavyAtoms`
/// (ROMol.cpp:187-196).
pub(crate) fn num_heavy_atoms_kernel(topology: &TopologyBlock) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumHeavyAtoms
    // RDKit✔️✔️: unsigned int calcNumHeavyAtoms(const ROMol &mol) {
    // RDKit✔️✔️:   return mol.getNumHeavyAtoms();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumHeavyAtoms
    // BEGIN RDKIT CPP FUNCTION: RDKit::ROMol::getNumHeavyAtoms
    // RDKit✔️✔️: unsigned int ROMol::getNumHeavyAtoms() const {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto atom : atoms()) {
    // RDKit✔️✔️:     if (atom->getAtomicNum() > 1) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::ROMol::getNumHeavyAtoms
    validate_topology(topology, "num_heavy_atoms")?;
    u32::try_from(
        topology
            .atoms
            .iter()
            .filter(|atom| atom.atomic_number() > 1)
            .count(),
    )
    .map_err(|_| DescriptorError::Unsupported {
        function: "num_heavy_atoms",
        detail: "heavy-atom count exceeds the RDKit unsigned result model".to_owned(),
    })
}

/// Prepared form of [`crate::num_heavy_atoms`].
///
/// Reads only the topology accessor: the source closure
/// (ROMol.cpp:187-196) touches no hydrogen-count, valence, or ring state,
/// so the remaining prepared inputs are intentionally unused.
pub fn num_heavy_atoms_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_heavy_atoms_kernel(input.topology())
}

/// Explicit atom rows plus attached implicit/explicit-property hydrogens.
///
/// Behavior owner for `RDKit::Descriptors::calcNumAtoms`
/// (MolDescriptors.cpp:33-37) and `RDKit::ROMol::getNumAtoms`
/// (ROMol.cpp:176-186). The borrowed assignment is shape-checked, never
/// recomputed. Source `getTotalNumHs()` keeps its default
/// `includeNeighbors=false`: explicit H-atom rows count once as rows and
/// are never summed again through neighbor state.
///
/// Complexity review: one O(n) row pass with O(1) hydrogen-row lookup per
/// atom through the pinned core owner; zero allocation beyond the result.
pub(crate) fn num_atoms_kernel(
    topology: &TopologyBlock,
    assignment: &ValenceAssignment,
) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAtoms
    // RDKit✔️✔️: unsigned int calcNumAtoms(const ROMol &mol) {
    // RDKit✔️✔️:   bool onlyExplicit = false;
    // RDKit✔️✔️:   return mol.getNumAtoms(onlyExplicit);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAtoms
    // BEGIN RDKIT CPP FUNCTION: RDKit::ROMol::getNumAtoms
    // RDKit✔️✔️: unsigned int ROMol::getNumAtoms(bool onlyExplicit) const {
    // RDKit✔️✔️:   int res = rdcast<int>(boost::num_vertices(d_graph));
    // RDKit✔️✔️:   if (!onlyExplicit) {
    // RDKit✔️✔️:     // if we are interested in hydrogens as well add them up from
    // RDKit✔️✔️:     // each
    // RDKit✔️✔️:     for (const auto atom : atoms()) {
    // RDKit✔️✔️:       res += atom->getTotalNumHs();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::ROMol::getNumAtoms
    validate_topology(topology, "num_atoms")?;
    let assignment = prepared_valence(topology, Some(assignment), "num_atoms")?;
    let mut result =
        u32::try_from(topology.atoms.len()).map_err(|_| DescriptorError::Unsupported {
            function: "num_atoms",
            detail: "explicit atom count exceeds the RDKit unsigned result model".to_owned(),
        })?;
    for atom in &topology.atoms {
        let hydrogens =
            total_hydrogen_count_from_validated(topology, assignment.as_ref(), atom.id(), false)
                .map_err(|source| DescriptorError::Valence {
                    function: "num_atoms",
                    source,
                })?;
        result = result
            .checked_add(hydrogens)
            .ok_or_else(|| DescriptorError::Unsupported {
                function: "num_atoms",
                detail: "total atom count exceeds the RDKit unsigned result model".to_owned(),
            })?;
    }
    Ok(result)
}

/// Prepared form of [`crate::num_atoms`].
///
/// Borrows the prepared valence rows; never reassigns valence.
pub fn num_atoms_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_atoms_kernel(input.topology(), input.valence())
}
