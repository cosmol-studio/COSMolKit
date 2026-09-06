//! Detached hydrogen transformations.
//!
//! The public functions in this module are the algorithm boundary that the
//! `cosmolkit` runtime will call after extracting authorized blocks. Runtime
//! bookkeeping, cache invalidation, and operation contracts deliberately do
//! not appear here.

use cosmolkit_model::{
    AdjacencyList, AtomId, AtomSpec, Bond, BondOrder, BondSpec, CoordinateBlock,
    CoordinateValidationError, Element, MoleculeProperties, TopologyBlock, TopologyValidationError,
};

/// Parameters corresponding to RDKit's `MolOps::AddHsParameters`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AddHsParams {
    pub explicit_only: bool,
    pub add_coords: bool,
    pub add_residue_info: bool,
    pub skip_queries: bool,
    pub only_on_atoms: Option<Vec<AtomId>>,
}

/// Parameters corresponding to RDKit's `MolOps::RemoveHsParameters`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemoveHsParams {
    pub remove_degree_zero: bool,
    pub remove_higher_degrees: bool,
    pub remove_only_h_neighbors: bool,
    pub remove_isotopes: bool,
    pub remove_and_track_isotopes: bool,
    pub remove_dummy_neighbors: bool,
    pub remove_defining_bond_stereo: bool,
    pub remove_with_wedged_bond: bool,
    pub remove_with_query: bool,
    pub remove_mapped: bool,
    pub remove_in_sgroups: bool,
    pub show_warnings: bool,
    pub remove_nonimplicit: bool,
    pub update_explicit_count: bool,
    pub remove_hydrides: bool,
    pub remove_nontetrahedral_neighbors: bool,
}

impl Default for RemoveHsParams {
    fn default() -> Self {
        Self {
            remove_degree_zero: false,
            remove_higher_degrees: false,
            remove_only_h_neighbors: false,
            remove_isotopes: false,
            remove_and_track_isotopes: false,
            remove_dummy_neighbors: false,
            remove_defining_bond_stereo: false,
            remove_with_wedged_bond: true,
            remove_with_query: false,
            remove_mapped: true,
            remove_in_sgroups: true,
            show_warnings: true,
            remove_nonimplicit: true,
            update_explicit_count: false,
            remove_hydrides: false,
            remove_nontetrahedral_neighbors: false,
        }
    }
}

/// Errors raised by detached core algorithms.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CoreOperationError {
    InvalidTopology(TopologyValidationError),
    InvalidCoordinates(CoordinateValidationError),
    Unsupported {
        operation: &'static str,
        reason: &'static str,
    },
    InvalidInput {
        operation: &'static str,
        reason: &'static str,
    },
}

impl std::fmt::Display for CoreOperationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidTopology(error) => write!(formatter, "invalid detached topology: {error}"),
            Self::InvalidCoordinates(error) => {
                write!(formatter, "invalid detached coordinates: {error}")
            }
            Self::Unsupported { operation, reason } => {
                write!(
                    formatter,
                    "{operation} is unsupported for this detached input: {reason}"
                )
            }
            Self::InvalidInput { operation, reason } => {
                write!(formatter, "{operation} received invalid input: {reason}")
            }
        }
    }
}

impl std::error::Error for CoreOperationError {}

impl From<TopologyValidationError> for CoreOperationError {
    fn from(error: TopologyValidationError) -> Self {
        Self::InvalidTopology(error)
    }
}

impl From<CoordinateValidationError> for CoreOperationError {
    fn from(error: CoordinateValidationError) -> Self {
        Self::InvalidCoordinates(error)
    }
}

pub type DetachedBlocks = (TopologyBlock, CoordinateBlock, MoleculeProperties);

/// Apply the default AddHs operation to detached blocks.
pub fn add_hydrogens_impl(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> Result<DetachedBlocks, CoreOperationError> {
    add_hydrogens_with_params(topology, coordinates, properties, &AddHsParams::default())
}

/// Apply AddHs to detached blocks with explicit source-shaped parameters.
pub fn add_hydrogens_with_params(
    mut topology: TopologyBlock,
    mut coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    params: &AddHsParams,
) -> Result<DetachedBlocks, CoreOperationError> {
    validate_blocks(&topology, &coordinates)?;
    if params.add_residue_info {
        return Err(CoreOperationError::Unsupported {
            operation: "add_hydrogens",
            reason: "PDB residue-info assignment remains source-port work",
        });
    }
    let selected = selected_atoms(&topology, params.only_on_atoms.as_deref())?;
    if !params.explicit_only {
        return Err(CoreOperationError::Unsupported {
            operation: "add_hydrogens",
            reason: "implicit-valence assignment has not yet moved below the runtime",
        });
    }

    let additions = topology
        .atoms
        .iter()
        .filter(|atom| selected[atom.id().index()])
        .flat_map(|atom| std::iter::repeat_n(atom.id(), usize::from(atom.explicit_hydrogens())))
        .collect::<Vec<_>>();
    for atom in &mut topology.atoms {
        if selected[atom.id().index()] {
            atom.set_explicit_hydrogens(0);
        }
    }
    for heavy_atom in additions {
        // RDKit✔️❌: newIdx = mol.addAtom(new Atom(1), false, true);
        // RDKit✔️❌: mol.addBond(aidx, newIdx, Bond::SINGLE);
        // The detached model uses owned rows and rebuilds adjacency after the
        // source-shaped append loop; coordinate/property installation remains
        // the caller's responsibility.
        let hydrogen_id = AtomId::new(topology.atoms.len());
        topology.atoms.push(cosmolkit_model::Atom::from_spec(
            hydrogen_id,
            AtomSpec::new(Element::H),
        ));
        let bond_id = cosmolkit_model::BondId::new(topology.bonds.len());
        topology.bonds.push(Bond::from_spec(
            bond_id,
            BondSpec::new(heavy_atom, hydrogen_id, BondOrder::Single),
        ));
        if params.add_coords {
            for conformer in &mut coordinates.conformers_2d {
                conformer.push_coord([0.0, 0.0]);
            }
            for conformer in &mut coordinates.conformers_3d {
                conformer.push_coord([0.0, 0.0, 0.0]);
            }
        }
    }
    topology.adjacency = AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
    validate_blocks(&topology, &coordinates)?;
    Ok((topology, coordinates, properties))
}

/// Apply the default RemoveHs operation to detached blocks.
pub fn remove_hydrogens_impl(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> Result<DetachedBlocks, CoreOperationError> {
    remove_hydrogens_with_params(
        topology,
        coordinates,
        properties,
        &RemoveHsParams::default(),
    )
}

/// Apply the detached RemoveHs candidate selection and topology compaction.
pub fn remove_hydrogens_with_params(
    mut topology: TopologyBlock,
    mut coordinates: CoordinateBlock,
    mut properties: MoleculeProperties,
    params: &RemoveHsParams,
) -> Result<DetachedBlocks, CoreOperationError> {
    validate_blocks(&topology, &coordinates)?;
    if params.remove_and_track_isotopes || params.remove_isotopes || params.remove_with_query {
        return Err(CoreOperationError::Unsupported {
            operation: "remove_hydrogens",
            reason: "isotope and query-state transitions remain source-port work",
        });
    }
    let atoms_to_remove = topology
        .atoms
        .iter()
        .filter_map(|atom| {
            // RDKit✔️❌: if (shouldRemoveH(mol, atom, ps)) {
            // RDKit✔️❌:   atomsToRemove.set(atom->getIdx());
            // RDKit✔️❌: }
            (atom.atomic_number() == 1 && removable_hydrogen(&topology, atom.id(), params))
                .then_some(atom.id())
        })
        .collect::<Vec<_>>();
    if !atoms_to_remove.is_empty() {
        // RDKit✔️❌: mol.removeAtom(atom, clearProps);
        // The model mapping performs the same detached compaction boundary;
        // stereo/cache side effects are intentionally still unsupported here.
        let mapping = topology.remove_atoms_with_mapping(&atoms_to_remove);
        coordinates.remap_topology(&mapping.retained_atom_indices());
        properties.remap_topology(&mapping.atoms.new_to_old, &mapping.bonds.new_to_old);
    }
    validate_blocks(&topology, &coordinates)?;
    Ok((topology, coordinates, properties))
}

fn validate_blocks(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
) -> Result<(), CoreOperationError> {
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    Ok(())
}

fn selected_atoms(
    topology: &TopologyBlock,
    only_on_atoms: Option<&[AtomId]>,
) -> Result<Vec<bool>, CoreOperationError> {
    let mut selected = vec![only_on_atoms.is_none(); topology.atoms.len()];
    if let Some(ids) = only_on_atoms {
        for id in ids {
            let Some(slot) = selected.get_mut(id.index()) else {
                return Err(CoreOperationError::InvalidInput {
                    operation: "add_hydrogens",
                    reason: "only_on_atoms contains an out-of-range atom id",
                });
            };
            *slot = true;
        }
    }
    Ok(selected)
}

fn removable_hydrogen(topology: &TopologyBlock, atom: AtomId, params: &RemoveHsParams) -> bool {
    let Some(atom_data) = topology.atoms.get(atom.index()) else {
        return false;
    };
    if atom_data.isotope().is_some() && !params.remove_isotopes {
        return false;
    }
    if atom_data.atom_map().is_some() && !params.remove_mapped {
        return false;
    }
    if !params.remove_nonimplicit && !atom_data.implicit_hydrogen() {
        return false;
    }
    let neighbors = topology.adjacency.neighbors_of(atom.index());
    if neighbors.is_empty() {
        return params.remove_degree_zero;
    }
    if neighbors.len() > 1 && !params.remove_higher_degrees {
        return false;
    }
    if params.remove_only_h_neighbors
        && neighbors.iter().any(|neighbor| {
            topology
                .atoms
                .get(neighbor.atom_index)
                .is_some_and(|neighbor_atom| neighbor_atom.atomic_number() != 1)
        })
    {
        return false;
    }
    if params.remove_dummy_neighbors
        && neighbors.iter().any(|neighbor| {
            topology
                .atoms
                .get(neighbor.atom_index)
                .is_some_and(|neighbor_atom| neighbor_atom.atomic_number() == 0)
        })
    {
        return false;
    }
    if params.remove_in_sgroups
        && topology.substance_groups.iter().any(|group| {
            group.atoms().contains(&atom)
                || group.parent_atoms().contains(&atom)
                || group
                    .attach_points()
                    .iter()
                    .any(|point| point.atom == atom || point.leaving_atom == Some(atom))
        })
    {
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn explicit_h_topology() -> TopologyBlock {
        let carbon = cosmolkit_model::Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_explicit_hydrogens(1),
        );
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms: vec![carbon],
            ..TopologyBlock::default()
        }
    }

    fn bonded_ch_topology() -> TopologyBlock {
        let atoms = vec![
            cosmolkit_model::Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            cosmolkit_model::Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::H)),
        ];
        let bonds = vec![Bond::from_spec(
            cosmolkit_model::BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    #[test]
    fn default_add_hydrogens_reports_unmigrated_valence_boundary() {
        let error = add_hydrogens_impl(
            TopologyBlock::default(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .expect_err("default AddHs still requires the valence port");
        assert!(matches!(error, CoreOperationError::Unsupported { .. }));
    }

    #[test]
    fn explicit_add_hydrogens_is_a_detached_topology_transform() {
        let (topology, coordinates, _) = add_hydrogens_with_params(
            explicit_h_topology(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &AddHsParams {
                explicit_only: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(topology.atoms.len(), 2);
        assert_eq!(topology.bonds.len(), 1);
        assert!(coordinates.conformers_2d.is_empty());
        assert_eq!(topology.atoms[0].explicit_hydrogens(), 0);
        assert!(topology.validate().is_ok());
    }

    #[test]
    fn remove_hydrogens_compacts_detached_blocks_and_remaps_coordinates() {
        let (topology, coordinates, _) = remove_hydrogens_impl(
            bonded_ch_topology(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .unwrap();
        assert_eq!(topology.atoms.len(), 1);
        assert!(topology.bonds.is_empty());
        assert!(coordinates.conformers_3d.is_empty());
        assert!(topology.validate().is_ok());
    }
}
