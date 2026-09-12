//! Checked construction of the single runtime-owned [`Molecule`](crate::Molecule).
//!
//! The builder owns detached model blocks only. It has no live-state, cache,
//! parser, chemistry, or operation-runtime authority.

use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, Conformer3D,
    CoordinateBlock, MoleculeProperties, StereoGroup, SubstanceGroup, SubstanceGroupId,
    TopologyBlock, TopologyEditError, TopologyMapping,
};

use crate::{Molecule, OperationError};

/// Detached construction state for one new [`Molecule`].
#[derive(Clone, Debug, PartialEq, Default)]
pub struct MoleculeBuilder {
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
}

impl MoleculeBuilder {
    /// Creates an empty detached builder.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Stores detached input blocks without claiming they are valid.
    #[must_use]
    pub fn from_parts(
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Self {
        Self {
            topology,
            coordinates,
            properties,
        }
    }

    /// Validates every detached block before publishing one live value.
    pub fn build(self) -> Result<Molecule, OperationError> {
        Molecule::from_validated_parts(self.topology, self.coordinates, self.properties)
    }

    /// Appends an isolated atom with its canonical row id.
    pub fn add_atom(&mut self, atom: AtomSpec) -> AtomId {
        let old_atom_count = self.topology.atoms.len();
        let old_bond_count = self.topology.bonds.len();
        let id = AtomId::new(old_atom_count);
        let atom = Atom::from_spec(id, atom);
        let mut atoms = self.topology.atoms.clone();
        atoms.push(atom.clone());
        match TopologyBlock::try_from_parts(
            atoms,
            self.topology.bonds.clone(),
            self.topology.substance_groups.clone(),
            self.topology.stereo_groups.clone(),
        ) {
            Ok(topology) => {
                let mapping = TopologyMapping::with_appended(old_atom_count, old_bond_count, 1, 0);
                self.properties
                    .remap_topology(mapping.atoms().new_to_old(), mapping.bonds().new_to_old());
                self.topology = topology;
            }
            Err(_) => {
                // `from_parts` deliberately accepts invalid detached input and
                // `add_atom` has an infallible canonical signature. Preserve
                // that invalid state for the structured final `build` error.
                self.topology.atoms.push(atom);
            }
        }
        id
    }

    /// Appends a checked bond with its canonical row id.
    pub fn add_bond(&mut self, bond: BondSpec) -> Result<BondId, OperationError> {
        let mut edit = self
            .topology
            .begin_batch_edit()
            .map_err(OperationError::InvalidTopologyEdit)?;
        let id = edit
            .add_bond(bond)
            .map_err(OperationError::InvalidTopologyEdit)?;
        let (topology, mapping) = edit.finish().map_err(OperationError::InvalidTopologyEdit)?;
        self.properties
            .remap_topology(mapping.atoms().new_to_old(), mapping.bonds().new_to_old());
        self.topology = topology;
        Ok(id)
    }

    /// Changes one detached atom's formal charge without exposing storage.
    pub fn set_atom_formal_charge(
        &mut self,
        atom_id: AtomId,
        formal_charge: i8,
    ) -> Result<(), OperationError> {
        self.validate_topology()?;
        let atom_count = self.topology.atoms.len();
        let Some(atom) = self.topology.atoms.get(atom_id.index()) else {
            return Err(OperationError::InvalidTopologyEdit(
                TopologyEditError::AtomOutOfRange {
                    atom: atom_id,
                    atom_count,
                },
            ));
        };
        let mut atoms = self.topology.atoms.clone();
        atoms[atom_id.index()] = atom.clone();
        atoms[atom_id.index()].set_formal_charge(formal_charge);
        self.topology = TopologyBlock::try_from_parts(
            atoms,
            self.topology.bonds.clone(),
            self.topology.substance_groups.clone(),
            self.topology.stereo_groups.clone(),
        )
        .map_err(OperationError::InvalidTopology)?;
        Ok(())
    }

    /// Changes one detached bond order without changing row identity.
    pub fn set_bond_order(
        &mut self,
        bond_id: BondId,
        order: BondOrder,
    ) -> Result<(), OperationError> {
        self.validate_topology()?;
        let bond_count = self.topology.bonds.len();
        let Some(bond) = self.topology.bonds.get(bond_id.index()) else {
            return Err(OperationError::InvalidTopologyEdit(
                TopologyEditError::BondOutOfRange {
                    bond: bond_id,
                    bond_count,
                },
            ));
        };
        let mut bonds = self.topology.bonds.clone();
        bonds[bond_id.index()] = bond.clone();
        bonds[bond_id.index()].set_order(order);
        self.topology = TopologyBlock::try_from_parts(
            self.topology.atoms.clone(),
            bonds,
            self.topology.substance_groups.clone(),
            self.topology.stereo_groups.clone(),
        )
        .map_err(OperationError::InvalidTopology)?;
        Ok(())
    }

    /// Removes the bond between two valid atoms, if one exists.
    pub fn remove_bond_between_atoms(
        &mut self,
        begin_atom: AtomId,
        end_atom: AtomId,
    ) -> Result<bool, OperationError> {
        self.validate_atom_id(begin_atom)?;
        self.validate_atom_id(end_atom)?;
        let Some(bond_id) = self
            .topology
            .bonds
            .iter()
            .find(|bond| {
                (bond.begin() == begin_atom && bond.end() == end_atom)
                    || (bond.begin() == end_atom && bond.end() == begin_atom)
            })
            .map(Bond::id)
        else {
            return Ok(false);
        };

        let mut edit = self
            .topology
            .begin_batch_edit()
            .map_err(OperationError::InvalidTopologyEdit)?;
        edit.remove_bond(bond_id)
            .map_err(OperationError::InvalidTopologyEdit)?;
        let (topology, mapping) = edit.finish().map_err(OperationError::InvalidTopologyEdit)?;
        self.properties
            .remap_topology(mapping.atoms().new_to_old(), mapping.bonds().new_to_old());
        self.topology = topology;
        Ok(true)
    }

    /// Returns the number of incident bonds for a valid atom.
    pub fn degree(&self, atom_id: AtomId) -> Result<usize, OperationError> {
        self.validate_atom_id(atom_id)?;
        Ok(self.topology.adjacency.neighbors_of(atom_id.index()).len())
    }

    /// Returns incident bond ids in topology input order.
    pub fn neighbor_bonds(&self, atom_id: AtomId) -> Result<Vec<BondId>, OperationError> {
        self.validate_atom_id(atom_id)?;
        Ok(self
            .topology
            .adjacency
            .neighbors_of(atom_id.index())
            .iter()
            .map(|neighbor| neighbor.bond)
            .collect())
    }

    /// Returns the bond joining two valid atoms, if present.
    pub fn bond_between_atoms(
        &self,
        begin_atom: AtomId,
        end_atom: AtomId,
    ) -> Result<Option<BondId>, OperationError> {
        self.validate_atom_id(begin_atom)?;
        self.validate_atom_id(end_atom)?;
        Ok(self
            .topology
            .adjacency
            .neighbors_of(begin_atom.index())
            .iter()
            .find(|neighbor| neighbor.atom_index == end_atom.index())
            .map(|neighbor| neighbor.bond))
    }

    #[must_use]
    pub fn atoms(&self) -> &[Atom] {
        &self.topology.atoms
    }

    #[must_use]
    pub fn bonds(&self) -> &[Bond] {
        &self.topology.bonds
    }

    #[must_use]
    pub fn substance_groups(&self) -> &[SubstanceGroup] {
        &self.topology.substance_groups
    }

    #[must_use]
    pub fn stereo_groups(&self) -> &[StereoGroup] {
        &self.topology.stereo_groups
    }

    #[must_use]
    pub const fn coordinates(&self) -> &CoordinateBlock {
        &self.coordinates
    }

    #[must_use]
    pub const fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }

    /// Replaces all 2D conformers with canonical conformer zero.
    pub fn set_2d_coordinates(&mut self, coordinates: Vec<[f64; 2]>) -> Result<(), OperationError> {
        let conformer = Conformer2D::new(0, coordinates);
        conformer
            .validate_for_atom_count(self.topology.atoms.len())
            .map_err(OperationError::InvalidCoordinates)?;
        self.coordinates.conformers_2d = vec![conformer];
        Ok(())
    }

    /// Adds one checked 2D conformer and returns its dimension-local id.
    pub fn add_2d_conformer(
        &mut self,
        coordinates: Vec<[f64; 2]>,
    ) -> Result<usize, OperationError> {
        self.add_conformer(DetachedConformer::TwoD(coordinates))
    }

    /// Adds one checked 3D conformer and returns its dimension-local id.
    pub fn add_3d_conformer(
        &mut self,
        coordinates: Vec<[f64; 3]>,
    ) -> Result<usize, OperationError> {
        self.add_conformer(DetachedConformer::ThreeD(coordinates))
    }

    /// Appends an SGroup after assigning its canonical row id.
    pub fn add_substance_group(
        &mut self,
        mut group: SubstanceGroup,
    ) -> Result<SubstanceGroupId, OperationError> {
        self.validate_topology()?;
        let id = SubstanceGroupId::new(self.topology.substance_groups.len());
        group.set_id(id);
        let mut groups = self.topology.substance_groups.clone();
        groups.push(group);
        self.topology = TopologyBlock::try_from_parts(
            self.topology.atoms.clone(),
            self.topology.bonds.clone(),
            groups,
            self.topology.stereo_groups.clone(),
        )
        .map_err(OperationError::InvalidTopology)?;
        Ok(id)
    }

    /// Appends a checked enhanced-stereo group.
    pub fn add_stereo_group(&mut self, group: StereoGroup) -> Result<usize, OperationError> {
        self.validate_topology()?;
        let id = self.topology.stereo_groups.len();
        let mut groups = self.topology.stereo_groups.clone();
        groups.push(group);
        self.topology = TopologyBlock::try_from_parts(
            self.topology.atoms.clone(),
            self.topology.bonds.clone(),
            self.topology.substance_groups.clone(),
            groups,
        )
        .map_err(OperationError::InvalidTopology)?;
        Ok(id)
    }

    #[must_use]
    pub fn with_name(mut self, name: String) -> Self {
        self.properties = self.properties.with_name(name);
        self
    }

    pub fn with_property(mut self, key: String, value: String) -> Result<Self, OperationError> {
        self.properties
            .set_prop(key, value)
            .map_err(OperationError::InvalidProperty)?;
        Ok(self)
    }

    #[must_use]
    pub fn with_sdf_data_field(mut self, key: String, value: String) -> Self {
        self.properties = self.properties.with_sdf_data_field(key, value);
        self
    }

    #[must_use]
    pub fn with_properties(mut self, properties: MoleculeProperties) -> Self {
        self.properties = properties;
        self
    }

    fn validate_topology(&self) -> Result<(), OperationError> {
        self.topology
            .validate()
            .map_err(OperationError::InvalidTopology)
    }

    fn validate_atom_id(&self, atom: AtomId) -> Result<(), OperationError> {
        self.validate_topology()?;
        if atom.index() >= self.topology.atoms.len() {
            return Err(OperationError::InvalidTopologyEdit(
                TopologyEditError::AtomOutOfRange {
                    atom,
                    atom_count: self.topology.atoms.len(),
                },
            ));
        }
        Ok(())
    }

    fn add_conformer(&mut self, input: DetachedConformer) -> Result<usize, OperationError> {
        match input {
            DetachedConformer::TwoD(coordinates) => {
                let id =
                    next_conformer_id(self.coordinates.conformers_2d.iter().map(Conformer2D::id));
                let conformer = Conformer2D::new(id, coordinates);
                conformer
                    .validate_for_atom_count(self.topology.atoms.len())
                    .map_err(OperationError::InvalidCoordinates)?;
                self.coordinates.conformers_2d.push(conformer);
                Ok(id)
            }
            DetachedConformer::ThreeD(coordinates) => {
                let id =
                    next_conformer_id(self.coordinates.conformers_3d.iter().map(Conformer3D::id));
                let conformer = Conformer3D::new(id, coordinates, true);
                conformer
                    .validate_for_atom_count(self.topology.atoms.len())
                    .map_err(OperationError::InvalidCoordinates)?;
                self.coordinates.conformers_3d.push(conformer);
                Ok(id)
            }
        }
    }
}

enum DetachedConformer {
    TwoD(Vec<[f64; 2]>),
    ThreeD(Vec<[f64; 3]>),
}

fn next_conformer_id(ids: impl Iterator<Item = usize> + Clone) -> usize {
    let count = ids.clone().count();
    (0..=count)
        .find(|candidate| !ids.clone().any(|id| id == *candidate))
        .expect("a finite conformer table always has a free id")
}
