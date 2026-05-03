use crate::{AdjacencyList, Atom, Bond};
use glam::{DVec2, DVec3};
use std::sync::Arc;

/// Errors returned by SMILES parsing routines.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SmilesParseError {
    /// Parsing is not implemented yet in bootstrap phase.
    #[error("SMILES parser is not implemented yet")]
    NotImplemented,
    /// Concrete parse error from the RDKit-aligned subset parser.
    #[error("{0}")]
    ParseError(String),
}

/// Errors returned by SMILES writing routines.
pub use crate::smiles_write::SmilesWriteError;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum CoordinateDimension {
    TwoD,
    ThreeD,
}

/// Persistent topology storage owned by a molecule value.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TopologyData {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub adjacency: Option<AdjacencyList>,
}

/// Coordinate storage owned separately from topology.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ConformerStore {
    pub coords_2d: Option<Vec<DVec2>>,
    pub conformers_3d: Vec<Vec<DVec3>>,
    pub source_coordinate_dim: Option<CoordinateDimension>,
}

/// Molecule-level metadata storage.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PropertyStore {
    pub name: Option<String>,
    pub sdf_data_fields: Vec<(String, String)>,
}

/// Molecular graph with split owned topology, coordinate, and metadata blocks.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Molecule {
    topology: Arc<TopologyData>,
    conformers: Arc<ConformerStore>,
    props: Arc<PropertyStore>,
}

impl Molecule {
    /// Create an empty molecule.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add one atom and return its assigned index.
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let index = self.atoms().len();
        let mut atom = atom;
        atom.index = index;
        self.atoms_mut().push(atom);
        self.clear_conformers();
        self.clear_adjacency_cache();
        index
    }

    /// Add one bond and return its assigned index.
    pub fn add_bond(&mut self, bond: Bond) -> usize {
        let index = self.bonds().len();
        let mut bond = bond;
        bond.index = index;
        self.bonds_mut().push(bond);
        self.clear_adjacency_cache();
        index
    }

    /// Rebuild adjacency representation from current topology.
    pub fn rebuild_adjacency(&mut self) {
        let adjacency = AdjacencyList::from_topology(self.atoms().len(), self.bonds());
        self.topology_mut().adjacency = Some(adjacency);
    }

    /// Construct one molecule from a SMILES string.
    pub fn from_smiles(smiles: &str) -> Result<Self, SmilesParseError> {
        crate::smiles::parse_smiles(smiles)
    }

    /// Serialize this molecule to a SMILES string.
    pub fn to_smiles(&self, isomeric_smiles: bool) -> Result<String, SmilesWriteError> {
        let params = crate::smiles_write::SmilesWriteParams {
            do_isomeric_smiles: isomeric_smiles,
            ..Default::default()
        };
        self.to_smiles_with_params(&params)
    }

    /// Serialize this molecule to a SMILES string using RDKit-like write params.
    pub fn to_smiles_with_params(
        &self,
        params: &crate::smiles_write::SmilesWriteParams,
    ) -> Result<String, SmilesWriteError> {
        crate::smiles_write::mol_to_smiles(self, params)
    }

    /// Compute RDKit-aligned 2D coordinates and store them on this molecule.
    pub fn compute_2d_coords(&mut self) -> Result<&mut Self, crate::io::molblock::MolWriteError> {
        let coords = crate::io::molblock::compute_2d_coords(self)?;
        let conformers = self.conformers_mut();
        conformers.coords_2d = Some(coords.into_iter().map(|(x, y)| DVec2::new(x, y)).collect());
        conformers.source_coordinate_dim = Some(CoordinateDimension::TwoD);
        Ok(self)
    }

    /// Return a new molecule with explicit hydrogens added.
    pub fn with_hydrogens(&self) -> Result<Self, crate::hydrogens::AddHydrogensError> {
        let mut out = self.clone();
        crate::hydrogens::add_hydrogens_in_place(&mut out)?;
        Ok(out)
    }

    /// Return a new molecule with explicit hydrogens removed.
    pub fn without_hydrogens(&self) -> Result<Self, crate::hydrogens::RemoveHydrogensError> {
        let mut out = self.clone();
        crate::hydrogens::remove_hydrogens_in_place(&mut out)?;
        Ok(out)
    }

    /// Return a new molecule with RDKit-aligned 2D coordinates.
    pub fn with_2d_coords(&self) -> Result<Self, crate::io::molblock::MolWriteError> {
        let mut out = self.clone();
        out.compute_2d_coords()?;
        Ok(out)
    }

    /// Return a new molecule with aromatic bonds converted to an explicit Kekule form.
    pub fn with_kekulized_bonds(
        &self,
        sanitize: bool,
    ) -> Result<Self, crate::kekulize::KekulizeError> {
        let mut out = self.clone();
        crate::kekulize::kekulize_in_place(&mut out, sanitize)?;
        Ok(out)
    }

    /// Return a new molecule with molecule-level name metadata replaced.
    #[must_use]
    pub fn with_name(&self, name: impl Into<String>) -> Self {
        let mut out = self.clone();
        out.props_mut().name = Some(name.into());
        out
    }

    /// Return persistent topology storage.
    #[must_use]
    pub fn topology(&self) -> &TopologyData {
        &self.topology
    }

    /// Return mutable topology storage.
    pub fn topology_mut(&mut self) -> &mut TopologyData {
        Arc::make_mut(&mut self.topology)
    }

    /// Return atom storage.
    #[must_use]
    pub fn atoms(&self) -> &[Atom] {
        &self.topology.atoms
    }

    /// Return mutable atom storage.
    pub fn atoms_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.topology_mut().atoms
    }

    /// Return bond storage.
    #[must_use]
    pub fn bonds(&self) -> &[Bond] {
        &self.topology.bonds
    }

    /// Return mutable bond storage.
    pub fn bonds_mut(&mut self) -> &mut Vec<Bond> {
        &mut self.topology_mut().bonds
    }

    /// Return cached adjacency, if present.
    #[must_use]
    pub fn adjacency(&self) -> Option<&AdjacencyList> {
        self.topology.adjacency.as_ref()
    }

    /// Return mutable cached adjacency slot.
    pub fn adjacency_mut(&mut self) -> &mut Option<AdjacencyList> {
        &mut self.topology_mut().adjacency
    }

    /// Clear cached topology-derived adjacency.
    pub fn clear_adjacency_cache(&mut self) {
        self.topology_mut().adjacency = None;
    }

    /// Return coordinate storage.
    #[must_use]
    pub fn conformers(&self) -> &ConformerStore {
        &self.conformers
    }

    /// Return mutable coordinate storage.
    pub fn conformers_mut(&mut self) -> &mut ConformerStore {
        Arc::make_mut(&mut self.conformers)
    }

    /// Return mutable 2D coordinate storage.
    pub fn coords_2d_mut(&mut self) -> &mut Option<Vec<DVec2>> {
        &mut self.conformers_mut().coords_2d
    }

    /// Replace stored 2D coordinates.
    pub fn set_coords_2d(&mut self, coords: Option<Vec<DVec2>>) {
        self.conformers_mut().coords_2d = coords;
    }

    /// Return stored 3D conformers.
    #[must_use]
    pub fn conformers_3d(&self) -> &[Vec<DVec3>] {
        &self.conformers.conformers_3d
    }

    /// Return mutable 3D conformer storage.
    pub fn conformers_3d_mut(&mut self) -> &mut Vec<Vec<DVec3>> {
        &mut self.conformers_mut().conformers_3d
    }

    /// Return the source coordinate dimensionality, if known.
    #[must_use]
    pub fn source_coordinate_dim(&self) -> Option<CoordinateDimension> {
        self.conformers.source_coordinate_dim
    }

    /// Set the source coordinate dimensionality.
    pub fn set_source_coordinate_dim(&mut self, coordinate_dim: Option<CoordinateDimension>) {
        self.conformers_mut().source_coordinate_dim = coordinate_dim;
    }

    /// Remove all stored conformer coordinates and dimensionality metadata.
    pub fn clear_conformers(&mut self) {
        let conformers = self.conformers_mut();
        conformers.coords_2d = None;
        conformers.conformers_3d.clear();
        conformers.source_coordinate_dim = None;
    }

    /// Return molecule-level metadata storage.
    #[must_use]
    pub fn props(&self) -> &PropertyStore {
        &self.props
    }

    /// Return mutable molecule-level metadata storage.
    pub fn props_mut(&mut self) -> &mut PropertyStore {
        Arc::make_mut(&mut self.props)
    }

    /// Return stored 2D coordinates, if present.
    #[must_use]
    pub fn coords_2d(&self) -> Option<&[DVec2]> {
        self.conformers.coords_2d.as_deref()
    }

    /// Return the default stored 3D conformer, if present.
    #[must_use]
    pub fn coords_3d(&self) -> Option<&[DVec3]> {
        self.conformers.conformers_3d.first().map(Vec::as_slice)
    }

    /// Return one stored 3D conformer by index.
    #[must_use]
    pub fn conformer_3d(&self, index: usize) -> Option<&[DVec3]> {
        self.conformers.conformers_3d.get(index).map(Vec::as_slice)
    }

    /// Return the number of stored 3D conformers.
    #[must_use]
    pub fn num_3d_conformers(&self) -> usize {
        self.conformers.conformers_3d.len()
    }

    /// Return atom atomic numbers in atom-index order.
    #[must_use]
    pub fn atomic_numbers(&self) -> Vec<u8> {
        self.atoms().iter().map(|atom| atom.atomic_num).collect()
    }

    /// Return the RDKit-style distance-geometry bounds matrix.
    pub fn dg_bounds_matrix(&self) -> Result<Vec<Vec<f64>>, crate::DgBoundsError> {
        crate::distgeom::dg_bounds_matrix(self)
    }

    /// Serialize this molecule to RDKit-style SVG.
    pub fn to_svg(&self, width: u32, height: u32) -> Result<String, crate::SvgDrawError> {
        crate::draw::mol_to_svg(self, width, height)
    }

    /// Rasterize this molecule's RDKit-style SVG drawing into PNG bytes.
    pub fn to_png(&self, width: u32, height: u32) -> Result<Vec<u8>, crate::SvgDrawError> {
        crate::draw::mol_to_png(self, width, height)
    }

    /// Return the RDKit `PrepareMolForDrawing()`-style prepared drawing snapshot.
    pub fn prepare_for_drawing_parity(
        &self,
    ) -> Result<crate::PreparedDrawMolecule, crate::SvgDrawError> {
        crate::draw::prepare_mol_for_drawing_parity(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondDirection, BondOrder, BondStereo, ChiralTag};

    fn carbon(index: usize) -> Atom {
        Atom {
            index,
            atomic_num: 6,
            is_aromatic: false,
            formal_charge: 0,
            explicit_hydrogens: 0,
            no_implicit: false,
            num_radical_electrons: 0,
            chiral_tag: ChiralTag::Unspecified,
            isotope: None,
            atom_map_num: None,
            props: Default::default(),
            rdkit_cip_rank: None,
        }
    }

    fn single_bond(begin_atom: usize, end_atom: usize) -> Bond {
        Bond {
            index: 0,
            begin_atom,
            end_atom,
            order: BondOrder::Single,
            is_aromatic: false,
            direction: BondDirection::None,
            stereo: BondStereo::None,
            stereo_atoms: Vec::new(),
        }
    }

    fn ethane_with_2d_coords() -> Molecule {
        let mut mol = Molecule::new();
        mol.add_atom(carbon(0));
        mol.add_atom(carbon(1));
        mol.add_bond(single_bond(0, 1));
        mol.set_coords_2d(Some(vec![DVec2::new(0.0, 0.0), DVec2::new(1.0, 0.0)]));
        mol
    }

    #[test]
    fn clone_shares_storage_blocks_until_touched_block_mutates() {
        let mol = ethane_with_2d_coords();
        let cloned = mol.clone();
        assert!(Arc::ptr_eq(&mol.topology, &cloned.topology));
        assert!(Arc::ptr_eq(&mol.conformers, &cloned.conformers));
        assert!(Arc::ptr_eq(&mol.props, &cloned.props));

        let mut topology_edit = mol.clone();
        topology_edit.atoms_mut()[0].formal_charge = 1;
        assert!(!Arc::ptr_eq(&mol.topology, &topology_edit.topology));
        assert!(Arc::ptr_eq(&mol.conformers, &topology_edit.conformers));
        assert!(Arc::ptr_eq(&mol.props, &topology_edit.props));

        let mut conformer_edit = mol.clone();
        conformer_edit.set_coords_2d(None);
        assert!(Arc::ptr_eq(&mol.topology, &conformer_edit.topology));
        assert!(!Arc::ptr_eq(&mol.conformers, &conformer_edit.conformers));
        assert!(Arc::ptr_eq(&mol.props, &conformer_edit.props));

        let mut props_edit = mol.clone();
        props_edit.props_mut().name = Some("ligand".to_owned());
        assert!(Arc::ptr_eq(&mol.topology, &props_edit.topology));
        assert!(Arc::ptr_eq(&mol.conformers, &props_edit.conformers));
        assert!(!Arc::ptr_eq(&mol.props, &props_edit.props));
    }

    #[test]
    fn value_coordinate_transform_detaches_only_conformers() {
        let mut mol = Molecule::new();
        mol.add_atom(carbon(0));

        let with_coords = mol
            .with_2d_coords()
            .expect("single atom 2D coordinates should compute");

        assert!(Arc::ptr_eq(&mol.topology, &with_coords.topology));
        assert!(!Arc::ptr_eq(&mol.conformers, &with_coords.conformers));
        assert!(Arc::ptr_eq(&mol.props, &with_coords.props));
        assert!(mol.coords_2d().is_none());
        assert!(with_coords.coords_2d().is_some());
    }

    #[test]
    fn value_metadata_transform_detaches_only_props() {
        let mol = ethane_with_2d_coords();
        let named = mol.with_name("ligand");

        assert!(Arc::ptr_eq(&mol.topology, &named.topology));
        assert!(Arc::ptr_eq(&mol.conformers, &named.conformers));
        assert!(!Arc::ptr_eq(&mol.props, &named.props));
        assert_eq!(mol.props().name, None);
        assert_eq!(named.props().name.as_deref(), Some("ligand"));
    }
}
