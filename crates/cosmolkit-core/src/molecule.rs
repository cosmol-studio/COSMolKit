use crate::{AdjacencyList, Atom, Bond};
use glam::DVec2;

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

/// Molecular graph with atom/bond tables, optional 2D coordinates, and adjacency cache.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub coords_2d: Option<Vec<DVec2>>,
    pub adjacency: Option<AdjacencyList>,
}

impl Molecule {
    /// Create an empty molecule.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add one atom and return its assigned index.
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let index = self.atoms.len();
        let mut atom = atom;
        atom.index = index;
        self.atoms.push(atom);
        self.coords_2d = None;
        self.adjacency = None;
        index
    }

    /// Add one bond and return its assigned index.
    pub fn add_bond(&mut self, bond: Bond) -> usize {
        let index = self.bonds.len();
        let mut bond = bond;
        bond.index = index;
        self.bonds.push(bond);
        self.adjacency = None;
        index
    }

    /// Rebuild adjacency representation from current topology.
    pub fn rebuild_adjacency(&mut self) {
        self.adjacency = Some(AdjacencyList::from_topology(self.atoms.len(), &self.bonds));
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
        let coords = crate::io::molblock::compute_2d_coords_minimal(self)?;
        self.coords_2d = Some(coords.into_iter().map(|(x, y)| DVec2::new(x, y)).collect());
        Ok(self)
    }

    /// Return stored 2D coordinates, if present.
    #[must_use]
    pub fn coords_2d(&self) -> Option<&[DVec2]> {
        self.coords_2d.as_deref()
    }

    /// Return atom atomic numbers in atom-index order.
    #[must_use]
    pub fn atomic_numbers(&self) -> Vec<u8> {
        self.atoms.iter().map(|atom| atom.atomic_num).collect()
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
