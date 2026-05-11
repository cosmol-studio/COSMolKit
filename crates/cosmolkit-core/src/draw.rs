use crate::{BondDirection, BondOrder, Molecule};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SvgDrawError {
    #[error("coordinate generation is not implemented: {0}")]
    CoordinateGeneration(String),
    #[error("drawing path is not implemented: {0}")]
    Unsupported(String),
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
    #[error("SVG parse failed: {0}")]
    SvgParse(String),
    #[error("PNG pixmap allocation failed for {width}x{height}")]
    PixmapAllocation { width: u32, height: u32 },
    #[error("PNG encoding failed: {0}")]
    PngEncode(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct PreparedDrawAtom {
    pub index: usize,
    pub atomic_number: u8,
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreparedDrawBond {
    pub index: usize,
    pub begin_atom: usize,
    pub end_atom: usize,
    pub bond_order: BondOrder,
    pub is_aromatic: bool,
    pub direction: BondDirection,
    pub rdkit_direction_name: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PreparedDrawMolecule {
    pub atoms: Vec<PreparedDrawAtom>,
    pub bonds: Vec<PreparedDrawBond>,
}

pub fn mol_to_svg(_molecule: &Molecule, _width: u32, _height: u32) -> Result<String, SvgDrawError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::DRAWING_FEATURE).into())
}

pub fn mol_to_png(
    _molecule: &Molecule,
    _width: u32,
    _height: u32,
) -> Result<Vec<u8>, SvgDrawError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::DRAWING_FEATURE).into())
}

pub fn prepare_mol_for_drawing_parity(
    _molecule: &Molecule,
) -> Result<PreparedDrawMolecule, SvgDrawError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::DRAWING_FEATURE).into())
}
