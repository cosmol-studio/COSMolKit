//! Chemistry algorithms over the core molecule model.
//!
//! These modules compute or transform chemical state, but public mutation of an
//! existing `Molecule` still routes through registered operations.

pub mod aromaticity;
pub(crate) mod atom_properties;
pub mod atropisomer;
pub(crate) mod ciplabeler;
pub mod coordinates;
pub mod distgeom;
pub mod forcefield;
pub mod hydrogens;
pub mod kekulize;
pub(crate) mod matrices;
pub mod mol_transforms;
pub mod rings;
pub mod stereo;
pub mod stereo_enumerate;
pub mod valence;
