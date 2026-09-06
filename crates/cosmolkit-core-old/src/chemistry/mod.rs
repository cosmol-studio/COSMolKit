//! Chemistry algorithms over the core molecule model.
//!
//! These modules compute or transform chemical state, but public mutation of an
//! existing `Molecule` still routes through registered operations.

pub mod aromaticity;
pub(crate) mod atom_properties;
pub mod atropisomer;
pub mod cip;
pub(crate) mod ciplabeler;
pub(crate) mod conformer_selection;
pub mod coordinates;
pub mod distgeom;
pub mod forcefield;
pub mod hydrogens;
pub mod kekulize;
pub(crate) mod matrices;
pub mod mol_align;
pub(crate) mod mol_align_support;
pub mod mol_transforms;
pub(crate) mod numerics;
pub mod potential_stereo;
pub mod rings;
pub mod sanitize;
pub mod stereo;
#[cfg(feature = "stereoisomers")]
pub mod stereo_enumerate;
pub(crate) mod subgraph;
pub mod tautomer;
mod tautomer_transforms;
pub mod valence;
