//! IO pipelines.
//!
//! IO must be staged as parse -> builder -> finalization policy -> molecule.
//! Agents must not reintroduce parser code that mutates `Molecule` directly or
//! silently sanitizes/finalizes without an explicit policy.
//!
//! RDKit marker convention is defined in `dev/source_reproduction_protocol.md`.
//!
//! PDB/mmCIF policy:
//!
//! - `io::bio` is the Gemmi-aligned structural reader path into `BioStructure`.
//! - RDKit-compatible molecule input is layered over explicit `BioStructure`
//!   conversion, not implemented as a second text parser.
//! - RDKit-derived PDB behavior belongs only to Molecule compatibility APIs.
//! - Do not expose parallel public parser modules for the same input format.

use crate::{MOLBLOCK_IO_FEATURE, UnsupportedFeatureError};

pub mod bio;
mod gemmi_spacegroup_table;
pub mod mol2;
pub mod molblock;
pub mod molfile;
pub mod pdb_molecule;
pub mod pdb_writer;
pub mod sdf;
pub mod xyz;

pub fn dependency_versions() -> Result<(&'static str, &'static str), UnsupportedFeatureError> {
    Err(UnsupportedFeatureError::from_spec(&MOLBLOCK_IO_FEATURE))
}
