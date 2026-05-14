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
//! - Future RDKit-compatible molecule input must be layered over `BioStructure`
//!   conversion, not implemented as a second text parser.
//! - RDKit-derived PDB behavior belongs only to Molecule compatibility APIs.
//! - Do not expose parallel public parser modules for the same input format.

use crate::{MOLBLOCK_IO_FEATURE, UnsupportedFeatureError};

pub mod bio;
pub mod molblock;
pub mod molfile;
pub mod pdb_writer;
pub mod sdf;

pub fn dependency_versions() -> Result<(&'static str, &'static str), UnsupportedFeatureError> {
    Err(UnsupportedFeatureError::from_spec(&MOLBLOCK_IO_FEATURE))
}
