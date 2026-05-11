//! Placeholder for future IO pipelines.
//!
//! IO must eventually be staged as parse -> builder -> finalization policy ->
//! molecule. Agents must not reintroduce parser code that mutates `Molecule`
//! directly or silently sanitizes/finalizes without an explicit policy.

use std::path::Path;

use crate::{MOLBLOCK_IO_FEATURE, Molecule, UnsupportedFeatureError};

pub mod bio;
pub mod sdf;

pub fn dependency_versions() -> Result<(&'static str, &'static str), UnsupportedFeatureError> {
    Err(UnsupportedFeatureError::from_spec(&MOLBLOCK_IO_FEATURE))
}

pub mod molfile {
    use super::*;
    use crate::io::sdf::SdfReadError;

    #[derive(Debug, Clone, PartialEq)]
    pub struct MolFileRecord {
        pub molecule: Molecule,
        pub name: Option<String>,
    }

    pub fn read_mol_file(path: impl AsRef<Path>) -> Result<MolFileRecord, SdfReadError> {
        let text =
            std::fs::read_to_string(path).map_err(|err| SdfReadError::Parse(err.to_string()))?;
        read_mol_record_from_str(&text)
    }

    pub fn read_mol_record_from_str(s: &str) -> Result<MolFileRecord, SdfReadError> {
        let record = crate::io::sdf::read_sdf_from_str_with_params(
            s,
            crate::io::sdf::SdfReadParams {
                process_property_lists: false,
                ..Default::default()
            },
        )?;
        let name = record.molecule.properties().name().map(str::to_string);
        Ok(MolFileRecord {
            molecule: record.molecule,
            name,
        })
    }
}

pub mod molblock {
    use super::*;

    #[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
    pub enum MolWriteError {
        #[error("MolBlock writing is not implemented")]
        NotImplemented,
        #[error(transparent)]
        UnsupportedFeature(#[from] UnsupportedFeatureError),
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum SdfFormat {
        V2000,
        V3000,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct MolBlockWriteParams {
        pub format: SdfFormat,
        pub force_2d: bool,
    }

    impl Default for MolBlockWriteParams {
        fn default() -> Self {
            Self {
                format: SdfFormat::V2000,
                force_2d: false,
            }
        }
    }

    pub fn mol_to_v2000_2d_block(_molecule: &Molecule) -> Result<String, MolWriteError> {
        Err(UnsupportedFeatureError::from_spec(&MOLBLOCK_IO_FEATURE).into())
    }

    pub fn mol_to_v2000_3d_block(_molecule: &Molecule) -> Result<String, MolWriteError> {
        Err(UnsupportedFeatureError::from_spec(&MOLBLOCK_IO_FEATURE).into())
    }

    pub fn mol_to_v2000_block(_molecule: &Molecule) -> Result<String, MolWriteError> {
        Err(UnsupportedFeatureError::from_spec(&MOLBLOCK_IO_FEATURE).into())
    }

    pub fn mol_to_2d_sdf_record(
        _molecule: &Molecule,
        _format: SdfFormat,
    ) -> Result<String, MolWriteError> {
        Err(UnsupportedFeatureError::from_spec(&MOLBLOCK_IO_FEATURE).into())
    }

    pub fn mol_to_3d_sdf_record(
        _molecule: &Molecule,
        _format: SdfFormat,
    ) -> Result<String, MolWriteError> {
        Err(UnsupportedFeatureError::from_spec(&MOLBLOCK_IO_FEATURE).into())
    }
}
