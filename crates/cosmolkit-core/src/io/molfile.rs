//! Molfile convenience reader layered over SDF parsing.

use std::path::Path;

use crate::Molecule;
use crate::io::sdf::SdfReadError;

#[derive(Debug, Clone, PartialEq)]
pub struct MolFileRecord {
    pub molecule: Molecule,
    pub name: Option<String>,
}

pub fn read_mol_file(path: impl AsRef<Path>) -> Result<MolFileRecord, SdfReadError> {
    let text = std::fs::read_to_string(path).map_err(|err| SdfReadError::Parse(err.to_string()))?;
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
