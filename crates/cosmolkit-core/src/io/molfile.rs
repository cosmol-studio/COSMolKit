use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::Molecule;

use super::sdf::{
    SdfCoordinateMode, SdfReadError, infer_header_coordinate_dimension, parse_mol_data_stream,
    resolve_coordinate_dimension,
};

#[derive(Debug, Clone, PartialEq)]
pub struct MolFileRecord {
    pub title: String,
    pub program_line: Option<String>,
    pub comment_line: Option<String>,
    pub molecule: Molecule,
    pub raw_molblock: String,
}

pub fn read_mol_file(path: impl AsRef<Path>) -> Result<MolFileRecord, SdfReadError> {
    read_mol_file_with_coordinate_mode(path, SdfCoordinateMode::Auto)
}

pub fn read_mol_file_with_coordinate_mode(
    path: impl AsRef<Path>,
    coordinate_mode: SdfCoordinateMode,
) -> Result<MolFileRecord, SdfReadError> {
    let mut text = String::new();
    File::open(path)?.read_to_string(&mut text)?;
    read_mol_record_from_str_with_coordinate_mode(&text, coordinate_mode)
}

pub fn read_mol_record_from_str(s: &str) -> Result<MolFileRecord, SdfReadError> {
    read_mol_record_from_str_with_coordinate_mode(s, SdfCoordinateMode::Auto)
}

pub fn read_mol_record_from_str_with_coordinate_mode(
    s: &str,
    coordinate_mode: SdfCoordinateMode,
) -> Result<MolFileRecord, SdfReadError> {
    let lines = s
        .lines()
        .map(|line| line.strip_suffix('\r').unwrap_or(line).to_owned())
        .collect::<Vec<_>>();
    if lines.iter().all(|line| line.trim().is_empty()) {
        return Err(SdfReadError::Parse(
            "molfile text did not contain a mol block".to_owned(),
        ));
    }

    let title = lines.first().cloned().unwrap_or_default();
    let program_line = lines.get(1).cloned();
    let comment_line = lines.get(2).cloned();
    let source_dim = resolve_coordinate_dimension(
        coordinate_mode,
        infer_header_coordinate_dimension(lines.get(1).map(String::as_str)),
    );

    // RDKit's MolFromMolBlock() is MolFromMolDataStream() over the string and
    // stops at the mol block terminator. SDF property parsing is intentionally
    // not part of this path.
    let (mut molecule, data_start) = parse_mol_data_stream(&lines, source_dim)?;
    molecule.rebuild_adjacency();

    let extra = lines[data_start..]
        .iter()
        .find(|line| !line.trim().is_empty());
    if let Some(line) = extra {
        return Err(SdfReadError::Parse(format!(
            "Extra non-molfile content after M  END is not implemented: '{line}'"
        )));
    }

    let raw_molblock = lines[..data_start].join("\n");
    Ok(MolFileRecord {
        title,
        program_line,
        comment_line,
        molecule,
        raw_molblock,
    })
}
