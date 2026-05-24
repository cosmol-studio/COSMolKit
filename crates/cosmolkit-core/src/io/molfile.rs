//! Molfile convenience reader layered over SDF parsing.

use std::path::Path;

use crate::Molecule;
use crate::io::sdf::{SdfReadError, SdfReadParams};

#[derive(Debug, Clone, PartialEq)]
pub struct MolFileRecord {
    pub molecule: Molecule,
    pub name: Option<String>,
}

pub fn read_mol_file(path: impl AsRef<Path>) -> Result<MolFileRecord, SdfReadError> {
    let text = std::fs::read_to_string(path).map_err(|err| SdfReadError::Parse(err.to_string()))?;
    read_mol_record_from_str(&text)
}

pub fn read_mol_file_with_params(
    path: impl AsRef<Path>,
    params: SdfReadParams,
) -> Result<MolFileRecord, SdfReadError> {
    let text = std::fs::read_to_string(path).map_err(|err| SdfReadError::Parse(err.to_string()))?;
    read_mol_record_from_str_with_params(&text, params)
}

pub fn read_mol_record_from_str(s: &str) -> Result<MolFileRecord, SdfReadError> {
    read_mol_record_from_str_with_params(
        s,
        SdfReadParams {
            process_property_lists: false,
            ..Default::default()
        },
    )
}

pub fn read_mol_record_from_str_with_params(
    s: &str,
    params: SdfReadParams,
) -> Result<MolFileRecord, SdfReadError> {
    reject_extra_molfile_content(s)?;
    let record = crate::io::sdf::read_sdf_from_str_with_params(
        s,
        SdfReadParams {
            process_property_lists: false,
            ..params
        },
    )?;
    let name = record.molecule.properties().name().map(str::to_string);
    Ok(MolFileRecord {
        molecule: record.molecule,
        name,
    })
}

fn reject_extra_molfile_content(s: &str) -> Result<(), SdfReadError> {
    let mut offset = 0usize;
    for line in s.split_inclusive('\n') {
        let line_without_newline = line.trim_end_matches(['\r', '\n']);
        let end_offset = offset + line.len();
        if line_without_newline == "M  END" {
            if s[end_offset..].trim().is_empty() {
                return Ok(());
            }
            return Err(SdfReadError::Parse(
                "Extra non-molfile content after M  END".to_string(),
            ));
        }
        offset = end_offset;
    }
    if s.lines()
        .last()
        .is_some_and(|line| line.trim_end_matches('\r') == "M  END")
    {
        return Ok(());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CoordinateDimension, io::sdf::SdfCoordinateMode};

    const FLAT_MOL: &str = r#"flat
  COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"#;

    #[test]
    fn molfile_reader_with_params_can_force_coordinate_dimension() {
        let as_2d = read_mol_record_from_str_with_params(
            FLAT_MOL,
            SdfReadParams {
                coordinate_mode: SdfCoordinateMode::Require2D,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            as_2d.molecule.source_coordinate_dim(),
            Some(CoordinateDimension::TwoD)
        );
        assert!(as_2d.molecule.coords_2d().is_some());

        let as_3d = read_mol_record_from_str_with_params(
            FLAT_MOL,
            SdfReadParams {
                coordinate_mode: SdfCoordinateMode::Require3D,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(as_3d.molecule.coords_2d().is_some());
        assert_eq!(as_3d.molecule.conformers_3d().len(), 1);
    }

    #[test]
    fn molfile_reader_rejects_sdf_record_separator_after_m_end() {
        let err = read_mol_record_from_str(&format!("{FLAT_MOL}$$$$\n")).unwrap_err();
        assert!(err.to_string().contains("Extra non-molfile content"));
    }
}
