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
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: MolFromMolBlock
    // RDKit✔️✔️: std::unique_ptr<RWMol> MolFromMolBlock(const std::string &molBlock,
    // RDKit✔️✔️:                                        const MolFileParserParams &params) {
    // RDKit✔️✔️:   std::istringstream inStream(molBlock);
    // RDKit✔️✔️:   unsigned int line = 0;
    // RDKit✔️✔️:   return MolFromMolDataStream(inStream, line, params);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: MolFromMolBlock
    //
    // RDKit's MolFromMolBlock does not parse unread text after the CTAB; it
    // lets MolFromMolDataStream return after M END.
    let mol_block = mol_block_through_m_end(s);
    let record = crate::io::sdf::read_sdf_from_str_with_params(
        mol_block,
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

fn mol_block_through_m_end(input: &str) -> &str {
    let mut end = input.len();
    let mut offset = 0;
    for line in input.split_inclusive('\n') {
        let content = line.trim_end_matches('\n').trim_end_matches('\r');
        offset += line.len();
        if content == "M  END" {
            end = offset;
            break;
        }
    }
    &input[..end]
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
        assert!(as_2d.molecule.coordinates_2d().is_some());

        let as_3d = read_mol_record_from_str_with_params(
            FLAT_MOL,
            SdfReadParams {
                coordinate_mode: SdfCoordinateMode::Require3D,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(as_3d.molecule.coordinates_2d().is_some());
        assert_eq!(as_3d.molecule.conformers_3d().len(), 1);
    }

    #[test]
    fn molfile_reader_accepts_sdf_record_separator_after_m_end_like_rdkit() {
        let record = read_mol_record_from_str(&format!("{FLAT_MOL}$$$$\n")).unwrap();
        assert_eq!(record.molecule.num_atoms(), 1);
    }
}
