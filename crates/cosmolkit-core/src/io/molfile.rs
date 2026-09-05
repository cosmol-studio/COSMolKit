//! Molfile convenience reader layered over SDF parsing.

use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use crate::Molecule;
use crate::io::sdf::{
    SdfReadError, SdfReadParams, read_mol_block_molecule_with_params, read_mol_data_stream_molecule_with_params,
};

#[derive(Debug, Clone, PartialEq)]
pub struct MolFileRecord {
    pub molecule: Molecule,
    pub name: Option<String>,
}

pub fn read_mol_file(path: impl AsRef<Path>) -> Result<MolFileRecord, SdfReadError> {
    read_mol_file_with_params(
        path,
        SdfReadParams {
            process_property_lists: false,
            ..Default::default()
        },
    )
}

pub fn read_mol_file_with_params(path: impl AsRef<Path>, params: SdfReadParams) -> Result<MolFileRecord, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: MolFromMolFile
    // RDKit✔️✔️: std::unique_ptr<RWMol> MolFromMolFile(const std::string &fName,
    // RDKit✔️✔️:                                       const MolFileParserParams &params) {
    // RDKit✔️✔️:   std::ifstream inStream(fName.c_str());
    // RDKit✔️✔️:   if (!inStream || (inStream.bad())) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Bad input file " << fName;
    // RDKit✔️✔️:     throw BadFileException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!inStream.eof()) {
    // RDKit✔️✔️:     unsigned int line = 0;
    // RDKit✔️✔️:     return MolFromMolDataStream(inStream, line, params);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return std::unique_ptr<RWMol>();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: MolFromMolFile
    //
    // COSMolKit's public molfile API cannot return a null molecule record.
    // The RDKit null unique_ptr branch is represented at this Result boundary
    // as a parse error, matching the generated RDKit failure golden rows.
    let path = path.as_ref();
    let file = File::open(path).map_err(|_| SdfReadError::Parse(format!("Bad input file {}", path.display())))?;
    let mut reader = BufReader::new(file);
    if reader
        .fill_buf()
        .map_err(|err| SdfReadError::Parse(err.to_string()))?
        .is_empty()
    {
        return Err(SdfReadError::Parse("empty molfile".to_string()));
    }
    let mut line_number = 0;
    let molecule = read_mol_data_stream_molecule_with_params(
        &mut reader,
        &mut line_number,
        SdfReadParams {
            process_property_lists: false,
            ..params
        },
    )?;
    let name = molecule.properties().name().map(str::to_string);
    Ok(MolFileRecord { molecule, name })
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

pub fn read_mol_record_from_str_with_params(s: &str, params: SdfReadParams) -> Result<MolFileRecord, SdfReadError> {
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
    let molecule = read_mol_block_molecule_with_params(
        s,
        SdfReadParams {
            process_property_lists: false,
            ..params
        },
    )?;
    let name = molecule.properties().name().map(str::to_string);
    Ok(MolFileRecord { molecule, name })
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

    const EXPLICIT_H_MOL: &str = r#"explicit-h
  COSMolKit      2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
"#;

    const INVALID_VERSION_MOL: &str = r#"invalid-version
  COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 X2000
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
        assert_eq!(as_2d.molecule.source_coordinate_dim(), Some(CoordinateDimension::TwoD));
        assert!(as_2d.molecule.coordinates_2d().is_some());
        assert!(as_2d.molecule.conformers_3d().is_empty());

        let as_3d = read_mol_record_from_str_with_params(
            FLAT_MOL,
            SdfReadParams {
                coordinate_mode: SdfCoordinateMode::Require3D,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(as_3d.molecule.coordinates_2d().is_none());
        assert_eq!(as_3d.molecule.conformers_3d().len(), 1);
        assert!(as_3d.molecule.conformers_3d()[0].is_3d());
    }

    #[test]
    fn molfile_reader_accepts_sdf_record_separator_after_m_end_like_rdkit() {
        let record = read_mol_record_from_str(&format!("{FLAT_MOL}$$$$\n")).unwrap();
        assert_eq!(record.molecule.num_atoms(), 1);
    }

    #[test]
    fn molfile_reader_ignores_unread_trailing_text_after_m_end() {
        let record = read_mol_record_from_str(&format!("{FLAT_MOL}>  <FIELD>\nvalue\n\n$$$$\n")).unwrap();
        assert_eq!(record.molecule.num_atoms(), 1);
        assert_eq!(record.molecule.prop("FIELD"), None);
    }

    #[test]
    fn molfile_reader_reports_empty_input_and_missing_counts_line() {
        assert!(read_mol_record_from_str("").is_err());
        assert!(read_mol_record_from_str("name\n  COSMolKit      2D\ncomment\n").is_err());
    }

    #[test]
    fn molfile_reader_honors_strict_parsing_for_ctab_version_string() {
        let strict = read_mol_record_from_str_with_params(
            INVALID_VERSION_MOL,
            SdfReadParams {
                strict_parsing: true,
                ..Default::default()
            },
        );
        assert!(strict.is_err());

        let lax = read_mol_record_from_str_with_params(
            INVALID_VERSION_MOL,
            SdfReadParams {
                strict_parsing: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(lax.molecule.num_atoms(), 1);
    }

    #[test]
    fn molfile_reader_only_removes_hs_when_sanitize_is_enabled() {
        let sanitized_remove_hs = read_mol_record_from_str_with_params(
            EXPLICIT_H_MOL,
            SdfReadParams {
                sanitize: true,
                remove_hs: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(sanitized_remove_hs.molecule.num_atoms(), 1);

        let sanitized_keep_hs = read_mol_record_from_str_with_params(
            EXPLICIT_H_MOL,
            SdfReadParams {
                sanitize: true,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(sanitized_keep_hs.molecule.num_atoms(), 2);

        let unsanitized_remove_hs = read_mol_record_from_str_with_params(
            EXPLICIT_H_MOL,
            SdfReadParams {
                sanitize: false,
                remove_hs: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(unsanitized_remove_hs.molecule.num_atoms(), 2);
    }

    #[test]
    fn molfile_path_reader_reports_file_open_failure() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("missing.mol");
        let err = read_mol_file_with_params(&missing, SdfReadParams::default()).unwrap_err();
        assert!(err.to_string().contains("Bad input file"));
    }

    #[test]
    fn molfile_path_reader_reports_file_read_failure() {
        let dir = tempfile::tempdir().unwrap();
        let err = read_mol_file_with_params(dir.path(), SdfReadParams::default()).unwrap_err();
        assert!(!err.to_string().is_empty());
    }

    #[test]
    fn molfile_path_reader_reports_empty_file_as_no_molecule() {
        let file = tempfile::NamedTempFile::new().unwrap();
        let err = read_mol_file_with_params(file.path(), SdfReadParams::default()).unwrap_err();
        assert_eq!(err.to_string(), "empty molfile");
    }

    #[test]
    fn molfile_path_reader_reads_single_record_molfile() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut file, FLAT_MOL.as_bytes()).unwrap();

        let record = read_mol_file_with_params(file.path(), SdfReadParams::default()).unwrap();
        assert_eq!(record.name.as_deref(), Some("flat"));
        assert_eq!(record.molecule.num_atoms(), 1);
        assert_eq!(record.molecule.num_bonds(), 0);
    }
}
