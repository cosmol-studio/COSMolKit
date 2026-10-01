#![cfg(feature = "cap-bio")]
use cosmolkit::{
    BINDING_CONTRACT, BindingDefault, BindingItem, BindingKind, BindingOwner, BindingTypeRole,
    BioCoordinateFormat as F, BioReadError, BioReadParams, BioStructure, FunctionStatus, Protein,
    ProteinReadError,
};
use std::error::Error;

#[test]
fn bio_read_biostructure_from_text_with_params_routes_and_signature() {
    let _: fn(&str, &BioReadParams) -> Result<BioStructure, BioReadError> =
        BioStructure::from_text_with_params;
    let row = BINDING_CONTRACT
        .iter()
        .find(|r| r.semantic_id == "BioStructure.from_text_with_params")
        .unwrap();
    assert_eq!(row.python_name, "from_text_with_params");
    assert_eq!(row.javascript_name, "fromTextWithParams");
    assert_eq!(row.callable.unwrap().parameters.len(), 2);
    let inputs = [
        (F::Pdb, "HEADER    TEST\n", F::Pdb),
        (F::Mmcif, "data_demo\n_entry.id DEMO\n", F::Mmcif),
        (
            F::Mmjson,
            r#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#,
            F::Mmcif,
        ),
        (F::Unknown, "data_demo\n_entry.id DEMO\n", F::Mmcif),
        (F::Detect, "HEADER    TEST\n", F::Pdb),
    ];
    for (format, text, expected) in inputs {
        let structure = BioStructure::from_text_with_params(
            text,
            &BioReadParams {
                format,
                source_name: "custom.source".into(),
            },
        )
        .unwrap();
        assert_eq!(structure.input_format(), expected);
        // Gemmi derives a structure name from the CIF data block, but PDB
        // retains the supplied source name.
        let expected_name = if expected == F::Pdb {
            "custom.source"
        } else if format == F::ChemComp {
            "comp_CMP"
        } else {
            "demo"
        };
        assert_eq!(structure.source_state().name, expected_name);
        structure.validate().unwrap();
    }
}

#[test]
fn bio_read_biostructure_from_text_detects_content_and_signature() {
    let _: fn(&str) -> Result<BioStructure, BioReadError> = BioStructure::from_text;
    let row = BINDING_CONTRACT
        .iter()
        .find(|r| r.semantic_id == "BioStructure.from_text")
        .unwrap();
    assert_eq!(row.python_name, "from_text");
    assert_eq!(row.javascript_name, "fromText");
    assert_eq!(row.callable.unwrap().parameters.len(), 1);
    for (text, expected) in [
        ("HEADER    TEST\n", F::Pdb),
        ("data_demo\n_entry.id DEMO\n", F::Mmcif),
        (r#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#, F::Mmcif),
        (
            "data_comp_CMP\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n",
            F::ChemComp,
        ),
    ] {
        let structure = BioStructure::from_text(text).unwrap();
        assert_eq!(structure.input_format(), expected);
        structure.validate().unwrap();
    }
    assert!(matches!(
        BioStructure::from_text(""),
        Err(BioReadError::WrongFormat { .. })
    ));
}

#[test]
fn bio_read_biostructure_read_with_format_dispatches_and_preserves_sources() {
    let _: fn(&std::path::Path, F) -> Result<BioStructure, BioReadError> =
        BioStructure::read_with_format;
    let row = BINDING_CONTRACT
        .iter()
        .find(|r| r.semantic_id == "BioStructure.read_with_format")
        .unwrap();
    assert_eq!(row.python_name, "read_with_format");
    assert_eq!(row.javascript_name, "readWithFormat");
    assert_eq!(row.callable.unwrap().parameters.len(), 2);
    let dir = std::env::temp_dir().join(format!(
        "cosmolkit-bio-read-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("sample.cif");
    std::fs::write(&path, "HEADER    TEST\n").unwrap();
    assert!(matches!(
        BioStructure::read_with_format(&path, F::Unknown),
        Err(BioReadError::Cif(_))
    ));
    assert_eq!(
        BioStructure::read_with_format(&path, F::Detect)
            .unwrap()
            .input_format(),
        F::Pdb
    );
    assert_eq!(
        BioStructure::read_with_format(&path, F::Pdb)
            .unwrap()
            .input_format(),
        F::Pdb
    );
    let cif = "data_demo\n_entry.id DEMO\n";
    std::fs::write(&path, cif).unwrap();
    assert_eq!(
        BioStructure::read_with_format(&path, F::Mmcif)
            .unwrap()
            .input_format(),
        F::Mmcif
    );
    assert!(matches!(
        BioStructure::read_with_format(&path, F::ChemComp),
        Err(BioReadError::ChemComp(_))
    ));
    assert_eq!(
        BioStructure::read_with_format(&path, F::Unknown)
            .unwrap()
            .input_format(),
        F::Mmcif
    );
    std::fs::write(&path, r#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#).unwrap();
    assert_eq!(
        BioStructure::read_with_format(&path, F::Mmjson)
            .unwrap()
            .input_format(),
        F::Mmjson
    );
    assert_eq!(
        BioStructure::from_text_with_params(
            r#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#,
            &BioReadParams {
                format: F::Mmjson,
                source_name: "memory".into()
            }
        )
        .unwrap()
        .input_format(),
        F::Mmcif
    );
    let missing = dir.join("missing.cif");
    assert!(
        BioStructure::read_with_format(&missing, F::Unknown)
            .unwrap_err()
            .source()
            .is_some()
    );
    let comp = dir.join("comp.cif");
    std::fs::write(&comp, "data_comp_CMP\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n").unwrap();
    assert_eq!(
        BioStructure::read_with_format(&comp, F::ChemComp)
            .unwrap()
            .input_format(),
        F::ChemComp
    );
    std::fs::remove_file(&comp).unwrap();
    std::fs::remove_file(&path).unwrap();
    std::fs::remove_dir(&dir).unwrap();
}

#[test]
fn bio_read_biostructure_read_extension_and_signature() {
    let _: fn(&std::path::Path) -> Result<BioStructure, BioReadError> = BioStructure::read;
    let row = BINDING_CONTRACT
        .iter()
        .find(|r| r.semantic_id == "BioStructure.read")
        .unwrap();
    assert_eq!(row.python_name, "read");
    assert_eq!(row.javascript_name, "read");
    assert_eq!(row.callable.unwrap().parameters.len(), 1);
    let dir = std::env::temp_dir().join(format!(
        "cosmolkit-bio-extensions-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&dir).unwrap();
    for (extension, text, expected) in [
        ("pdb", "HEADER    TEST\n", F::Pdb),
        ("ent", "HEADER    TEST\n", F::Pdb),
        ("cif", "data_demo\n_entry.id DEMO\n", F::Mmcif),
        ("mmcif", "data_demo\n_entry.id DEMO\n", F::Mmcif),
        (
            "json",
            r#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#,
            F::Mmjson,
        ),
    ] {
        let path = dir.join(format!("sample.{extension}"));
        std::fs::write(&path, text).unwrap();
        assert_eq!(BioStructure::read(&path).unwrap().input_format(), expected);
        std::fs::remove_file(path).unwrap();
    }
    let unknown = dir.join("sample.xyz");
    std::fs::write(&unknown, "HEADER    TEST\n").unwrap();
    assert!(
        matches!(BioStructure::read(&unknown), Err(BioReadError::UnknownFileFormat(ref path)) if path == &unknown)
    );
    std::fs::remove_file(&unknown).unwrap();
    let mismatched = dir.join("mismatch.cif");
    std::fs::write(&mismatched, "HEADER    TEST\n").unwrap();
    assert!(matches!(
        BioStructure::read(&mismatched),
        Err(BioReadError::Cif(_))
    ));
    std::fs::remove_file(&mismatched).unwrap();
    let missing = dir.join("missing.cif");
    assert!(BioStructure::read(&missing).unwrap_err().source().is_some());
    std::fs::remove_dir(dir).unwrap();
}

// Required atom-site columns and row order follow pinned mmcif.hpp's table.
const MIXED_CIF: &str = "data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\nATOM 1 C CA . ALA A 1 1 1.0 2.0 3.0 11 ALA X\nHETATM 2 O O . HOH B 2 . 4.0 5.0 6.0 12 HOH Y\n";
const MIXED_JSON: &str = r#"{"data_demo":{"atom_site":{"group_PDB":["ATOM","HETATM"],"id":[1,2],"type_symbol":["C","O"],"label_atom_id":["CA","O"],"label_alt_id":[null,null],"label_comp_id":["ALA","HOH"],"label_asym_id":["A","B"],"label_entity_id":["1","2"],"label_seq_id":[1,null],"Cartn_x":[1,4],"Cartn_y":[2,5],"Cartn_z":[3,6],"auth_seq_id":[11,12],"auth_comp_id":["ALA","HOH"],"auth_asym_id":["X","Y"]}}}"#;

#[test]
fn bio_read_protein_from_text_with_params_projects_and_retains_source() {
    let _: fn(&str, &BioReadParams) -> Result<Protein, ProteinReadError> =
        Protein::from_text_with_params;
    let row = BINDING_CONTRACT
        .iter()
        .find(|r| r.semantic_id == "Protein.from_text_with_params")
        .unwrap();
    assert_eq!(row.python_name, "from_text_with_params");
    assert_eq!(row.javascript_name, "fromTextWithParams");
    assert_eq!(row.callable.unwrap().parameters.len(), 2);
    let text = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nHETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O  \n";
    for source_name in ["custom.pdb", "custom.source"] {
        for format in [F::Unknown, F::Detect, F::Pdb] {
            let params = BioReadParams {
                format,
                source_name: source_name.into(),
            };
            let structure = BioStructure::from_text_with_params(text, &params).unwrap();
            let protein = Protein::from_text_with_params(text, &params).unwrap();
            assert_eq!(protein.input_format(), F::Pdb);
            assert_eq!(protein.num_atoms(), 1);
            assert_eq!(structure.num_atoms(), 2);
            assert_eq!(
                protein.as_bio_structure().coordinates().positions(),
                &[[1.0, 2.0, 3.0]]
            );
            assert_eq!(
                protein.as_bio_structure().source_state().name,
                structure.source_state().name
            );
            assert_eq!(
                structure.source_state().name,
                if source_name == "custom.pdb" {
                    "custom"
                } else {
                    "custom.source"
                }
            );
            assert_eq!(
                protein.as_bio_structure().atoms()[0]
                    .source()
                    .serial()
                    .unwrap()
                    .value(),
                1
            );
            assert_eq!(
                protein.as_bio_structure().residues()[0].name().as_str(),
                "ALA"
            );
            assert_eq!(
                protein.as_bio_structure().residues()[0]
                    .source()
                    .seq_id()
                    .unwrap()
                    .seq_num(),
                1
            );
        }
    }
    for (text, explicit, expected_format) in [
        (MIXED_CIF, F::Mmcif, F::Mmcif),
        (MIXED_JSON, F::Mmjson, F::Mmcif),
    ] {
        for source_name in ["custom.pdb", "custom.source"] {
            for format in [F::Unknown, F::Detect, explicit] {
                let params = BioReadParams {
                    format,
                    source_name: source_name.into(),
                };
                let structure = BioStructure::from_text_with_params(text, &params).unwrap();
                let protein = Protein::from_text_with_params(text, &params).unwrap();
                assert_eq!(structure.num_atoms(), 2, "{format:?}/{source_name}");
                assert_eq!(protein.num_atoms(), 1, "{format:?}/{source_name}");
                assert_eq!(protein.input_format(), expected_format);
                assert_eq!(protein.as_bio_structure().source_state().name, "demo");
                assert_eq!(
                    protein.as_bio_structure().coordinates().positions(),
                    &[[1.0, 2.0, 3.0]]
                );
                assert_eq!(protein.as_bio_structure().atoms()[0].name().as_str(), "CA");
                assert_eq!(
                    protein.as_bio_structure().atoms()[0]
                        .source()
                        .serial()
                        .unwrap()
                        .value(),
                    1
                );
                assert_eq!(
                    protein.as_bio_structure().residues()[0].name().as_str(),
                    "ALA"
                );
                assert_eq!(
                    protein.as_bio_structure().residues()[0]
                        .source()
                        .label_seq_id(),
                    Some(1)
                );
                assert_eq!(
                    protein.as_bio_structure().residues()[0]
                        .source()
                        .seq_id()
                        .unwrap()
                        .seq_num(),
                    11
                );
                assert_eq!(
                    protein.as_bio_structure().chains()[0]
                        .source()
                        .label_asym_id(),
                    Some("A")
                );
            }
        }
    }
    let parse_error = Protein::from_text_with_params("", &BioReadParams::default()).unwrap_err();
    assert!(matches!(
        parse_error,
        ProteinReadError::Structure(BioReadError::WrongFormat { .. })
    ));
    assert!(parse_error.source().is_some());
    let invalid = Protein::from_text_with_params(
        text,
        &BioReadParams {
            format: F::ChemComp,
            source_name: "special".into(),
        },
    )
    .unwrap_err();
    assert!(matches!(
        invalid,
        ProteinReadError::Structure(BioReadError::WrongFormat { .. })
    ));
}

#[test]
fn bio_read_protein_from_text_detects_and_projects() {
    let _: fn(&str) -> Result<Protein, ProteinReadError> = Protein::from_text;
    let row = BINDING_CONTRACT
        .iter()
        .find(|r| r.semantic_id == "Protein.from_text")
        .unwrap();
    assert_eq!(row.python_name, "from_text");
    assert_eq!(row.javascript_name, "fromText");
    assert_eq!(row.callable.unwrap().parameters.len(), 1);
    let aa = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n";
    let water =
        "HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O  \n";
    assert_eq!(Protein::from_text(aa).unwrap().num_atoms(), 1);
    assert_eq!(
        Protein::from_text(&format!("{aa}{water}"))
            .unwrap()
            .num_atoms(),
        1
    );
    let empty = Protein::from_text(water).unwrap();
    assert_eq!(empty.num_atoms(), 0);
    assert_eq!(empty.input_format(), F::Pdb);
    assert!(matches!(
        Protein::from_text(""),
        Err(ProteinReadError::Structure(
            BioReadError::WrongFormat { .. }
        ))
    ));
}

#[test]
fn bio_read_protein_read_with_format_dispatch_and_errors() {
    let _: fn(&std::path::Path, F) -> Result<Protein, ProteinReadError> = Protein::read_with_format;
    let row = BINDING_CONTRACT
        .iter()
        .find(|r| r.semantic_id == "Protein.read_with_format")
        .unwrap();
    assert_eq!(row.python_name, "read_with_format");
    assert_eq!(row.javascript_name, "readWithFormat");
    assert_eq!(row.callable.unwrap().parameters.len(), 2);
    let dir = std::env::temp_dir().join(format!(
        "cosmolkit-protein-format-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("sample.cif");
    let pdb = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n";
    std::fs::write(&path, pdb).unwrap();
    assert!(matches!(
        Protein::read_with_format(&path, F::Unknown),
        Err(ProteinReadError::Structure(BioReadError::Cif(_)))
    ));
    for format in [F::Detect, F::Pdb] {
        let protein = Protein::read_with_format(&path, format).unwrap();
        assert_eq!(protein.num_atoms(), 1);
        assert_eq!(protein.input_format(), F::Pdb);
        assert_eq!(
            protein.as_bio_structure().coordinates().positions(),
            &[[1.0, 2.0, 3.0]]
        );
    }
    std::fs::write(&path, "data_demo\n_entry.id DEMO\n").unwrap();
    assert_eq!(
        Protein::read_with_format(&path, F::Mmcif)
            .unwrap()
            .input_format(),
        F::Mmcif
    );
    assert!(matches!(
        Protein::read_with_format(&path, F::ChemComp),
        Err(ProteinReadError::Structure(BioReadError::ChemComp(_)))
    ));
    std::fs::write(&path, r#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#).unwrap();
    assert_eq!(
        Protein::read_with_format(&path, F::Mmjson)
            .unwrap()
            .input_format(),
        F::Mmjson
    );
    std::fs::remove_file(&path).unwrap();
    std::fs::remove_dir(dir).unwrap();
}

#[test]
fn bio_read_protein_read_extension_and_projection() {
    let _: fn(&std::path::Path) -> Result<Protein, ProteinReadError> = Protein::read;
    let row = BINDING_CONTRACT
        .iter()
        .find(|r| r.semantic_id == "Protein.read")
        .unwrap();
    assert_eq!(row.python_name, "read");
    assert_eq!(row.javascript_name, "read");
    assert_eq!(row.callable.unwrap().parameters.len(), 1);
    let dir = std::env::temp_dir().join(format!(
        "cosmolkit-protein-read-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("sample.pdb");
    std::fs::write(&path, "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nHETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O  \n").unwrap();
    let protein = Protein::read(&path).unwrap();
    assert_eq!(protein.num_atoms(), 1);
    assert_eq!(protein.input_format(), F::Pdb);
    std::fs::write(
        &path,
        "HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 20.00           O  \n",
    )
    .unwrap();
    assert_eq!(Protein::read(&path).unwrap().num_atoms(), 0);
    std::fs::remove_file(&path).unwrap();
    let missing = dir.join("missing.pdb");
    assert!(matches!(
        Protein::read(&missing),
        Err(ProteinReadError::Structure(BioReadError::Io { .. }))
    ));
    std::fs::remove_dir(dir).unwrap();
}

#[test]
fn bio_read_registry_all_eight_real_constructors_and_supporting_types() {
    for target in ["BioStructure", "Protein"] {
        let error_type = if target == "Protein" {
            "crate :: ProteinReadError"
        } else {
            "crate :: BioReadError"
        };
        for (name, js, argument_types) in [
            ("from_text", "fromText", vec!["& str"]),
            (
                "from_text_with_params",
                "fromTextWithParams",
                vec!["& str", "& crate :: BioReadParams"],
            ),
            ("read", "read", vec!["& std :: path :: Path"]),
            (
                "read_with_format",
                "readWithFormat",
                vec!["& std :: path :: Path", "crate :: BioCoordinateFormat"],
            ),
        ] {
            let id = format!("{target}.{name}");
            let rows: Vec<_> = BINDING_CONTRACT
                .iter()
                .filter(|r| r.semantic_id == id)
                .collect();
            assert_eq!(rows.len(), 1, "{id}");
            let row = rows[0];
            assert_eq!(row.item, BindingItem::Callable);
            assert_eq!(row.owner, BindingOwner::Type);
            assert_eq!(row.status, FunctionStatus::Experimental);
            assert_eq!(row.feature, "cap-bio");
            assert_eq!(row.status, FunctionStatus::Experimental);
            assert_eq!(row.python_name, name);
            assert_eq!(row.javascript_name, js);
            let callable = row.callable.unwrap();
            assert_eq!(callable.kind, BindingKind::Static);
            assert_eq!(callable.error_type, Some(error_type));
            assert_eq!(callable.parameters.len(), argument_types.len());
            for (parameter, expected) in callable.parameters.iter().zip(argument_types) {
                assert_eq!(parameter.type_name, expected);
                assert_eq!(parameter.default, BindingDefault::Required);
            }
        }
    }
    for (name, role) in [
        ("BioReadParams", BindingTypeRole::Parameter),
        ("BioReadError", BindingTypeRole::Error),
    ] {
        let id = format!("types.{name}");
        let rows: Vec<_> = BINDING_CONTRACT
            .iter()
            .filter(|r| r.semantic_id == id)
            .collect();
        assert_eq!(rows.len(), 1, "{id}");
        let row = rows[0];
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.type_role, Some(role));
        assert_eq!(row.feature, "cap-bio");
    }
    let _: fn(&str) -> Result<BioStructure, BioReadError> = BioStructure::from_text;
    let _: fn(&str) -> Result<Protein, ProteinReadError> = Protein::from_text;
}

#[test]
fn bio_read_biostructure_from_text_with_params_errors() {
    for (format, text) in [
        (F::Unknown, ""),
        (F::Detect, ""),
        (F::Mmjson, "{"),
        (F::Mmcif, "not a CIF document"),
        (F::ChemComp, "data_comp_CMP\n_chem_comp_atom.atom_id C1\n"),
    ] {
        let error = BioStructure::from_text_with_params(
            text,
            &BioReadParams {
                format,
                source_name: "bad.source".into(),
            },
        )
        .unwrap_err();
        assert!(
            error.source().is_some() || matches!(error, BioReadError::WrongFormat { .. }),
            "{format:?}: {error:?}"
        );
    }
}
