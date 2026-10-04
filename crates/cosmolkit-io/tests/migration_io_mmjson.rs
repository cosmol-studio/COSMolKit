use cosmolkit_bio::{BioCoordinateFormat, BioStructureData};
use cosmolkit_io::{BioReadError, BioReadParams, read_bio_structure, read_bio_structure_file};
use std::error::Error;

#[test]
fn mmjson_detached_entry_preserves_document_order_raw_numbers_and_source() {
    // Pinned json.cpp retains unsorted keys and number strings; mmread.hpp
    // populates an ordinary CIF-derived structure from the mmJSON document.
    let text =
        r#"{"data_second":{"entry":{"id":[1.2300E+04]}},"data_first":{"entry":{"id":["FIRST"]}}}"#;
    let data = read_bio_structure(
        text,
        &BioReadParams {
            format: BioCoordinateFormat::Mmjson,
            source_name: "named.mmjson".into(),
        },
    )
    .unwrap();
    assert_eq!(data.input_format(), BioCoordinateFormat::Mmcif);
    assert_eq!(data.source_state().name, "second");
    assert_eq!(data.source_state().info["_entry.id"], "1.2300E+04");
}

// All required atom-site columns come from mmcif.hpp::populate_structure_from_block.
// Deliberately place C2 before N1; source row order is not alphabetic order.
const POPULATED: &str = r#"{"data_demo":{"entry":{"id":["DEMO"]},"atom_site":{"id":[2,1],"group_PDB":["ATOM","ATOM"],"type_symbol":["C","N"],"label_atom_id":["C2","N1"],"label_alt_id":[null,null],"label_comp_id":["GLY","GLY"],"label_asym_id":["A","A"],"label_entity_id":["1","1"],"label_seq_id":[7,7],"pdbx_PDB_ins_code":[null,null],"Cartn_x":[3.25,-1.5],"Cartn_y":[4,8],"Cartn_z":[5,9],"auth_seq_id":[7,7],"auth_comp_id":["GLY","GLY"],"auth_asym_id":["X","X"],"pdbx_PDB_model_num":[1,1]}}}"#;

fn assert_populated(data: &BioStructureData) {
    data.validate().unwrap();
    assert_eq!(data.source_state().name, "demo");
    assert_eq!(data.source_state().info["_entry.id"], "DEMO");
    assert_eq!(
        (
            data.models().len(),
            data.chains().len(),
            data.residues().len(),
            data.atoms().len()
        ),
        (1, 1, 1, 2)
    );
    assert_eq!(
        (
            data.models()[0].chain_span().start(),
            data.models()[0].chain_span().len()
        ),
        (0, 1)
    );
    assert_eq!(
        (
            data.chains()[0].residue_span().start(),
            data.chains()[0].residue_span().len()
        ),
        (0, 1)
    );
    assert_eq!(
        (
            data.residues()[0].atom_span().start(),
            data.residues()[0].atom_span().len()
        ),
        (0, 2)
    );
    assert_eq!(data.models()[0].source_model_number(), Some(1));
    assert_eq!(data.chains()[0].source().label_asym_id(), Some("A"));
    assert_eq!(data.residues()[0].source().label_seq_id(), Some(7));
    assert_eq!(data.residues()[0].source().subchain_id(), Some("A"));
    assert_eq!(data.residues()[0].name().as_str(), "GLY");
    assert_eq!(
        data.atoms()
            .iter()
            .map(|atom| atom.name().as_str().to_owned())
            .collect::<Vec<_>>(),
        ["C2", "N1"]
    );
    assert_eq!(
        data.atoms()
            .iter()
            .map(|atom| atom.element().symbol())
            .collect::<Vec<_>>(),
        ["C", "N"]
    );
    assert_eq!(
        data.atoms()
            .iter()
            .map(|atom| atom.source().serial().map(|id| id.value()))
            .collect::<Vec<_>>(),
        [Some(2), Some(1)]
    );
    assert_eq!(
        data.coordinates().positions(),
        &[[3.25, 4.0, 5.0], [-1.5, 8.0, 9.0]]
    );
}

#[test]
fn mmjson_detached_populated_atoms_memory_and_file() {
    let memory = read_bio_structure(
        POPULATED,
        &BioReadParams {
            format: BioCoordinateFormat::Mmjson,
            source_name: "atoms.mmjson".into(),
        },
    )
    .unwrap();
    assert_populated(&memory);
    assert_eq!(memory.input_format(), BioCoordinateFormat::Mmcif);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("atoms.json");
    std::fs::write(&path, POPULATED).unwrap();
    let file = read_bio_structure_file(&path, BioCoordinateFormat::Unknown).unwrap();
    assert_populated(&file);
    assert_eq!(file.input_format(), BioCoordinateFormat::Mmjson);
}

#[test]
fn mmjson_detached_incomplete_atom_site_retains_zero_atom_source_result() {
    // The original Step212 counterexample has no required label_alt_id column.
    // Pinned mmcif.hpp's full table lookup returns no atom-site table.
    let incomplete = r#"{"data_demo":{"atom_site":{"id":["1"],"type_symbol":["C"],"label_atom_id":["C1"],"label_comp_id":["GLY"],"label_asym_id":["A"],"label_seq_id":[1],"Cartn_x":[1.2300E+00],"Cartn_y":[2],"Cartn_z":[3]}}}"#;
    let data = read_bio_structure(
        incomplete,
        &BioReadParams {
            format: BioCoordinateFormat::Mmjson,
            source_name: "atoms.mmjson".into(),
        },
    )
    .unwrap();
    assert!(data.atoms().is_empty());
    assert!(data.coordinates().positions().is_empty());
}

#[test]
fn mmjson_detached_file_and_memory_format_are_distinct() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("entry.json");
    let text = r#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#;
    std::fs::write(&path, text).unwrap();
    let memory = read_bio_structure(text, &BioReadParams::default()).unwrap();
    let file = read_bio_structure_file(&path, BioCoordinateFormat::Unknown).unwrap();
    let detect = read_bio_structure_file(&path, BioCoordinateFormat::Detect).unwrap();
    assert_eq!(memory.input_format(), BioCoordinateFormat::Mmcif);
    assert_eq!(detect.input_format(), BioCoordinateFormat::Mmcif);
    assert_eq!(file.input_format(), BioCoordinateFormat::Mmjson);
    assert_eq!(memory.source_state().name, file.source_state().name);
}

#[test]
fn mmjson_detached_errors_keep_parser_source_and_original_path() {
    let error = read_bio_structure(
        "{",
        &BioReadParams {
            format: BioCoordinateFormat::Mmjson,
            source_name: "broken.mmjson".into(),
        },
    )
    .unwrap_err();
    assert!(matches!(error, BioReadError::Mmjson(_)));
    assert!(error.to_string().contains("broken.mmjson"));
    assert!(error.source().is_some());
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("invalid.json");
    std::fs::write(&path, "{ bad").unwrap();
    let error = read_bio_structure_file(&path, BioCoordinateFormat::Unknown).unwrap_err();
    assert!(matches!(error, BioReadError::Mmjson(_)));
    assert!(error.to_string().contains("invalid.json"));
    assert!(error.source().is_some());
}
