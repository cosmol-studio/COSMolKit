use std::path::PathBuf;

use cosmolkit_core::{
    BioConnectionType, BioCoorFormat, EntityKind, read_bio_structure_from_str,
    read_bio_structure_from_str_with_format,
};

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/bio")
        .join(name)
}

fn fixture_text(name: &str) -> String {
    std::fs::read_to_string(fixture_path(name))
        .unwrap_or_else(|error| panic!("failed to read bio fixture {name}: {error}"))
}

#[test]
fn pdb_fixture_exercises_gemmi_aligned_structural_records() {
    let structure = read_bio_structure_from_str(
        &fixture_text("gemmi_full_feature_sample.pdb"),
        "gemmi_full_feature_sample.pdb",
    )
    .unwrap();

    assert_eq!(structure.input_format(), BioCoorFormat::Pdb);
    assert_eq!(structure.metadata().entry_id.as_deref(), Some("1ABC"));
    assert_eq!(structure.metadata().authors, vec!["DOE, A.-B."]);
    assert_eq!(
        structure.entities()[0].sequence,
        vec!["ALA", "GLY", "SER", "THR", "TYR", "VAL", "LEU", "ASP"]
    );
    assert_eq!(structure.ter_status(), 'e');
    assert_eq!(structure.connections().len(), 2);
    assert_eq!(structure.connections()[0].type_, BioConnectionType::Disulf);
    assert_eq!(structure.connections()[1].type_, BioConnectionType::MetalC);
    assert_eq!(structure.cispeps().len(), 1);
    assert_eq!(structure.cispeps()[0].reported_angle, Some(12.34));
    assert_eq!(structure.mod_residues().len(), 1);
    assert_eq!(structure.mod_residues()[0].parent_comp_id, "MET");
    assert_eq!(structure.entities()[0].dbrefs.len(), 2);
    assert_eq!(structure.entities()[0].dbrefs[0].accession_code, "Q12345");
    assert_eq!(
        structure.entities()[0].dbrefs[1].accession_code,
        "Q12345LONGACCESSION"
    );
    let crystal = structure.crystal().unwrap();
    assert_eq!(crystal.spacegroup_hm.as_deref(), Some("P 21 21 21"));
    assert!(structure.has_origx());
    assert_eq!(structure.origx().vec, [4.0, 5.0, 6.0]);
    assert_eq!(structure.ncs_oper_identity_id(), Some("1"));
    assert_eq!(structure.ncs_operators().len(), 1);
    assert_eq!(structure.ncs_operators()[0].id, "2");
}

#[test]
fn mmcif_fixture_exercises_gemmi_aligned_structural_categories() {
    let structure = read_bio_structure_from_str(
        &fixture_text("gemmi_full_feature_sample.cif"),
        "gemmi_full_feature_sample.cif",
    )
    .unwrap();

    assert_eq!(structure.input_format(), BioCoorFormat::Mmcif);
    assert_eq!(structure.name(), "demo");
    assert_eq!(structure.metadata().entry_id.as_deref(), Some("9XYZ"));
    assert_eq!(structure.metadata().authors, vec!["DOE, J.", "SMITH, A."]);
    assert_eq!(structure.entities().len(), 1);
    assert_eq!(structure.entities()[0].kind, EntityKind::Polymer);
    assert_eq!(structure.entities()[0].sequence, vec!["CYS", "CYS"]);
    assert_eq!(structure.entities()[0].dbrefs.len(), 1);
    assert_eq!(structure.entities()[0].dbrefs[0].accession_code, "P12345");
    assert_eq!(structure.connections().len(), 1);
    assert_eq!(structure.connections()[0].type_, BioConnectionType::Disulf);
    assert_eq!(structure.cispeps().len(), 1);
    assert_eq!(structure.cispeps()[0].reported_angle, Some(-12.5));
    assert_eq!(structure.mod_residues().len(), 1);
    assert_eq!(structure.mod_residues()[0].mod_id, "MOD1");
    assert_eq!(structure.assemblies().len(), 1);
    assert_eq!(structure.assemblies()[0].generators[0].operators.len(), 2);
    assert_eq!(structure.assemblies()[0].absa, Some(12.5));
    assert_eq!(structure.assemblies()[0].ssa, Some(8.0));
    assert_eq!(structure.ncs_oper_identity_id(), Some("I"));
    assert_eq!(structure.ncs_operators().len(), 1);
    assert_eq!(structure.ncs_operators()[0].id, "N");
    assert_eq!(structure.ncs_operators()[0].transform.vec, [1.0, 2.0, 3.0]);
    assert_eq!(structure.entities()[0].sifts_unp_acc, vec!["P12345"]);
    assert_eq!(structure.residues()[0].sifts_unp.unwrap().num, 15);
    let crystal = structure.crystal().unwrap();
    assert_eq!(crystal.spacegroup_hm.as_deref(), Some("P 1"));
    assert!(structure.has_origx());
    assert_eq!(structure.origx().vec, [4.0, 5.0, 6.0]);
}

#[test]
fn chemcomp_fixture_enters_public_structural_dispatch_surface() {
    let structure = read_bio_structure_from_str_with_format(
        &fixture_text("gemmi_chemcomp_sample.cif"),
        "gemmi_chemcomp_sample.cif",
        BioCoorFormat::Mmcif,
    )
    .unwrap();

    assert_eq!(structure.input_format(), BioCoorFormat::ChemComp);
    assert_eq!(structure.name(), "GLY");
    assert_eq!(structure.num_models(), 1);
    assert_eq!(structure.num_atoms(), 1);
}

#[test]
fn structural_dispatch_surface_reads_mmjson_through_public_surface() {
    let mmjson = r#"{
  "data_demo": {
    "atom_site": {
      "group_PDB": ["ATOM"],
      "id": ["1"],
      "type_symbol": ["C"],
      "label_atom_id": ["CA"],
      "label_alt_id": ["?"],
      "label_comp_id": ["ALA"],
      "label_asym_id": ["A"],
      "label_entity_id": ["1"],
      "label_seq_id": ["1"],
      "Cartn_x": [1.0],
      "Cartn_y": [2.0],
      "Cartn_z": [3.0]
    }
  }
}"#;
    let structure = read_bio_structure_from_str(mmjson, "sample.json").unwrap();
    assert_eq!(structure.input_format(), BioCoorFormat::Mmjson);
    assert_eq!(structure.name(), "demo");
    assert_eq!(structure.num_atoms(), 1);
    assert_eq!(
        structure.atom_position(cosmolkit_core::bio::AtomId::new(0)),
        Some([1.0, 2.0, 3.0])
    );
}
