use std::path::PathBuf;

use cosmolkit_core::io::bio::{BioWriteError, MmcifOutputGroups, MmcifWriteOptions};
use cosmolkit_core::{BioAtomId, BioStructure};

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../testdata/bio/fixtures")
        .join(name)
}

fn fixture_text(name: &str) -> String {
    std::fs::read_to_string(fixture_path(name))
        .unwrap_or_else(|error| panic!("failed to read bio fixture {name}: {error}"))
}

#[test]
fn pdb_to_mmcif_preserves_hierarchy_coordinates_identifiers_and_metadata() {
    let structure = BioStructure::from_pdb_str(&fixture_text("gemmi_full_feature_sample.pdb"))
        .expect("PDB fixture should parse");
    let source = structure.clone();

    let output = structure
        .to_mmcif()
        .expect("PDB structure should write as mmCIF");
    let reparsed =
        BioStructure::from_mmcif_str(&output, "roundtrip.cif").expect("written mmCIF should parse");

    assert_eq!(
        structure, source,
        "writing must not mutate the source structure"
    );
    assert_eq!(reparsed.num_models(), structure.num_models());
    assert_eq!(reparsed.num_chains(), structure.num_chains());
    assert_eq!(reparsed.num_residues(), structure.num_residues());
    assert_eq!(reparsed.num_atoms(), structure.num_atoms());
    assert_eq!(reparsed.metadata().entry_id, structure.metadata().entry_id);
    assert_eq!(reparsed.metadata().title, structure.metadata().title);
    assert_eq!(
        reparsed.atom_position(BioAtomId::new(0)),
        structure.atom_position(BioAtomId::new(0))
    );
    assert_eq!(reparsed.connections().len(), structure.connections().len());
    assert_eq!(reparsed.assemblies().len(), structure.assemblies().len());
    assert!(output.contains("_atom_site.Cartn_x"));
    assert!(output.contains("_pdbx_struct_mod_residue.auth_comp_id"));
}

#[test]
fn mmcif_roundtrip_preserves_assembly_connection_and_multimodel_state() {
    let structure = BioStructure::from_mmcif_str(
        &fixture_text("gemmi_full_feature_sample.cif"),
        "gemmi_full_feature_sample.cif",
    )
    .expect("mmCIF fixture should parse");
    let output = structure.to_mmcif().expect("mmCIF fixture should write");
    let reparsed = BioStructure::from_mmcif_str(&output, "rewritten.cif")
        .expect("rewritten fixture should parse");
    assert_eq!(reparsed.num_atoms(), structure.num_atoms());
    assert_eq!(reparsed.assemblies(), structure.assemblies());
    assert_eq!(reparsed.connections().len(), structure.connections().len());
    assert_eq!(reparsed.cispeps().len(), structure.cispeps().len());

    let multi = BioStructure::from_mmcif_str(
        r#"
data_multi
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 2 3 1 ALA X 1
ATOM 2 C CA . ALA A 1 4 5 6 1 ALA X 2
"#,
        "multi.cif",
    )
    .expect("multi-model fixture should parse");
    let multi_output = multi.to_mmcif().expect("multi-model fixture should write");
    let multi_reparsed = BioStructure::from_mmcif_str(&multi_output, "multi-out.cif")
        .expect("multi-model output should parse");
    assert_eq!(multi_reparsed.num_models(), 2);
    assert_eq!(multi_reparsed.num_atoms(), 2);
    assert_eq!(
        multi_reparsed.atom_position(BioAtomId::new(1)),
        Some([4.0, 5.0, 6.0])
    );
}

#[test]
fn output_groups_can_omit_atoms_without_weakening_metadata_output() {
    let structure = BioStructure::from_mmcif_str(
        &fixture_text("gemmi_full_feature_sample.cif"),
        "gemmi_full_feature_sample.cif",
    )
    .expect("fixture should parse");
    let mut groups = MmcifOutputGroups::all(true);
    groups.atoms = false;
    let output = structure
        .to_mmcif_with_options(MmcifWriteOptions {
            groups,
            ..MmcifWriteOptions::default()
        })
        .expect("metadata-only mmCIF should write");

    assert!(output.contains("_entry.id"));
    assert!(output.contains("_entity.id"));
    assert!(output.contains("_pdbx_struct_assembly.id"));
    assert!(!output.contains("_atom_site.id"));
}

#[test]
fn file_writer_matches_string_writer_and_reports_io_errors() {
    let structure = BioStructure::from_mmcif_str(
        &fixture_text("gemmi_full_feature_sample.cif"),
        "gemmi_full_feature_sample.cif",
    )
    .expect("fixture should parse");
    let expected = structure.to_mmcif().expect("string writer should succeed");
    let file = tempfile::NamedTempFile::new().expect("temporary file should be created");

    structure
        .write_mmcif(file.path())
        .expect("file writer should succeed");
    assert_eq!(std::fs::read_to_string(file.path()).unwrap(), expected);

    let options = MmcifWriteOptions {
        compact: true,
        ..MmcifWriteOptions::default()
    };
    let expected_with_options = structure
        .to_mmcif_with_options(options)
        .expect("configured string writer should succeed");
    structure
        .write_mmcif_with_options(file.path(), options)
        .expect("configured file writer should succeed");
    assert_eq!(
        std::fs::read_to_string(file.path()).unwrap(),
        expected_with_options
    );

    let error = structure
        .write_mmcif(file.path().parent().unwrap())
        .expect_err("writing to a directory should fail");
    assert!(matches!(error, BioWriteError::Io(_)));
}
