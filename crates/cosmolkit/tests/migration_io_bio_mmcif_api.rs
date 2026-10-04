#![cfg(feature = "cap-bio")]
use cosmolkit::{
    BINDING_CONTRACT, BioCoordinateFormat, BioMmcifReadStage, BioStructure, FunctionStatus,
    Protein, ProteinReadError,
};
use std::error::Error;
const CIF: &str = include_str!("../../../testdata/bio/fixtures/gemmi_full_feature_sample.cif");

#[test]
fn public_mmcif_preserves_hierarchy_metadata_and_short_names() {
    let structure = BioStructure::from_mmcif(CIF).unwrap();
    structure.validate().unwrap();
    assert_eq!(structure.input_format(), BioCoordinateFormat::Mmcif);
    assert_eq!(structure.atoms().len(), 2);
    assert_eq!(structure.atoms()[0].name().as_bytes(), b"SG");
    assert_eq!(structure.source_state().name, "demo");
    assert_eq!(structure.source_state().info["_entry.id"], "9XYZ");
    assert_eq!(structure.metadata().authors, ["DOE, J.", "SMITH, A."]);
    assert_eq!(
        structure.coordinates().positions(),
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
    );
    assert_eq!(structure.entities().len(), 1);
    assert_eq!(structure.assemblies().len(), 1);
    assert_eq!(
        Protein::from_mmcif(CIF).unwrap(),
        structure.protein().unwrap()
    );
}

#[test]
fn public_mmcif_errors_keep_stage_and_underlying_cause() {
    for (text, stage) in [
        ("not a CIF document", BioMmcifReadStage::CifDocument),
        (
            "data_first\ndata_second\n_atom_site.id 1\n",
            BioMmcifReadStage::CoordinateBlock,
        ),
    ] {
        let error = BioStructure::from_mmcif(text).unwrap_err();
        assert_eq!(error.stage(), stage);
        assert!(error.source().is_some());
        let error = Protein::from_mmcif(text).unwrap_err();
        assert!(error.source().is_some());
        assert!(matches!(error, ProteinReadError::Mmcif(ref e) if e.stage() == stage));
    }
}

#[test]
fn public_mmcif_registry_names_match_both_associated_constructors() {
    for target in ["BioStructure", "Protein"] {
        let id = format!("{target}.from_mmcif");
        let row = BINDING_CONTRACT
            .iter()
            .find(|r| r.semantic_id == id)
            .unwrap();
        assert_eq!(row.python_name, "from_mmcif");
        assert_eq!(row.javascript_name, "fromMmcif");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.callable.unwrap().parameters.len(), 1);
    }
    assert!(
        !BINDING_CONTRACT
            .iter()
            .any(|r| r.semantic_id.starts_with("bio.bio_structure_from_")
                || r.semantic_id.starts_with("bio.protein_from_"))
    );
}
