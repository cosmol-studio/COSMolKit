#![cfg(feature = "cap-bio")]

use cosmolkit::binding_contract::{BINDING_CONTRACT, FunctionStatus, StateModel};
use cosmolkit::{BioSelection, BioStructure, Protein};

#[test]
fn bio_selection_copy_public_methods_and_registry_contract() {
    fn assert_error<T: std::error::Error>() {}
    assert_error::<cosmolkit::BioSelectionCopyError>();
    assert_error::<cosmolkit::BioSelectionCopyCause>();
    assert_error::<cosmolkit::BioRowTraverseError>();
    assert_error::<cosmolkit::BioRowModelError>();
    assert_error::<cosmolkit::BioRowChainError>();
    let pdb = concat!(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n",
        "ATOM      2  N   GLY B   2       4.000   5.000   6.000  1.00 20.00           N  \n"
    );
    let source = BioStructure::from_pdb(pdb).unwrap();
    let protein = Protein::from_pdb(pdb).unwrap();
    for (cid, atoms, chains, residues) in [
        ("/", 2, 2, 2),
        ("//A", 1, 1, 1),
        ("//*//ZZ", 0, 2, 2),
        ("//Z", 0, 0, 0),
        ("/99", 0, 0, 0),
    ] {
        let selection = BioSelection::from_cid(cid).unwrap();
        let value = source.with_selection(&selection).unwrap();
        let mut inplace = source.clone();
        inplace.retain_selection_(&selection).unwrap();
        assert_eq!(
            (value.num_atoms(), value.num_chains(), value.num_residues()),
            (atoms, chains, residues)
        );
        assert_eq!(value, inplace);
        let value = protein.with_selection(&selection).unwrap();
        let mut inplace = protein.clone();
        inplace.retain_selection_(&selection).unwrap();
        assert_eq!(value, inplace);
        assert_eq!(value.as_bio_structure().num_atoms(), atoms);
        assert_eq!(source.num_atoms(), 2);
        assert_eq!(protein.num_atoms(), 2);
    }
    for target in ["BioStructure", "Protein"] {
        for (method, js, state) in [
            (
                "with_selection",
                "withSelection",
                StateModel::ValueReturning,
            ),
            ("retain_selection_", "retainSelection", StateModel::InPlace),
        ] {
            let id = format!("{target}.{method}");
            let entry = BINDING_CONTRACT
                .iter()
                .find(|entry| entry.semantic_id == id)
                .unwrap();
            assert_eq!(entry.status, FunctionStatus::Experimental);
            assert_eq!(entry.feature, "cap-bio");
            assert_eq!(entry.python_name, method);
            assert_eq!(entry.javascript_name, js);
            let contract = entry.callable.unwrap();
            assert_eq!(contract.state_model, state);
            assert_eq!(contract.operation_semantic_id, Some(method));
            assert!(contract.error_type.unwrap().ends_with("BioOperationError"));
        }
    }
    for name in [
        "BioSelectionCopyError",
        "BioSelectionCopyCause",
        "BioRowTraverseError",
        "BioRowModelError",
        "BioRowChainError",
    ] {
        assert!(
            BINDING_CONTRACT
                .iter()
                .any(|entry| entry.semantic_id == format!("types.{name}")
                    && entry.status == FunctionStatus::Experimental)
        );
    }
    assert!(
        !BINDING_CONTRACT
            .iter()
            .any(|entry| entry.semantic_id.contains("copy_selection_blocks"))
    );
}
