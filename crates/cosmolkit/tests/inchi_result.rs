#![cfg(feature = "cap-inchi")]

use cosmolkit::{InchiWriteParams, Molecule};

#[test]
fn complete_generation_preserves_identifier_and_native_details() {
    let mol = Molecule::from_smiles("CCO").unwrap();
    for options in ["", "/FixedH", "/SNon"] {
        let params = InchiWriteParams::new(options.into());
        let result = mol.to_inchi_result_with_params(&params).unwrap();
        assert_eq!(result.inchi(), mol.to_inchi_with_params(&params).unwrap());
        assert_eq!(result.return_code(), 0);
        assert!(result.aux_info().starts_with("AuxInfo="));
        assert!(result.message().is_empty());
        assert!(result.diagnostics().is_empty());
    }
}
