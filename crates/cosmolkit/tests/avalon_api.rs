use cosmolkit::{AvalonEngineError, AvalonFingerprintError, AvalonFingerprintParams, Molecule};

#[test]
fn avalon_validates_parameters_before_input_conversion() {
    use std::error::Error;
    let source = Molecule::from_smiles("CCO").unwrap();
    let before = source.clone();
    let error = source
        .fingerprint_avalon_with_params(&AvalonFingerprintParams {
            n_bits: 7,
            ..Default::default()
        })
        .unwrap_err();
    assert!(matches!(
        error,
        AvalonFingerprintError::Engine(AvalonEngineError::InvalidArguments { .. })
    ));
    assert!(error.source().is_some());
    assert_eq!(source, before);
}

#[cfg(feature = "cap-depict")]
#[test]
fn avalon_adapter_preserves_native_bits_and_source_state() {
    let source = Molecule::from_smiles("CCO").unwrap();
    let before = source.clone();
    let params = AvalonFingerprintParams {
        n_bits: 64,
        ..Default::default()
    };
    assert_eq!(
        source
            .fingerprint_avalon_with_params(&params)
            .unwrap()
            .on_bits(),
        [6, 14, 30, 31, 42]
    );
    assert_eq!(source.fingerprint_avalon().unwrap().n_bits(), 512);
    assert_eq!(source, before);
    assert!(source.coordinates_2d().is_none());
}

#[cfg(not(feature = "cap-depict"))]
#[test]
fn avalon_missing_coordinate_capability_is_an_input_error() {
    let source = Molecule::from_smiles("CCO").unwrap();
    assert!(matches!(
        source.fingerprint_avalon(),
        Err(AvalonFingerprintError::Input(_))
    ));
}
