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
    assert_eq!(
        source.fingerprint_avalon().unwrap().on_bits(),
        [42, 198, 222, 262, 334, 479, 490]
    );
    assert_eq!(source, before);
    assert!(source.coordinates_2d().is_none());
}

#[cfg(not(feature = "cap-depict"))]
#[test]
fn avalon_internal_coordinates_do_not_enable_public_mol_generation() {
    let source = Molecule::from_smiles("CCO").unwrap();
    assert!(source.to_mol().is_err());
    assert!(source.fingerprint_avalon().is_ok());
    assert!(source.coordinates_2d().is_none());
}

#[test]
fn avalon_stereo_conversion_preserves_native_bits_without_writeback() {
    // RDKit 2026.03.6 GetAvalonFP(mol, nBits=64, bitFlags=0x007fff).
    for (smiles, bits) in [
        (
            "N[C@@H](C)C(=O)O",
            &[
                0, 2, 3, 6, 9, 10, 12, 13, 14, 18, 20, 21, 22, 26, 28, 30, 31, 33, 34, 43, 45, 48,
                49, 51, 52, 55, 57, 58, 59, 60, 61,
            ][..],
        ),
        (
            "F/C=C/F",
            &[2, 8, 18, 22, 40, 42, 43, 44, 47, 53, 57, 60][..],
        ),
        (
            "F/C=C\\F",
            &[2, 8, 18, 22, 40, 42, 43, 44, 47, 53, 57, 60][..],
        ),
    ] {
        let source = Molecule::from_smiles(smiles).unwrap();
        let before = source.clone();
        let fingerprint = source
            .fingerprint_avalon_with_params(&AvalonFingerprintParams {
                n_bits: 64,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(fingerprint.on_bits(), bits, "{smiles}");
        assert_eq!(source, before, "{smiles}");
        assert!(source.coordinates_2d().is_none());
    }
}
