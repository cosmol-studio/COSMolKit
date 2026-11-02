#![cfg(all(
    feature = "cap-serialization",
    feature = "cap-smiles",
    feature = "cap-transforms"
))]

use cosmolkit::{Molecule, MoleculeTransferError};

#[test]
fn transfer_preserves_molecular_value_and_coordinate_bits_without_mutation() {
    let source = Molecule::from_smiles("N[C@@H](C)C(=O)O").unwrap();
    let rows: Vec<[f64; 3]> = (0..source.num_atoms())
        .map(|i| [i as f64 + 0.125, -0.0, i as f64 * -0.25])
        .collect();
    let source = source.with_only_3d_conformer(rows.clone()).unwrap();
    let source = source
        .to_builder()
        .with_name("visualization ligand".into())
        .build()
        .unwrap();
    let before = source.to_binary().unwrap();
    let transfer = source.to_transfer_binary().unwrap();
    let restored = Molecule::from_transfer_binary(&transfer).unwrap();
    assert_eq!(restored.topology(), source.topology());
    assert_eq!(restored.properties(), source.properties());
    assert_eq!(restored.conformers_3d(), source.conformers_3d());
    for (actual, expected) in restored.conformers_3d()[0].coordinates().iter().zip(rows) {
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
    }
    assert_eq!(source.to_binary().unwrap(), before);
    assert_eq!(restored.to_transfer_binary().unwrap(), transfer);
    assert!(Molecule::from_transfer_binary(&before).is_err());
    assert!(Molecule::from_binary(&transfer).is_err());
    assert_eq!(
        Molecule::from_binary(&before).unwrap().to_binary().unwrap(),
        before
    );
}

#[test]
fn transfer_rejects_mismatched_versions_before_payload_and_malformed_envelopes() {
    let source = Molecule::from_smiles("CCO").unwrap();
    let transfer = source.to_transfer_binary().unwrap();
    let mut mismatch = transfer.clone();
    mismatch[10] = b'9';
    // Even with the molecule payload absent, the sender version is rejected first.
    mismatch.truncate(10 + env!("CARGO_PKG_VERSION").len());
    match Molecule::from_transfer_binary(&mismatch).unwrap_err() {
        MoleculeTransferError::VersionMismatch { expected, actual } => {
            assert_eq!(expected, env!("CARGO_PKG_VERSION"));
            assert_ne!(actual, expected);
        }
        other => panic!("expected version mismatch, got {other:?}"),
    }
    for end in [0, 8, 10, transfer.len() - 1] {
        assert!(matches!(
            Molecule::from_transfer_binary(&transfer[..end]),
            Err(MoleculeTransferError::Codec(_))
        ));
    }
    let mut trailing = transfer;
    trailing.push(0);
    assert!(matches!(
        Molecule::from_transfer_binary(&trailing),
        Err(MoleculeTransferError::Codec(_))
    ));
}
