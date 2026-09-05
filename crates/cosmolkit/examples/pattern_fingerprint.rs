use cosmolkit::{BatchRecord, Molecule, MoleculeBatch, PatternFingerprintParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("c1ccccc1O")?;
    let ordinary = molecule.pattern_fingerprint(&PatternFingerprintParams::default())?;
    let tautomeric = molecule.pattern_fingerprint(&PatternFingerprintParams {
        n_bits: 2048,
        tautomeric: true,
    })?;

    println!("ordinary Pattern bits: {:?}", ordinary.on_bits());
    println!("tautomeric Pattern bits: {:?}", tautomeric.on_bits());

    let batch = MoleculeBatch::new(vec![
        BatchRecord::Molecule(molecule),
        BatchRecord::Molecule(Molecule::from_smiles("CCO")?),
    ]);
    let fingerprints =
        batch.pattern_fingerprint_list_with_options(&PatternFingerprintParams::default(), Some(2), Some(false))?;
    for (index, fingerprint) in fingerprints.into_iter().enumerate() {
        println!("batch {index}: {:?}", fingerprint.map(|fp| fp.on_bits()));
    }

    Ok(())
}
