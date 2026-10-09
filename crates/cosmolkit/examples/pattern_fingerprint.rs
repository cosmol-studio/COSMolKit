use cosmolkit::{Molecule, PatternFingerprintParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("c1ccccc1O")?;
    let ordinary = molecule.fingerprint_pattern()?;
    let tautomeric = molecule.fingerprint_pattern_with_params(&PatternFingerprintParams {
        n_bits: 2048,
        tautomeric: true,
    })?;

    println!("ordinary Pattern bits: {:?}", ordinary.on_bits());
    println!("tautomeric Pattern bits: {:?}", tautomeric.on_bits());

    let inputs = [molecule, Molecule::from_smiles("CCO")?];
    for (index, input) in inputs.iter().enumerate() {
        println!(
            "input {index}: {:?}",
            input.fingerprint_pattern()?.on_bits()
        );
    }

    Ok(())
}
