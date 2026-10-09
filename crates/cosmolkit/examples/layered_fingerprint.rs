use cosmolkit::{Fingerprint, LayeredFingerprintLayers, LayeredFingerprintParams, Molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("c1ccccc1O")?;
    let mask = Fingerprint::from_on_bits(257, (0..257).filter(|bit| bit % 2 == 0))?;
    let params = LayeredFingerprintParams {
        layers: LayeredFingerprintLayers::ACTIVE,
        min_path: 2,
        max_path: 4,
        fp_size: 257,
        atom_counts: Some(vec![10; molecule.num_atoms()]),
        set_only_bits: Some(mask),
        branched_paths: true,
        from_atoms: Some(vec![0]),
    };

    let before = molecule.clone();
    let result = molecule.fingerprint_layered_with_output_with_params(&params)?;
    println!("bits: {:?}", result.fingerprint.on_bits());
    println!("seeded atom counts: {:?}", result.atom_counts);
    assert_eq!(molecule, before);

    let inputs = [molecule.clone(), Molecule::from_smiles("CCCO")?];
    let fingerprints = inputs
        .iter()
        .map(|input| input.fingerprint_layered_with_params(&params))
        .collect::<Result<Vec<_>, _>>()?;
    assert_eq!(fingerprints[0], result.fingerprint);

    Ok(())
}
