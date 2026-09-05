use cosmolkit::{
    BatchRecord, Fingerprint, LayeredFingerprintLayers, LayeredFingerprintParams, Molecule, MoleculeBatch,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("c1ccccc1O")?;
    let mask = Fingerprint::from_on_bits(257, (0..257).filter(|bit| bit % 2 == 0));
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
    let result = molecule.layered_fingerprint_with_output(&params)?;
    println!("bits: {:?}", result.fingerprint.on_bits());
    println!("seeded atom counts: {:?}", result.atom_counts);
    assert_eq!(molecule, before);

    let batch = MoleculeBatch::new(vec![
        BatchRecord::Molecule(molecule.clone()),
        BatchRecord::Molecule(Molecule::from_smiles("CCCO")?),
    ]);
    let batch_result = batch.layered_fingerprint_list_with_options(&params, Some(2), Some(false))?;
    assert_eq!(batch_result[0], Some(result.fingerprint));

    Ok(())
}
