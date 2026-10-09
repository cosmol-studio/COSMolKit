use cosmolkit::{FingerprintAdditionalOutput, Molecule, TopologicalTorsionFingerprintParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CCCCO")?;
    let params = TopologicalTorsionFingerprintParams::default();

    let sparse_count =
        molecule.fingerprint_topological_torsion_sparse_count_with_params(&params, None)?;
    let sparse_bit = molecule.fingerprint_topological_torsion_sparse_with_params(&params, None)?;
    let count = molecule.fingerprint_topological_torsion_count_with_params(&params, None)?;
    let bit = molecule.fingerprint_topological_torsion_with_params(&params, None)?;
    println!("sparse count: {:?}", sparse_count.nonzero_elements());
    println!("sparse bit: {:?}", sparse_bit.on_bits());
    println!("folded count: {:?}", count.nonzero_elements());
    println!("explicit bit: {:?}", bit.on_bits());

    let mut output = FingerprintAdditionalOutput::default();
    output.allocate_atom_to_bits();
    output.allocate_atom_counts();
    output.allocate_bit_paths();
    output.allocate_atoms_per_bit();
    let with_output_count =
        molecule.fingerprint_topological_torsion_count_with_params(&params, Some(&mut output))?;
    assert_eq!(with_output_count, count);
    println!(
        "count with output: {:?}",
        with_output_count.nonzero_elements()
    );
    println!("provenance: {output:?}");

    let inputs = [molecule.clone(), Molecule::from_smiles("CCCCC")?];
    let fingerprints = inputs
        .iter()
        .map(|input| input.fingerprint_topological_torsion_with_params(&params, None))
        .collect::<Result<Vec<_>, _>>()?;
    assert_eq!(fingerprints[0], bit);
    Ok(())
}
