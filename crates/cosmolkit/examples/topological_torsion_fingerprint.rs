use cosmolkit::{
    BatchRecord, Molecule, MoleculeBatch, TopologicalTorsionFingerprintOutputRequest,
    TopologicalTorsionFingerprintParams, TopologicalTorsionFingerprintValue,
    TopologicalTorsionFingerprintVector, TopologicalTorsionLegacyKind,
    TopologicalTorsionLegacyParams, topological_torsion_count_fingerprint,
    topological_torsion_fingerprint, topological_torsion_fingerprint_with_output,
    topological_torsion_legacy_fingerprint, topological_torsion_sparse_count_fingerprint,
    topological_torsion_sparse_fingerprint,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CCCCO")?;
    let params = TopologicalTorsionFingerprintParams::default();

    let sparse_count = topological_torsion_sparse_count_fingerprint(&molecule, &params)?;
    let sparse_bit = topological_torsion_sparse_fingerprint(&molecule, &params)?;
    let count = topological_torsion_count_fingerprint(&molecule, &params)?;
    let bit = topological_torsion_fingerprint(&molecule, &params)?;
    println!("sparse count: {:?}", sparse_count.nonzero_elements());
    println!("sparse bit: {:?}", sparse_bit.on_bits());
    println!("folded count: {:?}", count.nonzero_elements());
    println!("explicit bit: {:?}", bit.on_bits());

    let with_output = topological_torsion_fingerprint_with_output(
        &molecule,
        &params,
        TopologicalTorsionFingerprintOutputRequest {
            vector: TopologicalTorsionFingerprintVector::Count,
            atom_to_bits: true,
            atom_counts: true,
            bit_paths: true,
            atoms_per_bit: true,
        },
    )?;
    let TopologicalTorsionFingerprintValue::Count(with_output_count) = with_output.fingerprint
    else {
        unreachable!("the requested vector variant is preserved")
    };
    println!(
        "count with output: {:?}",
        with_output_count.nonzero_elements()
    );
    println!("provenance: {:?}", with_output.additional_output);

    let legacy = topological_torsion_legacy_fingerprint(
        &molecule,
        &TopologicalTorsionLegacyParams {
            kind: TopologicalTorsionLegacyKind::HashedBit,
            ..Default::default()
        },
    )?;
    println!("legacy hashed bit: {legacy:?}");

    let batch = MoleculeBatch::new(vec![
        BatchRecord::Molecule(molecule.clone()),
        BatchRecord::Molecule(Molecule::from_smiles("CCCCC")?),
    ]);
    let batch_bits =
        batch.topological_torsion_fingerprint_list_with_options(&params, Some(2), Some(false))?;
    assert_eq!(batch_bits[0], Some(bit));

    Ok(())
}
