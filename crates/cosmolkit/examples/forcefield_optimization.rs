use cosmolkit::{MmffOptimizationParams, Molecule, UffOptimizationParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CCO")?.with_hydrogens()?.sanitize()?;

    let molecule = molecule.with_added_3d_conformer(vec![
        vec![0.000, 0.000, 0.000],
        vec![1.540, 0.000, 0.000],
        vec![2.100, 1.200, 0.000],
        vec![-0.600, 0.900, 0.000],
        vec![-0.600, -0.900, 0.000],
        vec![0.000, 0.000, 1.000],
        vec![1.900, -0.900, 0.000],
        vec![1.700, 0.000, 1.000],
        vec![2.900, 1.200, 0.000],
    ])?;

    if molecule.uff_has_all_molecule_params()? {
        let result = molecule.with_uff_optimized_with_params(&UffOptimizationParams {
            max_iterations: 200,
            vdw_threshold: 10.0,
            conformer_id: None,
            ignore_interfragment_interactions: true,
        })?;
        println!(
            "UFF needs_more={} energy={:.6}",
            result.status > 0,
            result.energy
        );
        println!(
            "optimized first atom: {:?}",
            result.molecule.conformers_3d()[0].coordinates()[0]
        );
    }

    if molecule.mmff_has_all_molecule_params()? {
        let result = molecule.with_mmff_optimized_with_params(&MmffOptimizationParams {
            mmff_variant: "MMFF94".into(),
            max_iterations: 200,
            non_bonded_threshold: 100.0,
            conformer_id: None,
            ignore_interfragment_interactions: true,
        })?;
        println!("MMFF94 needs_more={}", result.needs_more());
        println!(
            "optimized first atom: {:?}",
            result.molecule.conformers_3d()[0].coordinates()[0]
        );
    }

    Ok(())
}
