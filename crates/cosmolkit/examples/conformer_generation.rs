use cosmolkit::{EmbedParams, MmffOptimizationParams, Molecule, UffOptimizationParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CC(=O)NC")?
        .with_hydrogens()?
        .sanitize()?;

    let mut single_params = EmbedParams::etkdg_v3();
    single_params.random_seed = 0xF00D;
    single_params.num_threads = 1;
    single_params.track_failures = true;

    let embedded = molecule.with_3d_conformer_with_params(&single_params)?;
    println!("single conformers={}", embedded.conformers_3d().len());
    println!(
        "single first atom={:?}",
        embedded.conformers_3d()[0].coordinates()[0]
    );

    let mut multi_params = EmbedParams::etkdg();
    multi_params.random_seed = 123;
    multi_params.num_threads = 1;
    multi_params.prune_rms_thresh = 0.5;
    multi_params.enable_sequential_random_seeds = true;

    let multi = molecule.with_3d_conformers_with_params(5, &multi_params)?;
    println!("pruned conformers={}", multi.conformers_3d().len());

    if embedded.uff_has_all_molecule_params()? {
        let result = embedded.with_uff_optimized_with_params(&UffOptimizationParams {
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
    }

    if embedded.mmff_has_all_molecule_params()? {
        let result = embedded.with_mmff_optimized_with_params(&MmffOptimizationParams {
            mmff_variant: "MMFF94".into(),
            max_iterations: 200,
            non_bonded_threshold: 100.0,
            conformer_id: None,
            ignore_interfragment_interactions: true,
        })?;
        println!("MMFF94 needs_more={}", result.needs_more());
    }

    Ok(())
}
