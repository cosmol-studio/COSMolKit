use cosmolkit::{
    EmbedParameters, Molecule, mmff_has_all_molecule_params, mmff_optimize_molecule,
    uff_has_all_molecule_params, uff_optimize_molecule,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CC(=O)NC")?
        .with_hydrogens()?
        .sanitize()?;

    let mut single_params = EmbedParameters::etkdg_v3();
    single_params.random_seed = 0xF00D;
    single_params.num_threads = 1;
    single_params.track_failures = true;

    let embedded = molecule.with_3d_conformer_with_params(single_params.clone())?;
    println!("single conformers={}", embedded.conformers_3d().len());
    println!(
        "single first atom={:?}",
        embedded.conformers_3d()[0].coords()[0]
    );

    let mut multi_params = EmbedParameters::etkdg();
    multi_params.random_seed = 123;
    multi_params.num_threads = 1;
    multi_params.prune_rms_thresh = 0.5;
    multi_params.enable_sequential_random_seeds = true;

    let multi = molecule.with_3d_conformers_with_params(5, multi_params)?;
    println!("pruned conformers={}", multi.conformers_3d().len());

    if uff_has_all_molecule_params(&embedded)? {
        let result = uff_optimize_molecule(&embedded, 200, 10.0, -1, true)?;
        println!(
            "UFF needs_more={} energy={:.6}",
            result.needs_more, result.energy
        );
    }

    if mmff_has_all_molecule_params(&embedded)? {
        let result = mmff_optimize_molecule(&embedded, "MMFF94", 200, 100.0, -1, true)?;
        println!("MMFF94 needs_more={}", result.needs_more);
    }

    Ok(())
}
