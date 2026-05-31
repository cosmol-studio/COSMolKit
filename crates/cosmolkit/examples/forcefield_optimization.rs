use cosmolkit::{
    Molecule, mmff_has_all_molecule_params, mmff_optimize_molecule, uff_has_all_molecule_params,
    uff_optimize_molecule,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CCO")?.with_hydrogens()?.sanitize()?;

    let mut builder = molecule.to_builder();
    builder.add_3d_conformer(vec![
        [0.000, 0.000, 0.000],
        [1.540, 0.000, 0.000],
        [2.100, 1.200, 0.000],
        [-0.600, 0.900, 0.000],
        [-0.600, -0.900, 0.000],
        [0.000, 0.000, 1.000],
        [1.900, -0.900, 0.000],
        [1.700, 0.000, 1.000],
        [2.900, 1.200, 0.000],
    ])?;
    let molecule = builder.build()?;

    if uff_has_all_molecule_params(&molecule)? {
        let result = uff_optimize_molecule(&molecule, 200, 10.0, -1, true)?;
        println!(
            "UFF needs_more={} energy={:.6}",
            result.needs_more, result.energy
        );
        println!(
            "optimized first atom: {:?}",
            result.molecule.conformers_3d()[0].coords()[0]
        );
    }

    if mmff_has_all_molecule_params(&molecule)? {
        let result = mmff_optimize_molecule(&molecule, "MMFF94", 200, 100.0, -1, true)?;
        println!("MMFF94 needs_more={}", result.needs_more);
        println!(
            "optimized first atom: {:?}",
            result.molecule.conformers_3d()[0].coords()[0]
        );
    }

    Ok(())
}
