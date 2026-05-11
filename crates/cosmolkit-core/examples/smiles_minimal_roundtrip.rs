use cosmolkit_core::{Molecule, SmilesWriteParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CCO")?;
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: false,
        all_hydrogens_explicit: false,
        do_random: false,
        rooted_at_atom: None,
        include_dative_bonds: true,
        ignore_atom_map_numbers: false,
    };

    let smiles = molecule.to_smiles_with_params(&params)?;
    assert_eq!(smiles, "CCO");

    println!("{smiles}");
    Ok(())
}
