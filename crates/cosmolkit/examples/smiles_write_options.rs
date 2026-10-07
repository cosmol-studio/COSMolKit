use cosmolkit::{Molecule, SmilesWriteParams};

// Usage: cargo run -p cosmolkit --example smiles_write_options
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CCO")?;
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: false,
        all_hydrogens_explicit: false,
        rooted_at_atom: None,
        include_dative_bonds: true,
        ignore_atom_map_numbers: false,
    };

    let smiles = molecule.to_smiles_with_params(&params)?;
    assert_eq!(smiles.as_bytes(), b"CCO");

    use std::io::Write;
    let mut stdout = std::io::stdout().lock();
    stdout.write_all(smiles.as_bytes())?;
    stdout.write_all(b"\n")?;
    Ok(())
}
