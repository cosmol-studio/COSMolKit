mod common;

use std::io;

use common::MethaneToolkit;
use cosmolkit_inchi::mol_from_inchi;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output = mol_from_inchi(&mut MethaneToolkit, b"InChI=1S/CH4/h1H4", false, false)?;
    let molecule = output
        .molecule
        .ok_or_else(|| io::Error::other("InChI source returned no molecular graph"))?;

    println!("status: {}", output.return_values.return_code);
    println!("atoms: {}", molecule.atoms().len());
    println!("bonds: {}", molecule.bonds().len());
    for (index, atom) in molecule.atoms().iter().enumerate() {
        println!(
            "atom {index}: Z={} charge={} explicit_hydrogens={}",
            atom.atomic_number, atom.formal_charge, atom.num_explicit_hydrogens
        );
    }
    Ok(())
}
