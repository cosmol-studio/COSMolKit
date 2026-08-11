mod common;

use common::MethaneToolkit;
use cosmolkit_inchi::{InchiAtom, InchiMolecule, mol_to_inchi};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let methane = InchiMolecule::try_from_graph(
        vec![InchiAtom {
            atomic_number: 6,
            ..InchiAtom::default()
        }],
        Vec::new(),
        Vec::new(),
    )
    .expect("the static methane graph has valid atom and bond indices");

    let output = mol_to_inchi(&mut MethaneToolkit, &methane, None)?;
    println!("status: {}", output.return_values.return_code);
    println!("InChI: {}", String::from_utf8(output.inchi)?);
    Ok(())
}
