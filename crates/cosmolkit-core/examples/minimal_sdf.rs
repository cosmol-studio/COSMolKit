use std::env;
use std::fs;
use std::path::PathBuf;

use cosmolkit_core::Molecule;
use cosmolkit_core::io::molblock::mol_to_sdf_record_minimal;

fn main() {
    let args: Vec<String> = env::args().collect();

    // Usage:
    //   cargo run -p cosmolkit-core --example minimal_sdf
    //   cargo run -p cosmolkit-core --example minimal_sdf -- "CCO"
    //   cargo run -p cosmolkit-core --example minimal_sdf -- "CCO" /tmp/out.sdf
    let smiles_inputs: Vec<String> = if args.len() >= 2 {
        vec![args[1].clone()]
    } else {
        vec!["C".to_owned(), "CC".to_owned(), "CCO".to_owned()]
    };

    let mut output = String::new();
    for smiles in smiles_inputs {
        let mut mol = Molecule::from_smiles(&smiles)
            .unwrap_or_else(|err| panic!("failed to parse SMILES '{}': {}", smiles, err));
        mol.compute_2d_coords().unwrap_or_else(|err| {
            panic!("failed to compute 2D coordinates for '{}': {}", smiles, err)
        });
        let sdf = mol_to_sdf_record_minimal(&mol)
            .unwrap_or_else(|err| panic!("failed to write SDF for '{}': {}", smiles, err));
        output.push_str(&sdf);
    }

    if args.len() >= 3 {
        let out_path = PathBuf::from(&args[2]);
        fs::write(&out_path, output)
            .unwrap_or_else(|err| panic!("failed to write file '{}': {}", out_path.display(), err));
        println!("Wrote minimal SDF to {}", out_path.display());
    } else {
        print!("{}", output);
    }
}
