use std::env;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use cosmolkit_core::io::sdf::SdfReader;

fn main() {
    let path = env::args_os().nth(1).map(PathBuf::from).unwrap_or_else(|| {
        panic!("usage: cargo run -p cosmolkit-core --example sdf_to_smiles -- <file.sdf>")
    });

    let file =
        File::open(&path).unwrap_or_else(|err| panic!("failed to open {}: {err}", path.display()));
    let reader = BufReader::new(file);
    let mut sdf = SdfReader::new(reader);

    let mut found_any = false;
    while let Some(record) = sdf
        .next_record()
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()))
    {
        found_any = true;
        let smiles = record
            .molecule
            .to_smiles(true)
            .unwrap_or_else(|err| panic!("failed to write SMILES for {}: {err}", record.title));
        if record.title.is_empty() {
            println!("{smiles}");
        } else {
            println!("{}\t{}", record.title, smiles);
        }
    }

    if !found_any {
        panic!("no SDF records found in {}", path.display());
    }
}
