use std::env;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use cosmolkit_core::io::sdf::SdfReader;

fn main() {
    let path = env::args_os().nth(1).map(PathBuf::from).unwrap_or_else(|| {
        panic!("usage: cargo run -p cosmolkit-core --example sdf_count -- <file.sdf>")
    });

    let file =
        File::open(&path).unwrap_or_else(|err| panic!("failed to open {}: {err}", path.display()));
    let mut reader = SdfReader::new(BufReader::new(file));

    let mut records = 0usize;
    let mut atoms = 0usize;
    let mut bonds = 0usize;
    let mut fields = 0usize;
    while let Some(record) = reader
        .next_record()
        .unwrap_or_else(|err| panic!("failed at record {}: {err}", records + 1))
    {
        records += 1;
        atoms += record.molecule.atoms().len();
        bonds += record.molecule.bonds().len();
        fields += record.data_fields.len();
    }

    println!("records={records} atoms={atoms} bonds={bonds} fields={fields}");
}
