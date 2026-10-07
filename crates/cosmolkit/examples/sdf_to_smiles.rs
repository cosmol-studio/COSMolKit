use std::env;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use cosmolkit::SdfRecordStream;
use std::io::Write;

// Usage:
//   cargo run -p cosmolkit --example sdf_to_smiles -- path/to/input.sdf
fn main() {
    let path = env::args_os().nth(1).map(PathBuf::from).unwrap_or_else(|| {
        panic!("usage: cargo run -p cosmolkit --example sdf_to_smiles -- <file.sdf>")
    });

    let file =
        File::open(&path).unwrap_or_else(|err| panic!("failed to open {}: {err}", path.display()));
    let reader = BufReader::new(file);
    let mut sdf = SdfRecordStream::new(reader);

    let mut found_any = false;
    while let Some(record) = sdf
        .next_record()
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()))
    {
        found_any = true;
        let molecule = record
            .molecule()
            .unwrap_or_else(|err| panic!("failed to obtain concrete SDF record: {err}"));
        let smiles = molecule
            .to_smiles()
            .unwrap_or_else(|err| panic!("failed to write SMILES for record: {err}"));
        let mut stdout = std::io::stdout().lock();
        if let Some(name) = molecule.properties().name().filter(|name| !name.is_empty()) {
            stdout.write_all(name.as_bytes()).expect("write SDF title");
            stdout.write_all(b"\t").expect("write field separator");
        }
        stdout.write_all(smiles.as_bytes()).expect("write SMILES");
        stdout.write_all(b"\n").expect("write record separator");
    }

    if !found_any {
        panic!("no SDF records found in {}", path.display());
    }
}
