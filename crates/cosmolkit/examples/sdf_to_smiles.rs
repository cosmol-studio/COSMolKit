use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit::SdfRecord;

// Usage:
//   cargo run -p cosmolkit --example sdf_to_smiles -- path/to/input.sdf
fn print_record(text: &str) -> Result<(), Box<dyn std::error::Error>> {
    let record = SdfRecord::from_sdf(text)?;
    let molecule = record.molecule()?;
    let smiles = molecule.to_smiles()?;
    match molecule.properties().name().filter(|name| !name.is_empty()) {
        Some(name) => println!("{name}\t{smiles}"),
        None => println!("{smiles}"),
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .ok_or("usage: cargo run -p cosmolkit --example sdf_to_smiles -- <file.sdf>")?;
    let mut reader = BufReader::new(File::open(path)?);
    let mut text = String::new();
    let mut line = String::new();
    let mut found_any = false;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        text.push_str(&line);
        if line.trim_end_matches(['\r', '\n']) == "$$$$" {
            print_record(&text)?;
            found_any = true;
            text.clear();
        }
    }
    if !text.is_empty() {
        print_record(&text)?;
        found_any = true;
    }
    if !found_any {
        return Err("no SDF records found".into());
    }
    Ok(())
}
