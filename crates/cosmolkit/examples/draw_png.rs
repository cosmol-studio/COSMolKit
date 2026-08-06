use std::fs;
use std::path::PathBuf;

use cosmolkit::Molecule;

// Usage:
//   cargo run -p cosmolkit --example draw_png -- "CCO"
//   cargo run -p cosmolkit --example draw_png -- "CCO" tmp/cosmolkit_preview.png
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let smiles = args.next().unwrap_or_else(|| "N#C".to_string());
    let output = args
        .next()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("tmp/cosmolkit_preview.png"));

    let mol = Molecule::from_smiles(&smiles)?;
    let png = mol.to_png(1000, 1000)?;
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(&output, png)?;
    println!("wrote {}", output.display());
    Ok(())
}
