use std::fs;
use std::path::PathBuf;

use cosmolkit_core::Molecule;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let smiles = std::env::args().nth(1).unwrap_or_else(|| "N#C".to_string());

    // Usage:
    //   cargo run -p cosmolkit-core --example draw_png -- "CCO"
    let output = std::env::args()
        .nth(2)
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
