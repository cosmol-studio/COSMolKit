use std::fs;
use std::path::PathBuf;

use cosmolkit_core::Molecule;

// Usage:
//     cargo run -p cosmolkit-core --example draw_svg -- "CCO"
//     cargo run -p cosmolkit-core --example draw_svg -- "CCO" 400 300 -o tmp/cosmolkit_preview.svg
//     cargo run -p cosmolkit-core --example draw_svg -- "CCO" -o tmp/cosmolkit_preview.svg

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let mut smiles = "N#C".to_string();
    let mut width = 300;
    let mut height = 300;
    let mut output: Option<PathBuf> = None;

    if let Some(arg) = args.next() {
        smiles = arg;
    }
    if let Some(arg) = args.next() {
        if arg == "-o" {
            output = args.next().map(PathBuf::from);
        } else {
            width = arg.parse::<u32>().unwrap_or(300);
            if let Some(arg) = args.next() {
                if arg == "-o" {
                    output = args.next().map(PathBuf::from);
                } else {
                    height = arg.parse::<u32>().unwrap_or(300);
                    if let Some(arg) = args.next() {
                        if arg == "-o" {
                            output = args.next().map(PathBuf::from);
                        }
                    }
                }
            }
        }
    }

    let mol = Molecule::from_smiles(&smiles)?;
    let svg = mol.to_svg(width, height)?;
    if let Some(output) = output {
        if let Some(parent) = output.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(&output, svg)?;
        println!("wrote {}", output.display());
    } else {
        print!("{svg}");
    }
    Ok(())
}
