use std::env;
use std::fs;

use cosmolkit::Molecule;

const DEMO_XYZ: &str = "3
water
O 0.000 0.000 0.000
H 0.758 0.000 0.504
H -0.758 0.000 0.504
";

// Usage:
//   cargo run -p cosmolkit --example read_xyz
//   cargo run -p cosmolkit --example read_xyz -- path/to/molecule.xyz
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let xyz = match env::args().nth(1) {
        Some(path) => fs::read_to_string(path)?,
        None => DEMO_XYZ.to_string(),
    };

    let mol = Molecule::from_xyz_block(&xyz)?;
    println!("atoms: {}", mol.num_atoms());
    println!("bonds: {}", mol.num_bonds());
    println!("3d conformers: {}", mol.conformers_3d().len());

    if let Some(conformer) = mol.conformers_3d().first() {
        for (idx, (atom, coord)) in mol.atoms().iter().zip(conformer.coordinates()).enumerate() {
            println!(
                "{idx}: Z={} x={:.3} y={:.3} z={:.3}",
                atom.atomic_number(),
                coord[0],
                coord[1],
                coord[2]
            );
        }
    }

    Ok(())
}
