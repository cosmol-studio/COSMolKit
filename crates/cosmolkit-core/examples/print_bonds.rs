use cosmolkit_core::Molecule;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let smiles = std::env::args().nth(1).unwrap();
    let mut mol = Molecule::from_smiles(&smiles)?;
    let _ = cosmolkit_core::kekulize::kekulize_in_place(&mut mol, false);
    mol.compute_2d_coords()?;
    for (i, pt) in mol.coords_2d().unwrap().iter().enumerate() {
        if (10..=18).contains(&i) {
            eprintln!("coord {} {:.3} {:.3}", i, pt.x, pt.y);
        }
    }
    for bond in &mol.bonds {
        println!(
            "{} {} {} {:?} arom={}",
            bond.index, bond.begin_atom, bond.end_atom, bond.order, bond.is_aromatic
        );
    }
    Ok(())
}
