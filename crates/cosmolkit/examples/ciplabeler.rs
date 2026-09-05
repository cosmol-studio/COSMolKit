use cosmolkit::{CipLabelOptions, Molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("C[C@H](F)Cl")?;
    let labeled = molecule.with_cip_labels()?;
    println!(
        "full descriptor: {:?}",
        labeled.atoms()[1].cip_descriptor()?
    );

    let mut selected = Molecule::from_smiles("C[C@H](F)Cl")?;
    selected.assign_cip_labels_with_options_(
        CipLabelOptions::default().with_atoms([1].into_iter().map(cosmolkit::AtomId::new)),
    )?;
    println!(
        "selected descriptor: {:?}",
        selected.atoms()[1].cip_descriptor()?
    );

    let mut double_bond = Molecule::from_smiles("F/C=C/F")?;
    double_bond.assign_cip_labels_()?;
    println!(
        "double-bond descriptor: {:?}",
        double_bond.bonds()[1].cip_descriptor()?
    );
    Ok(())
}
