#[cfg(feature = "cap-smiles")]
#[test]
fn smiles_convenience_returns_the_original_value_and_error() {
    let source = "CCO";
    assert_eq!(
        cosmolkit::mol_from_smiles(source)
            .unwrap()
            .to_smiles()
            .unwrap(),
        cosmolkit::Molecule::from_smiles(source)
            .unwrap()
            .to_smiles()
            .unwrap()
    );
    assert_eq!(
        format!("{:?}", cosmolkit::mol_from_smiles("C1CC").unwrap_err()),
        format!(
            "{:?}",
            cosmolkit::Molecule::from_smiles("C1CC").unwrap_err()
        )
    );
}

#[cfg(feature = "cap-bio")]
#[test]
fn bio_conveniences_forward_to_the_same_factories() {
    let pdb =
        "ATOM      1  N   ALA A   1      11.104  13.207   9.900  1.00 20.00           N\nEND\n";
    let structure = cosmolkit::bio_from_pdb(pdb).unwrap();
    assert_eq!(
        structure.num_atoms(),
        cosmolkit::BioStructure::from_pdb(pdb).unwrap().num_atoms()
    );
    let cif = structure.to_mmcif().unwrap();
    assert_eq!(
        cosmolkit::bio_from_mmcif(&cif).unwrap().num_atoms(),
        structure.num_atoms()
    );
    assert_eq!(
        cosmolkit::protein_from_mmcif(&cif).unwrap().num_atoms(),
        cosmolkit::Protein::from_mmcif(&cif).unwrap().num_atoms()
    );
    assert_eq!(
        cosmolkit::protein_from_pdb(pdb).unwrap().num_atoms(),
        cosmolkit::Protein::from_pdb(pdb).unwrap().num_atoms()
    );
}
