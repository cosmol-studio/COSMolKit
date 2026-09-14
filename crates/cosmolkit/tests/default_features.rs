#[test]
fn default_selects_full() {
    assert!(include_str!("../Cargo.toml").contains("default = [\"full\"]"));
}

#[cfg(feature = "full")]
#[test]
fn full_enables_all_architectural_capabilities() {
    assert!(cfg!(all(
        feature = "alignment",
        feature = "batch",
        feature = "bio",
        feature = "conformer",
        feature = "confseq",
        feature = "depict",
        feature = "fingerprints",
        feature = "forcefields",
        feature = "hashing",
        feature = "inchi",
        feature = "io",
        feature = "search",
        feature = "serialization",
        feature = "smiles",
        feature = "stereoisomers",
        feature = "tautomer",
        feature = "hydrogens",
        feature = "descriptors",
        feature = "valence",
        feature = "radicals",
        feature = "rings",
        feature = "matrices",
        feature = "transforms",
        feature = "stereo",
        feature = "kekulize",
        feature = "aromaticity",
        feature = "sanitize",
    )));
}

#[cfg(feature = "full")]
#[test]
fn full_exposes_descriptor_methods_on_live_molecules() {
    use cosmolkit::{
        Atom, AtomId, AtomSpec, CoordinateBlock, Element, Molecule, MoleculeProperties,
        TopologyBlock,
    };

    let topology = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(4)
                .with_no_implicit(true),
        )],
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    let molecule = Molecule::from_parts(
        topology,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap();

    assert_eq!(molecule.molecular_formula().unwrap(), "CH4");
    assert_eq!(
        molecule
            .molecular_formula_with_options(false, false)
            .unwrap(),
        "CH4"
    );
    assert!((molecule.molecular_weight().unwrap() - 16.043).abs() < 0.001);
    assert!((molecule.molecular_weight_with_options(true).unwrap() - 12.011).abs() < 0.001);
    assert!((molecule.exact_molecular_weight().unwrap() - 16.0313).abs() < 0.0001);
    assert_eq!(
        molecule.exact_molecular_weight_with_options(true).unwrap(),
        12.0
    );
}
