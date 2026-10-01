#[test]
fn default_selects_full() {
    assert!(include_str!("../Cargo.toml").contains("default = [\"full\"]"));
}

#[cfg(feature = "full")]
#[test]
fn full_enables_all_architectural_capabilities() {
    assert!(cfg!(all(
        feature = "cap-alignment",
        feature = "cap-batch",
        feature = "cap-bio",
        feature = "cap-conformer",
        feature = "cap-confseq",
        feature = "cap-depict",
        feature = "cap-fingerprints",
        feature = "cap-forcefields",
        feature = "cap-hashing",
        feature = "cap-inchi",
        feature = "cap-io",
        feature = "cap-search",
        feature = "cap-serialization",
        feature = "cap-smiles",
        feature = "cap-stereoisomers",
        feature = "cap-tautomer",
        feature = "cap-hydrogens",
        feature = "cap-descriptors",
        feature = "cap-valence",
        feature = "cap-radicals",
        feature = "cap-rings",
        feature = "cap-matrices",
        feature = "cap-transforms",
        feature = "cap-stereo",
        feature = "cap-kekulize",
        feature = "cap-aromaticity",
        feature = "cap-sanitize",
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
