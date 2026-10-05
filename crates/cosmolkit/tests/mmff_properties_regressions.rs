use cosmolkit::{
    AtomSpec, Conformer3D, CoordinateBlock, Element, MmffPropertiesParams, Molecule,
    MoleculeBuilder,
};
fn mmff_has_all_molecule_params(mol: &Molecule) -> Result<bool, cosmolkit::MmffMolPropertiesError> {
    mol.mmff_has_all_molecule_params()
}
fn raw_smiles(input: &str) -> Molecule {
    Molecule::from_smiles_with_params(
        input,
        &cosmolkit::SmilesParseParams {
            sanitize: false,
            ..Default::default()
        },
    )
    .unwrap()
}
fn single_atom_with_3d_conformer(atom: AtomSpec, coords: [f64; 3]) -> Molecule {
    let mut b = MoleculeBuilder::new();
    b.add_atom(atom);
    b.add_3d_conformer(vec![coords]).unwrap();
    b.build().unwrap()
}
#[test]
fn mmff_mol_properties_constructor_assigns_source_atom_types_for_ethanol() {
    let molecule = raw_smiles("CCO");

    let props = molecule
        .mmff_properties()
        .expect("ethanol source atom typing is ported");

    assert!(props.is_valid());
    assert_eq!(props.atom_type(0).unwrap(), 1);
    assert_eq!(props.atom_type(1).unwrap(), 1);
    assert_eq!(props.atom_type(2).unwrap(), 6);
}

#[test]
fn mmff_mol_properties_constructor_types_all_explicit_ethene_hydrogens() {
    let molecule = Molecule::from_smiles("C=C")
        .expect("ethene parses")
        .with_hydrogens()
        .expect("ethene AddHs succeeds");
    let props = molecule
        .mmff_properties()
        .expect("explicit-H ethene MMFF typing succeeds");

    assert!(props.is_valid());
    assert_eq!(
        props
            .atoms()
            .iter()
            .map(|properties| properties.atom_type)
            .collect::<Vec<_>>(),
        vec![2, 2, 5, 5, 5, 5]
    );
    assert!(
        mmff_has_all_molecule_params(&molecule)
            .expect("explicit-H ethene parameter coverage computes")
    );
}

#[test]
fn mmff_public_api_mmff_has_all_molecule_params_returns_true_for_empty_molecule() {
    let molecule = Molecule::new();

    let found_all = mmff_has_all_molecule_params(&molecule)
        .expect("empty MMFF parameter coverage should compute");

    assert!(found_all);
}

#[test]
fn mmff_public_api_mmff_has_all_molecule_params_returns_false_for_missing_atom_params() {
    let molecule = single_atom_with_3d_conformer(AtomSpec::new(Element::HE), [0.0, 0.0, 0.0]);

    let found_all = mmff_has_all_molecule_params(&molecule)
        .expect("missing MMFF atom typing should be reported as false");

    assert!(!found_all);
}

#[test]
fn mmff_public_api_mmff_has_all_molecule_params_uses_rdkit_default_mmff94_constructor_path() {
    let molecule = Molecule::new();
    let props = molecule
        .mmff_properties()
        .expect("empty MMFF constructor path should succeed");

    let found_all = mmff_has_all_molecule_params(&molecule)
        .expect("empty MMFF parameter coverage should compute");

    assert_eq!(found_all, props.is_valid());
    assert_eq!(props.variant(), cosmolkit::MmffVariant::Mmff94);
}

#[test]
fn mmff_mol_properties_types_neutral_nitric_oxide_like_rdkit() {
    let molecule = Molecule::from_smiles("[N]=O").expect("neutral nitric oxide parses");

    assert!(
        mmff_has_all_molecule_params(&molecule)
            .expect("neutral nitric oxide MMFF coverage computes")
    );
    let props = molecule
        .mmff_properties()
        .expect("neutral nitric oxide MMFF properties compute");

    assert_eq!(props.atom_type(0).unwrap(), 8);
    assert_eq!(props.atom_type(1).unwrap(), 7);
    assert_eq!(props.formal_charge(0).unwrap(), 0.0);
    assert_eq!(props.formal_charge(1).unwrap(), 0.0);
    assert!((props.partial_charge(0).unwrap() - 0.43400000000000005).abs() <= 1.0e-12);
    assert!((props.partial_charge(1).unwrap() - -0.43400000000000005).abs() <= 1.0e-12);
}

#[test]
fn mmff_mol_properties_dearomatizes_complete_fused_quinone_ring_like_rdkit() {
    let molecule = Molecule::from_smiles(
        "CC1(O)O[C@@H]2c3c(ccc(=O)c(O)c3O)[C@H]1[C@]1(C)OCc3c(ccc(O)c3O)[C@H]21",
    )
    .expect("fused quinone regression molecule parses");
    let props = molecule
        .mmff_properties()
        .expect("fused quinone MMFF properties compute");

    assert_eq!(
        props
            .atoms()
            .iter()
            .map(|properties| properties.atom_type)
            .collect::<Vec<_>>(),
        vec![
            1, 1, 6, 6, 1, 2, 2, 2, 2, 3, 7, 2, 6, 2, 6, 1, 1, 1, 6, 1, 37, 37, 37, 37, 37, 6, 37,
            6, 1,
        ]
    );
}

#[test]
fn mmff_mol_properties_falls_through_aromatic_six_ring_sulfur_like_rdkit() {
    let molecule =
        Molecule::from_smiles("[s+]1ccccc1.[Cl-]").expect("thiopyrylium chloride parses");
    let props = molecule
        .mmff_properties()
        .expect("thiopyrylium chloride MMFF properties compute");

    assert!(props.is_valid());
    assert_eq!(
        props
            .atoms()
            .iter()
            .map(|properties| properties.atom_type)
            .collect::<Vec<_>>(),
        vec![15, 37, 37, 37, 37, 37, 90]
    );
    assert_eq!(
        props
            .atoms()
            .iter()
            .map(|properties| properties.formal_charge)
            .collect::<Vec<_>>(),
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    );
    assert_eq!(
        props
            .atoms()
            .iter()
            .map(|properties| properties.partial_charge)
            .collect::<Vec<_>>(),
        vec![-0.203, 0.1015, 0.0, 0.0, 0.0, 0.1015, -1.0]
    );
}

#[test]
fn mmff_mol_properties_types_chembl_aromatic_phosphorus_case_like_rdkit() {
    let molecule = Molecule::from_smiles(
        "CC(=O)C1(C(C)=O)C(Cl)C(=O)N1N(c1c(O)ccc2c(P(Cl)Cl)pc(C(=O)O)n12)[N+](=O)[O-]",
    )
    .expect("ChEMBL aromatic phosphorus regression molecule parses");
    let props = molecule
        .mmff_properties()
        .expect("ChEMBL aromatic phosphorus MMFF properties compute");

    assert!(props.is_valid());
    assert_eq!(
        props
            .atoms()
            .iter()
            .map(|properties| properties.atom_type)
            .collect::<Vec<_>>(),
        vec![
            1, 3, 7, 20, 3, 1, 7, 20, 12, 3, 7, 10, 40, 2, 2, 6, 2, 2, 63, 64, 26, 12, 12, 75, 63,
            3, 7, 6, 39, 45, 32, 32,
        ]
    );
    assert_eq!(props.partial_charge(20).unwrap(), 0.4614);
    assert_eq!(props.partial_charge(23).unwrap(), -0.14900000000000002);
}

#[test]
fn mmff_mol_properties_atom_types_terminal_s_double_bonded_to_carbon_as_source_type_16() {
    let molecule = Molecule::from_smiles("C=S").unwrap();

    let props = molecule
        .mmff_properties()
        .expect("thiocarbonyl sulfur source branch is ported");
    let sulfur_idx = molecule
        .atoms()
        .iter()
        .position(|atom| atom.atomic_number() == 16)
        .expect("test molecule has sulfur");

    assert_eq!(props.atom_type(sulfur_idx).unwrap(), 16);
}

#[test]
fn mmff_mol_properties_constructor_preserves_existing_sanitized_prop_on_empty_molecule() {
    let mut properties = cosmolkit::MoleculeProperties::default();
    properties.set_prop("_MMFFSanitized", "already").unwrap();
    let molecule =
        Molecule::from_parts(Default::default(), Default::default(), properties).unwrap();
    let original = molecule.clone();
    let props = molecule.mmff_properties().unwrap();
    assert!(props.is_valid());
    assert_eq!(props.variant(), cosmolkit::MmffVariant::Mmff94);
    assert!(props.atoms().is_empty());
    assert_eq!(
        molecule.properties().prop("_MMFFSanitized"),
        Some("already")
    );
    assert!(!molecule.properties().is_prop_computed("_MMFFSanitized"));
    assert_eq!(molecule, original);
}
