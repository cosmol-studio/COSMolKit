use cosmolkit_core::{
    AtomSpec, Element, EmbedParameters, Molecule, MoleculeBuilder, OperationError,
    SemanticPrecondition, SemanticPreconditionSet, TopologyTrust, WITH_3D_CONFORMER_SPEC,
    WITH_3D_CONFORMERS_SPEC, WITH_HYDROGENS_SPEC,
};

fn xyz_water() -> Molecule {
    Molecule::from_xyz_block(
        "3\nwater\nO 0.000 0.000 0.000\nH 0.958 0.000 0.000\nH -0.239 0.927 0.000\n",
    )
    .expect("xyz water should parse")
}

#[test]
fn with_hydrogens_registry_declares_semantic_preconditions() {
    let preconditions = WITH_HYDROGENS_SPEC.semantic_preconditions;

    assert!(preconditions.contains(SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY));
    assert!(preconditions.contains(SemanticPreconditionSet::HYDROGEN_OWNERSHIP_REPRESENTED));
    assert!(preconditions.contains(SemanticPreconditionSet::VALENCE_COMPUTABLE));
}

#[test]
fn conformer_generation_registry_declares_trusted_topology_precondition() {
    assert!(
        WITH_3D_CONFORMER_SPEC
            .semantic_preconditions
            .contains(SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY)
    );
    assert!(
        WITH_3D_CONFORMERS_SPEC
            .semantic_preconditions
            .contains(SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY)
    );
}

#[test]
fn xyz_coordinate_only_molecule_rejects_with_hydrogens_before_mutation() {
    let molecule = xyz_water();

    assert_eq!(molecule.topology_trust(), TopologyTrust::CoordinateOnly);
    let err = molecule
        .with_hydrogens()
        .expect_err("coordinate-only XYZ topology cannot support AddHs");

    assert!(matches!(
        err,
        OperationError::Precondition {
            requirement: SemanticPrecondition::TrustedBondTopology,
            ..
        }
    ));
    assert_eq!(molecule.num_atoms(), 3);
    assert_eq!(molecule.num_bonds(), 0);
}

#[test]
fn xyz_coordinate_only_molecule_rejects_3d_conformer_generation_before_mutation() {
    let molecule = xyz_water();
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 0xF00D;
    params.num_threads = 1;

    let err = molecule
        .with_3d_conformer_with_params(params)
        .expect_err("coordinate-only XYZ topology cannot support ETKDG");

    assert!(matches!(
        err,
        OperationError::Precondition {
            requirement: SemanticPrecondition::TrustedBondTopology,
            ..
        }
    ));
    assert_eq!(molecule.conformers_3d().len(), 1);
}

#[test]
fn trusted_smiles_graph_allows_with_hydrogens() {
    let molecule = Molecule::from_smiles("CCO").expect("smiles should parse");

    assert_eq!(molecule.topology_trust(), TopologyTrust::TrustedGraph);
    let with_hydrogens = molecule
        .with_hydrogens()
        .expect("trusted SMILES topology supports AddHs");

    assert_eq!(with_hydrogens.num_atoms(), 9);
    assert_eq!(with_hydrogens.num_bonds(), 8);
    assert_eq!(with_hydrogens.topology_trust(), TopologyTrust::TrustedGraph);
}

#[test]
fn trusted_smiles_graph_without_explicit_hydrogens_allows_3d_conformer_generation() {
    let molecule = Molecule::from_smiles("Oc1ccccc1-c2ccccc2O").expect("smiles should parse");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 0xF00D;
    params.num_threads = 1;

    assert_eq!(molecule.topology_trust(), TopologyTrust::TrustedGraph);
    assert_eq!(
        molecule
            .atoms()
            .iter()
            .filter(|atom| atom.atomic_number() == 1)
            .count(),
        0
    );

    let embedded = molecule
        .with_3d_conformer_with_params(params)
        .expect("RDKit allows heavy-atom-only ETKDG with a warning");

    assert_eq!(embedded.num_atoms(), molecule.num_atoms());
    assert_eq!(embedded.conformers_3d().len(), 1);
    assert_eq!(
        embedded.conformers_3d()[0].coordinates().len(),
        molecule.num_atoms()
    );
}

#[test]
fn trusted_graph_without_explicit_hydrogen_ownership_rejects_with_hydrogens() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::O));
    builder.add_atom(AtomSpec::new(Element::H));
    builder.add_atom(AtomSpec::new(Element::H));
    let molecule = builder
        .build()
        .expect("builder should produce a valid graph");

    assert_eq!(molecule.topology_trust(), TopologyTrust::TrustedGraph);
    let err = molecule
        .with_hydrogens()
        .expect_err("unowned explicit H atoms cannot support AddHs");

    assert!(matches!(
        err,
        OperationError::Precondition {
            requirement: SemanticPrecondition::HydrogenOwnershipRepresented,
            ..
        }
    ));
}

#[test]
fn pickle_roundtrip_preserves_coordinate_only_topology_capability() {
    let molecule = xyz_water();
    let data = cosmolkit_core::mol_to_binary(&molecule).expect("pickle should serialize");
    let roundtrip = cosmolkit_core::mol_from_binary(&data).expect("pickle should deserialize");

    assert_eq!(roundtrip.topology_trust(), TopologyTrust::CoordinateOnly);
    assert!(matches!(
        roundtrip.with_hydrogens(),
        Err(OperationError::Precondition {
            requirement: SemanticPrecondition::TrustedBondTopology,
            ..
        })
    ));
}
