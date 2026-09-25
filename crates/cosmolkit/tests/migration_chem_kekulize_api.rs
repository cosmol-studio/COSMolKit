#[path = "support/coordinate_views.rs"]
mod coordinate_views;

use std::error::Error as _;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec, Conformer2D,
    CoordinateBlock, Element, KekulizeError, KekulizeParams, Molecule, MoleculeOpKind,
    MoleculeOpOutput, MoleculeProperties, OperationDomain, OperationError, ParityPolicy,
    StateModel, StereoGroup, StereoGroupKind, SupportStatus, TopologyBlock, TopologyEditKind,
    feature_spec, operation_invariant, operation_parity, operation_spec, support_matrix,
};

fn aromatic_carbon_cycle(size: usize) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        (0..size)
            .map(|index| {
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C)
                        .with_aromatic(true)
                        .with_prop("atom-note", format!("atom-{index}"))
                        .unwrap()
                        .with_computed_prop("_CIPCode", "R")
                        .unwrap(),
                )
            })
            .collect(),
        (0..size)
            .map(|index| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(
                        AtomId::new(index),
                        AtomId::new((index + 1) % size),
                        BondOrder::Aromatic,
                    )
                    .with_aromatic(true)
                    .with_prop("bond-note", format!("bond-{index}"))
                    .unwrap()
                    .with_computed_prop("_CIPCode", "E")
                    .unwrap(),
                )
            })
            .collect(),
        Vec::new(),
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(0)],
            vec![BondId::new(0)],
        )],
    )
    .unwrap()
}

fn aromatic_molecule(size: usize) -> Molecule {
    Molecule::from_parts(
        aromatic_carbon_cycle(size),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                17,
                (0..size)
                    .map(|index| [index as f64, -(index as f64)])
                    .collect(),
            )],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default()
            .with_name("kekulize-public")
            .with_prop("source", "preserved")
            .unwrap()
            .with_computed_prop("_CIPComputed", "true")
            .unwrap(),
    )
    .unwrap()
}

#[test]
fn canonical_public_signatures_and_defaults_compile() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::with_kekulized_bonds;
    let _: for<'a, 'b> fn(&'a Molecule, &'b KekulizeParams) -> Result<Molecule, OperationError> =
        Molecule::with_kekulized_bonds_with_params;
    let _: fn(&mut Molecule) -> Result<(), OperationError> = Molecule::kekulize_bonds_;
    let _: for<'a, 'b> fn(&'a mut Molecule, &'b KekulizeParams) -> Result<(), OperationError> =
        Molecule::kekulize_bonds_with_params_;
    assert_eq!(
        KekulizeParams::default(),
        KekulizeParams {
            mark_atoms_bonds: true,
            canonical: true,
            max_backtracks: 100,
        }
    );
}

#[test]
fn binding_contract_exposes_exactly_the_frozen_six_entries() {
    let expected = [
        "types.KekulizeParams",
        "types.KekulizeError",
        "Molecule.with_kekulized_bonds",
        "Molecule.with_kekulized_bonds_with_params",
        "Molecule.kekulize_bonds_",
        "Molecule.kekulize_bonds_with_params_",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "kekulize")
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    for row in &rows[..2] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    for row in &rows[2..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    assert_eq!(
        rows[2].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(
        rows[3].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(rows[4].callable.unwrap().state_model, StateModel::InPlace);
    assert_eq!(rows[5].callable.unwrap().state_model, StateModel::InPlace);
}

#[test]
fn generated_registry_and_all_four_matrices_share_one_exact_operation() {
    let feature = feature_spec("kekulize").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_kekulized_bonds_with_params").unwrap();
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::Local);
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.access.read(), BlockSet::NONE);
    assert_eq!(
        spec.access.write(),
        BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(spec.may_mutate, spec.access.write());
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.derived_effects.recompute.bits(), 0);
    assert_eq!(
        spec.derived_effects.preserve.bits(),
        (1 << 0) | (1 << 1) | (1 << 5)
    );
    assert_eq!(
        spec.derived_effects.invalidate.bits(),
        (1 << 2) | (1 << 3) | (1 << 4) | (1 << 6) | (1 << 7)
    );
    assert_eq!(format!("{:?}", spec.cip_state), "Preserve");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "weak_kekulize_bond_assignment"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "kekulize_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| row.feature.name == "kekulize")
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));
}

#[test]
fn value_operation_preserves_identity_coordinates_stereo_and_ordinary_props() {
    let source = aromatic_molecule(6);
    let source_snapshot = source.clone();
    let output = source.with_kekulized_bonds().unwrap();

    assert_eq!(source, source_snapshot);
    assert_eq!(source.num_atoms(), output.num_atoms());
    assert_eq!(source.num_bonds(), output.num_bonds());
    assert!(
        source
            .atoms()
            .iter()
            .zip(output.atoms())
            .all(|(before, after)| before.id() == after.id())
    );
    assert!(
        source
            .bonds()
            .iter()
            .zip(output.bonds())
            .all(|(before, after)| before.id() == after.id())
    );
    assert_eq!(source.topology().adjacency, output.topology().adjacency);
    assert_eq!(
        source.topology().stereo_groups,
        output.topology().stereo_groups
    );
    coordinate_views::assert_shared_coordinates(&source, &output);
    assert_eq!(source.coordinates_2d(), output.coordinates_2d());
    assert_eq!(source.conformers_3d(), output.conformers_3d());
    assert_eq!(output.property("source"), Some("preserved"));
    assert_eq!(output.atoms()[0].prop("atom-note"), Some("atom-0"));
    assert_eq!(output.bonds()[0].prop("bond-note"), Some("bond-0"));
    assert_eq!(output.property("_CIPComputed"), Some("true"));
    assert_eq!(output.atoms()[0].prop("_CIPCode"), Some("R"));
    assert_eq!(output.bonds()[0].prop("_CIPCode"), Some("E"));
    assert!(output.atoms().iter().all(|atom| !atom.is_aromatic()));
    assert!(output.bonds().iter().all(|bond| !bond.is_aromatic()));
    assert_eq!(
        output
            .bonds()
            .iter()
            .filter(|bond| bond.order() == BondOrder::Double)
            .count(),
        3
    );
    assert!(format!("{output:?}").contains("derived_cache_is_empty: true"));

    assert!(source.atoms().iter().all(Atom::is_aromatic));
    assert!(source.bonds().iter().all(Bond::is_aromatic));
    assert_eq!(source.property("_CIPComputed"), Some("true"));
    assert!(!std::ptr::eq(source.topology(), output.topology()));
    assert!(std::ptr::eq(source.properties(), output.properties()));
}

#[test]
fn parameterized_and_inplace_forms_match_their_frozen_semantics() {
    let source = aromatic_molecule(6);
    let default_output = source.with_kekulized_bonds().unwrap();
    let explicit_default = source
        .with_kekulized_bonds_with_params(&KekulizeParams::default())
        .unwrap();
    assert_eq!(default_output, explicit_default);

    let keep_flags = KekulizeParams {
        mark_atoms_bonds: false,
        canonical: false,
        max_backtracks: 100,
    };
    let expected = source
        .with_kekulized_bonds_with_params(&keep_flags)
        .unwrap();
    assert!(expected.atoms().iter().all(Atom::is_aromatic));
    assert!(expected.bonds().iter().all(Bond::is_aromatic));
    assert_eq!(
        expected
            .bonds()
            .iter()
            .filter(|bond| bond.order() == BondOrder::Double)
            .count(),
        3
    );

    let mut short = source.clone();
    short.kekulize_bonds_().unwrap();
    assert_eq!(short, default_output);
    let mut explicit = source.clone();
    explicit.kekulize_bonds_with_params_(&keep_flags).unwrap();
    assert_eq!(explicit, expected);
    assert_eq!(source, aromatic_molecule(6));
}

#[test]
fn typed_algorithm_failure_is_atomic_for_value_and_inplace_wrappers() {
    let source = aromatic_molecule(5);
    let snapshot = source.clone();
    let value_error = source.with_kekulized_bonds().unwrap_err();
    assert!(matches!(
        &value_error,
        OperationError::Kekulize(KekulizeError::NotKekulizable { problem_atoms })
            if problem_atoms == &(0..5).map(AtomId::new).collect::<Vec<_>>()
    ));
    assert!(value_error.source().is_some());
    assert_eq!(source, snapshot);

    let mut target = source.clone();
    let observer = target.clone();
    let inplace_error = target.kekulize_bonds_().unwrap_err();
    assert_eq!(inplace_error, value_error);
    assert_eq!(target, source);
    assert!(std::ptr::eq(target.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&target, &observer);
    assert!(std::ptr::eq(target.properties(), observer.properties()));
    assert_eq!(target.property("_CIPComputed"), Some("true"));
}

#[test]
fn single_output_local_edit_makes_multi_candidate_ordering_not_applicable() {
    let spec = operation_spec("with_kekulized_bonds_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
    assert_eq!(spec.topology_edit, TopologyEditKind::Local);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.auto_remap, BlockSet::NONE);
}
