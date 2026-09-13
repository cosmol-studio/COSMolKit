use std::error::Error as _;

use cosmolkit::{
    AromaticityError, AromaticityModel, AromaticityParams, Atom, AtomId, AtomSpec,
    BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner, BindingParity, BindingSupport,
    BlockSet, Bond, BondId, BondOrder, BondSpec, Conformer2D, CoordinateBlock, Element, Molecule,
    MoleculeOpKind, MoleculeOpOutput, MoleculeProperties, OperationDomain, OperationError,
    ParityPolicy, StateModel, StereoGroup, StereoGroupKind, SupportStatus, TopologyBlock,
    TopologyEditKind, feature_spec, operation_invariant, operation_parity, operation_spec,
    support_matrix,
};

fn carbon_cycle(size: usize) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        (0..size)
            .map(|index| {
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C)
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
                        if index % 2 == 0 {
                            BondOrder::Double
                        } else {
                            BondOrder::Single
                        },
                    )
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

fn aromaticity_molecule(size: usize) -> Molecule {
    Molecule::from_parts(
        carbon_cycle(size),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                29,
                (0..size)
                    .map(|index| [index as f64, -(index as f64)])
                    .collect(),
            )],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default()
            .with_name("aromaticity-public")
            .with_prop("source", "preserved")
            .unwrap()
            .with_computed_prop("_CIPComputed", "true")
            .unwrap(),
    )
    .unwrap()
}

#[test]
fn canonical_public_signatures_and_defaults_compile() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::with_assigned_aromaticity;
    let _: for<'a, 'b> fn(&'a Molecule, &'b AromaticityParams) -> Result<Molecule, OperationError> =
        Molecule::with_assigned_aromaticity_with_params;
    let _: fn(&mut Molecule) -> Result<(), OperationError> = Molecule::assign_aromaticity_;
    let _: for<'a, 'b> fn(&'a mut Molecule, &'b AromaticityParams) -> Result<(), OperationError> =
        Molecule::assign_aromaticity_with_params_;
    assert_eq!(
        AromaticityParams::default(),
        AromaticityParams {
            model: AromaticityModel::Rdkit,
        }
    );
}

#[test]
fn binding_contract_exposes_exactly_the_frozen_seven_entries() {
    let expected = [
        "types.AromaticityModel",
        "types.AromaticityParams",
        "types.AromaticityError",
        "Molecule.with_assigned_aromaticity",
        "Molecule.with_assigned_aromaticity_with_params",
        "Molecule.assign_aromaticity_",
        "Molecule.assign_aromaticity_with_params_",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "aromaticity")
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    for row in &rows[..3] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    for row in &rows[3..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    assert_eq!(
        rows[3].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(
        rows[4].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(rows[5].callable.unwrap().state_model, StateModel::InPlace);
    assert_eq!(rows[6].callable.unwrap().state_model, StateModel::InPlace);
}

#[test]
fn generated_registry_and_all_four_matrices_share_one_exact_operation() {
    let feature = feature_spec("aromaticity").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_assigned_aromaticity_with_params").unwrap();
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
    assert_eq!(spec.derived_effects.recompute.bits(), 1 << 3);
    assert_eq!(
        spec.derived_effects.preserve.bits(),
        (1 << 0) | (1 << 1) | (1 << 5)
    );
    assert_eq!(
        spec.derived_effects.invalidate.bits(),
        (1 << 2) | (1 << 4) | (1 << 6) | (1 << 7)
    );
    assert_eq!(format!("{:?}", spec.cip_state), "ClearComputed");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "weak_aromaticity_assignment"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "assign_aromaticity_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| row.feature.name == "aromaticity")
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));
}

#[test]
fn value_operation_preserves_identity_coordinates_stereo_and_ordinary_props() {
    let source = aromaticity_molecule(6);
    let snapshot = source.clone();
    let output = source.with_assigned_aromaticity().unwrap();

    assert_eq!(source, snapshot);
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
    assert!(std::ptr::eq(source.coordinates(), output.coordinates()));
    assert_eq!(source.conformers(), output.conformers());
    assert_eq!(output.property("source"), Some("preserved"));
    assert_eq!(output.atoms()[0].prop("atom-note"), Some("atom-0"));
    assert_eq!(output.bonds()[0].prop("bond-note"), Some("bond-0"));
    assert_eq!(output.property("_CIPComputed"), None);
    assert_eq!(output.atoms()[0].prop("_CIPCode"), None);
    assert_eq!(output.bonds()[0].prop("_CIPCode"), None);
    assert!(output.atoms().iter().all(Atom::is_aromatic));
    assert!(output.bonds().iter().all(Bond::is_aromatic));
    assert!(
        output
            .bonds()
            .iter()
            .all(|bond| bond.order() == BondOrder::Aromatic)
    );
    assert!(format!("{output:?}").contains("derived_cache_is_empty: false"));
    assert!(source.atoms().iter().all(|atom| !atom.is_aromatic()));
    assert!(source.bonds().iter().all(|bond| !bond.is_aromatic()));
    assert_eq!(source.property("_CIPComputed"), Some("true"));
    assert!(!std::ptr::eq(source.topology(), output.topology()));
    assert!(!std::ptr::eq(source.properties(), output.properties()));
}

#[test]
fn explicit_model_and_both_inplace_forms_match_value_semantics() {
    let source = aromaticity_molecule(6);
    let default_output = source.with_assigned_aromaticity().unwrap();
    let explicit_default = source
        .with_assigned_aromaticity_with_params(&AromaticityParams::default())
        .unwrap();
    assert_eq!(default_output, explicit_default);

    let simple_params = AromaticityParams {
        model: AromaticityModel::Simple,
    };
    let simple = source
        .with_assigned_aromaticity_with_params(&simple_params)
        .unwrap();
    assert_eq!(simple, default_output);

    let mut short = source.clone();
    short.assign_aromaticity_().unwrap();
    assert_eq!(short, default_output);
    let mut explicit = source.clone();
    explicit
        .assign_aromaticity_with_params_(&simple_params)
        .unwrap();
    assert_eq!(explicit, simple);
    assert_eq!(source, aromaticity_molecule(6));
}

#[test]
fn custom_model_failure_is_structured_and_atomic_for_value_and_inplace() {
    let source = aromaticity_molecule(6);
    let snapshot = source.clone();
    let params = AromaticityParams {
        model: AromaticityModel::Custom,
    };
    let value_error = source
        .with_assigned_aromaticity_with_params(&params)
        .unwrap_err();
    assert!(matches!(
        &value_error,
        OperationError::Aromaticity(AromaticityError::UnsupportedModel {
            model: AromaticityModel::Custom,
            ..
        })
    ));
    assert!(value_error.source().is_some());
    assert_eq!(source, snapshot);

    let mut target = source.clone();
    let observer = target.clone();
    let inplace_error = target.assign_aromaticity_with_params_(&params).unwrap_err();
    assert_eq!(inplace_error, value_error);
    assert_eq!(target, source);
    assert!(std::ptr::eq(target.topology(), observer.topology()));
    assert!(std::ptr::eq(target.coordinates(), observer.coordinates()));
    assert!(std::ptr::eq(target.properties(), observer.properties()));
    assert_eq!(target.property("_CIPComputed"), Some("true"));
}

#[test]
fn single_output_local_edit_makes_multi_candidate_ordering_not_applicable() {
    let spec = operation_spec("with_assigned_aromaticity_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
    assert_eq!(spec.topology_edit, TopologyEditKind::Local);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.auto_remap, BlockSet::NONE);
}
