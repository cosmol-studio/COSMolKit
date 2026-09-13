use std::error::Error as _;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec, Conformer2D,
    CoordinateBlock, Element, Molecule, MoleculeOpKind, MoleculeOpOutput, MoleculeProperties,
    OperationDomain, OperationError, ParityPolicy, StateModel, StereoGroup, StereoGroupKind,
    SupportStatus, TopologyBlock, TopologyEditKind, feature_spec, operation_invariant,
    operation_parity, operation_spec, support_matrix,
};

fn radical_molecule() -> Molecule {
    let carbon = Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C)
            .with_no_implicit(true)
            .with_radical_electrons(7)
            .with_prop("atom-label", "carbon")
            .unwrap()
            .with_computed_prop("_CIPCode", "R")
            .unwrap(),
    );
    let oxygen = Atom::from_spec(
        AtomId::new(1),
        AtomSpec::new(Element::O).with_radical_electrons(2),
    );
    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
            .with_prop("bond-label", "single")
            .unwrap()
            .with_computed_prop("_CIPBondCode", "E")
            .unwrap(),
    );
    let topology = TopologyBlock::try_from_parts(
        vec![carbon, oxygen],
        vec![bond],
        Vec::new(),
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(0)],
            vec![BondId::new(0)],
        )],
    )
    .unwrap();
    let properties = MoleculeProperties::default()
        .with_name("radicals-public")
        .with_prop("source", "preserved")
        .unwrap()
        .with_computed_prop("_CIPComputed", "true")
        .unwrap();
    Molecule::from_parts(
        topology,
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(5, vec![[0.0, 0.0], [1.0, 0.0]])],
            ..CoordinateBlock::default()
        },
        properties,
    )
    .unwrap()
}

fn unsupported_bond_molecule() -> Molecule {
    let atoms = vec![
        Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_radical_electrons(6),
        ),
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
    ];
    let bonds = vec![Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::ThreeCenter),
    )];
    Molecule::from_parts(
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap(),
        CoordinateBlock::default(),
        MoleculeProperties::default()
            .with_prop("source", "failure-preserved")
            .unwrap(),
    )
    .unwrap()
}

#[test]
fn canonical_public_signatures_compile_without_an_invented_parameter_variant() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::with_assigned_radicals;
    let _: fn(&mut Molecule) -> Result<(), OperationError> = Molecule::assign_radicals_;
}

#[test]
fn binding_contract_exposes_exactly_the_frozen_two_callables() {
    let expected = [
        "Molecule.with_assigned_radicals",
        "Molecule.assign_radicals_",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "radicals")
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    for row in &rows {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
        assert_eq!(row.callable.unwrap().parameters.len(), 0);
    }
    assert_eq!(
        rows[0].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(rows[1].callable.unwrap().state_model, StateModel::InPlace);
}

#[test]
fn generated_registry_and_all_four_matrices_share_the_exact_operation() {
    let feature = feature_spec("radicals").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_assigned_radicals").unwrap();
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::None);
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
    assert_eq!(format!("{:?}", spec.cip_state), "ClearComputed");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "weak_radical_assignment"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "assign_radicals_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| row.feature.name == "radicals")
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));
}

#[test]
fn value_operation_assigns_rows_and_preserves_mapping_coordinates_stereo_and_props() {
    let source = radical_molecule();
    let output = source.with_assigned_radicals().unwrap();

    assert_eq!(source.num_atoms(), output.num_atoms());
    assert_eq!(source.num_bonds(), output.num_bonds());
    assert_eq!(output.atom(AtomId::new(0)).unwrap().radical_electrons(), 3);
    assert_eq!(output.atom(AtomId::new(1)).unwrap().radical_electrons(), 2);
    assert_eq!(output.atom(AtomId::new(0)).unwrap().id(), AtomId::new(0));
    assert_eq!(output.bond(BondId::new(0)).unwrap().id(), BondId::new(0));
    assert_eq!(source.topology().adjacency, output.topology().adjacency);
    assert_eq!(
        source.topology().stereo_groups,
        output.topology().stereo_groups
    );
    assert!(std::ptr::eq(source.coordinates(), output.coordinates()));
    assert_eq!(source.conformers(), output.conformers());
    assert_eq!(output.property("source"), Some("preserved"));
    assert_eq!(
        output.atom(AtomId::new(0)).unwrap().prop("atom-label"),
        Some("carbon")
    );
    assert_eq!(
        output.bond(BondId::new(0)).unwrap().prop("bond-label"),
        Some("single")
    );
    assert_eq!(output.property("_CIPComputed"), None);
    assert_eq!(output.atom(AtomId::new(0)).unwrap().prop("_CIPCode"), None);
    assert_eq!(
        output.bond(BondId::new(0)).unwrap().prop("_CIPBondCode"),
        None
    );
    assert!(format!("{output:?}").contains("derived_cache_is_empty: true"));

    assert_eq!(source.atom(AtomId::new(0)).unwrap().radical_electrons(), 7);
    assert_eq!(source.property("_CIPComputed"), Some("true"));
    assert_eq!(
        source.atom(AtomId::new(0)).unwrap().prop("_CIPCode"),
        Some("R")
    );
    assert_eq!(
        source.bond(BondId::new(0)).unwrap().prop("_CIPBondCode"),
        Some("E")
    );
    assert!(!std::ptr::eq(source.topology(), output.topology()));
    assert!(!std::ptr::eq(source.properties(), output.properties()));
}

#[test]
fn inplace_and_value_forms_commit_the_same_authoritative_result() {
    let source = radical_molecule();
    let expected = source.with_assigned_radicals().unwrap();
    let mut target = source.clone();
    target.assign_radicals_().unwrap();
    assert_eq!(target, expected);
    assert_eq!(source.atom(AtomId::new(0)).unwrap().radical_electrons(), 7);
}

#[test]
fn typed_algorithm_failure_is_atomic_for_value_and_inplace_wrappers() {
    let source = unsupported_bond_molecule();
    let value_error = source.with_assigned_radicals().unwrap_err();
    assert!(matches!(&value_error, OperationError::Radical(_)));
    assert_eq!(
        value_error.to_string(),
        "radical assignment failed: Bad bond type"
    );
    assert_eq!(
        format!("{value_error:?}"),
        "Radical(Valence(BadBondType { bond: Some(BondId(0)), order: ThreeCenter }))"
    );
    assert!(value_error.source().is_some());

    let mut target = source.clone();
    let observer = target.clone();
    let error = target.assign_radicals_().unwrap_err();
    assert_eq!(error, value_error);
    assert_eq!(target, source);
    assert_eq!(observer, source);
    assert!(std::ptr::eq(target.topology(), observer.topology()));
    assert!(std::ptr::eq(target.coordinates(), observer.coordinates()));
    assert!(std::ptr::eq(target.properties(), observer.properties()));
}

#[test]
fn single_output_no_mapping_contract_makes_candidate_ordering_not_applicable() {
    let spec = operation_spec("with_assigned_radicals").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.auto_remap, BlockSet::NONE);
}
