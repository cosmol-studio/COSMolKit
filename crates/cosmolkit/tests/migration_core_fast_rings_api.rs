#[path = "support/coordinate_views.rs"]
mod coordinate_views;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec, Conformer2D,
    CoordinateBlock, Element, Molecule, MoleculeOpKind, MoleculeOpOutput, MoleculeProperties,
    OperationDomain, OperationError, ParityPolicy, StateModel, StereoGroup, StereoGroupKind,
    SupportStatus, TopologyBlock, TopologyEditKind, feature_spec, operation_invariant,
    operation_parity, operation_spec, support_matrix,
};

fn ring_molecule() -> Molecule {
    let atoms = (0..5)
        .map(|index| {
            let spec = AtomSpec::new(Element::C)
                .with_prop("atom-row", index.to_string())
                .unwrap();
            let spec = if index == 0 {
                spec.with_computed_prop("_CIPCode", "R").unwrap()
            } else {
                spec
            };
            Atom::from_spec(AtomId::new(index), spec)
        })
        .collect();
    let endpoints = [(0, 1), (1, 2), (2, 3), (3, 0), (3, 4)];
    let bonds = endpoints
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            let spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single)
                .with_prop("bond-row", index.to_string())
                .unwrap();
            let spec = if index == 0 {
                spec.with_computed_prop("_CIPBondCode", "E").unwrap()
            } else {
                spec
            };
            Bond::from_spec(BondId::new(index), spec)
        })
        .collect();
    let topology = TopologyBlock::try_from_parts(
        atoms,
        bonds,
        Vec::new(),
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(0)],
            vec![BondId::new(0)],
        )],
    )
    .unwrap();
    Molecule::from_parts(
        topology,
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                17,
                vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [-1.0, 1.0]],
            )],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default()
            .with_name("fast-rings-public")
            .with_prop("source", "preserved")
            .unwrap()
            .with_computed_prop("_CIPComputed", "true")
            .unwrap(),
    )
    .unwrap()
}

#[test]
fn canonical_public_signatures_compile_without_an_invented_parameter_variant() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::with_assigned_rings;
    let _: fn(&mut Molecule) -> Result<(), OperationError> = Molecule::assign_rings_;
}

#[test]
fn binding_contract_exposes_exactly_the_frozen_two_callables() {
    let expected = ["Molecule.with_assigned_rings", "Molecule.assign_rings_"];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| expected.contains(&row.semantic_id))
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
    let feature = feature_spec("rings").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_assigned_rings").unwrap();
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::None);
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.access.read(), BlockSet::TOPOLOGY);
    assert_eq!(spec.access.write(), BlockSet::DERIVED_CACHE);
    assert_eq!(spec.may_mutate, BlockSet::DERIVED_CACHE);
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.derived_effects.recompute.bits(), 1 << 0);
    assert_eq!(
        spec.derived_effects.preserve.bits(),
        (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7)
    );
    assert_eq!(spec.derived_effects.invalidate.bits(), 1 << 1);
    assert_eq!(format!("{:?}", spec.cip_state), "Preserve");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "weak_ring_cache_assignment"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "fast_find_rings_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| row.feature.name == "rings")
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));
}

#[test]
fn value_operation_assigns_a_private_cache_and_preserves_every_input_block() {
    let source = ring_molecule();
    let output = source.with_assigned_rings().unwrap();

    assert!(std::ptr::eq(source.topology(), output.topology()));
    coordinate_views::assert_shared_coordinates(&source, &output);
    assert!(std::ptr::eq(source.properties(), output.properties()));
    assert_eq!(source.num_atoms(), output.num_atoms());
    assert_eq!(source.num_bonds(), output.num_bonds());
    assert_eq!(source.topology().adjacency, output.topology().adjacency);
    assert_eq!(
        source.topology().stereo_groups,
        output.topology().stereo_groups
    );
    assert_eq!(output.property("source"), Some("preserved"));
    assert_eq!(output.property("_CIPComputed"), Some("true"));
    assert_eq!(
        output.atom(AtomId::new(0)).unwrap().prop("_CIPCode"),
        Some("R")
    );
    assert_eq!(
        output.bond(BondId::new(0)).unwrap().prop("_CIPBondCode"),
        Some("E")
    );
    assert_eq!(
        output.atom(AtomId::new(4)).unwrap().prop("atom-row"),
        Some("4")
    );
    assert_eq!(
        output.bond(BondId::new(4)).unwrap().prop("bond-row"),
        Some("4")
    );
    assert!(format!("{source:?}").contains("derived_cache_is_empty: true"));
    assert!(format!("{output:?}").contains("derived_cache_is_empty: false"));
}

#[test]
fn repeated_value_and_inplace_forms_replace_the_assignment_idempotently() {
    let source = ring_molecule();
    let first = source.with_assigned_rings().unwrap();
    let second = first.with_assigned_rings().unwrap();
    assert_eq!(first, second);
    assert!(std::ptr::eq(first.topology(), second.topology()));
    coordinate_views::assert_shared_coordinates(&first, &second);
    assert!(std::ptr::eq(first.properties(), second.properties()));
    assert!(format!("{second:?}").contains("derived_cache_is_empty: false"));

    let mut target = source.clone();
    let observer = target.clone();
    target.assign_rings_().unwrap();
    assert_eq!(target, first);
    assert_eq!(observer, source);
    assert!(format!("{observer:?}").contains("derived_cache_is_empty: true"));
    assert!(format!("{target:?}").contains("derived_cache_is_empty: false"));
}

#[test]
fn invalid_topology_is_rejected_before_a_live_operation_can_observe_partial_state() {
    let source = ring_molecule();
    let observer = source.clone();
    let atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
    let self_loop = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(0), BondOrder::Single),
    );
    let error = TopologyBlock::try_from_parts(vec![atom], vec![self_loop], Vec::new(), Vec::new())
        .unwrap_err();
    assert!(format!("{error:?}").contains("SelfLoopBond"));
    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(std::ptr::eq(source.properties(), observer.properties()));
    assert!(format!("{source:?}").contains("derived_cache_is_empty: true"));
}

#[test]
fn single_output_no_mapping_contract_makes_candidate_ordering_not_applicable() {
    let spec = operation_spec("with_assigned_rings").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.auto_remap, BlockSet::NONE);
}
