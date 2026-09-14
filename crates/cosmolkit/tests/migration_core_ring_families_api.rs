#[path = "support/coordinate_views.rs"]
mod coordinate_views;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec, Conformer2D,
    CoordinateBlock, Element, Molecule, MoleculeOpKind, MoleculeOpOutput, MoleculeProperties,
    OperationDomain, OperationError, ParityPolicy, RingSearchParams, StateModel, StereoGroup,
    StereoGroupKind, SupportStatus, TopologyBlock, TopologyEditKind, feature_spec,
    operation_invariant, operation_parity, operation_spec, support_matrix,
};

fn molecule_with_bond_orders(orders: [BondOrder; 4]) -> Molecule {
    let atoms = (0..4)
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
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 0)]
        .into_iter()
        .zip(orders)
        .enumerate()
        .map(|(index, ((begin, end), order))| {
            let spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
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
                23,
                vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            )],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default()
            .with_name("ring-families-public")
            .with_prop("source", "preserved")
            .unwrap()
            .with_computed_prop("_CIPComputed", "true")
            .unwrap(),
    )
    .unwrap()
}

fn ordinary_ring_molecule() -> Molecule {
    molecule_with_bond_orders([BondOrder::Single; 4])
}

#[test]
fn canonical_public_signatures_and_source_defaults_are_exact() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> =
        Molecule::with_assigned_ring_families;
    let _: for<'a, 'b> fn(&'a Molecule, &'b RingSearchParams) -> Result<Molecule, OperationError> =
        Molecule::with_assigned_ring_families_with_params;
    let _: fn(&mut Molecule) -> Result<(), OperationError> = Molecule::assign_ring_families_;
    let _: for<'a, 'b> fn(&'a mut Molecule, &'b RingSearchParams) -> Result<(), OperationError> =
        Molecule::assign_ring_families_with_params_;
    assert_eq!(RingSearchParams::default().include_dative_bonds, false);
    assert_eq!(RingSearchParams::default().include_hydrogen_bonds, false);
}

#[test]
fn binding_contract_has_the_parameter_type_and_exact_four_family_callables() {
    let expected = [
        "types.RingSearchParams",
        "Molecule.with_assigned_ring_families",
        "Molecule.with_assigned_ring_families_with_params",
        "Molecule.assign_ring_families_",
        "Molecule.assign_ring_families_with_params_",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| expected.contains(&row.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(rows[0].item, BindingItem::Type);
    assert_eq!(rows[0].owner, BindingOwner::Type);
    assert_eq!(rows[0].exposure, BindingExposure::Public);
    assert_eq!(rows[0].support, BindingSupport::Supported);
    assert_eq!(rows[0].parity, BindingParity::NotApplicable);
    for row in &rows[1..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.feature, "rings");
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    assert_eq!(rows[1].callable.unwrap().parameters.len(), 0);
    assert_eq!(rows[2].callable.unwrap().parameters.len(), 1);
    assert_eq!(rows[3].callable.unwrap().parameters.len(), 0);
    assert_eq!(rows[4].callable.unwrap().parameters.len(), 1);
    assert_eq!(
        rows[1].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(
        rows[2].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(rows[3].callable.unwrap().state_model, StateModel::InPlace);
    assert_eq!(rows[4].callable.unwrap().state_model, StateModel::InPlace);
}

#[test]
fn generated_registry_and_all_four_matrices_share_the_family_operation() {
    let feature = feature_spec("rings").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_assigned_ring_families_with_params").unwrap();
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::None);
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.access.read(), BlockSet::TOPOLOGY);
    assert_eq!(spec.access.write(), BlockSet::DERIVED_CACHE);
    assert_eq!(spec.may_mutate, BlockSet::DERIVED_CACHE);
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.derived_effects.recompute.bits(), 1 << 1);
    assert_eq!(
        spec.derived_effects.preserve.bits(),
        (1 << 0) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7)
    );
    assert_eq!(spec.derived_effects.invalidate.bits(), 0);
    assert_eq!(format!("{:?}", spec.cip_state), "Preserve");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "weak_ring_family_cache_assignment"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "find_ring_families_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| {
            row.operation
                .is_some_and(|operation| std::ptr::eq(operation, spec))
        })
        .unwrap();
    assert!(std::ptr::eq(support.feature, feature));
}

#[test]
fn value_operation_installs_family_state_and_preserves_all_input_blocks() {
    let source = ordinary_ring_molecule();
    let output = source.with_assigned_ring_families().unwrap();

    assert!(std::ptr::eq(source.topology(), output.topology()));
    coordinate_views::assert_shared_coordinates(&source, &output);
    assert!(std::ptr::eq(source.properties(), output.properties()));
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
    assert!(format!("{source:?}").contains("derived_cache_is_empty: true"));
    assert!(format!("{output:?}").contains("derived_cache_is_empty: false"));
}

#[test]
fn all_four_parameter_combinations_cross_the_live_boundary() {
    let source = molecule_with_bond_orders([
        BondOrder::Single,
        BondOrder::Dative,
        BondOrder::Hydrogen,
        BondOrder::Single,
    ]);
    for include_dative_bonds in [false, true] {
        for include_hydrogen_bonds in [false, true] {
            let params = RingSearchParams {
                include_dative_bonds,
                include_hydrogen_bonds,
            };
            let output = source
                .with_assigned_ring_families_with_params(&params)
                .unwrap();
            assert!(std::ptr::eq(source.topology(), output.topology()));
            coordinate_views::assert_shared_coordinates(&source, &output);
            assert!(std::ptr::eq(source.properties(), output.properties()));
            assert!(format!("{output:?}").contains("derived_cache_is_empty: false"));
        }
    }
    assert!(format!("{source:?}").contains("derived_cache_is_empty: true"));
}

#[test]
fn family_assignment_preserves_an_existing_ordinary_ring_assignment() {
    let source = ordinary_ring_molecule();
    let rings = source.with_assigned_rings().unwrap();
    let families = rings.with_assigned_ring_families().unwrap();
    assert!(std::ptr::eq(rings.topology(), families.topology()));
    coordinate_views::assert_shared_coordinates(&rings, &families);
    assert!(std::ptr::eq(rings.properties(), families.properties()));
    assert!(format!("{rings:?}").contains("derived_cache_is_empty: false"));
    assert!(format!("{families:?}").contains("derived_cache_is_empty: false"));
    assert!(families.with_assigned_ring_families().is_ok());
}

#[test]
fn inplace_success_matches_value_semantics_and_keeps_observers_unchanged() {
    let source = ordinary_ring_molecule();
    let expected = source.with_assigned_ring_families().unwrap();
    let mut target = source.clone();
    let observer = target.clone();
    target.assign_ring_families_().unwrap();
    assert_eq!(target, expected);
    assert_eq!(observer, source);
    assert!(std::ptr::eq(observer.topology(), source.topology()));
    coordinate_views::assert_shared_coordinates(&observer, &source);
    assert!(std::ptr::eq(observer.properties(), source.properties()));
    assert!(format!("{observer:?}").contains("derived_cache_is_empty: true"));
    assert!(format!("{target:?}").contains("derived_cache_is_empty: false"));
}

#[test]
fn detached_failure_is_structured_and_inplace_failure_is_atomic() {
    let source = Molecule::new();
    let value_error = source.with_assigned_ring_families().unwrap_err();
    assert!(matches!(value_error, OperationError::Rings(_)));
    assert!(format!("{source:?}").contains("derived_cache_is_empty: true"));

    let mut target = source.clone();
    let observer = target.clone();
    let inplace_error = target.assign_ring_families_().unwrap_err();
    assert!(matches!(inplace_error, OperationError::Rings(_)));
    assert_eq!(target, observer);
    assert!(std::ptr::eq(target.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&target, &observer);
    assert!(std::ptr::eq(target.properties(), observer.properties()));
    assert!(format!("{target:?}").contains("derived_cache_is_empty: true"));
}

#[test]
fn single_output_no_mapping_contract_excludes_multi_output_ordering() {
    let spec = operation_spec("with_assigned_ring_families_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.auto_remap, BlockSet::NONE);
}
