#[path = "support/coordinate_views.rs"]
mod coordinate_views;

use std::error::Error as _;

use cosmolkit::{
    Atom, AtomId, AtomPositionParams, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem,
    BindingOwner, BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec,
    Conformer2D, Conformer3D, CoordinateBlock, Element, Molecule, MoleculeOpKind, MoleculeOpOutput,
    MoleculeProperties, OperationDomain, OperationError, ParityPolicy, StateModel, StereoGroup,
    StereoGroupKind, SupportStatus, TopologyBlock, TopologyEditKind, TransformError, feature_spec,
    operation_invariant, operation_parity, operation_spec, support_matrix,
};

fn transform_molecule() -> Molecule {
    let atoms = vec![
        Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("atom-label", "first")
                .unwrap()
                .with_computed_prop("_CIPCode", "R")
                .unwrap(),
        ),
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
    ];
    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
            .with_prop("bond-label", "single")
            .unwrap()
            .with_computed_prop("_CIPBondCode", "E")
            .unwrap(),
    );
    let topology = TopologyBlock::try_from_parts(
        atoms,
        vec![bond],
        Vec::new(),
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(0)],
            vec![BondId::new(0)],
        )],
    )
    .unwrap();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(2, vec![[0.0, 0.0], [1.0, 0.0]])],
        conformers_3d: vec![
            Conformer3D::new(3, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], false)
                .with_prop("kind", "not-3d"),
            Conformer3D::new(7, vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]], true)
                .with_prop("kind", "default"),
            Conformer3D::new(11, vec![[0.0, 2.0, 0.0], [1.0, 2.0, 0.0]], true)
                .with_prop("kind", "explicit"),
        ],
        source_coordinate_dim: None,
    };
    let properties = MoleculeProperties::default()
        .with_name("transforms-public")
        .with_prop("source", "preserved")
        .unwrap()
        .with_computed_prop("_CIPComputed", "true")
        .unwrap();
    Molecule::from_parts(topology, coordinates, properties).unwrap()
}

#[test]
fn canonical_public_signatures_and_default_are_exact() {
    let _: fn(&Molecule, AtomId, [f64; 3]) -> Result<Molecule, OperationError> =
        Molecule::with_atom_position;
    let _: for<'a, 'b> fn(
        &'a Molecule,
        AtomId,
        [f64; 3],
        &'b AtomPositionParams,
    ) -> Result<Molecule, OperationError> = Molecule::with_atom_position_with_params;
    let _: fn(&mut Molecule, AtomId, [f64; 3]) -> Result<(), OperationError> =
        Molecule::set_atom_position_;
    let _: for<'a, 'b> fn(
        &'a mut Molecule,
        AtomId,
        [f64; 3],
        &'b AtomPositionParams,
    ) -> Result<(), OperationError> = Molecule::set_atom_position_with_params_;
    assert_eq!(AtomPositionParams::default().conformer_id, None);
}

#[test]
fn binding_contract_exposes_exactly_two_types_and_four_callables() {
    let expected = [
        "types.AtomPositionParams",
        "types.TransformError",
        "Molecule.with_atom_position",
        "Molecule.with_atom_position_with_params",
        "Molecule.set_atom_position_",
        "Molecule.set_atom_position_with_params_",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "transforms")
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    for row in &rows {
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    for row in &rows[..2] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
    }
    for row in &rows[2..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
    }
    assert_eq!(rows[2].callable.unwrap().parameters.len(), 2);
    assert_eq!(rows[3].callable.unwrap().parameters.len(), 3);
    assert_eq!(rows[4].callable.unwrap().parameters.len(), 2);
    assert_eq!(rows[5].callable.unwrap().parameters.len(), 3);
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
    assert_eq!(rows[4].javascript_name, "setAtomPosition");
    assert_eq!(rows[5].javascript_name, "setAtomPositionWithParams");
}

#[test]
fn generated_registry_and_all_four_matrices_share_the_coordinate_operation() {
    let feature = feature_spec("transforms").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_atom_position_with_params").unwrap();
    assert_eq!(spec.domain, OperationDomain::Coordinate);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::None);
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.access.read(), BlockSet::NONE);
    assert_eq!(
        spec.access.write(),
        BlockSet::TOPOLOGY
            .union(BlockSet::COORDINATES)
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(spec.may_mutate, spec.access.write());
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.derived_effects.recompute.bits(), 0);
    assert_eq!(
        spec.derived_effects.preserve.bits(),
        (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 7)
    );
    assert_eq!(spec.derived_effects.invalidate.bits(), (1 << 4) | (1 << 6));
    assert_eq!(format!("{:?}", spec.cip_state), "ClearComputed");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert!(spec.io_roundtrip);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "coordinate_atom_position"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "set_atom_position_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| row.feature.name == "transforms")
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));
}

#[test]
fn default_value_operation_updates_first_true_3d_conformer_and_clears_cip() {
    let source = transform_molecule();
    let output = source
        .with_atom_position(AtomId::new(1), [9.0, 8.0, 7.0])
        .unwrap();

    assert_eq!(source.conformers_3d()[0], output.conformers_3d()[0]);
    assert_eq!(output.conformers_3d()[1].coordinates()[1], [9.0, 8.0, 7.0]);
    assert_eq!(source.conformers_3d()[2], output.conformers_3d()[2]);
    assert_eq!(
        source.to_builder().coordinates().conformers_2d,
        output.to_builder().coordinates().conformers_2d
    );
    assert_eq!(output.conformers_3d()[1].id(), 7);
    assert_eq!(
        output.conformers_3d()[1]
            .props()
            .get("kind")
            .map(String::as_str),
        Some("default")
    );
    assert_eq!(source.topology().adjacency, output.topology().adjacency);
    assert_eq!(
        source.topology().stereo_groups,
        output.topology().stereo_groups
    );
    assert_eq!(output.property("source"), Some("preserved"));
    assert_eq!(
        output.atom(AtomId::new(0)).unwrap().prop("atom-label"),
        Some("first")
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

    assert_eq!(source.conformers_3d()[1].coordinates()[1], [1.0, 1.0, 0.0]);
    assert_eq!(source.property("_CIPComputed"), Some("true"));
    assert_eq!(
        source.atom(AtomId::new(0)).unwrap().prop("_CIPCode"),
        Some("R")
    );
    assert!(!std::ptr::eq(source.topology(), output.topology()));
    coordinate_views::assert_detached_coordinates(&source, &output);
    assert!(!std::ptr::eq(source.properties(), output.properties()));
}

#[test]
fn explicit_conformer_id_updates_only_the_selected_row_set() {
    let source = transform_molecule();
    let params = AtomPositionParams {
        conformer_id: Some(11),
    };
    let output = source
        .with_atom_position_with_params(AtomId::new(0), [-1.0, -2.0, -3.0], &params)
        .unwrap();
    assert_eq!(source.conformers_3d()[0], output.conformers_3d()[0]);
    assert_eq!(source.conformers_3d()[1], output.conformers_3d()[1]);
    assert_eq!(
        output.conformers_3d()[2].coordinates()[0],
        [-1.0, -2.0, -3.0]
    );
    assert_eq!(output.conformers_3d()[2].id(), 11);
    assert!(output.conformers_3d()[2].is_3d());
    assert_eq!(
        output.conformers_3d()[2]
            .props()
            .get("kind")
            .map(String::as_str),
        Some("explicit")
    );
}

#[test]
fn inplace_and_value_forms_commit_the_same_result_and_keep_observers_unchanged() {
    let source = transform_molecule();
    let params = AtomPositionParams {
        conformer_id: Some(11),
    };
    let expected = source
        .with_atom_position_with_params(AtomId::new(0), [4.0, 5.0, 6.0], &params)
        .unwrap();
    let mut target = source.clone();
    let observer = target.clone();
    target
        .set_atom_position_with_params_(AtomId::new(0), [4.0, 5.0, 6.0], &params)
        .unwrap();
    assert_eq!(target, expected);
    assert_eq!(observer, source);
    assert!(std::ptr::eq(observer.topology(), source.topology()));
    coordinate_views::assert_shared_coordinates(&observer, &source);
    assert!(std::ptr::eq(observer.properties(), source.properties()));
}

#[test]
fn structured_algorithm_failures_are_atomic_for_value_and_inplace_forms() {
    let source = transform_molecule();
    let value_error = source
        .with_atom_position(AtomId::new(2), [1.0, 2.0, 3.0])
        .unwrap_err();
    assert_eq!(
        value_error,
        OperationError::Transform(TransformError::AtomOutOfRange {
            role: "atom",
            atom: AtomId::new(2),
            atom_count: 2,
        })
    );
    assert!(value_error.source().is_some());
    assert_eq!(source.conformers_3d()[1].coordinates()[0], [0.0, 1.0, 0.0]);

    let mut target = source.clone();
    let observer = target.clone();
    let inplace_error = target
        .set_atom_position_(AtomId::new(0), [f64::NAN, 2.0, 3.0])
        .unwrap_err();
    assert_eq!(
        inplace_error,
        OperationError::Transform(TransformError::NonFinitePoint {
            role: "position",
            axis: "x",
        })
    );
    assert_eq!(target, source);
    assert_eq!(observer, source);
    assert!(std::ptr::eq(target.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&target, &observer);
    assert!(std::ptr::eq(target.properties(), observer.properties()));
}

#[test]
fn missing_default_and_explicit_conformers_remain_distinct_errors() {
    let source = transform_molecule();
    let no_3d = Molecule::from_parts(
        source.topology().clone(),
        CoordinateBlock {
            conformers_2d: source.to_builder().coordinates().conformers_2d.to_vec(),
            conformers_3d: vec![Conformer3D::new(
                3,
                vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                false,
            )],
            source_coordinate_dim: None,
        },
        source.properties().clone(),
    )
    .unwrap();
    assert_eq!(
        no_3d
            .with_atom_position(AtomId::new(0), [1.0, 2.0, 3.0])
            .unwrap_err(),
        OperationError::Transform(TransformError::No3dConformer)
    );
    assert_eq!(
        source
            .with_atom_position_with_params(
                AtomId::new(0),
                [1.0, 2.0, 3.0],
                &AtomPositionParams {
                    conformer_id: Some(99),
                },
            )
            .unwrap_err(),
        OperationError::Transform(TransformError::ConformerNotFound { conformer_id: 99 })
    );
}

#[test]
fn single_output_no_mapping_contract_excludes_multi_output_ordering() {
    let spec = operation_spec("with_atom_position_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.auto_remap, BlockSet::NONE);
}
