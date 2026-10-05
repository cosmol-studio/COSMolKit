#[path = "support/coordinate_views.rs"]
mod coordinate_views;
use cosmolkit::*;
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

fn rows() -> Vec<Vec<f64>> {
    vec![vec![9., 8., 7.], vec![6., 5., 4.]]
}
fn unchanged_graph(source: &Molecule, result: &Molecule) {
    assert!(std::ptr::eq(source.topology(), result.topology()));
    assert!(std::ptr::eq(source.properties(), result.properties()));
    assert_eq!(result.property("_CIPComputed"), Some("true"));
    assert_eq!(
        result.atom(AtomId::new(0)).unwrap().prop("_CIPCode"),
        Some(&cosmolkit_model::PropertyValue::from("R"))
    );
    assert_eq!(
        result.bond(BondId::new(0)).unwrap().prop("_CIPBondCode"),
        Some(&cosmolkit_model::PropertyValue::from("E"))
    );
}
#[test]
fn two_dimensional_install_is_value_cow_with_original_mixed_provenance() {
    let source = transform_molecule();
    let before = source.to_builder();
    let result = source.with_2d_coordinate_block(rows()).unwrap();
    assert_eq!(source.to_builder().coordinates(), before.coordinates());
    assert_eq!(result.conformers_3d(), source.conformers_3d());
    assert_eq!(result.coordinates_2d().unwrap(), &[[9., 8.], [6., 5.]]);
    assert_eq!(
        result.to_builder().coordinates().source_coordinate_dim,
        Some(CoordinateDimension::TwoD)
    );
    unchanged_graph(&source, &result);
    coordinate_views::assert_detached_coordinates(&source, &result);
}
#[test]
fn replacement_by_actual_id_preserves_is3d_and_cip() {
    let source = transform_molecule();
    let result = source
        .with_3d_coordinates_with_params(rows(), &Replace3DCoordinatesParams { conformer_id: 7 })
        .unwrap();
    assert_eq!(result.conformers_3d()[1].id(), 7);
    assert!(result.conformers_3d()[1].is_3d());
    assert_eq!(
        result.conformers_3d()[1],
        Conformer3D::new(7, vec![[9., 8., 7.], [6., 5., 4.]], true)
    );
    assert_eq!(result.conformers_3d()[0], source.conformers_3d()[0]);
    assert_eq!(result.conformers_3d()[2], source.conformers_3d()[2]);
    assert_eq!(source.conformers_3d()[1].coordinates()[0], [0., 1., 0.]);
    unchanged_graph(&source, &result);
}
#[test]
fn append_inplace_commits_scalar_position_and_keeps_shared_observer() {
    let mut source = transform_molecule();
    let observer = source.clone();
    let value = source
        .with_added_3d_conformer_with_params(rows(), &Coordinate3DInputParams { is_3d: false })
        .unwrap();
    assert_eq!(
        source
            .add_3d_conformer_with_params_(rows(), &Coordinate3DInputParams { is_3d: false })
            .unwrap(),
        3
    );
    assert_eq!(source.conformers_3d()[3].id(), 12);
    assert!(!source.conformers_3d()[3].is_3d());
    assert_eq!(observer.conformers_3d().len(), 3);
    assert_eq!(
        source.to_builder().coordinates(),
        value.to_builder().coordinates()
    );
    unchanged_graph(&observer, &source);
}
#[test]
fn only_false_and_clear_preserve_xy_and_original_source_dimension() {
    let mut source = transform_molecule();
    let observer = source.clone();
    let params = Coordinate3DInputParams { is_3d: false };
    let value = source
        .with_only_3d_conformer_with_params(rows(), &params)
        .unwrap();
    assert_eq!(
        source
            .set_only_3d_conformer_with_params_(rows(), &params)
            .unwrap(),
        0
    );
    assert_eq!(
        source.to_builder().coordinates(),
        value.to_builder().coordinates()
    );
    assert_eq!(
        source.to_builder().coordinates().source_coordinate_dim,
        Some(CoordinateDimension::TwoD)
    );
    assert_eq!(source.coordinates_2d(), observer.coordinates_2d());
    let value = source.with_cleared_3d_conformers().unwrap();
    source.clear_3d_conformers_().unwrap();
    assert_eq!(
        source.to_builder().coordinates(),
        value.to_builder().coordinates()
    );
    assert!(source.conformers_3d().is_empty());
    assert_eq!(source.coordinates_2d(), observer.coordinates_2d());
    unchanged_graph(&observer, &source);
}
#[test]
fn constructor_retains_explicit_provenance_and_none_keeps_existing_fallback() {
    let source = transform_molecule();
    let builder = source.to_builder();
    for provenance in [
        Some(CoordinateDimension::TwoD),
        Some(CoordinateDimension::ThreeD),
        None,
    ] {
        let mut block = builder.coordinates().clone();
        block.source_coordinate_dim = provenance;
        let rebuilt = Molecule::from_parts(
            source.topology().clone(),
            block,
            source.properties().clone(),
        )
        .unwrap();
        assert_eq!(
            rebuilt.to_builder().coordinates().source_coordinate_dim,
            provenance.or(Some(CoordinateDimension::ThreeD))
        );
    }
}
#[test]
fn public_coordinate_errors_are_atomic_for_all_live_input_transitions() {
    let mut source = transform_molecule();
    let observer = source.clone();
    assert!(matches!(
        source.set_2d_coordinates_(vec![vec![0.; 2]]),
        Err(OperationError::CoordinateInput(
            CoordinateInputError::RowCount { .. }
        ))
    ));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(matches!(
        source.set_3d_coordinates_with_params_(
            rows(),
            &Replace3DCoordinatesParams { conformer_id: 99 }
        ),
        Err(OperationError::CoordinateInput(
            CoordinateInputError::ConformerNotFound {
                conformer_id: 99,
                count: 3
            }
        ))
    ));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(
        source
            .add_3d_conformer_(vec![vec![f64::NAN; 3]; 2])
            .is_err()
    );
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(source.set_only_3d_conformer_(vec![vec![0.; 2]; 2]).is_err());
    coordinate_views::assert_shared_coordinates(&source, &observer);
    unchanged_graph(&observer, &source);
}
#[test]
fn generated_contracts_keep_exact_coordinate_only_authority() {
    for method in [
        "with_2d_coordinate_block_with_params",
        "with_3d_coordinates_with_params",
        "with_added_3d_conformer_with_params",
        "with_only_3d_conformer_with_params",
        "with_cleared_3d_conformers",
    ] {
        let spec = operation_spec(method).unwrap();
        assert_eq!(spec.access.read(), BlockSet::TOPOLOGY);
        assert_eq!(
            spec.access.write(),
            BlockSet::COORDINATES.union(BlockSet::DERIVED_CACHE)
        );
        assert_eq!(spec.auto_remap, BlockSet::NONE);
        assert_eq!(format!("{:?}", spec.cip_state), "Preserve");
        assert_eq!(spec.status, FunctionStatus::Native);
        assert_eq!(spec.derived_effects.invalidate.bits(), 1 << 6);
    }
}

#[test]
fn xyz_query_resolves_noncontiguous_ids_preserves_borrowing_and_never_generates() {
    let source = transform_molecule();
    let observer = source.clone();
    assert!(std::ptr::eq(
        source.coordinates_3d(7).unwrap(),
        source.conformers_3d()[1].coordinates()
    ));
    assert_eq!(
        source.coordinates_3d(3).unwrap(),
        source.conformers_3d()[0].coordinates()
    );
    assert!(matches!(
        source.coordinates_3d(1),
        Err(Coordinate3DReadError::ConformerNotFound {
            conformer_id: 1,
            count: 3
        })
    ));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    unchanged_graph(&source, &observer);
    let empty = Molecule::new();
    assert!(matches!(
        empty.coordinates_3d(0),
        Err(Coordinate3DReadError::ConformerNotFound {
            conformer_id: 0,
            count: 0
        })
    ));
    assert!(!empty.has_2d_coordinates());
    assert!(empty.conformers_3d().is_empty());
}

#[test]
fn missing_replacement_ids_and_default_preserve_all_live_blocks_and_sharing() {
    let mut source = transform_molecule();
    let observer = source.clone();
    for conformer_id in [1, 0] {
        let params = Replace3DCoordinatesParams { conformer_id };
        for result in [
            source
                .with_3d_coordinates_with_params(rows(), &params)
                .map(|_| ()),
            source.set_3d_coordinates_with_params_(rows(), &params),
        ] {
            assert!(matches!(result, Err(OperationError::CoordinateInput(
                CoordinateInputError::ConformerNotFound { conformer_id: id, count: 3 }
            )) if id == conformer_id));
            coordinate_views::assert_shared_coordinates(&source, &observer);
            unchanged_graph(&observer, &source);
        }
    }
    assert!(matches!(
        source.set_3d_coordinates_(rows()),
        Err(OperationError::CoordinateInput(
            CoordinateInputError::ConformerNotFound {
                conformer_id: 0,
                count: 3
            }
        ))
    ));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    source
        .set_3d_coordinates_with_params_(rows(), &Replace3DCoordinatesParams { conformer_id: 7 })
        .unwrap();
    assert_eq!(
        source.coordinates_3d(7).unwrap(),
        &[[9., 8., 7.], [6., 5., 4.]]
    );
    assert_eq!(source.conformers_3d()[0], observer.conformers_3d()[0]);
    assert_eq!(source.conformers_3d()[2], observer.conformers_3d()[2]);
    assert_eq!(source.coordinates_2d(), observer.coordinates_2d());
    assert_eq!(
        source.to_builder().coordinates().source_coordinate_dim,
        observer.to_builder().coordinates().source_coordinate_dim
    );
    unchanged_graph(&observer, &source);
}
