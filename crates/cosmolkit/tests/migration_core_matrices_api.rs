#[path = "support/coordinate_views.rs"]
mod coordinate_views;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, Bond, BondId, BondOrder, BondSpec, Conformer3D, CoordinateBlock,
    DenseMatrix, DistanceMatrix3dParams, DistanceMatrixParams, Element, MatrixError, Molecule,
    MoleculeProperties, StateModel, StereoGroup, StereoGroupKind, TopologyBlock,
    operation_invariant, operation_parity, operation_spec, support_matrix,
};

fn molecule_with_coordinates(coordinates: CoordinateBlock) -> Molecule {
    let atoms = [Element::C, Element::O, Element::N]
        .into_iter()
        .enumerate()
        .map(|(index, element)| {
            let spec = AtomSpec::new(element)
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
    let bonds = [
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
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
        coordinates,
        MoleculeProperties::default()
            .with_name("matrix-public")
            .with_prop("source", "preserved")
            .unwrap(),
    )
    .unwrap()
}

fn molecule_3d() -> Molecule {
    molecule_with_coordinates(CoordinateBlock {
        conformers_3d: vec![
            Conformer3D::new(
                8,
                vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 4.0, 0.0]],
                true,
            ),
            Conformer3D::new(
                3,
                vec![[0.0, 0.0, 0.0], [0.0, 0.0, 12.0], [0.0, 5.0, 12.0]],
                true,
            ),
        ],
        ..CoordinateBlock::default()
    })
}

#[test]
fn canonical_signatures_and_source_defaults_are_exact() {
    let _: fn(&Molecule) -> Result<DenseMatrix, MatrixError> = Molecule::distance_matrix;
    let _: for<'a, 'b> fn(
        &'a Molecule,
        &'b DistanceMatrixParams,
    ) -> Result<DenseMatrix, MatrixError> = Molecule::distance_matrix_with_params;
    let _: fn(&Molecule) -> Result<DenseMatrix, MatrixError> = Molecule::distance_matrix_3d;
    let _: for<'a, 'b> fn(
        &'a Molecule,
        &'b DistanceMatrix3dParams,
    ) -> Result<DenseMatrix, MatrixError> = Molecule::distance_matrix_3d_with_params;

    assert_eq!(
        DistanceMatrixParams::default(),
        DistanceMatrixParams {
            use_bond_order: false,
            use_atom_weights: false,
        }
    );
    assert_eq!(
        DistanceMatrix3dParams::default(),
        DistanceMatrix3dParams {
            conformer_id: None,
            use_atom_weights: false,
        }
    );
}

#[test]
fn binding_contract_has_exact_types_and_four_read_only_callables() {
    let expected = [
        "types.DenseMatrix",
        "types.DistanceMatrixParams",
        "types.DistanceMatrix3dParams",
        "types.MatrixError",
        "Molecule.distance_matrix",
        "Molecule.distance_matrix_with_params",
        "Molecule.distance_matrix_3d",
        "Molecule.distance_matrix_3d_with_params",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| expected.contains(&row.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    for row in &rows[..4] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
        assert_eq!(row.feature, "matrices");
    }
    for row in &rows[4..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
        assert_eq!(row.feature, "matrices");
        assert_eq!(row.callable.unwrap().state_model, StateModel::ReadOnly);
        assert_eq!(row.callable.unwrap().operation_semantic_id, None);
    }
    assert_eq!(rows[4].callable.unwrap().parameters.len(), 0);
    assert_eq!(rows[5].callable.unwrap().parameters.len(), 1);
    assert_eq!(rows[6].callable.unwrap().parameters.len(), 0);
    assert_eq!(rows[7].callable.unwrap().parameters.len(), 1);
}

#[test]
fn read_only_queries_have_no_operation_or_generated_matrix_rows() {
    for method in [
        "distance_matrix",
        "distance_matrix_with_params",
        "distance_matrix_3d",
        "distance_matrix_3d_with_params",
    ] {
        assert!(operation_spec(method).is_none());
        assert!(operation_invariant(method).is_none());
        assert!(operation_parity(method).is_none());
        assert!(support_matrix().iter().all(|row| {
            row.operation
                .is_none_or(|operation| operation.method != method)
        }));
    }
}

#[test]
fn topological_queries_match_full_detached_rows_and_do_not_mutate_live_state() {
    let source = molecule_3d();
    let before = source.clone();

    let unweighted = source.distance_matrix().unwrap();
    assert_eq!(unweighted.dimension(), 3);
    assert_eq!(
        unweighted.values(),
        &[0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0]
    );

    let weighted = source
        .distance_matrix_with_params(&DistanceMatrixParams {
            use_bond_order: true,
            use_atom_weights: true,
        })
        .unwrap();
    assert_eq!(weighted.dimension(), 3);
    assert_eq!(weighted.get(0, 0), Some(1.0));
    assert_eq!(weighted.get(1, 1), Some(0.75));
    assert_eq!(weighted.get(2, 2), Some(6.0 / 7.0));
    assert_eq!(weighted.get(0, 1), Some(1.0));
    assert_eq!(weighted.get(1, 2), Some(0.5));
    assert_eq!(weighted.get(0, 2), Some(1.5));

    assert_eq!(source, before);
    assert_eq!(source.property("source"), Some("preserved"));
    assert_eq!(
        source.topology().stereo_groups,
        before.topology().stereo_groups
    );
    assert_eq!(
        source.atom(AtomId::new(0)).unwrap().prop("_CIPCode"),
        Some("R")
    );
    assert!(std::ptr::eq(source.topology(), before.topology()));
    coordinate_views::assert_shared_coordinates(&source, &before);
    assert!(std::ptr::eq(source.properties(), before.properties()));
}

#[test]
fn three_dimensional_queries_select_by_id_and_preserve_exact_rows() {
    let source = molecule_3d();
    let before = source.clone();
    assert_eq!(
        source.distance_matrix_3d().unwrap().values(),
        &[0.0, 3.0, 5.0, 3.0, 0.0, 4.0, 5.0, 4.0, 0.0]
    );

    let selected = source
        .distance_matrix_3d_with_params(&DistanceMatrix3dParams {
            conformer_id: Some(3),
            use_atom_weights: true,
        })
        .unwrap();
    assert_eq!(selected.get(0, 0), Some(1.0));
    assert_eq!(selected.get(1, 1), Some(0.75));
    assert_eq!(selected.get(2, 2), Some(6.0 / 7.0));
    assert_eq!(selected.get(0, 1), Some(12.0));
    assert_eq!(selected.get(1, 2), Some(5.0));
    assert_eq!(selected.get(0, 2), Some(13.0));
    assert_eq!(source, before);
}

#[test]
fn three_dimensional_failures_remain_structured_and_leave_source_unchanged() {
    let no_coordinates = molecule_with_coordinates(CoordinateBlock::default());
    let before = no_coordinates.clone();
    assert_eq!(
        no_coordinates.distance_matrix_3d(),
        Err(MatrixError::No3dConformer)
    );
    assert_eq!(no_coordinates, before);

    let source = molecule_3d();
    assert_eq!(
        source.distance_matrix_3d_with_params(&DistanceMatrix3dParams {
            conformer_id: Some(404),
            use_atom_weights: false,
        }),
        Err(MatrixError::ConformerNotFound { conformer_id: 404 })
    );
}

#[test]
fn empty_public_results_keep_dense_matrix_bounds_behavior() {
    let empty = Molecule::new();
    let matrix = empty.distance_matrix().unwrap();
    assert_eq!(matrix.dimension(), 0);
    assert!(matrix.values().is_empty());
    assert_eq!(matrix.get(0, 0), None);
}
