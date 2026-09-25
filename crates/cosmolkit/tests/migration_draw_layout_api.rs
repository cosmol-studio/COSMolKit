use std::collections::BTreeMap;
use std::error::Error as _;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec, Conformer2D,
    Conformer3D, Coordinate2DError, Coordinate2DLayoutError, Coordinate2DParams, CoordinateBlock,
    CoordinateDimension, Element, Molecule, MoleculeOpKind, MoleculeOpOutput, MoleculeProperties,
    OperationDomain, OperationError, ParityPolicy, StateModel, SupportStatus, TopologyBlock,
    TopologyEditKind, feature_spec, operation_invariant, operation_parity, operation_spec,
    support_matrix,
};

fn find_error_source<'a, T: std::error::Error + 'static>(
    error: &'a (dyn std::error::Error + 'static),
) -> Option<&'a T> {
    let mut current = Some(error);
    while let Some(source) = current {
        if let Some(typed) = source.downcast_ref::<T>() {
            return Some(typed);
        }
        current = source.source();
    }
    None
}

fn chain_molecule() -> Molecule {
    let atoms = (0..3)
        .map(|index| {
            Atom::from_spec(
                AtomId::new(index),
                AtomSpec::new(Element::C)
                    .with_prop("atom-label", format!("a{index}"))
                    .unwrap()
                    .with_computed_prop("_CIPCode", "R")
                    .unwrap(),
            )
        })
        .collect();
    let bonds = vec![
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        ),
        Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
        ),
    ];
    let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(7, vec![[10.0, 10.0], [11.0, 10.0], [12.0, 10.0]])
                .with_prop("source", "old-2d"),
        ],
        conformers_3d: vec![
            Conformer3D::new(
                4,
                vec![[0.0, 0.0, -0.0], [1.0, 0.0, 0.5], [2.0, 0.0, 1.0]],
                true,
            )
            .with_prop("source", "three-d"),
        ],
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    };
    let properties = MoleculeProperties::default()
        .with_name("draw-public")
        .with_prop("source", "preserved")
        .unwrap()
        .with_computed_prop("_CIPComputed", "true")
        .unwrap();
    Molecule::from_parts(topology, coordinates, properties).unwrap()
}

fn chain_molecule_with_conformer_ids(two_d_ids: &[usize], three_d_ids: &[usize]) -> Molecule {
    let source = chain_molecule();
    let atom_count = source.num_atoms();
    let coordinates = CoordinateBlock {
        conformers_2d: two_d_ids
            .iter()
            .map(|&id| Conformer2D::new(id, vec![[id as f64, 0.0]; atom_count]))
            .collect(),
        conformers_3d: three_d_ids
            .iter()
            .map(|&id| {
                Conformer3D::new(id, vec![[0.0, id as f64, -0.0]; atom_count], true)
                    .with_prop("source", format!("three-d-{id}"))
            })
            .collect(),
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    };
    Molecule::from_parts(
        source.topology().clone(),
        coordinates,
        source.properties().clone(),
    )
    .unwrap()
}

fn disconnected_molecule() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    Molecule::from_parts(
        topology,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap()
}

#[test]
fn short_and_parameterized_with_2d_coordinates() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::with_2d_coordinates;
    let _: for<'a, 'b> fn(
        &'a Molecule,
        &'b Coordinate2DParams,
    ) -> Result<Molecule, OperationError> = Molecule::with_2d_coordinates_with_params;
    let _: fn(&str) -> Option<&'static cosmolkit::FeatureSpec> = feature_spec;

    let source = chain_molecule();
    let short = source.with_2d_coordinates().unwrap();
    let configured = source
        .with_2d_coordinates_with_params(&Coordinate2DParams::default())
        .unwrap();
    assert_eq!(short.coordinates_2d(), configured.coordinates_2d());
    assert_eq!(
        short.coordinates_2d().unwrap(),
        &[
            [0.0, 0.0],
            [1.299038105676658, 0.7499999999999998],
            [2.598076211353316, -6.661338147750939e-16],
        ]
    );
    assert_eq!(
        Coordinate2DParams::default().coordinate_map,
        BTreeMap::new()
    );
}

#[test]
fn all_frozen_layout_options_and_constraint_errors() {
    let source = chain_molecule();
    let params = Coordinate2DParams {
        coordinate_map: BTreeMap::from([(0, [-2.0, 1.0]), (1, [0.0, 1.0]), (2, [2.0, 1.0])]),
        canonical_orientation: true,
        clear_existing_2d: false,
        flips_per_sample: 1,
        samples: 1,
        sample_seed: 17,
        permute_degree_four: true,
        force_rdkit: true,
        use_ring_templates: true,
    };
    let output = source.with_2d_coordinates_with_params(&params).unwrap();
    let stored = output.to_builder().coordinates().clone();
    assert_eq!(stored.conformers_2d.len(), 2);
    assert_eq!(stored.conformers_2d[0].id(), 7);
    assert_eq!(stored.conformers_2d[1].id(), 8);
    assert_eq!(
        stored.conformers_2d[1].coordinates(),
        &[[-2.0, 1.0], [0.0, 1.0], [2.0, 1.0]]
    );

    let mut invalid = Coordinate2DParams::default();
    invalid.coordinate_map.insert(3, [0.0, 0.0]);
    let error = source
        .with_2d_coordinates_with_params(&invalid)
        .unwrap_err();
    assert!(matches!(error, OperationError::Coordinate2D(_)));
    assert!(error.source().is_some());
}

#[test]
fn two_d_coordinate_rows_and_three_d_preservation() {
    let source = chain_molecule();
    let output = source.with_2d_coordinates().unwrap();
    let source_block = source.to_builder().coordinates().clone();
    let output_block = output.to_builder().coordinates().clone();

    assert_eq!(source_block.conformers_2d.len(), 1);
    assert_eq!(output_block.conformers_2d.len(), 1);
    assert_eq!(output_block.conformers_2d[0].id(), 0);
    assert_eq!(output_block.conformers_2d[0].coordinates().len(), 3);
    assert_eq!(&output_block.conformers_3d, &source_block.conformers_3d);
    assert_eq!(
        output_block.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(source.coordinates_2d().unwrap()[0], [10.0, 10.0]);
}

#[test]
fn conformer_id_empty_and_replace_paths_assign_zero() {
    let empty = chain_molecule_with_conformer_ids(&[], &[41]);
    let mut append = Coordinate2DParams::default();
    append.clear_existing_2d = false;
    let appended = empty.with_2d_coordinates_with_params(&append).unwrap();
    let appended_block = appended.to_builder().coordinates().clone();
    assert_eq!(
        appended_block
            .conformers_2d
            .iter()
            .map(Conformer2D::id)
            .collect::<Vec<_>>(),
        [0]
    );
    assert_eq!(appended_block.conformers_3d[0].id(), 41);

    let existing = chain_molecule_with_conformer_ids(&[9, 2], &[41]);
    let replaced = existing.with_2d_coordinates().unwrap();
    let replaced_block = replaced.to_builder().coordinates().clone();
    assert_eq!(
        replaced_block
            .conformers_2d
            .iter()
            .map(Conformer2D::id)
            .collect::<Vec<_>>(),
        [0]
    );
    assert_eq!(replaced_block.conformers_3d, existing.conformers_3d());
}

#[test]
fn conformer_id_append_uses_two_d_max_plus_one_and_preserves_three_d() {
    let source = chain_molecule_with_conformer_ids(&[9, 2, 7], &[usize::MAX, 3]);
    let source_block = source.to_builder().coordinates().clone();
    let mut params = Coordinate2DParams::default();
    params.clear_existing_2d = false;

    let output = source.with_2d_coordinates_with_params(&params).unwrap();
    let output_block = output.to_builder().coordinates().clone();
    assert_eq!(
        output_block
            .conformers_2d
            .iter()
            .map(Conformer2D::id)
            .collect::<Vec<_>>(),
        [9, 2, 7, 10]
    );
    assert_eq!(output_block.conformers_3d, source_block.conformers_3d);
    assert_eq!(source.to_builder().coordinates(), &source_block);
}

#[test]
fn conformer_id_overflow_is_structured_and_failure_atomic() {
    let source = chain_molecule_with_conformer_ids(&[usize::MAX], &[0]);
    let before = source.clone();
    let before_block = source.to_builder().coordinates().clone();
    let mut params = Coordinate2DParams::default();
    params.clear_existing_2d = false;

    let error = source
        .with_2d_coordinates_with_params(&params)
        .expect_err("max 2D conformer id must not wrap or panic");
    assert!(matches!(error, OperationError::Coordinate2D(_)));
    assert!(error.source().is_some());
    assert_eq!(source, before);
    assert_eq!(source.to_builder().coordinates(), &before_block);
}

#[test]
fn value_semantics_unchanged_blocks_and_failure_atomicity() {
    let source = chain_molecule();
    let output = source.with_2d_coordinates().unwrap();

    assert!(std::ptr::eq(source.topology(), output.topology()));
    assert!(std::ptr::eq(source.properties(), output.properties()));
    assert_eq!(output.property("source"), Some("preserved"));
    assert_eq!(output.property("_CIPComputed"), Some("true"));
    assert_eq!(
        output.atom(AtomId::new(0)).unwrap().prop("_CIPCode"),
        Some("R")
    );
    assert_eq!(source.coordinates_2d().unwrap()[0], [10.0, 10.0]);

    let before = source.clone();
    let mut invalid = Coordinate2DParams::default();
    invalid.coordinate_map.insert(99, [1.0, 2.0]);
    assert!(source.with_2d_coordinates_with_params(&invalid).is_err());
    assert_eq!(source, before);
    assert!(std::ptr::eq(source.topology(), before.topology()));
    assert!(std::ptr::eq(source.properties(), before.properties()));
    assert_eq!(source.coordinates_2d(), before.coordinates_2d());
}

#[test]
fn live_typed_error_chain_preserves_layout_category_and_source_object() {
    let source = chain_molecule();
    let before = source.clone();
    let mut params = Coordinate2DParams::default();
    params.coordinate_map.insert(3, [0.0, 0.0]);

    let error = source
        .with_2d_coordinates_with_params(&params)
        .expect_err("an out-of-range coordinate-map atom must remain an error");
    assert!(matches!(&error, OperationError::Coordinate2D(_)));
    let operation_source = error
        .source()
        .expect("OperationError must retain the coordinate-domain error");
    assert!(
        operation_source
            .downcast_ref::<Coordinate2DError>()
            .is_some()
    );
    assert!(matches!(
        find_error_source::<Coordinate2DLayoutError>(&error),
        Some(Coordinate2DLayoutError::AtomIndexOutOfRange {
            atom: 3,
            atom_count: 3
        })
    ));
    assert_eq!(source, before);
}

#[test]
fn live_typed_error_chain_preserves_sampling_safety_boundary() {
    let source = disconnected_molecule();
    let before = source.clone();
    let params = Coordinate2DParams {
        flips_per_sample: 1,
        samples: 1,
        sample_seed: 7,
        ..Default::default()
    };

    let error = source
        .with_2d_coordinates_with_params(&params)
        .expect_err("the CK boundary must not invent undefined sampling costs");
    assert!(matches!(&error, OperationError::Coordinate2D(_)));
    assert!(error.source().is_some());
    assert!(matches!(
        find_error_source::<Coordinate2DLayoutError>(&error),
        Some(Coordinate2DLayoutError::UndefinedSamplingDistance {
            first: 1,
            second: 0
        })
    ));
    assert_eq!(source, before);
}

#[test]
fn live_coordinate_failure_after_checkout_is_atomic_and_typed() {
    let source = chain_molecule_with_conformer_ids(&[usize::MAX], &[4]);
    let before = source.clone();
    let before_coordinates = source.to_builder().coordinates().clone();
    let mut params = Coordinate2DParams::default();
    params.clear_existing_2d = false;

    let error = source
        .with_2d_coordinates_with_params(&params)
        .expect_err("incrementing the maximum 2D conformer ID must fail without wrapping");
    assert!(matches!(
        &error,
        OperationError::Coordinate2D(Coordinate2DError::ConformerIdOverflow { max_id })
            if *max_id == usize::MAX
    ));
    assert!(
        error
            .source()
            .and_then(|source| source.downcast_ref::<Coordinate2DError>())
            .is_some()
    );
    assert_eq!(source, before);
    assert_eq!(source.to_builder().coordinates(), &before_coordinates);
}

#[test]
fn binding_registry_operation_contract_and_feature_isolation() {
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "depict")
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        [
            "types.Coordinate2DParams",
            "types.Coordinate2DError",
            "types.Coordinate2DTemplateError",
            "types.Coordinate2DLayoutError",
            "Molecule.with_2d_coordinates",
            "Molecule.with_2d_coordinates_with_params",
        ]
    );
    for row in &rows {
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::Experimental);
        assert_eq!(row.parity, BindingParity::RequiredWhenSupported);
    }
    for row in &rows[..4] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
    }
    for row in &rows[4..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(
            row.callable.unwrap().state_model,
            StateModel::ValueReturning
        );
    }
    assert_eq!(rows[4].callable.unwrap().parameters.len(), 0);
    assert_eq!(rows[5].callable.unwrap().parameters.len(), 1);

    let feature = feature_spec("depict").unwrap();
    assert_eq!(feature.status, SupportStatus::Experimental);
    let spec = operation_spec("with_2d_coordinates_with_params").unwrap();
    assert_eq!(spec.domain, OperationDomain::Coordinate);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::None);
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.access.read(), BlockSet::TOPOLOGY);
    assert_eq!(
        spec.access.write(),
        BlockSet::COORDINATES.union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(spec.may_mutate, spec.access.write());
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    assert_eq!(spec.derived_effects.invalidate.bits(), (1 << 4) | (1 << 6));
    assert_eq!(format!("{:?}", spec.cip_state), "Preserve");
    assert_eq!(spec.support, SupportStatus::Experimental);
    assert_eq!(spec.parity, ParityPolicy::RequiredWhenSupported);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "coordinate_2d_layout"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "compute_2d_coordinates_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| row.feature.name == "depict")
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));

    let _: Coordinate2DError = cosmolkit::Coordinate2DError::CoordGenUnavailable;
    assert!(
        BINDING_CONTRACT
            .iter()
            .filter(|row| row.feature == "depict")
            .all(|row| !row.semantic_id.contains("svg") && !row.semantic_id.contains("png"))
    );
}
