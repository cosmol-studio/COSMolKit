#[path = "support/coordinate_views.rs"]
mod coordinate_views;

use std::error::Error as _;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec, ChiralTag,
    Conformer3D, CoordinateBlock, Element, Molecule, MoleculeOpKind, MoleculeOpOutput,
    MoleculeProperties, OperationDomain, OperationError, ParityPolicy, StateModel, StereoError,
    StructureTagParams, SupportStatus, TopologyBlock, TopologyEditKind, feature_spec,
    operation_invariant, operation_parity, operation_spec, support_matrix,
};
use cosmolkit_core::{
    ValenceModel, ValenceParams, assign_chiral_tags_from_structure, assign_valence,
};

fn atom(index: usize, spec: AtomSpec) -> Atom {
    Atom::from_spec(AtomId::new(index), spec)
}

fn bond(index: usize, end: usize) -> Bond {
    Bond::from_spec(
        BondId::new(index),
        BondSpec::new(AtomId::new(0), AtomId::new(end), BondOrder::Single),
    )
}

fn tetra_topology(initial_tag: ChiralTag) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        vec![
            atom(
                0,
                AtomSpec::new(Element::C)
                    .with_explicit_hydrogens(1)
                    .with_chiral_tag(initial_tag)
                    .with_prop("center-label", "preserved")
                    .unwrap()
                    .with_computed_prop("_CIPCode", "R")
                    .unwrap(),
            ),
            atom(1, AtomSpec::new(Element::F)),
            atom(2, AtomSpec::new(Element::CL)),
            atom(3, AtomSpec::new(Element::BR)),
        ],
        vec![bond(0, 1), bond(1, 2), bond(2, 3)],
        Vec::new(),
        Vec::new(),
    )
    .unwrap()
}

fn coordinates(id: usize, is_3d: bool) -> CoordinateBlock {
    CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            id,
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            is_3d,
        )],
        ..CoordinateBlock::default()
    }
}

fn molecule(initial_tag: ChiralTag, conformer_id: usize, is_3d: bool) -> Molecule {
    Molecule::from_parts(
        tetra_topology(initial_tag),
        coordinates(conformer_id, is_3d),
        MoleculeProperties::default()
            .with_name("structure-tag-public")
            .with_prop("source", "preserved")
            .unwrap()
            .with_computed_prop("_StereochemDone", "1")
            .unwrap()
            .with_computed_prop("_CIPComputed", "true")
            .unwrap(),
    )
    .unwrap()
}

fn detached_expected(source: &Molecule, params: &StructureTagParams) -> TopologyBlock {
    let valence = assign_valence(
        source.topology(),
        &ValenceParams {
            model: ValenceModel::RdkitLike,
            strict: false,
        },
    )
    .unwrap();
    let mut topology = assign_chiral_tags_from_structure(
        source.topology(),
        source.to_builder().coordinates(),
        &valence,
        params,
    )
    .unwrap()
    .topology;
    for atom in &mut topology.atoms {
        atom.clear_computed_props();
    }
    for bond in &mut topology.bonds {
        bond.clear_computed_props();
    }
    topology
}

#[test]
fn canonical_signatures_and_source_defaults_are_exact() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> =
        Molecule::with_chiral_tags_from_structure;
    let _: for<'a, 'b> fn(
        &'a Molecule,
        &'b StructureTagParams,
    ) -> Result<Molecule, OperationError> = Molecule::with_chiral_tags_from_structure_with_params;
    let _: fn(&mut Molecule) -> Result<(), OperationError> =
        Molecule::assign_chiral_tags_from_structure_;
    let _: for<'a, 'b> fn(&'a mut Molecule, &'b StructureTagParams) -> Result<(), OperationError> =
        Molecule::assign_chiral_tags_from_structure_with_params_;
    assert_eq!(
        StructureTagParams::default(),
        StructureTagParams {
            conformer_id: -1,
            replace_existing_tags: true,
        }
    );
}

#[test]
fn binding_contract_has_two_types_and_four_canonical_callables() {
    let expected = [
        "types.StructureTagParams",
        "types.StereoError",
        "Molecule.with_chiral_tags_from_structure",
        "Molecule.with_chiral_tags_from_structure_with_params",
        "Molecule.assign_chiral_tags_from_structure_",
        "Molecule.assign_chiral_tags_from_structure_with_params_",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| expected.contains(&row.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    assert!(rows.iter().all(|row| row.feature == "stereo"));
    assert!(
        rows.iter()
            .all(|row| row.exposure == BindingExposure::Public)
    );
    assert_eq!(rows[0].item, BindingItem::Type);
    assert_eq!(rows[0].owner, BindingOwner::Type);
    assert_eq!(rows[0].support, BindingSupport::SupportedWithRdkitParity);
    assert_eq!(rows[0].parity, BindingParity::RequiredNow);
    assert_eq!(rows[1].item, BindingItem::Type);
    assert_eq!(rows[1].owner, BindingOwner::Type);
    assert_eq!(rows[1].support, BindingSupport::Supported);
    assert_eq!(rows[1].parity, BindingParity::NotApplicable);
    for row in &rows[2..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    assert_eq!(rows[2].javascript_name, "withChiralTagsFromStructure");
    assert_eq!(
        rows[3].javascript_name,
        "withChiralTagsFromStructureWithParams"
    );
    assert_eq!(rows[4].javascript_name, "assignChiralTagsFromStructure");
    assert_eq!(
        rows[5].javascript_name,
        "assignChiralTagsFromStructureWithParams"
    );
    assert_eq!(
        rows[2].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(rows[3].callable.unwrap().parameters.len(), 1);
    assert_eq!(rows[4].callable.unwrap().state_model, StateModel::InPlace);
    assert_eq!(rows[5].callable.unwrap().parameters.len(), 1);
}

#[test]
fn generated_registry_and_all_four_matrices_describe_one_weak_operation() {
    let feature = feature_spec("stereo").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_chiral_tags_from_structure_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::None);
    assert_eq!(spec.access.read(), BlockSet::COORDINATES);
    assert_eq!(
        spec.access.write(),
        BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(spec.may_mutate, spec.access.write());
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    assert_eq!(spec.semantic_preconditions.bits(), 0);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.derived_effects.recompute.bits(), 0);
    assert_eq!(
        spec.derived_effects.preserve.bits(),
        (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 5)
    );
    assert_eq!(
        spec.derived_effects.invalidate.bits(),
        (1 << 4) | (1 << 6) | (1 << 7)
    );
    assert_eq!(format!("{:?}", spec.cip_state), "ClearComputed");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert!(spec.io_roundtrip);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "weak_structure_tag_assignment"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "assign_chiral_tags_from_structure_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| {
            row.operation
                .is_some_and(|operation| std::ptr::eq(operation, spec))
        })
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));
}

#[test]
fn value_and_inplace_defaults_match_the_detached_owner_and_are_deterministic() {
    let source = molecule(ChiralTag::Unspecified, 7, true);
    let observer = source.clone();
    let expected = detached_expected(&source, &StructureTagParams::default());

    let short = source.with_chiral_tags_from_structure().unwrap();
    let explicit = source
        .with_chiral_tags_from_structure_with_params(&StructureTagParams::default())
        .unwrap();
    assert_eq!(short, explicit);
    assert_eq!(short.topology(), &expected);
    assert_eq!(
        short.atom(AtomId::new(0)).unwrap().chiral_tag(),
        ChiralTag::TetrahedralCcw
    );
    assert_eq!(
        short.to_builder().coordinates(),
        source.to_builder().coordinates()
    );
    coordinate_views::assert_shared_coordinates(&short, &source);
    assert_eq!(short.property("source"), Some("preserved"));
    assert_eq!(short.property("_StereochemDone"), None);
    assert_eq!(short.property("_CIPComputed"), None);
    assert_eq!(
        short.atom(AtomId::new(0)).unwrap().prop("center-label"),
        Some("preserved")
    );
    assert_eq!(short.atom(AtomId::new(0)).unwrap().prop("_CIPCode"), None);

    let mut default_inplace = source.clone();
    default_inplace
        .assign_chiral_tags_from_structure_()
        .unwrap();
    assert_eq!(default_inplace, short);
    let mut explicit_inplace = source.clone();
    explicit_inplace
        .assign_chiral_tags_from_structure_with_params_(&StructureTagParams::default())
        .unwrap();
    assert_eq!(explicit_inplace, short);

    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(std::ptr::eq(source.properties(), observer.properties()));
    assert_eq!(source.property("_StereochemDone"), Some("1"));
}

#[test]
fn explicit_conformer_selection_and_replace_flag_follow_source_order() {
    let mut source = molecule(ChiralTag::TetrahedralCw, 5, false);
    let mut coordinate_block = source.to_builder().coordinates().clone();
    coordinate_block.conformers_3d.push(Conformer3D::new(
        9,
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        true,
    ));
    source = Molecule::from_parts(
        source.topology().clone(),
        coordinate_block,
        source.properties().clone(),
    )
    .unwrap();

    let first_non_3d = source.with_chiral_tags_from_structure().unwrap();
    assert_eq!(
        first_non_3d.atom(AtomId::new(0)).unwrap().chiral_tag(),
        ChiralTag::TetrahedralCw
    );
    let keep = source
        .with_chiral_tags_from_structure_with_params(&StructureTagParams {
            conformer_id: 9,
            replace_existing_tags: false,
        })
        .unwrap();
    assert_eq!(
        keep.atom(AtomId::new(0)).unwrap().chiral_tag(),
        ChiralTag::TetrahedralCw
    );
    let replace = source
        .with_chiral_tags_from_structure_with_params(&StructureTagParams {
            conformer_id: 9,
            replace_existing_tags: true,
        })
        .unwrap();
    assert_eq!(
        replace.atom(AtomId::new(0)).unwrap().chiral_tag(),
        ChiralTag::TetrahedralCcw
    );
    assert_eq!(
        replace.to_builder().coordinates(),
        source.to_builder().coordinates()
    );
    assert_eq!(replace.property("source"), Some("preserved"));
}

#[test]
fn typed_stereo_failure_is_atomic_for_value_and_inplace_forms() {
    let source = molecule(ChiralTag::Unspecified, 7, true);
    let observer = source.clone();
    let params = StructureTagParams {
        conformer_id: 8,
        replace_existing_tags: true,
    };
    let error = source
        .with_chiral_tags_from_structure_with_params(&params)
        .unwrap_err();
    assert_eq!(
        error,
        OperationError::Stereo(StereoError::ConformerNotFound { requested: 8 })
    );
    assert!(error.source().is_some());
    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(std::ptr::eq(source.properties(), observer.properties()));

    let mut target = source.clone();
    let inplace_observer = target.clone();
    assert_eq!(
        target.assign_chiral_tags_from_structure_with_params_(&params),
        Err(OperationError::Stereo(StereoError::ConformerNotFound {
            requested: 8,
        }))
    );
    assert_eq!(target, inplace_observer);
    assert!(std::ptr::eq(target.topology(), inplace_observer.topology()));
    coordinate_views::assert_shared_coordinates(&target, &inplace_observer);
    assert!(std::ptr::eq(
        target.properties(),
        inplace_observer.properties()
    ));
}

#[test]
fn public_family_has_no_legacy_free_function_or_multiple_output_shape() {
    let spec = operation_spec("with_chiral_tags_from_structure_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert!(
        BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id.contains("chiral_tags_from_structure"))
            .all(|row| row.owner == BindingOwner::Molecule)
    );
    let first = molecule(ChiralTag::Unspecified, 7, true)
        .with_chiral_tags_from_structure()
        .unwrap();
    let second = molecule(ChiralTag::Unspecified, 7, true)
        .with_chiral_tags_from_structure()
        .unwrap();
    assert_eq!(first, second);
}
