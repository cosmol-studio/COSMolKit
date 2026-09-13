use std::error::Error as _;

use cosmolkit::{
    AddHsParams, Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem,
    BindingOwner, BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec,
    ChiralTag, Conformer2D, Conformer3D, CoordinateBlock, Element, HydrogenError, Molecule,
    MoleculeOpKind, MoleculeOpOutput, MoleculeProperties, OperationDomain, OperationError,
    ParityPolicy, SdfPropertyList, SdfPropertyListTarget, StateModel, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, SupportStatus, TopologyBlock,
    TopologyEditKind, feature_spec, operation_invariant, operation_parity, operation_spec,
    support_matrix,
};

fn atom(index: usize, spec: AtomSpec) -> Atom {
    Atom::from_spec(AtomId::new(index), spec)
}

fn bond(index: usize, begin: usize, end: usize, order: BondOrder) -> Bond {
    Bond::from_spec(
        BondId::new(index),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
    )
}

fn stateful_source() -> Molecule {
    let atoms = vec![
        atom(
            0,
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(1)
                .with_no_implicit(true)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_chiral_permutation(4)
                .with_prop("atom-note", "a0")
                .unwrap()
                .with_computed_prop("_CIPCode", "R")
                .unwrap(),
        ),
        atom(
            1,
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_prop("atom-note", "a1")
                .unwrap(),
        ),
    ];
    let bonds = vec![Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
            .with_prop("bond-note", "b0")
            .unwrap()
            .with_computed_prop("_CIPCode", "E")
            .unwrap(),
    )];
    let substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(0)])
            .with_bonds(vec![BondId::new(0)])
            .with_label("kept"),
    ];
    let stereo_groups = vec![
        StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(0)],
            vec![BondId::new(0)],
        )
        .with_id(19),
    ];
    let topology =
        TopologyBlock::try_from_parts(atoms, bonds, substance_groups, stereo_groups).unwrap();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(7, vec![[1.0, 2.0], [3.0, 4.0]]).with_prop("plane", "kept"),
        ],
        conformers_3d: vec![
            Conformer3D::new(8, vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], true)
                .with_prop("space", "kept"),
        ],
        ..CoordinateBlock::default()
    };
    let properties = MoleculeProperties::default()
        .with_name("hydrogens-public")
        .with_sdf_data_field("field", "kept")
        .with_prop("ordinary", "kept")
        .unwrap()
        .with_computed_prop("_CIPComputed", "true")
        .unwrap()
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atom-values",
            vec![Some("a0".into()), Some("a1".into())],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bond-values",
            vec![Some("b0".into())],
        ));
    Molecule::from_parts(topology, coordinates, properties).unwrap()
}

fn explicit_params() -> AddHsParams {
    AddHsParams {
        explicit_only: true,
        ..AddHsParams::default()
    }
}

fn ring_source_with_cache() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![
            atom(
                0,
                AtomSpec::new(Element::C)
                    .with_explicit_hydrogens(1)
                    .with_no_implicit(true),
            ),
            atom(1, AtomSpec::new(Element::C).with_no_implicit(true)),
            atom(2, AtomSpec::new(Element::C).with_no_implicit(true)),
        ],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 1, 2, BondOrder::Single),
            bond(2, 2, 0, BondOrder::Single),
        ],
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
    .with_assigned_rings()
    .unwrap()
    .with_assigned_valence()
    .unwrap()
}

#[test]
fn canonical_public_signatures_and_defaults_compile() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::with_hydrogens;
    let _: for<'a, 'b> fn(&'a Molecule, &'b AddHsParams) -> Result<Molecule, OperationError> =
        Molecule::with_hydrogens_with_params;
    let _: fn(&mut Molecule) -> Result<(), OperationError> = Molecule::add_hydrogens_;
    let _: for<'a, 'b> fn(&'a mut Molecule, &'b AddHsParams) -> Result<(), OperationError> =
        Molecule::add_hydrogens_with_params_;
    assert_eq!(
        AddHsParams::default(),
        AddHsParams {
            explicit_only: false,
            add_coords: false,
            add_residue_info: false,
            skip_queries: false,
            only_on_atoms: None,
        }
    );
}

#[test]
fn binding_contract_has_exact_add_hydrogen_types_and_callables() {
    let expected = [
        "types.AddHsParams",
        "types.HydrogenError",
        "Molecule.with_hydrogens",
        "Molecule.with_hydrogens_with_params",
        "Molecule.add_hydrogens_",
        "Molecule.add_hydrogens_with_params_",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| expected.contains(&row.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    for row in &rows[..2] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::Supported);
        assert_eq!(row.parity, BindingParity::NotApplicable);
    }
    for row in &rows[2..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.feature, "hydrogens");
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
fn generated_registry_and_four_matrices_share_the_strong_append_operation() {
    let feature = feature_spec("hydrogens").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_hydrogens_with_params").unwrap();
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Strong);
    assert_eq!(spec.topology_edit, TopologyEditKind::Appending);
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
    assert_eq!(
        spec.auto_remap,
        BlockSet::COORDINATES.union(BlockSet::PROPERTIES)
    );
    assert_eq!(format!("{:?}", spec.requires_mapping), "Required");
    assert_eq!(spec.derived_effects.recompute.bits(), 0);
    assert_eq!(spec.derived_effects.preserve.bits(), (1 << 0) | (1 << 1));
    assert_eq!(
        spec.derived_effects.invalidate.bits(),
        (1 << 2) | (1 << 3) | (1 << 4) | (1 << 6) | (1 << 7)
    );
    assert_eq!(format!("{:?}", spec.cip_state), "ClearComputed");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "strong_topology_with_coordinates"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "add_hydrogens_rdkit"
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
fn live_value_commit_preserves_state_projects_rows_and_detaches_written_blocks() {
    let source = stateful_source();
    let snapshot = source.clone();
    let output = source
        .with_hydrogens_with_params(&explicit_params())
        .unwrap();

    assert_eq!(source, snapshot);
    assert_eq!(output.num_atoms(), 3);
    assert_eq!(output.num_bonds(), 2);
    assert_eq!(output.atoms()[0].id(), AtomId::new(0));
    assert_eq!(output.atoms()[1].id(), AtomId::new(1));
    assert_eq!(output.atoms()[2].id(), AtomId::new(2));
    assert_eq!(output.bonds()[0].id(), BondId::new(0));
    assert_eq!(output.bonds()[1].id(), BondId::new(1));
    assert_eq!(output.bonds()[1].begin(), AtomId::new(0));
    assert_eq!(output.bonds()[1].end(), AtomId::new(2));
    assert_eq!(output.atoms()[0].chiral_tag(), ChiralTag::TetrahedralCw);
    assert_eq!(output.atoms()[0].chiral_permutation(), Some(4));
    assert_eq!(
        output.topology().substance_groups,
        source.topology().substance_groups
    );
    assert_eq!(
        output.topology().stereo_groups,
        source.topology().stereo_groups
    );
    assert_eq!(output.atoms()[0].prop("atom-note"), Some("a0"));
    assert_eq!(output.atoms()[1].prop("atom-note"), Some("a1"));
    assert_eq!(output.bonds()[0].prop("bond-note"), Some("b0"));
    assert_eq!(output.atoms()[0].prop("_CIPCode"), None);
    assert_eq!(output.bonds()[0].prop("_CIPCode"), None);
    assert_eq!(output.property("_CIPComputed"), None);
    assert_eq!(output.property("ordinary"), Some("kept"));
    assert_eq!(output.properties().name(), Some("hydrogens-public"));
    assert_eq!(
        output.properties().sdf_data_fields(),
        &[("field".into(), "kept".into())]
    );
    let lists = output.properties().sdf_property_lists();
    assert_eq!(
        lists[0].values(),
        &[Some("a0".into()), Some("a1".into()), None]
    );
    assert_eq!(lists[1].values(), &[Some("b0".into()), None]);

    assert_eq!(output.coordinates().conformers_2d[0].id(), 7);
    assert_eq!(
        output.coordinates().conformers_2d[0].coordinates(),
        &[[1.0, 2.0], [3.0, 4.0], [0.0, 0.0]]
    );
    assert_eq!(
        output.coordinates().conformers_2d[0]
            .props()
            .get("plane")
            .map(String::as_str),
        Some("kept")
    );
    assert_eq!(output.coordinates().conformers_3d[0].id(), 8);
    assert_eq!(
        output.coordinates().conformers_3d[0].coordinates(),
        &[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [0.0, 0.0, 0.0]]
    );
    assert!(!std::ptr::eq(source.topology(), output.topology()));
    assert!(!std::ptr::eq(source.coordinates(), output.coordinates()));
    assert!(!std::ptr::eq(source.properties(), output.properties()));
}

#[test]
fn strict_leaf_append_preserves_ring_cache_and_invalidates_other_declared_state() {
    let source = ring_source_with_cache();
    assert!(format!("{source:?}").contains("derived_cache_is_empty: false"));
    let output = source
        .with_hydrogens_with_params(&explicit_params())
        .unwrap();
    assert_eq!(output.num_atoms(), 4);
    assert_eq!(output.num_bonds(), 4);
    assert!(format!("{output:?}").contains("derived_cache_is_empty: false"));
    let ring_rechecked = output.with_assigned_rings().unwrap();
    assert_eq!(ring_rechecked, output);
    assert!(std::ptr::eq(output.topology(), ring_rechecked.topology()));
}

#[test]
fn default_parameterized_and_inplace_forms_have_identical_single_output_order() {
    let source = stateful_source();
    let default_output = source.with_hydrogens().unwrap();
    let explicit_default = source
        .with_hydrogens_with_params(&AddHsParams::default())
        .unwrap();
    assert_eq!(default_output, explicit_default);

    let expected = source
        .with_hydrogens_with_params(&explicit_params())
        .unwrap();
    let mut short = source.clone();
    short.add_hydrogens_().unwrap();
    assert_eq!(short, default_output);
    let mut explicit = source.clone();
    explicit
        .add_hydrogens_with_params_(&explicit_params())
        .unwrap();
    assert_eq!(explicit, expected);
    assert_eq!(
        explicit.atoms().iter().map(Atom::id).collect::<Vec<_>>(),
        vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
    );
}

#[test]
fn structured_parameter_failure_is_atomic_for_value_and_inplace_wrappers() {
    let source = stateful_source();
    let observer = source.clone();
    let params = AddHsParams {
        only_on_atoms: Some(vec![AtomId::new(99)]),
        ..AddHsParams::default()
    };
    let error = source.with_hydrogens_with_params(&params).unwrap_err();
    assert!(matches!(
        &error,
        OperationError::Hydrogen(HydrogenError::OnlyOnAtomOutOfRange {
            atom,
            atom_count: 2,
        }) if *atom == AtomId::new(99)
    ));
    assert!(error.source().is_some());
    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    assert!(std::ptr::eq(source.coordinates(), observer.coordinates()));
    assert!(std::ptr::eq(source.properties(), observer.properties()));

    let mut target = source.clone();
    let target_observer = target.clone();
    let error = target.add_hydrogens_with_params_(&params).unwrap_err();
    assert!(matches!(
        error,
        OperationError::Hydrogen(HydrogenError::OnlyOnAtomOutOfRange { .. })
    ));
    assert_eq!(target, source);
    assert!(std::ptr::eq(target.topology(), target_observer.topology()));
    assert!(std::ptr::eq(
        target.coordinates(),
        target_observer.coordinates()
    ));
    assert!(std::ptr::eq(
        target.properties(),
        target_observer.properties()
    ));
}
