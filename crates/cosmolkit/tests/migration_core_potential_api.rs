#[path = "support/coordinate_views.rs"]
mod coordinate_views;

use std::error::Error as _;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondDirection, BondId, BondOrder, BondSpec,
    BondStereo, ChiralTag, Conformer2D, CoordinateBlock, Element, Molecule, MoleculeOpKind,
    MoleculeOpOutput, MoleculeProperties, OperationDomain, OperationError, ParityPolicy,
    PotentialStereoCenter, PotentialStereoError, PotentialStereoParams, PotentialStereoResult,
    PotentialStereoSpecified, StateModel, StereoGroup, StereoGroupKind, SupportStatus,
    TopologyBlock, TopologyEditKind, feature_spec, operation_invariant, operation_parity,
    operation_spec, support_matrix,
};
use cosmolkit_core::{
    RingSearchParams, ValenceModel, ValenceParams, assign_valence, potential_stereo,
    symmetrized_sssr,
};

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn bond(index: usize, begin: usize, end: usize, order: BondOrder) -> Bond {
    Bond::from_spec(
        BondId::new(index),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
    )
}

fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}

fn four_distinct_ligands() -> Molecule {
    let topology = topology(
        vec![
            atom(0, Element::C),
            atom(1, Element::F),
            atom(2, Element::CL),
            atom(3, Element::BR),
            atom(4, Element::I),
        ],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 0, 2, BondOrder::Single),
            bond(2, 0, 3, BondOrder::Single),
            bond(3, 0, 4, BondOrder::Single),
        ],
    );
    Molecule::from_parts(
        topology,
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                7,
                vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]],
            )],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default()
            .with_name("potential-public")
            .with_prop("source", "preserved")
            .unwrap(),
    )
    .unwrap()
}

fn duplicate_ligands_with_stereo_state() -> Molecule {
    let center = AtomSpec::new(Element::C)
        .with_chiral_tag(ChiralTag::TetrahedralCw)
        .with_prop("atom-label", "center")
        .unwrap()
        .with_computed_prop("_CIPCode", "R")
        .unwrap();
    let wedge = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
        .with_direction(BondDirection::BeginWedge)
        .with_prop("bond-label", "wedge")
        .unwrap()
        .with_computed_prop("_CIPBondCode", "E")
        .unwrap();
    let topology = TopologyBlock::try_from_parts(
        vec![
            Atom::from_spec(AtomId::new(0), center),
            atom(1, Element::F),
            atom(2, Element::F),
            atom(3, Element::CL),
            atom(4, Element::BR),
        ],
        vec![
            Bond::from_spec(BondId::new(0), wedge),
            bond(1, 0, 2, BondOrder::Single),
            bond(2, 0, 3, BondOrder::Single),
            bond(3, 0, 4, BondOrder::Single),
        ],
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
                3,
                vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]],
            )],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default()
            .with_name("potential-clean")
            .with_prop("source", "preserved")
            .unwrap()
            .with_computed_prop("_CIPComputed", "true")
            .unwrap(),
    )
    .unwrap()
}

fn detached_result(
    topology: &TopologyBlock,
    params: &PotentialStereoParams,
) -> cosmolkit_core::PotentialStereoAssignment {
    let valence = assign_valence(
        topology,
        &ValenceParams {
            model: ValenceModel::RdkitLike,
            strict: false,
        },
    )
    .unwrap();
    let rings = symmetrized_sssr(topology, &RingSearchParams::default()).unwrap();
    potential_stereo(topology, &valence, &rings, params).unwrap()
}

#[test]
fn canonical_public_signatures_and_defaults_are_exact() {
    let _: fn(&Molecule) -> Result<PotentialStereoResult, OperationError> =
        Molecule::potential_stereo;
    let _: for<'a, 'b> fn(
        &'a Molecule,
        &'b PotentialStereoParams,
    ) -> Result<PotentialStereoResult, OperationError> = Molecule::potential_stereo_with_params;
    assert_eq!(
        PotentialStereoParams::default(),
        PotentialStereoParams {
            clean: false,
            flag_possible: true,
            allow_nontetrahedral: true,
        }
    );
}

#[test]
fn binding_contract_exposes_exactly_nine_types_and_two_callables() {
    let expected = [
        "types.PotentialStereoParams",
        "types.PotentialStereoType",
        "types.PotentialStereoSpecified",
        "types.PotentialStereoDescriptor",
        "types.PotentialStereoCenter",
        "types.PotentialStereoInfo",
        "types.RingStereoRelation",
        "types.PotentialStereoResult",
        "types.PotentialStereoError",
        "Molecule.potential_stereo",
        "Molecule.potential_stereo_with_params",
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
    for row in &rows[..8] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    assert_eq!(rows[8].item, BindingItem::Type);
    assert_eq!(rows[8].owner, BindingOwner::Type);
    assert_eq!(rows[8].support, BindingSupport::Supported);
    assert_eq!(rows[8].parity, BindingParity::NotApplicable);
    for row in &rows[9..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    assert_eq!(rows[9].javascript_name, "potentialStereo");
    assert_eq!(rows[10].javascript_name, "potentialStereoWithParams");
    assert_eq!(rows[9].callable.unwrap().parameters.len(), 0);
    assert_eq!(rows[10].callable.unwrap().parameters.len(), 1);
    assert_eq!(rows[9].callable.unwrap().state_model, StateModel::ReadOnly);
    assert_eq!(
        rows[10].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
}

#[test]
fn generated_registry_and_all_four_matrices_share_one_single_result_operation() {
    let feature = feature_spec("stereo").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("potential_stereo_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "crate :: PotentialStereoResult");
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::None);
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
        "weak_potential_stereo_cleanup"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "find_potential_stereo_rdkit"
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
fn both_public_callers_match_the_detached_oracle_and_preserve_input_order() {
    let source = four_distinct_ligands();
    let expected = detached_result(source.topology(), &PotentialStereoParams::default());
    let short = source.potential_stereo().unwrap();
    let explicit = source
        .potential_stereo_with_params(&PotentialStereoParams::default())
        .unwrap();
    assert_eq!(short.stereo, expected.stereo);
    assert_eq!(short.atom_ranks, expected.atom_ranks);
    assert_eq!(short.ring_relations, expected.ring_relations);
    assert!(short.cleaned_molecule.is_none());
    assert_eq!(short, explicit);
    assert_eq!(
        short.stereo[0].centered_on,
        PotentialStereoCenter::Atom(AtomId::new(0))
    );
    assert_eq!(
        short.stereo[0].specified,
        PotentialStereoSpecified::Unspecified
    );
    assert_eq!(short.stereo[0].controlling_atoms.len(), 4);
    assert_eq!(source.property("source"), Some("preserved"));
    assert_eq!(
        source.to_builder().coordinates().conformers_2d[0]
            .coordinates()
            .len(),
        5
    );
}

#[test]
fn clean_result_commits_only_source_defined_stereo_cleanup_and_cip_effects() {
    let source = duplicate_ligands_with_stereo_state();
    let observer = source.clone();
    let params = PotentialStereoParams {
        clean: true,
        ..PotentialStereoParams::default()
    };
    let mut expected = detached_result(source.topology(), &params);
    let result = source.potential_stereo_with_params(&params).unwrap();
    let cleaned = result.cleaned_molecule.as_ref().unwrap();

    let mut expected_topology = expected.cleaned_topology.take().unwrap();
    for atom in &mut expected_topology.atoms {
        atom.clear_computed_props();
    }
    for bond in &mut expected_topology.bonds {
        bond.clear_computed_props();
    }
    assert_eq!(cleaned.topology(), &expected_topology);
    assert_eq!(result.stereo, expected.stereo);
    assert_eq!(result.atom_ranks, expected.atom_ranks);
    assert_eq!(result.ring_relations, expected.ring_relations);
    assert_eq!(
        cleaned.to_builder().coordinates(),
        source.to_builder().coordinates()
    );
    coordinate_views::assert_shared_coordinates(&cleaned, &source);
    assert_eq!(cleaned.property("source"), Some("preserved"));
    assert_eq!(cleaned.property("_CIPComputed"), None);
    assert_eq!(
        cleaned.atom(AtomId::new(0)).unwrap().prop("atom-label"),
        Some("center")
    );
    assert_eq!(cleaned.atom(AtomId::new(0)).unwrap().prop("_CIPCode"), None);
    assert_eq!(
        cleaned.bond(BondId::new(0)).unwrap().prop("bond-label"),
        Some("wedge")
    );
    assert_eq!(
        cleaned.bond(BondId::new(0)).unwrap().prop("_CIPBondCode"),
        None
    );
    assert_eq!(
        cleaned.topology().stereo_groups,
        source.topology().stereo_groups
    );

    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(std::ptr::eq(source.properties(), observer.properties()));
    assert_eq!(
        source.atom(AtomId::new(0)).unwrap().chiral_tag(),
        ChiralTag::TetrahedralCw
    );
    assert_eq!(
        source.bond(BondId::new(0)).unwrap().direction(),
        BondDirection::BeginWedge
    );
    assert_eq!(source.property("_CIPComputed"), Some("true"));
}

#[test]
fn clean_and_flag_possible_combinations_have_exact_result_assembly() {
    let source = four_distinct_ligands();
    for clean in [false, true] {
        for flag_possible in [false, true] {
            let params = PotentialStereoParams {
                clean,
                flag_possible,
                allow_nontetrahedral: true,
            };
            let expected = detached_result(source.topology(), &params);
            let actual = source.potential_stereo_with_params(&params).unwrap();
            assert_eq!(actual.stereo, expected.stereo);
            assert_eq!(actual.atom_ranks, expected.atom_ranks);
            assert_eq!(actual.ring_relations, expected.ring_relations);
            assert_eq!(actual.cleaned_molecule.is_some(), clean);
        }
    }
}

#[test]
fn typed_valence_and_potential_stereo_failures_leave_the_source_unchanged() {
    let invalid_valence = Molecule::from_parts(
        topology(
            vec![atom(0, Element::C), atom(1, Element::C)],
            vec![bond(0, 0, 1, BondOrder::ThreeCenter)],
        ),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap();
    let observer = invalid_valence.clone();
    let error = invalid_valence.potential_stereo().unwrap_err();
    assert!(matches!(error, OperationError::Valence(_)));
    assert!(error.source().is_some());
    assert_eq!(invalid_valence, observer);
    assert!(std::ptr::eq(
        invalid_valence.topology(),
        observer.topology()
    ));

    let atrop_topology = topology(
        vec![atom(0, Element::C), atom(1, Element::C)],
        vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_stereo(BondStereo::AtropCw),
        )],
    );
    let atrop = Molecule::from_parts(
        atrop_topology,
        CoordinateBlock::default(),
        MoleculeProperties::default()
            .with_prop("source", "atrop")
            .unwrap(),
    )
    .unwrap();
    let observer = atrop.clone();
    assert_eq!(
        atrop.potential_stereo(),
        Err(OperationError::PotentialStereo(
            PotentialStereoError::AtropisomerDependencyUnavailable {
                bond: BondId::new(0),
            }
        ))
    );
    assert_eq!(atrop, observer);
    assert!(std::ptr::eq(atrop.topology(), observer.topology()));
    assert_eq!(atrop.property("source"), Some("atrop"));
}

#[test]
fn repeated_calls_are_deterministic_and_no_inplace_or_multi_output_surface_exists() {
    let source = four_distinct_ligands();
    let first = source.potential_stereo().unwrap();
    let second = source.potential_stereo().unwrap();
    assert_eq!(first, second);
    let spec = operation_spec("potential_stereo_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert!(
        BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id.starts_with("Molecule.potential_stereo"))
            .all(|row| row.callable.unwrap().state_model != StateModel::InPlace)
    );
}
