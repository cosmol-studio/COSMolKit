use std::error::Error as _;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec, Conformer2D,
    CoordinateBlock, Element, Molecule, MoleculeOpKind, MoleculeOpOutput, MoleculeProperties,
    OperationDomain, OperationError, ParityPolicy, StateModel, StereoGroup, StereoGroupKind,
    SupportStatus, TopologyBlock, TopologyEditKind, ValenceError, ValenceModel, ValenceParams,
    feature_spec, operation_invariant, operation_parity, operation_spec, support_matrix,
};

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn molecule() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![
            atom(0, Element::C),
            atom(1, Element::C),
            atom(2, Element::O),
        ],
        vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            ),
        ],
        Vec::new(),
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(1)],
            vec![BondId::new(0)],
        )],
    )
    .unwrap();
    Molecule::from_parts(
        topology,
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                7,
                vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
            )],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default()
            .with_name("valence-public")
            .with_prop("source", "preserved")
            .unwrap(),
    )
    .unwrap()
}

fn invalid_pentavalent_carbon() -> Molecule {
    let mut atoms = vec![atom(0, Element::C)];
    atoms.extend((1..=5).map(|index| {
        Atom::from_spec(
            AtomId::new(index),
            AtomSpec::new(Element::F).with_no_implicit(true),
        )
    }));
    let bonds = (1..=5)
        .map(|index| {
            Bond::from_spec(
                BondId::new(index - 1),
                BondSpec::new(AtomId::new(0), AtomId::new(index), BondOrder::Single),
            )
        })
        .collect();
    Molecule::from_parts(
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap(),
        CoordinateBlock::default(),
        MoleculeProperties::default().with_name("invalid-valence"),
    )
    .unwrap()
}

#[test]
fn canonical_public_signatures_and_defaults_compile() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::with_assigned_valence;
    let _: for<'a, 'b> fn(&'a Molecule, &'b ValenceParams) -> Result<Molecule, OperationError> =
        Molecule::with_assigned_valence_with_params;
    let _: fn(&mut Molecule) -> Result<(), OperationError> = Molecule::assign_valence_;
    let _: for<'a, 'b> fn(&'a mut Molecule, &'b ValenceParams) -> Result<(), OperationError> =
        Molecule::assign_valence_with_params_;
    let _: fn(&Molecule, AtomId) -> Result<bool, ValenceError> = Molecule::has_valence_violation;
    assert_eq!(
        ValenceParams::default(),
        ValenceParams {
            model: ValenceModel::RdkitLike,
            strict: true,
        }
    );
}

#[test]
fn binding_contract_exposes_the_frozen_eight_entries() {
    let expected = [
        "types.ValenceModel",
        "types.ValenceParams",
        "types.ValenceError",
        "Molecule.with_assigned_valence",
        "Molecule.with_assigned_valence_with_params",
        "Molecule.assign_valence_",
        "Molecule.assign_valence_with_params_",
        "Molecule.has_valence_violation",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "valence")
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    for row in &rows[..3] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::Supported);
        assert_eq!(row.parity, BindingParity::NotApplicable);
    }
    for row in &rows[3..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    assert_eq!(
        rows[3].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(rows[5].callable.unwrap().state_model, StateModel::InPlace);
    assert_eq!(rows[7].callable.unwrap().state_model, StateModel::ReadOnly);
}

#[test]
fn generated_registry_and_all_four_matrices_share_one_operation() {
    let feature = feature_spec("valence").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_assigned_valence_with_params").unwrap();
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::None);
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.access.read(), BlockSet::TOPOLOGY);
    assert_eq!(spec.access.write(), BlockSet::DERIVED_CACHE);
    assert_eq!(spec.may_mutate, BlockSet::DERIVED_CACHE);
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    // These metadata value types are intentionally not top-level API names. Their
    // public projections still expose exact discriminants/bits for matrix review.
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.derived_effects.recompute.bits(), 1 << 2);
    assert_eq!(spec.derived_effects.preserve.bits(), (1 << 0) | (1 << 4));
    assert_eq!(spec.derived_effects.invalidate.bits(), 0);
    assert_eq!(format!("{:?}", spec.cip_state), "Preserve");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "weak_valence_cache_assignment"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "assign_valence_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| row.feature.name == "valence")
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));
}

#[test]
fn value_operation_commits_valid_cache_and_shares_every_untouched_block() {
    let source = molecule();
    let output = source.with_assigned_valence().unwrap();

    assert!(std::ptr::eq(source.topology(), output.topology()));
    assert!(std::ptr::eq(source.coordinates(), output.coordinates()));
    assert!(std::ptr::eq(source.properties(), output.properties()));
    assert_eq!(source.property("source"), Some("preserved"));
    assert_eq!(output.property("source"), Some("preserved"));
    assert_eq!(
        source.topology().stereo_groups,
        output.topology().stereo_groups
    );
    assert!(format!("{source:?}").contains("derived_cache_is_empty: true"));
    assert!(format!("{output:?}").contains("derived_cache_is_empty: false"));
}

#[test]
fn short_and_parameterized_forms_agree_and_recompute_idempotently() {
    let source = molecule();
    let short = source.with_assigned_valence().unwrap();
    let explicit = source
        .with_assigned_valence_with_params(&ValenceParams::default())
        .unwrap();
    assert_eq!(short, explicit);
    assert!(format!("{short:?}").contains("derived_cache_is_empty: false"));
    assert!(format!("{explicit:?}").contains("derived_cache_is_empty: false"));

    let recomputed = short.with_assigned_valence().unwrap();
    assert_eq!(recomputed, short);
    assert!(std::ptr::eq(short.topology(), recomputed.topology()));
    assert!(format!("{recomputed:?}").contains("derived_cache_is_empty: false"));
}

#[test]
fn inplace_forms_commit_only_after_success() {
    let mut short = molecule();
    short.assign_valence_().unwrap();
    assert!(format!("{short:?}").contains("derived_cache_is_empty: false"));

    let mut explicit = molecule();
    explicit
        .assign_valence_with_params_(&ValenceParams::default())
        .unwrap();
    assert_eq!(explicit, short);
    assert!(format!("{explicit:?}").contains("derived_cache_is_empty: false"));
}

#[test]
fn typed_failure_is_atomic_and_read_only_predicate_keeps_exact_errors() {
    let source = invalid_pentavalent_carbon();
    let mut target = source.clone();
    let observer = target.clone();
    let error = target.assign_valence_().unwrap_err();
    assert!(matches!(
        &error,
        OperationError::Valence(ValenceError::InvalidValence {
            atom,
            atomic_number: 6,
            ..
        }) if *atom == AtomId::new(0)
    ));
    assert!(error.source().is_some());
    assert_eq!(target, source);
    assert!(std::ptr::eq(target.topology(), observer.topology()));
    assert!(std::ptr::eq(target.coordinates(), observer.coordinates()));
    assert!(std::ptr::eq(target.properties(), observer.properties()));
    assert!(format!("{target:?}").contains("derived_cache_is_empty: true"));

    assert_eq!(source.has_valence_violation(AtomId::new(0)), Ok(true));
    assert_eq!(
        source.has_valence_violation(AtomId::new(99)),
        Err(ValenceError::AtomOutOfRange {
            atom: AtomId::new(99),
            atom_count: 6,
        })
    );
}

#[test]
fn single_output_contract_makes_multi_candidate_ordering_not_applicable() {
    let spec = operation_spec("with_assigned_valence_with_params").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
}
