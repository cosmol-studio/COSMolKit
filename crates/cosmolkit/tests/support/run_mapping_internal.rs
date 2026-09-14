use cosmolkit_model::{
    Atom, AtomId, AtomMapping, AtomSpec, Bond, BondId, BondMapping, BondOrder, BondSpec,
    Conformer2D, Conformer3D, CoordinateBlock, Element, MappingValidationError, MoleculeProperties,
    SdfPropertyList, SdfPropertyListTarget, TopologyBlock, TopologyMapping,
};

use super::*;
use crate::ops::{
    BlockAccess, CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement, MoleculeOpKind,
    OperationDomain, ParityPolicy, SemanticPreconditionSet, SupportStatus, TopologyEditKind,
};

struct MappingAccess;

fn spec(
    method: &'static str,
    edit: TopologyEditKind,
    requirement: MappingRequirement,
    access: BlockAccess,
    auto_remap: BlockSet,
) -> &'static MoleculeOpSpec {
    spec_with_cip(
        method,
        edit,
        requirement,
        access,
        auto_remap,
        CipStatePolicy::Preserve,
    )
}

fn spec_with_cip(
    method: &'static str,
    edit: TopologyEditKind,
    requirement: MappingRequirement,
    access: BlockAccess,
    auto_remap: BlockSet,
    cip_state: CipStatePolicy,
) -> &'static MoleculeOpSpec {
    Box::leak(Box::new(MoleculeOpSpec {
        method,
        impl_fn: "mapping_test_impl",
        output: MoleculeOpOutput::Single,
        result_type: "Molecule",
        domain: OperationDomain::Topology,
        kind: if matches!(
            edit,
            TopologyEditKind::Compacting
                | TopologyEditKind::Appending
                | TopologyEditKind::Renumbering
        ) {
            MoleculeOpKind::Strong
        } else {
            MoleculeOpKind::Weak
        },
        topology_edit: edit,
        access,
        may_mutate: access.write(),
        auto_remap,
        derived_effects: DerivedEffects::new(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        cip_state,
        semantic_preconditions: SemanticPreconditionSet::NONE,
        requires_mapping: requirement,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    }))
}

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn topology(elements: &[Element]) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        elements
            .iter()
            .copied()
            .enumerate()
            .map(|(index, element)| atom(index, element))
            .collect(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap()
}

fn bonded_topology() -> TopologyBlock {
    TopologyBlock::try_from_parts(
        vec![atom(0, Element::C), atom(1, Element::N)],
        vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )],
        Vec::new(),
        Vec::new(),
    )
    .unwrap()
}

fn mapping(
    atom_old_to_new: Vec<Option<usize>>,
    atom_new_to_old: Vec<Option<usize>>,
    bond_old_to_new: Vec<Option<usize>>,
    bond_new_to_old: Vec<Option<usize>>,
) -> TopologyMapping {
    TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atom_old_to_new
                .into_iter()
                .map(|row| row.map(AtomId::new))
                .collect(),
            new_to_old: atom_new_to_old
                .into_iter()
                .map(|row| row.map(AtomId::new))
                .collect(),
        },
        bonds: BondMapping {
            old_to_new: bond_old_to_new
                .into_iter()
                .map(|row| row.map(BondId::new))
                .collect(),
            new_to_old: bond_new_to_old
                .into_iter()
                .map(|row| row.map(BondId::new))
                .collect(),
        },
    }
}

fn properties(
    atom_values: Vec<Option<&str>>,
    bond_values: Vec<Option<&str>>,
) -> MoleculeProperties {
    MoleculeProperties::default()
        .with_name("source")
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atoms",
            atom_values
                .into_iter()
                .map(|value| value.map(str::to_owned))
                .collect(),
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bonds",
            bond_values
                .into_iter()
                .map(|value| value.map(str::to_owned))
                .collect(),
        ))
}

fn molecule_with_rows(has_coordinates: bool) -> Molecule {
    let coordinates = if has_coordinates {
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                4,
                vec![[10.0, 11.0], [20.0, 21.0], [30.0, 31.0]],
            )],
            conformers_3d: vec![Conformer3D::new(
                9,
                vec![[10.0, 11.0, 12.0], [20.0, 21.0, 22.0], [30.0, 31.0, 32.0]],
                true,
            )],
            ..Default::default()
        }
    } else {
        CoordinateBlock::default()
    };
    Molecule::from_parts(
        topology(&[Element::C, Element::N, Element::O]),
        coordinates,
        properties(vec![Some("c"), None, Some("o")], Vec::new()),
    )
    .unwrap()
}

fn all_write() -> BlockAccess {
    BlockAccess::new(
        BlockSet::NONE,
        BlockSet::TOPOLOGY
            .union(BlockSet::COORDINATES)
            .union(BlockSet::PROPERTIES),
    )
}

#[test]
fn requirement_profiles_reject_missing_unexpected_and_nonidentity_mappings() {
    let one = topology(&[Element::C]);
    let none = spec(
        "none",
        TopologyEditKind::None,
        MappingRequirement::None,
        BlockAccess::new(BlockSet::TOPOLOGY, BlockSet::NONE),
        BlockSet::NONE,
    );
    assert_eq!(
        OpParts::<MappingAccess>::validate_mapping_obligation(none, &one, &one, None, None),
        Ok(())
    );
    assert!(matches!(
        OpParts::<MappingAccess>::validate_mapping_obligation(
            none,
            &one,
            &one,
            None,
            Some(&TopologyMapping::identity(1, 0)),
        ),
        Err(OperationError::MappingContract {
            operation: "none",
            issue: "mapping was not expected",
            requirement: MappingRequirement::None,
        })
    ));

    let identity = spec(
        "identity",
        TopologyEditKind::Local,
        MappingRequirement::Identity,
        all_write(),
        BlockSet::NONE,
    );
    assert!(matches!(
        OpParts::<MappingAccess>::validate_mapping_obligation(
            identity,
            &one,
            &one,
            Some(TopologyEditKind::Local),
            None,
        ),
        Err(OperationError::MappingContract {
            issue: "required mapping was not recorded",
            ..
        })
    ));
    assert_eq!(
        OpParts::<MappingAccess>::validate_mapping_obligation(
            identity,
            &one,
            &one,
            Some(TopologyEditKind::Local),
            Some(&TopologyMapping::identity(1, 0)),
        ),
        Ok(())
    );
    let nonidentity = mapping(vec![None], vec![None], Vec::new(), Vec::new());
    assert!(matches!(
        OpParts::<MappingAccess>::validate_mapping_obligation(
            identity,
            &one,
            &one,
            Some(TopologyEditKind::Local),
            Some(&nonidentity),
        ),
        Err(OperationError::MappingContract {
            issue: "mapping is not identity",
            ..
        })
    ));
}

#[test]
fn topology_edit_and_mapping_records_are_exact_and_at_most_once() {
    let source = molecule_with_rows(false);
    let compact = spec(
        "edit",
        TopologyEditKind::Compacting,
        MappingRequirement::Required,
        all_write(),
        BlockSet::NONE,
    );
    let mut parts = OpParts::<MappingAccess>::new(&source, compact).unwrap();
    assert!(matches!(
        parts.record_topology_edit_runtime(TopologyEditKind::Appending),
        Err(OperationError::TopologyEditContract {
            issue: "declaration-mismatch",
            expected: TopologyEditKind::Compacting,
            actual: Some(TopologyEditKind::Appending),
            ..
        })
    ));
    parts
        .record_topology_edit_runtime(TopologyEditKind::Compacting)
        .unwrap();
    assert!(matches!(
        parts.record_topology_edit_runtime(TopologyEditKind::Compacting),
        Err(OperationError::TopologyEditContract {
            issue: "duplicate-record",
            ..
        })
    ));
    parts
        .record_topology_mapping_runtime(TopologyMapping::identity(3, 0))
        .unwrap();
    assert!(matches!(
        parts.record_topology_mapping_runtime(TopologyMapping::identity(3, 0)),
        Err(OperationError::MappingContract {
            issue: "duplicate mapping record",
            ..
        })
    ));
}

#[test]
fn atom_mapping_failures_preserve_exact_direction_and_fields() {
    let topology = topology(&[Element::C]);
    let required = spec(
        "atom-errors",
        TopologyEditKind::Local,
        MappingRequirement::Required,
        all_write(),
        BlockSet::NONE,
    );
    let cases = [
        (
            mapping(Vec::new(), vec![Some(0)], Vec::new(), Vec::new()),
            MappingValidationError::Length {
                entity: "atom",
                direction: "old-to-new",
                actual: 0,
                expected: 1,
            },
        ),
        (
            mapping(vec![Some(0)], Vec::new(), Vec::new(), Vec::new()),
            MappingValidationError::Length {
                entity: "atom",
                direction: "new-to-old",
                actual: 0,
                expected: 1,
            },
        ),
        (
            mapping(vec![Some(2)], vec![None], Vec::new(), Vec::new()),
            MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 2,
                target_count: 1,
            },
        ),
        (
            mapping(vec![None], vec![Some(2)], Vec::new(), Vec::new()),
            MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "new-to-old",
                row: 0,
                mapped: 2,
                target_count: 1,
            },
        ),
        (
            mapping(vec![Some(0)], vec![None], Vec::new(), Vec::new()),
            MappingValidationError::InverseMismatch {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 0,
            },
        ),
        (
            mapping(vec![None], vec![Some(0)], Vec::new(), Vec::new()),
            MappingValidationError::InverseMismatch {
                entity: "atom",
                direction: "new-to-old",
                row: 0,
                mapped: 0,
            },
        ),
    ];
    for (candidate, expected) in cases {
        assert_eq!(
            OpParts::<MappingAccess>::validate_mapping_obligation(
                required,
                &topology,
                &topology,
                Some(TopologyEditKind::Local),
                Some(&candidate),
            ),
            Err(OperationError::InvalidTopologyMapping {
                operation: "atom-errors",
                source: expected,
            })
        );
    }
}

#[test]
fn bond_mapping_failures_preserve_exact_direction_and_fields() {
    let topology = bonded_topology();
    let required = spec(
        "bond-errors",
        TopologyEditKind::Local,
        MappingRequirement::Required,
        all_write(),
        BlockSet::NONE,
    );
    let atoms = || (vec![Some(0), Some(1)], vec![Some(0), Some(1)]);
    let cases = [
        {
            let (old, new) = atoms();
            (
                mapping(old, new, Vec::new(), vec![Some(0)]),
                MappingValidationError::Length {
                    entity: "bond",
                    direction: "old-to-new",
                    actual: 0,
                    expected: 1,
                },
            )
        },
        {
            let (old, new) = atoms();
            (
                mapping(old, new, vec![Some(0)], Vec::new()),
                MappingValidationError::Length {
                    entity: "bond",
                    direction: "new-to-old",
                    actual: 0,
                    expected: 1,
                },
            )
        },
        {
            let (old, new) = atoms();
            (
                mapping(old, new, vec![Some(3)], vec![None]),
                MappingValidationError::OutOfRange {
                    entity: "bond",
                    direction: "old-to-new",
                    row: 0,
                    mapped: 3,
                    target_count: 1,
                },
            )
        },
        {
            let (old, new) = atoms();
            (
                mapping(old, new, vec![None], vec![Some(3)]),
                MappingValidationError::OutOfRange {
                    entity: "bond",
                    direction: "new-to-old",
                    row: 0,
                    mapped: 3,
                    target_count: 1,
                },
            )
        },
        {
            let (old, new) = atoms();
            (
                mapping(old, new, vec![Some(0)], vec![None]),
                MappingValidationError::InverseMismatch {
                    entity: "bond",
                    direction: "old-to-new",
                    row: 0,
                    mapped: 0,
                },
            )
        },
        {
            let (old, new) = atoms();
            (
                mapping(old, new, vec![None], vec![Some(0)]),
                MappingValidationError::InverseMismatch {
                    entity: "bond",
                    direction: "new-to-old",
                    row: 0,
                    mapped: 0,
                },
            )
        },
    ];
    for (candidate, expected) in cases {
        assert_eq!(
            OpParts::<MappingAccess>::validate_mapping_obligation(
                required,
                &topology,
                &topology,
                Some(TopologyEditKind::Local),
                Some(&candidate),
            ),
            Err(OperationError::InvalidTopologyMapping {
                operation: "bond-errors",
                source: expected,
            })
        );
    }
}

#[test]
fn malformed_mapping_is_rejected_before_projection_or_any_installation() {
    let source = molecule_with_rows(true);
    let original_coordinates = source.coordinate_block_runtime().clone();
    let original_properties = source.properties().clone();
    let compact = spec(
        "invalid-before-helper",
        TopologyEditKind::Compacting,
        MappingRequirement::Required,
        all_write(),
        BlockSet::COORDINATES.union(BlockSet::PROPERTIES),
    );
    let mut parts = OpParts::<MappingAccess>::new(&source, compact).unwrap();
    let _ = parts.checkout_topology_runtime().unwrap();
    parts
        .install_topology_runtime(topology(&[Element::O, Element::C]))
        .unwrap();
    parts
        .record_topology_edit_runtime(TopologyEditKind::Compacting)
        .unwrap();
    parts
        .record_topology_mapping_runtime(mapping(
            vec![None, None, None],
            vec![None, Some(99)],
            Vec::new(),
            Vec::new(),
        ))
        .unwrap();

    assert!(matches!(
        parts.apply_runtime_remap_runtime(),
        Err(OperationError::InvalidTopologyMapping {
            source: MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "new-to-old",
                row: 1,
                mapped: 99,
                target_count: 3,
            },
            ..
        })
    ));
    assert!(matches!(parts.coordinates, WorkingBlock::Shared));
    assert!(matches!(parts.properties, WorkingBlock::Shared));
    assert_eq!(source.coordinate_block_runtime(), &original_coordinates);
    assert_eq!(source.properties(), &original_properties);
}

#[test]
fn deletion_and_reorder_remap_all_coordinate_and_property_rows_in_new_order() {
    let source = molecule_with_rows(true);
    let compact = spec_with_cip(
        "delete-reorder",
        TopologyEditKind::Compacting,
        MappingRequirement::Required,
        all_write(),
        BlockSet::COORDINATES.union(BlockSet::PROPERTIES),
        CipStatePolicy::ClearComputed,
    );
    let mut parts = OpParts::<MappingAccess>::new(&source, compact).unwrap();
    let _ = parts.checkout_topology_runtime().unwrap();
    parts
        .install_topology_runtime(topology(&[Element::O, Element::C]))
        .unwrap();
    parts
        .record_topology_edit_runtime(TopologyEditKind::Compacting)
        .unwrap();
    parts
        .record_topology_mapping_runtime(mapping(
            vec![Some(1), None, Some(0)],
            vec![Some(2), Some(0)],
            Vec::new(),
            Vec::new(),
        ))
        .unwrap();
    parts.apply_runtime_remap_runtime().unwrap();
    parts.apply_cip_policy_runtime().unwrap();
    let output = parts.finish().unwrap();

    assert_eq!(
        output.coordinate_block_runtime().conformers_2d[0].coordinates(),
        &[[30.0, 31.0], [10.0, 11.0]]
    );
    assert_eq!(
        output.coordinate_block_runtime().conformers_3d[0].coordinates(),
        &[[30.0, 31.0, 32.0], [10.0, 11.0, 12.0]]
    );
    assert_eq!(
        output.properties().sdf_property_lists()[0].values(),
        &[Some("o".to_owned()), Some("c".to_owned())]
    );
    assert_eq!(source.topology().atoms.len(), 3);
    assert_eq!(source.properties().name(), Some("source"));
}

#[test]
fn append_distinguishes_empty_coordinates_owner_values_and_property_none_rows() {
    let append_spec = spec_with_cip(
        "append",
        TopologyEditKind::Appending,
        MappingRequirement::Required,
        all_write(),
        BlockSet::COORDINATES.union(BlockSet::PROPERTIES),
        CipStatePolicy::ClearComputed,
    );
    let candidate = || topology(&[Element::C, Element::N, Element::O, Element::H]);

    let source = molecule_with_rows(true);
    let mut rejected = OpParts::<MappingAccess>::new(&source, append_spec).unwrap();
    let _ = rejected.checkout_topology_runtime().unwrap();
    rejected.install_topology_runtime(candidate()).unwrap();
    rejected
        .record_topology_edit_runtime(TopologyEditKind::Appending)
        .unwrap();
    rejected
        .record_topology_mapping_runtime(TopologyMapping::with_appended(3, 0, 1, 0))
        .unwrap();
    assert_eq!(
        rejected.apply_runtime_remap_runtime(),
        Err(OperationError::CoordinateAppendRequiresValues {
            operation: "append",
        })
    );
    assert!(matches!(rejected.properties, WorkingBlock::Shared));

    let source = molecule_with_rows(false);
    let mut empty = OpParts::<MappingAccess>::new(&source, append_spec).unwrap();
    let _ = empty.checkout_topology_runtime().unwrap();
    empty.install_topology_runtime(candidate()).unwrap();
    empty
        .record_topology_edit_runtime(TopologyEditKind::Appending)
        .unwrap();
    empty
        .record_topology_mapping_runtime(TopologyMapping::with_appended(3, 0, 1, 0))
        .unwrap();
    empty.apply_runtime_remap_runtime().unwrap();
    empty.apply_cip_policy_runtime().unwrap();
    let output = empty.finish().unwrap();
    assert!(output.coordinate_block_runtime().conformers_2d.is_empty());
    assert_eq!(
        output.properties().sdf_property_lists()[0].values(),
        &[Some("c".to_owned()), None, Some("o".to_owned()), None]
    );

    let source = molecule_with_rows(true);
    let mut supplied = OpParts::<MappingAccess>::new(&source, append_spec).unwrap();
    let _ = supplied.checkout_topology_runtime().unwrap();
    supplied.install_topology_runtime(candidate()).unwrap();
    let _ = supplied.checkout_coordinates_runtime().unwrap();
    supplied
        .install_coordinates_runtime(CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                0,
                vec![[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [4.0, 4.0]],
            )],
            ..Default::default()
        })
        .unwrap();
    supplied
        .record_topology_edit_runtime(TopologyEditKind::Appending)
        .unwrap();
    supplied
        .record_topology_mapping_runtime(TopologyMapping::with_appended(3, 0, 1, 0))
        .unwrap();
    supplied.apply_runtime_remap_runtime().unwrap();
    supplied.apply_cip_policy_runtime().unwrap();
    assert_eq!(
        supplied
            .finish()
            .unwrap()
            .coordinate_block_runtime()
            .conformers_2d[0]
            .coordinates()
            .len(),
        4
    );
}

#[test]
fn auto_remap_requires_write_authority_and_is_at_most_once() {
    let source = molecule_with_rows(false);
    let denied_spec = spec(
        "denied",
        TopologyEditKind::Local,
        MappingRequirement::Identity,
        BlockAccess::new(BlockSet::TOPOLOGY, BlockSet::NONE),
        BlockSet::COORDINATES,
    );
    let mut denied = OpParts::<MappingAccess>::new(&source, denied_spec).unwrap();
    denied
        .record_topology_edit_runtime(TopologyEditKind::Local)
        .unwrap();
    denied
        .record_topology_mapping_runtime(TopologyMapping::identity(3, 0))
        .unwrap();
    assert_eq!(
        denied.apply_runtime_remap_runtime(),
        Err(OperationError::AccessDenied {
            operation: "denied",
            block: "coordinates",
        })
    );

    let allowed_spec = spec(
        "once",
        TopologyEditKind::Local,
        MappingRequirement::Identity,
        all_write(),
        BlockSet::COORDINATES,
    );
    let mut allowed = OpParts::<MappingAccess>::new(&source, allowed_spec).unwrap();
    allowed
        .record_topology_edit_runtime(TopologyEditKind::Local)
        .unwrap();
    allowed
        .record_topology_mapping_runtime(TopologyMapping::identity(3, 0))
        .unwrap();
    allowed.apply_runtime_remap_runtime().unwrap();
    assert!(matches!(
        allowed.apply_runtime_remap_runtime(),
        Err(OperationError::AutoRemapContract {
            block: "coordinates",
            issue: "remap was already applied",
            ..
        })
    ));
}

#[test]
fn invalid_candidate_property_rows_fail_before_any_auto_remap_is_retained() {
    let source = molecule_with_rows(false);
    let local = spec(
        "bad-property-candidate",
        TopologyEditKind::Local,
        MappingRequirement::Identity,
        all_write(),
        BlockSet::COORDINATES.union(BlockSet::PROPERTIES),
    );
    let mut parts = OpParts::<MappingAccess>::new(&source, local).unwrap();
    let _ = parts.checkout_properties_runtime().unwrap();
    parts
        .install_properties_runtime(properties(vec![Some("only-one")], Vec::new()))
        .unwrap();
    parts
        .record_topology_edit_runtime(TopologyEditKind::Local)
        .unwrap();
    parts
        .record_topology_mapping_runtime(TopologyMapping::identity(3, 0))
        .unwrap();
    assert!(matches!(
        parts.apply_runtime_remap_runtime(),
        Err(OperationError::InvalidPropertyList {
            target: "atom",
            values: 1,
            expected: 3,
            ..
        })
    ));
    assert!(!parts.remapped_blocks.contains(BlockSet::COORDINATES));
    assert!(!parts.remapped_blocks.contains(BlockSet::PROPERTIES));
}

#[test]
fn source_guards_keep_mapping_runtime_private_and_algorithm_free() {
    let context = include_str!("../../src/ops/context.rs");
    let model_mapping = include_str!("../../../cosmolkit-model/src/mapping.rs");
    let registry = include_str!("../../src/ops/registry.rs");

    let obligation_start = context
        .find("fn validate_mapping_obligation(")
        .expect("mapping-obligation validator must remain present");
    let obligation_end = context[obligation_start..]
        .find("pub(super) fn apply_runtime_remap_runtime")
        .map(|offset| obligation_start + offset)
        .expect("mapping-obligation validator must end before remapping");
    let obligation = &context[obligation_start..obligation_end];
    let required_mapping = obligation
        .find("required mapping was not recorded")
        .expect("mapping obligation must first require the declared mapping");
    let obligation_validation = obligation
        .find(".validate_for_counts(")
        .expect("mapping obligation must validate mapping dimensions");
    let identity_validation = obligation
        .find("mapping != TopologyMapping::identity")
        .expect("identity mappings must receive the additional identity check");
    assert!(required_mapping < obligation_validation);
    assert!(obligation_validation < identity_validation);
    assert_eq!(obligation.match_indices(".validate_for_counts(").count(), 1);

    let leaf_start = context
        .find("PreservationProof::LeafAtomAppend => {")
        .expect("leaf-atom-append preservation proof must remain present");
    let leaf_end = context[leaf_start..]
        .find("PreservationProof::RadicalElectronAssignment => {")
        .map(|offset| leaf_start + offset)
        .expect("leaf-atom-append proof must end at the next proof arm");
    let leaf_proof = &context[leaf_start..leaf_end];
    let leaf_mapping = leaf_proof
        .find("leaf-atom-append proof requires a recorded topology mapping")
        .expect("leaf proof must require a recorded mapping");
    let leaf_validation = leaf_proof
        .find(".validate_for_counts(")
        .expect("leaf proof must validate mapping dimensions");
    let prefix_validation = leaf_proof
        .find("let atom_prefix_is_identity")
        .expect("leaf proof must validate preserved row identity");
    assert!(leaf_mapping < leaf_validation);
    assert!(leaf_validation < prefix_validation);
    assert_eq!(leaf_proof.match_indices(".validate_for_counts(").count(), 1);

    let authoritative_validations = context
        .match_indices(".validate_for_counts(")
        .map(|(position, _)| position)
        .collect::<Vec<_>>();
    assert_eq!(
        authoritative_validations,
        vec![
            obligation_start + obligation_validation,
            leaf_start + leaf_validation,
        ],
        "every mapping-dimension validation must belong to an audited authority boundary"
    );
    assert_eq!(context.matches(".remap_topology(").count(), 2);
    assert_eq!(
        model_mapping.matches("pub fn validate_for_counts(").count(),
        1
    );
    assert_eq!(registry.matches("molecule_ops!").count(), 1);
    for forbidden in [
        "pub topology_mapping:",
        "pub(crate) topology_mapping:",
        "pub fn record_topology_mapping",
        "pub(crate) fn record_topology_mapping_runtime",
        "cosmolkit_core",
        "RDKit",
    ] {
        assert!(
            !context.contains(forbidden),
            "found forbidden source: {forbidden}"
        );
    }
}
