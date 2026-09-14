use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomMapping, AtomSpec, Bond, BondId, BondMapping, BondOrder,
    BondSpec, Conformer2D, CoordinateBlock, CoordinateValidationError, Element,
    MappingValidationError, MoleculeProperties, SdfPropertyList, SdfPropertyListTarget,
    StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
    TopologyBlock, TopologyMapping, TopologyValidationError,
};

use super::*;
use crate::ops::{
    BlockAccess, CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement, MoleculeOpKind,
    OperationDomain, ParityPolicy, SemanticPreconditionSet, SupportStatus, TopologyEditKind,
};

struct CommitAccess;

#[derive(Clone, Copy)]
struct SpecFields {
    output: MoleculeOpOutput,
    access: BlockAccess,
    may_mutate: BlockSet,
    auto_remap: BlockSet,
    effects: DerivedEffects,
    cip: CipStatePolicy,
    preconditions: SemanticPreconditionSet,
    edit: TopologyEditKind,
    mapping: MappingRequirement,
}

fn no_effects() -> DerivedEffects {
    DerivedEffects::new(
        DerivedState::NONE,
        DerivedState::NONE,
        DerivedState::NONE,
        DerivedState::NONE,
    )
}

fn base_fields() -> SpecFields {
    SpecFields {
        output: MoleculeOpOutput::Single,
        access: BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
        may_mutate: BlockSet::NONE,
        auto_remap: BlockSet::NONE,
        effects: no_effects(),
        cip: CipStatePolicy::Preserve,
        preconditions: SemanticPreconditionSet::NONE,
        edit: TopologyEditKind::None,
        mapping: MappingRequirement::None,
    }
}

fn spec(method: &'static str, fields: SpecFields) -> &'static MoleculeOpSpec {
    Box::leak(Box::new(MoleculeOpSpec {
        method,
        impl_fn: "commit_test_impl",
        output: fields.output,
        result_type: "Molecule",
        domain: OperationDomain::Topology,
        kind: if matches!(
            fields.edit,
            TopologyEditKind::Compacting
                | TopologyEditKind::Appending
                | TopologyEditKind::Renumbering
        ) {
            MoleculeOpKind::Strong
        } else {
            MoleculeOpKind::Weak
        },
        topology_edit: fields.edit,
        access: fields.access,
        may_mutate: fields.may_mutate,
        auto_remap: fields.auto_remap,
        derived_effects: fields.effects,
        cip_state: fields.cip,
        semantic_preconditions: fields.preconditions,
        requires_mapping: fields.mapping,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    }))
}

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn topology() -> TopologyBlock {
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

fn properties() -> MoleculeProperties {
    MoleculeProperties::default()
        .with_name("source")
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atom_rows",
            vec![Some("c".to_owned()), Some("n".to_owned())],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bond_rows",
            vec![Some("cn".to_owned())],
        ))
}

fn coordinates() -> CoordinateBlock {
    CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(3, vec![[0.0, 0.0], [1.0, 0.0]])],
        ..Default::default()
    }
}

fn molecule() -> Molecule {
    Molecule::from_parts(topology(), coordinates(), properties()).unwrap()
}

fn all_blocks() -> BlockSet {
    BlockSet::TOPOLOGY
        .union(BlockSet::COORDINATES)
        .union(BlockSet::PROPERTIES)
        .union(BlockSet::DERIVED_CACHE)
}

fn construction_error(result: Result<OpParts<'_, CommitAccess>, OperationError>) -> OperationError {
    match result {
        Ok(_) => panic!("operation construction unexpectedly succeeded"),
        Err(error) => error,
    }
}

#[test]
fn semantic_preconditions_fail_closed_with_exact_missing_mask() {
    let source = molecule();
    assert!(OpParts::<CommitAccess>::new(&source, spec("empty", base_fields())).is_ok());
    for (method, missing) in [
        ("trusted", SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY),
        (
            "hydrogen-ownership",
            SemanticPreconditionSet::HYDROGEN_OWNERSHIP_REPRESENTED,
        ),
        (
            "combined",
            SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY
                .union(SemanticPreconditionSet::HYDROGEN_OWNERSHIP_REPRESENTED),
        ),
    ] {
        let mut fields = base_fields();
        fields.preconditions = missing;
        let operation = spec(method, fields);
        assert_eq!(
            construction_error(OpParts::<CommitAccess>::new(&source, operation)),
            OperationError::SemanticPreconditionContract {
                operation: method,
                missing,
                issue: "the live runtime has no authoritative provenance evidence for this precondition",
            }
        );
        let mut target = source.clone();
        let before = target.clone();
        assert_eq!(
            construction_error(OpParts::<CommitAccess>::new_in_place(
                &mut target,
                operation,
            )),
            OperationError::SemanticPreconditionContract {
                operation: method,
                missing,
                issue: "the live runtime has no authoritative provenance evidence for this precondition",
            }
        );
        assert_eq!(target, before);
    }
}

#[test]
fn output_mismatch_precedes_precondition_and_body_access() {
    let source = molecule();
    let mut fields = base_fields();
    fields.output = MoleculeOpOutput::Multiple;
    fields.preconditions = SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY;
    assert_eq!(
        construction_error(OpParts::<CommitAccess>::new(
            &source,
            spec("multiple", fields),
        )),
        OperationError::OutputMismatch {
            operation: "multiple",
            expected: MoleculeOpOutput::Single,
            actual: MoleculeOpOutput::Multiple,
        }
    );
}

#[test]
fn finish_rejects_each_general_declaration_error_with_exact_masks() {
    let source = molecule();
    let rows = [
        (
            "overlap",
            BlockAccess::new(BlockSet::TOPOLOGY, BlockSet::TOPOLOGY),
            BlockSet::TOPOLOGY,
            BlockSet::NONE,
            "access",
            0,
            BlockSet::TOPOLOGY.bits(),
        ),
        (
            "may-mutate-missing",
            BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY),
            BlockSet::NONE,
            BlockSet::NONE,
            "may_mutate",
            BlockSet::TOPOLOGY.bits(),
            0,
        ),
        (
            "may-mutate-extra",
            BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
            BlockSet::TOPOLOGY,
            BlockSet::NONE,
            "may_mutate",
            0,
            BlockSet::TOPOLOGY.bits(),
        ),
        (
            "unauthorized-remap",
            BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY),
            BlockSet::TOPOLOGY,
            BlockSet::COORDINATES,
            "auto_remap",
            BlockSet::TOPOLOGY.bits(),
            BlockSet::COORDINATES.bits(),
        ),
    ];
    for (method, access, may_mutate, auto_remap, field, expected, actual) in rows {
        let mut fields = base_fields();
        fields.access = access;
        fields.may_mutate = may_mutate;
        fields.auto_remap = auto_remap;
        assert_eq!(
            OpParts::<CommitAccess>::new(&source, spec(method, fields))
                .unwrap()
                .finish(),
            Err(OperationError::OperationContract {
                operation: method,
                field,
                issue: if field == "access" {
                    "read and write sets overlap"
                } else if field == "may_mutate" {
                    "must equal access.write"
                } else {
                    "contains a block without write authority"
                },
                expected,
                actual,
            })
        );
    }

    let mut valid = base_fields();
    valid.access = BlockAccess::new(BlockSet::TOPOLOGY, BlockSet::COORDINATES);
    valid.may_mutate = BlockSet::COORDINATES;
    assert!(
        OpParts::<CommitAccess>::new(&source, spec("valid", valid))
            .unwrap()
            .finish()
            .is_ok()
    );
}

#[test]
fn finish_closes_none_identity_and_required_mapping_profiles() {
    let source = molecule();
    assert!(
        OpParts::<CommitAccess>::new(&source, spec("none", base_fields()))
            .unwrap()
            .finish()
            .is_ok()
    );

    let mut identity_fields = base_fields();
    identity_fields.access = BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY);
    identity_fields.may_mutate = BlockSet::TOPOLOGY;
    identity_fields.edit = TopologyEditKind::Local;
    identity_fields.mapping = MappingRequirement::Identity;
    let identity_spec = spec("identity", identity_fields);
    assert!(matches!(
        OpParts::<CommitAccess>::new(&source, identity_spec)
            .unwrap()
            .finish(),
        Err(OperationError::TopologyEditContract {
            issue: "missing-record",
            ..
        })
    ));
    let mut identity = OpParts::<CommitAccess>::new(&source, identity_spec).unwrap();
    let detached = identity.checkout_topology_runtime().unwrap();
    identity.install_topology_runtime(detached).unwrap();
    identity
        .record_topology_edit_runtime(TopologyEditKind::Local)
        .unwrap();
    identity
        .record_topology_mapping_runtime(TopologyMapping::identity(2, 1))
        .unwrap();
    identity.apply_cip_policy_runtime().unwrap();
    assert!(identity.finish().is_ok());

    let mut required_fields = identity_fields;
    required_fields.edit = TopologyEditKind::Compacting;
    required_fields.mapping = MappingRequirement::Required;
    let required_spec = spec("required", required_fields);
    let mut required = OpParts::<CommitAccess>::new(&source, required_spec).unwrap();
    required.topology_edit = Some(TopologyEditKind::Compacting);
    assert!(matches!(
        required.finish(),
        Err(OperationError::MappingContract {
            issue: "required mapping was not recorded",
            ..
        })
    ));

    let mut malformed = OpParts::<CommitAccess>::new(&source, required_spec).unwrap();
    malformed.topology_edit = Some(TopologyEditKind::Compacting);
    malformed.topology_mapping = Some(TopologyMapping {
        atoms: AtomMapping {
            old_to_new: vec![Some(AtomId::new(9)), Some(AtomId::new(1))],
            new_to_old: vec![Some(AtomId::new(0)), Some(AtomId::new(1))],
        },
        bonds: BondMapping {
            old_to_new: vec![Some(BondId::new(0))],
            new_to_old: vec![Some(BondId::new(0))],
        },
    });
    assert!(matches!(
        malformed.finish(),
        Err(OperationError::InvalidTopologyMapping {
            operation: "required",
            source: MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 9,
                target_count: 2,
            },
        })
    ));
}

#[test]
fn finish_requires_exact_remap_completion_and_accepts_one_validated_application() {
    let source = molecule();
    let mut fields = base_fields();
    fields.access = BlockAccess::new(
        BlockSet::NONE,
        BlockSet::TOPOLOGY
            .union(BlockSet::COORDINATES)
            .union(BlockSet::PROPERTIES),
    );
    fields.may_mutate = fields.access.write();
    fields.auto_remap = BlockSet::COORDINATES.union(BlockSet::PROPERTIES);
    fields.edit = TopologyEditKind::Local;
    fields.mapping = MappingRequirement::Identity;
    let operation = spec("remap", fields);

    for completed in [BlockSet::NONE, BlockSet::COORDINATES, BlockSet::PROPERTIES] {
        let mut parts = OpParts::<CommitAccess>::new(&source, operation).unwrap();
        parts.topology_edit = Some(TopologyEditKind::Local);
        parts.topology_mapping = Some(TopologyMapping::identity(2, 1));
        parts.remapped_blocks = completed;
        assert_eq!(
            parts.finish(),
            Err(OperationError::OperationContract {
                operation: "remap",
                field: "auto_remap",
                issue: "declared remap was not completed",
                expected: fields.auto_remap.bits(),
                actual: completed.bits(),
            })
        );
    }

    let mut complete = OpParts::<CommitAccess>::new(&source, operation).unwrap();
    let detached = complete.checkout_topology_runtime().unwrap();
    complete.install_topology_runtime(detached).unwrap();
    complete
        .record_topology_edit_runtime(TopologyEditKind::Local)
        .unwrap();
    complete
        .record_topology_mapping_runtime(TopologyMapping::identity(2, 1))
        .unwrap();
    complete.apply_runtime_remap_runtime().unwrap();
    complete.apply_cip_policy_runtime().unwrap();
    let output = complete.finish().unwrap();
    assert_eq!(output.topology(), source.topology());
    assert_eq!(
        output.coordinate_block_runtime().conformers_2d[0].coordinates(),
        source.coordinate_block_runtime().conformers_2d[0].coordinates()
    );
    assert_eq!(output.coordinate_block_runtime().conformers_2d[0].id(), 0);
    assert_eq!(output.properties(), source.properties());
}

#[test]
fn finish_rejects_each_missing_effect_category_and_cip_then_accepts_full_trace() {
    let source = molecule();
    let categories = [
        (
            "missing-recompute",
            DerivedEffects::new(
                DerivedState::VALENCE,
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::NONE,
            ),
            "recompute",
            DerivedState::VALENCE,
        ),
        (
            "missing-preserve",
            DerivedEffects::new(
                DerivedState::NONE,
                DerivedState::RINGS,
                DerivedState::NONE,
                DerivedState::NONE,
            ),
            "preserve",
            DerivedState::RINGS,
        ),
        (
            "missing-invalidate",
            DerivedEffects::new(
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::STEREO,
                DerivedState::NONE,
            ),
            "invalidate",
            DerivedState::STEREO,
        ),
        (
            "without_hydrogens",
            DerivedEffects::new(
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::VALENCE,
            ),
            "operation_defined",
            DerivedState::VALENCE,
        ),
    ];
    for (method, effects, action, states) in categories {
        let mut fields = base_fields();
        fields.access = BlockAccess::new(
            BlockSet::NONE,
            BlockSet::PROPERTIES.union(BlockSet::DERIVED_CACHE),
        );
        fields.may_mutate = fields.access.write();
        fields.effects = effects;
        let mut parts = OpParts::<CommitAccess>::new(&source, spec(method, fields)).unwrap();
        parts.properties = WorkingBlock::Installed(source.properties().clone());
        assert!(matches!(
            parts.finish(),
            Err(OperationError::DerivedEffectContract {
                action: found_action,
                states: found_states,
                issue: "declared effect was not completed",
                ..
            }) if found_action == action && found_states == states
        ));
    }

    let mut cip_fields = base_fields();
    cip_fields.access = BlockAccess::new(BlockSet::NONE, BlockSet::PROPERTIES);
    cip_fields.may_mutate = BlockSet::PROPERTIES;
    let mut missing_cip =
        OpParts::<CommitAccess>::new(&source, spec("missing-cip", cip_fields)).unwrap();
    missing_cip.properties = WorkingBlock::Installed(source.properties().clone());
    assert!(matches!(
        missing_cip.finish(),
        Err(OperationError::CipStateContract {
            issue: "CIP policy was not applied",
            ..
        })
    ));

    let mut full_fields = base_fields();
    full_fields.access = BlockAccess::new(BlockSet::NONE, all_blocks());
    full_fields.may_mutate = all_blocks();
    full_fields.effects = DerivedEffects::new(
        DerivedState::VALENCE,
        DerivedState::RINGS,
        DerivedState::STEREO,
        DerivedState::NONE,
    );
    let mut complete = OpParts::<CommitAccess>::new(&source, spec("full", full_fields)).unwrap();
    let mut cache = complete.checkout_derived_cache_runtime().unwrap();
    cache.install_valence_assignment(
        cosmolkit_core::assign_valence(source.topology(), &Default::default()).unwrap(),
    );
    complete.install_derived_cache_runtime(cache).unwrap();
    complete
        .mark_cache_updated_runtime(DerivedState::VALENCE)
        .unwrap();
    complete
        .prove_preserved_runtime(DerivedState::RINGS, PreservationProof::UnchangedInput)
        .unwrap();
    complete.clear_cache_runtime(DerivedState::STEREO).unwrap();
    complete.apply_cip_policy_runtime().unwrap();
    assert!(complete.finish().is_ok());
}

#[test]
fn every_checked_out_block_is_rejected_before_materialization() {
    let source = molecule();
    let mut fields = base_fields();
    fields.access = BlockAccess::new(BlockSet::NONE, all_blocks());
    fields.may_mutate = all_blocks();
    let operation = spec("checked-out", fields);
    for block in ["topology", "coordinates", "properties", "derived_cache"] {
        let mut parts = OpParts::<CommitAccess>::new(&source, operation).unwrap();
        match block {
            "topology" => parts.topology = WorkingBlock::CheckedOut,
            "coordinates" => parts.coordinates = WorkingBlock::CheckedOut,
            "properties" => parts.properties = WorkingBlock::CheckedOut,
            "derived_cache" => parts.derived_cache = WorkingBlock::CheckedOut,
            _ => unreachable!(),
        }
        assert_eq!(
            parts.finish(),
            Err(OperationError::IncompleteCommit {
                operation: "checked-out",
                block,
            })
        );
    }
}

fn invariant_parts(source: &Molecule, method: &'static str) -> OpParts<'static, CommitAccess> {
    let leaked = Box::leak(Box::new(source.clone()));
    let mut fields = base_fields();
    fields.access = BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY);
    fields.may_mutate = BlockSet::TOPOLOGY;
    fields.edit = TopologyEditKind::Local;
    fields.mapping = MappingRequirement::Identity;
    let mut parts = OpParts::<CommitAccess>::new(leaked, spec(method, fields)).unwrap();
    parts.topology_edit = Some(TopologyEditKind::Local);
    parts.topology_mapping = Some(TopologyMapping::identity(2, 1));
    parts.effect_trace.cip_applied = true;
    parts
}

#[test]
fn finish_rejects_topology_adjacency_stereo_and_sgroup_invariants() {
    let source = molecule();

    let mut wrong_id = invariant_parts(&source, "atom-id");
    let mut candidate = source.topology().clone();
    candidate.atoms[0] = candidate.atoms[0].clone().with_id(AtomId::new(7));
    wrong_id.topology = WorkingBlock::Installed(candidate);
    assert!(matches!(
        wrong_id.finish(),
        Err(OperationError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch { .. }
        ))
    ));

    let mut adjacency = invariant_parts(&source, "adjacency");
    let mut candidate = source.topology().clone();
    candidate.adjacency = AdjacencyList::from_topology(2, &[]);
    adjacency.topology = WorkingBlock::Installed(candidate);
    assert!(matches!(
        adjacency.finish(),
        Err(OperationError::InvalidTopology(
            TopologyValidationError::AdjacencyMismatch
        ))
    ));

    let mut stereo = invariant_parts(&source, "stereo");
    let mut candidate = source.topology().clone();
    candidate.stereo_groups.push(StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomId::new(9)],
        Vec::new(),
    ));
    stereo.topology = WorkingBlock::Installed(candidate);
    assert!(matches!(
        stereo.finish(),
        Err(OperationError::InvalidTopology(
            TopologyValidationError::StereoGroupAtomOutOfRange { .. }
        ))
    ));

    let mut sgroup = invariant_parts(&source, "sgroup");
    let mut candidate = source.topology().clone();
    let mut group = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data);
    group.push_atom(AtomId::new(9));
    candidate.substance_groups.push(group);
    sgroup.topology = WorkingBlock::Installed(candidate);
    assert!(matches!(
        sgroup.finish(),
        Err(OperationError::InvalidTopology(
            TopologyValidationError::SubstanceGroupAtomOutOfRange { .. }
        ))
    ));
}

#[test]
fn finish_rejects_coordinate_and_typed_property_row_invariants() {
    let source = molecule();
    let mut bad_coordinates = invariant_parts(&source, "coordinates");
    bad_coordinates.topology = WorkingBlock::Shared;
    bad_coordinates.coordinates = WorkingBlock::Installed(CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(3, vec![[0.0, 0.0]])],
        ..Default::default()
    });
    assert!(matches!(
        bad_coordinates.finish(),
        Err(OperationError::InvalidCoordinates(
            CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: 3,
                rows: 1,
                atom_count: 2,
            }
        ))
    ));

    for (target, name, values, expected) in [
        (SdfPropertyListTarget::Atom, "bad_atoms", 1, 2),
        (SdfPropertyListTarget::Bond, "bad_bonds", 0, 1),
    ] {
        let mut bad_properties = invariant_parts(&source, name);
        bad_properties.topology = WorkingBlock::Shared;
        bad_properties.properties = WorkingBlock::Installed(
            MoleculeProperties::default().with_sdf_property_list(SdfPropertyList::new(
                target,
                name,
                vec![None; values],
            )),
        );
        assert_eq!(
            bad_properties.finish(),
            Err(OperationError::InvalidPropertyList {
                target: if target == SdfPropertyListTarget::Atom {
                    "atom"
                } else {
                    "bond"
                },
                name: name.to_owned(),
                values,
                expected,
            })
        );
    }
}

#[test]
fn untouched_and_staged_value_success_construct_once_and_preserve_source() {
    let source = molecule();
    let before = source.clone();
    let untouched = OpParts::<CommitAccess>::new(&source, spec("untouched", base_fields()))
        .unwrap()
        .finish()
        .unwrap();
    assert_eq!(untouched.topology(), source.topology());
    assert_eq!(
        untouched.coordinate_block_runtime(),
        source.coordinate_block_runtime()
    );
    assert_eq!(untouched.properties(), source.properties());
    assert_eq!(untouched.runtime_constructions(), 1);
    assert_eq!(source, before);

    let mut fields = base_fields();
    fields.access = BlockAccess::new(BlockSet::NONE, BlockSet::PROPERTIES);
    fields.may_mutate = BlockSet::PROPERTIES;
    let mut staged = OpParts::<CommitAccess>::new(&source, spec("staged", fields)).unwrap();
    let properties = staged
        .checkout_properties_runtime()
        .unwrap()
        .with_name("result");
    staged.install_properties_runtime(properties).unwrap();
    staged.apply_cip_policy_runtime().unwrap();
    let output = staged.finish().unwrap();
    assert_eq!(output.properties().name(), Some("result"));
    assert_eq!(output.runtime_constructions(), 1);
    assert_eq!(source.properties().name(), Some("source"));
    assert_eq!(source, before);
}

#[test]
fn in_place_install_is_atomic_for_failure_and_success() {
    let mut fields = base_fields();
    fields.access = BlockAccess::new(BlockSet::NONE, BlockSet::PROPERTIES);
    fields.may_mutate = BlockSet::PROPERTIES;
    let operation = spec("in-place", fields);

    let mut failed = molecule();
    let before = failed.clone();
    let mut transaction = OpParts::<CommitAccess>::new_in_place(&mut failed, operation).unwrap();
    transaction.properties = WorkingBlock::CheckedOut;
    assert_eq!(
        transaction.finish_in_place(),
        Err(OperationError::IncompleteCommit {
            operation: "in-place",
            block: "properties",
        })
    );
    assert_eq!(failed, before);

    let mut target = molecule();
    let mut transaction = OpParts::<CommitAccess>::new_in_place(&mut target, operation).unwrap();
    let properties = transaction
        .checkout_properties_runtime()
        .unwrap()
        .with_name("installed");
    transaction.install_properties_runtime(properties).unwrap();
    transaction.apply_cip_policy_runtime().unwrap();
    transaction.finish_in_place().unwrap();
    assert_eq!(target.properties().name(), Some("installed"));
    assert_eq!(target.runtime_constructions(), 1);
}

#[test]
fn source_guards_keep_one_private_non_domain_commit_owner() {
    let context = include_str!("../../src/ops/context.rs");
    let molecule = include_str!("../../src/molecule.rs");
    let lib = include_str!("../../src/lib.rs");
    assert_eq!(context.matches("fn validate_candidate(&self)").count(), 1);
    assert_eq!(context.matches("pub(super) fn finish(self)").count(), 1);
    assert!(!context.contains("pub(crate) fn finish(self)"));
    assert_eq!(context.matches("*target = replacement;").count(), 1);
    assert!(context.contains("self.validate_candidate()?;"));
    assert!(!context.contains("pub struct OpParts"));
    assert!(!lib.contains("pub use ops::OpParts"));
    assert!(!context.contains("cosmolkit_core"));
    assert!(!context.contains("MOLECULE_OPS"));
    assert_eq!(molecule.matches("struct MoleculeState").count(), 1);
}
