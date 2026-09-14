use std::panic::{AssertUnwindSafe, catch_unwind};

use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, CoordinateBlock,
    Element, MoleculeProperties, SdfPropertyList, SdfPropertyListTarget, TopologyBlock,
    TopologyMapping,
};

use super::*;
use crate::ops::{
    BlockAccess, CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement, MoleculeOpKind,
    OperationDomain, ParityPolicy, SemanticPreconditionSet, SupportStatus, TopologyEditKind,
};

struct FailureAccess;

fn no_effects() -> DerivedEffects {
    DerivedEffects::new(
        DerivedState::NONE,
        DerivedState::NONE,
        DerivedState::NONE,
        DerivedState::NONE,
    )
}

fn spec(
    method: &'static str,
    output: MoleculeOpOutput,
    access: BlockAccess,
    preconditions: SemanticPreconditionSet,
) -> &'static MoleculeOpSpec {
    Box::leak(Box::new(MoleculeOpSpec {
        method,
        impl_fn: "failure_test_impl",
        output,
        result_type: "Molecule",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::None,
        access,
        may_mutate: access.write(),
        auto_remap: BlockSet::NONE,
        derived_effects: no_effects(),
        cip_state: CipStatePolicy::Preserve,
        semantic_preconditions: preconditions,
        requires_mapping: MappingRequirement::None,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    }))
}

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn molecule() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![atom(0, Element::C), atom(1, Element::N)],
        vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )],
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(7, vec![[0.0, 0.0], [1.0, 0.0]])],
        ..Default::default()
    };
    let properties = MoleculeProperties::default()
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
        ));
    Molecule::from_parts(topology, coordinates, properties).unwrap()
}

fn property_write_spec(method: &'static str) -> &'static MoleculeOpSpec {
    spec(
        method,
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::PROPERTIES),
        SemanticPreconditionSet::NONE,
    )
}

fn algorithm_error(method: &'static str) -> OperationError {
    OperationError::Algorithm {
        operation: method,
        detail: "detached body rejected input".to_owned(),
    }
}

fn stage_property_then_fail(
    parts: &mut OpParts<'_, FailureAccess>,
    method: &'static str,
) -> Result<(), OperationError> {
    let properties = parts
        .checkout_properties_runtime()?
        .with_name("staged-not-authoritative");
    parts.install_properties_runtime(properties)?;
    Err(algorithm_error(method))
}

#[test]
fn value_body_error_preserves_the_complete_source_and_exact_error() {
    let source = molecule();
    let before = source.clone();
    let operation = property_write_spec("value-body-error");
    let mut parts = OpParts::<FailureAccess>::new(&source, operation).unwrap();
    let error = stage_property_then_fail(&mut parts, operation.method).unwrap_err();
    assert_eq!(error, algorithm_error("value-body-error"));
    drop(parts);
    assert_eq!(source, before);
    assert_eq!(source.properties().name(), Some("source"));
}

#[test]
fn value_body_panic_propagates_and_preserves_the_complete_source() {
    let source = molecule();
    let before = source.clone();
    let result = catch_unwind(AssertUnwindSafe(|| {
        let mut parts =
            OpParts::<FailureAccess>::new(&source, property_write_spec("value-panic")).unwrap();
        let properties = parts
            .checkout_properties_runtime()
            .unwrap()
            .with_name("staged-before-panic");
        parts.install_properties_runtime(properties).unwrap();
        panic!("body panic remains a panic");
    }));
    assert!(result.is_err());
    assert_eq!(source, before);
}

#[test]
fn in_place_constructor_rejections_leave_the_target_unchanged() {
    let mut target = molecule();
    let before = target.clone();
    let output_error = match OpParts::<FailureAccess>::new_in_place(
        &mut target,
        spec(
            "wrong-output",
            MoleculeOpOutput::Multiple,
            BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
            SemanticPreconditionSet::NONE,
        ),
    ) {
        Ok(_) => panic!("multiple-output spec unexpectedly constructed single-output parts"),
        Err(error) => error,
    };
    assert_eq!(
        output_error,
        OperationError::OutputMismatch {
            operation: "wrong-output",
            expected: MoleculeOpOutput::Single,
            actual: MoleculeOpOutput::Multiple,
        }
    );
    assert_eq!(target, before);

    let precondition_error = match OpParts::<FailureAccess>::new_in_place(
        &mut target,
        spec(
            "missing-precondition",
            MoleculeOpOutput::Single,
            BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
            SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY,
        ),
    ) {
        Ok(_) => panic!("missing precondition unexpectedly constructed parts"),
        Err(error) => error,
    };
    assert_eq!(
        precondition_error,
        OperationError::SemanticPreconditionContract {
            operation: "missing-precondition",
            missing: SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY,
            issue: "the live runtime has no authoritative provenance evidence for this precondition",
        }
    );
    assert_eq!(target, before);
}

#[test]
fn in_place_body_error_before_staging_returns_unchanged_target_and_error() {
    let mut target = molecule();
    let before = target.clone();
    let operation = property_write_spec("early-body-error");
    let returned = {
        let mut parts = OpParts::<FailureAccess>::new_in_place(&mut target, operation).unwrap();
        let error = algorithm_error(operation.method);
        parts.abort_in_place();
        error
    };
    assert_eq!(returned, algorithm_error("early-body-error"));
    assert_eq!(target, before);
}

#[test]
fn abort_discards_one_staged_block_and_returns_original_body_error() {
    let mut target = molecule();
    let before = target.clone();
    let operation = property_write_spec("staged-body-error");
    let returned = {
        let mut parts = OpParts::<FailureAccess>::new_in_place(&mut target, operation).unwrap();
        let error = stage_property_then_fail(&mut parts, operation.method).unwrap_err();
        parts.abort_in_place();
        assert!(matches!(parts.properties, WorkingBlock::Shared));
        error
    };
    assert_eq!(returned, algorithm_error("staged-body-error"));
    assert_eq!(target, before);
}

#[test]
fn abort_clears_all_staged_blocks_and_contract_evidence() {
    let mut target = molecule();
    let before = target.clone();
    {
        let access = BlockAccess::new(
            BlockSet::NONE,
            BlockSet::TOPOLOGY
                .union(BlockSet::COORDINATES)
                .union(BlockSet::PROPERTIES)
                .union(BlockSet::DERIVED_CACHE),
        );
        let mut parts = OpParts::<FailureAccess>::new_in_place(
            &mut target,
            spec(
                "all-staged",
                MoleculeOpOutput::Single,
                access,
                SemanticPreconditionSet::NONE,
            ),
        )
        .unwrap();

        let topology = parts.checkout_topology_runtime().unwrap();
        parts.install_topology_runtime(topology).unwrap();
        let coordinates = parts.checkout_coordinates_runtime().unwrap();
        parts.install_coordinates_runtime(coordinates).unwrap();
        let properties = parts
            .checkout_properties_runtime()
            .unwrap()
            .with_name("all-staged");
        parts.install_properties_runtime(properties).unwrap();
        let mut cache = parts.checkout_derived_cache_runtime().unwrap();
        cache.mark_valid(DerivedState::RINGS);
        parts.install_derived_cache_runtime(cache).unwrap();
        parts.topology_edit = Some(TopologyEditKind::Local);
        parts.topology_mapping = Some(TopologyMapping::identity(2, 1));
        parts.remapped_blocks = BlockSet::COORDINATES.union(BlockSet::PROPERTIES);
        parts.effect_trace = EffectTrace {
            updated: DerivedState::RINGS,
            cleared: DerivedState::VALENCE,
            preserved: DerivedState::STEREO,
            cip_applied: true,
        };

        parts.abort_in_place();
        assert!(matches!(parts.topology, WorkingBlock::Shared));
        assert!(matches!(parts.coordinates, WorkingBlock::Shared));
        assert!(matches!(parts.properties, WorkingBlock::Shared));
        assert!(matches!(parts.derived_cache, WorkingBlock::Shared));
        assert_eq!(parts.topology_edit, None);
        assert_eq!(parts.topology_mapping, None);
        assert!(parts.remapped_blocks.is_empty());
        assert_eq!(parts.effect_trace, EffectTrace::default());
    }
    assert_eq!(target, before);
}

#[test]
fn abort_recovers_a_checked_out_slot_without_exposing_a_placeholder() {
    let mut target = molecule();
    let before = target.clone();
    {
        let mut parts = OpParts::<FailureAccess>::new_in_place(
            &mut target,
            property_write_spec("checked-out-abort"),
        )
        .unwrap();
        let _detached = parts.checkout_properties_runtime().unwrap();
        assert!(matches!(parts.properties, WorkingBlock::CheckedOut));
        parts.abort_in_place();
        assert!(matches!(parts.properties, WorkingBlock::Shared));
    }
    assert_eq!(target, before);
}

#[test]
fn shared_arc_observer_and_unique_target_both_remain_isolated_on_abort() {
    let mut unique_target = molecule();
    let unique_before = unique_target.clone();
    {
        let mut parts = OpParts::<FailureAccess>::new_in_place(
            &mut unique_target,
            property_write_spec("unique-abort"),
        )
        .unwrap();
        let _ = stage_property_then_fail(&mut parts, "unique-abort").unwrap_err();
        parts.abort_in_place();
    }
    assert_eq!(unique_target, unique_before);

    let mut shared_target = molecule();
    let observer = shared_target.clone();
    let shared_before = shared_target.clone();
    {
        let mut parts = OpParts::<FailureAccess>::new_in_place(
            &mut shared_target,
            property_write_spec("shared-abort"),
        )
        .unwrap();
        let _ = stage_property_then_fail(&mut parts, "shared-abort").unwrap_err();
        parts.abort_in_place();
    }
    assert_eq!(shared_target, shared_before);
    assert_eq!(observer, shared_before);
}

#[test]
fn in_place_body_panic_propagates_and_target_remains_complete() {
    let mut target = molecule();
    let before = target.clone();
    let result = catch_unwind(AssertUnwindSafe(|| {
        let mut parts = OpParts::<FailureAccess>::new_in_place(
            &mut target,
            property_write_spec("in-place-panic"),
        )
        .unwrap();
        let properties = parts
            .checkout_properties_runtime()
            .unwrap()
            .with_name("staged-before-panic");
        parts.install_properties_runtime(properties).unwrap();
        panic!("in-place body panic remains a panic");
    }));
    assert!(result.is_err());
    assert_eq!(target, before);
    assert_eq!(target.topology().validate(), Ok(()));
    assert_eq!(
        target
            .coordinate_block_runtime()
            .validate_for_atom_count(target.num_atoms()),
        Ok(())
    );
}

#[test]
fn finish_error_is_atomic_after_staging() {
    let mut target = molecule();
    let before = target.clone();
    let error = {
        let mut parts = OpParts::<FailureAccess>::new_in_place(
            &mut target,
            property_write_spec("finish-error"),
        )
        .unwrap();
        let _detached = parts.checkout_properties_runtime().unwrap();
        parts.finish_in_place().unwrap_err()
    };
    assert_eq!(
        error,
        OperationError::IncompleteCommit {
            operation: "finish-error",
            block: "properties",
        }
    );
    assert_eq!(target, before);
}

#[test]
fn wrapper_and_owner_source_guards_preserve_failure_ordering_and_boundaries() {
    let wrappers = include_str!("../../../cosmolkit-macros/src/wrappers.rs");
    let context = include_str!("../../src/ops/context.rs");
    let molecule_source = include_str!("../../src/molecule.rs");
    let lib = include_str!("../../src/lib.rs");
    let core_manifest = include_str!("../../../cosmolkit-core/Cargo.toml");
    let model_manifest = include_str!("../../../cosmolkit-model/Cargo.toml");

    assert_eq!(wrappers.matches("parts.abort_in_place();").count(), 2);
    assert!(wrappers.contains("Err(error) => {\n                            parts.abort_in_place();\n                            return Err(error);"));
    assert!(wrappers.contains("if let Err(error) = #impl_fn(&mut parts, #(#call_args),*) {\n                        parts.abort_in_place();\n                        return Err(error);"));
    assert!(wrappers.contains("parts.finish_in_place()?;\n                    Ok(result)"));
    assert!(wrappers.contains("parts.finish_in_place()"));
    assert!(!wrappers.contains("catch_unwind"));

    assert_eq!(context.matches("pub(super) fn abort_in_place").count(), 1);
    assert_eq!(context.matches("pub(super) fn finish_in_place").count(), 1);
    assert!(!context.contains("pub(crate) fn abort_in_place"));
    assert!(!context.contains("pub(crate) fn finish_in_place"));
    assert_eq!(
        context
            .matches("\n    in_place_target: Option<&'a mut Molecule>,")
            .count(),
        1
    );
    assert_eq!(
        context
            .matches("\n        in_place_target: Option<&'a mut Molecule>,")
            .count(),
        1
    );
    assert!(!context.contains("rollback_target"));
    assert!(!context.contains("pub struct OpParts"));
    assert!(!lib.contains("pub use ops::OpParts"));
    assert!(!context.contains("cosmolkit_core"));
    assert_eq!(molecule_source.matches("struct MoleculeState").count(), 1);
    assert!(!core_manifest.contains("cosmolkit ="));
    assert!(!model_manifest.contains("cosmolkit ="));
}
