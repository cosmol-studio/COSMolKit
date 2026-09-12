use std::sync::Arc;

use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, CoordinateBlock, Element,
    MoleculeProperties, TopologyBlock,
};

use super::*;
use crate::molecule::DerivedCacheBlock;
use crate::ops::{
    BlockAccess, CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement, MoleculeOpKind,
    OperationDomain, ParityPolicy, SemanticPreconditionSet, SupportStatus, TopologyEditKind,
};

struct EffectsAccess;

fn all_effect_access() -> BlockAccess {
    BlockAccess::new(
        BlockSet::NONE,
        BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE),
    )
}

fn spec(
    method: &'static str,
    output: MoleculeOpOutput,
    access: BlockAccess,
    may_mutate: BlockSet,
    effects: DerivedEffects,
    cip_state: CipStatePolicy,
) -> &'static MoleculeOpSpec {
    Box::leak(Box::new(MoleculeOpSpec {
        method,
        impl_fn: "effects_test_impl",
        output,
        result_type: "Molecule",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::None,
        access,
        may_mutate,
        auto_remap: BlockSet::NONE,
        derived_effects: effects,
        cip_state,
        semantic_preconditions: SemanticPreconditionSet::NONE,
        requires_mapping: MappingRequirement::None,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    }))
}

fn effects(
    recompute: DerivedState,
    preserve: DerivedState,
    invalidate: DerivedState,
    operation_defined: DerivedState,
) -> DerivedEffects {
    DerivedEffects::new(recompute, preserve, invalidate, operation_defined)
}

fn atom(index: usize) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C))
}

fn molecule() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![atom(0), atom(1)],
        vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )],
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

fn molecule_with_valid(states: DerivedState) -> Molecule {
    let source = molecule();
    let mut cache = DerivedCacheBlock::default();
    cache.mark_valid(states);
    Molecule::from_runtime_parts(
        Arc::new(source.topology().clone()),
        Arc::new(source.coordinates().clone()),
        Arc::new(source.properties().clone()),
        Arc::new(cache),
    )
    .unwrap()
}

fn effect_error(result: Result<(), OperationError>) -> OperationError {
    result.expect_err("effect action unexpectedly succeeded")
}

#[test]
fn derived_state_set_algebra_preserves_exact_masks() {
    let a = DerivedState::RINGS.union(DerivedState::VALENCE);
    let b = DerivedState::VALENCE.union(DerivedState::STEREO);
    assert_eq!(a.intersection(b), DerivedState::VALENCE);
    assert_eq!(a.difference(b), DerivedState::RINGS);
    assert!(a.intersects(b));
    assert!(DerivedState::NONE.is_empty());
    assert_eq!(
        a.bits(),
        DerivedState::RINGS.bits() | DerivedState::VALENCE.bits()
    );
}

#[test]
fn empty_and_each_independent_category_validate() {
    for declared in [
        effects(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        effects(
            DerivedState::RINGS,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        effects(
            DerivedState::NONE,
            DerivedState::RINGS,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        effects(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::RINGS,
            DerivedState::NONE,
        ),
    ] {
        let operation = spec(
            "independent",
            MoleculeOpOutput::Single,
            all_effect_access(),
            all_effect_access().write(),
            declared,
            CipStatePolicy::Preserve,
        );
        assert_eq!(
            OpParts::<EffectsAccess>::validate_effect_contract(operation),
            Ok(())
        );
    }
}

#[test]
fn every_pairwise_overlap_reports_the_exact_state() {
    let ring = DerivedState::RINGS;
    let declarations = [
        effects(ring, ring, DerivedState::NONE, DerivedState::NONE),
        effects(ring, DerivedState::NONE, ring, DerivedState::NONE),
        effects(ring, DerivedState::NONE, DerivedState::NONE, ring),
        effects(DerivedState::NONE, ring, ring, DerivedState::NONE),
        effects(DerivedState::NONE, ring, DerivedState::NONE, ring),
        effects(DerivedState::NONE, DerivedState::NONE, ring, ring),
        effects(ring, ring, ring, ring),
    ];
    for declared in declarations {
        let operation = spec(
            "overlap",
            MoleculeOpOutput::Single,
            all_effect_access(),
            all_effect_access().write(),
            declared,
            CipStatePolicy::Preserve,
        );
        assert_eq!(
            OpParts::<EffectsAccess>::validate_effect_contract(operation),
            Err(OperationError::DerivedEffectContract {
                operation: "overlap",
                action: "declaration",
                states: ring,
                issue: "effect categories overlap",
            })
        );
    }
}

#[test]
fn effects_require_cache_write_and_may_mutate_authority() {
    let declared = effects(
        DerivedState::RINGS,
        DerivedState::NONE,
        DerivedState::NONE,
        DerivedState::NONE,
    );
    for (access, may_mutate) in [
        (
            BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
            BlockSet::NONE,
        ),
        (
            BlockAccess::new(BlockSet::NONE, BlockSet::DERIVED_CACHE),
            BlockSet::NONE,
        ),
    ] {
        let operation = spec(
            "authority",
            MoleculeOpOutput::Single,
            access,
            may_mutate,
            declared,
            CipStatePolicy::Preserve,
        );
        assert!(matches!(
            OpParts::<EffectsAccess>::validate_effect_contract(operation),
            Err(OperationError::DerivedEffectContract {
                operation: "authority",
                action: "declaration",
                states: DerivedState::RINGS,
                issue: "derived effects require derived_cache write and may_mutate authority",
            })
        ));
    }
}

#[test]
fn operation_defined_allow_list_is_exact() {
    for method in ["without_hydrogens", "without_hydrogens_with_params"] {
        let operation = spec(
            method,
            MoleculeOpOutput::Single,
            all_effect_access(),
            all_effect_access().write(),
            effects(
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::VALENCE,
            ),
            CipStatePolicy::Preserve,
        );
        assert_eq!(
            OpParts::<EffectsAccess>::validate_effect_contract(operation),
            Ok(())
        );
    }
    for (method, states) in [
        ("other", DerivedState::VALENCE),
        ("without_hydrogens", DerivedState::RINGS),
        (
            "without_hydrogens",
            DerivedState::VALENCE.union(DerivedState::RINGS),
        ),
    ] {
        let operation = spec(
            method,
            MoleculeOpOutput::Single,
            all_effect_access(),
            all_effect_access().write(),
            effects(
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::NONE,
                states,
            ),
            CipStatePolicy::Preserve,
        );
        assert!(matches!(
            OpParts::<EffectsAccess>::validate_effect_contract(operation),
            Err(OperationError::DerivedEffectContract {
                action: "operation_defined",
                issue: "only valence in the hydrogen-removal family is allow-listed",
                ..
            })
        ));
    }
}

#[test]
fn update_and_clear_change_only_declared_cache_bits() {
    let source = molecule_with_valid(DerivedState::RINGS.union(DerivedState::DRAWING));
    let operation = spec(
        "cache",
        MoleculeOpOutput::Single,
        all_effect_access(),
        all_effect_access().write(),
        effects(
            DerivedState::VALENCE,
            DerivedState::NONE,
            DerivedState::RINGS,
            DerivedState::NONE,
        ),
        CipStatePolicy::Preserve,
    );
    let mut parts = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    parts
        .mark_cache_updated_runtime(DerivedState::VALENCE)
        .unwrap();
    parts.clear_cache_runtime(DerivedState::RINGS).unwrap();
    parts.apply_cip_policy_runtime().unwrap();
    parts.validate_effect_completion().unwrap();
    let output = parts.finish().unwrap();
    let valid = output.derived_cache_runtime().valid_states();
    assert!(valid.contains(DerivedState::VALENCE));
    assert!(valid.contains(DerivedState::DRAWING));
    assert!(!valid.intersects(DerivedState::RINGS));
    assert!(
        source
            .derived_cache_runtime()
            .valid_states()
            .contains(DerivedState::RINGS)
    );
}

#[test]
fn empty_wrong_category_and_duplicate_actions_are_atomic() {
    let source = molecule_with_valid(DerivedState::DRAWING);
    let operation = spec(
        "bad-actions",
        MoleculeOpOutput::Single,
        all_effect_access(),
        all_effect_access().write(),
        effects(
            DerivedState::VALENCE,
            DerivedState::NONE,
            DerivedState::RINGS,
            DerivedState::NONE,
        ),
        CipStatePolicy::Preserve,
    );
    let mut parts = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    for error in [
        effect_error(parts.mark_cache_updated_runtime(DerivedState::NONE)),
        effect_error(parts.mark_cache_updated_runtime(DerivedState::RINGS)),
        effect_error(parts.clear_cache_runtime(DerivedState::DRAWING)),
    ] {
        assert!(matches!(
            error,
            OperationError::DerivedEffectContract { .. }
        ));
    }
    assert_eq!(parts.effect_trace, EffectTrace::default());
    assert_eq!(
        parts.current_cache_candidate().unwrap().valid_states(),
        DerivedState::DRAWING
    );
    parts
        .mark_cache_updated_runtime(DerivedState::VALENCE)
        .unwrap();
    let before = parts.current_cache_candidate().unwrap();
    assert!(matches!(
        parts.mark_cache_updated_runtime(DerivedState::VALENCE),
        Err(OperationError::DerivedEffectContract {
            issue: "state was already handled",
            ..
        })
    ));
    assert_eq!(parts.current_cache_candidate().unwrap(), before);
}

#[test]
fn preservation_requires_an_objective_supported_proof() {
    let source = molecule_with_valid(DerivedState::RINGS);
    let operation = spec(
        "preserve",
        MoleculeOpOutput::Single,
        all_effect_access(),
        all_effect_access().write(),
        effects(
            DerivedState::NONE,
            DerivedState::RINGS,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        CipStatePolicy::Preserve,
    );
    let mut parts = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    assert!(matches!(
        parts.prove_preserved_runtime(DerivedState::RINGS, PreservationProof::LeafAtomAppend),
        Err(OperationError::DerivedEffectContract {
            issue: "leaf-atom-append proof is not implemented by this runtime unit",
            ..
        })
    ));
    assert_eq!(parts.effect_trace, EffectTrace::default());
    parts
        .prove_preserved_runtime(DerivedState::RINGS, PreservationProof::UnchangedInput)
        .unwrap();

    let mut changed = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    let properties = changed.checkout_properties_runtime().unwrap();
    changed
        .install_properties_runtime(properties.with_name("changed"))
        .unwrap();
    assert!(matches!(
        changed.prove_preserved_runtime(DerivedState::RINGS, PreservationProof::UnchangedInput),
        Err(OperationError::DerivedEffectContract {
            issue: "unchanged-input proof failed",
            ..
        })
    ));
}

#[test]
fn completion_reports_each_missing_category_then_accepts_all() {
    let operation = spec(
        "complete",
        MoleculeOpOutput::Single,
        all_effect_access(),
        all_effect_access().write(),
        effects(
            DerivedState::VALENCE,
            DerivedState::RING_FAMILIES,
            DerivedState::RINGS,
            DerivedState::NONE,
        ),
        CipStatePolicy::Preserve,
    );
    let source = molecule();
    let mut parts = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    let properties = parts.checkout_properties_runtime().unwrap();
    parts.install_properties_runtime(properties).unwrap();
    assert!(matches!(
        parts.validate_effect_completion(),
        Err(OperationError::DerivedEffectContract {
            action: "recompute",
            states: DerivedState::VALENCE,
            ..
        })
    ));
    parts
        .mark_cache_updated_runtime(DerivedState::VALENCE)
        .unwrap();
    assert!(matches!(
        parts.validate_effect_completion(),
        Err(OperationError::DerivedEffectContract {
            action: "preserve",
            states: DerivedState::RING_FAMILIES,
            ..
        })
    ));
    parts
        .prove_preserved_runtime(
            DerivedState::RING_FAMILIES,
            PreservationProof::UnchangedInput,
        )
        .unwrap();
    assert!(matches!(
        parts.validate_effect_completion(),
        Err(OperationError::DerivedEffectContract {
            action: "invalidate",
            states: DerivedState::RINGS,
            ..
        })
    ));
    parts.clear_cache_runtime(DerivedState::RINGS).unwrap();
    assert!(matches!(
        parts.validate_effect_completion(),
        Err(OperationError::CipStateContract {
            issue: "CIP policy was not applied",
            ..
        })
    ));
    parts.apply_cip_policy_runtime().unwrap();
    assert_eq!(parts.validate_effect_completion(), Ok(()));
}

fn cip_source() -> Molecule {
    let mut source = molecule();
    let mut topology = source.topology().clone();
    topology.atoms[0].set_prop("_CIPCode", "R").unwrap();
    topology.atoms[0]
        .set_computed_prop("_CIPNeighborOrder", "1")
        .unwrap();
    topology.bonds[0]
        .set_computed_prop("_CIPNeighborOrder", "0")
        .unwrap();
    let mut properties = source.properties().clone();
    properties
        .set_computed_prop("_CIPComputed", "true")
        .unwrap();
    properties.set_prop("ordinary", "kept").unwrap();
    source = Molecule::from_parts(topology, CoordinateBlock::default(), properties).unwrap();
    source
}

#[test]
fn cip_preserve_accepts_equal_state_and_rejects_changed_observable_state() {
    let source = cip_source();
    let operation = spec(
        "cip-preserve",
        MoleculeOpOutput::Single,
        all_effect_access(),
        all_effect_access().write(),
        effects(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        CipStatePolicy::Preserve,
    );
    let mut unchanged = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    unchanged.apply_cip_policy_runtime().unwrap();

    let mut changed = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    let mut topology = changed.checkout_topology_runtime().unwrap();
    topology.atoms[0].set_prop("_CIPCode", "S").unwrap();
    changed.install_topology_runtime(topology).unwrap();
    assert_eq!(
        changed.apply_cip_policy_runtime(),
        Err(OperationError::CipStateContract {
            operation: "cip-preserve",
            policy: CipStatePolicy::Preserve,
            issue: "candidate changed CIP-observable state",
        })
    );
    assert!(!changed.effect_trace.cip_applied);
}

#[test]
fn cip_clear_removes_computed_members_only_and_is_atomic() {
    let source = cip_source();
    let operation = spec(
        "cip-clear",
        MoleculeOpOutput::Single,
        all_effect_access(),
        all_effect_access().write(),
        effects(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        CipStatePolicy::ClearComputed,
    );
    let mut parts = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    parts.apply_cip_policy_runtime().unwrap();
    let output = parts.finish().unwrap();
    assert_eq!(output.topology().atoms[0].prop("_CIPCode"), Some("R"));
    assert_eq!(output.topology().atoms[0].prop("_CIPNeighborOrder"), None);
    assert_eq!(output.topology().bonds[0].prop("_CIPNeighborOrder"), None);
    assert_eq!(output.properties().prop("_CIPComputed"), None);
    assert_eq!(output.properties().prop("ordinary"), Some("kept"));
    assert_eq!(
        source.topology().atoms[0].prop("_CIPNeighborOrder"),
        Some("1")
    );

    let denied = spec(
        "cip-clear-denied",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY),
        BlockSet::TOPOLOGY,
        effects(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        CipStatePolicy::ClearComputed,
    );
    assert!(matches!(
        OpParts::<EffectsAccess>::validate_effect_contract(denied),
        Err(OperationError::CipStateContract { .. })
    ));
}

#[test]
fn cip_assign_requires_exact_operation_authority_and_assignment_evidence() {
    for method in ["wrong", "with_cip_labels_with_options"] {
        let access = if method == "wrong" {
            all_effect_access()
        } else {
            BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY)
        };
        let operation = spec(
            method,
            MoleculeOpOutput::Single,
            access,
            access.write(),
            effects(
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::NONE,
            ),
            CipStatePolicy::Assign,
        );
        assert!(matches!(
            OpParts::<EffectsAccess>::validate_effect_contract(operation),
            Err(OperationError::CipStateContract { .. })
        ));
    }

    let operation = spec(
        "with_cip_labels_with_options",
        MoleculeOpOutput::Single,
        all_effect_access(),
        all_effect_access().write(),
        effects(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        CipStatePolicy::Assign,
    );
    let source = molecule();
    let mut missing = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    assert!(matches!(
        missing.apply_cip_policy_runtime(),
        Err(OperationError::CipStateContract {
            issue: "assignment did not install computed _CIPComputed evidence",
            ..
        })
    ));
    let mut assigned = OpParts::<EffectsAccess>::new(&source, operation).unwrap();
    let mut properties = assigned.checkout_properties_runtime().unwrap();
    properties
        .set_computed_prop("_CIPComputed", "true")
        .unwrap();
    assigned.install_properties_runtime(properties).unwrap();
    assigned.apply_cip_policy_runtime().unwrap();
}

#[test]
fn tautomer_transition_guard_is_exact_and_single_output_cannot_apply_it() {
    let valid = spec(
        "enumerate_tautomers_with_options",
        MoleculeOpOutput::Multiple,
        all_effect_access(),
        all_effect_access().write(),
        effects(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        CipStatePolicy::TautomerSourceTransition,
    );
    assert_eq!(
        OpParts::<EffectsAccess>::validate_effect_contract(valid),
        Ok(())
    );
    for (method, output, access) in [
        ("wrong", MoleculeOpOutput::Multiple, all_effect_access()),
        (
            "enumerate_tautomers_with_options",
            MoleculeOpOutput::Single,
            all_effect_access(),
        ),
        (
            "enumerate_tautomers_with_options",
            MoleculeOpOutput::Multiple,
            BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY),
        ),
    ] {
        let invalid = spec(
            method,
            output,
            access,
            access.write(),
            effects(
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::NONE,
                DerivedState::NONE,
            ),
            CipStatePolicy::TautomerSourceTransition,
        );
        assert!(matches!(
            OpParts::<EffectsAccess>::validate_effect_contract(invalid),
            Err(OperationError::CipStateContract { .. })
        ));
    }
}

#[test]
fn source_shape_keeps_effect_authority_private_and_domain_free() {
    let context = include_str!("../../src/ops/context.rs");
    let molecule = include_str!("../../src/molecule.rs");
    let manifest = include_str!("../../Cargo.toml");
    assert!(context.contains("struct EffectTrace"));
    assert!(context.contains("fn validate_effect_contract"));
    assert!(context.contains("fn clear_cache_runtime"));
    assert!(!context.contains("pub struct EffectTrace"));
    assert!(!context.contains("pub fn clear_cache"));
    assert!(!molecule.contains("pub struct DerivedCacheBlock"));
    assert!(!manifest.contains("cosmolkit ="));
    for forbidden in ["aromatize", "kekulize", "sanitize", "assign_valence"] {
        assert!(!context.contains(forbidden));
    }
}
