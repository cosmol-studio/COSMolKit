use std::sync::Arc;

use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, CoordinateBlock,
    Element, MoleculeProperties, SdfPropertyList, SdfPropertyListTarget, TopologyBlock,
};

use super::*;
use crate::molecule::DerivedCacheBlock;
use crate::ops::{
    BlockAccess, CipStatePolicy, DerivedEffects, DerivedState, FunctionStatus, MappingRequirement,
    MoleculeOpKind, OperationDomain, ParityPolicy, SemanticPreconditionSet, TopologyEditKind,
};

struct MultipleAccess;

#[derive(Clone, Copy)]
struct SpecFields {
    output: MoleculeOpOutput,
    access: BlockAccess,
    effects: DerivedEffects,
    cip: CipStatePolicy,
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

fn all_blocks() -> BlockSet {
    BlockSet::TOPOLOGY
        .union(BlockSet::COORDINATES)
        .union(BlockSet::PROPERTIES)
        .union(BlockSet::DERIVED_CACHE)
}

fn base_fields() -> SpecFields {
    SpecFields {
        output: MoleculeOpOutput::Multiple,
        access: BlockAccess::new(all_blocks(), BlockSet::NONE),
        effects: no_effects(),
        cip: CipStatePolicy::Preserve,
        edit: TopologyEditKind::None,
        mapping: MappingRequirement::None,
    }
}

fn spec(method: &'static str, fields: SpecFields) -> &'static MoleculeOpSpec {
    Box::leak(Box::new(MoleculeOpSpec {
        method,
        impl_fn: "multiple_test_impl",
        output: fields.output,
        result_type: "Vec<Molecule>",
        domain: OperationDomain::Topology,
        kind: if fields.edit == TopologyEditKind::None {
            MoleculeOpKind::Weak
        } else {
            MoleculeOpKind::Strong
        },
        topology_edit: fields.edit,
        access: fields.access,
        may_mutate: fields.access.write(),
        auto_remap: BlockSet::NONE,
        derived_effects: fields.effects,
        cip_state: fields.cip,
        semantic_preconditions: SemanticPreconditionSet::NONE,
        requires_mapping: fields.mapping,
        status: FunctionStatus::Experimental,
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

fn coordinates() -> CoordinateBlock {
    CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(9, vec![[0.0, 0.0], [1.0, 0.0]])],
        ..Default::default()
    }
}

fn properties(name: &str) -> MoleculeProperties {
    MoleculeProperties::default()
        .with_name(name)
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atom_rows",
            vec![
                Some(cosmolkit_model::PropertyValue::from("c")),
                Some(cosmolkit_model::PropertyValue::from("n")),
            ],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bond_rows",
            vec![Some(cosmolkit_model::PropertyValue::from("cn"))],
        ))
}

fn molecule() -> Molecule {
    Molecule::from_parts(topology(), coordinates(), properties("source")).unwrap()
}

fn tuple(source: &Molecule) -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
    (
        source.topology().clone(),
        source.coordinate_block_runtime().clone(),
        source.properties().clone(),
    )
}

fn construction_error(
    result: Result<MultiOutputOpParts<'_, MultipleAccess>, OperationError>,
) -> OperationError {
    match result {
        Ok(_) => panic!("multiple-output construction unexpectedly succeeded"),
        Err(error) => error,
    }
}

#[test]
fn constructor_rejects_single_output_before_preconditions_with_exact_fields() {
    let source = molecule();
    let mut fields = base_fields();
    fields.output = MoleculeOpOutput::Single;
    assert_eq!(
        construction_error(MultiOutputOpParts::<MultipleAccess>::new(
            &source,
            spec("single", fields),
        )),
        OperationError::OutputMismatch {
            operation: "single",
            expected: MoleculeOpOutput::Multiple,
            actual: MoleculeOpOutput::Single,
        }
    );
}

#[test]
fn source_reads_enforce_each_declared_block_independently() {
    let source = molecule();
    for (allowed, denied_name) in [
        (BlockSet::TOPOLOGY, "coordinates"),
        (BlockSet::COORDINATES, "topology"),
        (BlockSet::PROPERTIES, "topology"),
    ] {
        let mut fields = base_fields();
        fields.access = BlockAccess::new(allowed, BlockSet::NONE);
        let parts =
            MultiOutputOpParts::<MultipleAccess>::new(&source, spec(denied_name, fields)).unwrap();
        let (allowed_result, denied_result) = if allowed == BlockSet::TOPOLOGY {
            (
                parts.source_topology_runtime().map(|_| ()),
                parts.source_coordinates_runtime().map(|_| ()),
            )
        } else if allowed == BlockSet::COORDINATES {
            (
                parts.source_coordinates_runtime().map(|_| ()),
                parts.source_topology_runtime().map(|_| ()),
            )
        } else {
            (
                parts.source_properties_runtime().map(|_| ()),
                parts.source_topology_runtime().map(|_| ()),
            )
        };
        assert_eq!(allowed_result, Ok(()));
        assert_eq!(
            denied_result,
            Err(OperationError::AccessDenied {
                operation: denied_name,
                block: denied_name,
            })
        );
    }
}

#[test]
fn missing_empty_and_duplicate_emit_are_distinct() {
    let source = molecule();
    let operation = spec("emit-state", base_fields());
    assert_eq!(
        MultiOutputOpParts::<MultipleAccess>::new(&source, operation)
            .unwrap()
            .finish(),
        Err(OperationError::IncompleteCommit {
            operation: "emit-state",
            block: "outputs",
        })
    );

    let mut empty = MultiOutputOpParts::<MultipleAccess>::new(&source, operation).unwrap();
    empty.emit_all_runtime(Vec::new()).unwrap();
    assert!(empty.finish().unwrap().is_empty());

    let mut duplicate = MultiOutputOpParts::<MultipleAccess>::new(&source, operation).unwrap();
    duplicate.emit_all_runtime(vec![tuple(&source)]).unwrap();
    assert_eq!(
        duplicate.emit_all_runtime(Vec::new()),
        Err(OperationError::OperationContract {
            operation: "emit-state",
            field: "outputs",
            issue: "multiple-output operation emitted more than once",
            expected: 0,
            actual: 1,
        })
    );
}

#[test]
fn one_many_order_and_duplicate_candidates_are_preserved() {
    let source = molecule();
    let mut fields = base_fields();
    fields.access = BlockAccess::new(
        BlockSet::TOPOLOGY.union(BlockSet::COORDINATES),
        BlockSet::PROPERTIES,
    );
    let operation = spec("ordered", fields);

    let mut one = MultiOutputOpParts::<MultipleAccess>::new(&source, operation).unwrap();
    one.emit_all_runtime(vec![(
        source.topology().clone(),
        source.coordinate_block_runtime().clone(),
        properties("one"),
    )])
    .unwrap();
    let one = one.finish().unwrap();
    assert_eq!(one.len(), 1);
    assert_eq!(
        one[0].properties().name().map(|value| value.as_bytes()),
        Some(b"one".as_slice())
    );

    let candidates = ["first", "duplicate", "duplicate", "last"]
        .into_iter()
        .map(|name| {
            (
                source.topology().clone(),
                source.coordinate_block_runtime().clone(),
                properties(name),
            )
        })
        .collect();
    let mut many = MultiOutputOpParts::<MultipleAccess>::new(&source, operation).unwrap();
    many.emit_all_runtime(candidates).unwrap();
    let names = many
        .finish()
        .unwrap()
        .into_iter()
        .map(|candidate| {
            candidate
                .properties()
                .name()
                .map(|value| value.as_bytes())
                .unwrap()
                .to_owned()
        })
        .collect::<Vec<_>>();
    assert_eq!(
        names,
        [
            b"first".as_slice(),
            b"duplicate".as_slice(),
            b"duplicate".as_slice(),
            b"last".as_slice()
        ]
    );
}

#[test]
fn invalid_topology_coordinate_and_property_rows_remain_structured() {
    let source = molecule();
    let operation = spec("invalid-candidate", base_fields());

    let mut bad_topology = source.topology().clone();
    bad_topology.bonds[0] = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(9), BondOrder::Single),
    );
    let cases = [
        (
            bad_topology,
            source.coordinate_block_runtime().clone(),
            source.properties().clone(),
            "topology",
        ),
        (
            source.topology().clone(),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(1, vec![[0.0, 0.0]])],
                ..Default::default()
            },
            source.properties().clone(),
            "coordinates",
        ),
        (
            source.topology().clone(),
            source.coordinate_block_runtime().clone(),
            MoleculeProperties::default().with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "bad_bonds",
                Vec::new(),
            )),
            "properties",
        ),
    ];
    for (topology, coordinates, properties, kind) in cases {
        let mut parts = MultiOutputOpParts::<MultipleAccess>::new(&source, operation).unwrap();
        parts
            .emit_all_runtime(vec![(topology, coordinates, properties)])
            .unwrap();
        let error = parts.finish().unwrap_err();
        assert!(matches!(
            (kind, error),
            ("topology", OperationError::InvalidTopology(_))
                | ("coordinates", OperationError::InvalidCoordinates(_))
                | ("properties", OperationError::InvalidPropertyList { .. })
        ));
    }
}

#[test]
fn nonidentity_rows_fail_closed_without_mapping_payload() {
    let source = molecule();
    let mut fields = base_fields();
    fields.access = BlockAccess::new(
        BlockSet::COORDINATES.union(BlockSet::PROPERTIES),
        BlockSet::TOPOLOGY.union(BlockSet::DERIVED_CACHE),
    );
    fields.edit = TopologyEditKind::Renumbering;
    fields.mapping = MappingRequirement::Required;
    fields.effects = DerivedEffects::new(
        DerivedState::NONE,
        DerivedState::NONE,
        DerivedState::RINGS,
        DerivedState::NONE,
    );
    let operation = spec("nonidentity", fields);
    let mut reordered = source.topology().clone();
    reordered.atoms.swap(0, 1);
    let mut parts = MultiOutputOpParts::<MultipleAccess>::new(&source, operation).unwrap();
    parts
        .emit_all_runtime(vec![(
            reordered,
            source.coordinate_block_runtime().clone(),
            source.properties().clone(),
        )])
        .unwrap();
    assert_eq!(
        parts.finish(),
        Err(OperationError::MappingContract {
            operation: "nonidentity",
            issue: "multiple-output tuple has no non-identity mapping evidence",
            requirement: MappingRequirement::Required,
        })
    );
}

#[test]
fn changed_block_without_write_authority_is_rejected() {
    let source = molecule();
    let operation = spec("property-denied", base_fields());
    let mut parts = MultiOutputOpParts::<MultipleAccess>::new(&source, operation).unwrap();
    parts
        .emit_all_runtime(vec![(
            source.topology().clone(),
            source.coordinate_block_runtime().clone(),
            properties("changed"),
        )])
        .unwrap();
    assert_eq!(
        parts.finish(),
        Err(OperationError::AccessDenied {
            operation: "property-denied",
            block: "properties",
        })
    );
}

#[test]
#[cfg(feature = "cap-rings")]
fn invalidate_clears_each_candidate_cache_without_touching_source() {
    let base = molecule();
    let mut cache = DerivedCacheBlock::default();
    cache.install_ring_info(cosmolkit_core::RingInfo::new(
        cosmolkit_core::RingFindType::SymmSssr,
        base.topology().atoms.len(),
        base.topology().bonds.len(),
    ));
    cache.mark_valid(DerivedState::RINGS.union(DerivedState::STEREO));
    let source = Molecule::from_runtime_parts(
        Arc::new(base.topology().clone()),
        Arc::new(base.coordinate_block_runtime().clone()),
        Arc::new(base.properties().clone()),
        Arc::new(cache),
    )
    .unwrap();
    let before = source.clone();
    let mut fields = base_fields();
    fields.access = BlockAccess::new(
        BlockSet::TOPOLOGY
            .union(BlockSet::COORDINATES)
            .union(BlockSet::PROPERTIES),
        BlockSet::DERIVED_CACHE,
    );
    fields.effects = DerivedEffects::new(
        DerivedState::NONE,
        DerivedState::NONE,
        DerivedState::RINGS,
        DerivedState::NONE,
    );
    let mut parts =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("invalidate", fields)).unwrap();
    parts
        .emit_all_runtime(vec![tuple(&source), tuple(&source)])
        .unwrap();
    let outputs = parts.finish().unwrap();
    assert_eq!(outputs.len(), 2);
    for output in outputs {
        assert!(
            !output
                .derived_cache_runtime()
                .valid_states()
                .contains(DerivedState::RINGS)
        );
        assert!(
            output
                .derived_cache_runtime()
                .valid_states()
                .contains(DerivedState::STEREO)
        );
    }
    assert_eq!(source, before);
    assert!(
        source
            .derived_cache_runtime()
            .valid_states()
            .contains(DerivedState::RINGS)
    );
}

#[test]
fn preserve_requires_unchanged_input_proof_per_candidate() {
    let source = molecule();
    let mut fields = base_fields();
    fields.access = BlockAccess::new(
        BlockSet::TOPOLOGY.union(BlockSet::COORDINATES),
        BlockSet::PROPERTIES.union(BlockSet::DERIVED_CACHE),
    );
    fields.effects = DerivedEffects::new(
        DerivedState::NONE,
        DerivedState::RINGS,
        DerivedState::NONE,
        DerivedState::NONE,
    );
    let mut parts =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("preserve", fields)).unwrap();
    parts
        .emit_all_runtime(vec![(
            source.topology().clone(),
            source.coordinate_block_runtime().clone(),
            properties("changed"),
        )])
        .unwrap();
    assert!(matches!(
        parts.finish(),
        Err(OperationError::DerivedEffectContract {
            operation: "preserve",
            action: "preserve",
            states,
            issue: "unchanged-input proof failed",
        }) if states == DerivedState::RINGS
    ));
}

#[test]
fn cip_clear_and_tautomer_transition_are_applied_per_candidate() {
    let source = molecule();
    let mut computed = source.properties().clone();
    computed.set_computed_prop("_CIPComputed", "true").unwrap();
    let mut clear_fields = base_fields();
    clear_fields.access = BlockAccess::new(
        BlockSet::COORDINATES,
        BlockSet::TOPOLOGY.union(BlockSet::PROPERTIES),
    );
    clear_fields.cip = CipStatePolicy::ClearComputed;
    let mut clear =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("clear-cip", clear_fields))
            .unwrap();
    clear
        .emit_all_runtime(vec![(
            source.topology().clone(),
            source.coordinate_block_runtime().clone(),
            computed,
        )])
        .unwrap();
    let output = clear.finish().unwrap().pop().unwrap();
    assert_eq!(output.properties().prop("_CIPComputed"), None);

    let mut tautomer_fields = base_fields();
    tautomer_fields.access = BlockAccess::new(
        BlockSet::COORDINATES,
        BlockSet::TOPOLOGY.union(BlockSet::PROPERTIES),
    );
    tautomer_fields.cip = CipStatePolicy::TautomerSourceTransition;
    let mut tautomer = MultiOutputOpParts::<MultipleAccess>::new(
        &source,
        spec("enumerate_tautomers_with_params", tautomer_fields),
    )
    .unwrap();
    tautomer.emit_all_runtime(vec![tuple(&source)]).unwrap();
    assert_eq!(tautomer.finish().unwrap().len(), 1);
}

#[test]
fn invalid_candidate_at_each_position_rejects_all_and_preserves_source() {
    let source = molecule();
    let before = source.clone();
    let operation = spec("atomic", base_fields());
    for bad_index in 0..3 {
        let mut candidates = vec![tuple(&source), tuple(&source), tuple(&source)];
        candidates[bad_index].1 = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(3, vec![[0.0, 0.0]])],
            ..Default::default()
        };
        let mut parts = MultiOutputOpParts::<MultipleAccess>::new(&source, operation).unwrap();
        parts.emit_all_runtime(candidates).unwrap();
        assert!(matches!(
            parts.finish(),
            Err(OperationError::InvalidCoordinates(_))
        ));
        assert_eq!(source, before);
    }
}

#[test]
fn source_shape_has_one_private_owner_and_no_branch_runtime() {
    let multiple = include_str!("../../src/ops/multiple.rs");
    let context = include_str!("../../src/ops/context.rs");
    let ops = include_str!("../../src/ops/mod.rs");
    let runtime = include_str!("../../src/ops/runtime/mod.rs");
    let lib = include_str!("../../src/lib.rs");
    let manifest = include_str!("../../Cargo.toml");
    assert_eq!(multiple.matches("struct MultiOutputOpParts").count(), 1);
    assert!(multiple.contains("pub(crate) struct MultiOutputOpParts"));
    assert!(!multiple.contains("pub struct MultiOutputOpParts"));
    assert!(ops.contains("pub(crate) use runtime::multiple::MultiOutputOpParts"));
    assert!(ops.contains("mod runtime;"));
    assert!(!ops.contains("pub(crate) mod runtime;"));
    assert!(runtime.contains("pub(super) mod multiple;"));
    let runtime_reexport = lib
        .lines()
        .find(|line| {
            line.starts_with("pub(crate) use ops::{") && line.contains("MultiOutputOpParts")
        })
        .expect("lib.rs must keep one crate-private runtime re-export");
    assert!(runtime_reexport.contains("MultiOutputOpParts"));
    assert!(runtime_reexport.contains("OpParts"));
    assert!(!lib.lines().any(|line| line.starts_with("pub use ops::{")
        && (line.contains("MultiOutputOpParts") || line.contains("OpParts"))));
    assert!(context.contains("fn validate_multiple_candidate"));
    for forbidden in [
        "MultiMoleculeOpParts",
        "BranchHandle",
        "branch_id",
        "sanitize(",
        "kekulize(",
        "enumerate_tautomers(",
    ] {
        assert!(!multiple.contains(forbidden));
    }
    assert!(!manifest.contains("cosmolkit ="));
}

fn lazy_fields() -> SpecFields {
    SpecFields {
        output: MoleculeOpOutput::LazyMultiple,
        ..base_fields()
    }
}

fn lazy_tuple(source: &Molecule) -> (TopologyBlock, Option<CoordinateBlock>, MoleculeProperties) {
    (source.topology().clone(), None, source.properties().clone())
}

#[test]
fn lazy_missing_empty_duplicate_and_cardinality_mismatch_are_distinct() {
    let source = molecule();
    let operation = spec("lazy-state", lazy_fields());
    let missing = MultiOutputOpParts::<MultipleAccess>::new(&source, operation)
        .unwrap()
        .finish_lazy();
    assert!(matches!(
        missing,
        Err(OperationError::IncompleteCommit {
            operation: "lazy-state",
            block: "outputs"
        })
    ));
    let mut parts = MultiOutputOpParts::<MultipleAccess>::new(&source, operation).unwrap();
    assert_eq!(
        parts.emit_all_runtime(Vec::new()),
        Err(OperationError::OutputMismatch {
            operation: "lazy-state",
            expected: MoleculeOpOutput::Multiple,
            actual: MoleculeOpOutput::LazyMultiple
        })
    );
    parts.emit_lazy_runtime(std::iter::empty()).unwrap();
    assert!(matches!(
        parts.emit_lazy_runtime(std::iter::empty()),
        Err(OperationError::OperationContract {
            field: "outputs",
            ..
        })
    ));
    let mut empty = parts.finish_lazy().unwrap();
    assert!(empty.next().is_none());
    assert!(empty.next().is_none());
    assert_eq!(empty.yielded_count(), 0);
    let mut eager =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("eager-state", base_fields()))
            .unwrap();
    assert_eq!(
        eager.emit_lazy_runtime(std::iter::empty()),
        Err(OperationError::OutputMismatch {
            operation: "eager-state",
            expected: MoleculeOpOutput::LazyMultiple,
            actual: MoleculeOpOutput::Multiple
        })
    );
    assert!(matches!(
        eager.finish_lazy(),
        Err(OperationError::OutputMismatch {
            expected: MoleculeOpOutput::LazyMultiple,
            actual: MoleculeOpOutput::Multiple,
            ..
        })
    ));
}

#[test]
fn lazy_prefix_pulls_exactly_requested_candidates_and_keeps_unchanged_blocks_shared() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let source = molecule();
    let pulls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&pulls);
    let detached = lazy_tuple(&source);
    // An unbounded stream cannot be eagerly collected. Consume a small prefix.
    let stream = std::iter::from_fn(move || {
        observed.fetch_add(1, Ordering::SeqCst);
        Some(Ok(detached.clone()))
    });
    let mut parts =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("lazy-prefix", lazy_fields()))
            .unwrap();
    parts.emit_lazy_runtime(stream).unwrap();
    assert_eq!(pulls.load(Ordering::SeqCst), 0);
    let mut outputs = parts.finish_lazy().unwrap();
    assert_eq!(pulls.load(Ordering::SeqCst), 0);
    for expected in 1..=7 {
        let output = outputs.next().unwrap().unwrap();
        assert_eq!(pulls.load(Ordering::SeqCst), expected);
        assert_eq!(outputs.yielded_count(), expected);
        assert_eq!(output, source);
        assert!(Arc::ptr_eq(
            &output.topology_arc_runtime(),
            &source.topology_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &output.coordinates_arc_runtime(),
            &source.coordinates_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &output.properties_arc_runtime(),
            &source.properties_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &output.derived_cache_arc_runtime(),
            &source.derived_cache_arc_runtime()
        ));
    }
    drop(outputs);
    assert_eq!(pulls.load(Ordering::SeqCst), 7);
}

#[test]
fn lazy_snapshot_outlives_input_and_remains_independent_of_later_replacement() {
    let mut source = molecule();
    let before = source.clone();
    let expected_coordinates = Arc::clone(&source.coordinates_arc_runtime());
    let detached = lazy_tuple(&source);
    let mut parts =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("lazy-lifetime", lazy_fields()))
            .unwrap();
    parts
        .emit_lazy_runtime(std::iter::once(Ok(detached)))
        .unwrap();
    let mut outputs = parts.finish_lazy().unwrap();
    source = Molecule::from_parts(
        TopologyBlock::default(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap();
    assert_ne!(source, before);
    drop(source);
    let output = outputs.next().unwrap().unwrap();
    assert_eq!(output, before);
    assert!(Arc::ptr_eq(
        &output.coordinates_arc_runtime(),
        &expected_coordinates
    ));
    assert_eq!(outputs.yielded_count(), 1);
    assert!(outputs.next().is_none());
    assert!(outputs.next().is_none());
}

#[test]
fn lazy_callback_error_is_deferred_retains_success_count_and_fuses() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let source = molecule();
    let before = source.clone();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let candidate = lazy_tuple(&source);
    let error = OperationError::AccessDenied {
        operation: "callback",
        block: "user-source",
    };
    let stream = std::iter::from_fn(move || match observed.fetch_add(1, Ordering::SeqCst) {
        0 => Some(Ok(candidate.clone())),
        _ => Some(Err(error.clone())),
    });
    let mut parts =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("lazy-callback", lazy_fields()))
            .unwrap();
    parts.emit_lazy_runtime(stream).unwrap();
    let mut outputs = parts.finish_lazy().unwrap();
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    assert_eq!(outputs.next().unwrap().unwrap(), source);
    assert_eq!(outputs.yielded_count(), 1);
    assert_eq!(
        outputs.next(),
        Some(Err(OperationError::AccessDenied {
            operation: "callback",
            block: "user-source"
        }))
    );
    assert_eq!(outputs.yielded_count(), 1);
    for _ in 0..3 {
        assert!(outputs.next().is_none());
    }
    assert_eq!(calls.load(Ordering::SeqCst), 2);
    assert_eq!(source, before);
}

#[test]
fn lazy_invalid_candidate_at_each_position_rejects_before_counting_and_fuses() {
    let source = molecule();
    let before = source.clone();
    for bad_index in 0..3 {
        let mut candidates = vec![Ok(lazy_tuple(&source)); 3];
        candidates[bad_index].as_mut().unwrap().1 = Some(CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(3, vec![[0.0, 0.0]])],
            ..Default::default()
        });
        let mut parts =
            MultiOutputOpParts::<MultipleAccess>::new(&source, spec("lazy-invalid", lazy_fields()))
                .unwrap();
        parts.emit_lazy_runtime(candidates.into_iter()).unwrap();
        let mut outputs = parts.finish_lazy().unwrap();
        for _ in 0..bad_index {
            outputs.next().unwrap().unwrap();
        }
        assert!(matches!(
            outputs.next(),
            Some(Err(OperationError::InvalidCoordinates(_)))
        ));
        assert_eq!(outputs.yielded_count(), bad_index);
        assert!(outputs.next().is_none());
        assert!(outputs.next().is_none());
        assert_eq!(source, before);
    }
}

#[test]
fn lazy_declared_write_preserves_order_duplicates_and_shares_other_blocks() {
    let source = molecule();
    let mut fields = lazy_fields();
    fields.access = BlockAccess::new(
        BlockSet::TOPOLOGY.union(BlockSet::COORDINATES),
        BlockSet::PROPERTIES,
    );
    let topology = source.topology().clone();
    let candidates = ["first", "duplicate", "duplicate", "last"]
        .into_iter()
        .map(move |name| Ok((topology.clone(), None, properties(name))));
    let mut parts =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("lazy-order", fields)).unwrap();
    parts.emit_lazy_runtime(candidates).unwrap();
    let mut outputs = parts.finish_lazy().unwrap();
    for name in ["first", "duplicate", "duplicate", "last"] {
        let output = outputs.next().unwrap().unwrap();
        assert_eq!(
            output.properties().name().map(|value| value.as_bytes()),
            Some(name.as_bytes())
        );
        assert!(!Arc::ptr_eq(
            &output.properties_arc_runtime(),
            &source.properties_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &output.topology_arc_runtime(),
            &source.topology_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &output.coordinates_arc_runtime(),
            &source.coordinates_arc_runtime()
        ));
    }
    assert!(outputs.next().is_none());
    assert_eq!(outputs.yielded_count(), 4);
}

#[test]
fn lazy_access_mapping_and_preservation_contracts_fail_closed_per_candidate() {
    let source = molecule();
    let changed = (source.topology().clone(), None, properties("changed"));
    let mut fields = lazy_fields();
    let mut parts =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("lazy-access", fields)).unwrap();
    parts
        .emit_lazy_runtime(std::iter::once(Ok(changed.clone())))
        .unwrap();
    let mut outputs = parts.finish_lazy().unwrap();
    assert_eq!(
        outputs.next(),
        Some(Err(OperationError::AccessDenied {
            operation: "lazy-access",
            block: "properties"
        }))
    );
    assert_eq!(outputs.yielded_count(), 0);
    assert!(outputs.next().is_none());

    fields.access = BlockAccess::new(
        BlockSet::COORDINATES,
        BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE),
    );
    fields.edit = TopologyEditKind::Renumbering;
    fields.mapping = MappingRequirement::Required;
    let mut topology = source.topology().clone();
    topology.atoms.swap(0, 1);
    let mut parts =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("lazy-mapping", fields)).unwrap();
    parts
        .emit_lazy_runtime(std::iter::once(Ok((
            topology,
            None,
            source.properties().clone(),
        ))))
        .unwrap();
    let mut outputs = parts.finish_lazy().unwrap();
    assert!(matches!(
        outputs.next(),
        Some(Err(OperationError::MappingContract {
            operation: "lazy-mapping",
            ..
        }))
    ));
    assert_eq!(outputs.yielded_count(), 0);
    assert!(outputs.next().is_none());

    fields.edit = TopologyEditKind::None;
    fields.mapping = MappingRequirement::None;
    fields.effects = DerivedEffects::new(
        DerivedState::NONE,
        DerivedState::RINGS,
        DerivedState::NONE,
        DerivedState::NONE,
    );
    let mut parts =
        MultiOutputOpParts::<MultipleAccess>::new(&source, spec("lazy-preserve", fields)).unwrap();
    parts
        .emit_lazy_runtime(std::iter::once(Ok(changed)))
        .unwrap();
    let mut outputs = parts.finish_lazy().unwrap();
    assert!(matches!(
        outputs.next(),
        Some(Err(OperationError::DerivedEffectContract {
            operation: "lazy-preserve",
            action: "preserve",
            ..
        }))
    ));
    assert_eq!(outputs.yielded_count(), 0);
    assert!(outputs.next().is_none());
}

#[test]
fn lazy_stream_releases_captured_state_on_error_end_or_iterator_drop() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct DropProbe(Arc<AtomicUsize>);
    impl Drop for DropProbe {
        fn drop(&mut self) {
            self.0.fetch_add(1, Ordering::SeqCst);
        }
    }
    let source = molecule();
    for termination in ["error", "end", "drop"] {
        let drops = Arc::new(AtomicUsize::new(0));
        let probe = DropProbe(Arc::clone(&drops));
        let detached = lazy_tuple(&source);
        let stream = std::iter::from_fn(move || {
            let _keep_alive = &probe;
            match termination {
                "error" => Some(Err(OperationError::AccessDenied {
                    operation: "drop-test",
                    block: "callback",
                })),
                "end" => None,
                _ => Some(Ok(detached.clone())),
            }
        });
        let mut parts =
            MultiOutputOpParts::<MultipleAccess>::new(&source, spec("lazy-drop", lazy_fields()))
                .unwrap();
        parts.emit_lazy_runtime(stream).unwrap();
        let mut outputs = parts.finish_lazy().unwrap();
        assert_eq!(drops.load(Ordering::SeqCst), 0);
        if termination != "drop" {
            let _ = outputs.next();
            assert_eq!(drops.load(Ordering::SeqCst), 1);
        }
        drop(outputs);
        assert_eq!(drops.load(Ordering::SeqCst), 1);
    }
}

#[cfg(all(feature = "cap-stereoisomers", feature = "cap-smiles"))]
#[test]
fn source_stereoisomer_no_center_clears_atom_code_without_finalization_and_shares_unchanged_blocks()
{
    let base = Molecule::from_smiles("CC").unwrap();
    let mut topology = base.topology().clone();
    topology.atoms[0]
        .set_computed_prop("_CIPCode", "R")
        .unwrap();
    let mut properties = base.properties().clone();
    properties.set_computed_prop("_CIPComputed", true).unwrap();
    properties
        .set_computed_prop("preexisting_computed", "retained")
        .unwrap();
    let source = Molecule::from_parts(
        topology,
        base.coordinate_block_runtime().clone(),
        properties,
    )
    .unwrap();
    let before = source.clone();
    let mut stream = source.enumerate_stereoisomers().unwrap();
    assert_eq!(stream.yielded_count(), 0);
    let output = stream.next().unwrap().unwrap();
    // Pinned EnumerateStereoisomers.py: atom.ClearProp before no-center yield;
    // the no-center return bypasses ClearComputedProps and final assignment.
    assert!(output.topology().atoms[0].prop("_CIPCode").is_none());
    assert_eq!(
        output.properties().prop("_CIPComputed"),
        source.properties().prop("_CIPComputed")
    );
    assert_eq!(
        output.properties().prop("preexisting_computed"),
        source.properties().prop("preexisting_computed")
    );
    assert_eq!(
        output.properties().prop("_MolFileChiralFlag"),
        source.properties().prop("_MolFileChiralFlag")
    );
    assert!(!Arc::ptr_eq(
        &output.topology_arc_runtime(),
        &source.topology_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &output.coordinates_arc_runtime(),
        &source.coordinates_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &output.properties_arc_runtime(),
        &source.properties_arc_runtime()
    ));
    assert_eq!(stream.yielded_count(), 1);
    assert!(stream.next().is_none());
    assert!(stream.next().is_none());
    assert_eq!(source, before);
    assert!(Arc::ptr_eq(
        &source.topology_arc_runtime(),
        &before.topology_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.derived_cache_arc_runtime(),
        &before.derived_cache_arc_runtime()
    ));
}

#[cfg(all(feature = "cap-stereoisomers", feature = "cap-smiles"))]
#[test]
fn source_stereoisomer_candidates_clear_computed_then_keep_legacy_assigned_codes() {
    let base = Molecule::from_smiles("CC(F)Cl").unwrap();
    let mut properties = base.properties().clone();
    properties.set_computed_prop("_CIPComputed", true).unwrap();
    properties
        .set_computed_prop("preexisting_computed", "discarded")
        .unwrap();
    properties.set_prop("persistent", "retained").unwrap();
    let source = Molecule::from_parts(
        base.topology().clone(),
        base.coordinate_block_runtime().clone(),
        properties,
    )
    .unwrap();
    let before = source.clone();
    let outputs = source
        .enumerate_stereoisomers()
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(outputs.len(), 2);
    let mut codes = Vec::new();
    for output in outputs {
        let Some(cosmolkit_model::PropertyValue::String(code)) =
            output.topology().atoms[1].prop("_CIPCode")
        else {
            panic!("legacy candidate lost assigned CIP code")
        };
        codes.push(code.as_bytes().to_owned());
        // Native ClearComputedProps then legacy AssignStereochemistry does not
        // recreate the unrelated modern-CIP computed marker.
        assert!(output.properties().prop("_CIPComputed").is_none());
        assert!(output.properties().prop("preexisting_computed").is_none());
        assert_eq!(
            output.properties().prop("persistent"),
            source.properties().prop("persistent")
        );
        assert_eq!(
            output.properties().prop("_StereochemDone"),
            Some(&cosmolkit_model::PropertyValue::Bool(true))
        );
        assert!(Arc::ptr_eq(
            &output.coordinates_arc_runtime(),
            &source.coordinates_arc_runtime()
        ));
    }
    codes.sort();
    assert_eq!(codes, vec![b"R".to_vec(), b"S".to_vec()]);
    assert_eq!(source, before);
    assert!(Arc::ptr_eq(
        &source.topology_arc_runtime(),
        &before.topology_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.properties_arc_runtime(),
        &before.properties_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.derived_cache_arc_runtime(),
        &before.derived_cache_arc_runtime()
    ));
}

#[cfg(all(feature = "cap-stereoisomers", feature = "cap-smiles"))]
#[test]
fn source_stereoisomer_callback_failure_is_deferred_structured_atomic_and_fused() {
    use std::error::Error;
    use std::sync::atomic::{AtomicUsize, Ordering};
    let source = Molecule::from_smiles("CC(F)=CC(Cl)C").unwrap();
    let before = source.clone();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = calls.clone();
    let options = crate::StereoisomerOptions {
        max_isomers: 3,
        ..Default::default()
    };
    let mut stream = source
        .enumerate_stereoisomers_with_random_bits(
            &options,
            Box::new(move |width| {
                assert_eq!(width, 2);
                match observed.fetch_add(1, Ordering::SeqCst) {
                    0 => Ok(num_bigint::BigUint::from(0u8)),
                    _ => Err("injected-source-error".to_owned()),
                }
            }),
        )
        .unwrap();
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    let output = stream.next().unwrap().unwrap();
    assert_eq!(stream.yielded_count(), 1);
    let error = stream.next().unwrap().unwrap_err();
    let OperationError::Enumeration(detail) = error else {
        panic!("callback failure changed category")
    };
    assert!(
        matches!(detail.source().unwrap().downcast_ref::<crate::EnumerationError>(),Some(crate::EnumerationError::RandomBitsSource{bit_count:2,message}) if message=="injected-source-error")
    );
    assert_eq!(stream.yielded_count(), 1);
    for _ in 0..3 {
        assert!(stream.next().is_none());
    }
    assert_eq!(calls.load(Ordering::SeqCst), 2);
    assert_eq!(source, before);
    assert!(Arc::ptr_eq(
        &source.topology_arc_runtime(),
        &before.topology_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.coordinates_arc_runtime(),
        &before.coordinates_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.properties_arc_runtime(),
        &before.properties_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.derived_cache_arc_runtime(),
        &before.derived_cache_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &output.coordinates_arc_runtime(),
        &source.coordinates_arc_runtime()
    ));
}

#[test]
fn reconstruction_input_accepts_valid_blocks_without_manufacturing_cache_facts() {
    let source = molecule();
    let before = source.derived_cache_runtime().clone();
    assert_eq!(validate_reconstruction_input(source.topology(), source.coordinate_block_runtime(), source.properties(), source.derived_cache_runtime()), Ok(()));
    assert_eq!(source.derived_cache_runtime(), &before);
}

#[test]
fn reconstruction_input_rejects_topology_before_coordinates_or_property_lists() {
    let mut graph = topology();
    graph.atoms[0] = atom(7, Element::C);
    assert!(matches!(validate_reconstruction_input(&graph, &CoordinateBlock { conformers_2d: vec![Conformer2D::new(3, vec![])], ..Default::default() }, &properties("bad"), &DerivedCacheBlock::default()), Err(OperationError::InvalidTopology(_))));
}

#[test]
fn reconstruction_input_checks_every_coordinate_frame() {
    let mut frames = coordinates();
    frames.conformers_2d.push(Conformer2D::new(12, vec![[0.0, 0.0]]));
    assert_eq!(validate_reconstruction_input(&topology(), &frames, &properties("valid"), &DerivedCacheBlock::default()), Err(OperationError::InvalidCoordinates(cosmolkit_model::CoordinateValidationError::RowCount { dimension: "2D", conformer: 12, rows: 1, atom_count: 2 })));
}

#[test]
fn reconstruction_input_rejects_both_property_list_target_length_errors() {
    for (target, label, expected) in [(SdfPropertyListTarget::Atom, "atom", 2), (SdfPropertyListTarget::Bond, "bond", 1)] {
        let props = MoleculeProperties::default().with_sdf_property_list(SdfPropertyList::new(target, "bad_rows", vec![]));
        assert_eq!(validate_reconstruction_input(&topology(), &coordinates(), &props, &DerivedCacheBlock::default()), Err(OperationError::InvalidPropertyList { target: label, name: "bad_rows".into(), values: 0, expected }));
    }
}

#[cfg(feature = "cap-valence")]
#[test]
fn reconstruction_input_rejects_cache_validity_without_assignment() {
    let mut cache = DerivedCacheBlock::default();
    cache.mark_valid(DerivedState::VALENCE);
    assert_eq!(validate_reconstruction_input(&topology(), &coordinates(), &properties("valid"), &cache), Err(OperationError::InvalidDerivedCache { state: "valence", field: "assignment", actual: 0, expected: 1 }));
}

#[cfg(feature = "cap-valence")]
#[test]
fn reconstruction_input_rejects_assignment_without_cache_validity() {
    let mut cache = DerivedCacheBlock::default();
    cache.install_valence_assignment(cosmolkit_core::ValenceAssignment { explicit_valence: vec![1, 1], implicit_hydrogens: vec![3, 2] });
    assert_eq!(validate_reconstruction_input(&topology(), &coordinates(), &properties("valid"), &cache), Err(OperationError::InvalidDerivedCache { state: "valence", field: "validity_bit", actual: 0, expected: 1 }));
}

#[cfg(feature = "cap-valence")]
#[test]
fn reconstruction_input_checks_both_source_valence_row_counts() {
    for (explicit_valence, implicit_hydrogens, field) in [(vec![1], vec![3, 2], "explicit_valence"), (vec![1, 1], vec![3], "implicit_hydrogens")] {
        let mut cache = DerivedCacheBlock::default();
        cache.install_valence_assignment(cosmolkit_core::ValenceAssignment { explicit_valence, implicit_hydrogens });
        cache.mark_valid(DerivedState::VALENCE);
        assert_eq!(validate_reconstruction_input(&topology(), &coordinates(), &properties("valid"), &cache), Err(OperationError::InvalidDerivedCache { state: "valence", field, actual: 1, expected: 2 }));
    }
}

#[cfg(feature = "cap-rings")]
#[test]
fn reconstruction_input_rejects_missing_ring_storage() {
    let mut cache = DerivedCacheBlock::default();
    cache.mark_valid(DerivedState::RINGS);
    assert_eq!(validate_reconstruction_input(&topology(), &coordinates(), &properties("valid"), &cache), Err(OperationError::InvalidDerivedCache { state: "rings", field: "assignment", actual: 0, expected: 1 }));
}

#[cfg(feature = "cap-valence")]
#[test]
fn reconstruction_input_preserves_actual_valid_source_facts() {
    let mut cache = DerivedCacheBlock::default();
    cache.install_valence_assignment(cosmolkit_core::ValenceAssignment { explicit_valence: vec![1, 1], implicit_hydrogens: vec![3, 2] });
    cache.mark_valid(DerivedState::VALENCE);
    let before = cache.clone();
    assert_eq!(validate_reconstruction_input(&topology(), &coordinates(), &properties("valid"), &cache), Ok(()));
    assert_eq!(cache, before);
}

#[cfg(feature = "cap-reaction")]
fn typed_reconstruction_parts(source: &Molecule) -> MultiOutputOpParts<'_, crate::ReactionProductsFromInputsAccess> {
    MultiOutputOpParts::new(source, &super::super::registry::REACTION_PRODUCTS_FROM_INPUTS_SPEC).unwrap()
}

#[cfg(feature = "cap-reaction")]
fn typed_reconstruction_product() -> cosmolkit_reaction::ReactionProduct {
    use cosmolkit_reaction::{ReactionProduct, ReactionRowOrigin};
    ReactionProduct {
        topology: topology(), coordinates: coordinates(), properties: properties("product"),
        atom_origins: vec![Some(ReactionRowOrigin { input: 0, row: AtomId::new(0) }), Some(ReactionRowOrigin { input: 0, row: AtomId::new(1) })],
        bond_origins: vec![Some(ReactionRowOrigin { input: 0, row: BondId::new(0) })],
        valence: cosmolkit_core::ValenceAssignment { explicit_valence: vec![1, 1], implicit_hydrogens: vec![3, 2] },
        rings: None,
    }
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reconstruction_typed_inputs_retain_actual_order_duplicates_counts_and_borrows() {
    let source = molecule();
    let other = Molecule::from_parts(TopologyBlock::try_from_parts(vec![atom(0, Element::O)], vec![], vec![], vec![]).unwrap(), CoordinateBlock::default(), MoleculeProperties::default()).unwrap();
    let mut parts = typed_reconstruction_parts(&source);
    let inputs = parts.reconstruction_inputs_runtime(&[&source, &other, &source]).unwrap();
    assert_eq!(parts.reconstruction_inputs, [(2, 1), (1, 0), (2, 1)]);
    assert_eq!(inputs.len(), 3);
    assert!(std::ptr::eq(inputs[0].topology, source.topology()));
    assert!(std::ptr::eq(inputs[1].topology, other.topology()));
    assert!(std::ptr::eq(inputs[2].topology, source.topology()));
    assert!(std::ptr::eq(inputs[0].coordinates, source.coordinate_block_runtime()));
    assert!(std::ptr::eq(inputs[2].properties, source.properties()));
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reconstruction_typed_inputs_reject_second_declaration_without_replacing_evidence() {
    let source = molecule();
    let mut parts = typed_reconstruction_parts(&source);
    parts.reconstruction_inputs_runtime(&[&source, &source]).unwrap();
    assert!(matches!(parts.reconstruction_inputs_runtime(&[]), Err(OperationError::MappingContract { issue: "reconstruction input set must be declared exactly once", .. })));
    assert_eq!(parts.reconstruction_inputs, [(2, 1), (2, 1)]);
    assert!(parts.reconstruction_inputs_read);
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reconstruction_typed_source_uses_same_checked_single_input_boundary() {
    let source = molecule();
    let mut parts = MultiOutputOpParts::<crate::ReactionProductsAccess>::new(&source, &super::super::registry::REACTION_PRODUCTS_SPEC).unwrap();
    let input = parts.reconstruction_source_runtime().unwrap();
    assert!(std::ptr::eq(input.topology, source.topology()));
    assert!(std::ptr::eq(input.coordinates, source.coordinate_block_runtime()));
    assert_eq!(parts.reconstruction_inputs, [(2, 1)]);
    assert!(parts.reconstruction_inputs_read);
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reconstruction_typed_candidate_accepts_repeated_origins_and_none_rows() {
    let source = molecule();
    let mut product = typed_reconstruction_product();
    product.atom_origins[1] = product.atom_origins[0];
    product.bond_origins[0] = None;
    let mut parts = typed_reconstruction_parts(&source);
    parts.reconstruction_inputs_runtime(&[&source]).unwrap();
    parts.emit_reconstructed_runtime(vec![product]).unwrap();
    let outputs = parts.finish().unwrap();
    assert_eq!(outputs.len(), 1);
    assert_eq!(outputs[0].topology(), &topology());
    assert_eq!(outputs[0].coordinate_block_runtime(), &coordinates());
    assert_eq!(outputs[0].properties(), &properties("product"));
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reconstruction_typed_candidate_checks_later_product_before_returning_any_molecule() {
    let source = molecule();
    let observer = source.clone();
    let good = typed_reconstruction_product();
    let mut bad = typed_reconstruction_product();
    bad.atom_origins[0].as_mut().unwrap().input = 999;
    let mut parts = typed_reconstruction_parts(&source);
    parts.reconstruction_inputs_runtime(&[&source]).unwrap();
    parts.emit_reconstructed_runtime(vec![good, bad]).unwrap();
    assert!(matches!(parts.finish(), Err(OperationError::InvalidReconstructionOrigin { entity: "atom", destination: 0, input: 999, input_count: 1, row_count: None, .. })));
    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    assert!(std::ptr::eq(source.coordinate_block_runtime(), observer.coordinate_block_runtime()));
    assert!(std::ptr::eq(source.properties(), observer.properties()));
    assert!(Arc::ptr_eq(&source.derived_cache_arc_runtime(), &observer.derived_cache_arc_runtime()));
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reconstruction_typed_candidate_checks_actual_derived_fact_lengths() {
    let source = molecule();
    let mut bad = typed_reconstruction_product();
    bad.valence.implicit_hydrogens.pop();
    let mut parts = typed_reconstruction_parts(&source);
    parts.emit_reconstructed_runtime(vec![bad]).unwrap();
    assert!(matches!(parts.finish(), Err(OperationError::InvalidAlgorithmResult { field: "implicit hydrogens", actual: 1, expected: 2, .. })));
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reconstruction_typed_candidate_checks_actual_property_rows() {
    let source = molecule();
    let mut bad = typed_reconstruction_product();
    bad.properties = MoleculeProperties::default().with_sdf_property_list(SdfPropertyList::new(SdfPropertyListTarget::Bond, "bad_rows", vec![]));
    let mut parts = typed_reconstruction_parts(&source);
    parts.emit_reconstructed_runtime(vec![bad]).unwrap();
    assert!(matches!(parts.finish(), Err(OperationError::InvalidPropertyList { target: "bond", values: 0, expected: 1, .. })));
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reconstruction_typed_grouping_preserves_empty_sets_duplicate_values_and_order() {
    let source = molecule();
    let mut first = typed_reconstruction_product();
    first.properties = properties("first");
    let mut second = typed_reconstruction_product();
    second.properties = properties("second");
    let mut parts = typed_reconstruction_parts(&source);
    parts.emit_reconstructed_runtime(vec![first, second.clone(), second]).unwrap();
    let groups = crate::reaction::assemble_product_sets(parts.finish().unwrap(), vec![1, 0, 2]).unwrap();
    assert_eq!(groups.iter().map(Vec::len).collect::<Vec<_>>(), [1, 0, 2]);
    assert_eq!(groups[0][0].properties(), &properties("first"));
    assert_eq!(groups[2][0].properties(), &properties("second"));
    assert_eq!(groups[2][0], groups[2][1]);
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reconstruction_typed_grouping_rejects_incomplete_or_oversized_metadata() {
    for lengths in [vec![], vec![0], vec![2]] {
        assert!(matches!(crate::reaction::assemble_product_sets(vec![molecule()], lengths), Err(OperationError::InvalidAlgorithmResult { .. })));
    }
}

fn reconstruction_check(
    inputs: &[(usize, usize)],
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    atom_origins: &[Option<(usize, usize)>],
    bond_origins: &[Option<(usize, usize)>],
    valence_rows: (usize, usize),
    ring_rows: Option<(usize, usize)>,
) -> Result<(), OperationError> {
    let mut fields = base_fields();
    fields.mapping = MappingRequirement::Reconstruction;
    fields.edit = TopologyEditKind::Reconstruction;
    validate_reconstruction_rows(
        spec("reconstruction", fields),
        inputs,
        topology,
        coordinates,
        atom_origins.iter().copied(),
        bond_origins.iter().copied(),
        valence_rows,
        ring_rows,
    )
}

#[test]
fn reconstruction_origins_accept_one_to_many_and_independent_source_inputs() {
    assert_eq!(
        reconstruction_check(
            &[(2, 1), (1, 0)],
            &topology(),
            &coordinates(),
            &[Some((1, 0)), Some((1, 0))],
            &[Some((0, 0))],
            (2, 2),
            Some((2, 1)),
        ),
        Ok(()),
    );
}

#[test]
fn reconstruction_origins_accept_new_atom_and_bond_rows_without_sources() {
    assert_eq!(
        reconstruction_check(
            &[],
            &topology(),
            &coordinates(),
            &[None, None],
            &[None],
            (2, 2),
            None
        ),
        Ok(()),
    );
}

#[test]
fn reconstruction_origins_accept_empty_product_and_actual_zero_row_facts() {
    assert_eq!(
        reconstruction_check(
            &[],
            &TopologyBlock::default(),
            &CoordinateBlock::default(),
            &[],
            &[],
            (0, 0),
            Some((0, 0))
        ),
        Ok(()),
    );
}

#[test]
fn reconstruction_origins_reject_atom_input_index_with_destination_and_no_row_count() {
    assert_eq!(
        reconstruction_check(
            &[(2, 1)],
            &topology(),
            &coordinates(),
            &[None, Some((1, 0))],
            &[None],
            (2, 2),
            None
        ),
        Err(OperationError::InvalidReconstructionOrigin {
            operation: "reconstruction",
            entity: "atom",
            destination: 1,
            input: 1,
            row: 0,
            input_count: 1,
            row_count: None,
        }),
    );
}

#[test]
fn reconstruction_origins_reject_atom_row_at_exact_source_length() {
    assert_eq!(
        reconstruction_check(
            &[(2, 1)],
            &topology(),
            &coordinates(),
            &[None, Some((0, 2))],
            &[None],
            (2, 2),
            None
        ),
        Err(OperationError::InvalidReconstructionOrigin {
            operation: "reconstruction",
            entity: "atom",
            destination: 1,
            input: 0,
            row: 2,
            input_count: 1,
            row_count: Some(2),
        }),
    );
}

#[test]
fn reconstruction_origins_reject_bond_input_index_with_actual_input_count() {
    assert_eq!(
        reconstruction_check(
            &[(2, 1)],
            &topology(),
            &coordinates(),
            &[None, None],
            &[Some((3, 0))],
            (2, 2),
            None
        ),
        Err(OperationError::InvalidReconstructionOrigin {
            operation: "reconstruction",
            entity: "bond",
            destination: 0,
            input: 3,
            row: 0,
            input_count: 1,
            row_count: None,
        }),
    );
}

#[test]
fn reconstruction_origins_reject_bond_row_at_exact_source_length() {
    assert_eq!(
        reconstruction_check(
            &[(2, 1)],
            &topology(),
            &coordinates(),
            &[None, None],
            &[Some((0, 1))],
            (2, 2),
            None
        ),
        Err(OperationError::InvalidReconstructionOrigin {
            operation: "reconstruction",
            entity: "bond",
            destination: 0,
            input: 0,
            row: 1,
            input_count: 1,
            row_count: Some(1),
        }),
    );
}

#[test]
fn reconstruction_origins_check_source_bond_counts_independently_of_atom_counts() {
    assert_eq!(
        reconstruction_check(
            &[(2, 0)],
            &topology(),
            &coordinates(),
            &[Some((0, 0)), Some((0, 1))],
            &[Some((0, 0))],
            (2, 2),
            None
        ),
        Err(OperationError::InvalidReconstructionOrigin {
            operation: "reconstruction",
            entity: "bond",
            destination: 0,
            input: 0,
            row: 0,
            input_count: 1,
            row_count: Some(0),
        }),
    );
}

#[test]
fn reconstruction_origins_reject_all_four_destination_and_valence_length_mismatches() {
    for (atoms, bonds, valence, field, actual, expected) in [
        (vec![None], vec![None], (2, 2), "atom origins", 1, 2),
        (vec![None, None], vec![], (2, 2), "bond origins", 0, 1),
        (
            vec![None, None],
            vec![None],
            (1, 2),
            "explicit valence",
            1,
            2,
        ),
        (
            vec![None, None],
            vec![None],
            (2, 1),
            "implicit hydrogens",
            1,
            2,
        ),
    ] {
        assert_eq!(
            reconstruction_check(
                &[],
                &topology(),
                &coordinates(),
                &atoms,
                &bonds,
                valence,
                None
            ),
            Err(OperationError::InvalidAlgorithmResult {
                operation: "reconstruction",
                field,
                actual,
                expected
            }),
        );
    }
}

#[test]
fn reconstruction_origins_length_error_precedes_invalid_source_origin() {
    assert_eq!(
        reconstruction_check(
            &[],
            &topology(),
            &coordinates(),
            &[Some((usize::MAX, usize::MAX))],
            &[None],
            (2, 2),
            None
        ),
        Err(OperationError::InvalidAlgorithmResult {
            operation: "reconstruction",
            field: "atom origins",
            actual: 1,
            expected: 2
        }),
    );
}

#[test]
fn reconstruction_origins_first_atom_error_precedes_bond_error_and_ring_rows() {
    assert_eq!(
        reconstruction_check(
            &[],
            &topology(),
            &coordinates(),
            &[Some((4, 0)), Some((5, 0))],
            &[Some((6, 0))],
            (2, 2),
            Some((0, 0))
        ),
        Err(OperationError::InvalidReconstructionOrigin {
            operation: "reconstruction",
            entity: "atom",
            destination: 0,
            input: 4,
            row: 0,
            input_count: 0,
            row_count: None,
        }),
    );
}

#[test]
fn reconstruction_origins_reject_both_optional_ring_membership_row_mismatches() {
    for (rings, field, actual, expected) in [
        ((1, 1), "ring atom membership", 1, 2),
        ((2, 0), "ring bond membership", 0, 1),
    ] {
        assert_eq!(
            reconstruction_check(
                &[],
                &topology(),
                &coordinates(),
                &[None, None],
                &[None],
                (2, 2),
                Some(rings)
            ),
            Err(OperationError::InvalidAlgorithmResult {
                operation: "reconstruction",
                field,
                actual,
                expected
            }),
        );
    }
}

#[test]
fn reconstruction_origins_reject_invalid_physical_topology_before_coordinate_rows() {
    let mut broken = topology();
    broken.atoms[0] = atom(7, Element::C);
    let bad_coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(9, vec![[0.0, 0.0]])],
        ..Default::default()
    };
    assert_eq!(
        reconstruction_check(
            &[],
            &broken,
            &bad_coordinates,
            &[None, None],
            &[None],
            (2, 2),
            None
        ),
        Err(OperationError::InvalidTopology(
            cosmolkit_model::TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(7)
            }
        )),
    );
}

#[test]
fn reconstruction_origins_validate_later_physical_coordinate_frames() {
    let coordinates = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(9, vec![[0.0, 0.0], [1.0, 0.0]]),
            Conformer2D::new(12, vec![[0.0, 0.0]]),
        ],
        ..Default::default()
    };
    assert_eq!(
        reconstruction_check(
            &[],
            &topology(),
            &coordinates,
            &[None, None],
            &[None],
            (2, 2),
            None
        ),
        Err(OperationError::InvalidCoordinates(
            cosmolkit_model::CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: 12,
                rows: 1,
                atom_count: 2
            }
        )),
    );
}

#[test]
fn reconstruction_origins_allow_absent_rings_and_exact_present_rows() {
    for rings in [None, Some((2, 1))] {
        assert_eq!(
            reconstruction_check(
                &[(2, 1)],
                &topology(),
                &coordinates(),
                &[Some((0, 1)), Some((0, 0))],
                &[Some((0, 0))],
                (2, 2),
                rings
            ),
            Ok(()),
        );
    }
}

#[test]
fn reconstruction_origins_leave_all_borrowed_values_unchanged_after_success_and_failure() {
    let topology = topology();
    let coordinates = coordinates();
    let original = (topology.clone(), coordinates.clone());
    let inputs = [(2, 1)];
    let atoms = [Some((0, 0)), Some((0, 0))];
    let bonds = [Some((0, 0))];
    assert!(
        reconstruction_check(
            &inputs,
            &topology,
            &coordinates,
            &atoms,
            &bonds,
            (2, 2),
            None
        )
        .is_ok()
    );
    assert!(
        reconstruction_check(
            &inputs,
            &topology,
            &coordinates,
            &atoms,
            &bonds,
            (2, 1),
            None
        )
        .is_err()
    );
    assert_eq!((topology, coordinates), original);
    assert_eq!(inputs, [(2, 1)]);
    assert_eq!(atoms, [Some((0, 0)), Some((0, 0))]);
    assert_eq!(bonds, [Some((0, 0))]);
}
