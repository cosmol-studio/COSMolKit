use std::sync::Arc;

use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, CoordinateBlock,
    Element, MoleculeProperties, SdfPropertyList, SdfPropertyListTarget, TopologyBlock,
};

use super::*;
use crate::molecule::DerivedCacheBlock;
use crate::ops::{
    BlockAccess, CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement, MoleculeOpKind,
    OperationDomain, ParityPolicy, SemanticPreconditionSet, SupportStatus, TopologyEditKind,
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
            vec![Some("c".to_owned()), Some("n".to_owned())],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bond_rows",
            vec![Some("cn".to_owned())],
        ))
}

fn molecule() -> Molecule {
    Molecule::from_parts(topology(), coordinates(), properties("source")).unwrap()
}

fn tuple(source: &Molecule) -> DetachedCandidate {
    (
        source.topology().clone(),
        source.coordinates().clone(),
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
        let parts = MultiOutputOpParts::<MultipleAccess>::new(
            &source,
            spec(denied_name, fields),
        )
        .unwrap();
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
        source.coordinates().clone(),
        properties("one"),
    )])
    .unwrap();
    let one = one.finish().unwrap();
    assert_eq!(one.len(), 1);
    assert_eq!(one[0].properties().name(), Some("one"));

    let candidates = ["first", "duplicate", "duplicate", "last"]
        .into_iter()
        .map(|name| {
            (
                source.topology().clone(),
                source.coordinates().clone(),
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
        .map(|candidate| candidate.properties().name().unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(names, ["first", "duplicate", "duplicate", "last"]);
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
            source.coordinates().clone(),
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
            source.coordinates().clone(),
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
            source.coordinates().clone(),
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
            source.coordinates().clone(),
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
fn invalidate_clears_each_candidate_cache_without_touching_source() {
    let base = molecule();
    let mut cache = DerivedCacheBlock::default();
    cache.mark_valid(DerivedState::RINGS.union(DerivedState::STEREO));
    let source = Molecule::from_runtime_parts(
        Arc::new(base.topology().clone()),
        Arc::new(base.coordinates().clone()),
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
        assert!(!output
            .derived_cache_runtime()
            .valid_states()
            .contains(DerivedState::RINGS));
        assert!(output
            .derived_cache_runtime()
            .valid_states()
            .contains(DerivedState::STEREO));
    }
    assert_eq!(source, before);
    assert!(source
        .derived_cache_runtime()
        .valid_states()
        .contains(DerivedState::RINGS));
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
            source.coordinates().clone(),
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
    computed
        .set_computed_prop("_CIPComputed", "true")
        .unwrap();
    let mut clear_fields = base_fields();
    clear_fields.access = BlockAccess::new(
        BlockSet::COORDINATES,
        BlockSet::TOPOLOGY.union(BlockSet::PROPERTIES),
    );
    clear_fields.cip = CipStatePolicy::ClearComputed;
    let mut clear = MultiOutputOpParts::<MultipleAccess>::new(
        &source,
        spec("clear-cip", clear_fields),
    )
    .unwrap();
    clear
        .emit_all_runtime(vec![(
            source.topology().clone(),
            source.coordinates().clone(),
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
        spec("enumerate_tautomers_with_options", tautomer_fields),
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
