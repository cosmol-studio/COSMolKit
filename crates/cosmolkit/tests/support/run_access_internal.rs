use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Conformer2D, CoordinateBlock, CoordinateValidationError,
    Element, MoleculeProperties, TopologyBlock, TopologyValidationError,
};

use super::*;
use crate::ops::{
    BlockAccess, CipStatePolicy, DerivedEffects, DerivedState, FunctionStatus, MappingRequirement,
    MoleculeOpKind, OperationDomain, ParityPolicy, SemanticPreconditionSet, TopologyEditKind,
};

struct TestAccess;

#[derive(cosmolkit_macros::MoleculeResult)]
struct PendingReport<M = Molecule> {
    label: String,
    #[pending_molecule]
    molecule: Option<M>,
}

#[derive(cosmolkit_macros::MoleculeResult)]
struct RequiredPendingReport<M = Molecule> {
    #[pending_molecule]
    molecule: M,
}

fn pending_spec() -> &'static MoleculeOpSpec {
    spec(
        "pending_test",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::PROPERTIES),
    )
}

#[test]
fn pending_some_and_required_results_commit_detached_blocks_and_preserve_cow() {
    let source = molecule();
    let mut parts = OpParts::<TestAccess>::new(&source, pending_spec()).unwrap();
    let properties = parts
        .checkout_properties_runtime()
        .unwrap()
        .with_name("candidate");
    parts.install_properties_runtime(properties).unwrap();
    parts.apply_cip_policy_runtime().unwrap();
    let pending = parts.pending_molecule_runtime().unwrap();
    assert!(matches!(pending.topology, WorkingBlock::Shared));
    let result = parts
        .finish_result(PendingReport {
            label: "metadata".into(),
            molecule: Some(pending),
        })
        .unwrap();
    assert_eq!(result.label, "metadata");
    let finished = result.molecule.unwrap();
    assert_eq!(
        finished.properties().name().map(|value| value.as_bytes()),
        Some(b"candidate".as_slice())
    );
    assert_eq!(
        source.properties().name().map(|value| value.as_bytes()),
        Some(b"source".as_slice())
    );
    assert!(std::ptr::eq(finished.topology(), source.topology()));
    assert!(std::ptr::eq(
        finished.coordinate_block_runtime(),
        source.coordinate_block_runtime()
    ));

    let mut parts = OpParts::<TestAccess>::new(&source, pending_spec()).unwrap();
    parts.apply_cip_policy_runtime().unwrap();
    let molecule = parts.pending_molecule_runtime().unwrap();
    let required = parts
        .finish_result(RequiredPendingReport { molecule })
        .unwrap();
    assert!(std::ptr::eq(
        required.molecule.topology(),
        source.topology()
    ));
}

#[test]
fn pending_none_still_finishes_and_missing_blocks_cannot_be_hidden() {
    let source = molecule();
    let mut parts = OpParts::<TestAccess>::new(&source, pending_spec()).unwrap();
    parts.apply_cip_policy_runtime().unwrap();
    let result = parts
        .finish_result(PendingReport::<PendingMolecule<TestAccess>> {
            label: "analysis".into(),
            molecule: None,
        })
        .unwrap();
    assert_eq!(result.label, "analysis");
    assert!(result.molecule.is_none());

    let mut parts = OpParts::<TestAccess>::new(&source, pending_spec()).unwrap();
    let _detached = parts.checkout_properties_runtime().unwrap();
    assert!(parts.pending_molecule_runtime().is_err());
    assert!(
        parts
            .finish_result(PendingReport::<PendingMolecule<TestAccess>> {
                label: String::new(),
                molecule: None,
            })
            .is_err()
    );
}

#[test]
fn pending_duplicate_seal_dropped_candidate_and_foreign_transaction_are_rejected() {
    let source = molecule();
    let declaration = pending_spec();
    let mut parts = OpParts::<TestAccess>::new(&source, declaration).unwrap();
    let pending = parts.pending_molecule_runtime().unwrap();
    assert!(parts.pending_molecule_runtime().is_err());
    assert!(parts.ensure_unsealed_runtime().is_err());
    drop(pending);
    assert!(
        parts
            .finish_result(PendingReport::<PendingMolecule<TestAccess>> {
                label: String::new(),
                molecule: None,
            })
            .is_err()
    );

    let mut first = OpParts::<TestAccess>::new(&source, declaration).unwrap();
    let mut second = OpParts::<TestAccess>::new(&source, declaration).unwrap();
    let first_pending = first.pending_molecule_runtime().unwrap();
    let _second_pending = second.pending_molecule_runtime().unwrap();
    let error = operation_error(second.finish_result(RequiredPendingReport {
        molecule: first_pending,
    }));
    assert!(matches!(
        error,
        OperationError::IncompleteCommit {
            block: "foreign pending molecule",
            ..
        }
    ));
    assert_eq!(
        source.properties().name().map(|value| value.as_bytes()),
        Some(b"source".as_slice())
    );
}

#[test]
fn pending_body_error_drops_candidate_without_touching_source() {
    let source = molecule();
    let before = source.clone();
    let failed = (|| -> Result<(), OperationError> {
        let mut parts = OpParts::<TestAccess>::new(&source, pending_spec())?;
        let properties = parts
            .checkout_properties_runtime()?
            .with_name("not committed");
        parts.install_properties_runtime(properties)?;
        let _pending = parts.pending_molecule_runtime()?;
        Err(OperationError::IncompleteCommit {
            operation: "pending_test",
            block: "body failed",
        })
    })();
    assert!(failed.is_err());
    assert_eq!(source, before);
    assert!(std::ptr::eq(source.properties(), before.properties()));
}

#[cfg(feature = "op-contracts-strict")]
#[test]
fn pending_some_and_none_cannot_skip_strict_contract_validation() {
    let source = molecule();
    let mut parts = OpParts::<TestAccess>::new(&source, pending_spec()).unwrap();
    let properties = parts
        .checkout_properties_runtime()
        .unwrap()
        .with_name("staged");
    parts.install_properties_runtime(properties).unwrap();
    // Intentionally omit mandatory CIP bookkeeping. Sealing is not validation.
    let pending = parts.pending_molecule_runtime().unwrap();
    assert!(
        parts
            .finish_result(RequiredPendingReport { molecule: pending })
            .is_err()
    );
    let mut parts = OpParts::<TestAccess>::new(&source, pending_spec()).unwrap();
    let properties = parts
        .checkout_properties_runtime()
        .unwrap()
        .with_name("staged");
    parts.install_properties_runtime(properties).unwrap();
    assert!(
        parts
            .finish_result(PendingReport::<PendingMolecule<TestAccess>> {
                label: String::new(),
                molecule: None,
            })
            .is_err()
    );
}

#[cfg(feature = "cap-stereo")]
#[test]
fn pending_generated_capabilities_are_frozen_after_sealing() {
    let source = molecule();
    let mut parts =
        OpParts::<crate::ops::PotentialStereoAccess>::new(&source, pending_spec()).unwrap();
    let _pending = parts.pending_molecule().unwrap();
    assert!(parts.pending_molecule().is_err());
    assert!(parts.checkout_topology().is_err());
    assert!(parts.install_topology(TopologyBlock::default()).is_err());
    assert!(parts.checkout_properties().is_err());
    assert!(parts.checkout_derived_cache().is_err());
    assert!(parts.clear_cache(DerivedState::STEREO).is_err());
    assert!(parts.apply_cip_policy().is_err());
    assert!(
        parts
            .prove_preserved(DerivedState::RINGS, PreservationProof::UnchangedInput)
            .is_err()
    );
}

fn spec(
    method: &'static str,
    output: MoleculeOpOutput,
    access: BlockAccess,
) -> &'static MoleculeOpSpec {
    Box::leak(Box::new(MoleculeOpSpec {
        method,
        impl_fn: "test_impl",
        output,
        result_type: "Molecule",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::None,
        access,
        may_mutate: access.write(),
        auto_remap: BlockSet::NONE,
        derived_effects: DerivedEffects::new(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::NONE,
        ),
        cip_state: CipStatePolicy::Preserve,
        semantic_preconditions: SemanticPreconditionSet::NONE,
        requires_mapping: MappingRequirement::None,
        status: FunctionStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    }))
}

fn atom(index: usize) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C))
}

fn topology(atom_count: usize) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        (0..atom_count).map(atom).collect(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("test topology is valid")
}

fn molecule() -> Molecule {
    Molecule::from_parts(
        topology(1),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(7, vec![[1.0, 2.0]])],
            ..Default::default()
        },
        MoleculeProperties::default().with_name("source"),
    )
    .expect("test molecule is valid")
}

#[test]
fn conditional_cow_source_borrows_keep_shared_blocks() {
    let source = molecule();
    let operation = spec(
        "conditional-source-borrow",
        MoleculeOpOutput::Single,
        BlockAccess::new(
            BlockSet::NONE,
            BlockSet::TOPOLOGY
                .union(BlockSet::PROPERTIES)
                .union(BlockSet::DERIVED_CACHE),
        ),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, operation).unwrap();
    let (value, changed) = parts
        .stage_topology_properties_cow_runtime(|topology, properties, _| {
            Ok((17, Some((topology, properties))))
        })
        .unwrap();
    assert_eq!(value, 17);
    assert!(changed);
    assert!(matches!(parts.topology, WorkingBlock::Shared));
    assert!(matches!(parts.properties, WorkingBlock::Shared));
    assert!(std::ptr::eq(
        parts.current_topology_candidate().unwrap(),
        source.topology()
    ));
    assert!(std::ptr::eq(
        parts.current_properties_candidate().unwrap(),
        source.properties()
    ));
}

#[test]
fn conditional_cow_still_validates_owned_and_foreign_candidates() {
    let source = molecule();
    let operation = spec(
        "conditional-untrusted-candidate",
        MoleculeOpOutput::Single,
        BlockAccess::new(
            BlockSet::NONE,
            BlockSet::TOPOLOGY
                .union(BlockSet::PROPERTIES)
                .union(BlockSet::DERIVED_CACHE),
        ),
    );
    for owned in [false, true] {
        for malformed in [false, true] {
            let mut candidate = topology(1);
            if malformed {
                // Keep the old adjacency, breaking only its row alignment.
                candidate.atoms.push(atom(1));
            }
            let mut parts = OpParts::<TestAccess>::new(&source, operation).unwrap();
            let candidate = if owned {
                Cow::Owned(candidate)
            } else {
                // A static foreign borrow satisfies the callback's HRTB but
                // must still be rejected as unrelated to its source borrow.
                Cow::Borrowed(&*Box::leak(Box::new(candidate)))
            };
            let result = parts.stage_topology_properties_cow_runtime(|_, properties, _| {
                Ok(((), Some((candidate, properties))))
            });
            if malformed {
                // Structural errors retain priority over foreign-borrow errors.
                assert!(matches!(
                    result,
                    Err(OperationError::InvalidTopology(
                        TopologyValidationError::AdjacencyMismatch
                    ))
                ));
            } else if owned {
                assert_eq!(result.unwrap(), ((), true));
                assert!(matches!(parts.topology, WorkingBlock::Installed(_)));
            } else {
                assert!(matches!(
                    result,
                    Err(OperationError::IncompleteCommit {
                        block: "foreign borrowed topology candidate",
                        ..
                    })
                ));
            }
            if malformed || !owned {
                assert!(matches!(parts.topology, WorkingBlock::Shared));
                assert!(matches!(parts.properties, WorkingBlock::Shared));
            }
        }
    }
}

fn denied(operation: &'static str, block: &'static str) -> OperationError {
    OperationError::AccessDenied { operation, block }
}

fn checked_out(operation: &'static str, block: &'static str) -> OperationError {
    OperationError::BlockCheckedOut { operation, block }
}

fn not_checked_out(operation: &'static str, block: &'static str) -> OperationError {
    OperationError::BlockNotCheckedOut { operation, block }
}

fn construction_error(result: Result<OpParts<'_, TestAccess>, OperationError>) -> OperationError {
    match result {
        Ok(_) => panic!("operation context construction unexpectedly succeeded"),
        Err(error) => error,
    }
}

fn operation_error<T>(result: Result<T, OperationError>) -> OperationError {
    match result {
        Ok(_) => panic!("operation unexpectedly succeeded"),
        Err(error) => error,
    }
}

#[test]
fn constructors_are_lazy_and_reject_multiple_output_before_exposure() {
    let source = molecule();
    let single = spec(
        "lazy",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let parts = OpParts::<TestAccess>::new(&source, single).unwrap();
    assert!(std::ptr::eq(parts.source.topology(), source.topology()));
    assert!(std::ptr::eq(
        parts.source.coordinate_block_runtime(),
        source.coordinate_block_runtime()
    ));
    assert!(std::ptr::eq(parts.source.properties(), source.properties()));
    assert!(matches!(parts.topology, WorkingBlock::Shared));
    assert!(matches!(parts.coordinates, WorkingBlock::Shared));
    assert!(matches!(parts.properties, WorkingBlock::Shared));
    assert!(matches!(parts.derived_cache, WorkingBlock::Shared));

    let multiple = spec(
        "multiple",
        MoleculeOpOutput::Multiple,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    assert_eq!(
        construction_error(OpParts::<TestAccess>::new(&source, multiple)),
        OperationError::OutputMismatch {
            operation: "multiple",
            expected: MoleculeOpOutput::Single,
            actual: MoleculeOpOutput::Multiple,
        }
    );

    let mut target = molecule();
    let topology_ptr = target.topology() as *const TopologyBlock;
    {
        let parts = OpParts::<TestAccess>::new_in_place(&mut target, single).unwrap();
        assert_eq!(
            parts.source.topology() as *const TopologyBlock,
            topology_ptr
        );
    }
    assert_eq!(
        target.properties().name().map(|value| value.as_bytes()),
        Some(b"source".as_slice())
    );
    assert_eq!(
        construction_error(OpParts::<TestAccess>::new_in_place(&mut target, multiple,)),
        OperationError::OutputMismatch {
            operation: "multiple",
            expected: MoleculeOpOutput::Single,
            actual: MoleculeOpOutput::Multiple,
        }
    );
}

#[test]
fn topology_permissions_and_lifecycle_are_exact() {
    let source = molecule();
    let none = spec(
        "topology-none",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, none).unwrap();
    assert_eq!(
        parts.read_topology_runtime().unwrap_err(),
        denied("topology-none", "topology")
    );
    assert_eq!(
        parts.checkout_topology_runtime().unwrap_err(),
        denied("topology-none", "topology")
    );
    assert_eq!(
        parts.install_topology_runtime(topology(1)).unwrap_err(),
        denied("topology-none", "topology")
    );

    let read = spec(
        "topology-read",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::TOPOLOGY, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, read).unwrap();
    assert_eq!(parts.read_topology_runtime().unwrap().atoms.len(), 1);
    assert_eq!(
        parts.checkout_topology_runtime().unwrap_err(),
        denied("topology-read", "topology")
    );

    let write = spec(
        "topology-write",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, write).unwrap();
    assert_eq!(parts.read_topology_runtime().unwrap().atoms.len(), 1);
    assert_eq!(
        parts.install_topology_runtime(topology(1)).unwrap_err(),
        not_checked_out("topology-write", "topology")
    );
    let detached = parts.checkout_topology_runtime().unwrap();
    assert_eq!(detached.atoms.len(), 1);
    assert_eq!(
        parts.read_topology_runtime().unwrap_err(),
        checked_out("topology-write", "topology")
    );
    assert_eq!(
        parts.checkout_topology_runtime().unwrap_err(),
        checked_out("topology-write", "topology")
    );
    parts.install_topology_runtime(detached).unwrap();
    assert_eq!(parts.read_topology_runtime().unwrap().atoms.len(), 1);
    assert_eq!(
        parts.install_topology_runtime(topology(1)).unwrap_err(),
        not_checked_out("topology-write", "topology")
    );
}

#[test]
fn coordinate_permissions_and_lifecycle_are_exact() {
    let source = molecule();
    let none = spec(
        "coordinates-none",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, none).unwrap();
    assert_eq!(
        parts.read_coordinates_runtime().unwrap_err(),
        denied("coordinates-none", "coordinates")
    );
    assert_eq!(
        parts.checkout_coordinates_runtime().unwrap_err(),
        denied("coordinates-none", "coordinates")
    );
    assert_eq!(
        parts
            .install_coordinates_runtime(CoordinateBlock::default())
            .unwrap_err(),
        denied("coordinates-none", "coordinates")
    );

    let read = spec(
        "coordinates-read",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::COORDINATES, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, read).unwrap();
    assert_eq!(
        parts
            .read_coordinates_runtime()
            .unwrap()
            .conformers_2d
            .len(),
        1
    );
    assert_eq!(
        parts.checkout_coordinates_runtime().unwrap_err(),
        denied("coordinates-read", "coordinates")
    );

    let write = spec(
        "coordinates-write",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::COORDINATES),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, write).unwrap();
    assert_eq!(
        parts
            .read_coordinates_runtime()
            .unwrap()
            .conformers_2d
            .len(),
        1
    );
    assert_eq!(
        parts
            .install_coordinates_runtime(CoordinateBlock::default())
            .unwrap_err(),
        not_checked_out("coordinates-write", "coordinates")
    );
    let detached = parts.checkout_coordinates_runtime().unwrap();
    assert_eq!(
        parts.read_coordinates_runtime().unwrap_err(),
        checked_out("coordinates-write", "coordinates")
    );
    assert_eq!(
        parts.checkout_coordinates_runtime().unwrap_err(),
        checked_out("coordinates-write", "coordinates")
    );
    parts.install_coordinates_runtime(detached).unwrap();
    assert_eq!(
        parts
            .read_coordinates_runtime()
            .unwrap()
            .conformers_2d
            .len(),
        1
    );
    assert_eq!(
        parts
            .install_coordinates_runtime(CoordinateBlock::default())
            .unwrap_err(),
        not_checked_out("coordinates-write", "coordinates")
    );
}

#[test]
fn property_permissions_and_lifecycle_are_exact() {
    let source = molecule();
    let none = spec(
        "properties-none",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, none).unwrap();
    assert_eq!(
        parts.read_properties_runtime().unwrap_err(),
        denied("properties-none", "properties")
    );
    assert_eq!(
        parts.checkout_properties_runtime().unwrap_err(),
        denied("properties-none", "properties")
    );
    assert_eq!(
        parts
            .install_properties_runtime(MoleculeProperties::default())
            .unwrap_err(),
        denied("properties-none", "properties")
    );

    let read = spec(
        "properties-read",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::PROPERTIES, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, read).unwrap();
    assert_eq!(
        parts
            .read_properties_runtime()
            .unwrap()
            .name()
            .map(|value| value.as_bytes()),
        Some(b"source".as_slice())
    );
    assert_eq!(
        parts.checkout_properties_runtime().unwrap_err(),
        denied("properties-read", "properties")
    );

    let write = spec(
        "properties-write",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::PROPERTIES),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, write).unwrap();
    assert_eq!(
        parts
            .read_properties_runtime()
            .unwrap()
            .name()
            .map(|value| value.as_bytes()),
        Some(b"source".as_slice())
    );
    assert_eq!(
        parts
            .install_properties_runtime(MoleculeProperties::default())
            .unwrap_err(),
        not_checked_out("properties-write", "properties")
    );
    let detached = parts.checkout_properties_runtime().unwrap();
    assert_eq!(
        parts.read_properties_runtime().unwrap_err(),
        checked_out("properties-write", "properties")
    );
    assert_eq!(
        parts.checkout_properties_runtime().unwrap_err(),
        checked_out("properties-write", "properties")
    );
    parts.install_properties_runtime(detached).unwrap();
    assert_eq!(
        parts
            .read_properties_runtime()
            .unwrap()
            .name()
            .map(|value| value.as_bytes()),
        Some(b"source".as_slice())
    );
    assert_eq!(
        parts
            .install_properties_runtime(MoleculeProperties::default())
            .unwrap_err(),
        not_checked_out("properties-write", "properties")
    );
}

#[test]
fn derived_cache_permissions_and_lifecycle_are_exact() {
    let source = molecule();
    let none = spec(
        "cache-none",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, none).unwrap();
    assert_eq!(
        operation_error(parts.read_derived_cache_runtime()),
        denied("cache-none", "derived_cache")
    );
    assert_eq!(
        operation_error(parts.checkout_derived_cache_runtime()),
        denied("cache-none", "derived_cache")
    );
    assert_eq!(
        parts
            .install_derived_cache_runtime(DerivedCacheBlock::default())
            .unwrap_err(),
        denied("cache-none", "derived_cache")
    );

    let read = spec(
        "cache-read",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::DERIVED_CACHE, BlockSet::NONE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, read).unwrap();
    parts.read_derived_cache_runtime().unwrap();
    assert_eq!(
        operation_error(parts.checkout_derived_cache_runtime()),
        denied("cache-read", "derived_cache")
    );

    let write = spec(
        "cache-write",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::DERIVED_CACHE),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, write).unwrap();
    parts.read_derived_cache_runtime().unwrap();
    assert_eq!(
        parts
            .install_derived_cache_runtime(DerivedCacheBlock::default())
            .unwrap_err(),
        not_checked_out("cache-write", "derived_cache")
    );
    let detached = parts.checkout_derived_cache_runtime().unwrap();
    assert_eq!(
        operation_error(parts.read_derived_cache_runtime()),
        checked_out("cache-write", "derived_cache")
    );
    assert_eq!(
        operation_error(parts.checkout_derived_cache_runtime()),
        checked_out("cache-write", "derived_cache")
    );
    parts.install_derived_cache_runtime(detached).unwrap();
    parts.read_derived_cache_runtime().unwrap();
    assert_eq!(
        parts
            .install_derived_cache_runtime(DerivedCacheBlock::default())
            .unwrap_err(),
        not_checked_out("cache-write", "derived_cache")
    );
}

#[test]
fn checkout_bookkeeping_is_independent_and_sources_stay_unchanged() {
    let source = molecule();
    let all = BlockSet::TOPOLOGY
        .union(BlockSet::COORDINATES)
        .union(BlockSet::PROPERTIES)
        .union(BlockSet::DERIVED_CACHE);
    let access = spec(
        "independent",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, all),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, access).unwrap();
    let topology = parts.checkout_topology_runtime().unwrap();
    let coordinates = parts.checkout_coordinates_runtime().unwrap();
    assert_eq!(
        parts
            .read_properties_runtime()
            .unwrap()
            .name()
            .map(|value| value.as_bytes()),
        Some(b"source".as_slice())
    );
    parts.install_topology_runtime(topology).unwrap();
    assert_eq!(parts.read_topology_runtime().unwrap().atoms.len(), 1);
    assert_eq!(
        parts.read_coordinates_runtime().unwrap_err(),
        checked_out("independent", "coordinates")
    );
    parts.install_coordinates_runtime(coordinates).unwrap();
    assert_eq!(source.num_atoms(), 1);
    assert_eq!(
        source.properties().name().map(|value| value.as_bytes()),
        Some(b"source".as_slice())
    );

    let mut target = molecule();
    {
        let mut parts = OpParts::<TestAccess>::new_in_place(&mut target, access).unwrap();
        let properties = parts.checkout_properties_runtime().unwrap();
        parts
            .install_properties_runtime(properties.with_name("working"))
            .unwrap();
        assert_eq!(
            parts
                .read_properties_runtime()
                .unwrap()
                .name()
                .map(|value| value.as_bytes()),
            Some(b"working".as_slice())
        );
    }
    assert_eq!(
        target.properties().name().map(|value| value.as_bytes()),
        Some(b"source".as_slice())
    );
}

#[test]
fn invalid_replacements_are_rejected_without_replacing_working_or_live_state() {
    let source = molecule();
    let topology_write = spec(
        "invalid-topology",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, topology_write).unwrap();
    parts.checkout_topology_runtime().unwrap();
    let invalid_topology = TopologyBlock {
        atoms: vec![atom(1)],
        bonds: Vec::new(),
        adjacency: AdjacencyList::from_topology(1, &[]),
        substance_groups: Vec::new(),
        stereo_groups: Vec::new(),
    };
    assert_eq!(
        parts
            .install_topology_runtime(invalid_topology)
            .unwrap_err(),
        OperationError::InvalidTopology(TopologyValidationError::AtomIdMismatch {
            position: 0,
            id: AtomId::new(1),
        })
    );
    assert_eq!(
        parts.read_topology_runtime().unwrap_err(),
        checked_out("invalid-topology", "topology")
    );
    assert_eq!(source.num_atoms(), 1);

    let coordinate_write = spec(
        "invalid-coordinates",
        MoleculeOpOutput::Single,
        BlockAccess::new(BlockSet::NONE, BlockSet::COORDINATES),
    );
    let mut parts = OpParts::<TestAccess>::new(&source, coordinate_write).unwrap();
    parts.checkout_coordinates_runtime().unwrap();
    let invalid_coordinates = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(7, vec![[0.0, 0.0]]),
            Conformer2D::new(8, vec![[0.0, 0.0], [1.0, 1.0]]),
        ],
        ..Default::default()
    };
    assert_eq!(
        parts
            .install_coordinates_runtime(invalid_coordinates)
            .unwrap_err(),
        OperationError::InvalidCoordinates(CoordinateValidationError::RowCount {
            dimension: "2D",
            conformer: 8,
            rows: 2,
            atom_count: 1,
        })
    );
    assert_eq!(
        parts.read_coordinates_runtime().unwrap_err(),
        checked_out("invalid-coordinates", "coordinates")
    );
    assert_eq!(source.coordinate_block_runtime().conformers_2d.len(), 1);
}

#[cfg(all(feature = "cap-stereo", feature = "op-contracts-strict"))]
fn cip_metadata_preservation_source(with_existing_computed: bool) -> Molecule {
    let mut builder = crate::MoleculeBuilder::new();
    let mut first = AtomSpec::new(Element::C).with_prop("user", 7_i32).unwrap();
    if with_existing_computed {
        first = first.with_computed_prop("existing", "kept").unwrap();
    }
    let a = builder.add_atom(first);
    let b = builder.add_atom(AtomSpec::new(Element::C));
    let mut bond = crate::BondSpec::new(a, b, crate::BondOrder::Single)
        .with_prop("user", 9_i32)
        .unwrap();
    if with_existing_computed {
        bond = bond.with_computed_prop("existing", "kept").unwrap();
    }
    builder.add_bond(bond).unwrap();
    builder.build().unwrap()
}
#[cfg(all(feature = "cap-stereo", feature = "op-contracts-strict"))]
fn cip_metadata_preservation_check(
    source: &Molecule,
    candidate: TopologyBlock,
) -> Result<(), OperationError> {
    let declaration = &crate::ops::runtime::registry::WITH_CIP_LABELS_SPEC;
    let mut parts = OpParts::<TestAccess>::new(source, declaration)?;
    let _old = parts.checkout_topology_runtime()?;
    parts.install_topology_runtime(candidate)?;
    parts.prove_preserved_runtime(DerivedState::RINGS, PreservationProof::CipLabelAssignment)
}
#[cfg(all(feature = "cap-stereo", feature = "op-contracts-strict"))]
#[test]
fn cip_metadata_preservation_allows_only_owned_computed_lifecycle() {
    for existing in [false, true] {
        let source = cip_metadata_preservation_source(existing);
        let before = source.clone();
        let mut candidate = source.topology().clone();
        candidate.atoms[0].set_prop("_CIPCode", "R").unwrap();
        candidate.atoms[0]
            .set_computed_prop("_CIPNeighborOrder", "[1]")
            .unwrap();
        candidate.bonds[0].set_prop("_CIPCode", "E").unwrap();
        candidate.bonds[0]
            .set_computed_prop("_CIPNeighborOrder", "[0,1]")
            .unwrap();
        assert_eq!(cip_metadata_preservation_check(&source, candidate), Ok(()));
        assert_eq!(source, before);
        assert!(std::ptr::eq(source.topology(), before.topology()));
    }
}
#[cfg(all(feature = "cap-stereo", feature = "op-contracts-strict"))]
#[test]
fn cip_metadata_preservation_rejects_unowned_values_and_membership() {
    let source = cip_metadata_preservation_source(true);
    let before = source.clone();
    for change in 0..5 {
        let mut candidate = source.topology().clone();
        candidate.atoms[0]
            .set_computed_prop("_CIPNeighborOrder", "[1]")
            .unwrap();
        match change {
            0 => candidate.atoms[0].set_prop("user", 8_i32).unwrap(),
            1 => candidate.bonds[0].set_prop("user", 8_i32).unwrap(),
            2 => candidate.atoms[0].clear_prop("existing").unwrap(),
            3 => candidate.bonds[0].clear_prop("existing").unwrap(),
            4 => candidate.atoms[0].set_prop("_CIPRank", 3_i32).unwrap(),
            _ => unreachable!(),
        }
        assert!(
            matches!(
                cip_metadata_preservation_check(&source, candidate),
                Err(OperationError::DerivedEffectContract { .. })
            ),
            "unowned change {change}"
        );
        assert_eq!(source, before);
    }
}
#[cfg(all(feature = "cap-stereo", feature = "op-contracts-strict"))]
#[test]
fn cip_metadata_preservation_rejects_invented_or_malformed_metadata() {
    let source = cip_metadata_preservation_source(false);
    let mut invented = source.topology().clone();
    invented.atoms[0]
        .set_prop(
            "__computedProps",
            crate::PropertyValue::StringVector(vec![]),
        )
        .unwrap();
    assert!(matches!(
        cip_metadata_preservation_check(&source, invented),
        Err(OperationError::DerivedEffectContract { .. })
    ));
    let mut malformed = source.topology().clone();
    malformed.atoms[0]
        .set_prop("__computedProps", "wrong-tag")
        .unwrap();
    let error = cip_metadata_preservation_check(&source, malformed).unwrap_err();
    assert!(matches!(
        error,
        OperationError::AtomProperty(cosmolkit_model::AtomPropertyError::ComputedListKind(_))
    ));
    assert!(std::error::Error::source(&error).is_some());
}
