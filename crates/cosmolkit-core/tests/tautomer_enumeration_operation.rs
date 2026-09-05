use std::sync::Mutex;

use cosmolkit_core::{
    BlockSet, CipStatePolicy, DerivedState, InvariantCheckSet, MOLECULE_OPS, MappingRequirement, Molecule,
    MoleculeOpKind, MoleculeOpOutput, OPERATION_INVARIANT_MATRIX, OperationDomain, OperationError, PARITY_MATRIX,
    ParityPolicy, SUPPORT_MATRIX, SupportStatus, TopologyEditKind,
    chemistry::tautomer::{
        TautomerEnumeration, TautomerEnumerationCallback, TautomerEnumerationStatus, TautomerEnumerator,
        TautomerOptions, TautomerRunError,
    },
};

fn tautomer_operation() -> &'static cosmolkit_core::MoleculeOpSpec {
    MOLECULE_OPS
        .iter()
        .copied()
        .find(|operation| operation.method == "enumerate_tautomers_with_options")
        .expect("tautomer enumeration must be registered")
}

#[test]
fn registry_declares_the_complete_multiple_output_tautomer_contract() {
    let operation = tautomer_operation();

    assert_eq!(operation.impl_fn, "enumerate_tautomers_with_options_impl");
    assert_eq!(operation.output, MoleculeOpOutput::Multiple);
    assert_eq!(
        operation.result_type,
        "crate :: chemistry :: tautomer :: TautomerEnumeration"
    );
    assert_eq!(operation.domain, OperationDomain::Topology);
    assert_eq!(operation.kind, MoleculeOpKind::Weak);
    assert_eq!(operation.topology_edit, TopologyEditKind::Local);
    assert_eq!(operation.access.read(), BlockSet::COORDINATES);
    assert_eq!(
        operation.access.write(),
        BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(operation.may_mutate, operation.access.write());
    assert_eq!(operation.auto_remap, BlockSet::NONE);
    assert_eq!(operation.requires_mapping, MappingRequirement::None);
    assert_eq!(operation.cip_state, CipStatePolicy::TautomerSourceTransition);
    assert_eq!(
        operation.derived_effects.recompute(),
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
    );
    assert_eq!(operation.derived_effects.preserve(), DerivedState::NONE);
    assert_eq!(
        operation.derived_effects.invalidate(),
        DerivedState::DRAWING.union(DerivedState::FINGERPRINT)
    );
    assert_eq!(operation.derived_effects.operation_defined(), DerivedState::NONE);
    assert_eq!(operation.support, SupportStatus::Experimental);
    assert_eq!(operation.parity, ParityPolicy::RequiredWhenSupported);
    assert!(operation.io_roundtrip);

    let support = SUPPORT_MATRIX
        .iter()
        .find(|entry| {
            entry
                .operation
                .is_some_and(|candidate| candidate.method == operation.method)
        })
        .expect("tautomer operation must have a support entry");
    assert_eq!(support.feature.name, "tautomer.enumeration");

    let invariant = OPERATION_INVARIANT_MATRIX
        .iter()
        .find(|entry| entry.operation.method == operation.method)
        .expect("tautomer operation must have an invariant entry");
    assert_eq!(invariant.profile, "tautomer_enumeration_multi_output");
    assert_eq!(
        invariant.required_checks,
        InvariantCheckSet::GRAPH_INDEX_VALIDITY
            .union(InvariantCheckSet::COORDINATE_ROW_ALIGNMENT)
            .union(InvariantCheckSet::MAPPING_RECORDED_OR_EXPLICITLY_UNSUPPORTED)
            .union(InvariantCheckSet::STEREO_VALIDITY)
            .union(InvariantCheckSet::CACHE_INVALIDATION)
            .union(InvariantCheckSet::PROPERTY_POLICY)
            .union(InvariantCheckSet::SOURCE_UNCHANGED)
            .union(InvariantCheckSet::UNSUPPORTED_IS_EXPLICIT)
    );

    let parity = PARITY_MATRIX
        .iter()
        .find(|entry| entry.operation.method == operation.method)
        .expect("tautomer operation must have a parity entry");
    assert_eq!(parity.profile, "rdkit_tautomer_enumeration_full_plan");
    assert_eq!(parity.rdkit_version, None);
}

#[test]
fn enumeration_emits_canonical_key_order_and_aligned_rich_result_metadata() {
    let source = Molecule::from_smiles("CC(C)=O").expect("parse acetone");
    let result = TautomerEnumerator::new().enumerate(&source).expect("enumerate acetone");

    assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
    assert_eq!(result.len(), 2);
    assert_eq!(result.canonical_smiles(), ["C=C(C)O", "CC(C)=O"]);
    assert_eq!(result.get(result.len()), None);
    assert!(result.at(result.len()).is_err());

    let entries = result.iter_with_smiles().collect::<Vec<_>>();
    assert!(entries.windows(2).all(|pair| pair[0].0 < pair[1].0));
    for (index, (canonical_smiles, molecule)) in entries.iter().enumerate() {
        assert_eq!(molecule.to_smiles(true).unwrap(), *canonical_smiles);
        assert_eq!(result.get(index), Some(*molecule));
        assert_eq!(&result[index], *molecule);
    }
    assert!(!result.modified_atoms().is_empty());
    assert!(!result.modified_bonds().is_empty());
    assert!(
        result
            .modified_atoms()
            .iter()
            .all(|atom| atom.index() < source.num_atoms())
    );
    assert!(
        result
            .modified_bonds()
            .iter()
            .all(|bond| bond.index() < source.num_bonds())
    );

    let copied = result.molecules();
    assert_eq!(copied.len(), result.len());
    assert!(copied.iter().zip(result.iter()).all(|(left, right)| left == right));
    assert_eq!(result.clone(), result);
}

#[test]
fn enumeration_preserves_source_value_coordinates_properties_and_output_invariants() {
    let source = Molecule::from_smiles("CC(C)=O")
        .expect("parse coordinate fixture")
        .with_2d_coordinates()
        .expect("generate 2D coordinates");
    let coordinates_3d = (0..source.num_atoms())
        .map(|index| [index as f64, index as f64 + 0.25, -(index as f64)])
        .collect::<Vec<_>>();
    let source = source
        .with_only_3d_conformer(coordinates_3d, true)
        .expect("attach 3D coordinates")
        .with_name("coordinate-property-source")
        .with_prop("source_id", "tautomer-operation")
        .with_sdf_data_field("dataset", "operation-level");
    let before = source.clone();
    let expected_2d = source.conformers_2d().to_vec();
    let expected_3d = source.conformers_3d().to_vec();
    let expected_properties = source.properties().clone();

    let first = TautomerEnumerator::new()
        .enumerate(&source)
        .expect("enumerate decorated source");
    let second = TautomerEnumerator::new()
        .enumerate(&source)
        .expect("repeat decorated source enumeration");

    assert_eq!(source, before, "value-style enumeration mutated its source");
    assert_eq!(first, second, "repeated execution must be deterministic");
    for output in first.iter() {
        assert_eq!(output.num_atoms(), source.num_atoms());
        assert_eq!(output.num_bonds(), source.num_bonds());
        assert_eq!(output.conformers_2d(), expected_2d);
        assert_eq!(output.conformers_3d(), expected_3d);
        assert_eq!(output.properties(), &expected_properties);
        assert_eq!(output.prop("source_id"), Some("tautomer-operation"));
        assert_eq!(output.conformers_2d()[0].coordinates().len(), output.num_atoms());
        assert_eq!(output.conformers_3d()[0].coordinates().len(), output.num_atoms());
        output
            .with_assigned_valence()
            .expect("every emitted branch remains valid for a strict topology operation");
    }
}

#[test]
fn zero_and_exact_limits_preserve_source_defined_status_and_initial_branch() {
    let source = Molecule::from_smiles("CC(C)=O").expect("parse limit fixture");
    let source_smiles = source.to_smiles(true).unwrap();

    for max_transforms in [0, 1] {
        let enumerator =
            TautomerEnumerator::from_options(TautomerOptions::default().with_max_transforms(max_transforms));
        let result = enumerator.enumerate(&source).expect("transform-limited run");
        assert_eq!(result.status(), TautomerEnumerationStatus::MaxTransformsReached);
        assert_eq!(result.canonical_smiles(), [source_smiles.clone()]);
    }

    for max_tautomers in [0, 1] {
        let enumerator = TautomerEnumerator::from_options(TautomerOptions::default().with_max_tautomers(max_tautomers));
        let result = enumerator.enumerate(&source).expect("tautomer-limited run");
        assert_eq!(result.status(), TautomerEnumerationStatus::MaxTautomersReached);
        assert_eq!(result.canonical_smiles(), [source_smiles.clone()]);
    }
    assert_eq!(source.to_smiles(true).unwrap(), source_smiles);
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CallbackObservation {
    source_smiles: String,
    source_id: Option<String>,
    candidate_smiles: Vec<String>,
    status: TautomerEnumerationStatus,
    modified_atom_count: usize,
    modified_bond_count: usize,
}

#[derive(Debug, Default)]
struct CancelAtFirstMatch {
    observations: Mutex<Vec<CallbackObservation>>,
}

impl TautomerEnumerationCallback for CancelAtFirstMatch {
    fn should_continue(&self, molecule: &Molecule, result: &TautomerEnumeration) -> bool {
        self.observations.lock().unwrap().push(CallbackObservation {
            source_smiles: molecule.to_smiles(true).unwrap(),
            source_id: molecule.prop("source_id").map(str::to_owned),
            candidate_smiles: result.canonical_smiles(),
            status: result.status(),
            modified_atom_count: result.modified_atoms().len(),
            modified_bond_count: result.modified_bonds().len(),
        });
        false
    }
}

#[test]
fn callback_observes_the_pre_application_result_and_cancellation_status() {
    let callback = CancelAtFirstMatch::default();
    let source = Molecule::from_smiles("CC(C)=O")
        .expect("parse callback fixture")
        .with_prop("source_id", "callback-source");
    let before = source.clone();
    let mut enumerator = TautomerEnumerator::new();
    enumerator.set_callback(Some(&callback));

    let result = enumerator.enumerate(&source).expect("cancel enumeration");

    assert_eq!(source, before);
    assert_eq!(result.status(), TautomerEnumerationStatus::Canceled);
    assert_eq!(result.len(), 1);
    assert_eq!(
        callback.observations.into_inner().unwrap(),
        [CallbackObservation {
            source_smiles: "CC(C)=O".to_owned(),
            source_id: Some("callback-source".to_owned()),
            candidate_smiles: vec!["CC(C)=O".to_owned()],
            status: TautomerEnumerationStatus::Completed,
            modified_atom_count: 0,
            modified_bond_count: 0,
        }]
    );
}

#[test]
fn invalid_aromatic_input_returns_a_structured_error_without_mutating_source() {
    let source = Molecule::from_smiles_with_sanitize("c1cccc1", false)
        .expect("parse deliberately unkekulizable odd aromatic ring")
        .with_prop("source_id", "invalid-aromatic");
    let before = source.clone();

    let error = TautomerEnumerator::new()
        .enumerate(&source)
        .expect_err("invalid aromatic input must fail explicitly");

    assert!(matches!(
        error,
        OperationError::Tautomer {
            operation,
            source,
        } if operation.method == "enumerate_tautomers_with_options"
            && matches!(*source, TautomerRunError::Kekulize(_))
    ));
    assert_eq!(source, before, "failed enumeration mutated its source");
}
