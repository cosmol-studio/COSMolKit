//! Additional proposed runtime/COW cases; original effect negatives unchanged.
use super::*;
use crate::{Atom, AtomId, AtomSpec, BlockAccess, Element};

struct Probe;
fn operation() -> MoleculeOpSpec {
    *crate::operation_spec("with_atom_pair_atom_code").unwrap()
}
fn source() -> Molecule {
    Molecule::from_parts(
        TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap()
}

#[test]
fn new_exact_cip_method_rejects_each_missing_authority_and_unrelated_methods() {
    let valid = operation();
    OpParts::<Probe>::validate_effect_contract(Box::leak(Box::new(valid))).unwrap();
    for missing in [BlockSet::TOPOLOGY, BlockSet::PROPERTIES] {
        for write_permission in [true, false] {
            let mut invalid = valid;
            if write_permission {
                invalid.access =
                    BlockAccess::new(BlockSet::NONE, valid.access.write().difference(missing));
            } else {
                invalid.may_mutate = valid.may_mutate.difference(missing);
            }
            assert!(matches!(
                OpParts::<Probe>::validate_effect_contract(Box::leak(Box::new(invalid))),
                Err(OperationError::CipStateContract { .. })
            ));
        }
    }
    for method in [
        "wrong",
        "with_atom_pair_atom_code_with_params",
        "with_cip_labels",
        "atom_code",
    ] {
        let mut invalid = valid;
        invalid.method = method;
        assert!(matches!(
            OpParts::<Probe>::validate_effect_contract(Box::leak(Box::new(invalid))),
            Err(OperationError::CipStateContract { .. })
        ));
    }
}

#[test]
fn no_effect_scoped_callback_retains_every_arc_and_error_leaves_complete_slots() {
    let source = source();
    let spec = Box::leak(Box::new(operation()));
    let mut parts = OpParts::<Probe>::new(&source, spec).unwrap();
    let (answer, changed) = parts
        .stage_topology_properties_runtime(|topology, properties, _cache| {
            assert!(matches!(topology, Cow::Borrowed(_)));
            assert!(matches!(properties, Cow::Borrowed(_)));
            Ok((17, None))
        })
        .unwrap();
    assert_eq!(answer, 17);
    assert!(!changed);
    parts
        .record_topology_edit_runtime(TopologyEditKind::Local)
        .unwrap();
    let output = parts.finish().unwrap();
    assert!(Arc::ptr_eq(
        &source.topology_arc_runtime(),
        &output.topology_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.properties_arc_runtime(),
        &output.properties_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.coordinates_arc_runtime(),
        &output.coordinates_arc_runtime()
    ));
    assert!(Arc::ptr_eq(
        &source.derived_cache_arc_runtime(),
        &output.derived_cache_arc_runtime()
    ));
    let mut parts = OpParts::<Probe>::new(&source, spec).unwrap();
    let error = parts
        .stage_topology_properties_runtime::<()>(|topology, properties, _| {
            let mut candidate = topology.into_owned();
            candidate.atoms[0].set_prop("probe", "local").unwrap();
            let mut candidate_properties = properties.into_owned();
            candidate_properties.set_prop("probe", "local").unwrap();
            Err(OperationError::Algorithm {
                operation: "probe",
                detail: "source failure".into(),
            })
        })
        .unwrap_err();
    assert!(matches!(error, OperationError::Algorithm { .. }));
    parts.ensure_complete_blocks().unwrap();
    parts
        .record_topology_edit_runtime(TopologyEditKind::Local)
        .unwrap();
    let untouched = parts.finish().unwrap();
    assert_eq!(untouched, source);
    assert!(Arc::ptr_eq(
        &source.derived_cache_arc_runtime(),
        &untouched.derived_cache_arc_runtime()
    ));
}

#[test]
fn new_exact_method_still_requires_computed_assignment_evidence() {
    let source = source();
    let spec = Box::leak(Box::new(operation()));
    for state in [None, Some(false), Some(true)] {
        let mut parts = OpParts::<Probe>::new(&source, spec).unwrap();
        let mut properties = parts.checkout_properties_runtime().unwrap();
        if let Some(value) = state {
            properties
                .set_prop("_CIPComputed", if value { "1" } else { "0" })
                .unwrap();
        }
        parts.install_properties_runtime(properties).unwrap();
        assert!(matches!(
            parts.apply_cip_policy_runtime(),
            Err(OperationError::CipStateContract {
                issue: "assignment did not install computed _CIPComputed evidence",
                ..
            })
        ));
    }
    let mut parts = OpParts::<Probe>::new(&source, spec).unwrap();
    let mut properties = parts.checkout_properties_runtime().unwrap();
    properties.set_computed_prop("_CIPComputed", "1").unwrap();
    parts.install_properties_runtime(properties).unwrap();
    parts.apply_cip_policy_runtime().unwrap();
}
