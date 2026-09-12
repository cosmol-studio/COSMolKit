use std::sync::Arc;

use cosmolkit_macros::{mol_op_body, molecule_ops};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Conformer2D, CoordinateBlock, Element, MoleculeProperties,
    TopologyBlock,
};

use super::{FeatureSpec, OperationError, SupportStatus};
use crate::Molecule;

const COW_TEST_FEATURE: FeatureSpec = FeatureSpec {
    name: "cow-runtime-test",
    category: "internal-test",
    status: SupportStatus::Experimental,
    rdkit_parity_sensitive: false,
    docs: "Internal registered-operation coverage for block-level COW.",
};

#[mol_op_body(cow_coordinates_for_test, parts)]
fn cow_coordinates_for_test_impl() -> Result<(), OperationError> {
    let mut coordinates = parts.checkout_coordinates()?;
    coordinates.conformers_2d[0].coordinates_mut()[0] = [9.0, 8.0];
    parts.install_coordinates(coordinates)?;
    parts.apply_cip_policy()
}

#[mol_op_body(cow_coordinates_failure_for_test, parts)]
fn cow_coordinates_failure_for_test_impl() -> Result<(), OperationError> {
    let mut coordinates = parts.checkout_coordinates()?;
    coordinates.conformers_2d[0].coordinates_mut()[0] = [7.0, 6.0];
    Err(OperationError::Algorithm {
        operation: "cow_coordinates_failure_for_test",
        detail: "intentional failure after detachment".to_owned(),
    })
}

molecule_ops! {
    op cow_coordinates_for_test {
        method: cow_coordinates_for_test,
        impl_fn: crate::ops::cow_tests::cow_coordinates_for_test_impl,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        derived_effects: {
            recompute: [], preserve: [], invalidate: [], operation_defined: [],
        },
        cip_state: preserve,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::ops::cow_tests::COW_TEST_FEATURE,
        parity: not_applicable,
        io_roundtrip: false,
        invariant_profile: "cow-coordinate-write-test",
    }

    op cow_coordinates_failure_for_test {
        method: cow_coordinates_failure_for_test,
        impl_fn: crate::ops::cow_tests::cow_coordinates_failure_for_test_impl,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        derived_effects: {
            recompute: [], preserve: [], invalidate: [], operation_defined: [],
        },
        cip_state: preserve,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::ops::cow_tests::COW_TEST_FEATURE,
        parity: not_applicable,
        io_roundtrip: false,
        invariant_profile: "cow-coordinate-failure-test",
        inplace: true,
        inplace_method: cow_coordinates_failure_for_test_,
    }
}

fn molecule() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("test topology is valid");
    Molecule::from_parts(
        topology,
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(3, vec![[1.0, 2.0]])],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default().with_name("source"),
    )
    .expect("test molecule is valid")
}

#[test]
fn registered_value_operation_shares_untouched_blocks_and_detaches_only_the_write() {
    let source = molecule();
    let output = source.cow_coordinates_for_test().unwrap();

    assert!(std::ptr::eq(source.topology(), output.topology()));
    assert!(std::ptr::eq(source.properties(), output.properties()));
    assert!(Arc::ptr_eq(
        &source.derived_cache_arc_runtime(),
        &output.derived_cache_arc_runtime()
    ));
    assert!(!std::ptr::eq(source.coordinates(), output.coordinates()));
    assert_eq!(
        source.coordinates().conformers_2d[0].coordinates()[0],
        [1.0, 2.0]
    );
    assert_eq!(
        output.coordinates().conformers_2d[0].coordinates()[0],
        [9.0, 8.0]
    );
}

#[test]
fn registered_in_place_failure_keeps_every_original_block_allocation_and_value() {
    let mut target = molecule();
    let observer = target.clone();
    let before = target.clone();

    let error = target.cow_coordinates_failure_for_test_().unwrap_err();
    assert_eq!(
        error,
        OperationError::Algorithm {
            operation: "cow_coordinates_failure_for_test",
            detail: "intentional failure after detachment".to_owned(),
        }
    );
    assert_eq!(target, before);
    assert!(std::ptr::eq(target.topology(), observer.topology()));
    assert!(std::ptr::eq(target.coordinates(), observer.coordinates()));
    assert!(std::ptr::eq(target.properties(), observer.properties()));
    assert!(Arc::ptr_eq(
        &target.derived_cache_arc_runtime(),
        &observer.derived_cache_arc_runtime()
    ));
}
