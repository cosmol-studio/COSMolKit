use std::sync::Arc;

use cosmolkit_macros::mol_op_body;
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Conformer2D, CoordinateBlock, Element, MoleculeProperties,
    TopologyBlock,
};

use super::OperationError;
use crate::Molecule;

#[mol_op_body(cow_coordinates_for_test, parts)]
pub(crate) fn cow_coordinates_for_test_impl() -> Result<(), OperationError> {
    let mut coordinates = parts.checkout_coordinates()?;
    coordinates.conformers_2d[0].coordinates_mut()[0] = [9.0, 8.0];
    parts.install_coordinates(coordinates)?;
    parts.apply_cip_policy()
}

#[mol_op_body(cow_coordinates_failure_for_test, parts)]
pub(crate) fn cow_coordinates_failure_for_test_impl() -> Result<(), OperationError> {
    let mut coordinates = parts.checkout_coordinates()?;
    coordinates.conformers_2d[0].coordinates_mut()[0] = [7.0, 6.0];
    Err(OperationError::Algorithm {
        operation: "cow_coordinates_failure_for_test",
        detail: "intentional failure after detachment".to_owned(),
    })
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
    assert!(!std::ptr::eq(
        source.coordinate_block_runtime(),
        output.coordinate_block_runtime()
    ));
    assert_eq!(
        source.coordinate_block_runtime().conformers_2d[0].coordinates()[0],
        [1.0, 2.0]
    );
    assert_eq!(
        output.coordinate_block_runtime().conformers_2d[0].coordinates()[0],
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
    assert!(std::ptr::eq(
        target.coordinate_block_runtime(),
        observer.coordinate_block_runtime()
    ));
    assert!(std::ptr::eq(target.properties(), observer.properties()));
    assert!(Arc::ptr_eq(
        &target.derived_cache_arc_runtime(),
        &observer.derived_cache_arc_runtime()
    ));
}
