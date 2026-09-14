//! Public coordinate-view assertions; full metadata is inspected on detached
//! builder snapshots, not through an extra live Molecule accessor.
#![allow(dead_code)]

use cosmolkit::Molecule;

pub fn assert_shared_coordinates(left: &Molecule, right: &Molecule) {
    assert_eq!(
        left.to_builder().coordinates(),
        right.to_builder().coordinates()
    );
    match (left.coordinates_2d(), right.coordinates_2d()) {
        (Some(left), Some(right)) => assert!(std::ptr::eq(left, right)),
        (None, None) => {}
        _ => panic!("2D coordinate presence changed"),
    }
    assert!(std::ptr::eq(left.conformers_3d(), right.conformers_3d()));
    for (left, right) in left.conformers_3d().iter().zip(right.conformers_3d()) {
        assert!(std::ptr::eq(left.coordinates(), right.coordinates()));
    }
}

pub fn assert_detached_coordinates(left: &Molecule, right: &Molecule) {
    let detached_2d = match (left.coordinates_2d(), right.coordinates_2d()) {
        (Some(left), Some(right)) => !std::ptr::eq(left, right),
        (None, None) => false,
        _ => true,
    };
    let detached_3d = !std::ptr::eq(left.conformers_3d(), right.conformers_3d());
    assert!(
        detached_2d || detached_3d,
        "no visible coordinate storage detached"
    );
}
