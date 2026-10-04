use std::sync::Arc;

use cosmolkit_macros::mol_op_body;
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Conformer2D, CoordinateBlock, Element, MoleculeProperties,
    TopologyBlock,
};

use super::OperationError;
use crate::{DerivedState, Molecule};

/// Actual-site cfg(test) history for ring-cache installs: append-only,
/// never reset. Each entry records BOTH row-buffer addresses of the
/// finder-return RingInfo immediately before it is moved into the
/// derived cache.
#[cfg(any(
    feature = "cap-rings",
    feature = "cap-smiles",
    feature = "cap-sanitize",
    feature = "cap-hydrogens",
    feature = "cap-aromaticity",
    feature = "cap-kekulize",
    feature = "cap-descriptors"
))]
pub(crate) mod ring_live_probe {
    use cosmolkit_core::RingInfo;
    use std::cell::RefCell;

    thread_local! {
        static INSTALL_HISTORY: RefCell<Vec<(usize, usize)>> =
            const { RefCell::new(Vec::new()) };
    }

    pub(crate) fn record_install(rings: &RingInfo) {
        INSTALL_HISTORY.with(|history| {
            history.borrow_mut().push((
                rings.atom_rings().as_ptr() as usize,
                rings.bond_rings().as_ptr() as usize,
            ));
        });
    }

    pub(crate) fn history_len() -> usize {
        INSTALL_HISTORY.with(|history| history.borrow().len())
    }

    pub(crate) fn history_after(baseline: usize) -> Vec<(usize, usize)> {
        INSTALL_HISTORY.with(|history| {
            let history = history.borrow();
            history[baseline.min(history.len())..].to_vec()
        })
    }
}

/// Actual-site cfg(test) history for the Kekulize transaction: append-only.
/// Each entry records BOTH row-buffer addresses of the ordinary-ring value
/// at the moment the body holds it (the reused checked-out cache before
/// re-install, or the moved ring_update before installation).
#[cfg(feature = "cap-kekulize")]
pub(crate) mod ring_kekulize_probe {
    use cosmolkit_core::RingInfo;
    use std::cell::RefCell;

    thread_local! {
        static HISTORY: RefCell<Vec<(usize, usize)>> =
            const { RefCell::new(Vec::new()) };
    }

    pub(crate) fn record(rings: &RingInfo) {
        HISTORY.with(|history| {
            history.borrow_mut().push((
                rings.atom_rings().as_ptr() as usize,
                rings.bond_rings().as_ptr() as usize,
            ));
        });
    }

    pub(crate) fn len() -> usize {
        HISTORY.with(|history| history.borrow().len())
    }

    pub(crate) fn after(baseline: usize) -> Vec<(usize, usize)> {
        HISTORY.with(|history| {
            let history = history.borrow();
            history[baseline.min(history.len())..].to_vec()
        })
    }
}

/// Actual-site cfg(test) counter for root-level symmetrized_sssr
/// acquisitions in the aromaticity body: append-only, never reset.
#[cfg(feature = "cap-aromaticity")]
pub(crate) mod ring_aromaticity_probe {
    use std::cell::Cell;

    thread_local! {
        static ACQUISITIONS: Cell<u64> = const { Cell::new(0) };
    }

    pub(crate) fn record_acquisition() {
        ACQUISITIONS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn acquisitions() -> u64 {
        ACQUISITIONS.with(Cell::get)
    }
}

#[mol_op_body(ring_live_cow_checkout_conflict_for_test, parts)]
pub(crate) fn ring_live_cow_checkout_conflict_for_test_impl() -> Result<(), OperationError> {
    // The derived cache is CheckedOut for the rest of this body; a cache
    // mark/clear must keep returning the same structural failure.
    let _cache = parts.checkout_derived_cache()?;
    parts.clear_cache(DerivedState::RING_FAMILIES)
}

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

fn ring_molecule() -> Molecule {
    let topology = TopologyBlock::try_from_parts(
        (0..6)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(cosmolkit_model::Element::C)))
            .collect(),
        (0..6)
            .map(|id| {
                cosmolkit_model::Bond::from_spec(
                    cosmolkit_model::BondId::new(id),
                    cosmolkit_model::BondSpec::new(
                        AtomId::new(id),
                        AtomId::new((id + 1) % 6),
                        cosmolkit_model::BondOrder::Single,
                    ),
                )
            })
            .collect(),
        Vec::new(),
        Vec::new(),
    )
    .expect("cycle topology is valid");
    Molecule::from_parts(
        topology,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .expect("test molecule is valid")
}

#[test]
#[cfg(feature = "cap-rings")]
fn ring_live_cow_ring_row_buffers_survive_install_clear_mark_finish() {
    let source = ring_molecule();
    let peer = source.clone();
    // Baselines: no live ordinary rings are installed by from_parts; source
    // and peer still share ONE derived-cache Arc.
    let baseline_valid = source.derived_cache_runtime().valid_states();
    assert!(source.derived_cache_runtime().ring_info().is_none());
    assert!(Arc::ptr_eq(
        &source.derived_cache_arc_runtime(),
        &peer.derived_cache_arc_runtime()
    ));

    let history_before = ring_live_probe::history_len();
    let output = source.with_assigned_rings().expect("ring assignment");
    let observations = ring_live_probe::history_after(history_before);
    assert_eq!(observations.len(), 1, "one actual install observation");
    let (observed_atom_ptr, observed_bond_ptr) = observations[0];

    let cache = output.derived_cache_runtime();
    let rings = cache.ring_info().expect("rings installed");
    assert!(rings.is_initialized());
    assert_eq!(rings.find_type(), cosmolkit_core::RingFindType::Fast);
    assert_eq!(rings.atom_rings().len(), 1, "one ring row");
    assert!(!rings.atom_rings().is_empty());
    assert!(!rings.bond_rings().is_empty());
    // BOTH nonempty row-buffer allocations observed at the ACTUAL install
    // site survive install -> clear(RING_FAMILIES) -> mark(RINGS) ->
    // finish with no re-clone.
    assert_eq!(
        rings.atom_rings().as_ptr() as usize,
        observed_atom_ptr,
        "atom row buffer preserved"
    );
    assert_eq!(
        rings.bond_rings().as_ptr() as usize,
        observed_bond_ptr,
        "bond row buffer preserved"
    );
    assert!(cache.valid_states().contains(DerivedState::RINGS));

    // Shared source/peer values remain unchanged: neither gains rings, and
    // the peers still share the ORIGINAL cache Arc.
    assert!(source.derived_cache_runtime().ring_info().is_none());
    assert_eq!(
        source.derived_cache_runtime().valid_states(),
        baseline_valid
    );
    assert!(peer.derived_cache_runtime().ring_info().is_none());
    assert_eq!(peer.derived_cache_runtime().valid_states(), baseline_valid);
    assert!(Arc::ptr_eq(
        &source.derived_cache_arc_runtime(),
        &peer.derived_cache_arc_runtime()
    ));
}

#[test]
fn ring_live_cow_checked_out_clear_failure_remains() {
    let mut target = ring_molecule();
    let before = target.clone();
    let error = target
        .ring_live_cow_checkout_conflict_for_test_()
        .unwrap_err();
    assert!(
        matches!(
            &error,
            OperationError::BlockCheckedOut {
                block: "derived_cache",
                ..
            }
        ),
        "got {error:?}"
    );
    // Basic in-place failure safety: the receiver stays whole and equal.
    assert_eq!(target, before);
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
