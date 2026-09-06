//! Detached sanitization assignments.
//!
//! This module contains the value-level portion of sanitization.  It computes
//! assignments from detached model blocks and never installs derived state on
//! a live molecule.  The operation layer remains responsible for cache
//! invalidation, effect bookkeeping, and commit.

use cosmolkit_model::TopologyBlock;

/// Recomputes the RDKit-like property-cache assignment for detached topology
/// data.
///
/// The returned [`crate::ValenceAssignment`] is an ordinary value result.  It
/// is deliberately not written into a runtime cache here; callers at the
/// operation boundary decide whether the assignment is accepted and installed.
pub fn assign_property_cache_for_topology(
    topology: &TopologyBlock,
    strict: bool,
) -> Result<crate::ValenceAssignment, crate::ValenceError> {
    // Property-cache assignment is the valence/implicit-hydrogen calculation
    // performed by RDKit's updatePropertyCache() path.  Keep the source-backed
    // implementation in valence.rs so there is one algorithm owner.
    crate::valence::assign_valence_with_options_for_topology(
        topology,
        crate::ValenceModel::RdkitLike,
        strict,
    )
}
