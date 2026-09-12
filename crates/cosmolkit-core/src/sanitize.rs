//! Detached sanitization assignments.
//!
//! This module contains the value-level portion of RDKit's property-cache
//! stage.  It deliberately stops before cache installation: cache validity,
//! invalidation, and operation sequencing belong to the `cosmolkit` runtime.

use cosmolkit_model::TopologyBlock;

use crate::{
    ValenceAssignment, ValenceError, ValenceModel, assign_valence_with_options_for_topology,
};

/// Errors produced while evaluating detached sanitization properties.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SanitizeError {
    #[error("detached topology is invalid: {0}")]
    InvalidTopology(String),
    #[error("property-cache valence assignment failed: {0}")]
    Valence(#[from] ValenceError),
}

/// The value-level result of the RDKit property-cache sanitization stage.
///
/// This is intentionally an assignment rather than a runtime cache object.
/// The parent runtime decides whether and when these values become
/// authoritative derived state on a live molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PropertyCacheAssignment {
    pub valence: ValenceAssignment,
}

impl PropertyCacheAssignment {
    #[must_use]
    pub fn valence(&self) -> &ValenceAssignment {
        &self.valence
    }

    #[must_use]
    pub fn into_valence(self) -> ValenceAssignment {
        self.valence
    }
}

/// Assign RDKit-like explicit valence and implicit hydrogen facts for a
/// detached topology.
///
/// // RDKit❗✔️: if (atom->needsUpdatePropertyCache()) {
/// // RDKit❗✔️:   atom->updatePropertyCache(false);
/// // RDKit❗✔️: }
///
/// The source operation writes atom property-cache fields in-place.  At this
/// boundary we return the complete assignment as an owned value; installation
/// into a live molecule remains the runtime's responsibility.
pub fn assign_property_cache_for_topology(
    topology: &TopologyBlock,
    strict: bool,
) -> Result<PropertyCacheAssignment, SanitizeError> {
    topology
        .validate()
        .map_err(|error| SanitizeError::InvalidTopology(error.to_string()))?;
    let valence =
        assign_valence_with_options_for_topology(topology, ValenceModel::RdkitLike, strict)?;
    Ok(PropertyCacheAssignment { valence })
}

/// Convenience accessor used by callers that only need the assignment.
pub fn assign_valence_properties_for_topology(
    topology: &TopologyBlock,
    strict: bool,
) -> Result<ValenceAssignment, SanitizeError> {
    assign_property_cache_for_topology(topology, strict).map(PropertyCacheAssignment::into_valence)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec};
    use cosmolkit_types::{BondOrder, Element};

    fn ethanol_topology() -> TopologyBlock {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            ),
        ];
        let adjacency = AdjacencyList::from_topology(3, &bonds);
        TopologyBlock {
            atoms,
            bonds,
            adjacency,
            stereo_groups: Vec::new(),
            substance_groups: Vec::new(),
        }
    }

    #[test]
    fn detached_assignment_matches_formal_valence_boundary() {
        let topology = ethanol_topology();
        let assignment = assign_property_cache_for_topology(&topology, false).unwrap();
        assert_eq!(assignment.valence.explicit_valence, vec![1, 2, 1]);
        assert_eq!(assignment.valence.implicit_hydrogens, vec![3, 2, 1]);
    }
}
