//! Thin valence operation projection over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{DerivedState, Molecule, PreservationProof, ValenceError, ValenceParams};

#[mol_op_body(with_assigned_valence, parts)]
pub(crate) fn assign_valence_impl(params: &ValenceParams) -> Result<(), OperationError> {
    let assignment = cosmolkit_core::assign_valence(parts.topology()?, params)
        .map_err(OperationError::Valence)?;
    let mut cache = parts.checkout_derived_cache()?;
    cache.install_valence_assignment(assignment);
    parts.install_derived_cache(cache)?;
    parts.mark_cache_updated(DerivedState::VALENCE)?;
    parts.prove_preserved(
        DerivedState::RINGS.union(DerivedState::STEREO),
        PreservationProof::UnchangedInput,
    )?;
    parts.apply_cip_policy()
}

#[cfg(test)]
mod metadata_tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn atom_metadata_queries_preserve_both_receivers_cache_identity_and_values() {
        let source = Molecule::from_parts(
            crate::TopologyBlock::try_from_parts(
                vec![crate::Atom::from_spec(
                    crate::AtomId::new(0),
                    crate::AtomSpec::new(crate::Element::C),
                )],
                Vec::new(),
                Vec::new(),
                Vec::new(),
            )
            .unwrap(),
            crate::CoordinateBlock::default(),
            crate::MoleculeProperties::default(),
        )
        .unwrap();
        for molecule in [source.clone(), source.with_assigned_valence().unwrap()] {
            let peer = molecule.clone();
            let cache = molecule.derived_cache_arc_runtime();
            let values = (*cache).clone();
            for recalculate in [true, false, true, false] {
                assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
                assert!(Arc::ptr_eq(&cache, &peer.derived_cache_arc_runtime()));
                let result = molecule.atom_metadata(recalculate);
                assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
                assert!(Arc::ptr_eq(&cache, &peer.derived_cache_arc_runtime()));
                assert_eq!(molecule.derived_cache_runtime(), &values);
                assert_eq!(peer.derived_cache_runtime(), &values);
                assert!(std::ptr::eq(molecule.topology(), peer.topology()));
                assert!(Arc::ptr_eq(
                    &molecule.coordinates_arc_runtime(),
                    &peer.coordinates_arc_runtime()
                ));
                assert!(std::ptr::eq(molecule.properties(), peer.properties()));
                if recalculate || cache.valence_assignment().is_some() {
                    assert_eq!(result.unwrap()[0].total_valence, 4);
                } else {
                    assert_eq!(
                        result,
                        Err(ValenceError::ExplicitValenceCacheNotInitialized {
                            atom: crate::AtomId::new(0)
                        })
                    );
                }
            }
        }
    }
}

impl Molecule {
    /// Reports whether one atom violates the detached owner's valence rules.
    pub fn has_valence_violation(&self, atom_id: crate::AtomId) -> Result<bool, ValenceError> {
        cosmolkit_core::has_valence_violation(self.topology(), atom_id)
    }
}

impl Molecule {
    /// Read degree, valence and hydrogen metadata without installing cache state.
    ///
    /// With `recalculate = true` (the Python default), calculate strict valence
    /// from the current topology. With `false`, read only the existing valid
    /// valence assignment; missing or invalidated cache entries return a typed
    /// error, never trigger calculation. Neither branch changes molecule state.
    pub fn atom_metadata(
        &self,
        recalculate: bool,
    ) -> Result<Vec<crate::AtomMetadata>, ValenceError> {
        if recalculate {
            cosmolkit_core::atom_metadata(self.topology())
        } else {
            cosmolkit_core::atom_metadata_from_assignment(
                self.topology(),
                self.derived_cache_runtime().valence_assignment(),
            )
        }
    }
}
