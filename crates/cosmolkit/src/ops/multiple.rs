//! Private ordered multiple-output operation lifecycle.

use std::marker::PhantomData;

use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

use super::context::{OpParts, validate_multiple_candidate};
use super::{BlockSet, MoleculeOpOutput, MoleculeOpSpec, OperationError};
use crate::Molecule;

enum DetachedCandidate {
    Blocks(TopologyBlock, CoordinateBlock, MoleculeProperties),
    #[cfg(feature = "cap-tautomer")]
    Prepared(TopologyBlock, MoleculeProperties, PreparedCacheValues),
}

/// Typed detached facts; construction and cache authority stay in runtime.
pub(super) struct PreparedCacheValues {
    #[cfg(feature = "cap-tautomer")]
    pub(super) valence: cosmolkit_core::ValenceAssignment,
    #[cfg(feature = "cap-tautomer")]
    pub(super) rings: cosmolkit_core::RingInfo,
}

/// One private collection transaction for a generated multiple-output body.
///
/// The body can only emit detached tuples. No candidate becomes a public
/// molecule until every tuple has passed the shared runtime validators.
pub(crate) struct MultiOutputOpParts<'a, Access> {
    spec: &'static MoleculeOpSpec,
    source: &'a Molecule,
    emitted: Option<Vec<DetachedCandidate>>,
    access: PhantomData<Access>,
}

impl<'a, Access> MultiOutputOpParts<'a, Access> {
    pub(super) fn new(
        source: &'a Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        if spec.output != MoleculeOpOutput::Multiple {
            return Err(OperationError::OutputMismatch {
                operation: spec.method,
                expected: MoleculeOpOutput::Multiple,
                actual: spec.output,
            });
        }
        OpParts::<Access>::validate_semantic_preconditions(spec)?;
        Ok(Self {
            spec,
            source,
            emitted: None,
            access: PhantomData,
        })
    }

    fn ensure_read_access(
        &self,
        block: BlockSet,
        name: &'static str,
    ) -> Result<(), OperationError> {
        if self.spec.access.can_read(block) {
            Ok(())
        } else {
            Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block: name,
            })
        }
    }

    pub(super) fn source_topology_runtime(&self) -> Result<&TopologyBlock, OperationError> {
        self.ensure_read_access(BlockSet::TOPOLOGY, "topology")?;
        Ok(self.source.topology())
    }

    pub(super) fn source_coordinates_runtime(&self) -> Result<&CoordinateBlock, OperationError> {
        self.ensure_read_access(BlockSet::COORDINATES, "coordinates")?;
        Ok(self.source.coordinate_block_runtime())
    }

    pub(super) fn source_properties_runtime(&self) -> Result<&MoleculeProperties, OperationError> {
        self.ensure_read_access(BlockSet::PROPERTIES, "properties")?;
        Ok(self.source.properties())
    }

    pub(super) fn source_derived_cache_runtime(
        &self,
    ) -> Result<&crate::molecule::DerivedCacheBlock, OperationError> {
        self.ensure_read_access(BlockSet::DERIVED_CACHE, "derived_cache")?;
        Ok(self.source.derived_cache_runtime())
    }

    pub(super) fn emit_all_runtime(
        &mut self,
        candidates: Vec<(TopologyBlock, CoordinateBlock, MoleculeProperties)>,
    ) -> Result<(), OperationError> {
        if self.emitted.is_some() {
            return Err(OperationError::OperationContract {
                operation: self.spec.method,
                field: "outputs",
                issue: "multiple-output operation emitted more than once",
                expected: 0,
                actual: 1,
            });
        }
        self.emitted = Some(
            candidates
                .into_iter()
                .map(|(t, c, p)| DetachedCandidate::Blocks(t, c, p))
                .collect(),
        );
        Ok(())
    }

    #[cfg(feature = "cap-tautomer")]
    pub(super) fn emit_prepared_runtime(
        &mut self,
        candidates: Vec<(
            TopologyBlock,
            MoleculeProperties,
            cosmolkit_core::ValenceAssignment,
            cosmolkit_core::RingInfo,
        )>,
    ) -> Result<(), OperationError> {
        if !self.spec.access.can_write(BlockSet::DERIVED_CACHE) {
            return Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block: "derived_cache",
            });
        }
        if self.emitted.is_some() {
            return Err(OperationError::OperationContract {
                operation: self.spec.method,
                field: "outputs",
                issue: "multiple-output operation emitted more than once",
                expected: 0,
                actual: 1,
            });
        }
        self.emitted = Some(
            candidates
                .into_iter()
                .map(|(topology, properties, valence, rings)| {
                    DetachedCandidate::Prepared(
                        topology,
                        properties,
                        PreparedCacheValues { valence, rings },
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn finish_runtime(self) -> Result<Vec<Molecule>, OperationError> {
        let candidates = self.emitted.ok_or(OperationError::IncompleteCommit {
            operation: self.spec.method,
            block: "outputs",
        })?;

        let validated = candidates
            .into_iter()
            .map(|candidate| {
                let (topology, coordinates, properties, prepared) = match candidate {
                    DetachedCandidate::Blocks(t, c, p) => (t, Some(c), p, None),
                    #[cfg(feature = "cap-tautomer")]
                    DetachedCandidate::Prepared(t, p, facts) => (t, None, p, Some(facts)),
                };
                validate_multiple_candidate(
                    self.source,
                    self.spec,
                    topology,
                    coordinates,
                    properties,
                    prepared,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;

        validated
            .into_iter()
            .map(|(topology, coordinates, properties, cache)| {
                Molecule::from_runtime_parts(topology, coordinates, properties, cache)
            })
            .collect()
    }

    pub(super) fn finish(self) -> Result<Vec<Molecule>, OperationError> {
        self.finish_runtime()
    }
}

#[cfg(test)]
mod tests {
    include!("../../tests/support/run_multiple_internal.rs");
}
