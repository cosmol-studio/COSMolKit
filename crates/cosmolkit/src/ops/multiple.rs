//! Private ordered multiple-output operation lifecycle.

use std::marker::PhantomData;

use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

use super::context::{OpParts, validate_multiple_candidate};
use super::{BlockSet, MoleculeOpOutput, MoleculeOpSpec, OperationError};
use crate::Molecule;

enum DetachedCandidate {
    Blocks(TopologyBlock, CoordinateBlock, MoleculeProperties),
    SharedCoordinates(TopologyBlock, MoleculeProperties),
    #[cfg(feature = "cap-stereoisomers")]
    StereoPrepared(
        TopologyBlock,
        Option<CoordinateBlock>,
        MoleculeProperties,
        PreparedCacheValues,
    ),
    #[cfg(feature = "cap-tautomer")]
    Prepared(TopologyBlock, MoleculeProperties, PreparedCacheValues),
    #[cfg(feature = "cap-reaction")]
    Reconstructed(cosmolkit_reaction::ReactionProduct),
}

/// Typed detached facts; construction and cache authority stay in runtime.
pub(super) struct PreparedCacheValues {
    #[cfg(any(
        feature = "cap-tautomer",
        feature = "cap-reaction",
        feature = "cap-stereoisomers"
    ))]
    pub(super) valence: cosmolkit_core::ValenceAssignment,
    #[cfg(any(
        feature = "cap-tautomer",
        feature = "cap-reaction",
        feature = "cap-stereoisomers"
    ))]
    pub(super) rings: Option<cosmolkit_core::RingInfo>,
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
    lazy_emitted:
        Option<Box<dyn Iterator<Item = Result<DetachedCandidate, OperationError>> + Send>>,
    #[cfg(feature = "cap-reaction")]
    reconstruction_inputs: Vec<(usize, usize)>,
    #[cfg(feature = "cap-reaction")]
    reconstruction_inputs_read: bool,
}

impl<'a, Access> MultiOutputOpParts<'a, Access> {
    pub(super) fn new(
        source: &'a Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        if !matches!(
            spec.output,
            MoleculeOpOutput::Multiple | MoleculeOpOutput::LazyMultiple
        ) {
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
            lazy_emitted: None,
            #[cfg(feature = "cap-reaction")]
            reconstruction_inputs: vec![(
                source.topology().atoms.len(),
                source.topology().bonds.len(),
            )],
            #[cfg(feature = "cap-reaction")]
            reconstruction_inputs_read: false,
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

    #[cfg(feature = "cap-reaction")]
    pub(super) fn reconstruction_inputs_runtime<'b>(
        &mut self,
        inputs: &'b [&'b Molecule],
    ) -> Result<Vec<cosmolkit_reaction::ReactionInput<'b>>, OperationError> {
        if self.spec.requires_mapping != super::MappingRequirement::Reconstruction
            || self.reconstruction_inputs_read
        {
            return Err(OperationError::MappingContract {
                operation: self.spec.method,
                issue: "reconstruction input set must be declared exactly once",
                requirement: self.spec.requires_mapping,
            });
        }
        self.ensure_read_access(BlockSet::TOPOLOGY, "topology")?;
        self.ensure_read_access(BlockSet::COORDINATES, "coordinates")?;
        self.ensure_read_access(BlockSet::PROPERTIES, "properties")?;
        self.ensure_read_access(BlockSet::DERIVED_CACHE, "derived_cache")?;
        self.reconstruction_inputs = inputs
            .iter()
            .map(|input| (input.topology().atoms.len(), input.topology().bonds.len()))
            .collect();
        self.reconstruction_inputs_read = true;
        Ok(inputs
            .iter()
            .map(|input| cosmolkit_reaction::ReactionInput {
                topology: input.topology(),
                coordinates: input.coordinate_block_runtime(),
                properties: input.properties(),
                rings: input.derived_cache_runtime().valid_ring_info(),
                valence: input.derived_cache_runtime().valence_assignment(),
            })
            .collect())
    }

    #[cfg(feature = "cap-reaction")]
    pub(super) fn emit_reconstructed_runtime(
        &mut self,
        products: Vec<cosmolkit_reaction::ReactionProduct>,
    ) -> Result<(), OperationError> {
        if self.spec.requires_mapping != super::MappingRequirement::Reconstruction {
            return Err(OperationError::MappingContract {
                operation: self.spec.method,
                issue: "reconstructed products require a reconstruction declaration",
                requirement: self.spec.requires_mapping,
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
            products
                .into_iter()
                .map(DetachedCandidate::Reconstructed)
                .collect(),
        );
        Ok(())
    }

    fn require_output(&self, expected: MoleculeOpOutput) -> Result<(), OperationError> {
        if self.spec.output == expected {
            Ok(())
        } else {
            Err(OperationError::OutputMismatch {
                operation: self.spec.method,
                expected,
                actual: self.spec.output,
            })
        }
    }

    pub(super) fn emit_lazy_runtime<I>(&mut self, candidates: I) -> Result<(), OperationError>
    where
        I: Iterator<
                Item = Result<
                    (TopologyBlock, Option<CoordinateBlock>, MoleculeProperties),
                    OperationError,
                >,
            > + Send
            + 'static,
    {
        self.require_output(MoleculeOpOutput::LazyMultiple)?;
        if self.lazy_emitted.is_some() {
            return Err(OperationError::OperationContract {
                operation: self.spec.method,
                field: "outputs",
                issue: "multiple-output operation emitted more than once",
                expected: 0,
                actual: 1,
            });
        }
        // STEREO-ENUM-LAZY-001: storing the stream must not call next or collect.
        // The mapping closure carries detached values only; the runtime retains
        // every validation and live construction authority.
        self.lazy_emitted = Some(Box::new(candidates.map(|candidate| {
            candidate.map(|(t, c, p)| match c {
                Some(c) => DetachedCandidate::Blocks(t, c, p),
                None => DetachedCandidate::SharedCoordinates(t, p),
            })
        })));
        Ok(())
    }

    #[cfg(feature = "cap-stereoisomers")]
    pub(super) fn emit_lazy_prepared_runtime<I>(
        &mut self,
        candidates: I,
    ) -> Result<(), OperationError>
    where
        I: Iterator<
                Item = Result<
                    (
                        TopologyBlock,
                        Option<CoordinateBlock>,
                        MoleculeProperties,
                        cosmolkit_core::ValenceAssignment,
                        cosmolkit_core::RingInfo,
                    ),
                    OperationError,
                >,
            > + Send
            + 'static,
    {
        self.require_output(MoleculeOpOutput::LazyMultiple)?;
        if !self.spec.access.can_write(BlockSet::DERIVED_CACHE) {
            return Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block: "derived_cache",
            });
        }
        if self.lazy_emitted.is_some() {
            return Err(OperationError::OperationContract {
                operation: self.spec.method,
                field: "outputs",
                issue: "multiple-output operation emitted more than once",
                expected: 0,
                actual: 1,
            });
        }
        self.lazy_emitted = Some(Box::new(candidates.map(|candidate| {
            candidate.map(|(t, c, p, valence, rings)| {
                DetachedCandidate::StereoPrepared(
                    t,
                    c,
                    p,
                    PreparedCacheValues {
                        valence,
                        rings: Some(rings),
                    },
                )
            })
        })));
        Ok(())
    }

    pub(super) fn finish_lazy(self) -> Result<StereoisomerIterator, OperationError> {
        self.require_output(MoleculeOpOutput::LazyMultiple)?;
        let stream = self.lazy_emitted.ok_or(OperationError::IncompleteCommit {
            operation: self.spec.method,
            block: "outputs",
        })?;
        Ok(StereoisomerIterator {
            source: self.source.clone(),
            spec: self.spec,
            stream: Some(stream),
            yielded: 0,
        })
    }

    pub(super) fn emit_all_runtime(
        &mut self,
        candidates: Vec<(TopologyBlock, CoordinateBlock, MoleculeProperties)>,
    ) -> Result<(), OperationError> {
        self.require_output(MoleculeOpOutput::Multiple)?;
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
                        PreparedCacheValues {
                            valence,
                            rings: Some(rings),
                        },
                    )
                })
                .collect(),
        );
        Ok(())
    }

    fn finish_runtime(self) -> Result<Vec<Molecule>, OperationError> {
        self.require_output(MoleculeOpOutput::Multiple)?;
        let candidates = self.emitted.ok_or(OperationError::IncompleteCommit {
            operation: self.spec.method,
            block: "outputs",
        })?;

        let validated = candidates
            .into_iter()
            .map(|candidate| {
                let (topology, coordinates, properties, prepared, reconstruction_validated) =
                    match candidate {
                        DetachedCandidate::Blocks(t, c, p) => (t, Some(c), p, None, false),
                        DetachedCandidate::SharedCoordinates(t, p) => (t, None, p, None, false),
                        #[cfg(feature = "cap-stereoisomers")]
                        DetachedCandidate::StereoPrepared(t, c, p, facts) => {
                            (t, c, p, Some(facts), false)
                        }
                        #[cfg(feature = "cap-tautomer")]
                        DetachedCandidate::Prepared(t, p, facts) => {
                            (t, None, p, Some(facts), false)
                        }
                        #[cfg(feature = "cap-reaction")]
                        DetachedCandidate::Reconstructed(product) => {
                            validate_reconstruction_origins(
                                self.spec,
                                &self.reconstruction_inputs,
                                &product,
                            )?;
                            (
                                product.topology,
                                Some(product.coordinates),
                                product.properties,
                                Some(PreparedCacheValues {
                                    valence: product.valence,
                                    rings: product.rings,
                                }),
                                true,
                            )
                        }
                    };
                validate_multiple_candidate(
                    self.source,
                    self.spec,
                    topology,
                    coordinates,
                    properties,
                    prepared,
                    reconstruction_validated,
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

pub struct StereoisomerIterator {
    source: Molecule,
    spec: &'static MoleculeOpSpec,
    stream: Option<Box<dyn Iterator<Item = Result<DetachedCandidate, OperationError>> + Send>>,
    yielded: usize,
}

impl StereoisomerIterator {
    #[must_use]
    pub fn yielded_count(&self) -> usize {
        self.yielded
    }
}

impl Iterator for StereoisomerIterator {
    type Item = Result<Molecule, OperationError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Taking the stream before calling foreign/domain code also fuses a
        // panic if the caller catches it. Errors and exhaustion release it.
        let mut stream = self.stream.take()?;
        let candidate = match stream.next()? {
            Ok(candidate) => candidate,
            Err(error) => return Some(Err(error)),
        };
        let (topology, coordinates, properties, prepared) = match candidate {
            DetachedCandidate::Blocks(t, c, p) => (t, Some(c), p, None),
            DetachedCandidate::SharedCoordinates(t, p) => (t, None, p, None),
            #[cfg(feature = "cap-stereoisomers")]
            DetachedCandidate::StereoPrepared(t, c, p, facts) => (t, c, p, Some(facts)),
            #[cfg(feature = "cap-tautomer")]
            DetachedCandidate::Prepared(t, p, facts) => (t, None, p, Some(facts)),
            #[cfg(feature = "cap-reaction")]
            DetachedCandidate::Reconstructed(_) => {
                return Some(Err(OperationError::MappingContract {
                    operation: self.spec.method,
                    issue: "lazy output cannot carry reconstruction candidates",
                    requirement: self.spec.requires_mapping,
                }));
            }
        };
        // The same validator used by eager multiple outputs checks each
        // requested tuple. This path adds no operation-specific chemistry.
        let result = validate_multiple_candidate(
            &self.source,
            self.spec,
            topology,
            coordinates,
            properties,
            prepared,
            false,
        )
        .and_then(|(t, c, p, cache)| Molecule::from_runtime_parts(t, c, p, cache));
        match result {
            Err(error) => Some(Err(error)),
            Ok(molecule) => {
                let Some(yielded) = self.yielded.checked_add(1) else {
                    return Some(Err(OperationError::OperationContract {
                        operation: self.spec.method,
                        field: "yielded_count",
                        issue: "validated output counter exhausted",
                        expected: 0,
                        actual: 1,
                    }));
                };
                self.yielded = yielded;
                self.stream = Some(stream);
                Some(Ok(molecule))
            }
        }
    }
}
impl std::iter::FusedIterator for StereoisomerIterator {}

#[cfg(feature = "cap-reaction")]
fn validate_reconstruction_origins(
    spec: &'static MoleculeOpSpec,
    inputs: &[(usize, usize)],
    product: &cosmolkit_reaction::ReactionProduct,
) -> Result<(), OperationError> {
    // This is runtime evidence validation, not a chemistry algorithm.
    // Each destination is a vector row; repeated source origins are legal,
    // None means a newly constructed row, and no inverse map is fabricated.
    for (actual, expected, field) in [
        (
            product.atom_origins.len(),
            product.topology.atoms.len(),
            "atom origins",
        ),
        (
            product.bond_origins.len(),
            product.topology.bonds.len(),
            "bond origins",
        ),
        (
            product.valence.explicit_valence.len(),
            product.topology.atoms.len(),
            "explicit valence",
        ),
        (
            product.valence.implicit_hydrogens.len(),
            product.topology.atoms.len(),
            "implicit hydrogens",
        ),
    ] {
        if actual != expected {
            return Err(OperationError::InvalidAlgorithmResult {
                operation: spec.method,
                field,
                actual,
                expected,
            });
        }
    }
    let check = |entity: &'static str,
                 origins: &mut dyn Iterator<Item = Option<(usize, usize)>>| {
        for (destination, origin) in origins.enumerate() {
            if let Some((input, row)) = origin {
                let row_count = inputs
                    .get(input)
                    .map(|counts| if entity == "atom" { counts.0 } else { counts.1 });
                if row_count.is_none_or(|count| row >= count) {
                    return Err(OperationError::InvalidReconstructionOrigin {
                        operation: spec.method,
                        entity,
                        destination,
                        input,
                        row,
                        input_count: inputs.len(),
                        row_count,
                    });
                }
            }
        }
        Ok(())
    };
    check(
        "atom",
        &mut product
            .atom_origins
            .iter()
            .map(|origin| origin.as_ref().map(|o| (o.input, o.row.index()))),
    )?;
    check(
        "bond",
        &mut product
            .bond_origins
            .iter()
            .map(|origin| origin.as_ref().map(|o| (o.input, o.row.index()))),
    )?;
    if let Some(rings) = &product.rings {
        for (actual, expected, field) in [
            (
                rings.atom_row_count(),
                product.topology.atoms.len(),
                "ring atom membership",
            ),
            (
                rings.bond_row_count(),
                product.topology.bonds.len(),
                "ring bond membership",
            ),
        ] {
            if actual != expected {
                return Err(OperationError::InvalidAlgorithmResult {
                    operation: spec.method,
                    field,
                    actual,
                    expected,
                });
            }
        }
    }
    product
        .topology
        .validate()
        .map_err(OperationError::InvalidTopology)?;
    product
        .coordinates
        .validate_for_atom_count(product.topology.atoms.len())
        .map_err(OperationError::InvalidCoordinates)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("../../tests/support/run_multiple_internal.rs");
}
