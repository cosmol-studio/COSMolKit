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
    source: Option<&'a Molecule>,
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
            source: Some(source),
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

    #[cfg(feature = "cap-reaction")]
    pub(super) fn new_reconstruction(
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        OpParts::<Access>::validate_operation_spec(spec)?;
        OpParts::<Access>::validate_semantic_preconditions(spec)?;
        OpParts::<Access>::validate_effect_contract(spec)?;
        let all = BlockSet::TOPOLOGY
            .union(BlockSet::COORDINATES)
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE);
        if spec.output != MoleculeOpOutput::Multiple
            || spec.requires_mapping != super::MappingRequirement::Reconstruction
            || spec.access.write() != all
            || spec.may_mutate != all
            || !spec.auto_remap.is_empty()
            || !spec.derived_effects.preserve.is_empty()
            || !spec.derived_effects.operation_defined.is_empty()
            || spec.derived_effects.recompute
                != super::DerivedState::VALENCE.union(super::DerivedState::RINGS)
            || spec.cip_state != super::CipStatePolicy::ReactionSourceTransition
        {
            return Err(OperationError::MappingContract {
                operation: spec.method,
                issue: "source-less reconstruction requires complete detached output and explicit effects",
                requirement: spec.requires_mapping,
            });
        }
        Ok(Self {
            spec,
            source: None,
            emitted: None,
            access: PhantomData,
            lazy_emitted: None,
            reconstruction_inputs: Vec::new(),
            reconstruction_inputs_read: false,
        })
    }

    fn source(&self) -> Result<&'a Molecule, OperationError> {
        self.source.ok_or(OperationError::IncompleteCommit {
            operation: self.spec.method,
            block: "single source molecule",
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
        Ok(self.source()?.topology())
    }

    pub(super) fn source_coordinates_runtime(&self) -> Result<&CoordinateBlock, OperationError> {
        self.ensure_read_access(BlockSet::COORDINATES, "coordinates")?;
        Ok(self.source()?.coordinate_block_runtime())
    }

    pub(super) fn source_properties_runtime(&self) -> Result<&MoleculeProperties, OperationError> {
        self.ensure_read_access(BlockSet::PROPERTIES, "properties")?;
        Ok(self.source()?.properties())
    }

    pub(super) fn source_derived_cache_runtime(
        &self,
    ) -> Result<&crate::molecule::DerivedCacheBlock, OperationError> {
        self.ensure_read_access(BlockSet::DERIVED_CACHE, "derived_cache")?;
        Ok(self.source()?.derived_cache_runtime())
    }

    #[cfg(feature = "cap-reaction")]
    pub(super) fn reconstruction_inputs_runtime<'b>(
        &mut self,
        inputs: &[&'b Molecule],
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
        let mut counts = Vec::with_capacity(inputs.len());
        let mut detached = Vec::with_capacity(inputs.len());
        for input in inputs {
            let topology = input.topology();
            let coordinates = input.coordinate_block_runtime();
            let properties = input.properties();
            let cache = input.derived_cache_runtime();
            validate_reconstruction_input(topology, coordinates, properties, cache)?;
            counts.push((topology.atoms.len(), topology.bonds.len()));
            detached.push(cosmolkit_reaction::ReactionInput {
                topology,
                coordinates,
                properties,
                rings: cache.valid_ring_info(),
                valence: cache.valence_assignment(),
            });
        }
        // Publish input evidence only after every actual input is valid.
        // Iterating the caller's slice retains repeated molecules and order;
        // no input graph or valid cache facts are cloned or recomputed.
        self.reconstruction_inputs = counts;
        self.reconstruction_inputs_read = true;
        Ok(detached)
    }

    #[cfg(feature = "cap-reaction")]
    pub(super) fn reconstruction_source_runtime(
        &mut self,
    ) -> Result<cosmolkit_reaction::ReactionInput<'a>, OperationError> {
        let source = self.source()?;
        let mut inputs = self.reconstruction_inputs_runtime(&[source])?;
        Ok(inputs.remove(0))
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
        let source = self.source()?;
        let stream = self.lazy_emitted.ok_or(OperationError::IncompleteCommit {
            operation: self.spec.method,
            block: "outputs",
        })?;
        Ok(StereoisomerIterator {
            source: source.clone(),
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
        #[cfg(feature = "cap-reaction")]
        if self.source.is_none() && !self.reconstruction_inputs_read {
            return Err(OperationError::IncompleteCommit {
                operation: self.spec.method,
                block: "reconstruction inputs",
            });
        }
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
                            if self.source.is_none() {
                                // Complete detached products have no distinguished source.
                                // Runtime alone validates and installs their cache facts;
                                // all other derived state starts invalid, as declared.
                                let mut cache = crate::molecule::DerivedCacheBlock::default();
                                cache.install_valence_assignment(product.valence);
                                let mut prepared = super::DerivedState::VALENCE;
                                if let Some(rings) = product.rings {
                                    cache.install_ring_info(rings);
                                    prepared = prepared.union(super::DerivedState::RINGS);
                                }
                                cache.mark_valid(prepared);
                                validate_reconstruction_input(
                                    &product.topology,
                                    &product.coordinates,
                                    &product.properties,
                                    &cache,
                                )?;
                                return Ok((
                                    std::sync::Arc::new(product.topology),
                                    std::sync::Arc::new(product.coordinates),
                                    std::sync::Arc::new(product.properties),
                                    std::sync::Arc::new(cache),
                                ));
                            }
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
                    self.source.ok_or(OperationError::IncompleteCommit {
                        operation: self.spec.method,
                        block: "single source molecule",
                    })?,
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

#[cfg(any(feature = "cap-reaction", test))]
fn validate_reconstruction_input(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    cache: &crate::molecule::DerivedCacheBlock,
) -> Result<(), OperationError> {
    OpParts::<()>::validate_detached_candidate_invariants(topology, coordinates, properties)?;
    cache.validate_for_topology(topology)
}

#[cfg(feature = "cap-reaction")]
fn validate_reconstruction_origins(
    spec: &'static MoleculeOpSpec,
    inputs: &[(usize, usize)],
    product: &cosmolkit_reaction::ReactionProduct,
) -> Result<(), OperationError> {
    validate_reconstruction_rows(
        spec,
        inputs,
        &product.topology,
        &product.coordinates,
        product
            .atom_origins
            .iter()
            .map(|origin| origin.as_ref().map(|o| (o.input, o.row.index()))),
        product
            .bond_origins
            .iter()
            .map(|origin| origin.as_ref().map(|o| (o.input, o.row.index()))),
        (
            product.valence.explicit_valence.len(),
            product.valence.implicit_hydrogens.len(),
        ),
        product
            .rings
            .as_ref()
            .map(|rings| (rings.atom_row_count(), rings.bond_row_count())),
    )
}

#[cfg(any(feature = "cap-reaction", test))]
fn validate_reconstruction_rows(
    spec: &'static MoleculeOpSpec,
    inputs: &[(usize, usize)],
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    mut atom_origins: impl ExactSizeIterator<Item = Option<(usize, usize)>>,
    mut bond_origins: impl ExactSizeIterator<Item = Option<(usize, usize)>>,
    valence_rows: (usize, usize),
    ring_rows: Option<(usize, usize)>,
) -> Result<(), OperationError> {
    // This is runtime evidence validation, not a chemistry algorithm.
    // Each destination is a vector row; repeated source origins are legal,
    // None means a newly constructed row, and no inverse map is fabricated.
    // The typed product adapter borrows the original origin vectors and reads
    // actual fact row counts. Iteration allocates no replacement mapping or
    // graph and leaves every source and detached candidate unchanged.
    for (actual, expected, field) in [
        (atom_origins.len(), topology.atoms.len(), "atom origins"),
        (bond_origins.len(), topology.bonds.len(), "bond origins"),
        (valence_rows.0, topology.atoms.len(), "explicit valence"),
        (valence_rows.1, topology.atoms.len(), "implicit hydrogens"),
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
    check("atom", &mut atom_origins)?;
    check("bond", &mut bond_origins)?;
    if let Some((atom_rows, bond_rows)) = ring_rows {
        for (actual, expected, field) in [
            (atom_rows, topology.atoms.len(), "ring atom membership"),
            (bond_rows, topology.bonds.len(), "ring bond membership"),
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
    topology
        .validate()
        .map_err(OperationError::InvalidTopology)?;
    coordinates
        .validate_for_atom_count(topology.atoms.len())
        .map_err(OperationError::InvalidCoordinates)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("../../tests/support/run_multiple_internal.rs");
}
