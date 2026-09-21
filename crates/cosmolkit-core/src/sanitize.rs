//! RDKit-aligned detached sanitization orchestration and assignments.
//!
//! This module sequences source-backed algorithms over detached topology
//! values. It deliberately stops before live cache installation: cache
//! validity, invalidation, operation capabilities, and commit authority belong
//! exclusively to the parent `cosmolkit` runtime.

use std::ops::{BitAnd, BitOr, BitOrAssign};

use cosmolkit_model::{
    AtomId, QueryStateError, QueryStateRef, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::Hybridization;

use crate::{
    AromaticityError, AromaticityParams, AtropisomerError, CleanupError, CleanupParams,
    ConjugationError, HybridizationAssignment, HybridizationError, KekulizeError, KekulizeParams,
    RadicalError, RingFindingError, RingInfo, RingSearchParams, StereoError, ValenceAssignment,
    ValenceError, ValenceModel, ValenceParams, assign_aromaticity_with_query_state,
    assign_conjugation, assign_hybridization, assign_radicals, assign_valence,
    assign_valence_state_for_atom_from_parts, cleanup, find_sssr, kekulize,
    kekulize_with_query_state, symmetrized_sssr,
};

use crate::hcount::{AdjustHsError, adjust_hs};

const NAMED_SANITIZE_BITS: u32 = 0x0000_0fff;
const ALL_SANITIZE_BITS: u32 = 0x0fff_ffff;

/// Independently selectable RDKit sanitization stages.
///
/// `ALL` intentionally retains RDKit's reserved high bits. Raw construction
/// accepts either the exact `ALL` sentinel or a combination of named bits;
/// arbitrary unknown-bit patterns are rejected instead of silently doing
/// nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SanitizeOperations(u32);

impl SanitizeOperations {
    pub const NONE: Self = Self(0x000);
    pub const CLEANUP: Self = Self(0x001);
    pub const PROPERTIES: Self = Self(0x002);
    pub const SYMM_RINGS: Self = Self(0x004);
    pub const KEKULIZE: Self = Self(0x008);
    pub const FIND_RADICALS: Self = Self(0x010);
    pub const SET_AROMATICITY: Self = Self(0x020);
    pub const SET_CONJUGATION: Self = Self(0x040);
    pub const SET_HYBRIDIZATION: Self = Self(0x080);
    pub const CLEANUP_CHIRALITY: Self = Self(0x100);
    pub const ADJUST_HS: Self = Self(0x200);
    pub const CLEANUP_ORGANOMETALLICS: Self = Self(0x400);
    pub const CLEANUP_ATROPISOMERS: Self = Self(0x800);
    pub const ALL: Self = Self(ALL_SANITIZE_BITS);

    pub fn from_bits(bits: u32) -> Result<Self, SanitizeError> {
        let unknown_bits = bits & !NAMED_SANITIZE_BITS;
        if unknown_bits != 0 && bits != ALL_SANITIZE_BITS {
            return Err(SanitizeError::InvalidOperations { bits, unknown_bits });
        }
        Ok(Self(bits))
    }

    #[must_use]
    pub const fn bits(self) -> u32 {
        self.0
    }

    #[must_use]
    pub const fn contains(self, operation: Self) -> bool {
        operation.0 != 0 && self.0 & operation.0 == operation.0
    }

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }
}

impl Default for SanitizeOperations {
    fn default() -> Self {
        Self::ALL
    }
}

impl TryFrom<u32> for SanitizeOperations {
    type Error = SanitizeError;

    fn try_from(bits: u32) -> Result<Self, Self::Error> {
        Self::from_bits(bits)
    }
}

impl From<SanitizeOperations> for u32 {
    fn from(operations: SanitizeOperations) -> Self {
        operations.bits()
    }
}

impl BitOr for SanitizeOperations {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for SanitizeOperations {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

impl BitAnd for SanitizeOperations {
    type Output = Self;

    fn bitand(self, rhs: Self) -> Self::Output {
        Self(self.0 & rhs.0)
    }
}

/// Exact sanitization stage used for source-compatible failure reporting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum SanitizeStage {
    None = 0x000,
    Cleanup = 0x001,
    Properties = 0x002,
    SymmRings = 0x004,
    Kekulize = 0x008,
    FindRadicals = 0x010,
    SetAromaticity = 0x020,
    SetConjugation = 0x040,
    SetHybridization = 0x080,
    CleanupChirality = 0x100,
    AdjustHs = 0x200,
    CleanupOrganometallics = 0x400,
    CleanupAtropisomers = 0x800,
}

/// Parameters for the detached sanitization pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SanitizeParams {
    pub operations: SanitizeOperations,
}

impl Default for SanitizeParams {
    fn default() -> Self {
        Self {
            operations: SanitizeOperations::ALL,
        }
    }
}

/// The authoritative detached output of a successful sanitize operation.
#[derive(Debug, Clone, PartialEq)]
pub struct SanitizeAssignment {
    pub topology: TopologyBlock,
}

/// A sanitization failure with the exact active source stage and typed cause.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum SanitizeError {
    #[error("invalid sanitize operation bits 0x{bits:08x}; unknown bits 0x{unknown_bits:08x}")]
    InvalidOperations { bits: u32, unknown_bits: u32 },
    #[error("sanitization failed at {stage:?}: invalid topology: {source}")]
    InvalidTopology {
        stage: SanitizeStage,
        source: TopologyValidationError,
    },
    #[error("sanitization failed at {stage:?}: invalid query state: {source}")]
    InvalidQueryState {
        stage: SanitizeStage,
        source: QueryStateError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Cleanup {
        stage: SanitizeStage,
        source: CleanupError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Properties {
        stage: SanitizeStage,
        source: PropertyCacheError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Rings {
        stage: SanitizeStage,
        source: RingFindingError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Kekulize {
        stage: SanitizeStage,
        source: KekulizeError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Radicals {
        stage: SanitizeStage,
        source: RadicalError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Aromaticity {
        stage: SanitizeStage,
        source: AromaticityError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Conjugation {
        stage: SanitizeStage,
        source: ConjugationError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Hybridization {
        stage: SanitizeStage,
        source: HybridizationError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Atropisomers {
        stage: SanitizeStage,
        source: AtropisomerError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Chirality {
        stage: SanitizeStage,
        source: StereoError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    AdjustHs {
        stage: SanitizeStage,
        source: AdjustHsError,
    },
}

/// One source-equivalent chemistry problem found without stopping detection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChemistryProblem {
    pub operation: SanitizeStage,
    pub error: ChemistryProblemError,
}

/// Problems that RDKit's detector catches and returns to the caller.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum ChemistryProblemError {
    #[error(transparent)]
    Valence(ValenceError),
    #[error(transparent)]
    Kekulize(KekulizeError),
}

/// Ordered chemistry problems found while inspecting a detached topology.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ChemistryProblemReport {
    pub problems: Vec<ChemistryProblem>,
}

/// Configuration for detached property-cache calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PropertyCacheParams {
    pub strict: bool,
}

impl Default for PropertyCacheParams {
    fn default() -> Self {
        Self { strict: true }
    }
}

/// Errors produced while calculating detached property-cache values.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PropertyCacheError {
    #[error("property-cache valence assignment failed: {0}")]
    Valence(#[from] ValenceError),
}

/// The complete value-level result of the RDKit property-cache stage.
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

/// Calculate RDKit-like explicit valence and implicit-hydrogen facts for a
/// detached topology.
///
/// The source operation writes existing atom cache fields in place. This
/// boundary instead returns two newly allocated vectors so that no algorithm
/// crate can acquire live-cache authority. That preserves source behavior but
/// is a material allocation difference from the source implementation.
pub fn assign_property_cache(
    topology: &TopologyBlock,
    params: &PropertyCacheParams,
) -> Result<PropertyCacheAssignment, PropertyCacheError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::updatePropertyCache
    // RDKit✔️❌: void ROMol::updatePropertyCache(bool strict) {
    // RDKit✔️❌:   for (auto atom : atoms()) {
    // RDKit✔️❌:     atom->updatePropertyCache(strict);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto bond : bonds()) {
    // RDKit✔️❌:     bond->updatePropertyCache(strict);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION ROMol::updatePropertyCache

    // BEGIN RDKIT CPP FUNCTION Atom::updatePropertyCache
    // RDKit✔️❌: void Atom::updatePropertyCache(bool strict) {
    // RDKit✔️❌:   calcExplicitValence(strict);
    // RDKit✔️❌:   calcImplicitValence(strict);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Atom::updatePropertyCache

    // BEGIN RDKIT CPP FUNCTION Bond::updatePropertyCache
    // RDKit✔️❌: void updatePropertyCache(bool strict = true) { (void)strict; }
    // END RDKIT CPP FUNCTION Bond::updatePropertyCache

    let valence = assign_valence(
        topology,
        &ValenceParams {
            model: ValenceModel::RdkitLike,
            strict: params.strict,
        },
    )?;
    Ok(PropertyCacheAssignment { valence })
}

fn clear_topology_computed_properties(topology: &mut TopologyBlock) {
    for atom in &mut topology.atoms {
        atom.clear_computed_props();
    }
    for bond in &mut topology.bonds {
        bond.clear_computed_props();
    }
}

fn materialize_radicals(topology: &mut TopologyBlock, values: &[u8]) {
    for (atom, &value) in topology.atoms.iter_mut().zip(values) {
        atom.set_radical_electrons(value);
    }
}

fn materialize_hybridization(topology: &mut TopologyBlock, assignment: &HybridizationAssignment) {
    for (atom, &value) in topology.atoms.iter_mut().zip(&assignment.values) {
        atom.set_hybridization(value);
    }
}

fn current_hybridization(topology: &TopologyBlock) -> HybridizationAssignment {
    HybridizationAssignment {
        values: topology
            .atoms
            .iter()
            .map(cosmolkit_model::Atom::hybridization)
            .collect(),
    }
}

/// Run RDKit's selected sanitization stages in exact source order over a
/// detached topology value.
///
/// The borrowed input is never mutated and a failed stage returns no partial
/// topology. Intermediate assignments reproduce stage dependencies but do not
/// become a second live-cache authority.
pub fn sanitize_topology(
    topology: &TopologyBlock,
    params: &SanitizeParams,
) -> Result<SanitizeAssignment, SanitizeError> {
    sanitize_topology_with_query_state(topology, params, None)
}

#[doc(hidden)]
pub fn sanitize_topology_with_query_state(
    topology: &TopologyBlock,
    params: &SanitizeParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<SanitizeAssignment, SanitizeError> {
    topology
        .validate()
        .map_err(|source| SanitizeError::InvalidTopology {
            stage: SanitizeStage::None,
            source,
        })?;
    if let Some(state) = query_state {
        QueryStateRef::try_for_topology(state.atoms(), state.bonds(), topology).map_err(
            |source| SanitizeError::InvalidQueryState {
                stage: SanitizeStage::None,
                source,
            },
        )?;
    }

    // Complete pinned source: MolOps.cpp::sanitizeMol(RWMol &, unsigned int &, unsigned int).
    // RDKit✔️❌: void sanitizeMol(RWMol &mol, unsigned int &operationThatFailed,
    // RDKit✔️❌:                  unsigned int sanitizeOps) {
    // RDKit✔️❌:   // clear out any cached properties
    // RDKit✔️❌:   mol.clearComputedProps();
    // RDKit✔️❌:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_CLEANUP;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     // clean up things like nitro groups
    // RDKit✔️❌:     cleanUp(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // fix things like non-metal to metal bonds that should be dative.
    // RDKit✔️❌:   operationThatFailed = SANITIZE_CLEANUP_ORGANOMETALLICS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     cleanUpOrganometallics(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // update computed properties on atoms and bonds:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     mol.updatePropertyCache(true);
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     mol.updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_SYMMRINGS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     VECT_INT_VECT arings;
    // RDKit✔️❌:     MolOps::symmetrizeSSSR(mol, arings);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // kekulizations
    // RDKit✔️❌:   operationThatFailed = SANITIZE_KEKULIZE;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     Kekulize(mol, true, false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // look for radicals:
    // RDKit✔️❌:   // We do this now because we need to know
    // RDKit✔️❌:   // that the N in [N]1C=CC=C1 has a radical
    // RDKit✔️❌:   // before we move into setAromaticity().
    // RDKit✔️❌:   // It's important that this happen post-Kekulization
    // RDKit✔️❌:   // because there's no way of telling what to do
    // RDKit✔️❌:   // with the same molecule if it's in the form
    // RDKit✔️❌:   // [n]1cccc1
    // RDKit✔️❌:   operationThatFailed = SANITIZE_FINDRADICALS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     assignRadicals(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // then do aromaticity perception
    // RDKit✔️❌:   operationThatFailed = SANITIZE_SETAROMATICITY;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     setAromaticity(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // set conjugation
    // RDKit✔️❌:   operationThatFailed = SANITIZE_SETCONJUGATION;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     setConjugation(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // set hybridization
    // RDKit✔️❌:   operationThatFailed = SANITIZE_SETHYBRIDIZATION;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     setHybridization(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_CLEANUPATROPISOMERS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     cleanupAtropisomers(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // remove bogus chirality specs:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_CLEANUPCHIRALITY;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     cleanupChirality(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // adjust Hydrogen counts:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_ADJUSTHS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     adjustHs(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // now that everything has been cleaned up, go through and check/update the
    // RDKit✔️❌:   // computed valences on atoms and bonds one more time
    // RDKit✔️❌:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     mol.updatePropertyCache(true);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   operationThatFailed = 0;
    // RDKit✔️❌: }
    // The source mutates one `RWMol`. This detached boundary clones the full
    // topology and several owner stages return further owned clones, which is
    // materially more allocation despite preserving the source stage order
    // and each owner's asymptotic traversal shape.

    let operations = params.operations;
    let mut working = topology.clone();
    clear_topology_computed_properties(&mut working);

    if operations.contains(SanitizeOperations::CLEANUP) {
        working = cleanup(
            &working,
            &CleanupParams {
                charge_normalization: true,
                organometallics: false,
            },
        )
        .map_err(|source| SanitizeError::Cleanup {
            stage: SanitizeStage::Cleanup,
            source,
        })?;
    }

    if operations.contains(SanitizeOperations::CLEANUP_ORGANOMETALLICS) {
        working = cleanup(
            &working,
            &CleanupParams {
                charge_normalization: false,
                organometallics: true,
            },
        )
        .map_err(|source| SanitizeError::Cleanup {
            stage: SanitizeStage::CleanupOrganometallics,
            source,
        })?;
    }

    let valence = assign_property_cache(
        &working,
        &PropertyCacheParams {
            strict: operations.contains(SanitizeOperations::PROPERTIES),
        },
    )
    .map_err(|source| SanitizeError::Properties {
        stage: SanitizeStage::Properties,
        source,
    })?
    .into_valence();

    let mut rings: Option<RingInfo> = None;
    if operations.contains(SanitizeOperations::SYMM_RINGS) {
        rings = Some(
            symmetrized_sssr(&working, &RingSearchParams::default()).map_err(|source| {
                SanitizeError::Rings {
                    stage: SanitizeStage::SymmRings,
                    source,
                }
            })?,
        );
    }

    if operations.contains(SanitizeOperations::KEKULIZE) {
        working = kekulize_with_query_state(
            &working,
            &KekulizeParams {
                mark_atoms_bonds: true,
                canonical: false,
                max_backtracks: KekulizeParams::default().max_backtracks,
            },
            query_state,
        )
        .map_err(|source| SanitizeError::Kekulize {
            stage: SanitizeStage::Kekulize,
            source,
        })?
        .topology;
    }

    if operations.contains(SanitizeOperations::FIND_RADICALS) {
        let assignment = assign_radicals(&working).map_err(|source| SanitizeError::Radicals {
            stage: SanitizeStage::FindRadicals,
            source,
        })?;
        materialize_radicals(&mut working, &assignment.radical_electrons);
    }

    if operations.contains(SanitizeOperations::SET_AROMATICITY) {
        let ring_assignment = if let Some(existing) = &rings {
            existing.clone()
        } else {
            let calculated =
                symmetrized_sssr(&working, &RingSearchParams::default()).map_err(|source| {
                    SanitizeError::Rings {
                        stage: SanitizeStage::SetAromaticity,
                        source,
                    }
                })?;
            rings = Some(calculated.clone());
            calculated
        };
        working = assign_aromaticity_with_query_state(
            &working,
            &ring_assignment,
            &AromaticityParams::default(),
            query_state,
        )
        .map_err(|source| SanitizeError::Aromaticity {
            stage: SanitizeStage::SetAromaticity,
            source,
        })?
        .topology;
    }

    if operations.contains(SanitizeOperations::SET_CONJUGATION) {
        working = assign_conjugation(&working, &valence).map_err(|source| {
            SanitizeError::Conjugation {
                stage: SanitizeStage::SetConjugation,
                source,
            }
        })?;
    }

    let mut hybridization = None;
    if operations.contains(SanitizeOperations::SET_HYBRIDIZATION) {
        let assignment = assign_hybridization(&working, &valence).map_err(|source| {
            SanitizeError::Hybridization {
                stage: SanitizeStage::SetHybridization,
                source,
            }
        })?;
        materialize_hybridization(&mut working, &assignment);
        hybridization = Some(assignment);
    }

    if operations.contains(SanitizeOperations::CLEANUP_ATROPISOMERS) {
        let assignment = if let Some(assignment) = hybridization {
            assignment
        } else if working.atoms.first().map_or(true, |atom| {
            atom.hybridization() != Hybridization::Unspecified
        }) {
            current_hybridization(&working)
        } else {
            // `MolOps::Hybridizations(mol)` computes conjugation and
            // hybridization in a copy when the current atom state is absent.
            let nested_valence =
                assign_property_cache(&working, &PropertyCacheParams { strict: false })
                    .map_err(|source| SanitizeError::Properties {
                        stage: SanitizeStage::CleanupAtropisomers,
                        source,
                    })?
                    .into_valence();
            let nested_topology =
                assign_conjugation(&working, &nested_valence).map_err(|source| {
                    SanitizeError::Conjugation {
                        stage: SanitizeStage::CleanupAtropisomers,
                        source,
                    }
                })?;
            assign_hybridization(&nested_topology, &nested_valence).map_err(|source| {
                SanitizeError::Hybridization {
                    stage: SanitizeStage::CleanupAtropisomers,
                    source,
                }
            })?
        };
        let ring_assignment = if let Some(existing) = &rings {
            existing.clone()
        } else {
            let calculated =
                find_sssr(&working, &RingSearchParams::default()).map_err(|source| {
                    SanitizeError::Rings {
                        stage: SanitizeStage::CleanupAtropisomers,
                        source,
                    }
                })?;
            calculated
        };
        working = crate::atropisomer::cleanup_invalid_atropisomers(
            &working,
            &assignment,
            &ring_assignment,
        )
        .map_err(|source| SanitizeError::Atropisomers {
            stage: SanitizeStage::CleanupAtropisomers,
            source,
        })?;
    }

    if operations.contains(SanitizeOperations::CLEANUP_CHIRALITY) {
        working =
            crate::structure_tags::cleanup_chirality(&working, &valence).map_err(|source| {
                SanitizeError::Chirality {
                    stage: SanitizeStage::CleanupChirality,
                    source,
                }
            })?;
    }

    if operations.contains(SanitizeOperations::ADJUST_HS) {
        working = adjust_hs(&working, &valence)
            .map_err(|source| SanitizeError::AdjustHs {
                stage: SanitizeStage::AdjustHs,
                source,
            })?
            .topology;
    }

    if operations.contains(SanitizeOperations::PROPERTIES) {
        assign_property_cache(&working, &PropertyCacheParams { strict: true }).map_err(
            |source| SanitizeError::Properties {
                stage: SanitizeStage::Properties,
                source,
            },
        )?;
    }

    working
        .validate()
        .map_err(|source| SanitizeError::InvalidTopology {
            stage: SanitizeStage::None,
            source,
        })?;
    Ok(SanitizeAssignment { topology: working })
}

fn is_source_kekulize_problem(error: &KekulizeError) -> bool {
    match error {
        KekulizeError::AromaticAtomOutsideRing { .. } | KekulizeError::NotKekulizable { .. } => {
            true
        }
        KekulizeError::Valence(source) => {
            matches!(source, ValenceError::InvalidValence { .. })
        }
        _ => false,
    }
}

/// Detect all source-equivalent atom-valence problems followed by at most one
/// source-equivalent kekulization problem.
///
/// The borrowed input is never mutated. Errors outside RDKit's caught
/// `MolSanitizeException` boundary remain top-level [`SanitizeError`] values.
pub fn detect_chemistry_problems(
    topology: &TopologyBlock,
    params: &SanitizeParams,
) -> Result<ChemistryProblemReport, SanitizeError> {
    topology
        .validate()
        .map_err(|source| SanitizeError::InvalidTopology {
            stage: SanitizeStage::None,
            source,
        })?;

    // Complete pinned source: MolOps.cpp::detectChemistryProblems.
    // RDKit✔️❌: std::vector<std::unique_ptr<MolSanitizeException>> detectChemistryProblems(
    // RDKit✔️❌:     const ROMol &imol, unsigned int sanitizeOps) {
    // RDKit✔️❌:   RWMol mol(imol);
    // RDKit✔️❌:   std::vector<std::unique_ptr<MolSanitizeException>> res;
    // RDKit✔️❌:
    // RDKit✔️❌:   // clear out any cached properties
    // RDKit✔️❌:   mol.clearComputedProps();
    // RDKit✔️❌:
    // RDKit✔️❌:   int operation;
    // RDKit✔️❌:   operation = SANITIZE_CLEANUP;
    // RDKit✔️❌:   if (sanitizeOps & operation) {
    // RDKit✔️❌:     // clean up things like nitro groups
    // RDKit✔️❌:     cleanUp(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // update computed properties on atoms and bonds:
    // RDKit✔️❌:   operation = SANITIZE_PROPERTIES;
    // RDKit✔️❌:   if (sanitizeOps & operation) {
    // RDKit✔️❌:     for (auto &atom : mol.atoms()) {
    // RDKit✔️❌:       try {
    // RDKit✔️❌:         bool strict = true;
    // RDKit✔️❌:         atom->updatePropertyCache(strict);
    // RDKit✔️❌:       } catch (const MolSanitizeException &e) {
    // RDKit✔️❌:         res.emplace_back(e.copy());
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     mol.updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // kekulizations
    // RDKit✔️❌:   operation = SANITIZE_KEKULIZE;
    // RDKit✔️❌:   if (sanitizeOps & operation) {
    // RDKit✔️❌:     try {
    // RDKit✔️❌:       Kekulize(mol, true, false);
    // RDKit✔️❌:     } catch (const MolSanitizeException &e) {
    // RDKit✔️❌:       res.emplace_back(e.copy());
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // The detached owner APIs for cleanup and kekulization return owned
    // topologies, so selected stages add full-topology clones beyond the one
    // source-mandated `RWMol` copy. Atom traversal and problem collection
    // retain the source linear complexity and stable stored-row ordering.

    let operations = params.operations;
    let mut working = topology.clone();
    clear_topology_computed_properties(&mut working);
    let mut report = ChemistryProblemReport::default();

    if operations.contains(SanitizeOperations::CLEANUP) {
        working = cleanup(
            &working,
            &CleanupParams {
                charge_normalization: true,
                organometallics: false,
            },
        )
        .map_err(|source| SanitizeError::Cleanup {
            stage: SanitizeStage::Cleanup,
            source,
        })?;
    }

    if operations.contains(SanitizeOperations::PROPERTIES) {
        for atom_index in 0..working.atoms.len() {
            match assign_valence_state_for_atom_from_parts(
                &working.atoms,
                &working.bonds,
                &working.adjacency,
                AtomId::new(atom_index),
                true,
            ) {
                Ok(_) => {}
                Err(source @ ValenceError::InvalidValence { .. }) => {
                    report.problems.push(ChemistryProblem {
                        operation: SanitizeStage::Properties,
                        error: ChemistryProblemError::Valence(source),
                    });
                }
                Err(source) => {
                    return Err(SanitizeError::Properties {
                        stage: SanitizeStage::Properties,
                        source: PropertyCacheError::Valence(source),
                    });
                }
            }
        }
    } else {
        assign_property_cache(&working, &PropertyCacheParams { strict: false }).map_err(
            |source| SanitizeError::Properties {
                stage: SanitizeStage::Properties,
                source,
            },
        )?;
    }

    if operations.contains(SanitizeOperations::KEKULIZE) {
        if let Err(source) = kekulize(
            &working,
            &KekulizeParams {
                mark_atoms_bonds: true,
                canonical: false,
                max_backtracks: KekulizeParams::default().max_backtracks,
            },
        ) {
            if is_source_kekulize_problem(&source) {
                report.problems.push(ChemistryProblem {
                    operation: SanitizeStage::Kekulize,
                    error: ChemistryProblemError::Kekulize(source),
                });
            } else {
                return Err(SanitizeError::Kekulize {
                    stage: SanitizeStage::Kekulize,
                    source,
                });
            }
        }
    }

    Ok(report)
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
        let assignment =
            assign_property_cache(&topology, &PropertyCacheParams { strict: false }).unwrap();
        assert_eq!(assignment.valence.explicit_valence, vec![1, 2, 1]);
        assert_eq!(assignment.valence.implicit_hydrogens, vec![3, 2, 1]);
    }
}
