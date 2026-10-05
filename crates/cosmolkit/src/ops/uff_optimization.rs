//! Coordinate-only public projection of the detached prepared UFF owner.
use crate::{DerivedState, Molecule, OperationError, PreservationProof};
use cosmolkit_macros::{MoleculeResult, mol_op_body};
use std::{error::Error, fmt, sync::Arc};

/// Options matching UFFOptimizeMolecule; None selects the first stored 3D row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffOptimizationParams {
    pub max_iterations: i32,
    pub vdw_threshold: f64,
    pub ignore_interfragment_interactions: bool,
    pub conformer_id: Option<usize>,
}
impl Default for UffOptimizationParams {
    fn default() -> Self {
        // RDKit✔️✔️: inline std::pair<int, double> UFFOptimizeMolecule(
        // RDKit✔️✔️:     ROMol &mol, int maxIters = 1000, double vdwThresh = 10.0, int confId = -1,
        // RDKit✔️✔️:     bool ignoreInterfragInteractions = true) {
        // CK expresses the source default selector in the independent 3D table.
        // Constant options incur no chemistry work or allocation.
        Self {
            max_iterations: 1000,
            vdw_threshold: 10.0,
            ignore_interfragment_interactions: true,
            conformer_id: None,
        }
    }
}

/// A runtime-finalized molecule, convergence status and final UFF energy.
#[derive(Clone, Debug, PartialEq, MoleculeResult)]
pub struct UffOptimizationResult<M = Molecule> {
    #[pending_molecule]
    pub molecule: M,
    pub status: i32,
    pub energy: f64,
}

/// Options for source-ordered optimization of every stored 3D conformer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffConformerOptimizationParams {
    pub num_threads: i32,
    pub max_iterations: i32,
    pub vdw_threshold: f64,
    pub ignore_interfragment_interactions: bool,
}

impl Default for UffConformerOptimizationParams {
    fn default() -> Self {
        // RDKit✔️✔️: inline void UFFOptimizeMoleculeConfs(ROMol &mol,
        // RDKit✔️✔️:                                      std::vector<std::pair<int, double>> &res,
        // RDKit✔️✔️:                                      int numThreads = 1, int maxIters = 1000,
        // RDKit✔️✔️:                                      double vdwThresh = 10.0,
        // RDKit✔️✔️:                                      bool ignoreInterfragInteractions = true) {
        // Constructing the original options is constant time with no chemistry state.
        Self {
            num_threads: 1,
            max_iterations: 1000,
            vdw_threshold: 10.0,
            ignore_interfragment_interactions: true,
        }
    }
}

/// One source-ordered conformer status and final UFF energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffConformerResult {
    pub conformer_id: usize,
    pub status: i32,
    pub energy: f64,
}

/// A runtime-finalized molecule and its source-ordered conformer outcomes.
#[derive(Clone, Debug, PartialEq, MoleculeResult)]
pub struct UffConformerOptimizationResult<M = Molecule> {
    #[pending_molecule]
    pub molecule: M,
    pub conformers: Vec<UffConformerResult>,
}

/// Stable categories; unavailable geometry is not a guessed 2D fallback.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UffOptimizationErrorKind {
    MissingConformer { requested: Option<usize> },
    Rings,
    Optimization,
    ConformerOptimization,
    Evaluation,
}

#[derive(Debug)]
enum Failure {
    MissingConformer(Option<usize>),
    Rings(cosmolkit_core::RingFindingError),
    Optimization(cosmolkit_forcefields::UffSingleError),
    ConformerOptimization(cosmolkit_forcefields::UffConformerError),
    Evaluation(cosmolkit_forcefields::UffEvaluationError),
}

/// Clones retain the same concrete failure and borrowed source chain.
/// Equality is failure identity, not error-message or chemistry equivalence.
#[derive(Clone, Debug)]
pub struct UffOptimizationError(Arc<Failure>);
impl UffOptimizationError {
    pub fn kind(&self) -> UffOptimizationErrorKind {
        match self.0.as_ref() {
            Failure::MissingConformer(requested) => UffOptimizationErrorKind::MissingConformer {
                requested: *requested,
            },
            Failure::Rings(_) => UffOptimizationErrorKind::Rings,
            Failure::Optimization(_) => UffOptimizationErrorKind::Optimization,
            Failure::ConformerOptimization(_) => UffOptimizationErrorKind::ConformerOptimization,
            Failure::Evaluation(_) => UffOptimizationErrorKind::Evaluation,
        }
    }
    fn operation(cause: Failure) -> OperationError {
        OperationError::UffOptimization(Self(Arc::new(cause)))
    }
}
impl PartialEq for UffOptimizationError {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}
impl fmt::Display for UffOptimizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0.as_ref() {
            Failure::MissingConformer(id) => {
                write!(f, "no stored 3D conformer for selector {id:?}")
            }
            Failure::Rings(cause) => fmt::Display::fmt(cause, f),
            Failure::Optimization(cause) => fmt::Display::fmt(cause, f),
            Failure::ConformerOptimization(cause) => fmt::Display::fmt(cause, f),
            Failure::Evaluation(cause) => fmt::Display::fmt(cause, f),
        }
    }
}
impl Error for UffOptimizationError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self.0.as_ref() {
            Failure::MissingConformer(_) => None,
            Failure::Rings(cause) => Some(cause),
            Failure::Optimization(cause) => Some(cause),
            Failure::ConformerOptimization(cause) => Some(cause),
            Failure::Evaluation(cause) => Some(cause),
        }
    }
}

#[mol_op_body(with_uff_optimized, parts)]
pub(crate) fn with_uff_optimized_impl(
    params: &UffOptimizationParams,
) -> Result<
    UffOptimizationResult<crate::PendingMolecule<super::WithUffOptimizedAccess>>,
    OperationError,
> {
    let mut coordinates = parts.checkout_coordinates()?;
    let outcome = (|| {
        let rows = &coordinates.conformers_3d;
        let selected = match params.conformer_id {
            Some(id) => rows.iter().find(|row| row.id() == id),
            None => rows.first(),
        };
        let id = selected
            .ok_or_else(|| {
                UffOptimizationError::operation(Failure::MissingConformer(params.conformer_id))
            })?
            .id();
        let topology = parts.topology()?;
        let cache = parts.derived_cache()?;
        cache.validate_for_atom_count(topology.atoms.len())?;
        let valence = cache
            .valence_assignment()
            .ok_or(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "assignment",
                actual: 0,
                expected: 1,
            })?;
        #[cfg(feature = "cap-rings")]
        let cached_rings = if cache.valid_states().contains(DerivedState::RINGS) {
            cache.ring_info()
        } else {
            None
        };
        #[cfg(not(feature = "cap-rings"))]
        let cached_rings: Option<&cosmolkit_core::RingInfo> = None;
        let computed_rings;
        let rings = match cached_rings {
            Some(rings) => rings,
            None => {
                computed_rings = cosmolkit_core::symmetrized_sssr(
                    topology,
                    &cosmolkit_core::RingSearchParams::default(),
                )
                .map_err(|cause| UffOptimizationError::operation(Failure::Rings(cause)))?;
                &computed_rings
            }
        };
        // UFF.h:43-47 chemistry belongs exclusively to the detached owner.
        // Runtime borrows trusted caches; it neither reassigns valence nor
        // publishes the temporary ring assignment. Only coordinates use COW.
        cosmolkit_forcefields::optimize_uff_single_prepared(
            topology,
            &mut coordinates,
            valence,
            rings,
            parts.properties()?,
            cosmolkit_forcefields::UffSingleOptions {
                conformer_id: id,
                max_iterations: params.max_iterations,
                vdw_threshold: params.vdw_threshold,
                ignore_interfragment_interactions: params.ignore_interfragment_interactions,
            },
        )
        .map_err(|cause| UffOptimizationError::operation(Failure::Optimization(cause)))
    })();
    // Restore the checked-out block before propagating the domain failure.
    // An error aborts the value transaction; no candidate becomes live.
    parts.install_coordinates(coordinates)?;
    let outcome = outcome?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    Ok(UffOptimizationResult {
        molecule: parts.pending_molecule()?,
        status: outcome.status,
        energy: outcome.energy,
    })
}

#[mol_op_body(with_uff_optimized_confs, parts)]
pub(crate) fn with_uff_optimized_confs_impl(
    params: &UffConformerOptimizationParams,
) -> Result<
    UffConformerOptimizationResult<crate::PendingMolecule<super::WithUffOptimizedConfsAccess>>,
    OperationError,
> {
    // BEGIN RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs (UFF.h:69-78)
    // RDKit❗❌: inline void UFFOptimizeMoleculeConfs(ROMol &mol,
    // RDKit❗❌:                                      std::vector<std::pair<int, double>> &res,
    // RDKit❗❌:                                      int numThreads = 1, int maxIters = 1000,
    // RDKit❗❌:                                      double vdwThresh = 10.0,
    // RDKit❗❌:                                      bool ignoreInterfragInteractions = true) {
    // RDKit❗❌:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit❗❌:       mol, vdwThresh, -1, ignoreInterfragInteractions));
    // RDKit❗❌:   ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads, maxIters);
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs
    // Runtime supplies only trusted topology/cache/property borrows and a
    // coordinate checkout. The detached owner resolves confId=-1 to the
    // first stored 3D ID, constructs once, and returns rows in stored order.
    // This wrapper restores the checked-out coordinate block before surfacing
    // the typed owner failure; no cache or non-coordinate block is committed.
    // Complexity review — RDKit❗❌: the result path has three distinct O(C)
    // row buffers: the owner resizes `Vec<OptimizationOutcome>`, the detached
    // facade collects `Vec<UffConformerOutcome>`, and this wrapper collects
    // `Vec<UffConformerResult>`. This final projection is separate from the
    // runtime's coordinate-block write-checkout/COW clone above. The owner
    // constructs one field and reuses its O(A) position-handle Vec; this result
    // conversion does not copy chemistry, coordinate rows, or contributions.
    let mut coordinates = parts.checkout_coordinates()?;
    let conformers = (|| {
        let topology = parts.topology()?;
        let cache = parts.derived_cache()?;
        cache.validate_for_atom_count(topology.atoms.len())?;
        let valence = cache
            .valence_assignment()
            .ok_or(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "assignment",
                actual: 0,
                expected: 1,
            })?;
        #[cfg(feature = "cap-rings")]
        let cached_rings = if cache.valid_states().contains(DerivedState::RINGS) {
            cache.ring_info()
        } else {
            None
        };
        #[cfg(not(feature = "cap-rings"))]
        let cached_rings: Option<&cosmolkit_core::RingInfo> = None;
        let computed_rings;
        let rings = match cached_rings {
            Some(rings) => rings,
            None => {
                computed_rings = cosmolkit_core::symmetrized_sssr(
                    topology,
                    &cosmolkit_core::RingSearchParams::default(),
                )
                .map_err(|cause| UffOptimizationError::operation(Failure::Rings(cause)))?;
                &computed_rings
            }
        };
        // UFF.h:69-78 chemistry and first-3D selection belong to the existing
        // detached serial owner. Runtime borrows cached valence/rings and
        // publishes neither a recomputed assignment nor temporary cache state.
        cosmolkit_forcefields::optimize_uff_conformers_prepared(
            topology,
            &mut coordinates,
            valence,
            rings,
            parts.properties()?,
            cosmolkit_forcefields::UffConformerOptions {
                num_threads: params.num_threads,
                max_iterations: params.max_iterations,
                vdw_threshold: params.vdw_threshold,
                ignore_interfragment_interactions: params.ignore_interfragment_interactions,
            },
        )
        .map_err(|cause| UffOptimizationError::operation(Failure::ConformerOptimization(cause)))
    })();
    // Restore the checked-out block before propagating the domain failure.
    // An error aborts the value transaction; no candidate becomes live.
    parts.install_coordinates(coordinates)?;
    let conformers = conformers?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    Ok(UffConformerOptimizationResult {
        molecule: parts.pending_molecule()?,
        conformers: conformers
            .into_iter()
            .map(|outcome| UffConformerResult {
                conformer_id: outcome.conformer_id,
                status: outcome.status,
                energy: outcome.energy,
            })
            .collect(),
    })
}

#[cfg(all(
    test,
    feature = "cap-smiles",
    feature = "cap-valence",
    feature = "cap-rings"
))]
mod tests {
    use super::*;
    use crate::{Conformer2D, Conformer3D, CoordinateBlock};

    fn fixture() -> Molecule {
        let chemistry = Molecule::from_smiles("CC.CC").unwrap();
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.9, 0.2, 0.0],
            [5.0, 1.0, 0.0],
            [6.7, 1.1, 0.3],
        ];
        Molecule::from_parts(
            chemistry.topology().clone(),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(3, vec![[9.0, -0.0]; 4])],
                conformers_3d: vec![
                    Conformer3D::new(7, positions.clone(), true),
                    Conformer3D::new(3, positions, true),
                ],
                ..CoordinateBlock::default()
            },
            chemistry.properties().clone(),
        )
        .unwrap()
        .with_assigned_valence()
        .unwrap()
    }

    #[test]
    fn uff_public_option_selection_and_cache_product_preserves_storage() {
        let mut calls = 0;
        for warm_rings in [false, true] {
            let source = if warm_rings {
                fixture().with_assigned_rings().unwrap()
            } else {
                fixture()
            };
            let before = source.clone();
            for iterations in [0, 1, 50] {
                for threshold in [0.0, 10.0] {
                    for ignore in [false, true] {
                        for selector in [None, Some(7), Some(3)] {
                            let params = UffOptimizationParams {
                                max_iterations: iterations,
                                vdw_threshold: threshold,
                                ignore_interfragment_interactions: ignore,
                                conformer_id: selector,
                            };
                            let result = source.with_uff_optimized_with_params(&params).unwrap();
                            calls += 1;
                            assert!(matches!(result.status, 0 | 1));
                            assert!(result.energy.is_finite());
                            let output = &result.molecule;
                            assert!(std::ptr::eq(source.topology(), output.topology()));
                            assert!(std::ptr::eq(source.properties(), output.properties()));
                            assert!(Arc::ptr_eq(
                                &source.derived_cache_arc_runtime(),
                                &output.derived_cache_arc_runtime()
                            ));
                            assert_eq!(
                                output.coordinate_block_runtime().conformers_2d,
                                before.coordinate_block_runtime().conformers_2d
                            );
                            assert_eq!(
                                output
                                    .conformers_3d()
                                    .iter()
                                    .map(|c| c.id())
                                    .collect::<Vec<_>>(),
                                vec![7, 3]
                            );
                            let selected = selector.unwrap_or(7);
                            for row in output.conformers_3d() {
                                if row.id() != selected {
                                    assert_eq!(
                                        row,
                                        before
                                            .conformers_3d()
                                            .iter()
                                            .find(|c| c.id() == row.id())
                                            .unwrap()
                                    );
                                }
                            }
                            assert_eq!(source, before);
                            assert!(std::ptr::eq(
                                source.coordinate_block_runtime(),
                                before.coordinate_block_runtime()
                            ));
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 72);
        let source = fixture();
        assert_eq!(
            source.with_uff_optimized().unwrap(),
            source
                .with_uff_optimized_with_params(&UffOptimizationParams::default())
                .unwrap()
        );
    }

    #[test]
    fn uff_public_missing_geometry_and_invalid_cache_are_atomic_typed_errors() {
        let source = fixture();
        let before = source.clone();
        let error = source
            .with_uff_optimized_with_params(&UffOptimizationParams {
                conformer_id: Some(0),
                ..Default::default()
            })
            .unwrap_err();
        let OperationError::UffOptimization(ref cause) = error else {
            panic!("wrong error: {error}")
        };
        assert_eq!(
            cause.kind(),
            UffOptimizationErrorKind::MissingConformer { requested: Some(0) }
        );
        assert_eq!(error, error.clone());
        assert_eq!(source, before);
        let two_d_only = Molecule::from_parts(
            source.topology().clone(),
            CoordinateBlock {
                conformers_2d: source.coordinate_block_runtime().conformers_2d.clone(),
                ..Default::default()
            },
            source.properties().clone(),
        )
        .unwrap()
        .with_assigned_valence()
        .unwrap();
        let error = two_d_only.with_uff_optimized().unwrap_err();
        let OperationError::UffOptimization(cause) = error else {
            panic!("missing 3D must be typed")
        };
        assert_eq!(
            cause.kind(),
            UffOptimizationErrorKind::MissingConformer { requested: None }
        );
        let invalid_cache = Molecule::from_parts(
            source.topology().clone(),
            source.coordinate_block_runtime().clone(),
            source.properties().clone(),
        )
        .unwrap();
        assert!(matches!(
            invalid_cache.with_uff_optimized(),
            Err(OperationError::InvalidDerivedCache {
                state: "valence",
                ..
            })
        ));
        assert_eq!(
            invalid_cache.derived_cache_runtime().valid_states(),
            DerivedState::NONE
        );
    }

    #[test]
    fn uff_public_owner_failure_retains_borrowed_cause_and_source_storage() {
        // Builder.cpp:281 directionVector rejects zero and sub-tolerance
        // trigonal-bipyramid vectors. Both finite geometric inputs are valid
        // detached coordinates; the actual builder failure must cross runtime.
        for displacement in [0.0, 1.0e-17] {
            let chemistry = Molecule::from_smiles("P(F)(F)(F)(F)F").unwrap();
            let mut positions = vec![[0.0, 0.0, 1.0]; 6];
            positions[0] = [0.0; 3];
            positions[1] = [displacement, 0.0, 0.0];
            let source = Molecule::from_parts(
                chemistry.topology().clone(),
                CoordinateBlock {
                    conformers_3d: vec![Conformer3D::new(0, positions, true)],
                    ..Default::default()
                },
                chemistry.properties().clone(),
            )
            .unwrap()
            .with_assigned_valence()
            .unwrap();
            let before = source.clone();
            let error = source.with_uff_optimized().unwrap_err();
            let OperationError::UffOptimization(cause) = &error else {
                panic!("owner error flattened: {error}")
            };
            assert_eq!(cause.kind(), UffOptimizationErrorKind::Optimization);
            let opaque = cause.source().unwrap();
            assert!(
                opaque
                    .downcast_ref::<cosmolkit_forcefields::UffSingleError>()
                    .is_some()
            );
            let cloned = error.clone();
            assert_eq!(error, cloned);
            // OperationError stores distinct cloned wrapper values; their
            // concrete causes live in the one shared Arc, not the wrappers.
            assert!(std::ptr::eq(
                error.source().unwrap().source().unwrap(),
                cloned.source().unwrap().source().unwrap(),
            ));
            let mut depth = 0;
            let mut child: Option<&(dyn Error + 'static)> = Some(&error);
            while let Some(current) = child {
                depth += 1;
                child = current.source();
            }
            assert!(depth >= 5, "concrete builder cause chain was lost");
            assert_eq!(source, before);
            assert!(std::ptr::eq(
                source.coordinate_block_runtime(),
                before.coordinate_block_runtime()
            ));
            assert!(std::ptr::eq(source.topology(), before.topology()));
            assert!(std::ptr::eq(source.properties(), before.properties()));
            assert!(Arc::ptr_eq(
                &source.derived_cache_arc_runtime(),
                &before.derived_cache_arc_runtime()
            ));
        }
    }
}

impl UffOptimizationResult {
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
    pub fn status_code(&self) -> i32 {
        self.status
    }
    pub fn needs_more(&self) -> bool {
        self.status > 0
    }
    pub fn energy(&self) -> f64 {
        self.energy
    }
}
impl UffConformerResult {
    pub fn conformer_id(&self) -> usize {
        self.conformer_id
    }
    pub fn status_code(&self) -> i32 {
        self.status
    }
    pub fn needs_more(&self) -> bool {
        self.status > 0
    }
    pub fn energy(&self) -> f64 {
        self.energy
    }
}
impl UffConformerOptimizationResult {
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
    pub fn conformer_results(&self) -> &[UffConformerResult] {
        &self.conformers
    }
}

pub use cosmolkit_forcefields::{UffEnergyGradient, UffEvaluationParams};
impl Molecule {
    /// Read-only energy and gradient from the original UFF field constructor.
    pub fn uff_energy_gradient(&self) -> Result<UffEnergyGradient, OperationError> {
        self.uff_energy_gradient_with_params(&UffEvaluationParams::default())
    }
    pub fn uff_energy_gradient_with_params(
        &self,
        params: &UffEvaluationParams,
    ) -> Result<UffEnergyGradient, OperationError> {
        // Select only the independently stored 3D dimension, then borrow
        // prepared chemistry state through the existing public facade.
        let coordinates = self.coordinate_block_runtime();
        let row = match params.conformer_id {
            None => coordinates.conformers_3d.first(),
            Some(id) => coordinates.conformers_3d.iter().find(|row| row.id() == id),
        };
        row.ok_or_else(|| {
            UffOptimizationError::operation(Failure::MissingConformer(params.conformer_id))
        })?;
        let topology = self.topology();
        let cache = self.derived_cache_runtime();
        cache.validate_for_atom_count(topology.atoms.len())?;
        let valence = cache
            .valence_assignment()
            .ok_or(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "assignment",
                actual: 0,
                expected: 1,
            })?;
        let computed_rings;
        let rings = match cache.valid_ring_info() {
            Some(rings) => rings,
            None => {
                computed_rings = cosmolkit_core::symmetrized_sssr(
                    topology,
                    &cosmolkit_core::RingSearchParams::default(),
                )
                .map_err(|e| UffOptimizationError::operation(Failure::Rings(e)))?;
                &computed_rings
            }
        };
        cosmolkit_forcefields::evaluate_uff(
            topology,
            coordinates,
            valence,
            rings,
            self.properties(),
            params,
        )
        .map_err(|e| UffOptimizationError::operation(Failure::Evaluation(e)))
    }
}
