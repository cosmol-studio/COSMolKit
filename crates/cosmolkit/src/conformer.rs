//! Canonical conformer value/query projections; generation has one detached owner.
use crate::Molecule;
pub use cosmolkit_conformer::{ConformerError, EmbedFailureCause, EmbedParams};
use std::{error::Error, fmt, sync::Arc};
#[derive(Clone, Debug)]
pub struct ConformerRunError(Arc<ConformerRunFailure>);
#[derive(Debug)]
enum ConformerRunFailure {
    Bounds(cosmolkit_conformer::BoundsQueryError),
    Generation(cosmolkit_conformer::GenerationError),
}
impl PartialEq for ConformerRunError {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}
impl fmt::Display for ConformerRunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self.source().unwrap(), f)
    }
}
impl Error for ConformerRunError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(match self.0.as_ref() {
            ConformerRunFailure::Bounds(e) => e,
            ConformerRunFailure::Generation(e) => e,
        })
    }
}
impl ConformerRunError {
    pub(crate) fn generation(e: cosmolkit_conformer::GenerationError) -> crate::OperationError {
        crate::OperationError::Conformer(Self(Arc::new(ConformerRunFailure::Generation(e))))
    }
}
impl Molecule {
    pub fn num_3d_conformers(&self) -> usize {
        self.conformers_3d().len()
    }
    pub fn dg_bounds_matrix(&self) -> Result<Vec<Vec<f64>>, ConformerRunError> {
        // Source query calls the sole graph-bounds owner; model cache preparation
        // supplies the explicit detached state implicit in the source ROMol.
        // RDKit❗❌: PyObject *getMolBoundsMatrix(ROMol &mol, bool set15bounds = true,
        // RDKit❗❌:                              bool scaleVDW = false,
        // RDKit❗❌:                              bool doTriangleSmoothing = true,
        // RDKit❗❌:                              bool useMacrocycle14config = false) {
        // RDKit❗❌:   unsigned int nats = mol.getNumAtoms();
        // RDKit❗❌:   npy_intp dims[2];
        // RDKit❗❌:   dims[0] = nats;
        // RDKit❗❌:   dims[1] = nats;
        // RDKit❗❌:
        // RDKit❗❌:   DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(nats));
        // RDKit❗❌:   DGeomHelpers::initBoundsMat(mat);
        // RDKit❗❌:   DGeomHelpers::setTopolBounds(mol, mat, set15bounds, scaleVDW,
        // RDKit❗❌:                                useMacrocycle14config);
        // RDKit❗❌:   if (doTriangleSmoothing) {
        // RDKit❗❌:     DistGeom::triangleSmoothBounds(mat);
        // RDKit❗❌:   }
        // RDKit❗❌:   auto *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        // RDKit❗❌:   memcpy(static_cast<void *>(PyArray_DATA(res)),
        // RDKit❗❌:          static_cast<void *>(mat->getData()), nats * nats * sizeof(double));
        // RDKit❗❌:
        // RDKit❗❌:   return PyArray_Return(res);
        // RDKit❗❌: }
        // Canonical Rust materialization uses row vectors; the Python projection
        // returns the source-shaped dense array. Matrix materialization visits
        // O(N^2) entries, but this wrapper has known additional costs: the unique
        // public_bounds owner eagerly recomputes valence, SSSR, conjugation and
        // hybridization rather than reading an initialized source ROMol. It
        // copies the flat bounds into N separately allocated row buffers; the
        // Python projection then flattens them into another contiguous buffer.
        // The source wrapper allocates one array and directly memcpy's the flat
        // matrix. These preparation/allocation costs make wrapper equivalence
        // false; the unique dense build_bounds_matrix path is reviewed separately.
        cosmolkit_conformer::dg_bounds_matrix(self.topology())
            .map_err(|e| ConformerRunError(Arc::new(ConformerRunFailure::Bounds(e))))
    }
}

pub(crate) struct EmbedConformerReport {
    pub conf_ids: Vec<i32>,
    pub params: EmbedParams,
    pub requested_num_confs: u32,
}
#[derive(Clone)]
pub struct EmbedMoleculeResult {
    pub molecule: Molecule,
    pub conf_id: i32,
    pub params: EmbedParams,
}
#[derive(Clone)]
pub struct EmbedMultipleConfsResult {
    pub molecule: Molecule,
    pub conf_ids: Vec<i32>,
    pub requested_num_confs: u32,
    pub params: EmbedParams,
}

/// Invalid access to embedding diagnostics; not an embedding failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EmbedFailureCountError {
    InvalidCause,
    MissingCounter { index: usize, available: usize },
}
impl fmt::Display for EmbedFailureCountError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidCause => f.write_str("EndOfEnum is not an embedding failure cause"),
            Self::MissingCounter { index, available } => write!(
                f,
                "embedding failure counter {index} is absent (available: {available})"
            ),
        }
    }
}
impl std::error::Error for EmbedFailureCountError {}

fn failure_count(
    params: &EmbedParams,
    cause: EmbedFailureCause,
) -> Result<Option<u32>, EmbedFailureCountError> {
    if cause == EmbedFailureCause::EndOfEnum {
        return Err(EmbedFailureCountError::InvalidCause);
    }
    if !params.track_failures {
        return Ok(None);
    }
    let index = cause.rdkit_ordinal() as usize;
    params
        .failures
        .get(index)
        .copied()
        .map(Some)
        .ok_or(EmbedFailureCountError::MissingCounter {
            index,
            available: params.failures.len(),
        })
}

#[cfg(test)]
mod failure_diagnostic_tests {
    use super::*;

    #[test]
    fn direct_result_diagnostics_distinguish_untracked_zero_and_invalid_access() {
        let mut params = EmbedParams::etkdg_v3();
        params.failures = vec![0; EmbedFailureCause::EndOfEnum as usize];
        params.failures[EmbedFailureCause::InitialCoords as usize] = 3;
        let mut result = EmbedMultipleConfsResult {
            molecule: Molecule::default(),
            conf_ids: vec![],
            requested_num_confs: 8,
            params,
        };
        assert!(!result.failures_tracked());
        assert_eq!(result.failure_counts(), None);
        assert_eq!(
            result.failure_count(EmbedFailureCause::InitialCoords),
            Ok(None)
        );
        assert_eq!(
            result.failure_count(EmbedFailureCause::EndOfEnum),
            Err(EmbedFailureCountError::InvalidCause)
        );
        result.params.track_failures = true;
        assert!(std::ptr::eq(
            result.failure_counts().unwrap().as_ptr(),
            result.params.failures.as_ptr()
        ));
        assert_eq!(
            result.failure_count(EmbedFailureCause::InitialCoords),
            Ok(Some(3))
        );
        assert_eq!(result.failure_count(EmbedFailureCause::Clash), Ok(Some(0)));
        let single = EmbedMoleculeResult {
            molecule: result.molecule.clone(),
            conf_id: -1,
            params: result.params.clone(),
        };
        assert!(single.failures_tracked());
        assert_eq!(single.failure_counts(), result.failure_counts());
        assert_eq!(
            single.failure_count(EmbedFailureCause::InitialCoords),
            Ok(Some(3))
        );
        result.params.failures.clear();
        assert_eq!(
            result.failure_count(EmbedFailureCause::Clash),
            Err(EmbedFailureCountError::MissingCounter {
                index: 14,
                available: 0
            })
        );
    }
}
impl From<(Molecule, EmbedConformerReport)> for EmbedMoleculeResult {
    fn from((molecule, report): (Molecule, EmbedConformerReport)) -> Self {
        // RDKit❗✔️: inline int EmbedMolecule(ROMol &mol, EmbedParameters &params) {
        // RDKit❗✔️:   INT_VECT confIds;
        // RDKit❗✔️:   EmbedMultipleConfs(mol, confIds, 1, params);
        // RDKit❗✔️:   int res;
        // RDKit❗✔️:   if (confIds.size()) {
        // RDKit❗✔️:     res = confIds[0];
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     res = -1;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // The existing detached owner already performed the source generation.
        // Only its source ID selection and configuration move occur after checked finish.
        Self {
            molecule,
            conf_id: report.conf_ids.first().copied().unwrap_or(-1),
            params: report.params,
        }
    }
}
impl From<(Molecule, EmbedConformerReport)> for EmbedMultipleConfsResult {
    fn from((molecule, report): (Molecule, EmbedConformerReport)) -> Self {
        Self {
            molecule,
            conf_ids: report.conf_ids,
            requested_num_confs: report.requested_num_confs,
            params: report.params,
        }
    }
}
impl EmbedMoleculeResult {
    /// Whether failure events were collected for this invocation.
    pub fn failures_tracked(&self) -> bool {
        self.params.track_failures
    }
    /// Counts in EmbedFailureCause order, or None when tracking was disabled.
    /// These count intermediate failure events, not unsuccessful conformers.
    pub fn failure_counts(&self) -> Option<&[u32]> {
        self.failures_tracked()
            .then_some(self.params.failures.as_slice())
    }
    pub fn failure_count(
        &self,
        cause: EmbedFailureCause,
    ) -> Result<Option<u32>, EmbedFailureCountError> {
        failure_count(&self.params, cause)
    }
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
    pub fn params(&self) -> &EmbedParams {
        &self.params
    }
    pub fn conf_id(&self) -> i32 {
        self.conf_id
    }
    pub fn ok(&self) -> bool {
        self.conf_id >= 0
    }
}
impl fmt::Debug for EmbedMoleculeResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EmbedMoleculeResult")
            .field("conf_id", &self.conf_id)
            .field("ok", &self.ok())
            .finish()
    }
}
impl EmbedMultipleConfsResult {
    /// Whether failure events were collected for this invocation.
    pub fn failures_tracked(&self) -> bool {
        self.params.track_failures
    }
    /// Counts in EmbedFailureCause order, or None when tracking was disabled.
    /// These count intermediate failure events, not unsuccessful conformers.
    pub fn failure_counts(&self) -> Option<&[u32]> {
        self.failures_tracked()
            .then_some(self.params.failures.as_slice())
    }
    pub fn failure_count(
        &self,
        cause: EmbedFailureCause,
    ) -> Result<Option<u32>, EmbedFailureCountError> {
        failure_count(&self.params, cause)
    }
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
    pub fn params(&self) -> &EmbedParams {
        &self.params
    }
    pub fn conf_ids(&self) -> &[i32] {
        &self.conf_ids
    }
    pub fn requested_num_confs(&self) -> u32 {
        self.requested_num_confs
    }
    pub fn generated_count(&self) -> usize {
        self.conf_ids.len()
    }
}
impl fmt::Debug for EmbedMultipleConfsResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EmbedMultipleConfsResult")
            .field("requested_num_confs", &self.requested_num_confs)
            .field("generated_count", &self.generated_count())
            .finish()
    }
}
