use std::error::Error;
use std::fmt;

use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

use super::api::UffParameterError;
use super::atom_typer::UffTypingDiagnostic;
use super::convenience::{
    OptimizationOutcome, SerialUffOptimizationError, SingleConformerOptimizationError,
    SingleConformerOptions,
};

#[cfg(not(target_family = "wasm"))]
use super::convenience::{DispatchedUffOptimizationError, PreparedConformerDispatchOutcome};

#[derive(Debug)]
pub(super) enum UffPreparedOptimizationError {
    Preparation(UffParameterError),
    Single(SingleConformerOptimizationError),
    Serial(SerialUffOptimizationError),
    ConformerStage(super::convenience::SerialConformerOptimizationError),
    #[cfg(not(target_family = "wasm"))]
    Dispatch(DispatchedUffOptimizationError),
}

impl fmt::Display for UffPreparedOptimizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Preparation(source) => fmt::Display::fmt(source, formatter),
            Self::Single(source) => fmt::Display::fmt(source, formatter),
            Self::Serial(source) => fmt::Display::fmt(source, formatter),
            Self::ConformerStage(source) => fmt::Display::fmt(source, formatter),
            #[cfg(not(target_family = "wasm"))]
            Self::Dispatch(source) => fmt::Display::fmt(source, formatter),
        }
    }
}

impl Error for UffPreparedOptimizationError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Preparation(source) => Some(source),
            Self::Single(source) => Some(source),
            Self::Serial(source) => Some(source),
            Self::ConformerStage(source) => Some(source),
            #[cfg(not(target_family = "wasm"))]
            Self::Dispatch(source) => Some(source),
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn optimize_prepared_uff_single(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
) -> Result<OptimizationOutcome, UffPreparedOptimizationError> {
    // BEGIN RDKIT CPP FUNCTION UFF::UFFOptimizeMolecule (UFF.h:40-47)
    // RDKit❗❌: inline std::pair<int, double> UFFOptimizeMolecule(
    // RDKit❗❌:     ROMol &mol, int maxIters = 1000, double vdwThresh = 10.0, int confId = -1,
    // RDKit❗❌:     bool ignoreInterfragInteractions = true) {
    // RDKit❗❌:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit❗❌:       mol, vdwThresh, confId, ignoreInterfragInteractions));
    // RDKit❗❌:   std::pair<int, double> res =
    // RDKit❗❌:       ForceFieldsHelper::OptimizeMolecule(*ff, maxIters);
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION UFF::UFFOptimizeMolecule

    // BEGIN RDKIT CPP FUNCTION ForceFieldHelpers::OptimizeMolecule (FFConvenience.h:94-101)
    // RDKit❗✔️: inline std::pair<int, double> OptimizeMolecule(ForceFields::ForceField &ff,
    // RDKit❗✔️:                                                int maxIters = 1000) {
    // RDKit❗✔️:   ff.initialize();
    // RDKit❗✔️:   int res = ff.minimize(maxIters);
    // RDKit❗✔️:   double e = ff.calcEnergy();
    // RDKit❗✔️:   return std::make_pair(res, e);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION ForceFieldHelpers::OptimizeMolecule

    // Behavior review — RDKit❗❌: this boundary creates one borrowed
    // Cached { topology, assignment } input after the existing cache length,
    // signed-byte/noImplicit, then topology checks. The shared typer reads
    // total valence from those rows and requests ordered incident-bond state
    // only at the source SP2 expression after its aromatic short-circuit. The
    // selected conformer ID and existing construction/initialize/minimize/
    // energy order and typed outcomes remain unchanged.
    // Complexity review — RDKit❗❌: cache validation is O(V); topology
    // validation remains O(V+E) and constructs expected adjacency temporarily.
    // Neither the total-valence nor conjugation projection Vec is created.
    // Typing reads total valence in O(1) and scans incident bonds lazily only
    // when the source branch evaluates the conjugation predicate. The existing
    // builder/optimizer and numerical loop are reused once; their costs remain.
    let prepared = super::api::prepare_parameter_query(topology, valence)
        .map_err(UffPreparedOptimizationError::Preparation)?;

    super::convenience::optimize_single_conformer_with_state(
        topology,
        coordinates,
        prepared.typing_state,
        rings,
        valence,
        properties,
        diagnostics,
        options,
    )
    .map_err(UffPreparedOptimizationError::Single)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn optimize_prepared_uff_serial(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    results: &mut Vec<OptimizationOutcome>,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
) -> Result<(), UffPreparedOptimizationError> {
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

    // BEGIN RDKIT CPP FUNCTION ForceFieldsHelper::OptimizeMoleculeConfsST (FFConvenience.h:63-78)
    // RDKit❗❌: inline void OptimizeMoleculeConfsST(ROMol &mol, ForceFields::ForceField &ff,
    // RDKit❗❌:                                     std::vector<std::pair<int, double>> &res,
    // RDKit❗❌:                                     int maxIters) {
    // RDKit❗❌:   PRECONDITION(res.size() >= mol.getNumConformers(),
    // RDKit❗❌:                "res.size() must be >= mol.getNumConformers()");
    // RDKit❗❌:   unsigned int i = 0;
    // RDKit❗❌:   for (ROMol::ConformerIterator cit = mol.beginConformers();
    // RDKit❗❌:        cit != mol.endConformers(); ++cit, ++i) {
    // RDKit❗❌:     for (unsigned int aidx = 0; aidx < mol.getNumAtoms(); ++aidx) {
    // RDKit❗❌:       ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    // RDKit❗❌:     }
    // RDKit❗❌:     ff.initialize();
    // RDKit❗❌:     int needsMore = ff.minimize(maxIters);
    // RDKit❗❌:     double e = ff.calcEnergy();
    // RDKit❗❌:     res[i] = std::make_pair(needsMore, e);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION ForceFieldsHelper::OptimizeMoleculeConfsST

    // Behavior review — RDKit❗❌: the prepared input carries one validated
    // borrow of topology and cached assignment into the existing same-block
    // serial owner. Its shared typer preserves direct total-valence reads and
    // requests ordered incident-bond state only when the source SP2 expression
    // reaches the conjugation call after the aromatic short-circuit. The owner
    // still constructs once, resizes after success, visits stored 3D rows in
    // order, and runs initialize/minimize/final-energy with the same outcomes.
    // Complexity review — RDKit❗❌: cached-row validation remains O(V), and
    // topology validation remains O(V+E) with temporary expected adjacency;
    // neither chemistry projection Vec is allocated. Total valence is read in
    // O(1), while the source-triggered conjugation scan stays branch-lazy. The
    // serial row adapter is O(1), and the owner reuses its O(A) position-handle
    // Vec. Native worker dispatch separately retains O(C) borrowed-row
    // metadata; it is not removed by this serial adapter.
    let prepared = super::api::prepare_parameter_query(topology, valence)
        .map_err(UffPreparedOptimizationError::Preparation)?;

    super::convenience::optimize_serial_uff_coordinate_block(
        topology,
        coordinates,
        results,
        prepared.typing_state,
        rings,
        valence,
        properties,
        diagnostics,
        options,
    )
    .map_err(UffPreparedOptimizationError::Serial)
}

#[cfg(not(target_family = "wasm"))]
#[allow(clippy::too_many_arguments)]
pub(super) fn optimize_prepared_uff_dispatch(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    results: &mut Vec<OptimizationOutcome>,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
    requested_threads: i32,
    observed_hardware: u32,
    threadsafe: bool,
) -> Result<PreparedConformerDispatchOutcome, UffPreparedOptimizationError> {
    optimize_prepared_uff_dispatch_with_observer(
        topology,
        coordinates,
        results,
        valence,
        rings,
        properties,
        diagnostics,
        options,
        requested_threads,
        || Ok(observed_hardware),
        threadsafe,
    )
}

#[cfg(not(target_family = "wasm"))]
#[allow(clippy::too_many_arguments)]
pub(super) fn optimize_prepared_uff_dispatch_with_observer(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    results: &mut Vec<OptimizationOutcome>,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
    requested_threads: i32,
    observe_hardware: impl FnOnce() -> Result<u32, cosmolkit_core::ThreadCountError>,
    threadsafe: bool,
) -> Result<PreparedConformerDispatchOutcome, UffPreparedOptimizationError> {
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

    // BEGIN RDKIT CPP FUNCTION ForceFieldsHelper::OptimizeMoleculeConfs (FFConvenience.h:115-129)
    // RDKit❗✔️: inline void OptimizeMoleculeConfs(ROMol &mol, ForceFields::ForceField &ff,
    // RDKit❗✔️:                                   std::vector<std::pair<int, double>> &res,
    // RDKit❗✔️:                                   int numThreads = 1, int maxIters = 1000) {
    // RDKit❗✔️:   res.resize(mol.getNumConformers());
    // RDKit❗✔️:   numThreads = getNumThreadsToUse(numThreads);
    // RDKit❗✔️:   if (numThreads == 1) {
    // RDKit❗✔️:     detail::OptimizeMoleculeConfsST(mol, ff, res, maxIters);
    // RDKit❗✔️:   }
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:   else {
    // RDKit❗✔️:     detail::OptimizeMoleculeConfsMT(mol, ff, res, numThreads, maxIters);
    // RDKit❗✔️:   }
    // RDKit❗✔️: #endif
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION ForceFieldsHelper::OptimizeMoleculeConfs

    // Behavior review — RDKit❗❌: the caller-prepared boundary validates once
    // and passes one borrowed topology/assignment state to the existing
    // same-block builder and dispatcher. The shared typer reads cached total
    // valence directly and requests ordered incident-bond state only when the
    // source SP2 expression reaches the conjugation call after its aromatic
    // short-circuit. Requested threads, observed hardware, source-mode flag,
    // and every nested raw Serial/Workers result or panic payload are retained.
    // Complexity review — RDKit❗❌: cache validation remains O(V); topology
    // validation remains O(V+E) with temporary expected adjacency. No total-
    // valence or conjugation output Vec is allocated; typing uses O(1) valence
    // reads and a source-branch-lazy incident scan. Serial routing keeps the
    // O(1) row adapter and reused O(A) position handles. Worker routing keeps
    // its separate O(C) borrowed-row Vec. The unique builder, route, and
    // numerical loops are unchanged; no whole-operation allocation claim follows.
    let prepared = super::api::prepare_parameter_query(topology, valence)
        .map_err(UffPreparedOptimizationError::Preparation)?;

    super::convenience::optimize_dispatched_uff_coordinate_block_with_observer(
        topology,
        coordinates,
        results,
        prepared.typing_state,
        rings,
        valence,
        properties,
        diagnostics,
        options,
        requested_threads,
        observe_hardware,
        threadsafe,
    )
    .map_err(UffPreparedOptimizationError::Dispatch)
}
