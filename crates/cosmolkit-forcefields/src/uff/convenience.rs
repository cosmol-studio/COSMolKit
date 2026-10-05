use std::num::NonZeroU32;

use crate::kernel::{ForceField, ForceFieldKernelError};
use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

use super::atom_typer::{UffAtomStateRef, UffTypingDiagnostic};
use super::builder::{
    AutomaticForceFieldConstructionError, DEFAULT_TORSION_BOND_SMARTS,
    construct_force_field_with_automatic_typing_from_rows as construct_force_field_with_automatic_typing,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct SingleConformerOptions {
    pub(super) conformer_id: usize,
    pub(super) max_iterations: i32,
    pub(super) vdw_threshold: f64,
    pub(super) ignore_interfragment_interactions: bool,
}

impl SingleConformerOptions {
    pub(super) const fn for_conformer(conformer_id: usize) -> Self {
        // BEGIN RDKIT CPP FUNCTION UFF::UFFOptimizeMolecule option defaults (UFF.h:40-42)
        // RDKit✔️✔️: inline std::pair<int, double> UFFOptimizeMolecule(
        // RDKit✔️✔️:     ROMol &mol, int maxIters = 1000, double vdwThresh = 10.0, int confId = -1,
        // RDKit✔️✔️:     bool ignoreInterfragInteractions = true) {
        // END RDKIT CPP FUNCTION UFF::UFFOptimizeMolecule option defaults
        // confId is deliberately explicit and typed as usize at this private
        // selected-conformer boundary; the source sentinel is not a selection.
        Self {
            conformer_id,
            max_iterations: 1000,
            vdw_threshold: 10.0,
            ignore_interfragment_interactions: true,
        }
    }

    pub(super) const fn max_iterations_as_unsigned(&self) -> u32 {
        // RDKit's int maxIters is passed to this unsigned source parameter.
        // The cast preserves the language-defined modulo-2^32 conversion.
        // RDKit✔️✔️: int minimize(unsigned int maxIts = 200, double forceTol = 1e-4,
        // RDKit✔️✔️:              double energyTol = 1e-6);
        self.max_iterations as u32
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum UffThreadCountError {
    UndefinedSignedNegation,
}

fn resolve_uff_thread_count(
    requested: i32,
    observed_hardware: u32,
    threadsafe: bool,
) -> Result<NonZeroU32, UffThreadCountError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::getNumThreadsToUse (RDGeneral/RDThreads.h:17-42)
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️: #include <thread>
    // RDKit❗✔️:
    // RDKit❗✔️: namespace RDKit {
    // RDKit❗✔️: inline unsigned int getNumThreadsToUse(int target) {
    // RDKit❗✔️:   if (target >= 1) {
    // RDKit❗✔️:     return static_cast<unsigned int>(target);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int res = std::thread::hardware_concurrency();
    // RDKit❗✔️:   if (res > rdcast<unsigned int>(-target)) {
    // RDKit❗✔️:     return res + target;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return 1;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: }  // namespace RDKit
    // RDKit❗✔️:
    // RDKit❗✔️: #else
    // RDKit❗✔️:
    // RDKit❗✔️: namespace RDKit {
    // RDKit❗✔️: inline unsigned int getNumThreadsToUse(int target) {
    // RDKit❗✔️:   RDUNUSED_PARAM(target);
    // RDKit❗✔️:   return 1;
    // RDKit❗✔️: }
    // RDKit❗✔️: }  // namespace RDKit
    // RDKit❗✔️: #endif
    // END RDKIT CPP FUNCTION RDKit::getNumThreadsToUse
    // This is an O(1), allocation-free translation. The caller supplies the
    // source build mode and hardware count; this helper neither queries the
    // host nor selects a threaded build configuration.
    if !threadsafe {
        return Ok(NonZeroU32::new(1).expect("one is nonzero"));
    }
    if requested >= 1 {
        return Ok(NonZeroU32::new(requested as u32).expect("positive target"));
    }

    // RDThreads.h negates a signed int before its unsigned cast. Negation of
    // INT_MIN is undefined in C++; report that unrepresentable source boundary
    // rather than silently wrapping it in Rust.
    let magnitude = requested
        .checked_neg()
        .ok_or(UffThreadCountError::UndefinedSignedNegation)? as u32;
    if observed_hardware > magnitude {
        let count = observed_hardware - magnitude;
        Ok(NonZeroU32::new(count).expect("hardware count exceeds target magnitude"))
    } else {
        Ok(NonZeroU32::new(1).expect("one is nonzero"))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConformerDispatchRoute {
    Serial,
    Workers { source_count: i32 },
}

fn resolve_conformer_dispatch(
    requested: i32,
    observed_hardware: u32,
    threadsafe: bool,
) -> Result<ConformerDispatchRoute, UffThreadCountError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::getNumThreadsToUse (RDGeneral/RDThreads.h:17-42)
    // RDKit❗✔️: inline unsigned int getNumThreadsToUse(int target) {
    // RDKit❗✔️:   if (target >= 1) {
    // RDKit❗✔️:     return static_cast<unsigned int>(target);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int res = std::thread::hardware_concurrency();
    // RDKit❗✔️:   if (res > rdcast<unsigned int>(-target)) {
    // RDKit❗✔️:     return res + target;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return 1;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: inline unsigned int getNumThreadsToUse(int target) {
    // RDKit❗✔️:   RDUNUSED_PARAM(target);
    // RDKit❗✔️:   return 1;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::getNumThreadsToUse
    // BEGIN RDKIT CPP FUNCTION ForceFieldsHelper::OptimizeMoleculeConfs (FFConvenience.h:115-129)
    // RDKit❗✔️:   numThreads = getNumThreadsToUse(numThreads);
    // RDKit❗✔️:   if (numThreads == 1) {
    // RDKit❗✔️:     detail::OptimizeMoleculeConfsST(mol, ff, res, maxIters);
    // RDKit❗✔️:   }
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:   else {
    // RDKit❗✔️:     detail::OptimizeMoleculeConfsMT(mol, ff, res, numThreads, maxIters);
    // RDKit❗✔️:   }
    // RDKit❗✔️: #endif
    // END RDKIT CPP FUNCTION ForceFieldsHelper::OptimizeMoleculeConfs dispatch

    // RDThreads.h's comparison and unsigned addition determine an unsigned
    // source count. The assignment in FFConvenience.h is a separate
    // unsigned-to-signed conversion boundary. The vendored CMakeLists.txt
    // requires C++20; the current x86_64-linux-gnu reference target reports a
    // 32-bit int, so C++20 integral conversion maps this pair modulo 2^32 and
    // this Rust cast preserves the same low bits. The owner receives source
    // build mode and hardware count explicitly and performs no host query.
    // Resolution and route selection are O(1), allocation-free operations.
    let unsigned_count = resolve_uff_thread_count(requested, observed_hardware, threadsafe)?;
    let source_count = unsigned_count.get() as i32;
    if source_count == 1 {
        Ok(ConformerDispatchRoute::Serial)
    } else {
        Ok(ConformerDispatchRoute::Workers { source_count })
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct OptimizationOutcome {
    pub(crate) status: i32,
    pub(crate) energy: f64,
}

fn resize_conformer_results(results: &mut Vec<OptimizationOutcome>, conformer_count: usize) {
    // BEGIN RDKIT CPP FUNCTION ForceFieldsHelper::OptimizeMoleculeConfs (FFConvenience.h:115-129)
    // RDKit✔️✔️:   res.resize(mol.getNumConformers());
    // END RDKIT CPP FUNCTION ForceFieldsHelper::OptimizeMoleculeConfs result resize

    // Behavior review: Vec::resize retains the existing prefix, truncates any
    // excess tail, and initializes appended source pair values to (0, +0.0).
    // Allocation review: this mutates the same contiguous Vec, initializing
    // only newly appended rows; it adds no temporary result buffer or scan.
    results.resize(
        conformer_count,
        OptimizationOutcome {
            status: 0,
            energy: 0.0,
        },
    );
}

pub(crate) struct SerialConformer<'a> {
    pub(crate) id: usize,
    pub(crate) positions: &'a mut [[f64; 3]],
}

struct SourceIndexedWorkerRow<'rows, 'coordinates> {
    source_index: u32,
    conformer_id: usize,
    conformer: &'rows mut SerialConformer<'coordinates>,
}

struct PartitionedWorkerLane<'rows, 'coordinates> {
    rows: Vec<SourceIndexedWorkerRow<'rows, 'coordinates>>,
    source_result_slots: Vec<&'rows mut OptimizationOutcome>,
}

enum WorkerResultSlotView<'results> {
    Dense(&'results mut [OptimizationOutcome]),
    Strided {
        source_result_slots: Vec<&'results mut OptimizationOutcome>,
        lane_count: NonZeroU32,
    },
}

impl WorkerResultSlotView<'_> {
    fn get_mut(&mut self, source_index: u32) -> &mut OptimizationOutcome {
        // RDKit❗✔️: (*res)[i] = std::make_pair(needsMore, e);
        // FFConvenience.h:43. Dense callers use the wrapped source slot
        // directly; partitioned callers use the source quotient within the
        // lane's unique, source-ordered slot borrows. Both accesses are O(1)
        // and perform no allocation or scan.
        match self {
            Self::Dense(results) => results
                .get_mut(source_index as usize)
                .expect("the source result-capacity precondition covers this index"),
            Self::Strided {
                source_result_slots,
                lane_count,
            } => {
                let lane_slot_index = usize::try_from(source_index / lane_count.get())
                    .expect("a source lane slot index fits supported slice lengths");
                &mut **source_result_slots
                    .get_mut(lane_slot_index)
                    .expect("partitioned source slots contain each active lane index")
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct WorkerConformerAssignment {
    source_index: u32,
    thread_idx: u32,
    lane_count: NonZeroU32,
}

impl WorkerConformerAssignment {
    const fn new(thread_idx: u32, lane_count: NonZeroU32) -> Self {
        Self {
            source_index: 0,
            thread_idx,
            lane_count,
        }
    }

    fn selects_current(&self) -> bool {
        // RDKit❗✔️: if (i % numThreads != threadIdx) {
        // RDKit❗✔️:   continue;
        // RDKit❗✔️: }
        // FFConvenience.h:31-33. `NonZeroU32` excludes the source modulo-zero
        // case; modulo returns a lane below count, so an out-of-range thread
        // index selects no conformer. Constant-time arithmetic, no allocation.
        self.source_index % self.lane_count.get() == self.thread_idx
    }

    fn advance(&mut self) {
        // RDKit❗✔️: ++i
        // FFConvenience.h:29,32. Match unsigned int wrap explicitly.
        self.source_index = self.source_index.wrapping_add(1);
    }
}

fn wrapped_source_lane_cardinality(
    source_rows: u64,
    lane_count: NonZeroU32,
    thread_idx: u32,
) -> u64 {
    // RDKit✔️🔝: unsigned int i = 0;
    // RDKit✔️🔝: for (ROMol::ConformerIterator cit = mol->beginConformers();
    // RDKit✔️🔝:      cit != mol->endConformers(); ++cit, ++i) {
    // RDKit✔️🔝:   if (i % numThreads != threadIdx) {
    // RDKit✔️🔝:     continue;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝: }
    // FFConvenience.h:29-34. Count the exact source-selected row sequence for
    // lane-capacity planning. Each complete u32 counter cycle contributes
    // floor(2^32 / W) rows to every lane and one extra to the first remainder
    // lanes; the final prefix contributes its direct modulo counts. This is
    // O(1) checked integer arithmetic with no row scan or allocation, so
    // capacity planning does not add another C-dependent traversal.
    let lane_count = u64::from(lane_count.get());
    let thread_idx = u64::from(thread_idx);
    if thread_idx >= lane_count {
        return 0;
    }

    let counter_cycle = 1_u64 << u32::BITS;
    let complete_cycles = source_rows / counter_cycle;
    let remaining_rows = source_rows % counter_cycle;
    let rows_per_cycle_base = counter_cycle / lane_count;
    let extra_cycle_lanes = counter_cycle % lane_count;
    let rows_per_cycle = if thread_idx < extra_cycle_lanes {
        rows_per_cycle_base + 1
    } else {
        rows_per_cycle_base
    };
    let complete_cycle_rows = complete_cycles
        .checked_mul(rows_per_cycle)
        .expect("a lane count cannot exceed the total source row count");
    let remaining_lane_rows = if remaining_rows > thread_idx {
        (remaining_rows - 1 - thread_idx) / lane_count + 1
    } else {
        0
    };
    complete_cycle_rows
        .checked_add(remaining_lane_rows)
        .expect("a lane count cannot exceed the total source row count")
}

fn partition_worker_lanes<'rows, 'coordinates>(
    conformers: &'rows mut [SerialConformer<'coordinates>],
    results: &'rows mut [OptimizationOutcome],
    num_threads: i32,
) -> Result<Vec<PartitionedWorkerLane<'rows, 'coordinates>>, SerialConformerOptimizationError> {
    // RDKit✔️❌: PRECONDITION(res->size() >= mol->getNumConformers(),
    // RDKit✔️❌:              "res->size() must be >= mol->getNumConformers()");
    // RDKit✔️❌: unsigned int i = 0;
    // RDKit✔️❌: for (ROMol::ConformerIterator cit = mol->beginConformers();
    // RDKit✔️❌:      cit != mol->endConformers(); ++cit, ++i) {
    // RDKit✔️❌:   if (i % numThreads != threadIdx) {
    // RDKit✔️❌:     continue;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   (*res)[i] = std::make_pair(needsMore, e);
    // RDKit✔️❌: }
    // FFConvenience.h:27-44. Nonpositive thread counts create no source
    // workers, so they return before the worker capacity precondition. For a
    // positive count, validate the original lengths first, then retain each
    // row's wrapped source index/ID and borrow unique active result slots by
    // their original source index. Extra result suffix entries stay unborrowed.
    // The source workers each scan C rows; this one-time routing scans C rows
    // and stores O(C+W) borrow metadata, with no coordinate or field copies.
    if num_threads <= 0 {
        return Ok(Vec::new());
    }
    validate_result_capacity(conformers.len(), results.len())
        .map_err(SerialConformerOptimizationError::ResultCapacity)?;

    let worker_count = u32::try_from(num_threads)
        .expect("a positive i32 worker count fits the source unsigned count");
    let lane_count = NonZeroU32::new(worker_count).expect("positive workers are nonzero");
    let worker_count_usize = usize::try_from(worker_count)
        .expect("supported targets can represent a positive i32 lane count");
    let source_rows = u64::try_from(conformers.len())
        .expect("supported Rust slice lengths fit the source-row count");
    let active_result_rows = source_rows.min(1_u64 << u32::BITS);
    let active_result_slots = usize::try_from(active_result_rows)
        .expect("active source result indices fit supported slice lengths");

    let mut lanes: Vec<PartitionedWorkerLane<'rows, 'coordinates>> =
        Vec::with_capacity(worker_count_usize);
    for _ in 0..worker_count_usize {
        lanes.push(PartitionedWorkerLane {
            rows: Vec::new(),
            source_result_slots: Vec::new(),
        });
    }
    for (lane_index, lane) in lanes.iter_mut().enumerate() {
        let thread_idx = u32::try_from(lane_index)
            .expect("lane positions are below the positive i32 worker count");
        let row_capacity = usize::try_from(wrapped_source_lane_cardinality(
            source_rows,
            lane_count,
            thread_idx,
        ))
        .expect("each lane row count is bounded by the source slice length");
        let result_capacity = usize::try_from(wrapped_source_lane_cardinality(
            active_result_rows,
            lane_count,
            thread_idx,
        ))
        .expect("each lane result count is bounded by the active result length");
        lane.rows.reserve_exact(row_capacity);
        lane.source_result_slots.reserve_exact(result_capacity);
    }

    let mut assignment = WorkerConformerAssignment::new(0, lane_count);
    for conformer in conformers.iter_mut() {
        let source_index = assignment.source_index;
        let lane_index = usize::try_from(source_index % worker_count)
            .expect("source lane index fits supported slice lengths");
        let conformer_id = conformer.id;
        lanes[lane_index].rows.push(SourceIndexedWorkerRow {
            source_index,
            conformer_id,
            conformer,
        });
        assignment.advance();
    }

    for (source_slot, result) in results.iter_mut().take(active_result_slots).enumerate() {
        let source_slot = u64::try_from(source_slot)
            .expect("active result source indices fit the source-row count");
        let lane_index = usize::try_from(source_slot % u64::from(worker_count))
            .expect("source result lane index fits supported slice lengths");
        lanes[lane_index].source_result_slots.push(result);
    }
    Ok(lanes)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SerialResultCapacityError {
    pub(super) conformers: usize,
    pub(super) result_slots: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SerialCoordinateCountError {
    pub(super) conformer_id: usize,
    pub(super) atoms: usize,
    pub(super) coordinates: usize,
}

impl std::fmt::Display for SerialResultCapacityError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for SerialResultCapacityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

impl std::fmt::Display for SerialCoordinateCountError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for SerialCoordinateCountError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SerialConformerOptimizationError {
    ResultCapacity(SerialResultCapacityError),
    CoordinateCount {
        input_index: usize,
        source: SerialCoordinateCountError,
    },
    Optimization {
        input_index: usize,
        conformer_id: usize,
        source: OptimizationStageError,
    },
}

#[cfg(not(target_family = "wasm"))]
// Keep serial stage errors and every raw worker join result for the caller.
pub(crate) enum PreparedConformerDispatchOutcome {
    Serial(Result<(), SerialConformerOptimizationError>),
    Workers(
        Result<
            Vec<std::thread::Result<Result<(), SerialConformerOptimizationError>>>,
            SerialConformerOptimizationError,
        >,
    ),
}

#[derive(Debug)]
pub(super) enum SerialUffOptimizationError {
    Construction(AutomaticForceFieldConstructionError),
    Optimization(SerialConformerOptimizationError),
}

#[cfg(not(target_family = "wasm"))]
#[derive(Debug)]
pub(super) enum DispatchedUffOptimizationError {
    Construction(AutomaticForceFieldConstructionError),
    ThreadCount(UffThreadCountError),
}

impl std::fmt::Display for UffThreadCountError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for UffThreadCountError {}

#[cfg(not(target_family = "wasm"))]
impl std::fmt::Display for DispatchedUffOptimizationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Construction(source) => {
                write!(formatter, "UFF force-field construction failed: {source}")
            }
            Self::ThreadCount(source) => {
                write!(
                    formatter,
                    "UFF conformer thread-count resolution failed: {source}"
                )
            }
        }
    }
}

#[cfg(not(target_family = "wasm"))]
impl std::error::Error for DispatchedUffOptimizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Construction(source) => Some(source),
            Self::ThreadCount(source) => Some(source),
        }
    }
}

#[cfg(test)]
fn rebind_serial_conformer<'next>(
    field: ForceField<'_>,
    conformer: &'next mut SerialConformer<'_>,
    atom_count: usize,
) -> Result<ForceField<'next>, SerialCoordinateCountError> {
    // RDKit✔️❌: for (unsigned int aidx = 0; aidx < mol.getNumAtoms(); ++aidx) {
    // RDKit✔️❌:   ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    // RDKit✔️❌: }
    // RDKit✔️✔️: PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    // FFConvenience.h:71-73; Conformer.cpp:40. This one-shot constructor is
    // retained only for S04's isolated borrow-lifetime regressions. Production
    // serial traversal fills the existing PointPtrVect below.
    if atom_count != 0 && conformer.positions.len() != atom_count {
        return Err(SerialCoordinateCountError {
            conformer_id: conformer.id,
            atoms: atom_count,
            coordinates: conformer.positions.len(),
        });
    }
    let positions = conformer
        .positions
        .iter_mut()
        .take(atom_count)
        .map(|position| &mut position[..])
        .collect();
    Ok(field.rebind_positions(positions))
}

fn bind_serial_conformer_positions<'rows>(
    field: &mut ForceField<'rows>,
    conformer: &'rows mut SerialConformer<'_>,
    atom_count: usize,
) -> Result<(), SerialCoordinateCountError> {
    // RDKit✔️✔️: for (unsigned int aidx = 0; aidx < mol.getNumAtoms(); ++aidx) {
    // RDKit✔️✔️:   ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    // RDKit✔️✔️: }
    // FFConvenience.h:71-73; Conformer.cpp:40. The source-sized position
    // vector is allocated once before traversal; its first row fills existing
    // capacity and each later row replaces the same indexed mutable borrows.
    // This is O(atom_count) per source row with no per-row handle allocation.
    if atom_count != 0 && conformer.positions.len() != atom_count {
        return Err(SerialCoordinateCountError {
            conformer_id: conformer.id,
            atoms: atom_count,
            coordinates: conformer.positions.len(),
        });
    }

    let positions = field.positions_mut();
    if positions.len() != atom_count {
        positions.clear();
        positions.extend(
            conformer
                .positions
                .iter_mut()
                .take(atom_count)
                .map(|position| &mut position[..]),
        );
    } else {
        for (position, coordinate) in positions
            .iter_mut()
            .zip(conformer.positions.iter_mut().take(atom_count))
        {
            *position = &mut coordinate[..];
        }
    }
    Ok(())
}

fn validate_result_capacity(
    conformer_count: usize,
    result_slots: usize,
) -> Result<(), SerialResultCapacityError> {
    // RDKit✔️✔️: PRECONDITION(res->size() >= mol->getNumConformers(),
    // RDKit✔️✔️:              "res->size() must be >= mol->getNumConformers()");
    // RDKit✔️✔️: PRECONDITION(res.size() >= mol.getNumConformers(),
    // RDKit✔️✔️:              "res.size() must be >= mol.getNumConformers()");
    // FFConvenience.h:27-30, 66-67. This count owner is also used before a
    // lane is formed, when the complete borrowed slices are no longer present.
    // It preserves the same O(1), allocation-free source capacity precondition.
    if result_slots < conformer_count {
        return Err(SerialResultCapacityError {
            conformers: conformer_count,
            result_slots,
        });
    }
    Ok(())
}

#[cfg(test)]
fn validate_serial_result_capacity(
    conformers: &[SerialConformer<'_>],
    results: &[OptimizationOutcome],
) -> Result<(), SerialResultCapacityError> {
    // Supplied slice order replaces beginConformers iteration; IDs are
    // retained labels and never result-vector offsets.
    validate_result_capacity(conformers.len(), results.len())
}

fn prepare_worker_force_field<'field, 'rows>(
    field: ForceField<'field>,
    conformer_count: usize,
    result_slots: usize,
    atom_count: usize,
) -> Result<ForceField<'rows>, SerialConformerOptimizationError>
where
    'field: 'rows,
{
    // RDKit❗✔️: PRECONDITION(res->size() >= mol->getNumConformers(),
    // RDKit❗✔️:              "res->size() must be >= mol->getNumConformers()");
    // RDKit❗✔️: ff.positions().resize(mol->getNumAtoms());
    // FFConvenience.h:27-35. Reuse the serial capacity error before any
    // traversal, including empty or unassigned workers. The copied field's
    // handle vector stays empty while reserving source-sized capacity: Rust
    // references have no null value, and the selected row fills all handles
    // before initialize. This avoids default-initializing A null slots while
    // preserving O(A) storage and no per-row allocation.
    validate_result_capacity(conformer_count, result_slots)
        .map_err(SerialConformerOptimizationError::ResultCapacity)?;

    let mut field = field.clear_positions_and_shorten_lifetime::<'rows>();
    let positions = field.positions_mut();
    if positions.capacity() < atom_count {
        #[cfg(test)]
        crate::kernel::cf3d_uff_one_note_serial_position_buffer_growth();
        positions.reserve(atom_count);
    }
    #[cfg(test)]
    crate::kernel::cf3d_uff_one_note_serial_position_buffer_start(
        field.positions().as_ptr() as usize,
        field.positions_mut().capacity(),
    );
    Ok(field)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum OptimizationStageError {
    Initialize(ForceFieldKernelError),
    Minimize(ForceFieldKernelError),
    FinalEnergy(ForceFieldKernelError),
}

#[derive(Debug)]
pub(super) enum SingleConformerOptimizationError {
    Construction(AutomaticForceFieldConstructionError),
    Optimization(OptimizationStageError),
}

impl std::fmt::Display for SingleConformerOptimizationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for SingleConformerOptimizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Construction(source) => Some(source),
            Self::Optimization(source) => Some(source),
        }
    }
}

impl std::fmt::Display for SerialUffOptimizationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for SerialUffOptimizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Construction(source) => Some(source),
            Self::Optimization(source) => Some(source),
        }
    }
}

impl std::fmt::Display for OptimizationStageError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for OptimizationStageError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Initialize(source) => Some(source),
            Self::Minimize(source) => Some(source),
            Self::FinalEnergy(source) => Some(source),
        }
    }
}

impl std::fmt::Display for SerialConformerOptimizationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for SerialConformerOptimizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ResultCapacity(source) => Some(source),
            Self::CoordinateCount { source, .. } => Some(source),
            Self::Optimization { source, .. } => Some(source),
        }
    }
}

fn assign_serial_result(
    field: &mut ForceField<'_>,
    max_iterations: u32,
    result: &mut OptimizationOutcome,
) -> Result<(), OptimizationStageError> {
    // RDKit✔️✔️: ff.initialize();
    // RDKit✔️✔️: int needsMore = ff.minimize(maxIters);
    // RDKit✔️✔️: double e = ff.calcEnergy();
    // RDKit✔️✔️: res[i] = std::make_pair(needsMore, e);
    // FFConvenience.h:74-77. A failed stage leaves this slot untouched;
    // callers retain previously assigned slots and stop before later slots.
    // One constant-sized result assignment, with no allocation or cloning.
    let outcome = execute_serial_slot(field, max_iterations)?;
    *result = outcome;
    Ok(())
}

fn optimize_assigned_worker_rows<'field, 'rows, 'coordinates, 'results, Rows>(
    field: ForceField<'field>,
    rows: Rows,
    mut results: WorkerResultSlotView<'results>,
    conformer_count: usize,
    result_slot_count: usize,
    atom_count: usize,
    max_iterations: i32,
) -> Result<(), SerialConformerOptimizationError>
where
    'field: 'rows,
    'coordinates: 'rows,
    Rows: Iterator<Item = SourceIndexedWorkerRow<'rows, 'coordinates>>,
{
    // RDKit❗✔️: for (unsigned int aidx = 0; aidx < mol->getNumAtoms(); ++aidx) {
    // RDKit❗✔️:   ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    // RDKit❗✔️: }
    // RDKit❗✔️: ff.initialize();
    // RDKit❗✔️: int needsMore = ff.minimize(maxIters);
    // RDKit❗✔️: double e = ff.calcEnergy();
    // RDKit❗✔️: (*res)[i] = std::make_pair(needsMore, e);
    // FFConvenience.h:38-43. The caller supplies only rows selected by the
    // source u32 index modulo rule, or the already partitioned lane sequence.
    // Reuse the same capacity-before-reserve setup and complete serial stage
    // driver. Each row binds O(A) coordinate handles, then writes its wrapped
    // source result slot only after final energy succeeds; no row allocation,
    // coordinate copy, field copy, or result scan is performed here.
    let mut field =
        prepare_worker_force_field(field, conformer_count, result_slot_count, atom_count)?;
    let max_iterations = max_iterations as u32;

    for row in rows {
        let input_index = row.source_index as usize;
        let conformer_id = row.conformer_id;
        bind_serial_conformer_positions(&mut field, row.conformer, atom_count).map_err(
            |source| SerialConformerOptimizationError::CoordinateCount {
                input_index,
                source,
            },
        )?;
        if let Err(source) = assign_serial_result(
            &mut field,
            max_iterations,
            results.get_mut(row.source_index),
        ) {
            return Err(SerialConformerOptimizationError::Optimization {
                input_index,
                conformer_id,
                source,
            });
        }
    }
    #[cfg(test)]
    crate::kernel::cf3d_uff_one_note_serial_position_buffer_end(
        field.positions().as_ptr() as usize,
        field.positions_mut().capacity(),
    );
    Ok(())
}

fn optimize_partitioned_worker<'field, 'rows, 'coordinates>(
    field: ForceField<'field>,
    lane: PartitionedWorkerLane<'rows, 'coordinates>,
    conformer_count: usize,
    result_slot_count: usize,
    atom_count: usize,
    lane_count: NonZeroU32,
    max_iterations: i32,
) -> Result<(), SerialConformerOptimizationError>
where
    'field: 'rows,
    'coordinates: 'rows,
{
    // RDKit❗✔️:     ff.initialize();
    // RDKit❗✔️:     int needsMore = ff.minimize(maxIters);
    // RDKit❗✔️:     double e = ff.calcEnergy();
    // RDKit❗✔️:     (*res)[i] = std::make_pair(needsMore, e);
    // FFConvenience.h:38-43. M03 has already routed this lane's rows by the
    // source wrapped-index modulo rule. Consume its source-ordered row and
    // result-slot borrows, preserving the original global capacity counts;
    // M04 remains the single owner of coordinate binding and all result/error
    // stages. The caller moves one lane-local copied ForceField into this
    // adapter. This adds no source-row rescan, field/coordinate copy, or
    // per-row allocation; the lane's borrow metadata was reserved by M03.
    let PartitionedWorkerLane {
        rows,
        source_result_slots,
    } = lane;
    optimize_assigned_worker_rows(
        field,
        rows.into_iter(),
        WorkerResultSlotView::Strided {
            source_result_slots,
            lane_count,
        },
        conformer_count,
        result_slot_count,
        atom_count,
        max_iterations,
    )
}

#[cfg(all(test, not(target_family = "wasm")))]
type UffMtM10CopiedFieldState = (bool, usize, Vec<f64>, bool, Vec<i32>, u32, usize, usize);

#[cfg(all(test, not(target_family = "wasm")))]
#[derive(Debug)]
struct UffMtM10LaneEntry {
    lane_index: usize,
    row_len: usize,
    row_capacity: usize,
    result_slot_len: usize,
    result_slot_capacity: usize,
    copied_field_state: UffMtM10CopiedFieldState,
}

#[cfg(all(test, not(target_family = "wasm")))]
struct UffMtM10Observer {
    started_tx: std::sync::mpsc::Sender<UffMtM10LaneEntry>,
    release_receivers: Vec<std::sync::Mutex<std::sync::mpsc::Receiver<()>>>,
    worker_work_tx: std::sync::mpsc::Sender<(usize, crate::kernel::UffOneSerialWorkCounts)>,
    handle_storage_tx: std::sync::mpsc::Sender<(usize, usize)>,
}

#[cfg(all(test, not(target_family = "wasm")))]
impl UffMtM10Observer {
    fn new(
        worker_count: usize,
    ) -> (
        std::sync::Arc<Self>,
        std::sync::mpsc::Receiver<UffMtM10LaneEntry>,
        Vec<std::sync::mpsc::Sender<()>>,
        std::sync::mpsc::Receiver<(usize, crate::kernel::UffOneSerialWorkCounts)>,
        std::sync::mpsc::Receiver<(usize, usize)>,
    ) {
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (worker_work_tx, worker_work_rx) = std::sync::mpsc::channel();
        let (handle_storage_tx, handle_storage_rx) = std::sync::mpsc::channel();
        let mut release_receivers = Vec::with_capacity(worker_count);
        let mut release_senders = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let (release_tx, release_rx) = std::sync::mpsc::channel();
            release_receivers.push(std::sync::Mutex::new(release_rx));
            release_senders.push(release_tx);
        }
        (
            std::sync::Arc::new(Self {
                started_tx,
                release_receivers,
                worker_work_tx,
                handle_storage_tx,
            }),
            started_rx,
            release_senders,
            worker_work_rx,
            handle_storage_rx,
        )
    }

    fn enter_lane(&self, entry: UffMtM10LaneEntry) {
        let lane_index = entry.lane_index;
        let _ = self.started_tx.send(entry);
        let release = self
            .release_receivers
            .get(lane_index)
            .expect("M10 observer has one release receiver for each lane");
        let _ = release
            .lock()
            .expect("M10 release receiver mutex is not poisoned")
            .recv();
    }

    fn note_worker_work(&self, lane_index: usize) {
        let _ = self
            .worker_work_tx
            .send((lane_index, crate::kernel::cf3d_uff_one_serial_work_counts()));
    }

    fn note_handle_storage(&self, len: usize, capacity: usize) {
        let _ = self.handle_storage_tx.send((len, capacity));
    }
}

#[cfg(all(test, not(target_family = "wasm")))]
std::thread_local! {
    static UFF_MT_M10_OBSERVER: std::cell::RefCell<Option<std::sync::Arc<UffMtM10Observer>>> =
        const { std::cell::RefCell::new(None) };
}

#[cfg(all(test, not(target_family = "wasm")))]
struct UffMtM10ObserverScope;

#[cfg(all(test, not(target_family = "wasm")))]
impl UffMtM10ObserverScope {
    fn install(observer: std::sync::Arc<UffMtM10Observer>) -> Self {
        UFF_MT_M10_OBSERVER.with(|current| *current.borrow_mut() = Some(observer));
        Self
    }
}

#[cfg(all(test, not(target_family = "wasm")))]
impl Drop for UffMtM10ObserverScope {
    fn drop(&mut self) {
        UFF_MT_M10_OBSERVER.with(|current| {
            current.borrow_mut().take();
        });
    }
}

#[cfg(not(target_family = "wasm"))]
fn join_worker_handles<'scope>(
    handles: Vec<
        std::thread::ScopedJoinHandle<'scope, Result<(), SerialConformerOptimizationError>>,
    >,
) -> Vec<std::thread::Result<Result<(), SerialConformerOptimizationError>>> {
    // RDKit❗❌: for (auto &thread : tg) {
    // RDKit❗❌:   if (thread.joinable()) {
    // RDKit❗❌:     thread.join();
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // FFConvenience.h:56-60. Every scoped handle supplied by spawn is
    // joinable and consumed exactly once. Vec iteration keeps creation order;
    // collecting each raw join result retains both the worker's typed Result
    // and its original panic payload instead of choosing a public error policy.
    // This is one O(W) join pass, with one O(W) outcome vector allocation
    // beyond the source handle vector; the retained lossless outcomes require
    // that storage, so this helper does not claim source allocation parity.
    handles.into_iter().map(|handle| handle.join()).collect()
}

#[cfg(not(target_family = "wasm"))]
fn optimize_conformers_mt<'field, 'rows, 'coordinates>(
    source_field: &ForceField<'field>,
    conformers: &'rows mut [SerialConformer<'coordinates>],
    results: &'rows mut [OptimizationOutcome],
    num_threads: i32,
    atom_count: usize,
    max_iterations: i32,
) -> Result<
    Vec<std::thread::Result<Result<(), SerialConformerOptimizationError>>>,
    SerialConformerOptimizationError,
>
where
    'coordinates: 'rows,
{
    // RDKit❗❗: inline void OptimizeMoleculeConfsMT(ROMol &mol,
    // RDKit❗❗:                                     const ForceFields::ForceField &ff,
    // RDKit❗❗:                                     std::vector<std::pair<int, double>> &res,
    // RDKit❗❗:                                     int numThreads, int maxIters) {
    // RDKit❗❗:   std::vector<std::thread> tg;
    // RDKit❗❗:   for (int ti = 0; ti < numThreads; ++ti) {
    // RDKit❗❗:     tg.emplace_back(std::thread(detail::OptimizeMoleculeConfsHelper_, ff, &mol,
    // RDKit❗❗:                                 &res, ti, numThreads, maxIters));
    // RDKit❗❗:   }
    // RDKit❗❗:   for (auto &thread : tg) {
    // RDKit❗❗:     if (thread.joinable()) {
    // RDKit❗❗:       thread.join();
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗: }
    // FFConvenience.h:46-60. Preserve the signed source loop: a nonpositive
    // count creates no lanes, while each positive lane gets one independent
    // by-value field copy before its thread starts, even when it has no rows.
    // M03 validates and partitions original borrowed rows/result slots; each
    // child runs the shared M06 stage loop and M07 retains every join outcome.
    //
    // Complexity: source workers scan C rows per lane (O(W*C)); this
    // source-ordered partition scans C rows once and reserves O(C+W) borrow
    // metadata, then each child visits only its assigned rows. Each lane still
    // copies the field and its terms once as the source thread argument does.
    // The borrowed metadata and M07 raw-outcome vector are additional O(C+W)
    // storage, so this owner does not claim source allocation parity.
    let conformer_count = conformers.len();
    let result_slot_count = results.len();
    let lanes = partition_worker_lanes(conformers, results, num_threads)?;
    if lanes.is_empty() {
        return Ok(Vec::new());
    }

    let lane_count = NonZeroU32::new(
        u32::try_from(num_threads).expect("a positive signed worker count fits u32"),
    )
    .expect("positive source worker counts create nonempty lane sets");
    #[cfg(all(test, not(target_family = "wasm")))]
    let m10_observer = UFF_MT_M10_OBSERVER.with(|current| current.borrow().clone());
    Ok(std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(lanes.len());
        for (lane_index, lane) in lanes.into_iter().enumerate() {
            #[cfg(all(test, not(target_family = "wasm")))]
            let lane_observer = m10_observer.clone();
            let lane_field: ForceField<'rows> = source_field.copy();
            #[cfg(all(test, not(target_family = "wasm")))]
            let lane_entry = lane_observer.as_ref().map(|_| UffMtM10LaneEntry {
                lane_index,
                row_len: lane.rows.len(),
                row_capacity: lane.rows.capacity(),
                result_slot_len: lane.source_result_slots.len(),
                result_slot_capacity: lane.source_result_slots.capacity(),
                copied_field_state: crate::kernel::cf3d_uff_one_copy_state_for_worker_test(
                    &lane_field,
                ),
            });
            handles.push(scope.spawn(move || {
                #[cfg(all(test, not(target_family = "wasm")))]
                if let (Some(observer), Some(entry)) = (lane_observer.as_ref(), lane_entry) {
                    observer.enter_lane(entry);
                }
                let result = optimize_partitioned_worker(
                    lane_field,
                    lane,
                    conformer_count,
                    result_slot_count,
                    atom_count,
                    lane_count,
                    max_iterations,
                );
                #[cfg(all(test, not(target_family = "wasm")))]
                if let Some(observer) = lane_observer.as_ref() {
                    observer.note_worker_work(lane_index);
                }
                result
            }));
        }
        #[cfg(all(test, not(target_family = "wasm")))]
        if let Some(observer) = m10_observer.as_ref() {
            observer.note_handle_storage(handles.len(), handles.capacity());
        }
        join_worker_handles(handles)
    }))
}

#[cfg(not(target_family = "wasm"))]
pub(crate) fn optimize_prepared_conformers_dispatch<'field, 'rows, 'coordinates>(
    field: ForceField<'field>,
    conformers: &'rows mut [SerialConformer<'coordinates>],
    results: &'rows mut Vec<OptimizationOutcome>,
    atom_count: usize,
    max_iterations: i32,
    requested_threads: i32,
    observed_hardware: u32,
    threadsafe: bool,
) -> Result<PreparedConformerDispatchOutcome, UffThreadCountError>
where
    'field: 'rows,
    'coordinates: 'rows,
{
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

    // Behavior review: resize precedes exactly one explicit-input resolution;
    // the resolved signed count selects the same serial or worker owner. The
    // carrier retains the complete typed serial result or raw worker joins.
    // Complexity review: this dispatcher adds constant work, no result/field
    // copy and no numerical loop. D01 performs the source resize in place;
    // D02 is allocation-free; existing ST/MT owners retain their own costs.
    resize_conformer_results(results, conformers.len());
    let route = resolve_conformer_dispatch(requested_threads, observed_hardware, threadsafe)?;
    Ok(optimize_prepared_conformers_dispatch_resolved(
        field,
        conformers,
        results,
        atom_count,
        max_iterations,
        route,
    ))
}

#[cfg(not(target_family = "wasm"))]
fn optimize_prepared_conformers_dispatch_resolved<'field, 'rows, 'coordinates>(
    field: ForceField<'field>,
    conformers: &'rows mut [SerialConformer<'coordinates>],
    results: &'rows mut [OptimizationOutcome],
    atom_count: usize,
    max_iterations: i32,
    route: ConformerDispatchRoute,
) -> PreparedConformerDispatchOutcome
where
    'field: 'rows,
    'coordinates: 'rows,
{
    // RDKit❗✔️:   if (numThreads == 1) {
    // RDKit❗✔️:     detail::OptimizeMoleculeConfsST(mol, ff, res, maxIters);
    // RDKit❗✔️:   }
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:   else {
    // RDKit❗✔️:     detail::OptimizeMoleculeConfsMT(mol, ff, res, numThreads, maxIters);
    // RDKit❗✔️:   }
    // RDKit❗✔️: #endif
    // FFConvenience.h:121-128. This helper receives the already resolved
    // source count and reuses the sole serial stage owner or the existing raw
    // worker join owner without changing their result/error representations.
    match route {
        ConformerDispatchRoute::Serial => {
            PreparedConformerDispatchOutcome::Serial(optimize_serial_conformers(
                field,
                conformers,
                results,
                atom_count,
                max_iterations as u32,
            ))
        }
        ConformerDispatchRoute::Workers { source_count } => {
            PreparedConformerDispatchOutcome::Workers(optimize_conformers_mt(
                &field,
                conformers,
                results,
                source_count,
                atom_count,
                max_iterations,
            ))
        }
    }
}

pub(crate) fn optimize_serial_conformers<'field, 'rows>(
    field: ForceField<'field>,
    conformers: &'rows mut [SerialConformer<'_>],
    results: &mut [OptimizationOutcome],
    atom_count: usize,
    max_iterations: u32,
) -> Result<(), SerialConformerOptimizationError>
where
    'field: 'rows,
{
    optimize_serial_conformers_iter(
        field,
        conformers.iter_mut().map(|conformer| SerialConformer {
            id: conformer.id,
            positions: &mut *conformer.positions,
        }),
        results,
        atom_count,
        max_iterations,
    )
}

#[cfg(test)]
std::thread_local! {
    // Probe only when the adapter-specific regression arms it. The first
    // value counts rows materialized by the coordinate adapter; the second
    // records that count on entry to the shared stage owner, exposing eager
    // metadata collection without claiming the optimizer performs no other
    // allocations.
    static UFF_PUBLIC_PERF_SERIAL_ADAPTER_PROBE: std::cell::Cell<Option<(usize, Option<usize>)>> =
        const { std::cell::Cell::new(None) };
}

#[cfg(test)]
fn uff_public_perf_serial_adapter_probe_start() {
    UFF_PUBLIC_PERF_SERIAL_ADAPTER_PROBE.with(|probe| probe.set(Some((0, None))));
}

#[cfg(test)]
fn uff_public_perf_serial_adapter_probe_row() {
    UFF_PUBLIC_PERF_SERIAL_ADAPTER_PROBE.with(|probe| {
        if let Some((rows, stage_entry_rows)) = probe.get() {
            probe.set(Some((rows + 1, stage_entry_rows)));
        }
    });
}

#[cfg(test)]
fn uff_public_perf_serial_adapter_probe_stage_entry() {
    UFF_PUBLIC_PERF_SERIAL_ADAPTER_PROBE.with(|probe| {
        if let Some((rows, _)) = probe.get() {
            probe.set(Some((rows, Some(rows))));
        }
    });
}

#[cfg(test)]
fn uff_public_perf_serial_adapter_probe_finish() -> Option<(usize, Option<usize>)> {
    UFF_PUBLIC_PERF_SERIAL_ADAPTER_PROBE.with(|probe| {
        let value = probe.get();
        probe.set(None);
        value
    })
}

#[cfg(test)]
#[derive(Debug, Default)]
struct UffAllCostProbe {
    construction_calls: usize,
    rows_at_stage_entry: Option<usize>,
    visited_conformer_ids: Vec<usize>,
}

#[cfg(test)]
std::thread_local! {
    static UFF_ALL_COST_PROBE: std::cell::RefCell<Option<UffAllCostProbe>> =
        const { std::cell::RefCell::new(None) };
}

#[cfg(test)]
pub(super) fn uff_all_cost_probe_start() {
    UFF_ALL_COST_PROBE.with(|probe| {
        let previous = probe.borrow_mut().replace(UffAllCostProbe::default());
        assert!(previous.is_none(), "UFF all-cost probe was already active");
    });
}

#[cfg(test)]
fn uff_all_cost_note_construction_call() {
    UFF_ALL_COST_PROBE.with(|probe| {
        if let Some(probe) = probe.borrow_mut().as_mut() {
            probe.construction_calls += 1;
        }
    });
}

#[cfg(test)]
fn uff_all_cost_note_stage_entry() {
    UFF_ALL_COST_PROBE.with(|probe| {
        if let Some(probe) = probe.borrow_mut().as_mut() {
            probe.rows_at_stage_entry = Some(probe.visited_conformer_ids.len());
        }
    });
}

#[cfg(test)]
fn uff_all_cost_note_row(conformer_id: usize) {
    UFF_ALL_COST_PROBE.with(|probe| {
        if let Some(probe) = probe.borrow_mut().as_mut() {
            probe.visited_conformer_ids.push(conformer_id);
        }
    });
}

#[cfg(test)]
pub(super) fn uff_all_cost_probe_finish() -> Option<(usize, Option<usize>, Vec<usize>)> {
    UFF_ALL_COST_PROBE.with(|probe| {
        probe.borrow_mut().take().map(|probe| {
            (
                probe.construction_calls,
                probe.rows_at_stage_entry,
                probe.visited_conformer_ids,
            )
        })
    })
}

fn optimize_serial_conformers_iter<'field, 'rows, Conformers>(
    mut field: ForceField<'field>,
    conformers: Conformers,
    results: &mut [OptimizationOutcome],
    atom_count: usize,
    max_iterations: u32,
) -> Result<(), SerialConformerOptimizationError>
where
    'field: 'rows,
    Conformers: ExactSizeIterator<Item = SerialConformer<'rows>>,
{
    #[cfg(test)]
    uff_public_perf_serial_adapter_probe_stage_entry();
    #[cfg(test)]
    uff_all_cost_note_stage_entry();

    // RDKit✔️✔️: PRECONDITION(res.size() >= mol.getNumConformers(),
    // RDKit✔️✔️:              "res.size() must be >= mol.getNumConformers()");
    // RDKit✔️✔️: for (ROMol::ConformerIterator cit = mol.beginConformers();
    // RDKit✔️✔️:      cit != mol.endConformers(); ++cit, ++i) {
    // RDKit✔️✔️:   for (unsigned int aidx = 0; aidx < mol.getNumAtoms(); ++aidx) {
    // RDKit✔️✔️:     ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ff.initialize();
    // RDKit✔️✔️:   int needsMore = ff.minimize(maxIters);
    // RDKit✔️✔️:   double e = ff.calcEnergy();
    // RDKit✔️✔️:   res[i] = std::make_pair(needsMore, e);
    // RDKit✔️✔️: }
    // FFConvenience.h:66-77. Iterator order is the explicit source row order;
    // IDs remain labels. The source capacity precondition precedes all slots.
    validate_result_capacity(conformers.len(), results.len())
        .map_err(SerialConformerOptimizationError::ResultCapacity)?;

    // Consume the construction-reference borrows while retaining the existing
    // position vector allocation, then shorten its lifetime to this ordered
    // input iterator. The pinned source overwrites these same vector slots.
    let initial_position_buffer_address = field.positions().as_ptr() as usize;
    let initial_position_buffer_capacity = field.positions_mut().capacity();
    #[cfg(test)]
    crate::kernel::cf3d_uff_one_note_serial_position_buffer_start(
        initial_position_buffer_address,
        initial_position_buffer_capacity,
    );
    let mut field = field.clear_positions_and_shorten_lifetime::<'rows>();
    if field.positions_mut().capacity() < atom_count {
        #[cfg(test)]
        crate::kernel::cf3d_uff_one_note_serial_position_buffer_growth();
        field.positions_mut().reserve(atom_count);
    }
    for (input_index, conformer) in conformers.enumerate() {
        let SerialConformer {
            id: conformer_id,
            positions,
        } = conformer;
        if atom_count != 0 && positions.len() != atom_count {
            return Err(SerialConformerOptimizationError::CoordinateCount {
                input_index,
                source: SerialCoordinateCountError {
                    conformer_id,
                    atoms: atom_count,
                    coordinates: positions.len(),
                },
            });
        }
        let position_handles = field.positions_mut();
        position_handles.clear();
        position_handles.extend(
            positions
                .into_iter()
                .take(atom_count)
                .map(|position| &mut position[..]),
        );
        let slot_result =
            assign_serial_result(&mut field, max_iterations, &mut results[input_index]);
        if let Err(source) = slot_result {
            return Err(SerialConformerOptimizationError::Optimization {
                input_index,
                conformer_id,
                source,
            });
        }
    }
    #[cfg(test)]
    crate::kernel::cf3d_uff_one_note_serial_position_buffer_end(
        field.positions().as_ptr() as usize,
        field.positions_mut().capacity(),
    );
    Ok(())
}

fn optimize_worker_conformers<'field, 'rows>(
    field: ForceField<'field>,
    conformers: &'rows mut [SerialConformer<'_>],
    results: &mut [OptimizationOutcome],
    atom_count: usize,
    thread_idx: u32,
    lane_count: NonZeroU32,
    max_iterations: i32,
) -> Result<(), SerialConformerOptimizationError>
where
    'field: 'rows,
{
    // RDKit❗✔️: inline void OptimizeMoleculeConfsHelper_(
    // RDKit❗✔️:     ForceFields::ForceField ff, ROMol *mol,
    // RDKit❗✔️:     std::vector<std::pair<int, double>> *res, unsigned int threadIdx,
    // RDKit❗✔️:     unsigned int numThreads, int maxIters) {
    // RDKit❗✔️:   PRECONDITION(mol, "mol must not be nullptr");
    // RDKit❗✔️:   PRECONDITION(res, "res must not be nullptr");
    // RDKit❗✔️:   PRECONDITION(res->size() >= mol->getNumConformers(),
    // RDKit❗✔️:                "res->size() must be >= mol->getNumConformers()");
    // RDKit❗✔️:   unsigned int i = 0;
    // RDKit❗✔️:   ff.positions().resize(mol->getNumAtoms());
    // RDKit❗✔️:   for (ROMol::ConformerIterator cit = mol->beginConformers();
    // RDKit❗✔️:        cit != mol->endConformers(); ++cit, ++i) {
    // RDKit❗✔️:     if (i % numThreads != threadIdx) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (unsigned int aidx = 0; aidx < mol->getNumAtoms(); ++aidx) {
    // RDKit❗✔️:       ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ff.initialize();
    // RDKit❗✔️:     int needsMore = ff.minimize(maxIters);
    // RDKit❗✔️:     double e = ff.calcEnergy();
    // RDKit❗✔️:     (*res)[i] = std::make_pair(needsMore, e);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // FFConvenience.h:21-44. This is one borrowed synchronous worker; the
    // source MT creation/join layer is intentionally outside this function.
    // The assignment counter is source unsigned `i`, not the row ID; the
    // source-order iterator advances it for both selected and skipped rows.
    // M04 delegates preparation and the complete initialize/minimize/final
    // energy/result sequence to its shared loop, retaining source counts and
    // wrapped result identity without a second traversal-stage implementation.
    let conformer_count = conformers.len();
    let result_slot_count = results.len();
    let mut assignment = WorkerConformerAssignment::new(thread_idx, lane_count);
    let rows = conformers.iter_mut().filter_map(|conformer| {
        let source_index = assignment.source_index;
        let row = if assignment.selects_current() {
            Some(SourceIndexedWorkerRow {
                source_index,
                conformer_id: conformer.id,
                conformer,
            })
        } else {
            None
        };
        assignment.advance();
        row
    });
    optimize_assigned_worker_rows(
        field,
        rows,
        WorkerResultSlotView::Dense(results),
        conformer_count,
        result_slot_count,
        atom_count,
        max_iterations,
    )
}

fn execute_serial_slot(
    field: &mut ForceField<'_>,
    max_iterations: u32,
) -> Result<OptimizationOutcome, OptimizationStageError> {
    // RDKit✔️✔️: ff.initialize();
    // RDKit✔️✔️: int needsMore = ff.minimize(maxIters);
    // RDKit✔️✔️: double e = ff.calcEnergy();
    // FFConvenience.h:74-76 has exactly the same stages/default tolerances
    // as OptimizeMolecule. Reuse its complete driver and typed stage errors;
    // status1 remains an ordinary successful outcome with gathered positions.
    // Delegation introduces no field, contribution or coordinate copy.
    optimize_force_field(field, max_iterations)
}

pub(super) fn optimize_force_field(
    force_field: &mut ForceField<'_>,
    max_iterations: u32,
) -> Result<OptimizationOutcome, OptimizationStageError> {
    // BEGIN RDKIT CPP FUNCTION ForceFieldHelpers::OptimizeMolecule (FFConvenience.h:94-99)
    // RDKit✔️✔️: inline std::pair<int, double> OptimizeMolecule(ForceFields::ForceField &ff,
    // RDKit✔️✔️:                                                int maxIters = 1000) {
    // RDKit✔️✔️:   ff.initialize();
    // RDKit✔️✔️:   int res = ff.minimize(maxIters);
    // RDKit✔️✔️:   double e = ff.calcEnergy();
    // RDKit✔️✔️:   return std::make_pair(res, e);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ForceFieldHelpers::OptimizeMolecule
    // ForceField::minimize(maxIts) inherits these pinned default tolerances;
    // the selected path supplies its signed max-iteration conversion in U05.
    // RDKit✔️✔️: int minimize(unsigned int maxIts = 200, double forceTol = 1e-4,
    // RDKit✔️✔️:              double energyTol = 1e-6);
    force_field
        .initialize()
        .map_err(OptimizationStageError::Initialize)?;
    let status = force_field
        .minimize(max_iterations, 1e-4, 1e-6)
        .map_err(OptimizationStageError::Minimize)?;
    let energy = force_field
        .calc_energy_current(None)
        .map_err(OptimizationStageError::FinalEnergy)?;

    Ok(OptimizationOutcome { status, energy })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn optimize_single_conformer_with_state(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    typing_state: UffAtomStateRef<'_>,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
) -> Result<OptimizationOutcome, SingleConformerOptimizationError> {
    // BEGIN RDKIT CPP FUNCTION UFF::UFFOptimizeMolecule (UFF.h:40-47)
    // RDKit❗✔️: inline std::pair<int, double> UFFOptimizeMolecule(
    // RDKit❗✔️:     ROMol &mol, int maxIters = 1000, double vdwThresh = 10.0, int confId = -1,
    // RDKit❗✔️:     bool ignoreInterfragInteractions = true) {
    // RDKit❗✔️:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit❗✔️:       mol, vdwThresh, confId, ignoreInterfragInteractions));
    // RDKit❗✔️:   std::pair<int, double> res =
    // RDKit❗✔️:       ForceFieldsHelper::OptimizeMolecule(*ff, maxIters);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::UFFOptimizeMolecule
    let mut force_field = super::builder::construct_force_field_with_automatic_typing(
        topology,
        coordinates,
        options.conformer_id,
        typing_state,
        rings,
        valence,
        molecule_properties,
        diagnostics,
        DEFAULT_TORSION_BOND_SMARTS,
        options.vdw_threshold,
        options.ignore_interfragment_interactions,
    )
    .map_err(SingleConformerOptimizationError::Construction)?;

    let outcome = optimize_force_field(&mut force_field, options.max_iterations_as_unsigned())
        .map_err(SingleConformerOptimizationError::Optimization)?;

    // Behavior marker — RDKit❗✔️: the modeled entry accepts explicit 3D IDs
    // and caller-prepared chemistry state, then preserves the source builder
    // before initialize/minimize/final-energy order and returns its status and
    // energy. Constructor and kernel stage causes retain their original enums.
    // The source's negative conformer-ID default is intentionally outside this
    // explicit selected-conformer input contract.
    // Complexity marker — RDKit❗❌: this entry invokes automatic typing and
    // the source-ordered field builder once, borrowing topology, prepared
    // chemistry, and one mutable coordinate-row reference per selected atom.
    // The builder allocates one parameter row vector and builds its T accepted
    // contribution objects once (one boxed kernel term per stored term); the
    // optimizer reuses that same field and term vector in every callback. The
    // optimizer creates one D-coordinate point buffer, five additional
    // D-coordinate work vectors, one D-by-D inverse-Hessian vector, and the
    // field initializes one triangular distance cache. Final energy creates
    // one D-coordinate scatter buffer. Energy/gradient callbacks borrow these
    // arrays and the stored terms; they do not copy the field or rebuild terms.
    // The torsion caller now requests the atom-only Search projection: neither
    // top-level nor recursive hits materialize the unused bond-mapping vector.
    // Search still retains separate raw query/target row vectors per hit and
    // constructs a query-indexed optional mapping before the output atom row;
    // this differs from RDKit's raw pair row and reusable matchVect buffer.
    // The remaining result-buffer overhead keeps the performance gap explicit;
    // it is not attributed to optimizer callbacks. The source-required fragment
    // values and selected conformer metadata remain when that option is active.
    // The selected entry accepts no live Molecule and makes no direct clone
    // of its input blocks. Its source-required nonbonded fragment path still
    // materializes copied fragment values when interfragment interactions are
    // ignored; the pinned single-component getTheFrags branch copies the full
    // component.
    Ok(outcome)
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(super) fn optimize_single_conformer(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
) -> Result<OptimizationOutcome, SingleConformerOptimizationError> {
    optimize_single_conformer_with_state(
        topology,
        coordinates,
        UffAtomStateRef::SuppliedRows {
            total_valences,
            conjugated_presence: atom_has_conjugated_bond,
        },
        rings,
        valence,
        molecule_properties,
        diagnostics,
        options,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn optimize_serial_uff(
    topology: &TopologyBlock,
    reference_coordinates: &mut CoordinateBlock,
    conformers: &mut [SerialConformer<'_>],
    results: &mut [OptimizationOutcome],
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
) -> Result<(), SerialUffOptimizationError> {
    // BEGIN RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs (UFF.h:69-78)
    // RDKit✔️❌: inline void UFFOptimizeMoleculeConfs(ROMol &mol,
    // RDKit✔️❌:                                     std::vector<std::pair<int, double>> &res,
    // RDKit✔️❌:                                     int numThreads = 1, int maxIters = 1000,
    // RDKit✔️❌:                                     double vdwThresh = 10.0,
    // RDKit✔️❌:                                     bool ignoreInterfragInteractions = true) {
    // RDKit✔️❌:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit✔️❌:       mol, vdwThresh, -1, ignoreInterfragInteractions));
    // RDKit✔️❌:   ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads, maxIters);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs
    // Source constructs once before delegated result-capacity validation and
    // traversal. This private path receives the chosen 3D construction
    // reference explicitly; supplied mutable conformer rows define traversal
    // order. The negative source default and thread dispatch remain outside.
    let field = construct_force_field_with_automatic_typing(
        topology,
        reference_coordinates,
        options.conformer_id,
        total_valences,
        atom_has_conjugated_bond,
        rings,
        valence,
        molecule_properties,
        diagnostics,
        DEFAULT_TORSION_BOND_SMARTS,
        options.vdw_threshold,
        options.ignore_interfragment_interactions,
    )
    .map_err(SerialUffOptimizationError::Construction)?;

    optimize_serial_conformers(
        field,
        conformers,
        results,
        topology.atoms.len(),
        options.max_iterations_as_unsigned(),
    )
    .map_err(SerialUffOptimizationError::Optimization)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn optimize_serial_uff_coordinate_block(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    results: &mut Vec<OptimizationOutcome>,
    typing_state: UffAtomStateRef<'_>,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
) -> Result<(), SerialUffOptimizationError> {
    // BEGIN RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs (UFF.h:69-78)
    // RDKit❗❌: inline void UFFOptimizeMoleculeConfs(ROMol &mol,
    // RDKit❗❌:                                     std::vector<std::pair<int, double>> &res,
    // RDKit❗❌:                                     int numThreads = 1, int maxIters = 1000,
    // RDKit❗❌:                                     double vdwThresh = 10.0,
    // RDKit❗❌:                                     bool ignoreInterfragInteractions = true) {
    // RDKit❗❌:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit❗❌:       mol, vdwThresh, -1, ignoreInterfragInteractions));
    // RDKit❗❌:   ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads, maxIters);
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs
    // This selected-ID boundary reuses the source builder and serial stage
    // owner. Construction finishes and releases its reference-row borrows
    // before the same block's 3D rows are mutably borrowed again. The source
    // confId=-1 default remains outside this explicit input contract.
    // RDKit✔️✔️:   res.resize(mol.getNumConformers());
    // FFConvenience.h:115-129. Match the source boundary order: construct,
    // resize results, then enumerate and process source-ordered conformers.
    #[cfg(test)]
    uff_all_cost_note_construction_call();
    let field = construct_released_uff_field(
        topology,
        coordinates,
        typing_state,
        rings,
        valence,
        molecule_properties,
        diagnostics,
        options,
    )
    .map_err(SerialUffOptimizationError::Construction)?;

    resize_conformer_results(results, coordinates.conformers_3d.len());
    let conformers = coordinates.conformers_3d.iter_mut().map(|conformer| {
        #[cfg(test)]
        uff_public_perf_serial_adapter_probe_row();
        #[cfg(test)]
        uff_all_cost_note_row(conformer.id());
        SerialConformer {
            id: conformer.id(),
            positions: conformer.coordinates_mut(),
        }
    });

    // Complexity: the exact-size adapter visits each 3D row once and removes
    // the O(C) temporary borrowed-row Vec; the one O(A) position-handle
    // buffer and numerical stage loop are reused. No coordinate or
    // contribution is cloned, and 2D rows remain untouched.
    optimize_serial_conformers_iter(
        field,
        conformers,
        results,
        topology.atoms.len(),
        options.max_iterations_as_unsigned(),
    )
    .map_err(SerialUffOptimizationError::Optimization)
}

#[cfg(not(target_family = "wasm"))]
#[allow(clippy::too_many_arguments)]
fn optimize_dispatched_uff<'field, 'rows, 'positions>(
    topology: &TopologyBlock,
    reference_coordinates: &'field mut CoordinateBlock,
    conformers: &'rows mut [SerialConformer<'positions>],
    results: &'rows mut Vec<OptimizationOutcome>,
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
    requested_threads: i32,
    observed_hardware: u32,
    threadsafe: bool,
) -> Result<PreparedConformerDispatchOutcome, DispatchedUffOptimizationError>
where
    'field: 'rows,
    'positions: 'rows,
{
    // BEGIN RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs (UFF.h:69-78)
    // RDKit✔️❌: inline void UFFOptimizeMoleculeConfs(ROMol &mol,
    // RDKit✔️❌:                                     std::vector<std::pair<int, double>> &res,
    // RDKit✔️❌:                                     int numThreads = 1, int maxIters = 1000,
    // RDKit✔️❌:                                     double vdwThresh = 10.0,
    // RDKit✔️❌:                                     bool ignoreInterfragInteractions = true) {
    // RDKit✔️❌:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit✔️❌:       mol, vdwThresh, -1, ignoreInterfragInteractions));
    // RDKit✔️❌:   ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads, maxIters);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs
    // Source constructs one field before delegated result resizing and thread
    // resolution. This private entry takes the selected construction reference
    // and prepared chemistry explicitly; it does not select confId=-1 or
    // recompute those inputs. The existing builder owns the single construction
    // pass, and D04 owns resize, thread resolution and ST/MT stage delegation.
    // Complexity: this wrapper adds constant dispatch work and no field,
    // coordinate-block or contribution copy of its own. Construction allocates
    // the same parameter and term state as the delegated builder; MT retains
    // the already-reviewed per-lane source-shaped ForceField copies.
    let field = construct_force_field_with_automatic_typing(
        topology,
        reference_coordinates,
        options.conformer_id,
        total_valences,
        atom_has_conjugated_bond,
        rings,
        valence,
        molecule_properties,
        diagnostics,
        DEFAULT_TORSION_BOND_SMARTS,
        options.vdw_threshold,
        options.ignore_interfragment_interactions,
    )
    .map_err(DispatchedUffOptimizationError::Construction)?;

    optimize_prepared_conformers_dispatch(
        field,
        conformers,
        results,
        topology.atoms.len(),
        options.max_iterations,
        requested_threads,
        observed_hardware,
        threadsafe,
    )
    .map_err(DispatchedUffOptimizationError::ThreadCount)
}

#[cfg(not(target_family = "wasm"))]
#[allow(clippy::too_many_arguments)]
pub(super) fn optimize_dispatched_uff_coordinate_block<'rows>(
    topology: &TopologyBlock,
    coordinates: &'rows mut CoordinateBlock,
    results: &'rows mut Vec<OptimizationOutcome>,
    typing_state: UffAtomStateRef<'_>,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
    requested_threads: i32,
    observed_hardware: u32,
    threadsafe: bool,
) -> Result<PreparedConformerDispatchOutcome, DispatchedUffOptimizationError> {
    // BEGIN RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs (UFF.h:69-80)
    // RDKit✔️❌: inline void UFFOptimizeMoleculeConfs(ROMol &mol,
    // RDKit✔️❌:                                     std::vector<std::pair<int, double>> &res,
    // RDKit✔️❌:                                     int numThreads = 1, int maxIters = 1000,
    // RDKit✔️❌:                                     double vdwThresh = 10.0,
    // RDKit✔️❌:                                     bool ignoreInterfragInteractions = true) {
    // RDKit✔️❌:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit✔️❌:       mol, vdwThresh, -1, ignoreInterfragInteractions));
    // RDKit✔️❌:   ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads, maxIters);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs
    // Behavior review: the private input contract selects the construction
    // reference by explicit 3D ID, uses caller-prepared chemistry, and releases
    // that field's position borrows before borrowing every 3D row from this
    // same CoordinateBlock in stored order. The existing prepared dispatcher
    // retains source result-resize, explicit thread-resolution, and serial or
    // raw worker-outcome behavior. The source confId=-1 choice remains outside
    // this selected-ID boundary; 2D rows are never processed.
    // Complexity review: automatic construction runs once. The serial route
    // passes an O(1)-state exact-size iterator to the existing stage owner and
    // creates no borrowed-row Vec. The worker route retains one O(C) borrowed
    // SerialConformer Vec for the existing partitioner, without field, term,
    // or coordinate copies. No host query or error/panic aggregation is added.
    let field = construct_released_uff_field(
        topology,
        coordinates,
        typing_state,
        rings,
        valence,
        molecule_properties,
        diagnostics,
        options,
    )
    .map_err(DispatchedUffOptimizationError::Construction)?;

    // Match the source's construct, resize, resolve, then select-route order.
    // Keep the adapter lazy until the route is known: serial consumes it
    // directly; only the existing worker partition owner needs a row slice.
    resize_conformer_results(results, coordinates.conformers_3d.len());
    let route = resolve_conformer_dispatch(requested_threads, observed_hardware, threadsafe)
        .map_err(DispatchedUffOptimizationError::ThreadCount)?;
    let conformers = coordinates.conformers_3d.iter_mut().map(|conformer| {
        #[cfg(test)]
        uff_public_perf_serial_adapter_probe_row();
        SerialConformer {
            id: conformer.id(),
            positions: conformer.coordinates_mut(),
        }
    });

    match route {
        ConformerDispatchRoute::Serial => Ok(PreparedConformerDispatchOutcome::Serial(
            optimize_serial_conformers_iter(
                field,
                conformers,
                results,
                topology.atoms.len(),
                options.max_iterations_as_unsigned(),
            ),
        )),
        ConformerDispatchRoute::Workers { .. } => {
            let mut conformers = conformers.collect::<Vec<_>>();
            Ok(optimize_prepared_conformers_dispatch_resolved(
                field,
                &mut conformers,
                results,
                topology.atoms.len(),
                options.max_iterations,
                route,
            ))
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn construct_released_uff_field<'released>(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    typing_state: UffAtomStateRef<'_>,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    options: SingleConformerOptions,
) -> Result<ForceField<'released>, AutomaticForceFieldConstructionError> {
    // BEGIN RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs (UFF.h:78-80)
    // RDKit❗✔️:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit❗✔️:       mol, vdwThresh, -1, ignoreInterfragInteractions));
    // END RDKIT CPP FUNCTION UFF::UFFOptimizeMoleculeConfs construction
    // The pinned multi-conformer entry passes source confId -1; this private
    // boundary receives the caller's explicit selected 3D ID in the same
    // constructor slot. It reuses the automatic-typing constructor exactly
    // once, retaining diagnostic order and its typed construction error,
    // then consumes that same field to release its old mutable coordinate
    // borrows before the caller reborrows this block.
    // Complexity: the existing field, terms, coordinate block and position
    // Vec allocation are moved; no chemistry state or coordinates are copied.
    let field = super::builder::construct_force_field_with_automatic_typing(
        topology,
        coordinates,
        options.conformer_id,
        typing_state,
        rings,
        valence,
        molecule_properties,
        diagnostics,
        DEFAULT_TORSION_BOND_SMARTS,
        options.vdw_threshold,
        options.ignore_interfragment_interactions,
    )?;

    Ok(field.release_position_borrows())
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::error::Error as _;

    use cosmolkit_core::ValenceAssignment;
    use cosmolkit_core::fast_find_rings;
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer3D, CoordinateBlock,
        Element, Hybridization, MoleculeProperties, TopologyBlock,
    };

    #[test]
    fn uff_integrate_i02_outer_causes_borrow_actual_children() {
        use super::super::optimization::UffPreparedOptimizationError;

        fn assert_stored_child<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let reported = std::error::Error::source(parent)
                .expect("the outer error exposes its concrete stored child");
            let stored_error: &(dyn std::error::Error + 'static) = stored;
            assert!(std::ptr::eq(reported, stored_error));
            assert!(std::ptr::eq(
                reported
                    .downcast_ref::<Child>()
                    .expect("downcast keeps child type"),
                stored
            ));
            assert_eq!(parent.to_string(), stored.to_string());
        }

        // These outer wrappers are supplementary trait-dispatch checks. Their
        // children come from actual preparation/owner failures; I03-I05 also
        // exercise prepared single/serial/dispatch forwarding. The native raw
        // Serial/Workers outcome remains deliberately unconverted here.
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );

        let invalid_assignment = ValenceAssignment {
            explicit_valence: vec![0],
            implicit_hydrogens: vec![0, 0],
        };
        let preparation =
            super::super::api::prepare_parameter_query(&topology, &invalid_assignment)
                .err()
                .expect("real shared preparation reports the stored row mismatch");
        let wrapped = UffPreparedOptimizationError::Preparation(preparation);
        let stored = match &wrapped {
            UffPreparedOptimizationError::Preparation(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&wrapped, stored);

        let total_valences = [0, 0];
        let conjugated = [false, false];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring data");
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let mut options = super::SingleConformerOptions::for_conformer(99);
        options.ignore_interfragment_interactions = false;
        let mut actual_owner_error_calls = 1;

        let mut single_coordinates = worker_w06_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let single = super::optimize_single_conformer(
            &topology,
            &mut single_coordinates,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut Vec::new(),
            options,
        )
        .err()
        .expect("the real selected-ID constructor returns a typed failure");
        actual_owner_error_calls += 1;
        assert!(matches!(
            &single,
            super::SingleConformerOptimizationError::Construction(_)
        ));
        let wrapped = UffPreparedOptimizationError::Single(single);
        let stored = match &wrapped {
            UffPreparedOptimizationError::Single(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&wrapped, stored);

        let sentinel = super::OptimizationOutcome {
            status: -7,
            energy: 29.0,
        };
        let mut serial_coordinates = worker_w06_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let mut serial_results = vec![sentinel];
        let serial = super::optimize_serial_uff_coordinate_block(
            &topology,
            &mut serial_coordinates,
            &mut serial_results,
            super::UffAtomStateRef::SuppliedRows {
                total_valences: &total_valences,
                conjugated_presence: &conjugated,
            },
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut Vec::new(),
            options,
        )
        .err()
        .expect("the real serial entry preserves its construction failure");
        actual_owner_error_calls += 1;
        assert_eq!(serial_results, [sentinel]);
        assert!(matches!(
            &serial,
            super::SerialUffOptimizationError::Construction(_)
        ));
        let wrapped = UffPreparedOptimizationError::Serial(serial);
        let stored = match &wrapped {
            UffPreparedOptimizationError::Serial(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&wrapped, stored);

        #[cfg(not(target_family = "wasm"))]
        {
            let mut dispatch_coordinates =
                worker_w06_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
            let mut dispatch_results = vec![sentinel];
            let dispatch = super::optimize_dispatched_uff_coordinate_block(
                &topology,
                &mut dispatch_coordinates,
                &mut dispatch_results,
                super::UffAtomStateRef::SuppliedRows {
                    total_valences: &total_valences,
                    conjugated_presence: &conjugated,
                },
                &rings,
                &valence,
                &MoleculeProperties::default(),
                &mut Vec::new(),
                options,
                1,
                2,
                false,
            )
            .err()
            .expect("the real dispatch entry preserves construction errors");
            actual_owner_error_calls += 1;
            assert_eq!(dispatch_results, [sentinel]);
            assert!(matches!(
                &dispatch,
                super::DispatchedUffOptimizationError::Construction(_)
            ));
            let wrapped = UffPreparedOptimizationError::Dispatch(dispatch);
            let stored = match &wrapped {
                UffPreparedOptimizationError::Dispatch(stored) => stored,
                _ => unreachable!(),
            };
            assert_stored_child(&wrapped, stored);

            let mut dispatch_coordinates =
                worker_w06_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
            let mut dispatch_results = vec![sentinel];
            let dispatch = super::optimize_dispatched_uff_coordinate_block(
                &topology,
                &mut dispatch_coordinates,
                &mut dispatch_results,
                super::UffAtomStateRef::SuppliedRows {
                    total_valences: &total_valences,
                    conjugated_presence: &conjugated,
                },
                &rings,
                &valence,
                &MoleculeProperties::default(),
                &mut Vec::new(),
                super::SingleConformerOptions::for_conformer(31),
                i32::MIN,
                0,
                true,
            )
            .err()
            .expect("the real native dispatcher preserves thread-count failure");
            actual_owner_error_calls += 1;
            assert!(matches!(
                &dispatch,
                super::DispatchedUffOptimizationError::ThreadCount(_)
            ));
            let wrapped = UffPreparedOptimizationError::Dispatch(dispatch);
            let stored = match &wrapped {
                UffPreparedOptimizationError::Dispatch(stored) => stored,
                _ => unreachable!(),
            };
            assert_stored_child(&wrapped, stored);
        }

        #[cfg(not(target_family = "wasm"))]
        assert_eq!(actual_owner_error_calls, 5);
        #[cfg(target_family = "wasm")]
        assert_eq!(actual_owner_error_calls, 3);
    }

    #[test]
    fn uff_integrate_i03_prepared_single_matches_fixed_w06_source_matrix() {
        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let mut actual_calls = 0;

        for row_count in 1..=IDS.len() {
            for (selected_row, selected_id) in IDS[..row_count].iter().copied().enumerate() {
                for max_iterations in [0, 1] {
                    let mut coordinates = CoordinateBlock::default();
                    for row in 0..row_count {
                        coordinates.conformers_3d.push(
                            Conformer3D::new(
                                IDS[row],
                                [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]].to_vec(),
                                true,
                            )
                            .with_prop("source-row", format!("row-{row}")),
                        );
                    }
                    coordinates.conformers_2d.push(
                        cosmolkit_model::Conformer2D::new(selected_id, TWO_D_SENTINEL.to_vec())
                            .with_prop("layout", "preserve-this-overlapping-2d-row"),
                    );
                    let before_3d = coordinates.conformers_3d.clone();
                    let before_2d = coordinates.conformers_2d[0].clone();
                    let before_dimension = coordinates.source_coordinate_dim;
                    let mut diagnostics = Vec::new();
                    let mut options = super::SingleConformerOptions::for_conformer(selected_id);
                    options.max_iterations = max_iterations;
                    options.ignore_interfragment_interactions = false;

                    let outcome = super::super::optimization::optimize_prepared_uff_single(
                        &topology,
                        &mut coordinates,
                        &valence,
                        &rings,
                        &properties,
                        &mut diagnostics,
                        options,
                    )
                    .unwrap_or_else(|error| {
                        panic!("fixed W06 prepared single succeeds: rows={row_count}, selected={selected_id}, iterations={max_iterations}, error={error:?}")
                    });
                    actual_calls += 1;

                    assert_eq!(outcome.status, 1);
                    let (expected_xs, expected_energy) = if max_iterations == 0 {
                        (
                            [0.0, TARGET_DISTANCES[selected_row]],
                            source_energy(TARGET_DISTANCES[selected_row]),
                        )
                    } else {
                        source_one_step(TARGET_DISTANCES[selected_row])
                    };
                    assert!(
                        (outcome.energy - expected_energy).abs() <= 1.0e-12,
                        "fixed W06 source energy selected={selected_id}, iterations={max_iterations}"
                    );

                    assert_eq!(coordinates.conformers_3d.len(), row_count);
                    assert_eq!(
                        coordinates
                            .conformers_3d
                            .iter()
                            .map(Conformer3D::id)
                            .collect::<Vec<_>>(),
                        &IDS[..row_count]
                    );
                    for row in 0..row_count {
                        assert_eq!(
                            coordinates.conformers_3d[row].is_3d(),
                            before_3d[row].is_3d()
                        );
                        assert_eq!(
                            coordinates.conformers_3d[row].props(),
                            before_3d[row].props()
                        );
                        let expected_x = if row != selected_row || max_iterations == 0 {
                            [0.0, TARGET_DISTANCES[row]]
                        } else {
                            expected_xs
                        };
                        let expected_points =
                            [[expected_x[0], 0.0, 0.0], [expected_x[1], 0.0, 0.0]];
                        for (point, expected_point) in expected_points.into_iter().enumerate() {
                            for (axis, expected_component) in expected_point.into_iter().enumerate()
                            {
                                assert_eq!(
                                    coordinates.conformers_3d[row].coordinates()[point][axis]
                                        .to_bits(),
                                    expected_component.to_bits(),
                                    "fixed W06 selected-only coordinate row={row}, point={point}, axis={axis}, selected={selected_id}, iterations={max_iterations}"
                                );
                            }
                        }
                    }
                    assert_eq!(coordinates.conformers_2d[0], before_2d);
                    assert_eq!(coordinates.conformers_2d[0].id(), selected_id);
                    assert_eq!(coordinates.source_coordinate_dim, before_dimension);
                    assert!(diagnostics.is_empty());
                }
            }
        }

        assert_eq!(actual_calls, 12);
    }

    #[test]
    fn uff_prepare_p07_cached_and_supplied_single_match_fixed_w06_source_steps() {
        const DISTANCE: f64 = 4.0;
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let ratio = minimum / DISTANCE;
        let ratio3 = (ratio * ratio) * ratio;
        let ratio6 = ratio3 * ratio3;
        let ratio7 = (ratio3 * ratio3) * ratio;
        let ratio13 = (ratio6 * ratio6) * ratio;
        let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
        let first_gradient = (pre_factor * (0.0 - DISTANCE) / DISTANCE) * 0.1;
        let second_gradient = (pre_factor * (DISTANCE - 0.0) / DISTANCE) * 0.1;
        let expected_one_step = [
            0.0 + 1.0 * -first_gradient,
            DISTANCE + 1.0 * -second_gradient,
        ];
        let after_dx = expected_one_step[0] - expected_one_step[1];
        let after_distance = (after_dx * after_dx).sqrt();
        let expected = [
            ([0.0, DISTANCE], source_energy(DISTANCE)),
            (expected_one_step, source_energy(after_distance)),
        ];
        let supplied_valences = [0, 0];
        let supplied_conjugation = [false; 2];
        let mut actual_calls = 0;

        for cached in [true, false] {
            for (iteration_index, max_iterations) in [0, 1].into_iter().enumerate() {
                let mut coordinates =
                    worker_w06_coordinates(&[[0.0, 0.0, 0.0], [DISTANCE, 0.0, 0.0]]);
                let original_other_axes = coordinates.conformers_3d[0].coordinates().to_vec();
                let mut diagnostics = Vec::new();
                let mut options = super::SingleConformerOptions::for_conformer(31);
                options.max_iterations = max_iterations;
                options.ignore_interfragment_interactions = false;

                let outcome = if cached {
                    super::super::optimization::optimize_prepared_uff_single(
                        &topology,
                        &mut coordinates,
                        &valence,
                        &rings,
                        &properties,
                        &mut diagnostics,
                        options,
                    )
                    .unwrap_or_else(|error| {
                        panic!("fixed cached W06 p07 route succeeds: {error:?}")
                    })
                } else {
                    super::optimize_single_conformer(
                        &topology,
                        &mut coordinates,
                        &supplied_valences,
                        &supplied_conjugation,
                        &rings,
                        &valence,
                        &properties,
                        &mut diagnostics,
                        options,
                    )
                    .unwrap_or_else(|error| {
                        panic!("fixed supplied-row W06 p07 route succeeds: {error:?}")
                    })
                };
                actual_calls += 1;

                assert_eq!(
                    outcome.status, 1,
                    "cached={cached}, iterations={max_iterations}"
                );
                assert!(
                    (outcome.energy - expected[iteration_index].1).abs() <= 1.0e-12,
                    "fixed W06 source energy: cached={cached}, iterations={max_iterations}"
                );
                let actual = &coordinates.conformers_3d[0].coordinates();
                assert_eq!(
                    actual[0][0].to_bits(),
                    expected[iteration_index].0[0].to_bits()
                );
                assert_eq!(
                    actual[1][0].to_bits(),
                    expected[iteration_index].0[1].to_bits()
                );
                assert_eq!(actual[0][1].to_bits(), original_other_axes[0][1].to_bits());
                assert_eq!(actual[0][2].to_bits(), original_other_axes[0][2].to_bits());
                assert_eq!(actual[1][1].to_bits(), original_other_axes[1][1].to_bits());
                assert_eq!(actual[1][2].to_bits(), original_other_axes[1][2].to_bits());
                assert!(diagnostics.is_empty());
            }
        }
        assert_eq!(actual_calls, 4);
    }

    #[test]
    fn uff_prepare_p07_cached_and_supplied_single_keep_option_and_error_contracts() {
        use super::super::builder::{ForceFieldConstructionError, UffBuilderError};
        use super::super::optimization::UffPreparedOptimizationError;

        let defaults = super::SingleConformerOptions::for_conformer(31);
        assert_eq!(defaults.conformer_id, 31);
        assert_eq!(defaults.max_iterations, 1000);
        assert_eq!(defaults.vdw_threshold, 10.0);
        assert!(defaults.ignore_interfragment_interactions);
        for (signed, expected) in [
            (0, 0_u32),
            (23, 23_u32),
            (-1, u32::MAX),
            (-3, u32::MAX - 2),
            (i32::MIN, 1_u32 << 31),
        ] {
            let mut options = defaults;
            options.max_iterations = signed;
            assert_eq!(options.max_iterations_as_unsigned(), expected);
        }

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let properties = MoleculeProperties::default();
        let supplied_valences = [0, 0];
        let supplied_conjugation = [false; 2];

        for cached in [true, false] {
            let mut missing_id = worker_w06_coordinates(&[[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]]);
            let original_missing_id = missing_id.clone();
            let mut missing_options = super::SingleConformerOptions::for_conformer(99);
            missing_options.ignore_interfragment_interactions = false;
            let mut diagnostics = Vec::new();
            if cached {
                let error = super::super::optimization::optimize_prepared_uff_single(
                    &topology,
                    &mut missing_id,
                    &valence,
                    &rings,
                    &properties,
                    &mut diagnostics,
                    missing_options,
                )
                .expect_err("missing cached-route selected ID remains a typed construction error");
                assert!(matches!(
                    error,
                    UffPreparedOptimizationError::Single(
                        super::SingleConformerOptimizationError::Construction(
                            super::AutomaticForceFieldConstructionError::Construction(
                                ForceFieldConstructionError::Builder(
                                    UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                                        conformer_id: 99
                                    }
                                )
                            )
                        )
                    )
                ));
            } else {
                let error = super::optimize_single_conformer(
                    &topology,
                    &mut missing_id,
                    &supplied_valences,
                    &supplied_conjugation,
                    &rings,
                    &valence,
                    &properties,
                    &mut diagnostics,
                    missing_options,
                )
                .expect_err(
                    "missing supplied-route selected ID remains a typed construction error",
                );
                assert!(matches!(
                    error,
                    super::SingleConformerOptimizationError::Construction(
                        super::AutomaticForceFieldConstructionError::Construction(
                            ForceFieldConstructionError::Builder(
                                UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                                    conformer_id: 99
                                }
                            )
                        )
                    )
                ));
            }
            assert_eq!(missing_id, original_missing_id);

            let mut missing_coordinates = worker_w06_coordinates(&[[0.0, 0.0, 0.0]]);
            let original_missing_coordinates = missing_coordinates.clone();
            let mut options = super::SingleConformerOptions::for_conformer(31);
            options.ignore_interfragment_interactions = false;
            let mut diagnostics = Vec::new();
            if cached {
                let error = super::super::optimization::optimize_prepared_uff_single(
                    &topology,
                    &mut missing_coordinates,
                    &valence,
                    &rings,
                    &properties,
                    &mut diagnostics,
                    options,
                )
                .expect_err("cached route preserves selected coordinate-count failure");
                assert!(matches!(
                    error,
                    UffPreparedOptimizationError::Single(
                        super::SingleConformerOptimizationError::Construction(
                            super::AutomaticForceFieldConstructionError::Construction(
                                ForceFieldConstructionError::Builder(
                                    UffBuilderError::SelectedConformerCoordinateCountMismatch {
                                        conformer_id: 31,
                                        atoms: 2,
                                        coordinates: 1
                                    }
                                )
                            )
                        )
                    )
                ));
            } else {
                let error = super::optimize_single_conformer(
                    &topology,
                    &mut missing_coordinates,
                    &supplied_valences,
                    &supplied_conjugation,
                    &rings,
                    &valence,
                    &properties,
                    &mut diagnostics,
                    options,
                )
                .expect_err("supplied route preserves selected coordinate-count failure");
                assert!(matches!(
                    error,
                    super::SingleConformerOptimizationError::Construction(
                        super::AutomaticForceFieldConstructionError::Construction(
                            ForceFieldConstructionError::Builder(
                                UffBuilderError::SelectedConformerCoordinateCountMismatch {
                                    conformer_id: 31,
                                    atoms: 2,
                                    coordinates: 1
                                }
                            )
                        )
                    )
                ));
            }
            assert_eq!(missing_coordinates, original_missing_coordinates);
        }
    }

    #[test]
    fn uff_prepare_p07_cached_and_supplied_single_keep_typed_minimize_failure() {
        use super::super::optimization::UffPreparedOptimizationError;

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 6, Hybridization::Sp3, 0),
                worker_w06_atom(1, 6, Hybridization::Sp3, 0),
            ],
            &[(0, 1, BondOrder::Single)],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![1, 1],
            implicit_hydrogens: vec![0, 0],
        };
        let rings = fast_find_rings(&topology).expect("fixed p07 two-carbon topology has rings");
        let properties = MoleculeProperties::default();
        let supplied_valences = [4, 4];
        let supplied_conjugation = [false; 2];
        let mut actual_calls = 0;

        for cached in [true, false] {
            let mut coordinates = worker_w06_coordinates(&[[0.0, 0.0, 0.0], [1.514, 0.0, 0.0]]);
            let original = coordinates.clone();
            let mut options = super::SingleConformerOptions::for_conformer(31);
            options.ignore_interfragment_interactions = false;
            let mut diagnostics = Vec::new();

            if cached {
                let error = super::super::optimization::optimize_prepared_uff_single(
                    &topology,
                    &mut coordinates,
                    &valence,
                    &rings,
                    &properties,
                    &mut diagnostics,
                    options,
                )
                .expect_err("cached p07 source fixture reaches the typed minimizer failure");
                assert!(matches!(
                    error,
                    UffPreparedOptimizationError::Single(
                        super::SingleConformerOptimizationError::Optimization(
                            super::OptimizationStageError::Minimize(
                                super::ForceFieldKernelError::OptimizerBadDirection
                            )
                        )
                    )
                ));
            } else {
                let error = super::optimize_single_conformer(
                    &topology,
                    &mut coordinates,
                    &supplied_valences,
                    &supplied_conjugation,
                    &rings,
                    &valence,
                    &properties,
                    &mut diagnostics,
                    options,
                )
                .expect_err("supplied p07 source fixture reaches the typed minimizer failure");
                assert!(matches!(
                    error,
                    super::SingleConformerOptimizationError::Optimization(
                        super::OptimizationStageError::Minimize(
                            super::ForceFieldKernelError::OptimizerBadDirection
                        )
                    )
                ));
            }
            actual_calls += 1;
            assert_eq!(coordinates, original);
        }
        assert_eq!(actual_calls, 2);
    }

    #[test]
    fn uff_integrate_i04_prepared_serial_matches_fixed_w06_source_matrix() {
        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let mut actual_calls = 0;

        for row_count in 1..=IDS.len() {
            for (selected_row, selected_id) in IDS[..row_count].iter().copied().enumerate() {
                for max_iterations in [0, 1] {
                    let mut coordinates = CoordinateBlock::default();
                    for row in 0..row_count {
                        coordinates.conformers_3d.push(
                            Conformer3D::new(
                                IDS[row],
                                [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]].to_vec(),
                                true,
                            )
                            .with_prop("source-row", format!("row-{row}")),
                        );
                    }
                    coordinates.conformers_2d.push(
                        cosmolkit_model::Conformer2D::new(selected_id, TWO_D_SENTINEL.to_vec())
                            .with_prop("layout", "preserve-this-overlapping-2d-row"),
                    );
                    let before_3d = coordinates.conformers_3d.clone();
                    let before_2d = coordinates.conformers_2d[0].clone();
                    let before_dimension = coordinates.source_coordinate_dim;
                    let initial_result_len = if max_iterations == 0 && selected_row == 0 {
                        row_count - 1
                    } else if max_iterations == 1 && selected_row + 1 == row_count {
                        row_count + 1
                    } else {
                        row_count
                    };
                    let mut results = Vec::with_capacity(row_count + 1);
                    results.resize(initial_result_len, SENTINEL);
                    let result_storage = results.as_ptr();
                    let mut diagnostics = Vec::new();
                    let mut options = super::SingleConformerOptions::for_conformer(selected_id);
                    options.max_iterations = max_iterations;
                    options.ignore_interfragment_interactions = false;

                    super::super::optimization::optimize_prepared_uff_serial(
                        &topology,
                        &mut coordinates,
                        &mut results,
                        &valence,
                        &rings,
                        &properties,
                        &mut diagnostics,
                        options,
                    )
                    .unwrap_or_else(|error| {
                        panic!("fixed W06 prepared serial succeeds: rows={row_count}, construction-id={selected_id}, iterations={max_iterations}, error={error:?}")
                    });
                    actual_calls += 1;

                    assert_eq!(results.len(), row_count);
                    assert_eq!(results.as_ptr(), result_storage);
                    assert_eq!(coordinates.conformers_3d.len(), row_count);
                    assert_eq!(
                        coordinates
                            .conformers_3d
                            .iter()
                            .map(Conformer3D::id)
                            .collect::<Vec<_>>(),
                        &IDS[..row_count]
                    );
                    for row in 0..row_count {
                        assert_eq!(results[row].status, 1);
                        let (expected_xs, expected_energy) = if max_iterations == 0 {
                            (
                                [0.0, TARGET_DISTANCES[row]],
                                source_energy(TARGET_DISTANCES[row]),
                            )
                        } else {
                            source_one_step(TARGET_DISTANCES[row])
                        };
                        assert!(
                            (results[row].energy - expected_energy).abs() <= 1.0e-12,
                            "fixed W06 serial source energy row={row}, construction-id={selected_id}, iterations={max_iterations}"
                        );
                        assert_eq!(
                            coordinates.conformers_3d[row].is_3d(),
                            before_3d[row].is_3d()
                        );
                        assert_eq!(
                            coordinates.conformers_3d[row].props(),
                            before_3d[row].props()
                        );
                        let expected_points =
                            [[expected_xs[0], 0.0, 0.0], [expected_xs[1], 0.0, 0.0]];
                        for (point, expected_point) in expected_points.into_iter().enumerate() {
                            for (axis, expected_component) in expected_point.into_iter().enumerate()
                            {
                                assert_eq!(
                                    coordinates.conformers_3d[row].coordinates()[point][axis]
                                        .to_bits(),
                                    expected_component.to_bits(),
                                    "fixed W06 serial coordinate row={row}, point={point}, axis={axis}, construction-id={selected_id}, iterations={max_iterations}"
                                );
                            }
                        }
                    }
                    assert_eq!(coordinates.conformers_2d[0], before_2d);
                    assert_eq!(coordinates.conformers_2d[0].id(), selected_id);
                    assert_eq!(coordinates.source_coordinate_dim, before_dimension);
                    assert!(diagnostics.is_empty());
                }
            }
        }

        assert_eq!(actual_calls, 12);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_integrate_i05_prepared_dispatch_matches_fixed_w06_source_matrix() {
        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const REQUESTED: [i32; 5] = [-1, 0, 1, 2, 3];
        const HARDWARE: [u32; 2] = [0, 3];
        const SOURCE_THREAD_COUNTS: [[u32; 5]; 2] = [[1, 1, 1, 2, 3], [2, 3, 1, 2, 3]];
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let mut actual_calls = 0;

        for row_count in 1..=IDS.len() {
            for (selected_row, selected_id) in IDS[..row_count].iter().copied().enumerate() {
                for (requested_index, requested_threads) in REQUESTED.into_iter().enumerate() {
                    for (hardware_index, observed_hardware) in HARDWARE.into_iter().enumerate() {
                        for threadsafe in [false, true] {
                            let expected_source_count = if threadsafe {
                                SOURCE_THREAD_COUNTS[hardware_index][requested_index]
                            } else {
                                1
                            };
                            for max_iterations in [0_i32, 1] {
                                let mut coordinates = CoordinateBlock::default();
                                for row in 0..row_count {
                                    coordinates.conformers_3d.push(
                                        Conformer3D::new(
                                            IDS[row],
                                            [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]]
                                                .to_vec(),
                                            true,
                                        )
                                        .with_prop("source-row", format!("row-{row}")),
                                    );
                                }
                                coordinates.conformers_2d.push(
                                    cosmolkit_model::Conformer2D::new(
                                        selected_id,
                                        TWO_D_SENTINEL.to_vec(),
                                    )
                                    .with_prop("layout", "preserve-this-overlapping-2d-row"),
                                );
                                let before_3d = coordinates.conformers_3d.clone();
                                let before_2d = coordinates.conformers_2d[0].clone();
                                let before_dimension = coordinates.source_coordinate_dim;
                                let mut results = Vec::new();
                                let mut diagnostics = Vec::new();
                                let mut options =
                                    super::SingleConformerOptions::for_conformer(selected_id);
                                options.max_iterations = max_iterations;
                                options.ignore_interfragment_interactions = false;

                                let outcome = super::super::optimization::
                                    optimize_prepared_uff_dispatch(
                                        &topology,
                                        &mut coordinates,
                                        &mut results,
                                        &valence,
                                        &rings,
                                        &properties,
                                        &mut diagnostics,
                                        options,
                                        requested_threads,
                                        observed_hardware,
                                        threadsafe,
                                    )
                                    .unwrap_or_else(|error| {
                                        panic!("fixed W06 prepared dispatch succeeds: rows={row_count}, construction-id={selected_id}, requested={requested_threads}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}, error={error:?}")
                                    });
                                actual_calls += 1;

                                match outcome {
                                    super::PreparedConformerDispatchOutcome::Serial(result) => {
                                        assert_eq!(
                                            expected_source_count, 1,
                                            "W06 serial route: rows={row_count}, construction-id={selected_id}, requested={requested_threads}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                        );
                                        result.unwrap_or_else(|error| {
                                            panic!("fixed W06 serial result succeeds: {error:?}")
                                        });
                                    }
                                    super::PreparedConformerDispatchOutcome::Workers(result) => {
                                        assert_ne!(
                                            expected_source_count, 1,
                                            "W06 worker route: rows={row_count}, construction-id={selected_id}, requested={requested_threads}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                        );
                                        let joins = result.unwrap_or_else(|error| {
                                            panic!("fixed W06 worker setup succeeds: {error:?}")
                                        });
                                        assert_eq!(joins.len(), expected_source_count as usize);
                                        for (lane, join_result) in joins.into_iter().enumerate() {
                                            let row_result = join_result.unwrap_or_else(|panic| {
                                                panic!(
                                                    "fixed W06 worker lane {lane} joined: {panic:?}"
                                                )
                                            });
                                            row_result.unwrap_or_else(|error| {
                                                panic!("fixed W06 worker lane {lane} result succeeds: {error:?}")
                                            });
                                        }
                                    }
                                }

                                assert_eq!(results.len(), row_count);
                                assert_eq!(coordinates.conformers_3d.len(), row_count);
                                assert_eq!(
                                    coordinates
                                        .conformers_3d
                                        .iter()
                                        .map(Conformer3D::id)
                                        .collect::<Vec<_>>(),
                                    &IDS[..row_count]
                                );
                                for row in 0..row_count {
                                    assert_eq!(results[row].status, 1);
                                    let (expected_xs, expected_energy) = if max_iterations == 0 {
                                        (
                                            [0.0, TARGET_DISTANCES[row]],
                                            source_energy(TARGET_DISTANCES[row]),
                                        )
                                    } else {
                                        source_one_step(TARGET_DISTANCES[row])
                                    };
                                    assert!(
                                        (results[row].energy - expected_energy).abs() <= 1.0e-12,
                                        "fixed W06 source energy: row={row}, construction-id={selected_id}, requested={requested_threads}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                    );
                                    assert_eq!(
                                        coordinates.conformers_3d[row].is_3d(),
                                        before_3d[row].is_3d()
                                    );
                                    assert_eq!(
                                        coordinates.conformers_3d[row].props(),
                                        before_3d[row].props()
                                    );
                                    let expected_points =
                                        [[expected_xs[0], 0.0, 0.0], [expected_xs[1], 0.0, 0.0]];
                                    for (point, expected_point) in
                                        expected_points.into_iter().enumerate()
                                    {
                                        for (axis, expected_component) in
                                            expected_point.into_iter().enumerate()
                                        {
                                            assert_eq!(
                                                coordinates.conformers_3d[row].coordinates()[point]
                                                    [axis]
                                                    .to_bits(),
                                                expected_component.to_bits(),
                                                "fixed W06 source coordinate: row={row}, point={point}, axis={axis}, construction-id={selected_id}, requested={requested_threads}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                            );
                                        }
                                    }
                                }
                                assert_eq!(coordinates.conformers_2d[0], before_2d);
                                assert_eq!(coordinates.conformers_2d[0].id(), selected_id);
                                assert_eq!(coordinates.source_coordinate_dim, before_dimension);
                                assert!(diagnostics.is_empty());
                            }
                        }
                    }
                }
            }
        }

        assert_eq!(actual_calls, 240);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_integrate_i06_prepared_entries_preserve_borrowed_storage() {
        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let properties = MoleculeProperties::default();
        let topology_before = topology.clone();
        let valence_before = valence.clone();
        let rings_before = rings.clone();
        let properties_before = properties.clone();
        let topology_address = std::ptr::from_ref(&topology);
        let valence_address = std::ptr::from_ref(&valence);
        let rings_address = std::ptr::from_ref(&rings);
        let properties_address = std::ptr::from_ref(&properties);

        let assert_borrowed_inputs_unchanged = || {
            assert!(std::ptr::eq(
                topology_address,
                std::ptr::from_ref(&topology)
            ));
            assert!(std::ptr::eq(valence_address, std::ptr::from_ref(&valence)));
            assert!(std::ptr::eq(rings_address, std::ptr::from_ref(&rings)));
            assert!(std::ptr::eq(
                properties_address,
                std::ptr::from_ref(&properties)
            ));
            assert_eq!(topology, topology_before);
            assert_eq!(valence, valence_before);
            assert_eq!(rings, rings_before);
            assert_eq!(properties, properties_before);
        };

        let snapshot_3d = |coordinates: &CoordinateBlock| {
            coordinates
                .conformers_3d
                .iter()
                .map(|conformer| {
                    (
                        conformer.id(),
                        conformer.is_3d(),
                        conformer.props().clone(),
                        conformer
                            .coordinates()
                            .iter()
                            .map(|point| {
                                [point[0].to_bits(), point[1].to_bits(), point[2].to_bits()]
                            })
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let snapshot_2d = |coordinates: &CoordinateBlock| {
            coordinates
                .conformers_2d
                .iter()
                .map(|conformer| {
                    (
                        conformer.id(),
                        conformer.props().clone(),
                        conformer
                            .coordinates()
                            .iter()
                            .map(|point| [point[0].to_bits(), point[1].to_bits()])
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let coordinate_storage = |coordinates: &CoordinateBlock| {
            (
                coordinates.conformers_3d.as_ptr(),
                coordinates.conformers_2d.as_ptr(),
                coordinates
                    .conformers_3d
                    .iter()
                    .map(|conformer| conformer.coordinates().as_ptr())
                    .collect::<Vec<_>>(),
                coordinates
                    .conformers_2d
                    .iter()
                    .map(|conformer| conformer.coordinates().as_ptr())
                    .collect::<Vec<_>>(),
            )
        };

        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };

        let mut actual_calls = 0;
        for mode in 0..3 {
            for row_count in 1..=IDS.len() {
                for (selected_row, selected_id) in IDS[..row_count].iter().copied().enumerate() {
                    for max_iterations in [0, 1] {
                        assert_borrowed_inputs_unchanged();

                        let mut coordinates = CoordinateBlock::default();
                        for row in 0..row_count {
                            coordinates.conformers_3d.push(
                                Conformer3D::new(
                                    IDS[row],
                                    [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]].to_vec(),
                                    true,
                                )
                                .with_prop("source-row", format!("row-{row}")),
                            );
                        }
                        coordinates.conformers_2d.push(
                            cosmolkit_model::Conformer2D::new(selected_id, TWO_D_SENTINEL.to_vec())
                                .with_prop("layout", "preserve-this-overlapping-2d-row"),
                        );

                        let before_3d = snapshot_3d(&coordinates);
                        let before_2d = snapshot_2d(&coordinates);
                        let before_dimension = coordinates.source_coordinate_dim;
                        let before_coordinate_storage = coordinate_storage(&coordinates);
                        let mut diagnostics = Vec::new();
                        let mut options = super::SingleConformerOptions::for_conformer(selected_id);
                        options.max_iterations = max_iterations;
                        options.ignore_interfragment_interactions = false;

                        let expected_rows = (0..row_count)
                            .map(|row| {
                                let changes =
                                    max_iterations == 1 && (mode != 0 || row == selected_row);
                                let (xs, energy) = if changes {
                                    source_one_step(TARGET_DISTANCES[row])
                                } else {
                                    (
                                        [0.0, TARGET_DISTANCES[row]],
                                        source_energy(TARGET_DISTANCES[row]),
                                    )
                                };
                                ([[xs[0], 0.0, 0.0], [xs[1], 0.0, 0.0]], energy)
                            })
                            .collect::<Vec<_>>();

                        let outcomes = match mode {
                            0 => {
                                let outcome = super::super::optimization::
                                    optimize_prepared_uff_single(
                                        &topology,
                                        &mut coordinates,
                                        &valence,
                                        &rings,
                                        &properties,
                                        &mut diagnostics,
                                        options,
                                    )
                                    .unwrap_or_else(|error| {
                                        panic!("fixed W06 I06 single succeeds: rows={row_count}, selected={selected_id}, iterations={max_iterations}, error={error:?}")
                                    });
                                actual_calls += 1;
                                vec![outcome]
                            }
                            1 => {
                                let mut results = vec![SENTINEL; row_count];
                                let result_storage = results.as_ptr();
                                super::super::optimization::optimize_prepared_uff_serial(
                                    &topology,
                                    &mut coordinates,
                                    &mut results,
                                    &valence,
                                    &rings,
                                    &properties,
                                    &mut diagnostics,
                                    options,
                                )
                                .unwrap_or_else(|error| {
                                    panic!("fixed W06 I06 serial succeeds: rows={row_count}, selected={selected_id}, iterations={max_iterations}, error={error:?}")
                                });
                                actual_calls += 1;
                                assert_eq!(results.as_ptr(), result_storage);
                                results
                            }
                            _ => {
                                let mut results = vec![SENTINEL; row_count];
                                let result_storage = results.as_ptr();
                                let outcome = super::super::optimization::
                                    optimize_prepared_uff_dispatch(
                                        &topology,
                                        &mut coordinates,
                                        &mut results,
                                        &valence,
                                        &rings,
                                        &properties,
                                        &mut diagnostics,
                                        options,
                                        3,
                                        3,
                                        true,
                                    )
                                    .unwrap_or_else(|error| {
                                        panic!("fixed W06 I06 dispatch succeeds: rows={row_count}, selected={selected_id}, iterations={max_iterations}, error={error:?}")
                                    });
                                actual_calls += 1;
                                assert_eq!(results.as_ptr(), result_storage);
                                match outcome {
                                    super::PreparedConformerDispatchOutcome::Serial(result) => {
                                        panic!("fixed W06 I06 dispatch selected Serial: {result:?}")
                                    }
                                    super::PreparedConformerDispatchOutcome::Workers(result) => {
                                        let joins = result.unwrap_or_else(|error| {
                                            panic!("fixed W06 I06 worker setup succeeds: {error:?}")
                                        });
                                        assert_eq!(joins.len(), 3);
                                        for (lane, joined) in joins.into_iter().enumerate() {
                                            joined.unwrap_or_else(|panic| {
                                                panic!("fixed W06 I06 worker lane {lane} joins: {panic:?}")
                                            }).unwrap_or_else(|error| {
                                                panic!("fixed W06 I06 worker lane {lane} succeeds: {error:?}")
                                            });
                                        }
                                    }
                                }
                                results
                            }
                        };

                        assert_borrowed_inputs_unchanged();
                        assert_eq!(diagnostics, Vec::new());
                        assert_eq!(coordinates.source_coordinate_dim, before_dimension);
                        assert_eq!(coordinates.conformers_3d.len(), row_count);
                        assert_eq!(snapshot_2d(&coordinates), before_2d);
                        assert_eq!(
                            coordinate_storage(&coordinates),
                            before_coordinate_storage,
                            "W06 coordinate buffers remain borrowed in place: mode={mode}, rows={row_count}, selected={selected_id}, iterations={max_iterations}"
                        );

                        let after_3d = snapshot_3d(&coordinates);
                        assert_eq!(after_3d.len(), row_count);
                        for row in 0..row_count {
                            assert_eq!(after_3d[row].0, before_3d[row].0);
                            assert_eq!(after_3d[row].1, before_3d[row].1);
                            assert_eq!(after_3d[row].2, before_3d[row].2);
                            assert_eq!(after_3d[row].0, IDS[row]);
                            for point in 0..2 {
                                for axis in 0..3 {
                                    assert_eq!(
                                        after_3d[row].3[point][axis],
                                        expected_rows[row].0[point][axis].to_bits(),
                                        "fixed W06 source coordinate: mode={mode}, row={row}, point={point}, axis={axis}, selected={selected_id}, iterations={max_iterations}"
                                    );
                                }
                            }
                        }

                        if mode == 0 {
                            assert_eq!(outcomes.len(), 1);
                            assert_eq!(outcomes[0].status, 1);
                            assert!(
                                (outcomes[0].energy - expected_rows[selected_row].1).abs()
                                    <= 1.0e-12,
                                "fixed W06 source single energy: selected={selected_id}, iterations={max_iterations}"
                            );
                        } else {
                            assert_eq!(outcomes.len(), row_count);
                            for row in 0..row_count {
                                assert_eq!(outcomes[row].status, 1);
                                assert!(
                                    (outcomes[row].energy - expected_rows[row].1).abs() <= 1.0e-12,
                                    "fixed W06 source all-row energy: mode={mode}, row={row}, selected={selected_id}, iterations={max_iterations}"
                                );
                            }
                        }
                        assert!(diagnostics.is_empty());
                    }
                }
            }
        }

        assert_eq!(actual_calls, 36);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_integrate_i07_prepared_entries_preserve_errors_and_order() {
        use super::super::api::UffParameterError;
        use super::super::builder::{
            ForceFieldConstructionError, PreparedValenceField, UffBuilderError,
        };
        use super::super::optimization::UffPreparedOptimizationError;

        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        const RESULT_SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        fn require_error<T>(
            result: Result<T, UffPreparedOptimizationError>,
        ) -> UffPreparedOptimizationError {
            match result {
                Err(error) => error,
                Ok(_) => panic!("the fixed I07 source case must return its typed error"),
            }
        }

        fn assert_source_is_stored<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let source = std::error::Error::source(parent)
                .expect("the source-bearing error exposes its stored child");
            let stored_error: &(dyn std::error::Error + 'static) = stored;
            assert!(std::ptr::eq(source, stored_error));
            assert!(std::ptr::eq(
                source
                    .downcast_ref::<Child>()
                    .expect("the child keeps its type"),
                stored
            ));
        }

        fn assert_preparation_cause(
            error: &UffPreparedOptimizationError,
            expected: &UffBuilderError,
        ) {
            let parameter: &UffParameterError = match error {
                UffPreparedOptimizationError::Preparation(stored) => stored,
                _ => panic!("invalid assignment must fail at shared preparation: {error:?}"),
            };
            assert_source_is_stored(error, parameter);
            assert!(matches!(
                parameter.kind(),
                super::super::api::UffParameterErrorKind::Preparation
            ));
            let builder = std::error::Error::source(parameter)
                .expect("preparation retains the concrete builder cause")
                .downcast_ref::<UffBuilderError>()
                .expect("the preparation source remains UffBuilderError");
            assert_source_is_stored(parameter, builder);
            assert_eq!(builder, expected);
            assert!(std::error::Error::source(builder).is_none());
        }

        fn assert_selected_id_cause(error: &UffPreparedOptimizationError, selected_id: usize) {
            let automatic = match error {
                UffPreparedOptimizationError::Single(owner) => {
                    assert_source_is_stored(error, owner);
                    match owner {
                        super::SingleConformerOptimizationError::Construction(stored) => {
                            assert_source_is_stored(owner, stored);
                            stored
                        }
                        _ => panic!("single selected-ID failure is not construction"),
                    }
                }
                UffPreparedOptimizationError::Serial(owner) => {
                    assert_source_is_stored(error, owner);
                    match owner {
                        super::SerialUffOptimizationError::Construction(stored) => {
                            assert_source_is_stored(owner, stored);
                            stored
                        }
                        _ => panic!("serial selected-ID failure is not construction"),
                    }
                }
                UffPreparedOptimizationError::Dispatch(owner) => {
                    assert_source_is_stored(error, owner);
                    match owner {
                        super::DispatchedUffOptimizationError::Construction(stored) => {
                            assert_source_is_stored(owner, stored);
                            stored
                        }
                        _ => panic!("dispatch selected-ID failure is not construction"),
                    }
                }
                _ => panic!("selected-ID failure must be a delegated construction error"),
            };
            let construction = match automatic {
                super::AutomaticForceFieldConstructionError::Construction(stored) => stored,
                _ => panic!("selected-ID lookup follows table and typing stages"),
            };
            assert_source_is_stored(automatic, construction);
            let builder = match construction {
                ForceFieldConstructionError::Builder(stored) => stored,
                _ => panic!("selected-ID lookup is the builder construction cause"),
            };
            assert_source_is_stored(construction, builder);
            assert_eq!(
                builder,
                &UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                    conformer_id: selected_id,
                }
            );
            assert!(std::error::Error::source(builder).is_none());
        }

        fn assert_thread_count_cause(error: &UffPreparedOptimizationError) {
            let dispatch = match error {
                UffPreparedOptimizationError::Dispatch(stored) => stored,
                _ => panic!("INT_MIN thread resolution is an outer dispatch failure"),
            };
            assert_source_is_stored(error, dispatch);
            let thread_count = match dispatch {
                super::DispatchedUffOptimizationError::ThreadCount(stored) => stored,
                _ => panic!("INT_MIN reaches the thread-count stage after construction"),
            };
            assert_source_is_stored(dispatch, thread_count);
            assert_eq!(
                *thread_count,
                super::UffThreadCountError::UndefinedSignedNegation
            );
            assert!(std::error::Error::source(thread_count).is_none());
        }

        let snapshot_3d = |coordinates: &CoordinateBlock| {
            coordinates
                .conformers_3d
                .iter()
                .map(|conformer| {
                    (
                        conformer.id(),
                        conformer.is_3d(),
                        conformer.props().clone(),
                        conformer
                            .coordinates()
                            .iter()
                            .map(|point| {
                                [point[0].to_bits(), point[1].to_bits(), point[2].to_bits()]
                            })
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let snapshot_2d = |coordinates: &CoordinateBlock| {
            coordinates
                .conformers_2d
                .iter()
                .map(|conformer| {
                    (
                        conformer.id(),
                        conformer.props().clone(),
                        conformer
                            .coordinates()
                            .iter()
                            .map(|point| [point[0].to_bits(), point[1].to_bits()])
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let coordinate_storage = |coordinates: &CoordinateBlock| {
            (
                coordinates.conformers_3d.as_ptr(),
                coordinates.conformers_2d.as_ptr(),
                coordinates
                    .conformers_3d
                    .iter()
                    .map(|conformer| conformer.coordinates().as_ptr())
                    .collect::<Vec<_>>(),
                coordinates
                    .conformers_2d
                    .iter()
                    .map(|conformer| conformer.coordinates().as_ptr())
                    .collect::<Vec<_>>(),
            )
        };
        let coordinate_state = |coordinates: &CoordinateBlock| {
            (
                snapshot_3d(coordinates),
                snapshot_2d(coordinates),
                coordinates.source_coordinate_dim,
                coordinate_storage(coordinates),
            )
        };
        let make_coordinates = |stored_3d_id: Option<usize>, layout_id: usize| {
            let mut coordinates = CoordinateBlock::default();
            if let Some(id) = stored_3d_id {
                coordinates.conformers_3d.push(
                    Conformer3D::new(id, [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]].to_vec(), true)
                        .with_prop("source-row", "w06-row"),
                );
            }
            coordinates.conformers_2d.push(
                cosmolkit_model::Conformer2D::new(layout_id, TWO_D_SENTINEL.to_vec())
                    .with_prop("layout", "i07-preserved"),
            );
            coordinates
        };

        let atom = |row, atomic_number, charge| {
            Atom::from_spec(
                AtomId::new(row),
                AtomSpec::new(
                    Element::from_atomic_number(atomic_number)
                        .expect("fixed W06 element exists in the element table"),
                )
                .with_hybridization(Hybridization::Unspecified)
                .with_formal_charge(charge)
                .with_no_implicit(false),
            )
        };
        let preparation_topology = worker_w06_topology(vec![atom(0, 11, 1), atom(1, 17, -1)], &[]);
        assert!(
            preparation_topology
                .atoms
                .iter()
                .all(|atom| !atom.no_implicit())
        );
        let preparation_cases = vec![
            (
                "explicit-length-short",
                ValenceAssignment {
                    explicit_valence: vec![0],
                    implicit_hydrogens: vec![0, 0],
                },
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::Explicit,
                    expected: 2,
                    actual: 1,
                },
            ),
            (
                "explicit-length-long",
                ValenceAssignment {
                    explicit_valence: vec![0, 0, 0],
                    implicit_hydrogens: vec![0, 0],
                },
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::Explicit,
                    expected: 2,
                    actual: 3,
                },
            ),
            (
                "implicit-length-short",
                ValenceAssignment {
                    explicit_valence: vec![0, 0],
                    implicit_hydrogens: vec![0],
                },
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::ImplicitHydrogen,
                    expected: 2,
                    actual: 1,
                },
            ),
            (
                "implicit-length-long",
                ValenceAssignment {
                    explicit_valence: vec![0, 0],
                    implicit_hydrogens: vec![0, 0, 0],
                },
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::ImplicitHydrogen,
                    expected: 2,
                    actual: 3,
                },
            ),
            (
                "explicit-component-minus-one",
                ValenceAssignment {
                    explicit_valence: vec![-1, 0],
                    implicit_hydrogens: vec![0, 0],
                },
                UffBuilderError::SourceValencePrecondition {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::Explicit,
                    value: -1,
                },
            ),
            (
                "explicit-component-128",
                ValenceAssignment {
                    explicit_valence: vec![128, 0],
                    implicit_hydrogens: vec![0, 0],
                },
                UffBuilderError::SourceValenceOutOfRange {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::Explicit,
                    value: 128,
                },
            ),
            (
                "implicit-component-minus-one",
                ValenceAssignment {
                    explicit_valence: vec![0, 0],
                    implicit_hydrogens: vec![-1, 0],
                },
                UffBuilderError::SourceValencePrecondition {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::ImplicitHydrogen,
                    value: -1,
                },
            ),
            (
                "implicit-component-128",
                ValenceAssignment {
                    explicit_valence: vec![0, 0],
                    implicit_hydrogens: vec![128, 0],
                },
                UffBuilderError::SourceValenceOutOfRange {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::ImplicitHydrogen,
                    value: 128,
                },
            ),
        ];
        let mut preparation_calls = 0;
        for (case_name, assignment, expected_builder_error) in &preparation_cases {
            for mode in 0..3 {
                let mut coordinates = make_coordinates(Some(80), 80);
                let before_coordinates = coordinate_state(&coordinates);
                let mut diagnostics = Vec::new();
                let before_diagnostics = diagnostics.clone();
                let mut results = vec![RESULT_SENTINEL; 2];
                let before_results = results.clone();
                let results_address = results.as_ptr();
                let mut options = super::SingleConformerOptions::for_conformer(80);
                options.max_iterations = 0;
                options.ignore_interfragment_interactions = false;

                let error = match mode {
                    0 => require_error(super::super::optimization::optimize_prepared_uff_single(
                        &preparation_topology,
                        &mut coordinates,
                        assignment,
                        &fast_find_rings(&preparation_topology)
                            .expect("fixed W06 topology has ring information"),
                        &MoleculeProperties::default(),
                        &mut diagnostics,
                        options,
                    )),
                    1 => {
                        let rings = fast_find_rings(&preparation_topology)
                            .expect("fixed W06 topology has ring information");
                        super::super::optimization::optimize_prepared_uff_serial(
                            &preparation_topology,
                            &mut coordinates,
                            &mut results,
                            assignment,
                            &rings,
                            &MoleculeProperties::default(),
                            &mut diagnostics,
                            options,
                        )
                        .err()
                        .unwrap_or_else(|| panic!("I07 {case_name} serial must return Preparation"))
                    }
                    _ => {
                        let rings = fast_find_rings(&preparation_topology)
                            .expect("fixed W06 topology has ring information");
                        require_error(super::super::optimization::optimize_prepared_uff_dispatch(
                            &preparation_topology,
                            &mut coordinates,
                            &mut results,
                            assignment,
                            &rings,
                            &MoleculeProperties::default(),
                            &mut diagnostics,
                            options,
                            1,
                            3,
                            false,
                        ))
                    }
                };
                preparation_calls += 1;
                assert_preparation_cause(&error, expected_builder_error);
                assert_eq!(coordinate_state(&coordinates), before_coordinates);
                assert_eq!(diagnostics, before_diagnostics);
                assert_eq!(results.as_ptr(), results_address);
                assert_eq!(results, before_results);
            }
        }
        assert_eq!(preparation_calls, 24);

        let construction_topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valid_assignment = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let construction_rings = fast_find_rings(&construction_topology)
            .expect("fixed W06 topology has ring information");
        let construction_properties = MoleculeProperties::default();
        let mut construction_calls = 0;
        for (has_3d_row, mode) in [
            (true, 0),
            (true, 1),
            (true, 2),
            (false, 0),
            (false, 1),
            (false, 2),
        ] {
            let mut coordinates = make_coordinates(has_3d_row.then_some(80), 999);
            let before_coordinates = coordinate_state(&coordinates);
            let mut diagnostics = Vec::new();
            let before_diagnostics = diagnostics.clone();
            let mut results = vec![RESULT_SENTINEL; 2];
            let before_results = results.clone();
            let results_address = results.as_ptr();
            let mut options = super::SingleConformerOptions::for_conformer(999);
            options.max_iterations = 0;
            options.ignore_interfragment_interactions = false;

            let error = match mode {
                0 => require_error(super::super::optimization::optimize_prepared_uff_single(
                    &construction_topology,
                    &mut coordinates,
                    &valid_assignment,
                    &construction_rings,
                    &construction_properties,
                    &mut diagnostics,
                    options,
                )),
                1 => super::super::optimization::optimize_prepared_uff_serial(
                    &construction_topology,
                    &mut coordinates,
                    &mut results,
                    &valid_assignment,
                    &construction_rings,
                    &construction_properties,
                    &mut diagnostics,
                    options,
                )
                .err()
                .unwrap_or_else(|| panic!("I07 selected-ID case serial must fail construction")),
                _ => require_error(super::super::optimization::optimize_prepared_uff_dispatch(
                    &construction_topology,
                    &mut coordinates,
                    &mut results,
                    &valid_assignment,
                    &construction_rings,
                    &construction_properties,
                    &mut diagnostics,
                    options,
                    1,
                    3,
                    false,
                )),
            };
            construction_calls += 1;
            assert_selected_id_cause(&error, 999);
            assert_eq!(coordinate_state(&coordinates), before_coordinates);
            assert_eq!(diagnostics, before_diagnostics);
            assert_eq!(results.as_ptr(), results_address);
            assert_eq!(results, before_results);
        }
        assert_eq!(construction_calls, 6);

        let thread_topology = construction_topology;
        let thread_rings = construction_rings;
        let thread_properties = construction_properties;
        let thread_valence = valid_assignment;
        let source_energy = |distance: f64| {
            let minimum = (2.983_f64 * 3.947).sqrt();
            let well_depth = (0.03_f64 * 0.227).sqrt();
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let mut int_min_calls = 0;

        let mut true_coordinates = make_coordinates(Some(80), 80);
        let true_before = coordinate_state(&true_coordinates);
        let mut true_diagnostics = Vec::new();
        let mut true_results = Vec::new();
        let mut true_options = super::SingleConformerOptions::for_conformer(80);
        true_options.max_iterations = 0;
        true_options.ignore_interfragment_interactions = false;
        let true_error = require_error(super::super::optimization::optimize_prepared_uff_dispatch(
            &thread_topology,
            &mut true_coordinates,
            &mut true_results,
            &thread_valence,
            &thread_rings,
            &thread_properties,
            &mut true_diagnostics,
            true_options,
            i32::MIN,
            3,
            true,
        ));
        int_min_calls += 1;
        assert_thread_count_cause(&true_error);
        assert_eq!(true_results.len(), 1);
        assert_eq!(true_results[0].status, 0);
        assert_eq!(true_results[0].energy.to_bits(), 0.0_f64.to_bits());
        assert_eq!(coordinate_state(&true_coordinates), true_before);
        assert!(true_diagnostics.is_empty());

        let mut false_coordinates = make_coordinates(Some(80), 80);
        let false_before = coordinate_state(&false_coordinates);
        let mut false_diagnostics = Vec::new();
        let mut false_results = Vec::new();
        let mut false_options = super::SingleConformerOptions::for_conformer(80);
        false_options.max_iterations = 0;
        false_options.ignore_interfragment_interactions = false;
        let false_outcome = super::super::optimization::optimize_prepared_uff_dispatch(
            &thread_topology,
            &mut false_coordinates,
            &mut false_results,
            &thread_valence,
            &thread_rings,
            &thread_properties,
            &mut false_diagnostics,
            false_options,
            i32::MIN,
            3,
            false,
        )
        .unwrap_or_else(|error| panic!("source-mode false uses real ST execution: {error:?}"));
        int_min_calls += 1;
        match false_outcome {
            super::PreparedConformerDispatchOutcome::Serial(result) => result
                .unwrap_or_else(|error| panic!("source-mode false ST result succeeds: {error:?}")),
            super::PreparedConformerDispatchOutcome::Workers(_) => {
                panic!("source-mode false bypasses workers and uses the ST owner")
            }
        }
        assert_eq!(false_results.len(), 1);
        assert_eq!(false_results[0].status, 1);
        assert!(
            (false_results[0].energy - source_energy(4.0)).abs() <= 1.0e-12,
            "INT_MIN source-mode false has the independent zero-step W06 energy"
        );
        assert_eq!(coordinate_state(&false_coordinates), false_before);
        assert!(false_diagnostics.is_empty());

        assert_eq!(int_min_calls, 2);
        assert_eq!(preparation_calls + construction_calls + int_min_calls, 32);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_integrate_i08_prepared_entries_prepare_once_per_call() {
        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let properties = MoleculeProperties::default();
        let topology_before = topology.clone();
        let valence_before = valence.clone();
        let rings_before = rings.clone();
        let properties_before = properties.clone();
        let topology_address = std::ptr::from_ref(&topology);
        let valence_address = std::ptr::from_ref(&valence);
        let rings_address = std::ptr::from_ref(&rings);
        let properties_address = std::ptr::from_ref(&properties);
        let assert_borrowed_inputs_unchanged = || {
            assert!(std::ptr::eq(
                topology_address,
                std::ptr::from_ref(&topology)
            ));
            assert!(std::ptr::eq(valence_address, std::ptr::from_ref(&valence)));
            assert!(std::ptr::eq(rings_address, std::ptr::from_ref(&rings)));
            assert!(std::ptr::eq(
                properties_address,
                std::ptr::from_ref(&properties)
            ));
            assert_eq!(topology, topology_before);
            assert_eq!(valence, valence_before);
            assert_eq!(rings, rings_before);
            assert_eq!(properties, properties_before);
        };

        let snapshot_3d = |coordinates: &CoordinateBlock| {
            coordinates
                .conformers_3d
                .iter()
                .map(|conformer| {
                    (
                        conformer.id(),
                        conformer.is_3d(),
                        conformer.props().clone(),
                        conformer
                            .coordinates()
                            .iter()
                            .map(|point| {
                                [point[0].to_bits(), point[1].to_bits(), point[2].to_bits()]
                            })
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let snapshot_2d = |coordinates: &CoordinateBlock| {
            coordinates
                .conformers_2d
                .iter()
                .map(|conformer| {
                    (
                        conformer.id(),
                        conformer.props().clone(),
                        conformer
                            .coordinates()
                            .iter()
                            .map(|point| [point[0].to_bits(), point[1].to_bits()])
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let coordinate_storage = |coordinates: &CoordinateBlock| {
            (
                coordinates.conformers_3d.as_ptr(),
                coordinates.conformers_2d.as_ptr(),
                coordinates
                    .conformers_3d
                    .iter()
                    .map(|conformer| conformer.coordinates().as_ptr())
                    .collect::<Vec<_>>(),
                coordinates
                    .conformers_2d
                    .iter()
                    .map(|conformer| conformer.coordinates().as_ptr())
                    .collect::<Vec<_>>(),
            )
        };

        // W06's fixed reference uses the pinned UFF geometric means and
        // vdW energy/gradient arithmetic; it never reads an optimizer result.
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };

        let mut actual_calls = 0;
        for mode in 0..3 {
            for row_count in 1..=IDS.len() {
                for max_iterations in [0, 1] {
                    assert_borrowed_inputs_unchanged();
                    let selected_id = IDS[0];
                    let mut coordinates = CoordinateBlock::default();
                    for row in 0..row_count {
                        coordinates.conformers_3d.push(
                            Conformer3D::new(
                                IDS[row],
                                [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]].to_vec(),
                                true,
                            )
                            .with_prop("source-row", format!("row-{row}")),
                        );
                    }
                    coordinates.conformers_2d.push(
                        cosmolkit_model::Conformer2D::new(selected_id, TWO_D_SENTINEL.to_vec())
                            .with_prop("layout", "preserve-this-overlapping-2d-row"),
                    );

                    let before_3d = snapshot_3d(&coordinates);
                    let before_2d = snapshot_2d(&coordinates);
                    let before_dimension = coordinates.source_coordinate_dim;
                    let before_coordinate_storage = coordinate_storage(&coordinates);
                    let mut diagnostics = Vec::new();
                    let mut options = super::SingleConformerOptions::for_conformer(selected_id);
                    options.max_iterations = max_iterations;
                    options.ignore_interfragment_interactions = false;

                    let expected_rows = (0..row_count)
                        .map(|row| {
                            let changes = max_iterations == 1 && (mode != 0 || row == 0);
                            let (xs, energy) = if changes {
                                source_one_step(TARGET_DISTANCES[row])
                            } else {
                                (
                                    [0.0, TARGET_DISTANCES[row]],
                                    source_energy(TARGET_DISTANCES[row]),
                                )
                            };
                            ([[xs[0], 0.0, 0.0], [xs[1], 0.0, 0.0]], energy)
                        })
                        .collect::<Vec<_>>();

                    super::super::api::reset_prepare_parameter_query_calls();
                    assert_eq!(super::super::api::prepare_parameter_query_calls(), 0);
                    let outcomes = match mode {
                        0 => {
                            let outcome = super::super::optimization::optimize_prepared_uff_single(
                                &topology,
                                &mut coordinates,
                                &valence,
                                &rings,
                                &properties,
                                &mut diagnostics,
                                options,
                            )
                            .unwrap_or_else(|error| {
                                panic!("fixed W06 I08 single succeeds: rows={row_count}, iterations={max_iterations}, error={error:?}")
                            });
                            actual_calls += 1;
                            vec![outcome]
                        }
                        1 => {
                            let mut results = vec![SENTINEL; row_count];
                            super::super::optimization::optimize_prepared_uff_serial(
                                &topology,
                                &mut coordinates,
                                &mut results,
                                &valence,
                                &rings,
                                &properties,
                                &mut diagnostics,
                                options,
                            )
                            .unwrap_or_else(|error| {
                                panic!("fixed W06 I08 serial succeeds: rows={row_count}, iterations={max_iterations}, error={error:?}")
                            });
                            actual_calls += 1;
                            results
                        }
                        _ => {
                            let mut results = vec![SENTINEL; row_count];
                            let dispatch = super::super::optimization::optimize_prepared_uff_dispatch(
                                &topology,
                                &mut coordinates,
                                &mut results,
                                &valence,
                                &rings,
                                &properties,
                                &mut diagnostics,
                                options,
                                3,
                                3,
                                true,
                            )
                            .unwrap_or_else(|error| {
                                panic!("fixed W06 I08 dispatch succeeds: rows={row_count}, iterations={max_iterations}, error={error:?}")
                            });
                            actual_calls += 1;
                            match dispatch {
                                super::PreparedConformerDispatchOutcome::Serial(result) => {
                                    panic!("fixed W06 I08 dispatch selected Serial: {result:?}")
                                }
                                super::PreparedConformerDispatchOutcome::Workers(result) => {
                                    let joins = result.unwrap_or_else(|error| {
                                        panic!("fixed W06 I08 worker setup succeeds: {error:?}")
                                    });
                                    assert_eq!(joins.len(), 3);
                                    for (lane, joined) in joins.into_iter().enumerate() {
                                        joined.unwrap_or_else(|panic| {
                                            panic!("fixed W06 I08 worker lane {lane} joins: {panic:?}")
                                        })
                                        .unwrap_or_else(|error| {
                                            panic!("fixed W06 I08 worker lane {lane} succeeds: {error:?}")
                                        });
                                    }
                                }
                            }
                            results
                        }
                    };
                    assert_eq!(
                        super::super::api::prepare_parameter_query_calls(),
                        1,
                        "one shared preparation per entry: mode={mode}, rows={row_count}, iterations={max_iterations}"
                    );

                    assert_borrowed_inputs_unchanged();
                    assert!(diagnostics.is_empty());
                    assert_eq!(coordinates.source_coordinate_dim, before_dimension);
                    assert_eq!(snapshot_2d(&coordinates), before_2d);
                    assert_eq!(
                        coordinate_storage(&coordinates),
                        before_coordinate_storage,
                        "W06 coordinate buffers remain borrowed in place: mode={mode}, rows={row_count}, iterations={max_iterations}"
                    );
                    let after_3d = snapshot_3d(&coordinates);
                    assert_eq!(after_3d.len(), row_count);
                    for row in 0..row_count {
                        assert_eq!(after_3d[row].0, IDS[row]);
                        assert_eq!(after_3d[row].1, before_3d[row].1);
                        assert_eq!(after_3d[row].2, before_3d[row].2);
                        for point in 0..2 {
                            for axis in 0..3 {
                                assert_eq!(
                                    after_3d[row].3[point][axis],
                                    expected_rows[row].0[point][axis].to_bits(),
                                    "fixed W06 I08 source coordinate: mode={mode}, row={row}, point={point}, axis={axis}, iterations={max_iterations}"
                                );
                            }
                        }
                    }

                    if mode == 0 {
                        assert_eq!(outcomes.len(), 1);
                        assert_eq!(outcomes[0].status, 1);
                        assert!(
                            (outcomes[0].energy - expected_rows[0].1).abs() <= 1.0e-12,
                            "fixed W06 source single energy: rows={row_count}, iterations={max_iterations}"
                        );
                    } else {
                        assert_eq!(outcomes.len(), row_count);
                        for row in 0..row_count {
                            assert_eq!(outcomes[row].status, 1);
                            assert!(
                                (outcomes[row].energy - expected_rows[row].1).abs() <= 1.0e-12,
                                "fixed W06 source all-row energy: mode={mode}, row={row}, iterations={max_iterations}"
                            );
                        }
                    }
                }
            }
        }

        assert_eq!(actual_calls, 18);
    }

    #[test]
    fn uff_error_e03_serial_leaf_errors_preserve_fields() {
        // These manually supplied values verify Error trait dispatch only;
        // actual capacity and coordinate failures remain covered by S03/S04.
        let capacity = super::SerialResultCapacityError {
            conformers: 23,
            result_slots: 7,
        };
        let erased: &(dyn std::error::Error + 'static) = &capacity;
        let downcast = erased
            .downcast_ref::<super::SerialResultCapacityError>()
            .unwrap();
        assert!(std::ptr::eq(downcast, &capacity));
        assert_eq!((downcast.conformers, downcast.result_slots), (23, 7));
        assert!(erased.source().is_none());
        assert_eq!(capacity.to_string(), format!("{capacity:?}"));

        let coordinate = super::SerialCoordinateCountError {
            conformer_id: 907,
            atoms: 19,
            coordinates: 11,
        };
        let erased: &(dyn std::error::Error + 'static) = &coordinate;
        let downcast = erased
            .downcast_ref::<super::SerialCoordinateCountError>()
            .unwrap();
        assert!(std::ptr::eq(downcast, &coordinate));
        assert_eq!(
            (downcast.conformer_id, downcast.atoms, downcast.coordinates,),
            (907, 19, 11)
        );
        assert!(erased.source().is_none());
        assert_eq!(coordinate.to_string(), format!("{coordinate:?}"));
    }

    #[test]
    fn uff_error_e09_stage_and_serial_causes_borrow_stored_children() {
        fn assert_stored_child<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let exposed = std::error::Error::source(parent)
                .expect("the source-bearing variant exposes its stored child")
                .downcast_ref::<Child>()
                .expect("the source keeps the child's concrete type");
            assert!(std::ptr::eq(stored, exposed));
        }

        // These hand-built values prove trait dispatch only. Actual capacity,
        // coordinate, stage, and selected-conformer failures remain covered
        // by S03/S04, S11/T07, and U09 caller regressions.
        let stage_cases = [
            (
                super::OptimizationStageError::Initialize(super::ForceFieldKernelError::NoPoints),
                "Initialize(NoPoints)",
            ),
            (
                super::OptimizationStageError::Minimize(
                    super::ForceFieldKernelError::IndexOutOfRange {
                        argument: super::super::super::kernel::ForceFieldIndexArgument::J,
                        index: 17,
                        upper_bound: 23,
                    },
                ),
                "Minimize(IndexOutOfRange { argument: J, index: 17, upper_bound: 23 })",
            ),
            (
                super::OptimizationStageError::FinalEnergy(
                    super::ForceFieldKernelError::OptimizerBadDirection,
                ),
                "FinalEnergy(OptimizerBadDirection)",
            ),
        ];
        for (stage, expected_display) in stage_cases {
            assert_eq!(stage.to_string(), expected_display);
            let kernel_error = match &stage {
                super::OptimizationStageError::Initialize(stored)
                | super::OptimizationStageError::Minimize(stored)
                | super::OptimizationStageError::FinalEnergy(stored) => stored,
            };
            assert_stored_child(&stage, kernel_error);
            assert!(std::error::Error::source(kernel_error).is_none());
        }

        let capacity = super::SerialConformerOptimizationError::ResultCapacity(
            super::SerialResultCapacityError {
                conformers: 41,
                result_slots: 37,
            },
        );
        assert_eq!(
            capacity.to_string(),
            "ResultCapacity(SerialResultCapacityError { conformers: 41, result_slots: 37 })"
        );
        let capacity_child = match &capacity {
            super::SerialConformerOptimizationError::ResultCapacity(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&capacity, capacity_child);
        assert!(std::error::Error::source(capacity_child).is_none());

        let coordinate = super::SerialConformerOptimizationError::CoordinateCount {
            input_index: 29,
            source: super::SerialCoordinateCountError {
                conformer_id: 907,
                atoms: 17,
                coordinates: 13,
            },
        };
        assert_eq!(
            coordinate.to_string(),
            "CoordinateCount { input_index: 29, source: SerialCoordinateCountError { conformer_id: 907, atoms: 17, coordinates: 13 } }"
        );
        let (input_index, coordinate_child) = match &coordinate {
            super::SerialConformerOptimizationError::CoordinateCount {
                input_index,
                source,
            } => (input_index, source),
            _ => unreachable!(),
        };
        assert_eq!(*input_index, 29);
        assert_eq!(coordinate_child.conformer_id, 907);
        assert_stored_child(&coordinate, coordinate_child);
        assert!(std::error::Error::source(coordinate_child).is_none());

        let optimization = super::SerialConformerOptimizationError::Optimization {
            input_index: 43,
            conformer_id: 911,
            source: super::OptimizationStageError::FinalEnergy(
                super::ForceFieldKernelError::OptimizerBadDirection,
            ),
        };
        assert_eq!(
            optimization.to_string(),
            "Optimization { input_index: 43, conformer_id: 911, source: FinalEnergy(OptimizerBadDirection) }"
        );
        let (input_index, conformer_id, stage_child) = match &optimization {
            super::SerialConformerOptimizationError::Optimization {
                input_index,
                conformer_id,
                source,
            } => (input_index, conformer_id, source),
            _ => unreachable!(),
        };
        assert_eq!((*input_index, *conformer_id), (43, 911));
        assert_stored_child(&optimization, stage_child);
        let kernel_error = match stage_child {
            super::OptimizationStageError::FinalEnergy(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(stage_child, kernel_error);
        assert!(std::error::Error::source(kernel_error).is_none());
    }

    #[test]
    fn uff_error_e10_single_and_serial_wrappers_borrow_stored_causes() {
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        fn assert_stored_child<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let exposed = std::error::Error::source(parent)
                .expect("the source-bearing variant exposes its stored child")
                .downcast_ref::<Child>()
                .expect("the source keeps the child's concrete type");
            assert!(std::ptr::eq(stored, exposed));
        }

        // These hand-built values verify trait dispatch only. Actual selected
        // construction and optimizer failures remain covered by U09, while
        // the serial caller routes remain covered by S11/T07 and D06.
        let single_auto = super::AutomaticForceFieldConstructionError::Construction(
            ForceFieldConstructionError::Builder(
                UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 73 },
            ),
        );
        let single_construction =
            super::SingleConformerOptimizationError::Construction(single_auto);
        assert_eq!(
            single_construction.to_string(),
            "Construction(Construction(Builder(SelectedThreeDimensionalConformerNotFound { conformer_id: 73 })))"
        );
        let single_auto_child = match &single_construction {
            super::SingleConformerOptimizationError::Construction(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&single_construction, single_auto_child);
        let single_construction_child = match single_auto_child {
            super::AutomaticForceFieldConstructionError::Construction(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(single_auto_child, single_construction_child);
        let single_builder_child = match single_construction_child {
            ForceFieldConstructionError::Builder(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(single_construction_child, single_builder_child);
        assert!(matches!(
            single_builder_child,
            UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 73 }
        ));
        assert!(std::error::Error::source(single_builder_child).is_none());

        let single_kernel = super::ForceFieldKernelError::IndexOutOfRange {
            argument: super::super::super::kernel::ForceFieldIndexArgument::J,
            index: 17,
            upper_bound: 23,
        };
        let single_stage = super::OptimizationStageError::Minimize(single_kernel);
        let single_optimization =
            super::SingleConformerOptimizationError::Optimization(single_stage);
        assert_eq!(
            single_optimization.to_string(),
            "Optimization(Minimize(IndexOutOfRange { argument: J, index: 17, upper_bound: 23 }))"
        );
        let single_stage_child = match &single_optimization {
            super::SingleConformerOptimizationError::Optimization(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&single_optimization, single_stage_child);
        let single_kernel_child = match single_stage_child {
            super::OptimizationStageError::Minimize(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(single_stage_child, single_kernel_child);
        assert!(matches!(
            single_kernel_child,
            super::ForceFieldKernelError::IndexOutOfRange {
                argument: super::super::super::kernel::ForceFieldIndexArgument::J,
                index: 17,
                upper_bound: 23,
            }
        ));
        assert!(std::error::Error::source(single_kernel_child).is_none());

        let serial_auto = super::AutomaticForceFieldConstructionError::Construction(
            ForceFieldConstructionError::Builder(
                UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 991 },
            ),
        );
        let serial_construction = super::SerialUffOptimizationError::Construction(serial_auto);
        assert_eq!(
            serial_construction.to_string(),
            "Construction(Construction(Builder(SelectedThreeDimensionalConformerNotFound { conformer_id: 991 })))"
        );
        let serial_auto_child = match &serial_construction {
            super::SerialUffOptimizationError::Construction(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&serial_construction, serial_auto_child);
        let serial_construction_child = match serial_auto_child {
            super::AutomaticForceFieldConstructionError::Construction(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(serial_auto_child, serial_construction_child);
        let serial_builder_child = match serial_construction_child {
            ForceFieldConstructionError::Builder(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(serial_construction_child, serial_builder_child);
        assert!(matches!(
            serial_builder_child,
            UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 991 }
        ));
        assert!(std::error::Error::source(serial_builder_child).is_none());

        let serial_stage = super::OptimizationStageError::FinalEnergy(
            super::ForceFieldKernelError::OptimizerBadDirection,
        );
        let serial_optimization = super::SerialUffOptimizationError::Optimization(
            super::SerialConformerOptimizationError::Optimization {
                input_index: 43,
                conformer_id: 911,
                source: serial_stage,
            },
        );
        assert_eq!(
            serial_optimization.to_string(),
            "Optimization(Optimization { input_index: 43, conformer_id: 911, source: FinalEnergy(OptimizerBadDirection) })"
        );
        let serial_row_child = match &serial_optimization {
            super::SerialUffOptimizationError::Optimization(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&serial_optimization, serial_row_child);
        let (input_index, conformer_id, serial_stage_child) = match serial_row_child {
            super::SerialConformerOptimizationError::Optimization {
                input_index,
                conformer_id,
                source,
            } => (input_index, conformer_id, source),
            _ => unreachable!(),
        };
        assert_eq!((*input_index, *conformer_id), (43, 911));
        assert_stored_child(serial_row_child, serial_stage_child);
        let serial_kernel_child = match serial_stage_child {
            super::OptimizationStageError::FinalEnergy(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(serial_stage_child, serial_kernel_child);
        assert!(matches!(
            serial_kernel_child,
            super::ForceFieldKernelError::OptimizerBadDirection
        ));
        assert!(std::error::Error::source(serial_kernel_child).is_none());
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_error_e10_dispatched_wrapper_borrows_stored_causes() {
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        fn assert_stored_child<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let exposed = std::error::Error::source(parent)
                .expect("the source-bearing variant exposes its stored child")
                .downcast_ref::<Child>()
                .expect("the source keeps the child's concrete type");
            assert!(std::ptr::eq(stored, exposed));
        }

        // These hand-built values verify trait dispatch only. The real
        // constructor and resolver failures remain covered by dispatch D06.
        let automatic = super::AutomaticForceFieldConstructionError::Construction(
            ForceFieldConstructionError::Builder(
                UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 257 },
            ),
        );
        let construction = super::DispatchedUffOptimizationError::Construction(automatic);
        assert_eq!(
            construction.to_string(),
            "UFF force-field construction failed: Construction(Builder(SelectedThreeDimensionalConformerNotFound { conformer_id: 257 }))"
        );
        let automatic_child = match &construction {
            super::DispatchedUffOptimizationError::Construction(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&construction, automatic_child);
        let construction_child = match automatic_child {
            super::AutomaticForceFieldConstructionError::Construction(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(automatic_child, construction_child);
        let builder_child = match construction_child {
            ForceFieldConstructionError::Builder(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(construction_child, builder_child);
        assert!(matches!(
            builder_child,
            UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 257 }
        ));
        assert!(std::error::Error::source(builder_child).is_none());

        let thread_count = super::UffThreadCountError::UndefinedSignedNegation;
        let thread_failure = super::DispatchedUffOptimizationError::ThreadCount(thread_count);
        assert_eq!(
            thread_failure.to_string(),
            "UFF conformer thread-count resolution failed: UndefinedSignedNegation"
        );
        let thread_child = match &thread_failure {
            super::DispatchedUffOptimizationError::ThreadCount(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&thread_failure, thread_child);
        assert!(matches!(
            thread_child,
            super::UffThreadCountError::UndefinedSignedNegation
        ));
        assert!(std::error::Error::source(thread_child).is_none());
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_error_e11_same_block_selected_id_failure_keeps_the_real_chain() {
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        fn assert_stored_child<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let exposed = std::error::Error::source(parent)
                .expect("the source-bearing variant exposes its stored child")
                .downcast_ref::<Child>()
                .expect("the source keeps the child's concrete type");
            assert!(std::ptr::eq(stored, exposed));
        }

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let mut coordinates = worker_w06_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let coordinates_before = coordinates.clone();
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; 2],
        };
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut results = vec![sentinel];
        let mut diagnostics = Vec::new();
        let mut options = super::SingleConformerOptions::for_conformer(99);
        options.ignore_interfragment_interactions = false;

        let error = match super::optimize_dispatched_uff_coordinate_block(
            &topology,
            &mut coordinates,
            &mut results,
            super::UffAtomStateRef::SuppliedRows {
                total_valences: &total_valences,
                conjugated_presence: &conjugated,
            },
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            options,
            i32::MIN,
            0,
            true,
        ) {
            Err(error) => error,
            Ok(_) => panic!("the absent selected 3D ID fails in the real constructor"),
        };

        let automatic = match &error {
            super::DispatchedUffOptimizationError::Construction(stored) => stored,
            _ => unreachable!("construction error precedes dispatch"),
        };
        assert_stored_child(&error, automatic);
        let construction = match automatic {
            super::AutomaticForceFieldConstructionError::Construction(stored) => stored,
            _ => unreachable!("selected conformer lookup fails in field construction"),
        };
        assert_stored_child(automatic, construction);
        let builder = match construction {
            ForceFieldConstructionError::Builder(stored) => stored,
            _ => unreachable!("selected conformer lookup is a builder failure"),
        };
        assert_stored_child(construction, builder);
        assert!(matches!(
            builder,
            UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 99 }
        ));
        assert!(std::error::Error::source(builder).is_none());
        assert_eq!(results, [sentinel]);
        assert_eq!(coordinates, coordinates_before);
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn uff_error_e11_automatic_prepared_length_failure_keeps_the_real_child() {
        use super::super::atom_typer::{UffTypingError, UffTypingInput};

        fn assert_stored_child<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let exposed = std::error::Error::source(parent)
                .expect("the source-bearing variant exposes its stored child")
                .downcast_ref::<Child>()
                .expect("the source keeps the child's concrete type");
            assert!(std::ptr::eq(stored, exposed));
        }

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let mut coordinates = worker_w06_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let total_valences = [0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0; 2],
        };
        let error = match super::construct_force_field_with_automatic_typing(
            &topology,
            &mut coordinates,
            31,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut Vec::new(),
            super::DEFAULT_TORSION_BOND_SMARTS,
            10.0,
            false,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("the actual automatic constructor rejects the short prepared row")
            }
        };

        let typing = match &error {
            super::AutomaticForceFieldConstructionError::Typing(stored) => stored,
            _ => unreachable!("prepared row length is checked by atom typing"),
        };
        assert_stored_child(&error, typing);
        assert!(matches!(
            typing,
            UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 2,
                actual: 1,
            }
        ));
        assert!(std::error::Error::source(typing).is_none());
    }

    #[test]
    fn uff_error_e11_worker_stage_failures_keep_actual_context_and_causes() {
        fn assert_stored_child<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let exposed = std::error::Error::source(parent)
                .expect("the source-bearing variant exposes its stored child")
                .downcast_ref::<Child>()
                .expect("the source keeps the child's concrete type");
            assert!(std::ptr::eq(stored, exposed));
        }

        const IDS: [usize; 3] = [80, 7, 900];
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        // This is the retained W05 no-point initialization failure routed
        // through the real worker and serial-stage mapper.
        let mut empty_rows = vec![Vec::<[f64; 3]>::new(); IDS.len()];
        let mut initialization_conformers = empty_rows
            .iter_mut()
            .enumerate()
            .map(|(input_index, positions)| super::SerialConformer {
                id: IDS[input_index],
                positions: positions.as_mut_slice(),
            })
            .collect::<Vec<_>>();
        let mut initialization_results = vec![sentinel; IDS.len()];
        let initialization_error = match super::optimize_worker_conformers(
            super::ForceField::new(3),
            &mut initialization_conformers,
            &mut initialization_results,
            0,
            1,
            std::num::NonZeroU32::new(2).expect("fixed W05 lane count is nonzero"),
            0,
        ) {
            Err(error) => error,
            Ok(()) => panic!("the retained empty W05 row fails during initialization"),
        };
        let (input_index, conformer_id, initialization_stage) = match &initialization_error {
            super::SerialConformerOptimizationError::Optimization {
                input_index,
                conformer_id,
                source,
            } => (input_index, conformer_id, source),
            _ => unreachable!("the real worker reports a selected stage failure"),
        };
        assert_eq!((*input_index, *conformer_id), (1, IDS[1]));
        assert_stored_child(&initialization_error, initialization_stage);
        let initialization_kernel = match initialization_stage {
            super::OptimizationStageError::Initialize(stored) => stored,
            _ => unreachable!("the empty W05 row fails at initialize"),
        };
        assert_stored_child(initialization_stage, initialization_kernel);
        assert!(matches!(
            initialization_kernel,
            super::ForceFieldKernelError::NoPoints
        ));
        assert!(std::error::Error::source(initialization_kernel).is_none());
        drop(initialization_conformers);
        assert_eq!(initialization_results, [sentinel; IDS.len()]);

        // Reuse the exact W05 call-count rule: the first selected row consumes
        // two energy calls; the second selected row fails during minimize at
        // call3 or during final energy at call4.
        for (fail_on, expected_phase) in [(3_usize, 0_u8), (4, 1)] {
            let mut coordinate_rows = (0..4)
                .map(|source_index| {
                    if source_index % 2 == 0 {
                        let selected_ordinal = source_index / 2;
                        vec![[0.0, 0.0, 0.0], [2.0 + selected_ordinal as f64, 0.0, 0.0]]
                    } else {
                        Vec::new()
                    }
                })
                .collect::<Vec<_>>();
            let row_count = coordinate_rows.len();
            let mut conformers = coordinate_rows
                .iter_mut()
                .enumerate()
                .map(|(input_index, positions)| super::SerialConformer {
                    id: [IDS[0], IDS[1], IDS[2], 31][input_index],
                    positions: positions.as_mut_slice(),
                })
                .collect::<Vec<_>>();
            let mut results = vec![sentinel; row_count];
            let mut field = super::ForceField::new(3);
            field.add_contribution(Box::new(FailAtEnergyCallWithCoordinateNorm {
                calls: Cell::new(0),
                fail_on,
            }));
            let error = match super::optimize_worker_conformers(
                field,
                &mut conformers,
                &mut results,
                2,
                0,
                std::num::NonZeroU32::new(2).expect("fixed W05 lane count is nonzero"),
                0,
            ) {
                Err(error) => error,
                Ok(()) => panic!("the retained W05 energy-call failure must be reported"),
            };
            let (input_index, conformer_id, stage) = match &error {
                super::SerialConformerOptimizationError::Optimization {
                    input_index,
                    conformer_id,
                    source,
                } => (input_index, conformer_id, source),
                _ => unreachable!("the real worker reports a selected stage failure"),
            };
            assert_eq!((*input_index, *conformer_id), (2, IDS[2]));
            assert_stored_child(&error, stage);
            let kernel_error = match (expected_phase, stage) {
                (0, super::OptimizationStageError::Minimize(stored))
                | (1, super::OptimizationStageError::FinalEnergy(stored)) => stored,
                _ => unreachable!("W05 call3 is minimize and call4 is final energy"),
            };
            assert_stored_child(stage, kernel_error);
            assert_eq!(*kernel_error, expected_pair_index_error());
            assert!(std::error::Error::source(kernel_error).is_none());
            drop(conformers);
        }
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_dispatch_d03_raw_carrier_preserves_typed_errors_and_panic_payload() {
        let initialize_failure = super::SerialConformerOptimizationError::Optimization {
            input_index: 0,
            conformer_id: 101,
            source: super::OptimizationStageError::Initialize(
                super::ForceFieldKernelError::NoPoints,
            ),
        };
        let minimize_failure = super::SerialConformerOptimizationError::Optimization {
            input_index: 1,
            conformer_id: 203,
            source: super::OptimizationStageError::Minimize(
                super::ForceFieldKernelError::IndexOutOfRange {
                    argument: super::super::super::kernel::ForceFieldIndexArgument::J,
                    index: 2,
                    upper_bound: 2,
                },
            ),
        };
        let final_energy_failure = super::SerialConformerOptimizationError::Optimization {
            input_index: 4,
            conformer_id: 509,
            source: super::OptimizationStageError::FinalEnergy(
                super::ForceFieldKernelError::IndexOutOfRange {
                    argument: super::super::super::kernel::ForceFieldIndexArgument::I,
                    index: 2,
                    upper_bound: 2,
                },
            ),
        };

        for expected in [initialize_failure, minimize_failure, final_energy_failure] {
            match super::PreparedConformerDispatchOutcome::Serial(Err(expected)) {
                super::PreparedConformerDispatchOutcome::Serial(Err(actual)) => {
                    assert_eq!(actual, expected);
                }
                _ => unreachable!("serial carrier retains the original typed failure"),
            }
        }

        let setup_failure = super::SerialConformerOptimizationError::CoordinateCount {
            input_index: 6,
            source: super::SerialCoordinateCountError {
                conformer_id: 907,
                atoms: 3,
                coordinates: 2,
            },
        };
        match super::PreparedConformerDispatchOutcome::Workers(Err(setup_failure)) {
            super::PreparedConformerDispatchOutcome::Workers(Err(actual)) => {
                assert_eq!(actual, setup_failure);
            }
            _ => unreachable!("worker carrier retains the pre-spawn typed failure"),
        }

        const PANIC_MARKER: &str = "uff-dispatch-original-worker-panic";
        let worker_outcomes = std::thread::scope(|scope| {
            let handles = vec![
                scope.spawn(|| -> Result<(), super::SerialConformerOptimizationError> { Ok(()) }),
                scope.spawn(
                    move || -> Result<(), super::SerialConformerOptimizationError> {
                        Err(minimize_failure)
                    },
                ),
                scope.spawn(|| -> Result<(), super::SerialConformerOptimizationError> {
                    std::panic::panic_any(String::from(PANIC_MARKER))
                }),
            ];
            super::join_worker_handles(handles)
        });
        match super::PreparedConformerDispatchOutcome::Workers(Ok(worker_outcomes)) {
            super::PreparedConformerDispatchOutcome::Workers(Ok(outcomes)) => {
                assert_eq!(outcomes.len(), 3);
                assert!(matches!(outcomes.get(0), Some(Ok(Ok(())))));
                assert!(matches!(
                    outcomes.get(1),
                    Some(Ok(Err(actual))) if *actual == minimize_failure
                ));
                let panic_payload = outcomes
                    .get(2)
                    .expect("third creation-order handle carries the panic")
                    .as_ref()
                    .expect_err("third worker preserves its raw join panic");
                assert_eq!(
                    panic_payload.downcast_ref::<String>().map(String::as_str),
                    Some(PANIC_MARKER)
                );
            }
            _ => unreachable!("worker carrier retains the complete join vector"),
        }
    }

    #[test]
    fn uff_dispatch_d01_resize_preserves_prefix_truncates_and_appends_defaults() {
        const VECTOR_LENGTHS: [usize; 4] = [0, 1, 3, 5];
        const PREFIX_SENTINELS: [(i32, u64); 5] = [
            (19, 0x3ff8_0000_0000_0000),
            (-7, 0xc004_0000_0000_0000),
            (3, 0x4009_21fb_5444_2d18),
            (-11, 0x0000_0000_0000_0001),
            (28, 0xbff0_0000_0000_0000),
        ];

        let mut actual_resize_calls = 0;
        for initial_len in VECTOR_LENGTHS {
            for requested_len in VECTOR_LENGTHS {
                let mut results: Vec<_> = PREFIX_SENTINELS[..initial_len]
                    .iter()
                    .map(|&(status, energy_bits)| super::OptimizationOutcome {
                        status,
                        energy: f64::from_bits(energy_bits),
                    })
                    .collect();

                super::resize_conformer_results(&mut results, requested_len);
                assert_eq!(results.len(), requested_len);

                let preserved_prefix_len = initial_len.min(requested_len);
                for (index, &(expected_status, expected_energy_bits)) in
                    PREFIX_SENTINELS[..preserved_prefix_len].iter().enumerate()
                {
                    assert_eq!(results[index].status, expected_status, "row {index}");
                    assert_eq!(
                        results[index].energy.to_bits(),
                        expected_energy_bits,
                        "row {index} energy bits"
                    );
                }

                for index in preserved_prefix_len..requested_len {
                    assert_eq!(results[index].status, 0, "appended row {index}");
                    assert_eq!(
                        results[index].energy.to_bits(),
                        0.0_f64.to_bits(),
                        "appended row {index} must contain positive zero"
                    );
                }
                actual_resize_calls += 1;
            }
        }
        assert_eq!(actual_resize_calls, 16);
    }

    #[test]
    fn uff_dispatch_d02_thread_route_uses_literal_source_count_matrix() {
        #[derive(Clone, Copy)]
        enum ExpectedRoute {
            Serial,
            Workers(i32),
            UndefinedSignedNegation,
        }

        const REQUESTED: [i32; 10] = [
            -2_147_483_648,
            -2_147_483_647,
            -7,
            -2,
            -1,
            0,
            1,
            2,
            7,
            2_147_483_647,
        ];
        const HARDWARE: [u32; 8] = [0, 1, 2, 4, 8, 2_147_483_647, 2_147_483_648, 4_294_967_295];
        const EXPECTED_THREADSAFE: [[ExpectedRoute; 8]; 10] = [
            [ExpectedRoute::UndefinedSignedNegation; 8],
            [
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Workers(i32::MIN),
            ],
            [
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Workers(2_147_483_640),
                ExpectedRoute::Workers(2_147_483_641),
                ExpectedRoute::Workers(-8),
            ],
            [
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Workers(2),
                ExpectedRoute::Workers(6),
                ExpectedRoute::Workers(2_147_483_645),
                ExpectedRoute::Workers(2_147_483_646),
                ExpectedRoute::Workers(-3),
            ],
            [
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Workers(3),
                ExpectedRoute::Workers(7),
                ExpectedRoute::Workers(2_147_483_646),
                ExpectedRoute::Workers(2_147_483_647),
                ExpectedRoute::Workers(-2),
            ],
            [
                ExpectedRoute::Serial,
                ExpectedRoute::Serial,
                ExpectedRoute::Workers(2),
                ExpectedRoute::Workers(4),
                ExpectedRoute::Workers(8),
                ExpectedRoute::Workers(2_147_483_647),
                ExpectedRoute::Workers(i32::MIN),
                ExpectedRoute::Workers(-1),
            ],
            [ExpectedRoute::Serial; 8],
            [ExpectedRoute::Workers(2); 8],
            [ExpectedRoute::Workers(7); 8],
            [ExpectedRoute::Workers(i32::MAX); 8],
        ];

        // FFConvenience.h assigns the RDThreads.h unsigned return to signed
        // int before its ==1 branch. The vendored CMakeLists.txt requires
        // C++20; the current x86_64-linux-gnu target has 32-bit int, and
        // C++20 [conv.integral]/3 defines the signed conversion modulo 2^32.
        // These literal expectations freeze that source arithmetic without
        // asking the Rust resolver to produce the expected values.
        let mut actual_route_calls = 0;
        for (requested_index, &requested) in REQUESTED.iter().enumerate() {
            for (hardware_index, &observed_hardware) in HARDWARE.iter().enumerate() {
                let no_thread_route =
                    super::resolve_conformer_dispatch(requested, observed_hardware, false);
                assert_eq!(
                    no_thread_route,
                    Ok(super::ConformerDispatchRoute::Serial),
                    "non-threadsafe target={requested}, hardware={observed_hardware}"
                );
                actual_route_calls += 1;

                let threaded_route =
                    super::resolve_conformer_dispatch(requested, observed_hardware, true);
                match EXPECTED_THREADSAFE[requested_index][hardware_index] {
                    ExpectedRoute::Serial => assert_eq!(
                        threaded_route,
                        Ok(super::ConformerDispatchRoute::Serial),
                        "threadsafe serial target={requested}, hardware={observed_hardware}"
                    ),
                    ExpectedRoute::Workers(source_count) => assert_eq!(
                        threaded_route,
                        Ok(super::ConformerDispatchRoute::Workers { source_count }),
                        "threadsafe worker target={requested}, hardware={observed_hardware}"
                    ),
                    ExpectedRoute::UndefinedSignedNegation => assert_eq!(
                        threaded_route,
                        Err(super::UffThreadCountError::UndefinedSignedNegation),
                        "threadsafe INT_MIN target={requested}, hardware={observed_hardware}"
                    ),
                }
                actual_route_calls += 1;
            }
        }
        assert_eq!(actual_route_calls, 160);
    }

    #[test]
    fn uff_thread_t02_send_compile_proofs_cover_remaining_uff_terms() {
        fn assert_send<T: Send>() {}

        assert_send::<super::super::angle::AngleBendContrib>();
        assert_send::<super::super::inversion::InversionContrib>();
        assert_send::<super::super::nonbonded::VdwContrib>();
        assert_send::<super::super::torsion::TorsionAngleContrib>();
    }

    #[test]
    fn uff_thread_t03_thread_count_fixed_matrix_and_minimum_boundary() {
        const OBSERVED_HARDWARE: [u32; 6] = [0, 1, 2, 4, 8, u32::MAX];
        const CASES: [(i32, [u32; 6], [u32; 6]); 9] = [
            (1, [1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]),
            (2, [1, 1, 1, 1, 1, 1], [2, 2, 2, 2, 2, 2]),
            (7, [1, 1, 1, 1, 1, 1], [7, 7, 7, 7, 7, 7]),
            (
                i32::MAX,
                [1, 1, 1, 1, 1, 1],
                [
                    2_147_483_647,
                    2_147_483_647,
                    2_147_483_647,
                    2_147_483_647,
                    2_147_483_647,
                    2_147_483_647,
                ],
            ),
            (0, [1, 1, 1, 1, 1, 1], [1, 1, 2, 4, 8, 4_294_967_295]),
            (-1, [1, 1, 1, 1, 1, 1], [1, 1, 1, 3, 7, 4_294_967_294]),
            (-2, [1, 1, 1, 1, 1, 1], [1, 1, 1, 2, 6, 4_294_967_293]),
            (-7, [1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 4_294_967_288]),
            (
                -i32::MAX,
                [1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 2_147_483_648],
            ),
        ];

        let mut defined_input_calls = 0;
        for (requested, no_thread_expected, threaded_expected) in CASES {
            for (hardware_index, observed_hardware) in OBSERVED_HARDWARE.into_iter().enumerate() {
                let no_thread =
                    super::resolve_uff_thread_count(requested, observed_hardware, false)
                        .expect("no-thread source overload always returns one");
                assert_ne!(no_thread.get(), 0);
                assert_eq!(
                    no_thread.get(),
                    no_thread_expected[hardware_index],
                    "no-thread mode: target={requested}, hardware={observed_hardware}"
                );
                defined_input_calls += 1;

                let threaded = super::resolve_uff_thread_count(requested, observed_hardware, true)
                    .expect("all matrix targets are source-defined in threaded mode");
                assert_ne!(threaded.get(), 0);
                assert_eq!(
                    threaded.get(),
                    threaded_expected[hardware_index],
                    "threaded mode: target={requested}, hardware={observed_hardware}"
                );
                defined_input_calls += 1;
            }
        }
        assert_eq!(defined_input_calls, 108);

        let mut boundary_calls = 0;
        let no_thread_min = super::resolve_uff_thread_count(i32::MIN, 1234, false)
            .expect("no-thread mode ignores target");
        assert_eq!(no_thread_min.get(), 1);
        boundary_calls += 1;

        assert_eq!(
            super::resolve_uff_thread_count(i32::MIN, 0, true),
            Err(super::UffThreadCountError::UndefinedSignedNegation)
        );
        boundary_calls += 1;
        assert_eq!(boundary_calls, 2);
    }

    fn worker_w06_atom(
        row: usize,
        atomic_number: u8,
        hybridization: Hybridization,
        formal_charge: i8,
    ) -> Atom {
        Atom::from_spec(
            AtomId::new(row),
            AtomSpec::new(
                Element::from_atomic_number(atomic_number)
                    .expect("fixed W06 element exists in the element table"),
            )
            .with_hybridization(hybridization)
            .with_formal_charge(formal_charge)
            .with_no_implicit(true),
        )
    }

    fn worker_w06_topology(atoms: Vec<Atom>, edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(row, &(begin, end, order))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed W06 topology is structurally valid")
    }

    fn worker_w06_coordinates(points: &[[f64; 3]]) -> CoordinateBlock {
        let mut coordinates = CoordinateBlock::default();
        coordinates
            .conformers_3d
            .push(Conformer3D::new(31, points.to_vec(), true));
        coordinates
    }

    fn worker_w06_source_field<'a>(
        topology: &TopologyBlock,
        coordinates: &'a mut CoordinateBlock,
        total_valences: &[i32],
        vdw_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> ForceField<'a> {
        let atom_count = topology.atoms.len();
        let rings = fast_find_rings(topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; atom_count],
        };
        let mut diagnostics = Vec::new();
        super::construct_force_field_with_automatic_typing(
            topology,
            coordinates,
            31,
            total_valences,
            &vec![false; atom_count],
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            super::DEFAULT_TORSION_BOND_SMARTS,
            vdw_threshold,
            ignore_interfragment_interactions,
        )
        .unwrap_or_else(|error| panic!("fixed W06 reference field constructs: {error:?}"))
    }

    #[cfg(not(target_family = "wasm"))]
    fn run_worker_on_one_scoped_thread<'field, 'rows, 'coordinates>(
        field: ForceField<'field>,
        conformers: &'rows mut [super::SerialConformer<'coordinates>],
        results: &'rows mut [OptimizationOutcome],
        atom_count: usize,
        thread_idx: u32,
        lane_count: std::num::NonZeroU32,
        max_iterations: i32,
    ) -> (
        Result<(), super::SerialConformerOptimizationError>,
        Vec<usize>,
    )
    where
        'field: 'rows,
        'coordinates: 'rows,
    {
        std::thread::scope(|scope| {
            scope
                .spawn(move || {
                    let result = super::optimize_worker_conformers(
                        field,
                        &mut *conformers,
                        &mut *results,
                        atom_count,
                        thread_idx,
                        lane_count,
                        max_iterations,
                    );
                    let observed_ids = conformers
                        .iter()
                        .map(|conformer| conformer.id)
                        .collect::<Vec<_>>();
                    (result, observed_ids)
                })
                .join()
                .expect("scoped worker joins before its borrowed rows expire")
        })
    }

    #[cfg(not(target_family = "wasm"))]
    fn assert_scoped_worker_reference_matrix(
        source_field: &ForceField<'_>,
        row_ids: &[usize],
        initial_rows: &[Vec<[f64; 3]>],
        expected_statuses: &[i32],
        expected_energies: &[f64],
        atom_count: usize,
        energy_tolerance: f64,
    ) {
        assert_eq!(row_ids.len(), 3);
        assert_eq!(initial_rows.len(), 3);
        assert_eq!(expected_statuses.len(), 3);
        assert_eq!(expected_energies.len(), 3);
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        for row_count in 0..=3 {
            for lane_count_value in 1..=3 {
                let lane_count =
                    std::num::NonZeroU32::new(lane_count_value).expect("positive fixed lane count");
                for thread_idx in 0..lane_count_value {
                    let mut coordinate_rows = initial_rows[..row_count].to_vec();
                    let before_bits = coordinate_rows
                        .iter()
                        .flatten()
                        .flatten()
                        .map(|coordinate| coordinate.to_bits())
                        .collect::<Vec<_>>();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .zip(row_ids[..row_count].iter().copied())
                        .map(|(positions, id)| super::SerialConformer {
                            id,
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; row_count + 1];
                    let field_copy = source_field.copy();
                    let (worker_result, observed_ids) = run_worker_on_one_scoped_thread(
                        field_copy,
                        &mut conformers,
                        &mut results,
                        atom_count,
                        thread_idx,
                        lane_count,
                        0,
                    );
                    worker_result.expect("fixed selected source rows optimize successfully");
                    assert_eq!(observed_ids, row_ids[..row_count]);
                    drop(conformers);

                    for input_index in 0..row_count {
                        if input_index as u32 % lane_count_value == thread_idx {
                            assert_eq!(results[input_index].status, expected_statuses[input_index]);
                            assert!(
                                (results[input_index].energy - expected_energies[input_index])
                                    .abs()
                                    <= energy_tolerance,
                                "rows={row_count}, lanes={lane_count_value}, lane={thread_idx}, input={input_index}"
                            );
                        } else {
                            assert_eq!(results[input_index], sentinel);
                        }
                    }
                    assert_eq!(
                        results[row_count], sentinel,
                        "extra result slot stays untouched"
                    );
                    let after_bits = coordinate_rows
                        .iter()
                        .flatten()
                        .flatten()
                        .map(|coordinate| coordinate.to_bits())
                        .collect::<Vec<_>>();
                    assert_eq!(
                        after_bits, before_bits,
                        "zero-iteration rows stay bit-exact"
                    );
                }
            }
        }
    }

    #[test]
    fn uff_worker_w02_assignment_uses_source_indices_wraps_and_rejects_invalid_lanes() {
        let conformer_ids = [80, 2, 900, 17, 64];
        let lane_cases = [
            (1, vec![vec![80, 2, 900, 17, 64]]),
            (2, vec![vec![80, 900, 64], vec![2, 17]]),
            (3, vec![vec![80, 17], vec![2, 64], vec![900]]),
            (4, vec![vec![80, 64], vec![2], vec![900], vec![17]]),
            (
                7,
                vec![
                    vec![80],
                    vec![2],
                    vec![900],
                    vec![17],
                    vec![64],
                    vec![],
                    vec![],
                ],
            ),
        ];
        for (lane_count, expected_by_lane) in lane_cases {
            let lane_count = std::num::NonZeroU32::new(lane_count).unwrap();
            let mut assigned_per_row = [0; 5];
            for thread_idx in 0..lane_count.get() {
                let mut assignment = super::WorkerConformerAssignment::new(thread_idx, lane_count);
                let mut selected_ids = Vec::new();
                for (row_index, conformer_id) in conformer_ids.iter().copied().enumerate() {
                    if assignment.selects_current() {
                        selected_ids.push(conformer_id);
                        assigned_per_row[row_index] += 1;
                    }
                    assignment.advance();
                }
                assert_eq!(selected_ids, expected_by_lane[thread_idx as usize]);
                assert_eq!(assignment.source_index, conformer_ids.len() as u32);
            }
            assert_eq!(assigned_per_row, [1; 5]);

            for thread_idx in [lane_count.get(), lane_count.get() + 1, u32::MAX] {
                let mut assignment = super::WorkerConformerAssignment::new(thread_idx, lane_count);
                for _ in conformer_ids {
                    assert!(!assignment.selects_current());
                    assignment.advance();
                }
                assert_eq!(assignment.source_index, conformer_ids.len() as u32);
            }
        }

        let mut wrapping = super::WorkerConformerAssignment {
            source_index: u32::MAX,
            thread_idx: 3,
            lane_count: std::num::NonZeroU32::new(4).unwrap(),
        };
        assert!(wrapping.selects_current());
        wrapping.advance();
        assert_eq!(wrapping.source_index, 0);
        assert!(!wrapping.selects_current());
        wrapping.advance();
        assert_eq!(wrapping.source_index, 1);
        assert!(!wrapping.selects_current());
        wrapping.advance();
        assert_eq!(wrapping.source_index, 2);
        assert!(!wrapping.selects_current());
        wrapping.advance();
        assert_eq!(wrapping.source_index, 3);
        assert!(wrapping.selects_current());
    }

    #[test]
    fn uff_worker_w03_preflight_checks_slots_and_reserves_empty_or_unselected() {
        let empty: [super::SerialConformer<'_>; 0] = [];
        cf3d_uff_one_kernel_counts_reset();
        let mut prepared =
            super::prepare_worker_force_field(ForceField::new(3), empty.len(), 0, 3).unwrap();
        assert!(prepared.positions().is_empty());
        assert!(prepared.positions_mut().capacity() >= 3);
        drop(prepared);
        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 0);

        let mut coordinates = [[1.0, 2.0, 3.0]];
        let conformers = [super::SerialConformer {
            id: 800,
            positions: &mut coordinates,
        }];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let assignment =
            super::WorkerConformerAssignment::new(1, std::num::NonZeroU32::new(2).unwrap());
        assert!(!assignment.selects_current());

        for result_slots in [0, 1, 2] {
            let results = vec![sentinel; result_slots];
            let prepared = super::prepare_worker_force_field(
                ForceField::new(3),
                conformers.len(),
                results.len(),
                1,
            );
            if result_slots == 0 {
                assert!(matches!(
                    prepared,
                    Err(super::SerialConformerOptimizationError::ResultCapacity(
                        super::SerialResultCapacityError {
                            conformers: 1,
                            result_slots: 0
                        }
                    ))
                ));
            } else {
                let mut prepared = prepared.unwrap();
                assert!(prepared.positions().is_empty());
                assert!(prepared.positions_mut().capacity() >= 1);
                drop(prepared);
            }
            assert_eq!(results, vec![sentinel; result_slots]);
        }
        drop(conformers);
        assert_eq!(coordinates, [[1.0, 2.0, 3.0]]);
        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 0);
    }

    #[test]
    fn uff_mt_m01_capacity_count_matrix_preserves_rows_slots_and_reservation() {
        const CASES: [(usize, usize); 6] = [(0, 0), (0, 2), (1, 0), (1, 1), (3, 2), (3, 5)];
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        for (conformer_count, result_slots) in CASES {
            let mut coordinate_rows = (0..conformer_count)
                .map(|row| [[row as f64 + 1.0, 2.0, 3.0]])
                .collect::<Vec<_>>();
            let coordinate_bits_before = coordinate_rows
                .iter()
                .flat_map(|row| row[0].iter().copied())
                .map(f64::to_bits)
                .collect::<Vec<_>>();
            let mut conformers = coordinate_rows
                .iter_mut()
                .enumerate()
                .map(|(row, positions)| super::SerialConformer {
                    id: [800, 7, 900][row],
                    positions: positions.as_mut_slice(),
                })
                .collect::<Vec<_>>();
            let mut results = vec![sentinel; result_slots];
            let results_before = results.clone();
            let short_capacity = result_slots < conformer_count;

            cf3d_uff_one_kernel_counts_reset();
            let prepared = super::prepare_worker_force_field(
                ForceField::new(3),
                conformer_count,
                result_slots,
                3,
            );
            if short_capacity {
                assert!(matches!(
                    prepared,
                    Err(super::SerialConformerOptimizationError::ResultCapacity(
                        super::SerialResultCapacityError {
                            conformers,
                            result_slots: actual_slots,
                        }
                    )) if conformers == conformer_count && actual_slots == result_slots
                ));
            } else {
                let mut prepared = prepared.expect("sufficient count capacity is accepted");
                assert!(prepared.positions().is_empty());
                assert!(prepared.positions_mut().capacity() >= 3);
                drop(prepared);
            }

            drop(conformers);
            let coordinate_bits_after = coordinate_rows
                .iter()
                .flat_map(|row| row[0].iter().copied())
                .map(f64::to_bits)
                .collect::<Vec<_>>();
            assert_eq!(coordinate_bits_after, coordinate_bits_before);
            assert_eq!(results, results_before);

            let work = cf3d_uff_one_serial_work_counts();
            assert_eq!(work.initialize_calls, 0);
            assert_eq!(
                work.position_buffer_growths,
                usize::from(!short_capacity),
                "capacity failure precedes reserve for ({conformer_count}, {result_slots})"
            );
            if short_capacity {
                assert_eq!(work.initial_position_buffer_capacity, 0);
            } else {
                assert!(work.initial_position_buffer_capacity >= 3);
                assert_eq!(work.final_position_buffer_capacity, 0);
            }
            assert_eq!(cf3d_uff_one_kernel_counts(), (1, 0, 0));
        }
    }

    #[test]
    fn uff_mt_m02_wrapped_lane_cardinality_matches_literal_tables() {
        const SMALL_EXPECTED: [[[u64; 4]; 8]; 4] = [
            [
                [0, 0, 0, 0],
                [1, 0, 0, 0],
                [2, 0, 0, 0],
                [3, 0, 0, 0],
                [4, 0, 0, 0],
                [5, 0, 0, 0],
                [6, 0, 0, 0],
                [7, 0, 0, 0],
            ],
            [
                [0, 0, 0, 0],
                [1, 0, 0, 0],
                [1, 1, 0, 0],
                [2, 1, 0, 0],
                [2, 2, 0, 0],
                [3, 2, 0, 0],
                [3, 3, 0, 0],
                [4, 3, 0, 0],
            ],
            [
                [0, 0, 0, 0],
                [1, 0, 0, 0],
                [1, 1, 0, 0],
                [1, 1, 1, 0],
                [2, 1, 1, 0],
                [2, 2, 1, 0],
                [2, 2, 2, 0],
                [3, 2, 2, 0],
            ],
            [
                [0, 0, 0, 0],
                [1, 0, 0, 0],
                [1, 1, 0, 0],
                [1, 1, 1, 0],
                [1, 1, 1, 1],
                [2, 1, 1, 1],
                [2, 2, 1, 1],
                [2, 2, 2, 1],
            ],
        ];
        const LARGE_SOURCE_ROWS: [u64; 4] =
            [4_294_967_295, 4_294_967_296, 4_294_967_297, 8_589_934_595];
        const LARGE_EXPECTED: [[[u64; 4]; 4]; 4] = [
            [
                [4_294_967_295, 0, 0, 0],
                [2_147_483_648, 2_147_483_647, 0, 0],
                [1_431_655_765, 1_431_655_765, 1_431_655_765, 0],
                [1_073_741_824, 1_073_741_824, 1_073_741_824, 1_073_741_823],
            ],
            [
                [4_294_967_296, 0, 0, 0],
                [2_147_483_648, 2_147_483_648, 0, 0],
                [1_431_655_766, 1_431_655_765, 1_431_655_765, 0],
                [1_073_741_824, 1_073_741_824, 1_073_741_824, 1_073_741_824],
            ],
            [
                [4_294_967_297, 0, 0, 0],
                [2_147_483_649, 2_147_483_648, 0, 0],
                [1_431_655_767, 1_431_655_765, 1_431_655_765, 0],
                [1_073_741_825, 1_073_741_824, 1_073_741_824, 1_073_741_824],
            ],
            [
                [8_589_934_595, 0, 0, 0],
                [4_294_967_298, 4_294_967_297, 0, 0],
                [2_863_311_533, 2_863_311_531, 2_863_311_531, 0],
                [2_147_483_649, 2_147_483_649, 2_147_483_649, 2_147_483_648],
            ],
        ];
        const MAX_EXPECTED: [[u64; 4]; 4] = [
            [u64::MAX, 0, 0, 0],
            [9_223_372_036_854_775_808, 9_223_372_036_854_775_807, 0, 0],
            [
                6_148_914_694_099_828_735,
                6_148_914_689_804_861_440,
                6_148_914_689_804_861_440,
                0,
            ],
            [
                4_611_686_018_427_387_904,
                4_611_686_018_427_387_904,
                4_611_686_018_427_387_904,
                4_611_686_018_427_387_903,
            ],
        ];
        const WRAP_INDICES: [u32; 8] = [
            u32::MAX - 3,
            u32::MAX - 2,
            u32::MAX - 1,
            u32::MAX,
            0,
            1,
            2,
            3,
        ];
        const WRAP_EXPECTED_LANES: [(u32, [u32; 8]); 4] = [
            (1, [0, 0, 0, 0, 0, 0, 0, 0]),
            (2, [0, 1, 0, 1, 0, 1, 0, 1]),
            (3, [0, 1, 2, 0, 0, 1, 2, 0]),
            (4, [0, 1, 2, 3, 0, 1, 2, 3]),
        ];

        for source_rows in 0..=7_u64 {
            for lane_count_value in 1..=4_u32 {
                let lane_count = std::num::NonZeroU32::new(lane_count_value).unwrap();
                let expected =
                    SMALL_EXPECTED[(lane_count_value - 1) as usize][source_rows as usize];
                let mut total = 0_u64;
                for thread_idx in 0..lane_count_value {
                    let actual =
                        super::wrapped_source_lane_cardinality(source_rows, lane_count, thread_idx);
                    assert_eq!(
                        actual, expected[thread_idx as usize],
                        "small rows={source_rows}, lanes={lane_count_value}, lane={thread_idx}"
                    );
                    total = total.checked_add(actual).expect("small row sum fits u64");
                }
                assert_eq!(total, source_rows);
                assert_eq!(
                    super::wrapped_source_lane_cardinality(
                        source_rows,
                        lane_count,
                        lane_count_value,
                    ),
                    0,
                    "source modulo never selects an out-of-range lane"
                );
            }
        }

        for (count_index, source_rows) in LARGE_SOURCE_ROWS.into_iter().enumerate() {
            for lane_count_value in 1..=4_u32 {
                let lane_count = std::num::NonZeroU32::new(lane_count_value).unwrap();
                let expected = LARGE_EXPECTED[count_index][(lane_count_value - 1) as usize];
                let mut total = 0_u64;
                for thread_idx in 0..lane_count_value {
                    let actual =
                        super::wrapped_source_lane_cardinality(source_rows, lane_count, thread_idx);
                    assert_eq!(
                        actual, expected[thread_idx as usize],
                        "large rows={source_rows}, lanes={lane_count_value}, lane={thread_idx}"
                    );
                    total = total.checked_add(actual).expect("large row sum fits u64");
                }
                assert_eq!(total, source_rows);
            }
        }

        for lane_count_value in 1..=4_u32 {
            let lane_count = std::num::NonZeroU32::new(lane_count_value).unwrap();
            let expected = MAX_EXPECTED[(lane_count_value - 1) as usize];
            let mut total = 0_u64;
            for thread_idx in 0..lane_count_value {
                let actual =
                    super::wrapped_source_lane_cardinality(u64::MAX, lane_count, thread_idx);
                assert_eq!(
                    actual, expected[thread_idx as usize],
                    "maximum rows, lanes={lane_count_value}, lane={thread_idx}"
                );
                total = total.checked_add(actual).expect("maximum row sum fits u64");
            }
            assert_eq!(total, u64::MAX);
        }

        for (lane_count_value, expected_lanes) in WRAP_EXPECTED_LANES {
            let lane_count = std::num::NonZeroU32::new(lane_count_value).unwrap();
            let mut assignment = super::WorkerConformerAssignment {
                source_index: WRAP_INDICES[0],
                thread_idx: 0,
                lane_count,
            };
            for (wrap_offset, (source_index, expected_lane)) in
                WRAP_INDICES.into_iter().zip(expected_lanes).enumerate()
            {
                assert_eq!(assignment.source_index, source_index);
                for thread_idx in 0..lane_count_value {
                    let candidate = super::WorkerConformerAssignment {
                        source_index,
                        thread_idx,
                        lane_count,
                    };
                    assert_eq!(candidate.selects_current(), thread_idx == expected_lane);
                }
                assignment.advance();
                assert_eq!(
                    assignment.source_index,
                    source_index.wrapping_add(1),
                    "counter transition at wrap offset {wrap_offset}"
                );
            }
        }
    }

    #[test]
    fn uff_mt_m03_partition_preserves_source_rows_slots_and_preconditions() {
        const CONFORMER_IDS: [usize; 7] = [900, 2, 80, 17, 64, 7, 501];
        const WORKER_COUNTS: [i32; 6] = [-2, 0, 1, 2, 3, 4];
        const EXPECTED_SOURCE_INDICES: [[[&[usize]; 4]; 4]; 8] = [
            [
                [&[], &[], &[], &[]],
                [&[], &[], &[], &[]],
                [&[], &[], &[], &[]],
                [&[], &[], &[], &[]],
            ],
            [
                [&[0], &[], &[], &[]],
                [&[0], &[], &[], &[]],
                [&[0], &[], &[], &[]],
                [&[0], &[], &[], &[]],
            ],
            [
                [&[0, 1], &[], &[], &[]],
                [&[0], &[1], &[], &[]],
                [&[0], &[1], &[], &[]],
                [&[0], &[1], &[], &[]],
            ],
            [
                [&[0, 1, 2], &[], &[], &[]],
                [&[0, 2], &[1], &[], &[]],
                [&[0], &[1], &[2], &[]],
                [&[0], &[1], &[2], &[]],
            ],
            [
                [&[0, 1, 2, 3], &[], &[], &[]],
                [&[0, 2], &[1, 3], &[], &[]],
                [&[0, 3], &[1], &[2], &[]],
                [&[0], &[1], &[2], &[3]],
            ],
            [
                [&[0, 1, 2, 3, 4], &[], &[], &[]],
                [&[0, 2, 4], &[1, 3], &[], &[]],
                [&[0, 3], &[1, 4], &[2], &[]],
                [&[0, 4], &[1], &[2], &[3]],
            ],
            [
                [&[0, 1, 2, 3, 4, 5], &[], &[], &[]],
                [&[0, 2, 4], &[1, 3, 5], &[], &[]],
                [&[0, 3], &[1, 4], &[2, 5], &[]],
                [&[0, 4], &[1, 5], &[2], &[3]],
            ],
            [
                [&[0, 1, 2, 3, 4, 5, 6], &[], &[], &[]],
                [&[0, 2, 4, 6], &[1, 3, 5], &[], &[]],
                [&[0, 3, 6], &[1, 4], &[2, 5], &[]],
                [&[0, 4], &[1, 5], &[2, 6], &[3]],
            ],
        ];
        const EXPECTED_ROW_IDS: [[[&[usize]; 4]; 4]; 8] = [
            [
                [&[], &[], &[], &[]],
                [&[], &[], &[], &[]],
                [&[], &[], &[], &[]],
                [&[], &[], &[], &[]],
            ],
            [
                [&[900], &[], &[], &[]],
                [&[900], &[], &[], &[]],
                [&[900], &[], &[], &[]],
                [&[900], &[], &[], &[]],
            ],
            [
                [&[900, 2], &[], &[], &[]],
                [&[900], &[2], &[], &[]],
                [&[900], &[2], &[], &[]],
                [&[900], &[2], &[], &[]],
            ],
            [
                [&[900, 2, 80], &[], &[], &[]],
                [&[900, 80], &[2], &[], &[]],
                [&[900], &[2], &[80], &[]],
                [&[900], &[2], &[80], &[]],
            ],
            [
                [&[900, 2, 80, 17], &[], &[], &[]],
                [&[900, 80], &[2, 17], &[], &[]],
                [&[900, 17], &[2], &[80], &[]],
                [&[900], &[2], &[80], &[17]],
            ],
            [
                [&[900, 2, 80, 17, 64], &[], &[], &[]],
                [&[900, 80, 64], &[2, 17], &[], &[]],
                [&[900, 17], &[2, 64], &[80], &[]],
                [&[900, 64], &[2], &[80], &[17]],
            ],
            [
                [&[900, 2, 80, 17, 64, 7], &[], &[], &[]],
                [&[900, 80, 64], &[2, 17, 7], &[], &[]],
                [&[900, 17], &[2, 64], &[80, 7], &[]],
                [&[900, 64], &[2, 7], &[80], &[17]],
            ],
            [
                [&[900, 2, 80, 17, 64, 7, 501], &[], &[], &[]],
                [&[900, 80, 64, 501], &[2, 17, 7], &[], &[]],
                [&[900, 17, 501], &[2, 64], &[80, 7], &[]],
                [&[900, 64], &[2, 7], &[80, 501], &[17]],
            ],
        ];
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        for conformer_count in 0..=7_usize {
            for num_threads in WORKER_COUNTS {
                for extra_slots in [0, 2] {
                    let result_slots = conformer_count + extra_slots;
                    let mut coordinate_rows = (0..conformer_count)
                        .map(|row| [[row as f64 + 1.0, 2.0, 3.0]])
                        .collect::<Vec<_>>();
                    let before_coordinate_bits = coordinate_rows
                        .iter()
                        .flat_map(|row| row[0].iter().copied())
                        .map(f64::to_bits)
                        .collect::<Vec<_>>();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(row, positions)| super::SerialConformer {
                            id: CONFORMER_IDS[row],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; result_slots];

                    let mut lanes =
                        super::partition_worker_lanes(&mut conformers, &mut results, num_threads)
                            .expect("exact and extra capacity cases are accepted");
                    if num_threads <= 0 {
                        assert!(lanes.is_empty());
                    } else {
                        let lane_count = num_threads as usize;
                        assert_eq!(lanes.len(), lane_count);
                        let expected_indices =
                            EXPECTED_SOURCE_INDICES[conformer_count][lane_count - 1];
                        let expected_ids = EXPECTED_ROW_IDS[conformer_count][lane_count - 1];
                        let mut assigned_counts = [0_usize; 7];
                        for (lane_index, lane) in lanes.iter_mut().enumerate() {
                            let actual_indices = lane
                                .rows
                                .iter()
                                .map(|row| row.source_index as usize)
                                .collect::<Vec<_>>();
                            let actual_ids = lane
                                .rows
                                .iter()
                                .map(|row| row.conformer_id)
                                .collect::<Vec<_>>();
                            assert_eq!(actual_indices.as_slice(), expected_indices[lane_index]);
                            assert_eq!(actual_ids.as_slice(), expected_ids[lane_index]);
                            assert_eq!(
                                lane.source_result_slots.len(),
                                expected_indices[lane_index].len(),
                                "unique active slots follow their source slot modulo"
                            );

                            for (source_slot, result) in expected_indices[lane_index]
                                .iter()
                                .copied()
                                .zip(lane.source_result_slots.iter_mut())
                            {
                                result.status = 100 + lane_index as i32;
                                result.energy = source_slot as f64;
                            }
                            for row in &mut lane.rows {
                                assert_eq!(
                                    row.conformer_id,
                                    CONFORMER_IDS[row.source_index as usize]
                                );
                                row.conformer.positions[0][0] = 100.0 + row.source_index as f64;
                                assigned_counts[row.source_index as usize] += 1;
                            }
                        }
                        assert!(
                            assigned_counts[..conformer_count]
                                .iter()
                                .all(|count| *count == 1)
                        );
                    }

                    drop(lanes);
                    drop(conformers);
                    for (row, coordinates) in coordinate_rows.iter().enumerate() {
                        let expected_x = if num_threads > 0 {
                            100.0 + row as f64
                        } else {
                            row as f64 + 1.0
                        };
                        assert_eq!(coordinates[0], [expected_x, 2.0, 3.0]);
                    }
                    for source_slot in 0..conformer_count {
                        if num_threads > 0 {
                            assert_eq!(
                                results[source_slot].status,
                                100 + (source_slot % num_threads as usize) as i32
                            );
                            assert_eq!(results[source_slot].energy, source_slot as f64);
                        } else {
                            assert_eq!(results[source_slot], sentinel);
                        }
                    }
                    for result in results.iter().skip(conformer_count) {
                        assert_eq!(*result, sentinel, "extra result suffix stays untouched");
                    }
                    if num_threads <= 0 {
                        let after_coordinate_bits = coordinate_rows
                            .iter()
                            .flat_map(|row| row[0].iter().copied())
                            .map(f64::to_bits)
                            .collect::<Vec<_>>();
                        assert_eq!(after_coordinate_bits, before_coordinate_bits);
                    }
                }
            }
        }

        for conformer_count in 1..=7_usize {
            for num_threads in [-2, 0] {
                let mut coordinate_rows = (0..conformer_count)
                    .map(|row| [[row as f64 + 1.0, 2.0, 3.0]])
                    .collect::<Vec<_>>();
                let before_coordinate_bits = coordinate_rows
                    .iter()
                    .flat_map(|row| row[0].iter().copied())
                    .map(f64::to_bits)
                    .collect::<Vec<_>>();
                let mut conformers = coordinate_rows
                    .iter_mut()
                    .enumerate()
                    .map(|(row, positions)| super::SerialConformer {
                        id: CONFORMER_IDS[row],
                        positions: positions.as_mut_slice(),
                    })
                    .collect::<Vec<_>>();
                let mut results = vec![sentinel; conformer_count - 1];
                let results_before = results.clone();
                let lanes =
                    super::partition_worker_lanes(&mut conformers, &mut results, num_threads)
                        .expect("nonpositive thread counts skip the source precondition");
                assert!(lanes.is_empty());
                drop(lanes);
                drop(conformers);
                assert_eq!(results, results_before);
                let after_coordinate_bits = coordinate_rows
                    .iter()
                    .flat_map(|row| row[0].iter().copied())
                    .map(f64::to_bits)
                    .collect::<Vec<_>>();
                assert_eq!(after_coordinate_bits, before_coordinate_bits);
            }
        }

        for conformer_count in 1..=7_usize {
            for num_threads in 1..=4_i32 {
                let mut coordinate_rows = (0..conformer_count)
                    .map(|row| [[row as f64 + 1.0, 2.0, 3.0]])
                    .collect::<Vec<_>>();
                let before_coordinate_bits = coordinate_rows
                    .iter()
                    .flat_map(|row| row[0].iter().copied())
                    .map(f64::to_bits)
                    .collect::<Vec<_>>();
                let mut conformers = coordinate_rows
                    .iter_mut()
                    .enumerate()
                    .map(|(row, positions)| super::SerialConformer {
                        id: CONFORMER_IDS[row],
                        positions: positions.as_mut_slice(),
                    })
                    .collect::<Vec<_>>();
                let mut results = vec![sentinel; conformer_count - 1];
                let results_before = results.clone();
                assert!(matches!(
                    super::partition_worker_lanes(
                        &mut conformers,
                        &mut results,
                        num_threads,
                    ),
                    Err(super::SerialConformerOptimizationError::ResultCapacity(
                        super::SerialResultCapacityError {
                            conformers: actual_conformers,
                            result_slots,
                        }
                    )) if actual_conformers == conformer_count
                        && result_slots == conformer_count - 1
                ));
                drop(conformers);
                assert_eq!(results, results_before);
                let after_coordinate_bits = coordinate_rows
                    .iter()
                    .flat_map(|row| row[0].iter().copied())
                    .map(f64::to_bits)
                    .collect::<Vec<_>>();
                assert_eq!(after_coordinate_bits, before_coordinate_bits);
            }
        }
    }

    #[test]
    fn uff_mt_m04_dense_and_strided_views_share_the_w05_stage_loop() {
        const ID_LAYOUTS: [[usize; 3]; 2] = [[900, 2, 80], [80, 7, 900]];
        const W05_ENERGIES: [f64; 3] = [4.0, 9.0, 16.0];
        const W05_COORDINATES: [[[f64; 3]; 2]; 3] = [
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
        ];
        const W05_SOURCE_X: [f64; 3] = [2.0, 3.0, 4.0];
        let one_step_x = W05_SOURCE_X.map(|x| x + (-(2.0 * x * 0.1)));
        let one_step_energies = one_step_x.map(|x| x * x);
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let norm_reference = |coordinates: &[[f64; 3]]| {
            coordinates
                .iter()
                .flatten()
                .map(|coordinate| coordinate * coordinate)
                .sum::<f64>()
        };

        for ids in ID_LAYOUTS {
            for row_count in 0..=3_usize {
                for max_iterations in [0_i32, 1] {
                    for worker_count in 1..=3_i32 {
                        let mut dense_coordinates = W05_COORDINATES[..row_count]
                            .iter()
                            .copied()
                            .map(|row| row.to_vec())
                            .collect::<Vec<_>>();
                        let dense_before = dense_coordinates.clone();
                        let mut dense_conformers = dense_coordinates
                            .iter_mut()
                            .enumerate()
                            .map(|(row, positions)| super::SerialConformer {
                                id: ids[row],
                                positions: positions.as_mut_slice(),
                            })
                            .collect::<Vec<_>>();
                        let mut dense_results = vec![sentinel; row_count + 1];
                        let dense_result_slots = dense_results.len();
                        let dense_rows = dense_conformers.iter_mut().enumerate().map(
                            |(source_index, conformer)| super::SourceIndexedWorkerRow {
                                source_index: source_index as u32,
                                conformer_id: conformer.id,
                                conformer,
                            },
                        );
                        let mut dense_field = ForceField::new(3);
                        dense_field.add_contribution(Box::new(CoordinateNormContribution));
                        assert_eq!(
                            super::optimize_assigned_worker_rows(
                                dense_field,
                                dense_rows,
                                super::WorkerResultSlotView::Dense(&mut dense_results),
                                row_count,
                                dense_result_slots,
                                2,
                                max_iterations,
                            ),
                            Ok(())
                        );
                        drop(dense_conformers);

                        for source_index in 0..row_count {
                            let expected_coordinates = if max_iterations == 0 {
                                W05_COORDINATES[source_index]
                            } else {
                                [[0.0, 0.0, 0.0], [one_step_x[source_index], 0.0, 0.0]]
                            };
                            let expected_energy = if max_iterations == 0 {
                                W05_ENERGIES[source_index]
                            } else {
                                one_step_energies[source_index]
                            };
                            assert_eq!(
                                dense_results[source_index].status, 1,
                                "dense source slot {source_index}, iterations={max_iterations}, workers={worker_count}"
                            );
                            for point_index in 0..2 {
                                for axis in 0..3 {
                                    assert_eq!(
                                        dense_coordinates[source_index][point_index][axis]
                                            .to_bits(),
                                        expected_coordinates[point_index][axis].to_bits(),
                                        "dense fixed W05 coordinate: row={source_index}, point={point_index}, axis={axis}, iterations={max_iterations}, workers={worker_count}"
                                    );
                                }
                            }
                            assert_eq!(
                                dense_results[source_index].energy.to_bits(),
                                expected_energy.to_bits(),
                                "dense fixed source W05 energy: row={source_index}, iterations={max_iterations}, workers={worker_count}"
                            );
                            assert_eq!(
                                norm_reference(&dense_coordinates[source_index]).to_bits(),
                                expected_energy.to_bits(),
                                "dense independent coordinate-norm consistency"
                            );
                            if max_iterations == 0 {
                                assert_eq!(
                                    dense_coordinates[source_index],
                                    dense_before[source_index]
                                );
                            } else {
                                assert_ne!(
                                    dense_coordinates[source_index],
                                    dense_before[source_index]
                                );
                                assert!(
                                    dense_results[source_index].energy < W05_ENERGIES[source_index]
                                );
                            }
                        }
                        assert_eq!(dense_results[row_count], sentinel);

                        let mut strided_coordinates = W05_COORDINATES[..row_count]
                            .iter()
                            .copied()
                            .map(|row| row.to_vec())
                            .collect::<Vec<_>>();
                        let strided_before = strided_coordinates.clone();
                        let mut strided_conformers = strided_coordinates
                            .iter_mut()
                            .enumerate()
                            .map(|(row, positions)| super::SerialConformer {
                                id: ids[row],
                                positions: positions.as_mut_slice(),
                            })
                            .collect::<Vec<_>>();
                        let mut strided_results = vec![sentinel; row_count + 1];
                        let strided_result_slots = strided_results.len();
                        let lane_count = std::num::NonZeroU32::new(worker_count as u32).unwrap();
                        let lanes = super::partition_worker_lanes(
                            &mut strided_conformers,
                            &mut strided_results,
                            worker_count,
                        )
                        .expect("the source capacity precondition accepts the extra slot");
                        assert_eq!(lanes.len(), worker_count as usize);
                        for (lane_index, lane) in lanes.into_iter().enumerate() {
                            let expected_indices = (0..row_count)
                                .filter(|source_index| {
                                    source_index % worker_count as usize == lane_index
                                })
                                .collect::<Vec<_>>();
                            let actual_indices = lane
                                .rows
                                .iter()
                                .map(|row| row.source_index as usize)
                                .collect::<Vec<_>>();
                            assert_eq!(actual_indices, expected_indices);
                            for row in &lane.rows {
                                assert_eq!(row.conformer_id, ids[row.source_index as usize]);
                            }

                            let super::PartitionedWorkerLane {
                                rows,
                                source_result_slots,
                            } = lane;
                            let mut lane_field = ForceField::new(3);
                            lane_field.add_contribution(Box::new(CoordinateNormContribution));
                            assert_eq!(
                                super::optimize_assigned_worker_rows(
                                    lane_field,
                                    rows.into_iter(),
                                    super::WorkerResultSlotView::Strided {
                                        source_result_slots,
                                        lane_count,
                                    },
                                    row_count,
                                    strided_result_slots,
                                    2,
                                    max_iterations,
                                ),
                                Ok(())
                            );
                        }
                        drop(strided_conformers);

                        for source_index in 0..row_count {
                            let expected_coordinates = if max_iterations == 0 {
                                W05_COORDINATES[source_index]
                            } else {
                                [[0.0, 0.0, 0.0], [one_step_x[source_index], 0.0, 0.0]]
                            };
                            let expected_energy = if max_iterations == 0 {
                                W05_ENERGIES[source_index]
                            } else {
                                one_step_energies[source_index]
                            };
                            assert_eq!(
                                strided_results[source_index].status, 1,
                                "strided source slot {source_index}, lanes={worker_count}, iterations={max_iterations}"
                            );
                            for point_index in 0..2 {
                                for axis in 0..3 {
                                    assert_eq!(
                                        strided_coordinates[source_index][point_index][axis]
                                            .to_bits(),
                                        expected_coordinates[point_index][axis].to_bits(),
                                        "strided fixed W05 coordinate: row={source_index}, point={point_index}, axis={axis}, iterations={max_iterations}, workers={worker_count}"
                                    );
                                }
                            }
                            assert_eq!(
                                strided_results[source_index].energy.to_bits(),
                                expected_energy.to_bits(),
                                "strided fixed source W05 energy: row={source_index}, iterations={max_iterations}, workers={worker_count}"
                            );
                            assert_eq!(
                                norm_reference(&strided_coordinates[source_index]).to_bits(),
                                expected_energy.to_bits(),
                                "strided independent coordinate-norm consistency"
                            );
                            if max_iterations == 0 {
                                assert_eq!(
                                    strided_coordinates[source_index],
                                    strided_before[source_index]
                                );
                            } else {
                                assert_ne!(
                                    strided_coordinates[source_index],
                                    strided_before[source_index]
                                );
                                assert!(
                                    strided_results[source_index].energy
                                        < W05_ENERGIES[source_index]
                                );
                            }
                        }
                        assert_eq!(strided_results[row_count], sentinel);
                    }
                }
            }
        }

        // Synthetic wrap replay: the source u32 index returns to zero after
        // 2^32 rows. Repeated zero indices in one lane must overwrite its one
        // borrowed source result slot sequentially without duplicate borrows.
        let mut repeated_first_coordinates = [[1.0, 0.0, 0.0]];
        let mut repeated_second_coordinates = [[2.0, 0.0, 0.0]];
        let mut repeated_conformers = [
            super::SerialConformer {
                id: 10,
                positions: &mut repeated_first_coordinates,
            },
            super::SerialConformer {
                id: 20,
                positions: &mut repeated_second_coordinates,
            },
        ];
        let mut repeated_result = sentinel;
        {
            let (first, second) = repeated_conformers.split_at_mut(1);
            let repeated_rows = [
                super::SourceIndexedWorkerRow {
                    source_index: 0,
                    conformer_id: first[0].id,
                    conformer: &mut first[0],
                },
                super::SourceIndexedWorkerRow {
                    source_index: 0,
                    conformer_id: second[0].id,
                    conformer: &mut second[0],
                },
            ];
            let mut repeated_field = ForceField::new(3);
            repeated_field.add_contribution(Box::new(CoordinateNormContribution));
            assert_eq!(
                super::optimize_assigned_worker_rows(
                    repeated_field,
                    repeated_rows.into_iter(),
                    super::WorkerResultSlotView::Strided {
                        source_result_slots: vec![&mut repeated_result],
                        lane_count: std::num::NonZeroU32::new(2).unwrap(),
                    },
                    2,
                    2,
                    1,
                    0,
                ),
                Ok(())
            );
        }
        drop(repeated_conformers);
        assert_eq!(
            repeated_result,
            super::OptimizationOutcome {
                status: 1,
                energy: 4.0,
            }
        );
    }

    #[test]
    fn uff_mt_m05_dense_worker_matches_fixed_w05_effects_for_lane_matrix() {
        const SOURCE_X: [f64; 3] = [2.0, 3.0, 4.0];
        const IDS: [usize; 3] = [900, 2, 80];
        let one_step_x = SOURCE_X.map(|x| x + (-(2.0 * x * 0.1)));
        let one_step_energies = one_step_x.map(|x| x * x);
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        for row_count in 0..=SOURCE_X.len() {
            for worker_count in 1_u32..=3 {
                for lane in (0..=worker_count).chain([worker_count + 1, u32::MAX]) {
                    for max_iterations in [0_i32, 1] {
                        let mut coordinate_rows = SOURCE_X[..row_count]
                            .iter()
                            .copied()
                            .map(|x| vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]])
                            .collect::<Vec<_>>();
                        let before = coordinate_rows.clone();
                        let mut conformers = coordinate_rows
                            .iter_mut()
                            .enumerate()
                            .map(|(input_index, positions)| super::SerialConformer {
                                id: IDS[input_index],
                                positions: positions.as_mut_slice(),
                            })
                            .collect::<Vec<_>>();
                        let mut results = vec![sentinel; row_count + 1];
                        cf3d_uff_one_kernel_counts_reset();
                        let mut field = ForceField::new(3);
                        field.add_contribution(Box::new(CoordinateNormContribution));

                        assert_eq!(
                            super::optimize_worker_conformers(
                                field,
                                &mut conformers,
                                &mut results,
                                2,
                                lane,
                                std::num::NonZeroU32::new(worker_count).unwrap(),
                                max_iterations,
                            ),
                            Ok(()),
                            "rows={row_count}, workers={worker_count}, lane={lane}, iterations={max_iterations}"
                        );
                        drop(conformers);

                        let mut selected_count = 0;
                        for input_index in 0..row_count {
                            let selected = input_index as u32 % worker_count == lane;
                            let expected_coordinates = if selected && max_iterations == 1 {
                                selected_count += 1;
                                [[0.0, 0.0, 0.0], [one_step_x[input_index], 0.0, 0.0]]
                            } else {
                                if selected {
                                    selected_count += 1;
                                }
                                before[input_index].clone().try_into().unwrap()
                            };
                            for point_index in 0..2 {
                                for axis in 0..3 {
                                    assert_eq!(
                                        coordinate_rows[input_index][point_index][axis].to_bits(),
                                        expected_coordinates[point_index][axis].to_bits(),
                                        "W05 coordinate row={input_index}, point={point_index}, axis={axis}, rows={row_count}, workers={worker_count}, lane={lane}, iterations={max_iterations}"
                                    );
                                }
                            }

                            if selected {
                                let expected_energy = if max_iterations == 0 {
                                    SOURCE_X[input_index] * SOURCE_X[input_index]
                                } else {
                                    one_step_energies[input_index]
                                };
                                assert_eq!(
                                    results[input_index].status, 1,
                                    "W05 source status row={input_index}, workers={worker_count}, lane={lane}, iterations={max_iterations}"
                                );
                                assert_eq!(
                                    results[input_index].energy.to_bits(),
                                    expected_energy.to_bits(),
                                    "W05 fixed source energy row={input_index}, workers={worker_count}, lane={lane}, iterations={max_iterations}"
                                );
                            } else {
                                assert_eq!(results[input_index], sentinel);
                            }
                        }
                        assert_eq!(results[row_count], sentinel);
                        assert_eq!(
                            cf3d_uff_one_serial_work_counts().initialize_calls,
                            selected_count,
                            "rows={row_count}, workers={worker_count}, lane={lane}, iterations={max_iterations}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn uff_mt_m05_dense_worker_matches_fixed_w06_vdw_effects_for_lane_matrix() {
        use crate::kernel::Cf3dFragAcceptContributionIdentity as Identity;

        const TARGET_DISTANCES: [f64; 3] = [2.0, 3.0, 4.0];
        const IDS: [usize; 3] = [80, 7, 900];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let threshold = 10.0 * minimum;
        let source_energy = |distance: f64| {
            if distance <= 0.0 || distance > threshold {
                return 0.0;
            }
            let ratio = minimum / distance;
            let ratio3 = ratio * ratio * ratio;
            let ratio6 = ratio3 * ratio3;
            well_depth * (ratio6 * ratio6 - 2.0 * ratio6)
        };
        let reference_points = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut reference = worker_w06_coordinates(&reference_points);
        let source_field =
            worker_w06_source_field(&topology, &mut reference, &total_valences, 10.0, false);
        let expected_identity = vec![Identity::Vdw {
            at1_idx: 0,
            at2_idx: 1,
            x_ij: minimum,
            well_depth,
            threshold,
        }];
        assert_eq!(
            crate::kernel::cf3d_frag_accept_contribution_identities(&source_field),
            expected_identity
        );
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        for row_count in 0..=TARGET_DISTANCES.len() {
            for worker_count in 1_u32..=3 {
                for lane in (0..=worker_count).chain([worker_count + 1, u32::MAX]) {
                    let mut coordinate_rows = TARGET_DISTANCES[..row_count]
                        .iter()
                        .copied()
                        .map(|distance| vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
                        .collect::<Vec<_>>();
                    let before = coordinate_rows.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(input_index, positions)| super::SerialConformer {
                            id: IDS[input_index],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; row_count + 1];
                    let field = source_field.copy();
                    assert_eq!(
                        crate::kernel::cf3d_frag_accept_contribution_identities(&field),
                        expected_identity
                    );
                    cf3d_uff_one_kernel_counts_reset();
                    assert_eq!(
                        super::optimize_worker_conformers(
                            field,
                            &mut conformers,
                            &mut results,
                            2,
                            lane,
                            std::num::NonZeroU32::new(worker_count).unwrap(),
                            0,
                        ),
                        Ok(()),
                        "W06 rows={row_count}, workers={worker_count}, lane={lane}"
                    );
                    drop(conformers);

                    let mut selected_count = 0;
                    for input_index in 0..row_count {
                        let selected = input_index as u32 % worker_count == lane;
                        for point_index in 0..2 {
                            for axis in 0..3 {
                                assert_eq!(
                                    coordinate_rows[input_index][point_index][axis].to_bits(),
                                    before[input_index][point_index][axis].to_bits(),
                                    "W06 zero-iteration coordinate row={input_index}, point={point_index}, axis={axis}, rows={row_count}, workers={worker_count}, lane={lane}"
                                );
                            }
                        }
                        if selected {
                            selected_count += 1;
                            assert_eq!(results[input_index].status, 1);
                            assert!(
                                (results[input_index].energy
                                    - source_energy(TARGET_DISTANCES[input_index]))
                                .abs()
                                    <= 1.0e-12,
                                "W06 source VDW energy row={input_index}, workers={worker_count}, lane={lane}"
                            );
                        } else {
                            assert_eq!(results[input_index], sentinel);
                        }
                    }
                    assert_eq!(results[row_count], sentinel);
                    assert_eq!(
                        cf3d_uff_one_serial_work_counts().initialize_calls,
                        selected_count,
                        "W06 rows={row_count}, workers={worker_count}, lane={lane}"
                    );
                }
            }
        }
        assert_eq!(reference.conformers_3d[0].coordinates(), reference_points);
    }

    #[test]
    fn uff_mt_m05_dense_worker_preserves_fixed_stage_failure_effects_for_lane_matrix() {
        const IDS: [usize; 6] = [80, 7, 900, 31, 502, 8];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        for worker_count in 1_u32..=3 {
            let row_count = 2 * worker_count as usize;
            for lane in (0..=worker_count).chain([worker_count + 1, u32::MAX]) {
                for stage in [0_u8, 1, 2] {
                    let mut coordinate_rows = (0..row_count)
                        .map(|_| {
                            if stage == 0 {
                                Vec::<[f64; 3]>::new()
                            } else {
                                vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
                            }
                        })
                        .collect::<Vec<_>>();
                    let before = coordinate_rows.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(input_index, positions)| super::SerialConformer {
                            id: IDS[input_index],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; row_count + 1];
                    cf3d_uff_one_kernel_counts_reset();
                    let mut field = ForceField::new(3);
                    if stage != 0 {
                        field.add_contribution(Box::new(FailAtEnergyCallWithCoordinateNorm {
                            calls: Cell::new(0),
                            fail_on: if stage == 1 { 3 } else { 4 },
                        }));
                    }
                    let atom_count = if stage == 0 { 0 } else { 2 };
                    let result = super::optimize_worker_conformers(
                        field,
                        &mut conformers,
                        &mut results,
                        atom_count,
                        lane,
                        std::num::NonZeroU32::new(worker_count).unwrap(),
                        0,
                    );
                    drop(conformers);

                    if lane >= worker_count {
                        assert_eq!(result, Ok(()));
                        assert_eq!(results, vec![sentinel; row_count + 1]);
                        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 0);
                    } else if stage == 0 {
                        assert_eq!(
                            result,
                            Err(super::SerialConformerOptimizationError::Optimization {
                                input_index: lane as usize,
                                conformer_id: IDS[lane as usize],
                                source: super::OptimizationStageError::Initialize(
                                    ForceFieldKernelError::NoPoints
                                ),
                            })
                        );
                        assert_eq!(results, vec![sentinel; row_count + 1]);
                        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 1);
                    } else {
                        let failed_input_index = lane as usize + worker_count as usize;
                        let source = if stage == 1 {
                            super::OptimizationStageError::Minimize(expected_pair_index_error())
                        } else {
                            super::OptimizationStageError::FinalEnergy(expected_pair_index_error())
                        };
                        assert_eq!(
                            result,
                            Err(super::SerialConformerOptimizationError::Optimization {
                                input_index: failed_input_index,
                                conformer_id: IDS[failed_input_index],
                                source,
                            })
                        );
                        let mut expected_results = vec![sentinel; row_count + 1];
                        expected_results[lane as usize] = OptimizationOutcome {
                            status: 1,
                            energy: 4.0,
                        };
                        assert_eq!(results, expected_results);
                        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 2);
                    }

                    for (row_index, row) in coordinate_rows.iter().enumerate() {
                        assert_eq!(row.len(), before[row_index].len());
                        for (actual_point, expected_point) in row.iter().zip(&before[row_index]) {
                            for axis in 0..3 {
                                assert_eq!(
                                    actual_point[axis].to_bits(),
                                    expected_point[axis].to_bits(),
                                    "stage={stage}, rows={row_count}, workers={worker_count}, lane={lane}, row={row_index}, axis={axis}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn uff_mt_m06_partitioned_worker_matches_fixed_w05_reference_matrix() {
        const SOURCE_X: [f64; 3] = [2.0, 3.0, 4.0];
        const IDS: [usize; 3] = [900, 2, 80];
        let one_step_x = SOURCE_X.map(|x| x + (-(2.0 * x * 0.1)));
        let one_step_energies = one_step_x.map(|x| x * x);
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut source_field = ForceField::new(3);
        source_field.add_contribution(Box::new(CoordinateNormContribution));

        for row_count in 0..=SOURCE_X.len() {
            for worker_count in 1_u32..=3 {
                for max_iterations in [0_i32, 1] {
                    let mut coordinate_rows = SOURCE_X[..row_count]
                        .iter()
                        .copied()
                        .map(|x| vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]])
                        .collect::<Vec<_>>();
                    let before = coordinate_rows.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(input_index, positions)| super::SerialConformer {
                            id: IDS[input_index],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; row_count + 1];
                    let result_slot_count = results.len();
                    let lane_count = std::num::NonZeroU32::new(worker_count).unwrap();
                    let lanes = super::partition_worker_lanes(
                        &mut conformers,
                        &mut results,
                        worker_count as i32,
                    )
                    .expect("the original result capacity includes the extra slot");
                    assert_eq!(lanes.len(), worker_count as usize);
                    cf3d_uff_one_kernel_counts_reset();

                    for (lane_index, lane) in lanes.into_iter().enumerate() {
                        let expected_indices = (0..row_count)
                            .filter(|source_index| {
                                source_index % worker_count as usize == lane_index
                            })
                            .collect::<Vec<_>>();
                        let actual_indices = lane
                            .rows
                            .iter()
                            .map(|row| row.source_index as usize)
                            .collect::<Vec<_>>();
                        assert_eq!(actual_indices, expected_indices);
                        for row in &lane.rows {
                            assert_eq!(row.conformer_id, IDS[row.source_index as usize]);
                        }
                        assert_eq!(lane.source_result_slots.len(), expected_indices.len());

                        super::optimize_partitioned_worker(
                            source_field.copy(),
                            lane,
                            row_count,
                            result_slot_count,
                            2,
                            lane_count,
                            max_iterations,
                        )
                        .expect("each copied W05 lane processes its ordered source rows");
                    }
                    drop(conformers);

                    for source_index in 0..row_count {
                        let expected_coordinates = if max_iterations == 0 {
                            before[source_index].clone().try_into().unwrap()
                        } else {
                            [[0.0, 0.0, 0.0], [one_step_x[source_index], 0.0, 0.0]]
                        };
                        for point_index in 0..2 {
                            for axis in 0..3 {
                                assert_eq!(
                                    coordinate_rows[source_index][point_index][axis].to_bits(),
                                    expected_coordinates[point_index][axis].to_bits(),
                                    "W05 coordinate row={source_index}, point={point_index}, axis={axis}, rows={row_count}, workers={worker_count}, iterations={max_iterations}"
                                );
                            }
                        }
                        let expected_energy = if max_iterations == 0 {
                            SOURCE_X[source_index] * SOURCE_X[source_index]
                        } else {
                            one_step_energies[source_index]
                        };
                        assert_eq!(results[source_index].status, 1);
                        assert_eq!(
                            results[source_index].energy.to_bits(),
                            expected_energy.to_bits(),
                            "W05 fixed source result slot={source_index}, workers={worker_count}, iterations={max_iterations}"
                        );
                    }
                    assert_eq!(results[row_count], sentinel);
                    assert!(source_field.positions().is_empty());
                    assert_eq!(
                        cf3d_uff_one_serial_work_counts().initialize_calls,
                        row_count
                    );
                }
            }
        }
    }

    #[test]
    fn uff_mt_m06_partitioned_worker_matches_fixed_w06_vdw_reference_matrix() {
        use crate::kernel::Cf3dFragAcceptContributionIdentity as Identity;

        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const IDS: [usize; 3] = [80, 7, 900];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let threshold = 10.0 * minimum;
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = pre_factor * (0.0 - distance) / distance;
            let second_gradient = pre_factor * (distance - 0.0) / distance;
            let first_gradient = first_gradient * 0.1;
            let second_gradient = second_gradient * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let reference_points = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut reference = worker_w06_coordinates(&reference_points);
        let source_field =
            worker_w06_source_field(&topology, &mut reference, &total_valences, 10.0, false);
        let expected_identity = vec![Identity::Vdw {
            at1_idx: 0,
            at2_idx: 1,
            x_ij: minimum,
            well_depth,
            threshold,
        }];
        assert_eq!(
            crate::kernel::cf3d_frag_accept_contribution_identities(&source_field),
            expected_identity
        );
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        for row_count in 0..=TARGET_DISTANCES.len() {
            for worker_count in 1_u32..=3 {
                for max_iterations in [0_i32, 1] {
                    let mut coordinate_rows = TARGET_DISTANCES[..row_count]
                        .iter()
                        .copied()
                        .map(|distance| vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
                        .collect::<Vec<_>>();
                    let before = coordinate_rows.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(input_index, positions)| super::SerialConformer {
                            id: IDS[input_index],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; row_count + 1];
                    let result_slot_count = results.len();
                    let lane_count = std::num::NonZeroU32::new(worker_count).unwrap();
                    let lanes = super::partition_worker_lanes(
                        &mut conformers,
                        &mut results,
                        worker_count as i32,
                    )
                    .expect("the original result capacity includes the extra slot");
                    assert_eq!(lanes.len(), worker_count as usize);
                    cf3d_uff_one_kernel_counts_reset();

                    for (lane_index, lane) in lanes.into_iter().enumerate() {
                        let expected_indices = (0..row_count)
                            .filter(|source_index| {
                                source_index % worker_count as usize == lane_index
                            })
                            .collect::<Vec<_>>();
                        let actual_indices = lane
                            .rows
                            .iter()
                            .map(|row| row.source_index as usize)
                            .collect::<Vec<_>>();
                        assert_eq!(actual_indices, expected_indices);
                        for row in &lane.rows {
                            assert_eq!(row.conformer_id, IDS[row.source_index as usize]);
                        }
                        assert_eq!(lane.source_result_slots.len(), expected_indices.len());

                        let lane_field = source_field.copy();
                        assert_eq!(
                            crate::kernel::cf3d_frag_accept_contribution_identities(&lane_field),
                            expected_identity
                        );
                        super::optimize_partitioned_worker(
                            lane_field,
                            lane,
                            row_count,
                            result_slot_count,
                            2,
                            lane_count,
                            max_iterations,
                        )
                        .expect("each copied W06 lane processes its ordered source rows");
                    }
                    drop(conformers);

                    for source_index in 0..row_count {
                        let expected_coordinates = if max_iterations == 0 {
                            before[source_index].clone().try_into().unwrap()
                        } else {
                            let [first_x, second_x] =
                                source_one_step(TARGET_DISTANCES[source_index]).0;
                            [[first_x, 0.0, 0.0], [second_x, 0.0, 0.0]]
                        };
                        for point_index in 0..2 {
                            for axis in 0..3 {
                                assert_eq!(
                                    coordinate_rows[source_index][point_index][axis].to_bits(),
                                    expected_coordinates[point_index][axis].to_bits(),
                                    "W06 coordinate row={source_index}, point={point_index}, axis={axis}, rows={row_count}, workers={worker_count}, iterations={max_iterations}"
                                );
                            }
                        }
                        let expected_energy = if max_iterations == 0 {
                            source_energy(TARGET_DISTANCES[source_index])
                        } else {
                            source_one_step(TARGET_DISTANCES[source_index]).1
                        };
                        assert_eq!(results[source_index].status, 1);
                        assert!(
                            (results[source_index].energy - expected_energy).abs() <= 1.0e-12,
                            "W06 fixed VDW result slot={source_index}, workers={worker_count}, iterations={max_iterations}"
                        );
                    }
                    assert_eq!(results[row_count], sentinel);
                    assert_eq!(
                        crate::kernel::cf3d_frag_accept_contribution_identities(&source_field),
                        expected_identity
                    );
                    assert_eq!(
                        cf3d_uff_one_serial_work_counts().initialize_calls,
                        row_count
                    );
                }
            }
        }
        drop(source_field);
        assert_eq!(reference.conformers_3d[0].coordinates(), reference_points);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_mt_m08_real_transport_matches_fixed_w05_matrix() {
        const SOURCE_X: [f64; 3] = [2.0, 3.0, 4.0];
        const IDS: [usize; 3] = [900, 2, 80];
        const WORKER_COUNTS: [i32; 6] = [-2, 0, 1, 2, 3, 4];
        let one_step_x = SOURCE_X.map(|x| x + (-(2.0 * x * 0.1)));
        let one_step_energies = one_step_x.map(|x| x * x);
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut source_field = ForceField::new(3);
        source_field.add_contribution(Box::new(CoordinateNormContribution));
        let mut owner_calls = 0;

        for row_count in 0..=SOURCE_X.len() {
            for worker_count in WORKER_COUNTS {
                for max_iterations in [0_i32, 1] {
                    let mut coordinate_rows = SOURCE_X[..row_count]
                        .iter()
                        .copied()
                        .map(|x| vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]])
                        .collect::<Vec<_>>();
                    let before = coordinate_rows.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(input_index, positions)| super::SerialConformer {
                            id: IDS[input_index],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; row_count + 1];
                    cf3d_uff_one_kernel_counts_reset();
                    let outcomes = super::optimize_conformers_mt(
                        &source_field,
                        &mut conformers,
                        &mut results,
                        worker_count,
                        2,
                        max_iterations,
                    )
                    .expect("the fixed W05 result capacity is sufficient");
                    owner_calls += 1;

                    if worker_count <= 0 {
                        assert!(outcomes.is_empty());
                    } else {
                        assert_eq!(outcomes.len(), worker_count as usize);
                        for (lane, outcome) in outcomes.into_iter().enumerate() {
                            assert!(
                                matches!(outcome, Ok(Ok(()))),
                                "W05 lane={lane}, rows={row_count}, workers={worker_count}, iterations={max_iterations}"
                            );
                        }
                    }
                    assert_eq!(
                        conformers
                            .iter()
                            .map(|conformer| conformer.id)
                            .collect::<Vec<_>>(),
                        IDS[..row_count]
                    );
                    drop(conformers);

                    for source_index in 0..row_count {
                        let expected_coordinates = if worker_count > 0 && max_iterations == 1 {
                            [[0.0, 0.0, 0.0], [one_step_x[source_index], 0.0, 0.0]]
                        } else {
                            before[source_index].clone().try_into().unwrap()
                        };
                        for point_index in 0..2 {
                            for axis in 0..3 {
                                assert_eq!(
                                    coordinate_rows[source_index][point_index][axis].to_bits(),
                                    expected_coordinates[point_index][axis].to_bits(),
                                    "W05 coordinate row={source_index}, point={point_index}, axis={axis}, rows={row_count}, workers={worker_count}, iterations={max_iterations}"
                                );
                            }
                        }
                        if worker_count > 0 {
                            let expected_energy = if max_iterations == 0 {
                                SOURCE_X[source_index] * SOURCE_X[source_index]
                            } else {
                                one_step_energies[source_index]
                            };
                            assert_eq!(
                                results[source_index].status, 1,
                                "W05 source-index status row={source_index}, workers={worker_count}, iterations={max_iterations}"
                            );
                            assert_eq!(
                                results[source_index].energy.to_bits(),
                                expected_energy.to_bits(),
                                "W05 source-index energy row={source_index}, workers={worker_count}, iterations={max_iterations}"
                            );
                        } else {
                            assert_eq!(results[source_index], sentinel);
                        }
                    }
                    assert_eq!(results[row_count], sentinel);
                    assert!(source_field.positions().is_empty());
                    let spawned_lanes = worker_count.max(0) as usize;
                    assert_eq!(
                        cf3d_uff_one_kernel_counts(),
                        (spawned_lanes, 0, spawned_lanes),
                        "one source field and one term copy per spawned lane"
                    );
                }
            }
        }
        assert_eq!(owner_calls, 48);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_mt_m08_real_transport_matches_fixed_w06_vdw_matrix() {
        use crate::kernel::Cf3dFragAcceptContributionIdentity as Identity;

        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const IDS: [usize; 3] = [80, 7, 900];
        const WORKER_COUNTS: [i32; 6] = [-2, 0, 1, 2, 3, 4];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let threshold = 10.0 * minimum;
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = pre_factor * (0.0 - distance) / distance;
            let second_gradient = pre_factor * (distance - 0.0) / distance;
            let first_gradient = first_gradient * 0.1;
            let second_gradient = second_gradient * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let reference_points = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut reference = worker_w06_coordinates(&reference_points);
        let source_field =
            worker_w06_source_field(&topology, &mut reference, &total_valences, 10.0, false);
        let expected_identity = vec![Identity::Vdw {
            at1_idx: 0,
            at2_idx: 1,
            x_ij: minimum,
            well_depth,
            threshold,
        }];
        assert_eq!(
            crate::kernel::cf3d_frag_accept_contribution_identities(&source_field),
            expected_identity
        );
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut owner_calls = 0;

        for row_count in 0..=TARGET_DISTANCES.len() {
            for worker_count in WORKER_COUNTS {
                for max_iterations in [0_i32, 1] {
                    let mut coordinate_rows = TARGET_DISTANCES[..row_count]
                        .iter()
                        .copied()
                        .map(|distance| vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
                        .collect::<Vec<_>>();
                    let before = coordinate_rows.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(input_index, positions)| super::SerialConformer {
                            id: IDS[input_index],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; row_count + 1];
                    cf3d_uff_one_kernel_counts_reset();
                    let outcomes = super::optimize_conformers_mt(
                        &source_field,
                        &mut conformers,
                        &mut results,
                        worker_count,
                        2,
                        max_iterations,
                    )
                    .expect("the fixed W06 result capacity is sufficient");
                    owner_calls += 1;

                    if worker_count <= 0 {
                        assert!(outcomes.is_empty());
                    } else {
                        assert_eq!(outcomes.len(), worker_count as usize);
                        for (lane, outcome) in outcomes.into_iter().enumerate() {
                            assert!(
                                matches!(outcome, Ok(Ok(()))),
                                "W06 lane={lane}, rows={row_count}, workers={worker_count}, iterations={max_iterations}"
                            );
                        }
                    }
                    assert_eq!(
                        conformers
                            .iter()
                            .map(|conformer| conformer.id)
                            .collect::<Vec<_>>(),
                        IDS[..row_count]
                    );
                    drop(conformers);

                    for source_index in 0..row_count {
                        let expected_coordinates = if worker_count > 0 && max_iterations == 1 {
                            let [first_x, second_x] =
                                source_one_step(TARGET_DISTANCES[source_index]).0;
                            [[first_x, 0.0, 0.0], [second_x, 0.0, 0.0]]
                        } else {
                            before[source_index].clone().try_into().unwrap()
                        };
                        for point_index in 0..2 {
                            for axis in 0..3 {
                                assert_eq!(
                                    coordinate_rows[source_index][point_index][axis].to_bits(),
                                    expected_coordinates[point_index][axis].to_bits(),
                                    "W06 coordinate row={source_index}, point={point_index}, axis={axis}, rows={row_count}, workers={worker_count}, iterations={max_iterations}"
                                );
                            }
                        }
                        if worker_count > 0 {
                            let expected_energy = if max_iterations == 0 {
                                source_energy(TARGET_DISTANCES[source_index])
                            } else {
                                source_one_step(TARGET_DISTANCES[source_index]).1
                            };
                            assert_eq!(
                                results[source_index].status, 1,
                                "W06 source-index status row={source_index}, workers={worker_count}, iterations={max_iterations}"
                            );
                            assert!(
                                (results[source_index].energy - expected_energy).abs() <= 1.0e-12,
                                "W06 fixed VDW energy row={source_index}, workers={worker_count}, iterations={max_iterations}"
                            );
                        } else {
                            assert_eq!(results[source_index], sentinel);
                        }
                    }
                    assert_eq!(results[row_count], sentinel);
                    assert_eq!(source_field.positions().len(), 2);
                    assert_eq!(
                        crate::kernel::cf3d_frag_accept_contribution_identities(&source_field),
                        expected_identity
                    );
                    assert_eq!(source_field.positions().len(), reference_points.len());
                    for (point_index, expected_point) in reference_points.iter().enumerate() {
                        for axis in 0..3 {
                            assert_eq!(
                                source_field.positions()[point_index][axis].to_bits(),
                                expected_point[axis].to_bits(),
                                "worker field copies leave source position {point_index} axis {axis} unchanged"
                            );
                        }
                    }
                    let spawned_lanes = worker_count.max(0) as usize;
                    assert_eq!(
                        cf3d_uff_one_kernel_counts(),
                        (spawned_lanes, 0, spawned_lanes),
                        "one source field and one VDW term copy per spawned lane"
                    );
                }
            }
        }
        drop(source_field);
        assert_eq!(reference.conformers_3d[0].coordinates(), reference_points);
        assert_eq!(owner_calls, 48);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_dispatch_d04_dispatches_fixed_w05_success_matrix() {
        const SOURCE_X: [f64; 3] = [2.0, 3.0, 4.0];
        const IDS: [usize; 3] = [900, 2, 80];
        const REQUESTED: [i32; 6] = [-2, 0, 1, 2, 3, 4];
        const HARDWARE: [u32; 3] = [0, 1, 3];
        // Literal RDThreads.h results for thread-enabled builds; a disabled
        // build mode always selects one serial worker for every input above.
        const THREADED_COUNTS: [[i32; 6]; 3] =
            [[1, 1, 1, 2, 3, 4], [1, 1, 1, 2, 3, 4], [1, 3, 1, 2, 3, 4]];
        let one_step_x = SOURCE_X.map(|x| x + (-(2.0 * x * 0.1)));
        let one_step_energies = one_step_x.map(|x| x * x);
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut owner_calls = 0;
        let mut underfilled_result_vectors = 0;
        let mut oversized_result_vectors = 0;

        for row_count in 0..=SOURCE_X.len() {
            for (requested_index, requested) in REQUESTED.into_iter().enumerate() {
                for (hardware_index, observed_hardware) in HARDWARE.into_iter().enumerate() {
                    for threadsafe in [false, true] {
                        let expected_source_count = if threadsafe {
                            THREADED_COUNTS[hardware_index][requested_index]
                        } else {
                            1
                        };
                        for max_iterations in [0_i32, 1] {
                            let mut coordinate_rows = SOURCE_X[..row_count]
                                .iter()
                                .copied()
                                .map(|x| vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]])
                                .collect::<Vec<_>>();
                            let mut conformers = coordinate_rows
                                .iter_mut()
                                .enumerate()
                                .map(|(input_index, positions)| super::SerialConformer {
                                    id: IDS[input_index],
                                    positions: positions.as_mut_slice(),
                                })
                                .collect::<Vec<_>>();
                            let initial_result_len = if (row_count
                                + requested_index
                                + hardware_index
                                + usize::from(threadsafe)
                                + max_iterations as usize)
                                % 2
                                == 0
                            {
                                row_count + 1
                            } else {
                                row_count.saturating_sub(1)
                            };
                            underfilled_result_vectors +=
                                usize::from(initial_result_len < row_count);
                            oversized_result_vectors += usize::from(initial_result_len > row_count);
                            let mut results = vec![sentinel; initial_result_len];
                            let mut field = super::ForceField::new(2);
                            field.add_contribution(Box::new(CoordinateNormContribution));

                            let outcome = super::optimize_prepared_conformers_dispatch(
                                field,
                                &mut conformers,
                                &mut results,
                                2,
                                max_iterations,
                                requested,
                                observed_hardware,
                                threadsafe,
                            )
                            .expect("all D04 thread inputs have defined resolution");
                            owner_calls += 1;

                            match outcome {
                                super::PreparedConformerDispatchOutcome::Serial(result) => {
                                    assert_eq!(
                                        expected_source_count, 1,
                                        "W05 serial route: rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}"
                                    );
                                    result.unwrap_or_else(|error| {
                                        panic!("fixed W05 serial row matrix succeeds: {error:?}")
                                    });
                                }
                                super::PreparedConformerDispatchOutcome::Workers(result) => {
                                    assert_ne!(
                                        expected_source_count, 1,
                                        "W05 worker route: rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}"
                                    );
                                    let joins = result.unwrap_or_else(|error| {
                                        panic!("fixed W05 worker setup succeeds: {error:?}")
                                    });
                                    assert_eq!(joins.len(), expected_source_count as usize);
                                    for (lane, outcome) in joins.into_iter().enumerate() {
                                        assert!(
                                            matches!(outcome, Ok(Ok(()))),
                                            "W05 lane={lane}, rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                        );
                                    }
                                }
                            }

                            assert_eq!(
                                conformers
                                    .iter()
                                    .map(|conformer| conformer.id)
                                    .collect::<Vec<_>>(),
                                IDS[..row_count],
                                "W05 conformer IDs remain labels in supplied row order"
                            );
                            drop(conformers);

                            assert_eq!(results.len(), row_count);
                            for source_index in 0..row_count {
                                let expected_x = if max_iterations == 0 {
                                    SOURCE_X[source_index]
                                } else {
                                    one_step_x[source_index]
                                };
                                for (point_index, expected_point) in
                                    [[0.0, 0.0, 0.0], [expected_x, 0.0, 0.0]]
                                        .into_iter()
                                        .enumerate()
                                {
                                    for (axis, expected_component) in
                                        expected_point.into_iter().enumerate()
                                    {
                                        assert_eq!(
                                            coordinate_rows[source_index][point_index][axis]
                                                .to_bits(),
                                            expected_component.to_bits(),
                                            "W05 fixed source coordinate: row={source_index}, point={point_index}, axis={axis}, rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                        );
                                    }
                                }
                                let expected_energy = if max_iterations == 0 {
                                    SOURCE_X[source_index] * SOURCE_X[source_index]
                                } else {
                                    one_step_energies[source_index]
                                };
                                assert_eq!(results[source_index].status, 1);
                                assert_eq!(
                                    results[source_index].energy.to_bits(),
                                    expected_energy.to_bits(),
                                    "W05 fixed source energy: row={source_index}, rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                );
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(owner_calls, 288);
        assert!(underfilled_result_vectors > 0);
        assert!(oversized_result_vectors > 0);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_dispatch_d04_dispatches_fixed_w06_success_matrix() {
        use crate::kernel::Cf3dFragAcceptContributionIdentity as Identity;

        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const IDS: [usize; 3] = [80, 7, 900];
        const REQUESTED: [i32; 6] = [-2, 0, 1, 2, 3, 4];
        const HARDWARE: [u32; 3] = [0, 1, 3];
        const REFERENCE_POINTS: [[f64; 3]; 2] = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        const THREADED_COUNTS: [[i32; 6]; 3] =
            [[1, 1, 1, 2, 3, 4], [1, 1, 1, 2, 3, 4], [1, 3, 1, 2, 3, 4]];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let threshold = 10.0 * minimum;
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = pre_factor * (0.0 - distance) / distance;
            let second_gradient = pre_factor * (distance - 0.0) / distance;
            let first_gradient = first_gradient * 0.1;
            let second_gradient = second_gradient * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let expected_identity = vec![Identity::Vdw {
            at1_idx: 0,
            at2_idx: 1,
            x_ij: minimum,
            well_depth,
            threshold,
        }];
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut owner_calls = 0;
        let mut underfilled_result_vectors = 0;
        let mut oversized_result_vectors = 0;

        for row_count in 0..=TARGET_DISTANCES.len() {
            for (requested_index, requested) in REQUESTED.into_iter().enumerate() {
                for (hardware_index, observed_hardware) in HARDWARE.into_iter().enumerate() {
                    for threadsafe in [false, true] {
                        let expected_source_count = if threadsafe {
                            THREADED_COUNTS[hardware_index][requested_index]
                        } else {
                            1
                        };
                        for max_iterations in [0_i32, 1] {
                            let mut reference = worker_w06_coordinates(&REFERENCE_POINTS);
                            let field = worker_w06_source_field(
                                &topology,
                                &mut reference,
                                &total_valences,
                                10.0,
                                false,
                            );
                            assert_eq!(
                                crate::kernel::cf3d_frag_accept_contribution_identities(&field),
                                expected_identity,
                                "the fixed W06 matrix retains its pinned VDW source term"
                            );
                            let mut coordinate_rows = TARGET_DISTANCES[..row_count]
                                .iter()
                                .copied()
                                .map(|distance| vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
                                .collect::<Vec<_>>();
                            let mut conformers = coordinate_rows
                                .iter_mut()
                                .enumerate()
                                .map(|(input_index, positions)| super::SerialConformer {
                                    id: IDS[input_index],
                                    positions: positions.as_mut_slice(),
                                })
                                .collect::<Vec<_>>();
                            let initial_result_len = if (row_count
                                + requested_index
                                + hardware_index
                                + usize::from(threadsafe)
                                + max_iterations as usize)
                                % 2
                                == 0
                            {
                                row_count + 1
                            } else {
                                row_count.saturating_sub(1)
                            };
                            underfilled_result_vectors +=
                                usize::from(initial_result_len < row_count);
                            oversized_result_vectors += usize::from(initial_result_len > row_count);
                            let mut results = vec![sentinel; initial_result_len];

                            let outcome = super::optimize_prepared_conformers_dispatch(
                                field,
                                &mut conformers,
                                &mut results,
                                2,
                                max_iterations,
                                requested,
                                observed_hardware,
                                threadsafe,
                            )
                            .expect("all D04 thread inputs have defined resolution");
                            owner_calls += 1;

                            match outcome {
                                super::PreparedConformerDispatchOutcome::Serial(result) => {
                                    assert_eq!(
                                        expected_source_count, 1,
                                        "W06 serial route: rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}"
                                    );
                                    result.unwrap_or_else(|error| {
                                        panic!("fixed W06 serial row matrix succeeds: {error:?}")
                                    });
                                }
                                super::PreparedConformerDispatchOutcome::Workers(result) => {
                                    assert_ne!(
                                        expected_source_count, 1,
                                        "W06 worker route: rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}"
                                    );
                                    let joins = result.unwrap_or_else(|error| {
                                        panic!("fixed W06 worker setup succeeds: {error:?}")
                                    });
                                    assert_eq!(joins.len(), expected_source_count as usize);
                                    for (lane, outcome) in joins.into_iter().enumerate() {
                                        assert!(
                                            matches!(outcome, Ok(Ok(()))),
                                            "W06 lane={lane}, rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                        );
                                    }
                                }
                            }

                            assert_eq!(
                                conformers
                                    .iter()
                                    .map(|conformer| conformer.id)
                                    .collect::<Vec<_>>(),
                                IDS[..row_count],
                                "W06 conformer IDs remain labels in supplied row order"
                            );
                            drop(conformers);
                            assert_eq!(results.len(), row_count);
                            for source_index in 0..row_count {
                                let (expected_xs, expected_energy) = if max_iterations == 0 {
                                    (
                                        [0.0, TARGET_DISTANCES[source_index]],
                                        source_energy(TARGET_DISTANCES[source_index]),
                                    )
                                } else {
                                    source_one_step(TARGET_DISTANCES[source_index])
                                };
                                for (point_index, expected_point) in
                                    [[expected_xs[0], 0.0, 0.0], [expected_xs[1], 0.0, 0.0]]
                                        .into_iter()
                                        .enumerate()
                                {
                                    for (axis, expected_component) in
                                        expected_point.into_iter().enumerate()
                                    {
                                        assert_eq!(
                                            coordinate_rows[source_index][point_index][axis]
                                                .to_bits(),
                                            expected_component.to_bits(),
                                            "W06 fixed source coordinate: row={source_index}, point={point_index}, axis={axis}, rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                        );
                                    }
                                }
                                assert_eq!(results[source_index].status, 1);
                                assert!(
                                    (results[source_index].energy - expected_energy).abs()
                                        <= 1.0e-12,
                                    "W06 fixed source VDW energy: row={source_index}, rows={row_count}, requested={requested}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                );
                            }
                            for (point_index, expected_point) in
                                REFERENCE_POINTS.into_iter().enumerate()
                            {
                                for (axis, expected_component) in
                                    expected_point.into_iter().enumerate()
                                {
                                    assert_eq!(
                                        reference.conformers_3d[0].coordinates()[point_index][axis]
                                            .to_bits(),
                                        expected_component.to_bits(),
                                        "W06 construction reference remains unchanged"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(owner_calls, 288);
        assert!(underfilled_result_vectors > 0);
        assert!(oversized_result_vectors > 0);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_dispatch_d05_resizes_before_thread_error_and_honors_no_thread_mode() {
        const SOURCE_X: [f64; 3] = [2.0, 3.0, 4.0];
        const IDS: [usize; 3] = [900, 2, 80];
        const PREFIX_SENTINELS: [(i32, u64); 5] = [
            (17, 0x3ff8_0000_0000_0000),
            (-8, 0xc004_0000_0000_0000),
            (5, 0x4009_21fb_5444_2d18),
            (-13, 0x0000_0000_0000_0001),
            (29, 0xbff0_0000_0000_0000),
        ];
        let mut owner_calls = 0;

        for initial_len in [0_usize, 1, 5] {
            let mut coordinate_rows = SOURCE_X
                .iter()
                .copied()
                .map(|x| vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]])
                .collect::<Vec<_>>();
            let before = coordinate_rows.clone();
            let mut conformers = coordinate_rows
                .iter_mut()
                .enumerate()
                .map(|(row, positions)| super::SerialConformer {
                    id: IDS[row],
                    positions: positions.as_mut_slice(),
                })
                .collect::<Vec<_>>();
            let mut results = PREFIX_SENTINELS[..initial_len]
                .iter()
                .map(|&(status, energy_bits)| super::OptimizationOutcome {
                    status,
                    energy: f64::from_bits(energy_bits),
                })
                .collect::<Vec<_>>();

            let result = super::optimize_prepared_conformers_dispatch(
                super::ForceField::new(2),
                &mut conformers,
                &mut results,
                2,
                0,
                i32::MIN,
                0,
                true,
            );
            owner_calls += 1;
            assert!(matches!(
                result,
                Err(super::UffThreadCountError::UndefinedSignedNegation)
            ));
            drop(conformers);
            assert_eq!(
                coordinate_rows, before,
                "resolution failure precedes all row work"
            );

            assert_eq!(results.len(), SOURCE_X.len());
            for row in 0..SOURCE_X.len() {
                if row < initial_len {
                    assert_eq!(results[row].status, PREFIX_SENTINELS[row].0);
                    assert_eq!(
                        results[row].energy.to_bits(),
                        PREFIX_SENTINELS[row].1,
                        "resize preserves original pair bits at row {row}"
                    );
                } else {
                    assert_eq!(results[row].status, 0);
                    assert_eq!(
                        results[row].energy.to_bits(),
                        0.0_f64.to_bits(),
                        "resize appends source pair positive zero at row {row}"
                    );
                }
            }
        }

        let mut coordinate_rows = SOURCE_X
            .iter()
            .copied()
            .map(|x| vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]])
            .collect::<Vec<_>>();
        let before = coordinate_rows.clone();
        let mut conformers = coordinate_rows
            .iter_mut()
            .enumerate()
            .map(|(row, positions)| super::SerialConformer {
                id: IDS[row],
                positions: positions.as_mut_slice(),
            })
            .collect::<Vec<_>>();
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut results = vec![sentinel];
        let mut field = super::ForceField::new(2);
        field.add_contribution(Box::new(CoordinateNormContribution));
        let outcome = super::optimize_prepared_conformers_dispatch(
            field,
            &mut conformers,
            &mut results,
            2,
            0,
            i32::MIN,
            0,
            false,
        )
        .expect("source no-thread mode ignores INT_MIN before negation");
        owner_calls += 1;
        match outcome {
            super::PreparedConformerDispatchOutcome::Serial(Ok(())) => {}
            _ => panic!("source no-thread mode must execute the serial owner"),
        }
        drop(conformers);
        assert_eq!(
            coordinate_rows, before,
            "zero-iteration serial rows are unchanged"
        );
        assert_eq!(results.len(), SOURCE_X.len());
        for row in 0..SOURCE_X.len() {
            assert_eq!(results[row].status, 1);
            assert_eq!(
                results[row].energy.to_bits(),
                (SOURCE_X[row] * SOURCE_X[row]).to_bits()
            );
        }
        assert_eq!(owner_calls, 4);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_dispatch_d05_malformed_rows_preserve_serial_prefix_and_raw_lane_order() {
        const SOURCE_X: [f64; 3] = [2.0, 3.0, 4.0];
        const IDS: [usize; 3] = [900, 2, 80];
        const REQUESTED: [i32; 4] = [1, 2, 3, 4];
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut owner_calls = 0;

        for failed_index in 0..SOURCE_X.len() {
            for requested in REQUESTED {
                for threadsafe in [false, true] {
                    let mut coordinate_rows = SOURCE_X
                        .iter()
                        .enumerate()
                        .map(|(row, &x)| {
                            if row == failed_index {
                                vec![[x, -1.0, 0.5]]
                            } else {
                                vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]]
                            }
                        })
                        .collect::<Vec<_>>();
                    let before = coordinate_rows.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(row, positions)| super::SerialConformer {
                            id: IDS[row],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; SOURCE_X.len() + 2];
                    let mut field = super::ForceField::new(2);
                    field.add_contribution(Box::new(CoordinateNormContribution));
                    let outcome = super::optimize_prepared_conformers_dispatch(
                        field,
                        &mut conformers,
                        &mut results,
                        2,
                        0,
                        requested,
                        0,
                        threadsafe,
                    )
                    .expect("positive requested worker counts have defined routes");
                    owner_calls += 1;

                    let expected_failure =
                        super::SerialConformerOptimizationError::CoordinateCount {
                            input_index: failed_index,
                            source: super::SerialCoordinateCountError {
                                conformer_id: IDS[failed_index],
                                atoms: 2,
                                coordinates: 1,
                            },
                        };
                    let mut expected_results = vec![sentinel; SOURCE_X.len()];
                    let expected_source_count = if threadsafe && requested > 1 {
                        requested as usize
                    } else {
                        1
                    };

                    if expected_source_count == 1 {
                        for row in 0..failed_index {
                            expected_results[row] = super::OptimizationOutcome {
                                status: 1,
                                energy: SOURCE_X[row] * SOURCE_X[row],
                            };
                        }
                        match outcome {
                            super::PreparedConformerDispatchOutcome::Serial(Err(error)) => {
                                assert_eq!(error, expected_failure);
                            }
                            _ => panic!("serial dispatch must stop at the first malformed row"),
                        }
                    } else {
                        let joins = match outcome {
                            super::PreparedConformerDispatchOutcome::Workers(Ok(joins)) => joins,
                            _ => panic!("threaded dispatch must retain every raw lane outcome"),
                        };
                        assert_eq!(joins.len(), expected_source_count);
                        for (lane, outcome) in joins.iter().enumerate() {
                            let lane_indices =
                                (lane..SOURCE_X.len()).step_by(expected_source_count);
                            let mut expected_lane_error = false;
                            for row in lane_indices {
                                if row == failed_index {
                                    expected_lane_error = true;
                                    break;
                                }
                                expected_results[row] = super::OptimizationOutcome {
                                    status: 1,
                                    energy: SOURCE_X[row] * SOURCE_X[row],
                                };
                            }
                            if expected_lane_error {
                                match outcome {
                                    Ok(Err(error)) => assert_eq!(*error, expected_failure),
                                    _ => panic!(
                                        "creation-order lane {lane} must retain its malformed-row error"
                                    ),
                                }
                            } else {
                                assert!(
                                    matches!(outcome, Ok(Ok(()))),
                                    "creation-order lane {lane} completes its assigned source rows"
                                );
                            }
                        }
                    }

                    assert_eq!(results.len(), SOURCE_X.len());
                    assert_eq!(results, expected_results);
                    assert_eq!(
                        conformers
                            .iter()
                            .map(|conformer| conformer.id)
                            .collect::<Vec<_>>(),
                        IDS
                    );
                    drop(conformers);
                    assert_eq!(
                        coordinate_rows, before,
                        "failed rows and zero-iteration peers are unchanged"
                    );
                }
            }
        }
        assert_eq!(owner_calls, 24);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_dispatch_d05_retains_w05_w06_stage_failures_and_effects() {
        const IDS: [usize; 6] = [80, 7, 900, 31, 502, 8];
        const SOURCE_X: [f64; 6] = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        const W06_DISTANCES: [f64; 6] = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        const WORKER_COUNTS: [i32; 2] = [2, 3];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let w05_one_step_x = SOURCE_X.map(|x| x + (-(2.0 * x * 0.1)));
        let w05_one_step_energy = w05_one_step_x.map(|x| x * x);
        let mut owner_calls = [0_usize; 2];

        for worker_count in WORKER_COUNTS {
            // RDKit assigns source indices by modulo. One extra row gives
            // lane zero two rows and every other lane one, so the per-copy
            // failure counter reaches its fixed second-row stage only there.
            let row_count = worker_count as usize + 1;
            let failed_lane = 0_usize;
            let failed_index = failed_lane + worker_count as usize;
            for (max_iterations, fail_on) in [(0_i32, 3_usize), (1_i32, 6_usize)] {
                let mut coordinate_rows = SOURCE_X[..row_count]
                    .iter()
                    .copied()
                    .map(|x| vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]])
                    .collect::<Vec<_>>();
                let mut expected_coordinates = coordinate_rows.clone();
                let mut conformers = coordinate_rows
                    .iter_mut()
                    .enumerate()
                    .map(|(row, positions)| super::SerialConformer {
                        id: IDS[row],
                        positions: positions.as_mut_slice(),
                    })
                    .collect::<Vec<_>>();
                let mut results = vec![SENTINEL; row_count + 1];
                let mut field = super::ForceField::new(3);
                field.add_contribution(Box::new(FailAtEnergyCallWithCoordinateNorm {
                    calls: Cell::new(0),
                    fail_on,
                }));

                let outcome = super::optimize_prepared_conformers_dispatch(
                    field,
                    &mut conformers,
                    &mut results,
                    2,
                    max_iterations,
                    worker_count,
                    0,
                    true,
                )
                .expect("positive worker count selects the source threaded branch");
                owner_calls[0] += 1;
                let joins = match outcome {
                    super::PreparedConformerDispatchOutcome::Workers(Ok(joins)) => joins,
                    _ => panic!("W05 failure matrix retains actual threaded dispatch outcomes"),
                };
                assert_eq!(joins.len(), worker_count as usize);

                let expected_stage = if max_iterations == 0 {
                    super::OptimizationStageError::Minimize(expected_pair_index_error())
                } else {
                    super::OptimizationStageError::FinalEnergy(expected_pair_index_error())
                };
                let expected_error = super::SerialConformerOptimizationError::Optimization {
                    input_index: failed_index,
                    conformer_id: IDS[failed_index],
                    source: expected_stage,
                };
                let mut expected_results = vec![SENTINEL; row_count];
                for lane in 0..worker_count as usize {
                    let lane_indices = (lane..row_count).step_by(worker_count as usize);
                    let mut lane_failed = false;
                    for row in lane_indices {
                        if row == failed_index {
                            lane_failed = true;
                            if max_iterations == 1 {
                                expected_coordinates[row][1][0] = w05_one_step_x[row];
                            }
                            break;
                        }
                        if max_iterations == 1 {
                            expected_coordinates[row][1][0] = w05_one_step_x[row];
                        }
                        expected_results[row] = super::OptimizationOutcome {
                            status: 1,
                            energy: if max_iterations == 0 {
                                SOURCE_X[row] * SOURCE_X[row]
                            } else {
                                w05_one_step_energy[row]
                            },
                        };
                    }
                    if lane_failed {
                        assert!(matches!(
                            &joins[lane],
                            Ok(Err(error)) if *error == expected_error
                        ));
                    } else {
                        assert!(matches!(&joins[lane], Ok(Ok(()))));
                    }
                }
                assert_eq!(results.len(), row_count);
                assert_eq!(results, expected_results);
                assert_eq!(
                    conformers
                        .iter()
                        .map(|conformer| conformer.id)
                        .collect::<Vec<_>>(),
                    IDS[..row_count]
                );
                drop(conformers);
                for row in 0..row_count {
                    for point in 0..2 {
                        for axis in 0..3 {
                            assert_eq!(
                                coordinate_rows[row][point][axis].to_bits(),
                                expected_coordinates[row][point][axis].to_bits(),
                                "W05 detached effect row={row}, worker_count={worker_count}, iterations={max_iterations}, failed_index={failed_index}"
                            );
                        }
                    }
                }
            }
        }

        use crate::kernel::Cf3dFragAcceptContributionIdentity as Identity;

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let threshold = 10.0 * minimum;
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step_with_norm = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_vdw_gradient = pre_factor * (0.0 - distance) / distance;
            let second_vdw_gradient = pre_factor * (distance - 0.0) / distance;
            let first_gradient = (first_vdw_gradient + 2.0 * 0.0) * 0.1;
            let second_gradient = (second_vdw_gradient + 2.0 * distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            let norm_energy = first_x * first_x + second_x * second_x;
            (
                [first_x, second_x],
                source_energy(after_distance) + norm_energy,
            )
        };
        let expected_identity = [Identity::Vdw {
            at1_idx: 0,
            at2_idx: 1,
            x_ij: minimum,
            well_depth,
            threshold,
        }];

        for worker_count in WORKER_COUNTS {
            let row_count = worker_count as usize + 1;
            let failed_lane = 0_usize;
            let failed_index = failed_lane + worker_count as usize;
            for (max_iterations, fail_on) in [(0_i32, 3_usize), (1_i32, 6_usize)] {
                let mut reference = worker_w06_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
                let mut field = worker_w06_source_field(
                    &topology,
                    &mut reference,
                    &total_valences,
                    10.0,
                    false,
                );
                assert_eq!(
                    crate::kernel::cf3d_frag_accept_contribution_identities(&field),
                    expected_identity
                );
                field.add_contribution(Box::new(FailAtEnergyCallWithCoordinateNorm {
                    calls: Cell::new(0),
                    fail_on,
                }));
                let mut coordinate_rows = W06_DISTANCES[..row_count]
                    .iter()
                    .copied()
                    .map(|distance| vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
                    .collect::<Vec<_>>();
                let mut expected_coordinates = coordinate_rows.clone();
                let mut conformers = coordinate_rows
                    .iter_mut()
                    .enumerate()
                    .map(|(row, positions)| super::SerialConformer {
                        id: IDS[row],
                        positions: positions.as_mut_slice(),
                    })
                    .collect::<Vec<_>>();
                let mut results = vec![SENTINEL; row_count + 1];

                let outcome = super::optimize_prepared_conformers_dispatch(
                    field,
                    &mut conformers,
                    &mut results,
                    2,
                    max_iterations,
                    worker_count,
                    0,
                    true,
                )
                .expect("positive worker count selects the source threaded branch");
                owner_calls[1] += 1;
                let joins = match outcome {
                    super::PreparedConformerDispatchOutcome::Workers(Ok(joins)) => joins,
                    _ => panic!("W06 failure matrix retains actual threaded dispatch outcomes"),
                };
                assert_eq!(joins.len(), worker_count as usize);

                let expected_stage = if max_iterations == 0 {
                    super::OptimizationStageError::Minimize(expected_pair_index_error())
                } else {
                    super::OptimizationStageError::FinalEnergy(expected_pair_index_error())
                };
                let expected_error = super::SerialConformerOptimizationError::Optimization {
                    input_index: failed_index,
                    conformer_id: IDS[failed_index],
                    source: expected_stage,
                };
                let mut expected_results = vec![SENTINEL; row_count];
                for lane in 0..worker_count as usize {
                    let lane_indices = (lane..row_count).step_by(worker_count as usize);
                    let mut lane_failed = false;
                    for row in lane_indices {
                        if row == failed_index {
                            lane_failed = true;
                            if max_iterations == 1 {
                                expected_coordinates[row][0][0] =
                                    source_one_step_with_norm(W06_DISTANCES[row]).0[0];
                                expected_coordinates[row][1][0] =
                                    source_one_step_with_norm(W06_DISTANCES[row]).0[1];
                            }
                            break;
                        }
                        let (expected_xs, expected_energy) = if max_iterations == 0 {
                            let distance = W06_DISTANCES[row];
                            (
                                [0.0, distance],
                                source_energy(distance) + distance * distance,
                            )
                        } else {
                            source_one_step_with_norm(W06_DISTANCES[row])
                        };
                        if max_iterations == 1 {
                            expected_coordinates[row][0][0] = expected_xs[0];
                            expected_coordinates[row][1][0] = expected_xs[1];
                        }
                        expected_results[row] = super::OptimizationOutcome {
                            status: 1,
                            energy: expected_energy,
                        };
                    }
                    if lane_failed {
                        assert!(matches!(
                            &joins[lane],
                            Ok(Err(error)) if *error == expected_error
                        ));
                    } else {
                        assert!(matches!(&joins[lane], Ok(Ok(()))));
                    }
                }
                assert_eq!(results.len(), row_count);
                for row in 0..row_count {
                    assert_eq!(results[row].status, expected_results[row].status);
                    assert!(
                        (results[row].energy - expected_results[row].energy).abs() <= 1.0e-12,
                        "W06 typed failure lane preserves fixed source energy row={row}, worker_count={worker_count}, iterations={max_iterations}"
                    );
                }
                assert_eq!(
                    conformers
                        .iter()
                        .map(|conformer| conformer.id)
                        .collect::<Vec<_>>(),
                    IDS[..row_count]
                );
                drop(conformers);
                for row in 0..row_count {
                    for point in 0..2 {
                        for axis in 0..3 {
                            assert_eq!(
                                coordinate_rows[row][point][axis].to_bits(),
                                expected_coordinates[row][point][axis].to_bits(),
                                "W06 detached effect row={row}, worker_count={worker_count}, iterations={max_iterations}, failed_index={failed_index}"
                            );
                        }
                    }
                }
                for (point_index, expected_point) in [[0.0_f64, 0.0, 0.0], [2.0, 0.0, 0.0]]
                    .into_iter()
                    .enumerate()
                {
                    for (axis, expected_component) in expected_point.into_iter().enumerate() {
                        assert_eq!(
                            reference.conformers_3d[0].coordinates()[point_index][axis].to_bits(),
                            expected_component.to_bits(),
                            "W06 reference coordinates remain detached from processed rows"
                        );
                    }
                }
            }
        }

        assert_eq!(owner_calls, [4, 4]);
    }

    #[test]
    fn uff_block_l02_releases_full_uff_field_before_same_block_reborrow() {
        use crate::uff::atom_typer::{UffTypingDiagnostic, UffTypingDiagnosticKind};
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        const IDS: [usize; 3] = [80, 7, 900];
        const SOURCE_COORDINATES: [[[f64; 3]; 2]; 3] = [
            [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]],
        ];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; 2],
        };
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let expected_identity = [crate::kernel::Cf3dFragAcceptContributionIdentity::Vdw {
            at1_idx: 0,
            at2_idx: 1,
            x_ij: minimum,
            well_depth,
            threshold: 10.0 * minimum,
        }];
        let mut constructor_calls = 0;

        for (selected_row, selected_id) in IDS.into_iter().enumerate() {
            let mut coordinates = CoordinateBlock::default();
            for (row, id) in IDS.into_iter().enumerate() {
                coordinates.conformers_3d.push(Conformer3D::new(
                    id,
                    SOURCE_COORDINATES[row].to_vec(),
                    true,
                ));
            }
            let mut diagnostics = Vec::new();
            let mut options = super::SingleConformerOptions::for_conformer(selected_id);
            options.ignore_interfragment_interactions = false;

            crate::kernel::cf3d_uff_one_kernel_counts_reset();
            let field = super::construct_released_uff_field(
                &topology,
                &mut coordinates,
                super::UffAtomStateRef::SuppliedRows {
                    total_valences: &total_valences,
                    conjugated_presence: &conjugated,
                },
                &rings,
                &valence,
                &properties,
                &mut diagnostics,
                options,
            )
            .unwrap_or_else(|error| {
                panic!("fixed W06 selected id {selected_id} constructs: {error:?}")
            });
            constructor_calls += 1;

            assert!(field.positions().is_empty());
            assert_eq!(
                crate::kernel::cf3d_frag_accept_contribution_identities(&field),
                expected_identity,
                "W06 ordered source VDW identity for selected id {selected_id}"
            );
            assert!(diagnostics.is_empty());
            assert_eq!(
                crate::kernel::cf3d_uff_one_kernel_counts(),
                (1, 1, 0),
                "one builder field and one VDW contribution, no field copy"
            );
            assert_eq!(
                coordinates
                    .conformers_3d
                    .iter()
                    .map(Conformer3D::id)
                    .collect::<Vec<_>>(),
                IDS
            );

            // The force field stays live while the same block's source-ordered
            // rows are borrowed again and consumed by the existing serial loop.
            let mut conformers = coordinates
                .conformers_3d
                .iter_mut()
                .map(|conformer| super::SerialConformer {
                    id: conformer.id(),
                    positions: conformer.coordinates_mut(),
                })
                .collect::<Vec<_>>();
            assert_eq!(
                conformers
                    .iter()
                    .map(|conformer| conformer.id)
                    .collect::<Vec<_>>(),
                IDS
            );
            let mut results = [SENTINEL; 3];
            super::optimize_serial_conformers(field, &mut conformers, &mut results, 2, 0)
                .unwrap_or_else(|error| panic!("same-block W06 serial call succeeds: {error:?}"));
            drop(conformers);

            for row in 0..IDS.len() {
                assert_eq!(results[row].status, 1);
                assert!(
                    (results[row].energy - source_energy(SOURCE_COORDINATES[row][1][0])).abs()
                        <= 1.0e-12,
                    "fixed W06 source energy row={row}, selected id={selected_id}"
                );
                assert_eq!(
                    coordinates.conformers_3d[row].coordinates(),
                    SOURCE_COORDINATES[row],
                    "zero-iteration source coordinates row={row}, selected id={selected_id}"
                );
            }
            assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
            assert_eq!(
                selected_row,
                IDS.iter().position(|id| *id == selected_id).unwrap()
            );
        }
        assert_eq!(constructor_calls, 3);

        // Preserve the D06 source warning-before-missing-ID error path through
        // this released-field entry as well as the original dispatch test.
        const TBP_POINTS: [[f64; 3]; 6] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-0.6, 0.8, 0.0],
            [0.0, 0.6, 0.8],
            [0.0, -0.8, 0.6],
            [-0.8, 0.0, -0.6],
        ];
        let mut tbp_atoms = vec![worker_w06_atom(0, 15, Hybridization::Sp3d, 0)];
        tbp_atoms.extend((1..6).map(|row| worker_w06_atom(row, 6, Hybridization::Sp3, 0)));
        let tbp_topology = worker_w06_topology(
            tbp_atoms,
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
                (0, 5, BondOrder::Single),
            ],
        );
        let tbp_total_valences = [5, 1, 1, 1, 1, 1];
        let tbp_conjugated = [false; 6];
        let tbp_rings =
            fast_find_rings(&tbp_topology).expect("fixed TBP topology has ring information");
        let tbp_valence = ValenceAssignment {
            explicit_valence: tbp_total_valences.to_vec(),
            implicit_hydrogens: vec![0; 6],
        };
        let mut tbp_coordinates = CoordinateBlock::default();
        for id in IDS {
            tbp_coordinates
                .conformers_3d
                .push(Conformer3D::new(id, TBP_POINTS.to_vec(), true));
        }
        let mut diagnostics = Vec::new();
        let mut missing_id_options = super::SingleConformerOptions::for_conformer(99);
        missing_id_options.ignore_interfragment_interactions = false;
        crate::kernel::cf3d_uff_one_kernel_counts_reset();
        let error = match super::construct_released_uff_field(
            &tbp_topology,
            &mut tbp_coordinates,
            super::UffAtomStateRef::SuppliedRows {
                total_valences: &tbp_total_valences,
                conjugated_presence: &tbp_conjugated,
            },
            &tbp_rings,
            &tbp_valence,
            &properties,
            &mut diagnostics,
            missing_id_options,
        ) {
            Err(error) => error,
            Ok(_) => panic!("missing selected 3D ID fails during source construction"),
        };
        assert!(matches!(
            error,
            super::AutomaticForceFieldConstructionError::Construction(
                ForceFieldConstructionError::Builder(
                    UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 99 }
                )
            )
        ));
        assert_eq!(
            diagnostics,
            [UffTypingDiagnostic {
                atom_id: Some(AtomId::new(0)),
                kind: UffTypingDiagnosticKind::Warning,
                message_prefix: "UFFTYPER: Warning: hybridization set to SP3 for atom ",
            }]
        );
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts().0, 1);
        assert_eq!(
            tbp_coordinates
                .conformers_3d
                .iter()
                .map(Conformer3D::id)
                .collect::<Vec<_>>(),
            IDS
        );
    }

    #[test]
    fn uff_prepare_p08_released_cached_field_keeps_storage_and_error_order() {
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        const IDS: [usize; 3] = [80, 7, 900];
        const DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let prepared = super::super::api::prepare_parameter_query(&topology, &valence)
            .expect("fixed cached W06 typing input validates");
        let typing_state = prepared.typing_state;
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            well_depth * (ratio6 * ratio6 - 2.0 * ratio6)
        };
        let expected_identity = [crate::kernel::Cf3dFragAcceptContributionIdentity::Vdw {
            at1_idx: 0,
            at2_idx: 1,
            x_ij: minimum,
            well_depth,
            threshold: 10.0 * minimum,
        }];
        let mut actual_calls = 0;

        for row_count in [1_usize, 3] {
            let mut coordinates = CoordinateBlock::default();
            for row in 0..row_count {
                coordinates.conformers_3d.push(Conformer3D::new(
                    IDS[row],
                    [[0.0, 0.0, 0.0], [DISTANCES[row], 0.0, 0.0]].to_vec(),
                    true,
                ));
            }
            let before = coordinates.clone();
            let selected_id = IDS[row_count - 1];
            let mut options = super::SingleConformerOptions::for_conformer(selected_id);
            options.max_iterations = 0;
            options.ignore_interfragment_interactions = false;
            let mut diagnostics = Vec::new();

            crate::kernel::cf3d_uff_one_kernel_counts_reset();
            let mut field = super::construct_released_uff_field(
                &topology,
                &mut coordinates,
                typing_state,
                &rings,
                &valence,
                &properties,
                &mut diagnostics,
                options,
            )
            .unwrap_or_else(|error| {
                panic!("fixed cached p08 row_count={row_count} constructs: {error:?}")
            });
            actual_calls += 1;

            assert!(field.positions().is_empty());
            let released_position_address = field.positions().as_ptr() as usize;
            let released_position_capacity = field.positions_mut().capacity();
            assert_ne!(released_position_address, 0);
            assert!(released_position_capacity >= topology.atoms.len());
            assert_eq!(
                crate::kernel::cf3d_frag_accept_contribution_identities(&field),
                expected_identity,
                "fixed source term identity before same-block reborrow, rows={row_count}"
            );
            assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
            assert!(diagnostics.is_empty());

            let mut conformers = coordinates
                .conformers_3d
                .iter_mut()
                .map(|conformer| super::SerialConformer {
                    id: conformer.id(),
                    positions: conformer.coordinates_mut(),
                })
                .collect::<Vec<_>>();
            let mut results = vec![SENTINEL; row_count];
            super::optimize_serial_conformers(
                field,
                &mut conformers,
                &mut results,
                topology.atoms.len(),
                0,
            )
            .unwrap_or_else(|error| {
                panic!("fixed cached p08 same-block rows={row_count} optimize: {error:?}")
            });
            drop(conformers);
            actual_calls += 1;

            let position_work = crate::kernel::cf3d_uff_one_serial_work_counts();
            assert_eq!(position_work.position_buffer_growths, 0);
            assert_eq!(
                position_work.initial_position_buffer_address, released_position_address,
                "the released field's handle Vec is reused for rows={row_count}"
            );
            assert_eq!(
                position_work.final_position_buffer_address, released_position_address,
                "serial row rebinding keeps the same handle allocation, rows={row_count}"
            );
            assert_eq!(
                position_work.initial_position_buffer_capacity,
                released_position_capacity
            );
            assert_eq!(
                position_work.final_position_buffer_capacity,
                released_position_capacity
            );
            assert_eq!(position_work.initialize_calls, row_count);
            assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
            assert_eq!(results.len(), row_count);
            for row in 0..row_count {
                assert_eq!(results[row].status, 1);
                assert!(
                    (results[row].energy - source_energy(DISTANCES[row])).abs() <= 1.0e-12,
                    "fixed W06 source energy rows={row_count}, row={row}"
                );
                assert_eq!(coordinates.conformers_3d[row], before.conformers_3d[row]);
            }
        }
        assert_eq!(actual_calls, 4);

        // The source constructor/typing warning precedes selected-ID lookup;
        // failure remains before the caller resizes its result vector.
        let mut tbp_atoms = vec![worker_w06_atom(0, 15, Hybridization::Sp3d, 0)];
        tbp_atoms.extend((1..6).map(|row| worker_w06_atom(row, 6, Hybridization::Sp3, 0)));
        let tbp_topology = worker_w06_topology(
            tbp_atoms,
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
                (0, 5, BondOrder::Single),
            ],
        );
        let tbp_total_valences = [5, 1, 1, 1, 1, 1];
        let tbp_conjugation = [false; 6];
        let tbp_valence = ValenceAssignment {
            explicit_valence: tbp_total_valences.to_vec(),
            implicit_hydrogens: vec![0; 6],
        };
        let tbp_rings = fast_find_rings(&tbp_topology).expect("fixed TBP topology has rings");
        let mut zero_rows = CoordinateBlock::default();
        zero_rows.conformers_2d.push(
            cosmolkit_model::Conformer2D::new(99, vec![[17.0, -3.5]])
                .with_prop("layout", "2d-does-not-satisfy-3d-lookup"),
        );
        let zero_rows_before = zero_rows.clone();
        let mut results = vec![SENTINEL];
        let results_before = results.clone();
        let mut missing_options = super::SingleConformerOptions::for_conformer(99);
        missing_options.ignore_interfragment_interactions = false;
        let mut diagnostics = Vec::new();

        crate::kernel::cf3d_uff_one_kernel_counts_reset();
        let error = super::optimize_serial_uff_coordinate_block(
            &tbp_topology,
            &mut zero_rows,
            &mut results,
            super::UffAtomStateRef::SuppliedRows {
                total_valences: &tbp_total_valences,
                conjugated_presence: &tbp_conjugation,
            },
            &tbp_rings,
            &tbp_valence,
            &properties,
            &mut diagnostics,
            missing_options,
        )
        .expect_err("zero stored 3D rows leave selected-ID construction failure visible");
        assert!(matches!(
            error,
            super::SerialUffOptimizationError::Construction(
                super::AutomaticForceFieldConstructionError::Construction(
                    ForceFieldConstructionError::Builder(
                        UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                            conformer_id: 99
                        }
                    )
                )
            )
        ));
        assert_eq!(results, results_before, "construction fails before resize");
        assert_eq!(zero_rows, zero_rows_before);
        assert_eq!(
            diagnostics,
            [crate::uff::atom_typer::UffTypingDiagnostic {
                atom_id: Some(AtomId::new(0)),
                kind: crate::uff::atom_typer::UffTypingDiagnosticKind::Warning,
                message_prefix: "UFFTYPER: Warning: hybridization set to SP3 for atom ",
            }]
        );
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 0, 0));

        // A present selected row with a short atom-position vector fails at
        // the builder's exact coordinate-count check, also before result resize.
        let mut short_coordinates = worker_w06_coordinates(&[[0.0, 0.0, 0.0]]);
        let short_before = short_coordinates.clone();
        let mut short_results = vec![SENTINEL; 2];
        let short_results_before = short_results.clone();
        let mut short_options = super::SingleConformerOptions::for_conformer(31);
        short_options.ignore_interfragment_interactions = false;
        let mut short_diagnostics = Vec::new();

        crate::kernel::cf3d_uff_one_kernel_counts_reset();
        let short_error = super::optimize_serial_uff_coordinate_block(
            &topology,
            &mut short_coordinates,
            &mut short_results,
            super::UffAtomStateRef::SuppliedRows {
                total_valences: &[0, 0],
                conjugated_presence: &[false; 2],
            },
            &rings,
            &valence,
            &properties,
            &mut short_diagnostics,
            short_options,
        )
        .expect_err("short selected 3D row retains typed builder error");
        assert!(matches!(
            short_error,
            super::SerialUffOptimizationError::Construction(
                super::AutomaticForceFieldConstructionError::Construction(
                    ForceFieldConstructionError::Builder(
                        UffBuilderError::SelectedConformerCoordinateCountMismatch {
                            conformer_id: 31,
                            atoms: 2,
                            coordinates: 1
                        }
                    )
                )
            )
        ));
        assert_eq!(short_results, short_results_before);
        assert_eq!(short_coordinates, short_before);
        assert!(short_diagnostics.is_empty());
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 0, 0));
    }

    #[test]
    fn uff_prepare_p09_prepared_serial_uses_cached_rows_and_preserves_lazy_matrix() {
        use super::super::optimization::UffPreparedOptimizationError;
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        const IDS: [usize; 3] = [80, 7, 900];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let valence = ValenceAssignment {
            explicit_valence: vec![0, 0],
            implicit_hydrogens: vec![0, 0],
        };
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let mut actual_calls = 0;

        // A matching 2D ID is not a selected 3D row. Preparation succeeds,
        // then source construction fails before result resizing or row access.
        let mut zero_coordinates = CoordinateBlock::default();
        zero_coordinates.conformers_2d.push(
            cosmolkit_model::Conformer2D::new(IDS[0], TWO_D_SENTINEL.to_vec())
                .with_prop("layout", "preserve-this-2d-row"),
        );
        let zero_before = zero_coordinates.clone();
        let mut zero_results = vec![SENTINEL; 2];
        let mut zero_options = super::SingleConformerOptions::for_conformer(IDS[0]);
        zero_options.max_iterations = 0;
        zero_options.ignore_interfragment_interactions = false;
        super::super::api::reset_prepare_parameter_query_calls();
        super::uff_public_perf_serial_adapter_probe_start();
        let zero_error = super::super::optimization::optimize_prepared_uff_serial(
            &topology,
            &mut zero_coordinates,
            &mut zero_results,
            &valence,
            &rings,
            &properties,
            &mut Vec::new(),
            zero_options,
        )
        .expect_err("empty 3D storage retains the selected-ID source error");
        actual_calls += 1;
        let zero_probe = super::uff_public_perf_serial_adapter_probe_finish();
        assert!(matches!(
            zero_error,
            UffPreparedOptimizationError::Serial(super::SerialUffOptimizationError::Construction(
                super::AutomaticForceFieldConstructionError::Construction(
                    ForceFieldConstructionError::Builder(
                        UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                            conformer_id: 80
                        }
                    )
                )
            ))
        ));
        assert_eq!(super::super::api::prepare_parameter_query_calls(), 1);
        assert_eq!(zero_probe, Some((0, None)));
        assert_eq!(zero_results, [SENTINEL; 2]);
        assert_eq!(zero_coordinates, zero_before);

        // The fixed W06 source reference covers both no-iteration and one-step
        // results while caller excess capacity is retained through truncation.
        for max_iterations in [0, 1] {
            let mut coordinates = CoordinateBlock::default();
            coordinates.conformers_3d.push(
                Conformer3D::new(IDS[0], [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]].to_vec(), true)
                    .with_prop("source-row", "w06-row-0"),
            );
            coordinates.conformers_2d.push(
                cosmolkit_model::Conformer2D::new(IDS[0], TWO_D_SENTINEL.to_vec())
                    .with_prop("layout", "preserve-this-overlapping-2d-row"),
            );
            let before_3d = coordinates.conformers_3d.clone();
            let before_2d = coordinates.conformers_2d.clone();
            let before_dimension = coordinates.source_coordinate_dim;
            let mut results = Vec::with_capacity(2);
            results.resize(2, SENTINEL);
            let result_storage = results.as_ptr();
            let mut options = super::SingleConformerOptions::for_conformer(IDS[0]);
            options.max_iterations = max_iterations;
            options.ignore_interfragment_interactions = false;

            super::super::api::reset_prepare_parameter_query_calls();
            super::uff_public_perf_serial_adapter_probe_start();
            super::super::optimization::optimize_prepared_uff_serial(
                &topology,
                &mut coordinates,
                &mut results,
                &valence,
                &rings,
                &properties,
                &mut Vec::new(),
                options,
            )
            .unwrap_or_else(|error| {
                panic!("fixed prepared W06 row succeeds at {max_iterations} iterations: {error:?}")
            });
            actual_calls += 1;
            let probe = super::uff_public_perf_serial_adapter_probe_finish();

            assert_eq!(super::super::api::prepare_parameter_query_calls(), 1);
            assert_eq!(probe, Some((1, Some(0))));
            assert_eq!(results.as_ptr(), result_storage);
            assert_eq!(results.len(), 1);
            assert_eq!(results[0].status, 1);
            let (expected_xs, expected_energy) = if max_iterations == 0 {
                ([0.0, 4.0], source_energy(4.0))
            } else {
                source_one_step(4.0)
            };
            assert!((results[0].energy - expected_energy).abs() <= 1.0e-12);
            let expected_points = [[expected_xs[0], 0.0, 0.0], [expected_xs[1], 0.0, 0.0]];
            for (point, expected_point) in expected_points.into_iter().enumerate() {
                for (axis, expected_component) in expected_point.into_iter().enumerate() {
                    assert_eq!(
                        coordinates.conformers_3d[0].coordinates()[point][axis].to_bits(),
                        expected_component.to_bits(),
                        "fixed prepared W06 coordinate point={point}, axis={axis}, iterations={max_iterations}"
                    );
                }
            }
            assert_eq!(coordinates.conformers_3d[0].id(), before_3d[0].id());
            assert_eq!(coordinates.conformers_3d[0].is_3d(), before_3d[0].is_3d());
            assert_eq!(coordinates.conformers_3d[0].props(), before_3d[0].props());
            assert_eq!(coordinates.conformers_2d, before_2d);
            assert_eq!(coordinates.source_coordinate_dim, before_dimension);
        }

        // The first row completes; the malformed middle row returns its exact
        // indexed error and the final row remains unvisited with default slots.
        let mut three_coordinates = CoordinateBlock::default();
        three_coordinates.conformers_3d.extend([
            Conformer3D::new(IDS[0], [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]].to_vec(), true),
            Conformer3D::new(IDS[1], [[0.0, 0.0, 0.0]].to_vec(), true),
            Conformer3D::new(IDS[2], [[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]].to_vec(), true),
        ]);
        three_coordinates.conformers_2d.push(
            cosmolkit_model::Conformer2D::new(IDS[0], TWO_D_SENTINEL.to_vec())
                .with_prop("layout", "preserve-this-overlapping-2d-row"),
        );
        let three_before = three_coordinates.clone();
        let mut three_results = vec![SENTINEL];
        let mut three_options = super::SingleConformerOptions::for_conformer(IDS[0]);
        three_options.max_iterations = 0;
        three_options.ignore_interfragment_interactions = false;
        super::super::api::reset_prepare_parameter_query_calls();
        super::uff_public_perf_serial_adapter_probe_start();
        let three_error = super::super::optimization::optimize_prepared_uff_serial(
            &topology,
            &mut three_coordinates,
            &mut three_results,
            &valence,
            &rings,
            &properties,
            &mut Vec::new(),
            three_options,
        )
        .expect_err("the middle conformer preserves its source coordinate error");
        actual_calls += 1;
        let three_probe = super::uff_public_perf_serial_adapter_probe_finish();
        assert!(matches!(
            three_error,
            UffPreparedOptimizationError::Serial(super::SerialUffOptimizationError::Optimization(
                super::SerialConformerOptimizationError::CoordinateCount {
                    input_index: 1,
                    source: super::SerialCoordinateCountError {
                        conformer_id: 7,
                        atoms: 2,
                        coordinates: 1
                    }
                }
            ))
        ));
        assert_eq!(super::super::api::prepare_parameter_query_calls(), 1);
        assert_eq!(three_probe, Some((2, Some(0))));
        assert_eq!(three_results.len(), 3);
        assert_eq!(three_results[0].status, 1);
        assert!((three_results[0].energy - source_energy(4.0)).abs() <= 1.0e-12);
        assert_eq!(
            &three_results[1..],
            &[super::OptimizationOutcome {
                status: 0,
                energy: 0.0,
            }; 2]
        );
        assert_eq!(three_coordinates.conformers_3d, three_before.conformers_3d);
        assert_eq!(three_coordinates.conformers_2d, three_before.conformers_2d);
        assert_eq!(
            three_coordinates.source_coordinate_dim,
            three_before.source_coordinate_dim
        );
        assert!(three_results.capacity() >= 3);
        assert_eq!(actual_calls, 4);
    }

    #[test]
    fn uff_block_l03_same_block_serial_w06_reference_matrix() {
        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; 2],
        };
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let mut owner_calls = 0_usize;

        for row_count in 1..=IDS.len() {
            for selected_id in IDS[..row_count].iter().copied() {
                for max_iterations in [0, 1] {
                    let mut coordinates = CoordinateBlock::default();
                    for row in 0..row_count {
                        coordinates.conformers_3d.push(Conformer3D::new(
                            IDS[row],
                            [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]].to_vec(),
                            true,
                        ));
                    }
                    coordinates.conformers_2d.push(
                        cosmolkit_model::Conformer2D::new(selected_id, TWO_D_SENTINEL.to_vec())
                            .with_prop("layout", "preserve-this-2d-row"),
                    );
                    let two_d_before = coordinates.conformers_2d[0].clone();
                    let two_d_bits_before = two_d_before
                        .coordinates()
                        .iter()
                        .flatten()
                        .map(|coordinate| coordinate.to_bits())
                        .collect::<Vec<_>>();
                    let block_dimension_before = coordinates.source_coordinate_dim;
                    let mut options = super::SingleConformerOptions::for_conformer(selected_id);
                    options.max_iterations = max_iterations;
                    options.vdw_threshold = 10.0;
                    options.ignore_interfragment_interactions = false;
                    let mut diagnostics = Vec::new();
                    let mut results = vec![SENTINEL; row_count + 1];

                    crate::kernel::cf3d_uff_one_kernel_counts_reset();
                    super::optimize_serial_uff_coordinate_block(
                        &topology,
                        &mut coordinates,
                        &mut results,
                        super::UffAtomStateRef::SuppliedRows {
                            total_valences: &total_valences,
                            conjugated_presence: &conjugated,
                        },
                        &rings,
                        &valence,
                        &properties,
                        &mut diagnostics,
                        options,
                    )
                    .unwrap_or_else(|error| {
                        panic!("fixed W06 same-block case succeeds: rows={row_count}, selected={selected_id}, iterations={max_iterations}, error={error:?}")
                    });
                    owner_calls += 1;

                    assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
                    assert!(diagnostics.is_empty());
                    assert_eq!(results.len(), row_count);
                    assert_eq!(
                        coordinates
                            .conformers_3d
                            .iter()
                            .map(Conformer3D::id)
                            .collect::<Vec<_>>(),
                        &IDS[..row_count]
                    );
                    for row in 0..row_count {
                        let (expected_x, expected_energy) = if max_iterations == 0 {
                            (
                                [0.0, TARGET_DISTANCES[row]],
                                source_energy(TARGET_DISTANCES[row]),
                            )
                        } else {
                            source_one_step(TARGET_DISTANCES[row])
                        };
                        for (point, expected_point) in
                            [[expected_x[0], 0.0, 0.0], [expected_x[1], 0.0, 0.0]]
                                .into_iter()
                                .enumerate()
                        {
                            for (axis, expected_component) in expected_point.into_iter().enumerate()
                            {
                                assert_eq!(
                                    coordinates.conformers_3d[row].coordinates()[point][axis]
                                        .to_bits(),
                                    expected_component.to_bits(),
                                    "W06 same-block coordinate row={row}, point={point}, axis={axis}, selected={selected_id}, iterations={max_iterations}"
                                );
                            }
                        }
                        assert_eq!(results[row].status, 1);
                        assert!(
                            (results[row].energy - expected_energy).abs() <= 1.0e-12,
                            "W06 same-block source energy row={row}, selected={selected_id}, iterations={max_iterations}"
                        );
                    }
                    assert_eq!(coordinates.conformers_2d, [two_d_before.clone()]);
                    assert_eq!(
                        coordinates.conformers_2d[0]
                            .coordinates()
                            .iter()
                            .flatten()
                            .map(|coordinate| coordinate.to_bits())
                            .collect::<Vec<_>>(),
                        two_d_bits_before
                    );
                    assert_eq!(coordinates.conformers_2d[0].props(), two_d_before.props());
                    assert_eq!(
                        coordinates.conformers_2d[0].id(),
                        selected_id,
                        "2D sentinel deliberately overlaps the selected 3D ID"
                    );
                    assert_eq!(coordinates.source_coordinate_dim, block_dimension_before);
                }
            }
        }
        assert_eq!(owner_calls, 12);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_block_l04_same_block_dispatch_w06_reference_matrix() {
        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const REQUESTED: [i32; 4] = [-2, 0, 1, 3];
        const HARDWARE: [u32; 3] = [0, 1, 3];
        const THREADSAFE_COUNTS: [[i32; 4]; 3] = [[1, 1, 1, 3], [1, 1, 1, 3], [1, 3, 1, 3]];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; 2],
        };
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let mut owner_calls = 0_usize;
        let mut underfilled_results = 0_usize;
        let mut oversized_results = 0_usize;

        for row_count in 1..=IDS.len() {
            for (selected_index, selected_id) in IDS[..row_count].iter().copied().enumerate() {
                for (requested_index, requested_threads) in REQUESTED.into_iter().enumerate() {
                    for (hardware_index, observed_hardware) in HARDWARE.into_iter().enumerate() {
                        for threadsafe in [false, true] {
                            for max_iterations in [0_i32, 1] {
                                let expected_source_count = if threadsafe {
                                    THREADSAFE_COUNTS[hardware_index][requested_index]
                                } else {
                                    1
                                };
                                let expected_workers = if expected_source_count == 1 {
                                    0
                                } else {
                                    expected_source_count as usize
                                };
                                let mut coordinates = CoordinateBlock::default();
                                for row in 0..row_count {
                                    coordinates.conformers_3d.push(Conformer3D::new(
                                        IDS[row],
                                        [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]]
                                            .to_vec(),
                                        true,
                                    ));
                                }
                                coordinates.conformers_2d.push(
                                    cosmolkit_model::Conformer2D::new(
                                        selected_id,
                                        TWO_D_SENTINEL.to_vec(),
                                    )
                                    .with_prop("layout", "preserve-this-2d-row"),
                                );
                                let two_d_before = coordinates.conformers_2d[0].clone();
                                let two_d_bits_before = two_d_before
                                    .coordinates()
                                    .iter()
                                    .flatten()
                                    .map(|coordinate| coordinate.to_bits())
                                    .collect::<Vec<_>>();
                                let block_dimension_before = coordinates.source_coordinate_dim;
                                let initial_result_len = if (owner_calls
                                    + selected_index
                                    + requested_index
                                    + hardware_index
                                    + usize::from(threadsafe)
                                    + max_iterations as usize)
                                    % 2
                                    == 0
                                {
                                    row_count + 1
                                } else {
                                    row_count.saturating_sub(1)
                                };
                                underfilled_results += usize::from(initial_result_len < row_count);
                                oversized_results += usize::from(initial_result_len > row_count);
                                let mut results = vec![SENTINEL; initial_result_len];
                                let mut options =
                                    super::SingleConformerOptions::for_conformer(selected_id);
                                options.max_iterations = max_iterations;
                                options.vdw_threshold = 10.0;
                                options.ignore_interfragment_interactions = false;
                                let mut diagnostics = Vec::new();

                                crate::kernel::cf3d_uff_one_kernel_counts_reset();
                                let outcome = super::optimize_dispatched_uff_coordinate_block(
                                    &topology,
                                    &mut coordinates,
                                    &mut results,
                                    super::UffAtomStateRef::SuppliedRows {
                                        total_valences: &total_valences,
                                        conjugated_presence: &conjugated,
                                    },
                                    &rings,
                                    &valence,
                                    &properties,
                                    &mut diagnostics,
                                    options,
                                    requested_threads,
                                    observed_hardware,
                                    threadsafe,
                                )
                                .unwrap_or_else(|error| {
                                    panic!("fixed W06 same-block dispatch succeeds: rows={row_count}, selected={selected_id}, requested={requested_threads}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}, error={error:?}")
                                });
                                owner_calls += 1;

                                match outcome {
                                    super::PreparedConformerDispatchOutcome::Serial(result) => {
                                        assert_eq!(expected_source_count, 1);
                                        result.unwrap_or_else(|error| {
                                            panic!("fixed W06 same-block serial outcome: {error:?}")
                                        });
                                    }
                                    super::PreparedConformerDispatchOutcome::Workers(result) => {
                                        assert_eq!(expected_source_count, 3);
                                        let raw_joins = result.unwrap_or_else(|error| {
                                            panic!("fixed W06 same-block worker setup: {error:?}")
                                        });
                                        assert_eq!(raw_joins.len(), expected_workers);
                                        for (lane, raw_join) in raw_joins.iter().enumerate() {
                                            assert!(
                                                matches!(raw_join, Ok(Ok(()))),
                                                "ordered raw W06 join lane={lane}, rows={row_count}, selected={selected_id}, requested={requested_threads}, hardware={observed_hardware}, iterations={max_iterations}"
                                            );
                                        }
                                    }
                                }
                                assert_eq!(
                                    crate::kernel::cf3d_uff_one_kernel_counts(),
                                    (1 + expected_workers, 1, expected_workers),
                                    "one builder field plus exactly one existing copy per worker"
                                );
                                assert!(diagnostics.is_empty());
                                assert_eq!(results.len(), row_count);
                                assert_eq!(
                                    coordinates
                                        .conformers_3d
                                        .iter()
                                        .map(Conformer3D::id)
                                        .collect::<Vec<_>>(),
                                    &IDS[..row_count]
                                );
                                for row in 0..row_count {
                                    let (expected_x, expected_energy) = if max_iterations == 0 {
                                        (
                                            [0.0, TARGET_DISTANCES[row]],
                                            source_energy(TARGET_DISTANCES[row]),
                                        )
                                    } else {
                                        source_one_step(TARGET_DISTANCES[row])
                                    };
                                    for (point, expected_point) in
                                        [[expected_x[0], 0.0, 0.0], [expected_x[1], 0.0, 0.0]]
                                            .into_iter()
                                            .enumerate()
                                    {
                                        for (axis, expected_component) in
                                            expected_point.into_iter().enumerate()
                                        {
                                            assert_eq!(
                                                coordinates.conformers_3d[row].coordinates()[point]
                                                    [axis]
                                                    .to_bits(),
                                                expected_component.to_bits(),
                                                "W06 same-block dispatch coordinate row={row}, point={point}, axis={axis}, selected={selected_id}, requested={requested_threads}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                            );
                                        }
                                    }
                                    assert_eq!(results[row].status, 1);
                                    assert!(
                                        (results[row].energy - expected_energy).abs() <= 1.0e-12,
                                        "W06 same-block dispatch source energy row={row}, selected={selected_id}, requested={requested_threads}, hardware={observed_hardware}, threadsafe={threadsafe}, iterations={max_iterations}"
                                    );
                                }
                                assert_eq!(coordinates.conformers_2d, [two_d_before.clone()]);
                                assert_eq!(
                                    coordinates.conformers_2d[0]
                                        .coordinates()
                                        .iter()
                                        .flatten()
                                        .map(|coordinate| coordinate.to_bits())
                                        .collect::<Vec<_>>(),
                                    two_d_bits_before
                                );
                                assert_eq!(
                                    coordinates.conformers_2d[0].props(),
                                    two_d_before.props()
                                );
                                assert_eq!(coordinates.conformers_2d[0].id(), selected_id);
                                assert_eq!(
                                    coordinates.source_coordinate_dim,
                                    block_dimension_before
                                );
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(owner_calls, 288);
        assert!(underfilled_results > 0);
        assert!(oversized_results > 0);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_block_l06_same_block_storage_and_copy_matrix() {
        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const MODES: [(&str, i32, u32, bool, usize); 3] = [
            ("serial", 1, 3, false, 0),
            ("native1", 1, 3, true, 0),
            ("native3", 3, 3, true, 3),
        ];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; 2],
        };
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let mut owner_calls = 0_usize;

        for row_count in [1, 3] {
            for &(mode, requested_threads, observed_hardware, threadsafe, worker_count) in &MODES {
                for max_iterations in [0_i32, 1] {
                    let mut coordinates = CoordinateBlock::default();
                    for row in 0..row_count {
                        coordinates.conformers_3d.push(Conformer3D::new(
                            IDS[row],
                            [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]].to_vec(),
                            true,
                        ));
                    }
                    coordinates.conformers_2d.push(
                        cosmolkit_model::Conformer2D::new(IDS[0], TWO_D_SENTINEL.to_vec())
                            .with_prop("layout", "preserve-this-2d-row"),
                    );
                    let coordinate_row_addresses = coordinates
                        .conformers_3d
                        .iter()
                        .map(|conformer| conformer.coordinates().as_ptr() as usize)
                        .collect::<Vec<_>>();
                    let two_d_before = coordinates.conformers_2d[0].clone();
                    let two_d_bits_before = two_d_before
                        .coordinates()
                        .iter()
                        .flatten()
                        .map(|coordinate| coordinate.to_bits())
                        .collect::<Vec<_>>();
                    let block_dimension_before = coordinates.source_coordinate_dim;
                    let mut options = super::SingleConformerOptions::for_conformer(IDS[0]);
                    options.max_iterations = max_iterations;
                    options.vdw_threshold = 10.0;
                    options.ignore_interfragment_interactions = false;
                    let mut diagnostics = Vec::new();
                    let mut results = vec![SENTINEL; row_count + 1];

                    crate::kernel::cf3d_uff_one_kernel_counts_reset();
                    let outcome = super::optimize_dispatched_uff_coordinate_block(
                        &topology,
                        &mut coordinates,
                        &mut results,
                        super::UffAtomStateRef::SuppliedRows {
                            total_valences: &total_valences,
                            conjugated_presence: &conjugated,
                        },
                        &rings,
                        &valence,
                        &properties,
                        &mut diagnostics,
                        options,
                        requested_threads,
                        observed_hardware,
                        threadsafe,
                    )
                    .unwrap_or_else(|error| {
                        panic!("fixed W06 L06 same-block case succeeds: rows={row_count}, mode={mode}, iterations={max_iterations}, error={error:?}")
                    });
                    owner_calls += 1;

                    match outcome {
                        super::PreparedConformerDispatchOutcome::Serial(result) => {
                            assert_eq!(worker_count, 0, "L06 mode={mode}");
                            result.unwrap_or_else(|error| {
                                panic!("fixed W06 L06 serial outcome mode={mode}: {error:?}")
                            });
                        }
                        super::PreparedConformerDispatchOutcome::Workers(result) => {
                            assert_eq!(worker_count, 3, "L06 mode={mode}");
                            let raw_joins = result.unwrap_or_else(|error| {
                                panic!("fixed W06 L06 worker setup mode={mode}: {error:?}")
                            });
                            assert_eq!(raw_joins.len(), worker_count, "L06 mode={mode}");
                            for (lane, raw_join) in raw_joins.iter().enumerate() {
                                assert!(
                                    matches!(raw_join, Ok(Ok(()))),
                                    "ordered L06 raw join lane={lane}, rows={row_count}, mode={mode}, iterations={max_iterations}"
                                );
                            }
                        }
                    }

                    assert_eq!(
                        crate::kernel::cf3d_uff_one_kernel_counts(),
                        (1 + worker_count, 1, worker_count),
                        "one builder field and source-required per-worker field/term copies, including empty lanes; mode={mode}"
                    );
                    if worker_count == 0 {
                        let work = crate::kernel::cf3d_uff_one_serial_work_counts();
                        assert_eq!(work.position_buffer_growths, 0, "L06 mode={mode}");
                        assert_eq!(work.initialize_calls, row_count, "L06 mode={mode}");
                        assert_eq!(
                            work.initial_position_buffer_address,
                            work.final_position_buffer_address,
                            "released field retains and reuses its position allocation; mode={mode}"
                        );
                        assert_eq!(
                            work.initial_position_buffer_capacity,
                            work.final_position_buffer_capacity,
                            "released field retains its position capacity; mode={mode}"
                        );
                        assert!(
                            work.initial_position_buffer_capacity >= topology.atoms.len(),
                            "retained construction allocation avoids a first-row reserve; mode={mode}"
                        );
                    }
                    assert!(diagnostics.is_empty());
                    assert_eq!(results.len(), row_count);
                    assert_eq!(
                        coordinates
                            .conformers_3d
                            .iter()
                            .map(Conformer3D::id)
                            .collect::<Vec<_>>(),
                        &IDS[..row_count]
                    );
                    assert_eq!(
                        coordinates
                            .conformers_3d
                            .iter()
                            .map(|conformer| conformer.coordinates().as_ptr() as usize)
                            .collect::<Vec<_>>(),
                        coordinate_row_addresses,
                        "the original coordinate rows are optimized in place; mode={mode}"
                    );
                    for row in 0..row_count {
                        let (expected_x, expected_energy) = if max_iterations == 0 {
                            (
                                [0.0, TARGET_DISTANCES[row]],
                                source_energy(TARGET_DISTANCES[row]),
                            )
                        } else {
                            source_one_step(TARGET_DISTANCES[row])
                        };
                        for (point, expected_point) in
                            [[expected_x[0], 0.0, 0.0], [expected_x[1], 0.0, 0.0]]
                                .into_iter()
                                .enumerate()
                        {
                            for (axis, expected_component) in expected_point.into_iter().enumerate()
                            {
                                assert_eq!(
                                    coordinates.conformers_3d[row].coordinates()[point][axis]
                                        .to_bits(),
                                    expected_component.to_bits(),
                                    "W06 L06 same-block coordinate row={row}, point={point}, axis={axis}, mode={mode}, iterations={max_iterations}"
                                );
                            }
                        }
                        assert_eq!(results[row].status, 1);
                        assert!(
                            (results[row].energy - expected_energy).abs() <= 1.0e-12,
                            "W06 L06 source energy row={row}, mode={mode}, iterations={max_iterations}"
                        );
                    }
                    assert_eq!(coordinates.conformers_2d, [two_d_before.clone()]);
                    assert_eq!(
                        coordinates.conformers_2d[0]
                            .coordinates()
                            .iter()
                            .flatten()
                            .map(|coordinate| coordinate.to_bits())
                            .collect::<Vec<_>>(),
                        two_d_bits_before
                    );
                    assert_eq!(coordinates.conformers_2d[0].props(), two_d_before.props());
                    assert_eq!(coordinates.conformers_2d[0].id(), IDS[0]);
                    assert_eq!(coordinates.source_coordinate_dim, block_dimension_before);
                }
            }
        }
        assert_eq!(owner_calls, 12);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_block_l05_same_block_construction_and_thread_error_order() {
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const PREFIX: [(i32, u64); 5] = [
            (17, 0x3ff8_0000_0000_0000),
            (-8, 0xc004_0000_0000_0000),
            (5, 0x4009_21fb_5444_2d18),
            (-13, 0x0000_0000_0000_0001),
            (29, 0xbff0_0000_0000_0000),
        ];
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; 2],
        };
        let properties = MoleculeProperties::default();
        let make_results = |len: usize| {
            PREFIX[..len]
                .iter()
                .map(|&(status, energy_bits)| super::OptimizationOutcome {
                    status,
                    energy: f64::from_bits(energy_bits),
                })
                .collect::<Vec<_>>()
        };
        let result_bits = |results: &[super::OptimizationOutcome]| {
            results
                .iter()
                .map(|result| (result.status, result.energy.to_bits()))
                .collect::<Vec<_>>()
        };
        let coordinate_state = |coordinates: &CoordinateBlock| {
            let three_d = coordinates
                .conformers_3d
                .iter()
                .map(|conformer| {
                    (
                        conformer.id(),
                        conformer
                            .coordinates()
                            .iter()
                            .flatten()
                            .map(|coordinate| coordinate.to_bits())
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<Vec<_>>();
            let two_d = coordinates
                .conformers_2d
                .iter()
                .map(|conformer| {
                    (
                        conformer.id(),
                        conformer
                            .coordinates()
                            .iter()
                            .flatten()
                            .map(|coordinate| coordinate.to_bits())
                            .collect::<Vec<_>>(),
                        conformer.props().clone(),
                    )
                })
                .collect::<Vec<_>>();
            (three_d, two_d, coordinates.source_coordinate_dim)
        };
        let make_coordinates = |with_3d_rows: bool, two_d_id: usize| {
            let mut coordinates = CoordinateBlock::default();
            if with_3d_rows {
                for row in 0..IDS.len() {
                    coordinates.conformers_3d.push(Conformer3D::new(
                        IDS[row],
                        [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]].to_vec(),
                        true,
                    ));
                }
            }
            coordinates.conformers_2d.push(
                cosmolkit_model::Conformer2D::new(two_d_id, TWO_D_SENTINEL.to_vec())
                    .with_prop("layout", "preserve-this-2d-row"),
            );
            coordinates
        };
        let mut owner_calls = 0_usize;

        // A missing selected 3D ID (including a 2D-only block) fails in the
        // existing builder before the dispatcher can resize caller results.
        for with_3d_rows in [true, false] {
            for initial_len in [0_usize, 1, 5] {
                let mut coordinates = make_coordinates(with_3d_rows, 99);
                let coordinates_before = coordinate_state(&coordinates);
                let mut results = make_results(initial_len);
                let results_before = result_bits(&results);
                let mut options = super::SingleConformerOptions::for_conformer(99);
                options.ignore_interfragment_interactions = false;
                let mut diagnostics = Vec::new();

                crate::kernel::cf3d_uff_one_kernel_counts_reset();
                let error = match super::optimize_dispatched_uff_coordinate_block(
                    &topology,
                    &mut coordinates,
                    &mut results,
                    super::UffAtomStateRef::SuppliedRows {
                        total_valences: &total_valences,
                        conjugated_presence: &conjugated,
                    },
                    &rings,
                    &valence,
                    &properties,
                    &mut diagnostics,
                    options,
                    i32::MIN,
                    0,
                    true,
                ) {
                    Err(error) => error,
                    Ok(_) => panic!("missing selected 3D ID must fail during construction"),
                };
                owner_calls += 1;
                assert!(matches!(
                    error,
                    super::DispatchedUffOptimizationError::Construction(
                        super::AutomaticForceFieldConstructionError::Construction(
                            ForceFieldConstructionError::Builder(
                                UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                                    conformer_id: 99
                                }
                            )
                        )
                    )
                ));
                assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 0, 0));
                assert_eq!(result_bits(&results), results_before);
                assert_eq!(coordinate_state(&coordinates), coordinates_before);
                assert!(diagnostics.is_empty());
            }
        }

        // Construction succeeds, then source D04 resizes before its typed
        // INT_MIN thread-count error; no row stage may change coordinates.
        for initial_len in [0_usize, 1, 5] {
            let mut coordinates = make_coordinates(true, 7);
            let coordinates_before = coordinate_state(&coordinates);
            let mut results = make_results(initial_len);
            let mut options = super::SingleConformerOptions::for_conformer(7);
            options.max_iterations = 0;
            options.ignore_interfragment_interactions = false;
            let mut diagnostics = Vec::new();

            crate::kernel::cf3d_uff_one_kernel_counts_reset();
            let error = match super::optimize_dispatched_uff_coordinate_block(
                &topology,
                &mut coordinates,
                &mut results,
                super::UffAtomStateRef::SuppliedRows {
                    total_valences: &total_valences,
                    conjugated_presence: &conjugated,
                },
                &rings,
                &valence,
                &properties,
                &mut diagnostics,
                options,
                i32::MIN,
                0,
                true,
            ) {
                Err(error) => error,
                Ok(_) => panic!("threadsafe INT_MIN remains a typed source boundary"),
            };
            owner_calls += 1;
            assert!(matches!(
                error,
                super::DispatchedUffOptimizationError::ThreadCount(
                    super::UffThreadCountError::UndefinedSignedNegation
                )
            ));
            assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
            let resized_expectation = (0..IDS.len())
                .map(|row| {
                    if row < initial_len {
                        PREFIX[row]
                    } else {
                        (0, 0.0_f64.to_bits())
                    }
                })
                .collect::<Vec<_>>();
            assert_eq!(results.len(), IDS.len());
            assert_eq!(result_bits(&results), resized_expectation);
            assert_eq!(coordinate_state(&coordinates), coordinates_before);
            assert!(diagnostics.is_empty());
        }

        // The pinned non-threadsafe route returns one before negating the
        // requested count, so INT_MIN reaches the existing serial owner.
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        for initial_len in [0_usize, 1, 5] {
            let mut coordinates = make_coordinates(true, 7);
            let coordinates_before = coordinate_state(&coordinates);
            let mut results = make_results(initial_len);
            let mut options = super::SingleConformerOptions::for_conformer(7);
            options.max_iterations = 0;
            options.ignore_interfragment_interactions = false;
            let mut diagnostics = Vec::new();

            crate::kernel::cf3d_uff_one_kernel_counts_reset();
            let outcome = super::optimize_dispatched_uff_coordinate_block(
                &topology,
                &mut coordinates,
                &mut results,
                super::UffAtomStateRef::SuppliedRows {
                    total_valences: &total_valences,
                    conjugated_presence: &conjugated,
                },
                &rings,
                &valence,
                &properties,
                &mut diagnostics,
                options,
                i32::MIN,
                0,
                false,
            )
            .unwrap_or_else(|error| panic!("source non-threadsafe INT_MIN succeeds: {error:?}"));
            owner_calls += 1;
            match outcome {
                super::PreparedConformerDispatchOutcome::Serial(Ok(())) => {}
                _ => panic!("source non-threadsafe INT_MIN must use serial optimization"),
            }
            assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
            assert_eq!(results.len(), IDS.len());
            for row in 0..IDS.len() {
                assert_eq!(results[row].status, 1);
                assert!(
                    (results[row].energy - source_energy(TARGET_DISTANCES[row])).abs() <= 1.0e-12
                );
            }
            assert_eq!(coordinate_state(&coordinates), coordinates_before);
            assert!(diagnostics.is_empty());
        }

        // The full builder rejects a malformed selected atom-row count before
        // term assembly; retain its exact typed source failure and pre-dispatch
        // results for every caller-vector length.
        for initial_len in [0_usize, 1, 5] {
            let mut coordinates = CoordinateBlock::default();
            coordinates
                .conformers_3d
                .push(Conformer3D::new(7, vec![[0.0, 0.0, 0.0]], true));
            coordinates.conformers_2d.push(
                cosmolkit_model::Conformer2D::new(7, TWO_D_SENTINEL.to_vec())
                    .with_prop("layout", "preserve-this-2d-row"),
            );
            let coordinates_before = coordinate_state(&coordinates);
            let mut results = make_results(initial_len);
            let results_before = result_bits(&results);
            let mut options = super::SingleConformerOptions::for_conformer(7);
            options.ignore_interfragment_interactions = false;
            let mut diagnostics = Vec::new();

            crate::kernel::cf3d_uff_one_kernel_counts_reset();
            let error = match super::optimize_dispatched_uff_coordinate_block(
                &topology,
                &mut coordinates,
                &mut results,
                super::UffAtomStateRef::SuppliedRows {
                    total_valences: &total_valences,
                    conjugated_presence: &conjugated,
                },
                &rings,
                &valence,
                &properties,
                &mut diagnostics,
                options,
                i32::MIN,
                0,
                true,
            ) {
                Err(error) => error,
                Ok(_) => panic!("the full builder rejects the selected short row"),
            };
            owner_calls += 1;
            assert!(matches!(
                error,
                super::DispatchedUffOptimizationError::Construction(
                    super::AutomaticForceFieldConstructionError::Construction(
                        ForceFieldConstructionError::Builder(
                            UffBuilderError::SelectedConformerCoordinateCountMismatch {
                                conformer_id: 7,
                                atoms: 2,
                                coordinates: 1,
                            }
                        )
                    )
                )
            ));
            assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 0, 0));
            assert_eq!(result_bits(&results), results_before);
            assert_eq!(coordinate_state(&coordinates), coordinates_before);
            assert!(diagnostics.is_empty());
        }

        assert_eq!(owner_calls, 15);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_dispatch_d06_full_uff_builder_dispatch_matches_references_and_failures() {
        use crate::uff::atom_typer::{UffTypingDiagnostic, UffTypingDiagnosticKind};
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const INPUT_IDS: [usize; 3] = [80, 7, 900];
        const REFERENCE_12: [[f64; 3]; 2] = [[20.0, 0.0, 0.0], [22.0, 0.0, 0.0]];
        const REFERENCE_31: [[f64; 3]; 2] = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        const SENTINEL: super::OptimizationOutcome = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        const REQUESTED_THREADS: [i32; 3] = [1, 2, 3];
        const MAX_ITERATIONS: [i32; 2] = [0, 1];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; 2],
        };
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let source_one_step = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio7 = (ratio3 * ratio3) * ratio;
            let ratio13 = (ratio6 * ratio6) * ratio;
            let pre_factor = 12.0 * well_depth / minimum * (ratio7 - ratio13);
            let first_gradient = (pre_factor * (0.0 - distance) / distance) * 0.1;
            let second_gradient = (pre_factor * (distance - 0.0) / distance) * 0.1;
            let first_x = 0.0 + 1.0 * -first_gradient;
            let second_x = distance + 1.0 * -second_gradient;
            let dx = first_x - second_x;
            let after_distance = (dx * dx).sqrt();
            ([first_x, second_x], source_energy(after_distance))
        };
        let make_reference = || {
            let mut coordinates = CoordinateBlock::default();
            coordinates
                .conformers_3d
                .push(Conformer3D::new(12, REFERENCE_12.to_vec(), true));
            coordinates
                .conformers_3d
                .push(Conformer3D::new(31, REFERENCE_31.to_vec(), true));
            coordinates
        };
        let mut owner_calls = 0_usize;

        for source_mode in [false, true] {
            for requested_threads in REQUESTED_THREADS {
                for row_count in [0_usize, 1, 3] {
                    for max_iterations in MAX_ITERATIONS {
                        let mut reference_coordinates = make_reference();
                        let mut coordinate_rows = TARGET_DISTANCES[..row_count]
                            .iter()
                            .copied()
                            .map(|distance| vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
                            .collect::<Vec<_>>();
                        let mut conformers = coordinate_rows
                            .iter_mut()
                            .enumerate()
                            .map(|(row, positions)| super::SerialConformer {
                                id: INPUT_IDS[row],
                                positions: positions.as_mut_slice(),
                            })
                            .collect::<Vec<_>>();
                        let mut results = vec![SENTINEL; row_count + 1];
                        let mut diagnostics = Vec::new();
                        let mut options = super::SingleConformerOptions::for_conformer(31);
                        options.max_iterations = max_iterations;
                        options.vdw_threshold = 10.0;
                        options.ignore_interfragment_interactions = false;
                        let expected_workers = if source_mode && requested_threads > 1 {
                            requested_threads as usize
                        } else {
                            0
                        };

                        crate::kernel::cf3d_uff_one_kernel_counts_reset();
                        let outcome = super::optimize_dispatched_uff(
                            &topology,
                            &mut reference_coordinates,
                            &mut conformers,
                            &mut results,
                            &total_valences,
                            &conjugated,
                            &rings,
                            &valence,
                            &properties,
                            &mut diagnostics,
                            options,
                            requested_threads,
                            3,
                            source_mode,
                        )
                        .unwrap_or_else(|error| {
                            panic!(
                                "fixed W06 UFF builder-dispatch case succeeds: mode={source_mode}, requested={requested_threads}, rows={row_count}, iterations={max_iterations}, error={error:?}"
                            )
                        });
                        owner_calls += 1;

                        match outcome {
                            super::PreparedConformerDispatchOutcome::Serial(result) => {
                                assert_eq!(expected_workers, 0);
                                result.unwrap_or_else(|error| {
                                    panic!("fixed W06 serial case succeeds: {error:?}")
                                });
                            }
                            super::PreparedConformerDispatchOutcome::Workers(result) => {
                                assert_ne!(expected_workers, 0);
                                let joins = result.unwrap_or_else(|error| {
                                    panic!("fixed W06 worker setup succeeds: {error:?}")
                                });
                                assert_eq!(joins.len(), expected_workers);
                                for (lane, join) in joins.iter().enumerate() {
                                    assert!(
                                        matches!(join, Ok(Ok(()))),
                                        "W06 lane={lane}, mode={source_mode}, requested={requested_threads}, rows={row_count}, iterations={max_iterations}"
                                    );
                                }
                            }
                        }
                        assert_eq!(
                            crate::kernel::cf3d_uff_one_kernel_counts(),
                            (1 + expected_workers, 1, expected_workers),
                            "one automatic field build plus one field/term copy per source worker"
                        );
                        assert!(
                            diagnostics.is_empty(),
                            "valid Na/Cl fixture emits no diagnostics"
                        );
                        assert_eq!(
                            conformers
                                .iter()
                                .map(|conformer| conformer.id)
                                .collect::<Vec<_>>(),
                            INPUT_IDS[..row_count]
                        );
                        drop(conformers);

                        assert_eq!(results.len(), row_count);
                        for row in 0..row_count {
                            let (expected_x, expected_energy) = if max_iterations == 0 {
                                (
                                    [0.0, TARGET_DISTANCES[row]],
                                    source_energy(TARGET_DISTANCES[row]),
                                )
                            } else {
                                source_one_step(TARGET_DISTANCES[row])
                            };
                            for (point, expected_point) in
                                [[expected_x[0], 0.0, 0.0], [expected_x[1], 0.0, 0.0]]
                                    .into_iter()
                                    .enumerate()
                            {
                                for (axis, expected_component) in
                                    expected_point.into_iter().enumerate()
                                {
                                    assert_eq!(
                                        coordinate_rows[row][point][axis].to_bits(),
                                        expected_component.to_bits(),
                                        "W06 UFF builder-dispatch coordinate row={row}, point={point}, axis={axis}, mode={source_mode}, requested={requested_threads}, iterations={max_iterations}"
                                    );
                                }
                            }
                            assert_eq!(results[row].status, 1);
                            assert!(
                                (results[row].energy - expected_energy).abs() <= 1.0e-12,
                                "W06 UFF builder-dispatch source energy row={row}, mode={source_mode}, requested={requested_threads}, iterations={max_iterations}"
                            );
                        }
                        assert_eq!(
                            reference_coordinates
                                .conformers_3d
                                .iter()
                                .map(Conformer3D::id)
                                .collect::<Vec<_>>(),
                            [12, 31]
                        );
                        assert_eq!(
                            reference_coordinates.conformers_3d[0].coordinates(),
                            REFERENCE_12
                        );
                        assert_eq!(
                            reference_coordinates.conformers_3d[1].coordinates(),
                            REFERENCE_31
                        );
                    }
                }
            }
        }
        assert_eq!(owner_calls, 36);

        // The retained U10 TBP fixture's nonfatal typer warning is emitted
        // before a later missing-reference construction error, while D04
        // remains unentered and cannot resize or touch the result/row inputs.
        const TBP_POINTS: [[f64; 3]; 6] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-0.6, 0.8, 0.0],
            [0.0, 0.6, 0.8],
            [0.0, -0.8, 0.6],
            [-0.8, 0.0, -0.6],
        ];
        let mut tbp_atoms = vec![worker_w06_atom(0, 15, Hybridization::Sp3d, 0)];
        tbp_atoms.extend((1..6).map(|row| worker_w06_atom(row, 6, Hybridization::Sp3, 0)));
        let tbp_topology = worker_w06_topology(
            tbp_atoms,
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
                (0, 5, BondOrder::Single),
            ],
        );
        let tbp_total_valences = [5, 1, 1, 1, 1, 1];
        let tbp_conjugated = [false; 6];
        let tbp_rings =
            fast_find_rings(&tbp_topology).expect("fixed TBP topology has ring information");
        let tbp_valence = ValenceAssignment {
            explicit_valence: tbp_total_valences.to_vec(),
            implicit_hydrogens: vec![0; 6],
        };
        let mut tbp_reference = CoordinateBlock::default();
        tbp_reference
            .conformers_3d
            .push(Conformer3D::new(31, TBP_POINTS.to_vec(), true));
        let mut tbp_rows = vec![
            TBP_POINTS.to_vec(),
            TBP_POINTS.to_vec(),
            TBP_POINTS.to_vec(),
        ];
        let tbp_before = tbp_rows.clone();
        let mut tbp_conformers = tbp_rows
            .iter_mut()
            .enumerate()
            .map(|(row, positions)| super::SerialConformer {
                id: INPUT_IDS[row],
                positions: positions.as_mut_slice(),
            })
            .collect::<Vec<_>>();
        let mut sentinel_results = vec![SENTINEL];
        let mut tbp_diagnostics = Vec::new();
        let mut missing_reference_options = super::SingleConformerOptions::for_conformer(99);
        missing_reference_options.ignore_interfragment_interactions = false;
        crate::kernel::cf3d_uff_one_kernel_counts_reset();
        let construction_error = match super::optimize_dispatched_uff(
            &tbp_topology,
            &mut tbp_reference,
            &mut tbp_conformers,
            &mut sentinel_results,
            &tbp_total_valences,
            &tbp_conjugated,
            &tbp_rings,
            &tbp_valence,
            &properties,
            &mut tbp_diagnostics,
            missing_reference_options,
            i32::MIN,
            0,
            true,
        ) {
            Err(error) => error,
            Ok(_) => panic!("missing selected reference fails during construction"),
        };
        assert!(matches!(
            &construction_error,
            super::DispatchedUffOptimizationError::Construction(
                super::AutomaticForceFieldConstructionError::Construction(
                    ForceFieldConstructionError::Builder(
                        UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                            conformer_id: 99
                        }
                    )
                )
            )
        ));
        let construction_source = std::error::Error::source(&construction_error)
            .expect("construction variant remains the standard error source");
        assert!(construction_source.is::<super::AutomaticForceFieldConstructionError>());
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts().0, 1);
        assert_eq!(
            tbp_diagnostics,
            vec![UffTypingDiagnostic {
                atom_id: Some(cosmolkit_model::AtomId::new(0)),
                kind: UffTypingDiagnosticKind::Warning,
                message_prefix: "UFFTYPER: Warning: hybridization set to SP3 for atom ",
            }]
        );
        assert_eq!(sentinel_results, [SENTINEL]);
        drop(tbp_conformers);
        assert_eq!(tbp_rows, tbp_before);
        assert_eq!(tbp_reference.conformers_3d[0].id(), 31);
        assert_eq!(tbp_reference.conformers_3d[0].coordinates(), TBP_POINTS);

        // A successful build reaches D04: its resize happens before the
        // source-undefined INT_MIN resolver error, preserving the prefix and
        // source default pairs without entering any row stage.
        let mut reference_coordinates = make_reference();
        let mut coordinate_rows = TARGET_DISTANCES
            .iter()
            .copied()
            .map(|distance| vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
            .collect::<Vec<_>>();
        let before = coordinate_rows.clone();
        let mut conformers = coordinate_rows
            .iter_mut()
            .enumerate()
            .map(|(row, positions)| super::SerialConformer {
                id: INPUT_IDS[row],
                positions: positions.as_mut_slice(),
            })
            .collect::<Vec<_>>();
        let mut results = vec![SENTINEL];
        let mut diagnostics = Vec::new();
        let mut thread_options = super::SingleConformerOptions::for_conformer(31);
        thread_options.ignore_interfragment_interactions = false;
        crate::kernel::cf3d_uff_one_kernel_counts_reset();
        let thread_error = match super::optimize_dispatched_uff(
            &topology,
            &mut reference_coordinates,
            &mut conformers,
            &mut results,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            thread_options,
            i32::MIN,
            0,
            true,
        ) {
            Err(error) => error,
            Ok(_) => panic!("threadsafe INT_MIN is a typed source boundary"),
        };
        assert!(matches!(
            &thread_error,
            super::DispatchedUffOptimizationError::ThreadCount(
                super::UffThreadCountError::UndefinedSignedNegation
            )
        ));
        let thread_source = std::error::Error::source(&thread_error)
            .expect("thread-count variant remains the standard error source");
        assert!(
            thread_source
                .downcast_ref::<super::UffThreadCountError>()
                .is_some()
        );
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
        assert_eq!(results.len(), TARGET_DISTANCES.len());
        assert_eq!(results[0], SENTINEL);
        assert_eq!(results[1].status, 0);
        assert_eq!(results[1].energy.to_bits(), 0.0_f64.to_bits());
        assert_eq!(results[2].status, 0);
        assert_eq!(results[2].energy.to_bits(), 0.0_f64.to_bits());
        assert!(diagnostics.is_empty());
        assert_eq!(
            conformers
                .iter()
                .map(|conformer| conformer.id)
                .collect::<Vec<_>>(),
            INPUT_IDS
        );
        drop(conformers);
        assert_eq!(coordinate_rows, before);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_mt_m09_real_transport_preserves_typed_lane_failures_and_effects() {
        const IDS: [usize; 9] = [80, 7, 900, 31, 502, 8, 640, 91, 304];
        const SOURCE_X: [f64; 9] = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let coordinate_bits = |rows: &[Vec<[f64; 3]>]| {
            rows.iter()
                .flat_map(|row| row.iter())
                .flat_map(|point| point.iter())
                .map(|coordinate| coordinate.to_bits())
                .collect::<Vec<_>>()
        };
        let mut owner_calls = 0;

        // Capacity is checked before partition borrows or per-lane copies.
        {
            let mut coordinate_rows = SOURCE_X[..3]
                .iter()
                .copied()
                .map(|x| vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]])
                .collect::<Vec<_>>();
            let before = coordinate_bits(&coordinate_rows);
            let mut conformers = coordinate_rows
                .iter_mut()
                .enumerate()
                .map(|(row, positions)| super::SerialConformer {
                    id: IDS[row],
                    positions: positions.as_mut_slice(),
                })
                .collect::<Vec<_>>();
            let mut results = vec![sentinel; 2];
            let results_before = results.clone();
            let mut source_field = ForceField::new(3);
            source_field.add_contribution(Box::new(CoordinateNormContribution));
            cf3d_uff_one_kernel_counts_reset();

            let result = super::optimize_conformers_mt(
                &source_field,
                &mut conformers,
                &mut results,
                3,
                2,
                0,
            );
            owner_calls += 1;
            assert!(matches!(
                result,
                Err(super::SerialConformerOptimizationError::ResultCapacity(
                    super::SerialResultCapacityError {
                        conformers: 3,
                        result_slots: 2,
                    }
                ))
            ));
            drop(conformers);
            assert_eq!(coordinate_bits(&coordinate_rows), before);
            assert_eq!(results, results_before);
            assert_eq!(cf3d_uff_one_kernel_counts(), (0, 0, 0));
        }

        // Each lane's middle row is malformed. Its predecessor succeeds,
        // its successor remains untouched, and the other lanes finish.
        for failed_lane in 0..3_usize {
            let row_count = SOURCE_X.len();
            let failed_index = failed_lane + 3;
            let mut coordinate_rows = SOURCE_X
                .iter()
                .enumerate()
                .map(|(row, x)| {
                    if row == failed_index {
                        vec![[0.0, 0.0, 0.0]]
                    } else {
                        vec![[0.0, 0.0, 0.0], [*x, 0.0, 0.0]]
                    }
                })
                .collect::<Vec<_>>();
            let before = coordinate_bits(&coordinate_rows);
            let before_lengths = coordinate_rows.iter().map(Vec::len).collect::<Vec<_>>();
            let mut conformers = coordinate_rows
                .iter_mut()
                .enumerate()
                .map(|(row, positions)| super::SerialConformer {
                    id: IDS[row],
                    positions: positions.as_mut_slice(),
                })
                .collect::<Vec<_>>();
            let mut results = vec![sentinel; row_count + 1];
            let mut expected_results = vec![sentinel; row_count + 1];
            for lane in 0..3_usize {
                for source_index in (lane..row_count).step_by(3) {
                    if lane == failed_lane && source_index >= failed_index {
                        break;
                    }
                    expected_results[source_index] = OptimizationOutcome {
                        status: 1,
                        energy: SOURCE_X[source_index] * SOURCE_X[source_index],
                    };
                }
            }
            let mut source_field = ForceField::new(3);
            source_field.add_contribution(Box::new(CoordinateNormContribution));
            cf3d_uff_one_kernel_counts_reset();

            let outcomes = super::optimize_conformers_mt(
                &source_field,
                &mut conformers,
                &mut results,
                3,
                2,
                0,
            )
            .expect("the fixed malformed-row matrix has sufficient result capacity");
            owner_calls += 1;
            assert_eq!(outcomes.len(), 3);
            for lane in 0..3_usize {
                if lane == failed_lane {
                    match &outcomes[lane] {
                        Ok(Err(super::SerialConformerOptimizationError::CoordinateCount {
                            input_index,
                            source,
                        })) => {
                            assert_eq!(*input_index, failed_index);
                            assert_eq!(
                                *source,
                                super::SerialCoordinateCountError {
                                    conformer_id: IDS[failed_index],
                                    atoms: 2,
                                    coordinates: 1,
                                }
                            );
                        }
                        other => panic!("lane {lane} must retain its row error: {other:?}"),
                    }
                } else {
                    assert!(
                        matches!(&outcomes[lane], Ok(Ok(()))),
                        "independent valid lane {lane} must finish"
                    );
                }
            }
            assert_eq!(
                conformers
                    .iter()
                    .map(|conformer| conformer.id)
                    .collect::<Vec<_>>(),
                IDS
            );
            drop(conformers);
            assert_eq!(results, expected_results);
            assert_eq!(results[row_count], sentinel);
            assert_eq!(coordinate_bits(&coordinate_rows), before);
            assert_eq!(
                coordinate_rows.iter().map(Vec::len).collect::<Vec<_>>(),
                before_lengths
            );
            assert_eq!(cf3d_uff_one_kernel_counts(), (3, 0, 3));
        }

        // A zero-point field fails only when a lane reaches its first row.
        for row_count in [0_usize, 1, 7] {
            for worker_count in 1_i32..=4 {
                let mut coordinate_rows = vec![Vec::<[f64; 3]>::new(); row_count];
                let before = coordinate_bits(&coordinate_rows);
                let mut conformers = coordinate_rows
                    .iter_mut()
                    .enumerate()
                    .map(|(row, positions)| super::SerialConformer {
                        id: IDS[row],
                        positions: positions.as_mut_slice(),
                    })
                    .collect::<Vec<_>>();
                let mut results = vec![sentinel; row_count + 1];
                let mut source_field = ForceField::new(3);
                source_field.add_contribution(Box::new(CoordinateNormContribution));
                cf3d_uff_one_kernel_counts_reset();

                let outcomes = super::optimize_conformers_mt(
                    &source_field,
                    &mut conformers,
                    &mut results,
                    worker_count,
                    0,
                    0,
                )
                .expect("the zero-point fixture has sufficient result capacity");
                owner_calls += 1;
                assert_eq!(outcomes.len(), worker_count as usize);
                for lane in 0..worker_count as usize {
                    if lane < row_count {
                        match &outcomes[lane] {
                            Ok(Err(super::SerialConformerOptimizationError::Optimization {
                                input_index,
                                conformer_id,
                                source:
                                    super::OptimizationStageError::Initialize(
                                        ForceFieldKernelError::NoPoints,
                                    ),
                            })) => {
                                assert_eq!(*input_index, lane);
                                assert_eq!(*conformer_id, IDS[lane]);
                            }
                            other => panic!(
                                "nonempty lane {lane} must fail its first row at initialize: {other:?}"
                            ),
                        }
                    } else {
                        assert!(
                            matches!(&outcomes[lane], Ok(Ok(()))),
                            "empty lane {lane} has no initialize call"
                        );
                    }
                }
                assert!(results.iter().all(|result| *result == sentinel));
                assert!(
                    conformers
                        .iter()
                        .enumerate()
                        .all(|(row, conformer)| conformer.id == IDS[row])
                );
                drop(conformers);
                assert_eq!(coordinate_bits(&coordinate_rows), before);
                let lanes = worker_count as usize;
                assert_eq!(cf3d_uff_one_kernel_counts(), (lanes, 0, lanes));
            }
        }

        // Each copied field starts with an independent energy-call counter.
        // One completed row consumes two calls at minimize(0) and three at
        // minimize(1), so fail_on3 and fail_on6 both reach the second selected
        // row, at minimize and final energy respectively.
        const MATRIX_COUNTS: [usize; 4] = [1, 3, 6, 7];
        const WORKER_COUNTS: [i32; 3] = [1, 2, 3];
        let one_step_x = SOURCE_X.map(|x| x + (-(2.0 * x * 0.1)));
        let one_step_energies = one_step_x.map(|x| x * x);
        for row_count in MATRIX_COUNTS {
            for worker_count in WORKER_COUNTS {
                for (fail_on, max_iterations, expected_stage) in
                    [(3_usize, 0_i32, 0_u8), (6_usize, 1_i32, 1_u8)]
                {
                    let mut coordinate_rows = SOURCE_X[..row_count]
                        .iter()
                        .copied()
                        .map(|x| vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0]])
                        .collect::<Vec<_>>();
                    let mut expected_coordinates = coordinate_rows.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(row, positions)| super::SerialConformer {
                            id: IDS[row],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; row_count + 1];
                    let mut source_field = ForceField::new(3);
                    source_field.add_contribution(Box::new(FailAtEnergyCallWithCoordinateNorm {
                        calls: Cell::new(0),
                        fail_on,
                    }));
                    cf3d_uff_one_kernel_counts_reset();

                    let outcomes = super::optimize_conformers_mt(
                        &source_field,
                        &mut conformers,
                        &mut results,
                        worker_count,
                        2,
                        max_iterations,
                    )
                    .expect("the fixed failure matrix has sufficient result capacity");
                    owner_calls += 1;
                    assert_eq!(outcomes.len(), worker_count as usize);

                    let mut expected_results = vec![sentinel; row_count + 1];
                    for lane in 0..worker_count as usize {
                        let lane_indices = (lane..row_count)
                            .step_by(worker_count as usize)
                            .collect::<Vec<_>>();
                        let failed_index = lane_indices.get(1).copied();
                        if let Some(failed_index) = failed_index {
                            let expected_error = if expected_stage == 0 {
                                super::OptimizationStageError::Minimize(expected_pair_index_error())
                            } else {
                                let completed_index = lane_indices[0];
                                expected_coordinates[completed_index][1][0] =
                                    one_step_x[completed_index];
                                expected_coordinates[failed_index][1][0] = one_step_x[failed_index];
                                super::OptimizationStageError::FinalEnergy(
                                    expected_pair_index_error(),
                                )
                            };
                            match &outcomes[lane] {
                                Ok(Err(
                                    super::SerialConformerOptimizationError::Optimization {
                                        input_index,
                                        conformer_id,
                                        source,
                                    },
                                )) => {
                                    assert_eq!(*input_index, failed_index);
                                    assert_eq!(*conformer_id, IDS[failed_index]);
                                    assert_eq!(*source, expected_error);
                                }
                                other => panic!(
                                    "lane {lane} failed at source row {failed_index}: {other:?}"
                                ),
                            }
                            let completed_index = lane_indices[0];
                            expected_results[completed_index] = OptimizationOutcome {
                                status: 1,
                                energy: if max_iterations == 0 {
                                    SOURCE_X[completed_index] * SOURCE_X[completed_index]
                                } else {
                                    one_step_energies[completed_index]
                                },
                            };
                        } else {
                            assert!(
                                matches!(&outcomes[lane], Ok(Ok(()))),
                                "lane {lane} has no source row at its failure ordinal"
                            );
                            for source_index in lane_indices {
                                expected_results[source_index] = OptimizationOutcome {
                                    status: 1,
                                    energy: if max_iterations == 0 {
                                        SOURCE_X[source_index] * SOURCE_X[source_index]
                                    } else {
                                        one_step_energies[source_index]
                                    },
                                };
                                if expected_stage == 1 {
                                    expected_coordinates[source_index][1][0] =
                                        one_step_x[source_index];
                                }
                            }
                        }
                    }

                    assert_eq!(
                        conformers
                            .iter()
                            .map(|conformer| conformer.id)
                            .collect::<Vec<_>>(),
                        IDS[..row_count]
                    );
                    drop(conformers);
                    assert_eq!(results, expected_results);
                    assert_eq!(results[row_count], sentinel);
                    assert_eq!(
                        coordinate_bits(&coordinate_rows),
                        coordinate_bits(&expected_coordinates)
                    );
                    assert_eq!(
                        cf3d_uff_one_kernel_counts(),
                        (worker_count as usize, 0, worker_count as usize)
                    );
                }
            }
        }
        assert_eq!(owner_calls, 40);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_mt_m10_entry_copy_cache_and_storage_matrix() {
        use crate::kernel::Cf3dFragAcceptContributionIdentity as Identity;

        struct ReleaseAll(Vec<std::sync::mpsc::Sender<()>>);

        impl ReleaseAll {
            fn release_all(&mut self) {
                for release in std::mem::take(&mut self.0) {
                    let _ = release.send(());
                }
            }
        }

        impl Drop for ReleaseAll {
            fn drop(&mut self) {
                self.release_all();
            }
        }

        const IDS: [usize; 7] = [80, 7, 900, 31, 502, 8, 640];
        const TARGET_DISTANCES: [f64; 7] = [4.0, 5.0, 6.0, 4.25, 5.25, 6.25, 7.0];
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let reference_points = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut reference = worker_w06_coordinates(&reference_points);
        let mut source_field =
            worker_w06_source_field(&topology, &mut reference, &total_valences, 10.0, false);
        crate::kernel::cf3d_uff_one_set_fixed_points_for_serial_test(&mut source_field, &[1]);
        source_field
            .initialize()
            .expect("the fixed W06 source field initializes before transport");
        assert_eq!(
            crate::kernel::cf3d_uff_one_prefill_distance_cache_for_worker_test(&mut source_field),
            Ok(2.0)
        );
        let source_state = crate::kernel::cf3d_uff_one_copy_state_for_worker_test(&source_field);
        assert!(source_state.0);
        assert_eq!(source_state.1, 2);
        assert!(
            source_state.2[1] > 0.0,
            "the source pair distance is cached"
        );
        assert!(source_state.3);
        assert_eq!(source_state.4, [1]);
        assert_eq!(source_state.5, 2);
        assert_eq!(source_state.6, 3);
        assert_eq!(source_state.7, 1);
        let expected_identity = vec![Identity::Vdw {
            at1_idx: 0,
            at2_idx: 1,
            x_ij: minimum,
            well_depth,
            threshold: 10.0 * minimum,
        }];
        assert_eq!(
            crate::kernel::cf3d_frag_accept_contribution_identities(&source_field),
            expected_identity
        );

        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut owner_calls = 0;
        for row_count in [0_usize, 1, 7] {
            for worker_count in 1_usize..=4 {
                let worker_count_i32 =
                    i32::try_from(worker_count).expect("the fixed worker count fits i32");
                let mut coordinate_rows = TARGET_DISTANCES[..row_count]
                    .iter()
                    .copied()
                    .map(|distance| vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
                    .collect::<Vec<_>>();
                let rows_before = coordinate_rows.clone();
                let mut conformers = coordinate_rows
                    .iter_mut()
                    .enumerate()
                    .map(|(row, positions)| super::SerialConformer {
                        id: IDS[row],
                        positions: positions.as_mut_slice(),
                    })
                    .collect::<Vec<_>>();
                let mut results = vec![sentinel; row_count + 1];

                let (observer, started_rx, release_senders, worker_work_rx, handle_storage_rx) =
                    super::UffMtM10Observer::new(worker_count);
                let controller = std::thread::spawn(move || {
                    let mut release_all = ReleaseAll(release_senders);
                    let mut started = Vec::with_capacity(worker_count);
                    for _ in 0..worker_count {
                        match started_rx.recv_timeout(std::time::Duration::from_secs(30)) {
                            Ok(entry) => started.push(entry),
                            Err(error) => {
                                return Err(format!(
                                    "timed out before all M10 lanes entered: {error:?}"
                                ));
                            }
                        }
                    }
                    started.sort_unstable_by_key(|entry| entry.lane_index);
                    if started
                        .iter()
                        .enumerate()
                        .any(|(expected, entry)| entry.lane_index != expected)
                    {
                        return Err("M10 entry reports do not cover each lane once".to_owned());
                    }
                    release_all.release_all();
                    Ok(started)
                });

                cf3d_uff_one_kernel_counts_reset();
                let outcomes = {
                    let _observer_scope =
                        super::UffMtM10ObserverScope::install(std::sync::Arc::clone(&observer));
                    super::optimize_conformers_mt(
                        &source_field,
                        &mut conformers,
                        &mut results,
                        worker_count_i32,
                        2,
                        0,
                    )
                    .expect("the fixed M10 result vector covers all source rows")
                };
                owner_calls += 1;

                let started = controller
                    .join()
                    .expect("the M10 controller returns without panicking")
                    .unwrap_or_else(|error| panic!("M10 controller failed safely: {error}"));
                assert_eq!(started.len(), worker_count);
                for (lane, entry) in started.iter().enumerate() {
                    let expected_rows = (0..row_count)
                        .filter(|source_index| source_index % worker_count == lane)
                        .count();
                    assert_eq!(entry.lane_index, lane);
                    assert_eq!(entry.row_len, expected_rows);
                    assert_eq!(entry.row_capacity, expected_rows);
                    assert_eq!(entry.result_slot_len, expected_rows);
                    assert_eq!(entry.result_slot_capacity, expected_rows);
                    assert_eq!(
                        entry.copied_field_state,
                        (false, 0, Vec::new(), false, Vec::new(), 2, 0, 1),
                        "worker lane {lane} starts with source-copied terms and reset transient state"
                    );
                }

                assert_eq!(outcomes.len(), worker_count);
                assert!(outcomes.iter().all(|outcome| matches!(outcome, Ok(Ok(())))));
                let (handle_len, handle_capacity) = handle_storage_rx
                    .recv_timeout(std::time::Duration::from_secs(30))
                    .expect("the parent reports the joined-handle storage once");
                assert_eq!(handle_len, worker_count);
                assert_eq!(handle_capacity, worker_count);

                let mut worker_work = vec![None; worker_count];
                for _ in 0..worker_count {
                    let (lane, counts) = worker_work_rx
                        .recv_timeout(std::time::Duration::from_secs(30))
                        .expect("every released child reports its local work counters");
                    assert!(lane < worker_count);
                    assert!(worker_work[lane].replace(counts).is_none());
                }
                for lane in 0..worker_count {
                    let expected_rows = (0..row_count)
                        .filter(|source_index| source_index % worker_count == lane)
                        .count();
                    let counts = worker_work[lane]
                        .expect("each lane reports its private position-buffer observations");
                    assert_eq!(counts.position_buffer_growths, 1);
                    assert_eq!(counts.initialize_calls, expected_rows);
                    assert_ne!(counts.initial_position_buffer_address, 0);
                    assert_eq!(
                        counts.initial_position_buffer_address,
                        counts.final_position_buffer_address
                    );
                    assert!(counts.initial_position_buffer_capacity >= 2);
                    assert_eq!(
                        counts.initial_position_buffer_capacity,
                        counts.final_position_buffer_capacity
                    );
                }

                assert_eq!(
                    cf3d_uff_one_kernel_counts(),
                    (worker_count, 0, worker_count),
                    "one parent-local field and one ordered VDW term copy per spawned lane"
                );
                assert_eq!(
                    conformers
                        .iter()
                        .map(|conformer| conformer.id)
                        .collect::<Vec<_>>(),
                    IDS[..row_count]
                );
                drop(conformers);
                assert_eq!(coordinate_rows, rows_before);
                for source_index in 0..row_count {
                    assert_eq!(results[source_index].status, 1);
                    assert!(
                        (results[source_index].energy
                            - source_energy(TARGET_DISTANCES[source_index]))
                        .abs()
                            <= 1.0e-12,
                        "fixed W06 energy for source row {source_index}"
                    );
                }
                assert_eq!(results[row_count], sentinel);
                assert_eq!(
                    crate::kernel::cf3d_uff_one_copy_state_for_worker_test(&source_field),
                    source_state,
                    "worker copies preserve the initialized source cache and fixed points"
                );
                assert_eq!(
                    crate::kernel::cf3d_frag_accept_contribution_identities(&source_field),
                    expected_identity
                );
            }
        }
        assert_eq!(owner_calls, 12);
        drop(source_field);
        assert_eq!(reference.conformers_3d[0].coordinates(), reference_points);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_mt_m07_join_worker_handles_preserves_empty_and_success() {
        let empty = std::thread::scope(|_scope| super::join_worker_handles(Vec::new()));
        assert!(empty.is_empty());

        let outcomes = std::thread::scope(|scope| {
            let handles = (0..3)
                .map(|_| scope.spawn(|| Ok::<(), super::SerialConformerOptimizationError>(())))
                .collect();
            super::join_worker_handles(handles)
        });
        assert_eq!(outcomes.len(), 3);
        assert!(
            outcomes
                .into_iter()
                .all(|outcome| matches!(outcome, Ok(Ok(()))))
        );
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_mt_m07_join_worker_handles_preserves_typed_failure_positions() {
        fn fixed_error(lane: usize) -> super::SerialConformerOptimizationError {
            super::SerialConformerOptimizationError::CoordinateCount {
                input_index: lane,
                source: super::SerialCoordinateCountError {
                    conformer_id: 70 + lane,
                    atoms: 2,
                    coordinates: 1,
                },
            }
        }

        for failed_lane in 0..3 {
            let outcomes = std::thread::scope(|scope| {
                let handles = (0..3)
                    .map(|lane| {
                        scope.spawn(move || {
                            if lane == failed_lane {
                                Err(fixed_error(lane))
                            } else {
                                Ok(())
                            }
                        })
                    })
                    .collect();
                super::join_worker_handles(handles)
            });

            assert_eq!(outcomes.len(), 3);
            for (lane, outcome) in outcomes.into_iter().enumerate() {
                if lane == failed_lane {
                    match outcome {
                        Ok(Err(actual)) => assert_eq!(actual, fixed_error(lane)),
                        _ => panic!("typed failure for lane {lane} was not retained"),
                    }
                } else {
                    assert!(matches!(outcome, Ok(Ok(()))));
                }
            }
        }
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_mt_m07_join_worker_handles_preserves_panics_and_completion_order() {
        use std::sync::mpsc;

        struct ExitSignal {
            lane: usize,
            completed: mpsc::Sender<usize>,
            next_gate: Option<mpsc::Sender<()>>,
        }

        impl Drop for ExitSignal {
            fn drop(&mut self) {
                let _ = self.completed.send(self.lane);
                if let Some(next_gate) = self.next_gate.take() {
                    let _ = next_gate.send(());
                }
            }
        }

        for panic_lanes in [vec![1], vec![0, 2]] {
            let (completion_order, outcomes) = std::thread::scope(|scope| {
                let (gate0_tx, gate0_rx) = mpsc::channel();
                let (gate1_tx, gate1_rx) = mpsc::channel();
                let (gate2_tx, gate2_rx) = mpsc::channel();
                let (completed_tx, completed_rx) = mpsc::channel();
                let gates = [gate0_rx, gate1_rx, gate2_rx];
                let next_gates = [None, Some(gate0_tx), Some(gate1_tx)];
                let mut handles = Vec::with_capacity(3);

                for (lane, (gate, next_gate)) in gates.into_iter().zip(next_gates).enumerate() {
                    let completed = completed_tx.clone();
                    let should_panic = panic_lanes.contains(&lane);
                    handles.push(scope.spawn(move || {
                        let _ = gate.recv();
                        let exit_signal = ExitSignal {
                            lane,
                            completed,
                            next_gate,
                        };
                        if should_panic {
                            std::panic::panic_any(format!("lane-{lane}-panic-payload"));
                        }
                        drop(exit_signal);
                        Ok::<(), super::SerialConformerOptimizationError>(())
                    }));
                }
                drop(completed_tx);
                gate2_tx
                    .send(())
                    .expect("the third lane waits on the initial release gate");
                let completion_order = (0..3)
                    .map(|_| {
                        completed_rx
                            .recv()
                            .expect("every lane reports its exit before releasing its predecessor")
                    })
                    .collect::<Vec<_>>();
                let outcomes = super::join_worker_handles(handles);
                (completion_order, outcomes)
            });

            assert_eq!(completion_order, vec![2, 1, 0]);
            assert_eq!(outcomes.len(), 3);
            for (lane, outcome) in outcomes.into_iter().enumerate() {
                if panic_lanes.contains(&lane) {
                    let payload = match outcome {
                        Err(payload) => payload,
                        Ok(_) => panic!("child panic for lane {lane} was not retained"),
                    };
                    let expected = format!("lane-{lane}-panic-payload");
                    assert_eq!(
                        payload.downcast_ref::<String>(),
                        Some(&expected),
                        "panic payload remains at its creation-order lane"
                    );
                } else {
                    assert!(matches!(outcome, Ok(Ok(()))));
                }
            }
        }
    }

    #[test]
    fn uff_worker_w03_selected_bad_row_keeps_its_exact_conformer_id() {
        let mut skipped_coordinates = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut selected_coordinates = [[8.0, 9.0, 10.0]];
        let mut conformers = [
            super::SerialConformer {
                id: 800,
                positions: &mut skipped_coordinates,
            },
            super::SerialConformer {
                id: 7,
                positions: &mut selected_coordinates,
            },
        ];
        let mut field = ForceField::new(3);
        field.positions_mut().reserve(2);
        let mut assignment =
            super::WorkerConformerAssignment::new(1, std::num::NonZeroU32::new(2).unwrap());
        let mut selected_error = None;
        for conformer in &mut conformers {
            if assignment.selects_current() {
                selected_error =
                    super::bind_serial_conformer_positions(&mut field, conformer, 2).err();
            }
            assignment.advance();
        }
        assert_eq!(
            selected_error,
            Some(super::SerialCoordinateCountError {
                conformer_id: 7,
                atoms: 2,
                coordinates: 1,
            })
        );
        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 0);
        drop(field);
        drop(conformers);
        assert_eq!(skipped_coordinates, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        assert_eq!(selected_coordinates, [[8.0, 9.0, 10.0]]);
    }

    #[test]
    fn uff_worker_w04_only_selected_row_gets_status_and_result() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut skipped_first = [[1.0, 2.0, 3.0]];
        let mut selected = [[2.0, -1.0, 0.5]];
        let mut skipped_last = [[-3.0, 2.0, 1.0]];
        let before = [skipped_first, selected, skipped_last];
        let mut conformers = [
            super::SerialConformer {
                id: 80,
                positions: &mut skipped_first,
            },
            super::SerialConformer {
                id: 7,
                positions: &mut selected,
            },
            super::SerialConformer {
                id: 900,
                positions: &mut skipped_last,
            },
        ];
        let mut results = [sentinel; 3];
        cf3d_uff_one_kernel_counts_reset();
        let mut field = ForceField::new(3);
        field.add_contribution(Box::new(CoordinateNormContribution));

        assert_eq!(
            super::optimize_worker_conformers(
                field,
                &mut conformers,
                &mut results,
                1,
                1,
                std::num::NonZeroU32::new(2).unwrap(),
                0,
            ),
            Ok(())
        );
        drop(conformers);

        assert_eq!(skipped_first, before[0]);
        assert_eq!(selected, before[1]);
        assert_eq!(skipped_last, before[2]);
        assert_eq!(
            results,
            [
                sentinel,
                OptimizationOutcome {
                    status: 1,
                    energy: 5.25,
                },
                sentinel,
            ]
        );
        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 1);
    }

    #[test]
    fn uff_worker_w04_zero_and_negative_iterations_keep_source_statuses() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut zero_iteration_row = [[2.0, 0.0, 0.0]];
        let zero_iteration_before = zero_iteration_row;
        let mut zero_iteration_conformers = [super::SerialConformer {
            id: 32,
            positions: &mut zero_iteration_row,
        }];
        let mut zero_iteration_results = [sentinel];
        let mut zero_iteration_field = ForceField::new(3);
        zero_iteration_field.add_contribution(Box::new(CoordinateNormContribution));
        assert_eq!(
            super::optimize_worker_conformers(
                zero_iteration_field,
                &mut zero_iteration_conformers,
                &mut zero_iteration_results,
                1,
                0,
                std::num::NonZeroU32::new(1).unwrap(),
                0,
            ),
            Ok(())
        );
        drop(zero_iteration_conformers);
        assert_eq!(zero_iteration_row, zero_iteration_before);
        assert_eq!(
            zero_iteration_results,
            [OptimizationOutcome {
                status: 1,
                energy: 4.0,
            }]
        );

        let mut negative_iteration_row = [[2.0, 0.0, 0.0]];
        let negative_iteration_before = negative_iteration_row;
        let mut negative_iteration_conformers = [super::SerialConformer {
            id: 91,
            positions: &mut negative_iteration_row,
        }];
        let mut negative_iteration_results = [sentinel];
        let mut negative_iteration_field = ForceField::new(3);
        negative_iteration_field.add_contribution(Box::new(CoordinateNormContribution));
        assert_eq!(
            super::optimize_worker_conformers(
                negative_iteration_field,
                &mut negative_iteration_conformers,
                &mut negative_iteration_results,
                1,
                0,
                std::num::NonZeroU32::new(1).unwrap(),
                -1,
            ),
            Ok(())
        );
        drop(negative_iteration_conformers);
        let independent_energy = negative_iteration_row
            .iter()
            .flatten()
            .map(|coordinate| coordinate * coordinate)
            .sum::<f64>();
        assert_ne!(negative_iteration_row, negative_iteration_before);
        assert_eq!(negative_iteration_results[0].status, 0);
        assert!((negative_iteration_results[0].energy - independent_energy).abs() <= 1.0e-12);
        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 2);
    }

    #[test]
    fn uff_worker_w04_no_selected_lane_does_not_bind_or_initialize() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut short_row: [[f64; 3]; 0] = [];
        let mut conformers = [super::SerialConformer {
            id: 44,
            positions: &mut short_row,
        }];
        let mut results = [sentinel];
        cf3d_uff_one_kernel_counts_reset();
        let field = ForceField::new(3);

        assert_eq!(
            super::optimize_worker_conformers(
                field,
                &mut conformers,
                &mut results,
                1,
                2,
                std::num::NonZeroU32::new(2).unwrap(),
                0,
            ),
            Ok(())
        );
        assert_eq!(results, [sentinel]);
        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 0);
    }

    #[test]
    fn uff_worker_w04_final_energy_error_keeps_exact_index_id_and_call_order() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut skipped_coordinates: [[f64; 3]; 0] = [];
        let mut selected_coordinates = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let selected_before = selected_coordinates;
        let mut conformers = [
            super::SerialConformer {
                id: 800,
                positions: &mut skipped_coordinates,
            },
            super::SerialConformer {
                id: 7,
                positions: &mut selected_coordinates,
            },
        ];
        let mut results = [sentinel; 2];
        cf3d_uff_one_kernel_counts_reset();
        let mut field = ForceField::new(3);
        field.add_contribution(Box::new(FailAtEnergyCall {
            calls: Cell::new(0),
            fail_on: 2,
        }));

        assert_eq!(
            super::optimize_worker_conformers(
                field,
                &mut conformers,
                &mut results,
                2,
                1,
                std::num::NonZeroU32::new(2).unwrap(),
                0,
            ),
            Err(super::SerialConformerOptimizationError::Optimization {
                input_index: 1,
                conformer_id: 7,
                source: super::OptimizationStageError::FinalEnergy(expected_pair_index_error()),
            })
        );
        drop(conformers);
        assert_eq!(selected_coordinates, selected_before);
        assert_eq!(results, [sentinel; 2]);
        assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 1);
    }

    #[test]
    fn uff_worker_w05_initialize_failures_report_selected_input_positions() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let ids = [80, 7, 900];
        for row_count in 1..=3 {
            let lane_count = std::num::NonZeroU32::new(row_count as u32).unwrap();
            for failed_input_index in 0..row_count {
                let mut coordinate_rows = vec![Vec::<[f64; 3]>::new(); row_count];
                let before = coordinate_rows.clone();
                let mut conformers = coordinate_rows
                    .iter_mut()
                    .enumerate()
                    .map(|(input_index, positions)| super::SerialConformer {
                        id: ids[input_index],
                        positions: positions.as_mut_slice(),
                    })
                    .collect::<Vec<_>>();
                let mut results = vec![sentinel; row_count];
                cf3d_uff_one_kernel_counts_reset();

                assert_eq!(
                    super::optimize_worker_conformers(
                        ForceField::new(3),
                        &mut conformers,
                        &mut results,
                        0,
                        failed_input_index as u32,
                        lane_count,
                        0,
                    ),
                    Err(super::SerialConformerOptimizationError::Optimization {
                        input_index: failed_input_index,
                        conformer_id: ids[failed_input_index],
                        source: super::OptimizationStageError::Initialize(
                            ForceFieldKernelError::NoPoints
                        ),
                    })
                );
                drop(conformers);

                assert_eq!(coordinate_rows, before);
                assert_eq!(results, vec![sentinel; row_count]);
                assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 1);
            }
        }
    }

    #[test]
    fn uff_worker_w05_minimize_and_final_energy_failures_keep_result_order() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let ids = [80, 7, 900, 31, 502, 8];
        for selected_count in 1..=3 {
            let row_count = selected_count * 2;
            for failed_selected_ordinal in 0..selected_count {
                for stage in [
                    super::OptimizationStageError::Minimize(expected_pair_index_error()),
                    super::OptimizationStageError::FinalEnergy(expected_pair_index_error()),
                ] {
                    let mut coordinate_rows = (0..row_count)
                        .map(|input_index| {
                            if input_index % 2 == 0 {
                                let selected_ordinal = input_index / 2;
                                vec![[0.0, 0.0, 0.0], [2.0 + selected_ordinal as f64, 0.0, 0.0]]
                            } else {
                                Vec::new()
                            }
                        })
                        .collect::<Vec<_>>();
                    let before = coordinate_rows.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(input_index, positions)| super::SerialConformer {
                            id: ids[input_index],
                            positions: positions.as_mut_slice(),
                        })
                        .collect::<Vec<_>>();
                    let mut results = vec![sentinel; row_count];
                    let failed_input_index = failed_selected_ordinal * 2;
                    let fail_on = match stage {
                        super::OptimizationStageError::Minimize(_) => {
                            2 * failed_selected_ordinal + 1
                        }
                        super::OptimizationStageError::FinalEnergy(_) => {
                            2 * failed_selected_ordinal + 2
                        }
                        super::OptimizationStageError::Initialize(_) => unreachable!(),
                    };
                    cf3d_uff_one_kernel_counts_reset();
                    let mut field = ForceField::new(3);
                    field.add_contribution(Box::new(FailAtEnergyCallWithCoordinateNorm {
                        calls: Cell::new(0),
                        fail_on,
                    }));

                    assert_eq!(
                        super::optimize_worker_conformers(
                            field,
                            &mut conformers,
                            &mut results,
                            2,
                            0,
                            std::num::NonZeroU32::new(2).unwrap(),
                            0,
                        ),
                        Err(super::SerialConformerOptimizationError::Optimization {
                            input_index: failed_input_index,
                            conformer_id: ids[failed_input_index],
                            source: stage,
                        })
                    );
                    drop(conformers);

                    let mut expected_results = vec![sentinel; row_count];
                    for completed_ordinal in 0..failed_selected_ordinal {
                        expected_results[2 * completed_ordinal] = OptimizationOutcome {
                            status: 1,
                            energy: (2.0 + completed_ordinal as f64).powi(2),
                        };
                    }
                    assert_eq!(results, expected_results);
                    assert_eq!(coordinate_rows, before);
                    assert_eq!(
                        cf3d_uff_one_serial_work_counts().initialize_calls,
                        failed_selected_ordinal + 1,
                    );
                }
            }
        }
    }

    #[test]
    fn uff_worker_w05_final_energy_failure_keeps_minimized_coordinates() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let ids = [80, 7, 900, 31, 502, 8];
        for selected_count in 1..=3 {
            let row_count = selected_count * 2;
            for failed_selected_ordinal in 0..selected_count {
                let mut coordinate_rows = (0..row_count)
                    .map(|input_index| {
                        if input_index % 2 != 0 {
                            Vec::new()
                        } else {
                            vec![[1.0, 2.0, 0.5], [3.0, 4.0, 0.25]]
                        }
                    })
                    .collect::<Vec<_>>();
                let before = coordinate_rows.clone();
                let mut conformers = coordinate_rows
                    .iter_mut()
                    .enumerate()
                    .map(|(input_index, positions)| super::SerialConformer {
                        id: ids[input_index],
                        positions: positions.as_mut_slice(),
                    })
                    .collect::<Vec<_>>();
                let mut results = vec![sentinel; row_count];
                let failed_input_index = 2 * failed_selected_ordinal;
                cf3d_uff_one_kernel_counts_reset();
                let mut field = ForceField::new(3);
                field.add_contribution(Box::new(FailAtEnergyCallWithCoordinateNorm {
                    calls: Cell::new(0),
                    fail_on: 3 * failed_selected_ordinal + 3,
                }));

                assert_eq!(
                    super::optimize_worker_conformers(
                        field,
                        &mut conformers,
                        &mut results,
                        2,
                        0,
                        std::num::NonZeroU32::new(2).unwrap(),
                        1,
                    ),
                    Err(super::SerialConformerOptimizationError::Optimization {
                        input_index: failed_input_index,
                        conformer_id: ids[failed_input_index],
                        source: super::OptimizationStageError::FinalEnergy(
                            expected_pair_index_error()
                        ),
                    })
                );
                drop(conformers);

                let mut expected_results = vec![sentinel; row_count];
                for completed_ordinal in 0..failed_selected_ordinal {
                    let completed_coordinates = &coordinate_rows[2 * completed_ordinal];
                    let independent_energy = completed_coordinates
                        .iter()
                        .flatten()
                        .map(|coordinate| coordinate * coordinate)
                        .sum::<f64>();
                    expected_results[2 * completed_ordinal] = OptimizationOutcome {
                        status: 1,
                        energy: independent_energy,
                    };
                }
                assert_eq!(results, expected_results);
                assert_ne!(
                    coordinate_rows[failed_input_index],
                    before[failed_input_index]
                );
                for input_index in 0..row_count {
                    if input_index % 2 == 0 && input_index < failed_input_index {
                        assert_ne!(coordinate_rows[input_index], before[input_index]);
                    } else if input_index != failed_input_index {
                        assert_eq!(coordinate_rows[input_index], before[input_index]);
                    }
                }
                assert_eq!(
                    cf3d_uff_one_serial_work_counts().initialize_calls,
                    failed_selected_ordinal + 1,
                );
            }
        }
    }

    #[test]
    fn uff_serial_s07_single_conformer_reuses_one_field_and_term() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut only = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut conformers = [super::SerialConformer {
            id: 903,
            positions: &mut only,
        }];
        let mut results = [sentinel];
        cf3d_uff_one_kernel_counts_reset();
        let mut field = ForceField::new(1);
        field
            .positions_mut()
            .extend(reference.iter_mut().map(|row| &mut row[..]));
        field.add_contribution(Box::new(SquaredDistanceContribution {
            first: 0,
            second: 1,
            scale: 1.0,
        }));
        assert_eq!(
            super::optimize_serial_conformers(field, &mut conformers, &mut results, 2, 0),
            Ok(())
        );
        assert_eq!(conformers[0].id, 903);
        assert_eq!(
            results,
            [OptimizationOutcome {
                status: 1,
                energy: 1.0
            }]
        );
        assert_eq!(only, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        assert_eq!(cf3d_uff_one_kernel_counts(), (1, 1, 0));
    }

    #[test]
    fn uff_serial_s07_multiple_conformers_keep_input_order_and_term_identity() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut first = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut second = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let mut third = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let mut conformers = [
            super::SerialConformer {
                id: 80,
                positions: &mut first,
            },
            super::SerialConformer {
                id: 2,
                positions: &mut second,
            },
            super::SerialConformer {
                id: 500,
                positions: &mut third,
            },
        ];
        let mut results = [sentinel; 3];
        cf3d_uff_one_kernel_counts_reset();
        let mut field = ForceField::new(1);
        field
            .positions_mut()
            .extend(reference.iter_mut().map(|row| &mut row[..]));
        field.add_contribution(Box::new(SquaredDistanceContribution {
            first: 0,
            second: 1,
            scale: 1.0,
        }));
        assert_eq!(
            super::optimize_serial_conformers(field, &mut conformers, &mut results, 2, 0),
            Ok(())
        );
        assert_eq!(
            conformers
                .iter()
                .map(|conformer| conformer.id)
                .collect::<Vec<_>>(),
            [80, 2, 500]
        );
        assert_eq!(
            results,
            [
                OptimizationOutcome {
                    status: 1,
                    energy: 4.0
                },
                OptimizationOutcome {
                    status: 1,
                    energy: 9.0
                },
                OptimizationOutcome {
                    status: 1,
                    energy: 16.0
                },
            ]
        );
        assert_eq!(first, [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        assert_eq!(second, [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        assert_eq!(third, [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]]);
        assert_eq!(cf3d_uff_one_kernel_counts(), (1, 1, 0));
    }

    #[test]
    fn uff_serial_s11_fixed_points_survive_rebinding_across_conformers() {
        // ForceField.cpp:329-374 zeros every component of each fixed point
        // after contribution gradients; initialize and the serial positions
        // loop leave the separate fixed-point vector intact.
        let fixed_point_cases: [&[i32]; 4] = [&[], &[0], &[1], &[0, 1]];
        for fixed_points in fixed_point_cases {
            let mut reference = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
            let first_fixed = fixed_points.contains(&0);
            let second_fixed = fixed_points.contains(&1);
            let mut first = [[1.0, 2.0, 0.5], [3.0, 4.0, 0.25]];
            let mut second = [[10.0, 1.0, 2.0], [14.0, 5.0, 3.0]];
            let first_before = first;
            let second_before = second;
            let mut conformers = [
                super::SerialConformer {
                    id: 80,
                    positions: &mut first,
                },
                super::SerialConformer {
                    id: 2,
                    positions: &mut second,
                },
            ];
            let mut results = [
                OptimizationOutcome {
                    status: -9,
                    energy: 123.0,
                },
                OptimizationOutcome {
                    status: -9,
                    energy: 123.0,
                },
            ];
            cf3d_uff_one_kernel_counts_reset();
            let mut field = ForceField::new(3);
            field
                .positions_mut()
                .extend(reference.iter_mut().map(|row| &mut row[..]));
            field.add_contribution(Box::new(CoordinateNormContribution));
            crate::kernel::cf3d_uff_one_set_fixed_points_for_serial_test(&mut field, fixed_points);
            // With every point fixed, the source gradient wrapper zeros the
            // full vector and BFGS reports its source-defined bad direction.
            // Use zero iterations for that mask to exercise both serial
            // rebinds and final-energy writes; the nonempty zero-gradient
            // failure remains covered by its dedicated optimizer regression.
            let max_iterations = if first_fixed && second_fixed { 0 } else { 1 };
            assert_eq!(
                super::optimize_serial_conformers(
                    field,
                    &mut conformers,
                    &mut results,
                    2,
                    max_iterations,
                ),
                Ok(())
            );

            for (actual, before, outcome) in [
                (&first, &first_before, results[0]),
                (&second, &second_before, results[1]),
            ] {
                if first_fixed {
                    assert_eq!(actual[0].map(f64::to_bits), before[0].map(f64::to_bits),);
                } else {
                    assert_ne!(actual[0][0].to_bits(), before[0][0].to_bits());
                }
                if second_fixed {
                    assert_eq!(actual[1].map(f64::to_bits), before[1].map(f64::to_bits),);
                } else {
                    assert_ne!(actual[1][0].to_bits(), before[1][0].to_bits());
                }
                let independently_computed_energy = actual
                    .iter()
                    .flatten()
                    .map(|coordinate| coordinate * coordinate)
                    .sum::<f64>();
                assert!((outcome.energy - independently_computed_energy).abs() <= 1.0e-12);
                assert!(matches!(outcome.status, 0 | 1));
            }
            assert_eq!(cf3d_uff_one_kernel_counts(), (1, 1, 0));
        }
    }

    #[test]
    fn uff_serial_s11_mask_iteration_matrix() {
        // ForceField.cpp:353-375 applies every fixed-point mask to each
        // conformer's gradient. BFGSOpt::linearSearch rejects the all-zero
        // direction; a positive-iteration error leaves gather/result writes
        // unperformed, while zero iterations is a normal status-1 return.
        let fixed_point_cases: [&[i32]; 4] = [&[], &[0], &[1], &[0, 1]];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        for fixed_points in fixed_point_cases {
            let first_fixed = fixed_points.contains(&0);
            let second_fixed = fixed_points.contains(&1);
            for max_iterations in [0, 1] {
                let mut reference = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
                let mut first = [[1.0, 2.0, 0.5], [3.0, 4.0, 0.25]];
                let mut second = [[10.0, 1.0, 2.0], [14.0, 5.0, 3.0]];
                let first_before = first;
                let second_before = second;
                let mut conformers = [
                    super::SerialConformer {
                        id: 80,
                        positions: &mut first,
                    },
                    super::SerialConformer {
                        id: 2,
                        positions: &mut second,
                    },
                ];
                let mut results = [sentinel; 2];
                cf3d_uff_one_kernel_counts_reset();
                let mut field = ForceField::new(3);
                field
                    .positions_mut()
                    .extend(reference.iter_mut().map(|row| &mut row[..]));
                field.add_contribution(Box::new(CoordinateNormContribution));
                crate::kernel::cf3d_uff_one_set_fixed_points_for_serial_test(
                    &mut field,
                    fixed_points,
                );

                let outcome = super::optimize_serial_conformers(
                    field,
                    &mut conformers,
                    &mut results,
                    2,
                    max_iterations,
                );
                let all_fixed_positive_iterations =
                    first_fixed && second_fixed && max_iterations == 1;
                if all_fixed_positive_iterations {
                    assert!(matches!(
                        outcome,
                        Err(super::SerialConformerOptimizationError::Optimization {
                            input_index: 0,
                            conformer_id: 80,
                            source: OptimizationStageError::Minimize(
                                ForceFieldKernelError::OptimizerBadDirection
                            ),
                        })
                    ));
                    drop(conformers);
                    assert_eq!(
                        first.map(|row| row.map(f64::to_bits)),
                        first_before.map(|row| row.map(f64::to_bits)),
                    );
                    assert_eq!(
                        second.map(|row| row.map(f64::to_bits)),
                        second_before.map(|row| row.map(f64::to_bits)),
                    );
                    assert_eq!(results, [sentinel; 2]);
                    assert_eq!(cf3d_uff_one_kernel_counts(), (1, 1, 0));
                    continue;
                }

                assert_eq!(outcome, Ok(()));
                assert_eq!(
                    conformers.iter().map(|row| row.id).collect::<Vec<_>>(),
                    [80, 2]
                );
                drop(conformers);
                for (input_index, (actual, before)) in
                    [(&first, &first_before), (&second, &second_before)]
                        .into_iter()
                        .enumerate()
                {
                    for atom_index in 0..2 {
                        let fixed = if atom_index == 0 {
                            first_fixed
                        } else {
                            second_fixed
                        };
                        if max_iterations == 0 || fixed {
                            assert_eq!(
                                actual[atom_index].map(f64::to_bits),
                                before[atom_index].map(f64::to_bits),
                                "mask={fixed_points:?}, max_iterations={max_iterations}, input={input_index}, atom={atom_index}",
                            );
                        } else {
                            assert_ne!(
                                actual[atom_index][0].to_bits(),
                                before[atom_index][0].to_bits(),
                                "mask={fixed_points:?}, input={input_index}, atom={atom_index}",
                            );
                        }
                    }
                    let independent_energy = actual
                        .iter()
                        .flatten()
                        .map(|coordinate| coordinate * coordinate)
                        .sum::<f64>();
                    assert!(
                        (results[input_index].energy - independent_energy).abs() <= 1.0e-12,
                        "mask={fixed_points:?}, max_iterations={max_iterations}, input={input_index}",
                    );
                    if max_iterations == 0 {
                        assert_eq!(results[input_index].status, 1);
                    } else {
                        assert_eq!(results[input_index].status, 1);
                    }
                }
                assert_eq!(cf3d_uff_one_kernel_counts(), (1, 1, 0));
            }
        }
    }

    #[test]
    fn uff_serial_s06_failure_first_middle_last_preserves_result_slots() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        for failed_slot in 0..3 {
            let mut first = [0.0];
            let mut second = [2.0];
            let mut field = ForceField::new(1);
            field
                .positions_mut()
                .extend([first.as_mut_slice(), second.as_mut_slice()]);
            field.add_contribution(Box::new(FailAtEnergyCall {
                calls: Cell::new(0),
                fail_on: 2 * (failed_slot + 1),
            }));
            let mut results = [sentinel; 4];
            for (slot, result) in results[..3].iter_mut().enumerate() {
                let outcome = super::assign_serial_result(&mut field, 0, result);
                if slot == failed_slot {
                    assert_eq!(
                        outcome,
                        Err(OptimizationStageError::FinalEnergy(
                            expected_pair_index_error()
                        ))
                    );
                    break;
                }
                assert_eq!(outcome, Ok(()));
            }
            assert_eq!(
                &results[..failed_slot],
                vec![
                    OptimizationOutcome {
                        status: 1,
                        energy: 4.0
                    };
                    failed_slot
                ]
            );
            assert!(
                results[failed_slot..]
                    .iter()
                    .all(|result| *result == sentinel)
            );
            drop(field);
            assert_eq!(first, [0.0]);
            assert_eq!(second, [2.0]);
        }
    }

    #[test]
    fn uff_serial_s06_assigns_successful_statuses_in_call_order() {
        let mut point = [7.0];
        let mut field = ForceField::new(1);
        field.positions_mut().push(&mut point);
        let mut result = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        assert_eq!(
            super::assign_serial_result(&mut field, 0, &mut result),
            Ok(())
        );
        assert_eq!(
            result,
            OptimizationOutcome {
                status: 0,
                energy: 0.0
            }
        );
    }

    #[test]
    fn uff_serial_s05_status_zero_one_and_zero_iteration_coordinates() {
        let mut rows = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut field = ForceField::new(3);
        field
            .positions_mut()
            .extend(rows.iter_mut().map(|row| &mut row[..]));
        assert_eq!(
            super::execute_serial_slot(&mut field, 0),
            Ok(OptimizationOutcome {
                status: 0,
                energy: 0.0
            })
        );
        field.add_contribution(Box::new(SquaredDistanceContribution {
            first: 0,
            second: 1,
            scale: 1.0,
        }));
        assert_eq!(
            super::execute_serial_slot(&mut field, 0),
            Ok(OptimizationOutcome {
                status: 1,
                energy: 4.0
            })
        );
        drop(field);
        assert_eq!(rows, [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    }

    #[test]
    fn uff_serial_s05_errors_preserve_initialize_minimize_final_energy_stages() {
        let mut empty = ForceField::new(3);
        assert_eq!(
            super::execute_serial_slot(&mut empty, 0),
            Err(OptimizationStageError::Initialize(
                ForceFieldKernelError::NoPoints
            ))
        );
        for fail_on in [1, 2] {
            let mut first = [0.0];
            let mut second = [2.0];
            let mut field = ForceField::new(1);
            field
                .positions_mut()
                .extend([first.as_mut_slice(), second.as_mut_slice()]);
            field.add_contribution(Box::new(FailAtEnergyCall {
                calls: Cell::new(0),
                fail_on,
            }));
            let expected = if fail_on == 1 {
                OptimizationStageError::Minimize(expected_pair_index_error())
            } else {
                OptimizationStageError::FinalEnergy(expected_pair_index_error())
            };
            assert_eq!(super::execute_serial_slot(&mut field, 0), Err(expected));
            drop(field);
            assert_eq!(first, [0.0]);
            assert_eq!(second, [2.0]);
        }
    }

    #[test]
    fn uff_serial_s04_rebind_borrows_distinct_and_original_reference_rows() {
        let mut reference = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let mut other = [[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]];
        let reference_address = reference[0].as_ptr();
        let other_address = other[0].as_ptr();
        let mut field = ForceField::new(3);
        field
            .positions_mut()
            .extend(reference.iter_mut().map(|row| &mut row[..]));
        field.initialize().unwrap();
        let field = field.rebind_positions(Vec::new());
        let mut original = super::SerialConformer {
            id: 901,
            positions: &mut reference,
        };
        let field = super::rebind_serial_conformer(field, &mut original, 2).unwrap();
        assert_eq!(field.positions()[0].as_ptr(), reference_address);
        assert_eq!(field.positions()[1], [4.0, 5.0, 6.0]);
        let mut distinct = super::SerialConformer {
            id: 4,
            positions: &mut other,
        };
        let mut field = super::rebind_serial_conformer(field, &mut distinct, 2).unwrap();
        assert_eq!(field.positions()[0].as_ptr(), other_address);
        field.positions_mut()[0][0] = 99.0;
        drop(field);
        assert_eq!(reference, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert_eq!(other[0], [99.0, 8.0, 9.0]);
    }

    #[test]
    fn uff_serial_s04_coordinate_count_error_and_zero_atom_no_access() {
        let mut rows = [[1.0, 2.0, 3.0]];
        let mut input = super::SerialConformer {
            id: 500,
            positions: &mut rows,
        };
        assert_eq!(
            super::rebind_serial_conformer(ForceField::new(3), &mut input, 2).err(),
            Some(super::SerialCoordinateCountError {
                conformer_id: 500,
                atoms: 2,
                coordinates: 1
            })
        );
        let field = super::rebind_serial_conformer(ForceField::new(3), &mut input, 0).unwrap();
        assert!(field.positions().is_empty());
        drop(field);
        assert_eq!(rows, [[1.0, 2.0, 3.0]]);
    }

    #[test]
    fn uff_serial_s03_capacity_matrix_preserves_order_labels_and_storage() {
        let sentinel = super::OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut first = [[1.0, 2.0, 3.0]];
        let mut second = [[4.0, 5.0, 6.0]];
        let inputs = [
            super::SerialConformer {
                id: 800,
                positions: &mut first,
            },
            super::SerialConformer {
                id: 7,
                positions: &mut second,
            },
        ];
        for capacity in [0, 1, 2, 3] {
            let results = vec![sentinel; capacity];
            let expected = if capacity < 2 {
                Err(super::SerialResultCapacityError {
                    conformers: 2,
                    result_slots: capacity,
                })
            } else {
                Ok(())
            };
            assert_eq!(
                super::validate_serial_result_capacity(&inputs, &results),
                expected
            );
            assert_eq!(results, vec![sentinel; capacity]);
            assert_eq!(
                inputs.iter().map(|input| input.id).collect::<Vec<_>>(),
                [800, 7]
            );
            assert_eq!(inputs[0].positions, [[1.0, 2.0, 3.0]]);
            assert_eq!(inputs[1].positions, [[4.0, 5.0, 6.0]]);
        }
        assert_eq!(super::validate_serial_result_capacity(&[], &[]), Ok(()));
        assert_eq!(
            super::validate_serial_result_capacity(&[], &[sentinel]),
            Ok(())
        );
    }

    use super::{
        OptimizationOutcome, OptimizationStageError, SingleConformerOptions, optimize_force_field,
    };
    use crate::kernel::{
        EvaluationContext, ForceField, ForceFieldContribution, ForceFieldIndexArgument,
        ForceFieldKernelError, cf3d_uff_one_kernel_counts, cf3d_uff_one_kernel_counts_reset,
        cf3d_uff_one_serial_work_counts,
    };

    struct SquaredDistanceContribution {
        first: u32,
        second: u32,
        scale: f64,
    }

    struct CoordinateNormContribution;

    impl ForceFieldContribution for CoordinateNormContribution {
        fn get_energy(
            &self,
            context: &mut EvaluationContext<'_>,
        ) -> Result<f64, ForceFieldKernelError> {
            Ok(context
                .coordinates()
                .iter()
                .map(|coordinate| coordinate * coordinate)
                .sum())
        }

        fn get_grad(
            &self,
            context: &mut EvaluationContext<'_>,
            gradient: &mut [f64],
        ) -> Result<(), ForceFieldKernelError> {
            for (slot, coordinate) in gradient.iter_mut().zip(context.coordinates()) {
                *slot += 2.0 * coordinate;
            }
            Ok(())
        }

        fn copy(&self) -> Box<dyn ForceFieldContribution> {
            Box::new(Self)
        }
    }

    impl ForceFieldContribution for SquaredDistanceContribution {
        fn get_energy(
            &self,
            context: &mut EvaluationContext<'_>,
        ) -> Result<f64, ForceFieldKernelError> {
            let distance = context.distance(self.first, self.second)?;
            Ok(self.scale * distance * distance)
        }

        fn get_grad(
            &self,
            context: &mut EvaluationContext<'_>,
            gradient: &mut [f64],
        ) -> Result<(), ForceFieldKernelError> {
            let _ = context.distance(self.first, self.second)?;
            let first = self.first as usize;
            let second = self.second as usize;
            let derivative =
                2.0 * self.scale * (context.coordinates()[first] - context.coordinates()[second]);
            gradient[first] += derivative;
            gradient[second] -= derivative;
            Ok(())
        }

        fn copy(&self) -> Box<dyn ForceFieldContribution> {
            Box::new(Self {
                first: self.first,
                second: self.second,
                scale: self.scale,
            })
        }
    }

    struct FailAtEnergyCall {
        calls: Cell<usize>,
        fail_on: usize,
    }

    impl ForceFieldContribution for FailAtEnergyCall {
        fn get_energy(
            &self,
            context: &mut EvaluationContext<'_>,
        ) -> Result<f64, ForceFieldKernelError> {
            let call = self.calls.get() + 1;
            self.calls.set(call);
            if call == self.fail_on {
                let _ = context.distance(2, 0)?;
            }
            let distance = context.distance(0, 1)?;
            Ok(distance * distance)
        }

        fn get_grad(
            &self,
            _context: &mut EvaluationContext<'_>,
            _gradient: &mut [f64],
        ) -> Result<(), ForceFieldKernelError> {
            Ok(())
        }

        fn copy(&self) -> Box<dyn ForceFieldContribution> {
            Box::new(Self {
                calls: Cell::new(self.calls.get()),
                fail_on: self.fail_on,
            })
        }
    }

    struct FailAtEnergyCallWithCoordinateNorm {
        calls: Cell<usize>,
        fail_on: usize,
    }

    impl ForceFieldContribution for FailAtEnergyCallWithCoordinateNorm {
        fn get_energy(
            &self,
            context: &mut EvaluationContext<'_>,
        ) -> Result<f64, ForceFieldKernelError> {
            let call = self.calls.get() + 1;
            self.calls.set(call);
            if call == self.fail_on {
                let _ = context.distance(2, 0)?;
            }
            Ok(context
                .coordinates()
                .iter()
                .map(|coordinate| coordinate * coordinate)
                .sum())
        }

        fn get_grad(
            &self,
            context: &mut EvaluationContext<'_>,
            gradient: &mut [f64],
        ) -> Result<(), ForceFieldKernelError> {
            for (slot, coordinate) in gradient.iter_mut().zip(context.coordinates()) {
                *slot += 2.0 * coordinate;
            }
            Ok(())
        }

        fn copy(&self) -> Box<dyn ForceFieldContribution> {
            Box::new(Self {
                calls: Cell::new(self.calls.get()),
                fail_on: self.fail_on,
            })
        }
    }

    fn expected_pair_index_error() -> ForceFieldKernelError {
        ForceFieldKernelError::IndexOutOfRange {
            argument: ForceFieldIndexArgument::I,
            index: 2,
            upper_bound: 2,
        }
    }

    #[test]
    fn uff_public_perf_iterator_preserves_capacity_row_order_and_errors() {
        const IDS: [usize; 3] = [80, 7, 900];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let cases: [(usize, &[usize]); 3] = [(0, &[0, 1]), (1, &[0, 1, 3]), (3, &[2, 3, 5])];

        for (count, result_lengths) in cases {
            for max_iterations in [0, 1] {
                for &result_len in result_lengths {
                    let mut first = [[1.0, 2.0, 3.0]];
                    let mut second = [[4.0, 5.0, 6.0]];
                    let mut third = [[7.0, 8.0, 9.0]];
                    let before = [first, second, third];
                    let rows = [&mut first[..], &mut second[..], &mut third[..]];
                    let conformers = rows
                        .into_iter()
                        .zip(IDS)
                        .take(count)
                        .map(|(positions, id)| super::SerialConformer { id, positions });
                    let mut results = vec![sentinel; result_len];
                    cf3d_uff_one_kernel_counts_reset();

                    let outcome = super::optimize_serial_conformers_iter(
                        ForceField::new(3),
                        conformers,
                        &mut results,
                        1,
                        max_iterations,
                    );
                    let expected_outcome = if result_len < count {
                        Err(super::SerialConformerOptimizationError::ResultCapacity(
                            super::SerialResultCapacityError {
                                conformers: count,
                                result_slots: result_len,
                            },
                        ))
                    } else {
                        Ok(())
                    };
                    assert_eq!(outcome, expected_outcome);

                    let mut expected_results = vec![sentinel; result_len];
                    if result_len >= count {
                        for result in &mut expected_results[..count] {
                            *result = OptimizationOutcome {
                                status: 0,
                                energy: 0.0,
                            };
                        }
                    }
                    assert_eq!(results, expected_results);
                    assert_eq!(
                        cf3d_uff_one_serial_work_counts().initialize_calls,
                        if result_len < count { 0 } else { count },
                        "count={count}, result_len={result_len}, max_iterations={max_iterations}"
                    );
                    assert_eq!([first, second, third], before);
                }
            }
        }

        for bad_row in [1, 2] {
            let mut first = [[2.0, 0.0, 0.0]];
            let mut second = [[3.0, 0.0, 0.0]];
            let mut third = [[4.0, 0.0, 0.0]];
            let before = [first, second, third];
            let mut empty: [[f64; 3]; 0] = [];
            let rows = if bad_row == 1 {
                [&mut first[..], &mut empty[..], &mut third[..]]
            } else {
                [&mut first[..], &mut second[..], &mut empty[..]]
            };
            let conformers = rows
                .into_iter()
                .zip(IDS)
                .map(|(positions, id)| super::SerialConformer { id, positions });
            let mut results = [sentinel; 4];
            cf3d_uff_one_kernel_counts_reset();
            let mut field = ForceField::new(3);
            field.add_contribution(Box::new(CoordinateNormContribution));

            let outcome =
                super::optimize_serial_conformers_iter(field, conformers, &mut results, 1, 0);
            let expected_id = IDS[bad_row];
            assert_eq!(
                outcome,
                Err(super::SerialConformerOptimizationError::CoordinateCount {
                    input_index: bad_row,
                    source: super::SerialCoordinateCountError {
                        conformer_id: expected_id,
                        atoms: 1,
                        coordinates: 0,
                    },
                })
            );
            let expected_prefix = [
                OptimizationOutcome {
                    status: 1,
                    energy: 4.0,
                },
                OptimizationOutcome {
                    status: 1,
                    energy: 9.0,
                },
            ];
            assert_eq!(&results[..bad_row], &expected_prefix[..bad_row]);
            assert!(results[bad_row..].iter().all(|result| *result == sentinel));
            assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, bad_row);
            assert_eq!([first, second, third], before);
        }
    }

    #[test]
    fn uff_public_perf_block_keeps_same_block_rows_lazy_and_resizes_in_source_order() {
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        const IDS: [usize; 3] = [80, 7, 900];
        const SENTINEL: OptimizationOutcome = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; 2],
        };
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };

        // Source construction precedes result resize. With zero 3D rows the
        // selected ID is absent, so the typed builder failure preserves the
        // caller's result slots and never enters the lazy row adapter.
        let mut zero_coordinates = CoordinateBlock::default();
        zero_coordinates.conformers_2d.push(
            cosmolkit_model::Conformer2D::new(IDS[0], TWO_D_SENTINEL.to_vec())
                .with_prop("layout", "preserve-this-2d-row"),
        );
        let zero_2d_before = zero_coordinates.conformers_2d.clone();
        let mut zero_results = vec![SENTINEL; 2];
        let mut zero_options = super::SingleConformerOptions::for_conformer(IDS[0]);
        zero_options.vdw_threshold = 10.0;
        zero_options.ignore_interfragment_interactions = false;
        super::uff_public_perf_serial_adapter_probe_start();
        let zero_error = super::optimize_serial_uff_coordinate_block(
            &topology,
            &mut zero_coordinates,
            &mut zero_results,
            super::UffAtomStateRef::SuppliedRows {
                total_valences: &total_valences,
                conjugated_presence: &conjugated,
            },
            &rings,
            &valence,
            &properties,
            &mut Vec::new(),
            zero_options,
        )
        .expect_err("empty block has no selected 3D reference conformer");
        let zero_probe = super::uff_public_perf_serial_adapter_probe_finish();
        assert!(matches!(
            zero_error,
            super::SerialUffOptimizationError::Construction(
                super::AutomaticForceFieldConstructionError::Construction(
                    ForceFieldConstructionError::Builder(
                        UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                            conformer_id: 80
                        }
                    )
                )
            )
        ));
        assert_eq!(zero_probe, Some((0, None)));
        assert_eq!(zero_results, [SENTINEL; 2]);
        assert_eq!(zero_coordinates.conformers_2d, zero_2d_before);

        // One valid 3D row starts with excess caller storage. The block entry
        // truncates to the source conformer count, then assigns the fixed
        // zero-iteration W06 reference energy through the shared stage owner.
        let mut one_coordinates = CoordinateBlock::default();
        one_coordinates.conformers_3d.push(Conformer3D::new(
            IDS[0],
            [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]].to_vec(),
            true,
        ));
        one_coordinates.conformers_2d.push(
            cosmolkit_model::Conformer2D::new(IDS[0], TWO_D_SENTINEL.to_vec())
                .with_prop("layout", "preserve-this-2d-row"),
        );
        let one_2d_before = one_coordinates.conformers_2d.clone();
        let mut one_results = vec![SENTINEL; 2];
        let mut one_options = super::SingleConformerOptions::for_conformer(IDS[0]);
        one_options.max_iterations = 0;
        one_options.vdw_threshold = 10.0;
        one_options.ignore_interfragment_interactions = false;
        cf3d_uff_one_kernel_counts_reset();
        super::uff_public_perf_serial_adapter_probe_start();
        super::optimize_serial_uff_coordinate_block(
            &topology,
            &mut one_coordinates,
            &mut one_results,
            super::UffAtomStateRef::SuppliedRows {
                total_valences: &total_valences,
                conjugated_presence: &conjugated,
            },
            &rings,
            &valence,
            &properties,
            &mut Vec::new(),
            one_options,
        )
        .expect("one fixed W06 row succeeds");
        let one_probe = super::uff_public_perf_serial_adapter_probe_finish();
        assert_eq!(one_probe, Some((1, Some(0))));
        assert_eq!(one_results.len(), 1);
        assert_eq!(one_results[0].status, 1);
        assert!((one_results[0].energy - source_energy(4.0)).abs() <= 1.0e-12);
        assert_eq!(one_coordinates.conformers_2d, one_2d_before);

        // Three rows start with a short result vector. Row0 completes, row1
        // preserves its typed coordinate-count failure, and row2 is never
        // adapted or evaluated. The test-only probe observes zero rows at the
        // stage-owner entry and only the two rows demanded by its loop; this
        // detects eager metadata collection without claiming the optimizer
        // has no allocations elsewhere.
        let mut three_coordinates = CoordinateBlock::default();
        three_coordinates.conformers_3d.extend([
            Conformer3D::new(IDS[0], [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]].to_vec(), true),
            Conformer3D::new(IDS[1], [[0.0, 0.0, 0.0]].to_vec(), true),
            Conformer3D::new(IDS[2], [[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]].to_vec(), true),
        ]);
        three_coordinates.conformers_2d.push(
            cosmolkit_model::Conformer2D::new(IDS[0], TWO_D_SENTINEL.to_vec())
                .with_prop("layout", "preserve-this-2d-row"),
        );
        let three_3d_before = three_coordinates.conformers_3d.clone();
        let three_2d_before = three_coordinates.conformers_2d.clone();
        let mut three_results = vec![SENTINEL];
        let mut three_options = super::SingleConformerOptions::for_conformer(IDS[0]);
        three_options.max_iterations = 0;
        three_options.vdw_threshold = 10.0;
        three_options.ignore_interfragment_interactions = false;
        cf3d_uff_one_kernel_counts_reset();
        super::uff_public_perf_serial_adapter_probe_start();
        let three_error = super::optimize_serial_uff_coordinate_block(
            &topology,
            &mut three_coordinates,
            &mut three_results,
            super::UffAtomStateRef::SuppliedRows {
                total_valences: &total_valences,
                conjugated_presence: &conjugated,
            },
            &rings,
            &valence,
            &properties,
            &mut Vec::new(),
            three_options,
        )
        .expect_err("the second source row has too few atom coordinates");
        let three_probe = super::uff_public_perf_serial_adapter_probe_finish();
        assert!(matches!(
            three_error,
            super::SerialUffOptimizationError::Optimization(
                super::SerialConformerOptimizationError::CoordinateCount {
                    input_index: 1,
                    source: super::SerialCoordinateCountError {
                        conformer_id: 7,
                        atoms: 2,
                        coordinates: 1
                    }
                }
            )
        ));
        assert_eq!(three_probe, Some((2, Some(0))));
        assert_eq!(three_results.len(), 3);
        assert_eq!(three_results[0].status, 1);
        assert!((three_results[0].energy - source_energy(4.0)).abs() <= 1.0e-12);
        assert_eq!(
            &three_results[1..],
            &[OptimizationOutcome {
                status: 0,
                energy: 0.0,
            }; 2]
        );
        assert_eq!(three_coordinates.conformers_3d, three_3d_before);
        assert_eq!(three_coordinates.conformers_2d, three_2d_before);
        assert_eq!(
            crate::kernel::cf3d_uff_one_serial_work_counts().initialize_calls,
            1
        );
    }

    #[cfg(not(target_family = "wasm"))]
    fn run_uff_w06_dispatch_product(prepared: bool) {
        use crate::uff::builder::{ForceFieldConstructionError, UffBuilderError};

        fn assert_copy<T: Copy>() {}
        assert_copy::<super::UffAtomStateRef<'static>>();

        #[derive(Debug)]
        enum DispatchCallError {
            Direct(super::DispatchedUffOptimizationError),
            Prepared(super::super::optimization::UffPreparedOptimizationError),
        }

        const IDS: [usize; 3] = [80, 7, 900];
        const TARGET_DISTANCES: [f64; 3] = [4.0, 5.0, 6.0];
        const REQUESTED: [i32; 5] = [-1, 0, 1, 2, 3];
        const HARDWARE: [u32; 2] = [0, 3];
        const SOURCE_THREAD_COUNTS: [[u32; 5]; 2] = [[1, 1, 1, 2, 3], [2, 3, 1, 2, 3]];
        const SENTINEL: OptimizationOutcome = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        const TWO_D_SENTINEL: [[f64; 2]; 2] = [[17.0, -3.5], [-2.25, 18.125]];
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = fast_find_rings(&topology).expect("fixed W06 topology has ring information");
        let valence = ValenceAssignment {
            explicit_valence: total_valences.to_vec(),
            implicit_hydrogens: vec![0; 2],
        };
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let source_energy = |distance: f64| {
            let ratio = minimum / distance;
            let ratio3 = (ratio * ratio) * ratio;
            let ratio6 = ratio3 * ratio3;
            let ratio12 = ratio6 * ratio6;
            well_depth * (ratio12 - 2.0 * ratio6)
        };
        let mut matrix_calls = 0;

        for row_count in [0, 1, 3] {
            for (hardware_index, observed_hardware) in HARDWARE.into_iter().enumerate() {
                for (request_index, requested_threads) in REQUESTED.into_iter().enumerate() {
                    for threadsafe in [false, true] {
                        let mut coordinates = CoordinateBlock::default();
                        for row in 0..row_count {
                            coordinates.conformers_3d.push(Conformer3D::new(
                                IDS[row],
                                [[0.0, 0.0, 0.0], [TARGET_DISTANCES[row], 0.0, 0.0]].to_vec(),
                                true,
                            ));
                        }
                        coordinates.conformers_2d.push(
                            cosmolkit_model::Conformer2D::new(IDS[0], TWO_D_SENTINEL.to_vec())
                                .with_prop("layout", "preserve-this-2d-row"),
                        );
                        let before_3d = coordinates.conformers_3d.clone();
                        let before_2d = coordinates.conformers_2d.clone();
                        let initial_result_slots = match row_count {
                            0 | 1 => 2,
                            3 => 1,
                            _ => unreachable!("fixed matrix uses 0/1/3 rows"),
                        };
                        let mut results = vec![SENTINEL; initial_result_slots];
                        let mut options = super::SingleConformerOptions::for_conformer(IDS[0]);
                        options.max_iterations = 0;
                        options.vdw_threshold = 10.0;
                        options.ignore_interfragment_interactions = false;

                        cf3d_uff_one_kernel_counts_reset();
                        super::uff_public_perf_serial_adapter_probe_start();
                        let outcome = if prepared {
                            super::super::api::reset_prepare_parameter_query_calls();
                            let outcome =
                                super::super::optimization::optimize_prepared_uff_dispatch(
                                    &topology,
                                    &mut coordinates,
                                    &mut results,
                                    &valence,
                                    &rings,
                                    &properties,
                                    &mut Vec::new(),
                                    options,
                                    requested_threads,
                                    observed_hardware,
                                    threadsafe,
                                )
                                .map_err(DispatchCallError::Prepared);
                            assert_eq!(
                                super::super::api::prepare_parameter_query_calls(),
                                1,
                                "cached input prepares once before dispatch selection"
                            );
                            outcome
                        } else {
                            super::optimize_dispatched_uff_coordinate_block(
                                &topology,
                                &mut coordinates,
                                &mut results,
                                super::UffAtomStateRef::SuppliedRows {
                                    total_valences: &total_valences,
                                    conjugated_presence: &conjugated,
                                },
                                &rings,
                                &valence,
                                &properties,
                                &mut Vec::new(),
                                options,
                                requested_threads,
                                observed_hardware,
                                threadsafe,
                            )
                            .map_err(DispatchCallError::Direct)
                        };
                        let probe = super::uff_public_perf_serial_adapter_probe_finish();
                        matrix_calls += 1;

                        if row_count == 0 {
                            assert!(match (prepared, outcome) {
                                (
                                    false,
                                    Err(DispatchCallError::Direct(
                                        super::DispatchedUffOptimizationError::Construction(
                                            super::AutomaticForceFieldConstructionError::Construction(
                                                ForceFieldConstructionError::Builder(
                                                    UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                                                        conformer_id: 80
                                                    }
                                                )
                                            )
                                        ),
                                    )),
                                ) => true,
                                (
                                    true,
                                    Err(DispatchCallError::Prepared(
                                        super::super::optimization::UffPreparedOptimizationError::Dispatch(
                                            super::DispatchedUffOptimizationError::Construction(
                                                super::AutomaticForceFieldConstructionError::Construction(
                                                    ForceFieldConstructionError::Builder(
                                                        UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                                                            conformer_id: 80
                                                        }
                                                    )
                                                )
                                            )
                                        ),
                                    )),
                                ) => true,
                                _ => false,
                            });
                            assert_eq!(probe, Some((0, None)));
                            assert_eq!(results, vec![SENTINEL; 2]);
                        } else {
                            let expected_source_count = if threadsafe {
                                SOURCE_THREAD_COUNTS[hardware_index][request_index]
                            } else {
                                1
                            };
                            let outcome = match outcome {
                                Ok(outcome) => outcome,
                                Err(error) => {
                                    panic!("fixed W06 rows construct and dispatch: {error:?}")
                                }
                            };
                            if expected_source_count == 1 {
                                assert!(matches!(
                                    outcome,
                                    super::PreparedConformerDispatchOutcome::Serial(Ok(()))
                                ));
                                assert_eq!(probe, Some((row_count, Some(0))));
                            } else {
                                let worker_results = match outcome {
                                    super::PreparedConformerDispatchOutcome::Workers(Ok(
                                        results,
                                    )) => results,
                                    _ => panic!(
                                        "source count {expected_source_count} keeps raw workers"
                                    ),
                                };
                                assert_eq!(worker_results.len(), expected_source_count as usize);
                                assert!(
                                    worker_results
                                        .iter()
                                        .all(|result| matches!(result, Ok(Ok(()))))
                                );
                                assert_eq!(probe, Some((row_count, None)));
                            }
                            assert_eq!(results.len(), row_count);
                            for row in 0..row_count {
                                assert_eq!(results[row].status, 1);
                                assert!(
                                    (results[row].energy - source_energy(TARGET_DISTANCES[row]))
                                        .abs()
                                        <= 1.0e-12,
                                    "fixed source energy row={row}, count={row_count}, requested={requested_threads}, hardware={observed_hardware}, threadsafe={threadsafe}"
                                );
                            }
                        }
                        assert_eq!(coordinates.conformers_3d, before_3d);
                        assert_eq!(coordinates.conformers_2d, before_2d);
                    }
                }
            }
        }
        assert_eq!(matrix_calls, 60);

        // Preserve the worker engine's nested typed failure and every lane's
        // join result. With two source lanes, row1 fails in lane1 while lane0
        // still processes its later row2 and writes that original result slot.
        let mut failure_coordinates = CoordinateBlock::default();
        failure_coordinates.conformers_3d.extend([
            Conformer3D::new(IDS[0], [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]].to_vec(), true),
            Conformer3D::new(IDS[1], [[0.0, 0.0, 0.0]].to_vec(), true),
            Conformer3D::new(IDS[2], [[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]].to_vec(), true),
        ]);
        let mut failure_results = vec![SENTINEL];
        let mut failure_options = super::SingleConformerOptions::for_conformer(IDS[0]);
        failure_options.max_iterations = 0;
        failure_options.vdw_threshold = 10.0;
        failure_options.ignore_interfragment_interactions = false;
        super::uff_public_perf_serial_adapter_probe_start();
        let failure_outcome = if prepared {
            super::super::api::reset_prepare_parameter_query_calls();
            let outcome = super::super::optimization::optimize_prepared_uff_dispatch(
                &topology,
                &mut failure_coordinates,
                &mut failure_results,
                &valence,
                &rings,
                &properties,
                &mut Vec::new(),
                failure_options,
                2,
                0,
                true,
            )
            .map_err(DispatchCallError::Prepared);
            assert_eq!(super::super::api::prepare_parameter_query_calls(), 1);
            outcome.unwrap_or_else(|error| {
                panic!("fixed two-lane prepared worker route keeps raw results: {error:?}")
            })
        } else {
            super::optimize_dispatched_uff_coordinate_block(
                &topology,
                &mut failure_coordinates,
                &mut failure_results,
                super::UffAtomStateRef::SuppliedRows {
                    total_valences: &total_valences,
                    conjugated_presence: &conjugated,
                },
                &rings,
                &valence,
                &properties,
                &mut Vec::new(),
                failure_options,
                2,
                0,
                true,
            )
            .expect("fixed two-lane worker route keeps raw results")
        };
        let failure_probe = super::uff_public_perf_serial_adapter_probe_finish();
        assert_eq!(failure_probe, Some((3, None)));
        let worker_results = match failure_outcome {
            super::PreparedConformerDispatchOutcome::Workers(Ok(results)) => results,
            _ => panic!("two source lanes preserve the raw worker join vector"),
        };
        assert_eq!(worker_results.len(), 2);
        assert!(matches!(&worker_results[0], Ok(Ok(()))));
        assert!(matches!(
            &worker_results[1],
            Ok(Err(
                super::SerialConformerOptimizationError::CoordinateCount {
                    input_index: 1,
                    source: super::SerialCoordinateCountError {
                        conformer_id: 7,
                        atoms: 2,
                        coordinates: 1
                    }
                }
            ))
        ));
        assert_eq!(failure_results.len(), 3);
        assert_eq!(failure_results[0].status, 1);
        assert!((failure_results[0].energy - source_energy(4.0)).abs() <= 1.0e-12);
        assert_eq!(
            failure_results[1],
            OptimizationOutcome {
                status: 0,
                energy: 0.0,
            }
        );
        assert_eq!(failure_results[2].status, 1);
        assert!((failure_results[2].energy - source_energy(6.0)).abs() <= 1.0e-12);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_public_perf_dispatch_resolves_before_serial_metadata_and_keeps_raw_workers() {
        run_uff_w06_dispatch_product(false);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_prepare_p10_cached_dispatch_keeps_the_full_matrix_and_raw_worker_failure() {
        run_uff_w06_dispatch_product(true);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_prepare_p11_actual_adapters_keep_cache_and_topology_error_order() {
        use crate::uff::builder::{PreparedValenceField, UffBuilderError};
        use crate::uff::optimization::UffPreparedOptimizationError;

        const RESULT_SENTINEL: OptimizationOutcome = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        fn require_error<T>(
            result: Result<T, UffPreparedOptimizationError>,
        ) -> UffPreparedOptimizationError {
            match result {
                Err(error) => error,
                Ok(_) => panic!("fixed p11 preparation case must return its typed error"),
            }
        }

        fn assert_preparation_source(
            error: &UffPreparedOptimizationError,
            expected: &UffBuilderError,
        ) {
            let parameter = match error {
                UffPreparedOptimizationError::Preparation(parameter) => parameter,
                _ => panic!("p11 failure must come from cached preparation: {error:?}"),
            };
            let parameter_source =
                std::error::Error::source(error).expect("prepared error keeps its parameter error");
            let stored_parameter: &(dyn std::error::Error + 'static) = parameter;
            assert!(std::ptr::eq(parameter_source, stored_parameter));
            assert_eq!(
                parameter.kind(),
                super::super::api::UffParameterErrorKind::Preparation
            );

            let builder_source = std::error::Error::source(parameter)
                .expect("parameter preparation keeps its builder error");
            let builder = builder_source
                .downcast_ref::<UffBuilderError>()
                .expect("the concrete cache/topology cause remains UffBuilderError");
            assert!(std::ptr::eq(
                builder_source,
                builder as &(dyn std::error::Error + 'static)
            ));
            assert_eq!(builder, expected);

            match expected {
                UffBuilderError::TopologyValidation(expected_topology) => {
                    let topology_source = std::error::Error::source(builder)
                        .expect("topology validation retains its typed source");
                    let topology_error = topology_source
                        .downcast_ref::<cosmolkit_model::TopologyValidationError>()
                        .expect("the stored topology error keeps its concrete type");
                    assert!(std::ptr::eq(
                        topology_source,
                        topology_error as &(dyn std::error::Error + 'static)
                    ));
                    assert_eq!(topology_error, expected_topology);
                    assert!(std::error::Error::source(topology_error).is_none());
                }
                _ => assert!(std::error::Error::source(builder).is_none()),
            }
        }

        let valid_topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let rings = fast_find_rings(&valid_topology)
            .expect("the fixed W06 source topology has ring information");
        let mut invalid_topology = valid_topology;
        invalid_topology.adjacency = cosmolkit_model::AdjacencyList::default();
        let before_topology = invalid_topology.clone();
        let expected_topology_error = cosmolkit_model::TopologyValidationError::AdjacencyMismatch;
        assert_eq!(
            invalid_topology.validate(),
            Err(expected_topology_error.clone())
        );

        let cases = [
            (
                "explicit-row-length",
                ValenceAssignment {
                    explicit_valence: vec![0],
                    implicit_hydrogens: vec![0, 0],
                },
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::Explicit,
                    expected: 2,
                    actual: 1,
                },
            ),
            (
                "explicit-negative-source-cache",
                ValenceAssignment {
                    explicit_valence: vec![-1, 0],
                    implicit_hydrogens: vec![0, 0],
                },
                UffBuilderError::SourceValencePrecondition {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::Explicit,
                    value: -1,
                },
            ),
            (
                "explicit-above-source-cache-range",
                ValenceAssignment {
                    explicit_valence: vec![128, 0],
                    implicit_hydrogens: vec![0, 0],
                },
                UffBuilderError::SourceValenceOutOfRange {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::Explicit,
                    value: 128,
                },
            ),
            (
                "validated-cache-reaches-topology",
                ValenceAssignment {
                    explicit_valence: vec![0, 0],
                    implicit_hydrogens: vec![0, 0],
                },
                UffBuilderError::TopologyValidation(expected_topology_error),
            ),
        ];

        let properties = MoleculeProperties::default();
        let mut calls = 0;
        for (case_name, assignment, expected_error) in &cases {
            for mode in 0..3 {
                let mut coordinates = CoordinateBlock::default();
                coordinates.conformers_3d.push(
                    Conformer3D::new(80, [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]].to_vec(), true)
                        .with_prop("source-row", "p11-preparation-error"),
                );
                coordinates.conformers_2d.push(
                    cosmolkit_model::Conformer2D::new(80, [[17.0, -3.5], [-2.25, 18.125]].to_vec())
                        .with_prop("layout", "p11-preserved"),
                );
                let before_coordinates = coordinates.clone();
                let conformer_rows_address = coordinates.conformers_3d.as_ptr();
                let conformer_coords_address = coordinates.conformers_3d[0].coordinates().as_ptr();
                let layout_rows_address = coordinates.conformers_2d.as_ptr();
                let mut diagnostics = Vec::new();
                let mut results = vec![RESULT_SENTINEL; 2];
                let before_results = results.clone();
                let result_address = results.as_ptr();
                let result_capacity = results.capacity();
                let mut options = SingleConformerOptions::for_conformer(999);
                options.max_iterations = 0;

                super::super::api::reset_prepare_parameter_query_calls();
                super::super::builder::reset_typing_valence_projection_vec_constructions();
                super::super::builder::reset_conjugation_projection_vec_constructions();
                let error = match mode {
                    0 => require_error(super::super::optimization::optimize_prepared_uff_single(
                        &invalid_topology,
                        &mut coordinates,
                        assignment,
                        &rings,
                        &properties,
                        &mut diagnostics,
                        options,
                    )),
                    1 => super::super::optimization::optimize_prepared_uff_serial(
                        &invalid_topology,
                        &mut coordinates,
                        &mut results,
                        assignment,
                        &rings,
                        &properties,
                        &mut diagnostics,
                        options,
                    )
                    .err()
                    .unwrap_or_else(|| {
                        panic!("p11 {case_name} serial entry must return its typed error")
                    }),
                    _ => require_error(super::super::optimization::optimize_prepared_uff_dispatch(
                        &invalid_topology,
                        &mut coordinates,
                        &mut results,
                        assignment,
                        &rings,
                        &properties,
                        &mut diagnostics,
                        options,
                        1,
                        3,
                        false,
                    )),
                };
                calls += 1;

                assert_preparation_source(&error, expected_error);
                assert_eq!(
                    super::super::api::prepare_parameter_query_calls(),
                    1,
                    "one actual cached preparation: case={case_name}, adapter={mode}"
                );
                assert_eq!(
                    super::super::builder::typing_valence_projection_vec_constructions(),
                    0,
                    "cached adapter creates no total-valence projection: case={case_name}, adapter={mode}"
                );
                assert_eq!(
                    super::super::builder::conjugation_projection_vec_constructions(),
                    0,
                    "cached adapter creates no conjugation projection: case={case_name}, adapter={mode}"
                );
                assert_eq!(coordinates, before_coordinates);
                assert_eq!(coordinates.conformers_3d.as_ptr(), conformer_rows_address);
                assert_eq!(
                    coordinates.conformers_3d[0].coordinates().as_ptr(),
                    conformer_coords_address
                );
                assert_eq!(coordinates.conformers_2d.as_ptr(), layout_rows_address);
                assert!(diagnostics.is_empty());
                assert_eq!(results.as_ptr(), result_address);
                assert_eq!(results.capacity(), result_capacity);
                assert_eq!(results, before_results);
                assert_eq!(
                    invalid_topology, before_topology,
                    "invalid source topology remains borrowed: case={case_name}"
                );
            }
        }
        assert_eq!(calls, 12);
    }

    #[test]
    fn uff_one_u04_stage_failures_keep_precedence_and_typed_causes() {
        let mut uninitialized = ForceField::new(1);
        uninitialized.add_contribution(Box::new(FailAtEnergyCall {
            calls: Cell::new(0),
            fail_on: 1,
        }));
        assert_eq!(
            optimize_force_field(&mut uninitialized, 0),
            Err(OptimizationStageError::Initialize(
                ForceFieldKernelError::NoPoints
            ))
        );

        let mut first = [0.0];
        let mut second = [2.0];
        let mut minimization_error = ForceField::new(1);
        minimization_error
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        minimization_error.add_contribution(Box::new(FailAtEnergyCall {
            calls: Cell::new(0),
            fail_on: 1,
        }));
        assert_eq!(
            optimize_force_field(&mut minimization_error, 0),
            Err(OptimizationStageError::Minimize(expected_pair_index_error()))
        );
        assert_eq!(first, [0.0]);
        assert_eq!(second, [2.0]);

        let mut first = [0.0];
        let mut second = [2.0];
        let mut final_energy_error = ForceField::new(1);
        final_energy_error
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        final_energy_error.add_contribution(Box::new(FailAtEnergyCall {
            calls: Cell::new(0),
            fail_on: 2,
        }));
        assert_eq!(
            optimize_force_field(&mut final_energy_error, 0),
            Err(OptimizationStageError::FinalEnergy(
                expected_pair_index_error()
            ))
        );
        assert_eq!(first, [0.0]);
        assert_eq!(second, [2.0]);
    }

    #[test]
    fn uff_one_u04_zero_iteration_returns_status_one_coordinates_and_energy() {
        let mut first = [0.0];
        let mut second = [2.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.add_contribution(Box::new(SquaredDistanceContribution {
            first: 0,
            second: 1,
            scale: 1.0,
        }));

        assert_eq!(
            optimize_force_field(&mut force_field, 0),
            Ok(OptimizationOutcome {
                status: 1,
                energy: 4.0,
            })
        );
        assert_eq!(first, [0.0]);
        assert_eq!(second, [2.0]);
    }

    #[test]
    fn uff_one_u04_empty_contributions_return_source_status_zero_and_energy() {
        let mut position = [1.0, -2.0, 3.0];
        let mut force_field = ForceField::new(1);
        force_field.positions_mut().push(&mut position);

        assert_eq!(
            optimize_force_field(&mut force_field, 20),
            Ok(OptimizationOutcome {
                status: 0,
                energy: 0.0,
            })
        );
        assert_eq!(position, [1.0, -2.0, 3.0]);
    }

    #[test]
    fn uff_one_u04_nonempty_zero_gradient_keeps_source_bad_direction_error() {
        let mut first = [0.0];
        let mut second = [0.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.add_contribution(Box::new(SquaredDistanceContribution {
            first: 0,
            second: 1,
            scale: 1.0,
        }));

        assert_eq!(
            optimize_force_field(&mut force_field, 20),
            Err(OptimizationStageError::Minimize(
                ForceFieldKernelError::OptimizerBadDirection
            ))
        );
        assert_eq!(first, [0.0]);
        assert_eq!(second, [0.0]);
    }

    #[test]
    fn uff_one_u05_options_preserve_source_defaults_with_explicit_id() {
        let options = SingleConformerOptions::for_conformer(7);

        assert_eq!(options.conformer_id, 7);
        assert_eq!(options.max_iterations, 1000);
        assert_eq!(options.vdw_threshold, 10.0);
        assert!(options.ignore_interfragment_interactions);
        assert_eq!(options.max_iterations_as_unsigned(), 1000);
    }

    #[test]
    fn uff_one_u05_signed_iteration_conversion_matches_unsigned_source() {
        for (signed, expected) in [
            (0, 0_u32),
            (23, 23_u32),
            (-1, u32::MAX),
            (-3, u32::MAX - 2),
            (i32::MIN, 1_u32 << 31),
        ] {
            let mut options = SingleConformerOptions::for_conformer(1);
            options.max_iterations = signed;
            assert_eq!(options.max_iterations_as_unsigned(), expected);
        }
    }

    #[test]
    fn uff_one_u05_conformer_ids_remain_explicit_and_noncontiguous() {
        let options = [2_usize, 17, 1000].map(SingleConformerOptions::for_conformer);
        assert_eq!(
            options.map(|option| option.conformer_id),
            [2_usize, 17, 1000]
        );
    }

    #[test]
    fn uff_worker_w06_vdw_reference_and_processed_geometry_matrix() {
        use crate::kernel::Cf3dFragAcceptContributionIdentity as Identity;

        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let threshold = 10.0 * minimum;
        let source_energy = |distance: f64| {
            if distance <= 0.0 || distance > threshold {
                return 0.0;
            }
            let ratio = minimum / distance;
            let ratio3 = ratio * ratio * ratio;
            let ratio6 = ratio3 * ratio3;
            well_depth * (ratio6 * ratio6 - 2.0 * ratio6)
        };

        for reference_distance in [2.0, 40.0] {
            for ignore_interfragment_interactions in [false, true] {
                let reference_has_pair =
                    reference_distance < threshold && !ignore_interfragment_interactions;
                let reference_points = [[0.0, 0.0, 0.0], [reference_distance, 0.0, 0.0]];
                let mut reference = worker_w06_coordinates(&reference_points);
                let source_field = worker_w06_source_field(
                    &topology,
                    &mut reference,
                    &total_valences,
                    10.0,
                    ignore_interfragment_interactions,
                );
                let expected_identity = if reference_has_pair {
                    vec![Identity::Vdw {
                        at1_idx: 0,
                        at2_idx: 1,
                        x_ij: minimum,
                        well_depth,
                        threshold,
                    }]
                } else {
                    Vec::new()
                };
                let identities =
                    crate::kernel::cf3d_frag_accept_contribution_identities(&source_field);
                assert_eq!(identities, expected_identity);

                let expected_energies = if reference_has_pair {
                    [source_energy(2.0), source_energy(20.0)]
                } else {
                    [0.0, 0.0]
                };
                let worker_field = source_field.copy();
                assert_eq!(
                    crate::kernel::cf3d_frag_accept_contribution_identities(&worker_field),
                    identities
                );
                drop(source_field);

                let mut near_row = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
                let mut far_row = [[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]];
                let mut conformers = [
                    super::SerialConformer {
                        id: 80,
                        positions: &mut near_row,
                    },
                    super::SerialConformer {
                        id: 2,
                        positions: &mut far_row,
                    },
                ];
                let sentinel = OptimizationOutcome {
                    status: -9,
                    energy: 123.0,
                };
                let mut results = [sentinel; 2];
                super::optimize_worker_conformers(
                    worker_field,
                    &mut conformers,
                    &mut results,
                    2,
                    0,
                    std::num::NonZeroU32::new(1).unwrap(),
                    0,
                )
                .expect("the copied source VDW field evaluates both reordered rows");

                let expected_status = if reference_has_pair { 1 } else { 0 };
                for (result, expected_energy) in results.iter().zip(expected_energies) {
                    assert_eq!(result.status, expected_status);
                    assert!((result.energy - expected_energy).abs() <= 1.0e-12);
                }
                assert_eq!(near_row, [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
                assert_eq!(far_row, [[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]]);
                assert_eq!(reference.conformers_3d[0].coordinates(), reference_points);
            }
        }
    }

    #[test]
    fn uff_worker_w06_tbp_reference_terms_survive_copy_and_reordered_rows() {
        use crate::kernel::Cf3dFragAcceptContributionIdentity as Identity;

        const B20_TERMS: [(u32, u32, u32, u32); 10] = [
            (1, 0, 5, 2),
            (2, 0, 3, 3),
            (2, 0, 4, 3),
            (3, 0, 4, 3),
            (1, 0, 2, 0),
            (1, 0, 3, 0),
            (1, 0, 4, 0),
            (5, 0, 2, 0),
            (5, 0, 3, 0),
            (5, 0, 4, 0),
        ];
        const ALTERNATE_TERMS: [(u32, u32, u32, u32); 10] = [
            (2, 0, 4, 2),
            (1, 0, 3, 3),
            (1, 0, 5, 3),
            (3, 0, 5, 3),
            (2, 0, 1, 0),
            (2, 0, 3, 0),
            (2, 0, 5, 0),
            (4, 0, 1, 0),
            (4, 0, 3, 0),
            (4, 0, 5, 0),
        ];
        let diagonal = 1.0 / 3.0_f64.sqrt();
        let b20 = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-0.6, 0.8, 0.0],
            [0.0, 0.6, 0.8],
            [0.0, -0.8, 0.6],
            [-0.8, 0.0, -0.6],
        ];
        let alternate = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
            [diagonal, diagonal, diagonal],
        ];
        let topology = worker_w06_topology(
            std::iter::once(worker_w06_atom(0, 15, Hybridization::Sp3d, 0))
                .chain((1..6).map(|row| worker_w06_atom(row, 6, Hybridization::Sp3, 0)))
                .collect(),
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
                (0, 5, BondOrder::Single),
            ],
        );
        let total_valences = [5, 1, 1, 1, 1, 1];

        for (reference_label, reference_points, expected_terms) in [
            ("B20", b20, B20_TERMS),
            ("alternate", alternate, ALTERNATE_TERMS),
        ] {
            let mut reference = worker_w06_coordinates(&reference_points);
            let mut source_field =
                worker_w06_source_field(&topology, &mut reference, &total_valences, 0.0, false);
            let identities = crate::kernel::cf3d_frag_accept_contribution_identities(&source_field);
            assert_eq!(
                identities.len(),
                15,
                "{reference_label}: five bonds then ten angles"
            );
            for (bond_row, identity) in identities[..5].iter().enumerate() {
                assert!(matches!(
                    identity,
                    Identity::BondStretch {
                        end1_idx: 0,
                        end2_idx,
                        ..
                    } if *end2_idx == bond_row as u32 + 1
                ));
            }
            let actual_terms = identities[5..]
                .iter()
                .map(|identity| match identity {
                    Identity::AngleBend {
                        at1_idx,
                        at2_idx,
                        at3_idx,
                        order,
                        ..
                    } => (*at1_idx, *at2_idx, *at3_idx, *order),
                    other => panic!("{reference_label}: unexpected TBP source term {other:?}"),
                })
                .collect::<Vec<_>>();
            assert_eq!(
                actual_terms, expected_terms,
                "{reference_label}: source term order"
            );

            source_field
                .initialize()
                .expect("fixed source TBP field initializes before reference evaluation");
            let expected_energy = |points: &[[f64; 3]; 6], field: &mut ForceField<'_>| {
                let flat = points.iter().flatten().copied().collect::<Vec<_>>();
                crate::kernel::cf3d_bld_b05_calc_energy(field, &flat)
                    .expect("fixed source-listed TBP terms evaluate")
            };
            let expected_b20 = expected_energy(&b20, &mut source_field);
            let expected_alternate = expected_energy(&alternate, &mut source_field);
            if reference_label == "B20" {
                assert!(
                    (expected_b20 - 1142.363_631_778_066_3).abs() <= 1.0e-7,
                    "the pinned U10 B20 zero-iteration energy is an independent source oracle"
                );
            }
            for original_first in [false, true] {
                let (
                    first_points,
                    second_points,
                    first_expected,
                    second_expected,
                    first_id,
                    second_id,
                ) = if original_first {
                    (b20, alternate, expected_b20, expected_alternate, 2, 80)
                } else {
                    (alternate, b20, expected_alternate, expected_b20, 80, 2)
                };
                let worker_field = source_field.copy();
                assert_eq!(
                    crate::kernel::cf3d_frag_accept_contribution_identities(&worker_field),
                    identities,
                    "{reference_label}: copied worker retains reference term identities"
                );
                let mut first_row = first_points;
                let mut second_row = second_points;
                let mut conformers = [
                    super::SerialConformer {
                        id: first_id,
                        positions: &mut first_row,
                    },
                    super::SerialConformer {
                        id: second_id,
                        positions: &mut second_row,
                    },
                ];
                let sentinel = OptimizationOutcome {
                    status: -9,
                    energy: 123.0,
                };
                let mut results = [sentinel; 2];
                super::optimize_worker_conformers(
                    worker_field,
                    &mut conformers,
                    &mut results,
                    6,
                    0,
                    std::num::NonZeroU32::new(1).unwrap(),
                    0,
                )
                .expect("the copied source TBP field evaluates reordered rows");
                assert_eq!(results[0].status, 1, "{reference_label}: first row");
                assert_eq!(results[1].status, 1, "{reference_label}: second row");
                assert_eq!(results[0].energy.to_bits(), first_expected.to_bits());
                assert_eq!(results[1].energy.to_bits(), second_expected.to_bits());
                assert_eq!(first_row, first_points);
                assert_eq!(second_row, second_points);
            }
            drop(source_field);
            assert_eq!(reference.conformers_3d[0].coordinates(), reference_points);
        }
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_thread_t05_one_scoped_worker_matches_w05_and_w06_references() {
        let row_ids = [900, 2, 80];
        let w05_rows = [
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
        ];
        let mut w05_source_field = ForceField::new(3);
        w05_source_field.add_contribution(Box::new(CoordinateNormContribution));
        assert_scoped_worker_reference_matrix(
            &w05_source_field,
            &row_ids,
            &w05_rows,
            &[1, 1, 1],
            &[4.0, 9.0, 16.0],
            2,
            0.0,
        );

        use crate::kernel::Cf3dFragAcceptContributionIdentity as Identity;
        let topology = worker_w06_topology(
            vec![
                worker_w06_atom(0, 11, Hybridization::Unspecified, 1),
                worker_w06_atom(1, 17, Hybridization::Unspecified, -1),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let threshold = 10.0 * minimum;
        let mut reference = worker_w06_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let w06_source_field =
            worker_w06_source_field(&topology, &mut reference, &total_valences, 10.0, false);
        assert_eq!(
            crate::kernel::cf3d_frag_accept_contribution_identities(&w06_source_field),
            [Identity::Vdw {
                at1_idx: 0,
                at2_idx: 1,
                x_ij: minimum,
                well_depth,
                threshold,
            }]
        );
        let source_energy = |distance: f64| {
            if distance <= 0.0 || distance > threshold {
                return 0.0;
            }
            let ratio = minimum / distance;
            let ratio3 = ratio * ratio * ratio;
            let ratio6 = ratio3 * ratio3;
            well_depth * (ratio6 * ratio6 - 2.0 * ratio6)
        };
        let distances = [2.0, 5.0, 20.0];
        let w06_rows = distances
            .map(|distance| vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
            .to_vec();
        let w06_energies = distances.map(source_energy);
        assert_scoped_worker_reference_matrix(
            &w06_source_field,
            &row_ids,
            &w06_rows,
            &[1, 1, 1],
            &w06_energies,
            2,
            1.0e-12,
        );
        drop(w06_source_field);
        assert_eq!(
            reference.conformers_3d[0].coordinates(),
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
        );
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_thread_t07_scoped_worker_stage_failures_preserve_rows_and_results() {
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let lane_count = std::num::NonZeroU32::new(2).expect("two lanes are nonzero");

        let initialize_ids = [80, 7, 900];
        let mut initialize_rows = vec![Vec::<[f64; 3]>::new(); initialize_ids.len()];
        let initialize_before = initialize_rows.clone();
        let mut initialize_conformers = initialize_rows
            .iter_mut()
            .zip(initialize_ids)
            .map(|(positions, id)| super::SerialConformer {
                id,
                positions: positions.as_mut_slice(),
            })
            .collect::<Vec<_>>();
        let mut initialize_results = vec![sentinel; initialize_ids.len()];
        let (initialize_result, observed_ids) = run_worker_on_one_scoped_thread(
            ForceField::new(3),
            &mut initialize_conformers,
            &mut initialize_results,
            0,
            1,
            lane_count,
            0,
        );
        assert_eq!(observed_ids, initialize_ids);
        assert_eq!(
            initialize_result,
            Err(super::SerialConformerOptimizationError::Optimization {
                input_index: 1,
                conformer_id: 7,
                source: super::OptimizationStageError::Initialize(ForceFieldKernelError::NoPoints),
            })
        );
        drop(initialize_conformers);
        assert_eq!(initialize_rows, initialize_before);
        assert_eq!(initialize_results, vec![sentinel; initialize_ids.len()]);

        let ids = [80, 7, 900, 31, 502, 8];
        let mut minimize_rows = vec![
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[10.0, 0.0, 0.0], [12.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            vec![[20.0, 0.0, 0.0], [23.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
            vec![[30.0, 0.0, 0.0], [34.0, 0.0, 0.0]],
        ];
        let minimize_before = minimize_rows.clone();
        let mut minimize_conformers = minimize_rows
            .iter_mut()
            .zip(ids)
            .map(|(positions, id)| super::SerialConformer {
                id,
                positions: positions.as_mut_slice(),
            })
            .collect::<Vec<_>>();
        let mut minimize_results = vec![sentinel; ids.len()];
        let mut minimize_field = ForceField::new(3);
        minimize_field.add_contribution(Box::new(FailAtEnergyCallWithCoordinateNorm {
            calls: Cell::new(0),
            fail_on: 3,
        }));
        let (minimize_result, observed_ids) = run_worker_on_one_scoped_thread(
            minimize_field,
            &mut minimize_conformers,
            &mut minimize_results,
            2,
            0,
            lane_count,
            0,
        );
        assert_eq!(observed_ids, ids);
        assert_eq!(
            minimize_result,
            Err(super::SerialConformerOptimizationError::Optimization {
                input_index: 2,
                conformer_id: 900,
                source: super::OptimizationStageError::Minimize(expected_pair_index_error()),
            })
        );
        drop(minimize_conformers);
        let mut expected_minimize_results = vec![sentinel; ids.len()];
        expected_minimize_results[0] = OptimizationOutcome {
            status: 1,
            energy: 4.0,
        };
        assert_eq!(minimize_results, expected_minimize_results);
        assert_eq!(minimize_rows, minimize_before);

        let mut final_energy_rows = vec![
            vec![[1.0, 2.0, 0.5], [3.0, 4.0, 0.25]],
            vec![[11.0, 12.0, 0.5], [13.0, 14.0, 0.25]],
            vec![[2.0, 3.0, 0.75], [4.0, 5.0, 0.5]],
            vec![[21.0, 22.0, 0.5], [23.0, 24.0, 0.25]],
            vec![[5.0, 6.0, 0.75], [7.0, 8.0, 0.5]],
            vec![[31.0, 32.0, 0.5], [33.0, 34.0, 0.25]],
        ];
        let final_energy_before = final_energy_rows.clone();
        let mut final_energy_conformers = final_energy_rows
            .iter_mut()
            .zip(ids)
            .map(|(positions, id)| super::SerialConformer {
                id,
                positions: positions.as_mut_slice(),
            })
            .collect::<Vec<_>>();
        let mut final_energy_results = vec![sentinel; ids.len()];
        let mut final_energy_field = ForceField::new(3);
        final_energy_field.add_contribution(Box::new(FailAtEnergyCallWithCoordinateNorm {
            calls: Cell::new(0),
            fail_on: 6,
        }));
        let (final_energy_result, observed_ids) = run_worker_on_one_scoped_thread(
            final_energy_field,
            &mut final_energy_conformers,
            &mut final_energy_results,
            2,
            0,
            lane_count,
            1,
        );
        assert_eq!(observed_ids, ids);
        assert_eq!(
            final_energy_result,
            Err(super::SerialConformerOptimizationError::Optimization {
                input_index: 2,
                conformer_id: 900,
                source: super::OptimizationStageError::FinalEnergy(expected_pair_index_error()),
            })
        );
        drop(final_energy_conformers);
        let completed_energy = final_energy_rows[0]
            .iter()
            .flatten()
            .map(|coordinate| coordinate * coordinate)
            .sum::<f64>();
        assert_eq!(
            final_energy_results[0],
            OptimizationOutcome {
                status: 1,
                energy: completed_energy,
            }
        );
        assert_eq!(final_energy_results[2], sentinel);
        assert!(final_energy_rows[0] != final_energy_before[0]);
        assert!(final_energy_rows[2] != final_energy_before[2]);
        for input_index in 0..ids.len() {
            if input_index != 0 && input_index != 2 {
                assert_eq!(
                    final_energy_rows[input_index],
                    final_energy_before[input_index]
                );
                assert_eq!(final_energy_results[input_index], sentinel);
            }
        }
    }

    #[test]
    fn uff_worker_w07_copy_resets_fixed_points_and_cache_without_touching_source() {
        // ForceField.cpp:159-170 copies terms and point count, while omitted
        // positions, fixed points, distance storage and matrix size default empty.
        let fixed_point_sets: [&[i32]; 4] = [&[], &[0], &[1], &[0, 1]];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        for fixed_points in fixed_point_sets {
            for cache_prefilled in [false, true] {
                for selected_lane in [false, true] {
                    let mut reference_positions = [[1.0, 2.0, 0.5], [3.0, 4.0, 0.25]];
                    let reference_before = reference_positions;
                    let mut source_field = ForceField::new(3);
                    source_field
                        .positions_mut()
                        .extend(reference_positions.iter_mut().map(|row| &mut row[..]));
                    source_field.add_contribution(Box::new(CoordinateNormContribution));
                    crate::kernel::cf3d_uff_one_set_fixed_points_for_serial_test(
                        &mut source_field,
                        fixed_points,
                    );
                    if cache_prefilled {
                        source_field
                            .initialize()
                            .expect("fixed source field initializes for cache setup");
                        let cached_distance =
                            crate::kernel::cf3d_uff_one_prefill_distance_cache_for_worker_test(
                                &mut source_field,
                            )
                            .expect("source distance fills its initialized cache");
                        assert!(cached_distance > 0.0);
                    }

                    let source_state =
                        crate::kernel::cf3d_uff_one_copy_state_for_worker_test(&source_field);
                    assert_eq!(source_state.0, cache_prefilled);
                    assert_eq!(source_state.1, 2);
                    assert_eq!(source_state.2.len(), if cache_prefilled { 3 } else { 0 });
                    assert_eq!(source_state.3, cache_prefilled);
                    assert_eq!(source_state.4, fixed_points);
                    assert_eq!(source_state.5, if cache_prefilled { 2 } else { 0 });
                    assert_eq!(source_state.6, if cache_prefilled { 3 } else { 0 });
                    assert_eq!(source_state.7, 1);

                    let worker_field = source_field.copy();
                    let copied_state =
                        crate::kernel::cf3d_uff_one_copy_state_for_worker_test(&worker_field);
                    assert_eq!(
                        copied_state,
                        (
                            false,
                            0,
                            Vec::new(),
                            false,
                            Vec::new(),
                            if cache_prefilled { 2 } else { 0 },
                            0,
                            1,
                        )
                    );

                    let mut selected_row = [[2.0, 1.0, 0.5], [4.0, 3.0, 0.25]];
                    let mut skipped_row = [[10.0, 1.0, 2.0], [14.0, 5.0, 3.0]];
                    let selected_before = selected_row;
                    let skipped_before = skipped_row;
                    let mut conformers = [
                        super::SerialConformer {
                            id: 80,
                            positions: &mut selected_row,
                        },
                        super::SerialConformer {
                            id: 2,
                            positions: &mut skipped_row,
                        },
                    ];
                    let mut results = [sentinel; 2];
                    cf3d_uff_one_kernel_counts_reset();
                    super::optimize_worker_conformers(
                        worker_field,
                        &mut conformers,
                        &mut results,
                        2,
                        if selected_lane { 0 } else { 2 },
                        std::num::NonZeroU32::new(2).unwrap(),
                        1,
                    )
                    .expect("worker copy runs or skips without using source fixed/cache state");
                    assert_eq!(
                        cf3d_uff_one_serial_work_counts().initialize_calls,
                        usize::from(selected_lane)
                    );
                    assert_eq!(
                        crate::kernel::cf3d_uff_one_copy_state_for_worker_test(&source_field),
                        source_state
                    );
                    drop(conformers);

                    if selected_lane {
                        assert_ne!(
                            selected_row[0][0].to_bits(),
                            selected_before[0][0].to_bits()
                        );
                        assert_ne!(
                            selected_row[1][0].to_bits(),
                            selected_before[1][0].to_bits()
                        );
                        assert_eq!(skipped_row, skipped_before);
                        assert!(matches!(results[0].status, 0 | 1));
                        let fixed_source_energy = selected_row
                            .iter()
                            .flatten()
                            .map(|coordinate| coordinate * coordinate)
                            .sum::<f64>();
                        assert!((results[0].energy - fixed_source_energy).abs() <= 1.0e-12);
                        assert_eq!(results[1], sentinel);
                    } else {
                        assert_eq!(selected_row, selected_before);
                        assert_eq!(skipped_row, skipped_before);
                        assert_eq!(results, [sentinel; 2]);
                    }

                    drop(source_field);
                    assert_eq!(reference_positions, reference_before);
                }
            }
        }
    }

    #[test]
    fn uff_worker_w08_copies_terms_once_and_reuses_position_buffer() {
        // ForceField.cpp:159-170 copies each contribution once; the worker
        // reserves its source-sized position vector once before row traversal.
        let mut reference_positions = [[0.5, 1.0, 1.5], [2.0, 2.5, 3.0]];
        let reference_before = reference_positions;
        let mut source_field = ForceField::new(3);
        source_field
            .positions_mut()
            .extend(reference_positions.iter_mut().map(|row| &mut row[..]));
        source_field.add_contribution(Box::new(CoordinateNormContribution));
        source_field.add_contribution(Box::new(SquaredDistanceContribution {
            first: 0,
            second: 1,
            scale: 0.5,
        }));
        source_field
            .initialize()
            .expect("source field initializes before worker copy");

        let mut first_row = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let mut second_row = [[2.0, 0.0, 1.0], [3.0, 4.0, 5.0]];
        let mut third_row = [[8.0, 2.0, 1.0], [4.0, 3.0, 2.0]];
        let rows_before = [first_row, second_row, third_row];
        let mut conformers = [
            super::SerialConformer {
                id: 73,
                positions: &mut first_row,
            },
            super::SerialConformer {
                id: 4,
                positions: &mut second_row,
            },
            super::SerialConformer {
                id: 800,
                positions: &mut third_row,
            },
        ];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut results = [sentinel; 3];
        cf3d_uff_one_kernel_counts_reset();
        let worker_field = source_field.copy();
        super::optimize_worker_conformers(
            worker_field,
            &mut conformers,
            &mut results,
            2,
            0,
            std::num::NonZeroU32::new(1).unwrap(),
            0,
        )
        .expect("one copied worker evaluates all rows");
        drop(conformers);
        let kernel_counts = cf3d_uff_one_kernel_counts();
        let buffer_counts = cf3d_uff_one_serial_work_counts();
        assert_eq!(kernel_counts, (1, 0, 2));
        assert_eq!(buffer_counts.position_buffer_growths, 1);
        assert_eq!(buffer_counts.initialize_calls, 3);
        assert_eq!(
            buffer_counts.initial_position_buffer_address,
            buffer_counts.final_position_buffer_address
        );
        assert_eq!(
            buffer_counts.initial_position_buffer_capacity,
            buffer_counts.final_position_buffer_capacity
        );
        assert!(buffer_counts.final_position_buffer_capacity >= 2);
        assert_eq!(first_row, rows_before[0]);
        assert_eq!(second_row, rows_before[1]);
        assert_eq!(third_row, rows_before[2]);

        for (row_index, points) in rows_before.iter().enumerate() {
            let coordinate_norm = points
                .iter()
                .flatten()
                .map(|value| value * value)
                .sum::<f64>();
            let pair_distance_squared = (0..3)
                .map(|axis| (points[0][axis] - points[1][axis]).powi(2))
                .sum::<f64>();
            let expected_energy = coordinate_norm + 0.5 * pair_distance_squared;
            assert_eq!(results[row_index].status, 1);
            assert!((results[row_index].energy - expected_energy).abs() <= 1.0e-12);
        }
        drop(source_field);
        assert_eq!(reference_positions, reference_before);
    }

    #[test]
    fn uff_worker_w10_source_lanes_compose_in_input_order() {
        let conformer_ids = [80, 2, 900, 17, 64];
        let initial_rows = [
            [[1.0, 0.0, 0.5]],
            [[2.0, -1.0, 0.5]],
            [[3.0, -2.0, 0.5]],
            [[4.0, -3.0, 0.5]],
            [[5.0, -4.0, 0.5]],
        ];
        let expected_energies = [1.25, 5.25, 13.25, 25.25, 41.25];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut source_field = ForceField::new(3);
        source_field.add_contribution(Box::new(CoordinateNormContribution));

        for row_count in 0..=conformer_ids.len() {
            for lane_count_value in [1_u32, 2, 3, 4, 7] {
                let expected_by_lane: &[&[usize]] = match lane_count_value {
                    1 => &[&[0, 1, 2, 3, 4]],
                    2 => &[&[0, 2, 4], &[1, 3]],
                    3 => &[&[0, 3], &[1, 4], &[2]],
                    4 => &[&[0, 4], &[1], &[2], &[3]],
                    7 => &[&[0], &[1], &[2], &[3], &[4], &[], &[]],
                    _ => unreachable!(),
                };
                let lane_count = std::num::NonZeroU32::new(lane_count_value).unwrap();
                let mut coordinate_rows = initial_rows[..row_count].to_vec();
                let coordinates_before = coordinate_rows.clone();
                let mut results = vec![sentinel; row_count + 2];
                let mut assigned = vec![false; row_count];
                let mut total_initializations = 0;

                for thread_idx in 0..lane_count_value {
                    let expected_indices = expected_by_lane[thread_idx as usize]
                        .iter()
                        .copied()
                        .filter(|input_index| *input_index < row_count)
                        .collect::<Vec<_>>();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(input_index, positions)| super::SerialConformer {
                            id: conformer_ids[input_index],
                            positions: &mut positions[..],
                        })
                        .collect::<Vec<_>>();
                    assert_eq!(
                        conformers.iter().map(|row| row.id).collect::<Vec<_>>(),
                        conformer_ids[..row_count]
                    );

                    cf3d_uff_one_kernel_counts_reset();
                    super::optimize_worker_conformers(
                        source_field.copy(),
                        &mut conformers,
                        &mut results,
                        1,
                        thread_idx,
                        lane_count,
                        0,
                    )
                    .expect("fixed valid lane composition succeeds");
                    assert_eq!(
                        conformers.iter().map(|row| row.id).collect::<Vec<_>>(),
                        conformer_ids[..row_count]
                    );
                    drop(conformers);

                    let initialize_calls = cf3d_uff_one_serial_work_counts().initialize_calls;
                    assert_eq!(initialize_calls, expected_indices.len());
                    total_initializations += initialize_calls;
                    for input_index in 0..row_count {
                        if expected_indices.contains(&input_index) {
                            assigned[input_index] = true;
                        }
                        if assigned[input_index] {
                            assert_eq!(results[input_index].status, 1);
                            assert!(
                                (results[input_index].energy - expected_energies[input_index])
                                    .abs()
                                    <= 1.0e-12,
                                "rows={row_count}, lanes={lane_count_value}, lane={thread_idx}, input={input_index}",
                            );
                        } else {
                            assert_eq!(results[input_index], sentinel);
                        }
                    }
                    assert_eq!(&results[row_count..], &[sentinel; 2]);
                    assert_eq!(coordinate_rows, coordinates_before);
                }

                assert_eq!(total_initializations, row_count);
                assert!(assigned.iter().all(|was_assigned| *was_assigned));
                for invalid_thread_idx in [lane_count_value, lane_count_value + 1, u32::MAX] {
                    let results_before = results.clone();
                    let mut conformers = coordinate_rows
                        .iter_mut()
                        .enumerate()
                        .map(|(input_index, positions)| super::SerialConformer {
                            id: conformer_ids[input_index],
                            positions: &mut positions[..],
                        })
                        .collect::<Vec<_>>();
                    cf3d_uff_one_kernel_counts_reset();
                    super::optimize_worker_conformers(
                        source_field.copy(),
                        &mut conformers,
                        &mut results,
                        1,
                        invalid_thread_idx,
                        lane_count,
                        0,
                    )
                    .expect("out-of-range worker lane selects no rows");
                    drop(conformers);
                    assert_eq!(cf3d_uff_one_serial_work_counts().initialize_calls, 0);
                    assert_eq!(results, results_before);
                    assert_eq!(coordinate_rows, coordinates_before);
                }
            }
        }
    }
}
