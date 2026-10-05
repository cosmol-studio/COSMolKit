//! Private force-field storage and evaluation boundary.

#[cfg(test)]
std::thread_local! {
    static UFF_ONE_KERNEL_COUNTS: std::cell::Cell<(usize, usize, usize)> =
        const { std::cell::Cell::new((0, 0, 0)) };
}

#[cfg(test)]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) struct UffOneSerialWorkCounts {
    pub(super) position_buffer_growths: usize,
    pub(super) initialize_calls: usize,
    pub(super) initial_position_buffer_address: usize,
    pub(super) final_position_buffer_address: usize,
    pub(super) initial_position_buffer_capacity: usize,
    pub(super) final_position_buffer_capacity: usize,
}

#[cfg(test)]
std::thread_local! {
    static UFF_ONE_SERIAL_WORK_COUNTS: std::cell::Cell<UffOneSerialWorkCounts> =
        const { std::cell::Cell::new(UffOneSerialWorkCounts {
            position_buffer_growths: 0,
            initialize_calls: 0,
            initial_position_buffer_address: 0,
            final_position_buffer_address: 0,
            initial_position_buffer_capacity: 0,
            final_position_buffer_capacity: 0,
        }) };
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_kernel_counts_reset() {
    UFF_ONE_KERNEL_COUNTS.with(|counts| counts.set((0, 0, 0)));
    cf3d_uff_one_serial_work_counts_reset();
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_kernel_counts() -> (usize, usize, usize) {
    UFF_ONE_KERNEL_COUNTS.with(std::cell::Cell::get)
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_serial_work_counts_reset() {
    UFF_ONE_SERIAL_WORK_COUNTS.with(|counts| counts.set(UffOneSerialWorkCounts::default()));
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_serial_work_counts() -> UffOneSerialWorkCounts {
    UFF_ONE_SERIAL_WORK_COUNTS.with(std::cell::Cell::get)
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_note_serial_position_buffer_start(address: usize, capacity: usize) {
    UFF_ONE_SERIAL_WORK_COUNTS.with(|counts| {
        let mut current = counts.get();
        current.initial_position_buffer_address = address;
        current.initial_position_buffer_capacity = capacity;
        counts.set(current);
    });
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_note_serial_position_buffer_growth() {
    UFF_ONE_SERIAL_WORK_COUNTS.with(|counts| {
        let mut current = counts.get();
        current.position_buffer_growths += 1;
        counts.set(current);
    });
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_note_serial_position_buffer_end(address: usize, capacity: usize) {
    UFF_ONE_SERIAL_WORK_COUNTS.with(|counts| {
        let mut current = counts.get();
        current.final_position_buffer_address = address;
        current.final_position_buffer_capacity = capacity;
        counts.set(current);
    });
}

#[cfg(test)]
fn record_uff_one_initialize() {
    UFF_ONE_SERIAL_WORK_COUNTS.with(|counts| {
        let mut current = counts.get();
        current.initialize_calls += 1;
        counts.set(current);
    });
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_set_fixed_points_for_serial_test(
    field: &mut ForceField<'_>,
    fixed_points: &[i32],
) {
    field.fixed_points.extend_from_slice(fixed_points);
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_prefill_distance_cache_for_worker_test(
    field: &mut ForceField<'_>,
) -> Result<f64, ForceFieldKernelError> {
    field.distance(0, 1, None)
}

#[cfg(test)]
pub(super) fn cf3d_uff_one_copy_state_for_worker_test(
    field: &ForceField<'_>,
) -> (bool, usize, Vec<f64>, bool, Vec<i32>, u32, usize, usize) {
    (
        field.initialized,
        field.positions.len(),
        field.distance_matrix.clone(),
        field.matrix_allocated,
        field.fixed_points.clone(),
        field.num_points,
        field.matrix_size as usize,
        field.contributions.len(),
    )
}

#[cfg(test)]
fn record_uff_one_field() {
    UFF_ONE_KERNEL_COUNTS.with(|counts| {
        let (fields, terms, copies) = counts.get();
        counts.set((fields + 1, terms, copies));
    });
}

#[cfg(test)]
fn record_uff_one_term() {
    UFF_ONE_KERNEL_COUNTS.with(|counts| {
        let (fields, terms, copies) = counts.get();
        counts.set((fields, terms + 1, copies));
    });
}

#[cfg(test)]
fn record_uff_one_term_copy() {
    UFF_ONE_KERNEL_COUNTS.with(|counts| {
        let (fields, terms, copies) = counts.get();
        counts.set((fields, terms, copies + 1));
    });
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ForceFieldIndexArgument {
    I,
    J,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum BondIndexArgument {
    First,
    Second,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum AngleIndexArgument {
    First,
    Second,
    Third,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum TorsionIndexArgument {
    First,
    Second,
    Third,
    Fourth,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum AngleRangeBound {
    Minimum,
    Maximum,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ForceFieldKernelError {
    NoPoints,
    NoDistanceMatrix,
    MatrixSizeMismatch,
    NotInitialized,
    BadBounds,
    BadBondOrder,
    PositionCountMismatch {
        num_points: u32,
        position_count: usize,
    },
    IndexOutOfRange {
        argument: ForceFieldIndexArgument,
        index: u32,
        upper_bound: u32,
    },
    BondIndexOutOfRange {
        argument: BondIndexArgument,
        index: u32,
        upper_bound: usize,
    },
    AngleDegeneratePoints,
    AngleIndexOutOfRange {
        argument: AngleIndexArgument,
        index: u32,
        upper_bound: usize,
    },
    AngleBadOrder {
        order: u32,
    },
    TorsionIndexOutOfRange {
        argument: TorsionIndexArgument,
        index: u32,
        upper_bound: usize,
    },
    TorsionDegeneratePoints,
    TorsionBadHybridizations,
    TorsionBadOrder {
        order: u32,
    },
    AngleOutOfRange {
        bound: AngleRangeBound,
    },
    AngleBoundsOrder,
    PackedAngleBoundsOrder,
    TorsionBoundsOrder,
    BadIndex,
    OptimizerBadTolerance,
    OptimizerBadDirection,
    BadFixedPoint {
        index: i32,
        upper_bound: u32,
    },
    TransferPostcondition,
    // Proposed pinned OopBend addTerm nullable-parameter precondition.
    OopParametersMissing,
    // Proposed Native safety beyond source-defined nullable parameter/address inputs.
    // ROOT approval required; these are not C++ preconditions or Unsupported.
    TorsionParametersOutsideDefinedSource,
    TorsionAddressOutsideDefinedSource {
        index: i16,
        buffer_len: usize,
    },
    // Proposed Native safety policy for C++-undefined signed/buffer addresses;
    // ROOT authorized this Native safety policy; no C++ parity or Unsupported claim.
    OopAddressOutsideDefinedSource {
        index: i32,
        buffer_len: usize,
    },
}

impl ForceFieldKernelError {
    const fn source_category(self) -> &'static str {
        match self {
            Self::NoPoints
            | Self::NoDistanceMatrix
            | Self::MatrixSizeMismatch
            | Self::NotInitialized
            | Self::BadBounds
            | Self::BadBondOrder
            | Self::PositionCountMismatch { .. }
            | Self::OptimizerBadTolerance
            | Self::AngleBoundsOrder
            | Self::PackedAngleBoundsOrder
            | Self::AngleDegeneratePoints
            | Self::AngleBadOrder { .. }
            | Self::TorsionDegeneratePoints
            | Self::TorsionBadHybridizations
            | Self::TorsionBadOrder { .. }
            | Self::TorsionBoundsOrder
            | Self::OopParametersMissing => "Pre-condition Violation",
            Self::IndexOutOfRange { .. }
            | Self::BondIndexOutOfRange { .. }
            | Self::AngleIndexOutOfRange { .. }
            | Self::TorsionIndexOutOfRange { .. }
            | Self::AngleOutOfRange { .. } => "Range Error",
            Self::BadIndex
            | Self::BadFixedPoint { .. }
            | Self::OptimizerBadDirection
            | Self::OopAddressOutsideDefinedSource { .. }
            | Self::TorsionParametersOutsideDefinedSource
            | Self::TorsionAddressOutsideDefinedSource { .. } => "Invariant Violation",
            Self::TransferPostcondition => "Post-condition Violation",
        }
    }

    const fn source_message(self) -> &'static str {
        match self {
            Self::TorsionParametersOutsideDefinedSource => {
                "Native torsion parameter outside defined source"
            }
            Self::TorsionAddressOutsideDefinedSource { .. } => {
                "Native torsion address outside defined source"
            }
            Self::OopParametersMissing => "no OOP parameters",
            Self::OopAddressOutsideDefinedSource { .. } => {
                "OOP address outside defined source storage"
            }
            Self::NoPoints => "no points",
            Self::NoDistanceMatrix => "no distance matrix",
            Self::MatrixSizeMismatch => "matrix size mismatch",
            Self::NotInitialized => "not initialized",
            Self::BadBounds => "bad bounds",
            Self::BadBondOrder => "bad bond order",
            Self::AngleDegeneratePoints => "degenerate points",
            Self::AngleBadOrder { .. } => "bad order",
            Self::TorsionDegeneratePoints => "degenerate points",
            Self::TorsionBadHybridizations => "bad hybridizations",
            Self::TorsionBadOrder { .. } => "bad order",
            Self::PositionCountMismatch { .. } => "size mismatch",
            Self::OptimizerBadTolerance => "bad tolerance",
            Self::OptimizerBadDirection => "bad direction in linearSearch",
            Self::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                ..
            } => "i",
            Self::IndexOutOfRange {
                argument: ForceFieldIndexArgument::J,
                ..
            } => "j",
            Self::BondIndexOutOfRange {
                argument: BondIndexArgument::First,
                ..
            } => "idx1",
            Self::BondIndexOutOfRange {
                argument: BondIndexArgument::Second,
                ..
            } => "idx2",
            Self::AngleIndexOutOfRange {
                argument: AngleIndexArgument::First,
                ..
            } => "idx1",
            Self::AngleIndexOutOfRange {
                argument: AngleIndexArgument::Second,
                ..
            } => "idx2",
            Self::AngleIndexOutOfRange {
                argument: AngleIndexArgument::Third,
                ..
            } => "idx3",
            Self::TorsionIndexOutOfRange {
                argument: TorsionIndexArgument::First,
                ..
            } => "idx1",
            Self::TorsionIndexOutOfRange {
                argument: TorsionIndexArgument::Second,
                ..
            } => "idx2",
            Self::TorsionIndexOutOfRange {
                argument: TorsionIndexArgument::Third,
                ..
            } => "idx3",
            Self::TorsionIndexOutOfRange {
                argument: TorsionIndexArgument::Fourth,
                ..
            } => "idx4",
            Self::AngleOutOfRange {
                bound: AngleRangeBound::Minimum,
            } => "minAngleDeg",
            Self::AngleOutOfRange {
                bound: AngleRangeBound::Maximum,
            } => "maxAngleDeg",
            Self::AngleBoundsOrder => "minAngleDeg must be <= maxAngleDeg",
            Self::PackedAngleBoundsOrder => "minAngleDeg must be <= maxAngleDeg",
            Self::TorsionBoundsOrder => "minDihedralDeg must be <= maxDihedralDeg",
            Self::BadIndex => "Bad index",
            Self::BadFixedPoint { .. } => "bad fixed point index",
            Self::TransferPostcondition => "bad index",
        }
    }

    fn source_expression(self) -> Option<String> {
        match self {
            Self::OopParametersMissing => Some("mmffOopParams".to_owned()),
            // Native safe-boundary error: C++ defines no corresponding expression.
            Self::OopAddressOutsideDefinedSource { .. }
            | Self::TorsionParametersOutsideDefinedSource
            | Self::TorsionAddressOutsideDefinedSource { .. } => None,
            Self::NoPoints => Some("d_numPoints".to_owned()),
            Self::NoDistanceMatrix => Some("dp_distMat".to_owned()),
            Self::MatrixSizeMismatch => Some(
                "static_cast<unsigned int>(d_numPoints * (d_numPoints + 1) / 2) <= d_matSize"
                    .to_owned(),
            ),
            Self::NotInitialized => Some("df_init".to_owned()),
            Self::BadBounds => Some("maxLen >= minLen".to_owned()),
            Self::PositionCountMismatch { .. } => {
                Some("static_cast<unsigned int>(d_numPoints) == d_positions.size()".to_owned())
            }
            Self::OptimizerBadTolerance => Some("gradTol > 0".to_owned()),
            Self::OptimizerBadDirection => Some("status >= 0".to_owned()),
            Self::IndexOutOfRange {
                index, upper_bound, ..
            } => Some(format!("{index} < {upper_bound}")),
            Self::BondIndexOutOfRange {
                index, upper_bound, ..
            } => Some(format!("{index} < {upper_bound}")),
            Self::BadBondOrder => Some("bondOrder > 0".to_owned()),
            Self::AngleDegeneratePoints => {
                Some("(idx1 != idx2 && idx2 != idx3 && idx1 != idx3)".to_owned())
            }
            Self::AngleBadOrder { .. } => Some(
                "d_order == 0 || d_order == 1 || d_order == 2 || d_order == 3 || d_order == 4"
                    .to_owned(),
            ),
            Self::TorsionDegeneratePoints => Some(
                "(idx1 != idx2 && idx1 != idx3 && idx1 != idx4 && idx2 != idx3 && idx2 != idx4 && idx3 != idx4)"
                    .to_owned(),
            ),
            Self::TorsionBadHybridizations => Some(
                "(hyb2 == RDKit::Atom::SP2 || hyb2 == RDKit::Atom::SP3) && (hyb3 == RDKit::Atom::SP2 || hyb3 == RDKit::Atom::SP3)"
                    .to_owned(),
            ),
            Self::TorsionBadOrder { .. } => {
                Some("d_order == 2 || d_order == 3 || d_order == 6".to_owned())
            }
            Self::AngleIndexOutOfRange {
                index, upper_bound, ..
            } => Some(format!("{index} < {upper_bound}")),
            Self::TorsionIndexOutOfRange {
                index, upper_bound, ..
            } => Some(format!("{index} < {upper_bound}")),
            Self::AngleOutOfRange {
                bound: AngleRangeBound::Minimum,
            } => Some("minAngleDeg".to_owned()),
            Self::AngleOutOfRange {
                bound: AngleRangeBound::Maximum,
            } => Some("maxAngleDeg".to_owned()),
            Self::AngleBoundsOrder => Some("!(minAngleDeg > maxAngleDeg)".to_owned()),
            Self::PackedAngleBoundsOrder => Some("maxAngleDeg >= minAngleDeg".to_owned()),
            Self::TorsionBoundsOrder => Some("!(minDihedralDeg > maxDihedralDeg)".to_owned()),
            Self::BadIndex => Some("idx < d_matSize".to_owned()),
            Self::BadFixedPoint { index, .. } => {
                Some(format!("static_cast<unsigned int>({index}) < d_numPoints"))
            }
            Self::TransferPostcondition => {
                Some("tab == this->dimension() * d_positions.size()".to_owned())
            }
        }
    }
}

impl std::fmt::Display for ForceFieldKernelError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for ForceFieldKernelError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

#[path = "constraints/angle.rs"]
mod angle;
#[path = "constraints/angles.rs"]
mod angles;
#[path = "constraints/distance.rs"]
mod distance;
#[path = "constraints/distances.rs"]
mod distances;
#[path = "constraints/position.rs"]
mod position;
#[path = "constraints/torsion.rs"]
mod torsion;

pub(super) struct EvaluationContext<'a> {
    coordinates: &'a [f64],
    distance_matrix: &'a mut [f64],
    initialized: bool,
    dimension: u32,
    num_points: u32,
    matrix_size: u32,
}

#[cfg(test)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Cf3dFragAcceptContributionIdentity {
    BondStretch {
        end1_idx: u32,
        end2_idx: u32,
        rest_len: f64,
        force_constant: f64,
    },
    AngleBend {
        at1_idx: u32,
        at2_idx: u32,
        at3_idx: u32,
        order: u32,
        force_constant: f64,
        c0: f64,
        c1: f64,
        c2: f64,
        theta0: f64,
    },
    Vdw {
        at1_idx: u32,
        at2_idx: u32,
        x_ij: f64,
        well_depth: f64,
        threshold: f64,
    },
    TorsionAngle {
        at1_idx: u32,
        at2_idx: u32,
        at3_idx: u32,
        at4_idx: u32,
        order: u32,
        force_constant: f64,
        cos_term: f64,
    },
    Inversion {
        at1_idx: u32,
        at2_idx: u32,
        at3_idx: u32,
        at4_idx: u32,
        force_constant: f64,
        c0: f64,
        c1: f64,
        c2: f64,
    },
    Untracked,
}

pub(super) trait ForceFieldContribution: Send {
    // BEGIN RDKIT CPP INTERFACE ForceFields::ForceFieldContrib (Contrib.h:18-35)
    // RDKit❗✔️: class RDKIT_FORCEFIELD_EXPORT ForceFieldContrib {
    // RDKit❗✔️:  public:
    // RDKit❗✔️:   friend class ForceField;
    // RDKit❗✔️:   ForceFieldContrib() {}
    // RDKit❗✔️:   ForceFieldContrib(ForceFields::ForceField *owner) : dp_forceField(owner) {}
    // RDKit❗✔️:   virtual ~ForceFieldContrib() {}
    // RDKit❗✔️:   //! returns our contribution to the energy of a position
    // RDKit❗✔️:   virtual double getEnergy(double *pos) const = 0;
    // RDKit❗✔️:   //! calculates our contribution to the gradients of a position
    // RDKit❗✔️:   virtual void getGrad(double *pos, double *grad) const = 0;
    // RDKit❗✔️:   //! return a copy
    // RDKit❗✔️:   virtual ForceFieldContrib *copy() const = 0;
    // RDKit❗✔️:  protected:
    // RDKit❗✔️:   ForceField *dp_forceField{nullptr};  //!< our owning ForceField
    // RDKit❗✔️: };
    // END RDKIT CPP INTERFACE ForceFields::ForceFieldContrib
    // Rust passes owner-derived reads and cache access as a short-lived borrow;
    // contributions do not retain an owner pointer.
    #[cfg(test)]
    fn cf3d_frag_accept_test_identity(&self) -> Cf3dFragAcceptContributionIdentity {
        Cf3dFragAcceptContributionIdentity::Untracked
    }

    fn get_energy(&self, context: &mut EvaluationContext<'_>)
    -> Result<f64, ForceFieldKernelError>;
    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError>;
    fn copy(&self) -> Box<dyn ForceFieldContribution>;
}

impl ForceFieldContribution for crate::uff::bond::BondStretchContrib {
    #[cfg(test)]
    fn cf3d_frag_accept_test_identity(&self) -> Cf3dFragAcceptContributionIdentity {
        crate::uff::bond::BondStretchContrib::cf3d_frag_accept_stored_identity(self)
    }

    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        crate::uff::bond::BondStretchContrib::get_energy(self, context)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        crate::uff::bond::BondStretchContrib::get_grad(self, context, gradient)
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        Box::new(*self)
    }
}

pub(super) struct ForceField<'a> {
    dimension: u32,
    initialized: bool,
    num_points: u32,
    distance_matrix: Vec<f64>,
    matrix_allocated: bool,
    positions: Vec<&'a mut [f64]>,
    contributions: Vec<Box<dyn ForceFieldContribution>>,
    fixed_points: Vec<i32>,
    matrix_size: u32,
}

fn checked_force_field_pair(
    initialized: bool,
    num_points: u32,
    mut i: u32,
    mut j: u32,
) -> Result<(u32, u32), ForceFieldKernelError> {
    if !initialized {
        return Err(ForceFieldKernelError::NotInitialized);
    }
    if i >= num_points {
        return Err(ForceFieldKernelError::IndexOutOfRange {
            argument: ForceFieldIndexArgument::I,
            index: i,
            upper_bound: num_points,
        });
    }
    if j >= num_points {
        return Err(ForceFieldKernelError::IndexOutOfRange {
            argument: ForceFieldIndexArgument::J,
            index: j,
            upper_bound: num_points,
        });
    }
    if j < i {
        std::mem::swap(&mut i, &mut j);
    }
    Ok((i, j))
}

fn squared_coordinate_distance(
    dimension: u32,
    positions: &[&mut [f64]],
    i: u32,
    j: u32,
    pos: Option<&[f64]>,
) -> f64 {
    let mut result = 0.0;
    if let Some(pos) = pos {
        let mut pi = dimension.wrapping_mul(i) as usize;
        let mut pj = dimension.wrapping_mul(j) as usize;
        for _ in 0..dimension {
            let tmp = pos[pi] - pos[pj];
            result += tmp * tmp;
            pi += 1;
            pj += 1;
        }
    } else {
        for component in 0..dimension {
            let tmp = positions[i as usize][component as usize]
                - positions[j as usize][component as usize];
            result += tmp * tmp;
        }
    }
    result
}

fn source_distance2(
    initialized: bool,
    num_points: u32,
    dimension: u32,
    positions: &[&mut [f64]],
    i: u32,
    j: u32,
    pos: Option<&[f64]>,
) -> Result<f64, ForceFieldKernelError> {
    let (i, j) = checked_force_field_pair(initialized, num_points, i, j)?;
    Ok(squared_coordinate_distance(dimension, positions, i, j, pos))
}

fn source_distance(
    initialized: bool,
    num_points: u32,
    dimension: u32,
    matrix_size: u32,
    positions: &[&mut [f64]],
    distance_matrix: &mut [f64],
    i: u32,
    j: u32,
    pos: Option<&[f64]>,
) -> Result<f64, ForceFieldKernelError> {
    let (i, j) = checked_force_field_pair(initialized, num_points, i, j)?;
    let cache_index = i.wrapping_add(j.wrapping_mul(j.wrapping_add(1)) / 2);
    if cache_index >= matrix_size {
        return Err(ForceFieldKernelError::BadIndex);
    }

    let result = &mut distance_matrix[cache_index as usize];
    if *result < 0.0 {
        *result = squared_coordinate_distance(dimension, positions, i, j, pos).sqrt();
    }
    Ok(*result)
}

impl<'a> EvaluationContext<'a> {
    #[cfg(test)]
    pub(super) fn for_oop_cache_preservation_test(
        coordinates: &'a [f64],
        distance_matrix: &'a mut [f64],
        num_points: u32,
    ) -> Self {
        // Test fixture: preserve caller sentinel bytes before OOP evaluation.
        // Existing for_test resets caches for distance-using contribution tests.
        Self {
            coordinates,
            distance_matrix,
            initialized: true,
            dimension: 3,
            num_points,
            matrix_size: num_points * (num_points + 1) / 2,
        }
    }

    #[cfg(test)]
    pub(super) fn for_test(
        coordinates: &'a [f64],
        distance_matrix: &'a mut [f64],
        num_points: u32,
    ) -> Self {
        distance_matrix.fill(-1.0);
        Self {
            coordinates,
            distance_matrix,
            initialized: true,
            dimension: 3,
            num_points,
            matrix_size: num_points * (num_points + 1) / 2,
        }
    }

    pub(super) fn coordinates(&self) -> &[f64] {
        self.coordinates
    }

    pub(super) fn distance(&mut self, i: u32, j: u32) -> Result<f64, ForceFieldKernelError> {
        let (initialized, num_points, dimension, matrix_size, coordinates) = (
            self.initialized,
            self.num_points,
            self.dimension,
            self.matrix_size,
            self.coordinates,
        );
        source_distance(
            initialized,
            num_points,
            dimension,
            matrix_size,
            &[],
            &mut *self.distance_matrix,
            i,
            j,
            Some(coordinates),
        )
    }

    fn distance2(&self, i: u32, j: u32) -> Result<f64, ForceFieldKernelError> {
        source_distance2(
            self.initialized,
            self.num_points,
            self.dimension,
            &[],
            i,
            j,
            Some(self.coordinates),
        )
    }

    fn distance_const(&self, i: u32, j: u32) -> Result<f64, ForceFieldKernelError> {
        self.distance2(i, j).map(f64::sqrt)
    }
}

impl<'a> ForceField<'a> {
    pub(super) fn new(dimension: u32) -> Self {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::ForceField (ForceField.h:82)
        // RDKit✔️✔️: ForceField(unsigned int dimension = 3) : d_dimension(dimension) {}
        // END RDKIT CPP FUNCTION ForceFields::ForceField::ForceField
        #[cfg(test)]
        record_uff_one_field();
        Self {
            dimension,
            initialized: false,
            num_points: 0,
            distance_matrix: Vec::new(),
            matrix_allocated: false,
            positions: Vec::new(),
            contributions: Vec::new(),
            fixed_points: Vec::new(),
            matrix_size: 0,
        }
    }

    pub(super) fn copy<'copy>(&self) -> ForceField<'copy> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::ForceField(copy) (ForceField.cpp:159-170)
        // RDKit✔️✔️: ForceField::ForceField(const ForceField &other)
        // RDKit✔️✔️:     : d_dimension(other.d_dimension),
        // RDKit✔️✔️:       df_init(false),
        // RDKit✔️✔️:       d_numPoints(other.d_numPoints),
        // RDKit✔️✔️:       dp_distMat(nullptr) {
        // RDKit✔️✔️:   d_contribs.clear();
        // RDKit✔️✔️:   for (const auto &contrib : other.d_contribs) {
        // RDKit✔️✔️:     ForceFieldContrib *ncontrib = contrib->copy();
        // RDKit✔️✔️:     ncontrib->dp_forceField = this;
        // RDKit✔️✔️:     d_contribs.push_back(ContribPtr(ncontrib));
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: };
        // END RDKIT CPP FUNCTION ForceFields::ForceField::ForceField(copy)
        // The source leaves positions and fixed points at their default-empty
        // state. This independent empty borrow vector can later accept only
        // caller-supplied positions for the new kernel's chosen lifetime.
        let mut copied = ForceField::<'copy>::new(self.dimension);
        copied.initialized = false;
        copied.num_points = self.num_points;
        copied.contributions.reserve(self.contributions.len());
        for contribution in &self.contributions {
            // Source clones each term and rebinds its owner pointer. The Rust
            // term receives this copy's context when evaluated and stores no
            // pointer, so only the independent term copy is needed here.
            #[cfg(test)]
            record_uff_one_term_copy();
            copied.contributions.push(contribution.copy());
        }
        copied
    }

    fn dimension(&self) -> u32 {
        // RDKit✔️✔️: unsigned int dimension() const { return d_dimension; }
        self.dimension
    }

    fn num_points(&self) -> u32 {
        // RDKit✔️✔️: unsigned int numPoints() const { return d_numPoints; }
        self.num_points
    }

    pub(super) fn positions(&self) -> &[&'a mut [f64]] {
        // RDKit❗✔️: const RDGeom::PointPtrVect &positions() const { return d_positions; }
        &self.positions
    }

    pub(super) fn positions_mut(&mut self) -> &mut Vec<&'a mut [f64]> {
        // RDKit✔️✔️: RDGeom::PointPtrVect &positions() { return d_positions; }
        &mut self.positions
    }

    pub(super) fn rebind_positions<'next>(
        self,
        positions: Vec<&'next mut [f64]>,
    ) -> ForceField<'next> {
        // RDKit✔️✔️: RDGeom::PointPtrVect &positions() { return d_positions; }
        // RDKit✔️✔️: ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
        // ForceField.h:181; FFConvenience.h:72. Consuming the old borrowed
        // view releases its lifetime without copying any owned field state.
        // Source initialization remains a separate subsequent operation.
        let ForceField {
            dimension,
            initialized,
            num_points,
            distance_matrix,
            matrix_allocated,
            positions: old_positions,
            contributions,
            fixed_points,
            matrix_size,
        } = self;
        drop(old_positions);
        ForceField {
            dimension,
            initialized,
            num_points,
            distance_matrix,
            matrix_allocated,
            positions,
            contributions,
            fixed_points,
            matrix_size,
        }
    }

    pub(super) fn clear_positions_and_shorten_lifetime<'next>(mut self) -> ForceField<'next>
    where
        'a: 'next,
    {
        // RDKit✔️✔️: ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
        // FFConvenience.h:72 replaces every borrowed row in the existing
        // PointPtrVect. Drop the old row borrows, then use covariance to
        // shorten their now-empty Vec's lifetime while retaining its buffer.
        self.positions.clear();
        let ForceField {
            dimension,
            initialized,
            num_points,
            distance_matrix,
            matrix_allocated,
            positions,
            contributions,
            fixed_points,
            matrix_size,
        } = self;
        let positions: Vec<&'next mut [f64]> = positions;
        ForceField {
            dimension,
            initialized,
            num_points,
            distance_matrix,
            matrix_allocated,
            positions,
            contributions,
            fixed_points,
            matrix_size,
        }
    }

    pub(super) fn release_position_borrows<'next>(self) -> ForceField<'next> {
        // RDKit❗❗: ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
        // FFConvenience.h:72. This is the later source rebinding anchor; this
        // lifetime helper only releases the prior borrowed view.
        // Behavior review: consume every old mutable row reference before a
        // caller borrows those coordinates again. All non-position field
        // state moves unchanged; no chemistry stage runs here.
        // Allocation review: collect zero emitted references from the moved
        // Vec. L01 verifies address/capacity retention on this pinned build;
        // this test is not a portable std-library guarantee.
        let ForceField {
            dimension,
            initialized,
            num_points,
            distance_matrix,
            matrix_allocated,
            positions,
            contributions,
            fixed_points,
            matrix_size,
        } = self;
        let positions: Vec<&'next mut [f64]> = positions
            .into_iter()
            .filter_map(|_| None::<&'next mut [f64]>)
            .collect();
        ForceField {
            dimension,
            initialized,
            num_points,
            distance_matrix,
            matrix_allocated,
            positions,
            contributions,
            fixed_points,
            matrix_size,
        }
    }

    pub(super) fn add_contribution(&mut self, contribution: Box<dyn ForceFieldContribution>) {
        // BEGIN RDKIT CPP CALL UFF::Tools::addBonds (Builder.cpp:49)
        // RDKit❗✔️: field->contribs().push_back(ForceFields::ContribPtr(contrib));
        // END RDKIT CPP CALL UFF::Tools::addBonds
        #[cfg(test)]
        record_uff_one_term();
        self.contributions.push(contribution);
    }

    fn scatter_components<F>(&self, mut write_component: F) -> Result<(), ForceFieldKernelError>
    where
        F: FnMut(usize, f64),
    {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::scatter (ForceField.cpp:377-389)
        // RDKit❗✔️: void ForceField::scatter(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit❗✔️:   PRECONDITION(pos, "bad position vector");
        // RDKit❗✔️:   unsigned int tab = 0;
        // RDKit❗✔️:   for (auto d_position : d_positions) {
        // RDKit❗✔️:     for (unsigned int di = 0; di < this->dimension(); ++di) {
        // RDKit❗✔️:       pos[tab + di] = (*d_position)[di];  //->x;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     tab += this->dimension();
        // RDKit❗✔️:   }
        // RDKit❗✔️:   POSTCONDITION(tab == this->dimension() * d_positions.size(), "bad index");
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::scatter
        if !self.initialized {
            return Err(ForceFieldKernelError::NotInitialized);
        }
        // A Rust borrow is non-null; the source null-pointer condition is
        // represented structurally for this shared writer.
        let mut tab = 0_u32;
        for position in &self.positions {
            for di in 0..self.dimension {
                write_component(tab.wrapping_add(di) as usize, position[di as usize]);
            }
            tab = tab.wrapping_add(self.dimension);
        }
        if (tab as usize) != (self.dimension as usize).wrapping_mul(self.positions.len()) {
            return Err(ForceFieldKernelError::TransferPostcondition);
        }
        Ok(())
    }

    fn scatter(&self, pos: &mut [f64]) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::scatter (ForceField.cpp:377-389)
        // RDKit✔️✔️: void ForceField::scatter(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad position vector");
        // RDKit✔️✔️:   unsigned int tab = 0;
        // RDKit✔️✔️:   for (auto d_position : d_positions) {
        // RDKit✔️✔️:     for (unsigned int di = 0; di < this->dimension(); ++di) {
        // RDKit✔️✔️:       pos[tab + di] = (*d_position)[di];  //->x;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     tab += this->dimension();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   POSTCONDITION(tab == this->dimension() * d_positions.size(), "bad index");
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::scatter
        self.scatter_components(|index, value| pos[index] = value)
    }

    fn scatter_current_positions(&self) -> Result<Vec<f64>, ForceFieldKernelError> {
        let capacity = (self.dimension as usize).wrapping_mul(self.positions.len());
        let mut coordinates = Vec::with_capacity(capacity);
        self.scatter_components(|_, value| coordinates.push(value))?;
        Ok(coordinates)
    }

    fn gather(&mut self, pos: &[f64]) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::gather (ForceField.cpp:391-404)
        // RDKit✔️✔️: void ForceField::gather(double *pos) {
        // RDKit✔️✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad position vector");
        // RDKit✔️✔️:   unsigned int tab = 0;
        // RDKit✔️✔️:   for (auto &d_position : d_positions) {
        // RDKit✔️✔️:     for (unsigned int di = 0; di < this->dimension(); ++di) {
        // RDKit✔️✔️:       (*d_position)[di] = pos[tab + di];
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     tab += this->dimension();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::gather
        if !self.initialized {
            return Err(ForceFieldKernelError::NotInitialized);
        }

        // `&[f64]` cannot be null; preserve the source's no-copy input path.
        let mut tab = 0_u32;
        for position in &mut self.positions {
            for di in 0..self.dimension {
                position[di as usize] = pos[tab.wrapping_add(di) as usize];
            }
            tab = tab.wrapping_add(self.dimension);
        }
        Ok(())
    }

    pub(super) fn minimize(
        &mut self,
        max_its: u32,
        force_tol: f64,
        energy_tol: f64,
    ) -> Result<i32, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::minimize overload (ForceField.cpp:252-254)
        // RDKit✔️✔️: int ForceField::minimize(unsigned int maxIts, double forceTol,
        // RDKit✔️✔️:                          double energyTol) {
        // RDKit✔️✔️:   return minimize(0, nullptr, maxIts, forceTol, energyTol);
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::minimize overload
        self.minimize_with_snapshots(0, None, max_its, force_tol, energy_tol)
    }

    fn minimize_with_snapshots(
        &mut self,
        snapshot_freq: u32,
        snapshots: Option<&mut Vec<crate::optimizer::OptimizerSnapshot>>,
        max_its: u32,
        force_tol: f64,
        energy_tol: f64,
    ) -> Result<i32, ForceFieldKernelError> {
        use crate::optimizer::{OptimizerError, minimize};

        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::minimize (ForceField.cpp:257-282)
        // RDKit✔️✔️: int ForceField::minimize(unsigned int snapshotFreq,
        // RDKit✔️✔️:                          RDKit::SnapshotVect *snapshotVect, unsigned int maxIts,
        // RDKit✔️✔️:                          double forceTol, double energyTol) {
        // RDKit✔️✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️:   PRECONDITION(static_cast<unsigned int>(d_numPoints) == d_positions.size(),
        // RDKit✔️✔️:                "size mismatch");
        // RDKit✔️✔️:   if (d_contribs.empty()) {
        // RDKit✔️✔️:     return 0;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   unsigned int numIters = 0;
        // RDKit✔️✔️:   unsigned int dim = this->d_numPoints * d_dimension;
        // RDKit✔️✔️:   double finalForce = 0.0;
        // RDKit✔️✔️:   std::vector<double> points(dim);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   this->scatter(points.data());
        // RDKit✔️✔️:   ForceFieldsHelper::calcEnergy eCalc(this);
        // RDKit✔️✔️:   ForceFieldsHelper::calcGradient gCalc(this);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   int res =
        // RDKit✔️✔️:       BFGSOpt::minimize(dim, points.data(), forceTol, numIters, finalForce, eCalc,
        // RDKit✔️✔️:                         gCalc, snapshotFreq, snapshotVect, energyTol, maxIts);
        // RDKit✔️✔️:   this->gather(points.data());
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::minimize
        if !self.initialized {
            return Err(ForceFieldKernelError::NotInitialized);
        }
        if self.num_points as usize != self.positions.len() {
            return Err(ForceFieldKernelError::PositionCountMismatch {
                num_points: self.num_points,
                position_count: self.positions.len(),
            });
        }
        if self.contributions.is_empty() {
            return Ok(0);
        }

        let dimension = self.num_points.wrapping_mul(self.dimension) as usize;
        let mut points = vec![0.0; dimension];
        self.scatter(&mut points)?;
        let mut num_iters = 0;
        let mut final_force = 0.0;

        // BEGIN RDKIT CPP FUNCTION ForceFieldsHelper::calcEnergy::operator() (ForceField.cpp:97-104)
        // RDKit✔️✔️: class calcEnergy {
        // RDKit✔️✔️:  public:
        // RDKit✔️✔️:   calcEnergy(ForceFields::ForceField *ffHolder) : mp_ffHolder(ffHolder) {};
        // RDKit✔️✔️:   double operator()(double *pos) const { return mp_ffHolder->calcEnergy(pos); }
        // RDKit✔️✔️:
        // RDKit✔️✔️:  private:
        // RDKit✔️✔️:   ForceFields::ForceField *mp_ffHolder;
        // RDKit✔️✔️: };
        // END RDKIT CPP FUNCTION ForceFieldsHelper::calcEnergy::operator()
        let mut energy = |field: &mut Self, coordinates: &mut [f64]| field.calc_energy(coordinates);

        // BEGIN RDKIT CPP FUNCTION ForceFieldsHelper::calcGradient::operator() (ForceField.cpp:106-147)
        // RDKit✔️✔️: class calcGradient {
        // RDKit✔️✔️:  public:
        // RDKit✔️✔️:   calcGradient(ForceFields::ForceField *ffHolder) : mp_ffHolder(ffHolder) {};
        // RDKit✔️✔️:   double operator()(double *pos, double *grad) const {
        // RDKit✔️✔️:     double res = 1.0;
        // RDKit✔️✔️:     // the contribs to the gradient function use +=, so we need
        // RDKit✔️✔️:     // to zero the grad out before moving on:
        // RDKit✔️✔️:     for (unsigned int i = 0;
        // RDKit✔️✔️:          i < mp_ffHolder->numPoints() * mp_ffHolder->dimension(); i++) {
        // RDKit✔️✔️:       grad[i] = 0.0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     mp_ffHolder->calcGrad(pos, grad);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     // FIX: this hack reduces the gradients so that the
        // RDKit✔️✔️:     // minimizer is more efficient.
        // RDKit✔️✔️:     double maxGrad = -1e8;
        // RDKit✔️✔️:     double gradScale = 0.1;
        // RDKit✔️✔️:     for (unsigned int i = 0;
        // RDKit✔️✔️:          i < mp_ffHolder->numPoints() * mp_ffHolder->dimension(); i++) {
        // RDKit✔️✔️:       grad[i] *= gradScale;
        // RDKit✔️✔️:       if (fabs(grad[i]) > maxGrad) {
        // RDKit✔️✔️:         maxGrad = fabs(grad[i]);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     // this is a continuation of the same hack to avoid
        // RDKit✔️✔️:     // some potential numeric instabilities:
        // RDKit✔️✔️:     if (maxGrad > 10.0) {
        // RDKit✔️✔️:       while (maxGrad * gradScale > 10.0) {
        // RDKit✔️✔️:         gradScale *= .5;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       for (unsigned int i = 0;
        // RDKit✔️✔️:            i < mp_ffHolder->numPoints() * mp_ffHolder->dimension(); i++) {
        // RDKit✔️✔️:         grad[i] *= gradScale;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     res = gradScale;
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:  private:
        // RDKit✔️✔️:   ForceFields::ForceField *mp_ffHolder;
        // RDKit✔️✔️: };
        // END RDKIT CPP FUNCTION ForceFieldsHelper::calcGradient::operator()
        let mut gradient = |field: &mut Self, coordinates: &mut [f64], gradient: &mut [f64]| {
            gradient.fill(0.0);
            field.calc_grad(coordinates, gradient)?;

            let dimension = field.num_points.wrapping_mul(field.dimension) as usize;
            let mut max_grad = -1e8;
            let mut grad_scale = 0.1;
            for value in &mut gradient[..dimension] {
                *value *= grad_scale;
                if value.abs() > max_grad {
                    max_grad = value.abs();
                }
            }
            if max_grad > 10.0 {
                while max_grad * grad_scale > 10.0 {
                    grad_scale *= 0.5;
                }
                for value in &mut gradient[..dimension] {
                    *value *= grad_scale;
                }
            }
            Ok(grad_scale)
        };

        // The flat points vector is the source minimizer's sole owner-side
        // buffer. Both fallible callbacks borrow this same ForceField; there
        // is no owner clone, boxing, replay, or retained pointer.
        let optimizer_result = minimize(
            self,
            &mut points,
            force_tol,
            &mut num_iters,
            &mut final_force,
            &mut energy,
            &mut gradient,
            snapshot_freq,
            snapshots,
            energy_tol,
            max_its,
        );
        let status = match optimizer_result {
            Ok(status) => status,
            Err(OptimizerError::Evaluation(error)) => return Err(error),
            Err(OptimizerError::BadTolerance) => {
                return Err(ForceFieldKernelError::OptimizerBadTolerance);
            }
            Err(OptimizerError::BadDirection) => {
                return Err(ForceFieldKernelError::OptimizerBadDirection);
            }
        };

        // RDKit gathers after either normal optimizer status (0 or 1), but
        // never after a callback or invariant error.
        self.gather(&points)?;
        Ok(status)
    }

    pub(super) fn calc_energy_current(
        &mut self,
        mut contribution_energies: Option<&mut Vec<f64>>,
    ) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::calcEnergy(vector*) (ForceField.cpp:284-307)
        // RDKit✔️✔️: double ForceField::calcEnergy(std::vector<double> *contribs) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   if (d_contribs.empty()) {
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (contribs) {
        // RDKit✔️✔️:     contribs->clear();
        // RDKit✔️✔️:     contribs->reserve(d_contribs.size());
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   unsigned int N = d_positions.size();
        // RDKit✔️✔️:   auto *pos = new double[d_dimension * N];
        // RDKit✔️✔️:   this->scatter(pos);
        // RDKit✔️✔️:   // now loop over the contribs
        // RDKit✔️✔️:   for (const auto &d_contrib : d_contribs) {
        // RDKit✔️✔️:     double e = d_contrib->getEnergy(pos);
        // RDKit✔️✔️:     res += e;
        // RDKit✔️✔️:     if (contribs) {
        // RDKit✔️✔️:       contribs->push_back(e);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   delete[] pos;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::calcEnergy(vector*)
        if !self.initialized {
            return Err(ForceFieldKernelError::NotInitialized);
        }
        let mut energy = 0.0;
        if self.contributions.is_empty() {
            return Ok(energy);
        }
        if let Some(energies) = contribution_energies.as_deref_mut() {
            energies.clear();
            energies.reserve(self.contributions.len());
        }

        // A capacity-only Vec gives RDKit's single flat allocation without
        // zero-filling values that scatter immediately overwrites.
        let coordinates = self.scatter_current_positions()?;
        let mut context = EvaluationContext {
            coordinates: &coordinates,
            distance_matrix: &mut self.distance_matrix,
            initialized: self.initialized,
            dimension: self.dimension,
            num_points: self.num_points,
            matrix_size: self.matrix_size,
        };
        for contribution in &self.contributions {
            let term_energy = contribution.get_energy(&mut context)?;
            energy += term_energy;
            if let Some(energies) = contribution_energies.as_deref_mut() {
                energies.push(term_energy);
            }
        }
        Ok(energy)
    }

    fn calc_energy(&mut self, coordinates: &[f64]) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::calcEnergy(pos) (ForceField.cpp:310-325)
        // RDKit✔️✔️: double ForceField::calcEnergy(double *pos) {
        // RDKit✔️✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad position vector");
        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   this->initDistanceMatrix();
        // RDKit✔️✔️:   if (d_contribs.empty()) {
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   // now loop over the contribs
        // RDKit✔️✔️:   for (ContribPtrVect::const_iterator contrib = d_contribs.begin();
        // RDKit✔️✔️:        contrib != d_contribs.end(); contrib++) {
        // RDKit✔️✔️:     double E = (*contrib)->getEnergy(pos);
        // RDKit✔️✔️:     res += E;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::calcEnergy(pos)
        if !self.initialized {
            return Err(ForceFieldKernelError::NotInitialized);
        }
        // `&[f64]` cannot be null; an empty slice is not a null-pointer error.
        let mut energy = 0.0;
        self.init_distance_matrix()?;
        if self.contributions.is_empty() {
            return Ok(energy);
        }

        let mut context = EvaluationContext {
            coordinates,
            distance_matrix: &mut self.distance_matrix,
            initialized: self.initialized,
            dimension: self.dimension,
            num_points: self.num_points,
            matrix_size: self.matrix_size,
        };
        for contribution in &self.contributions {
            energy += contribution.get_energy(&mut context)?;
        }
        Ok(energy)
    }

    pub(super) fn calc_grad_current(
        &mut self,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::calcGrad(current) (ForceField.cpp:329-352)
        // RDKit✔️✔️: void ForceField::calcGrad(double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad gradient vector");
        // RDKit✔️✔️:   if (d_contribs.empty()) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   unsigned int N = d_positions.size();
        // RDKit✔️✔️:   auto *pos = new double[d_dimension * N];
        // RDKit✔️✔️:   this->scatter(pos);
        // RDKit✔️✔️:   for (const auto &d_contrib : d_contribs) {
        // RDKit✔️✔️:     d_contrib->getGrad(pos, grad);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   // zero out gradient values for any fixed points:
        // RDKit✔️✔️:   for (int d_fixedPoint : d_fixedPoints) {
        // RDKit✔️✔️:     CHECK_INVARIANT(static_cast<unsigned int>(d_fixedPoint) < d_numPoints,
        // RDKit✔️✔️:                     "bad fixed point index");
        // RDKit✔️✔️:     unsigned int idx = d_dimension * d_fixedPoint;
        // RDKit✔️✔️:     for (unsigned int di = 0; di < this->dimension(); ++di) {
        // RDKit✔️✔️:       grad[idx + di] = 0.0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   delete[] pos;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::calcGrad(current)
        if !self.initialized {
            return Err(ForceFieldKernelError::NotInitialized);
        }
        // A mutable Rust slice is non-null; preserve the source precondition
        // without rejecting an empty slice as a null pointer.
        if self.contributions.is_empty() {
            return Ok(());
        }

        let coordinates = self.scatter_current_positions()?;
        {
            let mut context = EvaluationContext {
                coordinates: &coordinates,
                distance_matrix: &mut self.distance_matrix,
                initialized: self.initialized,
                dimension: self.dimension,
                num_points: self.num_points,
                matrix_size: self.matrix_size,
            };
            for contribution in &self.contributions {
                contribution.get_grad(&mut context, gradient)?;
            }
        }
        self.zero_fixed_point_gradients(gradient)
    }

    fn calc_grad(
        &mut self,
        coordinates: &[f64],
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::calcGrad(pos) (ForceField.cpp:353-375)
        // RDKit✔️✔️: void ForceField::calcGrad(double *pos, double *grad) {
        // RDKit✔️✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad position vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad gradient vector");
        // RDKit✔️✔️:   if (d_contribs.empty()) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   for (ContribPtrVect::const_iterator contrib = d_contribs.begin();
        // RDKit✔️✔️:        contrib != d_contribs.end(); contrib++) {
        // RDKit✔️✔️:     (*contrib)->getGrad(pos, grad);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   for (INT_VECT::const_iterator it = d_fixedPoints.begin();
        // RDKit✔️✔️:        it != d_fixedPoints.end(); it++) {
        // RDKit✔️✔️:     CHECK_INVARIANT(static_cast<unsigned int>(*it) < d_numPoints,
        // RDKit✔️✔️:                     "bad fixed point index");
        // RDKit✔️✔️:     unsigned int idx = d_dimension * (*it);
        // RDKit✔️✔️:     for (unsigned int di = 0; di < this->dimension(); ++di) {
        // RDKit✔️✔️:       grad[idx + di] = 0.0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::calcGrad(pos)
        if !self.initialized {
            return Err(ForceFieldKernelError::NotInitialized);
        }
        // Borrowed slices encode both non-null preconditions structurally.
        if self.contributions.is_empty() {
            return Ok(());
        }

        {
            let mut context = EvaluationContext {
                coordinates,
                distance_matrix: &mut self.distance_matrix,
                initialized: self.initialized,
                dimension: self.dimension,
                num_points: self.num_points,
                matrix_size: self.matrix_size,
            };
            for contribution in &self.contributions {
                contribution.get_grad(&mut context, gradient)?;
            }
        }
        self.zero_fixed_point_gradients(gradient)
    }

    fn zero_fixed_point_gradients(
        &self,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FIXED-POINT BRANCHES ForceFields::ForceField::calcGrad
        // RDKit✔️✔️: for (int d_fixedPoint : d_fixedPoints) {
        // RDKit✔️✔️:   CHECK_INVARIANT(static_cast<unsigned int>(d_fixedPoint) < d_numPoints,
        // RDKit✔️✔️:                   "bad fixed point index");
        // RDKit✔️✔️:   unsigned int idx = d_dimension * d_fixedPoint;
        // RDKit✔️✔️:   for (unsigned int di = 0; di < this->dimension(); ++di) {
        // RDKit✔️✔️:     grad[idx + di] = 0.0;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FIXED-POINT BRANCHES
        for &fixed_point in &self.fixed_points {
            // Match the source's signed-to-unsigned comparison. Negative
            // fixed point IDs therefore fail the same invariant.
            let source_index = fixed_point as u32;
            if source_index >= self.num_points {
                return Err(ForceFieldKernelError::BadFixedPoint {
                    index: fixed_point,
                    upper_bound: self.num_points,
                });
            }
            let offset = self.dimension.wrapping_mul(source_index);
            for component in 0..self.dimension {
                gradient[offset.wrapping_add(component) as usize] = 0.0;
            }
        }
        Ok(())
    }

    fn fixed_points(&self) -> &[i32] {
        // RDKit✔️✔️: const INT_VECT &fixedPoints() const { return d_fixedPoints; }
        &self.fixed_points
    }

    fn fixed_points_mut(&mut self) -> &mut Vec<i32> {
        // RDKit✔️✔️: INT_VECT &fixedPoints() { return d_fixedPoints; }
        &mut self.fixed_points
    }

    fn distance(
        &mut self,
        i: u32,
        j: u32,
        pos: Option<&[f64]>,
    ) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::distance (ForceField.cpp:172-203)
        // RDKit✔️✔️: double ForceField::distance(unsigned int i, unsigned int j, double *pos) {
        // RDKit✔️✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️:   URANGE_CHECK(i, d_numPoints);
        // RDKit✔️✔️:   URANGE_CHECK(j, d_numPoints);
        // RDKit✔️✔️:   if (j < i) {
        // RDKit✔️✔️:     int tmp = j;
        // RDKit✔️✔️:     j = i;
        // RDKit✔️✔️:     i = tmp;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   unsigned int idx = i + j * (j + 1) / 2;
        // RDKit✔️✔️:   CHECK_INVARIANT(idx < d_matSize, "Bad index");
        // RDKit✔️✔️:   double &res = dp_distMat[idx];
        // RDKit✔️✔️:   if (res < 0.0) {
        // RDKit✔️✔️:     // we need to calculate this distance:
        // RDKit✔️✔️:     if (!pos) {
        // RDKit✔️✔️:       res = 0.0;
        // RDKit✔️✔️:       for (unsigned int idx = 0; idx < d_dimension; ++idx) {
        // RDKit✔️✔️:         double tmp =
        // RDKit✔️✔️:             (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
        // RDKit✔️✔️:         res += tmp * tmp;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       res = 0.0;
        // RDKit✔️✔️:       double *pi = &(pos[d_dimension * i]), *pj = &(pos[d_dimension * j]);
        // RDKit✔️✔️:       for (unsigned int idx = 0; idx < d_dimension; ++idx, ++pi, ++pj) {
        // RDKit✔️✔️:         double tmp = *pi - *pj;
        // RDKit✔️✔️:         res += tmp * tmp;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     res = sqrt(res);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::distance
        source_distance(
            self.initialized,
            self.num_points,
            self.dimension,
            self.matrix_size,
            &self.positions,
            &mut self.distance_matrix,
            i,
            j,
            pos,
        )
    }

    fn distance2(&self, i: u32, j: u32, pos: Option<&[f64]>) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::distance2 (ForceField.cpp:206-232)
        // RDKit✔️✔️: double ForceField::distance2(unsigned int i, unsigned int j,
        // RDKit✔️✔️:                              double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️:   URANGE_CHECK(i, d_numPoints);
        // RDKit✔️✔️:   URANGE_CHECK(j, d_numPoints);
        // RDKit✔️✔️:   if (j < i) {
        // RDKit✔️✔️:     int tmp = j;
        // RDKit✔️✔️:     j = i;
        // RDKit✔️✔️:     i = tmp;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   double res;
        // RDKit✔️✔️:   if (!pos) {
        // RDKit✔️✔️:     res = 0.0;
        // RDKit✔️✔️:     for (unsigned int idx = 0; idx < d_dimension; ++idx) {
        // RDKit✔️✔️:       double tmp =
        // RDKit✔️✔️:           (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
        // RDKit✔️✔️:       res += tmp * tmp;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     res = 0.0;
        // RDKit✔️✔️:     double *pi = &(pos[d_dimension * i]), *pj = &(pos[d_dimension * j]);
        // RDKit✔️✔️:     for (unsigned int idx = 0; idx < d_dimension; ++idx, ++pi, ++pj) {
        // RDKit✔️✔️:       double tmp = *pi - *pj;
        // RDKit✔️✔️:       res += tmp * tmp;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::distance2
        source_distance2(
            self.initialized,
            self.num_points,
            self.dimension,
            &self.positions,
            i,
            j,
            pos,
        )
    }

    fn distance_const(
        &self,
        i: u32,
        j: u32,
        pos: Option<&[f64]>,
    ) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::distance const overload (ForceField.cpp:234-237)
        // RDKit✔️✔️: double ForceField::distance(unsigned int i, unsigned int j, double *pos) const {
        // RDKit✔️✔️:   auto res = sqrt(distance2(i, j, pos));
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::distance const overload
        self.distance2(i, j, pos).map(f64::sqrt)
    }

    pub(super) fn initialize(&mut self) -> Result<(), ForceFieldKernelError> {
        // FFConvenience.h:74 calls this same initializer after replacing
        // position borrows. Preserve cleanup even for identical coordinates;
        // no row components are read until subsequent evaluation.
        #[cfg(test)]
        record_uff_one_initialize();
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::initialize (ForceField.cpp:239-250)
        // RDKit✔️✔️: void ForceField::initialize() {
        // RDKit✔️✔️:   df_init = false;
        // RDKit✔️✔️:   delete[] dp_distMat;
        // RDKit✔️✔️:   dp_distMat = nullptr;
        // RDKit✔️✔️:   d_numPoints = d_positions.size();
        // RDKit✔️✔️:   d_matSize = d_numPoints * (d_numPoints + 1) / 2;
        // RDKit✔️✔️:   dp_distMat = new double[d_matSize];
        // RDKit✔️✔️:   this->initDistanceMatrix();
        // RDKit✔️✔️:   df_init = true;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::initialize
        self.initialized = false;
        self.distance_matrix = Vec::new();
        self.matrix_allocated = false;
        self.num_points = self.positions.len() as u32;
        self.matrix_size = self
            .num_points
            .wrapping_mul(self.num_points.wrapping_add(1))
            / 2;
        self.distance_matrix
            .reserve_exact(self.matrix_size as usize);
        self.matrix_allocated = true;
        self.init_distance_matrix()?;
        self.initialized = true;
        Ok(())
    }

    fn init_distance_matrix(&mut self) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::ForceField::initDistanceMatrix (ForceField.cpp:404-415)
        // RDKit✔️✔️: void ForceField::initDistanceMatrix() {
        // RDKit✔️✔️:   PRECONDITION(d_numPoints, "no points");
        // RDKit✔️✔️:   PRECONDITION(dp_distMat, "no distance matrix");
        // RDKit✔️✔️:   PRECONDITION(static_cast<unsigned int>(d_numPoints * (d_numPoints + 1) / 2) <=
        // RDKit✔️✔️:                    d_matSize,
        // RDKit✔️✔️:                 "matrix size mismatch");
        // RDKit✔️✔️:   for (unsigned int i = 0; i < d_numPoints * (d_numPoints + 1) / 2; i++) {
        // RDKit✔️✔️:     dp_distMat[i] = -1.0;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::ForceField::initDistanceMatrix
        if self.num_points == 0 {
            return Err(ForceFieldKernelError::NoPoints);
        }
        if !self.matrix_allocated {
            return Err(ForceFieldKernelError::NoDistanceMatrix);
        }
        let expected_size = self
            .num_points
            .wrapping_mul(self.num_points.wrapping_add(1))
            / 2;
        if expected_size > self.matrix_size {
            return Err(ForceFieldKernelError::MatrixSizeMismatch);
        }

        let expected_len = expected_size as usize;
        if self.distance_matrix.len() < expected_len {
            self.distance_matrix.resize(expected_len, -1.0);
        } else {
            self.distance_matrix[..expected_len].fill(-1.0);
        }
        Ok(())
    }

    fn evaluation_context<'ctx>(
        &'ctx mut self,
        coordinates: &'ctx [f64],
    ) -> EvaluationContext<'ctx> {
        EvaluationContext {
            coordinates,
            distance_matrix: &mut self.distance_matrix,
            initialized: self.initialized,
            dimension: self.dimension,
            num_points: self.num_points,
            matrix_size: self.matrix_size,
        }
    }
}

impl Default for ForceField<'_> {
    fn default() -> Self {
        Self::new(3)
    }
}

#[cfg(test)]
pub(super) fn cf3d_bld_b05_copy_force_field<'copy>(
    force_field: &ForceField<'_>,
) -> ForceField<'copy> {
    force_field.copy()
}

#[cfg(test)]
pub(super) fn cf3d_bld_b05_calc_energy(
    force_field: &mut ForceField<'_>,
    coordinates: &[f64],
) -> Result<f64, ForceFieldKernelError> {
    force_field.calc_energy(coordinates)
}

#[cfg(test)]
pub(super) fn cf3d_bld_b05_calc_grad(
    force_field: &mut ForceField<'_>,
    coordinates: &[f64],
    gradient: &mut [f64],
) -> Result<(), ForceFieldKernelError> {
    force_field.calc_grad(coordinates, gradient)
}

#[cfg(test)]
pub(super) fn cf3d_frag_f24_contribution_energies(
    force_field: &mut ForceField<'_>,
) -> Result<Vec<f64>, ForceFieldKernelError> {
    let mut energies = Vec::new();
    force_field.calc_energy_current(Some(&mut energies))?;
    Ok(energies)
}

#[cfg(test)]
pub(super) fn cf3d_frag_accept_contribution_identities(
    force_field: &ForceField<'_>,
) -> Vec<Cf3dFragAcceptContributionIdentity> {
    force_field
        .contributions
        .iter()
        .map(|contribution| contribution.cf3d_frag_accept_test_identity())
        .collect()
}

#[cfg(test)]
pub(super) fn cf3d_frag_accept_current_energy_and_gradient(
    force_field: &mut ForceField<'_>,
) -> Result<(f64, Vec<f64>), ForceFieldKernelError> {
    let mut gradient = vec![0.0; force_field.dimension as usize * force_field.positions.len()];
    let energy = force_field.calc_energy_current(None)?;
    force_field.calc_grad_current(&mut gradient)?;
    Ok((energy, gradient))
}

#[cfg(test)]
pub(super) fn cf3d_bld_integration_minimize(
    force_field: &mut ForceField<'_>,
    max_its: u32,
    force_tol: f64,
    energy_tol: f64,
) -> Result<i32, ForceFieldKernelError> {
    force_field.minimize(max_its, force_tol, energy_tol)
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use super::{
        AngleIndexArgument, AngleRangeBound, BondIndexArgument, EvaluationContext, ForceField,
        ForceFieldContribution, ForceFieldIndexArgument, ForceFieldKernelError,
        TorsionIndexArgument, distance::DistanceConstraintContrib,
        distances::DistanceConstraintContribs,
    };
    use crate::optimizer::OptimizerSnapshot;
    use crate::uff::{bond::BondStretchContrib, params::AtomicParams};

    #[test]
    fn uff_error_e01_kernel_variants_are_leaf_errors() {
        // These manually supplied variants verify Error trait dispatch only;
        // the source-driven optimizer failure remains covered separately.
        let errors = [
            ForceFieldKernelError::NoPoints,
            ForceFieldKernelError::NoDistanceMatrix,
            ForceFieldKernelError::MatrixSizeMismatch,
            ForceFieldKernelError::NotInitialized,
            ForceFieldKernelError::BadBounds,
            ForceFieldKernelError::BadBondOrder,
            ForceFieldKernelError::PositionCountMismatch {
                num_points: 17,
                position_count: 29,
            },
            ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: 31,
                upper_bound: 47,
            },
            ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::Second,
                index: 37,
                upper_bound: 53,
            },
            ForceFieldKernelError::AngleDegeneratePoints,
            ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::Third,
                index: 41,
                upper_bound: 59,
            },
            ForceFieldKernelError::AngleBadOrder { order: 17 },
            ForceFieldKernelError::TorsionIndexOutOfRange {
                argument: TorsionIndexArgument::Fourth,
                index: 43,
                upper_bound: 61,
            },
            ForceFieldKernelError::TorsionDegeneratePoints,
            ForceFieldKernelError::TorsionBadHybridizations,
            ForceFieldKernelError::TorsionBadOrder { order: 19 },
            ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Maximum,
            },
            ForceFieldKernelError::AngleBoundsOrder,
            ForceFieldKernelError::PackedAngleBoundsOrder,
            ForceFieldKernelError::TorsionBoundsOrder,
            ForceFieldKernelError::BadIndex,
            ForceFieldKernelError::OptimizerBadTolerance,
            ForceFieldKernelError::OptimizerBadDirection,
            ForceFieldKernelError::BadFixedPoint {
                index: -23,
                upper_bound: 67,
            },
            ForceFieldKernelError::TransferPostcondition,
        ];

        assert_eq!(errors.len(), 25);
        for error in &errors {
            let erased: &(dyn std::error::Error + 'static) = error;
            assert_eq!(erased.downcast_ref::<ForceFieldKernelError>(), Some(error));
            assert!(erased.source().is_none());
            assert_eq!(error.to_string(), format!("{error:?}"));
        }
    }

    #[test]
    fn uff_thread_t02_send_compile_proofs_cover_field_and_contributions() {
        fn assert_send<T: Send>() {}

        assert_send::<ForceField<'static>>();
        assert_send::<Box<dyn ForceFieldContribution>>();
        assert_send::<BondStretchContrib>();
        assert_send::<super::angles::AngleConstraintContribs>();
        assert_send::<super::distance::DistanceConstraintContrib>();
        assert_send::<super::distances::DistanceConstraintContribs>();
        assert_send::<super::position::PositionConstraintContrib>();
        assert_send::<super::torsion::TorsionConstraintContrib>();
        assert_send::<Cell<usize>>();
        assert_send::<FailOnEnergyCallContribution>();

        // The trait supertrait is checked at each concrete production impl,
        // including UFF terms that remain private to their owning modules.
    }

    #[test]
    fn cf3d_bld_b06_torsion_error_metadata_matches_source() {
        assert_eq!(
            ForceFieldKernelError::TorsionDegeneratePoints.source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            ForceFieldKernelError::TorsionDegeneratePoints.source_message(),
            "degenerate points"
        );
        assert_eq!(
            ForceFieldKernelError::TorsionDegeneratePoints.source_expression(),
            Some(
                "(idx1 != idx2 && idx1 != idx3 && idx1 != idx4 && idx2 != idx3 && idx2 != idx4 && idx3 != idx4)"
                    .to_owned()
            )
        );

        assert_eq!(
            ForceFieldKernelError::TorsionBadHybridizations.source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            ForceFieldKernelError::TorsionBadHybridizations.source_message(),
            "bad hybridizations"
        );
        assert_eq!(
            ForceFieldKernelError::TorsionBadHybridizations.source_expression(),
            Some(
                "(hyb2 == RDKit::Atom::SP2 || hyb2 == RDKit::Atom::SP3) && (hyb3 == RDKit::Atom::SP2 || hyb3 == RDKit::Atom::SP3)"
                    .to_owned()
            )
        );

        let bad_order = ForceFieldKernelError::TorsionBadOrder { order: 4 };
        assert_eq!(bad_order.source_category(), "Pre-condition Violation");
        assert_eq!(bad_order.source_message(), "bad order");
        assert_eq!(
            bad_order.source_expression(),
            Some("d_order == 2 || d_order == 3 || d_order == 6".to_owned())
        );
    }

    fn initialized_bond_test_field<'a>(
        first: &'a mut Vec<f64>,
        second: &'a mut Vec<f64>,
    ) -> ForceField<'a> {
        let mut force_field = ForceField::new(3);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field
            .initialize()
            .expect("two-position source fixture initializes");
        force_field
    }

    #[test]
    fn uff_serial_s01_moves_owned_state_and_releases_old_borrow() {
        let mut old_first = vec![0.0, 0.0, 0.0];
        let mut old_second = vec![2.0, 0.0, 0.0];
        let mut next_first = vec![4.0, 1.0, 0.0];
        let mut next_second = vec![7.0, 1.0, 0.0];
        let next_address = next_first.as_ptr();
        let mut field = initialized_bond_test_field(&mut old_first, &mut old_second);
        attach_u06_bond_contribution(&mut field);
        field.fixed_points.extend([0, 1]);
        assert_eq!(field.distance(0, 1, None), Ok(2.0));
        let contribution = &*field.contributions[0] as *const dyn ForceFieldContribution;
        let contributions = field.contributions.as_ptr();
        let matrix = field.distance_matrix.as_ptr();
        let fixed = field.fixed_points.as_ptr();
        super::cf3d_uff_one_kernel_counts_reset();
        let rebound =
            field.rebind_positions(vec![next_first.as_mut_slice(), next_second.as_mut_slice()]);
        old_first[0] = -9.0;
        old_second[0] = -8.0;
        assert_eq!(old_first[0], -9.0);
        assert_eq!(old_second[0], -8.0);
        assert!(std::ptr::eq(contribution, &*rebound.contributions[0]));
        assert_eq!(rebound.contributions.as_ptr(), contributions);
        assert_eq!(rebound.distance_matrix.as_ptr(), matrix);
        assert_eq!(rebound.fixed_points.as_ptr(), fixed);
        assert_eq!(rebound.fixed_points, [0, 1]);
        assert_eq!(rebound.positions()[0].as_ptr(), next_address);
        assert_eq!(rebound.dimension, 3);
        assert_eq!(rebound.num_points, 2);
        assert_eq!(rebound.matrix_size, 3);
        assert!(rebound.initialized && rebound.matrix_allocated);
        assert_eq!(super::cf3d_uff_one_kernel_counts(), (0, 0, 0));
    }

    #[test]
    fn uff_block_l01_release_position_borrows_preserves_field_storage() {
        const POSITION_LENGTHS: [usize; 4] = [0, 1, 2, 7];
        let mut release_calls = 0;

        for position_length in POSITION_LENGTHS {
            for initialized in [false, true] {
                // An initialized field with no position entries is a real
                // rebinding intermediate: initialize one point, then clear
                // its old borrowed view before the new rows are installed.
                let binding_length = if initialized && position_length == 0 {
                    1
                } else {
                    position_length
                };
                let (
                    released,
                    scalar_state,
                    cache_address,
                    cache_capacity,
                    cache_values,
                    fixed_address,
                    fixed_capacity,
                    fixed_values,
                    contributions_address,
                    contributions_capacity,
                    contribution_address,
                ) = {
                    let mut old_rows = (0..binding_length)
                        .map(|row| vec![row as f64 + 0.5, 1.25, -2.5])
                        .collect::<Vec<_>>();
                    let mut field = ForceField::new(3);
                    field
                        .positions_mut()
                        .extend(old_rows.iter_mut().map(|row| row.as_mut_slice()));
                    field.add_contribution(Box::new(FailOnGradientContribution));
                    field.fixed_points.extend(0..binding_length as i32);
                    if initialized {
                        field
                            .initialize()
                            .expect("every initialized fixture has at least one real row");
                    }

                    let scalar_state = (
                        field.dimension,
                        field.initialized,
                        field.num_points,
                        field.matrix_allocated,
                        field.matrix_size,
                    );
                    let position_address = field.positions().as_ptr() as usize;
                    let position_capacity = field.positions_mut().capacity();
                    let cache_address = field.distance_matrix.as_ptr() as usize;
                    let cache_capacity = field.distance_matrix.capacity();
                    let cache_values = field.distance_matrix.clone();
                    let fixed_address = field.fixed_points.as_ptr() as usize;
                    let fixed_capacity = field.fixed_points.capacity();
                    let fixed_values = field.fixed_points.clone();
                    let contributions_address = field.contributions.as_ptr() as usize;
                    let contributions_capacity = field.contributions.capacity();
                    let contribution_address =
                        &*field.contributions[0] as *const dyn ForceFieldContribution;

                    super::cf3d_uff_one_kernel_counts_reset();
                    let mut released = if initialized && position_length == 0 {
                        field
                            .clear_positions_and_shorten_lifetime()
                            .release_position_borrows()
                    } else {
                        field.release_position_borrows()
                    };
                    release_calls += 1;

                    assert!(released.positions().is_empty());
                    assert_eq!(released.positions_mut().capacity(), position_capacity);
                    if position_capacity > 0 {
                        assert_eq!(
                            released.positions().as_ptr() as usize,
                            position_address,
                            "positions buffer remains allocated at length={position_length}, initialized={initialized}"
                        );
                    }

                    // These writes compile and run while the released field
                    // is alive, proving its type no longer borrows old_rows.
                    for (row_index, row) in old_rows.iter_mut().enumerate() {
                        row[0] = -10.0 - row_index as f64;
                    }
                    for (row_index, row) in old_rows.iter().enumerate() {
                        assert_eq!(row[0], -10.0 - row_index as f64);
                    }
                    assert_eq!(super::cf3d_uff_one_kernel_counts(), (0, 0, 0));

                    (
                        released,
                        scalar_state,
                        cache_address,
                        cache_capacity,
                        cache_values,
                        fixed_address,
                        fixed_capacity,
                        fixed_values,
                        contributions_address,
                        contributions_capacity,
                        contribution_address,
                    )
                };

                // The old coordinate rows have expired at this point; the
                // released ForceField remains usable without that lifetime.
                assert_eq!(
                    (
                        released.dimension,
                        released.initialized,
                        released.num_points,
                        released.matrix_allocated,
                        released.matrix_size,
                    ),
                    scalar_state
                );
                assert_eq!(released.distance_matrix.as_ptr() as usize, cache_address);
                assert_eq!(released.distance_matrix.capacity(), cache_capacity);
                assert_eq!(released.distance_matrix, cache_values);
                assert_eq!(released.fixed_points.as_ptr() as usize, fixed_address);
                assert_eq!(released.fixed_points.capacity(), fixed_capacity);
                assert_eq!(released.fixed_points, fixed_values);
                assert_eq!(
                    released.contributions.as_ptr() as usize,
                    contributions_address
                );
                assert_eq!(released.contributions.capacity(), contributions_capacity);
                assert_eq!(released.contributions.len(), 1);
                assert!(std::ptr::eq(
                    contribution_address,
                    &*released.contributions[0]
                ));
            }
        }

        assert_eq!(release_calls, 8);
    }

    #[test]
    fn uff_serial_s01_keeps_dimension_and_defers_row_count_to_initialize() {
        for dimension in [2, 3] {
            let mut old = vec![0.0; dimension as usize];
            let mut first = vec![1.0; dimension as usize];
            let mut second = vec![2.0; dimension as usize];
            let mut field = ForceField::new(dimension);
            field.positions_mut().push(old.as_mut_slice());
            field.initialize().unwrap();
            let mut rebound =
                field.rebind_positions(vec![first.as_mut_slice(), second.as_mut_slice()]);
            assert_eq!(rebound.dimension(), dimension);
            assert_eq!(rebound.positions().len(), 2);
            assert_eq!(rebound.num_points(), 1);
            rebound.initialize().unwrap();
            assert_eq!(rebound.num_points(), 2);
        }
    }

    #[test]
    fn uff_serial_s02_reinitialize_invalidates_same_and_changed_geometry_cache() {
        for separation in [2.0, 5.0] {
            let mut old_first = vec![0.0, 0.0, 0.0];
            let mut old_second = vec![2.0, 0.0, 0.0];
            let mut first = vec![0.0, 0.0, 0.0];
            let mut second = vec![separation, 0.0, 0.0];
            let mut field = initialized_bond_test_field(&mut old_first, &mut old_second);
            assert_eq!(field.distance(0, 1, None), Ok(2.0));
            let mut rebound =
                field.rebind_positions(vec![first.as_mut_slice(), second.as_mut_slice()]);
            assert_eq!(rebound.distance(0, 1, None), Ok(2.0));
            rebound.initialize().unwrap();
            assert_eq!(rebound.distance_matrix, [-1.0; 3]);
            assert_eq!(rebound.distance(0, 1, None), Ok(separation));
            assert!(rebound.initialized);
        }
    }

    #[test]
    fn uff_serial_s02_initialize_defers_component_access_and_zero_points_fails() {
        // ForceField.cpp:239-250 does not access Point components, even if a
        // Point2D has been installed in a dimension-3 field. Do not invent an
        // eager component error here. Complete serial inputs use [f64; 3].
        let mut short_row = [1.0, 2.0];
        let mut field = ForceField::new(3).rebind_positions(vec![&mut short_row]);
        assert_eq!(field.initialize(), Ok(()));
        assert_eq!(field.dimension(), 3);
        assert_eq!(field.num_points(), 1);
        assert_eq!(field.positions()[0].len(), 2);
        let mut empty = field.rebind_positions(Vec::new());
        assert_eq!(empty.initialize(), Err(ForceFieldKernelError::NoPoints));
        assert!(!empty.initialized);
        assert_eq!(empty.num_points(), 0);
        assert_eq!(empty.matrix_size, 0);
        assert!(empty.distance_matrix.is_empty());
        assert!(empty.matrix_allocated);
    }

    fn uff_sp3_carbon_params() -> AtomicParams {
        AtomicParams {
            r1: 0.757,
            theta0: 0.0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1: 1.912,
            v1: 0.0,
            u1: 0.0,
            gmp_xi: 5.343,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn assert_u05_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:.17}, got {actual:.17} (tolerance {tolerance})"
        );
    }

    fn attach_u06_bond_contribution(force_field: &mut ForceField<'_>) -> BondStretchContrib {
        let carbon = uff_sp3_carbon_params();
        let contribution =
            BondStretchContrib::new(force_field.positions(), 0, 1, 1.0, &carbon, &carbon)
                .expect("valid pinned sp3 C-C contribution");
        force_field.contributions.push(Box::new(contribution));
        contribution
    }

    fn evaluate_u06_gradient(
        force_field: &mut ForceField<'_>,
        coordinates: &[f64],
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // The pinned test reinitializes the force field cache before each new
        // coordinate state; calcGrad itself retains the source cache.
        force_field.init_distance_matrix()?;
        force_field.calc_grad(coordinates, gradient)
    }

    fn assert_u06_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:.17}, got {actual:.17} (tolerance {tolerance})"
        );
    }

    #[test]
    fn cf3d_u05_constructor_preserves_source_check_order_and_errors() {
        // Pinned RDKit Code/ForceField/UFF/BondStretch.cpp constructor and
        // Code/RDGeneral/Invariant.h PRECONDITION/URANGE_CHECK definitions.
        let mut first = vec![0.0; 3];
        let mut second = vec![0.0; 3];
        let force_field = initialized_bond_test_field(&mut first, &mut second);
        let carbon = uff_sp3_carbon_params();

        // Both indices and bond order fail; the constructor checks idx1 first.
        let first_error =
            BondStretchContrib::new(force_field.positions(), 2, 2, 0.0, &carbon, &carbon)
                .expect_err("source checks first endpoint range before later inputs");
        assert_eq!(
            first_error,
            ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::First,
                index: 2,
                upper_bound: 2,
            }
        );
        assert_eq!(first_error.source_category(), "Range Error");
        assert_eq!(first_error.source_message(), "idx1");
        assert_eq!(first_error.source_expression().as_deref(), Some("2 < 2"));

        let second_error =
            BondStretchContrib::new(force_field.positions(), 0, 2, 0.0, &carbon, &carbon)
                .expect_err("source checks second endpoint before bond-order calculation");
        assert_eq!(
            second_error,
            ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::Second,
                index: 2,
                upper_bound: 2,
            }
        );
        assert_eq!(second_error.source_category(), "Range Error");
        assert_eq!(second_error.source_message(), "idx2");
        assert_eq!(second_error.source_expression().as_deref(), Some("2 < 2"));

        // PRECONDITION(bondOrder > 0) rejects zero, negative and unordered NaN.
        for bond_order in [0.0, -0.0, -1.0, f64::NAN] {
            let error = BondStretchContrib::new(
                force_field.positions(),
                0,
                1,
                bond_order,
                &carbon,
                &carbon,
            )
            .expect_err("nonpositive or unordered source bond order");
            assert_eq!(error, ForceFieldKernelError::BadBondOrder);
            assert_eq!(error.source_category(), "Pre-condition Violation");
            assert_eq!(error.source_message(), "bad bond order");
            assert_eq!(error.source_expression().as_deref(), Some("bondOrder > 0"));
        }

        let empty_field = ForceField::new(3);
        assert_eq!(
            BondStretchContrib::new(empty_field.positions(), 0, 0, 1.0, &carbon, &carbon),
            Err(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::First,
                index: 0,
                upper_bound: 0,
            })
        );
    }

    #[test]
    fn cf3d_u05_energy_matches_fixed_equilibrium_compressed_and_stretched_values() {
        // Pinned RDKit testUFFForceField.cpp::testUFF1 fixes the sp3 C-C
        // rest length at 1.514 and force constant at 699.5918; testUFF2
        // checks zero energy at rest and positive energy at zero distance.
        let mut first = vec![0.0; 3];
        let mut second = vec![0.0; 3];
        let mut force_field = initialized_bond_test_field(&mut first, &mut second);
        let carbon = uff_sp3_carbon_params();
        let contribution =
            BondStretchContrib::new(force_field.positions(), 0, 1, 1.0, &carbon, &carbon)
                .expect("valid pinned sp3 C-C contribution");

        const SOURCE_REST_LENGTH: f64 = 1.514;
        const SOURCE_FORCE_CONSTANT: f64 = 699.5918;
        const SOURCE_QUADRATIC_ENERGY: f64 = 0.5 * SOURCE_FORCE_CONSTANT * 0.25 * 0.25;

        for (distance, expected_energy) in [
            (0.0, None),
            (SOURCE_REST_LENGTH - 0.25, Some(SOURCE_QUADRATIC_ENERGY)),
            (SOURCE_REST_LENGTH, Some(0.0)),
            (SOURCE_REST_LENGTH + 0.25, Some(SOURCE_QUADRATIC_ENERGY)),
            (1.814, Some(31.4816)),
        ] {
            // Source calcEnergy(pos) resets the one triangular distance cache
            // before evaluating the contribution.
            force_field
                .init_distance_matrix()
                .expect("initialized source cache resets");
            let coordinates = [0.0, 0.0, 0.0, distance, 0.0, 0.0];
            let energy = {
                let mut context = force_field.evaluation_context(&coordinates);
                contribution
                    .get_energy(&mut context)
                    .expect("valid source bond energy")
            };

            assert_u05_close(force_field.distance_matrix[1], distance, 1.0e-12);
            match expected_energy {
                Some(expected) => assert_u05_close(energy, expected, 1.0e-4),
                None => assert!(energy > 0.0, "source zero-distance energy is positive"),
            }
        }

        // The contribution constructor accepts an owner with positions before
        // initialization; source getEnergy then propagates distance's exact
        // not-initialized precondition.
        let mut uninitialized_first = vec![0.0; 3];
        let mut uninitialized_second = vec![0.0; 3];
        let mut uninitialized_field = ForceField::new(3);
        uninitialized_field.positions_mut().extend([
            uninitialized_first.as_mut_slice(),
            uninitialized_second.as_mut_slice(),
        ]);
        let uninitialized_contribution =
            BondStretchContrib::new(uninitialized_field.positions(), 0, 1, 1.0, &carbon, &carbon)
                .expect("source constructor does not require initialized owner");
        let mut context = uninitialized_field.evaluation_context(&[0.0; 6]);
        assert_eq!(
            uninitialized_contribution.get_energy(&mut context),
            Err(ForceFieldKernelError::NotInitialized)
        );
    }

    #[test]
    fn cf3d_u06_zero_and_nan_distance_use_source_floor_additively() {
        // Pinned RDKit Code/ForceField/UFF/testUFFForceField.cpp::testUFF2
        // exercises coincident endpoints and requires all six gradient values
        // to be nonzero. The source branch is `if (dist > 0.0)`, so unordered
        // NaN distance takes the same fixed fallback.
        let mut first = vec![0.0; 3];
        let mut second = vec![0.0; 3];
        let mut force_field = initialized_bond_test_field(&mut first, &mut second);
        let _contribution = attach_u06_bond_contribution(&mut force_field);

        let mut gradient = [10.0, -1.0, 2.0, -2.0, 5.0, 0.5];
        evaluate_u06_gradient(&mut force_field, &[0.0; 6], &mut gradient)
            .expect("initialized source context");
        for (actual, expected) in gradient.into_iter().zip([
            16.995918, 5.995918, 8.995918, -8.995918, -1.995918, -6.495918,
        ]) {
            assert_u06_close(actual, expected, 1.0e-4);
        }

        let mut nan_gradient = [0.0; 6];
        evaluate_u06_gradient(
            &mut force_field,
            &[0.0, 0.0, 0.0, f64::NAN, 0.0, 0.0],
            &mut nan_gradient,
        )
        .expect("NaN coordinates still follow source comparison branch");
        for (actual, expected) in nan_gradient.into_iter().zip([
            6.995918, 6.995918, 6.995918, -6.995918, -6.995918, -6.995918,
        ]) {
            assert_u06_close(actual, expected, 1.0e-4);
        }
    }

    #[test]
    fn cf3d_u06_equilibrium_compression_and_stretch_match_pinned_values() {
        // RDKit testUFF2 fixes r0=1.514, k01=699.5918, the 1.814 stretch
        // gradient at +/-209.8775, and repeats the case on the z axis.
        let mut first = vec![0.0; 3];
        let mut second = vec![0.0; 3];
        let mut force_field = initialized_bond_test_field(&mut first, &mut second);
        let _contribution = attach_u06_bond_contribution(&mut force_field);

        let mut equilibrium_gradient = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let equilibrium_before = equilibrium_gradient;
        evaluate_u06_gradient(
            &mut force_field,
            &[0.0, 0.0, 0.0, 1.514, 0.0, 0.0],
            &mut equilibrium_gradient,
        )
        .expect("equilibrium source gradient");
        assert_eq!(equilibrium_gradient, equilibrium_before);

        for (distance, expected_first, expected_second) in
            [(1.214, 209.8775, -209.8775), (1.814, -209.8775, 209.8775)]
        {
            let mut gradient = [0.0; 6];
            evaluate_u06_gradient(
                &mut force_field,
                &[0.0, 0.0, 0.0, distance, 0.0, 0.0],
                &mut gradient,
            )
            .expect("source x-axis gradient");
            assert_u06_close(gradient[0], expected_first, 1.0e-4);
            assert_u06_close(gradient[3], expected_second, 1.0e-4);
            assert_eq!(
                [gradient[1], gradient[2], gradient[4], gradient[5]],
                [0.0; 4]
            );
            assert_u06_close(gradient[0] + gradient[3], 0.0, 1.0e-12);
        }

        let mut z_gradient = [0.0; 6];
        evaluate_u06_gradient(
            &mut force_field,
            &[0.0, 0.0, 0.0, 0.0, 0.0, 1.814],
            &mut z_gradient,
        )
        .expect("source z-axis gradient");
        assert_u06_close(z_gradient[2], -209.8775, 1.0e-4);
        assert_u06_close(z_gradient[5], 209.8775, 1.0e-4);
        assert_eq!(
            [z_gradient[0], z_gradient[1], z_gradient[3], z_gradient[4]],
            [0.0; 4]
        );
    }

    #[test]
    fn cf3d_u06_gradient_uses_each_source_component_and_preserves_additive_values() {
        // A fixed non-axis-aligned source input makes the three-coordinate
        // order observable independently of the x/z fixtures in testUFF2.
        let mut first = vec![0.0; 3];
        let mut second = vec![0.0; 3];
        let mut force_field = initialized_bond_test_field(&mut first, &mut second);
        let _contribution = attach_u06_bond_contribution(&mut force_field);
        let mut gradient = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0];

        evaluate_u06_gradient(
            &mut force_field,
            &[0.0, 0.0, 0.0, 1.0, 2.0, 2.0],
            &mut gradient,
        )
        .expect("source three-component gradient");

        let component = 346.5311382666667;
        for (actual, expected) in gradient.into_iter().zip([
            10.0 - component,
            20.0 - 2.0 * component,
            30.0 - 2.0 * component,
            40.0 + component,
            50.0 + 2.0 * component,
            60.0 + 2.0 * component,
        ]) {
            assert_u06_close(actual, expected, 1.0e-4);
        }
        for axis in 0..3 {
            assert_u06_close(
                (gradient[axis] - [10.0, 20.0, 30.0][axis])
                    + (gradient[axis + 3] - [40.0, 50.0, 60.0][axis]),
                0.0,
                1.0e-10,
            );
        }
    }

    #[test]
    fn cf3d_u06_distance_error_precedes_gradient_mutation() {
        let mut first = vec![0.0; 3];
        let mut second = vec![0.0; 3];
        let mut force_field = ForceField::new(3);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        let carbon = uff_sp3_carbon_params();
        let contribution =
            BondStretchContrib::new(force_field.positions(), 0, 1, 1.0, &carbon, &carbon)
                .expect("constructor precedes source initialization");
        let mut gradient = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let before = gradient;
        let mut context = force_field.evaluation_context(&[0.0; 6]);

        assert_eq!(
            contribution.get_grad(&mut context, &mut gradient),
            Err(ForceFieldKernelError::NotInitialized)
        );
        assert_eq!(gradient, before);
    }

    #[derive(Clone)]
    struct SquaredDistanceContribution {
        first: u32,
        second: u32,
        scale: f64,
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
            for component in 0..context.dimension as usize {
                let first_index = context.dimension as usize * self.first as usize + component;
                let second_index = context.dimension as usize * self.second as usize + component;
                let derivative = 2.0
                    * self.scale
                    * (context.coordinates[first_index] - context.coordinates[second_index]);
                gradient[first_index] += derivative;
                gradient[second_index] -= derivative;
            }
            Ok(())
        }

        fn copy(&self) -> Box<dyn ForceFieldContribution> {
            Box::new(self.clone())
        }
    }

    struct FailOnEnergyCallContribution {
        calls: Cell<usize>,
        fail_on: usize,
    }

    impl ForceFieldContribution for FailOnEnergyCallContribution {
        fn get_energy(
            &self,
            context: &mut EvaluationContext<'_>,
        ) -> Result<f64, ForceFieldKernelError> {
            let call = self.calls.get() + 1;
            self.calls.set(call);
            if call == self.fail_on {
                // Trigger the kernel's real source-index error at a fixed
                // callback count; this helper only injects the failure point.
                let _ = context.distance(2, 0)?;
            }
            SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }
            .get_energy(context)
        }

        fn get_grad(
            &self,
            context: &mut EvaluationContext<'_>,
            gradient: &mut [f64],
        ) -> Result<(), ForceFieldKernelError> {
            SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }
            .get_grad(context, gradient)
        }

        fn copy(&self) -> Box<dyn ForceFieldContribution> {
            Box::new(Self {
                calls: self.calls.clone(),
                fail_on: self.fail_on,
            })
        }
    }

    struct FailOnGradientContribution;

    impl ForceFieldContribution for FailOnGradientContribution {
        fn get_energy(
            &self,
            context: &mut EvaluationContext<'_>,
        ) -> Result<f64, ForceFieldKernelError> {
            SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }
            .get_energy(context)
        }

        fn get_grad(
            &self,
            context: &mut EvaluationContext<'_>,
            _gradient: &mut [f64],
        ) -> Result<(), ForceFieldKernelError> {
            // Inject a genuine kernel range error from the gradient callback.
            let _ = context.distance(2, 0)?;
            Ok(())
        }

        fn copy(&self) -> Box<dyn ForceFieldContribution> {
            Box::new(Self)
        }
    }

    #[test]
    fn uff_one_u02_iteration_limits_keep_source_status_and_coordinates() {
        // RDKit ForceField.cpp:257-282 and BFGSOpt.h:184-327: both callbacks
        // run before the maxIts loop; each normal 0/1 status gathers positions.
        let mut zero_first = [0.0];
        let mut zero_second = [2.0];
        let mut zero_field = ForceField::new(1);
        zero_field
            .positions_mut()
            .extend([zero_first.as_mut_slice(), zero_second.as_mut_slice()]);
        zero_field.initialize().unwrap();
        zero_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        assert_eq!(zero_field.minimize(0, 0.1, 0.0), Ok(1));
        assert_eq!(zero_first, [0.0]);
        assert_eq!(zero_second, [2.0]);

        // This is the pinned one-iteration source case: BFGS accepts its
        // 0.4-unit first step but reports that more work is required.
        let mut one_first = [0.0];
        let mut one_second = [2.0];
        let mut one_field = ForceField::new(1);
        one_field
            .positions_mut()
            .extend([one_first.as_mut_slice(), one_second.as_mut_slice()]);
        one_field.initialize().unwrap();
        one_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        assert_eq!(one_field.minimize(1, 0.1, 0.0), Ok(1));
        assert_eq!(one_first, [0.4]);
        assert_eq!(one_second, [1.6]);

        // A longer source run reaches status 0 and records at least two
        // distinct iteration snapshots; final coordinates satisfy the
        // zero-separation minimum of this fixed squared-distance contribution.
        let mut many_first = [0.0];
        let mut many_second = [2.0];
        let mut many_field = ForceField::new(1);
        many_field
            .positions_mut()
            .extend([many_first.as_mut_slice(), many_second.as_mut_slice()]);
        many_field.initialize().unwrap();
        many_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        let mut snapshots = Vec::new();
        assert_eq!(
            many_field.minimize_with_snapshots(1, Some(&mut snapshots), 20, 1.0e-4, 1.0e-6,),
            Ok(0)
        );
        assert!(snapshots.len() >= 2, "snapshots={snapshots:?}");
        assert!((many_first[0] - many_second[0]).abs() < 1.0e-4);
    }

    #[test]
    fn uff_one_u02_energy_callback_errors_keep_type_and_do_not_gather() {
        let mut first = [0.0];
        let mut second = [2.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(FailOnEnergyCallContribution {
                calls: Cell::new(0),
                fail_on: 1,
            }));

        assert_eq!(
            force_field.minimize(4, 0.1, 0.0),
            Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(first, [0.0]);
        assert_eq!(second, [2.0]);
    }

    #[test]
    fn uff_one_u02_gradient_callback_errors_keep_type_and_do_not_gather() {
        let mut first = [0.0];
        let mut second = [2.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(FailOnGradientContribution));

        assert_eq!(
            force_field.minimize(4, 0.1, 0.0),
            Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(first, [0.0]);
        assert_eq!(second, [2.0]);
    }

    #[test]
    fn cf3d_opt_bridge_visibility_sibling_driver_snapshot() {
        use crate::optimizer::{OptimizerSnapshot, minimize};

        // RDKit BFGSOpt.h:184-327: this source-fixed displacement is above
        // MOVETOL and below TOLX, so the driver records its terminal snapshot.
        let mut positions = [1.1e-7];
        let mut iterations = 0;
        let mut final_energy = 0.0;
        let mut energy_calls = 0;
        let mut energy = |_: &mut (), point: &mut [f64]| -> Result<f64, std::convert::Infallible> {
            energy_calls += 1;
            Ok(0.5 * point[0] * point[0])
        };
        let mut gradient_calls = 0;
        let mut gradient = |_: &mut (),
                            point: &mut [f64],
                            grad: &mut [f64]|
         -> Result<f64, std::convert::Infallible> {
            gradient_calls += 1;
            grad[0] = point[0];
            Ok(1.0)
        };
        let mut snapshots: Vec<OptimizerSnapshot> = Vec::new();

        let status = minimize(
            &mut (),
            &mut positions,
            1.0e-6,
            &mut iterations,
            &mut final_energy,
            &mut energy,
            &mut gradient,
            2,
            Some(&mut snapshots),
            0.0,
            4,
        )
        .unwrap();

        assert_eq!(status, 0);
        assert_eq!(energy_calls, 2);
        assert_eq!(gradient_calls, 1);
        assert_eq!(iterations, 1);
        assert_eq!(final_energy, 0.0);
        assert_eq!(positions, [0.0]);
        assert_eq!(snapshots.len(), 1);
        assert_eq!(snapshots[0].positions, [0.0]);
        assert_eq!(snapshots[0].energy, 0.0);
    }

    #[test]
    fn cf3d_f13_minimize_source_preconditions_precede_empty_term_return() {
        // RDKit source: ForceField.cpp:257-265.
        let mut uninitialized = ForceField::new(1);
        assert_eq!(
            uninitialized.minimize(200, 1.0e-4, 1.0e-6),
            Err(ForceFieldKernelError::NotInitialized)
        );

        let mut point = vec![7.0];
        let mut force_field = ForceField::new(1);
        force_field.positions_mut().push(point.as_mut_slice());
        force_field.initialize().unwrap();
        force_field.num_points = 2;

        assert_eq!(
            force_field.minimize(200, 1.0e-4, 1.0e-6),
            Err(ForceFieldKernelError::PositionCountMismatch {
                num_points: 2,
                position_count: 1,
            })
        );
        assert_eq!(
            ForceFieldKernelError::PositionCountMismatch {
                num_points: 2,
                position_count: 1,
            }
            .source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            ForceFieldKernelError::PositionCountMismatch {
                num_points: 2,
                position_count: 1,
            }
            .source_message(),
            "size mismatch"
        );
        assert_eq!(point, [7.0]);
    }

    #[test]
    fn cf3d_f13_minimize_empty_contributions_keep_source_noop_and_snapshots() {
        // RDKit source: ForceField.cpp:257-282; BFGSOpt.h:184-327.
        let mut point = vec![7.0];
        let mut force_field = ForceField::new(1);
        force_field.positions_mut().push(point.as_mut_slice());
        force_field.initialize().unwrap();
        let mut snapshots = vec![OptimizerSnapshot {
            positions: vec![-1.0],
            energy: -1.0,
        }];

        assert_eq!(
            force_field.minimize_with_snapshots(1, Some(&mut snapshots), 1, 1.0e-4, 1.0e-6,),
            Ok(0)
        );
        assert_eq!(point, [7.0]);
        assert_eq!(
            snapshots,
            [OptimizerSnapshot {
                positions: vec![-1.0],
                energy: -1.0,
            }]
        );
    }

    #[test]
    fn cf3d_f13_minimize_default_overload_preserves_fixed_point_and_gathers_success() {
        // RDKit source: ForceField.cpp:252-282,106-145; ForceField.h:153-175.
        let mut fixed = vec![0.0];
        let mut moving = vec![2.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([fixed.as_mut_slice(), moving.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        force_field.fixed_points_mut().push(0);

        // These are the pinned ForceField.h defaults; the no-snapshot
        // overload also supplies snapshotFreq=0 and a null snapshot vector.
        assert_eq!(force_field.minimize(200, 1.0e-4, 1.0e-6), Ok(0));
        assert_eq!(fixed, [0.0]);
        assert!(moving[0].abs() < 1.0e-3, "moving={:?}", moving[0]);
    }

    #[test]
    fn cf3d_f13_minimize_status_one_gathers_and_keeps_periodic_snapshot() {
        // RDKit source: ForceField.cpp:257-282; BFGSOpt.h:184-327.
        let mut first = vec![0.0];
        let mut second = vec![2.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        let mut snapshots = Vec::new();

        assert_eq!(
            force_field.minimize_with_snapshots(1, Some(&mut snapshots), 1, 0.1, 0.0,),
            Ok(1)
        );
        assert!((first[0] - 0.4).abs() < 1.0e-12);
        assert!((second[0] - 1.6).abs() < 1.0e-12);
        assert_eq!(snapshots.len(), 1);
        assert!((snapshots[0].positions[0] - first[0]).abs() < 1.0e-12);
        assert!((snapshots[0].positions[1] - second[0]).abs() < 1.0e-12);
        assert!((snapshots[0].energy - 1.44).abs() < 1.0e-12);
    }

    #[test]
    fn cf3d_f13_minimize_callback_error_preserves_positions_and_prior_snapshot() {
        // RDKit source: ForceField.cpp:257-282; BFGSOpt.h:184-327.
        let mut first = vec![0.0];
        let mut second = vec![2.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(FailOnEnergyCallContribution {
                calls: Cell::new(0),
                fail_on: 3,
            }));
        let mut snapshots = Vec::new();

        assert_eq!(
            force_field.minimize_with_snapshots(1, Some(&mut snapshots), 4, 0.1, 0.0,),
            Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(first, [0.0]);
        assert_eq!(second, [2.0]);
        assert_eq!(snapshots.len(), 1);
        assert!((snapshots[0].positions[0] - 0.4).abs() < 1.0e-12);
        assert!((snapshots[0].positions[1] - 1.6).abs() < 1.0e-12);
        assert!((snapshots[0].energy - 1.44).abs() < 1.0e-12);
    }

    #[test]
    fn cf3d_f13_minimize_typed_optimizer_invariants_keep_source_categories() {
        // RDKit source: ForceField.cpp:257-282; BFGSOpt.h:184-327.
        let mut first = vec![0.0];
        let mut second = vec![0.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));

        assert_eq!(
            force_field.minimize(200, 0.0, 1.0e-6),
            Err(ForceFieldKernelError::OptimizerBadTolerance)
        );
        assert_eq!(
            ForceFieldKernelError::OptimizerBadTolerance.source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            ForceFieldKernelError::OptimizerBadTolerance.source_message(),
            "bad tolerance"
        );
        assert_eq!(
            ForceFieldKernelError::OptimizerBadTolerance
                .source_expression()
                .as_deref(),
            Some("gradTol > 0")
        );

        assert_eq!(
            force_field.minimize(200, 1.0e-4, 1.0e-6),
            Err(ForceFieldKernelError::OptimizerBadDirection)
        );
        assert_eq!(
            ForceFieldKernelError::OptimizerBadDirection.source_category(),
            "Invariant Violation"
        );
        assert_eq!(
            ForceFieldKernelError::OptimizerBadDirection.source_message(),
            "bad direction in linearSearch"
        );
        assert_eq!(
            ForceFieldKernelError::OptimizerBadDirection
                .source_expression()
                .as_deref(),
            Some("status >= 0")
        );
        assert_eq!(first, [0.0]);
        assert_eq!(second, [0.0]);
    }

    #[test]
    fn cf3d_f07_ordered_sums_copy_and_both_coordinate_overloads() {
        // RDKit source: ForceField.cpp:284-325 and Contrib.h:18-35.
        let mut first = vec![0.0, 0.0];
        let mut second = vec![3.0, 4.0];
        let mut force_field = ForceField::new(2);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 2.0,
            }));
        let copied = force_field.contributions[0].copy();
        force_field.contributions.push(copied);

        // Source order is left-to-right; the requested terms are cleared then
        // appended in contribution order, while an omitted vector is unused.
        let mut terms = vec![999.0];
        assert_eq!(force_field.calc_energy_current(Some(&mut terms)), Ok(100.0));
        assert_eq!(terms, [25.0, 50.0, 25.0]);
        assert_eq!(force_field.distance_matrix, [-1.0, 5.0, -1.0]);
        assert_eq!(force_field.calc_energy_current(None), Ok(100.0));

        // The explicit overload resets the cache and consumes its borrowed
        // coordinates directly; the next current call reuses that cache.
        assert_eq!(force_field.calc_energy(&[0.0, 0.0, 6.0, 8.0]), Ok(400.0));
        assert_eq!(force_field.distance_matrix, [-1.0, 10.0, -1.0]);
        assert_eq!(force_field.calc_energy_current(None), Ok(400.0));

        let coordinates = [0.0, 0.0, 3.0, 4.0];
        let contribution = SquaredDistanceContribution {
            first: 0,
            second: 1,
            scale: 1.0,
        };
        let mut context = force_field.evaluation_context(&coordinates);
        let mut gradient = [0.0; 4];
        contribution.get_grad(&mut context, &mut gradient).unwrap();
        assert_eq!(gradient, [-6.0, -8.0, 6.0, 8.0]);
    }

    #[test]
    fn cf3d_f07_empty_contributions_preserve_current_output_and_reset_explicit_cache() {
        // RDKit source: ForceField.cpp:284-325; empty current calls return
        // before output mutation, while explicit calls reset before returning.
        let mut force_field = ForceField::new(2);
        let mut terms = vec![7.0];
        assert_eq!(
            force_field.calc_energy_current(Some(&mut terms)),
            Err(ForceFieldKernelError::NotInitialized)
        );
        assert_eq!(
            force_field.calc_energy(&[]),
            Err(ForceFieldKernelError::NotInitialized)
        );
        assert_eq!(terms, [7.0]);

        let mut first = vec![1.0, 2.0];
        let mut second = vec![4.0, 6.0];
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field.distance_matrix[1] = 9.0;

        assert_eq!(force_field.calc_energy_current(Some(&mut terms)), Ok(0.0));
        assert_eq!(terms, [7.0]);
        assert_eq!(force_field.distance_matrix[1], 9.0);
        assert_eq!(force_field.calc_energy(&[]), Ok(0.0));
        assert_eq!(force_field.distance_matrix, [-1.0, -1.0, -1.0]);
    }

    #[test]
    fn cf3d_f07_callback_error_preserves_source_order_and_prior_terms() {
        // RDKit source: ForceField.cpp:284-307 accumulates and appends each
        // contribution before advancing to the next one.
        let mut first = vec![0.0, 0.0];
        let mut second = vec![3.0, 4.0];
        let mut force_field = ForceField::new(2);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 2,
                second: 0,
                scale: 1.0,
            }));

        let mut terms = vec![99.0];
        assert_eq!(
            force_field.calc_energy_current(Some(&mut terms)),
            Err(ForceFieldKernelError::IndexOutOfRange {
                argument: super::ForceFieldIndexArgument::I,
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(terms, [25.0]);
        assert_eq!(force_field.distance_matrix[1], 5.0);
    }

    #[test]
    fn cf3d_f08_accumulates_current_and_explicit_gradients_before_fixed_zeroing() {
        // RDKit source: ForceField.cpp:329-375. calcGrad accumulates into the
        // incoming vector, and fixed points 0, 2, and 4 are zeroed afterward.
        let mut points = [
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![2.0, 0.0],
            vec![3.0, 0.0],
            vec![4.0, 0.0],
        ];
        let mut force_field = ForceField::new(2);
        force_field
            .positions_mut()
            .extend(points.iter_mut().map(Vec::as_mut_slice));
        force_field.initialize().unwrap();
        force_field.fixed_points_mut().extend([0, 2, 4]);
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 1,
                second: 3,
                scale: 1.0,
            }));

        let mut current_gradient = [1.0, 2.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];
        force_field
            .calc_grad_current(&mut current_gradient)
            .unwrap();
        assert_eq!(
            current_gradient,
            [0.0, 0.0, 6.0, 20.0, 0.0, 0.0, 54.0, 60.0, 0.0, 0.0]
        );
        assert_eq!(force_field.distance_matrix[7], 2.0);

        // Repeated calls accumulate into free-point components; fixed points
        // remain zero because source zeroing follows each successful pass.
        force_field
            .calc_grad_current(&mut current_gradient)
            .unwrap();
        assert_eq!(
            current_gradient,
            [0.0, 0.0, 2.0, 20.0, 0.0, 0.0, 58.0, 60.0, 0.0, 0.0]
        );

        // The explicit overload uses its coordinates directly and does not
        // reset the distance cache before evaluating the contribution.
        let explicit_coordinates = [0.0, 0.0, 5.0, 0.0, 2.0, 0.0, 0.0, 0.0, 4.0, 0.0];
        let mut explicit_gradient = [1.0; 10];
        force_field
            .calc_grad(&explicit_coordinates, &mut explicit_gradient)
            .unwrap();
        assert_eq!(
            explicit_gradient,
            [0.0, 0.0, 11.0, 1.0, 0.0, 0.0, -9.0, 1.0, 0.0, 0.0]
        );
        assert_eq!(force_field.distance_matrix[7], 2.0);
    }

    #[test]
    fn cf3d_f08_empty_contributions_leave_gradient_and_cache_untouched() {
        // RDKit source: ForceField.cpp:329-375 returns before contribution
        // work and fixed-point zeroing when the source list is empty.
        let mut force_field = ForceField::new(2);
        let mut gradient = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(
            force_field.calc_grad_current(&mut gradient),
            Err(ForceFieldKernelError::NotInitialized)
        );
        assert_eq!(
            force_field.calc_grad(&[], &mut gradient),
            Err(ForceFieldKernelError::NotInitialized)
        );

        let mut first = vec![0.0, 0.0];
        let mut second = vec![1.0, 0.0];
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field.fixed_points_mut().push(0);
        force_field.distance_matrix[1] = 17.0;

        assert_eq!(force_field.calc_grad_current(&mut gradient), Ok(()));
        assert_eq!(force_field.calc_grad(&[], &mut gradient), Ok(()));
        assert_eq!(gradient, [1.0, 2.0, 3.0, 4.0]);
        assert_eq!(force_field.distance_matrix[1], 17.0);
    }

    #[test]
    fn cf3d_f08_fixed_point_errors_preserve_source_cast_and_iteration_order() {
        // RDKit source: ForceField.cpp:343-350 and 366-374 checks the signed
        // point after conversion to unsigned and zeroes IDs in source order.
        let mut first = vec![0.0, 0.0];
        let mut second = vec![3.0, 4.0];
        let mut force_field = ForceField::new(2);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field.fixed_points_mut().extend([0, -1]);
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));

        let mut current_gradient = [9.0; 4];
        assert_eq!(
            force_field.calc_grad_current(&mut current_gradient),
            Err(ForceFieldKernelError::BadFixedPoint {
                index: -1,
                upper_bound: 2,
            })
        );
        assert_eq!(current_gradient, [0.0, 0.0, 15.0, 17.0]);
        assert_eq!(
            ForceFieldKernelError::BadFixedPoint {
                index: -1,
                upper_bound: 2,
            }
            .source_category(),
            "Invariant Violation"
        );
        assert_eq!(
            ForceFieldKernelError::BadFixedPoint {
                index: -1,
                upper_bound: 2,
            }
            .source_message(),
            "bad fixed point index"
        );
        assert_eq!(
            ForceFieldKernelError::BadFixedPoint {
                index: -1,
                upper_bound: 2,
            }
            .source_expression()
            .as_deref(),
            Some("static_cast<unsigned int>(-1) < d_numPoints")
        );

        force_field.fixed_points_mut().clear();
        force_field.fixed_points_mut().push(2);
        let mut explicit_gradient = [9.0; 4];
        assert_eq!(
            force_field.calc_grad(&[0.0, 0.0, 3.0, 4.0], &mut explicit_gradient),
            Err(ForceFieldKernelError::BadFixedPoint {
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(explicit_gradient, [3.0, 1.0, 15.0, 17.0]);
    }

    #[test]
    fn cf3d_f08_callback_error_returns_before_fixed_point_zeroing() {
        // RDKit source: ForceField.cpp:339-350 calls gradients before applying
        // fixed-point zeroing; a failing earlier callback stops that sequence.
        let mut first = vec![0.0, 0.0];
        let mut second = vec![3.0, 4.0];
        let mut force_field = ForceField::new(2);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field.fixed_points_mut().push(0);
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 2,
                second: 0,
                scale: 1.0,
            }));

        let mut gradient = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(
            force_field.calc_grad_current(&mut gradient),
            Err(ForceFieldKernelError::IndexOutOfRange {
                argument: super::ForceFieldIndexArgument::I,
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(gradient, [-5.0, -6.0, 9.0, 12.0]);
    }

    #[test]
    fn cf3d_f09_copy_moves_safely_and_keeps_kernel_and_term_state_independent() {
        // RDKit source: ForceField.cpp:159-170 and ForceField.h:244-251.
        let mut original_first = vec![0.0, 0.0];
        let mut original_second = vec![3.0, 4.0];
        let mut force_field = ForceField::new(2);
        force_field.positions_mut().extend([
            original_first.as_mut_slice(),
            original_second.as_mut_slice(),
        ]);
        force_field.fixed_points_mut().push(1);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 2.0,
            }));

        // Moving the initialized kernel before any evaluation cannot leave a
        // stale contribution owner reference.
        let mut moved_before_evaluation = force_field;
        let mut original_terms = Vec::new();
        assert_eq!(
            moved_before_evaluation.calc_energy_current(Some(&mut original_terms)),
            Ok(75.0)
        );
        assert_eq!(original_terms, [25.0, 50.0]);
        assert_eq!(moved_before_evaluation.distance_matrix, [-1.0, 5.0, -1.0]);

        let mut copied = moved_before_evaluation.copy();
        assert_eq!(copied.dimension(), 2);
        assert_eq!(copied.num_points(), 2);
        assert!(!copied.initialized);
        assert!(copied.positions().is_empty());
        assert!(copied.fixed_points().is_empty());
        assert!(copied.distance_matrix.is_empty());
        assert!(!copied.matrix_allocated);
        assert_eq!(copied.matrix_size, 0);
        assert_eq!(copied.contributions.len(), 2);

        // The source copy starts without point references. Supplying distinct
        // points and mutating them demonstrates independence from the original.
        let mut copied_first = vec![0.0, 0.0];
        let mut copied_second = vec![6.0, 8.0];
        copied
            .positions_mut()
            .extend([copied_first.as_mut_slice(), copied_second.as_mut_slice()]);
        copied.initialize().unwrap();
        copied.fixed_points_mut().push(0);
        let mut copied_terms = vec![999.0];
        assert_eq!(
            copied.calc_energy_current(Some(&mut copied_terms)),
            Ok(300.0)
        );
        assert_eq!(copied_terms, [100.0, 200.0]);
        assert_eq!(moved_before_evaluation.fixed_points(), [1]);
        assert_eq!(moved_before_evaluation.positions().len(), 2);
        assert_eq!(moved_before_evaluation.distance_matrix[1], 5.0);

        copied.positions[1][0] = 0.0;
        copied.positions[1][1] = 1.0;
        assert_eq!(copied.calc_energy(&[0.0, 0.0, 0.0, 1.0]), Ok(3.0));
        let mut copied_after_mutation_terms = Vec::new();
        assert_eq!(
            copied.calc_energy_current(Some(&mut copied_after_mutation_terms)),
            Ok(3.0)
        );
        assert_eq!(copied_after_mutation_terms, [1.0, 2.0]);
        assert_eq!(moved_before_evaluation.calc_energy_current(None), Ok(75.0));

        // Moving the original after evaluation preserves its cache and term
        // order. Dropping it then leaves the copied contributions callable.
        let mut moved_after_evaluation = moved_before_evaluation;
        assert_eq!(moved_after_evaluation.calc_energy_current(None), Ok(75.0));
        drop(moved_after_evaluation);
        assert_eq!(copied.calc_energy_current(None), Ok(3.0));
    }

    #[test]
    fn uff_worker_w01_copy_preserves_source_empty_state_and_independence() {
        let mut original_first = vec![0.0, 0.0];
        let mut original_second = vec![3.0, 4.0];
        let mut original = ForceField::new(2);
        original.positions_mut().extend([
            original_first.as_mut_slice(),
            original_second.as_mut_slice(),
        ]);
        original.fixed_points_mut().push(1);
        original
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        original.initialize().unwrap();
        assert_eq!(original.calc_energy_current(None), Ok(25.0));
        assert_eq!(original.distance_matrix, [-1.0, 5.0, -1.0]);

        super::cf3d_uff_one_kernel_counts_reset();
        let mut copied = super::cf3d_bld_b05_copy_force_field(&original);
        assert_eq!(super::cf3d_uff_one_kernel_counts(), (1, 0, 1));
        assert_eq!(copied.dimension(), 2);
        assert_eq!(copied.num_points(), 2);
        assert!(!copied.initialized);
        assert!(copied.positions().is_empty());
        assert!(copied.fixed_points().is_empty());
        assert!(copied.distance_matrix.is_empty());
        assert!(!copied.matrix_allocated);
        assert_eq!(copied.matrix_size, 0);
        assert_eq!(copied.contributions.len(), 1);
        assert!(!std::ptr::eq(
            original.contributions[0].as_ref(),
            copied.contributions[0].as_ref()
        ));

        let mut copied_first = vec![0.0, 0.0];
        let mut copied_second = vec![6.0, 8.0];
        copied
            .positions_mut()
            .extend([copied_first.as_mut_slice(), copied_second.as_mut_slice()]);
        copied.initialize().unwrap();
        assert_eq!(
            super::cf3d_bld_b05_calc_energy(&mut copied, &[0.0, 0.0, 6.0, 8.0]),
            Ok(100.0)
        );
        copied.fixed_points_mut().push(0);

        assert_eq!(original.fixed_points(), [1]);
        assert_eq!(original.positions()[0], [0.0, 0.0]);
        assert_eq!(original.positions()[1], [3.0, 4.0]);
        assert_eq!(original.distance_matrix, [-1.0, 5.0, -1.0]);
        assert_eq!(original.calc_energy_current(None), Ok(25.0));

        drop(copied);
        drop(original);
        assert_eq!(original_first, [0.0, 0.0]);
        assert_eq!(original_second, [3.0, 4.0]);
        assert_eq!(copied_first, [0.0, 0.0]);
        assert_eq!(copied_second, [6.0, 8.0]);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn uff_thread_t04_copy_outlives_source_borrows_on_scoped_thread() {
        // RDKit ForceField.cpp:159-170 copies metadata/terms but leaves
        // positions and cache empty; ForceField.h:245-251 supplies defaults.
        let mut source_first = [0.0, 0.0, 0.0];
        let mut source_second = [1.814, 0.0, 0.0];
        let source_rows_before = [
            source_first.map(f64::to_bits),
            source_second.map(f64::to_bits),
        ];

        let (copied, source_identities, source_energy) = {
            let mut original = ForceField::new(3);
            original
                .positions_mut()
                .extend([source_first.as_mut_slice(), source_second.as_mut_slice()]);
            original.fixed_points_mut().push(1);
            let _forward = attach_u06_bond_contribution(&mut original);
            let carbon = uff_sp3_carbon_params();
            let reverse =
                BondStretchContrib::new(original.positions(), 1, 0, 1.0, &carbon, &carbon)
                    .expect("valid reversed source bond contribution");
            original.add_contribution(Box::new(reverse));
            original.initialize().unwrap();

            let source_energy = original.calc_energy_current(None).unwrap();
            assert_u05_close(source_energy, 62.963262, 1.0e-3);
            let source_identities = super::cf3d_frag_accept_contribution_identities(&original);
            assert!(matches!(
                source_identities.as_slice(),
                [
                    super::Cf3dFragAcceptContributionIdentity::BondStretch {
                        end1_idx: 0,
                        end2_idx: 1,
                        ..
                    },
                    super::Cf3dFragAcceptContributionIdentity::BondStretch {
                        end1_idx: 1,
                        end2_idx: 0,
                        ..
                    }
                ]
            ));
            let source_cache = original.distance_matrix.clone();
            let source_fixed_points = original.fixed_points().to_vec();

            let copied = original.copy();
            assert_eq!(copied.dimension(), 3);
            assert_eq!(copied.num_points(), 2);
            assert!(!copied.initialized);
            assert!(copied.positions().is_empty());
            assert!(copied.fixed_points().is_empty());
            assert!(copied.distance_matrix.is_empty());
            assert!(!copied.matrix_allocated);
            assert_eq!(copied.matrix_size, 0);
            assert_eq!(
                super::cf3d_frag_accept_contribution_identities(&copied),
                source_identities
            );
            assert!(
                original
                    .contributions
                    .iter()
                    .zip(&copied.contributions)
                    .all(|(source, copy)| !std::ptr::eq(source.as_ref(), copy.as_ref()))
            );

            assert_eq!(original.fixed_points(), source_fixed_points);
            assert_eq!(original.distance_matrix, source_cache);
            assert_eq!(
                original.calc_energy_current(None).unwrap().to_bits(),
                source_energy.to_bits()
            );

            // Dropping `original` ends both old coordinate borrows before the
            // independent copy is moved to the scoped worker below.
            drop(original);
            assert_eq!(source_first.map(f64::to_bits), source_rows_before[0]);
            assert_eq!(source_second.map(f64::to_bits), source_rows_before[1]);
            (copied, source_identities, source_energy)
        };

        let expected_identities = source_identities.clone();
        let (worker_energy, worker_identities, worker_rows) = std::thread::scope(|scope| {
            scope
                .spawn(move || {
                    let mut copied = copied;
                    assert_eq!(copied.dimension(), 3);
                    assert_eq!(copied.num_points(), 2);
                    assert!(!copied.initialized);
                    assert!(copied.positions().is_empty());
                    assert!(copied.fixed_points().is_empty());
                    assert!(copied.distance_matrix.is_empty());
                    assert_eq!(
                        super::cf3d_frag_accept_contribution_identities(&copied),
                        expected_identities
                    );

                    let mut worker_first = [0.0, 0.0, 0.0];
                    let mut worker_second = [1.814, 0.0, 0.0];
                    copied
                        .positions_mut()
                        .extend([worker_first.as_mut_slice(), worker_second.as_mut_slice()]);
                    copied.initialize().unwrap();
                    assert_eq!(copied.positions().len(), 2);
                    let worker_energy = copied.calc_energy_current(None).unwrap();
                    assert_u05_close(worker_energy, 62.963262, 1.0e-3);
                    let worker_identities =
                        super::cf3d_frag_accept_contribution_identities(&copied);
                    drop(copied);
                    (
                        worker_energy,
                        worker_identities,
                        [worker_first, worker_second],
                    )
                })
                .join()
                .expect("scoped copy worker completes")
        });

        assert_eq!(worker_energy.to_bits(), source_energy.to_bits());
        assert_eq!(worker_identities, source_identities);
        assert_eq!(source_first.map(f64::to_bits), source_rows_before[0]);
        assert_eq!(source_second.map(f64::to_bits), source_rows_before[1]);
        assert_eq!(
            worker_rows.map(|row| row.map(f64::to_bits)),
            source_rows_before
        );
    }

    #[test]
    fn cf3d_f04_empty_initialize_matches_source_precondition() {
        // RDKit source: ForceField.cpp:239-250 and ForceField.cpp:404-415.
        let mut force_field = ForceField::default();

        assert_eq!(force_field.dimension(), 3);
        assert!(!force_field.initialized);
        assert_eq!(
            force_field.initialize(),
            Err(ForceFieldKernelError::NoPoints)
        );
        assert_eq!(force_field.num_points(), 0);
        assert_eq!(force_field.matrix_size, 0);
        assert!(force_field.matrix_allocated);
        assert!(force_field.distance_matrix.is_empty());
        assert!(!force_field.initialized);
        assert_eq!(
            ForceFieldKernelError::NoPoints.source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            ForceFieldKernelError::NoPoints.source_message(),
            "no points"
        );
    }

    #[test]
    fn cf3d_f04_cache_preconditions_keep_source_order_and_messages() {
        // RDKit source: ForceField.cpp:404-415.
        let mut force_field = ForceField::default();

        // The first source precondition wins when points and matrix are absent.
        assert_eq!(
            force_field.init_distance_matrix(),
            Err(ForceFieldKernelError::NoPoints)
        );

        force_field.num_points = 1;
        assert_eq!(
            force_field.init_distance_matrix(),
            Err(ForceFieldKernelError::NoDistanceMatrix)
        );
        assert_eq!(
            ForceFieldKernelError::NoDistanceMatrix.source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            ForceFieldKernelError::NoDistanceMatrix.source_message(),
            "no distance matrix"
        );

        force_field.num_points = 2;
        force_field.matrix_allocated = true;
        force_field.matrix_size = 2;
        force_field.distance_matrix = vec![17.0, 19.0];
        assert_eq!(
            force_field.init_distance_matrix(),
            Err(ForceFieldKernelError::MatrixSizeMismatch)
        );
        assert_eq!(
            ForceFieldKernelError::MatrixSizeMismatch.source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            ForceFieldKernelError::MatrixSizeMismatch.source_message(),
            "matrix size mismatch"
        );
        assert_eq!(force_field.distance_matrix, [17.0, 19.0]);

        force_field.matrix_size = 4;
        force_field.distance_matrix = vec![17.0, 19.0, 23.0, 29.0];
        force_field.init_distance_matrix().unwrap();
        assert_eq!(force_field.distance_matrix, [-1.0, -1.0, -1.0, 29.0]);
    }

    #[test]
    fn cf3d_f04_nonempty_initialization_tracks_dimension_and_point_count() {
        // RDKit source: ForceField.h:82, 244-250 and ForceField.cpp:404-415.
        let mut first = vec![1.0, 2.0];
        let mut second = vec![3.0, 5.0];
        let mut force_field = ForceField::new(2);
        force_field.fixed_points_mut().extend([1, -2]);

        force_field.positions_mut().push(first.as_mut_slice());
        force_field.initialize().unwrap();
        assert_eq!(force_field.dimension(), 2);
        assert_eq!(force_field.num_points(), 1);
        assert_eq!(force_field.matrix_size, 1);
        assert_eq!(force_field.distance_matrix, [-1.0]);
        assert!(force_field.initialized);

        force_field.distance_matrix[0] = 13.0;
        force_field.positions_mut().push(second.as_mut_slice());
        force_field.initialize().unwrap();
        assert_eq!(force_field.dimension(), 2);
        assert_eq!(force_field.num_points(), 2);
        assert_eq!(force_field.matrix_size, 3);
        assert_eq!(force_field.distance_matrix, [-1.0, -1.0, -1.0]);
        assert_eq!(force_field.fixed_points(), [1, -2]);
        assert!(force_field.initialized);

        let coordinates = [7.0, 11.0, 13.0, 17.0];
        let context = force_field.evaluation_context(&coordinates);
        assert_eq!(context.coordinates, coordinates);
        assert_eq!(context.dimension, 2);
        assert_eq!(context.num_points, 2);
        assert_eq!(context.matrix_size, 3);
        assert_eq!(context.distance_matrix, [-1.0, -1.0, -1.0]);
        context.distance_matrix[0] = 29.0;
        drop(context);

        drop(force_field.positions_mut().pop());
        force_field.initialize().unwrap();
        assert_eq!(force_field.num_points(), 1);
        assert_eq!(force_field.matrix_size, 1);
        assert_eq!(force_field.distance_matrix, [-1.0]);
        assert_eq!(force_field.fixed_points(), [1, -2]);
    }

    #[test]
    fn cf3d_f04_nondefault_dimensions_are_preserved() {
        // RDKit source: ForceField.h:82; the constructor stores its argument.
        assert_eq!(ForceField::new(2).dimension(), 2);
        assert_eq!(ForceField::new(4).dimension(), 4);
    }

    #[test]
    fn cf3d_f05_current_positions_reverse_pair_and_diagonal_use_cache() {
        // RDKit source: ForceField.cpp:172-203, including triangular cache indexing.
        let mut first = vec![1.0, 2.0, 3.0];
        let mut second = vec![4.0, 6.0, 15.0];
        let mut third = vec![-2.0, 2.0, 3.0];
        let mut force_field = ForceField::default();
        force_field.positions_mut().extend([
            first.as_mut_slice(),
            second.as_mut_slice(),
            third.as_mut_slice(),
        ]);
        force_field.initialize().unwrap();

        assert_eq!(force_field.distance(1, 0, None).unwrap(), 13.0);
        assert_eq!(force_field.distance_matrix[1], 13.0);
        assert_eq!(force_field.distance(0, 1, None).unwrap(), 13.0);
        assert_eq!(force_field.distance(0, 2, None).unwrap(), 3.0);
        assert_eq!(force_field.distance_matrix[3], 3.0);
        assert_eq!(force_field.distance(2, 2, None).unwrap(), 0.0);
        assert_eq!(force_field.distance_matrix[5], 0.0);
    }

    #[test]
    fn cf3d_f05_squared_and_const_paths_use_current_or_explicit_coordinates_without_cache() {
        // RDKit source: ForceField.cpp:206-237; distance2 and const distance are pure.
        let mut first = vec![1.0, 2.0];
        let mut second = vec![4.0, 6.0];
        let mut force_field = ForceField::new(2);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();

        let alternate = [0.0, 0.0, 0.0, 3.0];
        let other_alternate = [0.0, 0.0, 0.0, 5.0];
        assert_eq!(force_field.distance2(1, 0, None).unwrap(), 25.0);
        assert_eq!(force_field.distance_const(0, 1, None).unwrap(), 5.0);
        assert_eq!(force_field.distance2(0, 1, Some(&alternate)).unwrap(), 9.0);
        assert_eq!(
            force_field.distance_const(1, 0, Some(&alternate)).unwrap(),
            3.0
        );
        assert_eq!(force_field.distance_matrix, [-1.0, -1.0, -1.0]);

        assert_eq!(force_field.distance(0, 1, Some(&alternate)).unwrap(), 3.0);
        assert_eq!(force_field.distance_matrix[1], 3.0);
        assert_eq!(
            force_field.distance(1, 0, Some(&other_alternate)).unwrap(),
            3.0
        );
        assert_eq!(
            force_field.distance2(0, 1, Some(&other_alternate)).unwrap(),
            25.0
        );
        assert_eq!(
            force_field
                .distance_const(0, 1, Some(&other_alternate))
                .unwrap(),
            5.0
        );
        assert_eq!(force_field.distance_matrix[1], 3.0);

        force_field.distance_matrix[1] = -2.0;
        assert_eq!(force_field.distance(1, 0, None).unwrap(), 5.0);
        assert_eq!(force_field.distance_matrix[1], 5.0);

        force_field.initialize().unwrap();
        let nan_coordinates = [f64::NAN, 0.0, 0.0, 0.0];
        assert!(
            force_field
                .distance(0, 1, Some(&nan_coordinates))
                .unwrap()
                .is_nan()
        );
        assert!(force_field.distance_matrix[1].is_nan());
        assert!(
            force_field
                .distance(0, 1, Some(&alternate))
                .unwrap()
                .is_nan()
        );
        assert!(force_field.distance_matrix[1].is_nan());
    }

    #[test]
    fn cf3d_f05_four_dimensional_paths_sum_each_source_component() {
        // RDKit source: ForceField.cpp:172-203 and :206-237 loop to d_dimension.
        let mut first = vec![1.0, 2.0, 3.0, 4.0];
        let mut second = vec![2.0, 4.0, 6.0, 8.0];
        let mut force_field = ForceField::new(4);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();

        assert_eq!(force_field.distance2(0, 1, None).unwrap(), 30.0);
        assert_eq!(force_field.distance(0, 1, None).unwrap(), 30.0_f64.sqrt());
        let alternate = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 4.0];
        assert_eq!(force_field.distance2(0, 1, Some(&alternate)).unwrap(), 25.0);
        assert_eq!(
            force_field.distance_const(0, 1, Some(&alternate)).unwrap(),
            5.0
        );
        assert_eq!(force_field.distance_matrix[1], 30.0_f64.sqrt());
    }

    #[test]
    fn cf3d_f05_reinitialize_resets_cached_current_coordinate_distance() {
        // RDKit source: ForceField.cpp:172-203 and :239-250.
        let mut first = vec![1.0, 2.0];
        let mut second = vec![4.0, 6.0];
        let mut force_field = ForceField::new(2);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();

        assert_eq!(force_field.distance(0, 1, None).unwrap(), 5.0);
        assert_eq!(force_field.distance_matrix[1], 5.0);
        force_field.positions_mut()[1][0] = 7.0;
        force_field.initialize().unwrap();
        assert_eq!(force_field.distance_matrix, [-1.0, -1.0, -1.0]);
        assert_eq!(force_field.distance2(0, 1, None).unwrap(), 52.0);
        assert_eq!(force_field.distance(0, 1, None).unwrap(), 52.0_f64.sqrt());
        assert_eq!(force_field.distance_matrix[1], 52.0_f64.sqrt());
    }

    #[test]
    fn cf3d_f05_validation_order_and_source_error_context_are_preserved() {
        // RDKit source: ForceField.cpp:172-232 and RDGeneral/Invariant.h:108-150.
        let mut uninitialized = ForceField::new(2);
        let not_initialized = uninitialized.distance(0, 0, None).unwrap_err();
        assert_eq!(not_initialized.source_category(), "Pre-condition Violation");
        assert_eq!(not_initialized.source_message(), "not initialized");
        assert_eq!(
            not_initialized.source_expression().as_deref(),
            Some("df_init")
        );
        assert_eq!(
            uninitialized.distance2(0, 0, None),
            Err(ForceFieldKernelError::NotInitialized)
        );

        let mut first = vec![0.0, 0.0];
        let mut second = vec![3.0, 4.0];
        let mut force_field = ForceField::new(2);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();

        let bad_i = force_field.distance(2, 2, None).unwrap_err();
        assert_eq!(bad_i.source_category(), "Range Error");
        assert_eq!(bad_i.source_message(), "i");
        assert_eq!(bad_i.source_expression().as_deref(), Some("2 < 2"));
        let bad_j_distance = force_field.distance(1, 2, None).unwrap_err();
        assert_eq!(bad_j_distance.source_category(), "Range Error");
        assert_eq!(bad_j_distance.source_message(), "j");
        assert_eq!(bad_j_distance.source_expression().as_deref(), Some("2 < 2"));
        let bad_i_distance2 = force_field.distance2(2, 2, None).unwrap_err();
        assert_eq!(bad_i_distance2.source_category(), "Range Error");
        assert_eq!(bad_i_distance2.source_message(), "i");
        assert_eq!(
            bad_i_distance2.source_expression().as_deref(),
            Some("2 < 2")
        );
        let bad_j = force_field.distance2(1, 2, None).unwrap_err();
        assert_eq!(bad_j.source_category(), "Range Error");
        assert_eq!(bad_j.source_message(), "j");
        assert_eq!(bad_j.source_expression().as_deref(), Some("2 < 2"));
        assert_eq!(force_field.distance_const(1, 2, None), Err(bad_j));

        force_field.matrix_size = 1;
        let bad_cache_index = force_field.distance(0, 1, None).unwrap_err();
        assert_eq!(bad_cache_index.source_category(), "Invariant Violation");
        assert_eq!(bad_cache_index.source_message(), "Bad index");
        assert_eq!(
            bad_cache_index.source_expression().as_deref(),
            Some("idx < d_matSize")
        );
        assert_eq!(force_field.distance2(0, 1, None).unwrap(), 25.0);
        assert_eq!(force_field.distance_matrix[1], -1.0);
    }

    #[test]
    fn cf3d_f05_borrowed_evaluation_context_reuses_distance_semantics() {
        // RDKit source: ForceField.cpp:172-237; Contrib.h:16-34 uses its owner for these reads.
        let mut first = vec![1.0, 2.0];
        let mut second = vec![4.0, 6.0];
        let mut force_field = ForceField::new(2);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();

        let explicit_coordinates = [0.0, 0.0, 0.0, 3.0];
        let mut context = force_field.evaluation_context(&explicit_coordinates);
        assert_eq!(context.distance2(1, 0), Ok(9.0));
        assert_eq!(context.distance_const(0, 1), Ok(3.0));
        assert_eq!(context.distance_matrix, [-1.0, -1.0, -1.0]);
        assert_eq!(context.distance(0, 1), Ok(3.0));
        assert_eq!(context.distance(1, 0), Ok(3.0));
        assert_eq!(context.distance_matrix[1], 3.0);

        let mut uninitialized = ForceField::new(2);
        let mut context = uninitialized.evaluation_context(&[]);
        assert_eq!(
            context.distance(0, 0),
            Err(ForceFieldKernelError::NotInitialized)
        );
    }

    #[test]
    fn cf3d_f06_scatter_gather_preserve_point_component_order_for_dimensions_2_3_4() {
        // RDKit source: ForceField.cpp:377-404; each point's components are contiguous.
        let cases: [(u32, &[f64], &[f64], &[f64], &[f64]); 3] = [
            (2, &[1.0, -2.0], &[3.5, 4.25], &[-5.0, 6.0], &[7.25, -8.5]),
            (
                3,
                &[1.0, -2.0, 3.5],
                &[4.25, -5.0, 6.0],
                &[-7.25, 8.5, 9.0],
                &[10.5, -11.0, 12.25],
            ),
            (
                4,
                &[1.0, -2.0, 3.5, -4.25],
                &[5.0, -6.0, 7.25, -8.5],
                &[-9.0, 10.5, -11.0, 12.25],
                &[13.5, -14.0, 15.25, -16.5],
            ),
        ];

        for (dimension, first_start, second_start, first_update, second_update) in cases {
            let mut first = first_start.to_vec();
            let mut second = second_start.to_vec();
            let mut force_field = ForceField::new(dimension);
            force_field
                .positions_mut()
                .extend([first.as_mut_slice(), second.as_mut_slice()]);
            force_field.initialize().unwrap();
            force_field.distance_matrix[0] = 41.5;
            let cache_before = force_field.distance_matrix.clone();

            let transfer_len = (dimension * 2) as usize;
            let mut flat = vec![-91.25; transfer_len + 1];
            assert_eq!(force_field.scatter(&mut flat), Ok(()));
            let expected_scatter: Vec<_> = first_start
                .iter()
                .chain(second_start.iter())
                .copied()
                .collect();
            assert_eq!(&flat[..transfer_len], expected_scatter);
            assert_eq!(flat[transfer_len], -91.25);
            assert_eq!(force_field.distance_matrix, cache_before);

            let mut gathered: Vec<_> = first_update
                .iter()
                .chain(second_update.iter())
                .copied()
                .collect();
            gathered.push(123.75);
            force_field.gather(&gathered).unwrap();
            assert_eq!(&force_field.positions()[0][..], first_update);
            assert_eq!(&force_field.positions()[1][..], second_update);
            assert_eq!(force_field.distance_matrix, cache_before);
        }
    }

    #[test]
    fn cf3d_f06_empty_point_transfer_accepts_empty_and_ignores_unused_flat_tail() {
        // RDKit source: ForceField.cpp:377-404. These calls isolate the source
        // zero-iteration transfer path; initialize() itself rejects zero points.
        let mut force_field = ForceField::new(4);
        force_field.initialized = true;
        let mut empty: [f64; 0] = [];

        assert_eq!(force_field.scatter(&mut empty), Ok(()));
        assert_eq!(force_field.gather(&empty), Ok(()));

        let mut unused_output = [23.0, -17.0];
        assert_eq!(force_field.scatter(&mut unused_output), Ok(()));
        assert_eq!(unused_output, [23.0, -17.0]);
        assert_eq!(force_field.gather(&unused_output), Ok(()));
    }

    #[test]
    fn cf3d_f06_transfer_preconditions_and_postcondition_metadata_match_source() {
        // RDKit source: ForceField.cpp:377-404 and Invariant.h:108-120.
        let mut force_field = ForceField::new(3);
        let mut empty: [f64; 0] = [];
        let scatter_error = force_field.scatter(&mut empty).unwrap_err();
        assert_eq!(scatter_error, ForceFieldKernelError::NotInitialized);
        assert_eq!(scatter_error.source_category(), "Pre-condition Violation");
        assert_eq!(scatter_error.source_message(), "not initialized");
        assert_eq!(
            scatter_error.source_expression().as_deref(),
            Some("df_init")
        );
        assert_eq!(
            force_field.gather(&empty),
            Err(ForceFieldKernelError::NotInitialized)
        );

        let postcondition = ForceFieldKernelError::TransferPostcondition;
        assert_eq!(postcondition.source_category(), "Post-condition Violation");
        assert_eq!(postcondition.source_message(), "bad index");
        assert_eq!(
            postcondition.source_expression().as_deref(),
            Some("tab == this->dimension() * d_positions.size()")
        );
    }

    #[test]
    fn cf3d_bld_b04_vec_rows_alias_through_scatter_and_gather() {
        // RDKit source: ForceField.cpp:377-404; source order is point, then component.
        let mut first = vec![1.0, -2.0, 3.5];
        let mut second = vec![4.25, -5.0, 6.0];
        let mut force_field = ForceField::new(3);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();

        let mut scattered = [-91.25; 7];
        assert_eq!(force_field.scatter(&mut scattered), Ok(()));
        assert_eq!(scattered, [1.0, -2.0, 3.5, 4.25, -5.0, 6.0, -91.25]);

        let gathered = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0];
        assert_eq!(force_field.gather(&gathered), Ok(()));
        assert_eq!(&force_field.positions()[0][..], &[10.0, 11.0, 12.0]);
        assert_eq!(&force_field.positions()[1][..], &[13.0, 14.0, 15.0]);
        drop(force_field);
        assert_eq!(first, vec![10.0, 11.0, 12.0]);
        assert_eq!(second, vec![13.0, 14.0, 15.0]);
    }

    #[test]
    fn cf3d_bld_b04_fixed_array_rows_are_borrowed_without_copying() {
        // RDKit's positions() stores point references; [f64; 3] rows are lent directly.
        let mut first = [1.0, 2.0, 3.0];
        let mut second = [4.0, 5.0, 6.0];
        let mut force_field = ForceField::new(3);
        force_field.positions_mut().push(first.as_mut_slice());
        force_field.positions_mut().push(second.as_mut_slice());
        force_field.initialize().unwrap();

        let mut scattered = [0.0; 6];
        assert_eq!(force_field.scatter(&mut scattered), Ok(()));
        assert_eq!(scattered, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(
            force_field.gather(&[7.0, 8.0, 9.0, 10.0, 11.0, 12.0]),
            Ok(())
        );
        assert_eq!(&force_field.positions()[0][..], &[7.0, 8.0, 9.0]);
        assert_eq!(&force_field.positions()[1][..], &[10.0, 11.0, 12.0]);
        drop(force_field);
        assert_eq!(first, [7.0, 8.0, 9.0]);
        assert_eq!(second, [10.0, 11.0, 12.0]);
    }

    #[test]
    fn cf3d_bld_b04_empty_and_malformed_dimensions_keep_source_indexing() {
        // ForceField::initialize checks the empty point count; scatter/gather
        // otherwise visit exactly dimension components and ignore row tails.
        let mut empty = ForceField::new(3);
        assert_eq!(empty.initialize(), Err(ForceFieldKernelError::NoPoints));
        assert_eq!(
            empty.scatter(&mut []),
            Err(ForceFieldKernelError::NotInitialized)
        );

        let mut longer_row = [1.0, 2.0, 99.0];
        let mut two_dimensional = ForceField::new(2);
        two_dimensional
            .positions_mut()
            .push(longer_row.as_mut_slice());
        two_dimensional.initialize().unwrap();
        let mut flat = [-7.0, -8.0, -9.0];
        assert_eq!(two_dimensional.scatter(&mut flat), Ok(()));
        assert_eq!(flat, [1.0, 2.0, -9.0]);
        assert_eq!(two_dimensional.gather(&[4.0, 5.0]), Ok(()));
        assert_eq!(&two_dimensional.positions()[0][..], &[4.0, 5.0, 99.0]);
        drop(two_dimensional);
        assert_eq!(longer_row, [4.0, 5.0, 99.0]);

        // Source initialize() derives only the row count; it does not inspect
        // row width. Do not dereference this short row: C++ has no defined
        // behavior for the corresponding out-of-range Point3D component.
        let mut short_row = [3.0, 4.0];
        let mut three_dimensional = ForceField::new(3);
        three_dimensional
            .positions_mut()
            .push(short_row.as_mut_slice());
        assert_eq!(three_dimensional.initialize(), Ok(()));
        assert_eq!(three_dimensional.num_points(), 1);
    }

    #[test]
    fn uff_one_u01_empty_contributions_initialize_and_have_zero_energy() {
        // Pinned ForceField.cpp::initialize accepts any nonzero point count,
        // initializes the triangular cache to -1.0, and does not require terms.
        // calcEnergy(vector*) returns 0.0 after its initialized precondition
        // when the contribution list is empty.
        let mut position = [1.25, -2.5, 0.5];
        let mut force_field = ForceField::new(3);
        force_field.positions_mut().push(&mut position);

        assert_eq!(force_field.initialize(), Ok(()));
        assert!(force_field.initialized);
        assert!(force_field.contributions.is_empty());
        assert_eq!(force_field.distance_matrix, [-1.0]);
        assert_eq!(force_field.calc_energy_current(None), Ok(0.0));
    }

    #[test]
    fn uff_one_u01_empty_positions_preserve_typed_initialization_failure() {
        // initDistanceMatrix's first source precondition is the only ordinary
        // input-reachable initialize failure: initialization derives a nonnull
        // cache and a matching triangular size itself.
        let mut force_field = ForceField::new(3);

        assert_eq!(
            force_field.initialize(),
            Err(ForceFieldKernelError::NoPoints)
        );
        assert!(!force_field.initialized);
        assert_eq!(
            ForceFieldKernelError::NoPoints.source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            ForceFieldKernelError::NoPoints.source_message(),
            "no points"
        );
        assert_eq!(
            force_field.calc_energy_current(None),
            Err(ForceFieldKernelError::NotInitialized)
        );
    }

    #[test]
    fn uff_one_u03_current_energy_sums_fixed_contributions_in_order() {
        let mut first = [0.0];
        let mut second = [2.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 1.0,
            }));
        force_field
            .contributions
            .push(Box::new(SquaredDistanceContribution {
                first: 0,
                second: 1,
                scale: 0.5,
            }));

        let mut contribution_energies = vec![-11.0];
        assert_eq!(
            force_field.calc_energy_current(Some(&mut contribution_energies)),
            Ok(6.0)
        );
        assert_eq!(contribution_energies, [4.0, 2.0]);
        assert_eq!(first, [0.0]);
        assert_eq!(second, [2.0]);
    }

    #[test]
    fn uff_one_u03_empty_contributions_keep_source_output_and_precondition_order() {
        let mut contribution_energies = vec![17.5];
        let mut uninitialized = ForceField::new(3);
        assert_eq!(
            uninitialized.calc_energy_current(Some(&mut contribution_energies)),
            Err(ForceFieldKernelError::NotInitialized)
        );
        assert_eq!(contribution_energies, [17.5]);

        let mut position = [1.0, -2.0, 3.0];
        let mut initialized = ForceField::new(3);
        initialized.positions_mut().push(&mut position);
        initialized.initialize().unwrap();
        assert!(initialized.contributions.is_empty());
        assert_eq!(
            initialized.calc_energy_current(Some(&mut contribution_energies)),
            Ok(0.0)
        );
        // ForceField.cpp returns before clearing `contribs` when there are no
        // terms; the preexisting value is therefore retained.
        assert_eq!(contribution_energies, [17.5]);
    }

    #[test]
    fn uff_one_u03_current_energy_propagates_typed_contribution_error() {
        let mut first = [0.0];
        let mut second = [2.0];
        let mut force_field = ForceField::new(1);
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
            .contributions
            .push(Box::new(FailOnEnergyCallContribution {
                calls: Cell::new(0),
                fail_on: 1,
            }));
        let mut contribution_energies = vec![99.0];

        assert_eq!(
            force_field.calc_energy_current(Some(&mut contribution_energies)),
            Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: 2,
                upper_bound: 2,
            })
        );
        assert!(contribution_energies.is_empty());
        assert_eq!(first, [0.0]);
        assert_eq!(second, [2.0]);
    }

    fn cf3d_bld_b04_field_borrowing_array_row<'a>(row: &'a mut [f64; 3]) -> ForceField<'a> {
        let mut force_field = ForceField::new(3);
        let row_slice: &'a mut [f64] = row;
        force_field.positions_mut().push(row_slice);
        force_field
    }

    #[test]
    fn cf3d_bld_b04_mutable_alias_ends_when_field_is_dropped() {
        let mut coordinate = [2.0, 3.0, 4.0];
        {
            let mut force_field = cf3d_bld_b04_field_borrowing_array_row(&mut coordinate);
            force_field.positions_mut()[0][0] = 17.0;
            assert_eq!(&force_field.positions()[0][..], &[17.0, 3.0, 4.0]);
        }
        assert_eq!(coordinate, [17.0, 3.0, 4.0]);
    }

    #[test]
    fn cf3d_bld_b04_copy_clones_terms_without_copying_position_rows() {
        let mut first = [0.0, 0.0, 0.0];
        let mut second = [2.0, 0.0, 0.0];
        let mut force_field = ForceField::new(3);
        force_field.positions_mut().push(first.as_mut_slice());
        force_field.positions_mut().push(second.as_mut_slice());
        let carbon = uff_sp3_carbon_params();
        let bond = BondStretchContrib::new(force_field.positions(), 0, 1, 1.0, &carbon, &carbon)
            .expect("fixed source indices form a valid bond term");
        force_field.add_contribution(Box::new(bond));
        force_field.initialize().unwrap();

        let copied = force_field.copy();
        assert_eq!(copied.dimension(), 3);
        assert_eq!(copied.num_points(), 2);
        assert!(!copied.initialized);
        assert!(copied.positions().is_empty());
        assert_eq!(copied.contributions.len(), 1);

        force_field.positions_mut()[0][0] = 1.5;
        let mut scattered = [0.0; 6];
        assert_eq!(force_field.scatter(&mut scattered), Ok(()));
        assert_eq!(scattered, [1.5, 0.0, 0.0, 2.0, 0.0, 0.0]);
        drop(copied);
        drop(force_field);
        assert_eq!(first, [1.5, 0.0, 0.0]);
        assert_eq!(second, [2.0, 0.0, 0.0]);
    }

    fn cf3d_f14_force_field<'a>(
        first: &'a mut Vec<f64>,
        second: &'a mut Vec<f64>,
    ) -> ForceField<'a> {
        let mut force_field = ForceField::default();
        force_field
            .positions_mut()
            .extend([first.as_mut_slice(), second.as_mut_slice()]);
        force_field.initialize().unwrap();
        force_field
    }

    fn cf3d_f14_energy(
        force_field: &mut ForceField<'_>,
        contribution: &DistanceConstraintContrib,
        coordinates: &[f64],
    ) -> Result<f64, ForceFieldKernelError> {
        // RDKit source: ForceField.cpp:310-325 calls initDistanceMatrix before
        // evaluating each explicit trial-position contribution.
        force_field.init_distance_matrix()?;
        let mut context = force_field.evaluation_context(coordinates);
        contribution.get_energy(&mut context)
    }

    #[test]
    fn cf3d_f14_constructors_keep_absolute_relative_offsets_and_source_clamps() {
        // RDKit source: DistanceConstraint.cpp:16-53,
        // Geometry/point.cpp:64-70, and Geometry/point.h:158-161.
        let mut first = vec![0.0, 0.0, 0.0];
        let mut second = vec![3.0, 4.0, 0.0];
        let force_field = cf3d_f14_force_field(&mut first, &mut second);

        let absolute =
            DistanceConstraintContrib::new_absolute(&force_field, 0, 1, 2.0, 6.0, 4.0).unwrap();
        assert_eq!(absolute.min_len(), 2.0);
        assert_eq!(absolute.max_len(), 6.0);

        let relative =
            DistanceConstraintContrib::new_relative(&force_field, 0, 1, true, 2.0, 4.0, 4.0)
                .unwrap();
        assert_eq!(relative.min_len(), 7.0);
        assert_eq!(relative.max_len(), 9.0);

        let independently_clamped =
            DistanceConstraintContrib::new_relative(&force_field, 0, 1, true, -10.0, -3.0, 4.0)
                .unwrap();
        assert_eq!(independently_clamped.min_len(), 0.0);
        assert_eq!(independently_clamped.max_len(), 2.0);

        let not_relative =
            DistanceConstraintContrib::new_relative(&force_field, 0, 1, false, 1.0, 2.0, 4.0)
                .unwrap();
        assert_eq!(not_relative.min_len(), 1.0);
        assert_eq!(not_relative.max_len(), 2.0);
    }

    #[test]
    fn cf3d_f14_relative_std_max_preserves_nan_from_initial_point_distance() {
        // RDKit source: DistanceConstraint.cpp:40-45 and Point3D subtraction/
        // length. std::max(NaN, 0.0) returns its first operand by comparison.
        let mut first = vec![f64::NAN, 0.0, 0.0];
        let mut second = vec![0.0, 0.0, 0.0];
        let force_field = cf3d_f14_force_field(&mut first, &mut second);

        let relative =
            DistanceConstraintContrib::new_relative(&force_field, 0, 1, true, 1.0, 2.0, 4.0)
                .unwrap();
        assert!(relative.min_len().is_nan());
        assert!(relative.max_len().is_nan());
    }

    #[test]
    fn cf3d_f14_constructor_errors_keep_endpoint_and_bound_order() {
        // RDKit source: DistanceConstraint.cpp:19-22 and 35-39.
        let mut first = vec![0.0, 0.0, 0.0];
        let mut second = vec![3.0, 4.0, 0.0];
        let force_field = cf3d_f14_force_field(&mut first, &mut second);

        assert_eq!(
            DistanceConstraintContrib::new_absolute(&force_field, 2, 2, 5.0, 1.0, 4.0).unwrap_err(),
            ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: 2,
                upper_bound: 2,
            }
        );
        assert_eq!(
            DistanceConstraintContrib::new_absolute(&force_field, 0, 2, 5.0, 1.0, 4.0).unwrap_err(),
            ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::J,
                index: 2,
                upper_bound: 2,
            }
        );
        assert_eq!(
            DistanceConstraintContrib::new_relative(&force_field, 2, 2, false, 5.0, 1.0, 4.0)
                .unwrap_err(),
            ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: 2,
                upper_bound: 2,
            }
        );
        assert_eq!(
            DistanceConstraintContrib::new_relative(&force_field, 0, 2, false, 5.0, 1.0, 4.0)
                .unwrap_err(),
            ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::J,
                index: 2,
                upper_bound: 2,
            }
        );
        assert_eq!(
            DistanceConstraintContrib::new_absolute(&force_field, 0, 1, 5.0, 1.0, 4.0).unwrap_err(),
            ForceFieldKernelError::BadBounds
        );
        assert_eq!(
            DistanceConstraintContrib::new_relative(&force_field, 0, 1, true, 1.0, f64::NAN, 4.0)
                .unwrap_err(),
            ForceFieldKernelError::BadBounds
        );
        assert_eq!(
            ForceFieldKernelError::BadBounds.source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            ForceFieldKernelError::BadBounds.source_message(),
            "bad bounds"
        );
        assert_eq!(
            ForceFieldKernelError::BadBounds
                .source_expression()
                .as_deref(),
            Some("maxLen >= minLen")
        );
    }

    #[test]
    fn cf3d_f14_energy_keeps_strict_interval_branches_and_degenerate_bounds() {
        // RDKit source: DistanceConstraint.cpp:55-69 and ForceField.cpp:172-203.
        let mut first = vec![0.0, 0.0, 0.0];
        let mut second = vec![1.0, 0.0, 0.0];
        let mut force_field = cf3d_f14_force_field(&mut first, &mut second);
        let contribution =
            DistanceConstraintContrib::new_absolute(&force_field, 0, 1, 2.0, 6.0, 4.0).unwrap();

        assert_eq!(
            cf3d_f14_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            ),
            Ok(2.0)
        );
        assert_eq!(
            cf3d_f14_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, 2.0, 0.0, 0.0]
            ),
            Ok(0.0)
        );
        assert_eq!(
            cf3d_f14_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, 4.0, 0.0, 0.0]
            ),
            Ok(0.0)
        );
        assert_eq!(
            cf3d_f14_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, 6.0, 0.0, 0.0]
            ),
            Ok(0.0)
        );
        assert_eq!(
            cf3d_f14_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, 9.0, 0.0, 0.0]
            ),
            Ok(18.0)
        );
        assert_eq!(
            cf3d_f14_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, f64::NAN, 0.0, 0.0]
            ),
            Ok(0.0)
        );

        let point =
            DistanceConstraintContrib::new_absolute(&force_field, 0, 1, 3.0, 3.0, 4.0).unwrap();
        assert_eq!(
            cf3d_f14_energy(&mut force_field, &point, &[0.0, 0.0, 0.0, 3.0, 0.0, 0.0]),
            Ok(0.0)
        );
        assert_eq!(
            cf3d_f14_energy(&mut force_field, &point, &[0.0, 0.0, 0.0, 2.0, 0.0, 0.0]),
            Ok(2.0)
        );
        assert_eq!(
            cf3d_f14_energy(&mut force_field, &point, &[0.0, 0.0, 0.0, 4.0, 0.0, 0.0]),
            Ok(2.0)
        );
    }

    fn cf3d_f15_add_constraint(
        force_field: &mut ForceField<'_>,
        first: u32,
        second: u32,
        min_len: f64,
        max_len: f64,
        force_constant: f64,
    ) {
        let contribution = DistanceConstraintContrib::new_absolute(
            force_field,
            first,
            second,
            min_len,
            max_len,
            force_constant,
        )
        .unwrap();
        force_field.contributions.push(Box::new(contribution));
    }

    #[test]
    fn cf3d_f15_gradient_violations_write_additive_equal_and_opposite_terms() {
        // RDKit source: DistanceConstraint.cpp:78-94. Lower-side writes add
        // to both existing endpoint values, with opposite signs.
        let mut first = vec![0.0, 0.0, 0.0];
        let mut second = vec![1.0, 0.0, 0.0];
        let mut lower = cf3d_f14_force_field(&mut first, &mut second);
        cf3d_f15_add_constraint(&mut lower, 0, 1, 2.0, 6.0, 4.0);
        let mut lower_gradient = [1.0; 6];
        lower
            .calc_grad(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], &mut lower_gradient)
            .unwrap();
        assert_eq!(lower_gradient, [5.0, 1.0, 1.0, -3.0, 1.0, 1.0]);
        assert_eq!(lower_gradient[0] - 1.0, -(lower_gradient[3] - 1.0));

        // The upper violation uses dist-maxLen and preserves endpoint order.
        let mut first = vec![0.0, 0.0, 0.0];
        let mut second = vec![9.0, 0.0, 0.0];
        let mut upper = cf3d_f14_force_field(&mut first, &mut second);
        cf3d_f15_add_constraint(&mut upper, 0, 1, 2.0, 6.0, 4.0);
        let mut upper_gradient = [0.0; 6];
        upper
            .calc_grad(&[0.0, 0.0, 0.0, 9.0, 0.0, 0.0], &mut upper_gradient)
            .unwrap();
        assert_eq!(upper_gradient, [-12.0, 0.0, 0.0, 12.0, 0.0, 0.0]);
        assert_eq!(upper_gradient[0], -upper_gradient[3]);

        // Reversing the stored endpoints reverses which point receives each
        // additive term while retaining equal and opposite values.
        let mut first = vec![0.0, 0.0, 0.0];
        let mut second = vec![9.0, 0.0, 0.0];
        let mut reversed = cf3d_f14_force_field(&mut first, &mut second);
        cf3d_f15_add_constraint(&mut reversed, 1, 0, 2.0, 6.0, 4.0);
        let mut reversed_gradient = [0.0; 6];
        reversed
            .calc_grad(&[0.0, 0.0, 0.0, 9.0, 0.0, 0.0], &mut reversed_gradient)
            .unwrap();
        assert_eq!(reversed_gradient, [-12.0, 0.0, 0.0, 12.0, 0.0, 0.0]);
    }

    #[test]
    fn cf3d_f15_gradient_inside_and_both_endpoints_preserve_the_input() {
        // RDKit source: DistanceConstraint.cpp:78-85 returns before writes
        // unless one of the strict violation comparisons succeeds.
        for distance in [2.0, 4.0, 6.0] {
            let mut first = vec![0.0, 0.0, 0.0];
            let mut second = vec![distance, 0.0, 0.0];
            let mut force_field = cf3d_f14_force_field(&mut first, &mut second);
            cf3d_f15_add_constraint(&mut force_field, 0, 1, 2.0, 6.0, 4.0);
            let coordinates = [0.0, 0.0, 0.0, distance, 0.0, 0.0];
            let mut gradient = [1.0; 6];
            force_field.calc_grad(&coordinates, &mut gradient).unwrap();
            assert_eq!(gradient, [1.0; 6], "distance={distance}");
        }
    }

    #[test]
    fn cf3d_f15_gradient_zero_distance_uses_the_source_floor() {
        // RDKit source: DistanceConstraint.cpp:87-92 uses max(dist, 1e-8).
        // At zero distance the source still writes zero, not 0/0 NaNs.
        let mut first = vec![0.0, 0.0, 0.0];
        let mut second = vec![0.0, 0.0, 0.0];
        let mut force_field = cf3d_f14_force_field(&mut first, &mut second);
        cf3d_f15_add_constraint(&mut force_field, 0, 1, 1.0, 2.0, 4.0);
        let mut gradient = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

        force_field.calc_grad(&[0.0; 6], &mut gradient).unwrap();
        assert_eq!(gradient, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }

    #[test]
    fn cf3d_f15_forcefield_dispatch_preserves_the_shared_distance_cache() {
        // RDKit source: ForceField.cpp:353-375 does not reset the cache before
        // calcGrad(pos); DistanceConstraint::getGrad uses ForceField::distance.
        // This distinguishes the source path from the legacy distance2 call.
        let mut first = vec![0.0, 0.0, 0.0];
        let mut second = vec![5.0, 0.0, 0.0];
        let mut force_field = cf3d_f14_force_field(&mut first, &mut second);
        cf3d_f15_add_constraint(&mut force_field, 0, 1, 2.0, 6.0, 4.0);
        let mut gradient = [1.0; 6];

        force_field
            .calc_grad(&[0.0, 0.0, 0.0, 5.0, 0.0, 0.0], &mut gradient)
            .unwrap();
        assert_eq!(force_field.distance_matrix[1], 5.0);
        force_field
            .calc_grad(&[0.0, 0.0, 0.0, 9.0, 0.0, 0.0], &mut gradient)
            .unwrap();
        assert_eq!(force_field.distance_matrix[1], 5.0);
        assert_eq!(gradient, [1.0; 6]);
    }

    fn cf3d_f16_force_field<'a>(positions: &'a mut [Vec<f64>]) -> ForceField<'a> {
        let mut force_field = ForceField::default();
        force_field
            .positions_mut()
            .extend(positions.iter_mut().map(Vec::as_mut_slice));
        force_field.initialize().unwrap();
        force_field
    }

    fn cf3d_f16_energy(
        force_field: &mut ForceField<'_>,
        contribution: &DistanceConstraintContribs,
        coordinates: &[f64],
    ) -> Result<f64, ForceFieldKernelError> {
        let context = force_field.evaluation_context(coordinates);
        contribution.get_energy(&context)
    }

    #[test]
    fn cf3d_f16_packed_empty_size_and_add_ordered_errors() {
        // RDKit source: DistanceConstraints.cpp:22-31 validates idx1, idx2,
        // then bounds, and appends only after every check succeeds.
        let mut positions = vec![vec![0.0; 3], vec![1.0, 0.0, 0.0]];
        let mut force_field = cf3d_f16_force_field(&mut positions);
        let mut contribution = DistanceConstraintContribs::new(&force_field);

        assert!(contribution.empty());
        assert_eq!(contribution.size(), 0);
        assert_eq!(
            cf3d_f16_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            ),
            Ok(0.0)
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 2, 3, 2.0, 1.0, 4.0),
            Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 0, 2, 2.0, 1.0, 4.0),
            Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::J,
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 0, 1, 2.0, 1.0, 4.0),
            Err(ForceFieldKernelError::BadBounds)
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 0, 1, 1.0, f64::NAN, 4.0),
            Err(ForceFieldKernelError::BadBounds)
        );
        assert!(contribution.empty());
        assert_eq!(contribution.size(), 0);
    }

    #[test]
    fn cf3d_f16_packed_energy_keeps_multiple_repeated_pairs_and_intervals() {
        // RDKit source: DistanceConstraints.cpp:22-31, 33-52, 54-72.
        // Terms remain distinct and accumulate in insertion order; lower,
        // upper, strict-endpoint and inside branches use independent bounds.
        let mut positions = vec![
            vec![0.0, 0.0, 0.0],
            vec![2.0, 0.0, 0.0],
            vec![0.0, 3.0, 0.0],
        ];
        let mut force_field = cf3d_f16_force_field(&mut positions);
        let mut contribution = DistanceConstraintContribs::new(&force_field);
        contribution
            .add_contrib(&force_field, 0, 1, 1.0, 3.0, 5.0)
            .unwrap();
        contribution
            .add_contrib_relative(&force_field, 0, 2, true, 0.0, 1.0, 2.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 0, 1, 0.0, 1.0, 1.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 1, 2, 4.0, 6.0, 4.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 0, 2, 7.0, 9.0, 10.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 0, 1, 5.0, 5.0, 9.0)
            .unwrap();
        assert!(!contribution.empty());
        assert_eq!(contribution.size(), 6);

        let coordinates = [
            0.0, 0.0, 0.0, // point 0
            5.0, 0.0, 0.0, // point 1: distance(0,1)=5
            8.0, 0.0, 0.0, // point 2: distance(0,2)=8, distance(1,2)=3
        ];
        // Energies in insertion order: 10 + 16 + 8 + 2 + 0 + 0.
        assert_eq!(
            cf3d_f16_energy(&mut force_field, &contribution, &coordinates),
            Ok(36.0)
        );
    }

    #[test]
    fn cf3d_f16_relative_insertion_preserves_clamps_false_mode_and_nan() {
        // RDKit source: DistanceConstraints.cpp:33-52 and Point3D
        // subtraction/length. Relative=true offsets each bound independently;
        // relative=false leaves both supplied values unchanged.
        let mut positions = vec![vec![0.0, 0.0, 0.0], vec![3.0, 4.0, 0.0]];
        let mut force_field = cf3d_f16_force_field(&mut positions);
        let mut contribution = DistanceConstraintContribs::new(&force_field);
        contribution
            .add_contrib_relative(&force_field, 0, 1, true, -10.0, -3.0, 2.0)
            .unwrap(); // [0, 2]
        contribution
            .add_contrib_relative(&force_field, 0, 1, true, -10.0, 1.0, 2.0)
            .unwrap(); // [0, 6]
        contribution
            .add_contrib_relative(&force_field, 0, 1, false, -10.0, -3.0, 2.0)
            .unwrap(); // [-10, -3]
        assert_eq!(
            cf3d_f16_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, 5.0, 0.0, 0.0]
            ),
            Ok(234.0)
        );

        let mut nan_positions = vec![vec![0.0, 0.0, 0.0], vec![f64::NAN, 0.0, 0.0]];
        let mut nan_force_field = cf3d_f16_force_field(&mut nan_positions);
        let mut nan_contribution = DistanceConstraintContribs::new(&nan_force_field);
        nan_contribution
            .add_contrib_relative(&nan_force_field, 0, 1, true, 1.0, 2.0, 1.0)
            .unwrap();
        assert_eq!(
            cf3d_f16_energy(
                &mut nan_force_field,
                &nan_contribution,
                &[0.0, 0.0, 0.0, 5.0, 0.0, 0.0]
            ),
            Ok(0.0)
        );
    }

    #[test]
    fn cf3d_f16_energy_keeps_squared_bounds_and_does_not_fill_distance_cache() {
        // RDKit source: DistanceConstraints.cpp:54-72 uses distance2 and
        // squared thresholds, unlike DistanceConstraint.cpp's distance path.
        let mut positions = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]];
        let mut force_field = cf3d_f16_force_field(&mut positions);
        let mut contribution = DistanceConstraintContribs::new(&force_field);
        contribution
            .add_contrib(&force_field, 0, 1, -5.0, -3.0, 1.0)
            .unwrap();
        assert_eq!(force_field.distance_matrix[1], -1.0);
        assert_eq!(
            cf3d_f16_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, 2.0, 0.0, 0.0]
            ),
            Ok(24.5)
        );
        assert_eq!(force_field.distance_matrix[1], -1.0);
        assert_eq!(
            cf3d_f16_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, f64::NAN, 0.0, 0.0]
            ),
            Ok(0.0)
        );
        assert_eq!(force_field.distance_matrix[1], -1.0);
    }

    #[test]
    fn cf3d_f16_energy_accumulates_repeated_terms_in_source_order() {
        // RDKit source: DistanceConstraints.cpp:59-71 accumulates each vector
        // element directly; repeated pairs are not deduplicated or reordered.
        let mut positions = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]];
        let mut force_field = cf3d_f16_force_field(&mut positions);
        let mut contribution = DistanceConstraintContribs::new(&force_field);
        contribution
            .add_contrib(&force_field, 0, 1, 1.0, 1.0, 2.0e16)
            .unwrap();
        contribution
            .add_contrib(&force_field, 0, 1, 1.0, 1.0, 2.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 0, 1, 1.0, 1.0, -2.0e16)
            .unwrap();
        assert_eq!(contribution.size(), 3);
        // At distance zero these terms are 1e16, 1, -1e16 in source order.
        assert_eq!(
            cf3d_f16_energy(
                &mut force_field,
                &contribution,
                &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            ),
            Ok(0.0)
        );
    }

    fn cf3d_f17_add_packed_constraint(
        force_field: &mut ForceField<'_>,
        idx1: u32,
        idx2: u32,
        min_len: f64,
        max_len: f64,
        force_constant: f64,
    ) {
        let mut contribution = DistanceConstraintContribs::new(force_field);
        contribution
            .add_contrib(force_field, idx1, idx2, min_len, max_len, force_constant)
            .unwrap();
        force_field.contributions.push(Box::new(contribution));
    }

    #[test]
    fn cf3d_f17_packed_gradient_matches_fixed_single_term_source_values() {
        // RDKit source: DistanceConstraints.cpp:71-96 uses strict squared
        // lower/upper comparisons and accumulates additive endpoint updates.
        // Fixed source-derived results use p0=(0,0,0), p1=(4,0,0), k=2.
        let coordinates = [0.0, 0.0, 0.0, 4.0, 0.0, 0.0];
        let cases = [
            (5.0, 7.0, [3.0, 1.0, 1.0, -1.0, 1.0, 1.0]), // lower violation
            (1.0, 3.0, [-1.0, 1.0, 1.0, 3.0, 1.0, 1.0]), // upper violation
            (3.0, 5.0, [1.0; 6]),                        // interval interior
            (4.0, 4.0, [1.0; 6]),                        // both strict endpoints
        ];

        for (min_len, max_len, expected) in cases {
            let mut packed_positions = vec![vec![0.0; 3], vec![4.0, 0.0, 0.0]];
            let mut packed_field = cf3d_f16_force_field(&mut packed_positions);
            cf3d_f17_add_packed_constraint(&mut packed_field, 0, 1, min_len, max_len, 2.0);
            let mut packed_gradient = [1.0; 6];
            packed_field
                .calc_grad(&coordinates, &mut packed_gradient)
                .unwrap();
            assert_eq!(packed_gradient, expected, "packed [{min_len}, {max_len}]");
            assert_eq!(packed_field.distance_matrix[1], -1.0);

            let mut single_positions = vec![vec![0.0; 3], vec![4.0, 0.0, 0.0]];
            let mut single_field = cf3d_f16_force_field(&mut single_positions);
            let single =
                DistanceConstraintContrib::new_absolute(&single_field, 0, 1, min_len, max_len, 2.0)
                    .unwrap();
            single_field.contributions.push(Box::new(single));
            let mut single_gradient = [1.0; 6];
            single_field
                .calc_grad(&coordinates, &mut single_gradient)
                .unwrap();
            assert_eq!(single_gradient, expected, "single [{min_len}, {max_len}]");
            assert_eq!(packed_gradient, single_gradient);
        }
    }

    #[test]
    fn cf3d_f17_overlapping_packed_pairs_accumulate_source_endpoint_updates() {
        // RDKit source: DistanceConstraints.cpp:75-96 visits each packed term
        // in insertion order and adds both endpoint updates to the same buffer.
        let mut positions = vec![
            vec![0.0, 0.0, 0.0],
            vec![4.0, 0.0, 0.0],
            vec![4.0, 3.0, 0.0],
        ];
        let mut force_field = cf3d_f16_force_field(&mut positions);
        let mut contribution = DistanceConstraintContribs::new(&force_field);
        contribution
            .add_contrib(&force_field, 0, 1, 5.0, 9.0, 2.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 1, 2, 0.0, 2.0, 3.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 0, 2, 0.0, 4.0, 5.0)
            .unwrap();
        force_field.contributions.push(Box::new(contribution));

        let coordinates = [0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 4.0, 3.0, 0.0];
        let mut gradient = [0.0; 9];
        force_field.calc_grad(&coordinates, &mut gradient).unwrap();
        assert_eq!(gradient, [-2.0, -3.0, 0.0, -2.0, -3.0, 0.0, 4.0, 6.0, 0.0]);
        assert_eq!(force_field.distance_matrix[1], -1.0);
        assert_eq!(force_field.distance_matrix[3], -1.0);
        assert_eq!(force_field.distance_matrix[4], -1.0);
    }

    #[test]
    fn cf3d_f17_repeated_gradient_terms_keep_source_accumulation_order() {
        // RDKit source: DistanceConstraints.cpp:76-96 applies each dGrad
        // immediately; 1e16 + 1 - 1e16 therefore rounds to zero in this order.
        let mut positions = vec![vec![0.0; 3], vec![1.0, 0.0, 0.0]];
        let mut force_field = cf3d_f16_force_field(&mut positions);
        let mut contribution = DistanceConstraintContribs::new(&force_field);
        for force_constant in [1.0e16, 1.0, -1.0e16] {
            contribution
                .add_contrib(&force_field, 0, 1, 2.0, 2.0, force_constant)
                .unwrap();
        }
        force_field.contributions.push(Box::new(contribution));

        let coordinates = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let mut gradient = [0.0; 6];
        force_field.calc_grad(&coordinates, &mut gradient).unwrap();
        assert_eq!(gradient, [0.0; 6]);
    }

    #[test]
    fn cf3d_f17_zero_distance_uses_floor_without_changing_gradient_or_cache() {
        // RDKit source: DistanceConstraints.cpp:87-92 divides by
        // max(1e-8, distance); zero coordinate deltas still produce zero.
        let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];
        let mut force_field = cf3d_f16_force_field(&mut positions);
        cf3d_f17_add_packed_constraint(&mut force_field, 0, 1, 1.0, 2.0, 4.0);
        let mut gradient = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

        force_field.calc_grad(&[0.0; 6], &mut gradient).unwrap();
        assert_eq!(gradient, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(force_field.distance_matrix[1], -1.0);
    }

    #[test]
    fn cf3d_f17_propagates_typed_context_error_before_gradient_writes() {
        // RDKit source: ForceField.cpp:206-213 rejects an uninitialized
        // owner before distance evaluation or this contribution's writes.
        let mut positions = vec![vec![0.0; 3], vec![1.0, 0.0, 0.0]];
        let mut force_field = ForceField::default();
        force_field
            .positions_mut()
            .extend(positions.iter_mut().map(Vec::as_mut_slice));
        let mut contribution = DistanceConstraintContribs::new(&force_field);
        contribution
            .add_contrib(&force_field, 0, 1, 2.0, 3.0, 4.0)
            .unwrap();
        let coordinates = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let context = force_field.evaluation_context(&coordinates);
        let mut gradient = [1.0; 6];

        assert_eq!(
            contribution.get_grad(&context, &mut gradient),
            Err(ForceFieldKernelError::NotInitialized)
        );
        assert_eq!(gradient, [1.0; 6]);
    }
}
