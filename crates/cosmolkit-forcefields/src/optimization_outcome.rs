//! Typed outcomes of a completed force-field optimization attempt.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OptimizationStatus {
    Converged,
    NotConverged,
    MissingParameters,
}

/// Nonconvergence retains the final energy; unavailable parameters have no energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OptimizationOutcome {
    Converged { energy: f64 },
    NotConverged { energy: f64 },
    MissingParameters,
}

impl OptimizationOutcome {
    pub const fn converged_with_energy(energy: f64) -> Self {
        Self::Converged { energy }
    }
    pub const fn not_converged_with_energy(energy: f64) -> Self {
        Self::NotConverged { energy }
    }
    pub const fn missing_parameters() -> Self {
        Self::MissingParameters
    }
    pub const fn status(&self) -> OptimizationStatus {
        match self {
            Self::Converged { .. } => OptimizationStatus::Converged,
            Self::NotConverged { .. } => OptimizationStatus::NotConverged,
            Self::MissingParameters => OptimizationStatus::MissingParameters,
        }
    }
    pub const fn energy(&self) -> Option<f64> {
        match *self {
            Self::Converged { energy } | Self::NotConverged { energy } => Some(energy),
            Self::MissingParameters => None,
        }
    }
    pub const fn converged(&self) -> bool {
        matches!(self, Self::Converged { .. })
    }
    pub const fn needs_more(&self) -> bool {
        matches!(self, Self::NotConverged { .. })
    }
    /// RDKit status for source comparisons, not an energy sentinel.
    pub const fn status_code(&self) -> i32 {
        match self {
            Self::Converged { .. } => 0,
            Self::NotConverged { .. } => 1,
            Self::MissingParameters => -1,
        }
    }
}

impl From<crate::uff::convenience::OptimizationOutcome> for OptimizationOutcome {
    fn from(value: crate::uff::convenience::OptimizationOutcome) -> Self {
        // The private minimizer produces exactly 0 or 1; missing parameters
        // are handled before building a field, never inferred from energy.
        match value.status {
            0 => Self::Converged {
                energy: value.energy,
            },
            1 => Self::NotConverged {
                energy: value.energy,
            },
            _ => unreachable!("force-field minimizer returned an invalid status"),
        }
    }
}
