use std::ops::{BitOr, BitOrAssign};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SanitizeOps(u32);

impl SanitizeOps {
    pub const NONE: Self = Self(0);
    pub const CLEANUP: Self = Self(1 << 0);
    pub const PROPERTIES: Self = Self(1 << 1);
    pub const SYMMRINGS: Self = Self(1 << 2);
    pub const KEKULIZE: Self = Self(1 << 3);
    pub const FIND_RADICALS: Self = Self(1 << 4);
    pub const SET_AROMATICITY: Self = Self(1 << 5);
    pub const SET_CONJUGATION: Self = Self(1 << 6);
    pub const SET_HYBRIDIZATION: Self = Self(1 << 7);
    pub const CLEANUP_CHIRALITY: Self = Self(1 << 8);
    pub const ADJUST_HYDROGENS: Self = Self(1 << 9);
    pub const CLEANUP_ORGANOMETALLICS: Self = Self(1 << 10);
    pub const CLEANUP_ATROPISOMERS: Self = Self(1 << 11);
    pub const ALL: Self = Self((1 << 12) - 1);

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }
}

impl BitOr for SanitizeOps {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for SanitizeOps {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SanitizeStep {
    Cleanup,
    CleanupOrganometallics,
    Properties,
    SymmRings,
    Kekulize,
    FindRadicals,
    SetAromaticity,
    SetConjugation,
    SetHybridization,
    CleanupAtropisomers,
    CleanupChirality,
    AdjustHydrogens,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("sanitize failed at {step:?}: {message}")]
pub struct SanitizeError {
    pub step: SanitizeStep,
    pub message: String,
    pub unsupported: Option<crate::UnsupportedFeatureError>,
}
