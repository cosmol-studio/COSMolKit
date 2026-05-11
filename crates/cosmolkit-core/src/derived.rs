use std::ops::{BitOr, BitOrAssign};

/// Derived state invalidated by molecule operations.
///
/// This is a contract surface. Future agents must add a new bit here before
/// adding a cache or derived perception result that can become stale.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct DerivedState(u32);

impl DerivedState {
    pub const NONE: Self = Self(0);
    pub const ADJACENCY: Self = Self(1 << 0);
    pub const RINGS: Self = Self(1 << 1);
    pub const RING_FAMILIES: Self = Self(1 << 2);
    pub const VALENCE: Self = Self(1 << 3);
    pub const AROMATICITY: Self = Self(1 << 4);
    pub const STEREO: Self = Self(1 << 5);
    pub const COORDINATES: Self = Self(1 << 6);
    pub const DRAWING: Self = Self(1 << 7);
    pub const FINGERPRINT: Self = Self(1 << 8);

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }

    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    #[must_use]
    pub const fn touches_cache(self) -> bool {
        self.contains(Self::ADJACENCY)
            || self.contains(Self::RINGS)
            || self.contains(Self::RING_FAMILIES)
            || self.contains(Self::VALENCE)
            || self.contains(Self::AROMATICITY)
            || self.contains(Self::STEREO)
    }
}

impl BitOr for DerivedState {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for DerivedState {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}
