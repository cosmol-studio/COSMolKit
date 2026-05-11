use std::ops::{BitOr, BitOrAssign};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TopologyOpKind {
    Strong,
    Weak,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DependentData {
    Coordinates,
    Stereo,
    MoleculeProps,
    AtomProps,
    BondProps,
    AdjacencyCache,
    RingCache,
    ValenceCache,
    AromaticityState,
    DrawingState,
    FingerprintState,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct DerivedState(u32);

impl DerivedState {
    pub(crate) const NONE: Self = Self(0);
    pub(crate) const ADJACENCY: Self = Self(1 << 0);
    pub(crate) const RINGS: Self = Self(1 << 1);
    pub(crate) const VALENCE: Self = Self(1 << 2);
    pub(crate) const AROMATICITY: Self = Self(1 << 3);
    pub(crate) const STEREO: Self = Self(1 << 4);
    pub(crate) const DRAWING: Self = Self(1 << 5);
    pub(crate) const FINGERPRINT: Self = Self(1 << 6);

    pub(crate) const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }

    pub(crate) const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
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

#[derive(Debug)]
pub(crate) struct TopologyOpSpec {
    pub(crate) name: &'static str,
    pub(crate) kind: TopologyOpKind,
    pub(crate) affects: &'static [DependentData],
    pub(crate) invalidates: DerivedState,
    pub(crate) rdkit_parity: bool,
    pub(crate) io_roundtrip: bool,
}

pub(crate) const TOPOLOGY_STATE_INVALIDATES: DerivedState = DerivedState::ADJACENCY
    .union(DerivedState::RINGS)
    .union(DerivedState::VALENCE)
    .union(DerivedState::STEREO)
    .union(DerivedState::DRAWING)
    .union(DerivedState::FINGERPRINT);

pub(crate) static WITH_HYDROGENS_SPEC: TopologyOpSpec = TopologyOpSpec {
    name: "with_hydrogens",
    kind: TopologyOpKind::Strong,
    affects: &[
        DependentData::Coordinates,
        DependentData::Stereo,
        DependentData::MoleculeProps,
        DependentData::AtomProps,
        DependentData::BondProps,
        DependentData::AdjacencyCache,
        DependentData::RingCache,
        DependentData::ValenceCache,
    ],
    invalidates: TOPOLOGY_STATE_INVALIDATES,
    rdkit_parity: true,
    io_roundtrip: true,
};

pub(crate) static WITHOUT_HYDROGENS_SPEC: TopologyOpSpec = TopologyOpSpec {
    name: "without_hydrogens",
    kind: TopologyOpKind::Strong,
    affects: &[
        DependentData::Coordinates,
        DependentData::Stereo,
        DependentData::MoleculeProps,
        DependentData::AtomProps,
        DependentData::BondProps,
        DependentData::AdjacencyCache,
        DependentData::RingCache,
        DependentData::ValenceCache,
    ],
    invalidates: TOPOLOGY_STATE_INVALIDATES,
    rdkit_parity: true,
    io_roundtrip: true,
};

pub(crate) static KEKULIZE_SPEC: TopologyOpSpec = TopologyOpSpec {
    name: "kekulize",
    kind: TopologyOpKind::Weak,
    affects: &[
        DependentData::AromaticityState,
        DependentData::ValenceCache,
        DependentData::DrawingState,
        DependentData::FingerprintState,
    ],
    invalidates: DerivedState::AROMATICITY
        .union(DerivedState::VALENCE)
        .union(DerivedState::DRAWING)
        .union(DerivedState::FINGERPRINT),
    rdkit_parity: true,
    io_roundtrip: true,
};

pub(crate) static SANITIZE_SPEC: TopologyOpSpec = TopologyOpSpec {
    name: "sanitize",
    kind: TopologyOpKind::Weak,
    affects: &[
        DependentData::Stereo,
        DependentData::ValenceCache,
        DependentData::AromaticityState,
        DependentData::DrawingState,
        DependentData::FingerprintState,
    ],
    invalidates: TOPOLOGY_STATE_INVALIDATES,
    rdkit_parity: true,
    io_roundtrip: true,
};
