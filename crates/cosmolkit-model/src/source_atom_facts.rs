//! Detached source facts, without chemistry computation or runtime cache authority.

/// The actual signed-byte values stored by the pinned source atom.
///
/// Negative values retain the source's uninitialized/failed getter state. A
/// producer copies genuine calculation results with the source storage width;
/// this value neither computes valence nor decides when a getter may run.
#[doc(hidden)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SourceAtomValenceFacts {
    pub explicit_valence: i8,
    pub implicit_valence: i8,
}

impl SourceAtomValenceFacts {
    pub const UNINITIALIZED: Self = Self {
        // RDKit✔️✔️:   d_implicitValence = -1;
        // RDKit✔️✔️:   d_explicitValence = -1;
        explicit_valence: -1,
        implicit_valence: -1,
    };
}

impl Default for SourceAtomValenceFacts {
    fn default() -> Self {
        // RDKit✔️✔️:   d_implicitValence = -1;
        // RDKit✔️✔️:   d_explicitValence = -1;
        // Behavior: only the two source initialization facts are represented;
        // no source getter or preparation algorithm executes in MODEL.
        // Complexity: two scalar initializations, O(1), no allocation.
        Self::UNINITIALIZED
    }
}
