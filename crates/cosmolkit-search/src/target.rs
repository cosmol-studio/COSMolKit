//! Explicit detached target data consumed by search algorithms.

use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{
    AdjacencyList, Atom, Bond, Conformer3D, CoordinateBlock, StereoGroup, TopologyBlock,
};

/// Read-only target capability required by query evaluation.
///
/// This trait deliberately exposes model blocks and value assignments rather
/// than a live `Molecule`. Implementations cannot provide mutation or runtime
/// cache authority to the search crate.
pub trait SearchTargetAccess {
    fn topology_block(&self) -> &TopologyBlock;
    fn coordinate_block(&self) -> &CoordinateBlock;
    fn ring_info(&self) -> Option<&RingInfo>;
    fn valence(&self) -> Option<&ValenceAssignment>;

    /// Atomic number visible to query evaluation. Ordinary targets use the
    /// validated model value; a detached depiction target may carry RDKit's
    /// temporary non-element sentinel without changing its topology.
    fn query_atomic_number(&self, atom: &Atom) -> u8 {
        atom.atomic_number()
    }

    fn stereo_groups(&self) -> &[StereoGroup] {
        &self.topology_block().stereo_groups
    }

    fn atoms(&self) -> &[Atom] {
        &self.topology_block().atoms
    }

    fn bonds(&self) -> &[Bond] {
        &self.topology_block().bonds
    }

    fn adjacency(&self) -> &AdjacencyList {
        &self.topology_block().adjacency
    }

    fn num_atoms(&self) -> usize {
        self.atoms().len()
    }

    fn num_bonds(&self) -> usize {
        self.bonds().len()
    }

    fn conformers_3d(&self) -> &[Conformer3D] {
        &self.coordinate_block().conformers_3d
    }
}

/// Borrowed detached search input assembled by the facade.
#[derive(Debug, Clone, Copy)]
pub struct SearchTarget<'a> {
    topology: &'a TopologyBlock,
    coordinates: &'a CoordinateBlock,
    stereo_groups: &'a [StereoGroup],
    ring_info: Option<&'a RingInfo>,
    valence: Option<&'a ValenceAssignment>,
    atomic_number_overrides: Option<&'a [Option<u8>]>,
}

impl<'a> SearchTarget<'a> {
    #[must_use]
    pub const fn new(
        topology: &'a TopologyBlock,
        coordinates: &'a CoordinateBlock,
        stereo_groups: &'a [StereoGroup],
        ring_info: Option<&'a RingInfo>,
        valence: Option<&'a ValenceAssignment>,
    ) -> Self {
        Self {
            topology,
            coordinates,
            stereo_groups,
            ring_info,
            valence,
            atomic_number_overrides: None,
        }
    }

    /// Attach temporary query-visible atomic numbers aligned with target atom
    /// indices. This does not modify canonical `Element` or the target block.
    pub fn with_atomic_number_overrides(mut self, overrides: &'a [Option<u8>]) -> Self {
        assert_eq!(overrides.len(), self.topology.atoms.len());
        self.atomic_number_overrides = Some(overrides);
        self
    }
}

impl SearchTargetAccess for SearchTarget<'_> {
    fn topology_block(&self) -> &TopologyBlock {
        self.topology
    }

    fn coordinate_block(&self) -> &CoordinateBlock {
        self.coordinates
    }

    fn ring_info(&self) -> Option<&RingInfo> {
        self.ring_info
    }

    fn valence(&self) -> Option<&ValenceAssignment> {
        self.valence
    }

    fn query_atomic_number(&self, atom: &Atom) -> u8 {
        // RDKit❗✔️: constexpr int DUMMY_ATOMIC_NUM = 200;
        // RDKit❗✔️: for (auto &at : rs_mol.atoms()) {
        // RDKit❗✔️:   if (!rs_atoms.test(at->getIdx())) {
        // RDKit❗✔️:     at->setAtomicNum(DUMMY_ATOMIC_NUM);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: only the query-visible number changes; unlike the source
        // clone, the typed model atom remains a valid element. Other source
        // effects of setAtomicNum must be audited by the owning caller.
        // Complexity: one indexed optional read per query atom, with no clone
        // or allocation in this hot path, versus an O(V+E) source mol clone.
        self.atomic_number_overrides
            .and_then(|overrides| overrides[atom.id().index()])
            .unwrap_or_else(|| atom.atomic_number())
    }

    fn stereo_groups(&self) -> &[StereoGroup] {
        self.stereo_groups
    }
}
