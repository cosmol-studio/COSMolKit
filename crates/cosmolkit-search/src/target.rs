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
        }
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

    fn stereo_groups(&self) -> &[StereoGroup] {
        self.stereo_groups
    }
}
