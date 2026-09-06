//! Detached target data used by search algorithms.
//!
//! This view deliberately contains no live `Molecule` or operation-runtime
//! state. Runtime adapters may construct it from authorized blocks, while the
//! matcher and predicate implementations consume only this explicit value.

use crate::{AdjacencyList, Atom, Bond, RingInfo, ValenceAssignment};
use cosmolkit_model::{CoordinateBlock, StereoGroup, TopologyBlock};

#[derive(Debug, Clone, Copy)]
pub(crate) struct SearchTarget<'a> {
    topology: &'a TopologyBlock,
    coordinates: &'a CoordinateBlock,
    stereo_groups: &'a [StereoGroup],
    ring_info: Option<&'a RingInfo>,
    valence: Option<&'a ValenceAssignment>,
}

/// Read-only model vocabulary required by search implementations.
///
/// Algorithms are generic over this narrow capability instead of accepting a
/// live runtime object. The `Molecule` implementation is a temporary adapter
/// and can move to `cosmolkit` unchanged when runtime ownership is migrated.
pub(crate) trait SearchTargetAccess {
    fn topology_block(&self) -> &TopologyBlock;
    fn coordinate_block(&self) -> &CoordinateBlock;
    fn ring_info(&self) -> Option<&RingInfo>;
    fn valence(&self) -> Option<&ValenceAssignment>;

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

    fn stereo_groups(&self) -> &[StereoGroup] {
        &self.topology_block().stereo_groups
    }

    fn conformers_3d(&self) -> &[cosmolkit_model::Conformer3D] {
        &self.coordinate_block().conformers_3d
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
}

impl<'a> SearchTarget<'a> {
    #[must_use]
    pub(crate) const fn new(
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

    #[must_use]
    pub(crate) fn atoms(self) -> &'a [Atom] {
        &self.topology.atoms
    }

    #[must_use]
    pub(crate) fn bonds(self) -> &'a [Bond] {
        &self.topology.bonds
    }

    #[must_use]
    pub(crate) const fn topology_block(self) -> &'a TopologyBlock {
        self.topology
    }

    #[must_use]
    pub(crate) const fn coordinate_block(self) -> &'a CoordinateBlock {
        self.coordinates
    }

    #[must_use]
    pub(crate) const fn adjacency(self) -> &'a AdjacencyList {
        &self.topology.adjacency
    }

    #[must_use]
    pub(crate) const fn num_atoms(self) -> usize {
        self.topology.atoms.len()
    }

    #[must_use]
    pub(crate) const fn num_bonds(self) -> usize {
        self.topology.bonds.len()
    }

    #[must_use]
    pub(crate) const fn stereo_groups(self) -> &'a [StereoGroup] {
        self.stereo_groups
    }

    #[must_use]
    pub(crate) const fn ring_info(self) -> Option<&'a RingInfo> {
        self.ring_info
    }

    #[must_use]
    pub(crate) const fn valence(self) -> Option<&'a ValenceAssignment> {
        self.valence
    }

    #[must_use]
    pub(crate) fn conformers_3d(self) -> &'a [cosmolkit_model::Conformer3D] {
        &self.coordinates.conformers_3d
    }

    #[must_use]
    pub(crate) fn vf2_num_atoms(self) -> usize {
        self.num_atoms()
    }

    #[must_use]
    pub(crate) fn vf2_num_bonds(self) -> usize {
        self.num_bonds()
    }

    #[must_use]
    pub(crate) fn vf2_bond_endpoints(self, bond_idx: usize) -> (usize, usize) {
        let bond = &self.topology.bonds[bond_idx];
        (bond.begin().index(), bond.end().index())
    }
}
