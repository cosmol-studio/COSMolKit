use crate::{
    AdjacencyList, Atom, AtomId, Bond, Molecule, MoleculeCapabilities, MoleculeProperties,
    molecule::{CoordinateBlock, DerivedCacheBlock, TopologyBlock},
    sgroup::SubstanceGroup,
    stereo::StereoGroup,
};

#[derive(Clone, Copy)]
pub(crate) struct MoleculeReadParts<'a> {
    topology: &'a TopologyBlock,
    coordinates: &'a CoordinateBlock,
    properties: &'a MoleculeProperties,
    derived_cache: &'a DerivedCacheBlock,
    capabilities: MoleculeCapabilities,
}

impl<'a> MoleculeReadParts<'a> {
    // Agent guardrail:
    // MoleculeReadParts is the read capability surface for operation bodies.
    // Agents must not extend, weaken, or bypass this type to regain raw
    // Molecule access from registered operation code. If a helper needs more
    // read capability, add a narrow accessor only after confirming it matches
    // operation policy; do not expose the underlying Molecule.
    #[must_use]
    pub(crate) fn from_molecule(molecule: &'a Molecule) -> Self {
        Self::from_blocks(
            molecule.topology_block(),
            molecule.coordinate_block(),
            molecule.properties(),
            molecule.derived_cache(),
            molecule.capabilities(),
        )
    }

    #[must_use]
    pub(crate) const fn from_blocks(
        topology: &'a TopologyBlock,
        coordinates: &'a CoordinateBlock,
        properties: &'a MoleculeProperties,
        derived_cache: &'a DerivedCacheBlock,
        capabilities: MoleculeCapabilities,
    ) -> Self {
        Self {
            topology,
            coordinates,
            properties,
            derived_cache,
            capabilities,
        }
    }

    #[must_use]
    pub(crate) fn atoms(self) -> &'a [Atom] {
        &self.topology.atoms
    }

    #[must_use]
    pub(crate) fn topology(self) -> &'a TopologyBlock {
        self.topology
    }

    #[must_use]
    pub(crate) fn coordinates(self) -> &'a CoordinateBlock {
        self.coordinates
    }

    #[must_use]
    pub(crate) fn bonds(self) -> &'a [Bond] {
        &self.topology.bonds
    }

    #[must_use]
    pub(crate) fn adjacency(self) -> &'a AdjacencyList {
        &self.topology.adjacency
    }

    #[must_use]
    pub(crate) fn atom(self, atom: AtomId) -> Option<&'a Atom> {
        self.topology.atoms.get(atom.index())
    }

    #[must_use]
    pub(crate) fn num_atoms(self) -> usize {
        self.topology.atoms.len()
    }

    #[must_use]
    pub(crate) fn num_bonds(self) -> usize {
        self.topology.bonds.len()
    }

    #[must_use]
    pub(crate) fn coordinates_2d(self) -> Option<&'a [[f64; 2]]> {
        self.coordinates
            .conformers_2d
            .first()
            .map(crate::Conformer2D::coordinates)
    }

    #[must_use]
    pub(crate) fn conformers_3d(self) -> &'a [crate::Conformer3D] {
        &self.coordinates.conformers_3d
    }

    #[must_use]
    pub(crate) fn substance_groups(self) -> &'a [SubstanceGroup] {
        &self.topology.substance_groups
    }

    #[allow(dead_code)]
    #[must_use]
    pub(crate) fn stereo_groups(self) -> &'a [StereoGroup] {
        &self.topology.stereo_groups
    }

    #[must_use]
    pub(crate) fn properties(self) -> &'a MoleculeProperties {
        self.properties
    }

    #[must_use]
    pub(crate) fn derived_cache(self) -> &'a DerivedCacheBlock {
        self.derived_cache
    }

    #[must_use]
    pub(crate) fn capabilities(self) -> MoleculeCapabilities {
        self.capabilities
    }

    pub(crate) fn add_hs_assignment(
        self,
        params: &crate::hydrogens::AddHsParams,
    ) -> Result<crate::hydrogens::AddHsAssignment, crate::AddHydrogensError> {
        crate::hydrogens::add_hs_assignment(self, params)
    }

    pub(crate) fn remove_hs_assignment(
        self,
        params: &crate::hydrogens::RemoveHsParams,
        sanitize: bool,
    ) -> Result<crate::hydrogens::RemoveHsAssignment, crate::RemoveHydrogensError> {
        crate::hydrogens::remove_hs_assignment(self, params, sanitize)
    }

    pub(crate) fn assign_valence_with_options(
        self,
        model: crate::ValenceModel,
        strict: bool,
    ) -> Result<crate::ValenceAssignment, crate::ValenceError> {
        crate::valence::assign_valence_with_options_from_read_parts(self, model, strict)
    }

    pub(crate) fn assign_radicals(self) -> Result<Vec<u8>, crate::ValenceError> {
        crate::valence::assign_radicals_from_read_parts(self)
    }

    pub(crate) fn symmetrize_sssr(self) -> Result<crate::RingInfo, crate::RingFindingError> {
        crate::rings::symmetrize_sssr_from_read_parts(self)
    }

    pub(crate) fn find_ring_families(
        self,
        include_dative_bonds: bool,
        include_hydrogen_bonds: bool,
    ) -> Result<crate::RingInfo, crate::RingFindingError> {
        crate::rings::find_ring_families_from_read_parts(
            self,
            include_dative_bonds,
            include_hydrogen_bonds,
        )
    }

    pub(crate) fn set_aromaticity(
        self,
        model: crate::AromaticityModel,
    ) -> Result<crate::AromaticityAssignment, crate::AromaticityError> {
        crate::aromaticity::set_aromaticity_from_read_parts(self, model)
    }

    pub(crate) fn kekulize_assignment(
        self,
        ring_info: Option<&crate::RingInfo>,
        clear_aromatic_flags: bool,
        mark_atoms_bonds: bool,
        max_backtracks: usize,
    ) -> Result<crate::kekulize::KekulizeAssignment, crate::KekulizeError> {
        crate::kekulize::kekulize_assignment_from_read_parts(
            self,
            ring_info,
            clear_aromatic_flags,
            mark_atoms_bonds,
            max_backtracks,
        )
    }

    pub(crate) fn rank_mol_atoms(self) -> Result<Vec<usize>, crate::KekulizeError> {
        crate::canon_rank::rank_mol_atoms_from_read_parts(self)
    }
}
