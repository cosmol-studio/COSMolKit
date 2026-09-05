use crate::{
    AdjacencyList, Atom, AtomId, Bond, Molecule, MoleculeCapabilities, MoleculeProperties,
    molecule::{CoordinateBlock, DerivedCacheBlock, TopologyBlock},
    sgroup::SubstanceGroup,
    stereo::StereoGroup,
};

#[derive(Clone, Copy)]
pub(crate) struct MoleculeReadAccess(u8);

impl MoleculeReadAccess {
    pub(crate) const NONE: Self = Self(0);
    pub(crate) const TOPOLOGY: Self = Self(1 << 0);
    pub(crate) const COORDINATES: Self = Self(1 << 1);
    pub(crate) const PROPERTIES: Self = Self(1 << 2);
    pub(crate) const DERIVED_CACHE: Self = Self(1 << 3);
    pub(crate) const ALL: Self =
        Self(Self::TOPOLOGY.0 | Self::COORDINATES.0 | Self::PROPERTIES.0 | Self::DERIVED_CACHE.0);

    #[must_use]
    pub(crate) const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[must_use]
    const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }
}

#[derive(Clone, Copy)]
pub(crate) struct MoleculeReadParts<'a> {
    source: Option<&'a Molecule>,
    topology: &'a TopologyBlock,
    coordinates: &'a CoordinateBlock,
    properties: &'a MoleculeProperties,
    derived_cache: &'a DerivedCacheBlock,
    capabilities: MoleculeCapabilities,
    #[cfg(feature = "op-contracts")]
    access: MoleculeReadAccess,
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
        Self::from_molecule_with_access(molecule, MoleculeReadAccess::ALL)
    }

    #[must_use]
    pub(crate) fn from_molecule_with_access(molecule: &'a Molecule, access: MoleculeReadAccess) -> Self {
        Self {
            source: Some(molecule),
            topology: molecule.topology_block(),
            coordinates: molecule.coordinate_block(),
            properties: molecule.properties(),
            derived_cache: molecule.derived_cache(),
            capabilities: molecule.capabilities(),
            #[cfg(feature = "op-contracts")]
            access,
        }
    }

    #[must_use]
    pub(crate) const fn from_blocks(
        topology: &'a TopologyBlock,
        coordinates: &'a CoordinateBlock,
        properties: &'a MoleculeProperties,
        derived_cache: &'a DerivedCacheBlock,
        capabilities: MoleculeCapabilities,
    ) -> Self {
        Self::from_blocks_with_access(
            topology,
            coordinates,
            properties,
            derived_cache,
            capabilities,
            MoleculeReadAccess::ALL,
        )
    }

    #[must_use]
    pub(crate) const fn from_blocks_with_access(
        topology: &'a TopologyBlock,
        coordinates: &'a CoordinateBlock,
        properties: &'a MoleculeProperties,
        derived_cache: &'a DerivedCacheBlock,
        capabilities: MoleculeCapabilities,
        access: MoleculeReadAccess,
    ) -> Self {
        #[cfg(not(feature = "op-contracts"))]
        let _ = access;
        Self {
            source: None,
            topology,
            coordinates,
            properties,
            derived_cache,
            capabilities,
            #[cfg(feature = "op-contracts")]
            access,
        }
    }

    #[inline]
    fn assert_access(self, required: MoleculeReadAccess, block: &'static str) {
        #[cfg(feature = "op-contracts")]
        assert!(
            self.access.contains(required),
            "operation read capability does not permit access to the {block} block"
        );
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = required;
            let _ = block;
        }
    }

    #[must_use]
    pub(crate) fn atoms(self) -> &'a [Atom] {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        &self.topology.atoms
    }

    #[must_use]
    pub(crate) fn topology(self) -> &'a TopologyBlock {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        self.topology
    }

    #[must_use]
    pub(crate) fn coordinates(self) -> &'a CoordinateBlock {
        self.assert_access(MoleculeReadAccess::COORDINATES, "coordinates");
        self.coordinates
    }

    #[must_use]
    pub(crate) fn bonds(self) -> &'a [Bond] {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        &self.topology.bonds
    }

    #[must_use]
    pub(crate) fn adjacency(self) -> &'a AdjacencyList {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        &self.topology.adjacency
    }

    #[must_use]
    pub(crate) fn atom(self, atom: AtomId) -> Option<&'a Atom> {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        self.topology.atoms.get(atom.index())
    }

    #[must_use]
    pub(crate) fn num_atoms(self) -> usize {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        self.topology.atoms.len()
    }

    #[must_use]
    pub(crate) fn num_bonds(self) -> usize {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        self.topology.bonds.len()
    }

    #[must_use]
    pub(crate) fn coordinates_2d(self) -> Option<&'a [[f64; 2]]> {
        self.assert_access(MoleculeReadAccess::COORDINATES, "coordinates");
        self.coordinates
            .conformers_2d
            .first()
            .map(crate::Conformer2D::coordinates)
    }

    #[must_use]
    pub(crate) fn conformers_3d(self) -> &'a [crate::Conformer3D] {
        self.assert_access(MoleculeReadAccess::COORDINATES, "coordinates");
        &self.coordinates.conformers_3d
    }

    #[must_use]
    pub(crate) fn substance_groups(self) -> &'a [SubstanceGroup] {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        &self.topology.substance_groups
    }

    #[allow(dead_code)]
    #[must_use]
    pub(crate) fn stereo_groups(self) -> &'a [StereoGroup] {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        &self.topology.stereo_groups
    }

    #[must_use]
    pub(crate) fn properties(self) -> &'a MoleculeProperties {
        self.assert_access(MoleculeReadAccess::PROPERTIES, "properties");
        self.properties
    }

    #[must_use]
    pub(crate) fn derived_cache(self) -> &'a DerivedCacheBlock {
        self.assert_access(MoleculeReadAccess::DERIVED_CACHE, "derived cache");
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

    pub(crate) fn find_sssr(self) -> Result<crate::RingInfo, crate::RingFindingError> {
        crate::rings::find_sssr_from_parts(self.num_atoms(), self.bonds(), self.adjacency())
    }

    pub(crate) fn find_ring_families(
        self,
        include_dative_bonds: bool,
        include_hydrogen_bonds: bool,
    ) -> Result<crate::RingInfo, crate::RingFindingError> {
        crate::rings::find_ring_families_from_read_parts(self, include_dative_bonds, include_hydrogen_bonds)
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

    pub(crate) fn canonical_isomeric_smiles(self) -> Result<String, crate::SmilesWriteError> {
        self.assert_access(MoleculeReadAccess::TOPOLOGY, "topology");
        self.assert_access(MoleculeReadAccess::COORDINATES, "coordinates");
        self.assert_access(MoleculeReadAccess::PROPERTIES, "properties");
        self.assert_access(MoleculeReadAccess::DERIVED_CACHE, "derived cache");

        let owned;
        let molecule = if let Some(source) = self.source {
            source
        } else {
            owned = Molecule::from_operation_blocks(
                self.topology.clone(),
                self.coordinates.clone(),
                self.properties.clone(),
                self.derived_cache.clone(),
                self.capabilities,
            )
            .map_err(|source| crate::SmilesWriteError::MoleculeInvariant { source })?;
            &owned
        };
        crate::smiles_write::mol_to_smiles(molecule, &crate::SmilesWriteParams::default())
    }

    pub(crate) fn kekulize_assignment_with_valence(
        self,
        ring_info: Option<&crate::RingInfo>,
        valence: Option<&crate::ValenceAssignment>,
        clear_aromatic_flags: bool,
        canonical: bool,
        max_backtracks: usize,
    ) -> Result<crate::kekulize::KekulizeAssignment, crate::KekulizeError> {
        crate::kekulize::kekulize_assignment_from_read_parts_with_valence(
            self,
            ring_info,
            valence,
            clear_aromatic_flags,
            canonical,
            max_backtracks,
        )
    }
}
