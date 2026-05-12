use crate::{
    Atom, AtomId, Bond, Molecule, MoleculeProperties,
    molecule::{CoordinateBlock, DerivedCacheBlock, TopologyBlock},
    sgroup::SubstanceGroup,
    stereo::StereoGroup,
};

#[derive(Clone, Copy)]
pub(crate) struct MoleculeReadParts<'a> {
    molecule: &'a Molecule,
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
        Self { molecule }
    }

    #[must_use]
    pub(crate) fn atoms(self) -> &'a [Atom] {
        self.molecule.atoms()
    }

    #[must_use]
    pub(crate) fn topology(self) -> &'a TopologyBlock {
        self.molecule.topology_block()
    }

    #[must_use]
    pub(crate) fn coordinates(self) -> &'a CoordinateBlock {
        self.molecule.coordinate_block()
    }

    #[must_use]
    pub(crate) fn bonds(self) -> &'a [Bond] {
        self.molecule.bonds()
    }

    #[must_use]
    pub(crate) fn atom(self, atom: AtomId) -> Option<&'a Atom> {
        self.molecule.atom(atom)
    }

    #[must_use]
    pub(crate) fn num_atoms(self) -> usize {
        self.molecule.num_atoms()
    }

    #[must_use]
    pub(crate) fn num_bonds(self) -> usize {
        self.molecule.num_bonds()
    }

    #[must_use]
    pub(crate) fn coords_2d(self) -> Option<&'a [[f64; 2]]> {
        self.molecule.coords_2d()
    }

    #[must_use]
    pub(crate) fn conformers_3d(self) -> &'a [crate::Conformer3D] {
        self.molecule.conformers_3d()
    }

    #[must_use]
    pub(crate) fn substance_groups(self) -> &'a [SubstanceGroup] {
        self.molecule.substance_groups()
    }

    #[allow(dead_code)]
    #[must_use]
    pub(crate) fn stereo_groups(self) -> &'a [StereoGroup] {
        self.molecule.stereo_groups()
    }

    #[must_use]
    pub(crate) fn properties(self) -> &'a MoleculeProperties {
        self.molecule.properties()
    }

    #[must_use]
    pub(crate) fn derived_cache(self) -> &'a DerivedCacheBlock {
        self.molecule.derived_cache()
    }

    fn with_molecule_read<R>(self, f: impl FnOnce(&Molecule) -> R) -> R {
        f(self.molecule)
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
        self.with_molecule_read(|molecule| {
            crate::assign_valence_with_options(molecule, model, strict)
        })
    }

    pub(crate) fn assign_radicals(self) -> Result<Vec<u8>, crate::ValenceError> {
        self.with_molecule_read(crate::assign_radicals)
    }

    pub(crate) fn symmetrize_sssr(self) -> Result<crate::RingInfo, crate::RingFindingError> {
        self.with_molecule_read(crate::symmetrize_sssr)
    }

    pub(crate) fn find_ring_families(
        self,
        include_dative_bonds: bool,
        include_hydrogen_bonds: bool,
    ) -> Result<crate::RingInfo, crate::RingFindingError> {
        self.with_molecule_read(|molecule| {
            crate::find_ring_families(molecule, include_dative_bonds, include_hydrogen_bonds)
        })
    }

    pub(crate) fn set_aromaticity(
        self,
        model: crate::AromaticityModel,
    ) -> Result<crate::AromaticityAssignment, crate::AromaticityError> {
        self.with_molecule_read(|molecule| crate::set_aromaticity(molecule, model))
    }

    pub(crate) fn kekulize_assignment(
        self,
        ring_info: Option<&crate::RingInfo>,
        clear_aromatic_flags: bool,
        mark_atoms_bonds: bool,
        max_backtracks: usize,
    ) -> Result<crate::kekulize::KekulizeAssignment, crate::KekulizeError> {
        self.with_molecule_read(|molecule| {
            crate::kekulize::kekulize_assignment(
                molecule,
                ring_info,
                clear_aromatic_flags,
                mark_atoms_bonds,
                max_backtracks,
            )
        })
    }

    pub(crate) fn rank_mol_atoms(self) -> Result<Vec<usize>, crate::KekulizeError> {
        self.with_molecule_read(crate::canon_rank::rank_mol_atoms)
    }
}
