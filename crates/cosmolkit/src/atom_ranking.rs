use crate::{AtomId, BondId, CanonicalRankError, Molecule};

/// Comparisons used to compute canonical ranks in the source atom order.
#[derive(Debug, Clone, PartialEq, Eq, cosmolkit_macros::Configuration)]
#[configuration(scalar)]
pub struct AtomRankParams {
    /// Assign distinct ranks to symmetry-equivalent atoms.
    pub break_ties: bool,
    /// Compare atom and bond stereochemistry.
    pub include_chirality: bool,
    /// Include isotope mass numbers in atom comparisons.
    pub include_isotopes: bool,
    /// Use atom-map labels in ordinary atom comparisons.
    pub include_atom_maps: bool,
    /// Compare the presence of specified chirality independently of its orientation.
    pub include_chiral_presence: bool,
}

impl Default for AtomRankParams {
    fn default() -> Self {
        let source = cosmolkit_core::CanonicalRankParams::default();
        Self {
            break_ties: source.break_ties,
            include_chirality: source.include_chirality,
            include_isotopes: source.include_isotopes,
            include_atom_maps: source.include_atom_maps,
            include_chiral_presence: source.include_chiral_presence,
        }
    }
}

impl AtomRankParams {
    fn comparison(&self) -> cosmolkit_core::CanonicalRankParams {
        let mut params = cosmolkit_core::CanonicalRankParams::default();
        params.break_ties = self.break_ties;
        params.include_chirality = self.include_chirality;
        params.include_isotopes = self.include_isotopes;
        params.include_atom_maps = self.include_atom_maps;
        params.include_chiral_presence = self.include_chiral_presence;
        params
    }
}

impl Molecule {
    /// Returns ranks in source atom order. Sort atom indices by these values to obtain a permutation.
    pub fn canonical_atom_ranks(&self) -> Result<Vec<usize>, CanonicalRankError> {
        self.canonical_atom_ranks_with_params(&AtomRankParams::default())
    }

    pub fn canonical_atom_ranks_with_params(
        &self,
        params: &AtomRankParams,
    ) -> Result<Vec<usize>, CanonicalRankError> {
        if self.num_atoms() == 0 {
            return Ok(Vec::new());
        }
        let valence = self.ranking_valence()?;
        cosmolkit_core::rank_mol_atoms_with_prepared_state(
            self.topology(),
            &params.comparison(),
            &valence,
            self.derived_cache_runtime().valid_ring_info(),
        )
    }

    fn ranking_valence(
        &self,
    ) -> Result<std::borrow::Cow<'_, cosmolkit_core::ValenceAssignment>, CanonicalRankError> {
        if let Some(valence) = self.derived_cache_runtime().valence_assignment() {
            return Ok(std::borrow::Cow::Borrowed(valence));
        }
        Ok(std::borrow::Cow::Owned(
            cosmolkit_core::assign_valence_with_options_for_topology(
                self.topology(),
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )?,
        ))
    }

    /// Ranks just the specified atoms and bonds; an omitted bond selection includes no bonds.
    /// The output has one entry per source atom, with -1 for unselected atoms.
    pub fn canonical_atom_ranks_in_fragment(
        &self,
        atoms: &[AtomId],
        bonds: &[BondId],
    ) -> Result<Vec<i64>, CanonicalRankError> {
        self.canonical_atom_ranks_in_fragment_with_params(
            atoms,
            bonds,
            None,
            &AtomRankParams::default(),
        )
    }

    pub fn canonical_atom_ranks_in_fragment_with_params(
        &self,
        atoms: &[AtomId],
        bonds: &[BondId],
        atom_symbols: Option<&[String]>,
        params: &AtomRankParams,
    ) -> Result<Vec<i64>, CanonicalRankError> {
        if atoms.is_empty() {
            return Err(CanonicalRankError::EmptyFragment);
        }
        let mut selected = vec![false; self.num_atoms()];
        for atom in atoms {
            *selected.get_mut(atom.index()).ok_or(
                CanonicalRankError::ComponentAtomOutOfRange {
                    atom_index: atom.index(),
                    atom_count: self.num_atoms(),
                },
            )? = true;
        }
        let mut edges = vec![false; self.num_bonds()];
        for bond in bonds {
            *edges
                .get_mut(bond.index())
                .ok_or(CanonicalRankError::FragmentBondOutOfRange {
                    bond_index: bond.index(),
                    bond_count: self.num_bonds(),
                })? = true;
        }
        let valence = self.ranking_valence()?;
        let ranks = cosmolkit_core::rank_fragment_atoms_with_prepared_state(
            self.topology(),
            &valence,
            self.derived_cache_runtime().valid_ring_info(),
            &selected,
            &edges,
            atom_symbols,
            None,
            &params.comparison(),
        )?;
        // Ranks are bounded by the length of an allocated atom vector (at most isize::MAX).
        Ok(ranks
            .into_iter()
            .zip(selected)
            .map(|(rank, used)| if used { rank as i64 } else { -1 })
            .collect())
    }
}
