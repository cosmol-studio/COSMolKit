use crate::{Molecule, SmilesParseError, ValenceModel, assign_radicals_rdkit_2025, assign_valence};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct SanitizeOps(u32);

impl SanitizeOps {
    pub const NONE: Self = Self(0);
    pub const CLEANUP: Self = Self(1 << 0);
    pub const PROPERTIES: Self = Self(1 << 1);
    pub const SYMMRINGS: Self = Self(1 << 2);
    pub const KEKULIZE: Self = Self(1 << 3);
    pub const FINDRADICALS: Self = Self(1 << 4);
    pub const SETAROMATICITY: Self = Self(1 << 5);
    pub const SETCONJUGATION: Self = Self(1 << 6);
    pub const SETHYBRIDIZATION: Self = Self(1 << 7);
    pub const CLEANUPCHIRALITY: Self = Self(1 << 8);
    pub const ADJUSTHS: Self = Self(1 << 9);
    pub const CLEANUP_ORGANOMETALLICS: Self = Self(1 << 10);
    pub const CLEANUPATROPISOMERS: Self = Self(1 << 11);

    pub const ALL: Self = Self(
        Self::CLEANUP.0
            | Self::CLEANUP_ORGANOMETALLICS.0
            | Self::PROPERTIES.0
            | Self::SYMMRINGS.0
            | Self::KEKULIZE.0
            | Self::FINDRADICALS.0
            | Self::SETAROMATICITY.0
            | Self::SETCONJUGATION.0
            | Self::SETHYBRIDIZATION.0
            | Self::CLEANUPATROPISOMERS.0
            | Self::CLEANUPCHIRALITY.0
            | Self::ADJUSTHS.0,
    );

    pub const SUPPORTED_ALL: Self = Self(
        Self::CLEANUP.0
            | Self::CLEANUP_ORGANOMETALLICS.0
            | Self::PROPERTIES.0
            | Self::SYMMRINGS.0
            | Self::KEKULIZE.0
            | Self::FINDRADICALS.0
            | Self::SETAROMATICITY.0
            | Self::CLEANUPCHIRALITY.0
            | Self::ADJUSTHS.0,
    );

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SanitizeStep {
    Cleanup,
    CleanupOrganometallics,
    Properties,
    SymmRings,
    Kekulize,
    FindRadicals,
    SetAromaticity,
    SetConjugation,
    SetHybridization,
    CleanupAtropisomers,
    CleanupChirality,
    AdjustHs,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("sanitize failed at {step:?}: {message}")]
pub struct SanitizeError {
    pub step: SanitizeStep,
    pub message: String,
}

impl SanitizeError {
    fn new(step: SanitizeStep, message: impl Into<String>) -> Self {
        Self {
            step,
            message: message.into(),
        }
    }
}

impl From<SanitizeError> for SmilesParseError {
    fn from(error: SanitizeError) -> Self {
        Self::ParseError(error.to_string())
    }
}

pub(crate) fn apply_sanitize_pipeline(
    mol: &mut Molecule,
    ops: SanitizeOps,
) -> Result<(), SanitizeError> {
    let mut original_implicit_hydrogens = None;

    if ops.contains(SanitizeOps::CLEANUP) {
        crate::smiles::cleanup_neutral_five_coordinate_nitrogens(mol)
            .map_err(|error| SanitizeError::new(SanitizeStep::Cleanup, error.to_string()))?;
    }

    if ops.contains(SanitizeOps::CLEANUP_ORGANOMETALLICS) {
        crate::smiles::cleanup_organometallic_single_bonds(mol).map_err(|error| {
            SanitizeError::new(SanitizeStep::CleanupOrganometallics, error.to_string())
        })?;
    }

    if ops.contains(SanitizeOps::PROPERTIES) {
        original_implicit_hydrogens = Some(
            assign_valence(mol, ValenceModel::RdkitLike)
                .map_err(|error| {
                    SanitizeError::new(SanitizeStep::Properties, format!("{error:?}"))
                })?
                .implicit_hydrogens,
        );
    } else if ops.contains(SanitizeOps::ADJUSTHS) {
        original_implicit_hydrogens = assign_valence(mol, ValenceModel::RdkitLike)
            .ok()
            .map(|assignment| assignment.implicit_hydrogens);
    }

    if ops.contains(SanitizeOps::SYMMRINGS) {
        // COSMolKit currently computes ring-derived facts on demand instead of
        // storing an RDKit-style symmetrized SSSR cache.
        mol.rebuild_adjacency();
    }

    if ops.contains(SanitizeOps::KEKULIZE) {
        crate::kekulize::kekulize_in_place(mol, true)
            .map_err(|error| SanitizeError::new(SanitizeStep::Kekulize, error.to_string()))?;
    }

    if ops.contains(SanitizeOps::FINDRADICALS) {
        assign_sanitized_radicals(mol)
            .map_err(|error| SanitizeError::new(SanitizeStep::FindRadicals, error))?;
    }

    if ops.contains(SanitizeOps::SETAROMATICITY) {
        crate::smiles::perceive_aromaticity(mol, &[])
            .map_err(|error| SanitizeError::new(SanitizeStep::SetAromaticity, error.to_string()))?;
        crate::smiles::prune_noncyclic_aromatic_bonds(mol);
    }

    if ops.contains(SanitizeOps::SETCONJUGATION) {
        return Err(SanitizeError::new(
            SanitizeStep::SetConjugation,
            "SANITIZE_SETCONJUGATION is not implemented because Molecule does not store conjugation flags yet",
        ));
    }

    if ops.contains(SanitizeOps::SETHYBRIDIZATION) {
        return Err(SanitizeError::new(
            SanitizeStep::SetHybridization,
            "SANITIZE_SETHYBRIDIZATION is not implemented because Molecule does not store hybridization flags yet",
        ));
    }

    if ops.contains(SanitizeOps::CLEANUPATROPISOMERS) {
        return Err(SanitizeError::new(
            SanitizeStep::CleanupAtropisomers,
            "SANITIZE_CLEANUPATROPISOMERS is not implemented because atropisomer stereo is not represented yet",
        ));
    }

    if ops.contains(SanitizeOps::CLEANUPCHIRALITY) {
        crate::smiles::assign_double_bond_stereo_from_directions(mol);
        crate::smiles::cleanup_nonstereo_double_bond_dirs(mol);
    }

    if ops.contains(SanitizeOps::ADJUSTHS) {
        if let Some(original_implicit_hydrogens) = original_implicit_hydrogens {
            crate::hydrogens::adjust_hydrogens_after_aromaticity_in_place(
                mol,
                &original_implicit_hydrogens,
            );
        }
        crate::hydrogens::remove_hydrogens_after_smiles_parse_in_place(mol)
            .map_err(|error| SanitizeError::new(SanitizeStep::AdjustHs, error.to_string()))?;
    }

    if ops.contains(SanitizeOps::CLEANUPCHIRALITY) {
        // COSMolKit's tetrahedral cleanup examines explicit-H bookkeeping, so
        // it must run after the RDKit-aligned H adjustment for current storage.
        crate::smiles::cleanup_invalid_tetrahedral_stereo(mol);
    }

    if ops.contains(SanitizeOps::PROPERTIES) {
        assign_valence(mol, ValenceModel::RdkitLike)
            .map_err(|error| SanitizeError::new(SanitizeStep::Properties, format!("{error:?}")))?;
    }

    crate::stereo::cache_rdkit_legacy_cip_ranks(mol);
    mol.rebuild_adjacency();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{SanitizeOps, SanitizeStep};
    use crate::Molecule;

    #[test]
    fn cleanup_step_cleans_nitro_without_full_pipeline() {
        let raw = Molecule::from_smiles_with_sanitize("CN(=O)=O", false)
            .expect("unsanitized nitro SMILES should parse");

        let cleaned = raw
            .sanitize_with_ops(SanitizeOps::CLEANUP)
            .expect("cleanup should normalize nitro charge form");

        assert_eq!(
            cleaned
                .atoms()
                .iter()
                .map(|atom| atom.formal_charge)
                .collect::<Vec<_>>(),
            vec![0, 1, -1, 0]
        );
        assert_eq!(
            raw.atoms()
                .iter()
                .map(|atom| atom.formal_charge)
                .collect::<Vec<_>>(),
            vec![0, 0, 0, 0],
            "sanitize_with_ops must preserve COW value semantics"
        );
    }

    #[test]
    fn unsupported_sanitize_steps_return_explicit_step_errors() {
        let mol = Molecule::from_smiles("CCO").expect("SMILES should parse");

        let conjugation = mol
            .sanitize_with_ops(SanitizeOps::SETCONJUGATION)
            .expect_err("conjugation storage is not implemented");
        assert_eq!(conjugation.step, SanitizeStep::SetConjugation);

        let hybridization = mol
            .sanitize_with_ops(SanitizeOps::SETHYBRIDIZATION)
            .expect_err("hybridization storage is not implemented");
        assert_eq!(hybridization.step, SanitizeStep::SetHybridization);

        let atropisomers = mol
            .sanitize_with_ops(SanitizeOps::CLEANUPATROPISOMERS)
            .expect_err("atropisomer cleanup is not implemented");
        assert_eq!(atropisomers.step, SanitizeStep::CleanupAtropisomers);
    }

    #[test]
    fn rdkit_all_flags_do_not_silently_skip_unsupported_steps() {
        let mol = Molecule::from_smiles("CCO").expect("SMILES should parse");

        let error = mol
            .sanitize_with_ops(SanitizeOps::ALL)
            .expect_err("full RDKit sanitize flags should expose unsupported steps");

        assert_eq!(error.step, SanitizeStep::SetConjugation);
    }
}

fn assign_sanitized_radicals(mol: &mut Molecule) -> Result<(), String> {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike)
        .map_err(|error| format!("valence assignment failed: {error:?}"))?;
    let radicals = assign_radicals_rdkit_2025(mol, &assignment.explicit_valence)
        .map_err(|error| format!("radical assignment failed: {error:?}"))?;
    for (atom, radical) in mol.atoms_mut().iter_mut().zip(radicals) {
        atom.num_radical_electrons = radical;
    }
    Ok(())
}
