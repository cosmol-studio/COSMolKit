//! Thin read-only path-score facade; no live-state mutation or ring preparation.
use crate::{Molecule, TopologicalTorsionPathScoreError};
impl Molecule {
    pub fn topological_torsion_path_score(
        &self,
        path: &[usize],
        size: usize,
        atom_codes: Option<&[u32]>,
    ) -> Result<u64, TopologicalTorsionPathScoreError> {
        cosmolkit_fingerprints::topological_torsion_path_score(
            self.topology(),
            self.properties(),
            self.derived_cache_runtime().valence_assignment(),
            path,
            size,
            atom_codes,
        )
    }
}
