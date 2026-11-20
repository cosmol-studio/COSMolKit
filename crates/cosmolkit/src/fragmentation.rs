use crate::{BricsBondMatch, Molecule, SubstructMatchError};

impl Molecule {
    /// Return unique terminal canonical isomeric SMILES from standard BRICS decomposition.
    /// Order is unspecified. An uncut molecule contributes its own SMILES.
    pub fn decompose_brics(&self) -> Result<Vec<String>, crate::BricsDecompositionError> {
        cosmolkit_reaction::decompose_brics(self.brics_input())
    }

    /// Find standard BRICS sites in RDKit rule and match order, leaving the molecule unchanged.
    /// Each bond is returned once, with the labels of the first matching rule.
    pub fn find_brics_bonds(&self) -> Result<Vec<BricsBondMatch>, SubstructMatchError> {
        cosmolkit_search::find_brics_bonds(&self.detached_search_target())
    }
    /// Enumerate unique terminal SMILES with regular r-BRICS, without the reBRICS post-pass.
    pub fn decompose_r_brics(&self) -> Result<Vec<String>, crate::BricsDecompositionError> {
        cosmolkit_reaction::decompose_r_brics(self.brics_input())
    }

    /// Return regular r-BRICS sites and their oriented environment labels in upstream order.
    pub fn find_r_brics_bonds(&self) -> Result<Vec<crate::RBricsBondMatch>, SubstructMatchError> {
        cosmolkit_search::find_r_brics_bonds(&self.detached_search_target())
    }

    fn brics_input(&self) -> cosmolkit_reaction::ReactionInput<'_> {
        let cache = self.derived_cache_runtime();
        cosmolkit_reaction::ReactionInput {
            topology: self.topology(),
            coordinates: self.coordinate_block_runtime(),
            properties: self.properties(),
            rings: cache.valid_ring_info(),
            valence: cache.valence_assignment(),
        }
    }
}
