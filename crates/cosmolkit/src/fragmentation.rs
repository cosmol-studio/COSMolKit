use crate::{BricsBondMatch, Molecule, SubstructMatchError};

impl Molecule {
    /// Find standard BRICS sites in RDKit rule and match order, leaving the molecule unchanged.
    /// Each bond is returned once, with the labels of the first matching rule.
    pub fn find_brics_bonds(&self) -> Result<Vec<BricsBondMatch>, SubstructMatchError> {
        cosmolkit_search::find_brics_bonds(&self.detached_search_target())
    }
}
