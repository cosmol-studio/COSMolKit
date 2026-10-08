//! SMILES serialization delegates to the canonical public facade.
impl crate::Molecule {
    pub fn from_smiles(smiles: &str) -> Result<Self, cosmolkit::SmilesError> {
        cosmolkit::Molecule::from_smiles(smiles).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_smiles_with_params(
        smiles: &str,
        params: &cosmolkit::SmilesParseParams,
    ) -> Result<Self, cosmolkit::SmilesError> {
        cosmolkit::Molecule::from_smiles_with_params(smiles, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn to_smiles(&self) -> Result<cosmolkit::PropertyText, cosmolkit::SmilesWriteError> {
        self.inner.borrow().to_smiles()
    }
}
