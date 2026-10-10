//! SMILES serialization delegates to the canonical public facade.
impl crate::Molecule {
    pub fn to_cx_smiles_with_params(
        &self,
        params: &cosmolkit::CxSmilesWriteParams,
    ) -> Result<cosmolkit::PropertyText, cosmolkit::SmilesWriteError> {
        self.inner.borrow().to_cx_smiles_with_params(params)
    }
    pub fn to_fragment_smiles(
        &self,
        atoms: &[cosmolkit::AtomId],
    ) -> Result<cosmolkit::PropertyText, cosmolkit::SmilesWriteError> {
        self.inner.borrow().to_fragment_smiles(atoms)
    }
    pub fn to_fragment_cx_smiles(
        &self,
        atoms: &[cosmolkit::AtomId],
    ) -> Result<cosmolkit::PropertyText, cosmolkit::SmilesWriteError> {
        self.inner.borrow().to_fragment_cx_smiles(atoms)
    }
    pub fn to_fragment_smiles_with_params(
        &self,
        params: &cosmolkit::FragmentSmilesWriteParams,
    ) -> Result<cosmolkit::PropertyText, cosmolkit::SmilesWriteError> {
        self.inner.borrow().to_fragment_smiles_with_params(params)
    }
    pub fn to_fragment_cx_smiles_with_params(
        &self,
        params: &cosmolkit::FragmentCxSmilesWriteParams,
    ) -> Result<cosmolkit::PropertyText, cosmolkit::SmilesWriteError> {
        self.inner
            .borrow()
            .to_fragment_cx_smiles_with_params(params)
    }
    pub fn to_random_smiles_with_params(
        &self,
        count: u32,
        seed: u32,
        params: &cosmolkit::RandomSmilesWriteParams,
    ) -> Result<Vec<cosmolkit::PropertyText>, cosmolkit::SmilesWriteError> {
        self.inner
            .borrow()
            .to_random_smiles_with_params(count, seed, params)
    }
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
    pub fn to_smiles_with_params(
        &self,
        params: &cosmolkit::SmilesWriteParams,
    ) -> Result<cosmolkit::PropertyText, cosmolkit::SmilesWriteError> {
        self.inner.borrow().to_smiles_with_params(params)
    }
}
