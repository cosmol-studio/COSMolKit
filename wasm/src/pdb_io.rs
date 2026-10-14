//! Thin molecular PDB/mmCIF projections of the canonical Rust API.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn to_pdb(&self) -> Result<String, ck::MolecularIoError> {
        self.inner.borrow().to_pdb()
    }
    pub fn to_pdb_with_params(
        &self,
        params: &ck::PdbWriteParams,
    ) -> Result<String, ck::MolecularIoError> {
        self.inner.borrow().to_pdb_with_params(params)
    }
}
