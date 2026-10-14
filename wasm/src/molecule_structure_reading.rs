//! Canonical structural text readers; enabled only with BIO and molecular IO.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn from_pdb(text: &str) -> Result<Self, ck::MolecularIoError> {
        ck::Molecule::from_pdb(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_pdb_with_params(
        text: &str,
        params: &ck::BioMoleculeParams,
    ) -> Result<Self, ck::MolecularIoError> {
        ck::Molecule::from_pdb_with_params(text, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_mmcif(text: &str) -> Result<Self, ck::MolecularIoError> {
        ck::Molecule::from_mmcif(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_mmcif_with_params(
        text: &str,
        params: &ck::BioMoleculeParams,
    ) -> Result<Self, ck::MolecularIoError> {
        ck::Molecule::from_mmcif_with_params(text, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
}
