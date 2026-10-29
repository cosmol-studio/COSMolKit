//! Avalon projection delegates only to the canonical molecule facade.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn fingerprint_avalon(&self) -> Result<ck::Fingerprint, ck::AvalonFingerprintError> {
        self.inner.borrow().fingerprint_avalon()
    }
    pub fn fingerprint_avalon_with_params(
        &self,
        params: &ck::AvalonFingerprintParams,
    ) -> Result<ck::Fingerprint, ck::AvalonFingerprintError> {
        self.inner.borrow().fingerprint_avalon_with_params(params)
    }
}
