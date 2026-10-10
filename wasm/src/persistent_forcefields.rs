//! Owned force-field factories delegated to the public facade.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn mmff_force_field(&self) -> Result<ck::MolecularForceField, ck::MmffForceFieldError> {
        self.inner.borrow().mmff_force_field()
    }
    pub fn mmff_force_field_with_params(
        &self,
        params: &ck::MmffForceFieldParams,
    ) -> Result<ck::MolecularForceField, ck::MmffForceFieldError> {
        self.inner.borrow().mmff_force_field_with_params(params)
    }
    pub fn uff_force_field(&self) -> Result<ck::MolecularForceField, ck::UffForceFieldError> {
        self.inner.borrow().uff_force_field()
    }
    pub fn uff_force_field_with_params(
        &self,
        params: &ck::UffForceFieldParams,
    ) -> Result<ck::MolecularForceField, ck::UffForceFieldError> {
        self.inner.borrow().uff_force_field_with_params(params)
    }
}
