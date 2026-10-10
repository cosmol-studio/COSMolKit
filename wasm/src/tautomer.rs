//! Thin canonical tautomer entry points; no enumeration algorithm here.
use crate::Molecule;
use cosmolkit as ck;
pub fn canonical_tautomer_from_molecules_with_params(
    molecules: &[std::sync::Arc<Molecule>],
    params: &ck::TautomerParams,
) -> Result<Molecule, ck::OperationError> {
    ck::canonical_tautomer_from_molecule_hosts_with_params(
        molecules.len(),
        |index, visit| visit(&molecules[index].inner.borrow()),
        params,
    )
    .map(Into::into)
}
impl Molecule {
    pub fn enumerate_tautomers_with_params(
        &self,
        params: &ck::TautomerParams,
    ) -> Result<ck::TautomerEnumeration, ck::OperationError> {
        self.inner.borrow().enumerate_tautomers_with_params(params)
    }
    pub fn canonical_tautomer_with_params(
        &self,
        params: &ck::TautomerParams,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .canonical_tautomer_with_params(params)
            .map(Into::into)
    }
    pub fn tautomer_score_with_params(
        &self,
        params: &ck::TautomerScoreParams,
    ) -> Result<ck::TautomerScore, ck::OperationError> {
        self.inner.borrow().tautomer_score_with_params(params)
    }
}
