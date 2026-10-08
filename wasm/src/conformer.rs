//! Canonical conformer operations and report ownership; no embedding algorithms.
use crate::Molecule;
use cosmolkit as ck;
#[derive(Clone)]
pub struct EmbedMoleculeResult {
    inner: ck::EmbedMoleculeResult,
}
impl EmbedMoleculeResult {
    pub fn molecule(&self) -> Molecule {
        // COSMolKit❗✔️: inner: self.inner.molecule().clone(),
        Molecule {
            inner: self.inner.molecule().clone().into(),
        }
    }
    pub fn params(&self) -> &ck::EmbedParams {
        self.inner.params()
    }
    pub fn conf_id(&self) -> i32 {
        self.inner.conf_id()
    }
    pub fn ok(&self) -> bool {
        self.inner.ok()
    }
}
#[derive(Clone)]
pub struct EmbedMultipleConfsResult {
    inner: ck::EmbedMultipleConfsResult,
}
impl EmbedMultipleConfsResult {
    pub fn molecule(&self) -> Molecule {
        // COSMolKit❗✔️: inner: self.inner.molecule().clone(),
        Molecule {
            inner: self.inner.molecule().clone().into(),
        }
    }
    pub fn params(&self) -> &ck::EmbedParams {
        self.inner.params()
    }
    pub fn conf_ids(&self) -> &[i32] {
        self.inner.conf_ids()
    }
    pub fn requested_num_confs(&self) -> u32 {
        self.inner.requested_num_confs()
    }
    pub fn generated_count(&self) -> usize {
        self.inner.generated_count()
    }
}
impl Molecule {
    pub fn with_3d_conformer(&self) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: .with_3d_conformer()
        self.inner.borrow().with_3d_conformer().map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn with_3d_conformer_with_params(
        &self,
        params: &ck::EmbedParams,
    ) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: .with_3d_conformer_with_params(params)
        self.inner
            .borrow()
            .with_3d_conformer_with_params(params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn embed_3d_conformer_(&self) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: .embed_3d_conformer_()
        self.inner.borrow_mut().embed_3d_conformer_()
    }
    pub fn embed_3d_conformer_with_params_(
        &self,
        params: &ck::EmbedParams,
    ) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: .embed_3d_conformer_with_params_(params)
        self.inner
            .borrow_mut()
            .embed_3d_conformer_with_params_(params)
    }
    pub fn with_3d_conformer_result(&self) -> Result<EmbedMoleculeResult, ck::OperationError> {
        // COSMolKit❗✔️: .with_3d_conformer_result()
        self.inner
            .borrow()
            .with_3d_conformer_result()
            .map(|inner| EmbedMoleculeResult { inner })
    }
    pub fn with_3d_conformer_result_with_params(
        &self,
        params: &ck::EmbedParams,
    ) -> Result<EmbedMoleculeResult, ck::OperationError> {
        // COSMolKit❗✔️: .with_3d_conformer_result_with_params(params)
        self.inner
            .borrow()
            .with_3d_conformer_result_with_params(params)
            .map(|inner| EmbedMoleculeResult { inner })
    }
    pub fn embed_3d_conformer_result_(&self) -> Result<EmbedMoleculeResult, ck::OperationError> {
        // COSMolKit❗✔️: .embed_3d_conformer_result_()
        self.inner
            .borrow_mut()
            .embed_3d_conformer_result_()
            .map(|inner| EmbedMoleculeResult { inner })
    }
    pub fn embed_3d_conformer_result_with_params_(
        &self,
        params: &ck::EmbedParams,
    ) -> Result<EmbedMoleculeResult, ck::OperationError> {
        // COSMolKit❗✔️: .embed_3d_conformer_result_with_params_(params)
        self.inner
            .borrow_mut()
            .embed_3d_conformer_result_with_params_(params)
            .map(|inner| EmbedMoleculeResult { inner })
    }
    pub fn with_3d_conformers(&self, num_confs: u32) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: .with_3d_conformers(num_confs)
        self.inner
            .borrow()
            .with_3d_conformers(num_confs)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn with_3d_conformers_with_params(
        &self,
        num_confs: u32,
        params: &ck::EmbedParams,
    ) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: .with_3d_conformers_with_params(num_confs, params)
        self.inner
            .borrow()
            .with_3d_conformers_with_params(num_confs, params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn embed_3d_conformers_(&self, num_confs: u32) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: .embed_3d_conformers_(num_confs)
        self.inner.borrow_mut().embed_3d_conformers_(num_confs)
    }
    pub fn embed_3d_conformers_with_params_(
        &self,
        num_confs: u32,
        params: &ck::EmbedParams,
    ) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: .embed_3d_conformers_with_params_(num_confs, params)
        self.inner
            .borrow_mut()
            .embed_3d_conformers_with_params_(num_confs, params)
    }
    pub fn with_3d_conformers_result(
        &self,
        num_confs: u32,
    ) -> Result<EmbedMultipleConfsResult, ck::OperationError> {
        // COSMolKit❗✔️: .with_3d_conformers_result(num_confs)
        self.inner
            .borrow()
            .with_3d_conformers_result(num_confs)
            .map(|inner| EmbedMultipleConfsResult { inner })
    }
    pub fn with_3d_conformers_result_with_params(
        &self,
        num_confs: u32,
        params: &ck::EmbedParams,
    ) -> Result<EmbedMultipleConfsResult, ck::OperationError> {
        // COSMolKit❗✔️: .with_3d_conformers_result_with_params(num_confs, params)
        self.inner
            .borrow()
            .with_3d_conformers_result_with_params(num_confs, params)
            .map(|inner| EmbedMultipleConfsResult { inner })
    }
    pub fn embed_3d_conformers_result_(
        &self,
        num_confs: u32,
    ) -> Result<EmbedMultipleConfsResult, ck::OperationError> {
        // COSMolKit❗✔️: .embed_3d_conformers_result_(num_confs)
        self.inner
            .borrow_mut()
            .embed_3d_conformers_result_(num_confs)
            .map(|inner| EmbedMultipleConfsResult { inner })
    }
    pub fn embed_3d_conformers_result_with_params_(
        &self,
        num_confs: u32,
        params: &ck::EmbedParams,
    ) -> Result<EmbedMultipleConfsResult, ck::OperationError> {
        // COSMolKit❗✔️: .embed_3d_conformers_result_with_params_(num_confs, params)
        self.inner
            .borrow_mut()
            .embed_3d_conformers_result_with_params_(num_confs, params)
            .map(|inner| EmbedMultipleConfsResult { inner })
    }
    pub fn num_3d_conformers(&self) -> usize {
        self.inner.borrow().num_3d_conformers()
    }
    pub fn dg_bounds_matrix(&self) -> Result<Vec<Vec<f64>>, ck::ConformerRunError> {
        self.inner.borrow().dg_bounds_matrix()
    }
}
