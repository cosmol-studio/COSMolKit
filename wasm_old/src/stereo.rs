//! Canonical stereo operations and read-only queries.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn with_chiral_tags_from_structure_with_params(
        &self,
        params: &ck::StructureTagParams,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_chiral_tags_from_structure_with_params(params)
            .map(Into::into)
    }
    pub fn assign_chiral_tags_from_structure_with_params_(
        &self,
        params: &ck::StructureTagParams,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .assign_chiral_tags_from_structure_with_params_(params)
    }
    pub fn potential_stereo_with_params(
        &self,
        params: &ck::PotentialStereoParams,
    ) -> Result<ck::PotentialStereoResult, ck::OperationError> {
        self.inner.borrow().potential_stereo_with_params(params)
    }
    pub fn with_cip_labels_with_options(
        &self,
        options: &ck::CipLabelOptions,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_cip_labels_with_options(options)
            .map(Into::into)
    }
    pub fn assign_cip_labels_with_options_(
        &self,
        options: &ck::CipLabelOptions,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .assign_cip_labels_with_options_(options)
    }
    pub fn cip_computed(&self) -> Result<bool, ck::PropertyValueError> {
        self.inner.borrow().cip_computed()
    }
    pub fn tetrahedral_stereo(&self) -> Result<Vec<ck::TetrahedralStereo>, ck::StereoReadError> {
        self.inner.borrow().tetrahedral_stereo()
    }
    pub fn perceive_stereochemistry(&self) -> Result<(), ck::StereoReadError> {
        self.inner.borrow().perceive_stereochemistry()
    }
}
