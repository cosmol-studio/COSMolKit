//! Manual coordinate calls delegate to the canonical operation runtime.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn with_2d_coordinate_block(
        &self,
        coordinates: Vec<Vec<f64>>,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_2d_coordinate_block(coordinates)
            .map(Into::into)
    }
    pub fn with_2d_coordinate_block_with_params(
        &self,
        coordinates: Vec<Vec<f64>>,
        params: &ck::Coordinate2DInputParams,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_2d_coordinate_block_with_params(coordinates, params)
            .map(Into::into)
    }
    pub fn set_2d_coordinates_(
        &self,
        coordinates: Vec<Vec<f64>>,
    ) -> Result<(), ck::OperationError> {
        self.inner.borrow_mut().set_2d_coordinates_(coordinates)
    }
    pub fn set_2d_coordinates_with_params_(
        &self,
        coordinates: Vec<Vec<f64>>,
        params: &ck::Coordinate2DInputParams,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_2d_coordinates_with_params_(coordinates, params)
    }
    pub fn with_3d_coordinates(
        &self,
        coordinates: Vec<Vec<f64>>,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_3d_coordinates(coordinates)
            .map(Into::into)
    }
    pub fn with_3d_coordinates_with_params(
        &self,
        coordinates: Vec<Vec<f64>>,
        params: &ck::Replace3DCoordinatesParams,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_3d_coordinates_with_params(coordinates, params)
            .map(Into::into)
    }
    pub fn set_3d_coordinates_(
        &self,
        coordinates: Vec<Vec<f64>>,
    ) -> Result<(), ck::OperationError> {
        self.inner.borrow_mut().set_3d_coordinates_(coordinates)
    }
    pub fn set_3d_coordinates_with_params_(
        &self,
        coordinates: Vec<Vec<f64>>,
        params: &ck::Replace3DCoordinatesParams,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_3d_coordinates_with_params_(coordinates, params)
    }
    pub fn with_added_3d_conformer(
        &self,
        coordinates: Vec<Vec<f64>>,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_added_3d_conformer(coordinates)
            .map(Into::into)
    }
    pub fn with_added_3d_conformer_with_params(
        &self,
        coordinates: Vec<Vec<f64>>,
        params: &ck::Coordinate3DInputParams,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_added_3d_conformer_with_params(coordinates, params)
            .map(Into::into)
    }
    pub fn add_3d_conformer_(
        &self,
        coordinates: Vec<Vec<f64>>,
    ) -> Result<usize, ck::OperationError> {
        self.inner.borrow_mut().add_3d_conformer_(coordinates)
    }
    pub fn add_3d_conformer_with_params_(
        &self,
        coordinates: Vec<Vec<f64>>,
        params: &ck::Coordinate3DInputParams,
    ) -> Result<usize, ck::OperationError> {
        self.inner
            .borrow_mut()
            .add_3d_conformer_with_params_(coordinates, params)
    }
    pub fn with_only_3d_conformer(
        &self,
        coordinates: Vec<Vec<f64>>,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_only_3d_conformer(coordinates)
            .map(Into::into)
    }
    pub fn with_only_3d_conformer_with_params(
        &self,
        coordinates: Vec<Vec<f64>>,
        params: &ck::Coordinate3DInputParams,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_only_3d_conformer_with_params(coordinates, params)
            .map(Into::into)
    }
    pub fn set_only_3d_conformer_(
        &self,
        coordinates: Vec<Vec<f64>>,
    ) -> Result<usize, ck::OperationError> {
        self.inner.borrow_mut().set_only_3d_conformer_(coordinates)
    }
    pub fn set_only_3d_conformer_with_params_(
        &self,
        coordinates: Vec<Vec<f64>>,
        params: &ck::Coordinate3DInputParams,
    ) -> Result<usize, ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_only_3d_conformer_with_params_(coordinates, params)
    }
    pub fn with_cleared_3d_conformers(&self) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_cleared_3d_conformers()
            .map(Into::into)
    }
    pub fn clear_3d_conformers_(&self) -> Result<(), ck::OperationError> {
        self.inner.borrow_mut().clear_3d_conformers_()
    }
    pub fn with_atom_position(
        &self,
        atom: ck::AtomId,
        position: [f64; 3],
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_position(atom, position)
            .map(Into::into)
    }
    pub fn with_atom_position_with_params(
        &self,
        atom: ck::AtomId,
        position: [f64; 3],
        params: &ck::AtomPositionParams,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_position_with_params(atom, position, params)
            .map(Into::into)
    }
    pub fn set_atom_position_(
        &self,
        atom: ck::AtomId,
        position: [f64; 3],
    ) -> Result<(), ck::OperationError> {
        self.inner.borrow_mut().set_atom_position_(atom, position)
    }
    pub fn set_atom_position_with_params_(
        &self,
        atom: ck::AtomId,
        position: [f64; 3],
        params: &ck::AtomPositionParams,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_atom_position_with_params_(atom, position, params)
    }
}
