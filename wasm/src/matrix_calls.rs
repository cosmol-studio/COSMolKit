//! Dense matrices retain the native row-major shape and checked lookup.
use crate::{Molecule, input, matrix_errors::MatrixError};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
ck::__ck_configuration_DistanceMatrixParams!(
    crate::project_parameters,
    ck::DistanceMatrixParams,
    {},
    { DistanceMatrixParams }
);
ck::__ck_configuration_DistanceMatrix3dParams!(
    crate::project_parameters,
    ck::DistanceMatrix3dParams,
    {},
    { DistanceMatrix3dParams }
);
#[wasm_bindgen]
pub struct DenseMatrix(ck::DenseMatrix);
#[wasm_bindgen]
impl DenseMatrix {
    #[wasm_bindgen(getter)]
    pub fn dimension(&self) -> usize {
        self.0.dimension()
    }
    #[wasm_bindgen(getter)]
    pub fn values(&self) -> Vec<f64> {
        self.0.values().to_vec()
    }
    pub fn get(
        &self,
        row: tsify::Ts<input::IndexInput>,
        column: tsify::Ts<input::IndexInput>,
    ) -> Result<Option<f64>, JsValue> {
        Ok(self.0.get(input::read(row)?, input::read(column)?))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=distanceMatrix)]
    pub fn distance_matrix(
        &self,
        params: Option<DistanceMatrixParamsInput>,
    ) -> Result<DenseMatrix, JsValue> {
        DistanceMatrixParams::with_input(params.as_ref(), |p| {
            self.0
                .distance_matrix_with_params(p)
                .map(DenseMatrix)
                .map_err(MatrixError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=distanceMatrix3d)]
    pub fn distance_matrix_3d(
        &self,
        params: Option<DistanceMatrix3dParamsInput>,
    ) -> Result<DenseMatrix, JsValue> {
        DistanceMatrix3dParams::with_input(params.as_ref(), |p| {
            self.0
                .distance_matrix_3d_with_params(p)
                .map(DenseMatrix)
                .map_err(MatrixError::from_rust)
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = distanceMatrixWithParams)]
    pub fn distance_matrix_with_params(
        &self,
        params: DistanceMatrixParamsInput,
    ) -> Result<DenseMatrix, JsValue> {
        self.distance_matrix(Some(params))
    }
    #[wasm_bindgen(js_name = distanceMatrix3dWithParams)]
    pub fn distance_matrix_3d_with_params(
        &self,
        params: DistanceMatrix3dParamsInput,
    ) -> Result<DenseMatrix, JsValue> {
        self.distance_matrix_3d(Some(params))
    }
}
