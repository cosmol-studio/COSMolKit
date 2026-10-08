//! Complete immutable distance-matrix values, policies and canonical query projection.
use crate::Molecule;
use crate::host_values::*;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct DenseMatrix {
    inner: ck::DenseMatrix,
}
#[wasm_bindgen]
impl DenseMatrix {
    pub fn dimension(&self) -> usize {
        self.inner.dimension()
    }
    pub fn values(&self) -> Vec<f64> {
        self.inner.values().to_vec()
    }
    #[wasm_bindgen(unchecked_return_type = "number | null")]
    pub fn get(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] row: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] column: JsValue,
    ) -> Result<JsValue, JsValue> {
        Ok(self
            .inner
            .get(usize_value(&row, "row")?, usize_value(&column, "column")?)
            .map_or(JsValue::NULL, JsValue::from_f64))
    }
}
#[wasm_bindgen]
pub struct DistanceMatrixParams {
    inner: ck::DistanceMatrixParams,
}
#[wasm_bindgen]
impl DistanceMatrixParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_bond_order: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_atom_weights: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::DistanceMatrixParams::default();
        if !use_bond_order.is_undefined() {
            inner.use_bond_order = bool_value(&use_bond_order, "useBondOrder")?;
        }
        if !use_atom_weights.is_undefined() {
            inner.use_atom_weights = bool_value(&use_atom_weights, "useAtomWeights")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=useBondOrder)]
    pub fn use_bond_order(&self) -> bool {
        self.inner.use_bond_order
    }
    #[wasm_bindgen(getter,js_name=useAtomWeights)]
    pub fn use_atom_weights(&self) -> bool {
        self.inner.use_atom_weights
    }
}
#[wasm_bindgen]
pub struct DistanceMatrix3dParams {
    inner: ck::DistanceMatrix3dParams,
}
#[wasm_bindgen]
impl DistanceMatrix3dParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_atom_weights: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::DistanceMatrix3dParams::default();
        if !conformer_id.is_undefined() && !conformer_id.is_null() {
            inner.conformer_id = Some(usize_value(&conformer_id, "conformerId")?);
        }
        if !use_atom_weights.is_undefined() {
            inner.use_atom_weights = bool_value(&use_atom_weights, "useAtomWeights")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=conformerId,unchecked_return_type="number | null")]
    pub fn conformer_id(&self) -> JsValue {
        self.inner
            .conformer_id
            .map_or(JsValue::NULL, |n| JsValue::from(n as u32))
    }
    #[wasm_bindgen(getter,js_name=useAtomWeights)]
    pub fn use_atom_weights(&self) -> bool {
        self.inner.use_atom_weights
    }
}
#[wasm_bindgen(
    inline_js = "export function visitDistanceMatrixParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid DistanceMatrixParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitDistanceMatrixParams)]
    fn visit_distance_matrix(
        v: &JsValue,
        f: &mut dyn FnMut(&DistanceMatrixParams),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen(
    inline_js = "export function visitDistanceMatrix3dParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid DistanceMatrix3dParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitDistanceMatrix3dParams)]
    fn visit_distance_matrix_3d(
        v: &JsValue,
        f: &mut dyn FnMut(&DistanceMatrix3dParams),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=distanceMatrix)]
    pub fn distance_matrix(&self) -> Result<DenseMatrix, JsValue> {
        // COSMolKit❗✔️: self.inner.distance_matrix()
        self.inner
            .distance_matrix()
            .map(|inner| DenseMatrix { inner })
            .map_err(|e| crate::matrix_errors::matrix_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=distanceMatrixWithParams)]
    pub fn distance_matrix_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "DistanceMatrixParams")] params: JsValue,
    ) -> Result<DenseMatrix, JsValue> {
        // COSMolKit❗✔️: self.inner.distance_matrix_with_params(&p.inner)
        let mut result = None;
        visit_distance_matrix(&params, &mut |p: &DistanceMatrixParams| {
            result = Some(
                self.inner
                    .distance_matrix_with_params(&p.inner)
                    .map(|inner| DenseMatrix { inner })
                    .map_err(|e| crate::matrix_errors::matrix_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=distanceMatrix3d)]
    pub fn distance_matrix_3d(&self) -> Result<DenseMatrix, JsValue> {
        // COSMolKit❗✔️: self.inner.distance_matrix_3d()
        self.inner
            .distance_matrix_3d()
            .map(|inner| DenseMatrix { inner })
            .map_err(|e| crate::matrix_errors::matrix_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=distanceMatrix3dWithParams)]
    pub fn distance_matrix_3d_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "DistanceMatrix3dParams")] params: JsValue,
    ) -> Result<DenseMatrix, JsValue> {
        // COSMolKit❗✔️: self.inner.distance_matrix_3d_with_params(&p.inner)
        let mut result = None;
        visit_distance_matrix_3d(&params, &mut |p: &DistanceMatrix3dParams| {
            result = Some(
                self.inner
                    .distance_matrix_3d_with_params(&p.inner)
                    .map(|inner| DenseMatrix { inner })
                    .map_err(|e| crate::matrix_errors::matrix_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
