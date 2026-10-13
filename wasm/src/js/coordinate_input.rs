//! Checked numeric transport; interpretation remains in canonical coordinate APIs.
use crate::construction::matrix;
use crate::{Molecule, host_values::*};
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
fn operation<T>(result: Result<T, ck::OperationError>) -> Result<T, JsValue> {
    result.map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen]
pub struct CoordinateZPolicy {
    inner: ck::CoordinateZPolicy,
}
#[wasm_bindgen]
impl CoordinateZPolicy {
    #[wasm_bindgen(getter,js_name=Ignore)]
    pub fn ignore() -> Self {
        Self {
            inner: ck::CoordinateZPolicy::Ignore,
        }
    }
    #[wasm_bindgen(getter,js_name=RequireZero)]
    pub fn require_zero() -> Self {
        Self {
            inner: ck::CoordinateZPolicy::RequireZero,
        }
    }
    #[wasm_bindgen(getter,js_name=Error)]
    pub fn error() -> Self {
        Self {
            inner: ck::CoordinateZPolicy::Error,
        }
    }
    #[wasm_bindgen(js_name=fromName)]
    pub fn from_name(value: &str) -> Result<Self, JsValue> {
        ck::CoordinateZPolicy::from_name(value)
            .map(|inner| Self { inner })
            .map_err(|e| crate::coordinate_errors::input_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen(
    inline_js = "export function visitCoordinateZPolicy(v,f){try{f(v);}catch(cause){throw new TypeError('invalid CoordinateZPolicy',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitCoordinateZPolicy)]
    fn visit_z_policy(v: &JsValue, f: &mut dyn FnMut(&CoordinateZPolicy)) -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct Coordinate2DInputParams {
    inner: ck::Coordinate2DInputParams,
}
#[cosmolkit_wasm::javascript_options("Coordinate2DInputParamsOptions")]
#[wasm_bindgen]
impl Coordinate2DInputParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "CoordinateZPolicy | string")]
        z_policy: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::Coordinate2DInputParams::default();
        if let Some(s) = z_policy.as_string() {
            inner.z_policy = ck::CoordinateZPolicy::from_name(&s)
                .map_err(|e| crate::coordinate_errors::input_error(&e).unwrap_or_else(|e| e))?;
        } else if !z_policy.is_undefined() {
            visit_z_policy(&z_policy, &mut |v| inner.z_policy = v.inner)?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=zPolicy)]
    pub fn z_policy(&self) -> CoordinateZPolicy {
        CoordinateZPolicy {
            inner: self.inner.z_policy,
        }
    }
}
#[wasm_bindgen]
pub struct Coordinate3DInputParams {
    inner: ck::Coordinate3DInputParams,
}
#[cosmolkit_wasm::javascript_options("Coordinate3DInputParamsOptions")]
#[wasm_bindgen]
impl Coordinate3DInputParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] is_3d: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::Coordinate3DInputParams::default();
        if !is_3d.is_undefined() {
            inner.is_3d = bool_value(&is_3d, "is3d")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=is3d)]
    pub fn is_3d(&self) -> bool {
        self.inner.is_3d
    }
}
#[wasm_bindgen]
pub struct Replace3DCoordinatesParams {
    inner: ck::Replace3DCoordinatesParams,
}
#[cosmolkit_wasm::javascript_options("Replace3DCoordinatesParamsOptions")]
#[wasm_bindgen]
impl Replace3DCoordinatesParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] conformer_id: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::Replace3DCoordinatesParams::default();
        if !conformer_id.is_undefined() {
            inner.conformer_id = usize_value(&conformer_id, "conformerId")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=conformerId)]
    pub fn conformer_id(&self) -> usize {
        self.inner.conformer_id
    }
}
#[wasm_bindgen]
pub struct AtomPositionParams {
    inner: ck::AtomPositionParams,
}
#[cosmolkit_wasm::javascript_options("AtomPositionParamsOptions")]
#[wasm_bindgen]
impl AtomPositionParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] conformer_id: JsValue,
    ) -> Result<Self, JsValue> {
        let conformer_id = if conformer_id.is_null() || conformer_id.is_undefined() {
            None
        } else {
            Some(usize_value(&conformer_id, "conformerId")?)
        };
        Ok(Self {
            inner: ck::AtomPositionParams { conformer_id },
        })
    }
    #[wasm_bindgen(getter,js_name=conformerId,unchecked_return_type="number | null")]
    pub fn conformer_id(&self) -> JsValue {
        self.inner
            .conformer_id
            .map_or(JsValue::NULL, |v| JsValue::from(v as u32))
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=coordinates3d)]
    pub fn coordinates_3d(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] conformer_id: JsValue,
    ) -> Result<Vec<f64>, JsValue> {
        let id = if conformer_id.is_undefined() {
            0
        } else {
            usize_value(&conformer_id, "conformerId")?
        };
        self.inner
            .coordinates_3d(id)
            .map_err(|e| crate::coordinate_errors::read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=with2dCoordinateBlock)]
    pub fn with_2d_coordinate_block(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Coordinate2DInputParams | Coordinate2DInputParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Molecule, JsValue> {
        self.with_2d_coordinate_block_with_params(
            coordinates,
            &Coordinate2DInputParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=with2dCoordinateBlockWithParams)]
    pub fn with_2d_coordinate_block_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        params: &Coordinate2DInputParams,
    ) -> Result<Molecule, JsValue> {
        operation(
            self.inner
                .with_2d_coordinate_block_with_params(matrix(&coordinates)?, &params.inner),
        )
        .map(|inner| Molecule {
            inner: inner.into(),
        })
    }
    #[wasm_bindgen(js_name=set2dCoordinates_)]
    pub fn set_2d_coordinates_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Coordinate2DInputParams | Coordinate2DInputParamsOptions"
        )]
        params: JsValue,
    ) -> Result<(), JsValue> {
        self.set_2d_coordinates_with_params_(
            coordinates,
            &Coordinate2DInputParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=set2dCoordinatesWithParams_)]
    pub fn set_2d_coordinates_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        params: &Coordinate2DInputParams,
    ) -> Result<(), JsValue> {
        operation(
            self.inner
                .set_2d_coordinates_with_params_(matrix(&coordinates)?, &params.inner),
        )
    }
    #[wasm_bindgen(js_name=with3dCoordinates)]
    pub fn with_3d_coordinates(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Replace3DCoordinatesParams | Replace3DCoordinatesParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Molecule, JsValue> {
        self.with_3d_coordinates_with_params(
            coordinates,
            &Replace3DCoordinatesParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=with3dCoordinatesWithParams)]
    pub fn with_3d_coordinates_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        params: &Replace3DCoordinatesParams,
    ) -> Result<Molecule, JsValue> {
        operation(
            self.inner
                .with_3d_coordinates_with_params(matrix(&coordinates)?, &params.inner),
        )
        .map(|inner| Molecule {
            inner: inner.into(),
        })
    }
    #[wasm_bindgen(js_name=set3dCoordinates_)]
    pub fn set_3d_coordinates_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Replace3DCoordinatesParams | Replace3DCoordinatesParamsOptions"
        )]
        params: JsValue,
    ) -> Result<(), JsValue> {
        self.set_3d_coordinates_with_params_(
            coordinates,
            &Replace3DCoordinatesParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=set3dCoordinatesWithParams_)]
    pub fn set_3d_coordinates_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        params: &Replace3DCoordinatesParams,
    ) -> Result<(), JsValue> {
        operation(
            self.inner
                .set_3d_coordinates_with_params_(matrix(&coordinates)?, &params.inner),
        )
    }
    #[wasm_bindgen(js_name=withAdded3dConformer)]
    pub fn with_added_3d_conformer(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Coordinate3DInputParams | Coordinate3DInputParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Molecule, JsValue> {
        self.with_added_3d_conformer_with_params(
            coordinates,
            &Coordinate3DInputParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=withAdded3dConformerWithParams)]
    pub fn with_added_3d_conformer_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        params: &Coordinate3DInputParams,
    ) -> Result<Molecule, JsValue> {
        operation(
            self.inner
                .with_added_3d_conformer_with_params(matrix(&coordinates)?, &params.inner),
        )
        .map(|inner| Molecule {
            inner: inner.into(),
        })
    }
    #[wasm_bindgen(js_name=add3dConformer_)]
    pub fn add_3d_conformer_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Coordinate3DInputParams | Coordinate3DInputParamsOptions"
        )]
        params: JsValue,
    ) -> Result<usize, JsValue> {
        self.add_3d_conformer_with_params_(
            coordinates,
            &Coordinate3DInputParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=add3dConformerWithParams_)]
    pub fn add_3d_conformer_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        params: &Coordinate3DInputParams,
    ) -> Result<usize, JsValue> {
        operation(
            self.inner
                .add_3d_conformer_with_params_(matrix(&coordinates)?, &params.inner),
        )
    }
    #[wasm_bindgen(js_name=withOnly3dConformer)]
    pub fn with_only_3d_conformer(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Coordinate3DInputParams | Coordinate3DInputParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Molecule, JsValue> {
        self.with_only_3d_conformer_with_params(
            coordinates,
            &Coordinate3DInputParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=withOnly3dConformerWithParams)]
    pub fn with_only_3d_conformer_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        params: &Coordinate3DInputParams,
    ) -> Result<Molecule, JsValue> {
        operation(
            self.inner
                .with_only_3d_conformer_with_params(matrix(&coordinates)?, &params.inner),
        )
        .map(|inner| Molecule {
            inner: inner.into(),
        })
    }
    #[wasm_bindgen(js_name=setOnly3dConformer_)]
    pub fn set_only_3d_conformer_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Coordinate3DInputParams | Coordinate3DInputParamsOptions"
        )]
        params: JsValue,
    ) -> Result<usize, JsValue> {
        self.set_only_3d_conformer_with_params_(
            coordinates,
            &Coordinate3DInputParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=setOnly3dConformerWithParams_)]
    pub fn set_only_3d_conformer_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[][]")] coordinates: JsValue,
        params: &Coordinate3DInputParams,
    ) -> Result<usize, JsValue> {
        operation(
            self.inner
                .set_only_3d_conformer_with_params_(matrix(&coordinates)?, &params.inner),
        )
    }
    #[wasm_bindgen(js_name=withCleared3dConformers)]
    pub fn with_cleared_3d_conformers(&self) -> Result<Molecule, JsValue> {
        operation(self.inner.with_cleared_3d_conformers()).map(|inner| Molecule {
            inner: inner.into(),
        })
    }
    #[wasm_bindgen(js_name=clear3dConformers_)]
    pub fn clear_3d_conformers_(&self) -> Result<(), JsValue> {
        operation(self.inner.clear_3d_conformers_())
    }
    #[wasm_bindgen(js_name=withAtomPosition)]
    pub fn with_atom_position(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        #[wasm_bindgen(unchecked_param_type = "[number, number, number]")] position: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "AtomPositionParams | AtomPositionParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Molecule, JsValue> {
        self.with_atom_position_with_params(
            atom,
            position,
            &AtomPositionParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=withAtomPositionWithParams)]
    pub fn with_atom_position_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        #[wasm_bindgen(unchecked_param_type = "[number, number, number]")] position: JsValue,
        params: &AtomPositionParams,
    ) -> Result<Molecule, JsValue> {
        operation(self.inner.with_atom_position_with_params(
            ck::AtomId::new(usize_value(&atom, "atom")?),
            point(&position)?,
            &params.inner,
        ))
        .map(|inner| Molecule {
            inner: inner.into(),
        })
    }
    #[wasm_bindgen(js_name=setAtomPosition_)]
    pub fn set_atom_position_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        #[wasm_bindgen(unchecked_param_type = "[number, number, number]")] position: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "AtomPositionParams | AtomPositionParamsOptions"
        )]
        params: JsValue,
    ) -> Result<(), JsValue> {
        self.set_atom_position_with_params_(
            atom,
            position,
            &AtomPositionParams::from_js_configuration(&params)?,
        )
    }
    #[wasm_bindgen(js_name=setAtomPositionWithParams_)]
    pub fn set_atom_position_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom: JsValue,
        #[wasm_bindgen(unchecked_param_type = "[number, number, number]")] position: JsValue,
        params: &AtomPositionParams,
    ) -> Result<(), JsValue> {
        operation(self.inner.set_atom_position_with_params_(
            ck::AtomId::new(usize_value(&atom, "atom")?),
            point(&position)?,
            &params.inner,
        ))
    }
}
