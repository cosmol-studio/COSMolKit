//! Coordinate layout remains a native operation; configurations are generated.
use crate::{Molecule, operation_errors::OperationError};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ReadonlyArray<readonly [number, number]>")]
    pub type Coordinates2D;
}
#[cfg(feature = "cap-depict")]
ck::__ck_configuration_Coordinate2DParams!(
    crate::project_parameters,
    ck::Coordinate2DParams,
    {},
    { Coordinate2DParams }
);
#[wasm_bindgen]
impl Molecule {
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=with2dCoordinates)]
    pub fn with_2d_coordinates(
        &self,
        params: Option<Coordinate2DParamsInput>,
    ) -> Result<Self, JsValue> {
        Coordinate2DParams::with_input(params.as_ref(), |p| {
            self.0
                .with_2d_coordinates_with_params(p)
                .map(Self)
                .map_err(OperationError::from_rust)
        })
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=add2dCoordinates_)]
    pub fn add_2d_coordinates(
        &mut self,
        params: Option<Coordinate2DParamsInput>,
    ) -> Result<(), JsValue> {
        Coordinate2DParams::with_input(params.as_ref(), |p| {
            self.0
                .add_2d_coordinates_with_params_(p)
                .map_err(OperationError::from_rust)
        })
    }
    #[wasm_bindgen(getter,js_name=has2dCoordinates)]
    pub fn has_2d_coordinates(&self) -> bool {
        self.0.has_2d_coordinates()
    }
    #[wasm_bindgen(getter,js_name=coordinates2d)]
    pub fn coordinates_2d(&self) -> Option<Coordinates2D> {
        self.0.coordinates_2d().map(|rows| {
            let values = js_sys::Array::new();
            for row in rows {
                let value = row
                    .iter()
                    .copied()
                    .map(JsValue::from)
                    .collect::<js_sys::Array>();
                js_sys::Object::freeze(&value);
                values.push(&value);
            }
            js_sys::Object::freeze(&values);
            JsValue::from(values).unchecked_into()
        })
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=with2dCoordinatesWithParams)]
    pub fn with_2d_coordinates_with_params(
        &self,
        params: Coordinate2DParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_2d_coordinates(Some(params))
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=add2dCoordinatesWithParams_)]
    pub fn add_2d_coordinates_with_params(
        &mut self,
        params: Coordinate2DParamsInput,
    ) -> Result<(), JsValue> {
        self.add_2d_coordinates(Some(params))
    }
}
