//! Complete depiction registry projection, with checked dimensions and borrowed parameters.
use crate::{
    Molecule, alignment_values::operation_error, drawing_errors::drawing_error,
    host_values::u32_value, image_errors::drawing_write_error,
    transform_parameters::Coordinate2DParams,
};
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=add2dCoordinates)]
    pub fn add_2d_coordinates(&self, #[wasm_bindgen(unchecked_optional_param_type = "Coordinate2DParams | Coordinate2DParamsOptions")] params: JsValue) -> Result<Molecule, JsValue> {
        if params.is_undefined() {
            return {
        // COSMolKit❗✔️: .add_2d_coordinates()
        self.inner
            .add_2d_coordinates()
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))

            };
        }
        let configuration_0 = crate::transform_parameters::Coordinate2DParams::from_js_configuration(&params)?;
        self.add_2d_coordinates_with_params(&configuration_0)
    }
    #[wasm_bindgen(js_name=add2dCoordinatesWithParams)]
    pub fn add_2d_coordinates_with_params(
        &self,
        params: &Coordinate2DParams,
    ) -> Result<Molecule, JsValue> {
        // COSMolKit❗✔️: .add_2d_coordinates_with_params(&params.inner)
        self.inner
            .add_2d_coordinates_with_params(&params.inner)
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toSvg)]
    pub fn to_svg(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] width: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] height: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: .to_svg(width, height)
        self.inner
            .to_svg(u32_value(&width, "width")?, u32_value(&height, "height")?)
            .map_err(|e| drawing_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toPng)]
    pub fn to_png(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] width: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] height: JsValue,
    ) -> Result<Vec<u8>, JsValue> {
        // COSMolKit❗✔️: .to_png(width, height)
        self.inner
            .to_png(u32_value(&width, "width")?, u32_value(&height, "height")?)
            .map_err(|e| drawing_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=add2dCoordinates_)]
    pub fn add_2d_coordinates_(&self, #[wasm_bindgen(unchecked_optional_param_type = "Coordinate2DParams | Coordinate2DParamsOptions")] params: JsValue) -> Result<(), JsValue> {
        if params.is_undefined() {
            return {
        // COSMolKit❗✔️: .add_2d_coordinates_()
        self.inner
            .add_2d_coordinates_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))

            };
        }
        let configuration_0 = crate::transform_parameters::Coordinate2DParams::from_js_configuration(&params)?;
        self.add_2d_coordinates_with_params_(&configuration_0)
    }
    #[wasm_bindgen(js_name=add2dCoordinatesWithParams_)]
    pub fn add_2d_coordinates_with_params_(
        &self,
        params: &Coordinate2DParams,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .add_2d_coordinates_with_params_(&params.inner)
        self.inner
            .add_2d_coordinates_with_params_(&params.inner)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writeSvg)]
    pub fn write_svg(
        &self,
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "number")] width: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] height: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .write_svg(std::path::Path::new(path), width, height)
        self.inner
            .write_svg(
                path,
                u32_value(&width, "width")?,
                u32_value(&height, "height")?,
            )
            .map_err(|e| drawing_write_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writePng)]
    pub fn write_png(
        &self,
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "number")] width: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] height: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .write_png(std::path::Path::new(path), width, height)
        self.inner
            .write_png(
                path,
                u32_value(&width, "width")?,
                u32_value(&height, "height")?,
            )
            .map_err(|e| drawing_write_error(&e).unwrap_or_else(|e| e))
    }
}
