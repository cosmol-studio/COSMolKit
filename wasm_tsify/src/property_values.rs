//! Exact native property alternatives; no conversion between numeric variants.
use crate::{
    fields::Field,
    input,
    model_errors::{PropertyValueError, PropertyValueKind},
    property_errors::text,
};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen(js_name=propertyValueToText)]
pub fn property_value_to_text(value: PropertyValueInput) -> Result<String, JsValue> {
    PropertyValue::with_ref(&value, |value| {
        ck::property_value_to_text(value)
            .map_err(crate::property_errors::PropertyStringError::from_rust)
            .and_then(|value| text(&value))
    })
}
#[wasm_bindgen]
pub struct PropertyValue(pub(crate) ck::PropertyValue);
crate::fields::value_field!(
    ck::PropertyValue,
    PropertyValue,
    PropertyValueInput,
    "PropertyValue"
);
crate::fields::optional_value_array!(
    ck::PropertyValue,
    PropertyValue,
    PropertyValueArray,
    "ReadonlyArray<PropertyValue | null | undefined>"
);
#[wasm_bindgen]
impl PropertyValue {
    #[wasm_bindgen(js_name=String)]
    pub fn string(value: tsify::Ts<input::TextInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::PropertyValue::String(
            input::read::<String>(value)?.into(),
        )))
    }
    #[wasm_bindgen(js_name=Int)]
    pub fn int(value: tsify::Ts<input::I32Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::PropertyValue::Int(input::read(value)?)))
    }
    #[wasm_bindgen(js_name=UInt)]
    pub fn uint(value: tsify::Ts<input::U32Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::PropertyValue::UInt(input::read(value)?)))
    }
    #[wasm_bindgen(js_name=Double)]
    pub fn double(value: tsify::Ts<input::F64Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::PropertyValue::Double(input::read(value)?)))
    }
    #[wasm_bindgen(js_name=Bool)]
    pub fn boolean(value: tsify::Ts<input::BoolInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::PropertyValue::Bool(input::read(value)?)))
    }
    #[wasm_bindgen(js_name=IntVector)]
    pub fn int_vector(value: <Vec<i32> as Field>::Input) -> Result<Self, JsValue> {
        Ok(Self(ck::PropertyValue::IntVector(Vec::<i32>::input(
            value,
        )?)))
    }
    #[wasm_bindgen(js_name=StringVector)]
    pub fn string_vector(value: <Vec<String> as Field>::Input) -> Result<Self, JsValue> {
        Ok(Self(ck::PropertyValue::StringVector(
            Vec::<String>::input(value)?
                .into_iter()
                .map(Into::into)
                .collect(),
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> PropertyValueKind {
        self.0.kind().into()
    }
    #[wasm_bindgen(js_name=asString)]
    pub fn as_string(&self) -> Result<String, JsValue> {
        text(self.0.as_string().map_err(PropertyValueError::from_rust)?)
    }
    #[wasm_bindgen(js_name=asInt)]
    pub fn as_int(&self) -> Result<i32, JsValue> {
        self.0.as_int().map_err(PropertyValueError::from_rust)
    }
    #[wasm_bindgen(js_name=asUint)]
    pub fn as_uint(&self) -> Result<u32, JsValue> {
        self.0.as_uint().map_err(PropertyValueError::from_rust)
    }
    #[wasm_bindgen(js_name=asDouble)]
    pub fn as_double(&self) -> Result<f64, JsValue> {
        self.0.as_double().map_err(PropertyValueError::from_rust)
    }
    #[wasm_bindgen(js_name=asBool)]
    pub fn as_bool(&self) -> Result<bool, JsValue> {
        self.0.as_bool().map_err(PropertyValueError::from_rust)
    }
    #[wasm_bindgen(js_name=asIntVector)]
    pub fn as_int_vector(&self) -> Result<Vec<i32>, JsValue> {
        self.0
            .as_int_vector()
            .map(<[i32]>::to_vec)
            .map_err(PropertyValueError::from_rust)
    }
}
