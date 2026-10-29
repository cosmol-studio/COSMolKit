//! Detached metadata values use the public model's typed property accessors.
use crate::{
    enums::{SdfPropertyListTarget, SdfPropertyListTargetInput},
    fields::Field,
    input,
    model_errors::PropertyValueError,
    property_errors::text,
    property_values::{PropertyValue, PropertyValueArray},
    sdf_records::SdfDataFields,
};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub struct SdfPropertyList(pub(crate) ck::SdfPropertyList);
crate::fields::value_field!(
    ck::SdfPropertyList,
    SdfPropertyList,
    SdfPropertyListInput,
    "SdfPropertyList"
);
#[wasm_bindgen]
impl SdfPropertyList {
    #[wasm_bindgen(constructor)]
    pub fn new(
        target: SdfPropertyListTargetInput,
        name: Ts<input::TextInput>,
        values: PropertyValueArray,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::SdfPropertyList::new(
            ck::SdfPropertyListTarget::input(target)?,
            input::read::<String>(name)?,
            values.native_values()?,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn target(&self) -> SdfPropertyListTarget {
        self.0.target().into()
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> Result<String, JsValue> {
        text(self.0.name())
    }
    #[wasm_bindgen(getter)]
    pub fn values(&self) -> PropertyValueArray {
        PropertyValueArray::from_native(self.0.values().to_vec())
    }
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ReadonlyArray<SdfPropertyList>")]
    pub type SdfPropertyLists;
    #[wasm_bindgen(typescript_type = "ReadonlyMap<string, PropertyValue>")]
    pub type PropertyMap;
}
#[wasm_bindgen]
pub struct MoleculeProperties(pub(crate) ck::MoleculeProperties);
crate::fields::value_field!(
    ck::MoleculeProperties,
    MoleculeProperties,
    MoleculePropertiesInput,
    "MoleculeProperties"
);
#[wasm_bindgen]
impl MoleculeProperties {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self(ck::MoleculeProperties::default())
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> Result<Option<String>, JsValue> {
        self.0.name().map(text).transpose()
    }
    #[wasm_bindgen(getter,js_name=sdfDataFields)]
    pub fn sdf_data_fields(&self) -> Result<Ts<SdfDataFields>, JsValue> {
        SdfDataFields(
            self.0
                .sdf_data_fields()
                .iter()
                .map(|(k, v)| Ok((text(k)?, text(v)?)))
                .collect::<Result<Vec<_>, JsValue>>()?,
        )
        .into_ts()
        .map_err(wasm_bindgen::JsError::from)
        .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=sdfPropertyLists)]
    pub fn sdf_property_lists(&self) -> SdfPropertyLists {
        let a = self
            .0
            .sdf_property_lists()
            .iter()
            .map(|v| JsValue::from(SdfPropertyList(v.clone())))
            .collect::<js_sys::Array>();
        js_sys::Object::freeze(&a).unchecked_into()
    }
    #[wasm_bindgen(getter)]
    pub fn props(&self) -> Result<PropertyMap, JsValue> {
        let map = js_sys::Map::new();
        for (k, v) in self.0.props() {
            map.set(&text(k)?.into(), &PropertyValue(v.clone()).into());
        }
        Ok(crate::fields::readonly_map(&map)?.unchecked_into())
    }
    pub fn prop(&self, key: Ts<input::TextInput>) -> Result<Option<PropertyValue>, JsValue> {
        Ok(self
            .0
            .prop(input::read::<String>(key)?)
            .cloned()
            .map(PropertyValue))
    }
    #[wasm_bindgen(js_name=isPropComputed)]
    pub fn is_prop_computed(&self, key: Ts<input::TextInput>) -> Result<bool, JsValue> {
        self.0
            .is_prop_computed(input::read::<String>(key)?)
            .map_err(PropertyValueError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=computedPropNames)]
    pub fn computed_prop_names(&self) -> Result<Option<Vec<String>>, JsValue> {
        self.0
            .computed_prop_names()
            .map_err(PropertyValueError::from_rust)?
            .map(|names| names.iter().map(text).collect())
            .transpose()
    }
}
