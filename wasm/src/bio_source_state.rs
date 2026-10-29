//! Native source metadata, preserving original numeric sentinels and records.
use crate::{
    bio_geometry::{BioTransform, BioTransformInput},
    bio_hierarchy::BioStructure,
    fields::{Field, I32ListMap, TextArray, TextMap},
    input,
    storage::NativeField,
};
use cosmolkit as ck;
use std::collections::BTreeMap;
use tsify::Ts;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct BioStructureSourceState(ck::BioStructureSourceState);
crate::fields::value_field!(
    ck::BioStructureSourceState,
    BioStructureSourceState,
    BioStructureSourceStateInput,
    "BioStructureSourceState"
);
#[wasm_bindgen]
impl BioStructureSourceState {
    #[wasm_bindgen(constructor)]
    pub fn new(
        name: Option<Ts<input::TextInput>>,
        resolution: Option<Ts<input::F64Input>>,
        conect_map: Option<Ts<I32ListMap>>,
        has_d_fraction: Option<Ts<input::BoolInput>>,
        non_ascii_line: Option<Ts<input::I32Input>>,
        ter_status: Option<Ts<input::U8Input>>,
        has_origx: Option<Ts<input::BoolInput>>,
        origx: Option<BioTransformInput>,
        info: Option<Ts<TextMap>>,
        raw_remarks: Option<Ts<TextArray>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioStructureSourceState::default();
        Ok(Self(ck::BioStructureSourceState {
            name: name
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.name),
            resolution: resolution
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.resolution),
            conect_map: conect_map
                .map(<BTreeMap<i32, Vec<i32>> as Field>::input)
                .transpose()?
                .unwrap_or(d.conect_map),
            has_d_fraction: has_d_fraction
                .map(input::read::<bool>)
                .transpose()?
                .unwrap_or(d.has_d_fraction),
            non_ascii_line: non_ascii_line
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.non_ascii_line),
            ter_status: ter_status
                .map(input::read::<u8>)
                .transpose()?
                .unwrap_or(d.ter_status),
            has_origx: has_origx
                .map(input::read::<bool>)
                .transpose()?
                .unwrap_or(d.has_origx),
            origx: origx
                .map(<ck::BioTransform as Field>::input)
                .transpose()?
                .unwrap_or(d.origx),
            info: info
                .map(<BTreeMap<String, String> as Field>::input)
                .transpose()?
                .unwrap_or(d.info),
            raw_remarks: raw_remarks
                .map(<Vec<String> as Field>::input)
                .transpose()?
                .unwrap_or(d.raw_remarks),
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn resolution(&self) -> f64 {
        self.0.resolution
    }
    #[wasm_bindgen(getter,js_name=conectMap)]
    pub fn conect_map(&self) -> Result<Ts<I32ListMap>, JsValue> {
        <BTreeMap<i32, Vec<i32>> as Field>::get(NativeField::owned(self.0.conect_map.clone()))
    }
    #[wasm_bindgen(getter,js_name=hasDFraction)]
    pub fn has_d_fraction(&self) -> bool {
        self.0.has_d_fraction
    }
    #[wasm_bindgen(getter,js_name=nonAsciiLine)]
    pub fn non_ascii_line(&self) -> i32 {
        self.0.non_ascii_line
    }
    #[wasm_bindgen(getter,js_name=terStatus)]
    pub fn ter_status(&self) -> u8 {
        self.0.ter_status
    }
    #[wasm_bindgen(getter,js_name=hasOrigx)]
    pub fn has_origx(&self) -> bool {
        self.0.has_origx
    }
    #[wasm_bindgen(getter)]
    pub fn origx(&self) -> BioTransform {
        BioTransform(self.0.origx)
    }
    #[wasm_bindgen(getter)]
    pub fn info(&self) -> Result<Ts<TextMap>, JsValue> {
        <BTreeMap<String, String> as Field>::get(NativeField::owned(self.0.info.clone()))
    }
    #[wasm_bindgen(getter,js_name=rawRemarks)]
    pub fn raw_remarks(&self) -> Result<Ts<TextArray>, JsValue> {
        <Vec<String> as Field>::get(NativeField::owned(self.0.raw_remarks.clone()))
    }
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(getter,js_name=sourceState)]
    pub fn source_state(&self) -> BioStructureSourceState {
        BioStructureSourceState(self.0.source_state().clone())
    }
    #[wasm_bindgen(getter,js_name=ncsOperIdentityId)]
    pub fn ncs_oper_identity_id(&self) -> Option<String> {
        self.0.ncs_oper_identity_id().map(str::to_owned)
    }
}
