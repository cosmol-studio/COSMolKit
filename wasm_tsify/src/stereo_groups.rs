//! Detached enhanced-stereo membership delegates native validation.
use crate::{
    enums::{StereoGroupKind, StereoGroupKindInput},
    fields::Field,
    input,
    model_errors::StereoGroupError,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct StereoGroup(pub(crate) ck::StereoGroup);
crate::fields::value_field!(
    ck::StereoGroup,
    StereoGroup,
    StereoGroupInput,
    "StereoGroup"
);
crate::fields::value_array!(ck::StereoGroup, StereoGroups, "ReadonlyArray<StereoGroup>");
#[wasm_bindgen]
impl StereoGroup {
    #[wasm_bindgen(constructor)]
    pub fn new(
        kind: StereoGroupKindInput,
        atoms: <Vec<ck::AtomId> as Field>::Input,
        bonds: <Vec<ck::BondId> as Field>::Input,
    ) -> Result<Self, JsValue> {
        ck::StereoGroup::new(
            ck::StereoGroupKind::input(kind)?,
            Vec::<ck::AtomId>::input(atoms)?,
            Vec::<ck::BondId>::input(bonds)?,
        )
        .map(Self)
        .map_err(StereoGroupError::from_rust)
    }
    #[wasm_bindgen(js_name=withId)]
    pub fn with_id(&self, id: Ts<input::U32Input>) -> Result<Self, JsValue> {
        Ok(Self(self.0.clone().with_id(input::read(id)?)))
    }
    #[wasm_bindgen(js_name=withWriteId)]
    pub fn with_write_id(&self, id: Ts<input::U32Input>) -> Result<Self, JsValue> {
        Ok(Self(self.0.clone().with_write_id(input::read(id)?)))
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> Option<u32> {
        self.0.id()
    }
    #[wasm_bindgen(getter,js_name=writeId)]
    pub fn write_id(&self) -> u32 {
        self.0.write_id()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> StereoGroupKind {
        self.0.kind().into()
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Result<<Vec<ck::AtomId> as Field>::Output, JsValue> {
        Vec::<ck::AtomId>::get(NativeField::owned(self.0.atoms().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn bonds(&self) -> Result<<Vec<ck::BondId> as Field>::Output, JsValue> {
        Vec::<ck::BondId>::get(NativeField::owned(self.0.bonds().to_vec()))
    }
}
