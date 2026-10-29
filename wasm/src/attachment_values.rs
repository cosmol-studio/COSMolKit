//! Detached template attachment values use the native validation entrypoints.
use crate::{fields::Field, input, storage::NativeField};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct TemplateAttachment(ck::TemplateAttachment);
crate::fields::value_field!(
    ck::TemplateAttachment,
    TemplateAttachment,
    TemplateAttachmentInput,
    "TemplateAttachment"
);
crate::fields::value_array!(
    ck::TemplateAttachment,
    TemplateAttachments,
    "ReadonlyArray<TemplateAttachment>"
);
#[wasm_bindgen]
impl TemplateAttachment {
    #[wasm_bindgen(constructor)]
    pub fn new(
        target: Ts<input::IndexInput>,
        label: Ts<input::TextInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::TemplateAttachment::new(
            ck::AtomId::new(input::read::<usize>(target)?),
            input::read::<String>(label)?,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn target(&self) -> usize {
        self.0.target().index()
    }
    #[wasm_bindgen(getter)]
    pub fn label(&self) -> String {
        self.0.label().into()
    }
}
#[wasm_bindgen]
pub struct TemplateAttachmentOrder(ck::TemplateAttachmentOrder);
crate::fields::value_field!(
    ck::TemplateAttachmentOrder,
    TemplateAttachmentOrder,
    TemplateAttachmentOrderInput,
    "TemplateAttachmentOrder"
);
#[wasm_bindgen]
impl TemplateAttachmentOrder {
    #[wasm_bindgen(constructor)]
    pub fn new(entries: TemplateAttachments) -> Result<Self, JsValue> {
        ck::TemplateAttachmentOrder::new(Vec::<ck::TemplateAttachment>::input(entries)?)
            .map(Self)
            .map_err(crate::model_errors::TemplateAttachmentOrderError::from_rust)
    }
    #[wasm_bindgen(getter)]
    pub fn entries(&self) -> Result<TemplateAttachments, JsValue> {
        Vec::<ck::TemplateAttachment>::get(NativeField::owned(self.0.entries().to_vec()))
    }
    #[wasm_bindgen(js_name=validateForAtomCount)]
    pub fn validate_for_atom_count(&self, count: Ts<input::IndexInput>) -> Result<(), JsValue> {
        self.0
            .validate_for_atom_count(input::read::<usize>(count)?)
            .map_err(crate::model_errors::TemplateAttachmentOrderError::from_rust)
    }
}
