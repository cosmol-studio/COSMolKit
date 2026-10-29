//! Explicit detached parts conversion through the public BIO facade.
use crate::bio_readers::BioStructure;
use cosmolkit_wasm::rust as ck;
use std::rc::Rc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct BioStructureParts {
    inner: ck::BioStructureParts,
}
#[wasm_bindgen]
impl BioStructureParts {
    #[wasm_bindgen(getter, js_name=inputFormat, unchecked_return_type="BioCoordinateFormat")]
    pub fn input_format(&self) -> u32 {
        self.inner.input_format as u32
    }
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=fromParts)]
    pub fn from_parts(parts: &BioStructureParts) -> Result<BioStructure, JsValue> {
        // COSMolKit❗✔️: ck::BioStructure::from_parts(parts.inner.clone())
        // Detached input is retained, matching Python's public value contract.
        // Only explicit detached conversion copies payloads; traversal shares the owner.
        ck::BioStructure::from_parts(parts.inner.clone())
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_hierarchy_errors::structure_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=validateParts)]
    pub fn validate_parts(parts: &BioStructureParts) -> Result<(), JsValue> {
        // COSMolKit❗✔️: ck::BioStructure::validate_parts(&parts.inner)
        ck::BioStructure::validate_parts(&parts.inner)
            .map_err(|e| crate::bio_hierarchy_errors::structure_error(&e).unwrap_or_else(|e| e))
    }
    pub fn validate(&self) -> Result<(), JsValue> {
        self.inner
            .validate()
            .map_err(|e| crate::bio_hierarchy_errors::structure_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=intoParts)]
    pub fn into_parts(&self) -> BioStructureParts {
        // COSMolKit❗✔️: inner: self.inner.as_ref().clone().into_parts(),
        // The source binding retains its receiver and materializes the explicit
        // detached result through the sole public facade, including all blocks.
        BioStructureParts {
            inner: self.inner.as_ref().clone().into_parts(),
        }
    }
}
