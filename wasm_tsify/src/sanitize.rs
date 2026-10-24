//! Native sanitization stage selection, not independently interpreted flags.
use crate::{Molecule, fields::Field, input, sanitize_errors::SanitizeError};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct SanitizeOperations(pub(crate) ck::SanitizeOperations);
crate::fields::value_field!(
    ck::SanitizeOperations,
    SanitizeOperations,
    SanitizeOperationsInput,
    "SanitizeOperations"
);
#[wasm_bindgen]
impl SanitizeOperations {
    #[wasm_bindgen(getter,js_name=NONE)]
    pub fn none() -> Self {
        Self(ck::SanitizeOperations::NONE)
    }
    #[wasm_bindgen(getter,js_name=CLEANUP)]
    pub fn cleanup() -> Self {
        Self(ck::SanitizeOperations::CLEANUP)
    }
    #[wasm_bindgen(getter,js_name=PROPERTIES)]
    pub fn properties() -> Self {
        Self(ck::SanitizeOperations::PROPERTIES)
    }
    #[wasm_bindgen(getter,js_name=SYMM_RINGS)]
    pub fn symm_rings() -> Self {
        Self(ck::SanitizeOperations::SYMM_RINGS)
    }
    #[wasm_bindgen(getter,js_name=KEKULIZE)]
    pub fn kekulize() -> Self {
        Self(ck::SanitizeOperations::KEKULIZE)
    }
    #[wasm_bindgen(getter,js_name=FIND_RADICALS)]
    pub fn find_radicals() -> Self {
        Self(ck::SanitizeOperations::FIND_RADICALS)
    }
    #[wasm_bindgen(getter,js_name=SET_AROMATICITY)]
    pub fn set_aromaticity() -> Self {
        Self(ck::SanitizeOperations::SET_AROMATICITY)
    }
    #[wasm_bindgen(getter,js_name=SET_CONJUGATION)]
    pub fn set_conjugation() -> Self {
        Self(ck::SanitizeOperations::SET_CONJUGATION)
    }
    #[wasm_bindgen(getter,js_name=SET_HYBRIDIZATION)]
    pub fn set_hybridization() -> Self {
        Self(ck::SanitizeOperations::SET_HYBRIDIZATION)
    }
    #[wasm_bindgen(getter,js_name=CLEANUP_CHIRALITY)]
    pub fn cleanup_chirality() -> Self {
        Self(ck::SanitizeOperations::CLEANUP_CHIRALITY)
    }
    #[wasm_bindgen(getter,js_name=ADJUST_HS)]
    pub fn adjust_hs() -> Self {
        Self(ck::SanitizeOperations::ADJUST_HS)
    }
    #[wasm_bindgen(getter,js_name=CLEANUP_ORGANOMETALLICS)]
    pub fn cleanup_organometallics() -> Self {
        Self(ck::SanitizeOperations::CLEANUP_ORGANOMETALLICS)
    }
    #[wasm_bindgen(getter,js_name=CLEANUP_ATROPISOMERS)]
    pub fn cleanup_atropisomers() -> Self {
        Self(ck::SanitizeOperations::CLEANUP_ATROPISOMERS)
    }
    #[wasm_bindgen(getter,js_name=ALL)]
    pub fn all() -> Self {
        Self(ck::SanitizeOperations::ALL)
    }
    #[wasm_bindgen(js_name=fromBits)]
    pub fn from_bits(bits: tsify::Ts<input::U32Input>) -> Result<Self, JsValue> {
        ck::SanitizeOperations::from_bits(input::read::<u32>(bits)?)
            .map(Self)
            .map_err(SanitizeError::from_rust)
    }
    #[wasm_bindgen(getter)]
    pub fn bits(&self) -> u32 {
        self.0.bits()
    }
    #[wasm_bindgen(getter,js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    pub fn contains(&self, other: SanitizeOperationsInput) -> Result<bool, JsValue> {
        Ok(self.0.contains(ck::SanitizeOperations::input(other)?))
    }
    pub fn union(&self, other: SanitizeOperationsInput) -> Result<Self, JsValue> {
        Ok(Self(self.0 | ck::SanitizeOperations::input(other)?))
    }
}
ck::__ck_configuration_SanitizeParams!(crate::project_parameters, ck::SanitizeParams, {}, {
    SanitizeParams
});
#[wasm_bindgen]
impl Molecule {
    pub fn sanitize(&self, params: Option<SanitizeParamsInput>) -> Result<Molecule, JsValue> {
        SanitizeParams::with_input(params.as_ref(), |params| {
            self.0
                .sanitize_with_params(params)
                .map(Molecule)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    pub fn sanitize_(&mut self, params: Option<SanitizeParamsInput>) -> Result<(), JsValue> {
        SanitizeParams::with_input(params.as_ref(), |params| {
            self.0
                .sanitize_with_params_(params)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = sanitizeWithParams)]
    pub fn sanitize_with_params(&self, params: SanitizeParamsInput) -> Result<Molecule, JsValue> {
        self.sanitize(Some(params))
    }
    #[wasm_bindgen(js_name = sanitizeWithParams_)]
    pub fn sanitize_with_params_(&mut self, params: SanitizeParamsInput) -> Result<(), JsValue> {
        self.sanitize_(Some(params))
    }
}
