//! Native identifier generation and complete diagnostic output.
use crate::{
    Molecule,
    enums::{InchiDiagnosticLevel, InchiErrorKind},
    errors::{error_projection, exception},
    fields::Field,
    input,
    storage::NativeField,
};
use cosmolkit as ck;

#[wasm_bindgen]
impl InchiWriteParams {
    #[wasm_bindgen(js_name = setOptions)]
    pub fn set_options_method(&self, value: <String as Field>::Input) -> Result<(), JsValue> {
        self.set_options(value)
    }
}

#[wasm_bindgen]
impl InchiReadParams {
    #[wasm_bindgen(js_name = setSanitize)]
    pub fn set_sanitize_method(&self, value: <bool as Field>::Input) -> Result<(), JsValue> {
        self.set_sanitize(value)
    }
    #[wasm_bindgen(js_name = setRemoveHs)]
    pub fn set_remove_hs_method(&self, value: <bool as Field>::Input) -> Result<(), JsValue> {
        self.set_remove_hs(value)
    }
}

use tsify::Ts;
use wasm_bindgen::prelude::*;
error_projection!(InchiError, ck::InchiError, {
    #[wasm_bindgen(getter)]
    pub fn operation(&self) -> String {
        self.inner.operation.into()
    }
    #[wasm_bindgen(getter)]
    pub fn detail(&self) -> String {
        self.inner.detail.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> InchiErrorKind {
        self.inner.kind.into()
    }
});
#[wasm_bindgen]
impl InchiError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> InchiErrorKind {
        self.inner.kind.into()
    }
}
ck::__ck_configuration_InchiReadParams!(crate::project_parameters, ck::InchiReadParams, {}, {
    InchiReadParams
});
ck::__ck_configuration_InchiWriteParams!(crate::project_parameters, ck::InchiWriteParams, {}, {
    InchiWriteParams
});
#[wasm_bindgen]
pub struct InchiDiagnostic(pub(crate) ck::InchiDiagnostic);
crate::fields::value_field!(
    ck::InchiDiagnostic,
    InchiDiagnostic,
    InchiDiagnosticInput,
    "InchiDiagnostic"
);
crate::fields::value_array!(
    ck::InchiDiagnostic,
    InchiDiagnostics,
    "ReadonlyArray<InchiDiagnostic>"
);
#[wasm_bindgen]
impl InchiDiagnostic {
    #[wasm_bindgen(getter)]
    pub fn level(&self) -> InchiDiagnosticLevel {
        self.0.level().into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.0.message().into()
    }
}
#[wasm_bindgen]
pub struct InchiResult(ck::InchiResult);
#[wasm_bindgen]
impl InchiResult {
    #[wasm_bindgen(getter)]
    pub fn inchi(&self) -> String {
        self.0.inchi().into()
    }
    #[wasm_bindgen(getter,js_name=returnCode)]
    pub fn return_code(&self) -> i32 {
        self.0.return_code()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.0.message().into()
    }
    #[wasm_bindgen(getter)]
    pub fn log(&self) -> String {
        self.0.log().into()
    }
    #[wasm_bindgen(getter,js_name=auxInfo)]
    pub fn aux_info(&self) -> String {
        self.0.aux_info().into()
    }
    #[wasm_bindgen(getter)]
    pub fn diagnostics(&self) -> Result<InchiDiagnostics, JsValue> {
        Vec::<ck::InchiDiagnostic>::get(NativeField::owned(self.0.diagnostics().to_vec()))
    }
}
#[wasm_bindgen(js_name=inchiToKey)]
pub fn inchi_to_key(text: Ts<input::TextInput>) -> Result<String, JsValue> {
    ck::inchi_to_key(&input::read::<String>(text)?).map_err(InchiError::from_rust)
}
#[wasm_bindgen(js_name=molFromInchi)]
pub fn mol_from_inchi(
    text: Ts<input::TextInput>,
    params: Option<InchiReadParamsInput>,
) -> Result<Molecule, JsValue> {
    Molecule::from_inchi(text, params)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fromInchi)]
    pub fn from_inchi(
        text: Ts<input::TextInput>,
        params: Option<InchiReadParamsInput>,
    ) -> Result<Self, JsValue> {
        let text = input::read::<String>(text)?;
        InchiReadParams::with_input(params.as_ref(), |p| {
            ck::Molecule::from_inchi_with_params(&text, p)
                .map(Self)
                .map_err(InchiError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toInchi)]
    pub fn to_inchi(&self, params: Option<InchiWriteParamsInput>) -> Result<String, JsValue> {
        InchiWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_inchi_with_params(p)
                .map_err(InchiError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toInchiResult)]
    pub fn to_inchi_result(
        &self,
        params: Option<InchiWriteParamsInput>,
    ) -> Result<InchiResult, JsValue> {
        InchiWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_inchi_result_with_params(p)
                .map(InchiResult)
                .map_err(InchiError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toInchiKey)]
    pub fn to_inchi_key(&self, params: Option<InchiWriteParamsInput>) -> Result<String, JsValue> {
        InchiWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_inchi_key_with_params(p)
                .map_err(InchiError::from_rust)
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = toInchiResultWithParams)]
    pub fn to_inchi_result_with_params(
        &self,
        params: InchiWriteParamsInput,
    ) -> Result<InchiResult, JsValue> {
        self.to_inchi_result(Some(params))
    }
    #[wasm_bindgen(js_name = fromInchiWithParams)]
    pub fn from_inchi_with_params(
        text: Ts<input::TextInput>,
        params: InchiReadParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_inchi(text, Some(params))
    }
    #[wasm_bindgen(js_name = toInchiWithParams)]
    pub fn to_inchi_with_params(&self, params: InchiWriteParamsInput) -> Result<String, JsValue> {
        self.to_inchi(Some(params))
    }
    #[wasm_bindgen(js_name = toInchiKeyWithParams)]
    pub fn to_inchi_key_with_params(
        &self,
        params: InchiWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_inchi_key(Some(params))
    }
}
