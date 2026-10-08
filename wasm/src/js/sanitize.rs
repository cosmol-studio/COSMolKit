//! Canonical sanitization and chemistry problems with complete source payloads.
use crate::Molecule;
use crate::alignment_values::{operation_error, set};
use crate::host_values::*;
use crate::transform_errors::{kekulize_error, sanitize_error};
use crate::transform_parameters::{SanitizeParams, SanitizeStage, stage};
use cosmolkit_wasm::rust as ck;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen(
    inline_js = "export function visitSanitizeParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid SanitizeParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitSanitizeParams)]
    fn visit_params(v: &JsValue, f: &mut dyn FnMut(&SanitizeParams)) -> Result<(), JsValue>;
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=sanitize)]
    pub fn sanitize(&self) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.sanitize()
        self.inner
            .sanitize()
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=sanitizeWithParams)]
    pub fn sanitize_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "SanitizeParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.sanitize_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SanitizeParams| {
            result = Some(
                self.inner
                    .sanitize_with_params(&p.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=sanitize_)]
    pub fn sanitize_(&self) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.sanitize_()
        self.inner
            .sanitize_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=sanitizeWithParams_)]
    pub fn sanitize_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "SanitizeParams")] params: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.sanitize_with_params_(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SanitizeParams| {
            result = Some(
                self.inner
                    .sanitize_with_params_(&p.inner)
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=detectChemistryProblems)]
    pub fn detect_chemistry_problems(&self) -> Result<ChemistryProblemReport, JsValue> {
        // COSMolKit❗✔️: self.inner.detect_chemistry_problems()
        self.inner
            .detect_chemistry_problems()
            .map(|inner| ChemistryProblemReport { inner })
            .map_err(|e| sanitize_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=detectChemistryProblemsWithParams)]
    pub fn detect_chemistry_problems_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "SanitizeParams")] params: JsValue,
    ) -> Result<ChemistryProblemReport, JsValue> {
        // COSMolKit❗✔️: self.inner.detect_chemistry_problems_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SanitizeParams| {
            result = Some(
                self.inner
                    .detect_chemistry_problems_with_params(&p.inner)
                    .map(|inner| ChemistryProblemReport { inner })
                    .map_err(|e| sanitize_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
#[wasm_bindgen]
pub struct ChemistryProblemReport {
    inner: ck::ChemistryProblemReport,
}
#[wasm_bindgen]
impl ChemistryProblemReport {
    #[wasm_bindgen(getter)]
    pub fn problems(&self) -> Vec<ChemistryProblem> {
        self.inner
            .problems
            .iter()
            .cloned()
            .map(|inner| ChemistryProblem { inner })
            .collect()
    }
}
#[wasm_bindgen]
pub struct ChemistryProblem {
    inner: ck::ChemistryProblem,
}
#[wasm_bindgen]
impl ChemistryProblem {
    #[wasm_bindgen(getter)]
    pub fn operation(&self) -> SanitizeStage {
        stage(self.inner.operation)
    }
    #[wasm_bindgen(
        getter,
        unchecked_return_type = "Error & {readonly domain: string; readonly kind: string; readonly cause: Error; readonly detail: ChemistryProblemError}"
    )]
    pub fn error(&self) -> Result<JsValue, JsValue> {
        problem_error(&self.inner.error)
    }
}
#[wasm_bindgen]
pub struct ChemistryProblemError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ChemistryProblemError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "sanitize".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
pub(crate) fn problem_error(source: &ck::ChemistryProblemError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: let cause = match &self.inner.error {
    // COSMolKit❗✔️:     ck::ChemistryProblemError::Valence(source) => {
    // COSMolKit❗✔️:         crate::canonical_atom_bond::valence_pyerr(py, source.clone())
    // COSMolKit❗✔️:     }
    // COSMolKit❗✔️:     ck::ChemistryProblemError::Kekulize(source) => kekulize_pyerr(py, source.clone()),
    // COSMolKit❗✔️: };
    let (kind, cause) = match source {
        ck::ChemistryProblemError::Valence(e) => {
            ("Valence", crate::valence_errors::valence_error(e)?)
        }
        ck::ChemistryProblemError::Kekulize(e) => ("Kekulize", kekulize_error(e)?),
    };
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("ChemistryProblemError");
    set(&error, "domain", "sanitize".into())?;
    set(&error, "kind", kind.into())?;
    set(&error, "cause", cause.clone())?;
    set(
        &error,
        "detail",
        ChemistryProblemError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(error.into())
}
