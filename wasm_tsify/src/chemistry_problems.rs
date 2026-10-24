//! Ordered native sanitization diagnostics, including their precise source graph.
use crate::{
    Molecule,
    errors::{enum_error, exception},
    fields::Field,
    sanitize::{SanitizeParams, SanitizeParamsInput},
    storage::NativeField,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | NonNullable<KekulizeError['cause']>"
    )]
    #[derive(Clone)]
    pub type ChemistryCause;
    #[wasm_bindgen(typescript_type = "ValenceError | KekulizeError")]
    pub type ChemistryPayload;
    #[wasm_bindgen(typescript_type = "ChemistryProblemError")]
    pub type ChemistryError;
}
enum_error!(ChemistryProblemError,ChemistryProblemErrorKind,ck::ChemistryProblemError,{Valence=>ck::ChemistryProblemError::Valence(_),Kekulize=>ck::ChemistryProblemError::Kekulize(_)},{
#[wasm_bindgen(getter)]pub fn source(&self)->ChemistryPayload{match self.inner.as_ref(){ck::ChemistryProblemError::Valence(v)=>crate::valence_errors::ValenceError::from_rust(v.clone()),ck::ChemistryProblemError::Kekulize(v)=>crate::kekulize_errors::KekulizeError::from_rust(v.clone())}.unchecked_into()}
},source=ChemistryCause,cause=|e:&ck::ChemistryProblemError|match e{ck::ChemistryProblemError::Valence(v)=>crate::valence_errors::cause(v).map(JsCast::unchecked_into),ck::ChemistryProblemError::Kekulize(v)=>crate::kekulize_errors::cause(v).map(JsCast::unchecked_into)});
#[wasm_bindgen]
pub struct ChemistryProblem(pub(crate) ck::ChemistryProblem);
crate::fields::value_field!(
    ck::ChemistryProblem,
    ChemistryProblem,
    ChemistryProblemInput,
    "ChemistryProblem"
);
crate::fields::value_array!(
    ck::ChemistryProblem,
    ChemistryProblems,
    "ReadonlyArray<ChemistryProblem>"
);
#[wasm_bindgen]
impl ChemistryProblem {
    #[wasm_bindgen(constructor)]
    pub fn new(
        operation: crate::enums::SanitizeStageInput,
        error: ChemistryError,
    ) -> Result<Self, JsValue> {
        let operation = ck::SanitizeStage::input(operation)?;
        ChemistryProblemError::with_native(error.as_ref(), |error| {
            Self(ck::ChemistryProblem {
                operation,
                error: error.clone(),
            })
        })
    }
    #[wasm_bindgen(getter)]
    pub fn operation(&self) -> crate::sanitize_errors::SanitizeStage {
        self.0.operation.into()
    }
    #[wasm_bindgen(getter)]
    pub fn error(&self) -> ChemistryError {
        ChemistryProblemError::from_rust(self.0.error.clone()).unchecked_into()
    }
}
#[wasm_bindgen]
pub struct ChemistryProblemReport(ck::ChemistryProblemReport);
#[wasm_bindgen]
impl ChemistryProblemReport {
    #[wasm_bindgen(constructor)]
    pub fn new(problems: Option<ChemistryProblems>) -> Result<Self, JsValue> {
        Ok(Self(match problems {
            Some(problems) => ck::ChemistryProblemReport {
                problems: Vec::<ck::ChemistryProblem>::input(problems)?,
            },
            None => ck::ChemistryProblemReport::default(),
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn problems(&self) -> Result<ChemistryProblems, JsValue> {
        Vec::<ck::ChemistryProblem>::get(NativeField::owned(self.0.problems.clone()))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=detectChemistryProblems)]
    pub fn detect_chemistry_problems(
        &self,
        params: Option<SanitizeParamsInput>,
    ) -> Result<ChemistryProblemReport, JsValue> {
        SanitizeParams::with_input(params.as_ref(), |p| {
            self.0
                .detect_chemistry_problems_with_params(p)
                .map(ChemistryProblemReport)
                .map_err(crate::sanitize_errors::SanitizeError::from_rust)
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = detectChemistryProblemsWithParams)]
    pub fn detect_chemistry_problems_with_params(
        &self,
        params: SanitizeParamsInput,
    ) -> Result<ChemistryProblemReport, JsValue> {
        self.detect_chemistry_problems(Some(params))
    }
}
