//! Native tautomer scoring values; enumeration callbacks are a separate capability.
use crate::{Molecule, fields::Field, input, storage::NativeField};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct TautomerScoreTerm(pub(crate) ck::TautomerScoreTerm);
crate::fields::value_field!(
    ck::TautomerScoreTerm,
    TautomerScoreTerm,
    TautomerScoreTermInput,
    "TautomerScoreTerm"
);
crate::fields::value_array!(
    ck::TautomerScoreTerm,
    TautomerScoreTerms,
    "ReadonlyArray<TautomerScoreTerm>"
);
#[wasm_bindgen]
impl TautomerScoreTerm {
    #[wasm_bindgen(constructor)]
    pub fn new(
        name: Ts<input::TextInput>,
        smarts: Ts<input::TextInput>,
        score: Ts<input::I32Input>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::TautomerScoreTerm::new(
            input::read::<String>(name)?,
            input::read::<String>(smarts)?,
            input::read::<i32>(score)?,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name().into()
    }
    #[wasm_bindgen(getter)]
    pub fn smarts(&self) -> String {
        self.0.smarts().into()
    }
    #[wasm_bindgen(getter)]
    pub fn score(&self) -> i32 {
        self.0.score()
    }
}
#[wasm_bindgen]
pub struct TautomerScore(pub(crate) ck::TautomerScore);
#[wasm_bindgen]
impl TautomerScore {
    #[wasm_bindgen(getter)]
    pub fn ring(&self) -> i32 {
        self.0.ring()
    }
    #[wasm_bindgen(getter)]
    pub fn substructure(&self) -> i32 {
        self.0.substructure()
    }
    #[wasm_bindgen(getter,js_name=heteroHydrogen)]
    pub fn hetero_hydrogen(&self) -> i32 {
        self.0.hetero_hydrogen()
    }
    #[wasm_bindgen(getter)]
    pub fn total(&self) -> i32 {
        self.0.total()
    }
}
ck::__ck_configuration_TautomerScoreParams!(
    crate::project_parameters,
    ck::TautomerScoreParams,
    {},
    { TautomerScoreParams }
);
#[wasm_bindgen(js_name=defaultTautomerScoreTerms)]
pub fn default_tautomer_score_terms() -> Result<TautomerScoreTerms, JsValue> {
    Vec::<ck::TautomerScoreTerm>::get(NativeField::owned(
        ck::default_tautomer_score_terms().to_vec(),
    ))
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=tautomerScore)]
    pub fn tautomer_score(
        &self,
        params: Option<TautomerScoreParamsInput>,
    ) -> Result<TautomerScore, JsValue> {
        TautomerScoreParams::with_input(params.as_ref(), |p| {
            self.0
                .tautomer_score_with_params(p)
                .map(TautomerScore)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=tautomerScoreWithParams)]
    pub fn tautomer_score_with_params(
        &self,
        params: TautomerScoreParamsInput,
    ) -> Result<TautomerScore, JsValue> {
        self.tautomer_score(Some(params))
    }
}
