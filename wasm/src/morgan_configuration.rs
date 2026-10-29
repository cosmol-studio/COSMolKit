//! Morgan input policies and canonical flattened configuration metadata.
use crate::{fields::Field, storage::NativeField};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct MorganInvariants(pub(crate) ck::MorganInvariants);
crate::fields::value_field!(
    ck::MorganInvariants,
    MorganInvariants,
    MorganInvariantsInput,
    "MorganInvariants"
);
#[wasm_bindgen]
pub enum MorganInvariantsKind {
    Connectivity,
    Features,
    FeaturePatterns,
}
#[wasm_bindgen]
impl MorganInvariants {
    #[wasm_bindgen(getter,js_name=Connectivity)]
    pub fn connectivity() -> Self {
        Self(ck::MorganInvariants::Connectivity)
    }
    #[wasm_bindgen(getter,js_name=Features)]
    pub fn features() -> Self {
        Self(ck::MorganInvariants::Features)
    }
    #[wasm_bindgen(js_name=featurePatterns)]
    pub fn feature_patterns(
        patterns: crate::query_values::QueryGraphArray,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::MorganInvariants::FeaturePatterns(Vec::<
            ck::QueryGraph,
        >::input(
            patterns
        )?)))
    }
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> MorganInvariantsKind {
        match self.0 {
            ck::MorganInvariants::Connectivity => MorganInvariantsKind::Connectivity,
            ck::MorganInvariants::Features => MorganInvariantsKind::Features,
            ck::MorganInvariants::FeaturePatterns(_) => MorganInvariantsKind::FeaturePatterns,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn patterns(&self) -> Result<Option<crate::query_values::QueryGraphArray>, JsValue> {
        match &self.0 {
            ck::MorganInvariants::FeaturePatterns(v) => {
                <Vec<ck::QueryGraph> as Field>::get(NativeField::owned(v.clone())).map(Some)
            }
            _ => Ok(None),
        }
    }
}
ck::__ck_configuration_MorganGeneratorParams!(
    crate::project_parameters,
    ck::MorganGeneratorParams,
    {},
    { MorganGeneratorParams }
);
ck::__ck_configuration_MorganFingerprintParams!(
    crate::project_parameters,
    ck::MorganFingerprintParams,
    {},
    { MorganFingerprintParams }
);
ck::__ck_configuration_MorganCallParams!(crate::project_parameters, ck::MorganCallParams, {}, {
    MorganCallParams
});
#[wasm_bindgen]
impl MorganGeneratorParams {
    #[wasm_bindgen(js_name = infoString)]
    pub fn info_string(&self) -> Result<String, JsValue> {
        self.inner.read(ck::MorganGeneratorParams::info_string)
    }
    #[wasm_bindgen(js_name = toJson)]
    pub fn to_json(&self) -> Result<String, JsValue> {
        self.inner.read(ck::MorganGeneratorParams::to_json)
    }
    #[wasm_bindgen(js_name = withJson)]
    pub fn with_json(&self, json: tsify::Ts<crate::input::TextInput>) -> Result<Self, JsValue> {
        let json = crate::input::read::<String>(json)?;
        let native = self
            .inner
            .read(|value| value.with_json(&json))?
            .map_err(crate::morgan_errors::FingerprintJsonError::from_rust)?;
        Ok(Self {
            inner: crate::storage::NativeField::owned(native),
        })
    }
}
