//! Generated facade tables, in native order, without a second declaration list.
use crate::{input, operation_errors::MoleculeOpSpec, storage::NativeField};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub struct FeatureSpec(&'static ck::FeatureSpec);
#[wasm_bindgen]
impl FeatureSpec {
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name.into()
    }
    #[wasm_bindgen(getter)]
    pub fn category(&self) -> String {
        self.0.category.into()
    }
    #[wasm_bindgen(getter)]
    pub fn docs(&self) -> String {
        self.0.docs.into()
    }
}
#[wasm_bindgen]
pub struct FeatureSpecIter(NativeField<ck::FeatureSpecIter>);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "IteratorResult<FeatureSpec, undefined>")]
    pub type FeatureIteration;
}
#[wasm_bindgen]
impl FeatureSpecIter {
    #[wasm_bindgen(js_name = "[Symbol.iterator]")]
    pub fn iterator(&self) -> FeatureSpecIter {
        Self(self.0.clone())
    }
    #[wasm_bindgen]
    pub fn next(&self) -> Result<FeatureIteration, JsValue> {
        let value = self
            .0
            .write(|v| v.next())?
            .map(|v| JsValue::from(FeatureSpec(v)));
        Ok(crate::iterator::result(value)?.unchecked_into())
    }
}
#[wasm_bindgen(js_name=featureSpecs)]
pub fn feature_specs() -> FeatureSpecIter {
    FeatureSpecIter(NativeField::owned(ck::feature_specs()))
}
#[wasm_bindgen(js_name=featureSpec)]
pub fn feature_spec(name: Ts<input::TextInput>) -> Result<Option<FeatureSpec>, JsValue> {
    Ok(ck::feature_spec(&input::read::<String>(name)?).map(FeatureSpec))
}
#[wasm_bindgen]
pub struct SupportMatrixEntry(&'static ck::SupportMatrixEntry);
#[wasm_bindgen]
impl SupportMatrixEntry {
    #[wasm_bindgen(getter)]
    pub fn feature(&self) -> FeatureSpec {
        FeatureSpec(self.0.feature)
    }
    #[wasm_bindgen(getter)]
    pub fn operation(&self) -> Option<MoleculeOpSpec> {
        self.0.operation.map(MoleculeOpSpec)
    }
}
#[wasm_bindgen]
pub struct OperationInvariantEntry(&'static ck::OperationInvariantEntry);
#[wasm_bindgen]
impl OperationInvariantEntry {
    #[wasm_bindgen(getter)]
    pub fn operation(&self) -> MoleculeOpSpec {
        MoleculeOpSpec(self.0.operation)
    }
    #[wasm_bindgen(getter)]
    pub fn profile(&self) -> String {
        self.0.profile.into()
    }
}
#[wasm_bindgen]
pub struct ParityMatrixEntry(&'static ck::ParityMatrixEntry);
#[wasm_bindgen]
impl ParityMatrixEntry {
    #[wasm_bindgen(getter)]
    pub fn operation(&self) -> MoleculeOpSpec {
        MoleculeOpSpec(self.0.operation)
    }
    #[wasm_bindgen(getter)]
    pub fn feature(&self) -> FeatureSpec {
        FeatureSpec(self.0.feature)
    }
    #[wasm_bindgen(getter)]
    pub fn profile(&self) -> String {
        self.0.profile.into()
    }
    #[wasm_bindgen(getter,js_name=rdkitVersion)]
    pub fn rdkit_version(&self) -> Option<String> {
        self.0.rdkit_version.map(Into::into)
    }
}
#[wasm_bindgen(js_name=operationSpecs)]
pub fn operation_specs() -> Vec<MoleculeOpSpec> {
    ck::operation_specs()
        .iter()
        .map(|&v| MoleculeOpSpec(v))
        .collect()
}
#[wasm_bindgen(js_name=operationSpec)]
pub fn operation_spec(name: Ts<input::TextInput>) -> Result<Option<MoleculeOpSpec>, JsValue> {
    Ok(ck::operation_spec(&input::read::<String>(name)?).map(MoleculeOpSpec))
}
#[wasm_bindgen(js_name=supportMatrix)]
pub fn support_matrix() -> Vec<SupportMatrixEntry> {
    ck::support_matrix()
        .iter()
        .map(SupportMatrixEntry)
        .collect()
}
#[wasm_bindgen(js_name=operationInvariantMatrix)]
pub fn operation_invariant_matrix() -> Vec<OperationInvariantEntry> {
    ck::operation_invariant_matrix()
        .iter()
        .map(OperationInvariantEntry)
        .collect()
}
#[wasm_bindgen(js_name=parityMatrix)]
pub fn parity_matrix() -> Vec<ParityMatrixEntry> {
    ck::parity_matrix().iter().map(ParityMatrixEntry).collect()
}
#[wasm_bindgen(js_name=operationInvariant)]
pub fn operation_invariant(
    name: Ts<input::TextInput>,
) -> Result<Option<OperationInvariantEntry>, JsValue> {
    Ok(ck::operation_invariant(&input::read::<String>(name)?).map(OperationInvariantEntry))
}
#[wasm_bindgen(js_name=operationParity)]
pub fn operation_parity(name: Ts<input::TextInput>) -> Result<Option<ParityMatrixEntry>, JsValue> {
    Ok(ck::operation_parity(&input::read::<String>(name)?).map(ParityMatrixEntry))
}
