//! Read-only projections of canonical generated declarations, not runtime authority.
use crate::alignment_values::set;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct FeatureSpec {
    inner: &'static ck::FeatureSpec,
}
#[wasm_bindgen]
impl FeatureSpec {
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.inner.name.into()
    }
    #[wasm_bindgen(getter)]
    pub fn category(&self) -> String {
        self.inner.category.into()
    }
    #[wasm_bindgen(getter)]
    pub fn docs(&self) -> String {
        self.inner.docs.into()
    }
}
#[wasm_bindgen]
pub struct FunctionStatus {
    inner: ck::FunctionStatus,
}
#[wasm_bindgen]
impl FunctionStatus {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        match self.inner {
            ck::FunctionStatus::Parity { .. } => "Parity",
            ck::FunctionStatus::ParityWithDifferences { .. } => "ParityWithDifferences",
            ck::FunctionStatus::Native => "Native",
            ck::FunctionStatus::Experimental => "Experimental",
        }
        .into()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn reference(&self) -> JsValue {
        match self.inner {
            ck::FunctionStatus::Parity { reference }
            | ck::FunctionStatus::ParityWithDifferences { reference, .. } => reference.into(),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn explanation(&self) -> JsValue {
        match self.inner {
            ck::FunctionStatus::ParityWithDifferences { explanation, .. } => explanation.into(),
            _ => JsValue::NULL,
        }
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum ParityPolicy {
    NotApplicable,
    RequiredWhenSupported,
    RequiredNow,
}
#[wasm_bindgen]
pub struct MoleculeOpSpec {
    inner: &'static ck::MoleculeOpSpec,
}
#[wasm_bindgen]
impl MoleculeOpSpec {
    #[wasm_bindgen(getter)]
    pub fn method(&self) -> String {
        self.inner.method.into()
    }
    #[wasm_bindgen(getter, js_name=implFn)]
    pub fn impl_fn(&self) -> String {
        self.inner.impl_fn.into()
    }
    #[wasm_bindgen(getter, js_name=resultType)]
    pub fn result_type(&self) -> String {
        self.inner.result_type.into()
    }
    #[wasm_bindgen(getter)]
    pub fn output(&self) -> String {
        format!("{:?}", self.inner.output)
    }
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        format!("{:?}", self.inner.domain)
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        format!("{:?}", self.inner.kind)
    }
    #[wasm_bindgen(getter, js_name=topologyEdit)]
    pub fn topology_edit(&self) -> String {
        format!("{:?}", self.inner.topology_edit)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "{ read: number; write: number }")]
    pub fn access(&self) -> Result<JsValue, JsValue> {
        let value = js_sys::Object::new();
        set(&value, "read", self.inner.access.read().bits().into())?;
        set(&value, "write", self.inner.access.write().bits().into())?;
        Ok(value.into())
    }
    #[wasm_bindgen(getter, js_name=mayMutate)]
    pub fn may_mutate(&self) -> u8 {
        self.inner.may_mutate.bits()
    }
    #[wasm_bindgen(getter, js_name=autoRemap)]
    pub fn auto_remap(&self) -> u8 {
        self.inner.auto_remap.bits()
    }
    #[wasm_bindgen(getter, js_name=derivedEffects, unchecked_return_type="{ recompute: number; preserve: number; invalidate: number; operationDefined: number }")]
    pub fn derived_effects(&self) -> Result<JsValue, JsValue> {
        let value = js_sys::Object::new();
        let effects = self.inner.derived_effects;
        set(&value, "recompute", effects.recompute.bits().into())?;
        set(&value, "preserve", effects.preserve.bits().into())?;
        set(&value, "invalidate", effects.invalidate.bits().into())?;
        set(
            &value,
            "operationDefined",
            effects.operation_defined.bits().into(),
        )?;
        Ok(value.into())
    }
    #[wasm_bindgen(getter, js_name=cipState)]
    pub fn cip_state(&self) -> String {
        format!("{:?}", self.inner.cip_state)
    }
    #[wasm_bindgen(getter, js_name=semanticPreconditions)]
    pub fn semantic_preconditions(&self) -> u8 {
        self.inner.semantic_preconditions.bits()
    }
    #[wasm_bindgen(getter, js_name=requiresMapping)]
    pub fn requires_mapping(&self) -> String {
        format!("{:?}", self.inner.requires_mapping)
    }
    #[wasm_bindgen(getter)]
    pub fn status(&self) -> FunctionStatus {
        FunctionStatus {
            inner: self.inner.status,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn parity(&self) -> ParityPolicy {
        match self.inner.parity {
            ck::ParityPolicy::NotApplicable => ParityPolicy::NotApplicable,
            ck::ParityPolicy::RequiredWhenSupported => ParityPolicy::RequiredWhenSupported,
            ck::ParityPolicy::RequiredNow => ParityPolicy::RequiredNow,
        }
    }
    #[wasm_bindgen(getter, js_name=ioRoundtrip)]
    pub fn io_roundtrip(&self) -> bool {
        self.inner.io_roundtrip
    }
}
#[wasm_bindgen]
pub struct SupportMatrixEntry {
    inner: &'static ck::SupportMatrixEntry,
}
#[wasm_bindgen]
impl SupportMatrixEntry {
    #[wasm_bindgen(getter)]
    pub fn feature(&self) -> FeatureSpec {
        FeatureSpec {
            inner: self.inner.feature,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "MoleculeOpSpec | null")]
    pub fn operation(&self) -> JsValue {
        self.inner
            .operation
            .map_or(JsValue::NULL, |inner| MoleculeOpSpec { inner }.into())
    }
}
#[wasm_bindgen]
pub struct OperationInvariantEntry {
    inner: &'static ck::OperationInvariantEntry,
}
#[wasm_bindgen]
impl OperationInvariantEntry {
    #[wasm_bindgen(getter)]
    pub fn operation(&self) -> MoleculeOpSpec {
        MoleculeOpSpec {
            inner: self.inner.operation,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn profile(&self) -> String {
        self.inner.profile.into()
    }
}
#[wasm_bindgen]
pub struct ParityMatrixEntry {
    inner: &'static ck::ParityMatrixEntry,
}
#[wasm_bindgen]
impl ParityMatrixEntry {
    #[wasm_bindgen(getter)]
    pub fn operation(&self) -> MoleculeOpSpec {
        MoleculeOpSpec {
            inner: self.inner.operation,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn feature(&self) -> FeatureSpec {
        FeatureSpec {
            inner: self.inner.feature,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn profile(&self) -> String {
        self.inner.profile.into()
    }
    #[wasm_bindgen(getter, js_name=rdkitVersion, unchecked_return_type="string | null")]
    pub fn rdkit_version(&self) -> JsValue {
        self.inner.rdkit_version.map_or(JsValue::NULL, Into::into)
    }
}
#[wasm_bindgen]
pub struct FeatureSpecIter {
    inner: ck::FeatureSpecIter,
}
#[wasm_bindgen]
impl FeatureSpecIter {
    #[wasm_bindgen(unchecked_return_type = "FeatureSpec | null")]
    pub fn next(&mut self) -> JsValue {
        self.inner
            .next()
            .map_or(JsValue::NULL, |inner| FeatureSpec { inner }.into())
    }
}
#[wasm_bindgen(js_name=featureSpecs)]
pub fn feature_specs() -> FeatureSpecIter {
    FeatureSpecIter {
        inner: ck::feature_specs(),
    }
}
#[wasm_bindgen(js_name=featureSpec, unchecked_return_type="FeatureSpec | null")]
pub fn feature_spec(name: &str) -> JsValue {
    ck::feature_spec(name).map_or(JsValue::NULL, |inner| FeatureSpec { inner }.into())
}
#[wasm_bindgen(js_name=operationSpecs, unchecked_return_type="MoleculeOpSpec[]")]
pub fn operation_specs() -> JsValue {
    ck::operation_specs()
        .iter()
        .map(|&inner| JsValue::from(MoleculeOpSpec { inner }))
        .collect::<js_sys::Array>()
        .into()
}
#[wasm_bindgen(js_name=operationSpec, unchecked_return_type="MoleculeOpSpec | null")]
pub fn operation_spec(method: &str) -> JsValue {
    ck::operation_spec(method).map_or(JsValue::NULL, |inner| MoleculeOpSpec { inner }.into())
}
#[wasm_bindgen(js_name=supportMatrix, unchecked_return_type="SupportMatrixEntry[]")]
pub fn support_matrix() -> JsValue {
    ck::support_matrix()
        .iter()
        .map(|inner| JsValue::from(SupportMatrixEntry { inner }))
        .collect::<js_sys::Array>()
        .into()
}
#[wasm_bindgen(js_name=operationInvariantMatrix, unchecked_return_type="OperationInvariantEntry[]")]
pub fn operation_invariant_matrix() -> JsValue {
    ck::operation_invariant_matrix()
        .iter()
        .map(|inner| JsValue::from(OperationInvariantEntry { inner }))
        .collect::<js_sys::Array>()
        .into()
}
#[wasm_bindgen(js_name=parityMatrix, unchecked_return_type="ParityMatrixEntry[]")]
pub fn parity_matrix() -> JsValue {
    ck::parity_matrix()
        .iter()
        .map(|inner| JsValue::from(ParityMatrixEntry { inner }))
        .collect::<js_sys::Array>()
        .into()
}
#[wasm_bindgen(js_name=operationInvariant, unchecked_return_type="OperationInvariantEntry | null")]
pub fn operation_invariant(method: &str) -> JsValue {
    ck::operation_invariant(method).map_or(JsValue::NULL, |inner| {
        OperationInvariantEntry { inner }.into()
    })
}
#[wasm_bindgen(js_name=operationParity, unchecked_return_type="ParityMatrixEntry | null")]
pub fn operation_parity(method: &str) -> JsValue {
    ck::operation_parity(method).map_or(JsValue::NULL, |inner| ParityMatrixEntry { inner }.into())
}
