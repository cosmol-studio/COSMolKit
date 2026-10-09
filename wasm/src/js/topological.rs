//! Canonical topological fingerprint configuration, provenance and calls.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use crate::atom_pair_parameters::{optional_array, optional_u32_array};
use crate::fingerprint_values::Fingerprint;
use crate::host_values::{bool_value, type_error, u32_value};
use crate::query_values::QueryGraph;
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Map};
use std::collections::BTreeMap;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct TopologicalFingerprintParams {
    inner: ck::TopologicalFingerprintParams,
}
#[wasm_bindgen]
impl TopologicalFingerprintParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] min_path: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_path: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] fp_size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_bits_per_feature: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_hs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] target_density: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] min_size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] branched_paths: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_bond_order: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")]from_atoms:JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::TopologicalFingerprintParams {
        let mut inner = ck::TopologicalFingerprintParams::default();
        if !min_path.is_undefined() {
            inner.min_path = u32_value(&min_path, "minPath")?;
        }
        if !max_path.is_undefined() {
            inner.max_path = u32_value(&max_path, "maxPath")?;
        }
        if !fp_size.is_undefined() {
            inner.fp_size = u32_value(&fp_size, "fpSize")?;
        }
        if !num_bits_per_feature.is_undefined() {
            inner.num_bits_per_feature = u32_value(&num_bits_per_feature, "numBitsPerFeature")?;
        }
        if !use_hs.is_undefined() {
            inner.use_hs = bool_value(&use_hs, "useHs")?;
        }
        if !target_density.is_undefined() {
            inner.target_density = target_density
                .as_f64()
                .ok_or_else(|| type_error("targetDensity"))?;
        }
        if !min_size.is_undefined() {
            inner.min_size = u32_value(&min_size, "minSize")?;
        }
        if !branched_paths.is_undefined() {
            inner.branched_paths = bool_value(&branched_paths, "branchedPaths")?;
        }
        if !use_bond_order.is_undefined() {
            inner.use_bond_order = bool_value(&use_bond_order, "useBondOrder")?;
        }
        inner.atom_invariants = optional_u32_array(&atom_invariants, "atomInvariants")?;
        inner.from_atoms = optional_u32_array(&from_atoms, "fromAtoms")?;
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=minPath)]
    pub fn min_path(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.min_path
        self.inner.min_path
    }
    #[wasm_bindgen(getter,js_name=maxPath)]
    pub fn max_path(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.max_path
        self.inner.max_path
    }
    #[wasm_bindgen(getter,js_name=fpSize)]
    pub fn fp_size(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.fp_size
        self.inner.fp_size
    }
    #[wasm_bindgen(getter,js_name=numBitsPerFeature)]
    pub fn num_bits_per_feature(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.num_bits_per_feature
        self.inner.num_bits_per_feature
    }
    #[wasm_bindgen(getter,js_name=useHs)]
    pub fn use_hs(&self) -> bool {
        // COSMolKit❗✔️: self.inner.use_hs
        self.inner.use_hs
    }
    #[wasm_bindgen(getter,js_name=targetDensity)]
    pub fn target_density(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.target_density
        self.inner.target_density
    }
    #[wasm_bindgen(getter,js_name=minSize)]
    pub fn min_size(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.min_size
        self.inner.min_size
    }
    #[wasm_bindgen(getter,js_name=branchedPaths)]
    pub fn branched_paths(&self) -> bool {
        // COSMolKit❗✔️: self.inner.branched_paths
        self.inner.branched_paths
    }
    #[wasm_bindgen(getter,js_name=useBondOrder)]
    pub fn use_bond_order(&self) -> bool {
        // COSMolKit❗✔️: self.inner.use_bond_order
        self.inner.use_bond_order
    }
    #[wasm_bindgen(getter,js_name=atomInvariants,unchecked_return_type="number[] | null")]
    pub fn atom_invariants(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.atom_invariants
        optional_array(self.inner.atom_invariants.as_deref())
    }
    #[wasm_bindgen(getter,js_name=fromAtoms,unchecked_return_type="number[] | null")]
    pub fn from_atoms(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.from_atoms
        optional_array(self.inner.from_atoms.as_deref())
    }
}
#[wasm_bindgen]
pub struct TopologicalFingerprintOutputRequest {
    inner: ck::TopologicalFingerprintOutputRequest,
}
#[wasm_bindgen]
impl TopologicalFingerprintOutputRequest {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] atom_bits: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] bit_info: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::TopologicalFingerprintOutputRequest {
        let mut inner = ck::TopologicalFingerprintOutputRequest::default();
        if !atom_bits.is_undefined() {
            inner.atom_bits = bool_value(&atom_bits, "atomBits")?;
        }
        if !bit_info.is_undefined() {
            inner.bit_info = bool_value(&bit_info, "bitInfo")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=atomBits)]
    pub fn atom_bits(&self) -> bool {
        // COSMolKit❗✔️: self.inner.atom_bits
        self.inner.atom_bits
    }
    #[wasm_bindgen(getter,js_name=bitInfo)]
    pub fn bit_info(&self) -> bool {
        // COSMolKit❗✔️: self.inner.bit_info
        self.inner.bit_info
    }
}
fn atom_rows(rows: &[Vec<u32>]) -> Array {
    rows.iter()
        .map(|r| r.iter().copied().map(JsValue::from).collect::<Array>())
        .collect()
}
fn info_rows(rows: &BTreeMap<u32, Vec<Vec<i32>>>) -> Map {
    let map = Map::new();
    for (k, paths) in rows {
        let p: Array = paths
            .iter()
            .map(|p| p.iter().copied().map(JsValue::from).collect::<Array>())
            .collect();
        map.set(&JsValue::from(*k), &p);
    }
    map
}
#[wasm_bindgen]
pub struct TopologicalFingerprintOutput {
    inner: ck::TopologicalFingerprintOutput,
}
#[wasm_bindgen]
impl TopologicalFingerprintOutput {
    #[wasm_bindgen(getter,js_name=atomBits,unchecked_return_type="number[][] | null")]
    pub fn atom_bits(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.atom_bits.clone()
        self.inner
            .atom_bits
            .as_deref()
            .map_or(JsValue::NULL, |v| atom_rows(v).into())
    }
    #[wasm_bindgen(getter,js_name=bitInfo,unchecked_return_type="Map<number, number[][]> | null")]
    pub fn bit_info(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.bit_info.clone()
        self.inner
            .bit_info
            .as_ref()
            .map_or(JsValue::NULL, |v| info_rows(v).into())
    }
}
#[wasm_bindgen]
pub struct TopologicalFingerprintResult {
    inner: ck::TopologicalFingerprintResult,
}
#[wasm_bindgen]
impl TopologicalFingerprintResult {
    pub fn fingerprint(&self) -> Fingerprint {
        // COSMolKit❗✔️: inner: self.inner.fingerprint().clone(),
        Fingerprint {
            inner: self.inner.fingerprint().clone(),
        }
    }
    #[wasm_bindgen(js_name=atomBits,unchecked_return_type="number[][]")]
    pub fn atom_bits(&self) -> Result<Array, JsValue> {
        // COSMolKit❗✔️: .atom_bits()
        self.inner.atom_bits().map(atom_rows).map_err(error)
    }
    #[wasm_bindgen(js_name=bitInfo,unchecked_return_type="Map<number, number[][]>")]
    pub fn bit_info(&self) -> Result<Map, JsValue> {
        // COSMolKit❗✔️: .bit_info()
        self.inner.bit_info().map(info_rows).map_err(error)
    }
}
#[wasm_bindgen]
pub struct TopologicalFingerprintError {
    kind: String,
    message: String,
    cause: JsValue,
    reason: Option<String>,
    field: Option<String>,
}
#[wasm_bindgen]
impl TopologicalFingerprintError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "Fingerprint".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn reason(&self) -> JsValue {
        self.reason
            .as_deref()
            .map_or(JsValue::NULL, JsValue::from_str)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn field(&self) -> JsValue {
        self.field
            .as_deref()
            .map_or(JsValue::NULL, JsValue::from_str)
    }
}
pub(crate) fn topological_error(
    source: &ck::TopologicalFingerprintError,
) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: object.setattr("domain", "Fingerprint")?;
    use ck::TopologicalFingerprintError as E;
    let kind = match source {
        E::InvalidArguments { .. } => "InvalidArguments",
        E::OutputNotRequested { .. } => "OutputNotRequested",
        E::Topology(_) => "Topology",
        E::Query(_) => "Query",
        E::Paths(_) => "Paths",
        E::Value(_) => "Value",
    };
    let (reason, field) = match source {
        E::InvalidArguments { reason } => (Some((*reason).to_owned()), None),
        E::OutputNotRequested { field } => (None, Some((*field).to_owned())),
        _ => (None, None),
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("TopologicalFingerprintError");
    let e: JsValue = e.into();
    set(&e, "domain", "Fingerprint".into())?;
    set(&e, "kind", kind.into())?;
    if let Some(v) = &reason {
        set(&e, "reason", v.as_str().into())?;
    }
    if let Some(v) = &field {
        set(&e, "field", v.as_str().into())?;
    }
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        TopologicalFingerprintError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            reason,
            field,
        }
        .into(),
    )?;
    Ok(e)
}
fn error(e: ck::TopologicalFingerprintError) -> JsValue {
    topological_error(&e).unwrap_or_else(|e| e)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintTopological)]
    pub fn fingerprint_topological(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_fingerprint(
        self.inner
            .fingerprint_topological()
            .map(|inner| Fingerprint { inner })
            .map_err(error)
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintTopologicalWithParams)]
    pub fn fingerprint_topological_with_params(
        &self,
        params: &TopologicalFingerprintParams,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_fingerprint_with_params(
        self.inner
            .fingerprint_topological_with_params(&params.inner)
            .map(|inner| Fingerprint { inner })
            .map_err(error)
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintTopologicalWithOutput)]
    pub fn fingerprint_topological_with_output(
        &self,
    ) -> Result<TopologicalFingerprintResult, JsValue> {
        // COSMolKit❗✔️: .topological_fingerprint_with_output(
        self.inner
            .fingerprint_topological_with_output()
            .map(|inner| TopologicalFingerprintResult { inner })
            .map_err(error)
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintTopologicalWithOutputWithParams)]
    pub fn fingerprint_topological_with_output_with_params(
        &self,
        params: &TopologicalFingerprintParams,
        request: &TopologicalFingerprintOutputRequest,
    ) -> Result<TopologicalFingerprintResult, JsValue> {
        // COSMolKit❗✔️: .topological_fingerprint_with_output_with_params(
        self.inner
            .fingerprint_topological_with_output_with_params(&params.inner, request.inner)
            .map(|inner| TopologicalFingerprintResult { inner })
            .map_err(error)
    }
}
#[wasm_bindgen(js_name=fingerprintTopologicalQueryWithParams)]
pub fn fingerprint_topological_query_with_params(
    query: &QueryGraph,
    params: &TopologicalFingerprintParams,
) -> Result<Fingerprint, JsValue> {
    // COSMolKit❗✔️: ck::topological_query_fingerprint_with_params(
    ck::fingerprint_topological_query_with_params(&query.inner, &params.inner)
        .map(|inner| Fingerprint { inner })
        .map_err(error)
}
#[wasm_bindgen(js_name=fingerprintTopologicalQueryWithOutputWithParams)]
pub fn fingerprint_topological_query_with_output_with_params(
    query: &QueryGraph,
    params: &TopologicalFingerprintParams,
    request: &TopologicalFingerprintOutputRequest,
) -> Result<TopologicalFingerprintResult, JsValue> {
    // COSMolKit❗✔️: ck::topological_query_fingerprint_with_output_with_params(
    ck::fingerprint_topological_query_with_output_with_params(
        &query.inner,
        &params.inner,
        request.inner,
    )
    .map(|inner| TopologicalFingerprintResult { inner })
    .map_err(error)
}
