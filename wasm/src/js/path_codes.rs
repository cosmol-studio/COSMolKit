//! Lossless public atom-code and path-score transport; no chemistry here.
use crate::Molecule;
use crate::alignment_values::{operation_error, set, source_error};
use crate::fingerprint_values::u64_value;
use crate::host_values::{bool_value, sequence, type_error, u32_value, usize_value};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::{error::Error as RustError, sync::Arc};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct AtomPairsParameters;
#[wasm_bindgen]
impl AtomPairsParameters {
    #[wasm_bindgen(js_name=version)]
    pub fn version() -> String {
        // COSMolKit❗✔️: ck::AtomPairsParameters::version()
        ck::AtomPairsParameters::version().to_owned()
    }
    #[wasm_bindgen(js_name=numTypeBits)]
    pub fn num_type_bits() -> u32 {
        // COSMolKit❗✔️: ck::AtomPairsParameters::num_type_bits()
        ck::AtomPairsParameters::num_type_bits()
    }
    #[wasm_bindgen(js_name=numPiBits)]
    pub fn num_pi_bits() -> u32 {
        // COSMolKit❗✔️: ck::AtomPairsParameters::num_pi_bits()
        ck::AtomPairsParameters::num_pi_bits()
    }
    #[wasm_bindgen(js_name=numBranchBits)]
    pub fn num_branch_bits() -> u32 {
        // COSMolKit❗✔️: ck::AtomPairsParameters::num_branch_bits()
        ck::AtomPairsParameters::num_branch_bits()
    }
    #[wasm_bindgen(js_name=numChiralBits)]
    pub fn num_chiral_bits() -> u32 {
        // COSMolKit❗✔️: ck::AtomPairsParameters::num_chiral_bits()
        ck::AtomPairsParameters::num_chiral_bits()
    }
    #[wasm_bindgen(js_name=codeSize)]
    pub fn code_size() -> u32 {
        // COSMolKit❗✔️: ck::AtomPairsParameters::code_size()
        ck::AtomPairsParameters::code_size()
    }
    #[wasm_bindgen(js_name=numPathBits)]
    pub fn num_path_bits() -> u32 {
        // COSMolKit❗✔️: ck::AtomPairsParameters::num_path_bits()
        ck::AtomPairsParameters::num_path_bits()
    }
    #[wasm_bindgen(js_name=maxPathLength)]
    pub fn max_path_length() -> u32 {
        // COSMolKit❗✔️: ck::AtomPairsParameters::max_path_length()
        ck::AtomPairsParameters::max_path_length()
    }
    #[wasm_bindgen(js_name=numAtomPairFingerprintBits)]
    pub fn num_atom_pair_fingerprint_bits() -> u32 {
        // COSMolKit❗✔️: ck::AtomPairsParameters::num_atom_pair_fingerprint_bits()
        ck::AtomPairsParameters::num_atom_pair_fingerprint_bits()
    }
    #[wasm_bindgen(js_name=atomTypes)]
    pub fn atom_types() -> Vec<u32> {
        // COSMolKit❗✔️: ck::AtomPairsParameters::atom_types()
        ck::AtomPairsParameters::atom_types()
    }
}
#[wasm_bindgen]
pub struct AtomCodeExplanation {
    inner: ck::AtomCodeExplanation,
}
#[wasm_bindgen]
impl AtomCodeExplanation {
    #[wasm_bindgen(js_name=fromCode)]
    pub fn from_code(
        #[wasm_bindgen(unchecked_param_type = "bigint")] code: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "bigint")] branch_subtract: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::AtomCodeExplanation::from_code(code, branch_subtract, include_chirality)
        let branch = if branch_subtract.is_undefined() {
            0
        } else {
            if !branch_subtract.is_bigint() {
                return Err(type_error("branchSubtract"));
            }
            i64::try_from(branch_subtract)
                .map_err(|_| js_sys::RangeError::new("branchSubtract is outside i64 range"))?
        };
        ck::AtomCodeExplanation::from_code(
            u64_value(&code, "code")?,
            branch,
            if include_chirality.is_undefined() {
                false
            } else {
                bool_value(&include_chirality, "includeChirality")?
            },
        )
        .map(|inner| Self { inner })
        .map_err(|e| explanation_error(&e).unwrap_or_else(|e| e))
    }
    pub fn symbol(&self) -> String {
        // COSMolKit❗✔️: self.inner.symbol()
        self.inner.symbol().into()
    }
    #[wasm_bindgen(js_name=branchCount)]
    pub fn branch_count(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.branch_count()
        self.inner.branch_count()
    }
    #[wasm_bindgen(js_name=piElectrons)]
    pub fn pi_electrons(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.pi_electrons()
        self.inner.pi_electrons()
    }
    #[wasm_bindgen(unchecked_return_type = "string | null")]
    pub fn chirality(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.chirality()
        self.inner.chirality().map_or(JsValue::NULL, JsValue::from)
    }
}
#[wasm_bindgen]
pub struct AtomCodeExplanationError {
    code: u8,
    message: String,
}
#[wasm_bindgen]
impl AtomCodeExplanationError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "Fingerprint".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        "UnknownChirality".into()
    }
    #[wasm_bindgen(getter)]
    pub fn code(&self) -> u8 {
        self.code
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "never | null")]
    pub fn cause(&self) -> JsValue {
        JsValue::NULL
    }
}
pub(crate) fn explanation_error(source: &ck::AtomCodeExplanationError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: let ck::AtomCodeExplanationError::UnknownChirality { code } = error;
    let ck::AtomCodeExplanationError::UnknownChirality { code } = *source;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("AtomCodeExplanationError");
    let e: JsValue = e.into();
    set(&e, "domain", "Fingerprint".into())?;
    set(&e, "kind", "UnknownChirality".into())?;
    set(&e, "code", code.into())?;
    set(
        &e,
        "detail",
        AtomCodeExplanationError {
            code,
            message: source.to_string(),
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
pub struct TopologicalTorsionPathScoreError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: JsValue,
}
#[wasm_bindgen]
impl TopologicalTorsionPathScoreError {
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
    #[wasm_bindgen(getter,js_name=actual,unchecked_return_type="number | null")]
    pub fn actual(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"actual".into())
    }
    #[wasm_bindgen(getter,js_name=required,unchecked_return_type="number | null")]
    pub fn required(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"required".into())
    }
    #[wasm_bindgen(getter,js_name=index,unchecked_return_type="number | null")]
    pub fn index(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"index".into())
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn atom_count(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"atomCount".into())
    }
    #[wasm_bindgen(getter,js_name=code,unchecked_return_type="number | null")]
    pub fn code(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"code".into())
    }
    #[wasm_bindgen(getter,js_name=subtract,unchecked_return_type="number | null")]
    pub fn subtract(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"subtract".into())
    }
}
pub(crate) fn score_error(
    source: &ck::TopologicalTorsionPathScoreError,
) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: E::AtomCodeUnderflow { .. } => "AtomCodeUnderflow",
    use ck::TopologicalTorsionPathScoreError as E;
    let kind = match source {
        E::ZeroSize => "ZeroSize",
        E::ShortPath { .. } => "ShortPath",
        E::ShortAtomCodes { .. } => "ShortAtomCodes",
        E::AtomIndexOutOfRange { .. } => "AtomIndexOutOfRange",
        E::AtomCodeUnderflow { .. } => "AtomCodeUnderflow",
        E::InvalidTopology(_) => "InvalidTopology",
        E::AtomCode(_) => "AtomCode",
        E::PackedCode(_) => "PackedCode",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("TopologicalTorsionPathScoreError");
    let e: JsValue = e.into();
    let fields: JsValue = js_sys::Object::new().into();
    for f in [
        "actual",
        "required",
        "index",
        "atomCount",
        "code",
        "subtract",
    ] {
        set(&fields, f, JsValue::NULL)?;
    }
    let mut field = |name: &str, v: JsValue| -> Result<(), JsValue> {
        set(&e, name, v.clone())?;
        set(&fields, name, v)
    };
    match *source {
        E::ShortPath { actual, required } | E::ShortAtomCodes { actual, required } => {
            field("actual", actual.into())?;
            field("required", required.into())?;
        }
        E::AtomIndexOutOfRange { index, atom_count } => {
            field("index", index.into())?;
            field("atomCount", atom_count.into())?;
        }
        E::AtomCodeUnderflow {
            index,
            code,
            subtract,
        } => {
            field("index", index.into())?;
            field("code", code.into())?;
            field("subtract", subtract.into())?;
        }
        _ => {}
    }
    set(&e, "domain", "Fingerprint".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        TopologicalTorsionPathScoreError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen(js_name=explainPathScore,unchecked_return_type="[string, number, number][]")]
pub fn explain_path_score(
    #[wasm_bindgen(unchecked_param_type = "bigint")] score: JsValue,
    #[wasm_bindgen(unchecked_optional_param_type = "number")] size: JsValue,
) -> Result<Array, JsValue> {
    // COSMolKit❗✔️: PyTuple::new(py, ck::explain_path_score(score, size))
    Ok(ck::explain_path_score(
        u64_value(&score, "score")?,
        if size.is_undefined() {
            4
        } else {
            usize_value(&size, "size")?
        },
    )
    .into_iter()
    .map(|(s, b, p)| -> JsValue { Array::of3(&s.into(), &b.into(), &p.into()).into() })
    .collect())
}
#[wasm_bindgen]
pub struct AtomPairAtomCodeResult {
    inner: cosmolkit_wasm::AtomPairAtomCodeResult,
}
#[wasm_bindgen]
impl AtomPairAtomCodeResult {
    #[wasm_bindgen(getter)]
    pub fn code(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.code
        self.inner.code()
    }
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> Molecule {
        // COSMolKit❗✔️: inner: self.inner.molecule.clone(),
        Molecule {
            inner: Arc::new(self.inner.molecule()),
        }
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=topologicalTorsionPathScore)]
    pub fn topological_torsion_path_score(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Uint32Array")] path: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")]atom_codes:JsValue,
    ) -> Result<u64, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_path_score(&path, size, atom_codes.as_deref())
        let path = sequence(&path, "path")?
            .iter()
            .map(|v| usize_value(&v, "path"))
            .collect::<Result<Vec<_>, _>>()?;
        let codes = crate::atom_pair_parameters::optional_u32_array(&atom_codes, "atomCodes")?;
        self.inner
            .topological_torsion_path_score(&path, usize_value(&size, "size")?, codes.as_deref())
            .map_err(|e| score_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withAtomPairAtomCode)]
    pub fn with_atom_pair_atom_code(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] branch_subtract: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        use_legacy_stereo_perception: JsValue,
    ) -> Result<AtomPairAtomCodeResult, JsValue> {
        // COSMolKit❗✔️: .with_atom_pair_atom_code(
        self.inner
            .with_atom_pair_atom_code(
                usize_value(&atom_id, "atomId")?,
                if branch_subtract.is_undefined() {
                    0
                } else {
                    u32_value(&branch_subtract, "branchSubtract")?
                },
                if include_chirality.is_undefined() {
                    false
                } else {
                    bool_value(&include_chirality, "includeChirality")?
                },
                if use_legacy_stereo_perception.is_undefined() {
                    true
                } else {
                    bool_value(&use_legacy_stereo_perception, "useLegacyStereoPerception")?
                },
            )
            .map(|inner| AtomPairAtomCodeResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
}
