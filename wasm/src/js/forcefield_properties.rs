//! UFF coverage and full MMFF atom properties via the canonical facade.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use crate::host_values::{type_error, usize_value};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum MmffVariant {
    Mmff94,
    Mmff94s,
}
#[wasm_bindgen]
pub struct MmffPropertiesParams {
    inner: ck::MmffPropertiesParams,
}
#[wasm_bindgen]
impl MmffPropertiesParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "string")] mmff_variant: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: mmff_variant: mmff_variant.into(),
        let mut inner = ck::MmffPropertiesParams::default();
        if !mmff_variant.is_undefined() {
            inner.mmff_variant = mmff_variant
                .as_string()
                .ok_or_else(|| type_error("mmffVariant"))?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=mmffVariant)]
    pub fn mmff_variant(&self) -> String {
        self.inner.mmff_variant.clone()
    }
}
#[wasm_bindgen]
pub struct MmffAtomProperties {
    inner: ck::MmffAtomProperties,
}
#[wasm_bindgen]
impl MmffAtomProperties {
    #[wasm_bindgen(js_name=atomType)]
    pub fn atom_type(&self) -> u8 {
        // COSMolKit❗✔️: self.inner.atom_type
        self.inner.atom_type()
    }
    #[wasm_bindgen(js_name=formalCharge)]
    pub fn formal_charge(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.formal_charge
        self.inner.formal_charge()
    }
    #[wasm_bindgen(js_name=partialCharge)]
    pub fn partial_charge(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.partial_charge
        self.inner.partial_charge()
    }
}
#[wasm_bindgen]
pub struct MmffProperties {
    inner: ck::MmffProperties,
}
#[wasm_bindgen]
impl MmffProperties {
    #[wasm_bindgen(js_name=isValid)]
    pub fn is_valid(&self) -> bool {
        // COSMolKit❗✔️: self.inner.is_valid()
        self.inner.is_valid()
    }
    pub fn variant(&self) -> MmffVariant {
        // COSMolKit❗✔️: match self.inner.variant() {
        match self.inner.variant() {
            ck::MmffVariant::Mmff94 => MmffVariant::Mmff94,
            ck::MmffVariant::Mmff94s => MmffVariant::Mmff94s,
        }
    }
    #[wasm_bindgen(unchecked_return_type = "MmffAtomProperties[]")]
    pub fn atoms(&self) -> Array {
        // COSMolKit❗✔️: .map(|inner| MmffAtomProperties { inner })
        self.inner
            .atoms()
            .iter()
            .copied()
            .map(|inner| JsValue::from(MmffAtomProperties { inner }))
            .collect()
    }
    #[wasm_bindgen(js_name=atomType)]
    pub fn atom_type(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom_index: JsValue,
    ) -> Result<u8, JsValue> {
        // COSMolKit❗✔️: .atom_type(atom_index)
        self.inner
            .atom_type(usize_value(&atom_index, "atomIndex")?)
            .map_err(mmff_error)
    }
    #[wasm_bindgen(js_name=formalCharge)]
    pub fn formal_charge(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom_index: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .formal_charge(atom_index)
        self.inner
            .formal_charge(usize_value(&atom_index, "atomIndex")?)
            .map_err(mmff_error)
    }
    #[wasm_bindgen(js_name=partialCharge)]
    pub fn partial_charge(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] atom_index: JsValue,
    ) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: .partial_charge(atom_index)
        self.inner
            .partial_charge(usize_value(&atom_index, "atomIndex")?)
            .map_err(mmff_error)
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum UffParameterErrorKind {
    Preparation,
    ParameterTable,
    Typing,
}
#[wasm_bindgen]
pub struct UffParameterError {
    kind: UffParameterErrorKind,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl UffParameterError {
    pub fn kind(&self) -> UffParameterErrorKind {
        self.kind
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
pub(crate) fn uff_parameter_error(source: &ck::UffParameterError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: crate::canonical_error_values::UffParameterErrorKind::from(cause.kind()),
    let (kind, tag) = match source.kind() {
        ck::UffParameterErrorKind::Preparation => {
            (UffParameterErrorKind::Preparation, "Preparation")
        }
        ck::UffParameterErrorKind::ParameterTable => {
            (UffParameterErrorKind::ParameterTable, "ParameterTable")
        }
        ck::UffParameterErrorKind::Typing => (UffParameterErrorKind::Typing, "Typing"),
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("UffParameterError");
    let e: JsValue = e.into();
    set(&e, "kind", tag.into())?;
    set(&e, "_kind", JsValue::from(kind as u32))?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        UffParameterError {
            kind,
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(e)
}
fn mmff_error(e: ck::MmffMolPropertiesError) -> JsValue {
    mmff_properties_error(&e).unwrap_or_else(|e| e)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=uffHasAllMoleculeParams)]
    pub fn uff_has_all_molecule_params(&self) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: fn uff_has_all_molecule_params(&self)
        self.inner
            .uff_has_all_molecule_params()
            .map_err(|e| uff_query_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=mmffHasAllMoleculeParams)]
    pub fn mmff_has_all_molecule_params(&self) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: fn mmff_has_all_molecule_params(&self)
        self.inner
            .mmff_has_all_molecule_params()
            .map_err(mmff_error)
    }
    #[wasm_bindgen(js_name=mmffProperties)]
    pub fn mmff_properties(&self) -> Result<MmffProperties, JsValue> {
        // COSMolKit❗✔️: fn mmff_properties(&self)
        self.inner
            .mmff_properties()
            .map(|inner| MmffProperties { inner })
            .map_err(mmff_error)
    }
    #[wasm_bindgen(js_name=mmffPropertiesWithParams)]
    pub fn mmff_properties_with_params(
        &self,
        params: &MmffPropertiesParams,
    ) -> Result<MmffProperties, JsValue> {
        // COSMolKit❗✔️: .mmff_properties_with_params(&params.inner)
        self.inner
            .mmff_properties_with_params(&params.inner)
            .map(|inner| MmffProperties { inner })
            .map_err(mmff_error)
    }
}
#[wasm_bindgen]
pub struct UffParameterQueryError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl UffParameterQueryError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "uff_parameters".into()
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
}
pub(crate) fn uff_query_error(source: &ck::UffParameterQueryError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: .setattr("domain", "uff_parameters")
    use ck::UffParameterQueryError as E;
    let kind = match source {
        E::Cache(_) => "Cache",
        E::Parameters(_) => "Parameters",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("UffParameterQueryError");
    let e: JsValue = e.into();
    set(&e, "domain", "uff_parameters".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        UffParameterQueryError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
pub struct MmffMolPropertiesError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl MmffMolPropertiesError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "mmff_properties".into()
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
}
pub(crate) fn mmff_properties_error(
    source: &ck::MmffMolPropertiesError,
) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: .setattr("domain", "mmff_properties")
    use ck::MmffMolPropertiesError as E;
    let kind = match source {
        E::Params(_) => "Params",
        E::Kekulize(_) => "Kekulize",
        E::Aromaticity(_) => "Aromaticity",
        E::RingFinding(_) => "RingFinding",
        E::Valence(_) => "Valence",
        E::AtomIndexOutOfRange { .. } => "AtomIndexOutOfRange",
        E::AtomTypePropertiesMissing { .. } => "AtomTypePropertiesMissing",
        E::AtomTypePbciMissing { .. } => "AtomTypePbciMissing",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("MmffMolPropertiesError");
    let e: JsValue = e.into();
    set(&e, "domain", "mmff_properties".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        MmffMolPropertiesError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(e)
}
