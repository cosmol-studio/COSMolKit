//! Molecular structural IO projections; parsing and chemistry belong to CK.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use crate::host_values::{bool_value, type_error, u32_value};
use crate::io_errors::io_error;
use crate::mol_sdf_writing::{MolCoordinateSelection, visit_selection};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;

/// PDB writer options. No coordinates are generated or written back.
#[wasm_bindgen]
pub struct PdbWriteParams {
    inner: ck::PdbWriteParams,
}
#[cosmolkit_wasm::javascript_options("PdbWriteParamsOptions")]
#[wasm_bindgen]
impl PdbWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "MolCoordinateSelection | null")]
        coordinate_selection: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] all_conformers: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] flavor: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::PdbWriteParams::default();
        if !coordinate_selection.is_undefined() && !coordinate_selection.is_null() {
            visit_selection(&coordinate_selection, &mut |s: &MolCoordinateSelection| {
                inner.coordinate_selection = s.inner;
            })?;
        }
        if !all_conformers.is_undefined() {
            inner.all_conformers = bool_value(&all_conformers, "allConformers")?;
        }
        if !flavor.is_undefined() {
            inner.flavor = u32_value(&flavor, "flavor")?;
        }
        Ok(Self { inner })
    }
    /// Unique-or-error automatic resolution or an explicit dimension and stored ID.
    #[wasm_bindgen(getter, js_name=coordinateSelection)]
    pub fn coordinate_selection(&self) -> MolCoordinateSelection {
        MolCoordinateSelection {
            inner: self.inner.coordinate_selection,
        }
    }
    /// Write all stored sets: 2D then 3D, each in storage order. Mutually
    /// exclusive with explicit single-set selection. PDB loses CK IDs/dimensions.
    #[wasm_bindgen(getter, js_name=allConformers)]
    pub fn all_conformers(&self) -> bool {
        self.inner.all_conformers
    }
    /// RDKit PDB writer flavor bit mask.
    #[wasm_bindgen(getter)]
    pub fn flavor(&self) -> u32 {
        self.inner.flavor
    }
}
#[wasm_bindgen(
    inline_js = "export function visitPdbWriteParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid PdbWriteParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch, js_name=visitPdbWriteParams)]
    fn visit_params(value: &JsValue, visit: &mut dyn FnMut(&PdbWriteParams))
    -> Result<(), JsValue>;
}

#[wasm_bindgen]
impl Molecule {
    /// Return PDB text, never mutate the molecule. Multiple coordinate sets
    /// require dimension/ID selection or explicit allConformers.
    #[wasm_bindgen(js_name=toPdb)]
    pub fn to_pdb(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "PdbWriteParams | PdbWriteParamsOptions")]
        params: JsValue,
    ) -> Result<String, JsValue> {
        if params.is_undefined() {
            return self
                .inner
                .to_pdb()
                .map_err(|e| io_error(&e).unwrap_or_else(|e| e));
        }
        let configuration = PdbWriteParams::from_js_configuration(&params)?;
        self.to_pdb_with_params(configuration.into())
    }
    #[wasm_bindgen(js_name=toPdbWithParams)]
    pub fn to_pdb_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "PdbWriteParams")] params: JsValue,
    ) -> Result<String, JsValue> {
        let mut result = None;
        visit_params(&params, &mut |p: &PdbWriteParams| {
            result = Some(
                self.inner
                    .to_pdb_with_params(&p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}

/// Structured PDB writer failure. The thrown Error exposes this value as detail.
#[wasm_bindgen]
pub struct PdbWriteError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl PdbWriteError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "io".into()
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
    #[wasm_bindgen(getter, js_name=twoD, unchecked_return_type="number | null")]
    pub fn two_d(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"twoD".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter, js_name=threeD, unchecked_return_type="number | null")]
    pub fn three_d(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"threeD".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "CoordinateDimension | null")]
    pub fn dimension(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"dimension".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter, js_name=coordinateId, unchecked_return_type="number | null")]
    pub fn coordinate_id(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"coordinateId".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn feature(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"feature".into()).unwrap_or(JsValue::NULL)
    }
}
pub(crate) fn pdb_write_error(source: &ck::PdbWriteError) -> Result<JsValue, JsValue> {
    use ck::PdbWriteError as E;
    let fields = js_sys::Object::new();
    for field in ["twoD", "threeD", "dimension", "coordinateId", "feature"] {
        set(&fields, field, JsValue::NULL)?;
    }
    let kind = match source {
        E::Property(_) => "Property",
        E::ConformerNotFound { id } => {
            set(&fields, "coordinateId", (*id as u32).into())?;
            "ConformerNotFound"
        }
        E::Unsupported { feature } => {
            set(&fields, "feature", (*feature).into())?;
            "Unsupported"
        }
        E::Topology(_) => "Topology",
        E::Coordinates(_) => "Coordinates",
        E::Kekulize(_) => "Kekulize",
        E::ConflictingCoordinateSelection => "ConflictingCoordinateSelection",
        E::AmbiguousCoordinates { two_d, three_d } => {
            set(&fields, "twoD", (*two_d as u32).into())?;
            set(&fields, "threeD", (*three_d as u32).into())?;
            "AmbiguousCoordinates"
        }
        E::MissingCoordinate { dimension, id } => {
            set(
                &fields,
                "dimension",
                JsValue::from(match dimension {
                    ck::CoordinateDimension::TwoD => 0u32,
                    ck::CoordinateDimension::ThreeD => 1u32,
                }),
            )?;
            set(&fields, "coordinateId", (*id as u32).into())?;
            "MissingCoordinate"
        }
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("PdbWriteError");
    let e: JsValue = error.into();
    set(&e, "domain", "io".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    for key in js_sys::Object::keys(&fields).iter() {
        let value = js_sys::Reflect::get(&fields, &key)?;
        if !value.is_null() {
            set(&e, &key.as_string().unwrap(), value)?;
        }
    }
    set(
        &e,
        "detail",
        PdbWriteError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
