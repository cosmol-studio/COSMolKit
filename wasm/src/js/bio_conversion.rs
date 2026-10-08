//! BIO chemistry conversion is delegated to canonical root APIs, preserving typed failures.
use crate::alignment_values::{set, source_error};
use crate::host_values::{bool_value, u32_value};
use crate::{Molecule, bio_protein::Protein, bio_readers::BioStructure};
use cosmolkit_wasm::rust as ck;
use std::{error::Error as RustError, sync::Arc};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct BioMoleculeParams {
    inner: ck::BioMoleculeParams,
}
#[wasm_bindgen]
impl BioMoleculeParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] sanitize: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_hs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] flavor: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] proximity_bonding: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::BioMoleculeParams::default();
        if !sanitize.is_undefined() {
            inner.sanitize = bool_value(&sanitize, "sanitize")?;
        }
        if !remove_hs.is_undefined() {
            inner.remove_hs = bool_value(&remove_hs, "removeHs")?;
        }
        if !flavor.is_undefined() {
            inner.flavor = u32_value(&flavor, "flavor")?;
        }
        if !proximity_bonding.is_undefined() {
            inner.proximity_bonding = bool_value(&proximity_bonding, "proximityBonding")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter)]
    pub fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[wasm_bindgen(getter,js_name=removeHs)]
    pub fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    #[wasm_bindgen(getter)]
    pub fn flavor(&self) -> u32 {
        self.inner.flavor
    }
    #[wasm_bindgen(getter,js_name=proximityBonding)]
    pub fn proximity_bonding(&self) -> bool {
        self.inner.proximity_bonding
    }
}
#[wasm_bindgen(
    inline_js = "export function visitBioMoleculeParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid BioMoleculeParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitBioMoleculeParams)]
    fn visit_params(
        value: &JsValue,
        visit: &mut dyn FnMut(&BioMoleculeParams),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct BioMoleculeError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl BioMoleculeError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
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
pub(crate) fn molecule_error(source: &ck::BioMoleculeError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: use ck::BioMoleculeError as E;
    use ck::BioMoleculeError as E;
    let kind = match source {
        E::Conversion(_) => "Conversion",
        E::Construction(_) => "Construction",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("BioMoleculeError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        BioMoleculeError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
pub struct BioMoleculeConversionError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl BioMoleculeConversionError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
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
pub(crate) fn conversion_error(
    source: &ck::BioMoleculeConversionError,
) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: use ck::BioMoleculeConversionError as E;
    use ck::BioMoleculeConversionError as E;
    let kind = match source {
        E::Structure(_) => "Structure",
        E::Topology(_) => "Topology",
        E::Chemistry(_) => "Chemistry",
        E::Sanitize(_) => "Sanitize",
        E::Hydrogens(_) => "Hydrogens",
        E::Valence(_) => "Valence",
        E::Stereo(_) => "Stereo",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("BioMoleculeConversionError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        BioMoleculeConversionError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=toMolecule)]
    pub fn to_molecule(&self) -> Result<Molecule, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::bio_structure_to_molecule(&self.inner)
        cosmolkit_wasm::bio_structure_to_molecule(&self.inner)
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| molecule_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toMoleculeWithParams)]
    pub fn to_molecule_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "BioMoleculeParams")] params: JsValue,
    ) -> Result<Molecule, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::bio_structure_to_molecule_with_params(&self.inner,&params.inner)
        let mut result = None;
        visit_params(&params, &mut |params: &BioMoleculeParams| {
            result = Some(
                cosmolkit_wasm::bio_structure_to_molecule_with_params(&self.inner, &params.inner)
                    .map(|inner| Molecule {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| molecule_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| crate::host_values::type_error("params"))?
    }
}
#[wasm_bindgen]
impl Protein {
    #[wasm_bindgen(js_name=toMolecule)]
    pub fn to_molecule(&self) -> Result<Molecule, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::protein_to_molecule(&self.inner)
        cosmolkit_wasm::protein_to_molecule(&self.inner)
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| molecule_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toMoleculeWithParams)]
    pub fn to_molecule_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "BioMoleculeParams")] params: JsValue,
    ) -> Result<Molecule, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::protein_to_molecule_with_params(&self.inner,&params.inner)
        let mut result = None;
        visit_params(&params, &mut |params: &BioMoleculeParams| {
            result = Some(
                cosmolkit_wasm::protein_to_molecule_with_params(&self.inner, &params.inner)
                    .map(|inner| Molecule {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| molecule_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| crate::host_values::type_error("params"))?
    }
}
