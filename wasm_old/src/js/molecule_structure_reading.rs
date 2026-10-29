//! Structural text reader projections using the canonical BIO conversion.
use crate::Molecule;
use crate::host_values::type_error;
use crate::io_errors::io_error;
use cosmolkit_wasm::rust as ck;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl Molecule {
    /// Parse PDB text through the shared BIO-to-molecule conversion.
    #[wasm_bindgen(js_name=fromPdb)]
    pub fn from_pdb(
        text: &str,
        #[wasm_bindgen(
            unchecked_optional_param_type = "BioMoleculeParams | BioMoleculeParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Self, JsValue> {
        if params.is_undefined() {
            return ck::Molecule::from_pdb(text)
                .map(|inner| Self {
                    inner: Arc::new(inner.into()),
                })
                .map_err(|e| io_error(&e).unwrap_or_else(|e| e));
        }
        let configuration =
            crate::bio_conversion::BioMoleculeParams::from_js_configuration(&params)?;
        Self::from_pdb_with_params(text, configuration.into())
    }
    #[wasm_bindgen(js_name=fromPdbWithParams)]
    pub fn from_pdb_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "BioMoleculeParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        let mut result = None;
        crate::bio_conversion::visit_params(
            &params,
            &mut |p: &crate::bio_conversion::BioMoleculeParams| {
                result = Some(
                    ck::Molecule::from_pdb_with_params(text, &p.inner)
                        .map(|inner| Self {
                            inner: Arc::new(inner.into()),
                        })
                        .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
                );
            },
        )?;
        result.ok_or_else(|| type_error("params"))?
    }
    /// Parse mmCIF text through the shared BIO-to-molecule conversion.
    #[wasm_bindgen(js_name=fromMmcif)]
    pub fn from_mmcif(
        text: &str,
        #[wasm_bindgen(
            unchecked_optional_param_type = "BioMoleculeParams | BioMoleculeParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Self, JsValue> {
        if params.is_undefined() {
            return ck::Molecule::from_mmcif(text)
                .map(|inner| Self {
                    inner: Arc::new(inner.into()),
                })
                .map_err(|e| io_error(&e).unwrap_or_else(|e| e));
        }
        let configuration =
            crate::bio_conversion::BioMoleculeParams::from_js_configuration(&params)?;
        Self::from_mmcif_with_params(text, configuration.into())
    }
    #[wasm_bindgen(js_name=fromMmcifWithParams)]
    pub fn from_mmcif_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "BioMoleculeParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        let mut result = None;
        crate::bio_conversion::visit_params(
            &params,
            &mut |p: &crate::bio_conversion::BioMoleculeParams| {
                result = Some(
                    ck::Molecule::from_mmcif_with_params(text, &p.inner)
                        .map(|inner| Self {
                            inner: Arc::new(inner.into()),
                        })
                        .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
                );
            },
        )?;
        result.ok_or_else(|| type_error("params"))?
    }
}
