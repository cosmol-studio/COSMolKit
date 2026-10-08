//! All sixteen registry embedding operations and complete reports.
use crate::{
    Molecule,
    alignment_values::{operation_error, set, source_error},
    embed_parameters::EmbedParams,
    host_values::u32_value,
};
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct EmbedMoleculeResult {
    inner: cosmolkit_wasm::EmbedMoleculeResult,
}
#[wasm_bindgen]
impl EmbedMoleculeResult {
    pub fn molecule(&self) -> Molecule {
        Molecule {
            inner: Arc::new(self.inner.molecule()),
        }
    }
    pub fn params(&self) -> EmbedParams {
        EmbedParams {
            inner: self.inner.params().clone(),
        }
    }
    #[wasm_bindgen(js_name=confId)]
    pub fn conf_id(&self) -> i32 {
        self.inner.conf_id()
    }
    pub fn ok(&self) -> bool {
        self.inner.ok()
    }
}
#[wasm_bindgen]
pub struct EmbedMultipleConfsResult {
    inner: cosmolkit_wasm::EmbedMultipleConfsResult,
}
#[wasm_bindgen]
impl EmbedMultipleConfsResult {
    pub fn molecule(&self) -> Molecule {
        Molecule {
            inner: Arc::new(self.inner.molecule()),
        }
    }
    pub fn params(&self) -> EmbedParams {
        EmbedParams {
            inner: self.inner.params().clone(),
        }
    }
    #[wasm_bindgen(js_name=confIds)]
    pub fn conf_ids(&self) -> Vec<i32> {
        self.inner.conf_ids().to_vec()
    }
    #[wasm_bindgen(js_name=requestedNumConfs)]
    pub fn requested_num_confs(&self) -> u32 {
        self.inner.requested_num_confs()
    }
    #[wasm_bindgen(js_name=generatedCount)]
    pub fn generated_count(&self) -> usize {
        self.inner.generated_count()
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=with3dConformer)]
    pub fn with_3d_conformer(&self) -> Result<Molecule, JsValue> {
        // COSMolKit❗✔️: .with_3d_conformer()
        self.inner
            .with_3d_conformer()
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=with3dConformerWithParams)]
    pub fn with_3d_conformer_with_params(&self, params: &EmbedParams) -> Result<Molecule, JsValue> {
        // COSMolKit❗✔️: .with_3d_conformer_with_params(&params.inner)
        self.inner
            .with_3d_conformer_with_params(&params.inner)
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=embed3dConformer)]
    pub fn embed_3d_conformer_(&self) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .embed_3d_conformer_()
        self.inner
            .embed_3d_conformer_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=embed3dConformerWithParams)]
    pub fn embed_3d_conformer_with_params_(&self, params: &EmbedParams) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .embed_3d_conformer_with_params_(&params.inner)
        self.inner
            .embed_3d_conformer_with_params_(&params.inner)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=with3dConformerResult)]
    pub fn with_3d_conformer_result(&self) -> Result<EmbedMoleculeResult, JsValue> {
        // COSMolKit❗✔️: .with_3d_conformer_result()
        self.inner
            .with_3d_conformer_result()
            .map(|inner| EmbedMoleculeResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=with3dConformerResultWithParams)]
    pub fn with_3d_conformer_result_with_params(
        &self,
        params: &EmbedParams,
    ) -> Result<EmbedMoleculeResult, JsValue> {
        // COSMolKit❗✔️: .with_3d_conformer_result_with_params(&params.inner)
        self.inner
            .with_3d_conformer_result_with_params(&params.inner)
            .map(|inner| EmbedMoleculeResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=embed3dConformerResult)]
    pub fn embed_3d_conformer_result_(&self) -> Result<EmbedMoleculeResult, JsValue> {
        // COSMolKit❗✔️: .embed_3d_conformer_result_()
        self.inner
            .embed_3d_conformer_result_()
            .map(|inner| EmbedMoleculeResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=embed3dConformerResultWithParams)]
    pub fn embed_3d_conformer_result_with_params_(
        &self,
        params: &EmbedParams,
    ) -> Result<EmbedMoleculeResult, JsValue> {
        // COSMolKit❗✔️: .embed_3d_conformer_result_with_params_(&params.inner)
        self.inner
            .embed_3d_conformer_result_with_params_(&params.inner)
            .map(|inner| EmbedMoleculeResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=with3dConformers)]
    pub fn with_3d_conformers(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_confs: JsValue,
    ) -> Result<Molecule, JsValue> {
        // COSMolKit❗✔️: .with_3d_conformers(num_confs)
        self.inner
            .with_3d_conformers(u32_value(&num_confs, "numConfs")?)
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=with3dConformersWithParams)]
    pub fn with_3d_conformers_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_confs: JsValue,
        params: &EmbedParams,
    ) -> Result<Molecule, JsValue> {
        // COSMolKit❗✔️: .with_3d_conformers_with_params(num_confs, &params.inner)
        self.inner
            .with_3d_conformers_with_params(u32_value(&num_confs, "numConfs")?, &params.inner)
            .map(|inner| Molecule {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=embed3dConformers)]
    pub fn embed_3d_conformers_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_confs: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .embed_3d_conformers_(num_confs)
        self.inner
            .embed_3d_conformers_(u32_value(&num_confs, "numConfs")?)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=embed3dConformersWithParams)]
    pub fn embed_3d_conformers_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_confs: JsValue,
        params: &EmbedParams,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .embed_3d_conformers_with_params_(num_confs, &params.inner)
        self.inner
            .embed_3d_conformers_with_params_(u32_value(&num_confs, "numConfs")?, &params.inner)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=with3dConformersResult)]
    pub fn with_3d_conformers_result(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_confs: JsValue,
    ) -> Result<EmbedMultipleConfsResult, JsValue> {
        // COSMolKit❗✔️: .with_3d_conformers_result(num_confs)
        self.inner
            .with_3d_conformers_result(u32_value(&num_confs, "numConfs")?)
            .map(|inner| EmbedMultipleConfsResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=with3dConformersResultWithParams)]
    pub fn with_3d_conformers_result_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_confs: JsValue,
        params: &EmbedParams,
    ) -> Result<EmbedMultipleConfsResult, JsValue> {
        // COSMolKit❗✔️: .with_3d_conformers_result_with_params(num_confs, &params.inner)
        self.inner
            .with_3d_conformers_result_with_params(
                u32_value(&num_confs, "numConfs")?,
                &params.inner,
            )
            .map(|inner| EmbedMultipleConfsResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=embed3dConformersResult)]
    pub fn embed_3d_conformers_result_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_confs: JsValue,
    ) -> Result<EmbedMultipleConfsResult, JsValue> {
        // COSMolKit❗✔️: .embed_3d_conformers_result_(num_confs)
        self.inner
            .embed_3d_conformers_result_(u32_value(&num_confs, "numConfs")?)
            .map(|inner| EmbedMultipleConfsResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=embed3dConformersResultWithParams)]
    pub fn embed_3d_conformers_result_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_confs: JsValue,
        params: &EmbedParams,
    ) -> Result<EmbedMultipleConfsResult, JsValue> {
        // COSMolKit❗✔️: .embed_3d_conformers_result_with_params_(num_confs, &params.inner)
        self.inner
            .embed_3d_conformers_result_with_params_(
                u32_value(&num_confs, "numConfs")?,
                &params.inner,
            )
            .map(|inner| EmbedMultipleConfsResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=num3dConformers)]
    pub fn num_3d_conformers(&self) -> usize {
        self.inner.num_3d_conformers()
    }
    #[wasm_bindgen(js_name=dgBoundsMatrix,unchecked_return_type="number[][]")]
    pub fn dg_bounds_matrix(&self) -> Result<JsValue, JsValue> {
        self.inner
            .dg_bounds_matrix()
            .map(|rows| {
                rows.into_iter()
                    .map(|row| -> JsValue {
                        row.into_iter()
                            .map(JsValue::from)
                            .collect::<js_sys::Array>()
                            .into()
                    })
                    .collect::<js_sys::Array>()
                    .into()
            })
            .map_err(|e| crate::conformer_errors::run_error(&e).unwrap_or_else(|e| e))
    }
}
