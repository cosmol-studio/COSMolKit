//! Canonical alignment operation projections. Native methods own host conversion.
use crate::Molecule;
use crate::alignment_parameters::*;
use crate::alignment_values::*;
use js_sys::Array;
use std::sync::Arc;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = alignmentTransformTo)]
    pub fn alignment_transform_to(
        &self,
        reference: &Molecule,
        #[wasm_bindgen(
            unchecked_optional_param_type = "AlignmentParameters | AlignmentParametersOptions"
        )]
        params: JsValue,
    ) -> Result<AlignmentResult, JsValue> {
        if params.is_undefined() {
            return {
                let result = self
                    .inner
                    .alignment_transform_to(&reference.inner)
                    .map_err(|error| {
                        alignment_error(&error).unwrap_or_else(|conversion| conversion)
                    })?;
                Ok(AlignmentResult { inner: result })
            };
        }
        let configuration_0 =
            crate::alignment_parameters::AlignmentParameters::from_js_configuration(&params)?;
        self.alignment_transform_to_with_params(reference, &configuration_0)
    }
    #[wasm_bindgen(js_name = alignmentTransformToWithParams)]
    pub fn alignment_transform_to_with_params(
        &self,
        reference: &Molecule,
        params: &AlignmentParameters,
    ) -> Result<AlignmentResult, JsValue> {
        let result = self
            .inner
            .alignment_transform_to_with_params(&reference.inner, &params.inner)
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(AlignmentResult { inner: result })
    }
    #[wasm_bindgen(js_name = bestAlignmentTo)]
    pub fn best_alignment_to(
        &self,
        reference: &Molecule,
        #[wasm_bindgen(
            unchecked_optional_param_type = "BestAlignmentParameters | BestAlignmentParametersOptions"
        )]
        params: JsValue,
    ) -> Result<AlignmentResult, JsValue> {
        if params.is_undefined() {
            return {
                let result = self
                    .inner
                    .best_alignment_to(&reference.inner)
                    .map_err(|error| {
                        alignment_error(&error).unwrap_or_else(|conversion| conversion)
                    })?;
                Ok(AlignmentResult { inner: result })
            };
        }
        let configuration_0 =
            crate::alignment_parameters::BestAlignmentParameters::from_js_configuration(&params)?;
        self.best_alignment_to_with_params(reference, &configuration_0)
    }
    #[wasm_bindgen(js_name = bestAlignmentToWithParams)]
    pub fn best_alignment_to_with_params(
        &self,
        reference: &Molecule,
        params: &BestAlignmentParameters,
    ) -> Result<AlignmentResult, JsValue> {
        let result = self
            .inner
            .best_alignment_to_with_params(&reference.inner, &params.inner)
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(AlignmentResult { inner: result })
    }
    #[wasm_bindgen(js_name = bestRmsdTo)]
    pub fn best_rmsd_to(
        &self,
        reference: &Molecule,
        #[wasm_bindgen(
            unchecked_optional_param_type = "BestAlignmentParameters | BestAlignmentParametersOptions"
        )]
        params: JsValue,
    ) -> Result<f64, JsValue> {
        if params.is_undefined() {
            return {
                let result = self.inner.best_rmsd_to(&reference.inner).map_err(|error| {
                    alignment_error(&error).unwrap_or_else(|conversion| conversion)
                })?;
                Ok(result)
            };
        }
        let configuration_0 =
            crate::alignment_parameters::BestAlignmentParameters::from_js_configuration(&params)?;
        self.best_rmsd_to_with_params(reference, &configuration_0)
    }
    #[wasm_bindgen(js_name = bestRmsdToWithParams)]
    pub fn best_rmsd_to_with_params(
        &self,
        reference: &Molecule,
        params: &BestAlignmentParameters,
    ) -> Result<f64, JsValue> {
        let result = self
            .inner
            .best_rmsd_to_with_params(&reference.inner, &params.inner)
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(result)
    }
    #[wasm_bindgen(js_name = coordinateRmsdTo)]
    pub fn coordinate_rmsd_to(
        &self,
        reference: &Molecule,
        #[wasm_bindgen(
            unchecked_optional_param_type = "CoordinateRmsdParameters | CoordinateRmsdParametersOptions"
        )]
        params: JsValue,
    ) -> Result<f64, JsValue> {
        if params.is_undefined() {
            return {
                let result = self
                    .inner
                    .coordinate_rmsd_to(&reference.inner)
                    .map_err(|error| {
                        alignment_error(&error).unwrap_or_else(|conversion| conversion)
                    })?;
                Ok(result)
            };
        }
        let configuration_0 =
            crate::alignment_parameters::CoordinateRmsdParameters::from_js_configuration(&params)?;
        self.coordinate_rmsd_to_with_params(reference, &configuration_0)
    }
    #[wasm_bindgen(js_name = coordinateRmsdToWithParams)]
    pub fn coordinate_rmsd_to_with_params(
        &self,
        reference: &Molecule,
        params: &CoordinateRmsdParameters,
    ) -> Result<f64, JsValue> {
        let result = self
            .inner
            .coordinate_rmsd_to_with_params(&reference.inner, &params.inner)
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(result)
    }
    #[wasm_bindgen(js_name = allConformerBestRmsds)]
    pub fn all_conformer_best_rmsds(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "AllConformerRmsdParameters | AllConformerRmsdParametersOptions"
        )]
        params: JsValue,
    ) -> Result<Vec<ConformerRmsd>, JsValue> {
        if params.is_undefined() {
            return {
                let result = self.inner.all_conformer_best_rmsds().map_err(|error| {
                    alignment_error(&error).unwrap_or_else(|conversion| conversion)
                })?;
                Ok(result
                    .into_iter()
                    .map(|inner| ConformerRmsd { inner })
                    .collect())
            };
        }
        let configuration_0 =
            crate::alignment_parameters::AllConformerRmsdParameters::from_js_configuration(
                &params,
            )?;
        self.all_conformer_best_rmsds_with_params(&configuration_0)
    }
    #[wasm_bindgen(js_name = allConformerBestRmsdsWithParams)]
    pub fn all_conformer_best_rmsds_with_params(
        &self,
        params: &AllConformerRmsdParameters,
    ) -> Result<Vec<ConformerRmsd>, JsValue> {
        let result = self
            .inner
            .all_conformer_best_rmsds_with_params(&params.inner)
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(result
            .into_iter()
            .map(|inner| ConformerRmsd { inner })
            .collect())
    }
    #[wasm_bindgen(js_name = withAlignmentTo, unchecked_return_type = "[Molecule, AlignmentResult]")]
    pub fn with_alignment_to(
        &self,
        reference: &Molecule,
        #[wasm_bindgen(
            unchecked_optional_param_type = "AlignmentParameters | AlignmentParametersOptions"
        )]
        params: JsValue,
    ) -> Result<JsValue, JsValue> {
        if params.is_undefined() {
            return {
                let result = self
                    .inner
                    .with_alignment_to(&reference.inner)
                    .map_err(|error| {
                        operation_error(&error).unwrap_or_else(|conversion| conversion)
                    })?;
                let values = Array::new();
                values.push(
                    &Molecule {
                        inner: Arc::new(result.0),
                    }
                    .into(),
                );
                values.push(&AlignmentResult { inner: result.1 }.into());
                Ok(values.into())
            };
        }
        let configuration_0 =
            crate::alignment_parameters::AlignmentParameters::from_js_configuration(&params)?;
        self.with_alignment_to_with_params(reference, &configuration_0)
    }
    #[wasm_bindgen(js_name = withAlignmentToWithParams, unchecked_return_type = "[Molecule, AlignmentResult]")]
    pub fn with_alignment_to_with_params(
        &self,
        reference: &Molecule,
        params: &AlignmentParameters,
    ) -> Result<JsValue, JsValue> {
        let result = self
            .inner
            .with_alignment_to_with_params(&reference.inner, &params.inner)
            .map_err(|error| operation_error(&error).unwrap_or_else(|conversion| conversion))?;
        let values = Array::new();
        values.push(
            &Molecule {
                inner: Arc::new(result.0),
            }
            .into(),
        );
        values.push(&AlignmentResult { inner: result.1 }.into());
        Ok(values.into())
    }
    #[wasm_bindgen(js_name = alignTo)]
    pub fn align_to_(
        &self,
        reference: &Molecule,
        #[wasm_bindgen(
            unchecked_optional_param_type = "AlignmentParameters | AlignmentParametersOptions"
        )]
        params: JsValue,
    ) -> Result<AlignmentResult, JsValue> {
        if params.is_undefined() {
            return {
                let result = self.inner.align_to_(&reference.inner).map_err(|error| {
                    operation_error(&error).unwrap_or_else(|conversion| conversion)
                })?;
                Ok(AlignmentResult { inner: result })
            };
        }
        let configuration_0 =
            crate::alignment_parameters::AlignmentParameters::from_js_configuration(&params)?;
        self.align_to_with_params_(reference, &configuration_0)
    }
    #[wasm_bindgen(js_name = alignToWithParams)]
    pub fn align_to_with_params_(
        &self,
        reference: &Molecule,
        params: &AlignmentParameters,
    ) -> Result<AlignmentResult, JsValue> {
        let result = self
            .inner
            .align_to_with_params_(&reference.inner, &params.inner)
            .map_err(|error| operation_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(AlignmentResult { inner: result })
    }
    #[wasm_bindgen(js_name = withAlignedConformers, unchecked_return_type = "[Molecule, ConformerAlignmentReport]")]
    pub fn with_aligned_conformers(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "ConformerAlignmentParameters | ConformerAlignmentParametersOptions"
        )]
        params: JsValue,
    ) -> Result<JsValue, JsValue> {
        if params.is_undefined() {
            return {
                let result = self.inner.with_aligned_conformers().map_err(|error| {
                    operation_error(&error).unwrap_or_else(|conversion| conversion)
                })?;
                let values = Array::new();
                values.push(
                    &Molecule {
                        inner: Arc::new(result.0),
                    }
                    .into(),
                );
                values.push(&ConformerAlignmentReport { inner: result.1 }.into());
                Ok(values.into())
            };
        }
        let configuration_0 =
            crate::alignment_parameters::ConformerAlignmentParameters::from_js_configuration(
                &params,
            )?;
        self.with_aligned_conformers_with_params(&configuration_0)
    }
    #[wasm_bindgen(js_name = withAlignedConformersWithParams, unchecked_return_type = "[Molecule, ConformerAlignmentReport]")]
    pub fn with_aligned_conformers_with_params(
        &self,
        params: &ConformerAlignmentParameters,
    ) -> Result<JsValue, JsValue> {
        let result = self
            .inner
            .with_aligned_conformers_with_params(&params.inner)
            .map_err(|error| operation_error(&error).unwrap_or_else(|conversion| conversion))?;
        let values = Array::new();
        values.push(
            &Molecule {
                inner: Arc::new(result.0),
            }
            .into(),
        );
        values.push(&ConformerAlignmentReport { inner: result.1 }.into());
        Ok(values.into())
    }
    #[wasm_bindgen(js_name = alignConformers)]
    pub fn align_conformers_(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "ConformerAlignmentParameters | ConformerAlignmentParametersOptions"
        )]
        params: JsValue,
    ) -> Result<ConformerAlignmentReport, JsValue> {
        if params.is_undefined() {
            return {
                let result = self.inner.align_conformers_().map_err(|error| {
                    operation_error(&error).unwrap_or_else(|conversion| conversion)
                })?;
                Ok(ConformerAlignmentReport { inner: result })
            };
        }
        let configuration_0 =
            crate::alignment_parameters::ConformerAlignmentParameters::from_js_configuration(
                &params,
            )?;
        self.align_conformers_with_params_(&configuration_0)
    }
    #[wasm_bindgen(js_name = alignConformersWithParams)]
    pub fn align_conformers_with_params_(
        &self,
        params: &ConformerAlignmentParameters,
    ) -> Result<ConformerAlignmentReport, JsValue> {
        let result = self
            .inner
            .align_conformers_with_params_(&params.inner)
            .map_err(|error| operation_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(ConformerAlignmentReport { inner: result })
    }
}
