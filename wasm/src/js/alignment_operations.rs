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
    pub fn alignment_transform_to(&self, reference: &Molecule) -> Result<AlignmentResult, JsValue> {
        let result = self
            .inner
            .alignment_transform_to(&reference.inner)
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(AlignmentResult { inner: result })
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
    pub fn best_alignment_to(&self, reference: &Molecule) -> Result<AlignmentResult, JsValue> {
        let result = self
            .inner
            .best_alignment_to(&reference.inner)
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(AlignmentResult { inner: result })
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
    pub fn best_rmsd_to(&self, reference: &Molecule) -> Result<f64, JsValue> {
        let result = self
            .inner
            .best_rmsd_to(&reference.inner)
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(result)
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
    pub fn coordinate_rmsd_to(&self, reference: &Molecule) -> Result<f64, JsValue> {
        let result = self
            .inner
            .coordinate_rmsd_to(&reference.inner)
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(result)
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
    pub fn all_conformer_best_rmsds(&self) -> Result<Vec<ConformerRmsd>, JsValue> {
        let result = self
            .inner
            .all_conformer_best_rmsds()
            .map_err(|error| alignment_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(result
            .into_iter()
            .map(|inner| ConformerRmsd { inner })
            .collect())
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
    pub fn with_alignment_to(&self, reference: &Molecule) -> Result<JsValue, JsValue> {
        let result = self
            .inner
            .with_alignment_to(&reference.inner)
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
    pub fn align_to_(&self, reference: &Molecule) -> Result<AlignmentResult, JsValue> {
        let result = self
            .inner
            .align_to_(&reference.inner)
            .map_err(|error| operation_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(AlignmentResult { inner: result })
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
    pub fn with_aligned_conformers(&self) -> Result<JsValue, JsValue> {
        let result = self
            .inner
            .with_aligned_conformers()
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
    pub fn align_conformers_(&self) -> Result<ConformerAlignmentReport, JsValue> {
        let result = self
            .inner
            .align_conformers_()
            .map_err(|error| operation_error(&error).unwrap_or_else(|conversion| conversion))?;
        Ok(ConformerAlignmentReport { inner: result })
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
