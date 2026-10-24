//! Public facade alignment calls; no geometry is computed in the projection.
use crate::{
    Molecule, MoleculeInput, alignment_errors::AlignmentError, alignment_values::*, fields::Field,
    operation_errors::OperationError, storage::NativeField,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "readonly [Molecule, AlignmentResult]")]
    pub type AlignedMolecule;
    #[wasm_bindgen(typescript_type = "readonly [Molecule, ConformerAlignmentReport]")]
    pub type AlignedConformers;
}
macro_rules! queries {($(($fn:ident,$js:ident,$params:ident,$input:ident,$native:ident,$out:ty,$map:expr)),+)=>{#[wasm_bindgen]impl Molecule{$(
 #[wasm_bindgen(js_name=$js)]pub fn $fn(&self,reference:MoleculeInput,params:Option<$input>)->Result<$out,JsValue>{Molecule::with_ref(&reference,|reference|$params::with_input(params.as_ref(),|p|self.0.$native(reference,p).map($map).map_err(AlignmentError::from_rust)))}
)+}};}
queries!(
    (
        alignment_transform_to,
        alignmentTransformTo,
        AlignmentParameters,
        AlignmentParametersInput,
        alignment_transform_to_with_params,
        AlignmentResult,
        AlignmentResult
    ),
    (
        best_alignment_to,
        bestAlignmentTo,
        BestAlignmentParameters,
        BestAlignmentParametersInput,
        best_alignment_to_with_params,
        AlignmentResult,
        AlignmentResult
    ),
    (
        best_rmsd_to,
        bestRmsdTo,
        BestAlignmentParameters,
        BestAlignmentParametersInput,
        best_rmsd_to_with_params,
        f64,
        std::convert::identity
    ),
    (
        coordinate_rmsd_to,
        coordinateRmsdTo,
        CoordinateRmsdParameters,
        CoordinateRmsdParametersInput,
        coordinate_rmsd_to_with_params,
        f64,
        std::convert::identity
    )
);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=allConformerBestRmsds)]
    pub fn all_conformer_best_rmsds(
        &self,
        params: Option<AllConformerRmsdParametersInput>,
    ) -> Result<ConformerRmsds, JsValue> {
        AllConformerRmsdParameters::with_input(params.as_ref(), |p| {
            let v = self
                .0
                .all_conformer_best_rmsds_with_params(p)
                .map_err(AlignmentError::from_rust)?;
            Vec::<ck::ConformerRmsd>::get(NativeField::owned(v))
        })
    }
    #[wasm_bindgen(js_name=withAlignmentTo)]
    pub fn with_alignment_to(
        &self,
        reference: MoleculeInput,
        params: Option<AlignmentParametersInput>,
    ) -> Result<AlignedMolecule, JsValue> {
        Molecule::with_ref(&reference, |r| {
            AlignmentParameters::with_input(params.as_ref(), |p| {
                let (m, report) = self
                    .0
                    .with_alignment_to_with_params(r, p)
                    .map_err(OperationError::from_rust)?;
                let pair = js_sys::Array::of2(&Molecule(m).into(), &AlignmentResult(report).into());
                js_sys::Object::freeze(&pair);
                Ok(pair.unchecked_into())
            })
        })
    }
    #[wasm_bindgen(js_name=alignTo_)]
    pub fn align_to(
        &mut self,
        reference: MoleculeInput,
        params: Option<AlignmentParametersInput>,
    ) -> Result<AlignmentResult, JsValue> {
        Molecule::with_ref(&reference, |r| {
            AlignmentParameters::with_input(params.as_ref(), |p| {
                self.0
                    .align_to_with_params_(r, p)
                    .map(AlignmentResult)
                    .map_err(OperationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=withAlignedConformers)]
    pub fn with_aligned_conformers(
        &self,
        params: Option<ConformerAlignmentParametersInput>,
    ) -> Result<AlignedConformers, JsValue> {
        ConformerAlignmentParameters::with_input(params.as_ref(), |p| {
            let (m, report) = self
                .0
                .with_aligned_conformers_with_params(p)
                .map_err(OperationError::from_rust)?;
            let pair = js_sys::Array::of2(
                &Molecule(m).into(),
                &ConformerAlignmentReport(report).into(),
            );
            js_sys::Object::freeze(&pair);
            Ok(pair.unchecked_into())
        })
    }
    #[wasm_bindgen(js_name=alignConformers_)]
    pub fn align_conformers(
        &mut self,
        params: Option<ConformerAlignmentParametersInput>,
    ) -> Result<ConformerAlignmentReport, JsValue> {
        ConformerAlignmentParameters::with_input(params.as_ref(), |p| {
            self.0
                .align_conformers_with_params_(p)
                .map(ConformerAlignmentReport)
                .map_err(OperationError::from_rust)
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = allConformerBestRmsdsWithParams)]
    pub fn all_conformer_best_rmsds_with_params(
        &self,
        params: AllConformerRmsdParametersInput,
    ) -> Result<ConformerRmsds, JsValue> {
        self.all_conformer_best_rmsds(Some(params))
    }
    #[wasm_bindgen(js_name = withAlignmentToWithParams)]
    pub fn with_alignment_to_with_params(
        &self,
        reference: MoleculeInput,
        params: AlignmentParametersInput,
    ) -> Result<AlignedMolecule, JsValue> {
        self.with_alignment_to(reference, Some(params))
    }
    #[wasm_bindgen(js_name = withAlignedConformersWithParams)]
    pub fn with_aligned_conformers_with_params(
        &self,
        params: ConformerAlignmentParametersInput,
    ) -> Result<AlignedConformers, JsValue> {
        self.with_aligned_conformers(Some(params))
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=alignmentTransformToWithParams)]
    pub fn alignment_transform_to_with_params(
        &self,
        reference: MoleculeInput,
        params: AlignmentParametersInput,
    ) -> Result<AlignmentResult, JsValue> {
        self.alignment_transform_to(reference, Some(params))
    }
    #[wasm_bindgen(js_name=bestAlignmentToWithParams)]
    pub fn best_alignment_to_with_params(
        &self,
        reference: MoleculeInput,
        params: BestAlignmentParametersInput,
    ) -> Result<AlignmentResult, JsValue> {
        self.best_alignment_to(reference, Some(params))
    }
    #[wasm_bindgen(js_name=bestRmsdToWithParams)]
    pub fn best_rmsd_to_with_params(
        &self,
        reference: MoleculeInput,
        params: BestAlignmentParametersInput,
    ) -> Result<f64, JsValue> {
        self.best_rmsd_to(reference, Some(params))
    }
    #[wasm_bindgen(js_name=coordinateRmsdToWithParams)]
    pub fn coordinate_rmsd_to_with_params(
        &self,
        reference: MoleculeInput,
        params: CoordinateRmsdParametersInput,
    ) -> Result<f64, JsValue> {
        self.coordinate_rmsd_to(reference, Some(params))
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=alignToWithParams_)]
    pub fn align_to_with_params_(
        &mut self,
        reference: MoleculeInput,
        params: AlignmentParametersInput,
    ) -> Result<AlignmentResult, JsValue> {
        self.align_to(reference, Some(params))
    }
    #[wasm_bindgen(js_name=alignConformersWithParams_)]
    pub fn align_conformers_with_params_(
        &mut self,
        params: ConformerAlignmentParametersInput,
    ) -> Result<ConformerAlignmentReport, JsValue> {
        self.align_conformers(Some(params))
    }
}
