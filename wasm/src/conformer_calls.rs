//! Embedding calls borrow canonical configuration, including its native state.
use crate::conformer_errors::{ConformerError, ConformerRunError};
use crate::{Molecule, fields::Field, input, storage::NativeField};
use cosmolkit as ck;
use serde::Deserialize;
use tsify::{Ts, Tsify};
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ReadonlyArray<ReadonlyArray<number>>")]
    pub type BoundsMatrixRows;
}

ck::__ck_configuration_EmbedParams!(crate::project_parameters, ck::EmbedParams, {}, {
    EmbedParams
});
crate::enums::enum_value!(EmbedTermination,"EmbedTermination"{Completed="completed",TimedOut="timed_out",Interrupted="interrupted"});
crate::enums::enum_value!(EmbedFailureCause,"EmbedFailureCause"{
    InitialCoords="initial_coords",FirstMinimization="first_minimization",CheckTetrahedralCenters="check_tetrahedral_centers",CheckChiralCenters="check_chiral_centers",MinimizeFourthDimension="minimize_fourth_dimension",EtkMinimization="etk_minimization",FinalChiralBounds="final_chiral_bounds",FinalCenterInVolume="final_center_in_volume",LinearDoubleBond="linear_double_bond",BadDoubleBondStereo="bad_double_bond_stereo",CheckChiralCenters2="check_chiral_centers2",ExceededTimeout="exceeded_timeout",Minimization="minimization",KtermViolation="kterm_violation",Clash="clash",EndOfEnum="end_of_enum"
});
#[derive(Deserialize, Tsify)]
#[serde(transparent)]
pub struct BoundsRows(#[tsify(type = "ReadonlyArray<ReadonlyArray<number>>")] Vec<Vec<f64>>);
pub(crate) fn matrix_rows(rows: Vec<Vec<f64>>) -> BoundsMatrixRows {
    let array = rows
        .into_iter()
        .map(|row| {
            js_sys::Object::freeze(
                &row.into_iter()
                    .map(JsValue::from)
                    .collect::<js_sys::Array>(),
            )
        })
        .collect::<js_sys::Array>();
    js_sys::Object::freeze(&array).unchecked_into()
}
#[wasm_bindgen]
impl EmbedParams {
    pub fn dg() -> Self {
        Self {
            inner: NativeField::owned(ck::EmbedParams::dg()),
        }
    }
    pub fn kdg() -> Self {
        Self {
            inner: NativeField::owned(ck::EmbedParams::kdg()),
        }
    }
    pub fn etdg() -> Self {
        Self {
            inner: NativeField::owned(ck::EmbedParams::etdg()),
        }
    }
    #[wasm_bindgen(js_name=etdgV2)]
    pub fn etdg_v2() -> Self {
        Self {
            inner: NativeField::owned(ck::EmbedParams::etdg_v2()),
        }
    }
    pub fn etkdg() -> Self {
        Self {
            inner: NativeField::owned(ck::EmbedParams::etkdg()),
        }
    }
    #[wasm_bindgen(js_name=etkdgV2)]
    pub fn etkdg_v2() -> Self {
        Self {
            inner: NativeField::owned(ck::EmbedParams::etkdg_v2()),
        }
    }
    #[wasm_bindgen(js_name=etkdgV3)]
    pub fn etkdg_v3() -> Self {
        Self {
            inner: NativeField::owned(ck::EmbedParams::etkdg_v3()),
        }
    }
    #[wasm_bindgen(js_name=srEtkdgV3)]
    pub fn sr_etkdg_v3() -> Self {
        Self {
            inner: NativeField::owned(ck::EmbedParams::sr_etkdg_v3()),
        }
    }
    #[wasm_bindgen(js_name=withJson)]
    pub fn with_json(&self, json: Ts<input::TextInput>) -> Result<Self, JsValue> {
        let json = input::read::<String>(json)?;
        self.inner.read(|p| {
            p.with_json(&json)
                .map(|p| Self {
                    inner: NativeField::owned(p),
                })
                .map_err(ConformerError::from_rust)
        })?
    }
    #[wasm_bindgen(js_name=toJson)]
    pub fn to_json(&self) -> Result<String, JsValue> {
        self.inner.read(ck::EmbedParams::to_json)
    }
    #[wasm_bindgen(js_name=setBoundsMatrix)]
    pub fn set_bounds_matrix(&self, rows: Ts<BoundsRows>) -> Result<(), JsValue> {
        let rows = input::decode(rows)?.0;
        self.inner
            .write(|p| p.set_bounds_matrix(rows).map_err(ConformerError::from_rust))?
    }
    #[wasm_bindgen(getter,js_name=hasBoundsMatrix)]
    pub fn has_bounds_matrix(&self) -> Result<bool, JsValue> {
        self.inner.read(ck::EmbedParams::has_bounds_matrix)
    }
    #[wasm_bindgen(getter)]
    pub fn failures(&self) -> Result<Ts<crate::fields::U32Array>, JsValue> {
        self.inner
            .read(|p| Vec::<u32>::get(NativeField::owned(p.failures.clone())))?
    }
}

fn options<R>(
    params: Option<&EmbedParamsInput>,
    call: impl FnOnce(&ck::EmbedParams) -> Result<R, JsValue>,
) -> Result<R, JsValue> {
    match params {
        Some(p) if input::is_plain_options(p.as_ref()) => call(&EmbedParams::options_on(
            ck::EmbedParams::etkdg_v3(),
            Some(p.as_ref()),
        )?),
        Some(p) => EmbedParams::with_input(Some(p), call),
        None => call(&ck::EmbedParams::etkdg_v3()),
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(getter,js_name=num3dConformers)]
    pub fn num_3d_conformers(&self) -> usize {
        self.0.num_3d_conformers()
    }
    #[wasm_bindgen(js_name=dgBoundsMatrix)]
    pub fn dg_bounds_matrix(&self) -> Result<BoundsMatrixRows, JsValue> {
        self.0
            .dg_bounds_matrix()
            .map(matrix_rows)
            .map_err(ConformerRunError::from_rust)
    }
    #[wasm_bindgen(js_name=with3dConformer)]
    pub fn with_3d_conformer(&self, params: Option<EmbedParamsInput>) -> Result<Molecule, JsValue> {
        options(params.as_ref(), |p| {
            self.0
                .with_3d_conformer_with_params(p)
                .map(Molecule)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=embed3dConformer_)]
    pub fn embed_3d_conformer(&mut self, params: Option<EmbedParamsInput>) -> Result<(), JsValue> {
        options(params.as_ref(), |p| {
            self.0
                .embed_3d_conformer_with_params_(p)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=with3dConformerResult)]
    pub fn with_3d_conformer_result(
        &self,
        params: Option<EmbedParamsInput>,
    ) -> Result<EmbedMoleculeResult, JsValue> {
        options(params.as_ref(), |p| {
            self.0
                .with_3d_conformer_result_with_params(p)
                .map(EmbedMoleculeResult)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=embed3dConformerResult_)]
    pub fn embed_3d_conformer_result(
        &mut self,
        params: Option<EmbedParamsInput>,
    ) -> Result<EmbedMoleculeResult, JsValue> {
        options(params.as_ref(), |p| {
            self.0
                .embed_3d_conformer_result_with_params_(p)
                .map(EmbedMoleculeResult)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=with3dConformers)]
    pub fn with_3d_conformers(
        &self,
        count: Ts<input::U32Input>,
        params: Option<EmbedParamsInput>,
    ) -> Result<Molecule, JsValue> {
        let count = input::read::<u32>(count)?;
        options(params.as_ref(), |p| {
            self.0
                .with_3d_conformers_with_params(count, p)
                .map(Molecule)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=embed3dConformers_)]
    pub fn embed_3d_conformers(
        &mut self,
        count: Ts<input::U32Input>,
        params: Option<EmbedParamsInput>,
    ) -> Result<(), JsValue> {
        let count = input::read::<u32>(count)?;
        options(params.as_ref(), |p| {
            self.0
                .embed_3d_conformers_with_params_(count, p)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=with3dConformersResult)]
    pub fn with_3d_conformers_result(
        &self,
        count: Ts<input::U32Input>,
        params: Option<EmbedParamsInput>,
    ) -> Result<EmbedMultipleConfsResult, JsValue> {
        let count = input::read::<u32>(count)?;
        options(params.as_ref(), |p| {
            self.0
                .with_3d_conformers_result_with_params(count, p)
                .map(EmbedMultipleConfsResult)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=embed3dConformersResult_)]
    pub fn embed_3d_conformers_result(
        &mut self,
        count: Ts<input::U32Input>,
        params: Option<EmbedParamsInput>,
    ) -> Result<EmbedMultipleConfsResult, JsValue> {
        let count = input::read::<u32>(count)?;
        options(params.as_ref(), |p| {
            self.0
                .embed_3d_conformers_result_with_params_(count, p)
                .map(EmbedMultipleConfsResult)
                .map_err(crate::operation_errors::OperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=with3dConformerWithParams)]
    pub fn with_3d_conformer_with_params(
        &self,
        params: EmbedParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.with_3d_conformer(Some(params))
    }
    #[wasm_bindgen(js_name=with3dConformerResultWithParams)]
    pub fn with_3d_conformer_result_with_params(
        &self,
        params: EmbedParamsInput,
    ) -> Result<EmbedMoleculeResult, JsValue> {
        self.with_3d_conformer_result(Some(params))
    }
    #[wasm_bindgen(js_name=embed3dConformerWithParams_)]
    pub fn embed_3d_conformer_with_params(
        &mut self,
        params: EmbedParamsInput,
    ) -> Result<(), JsValue> {
        self.embed_3d_conformer(Some(params))
    }
    #[wasm_bindgen(js_name=embed3dConformerResultWithParams_)]
    pub fn embed_3d_conformer_result_with_params(
        &mut self,
        params: EmbedParamsInput,
    ) -> Result<EmbedMoleculeResult, JsValue> {
        self.embed_3d_conformer_result(Some(params))
    }
    #[wasm_bindgen(js_name=with3dConformersWithParams)]
    pub fn with_3d_conformers_with_params(
        &self,
        count: Ts<input::U32Input>,
        params: EmbedParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.with_3d_conformers(count, Some(params))
    }
    #[wasm_bindgen(js_name=with3dConformersResultWithParams)]
    pub fn with_3d_conformers_result_with_params(
        &self,
        count: Ts<input::U32Input>,
        params: EmbedParamsInput,
    ) -> Result<EmbedMultipleConfsResult, JsValue> {
        self.with_3d_conformers_result(count, Some(params))
    }
    #[wasm_bindgen(js_name=embed3dConformersWithParams_)]
    pub fn embed_3d_conformers_with_params(
        &mut self,
        count: Ts<input::U32Input>,
        params: EmbedParamsInput,
    ) -> Result<(), JsValue> {
        self.embed_3d_conformers(count, Some(params))
    }
    #[wasm_bindgen(js_name=embed3dConformersResultWithParams_)]
    pub fn embed_3d_conformers_result_with_params(
        &mut self,
        count: Ts<input::U32Input>,
        params: EmbedParamsInput,
    ) -> Result<EmbedMultipleConfsResult, JsValue> {
        self.embed_3d_conformers_result(count, Some(params))
    }
}
#[wasm_bindgen]
pub struct EmbedMoleculeResult(ck::EmbedMoleculeResult);
#[wasm_bindgen]
pub struct EmbedMultipleConfsResult(ck::EmbedMultipleConfsResult);
macro_rules! result_fields {
    ($ty:ident) => {
        #[wasm_bindgen]
        impl $ty {
            #[wasm_bindgen(getter)]
            pub fn molecule(&self) -> Molecule {
                Molecule(self.0.molecule.clone())
            }
            #[wasm_bindgen(getter)]
            pub fn params(&self) -> EmbedParams {
                EmbedParams {
                    inner: NativeField::owned(self.0.params.clone()),
                }
            }
            #[wasm_bindgen(getter)]
            pub fn termination(&self) -> EmbedTermination {
                self.0.termination.into()
            }
            #[wasm_bindgen(getter,js_name=failuresTracked)]
            pub fn failures_tracked(&self) -> bool {
                self.0.failures_tracked()
            }
            #[wasm_bindgen(js_name=failureCounts)]
            pub fn failure_counts(&self) -> Result<Option<Ts<crate::fields::U32Array>>, JsValue> {
                self.0
                    .failure_counts()
                    .map(|v| Vec::<u32>::get(NativeField::owned(v.to_vec())))
                    .transpose()
            }
            #[wasm_bindgen(js_name=failureCount)]
            pub fn failure_count(
                &self,
                cause: EmbedFailureCauseInput,
            ) -> Result<Option<u32>, JsValue> {
                self.0
                    .failure_count(ck::EmbedFailureCause::input(cause)?)
                    .map_err(EmbedFailureCountError::from_rust)
            }
        }
    };
}
result_fields!(EmbedMoleculeResult);
result_fields!(EmbedMultipleConfsResult);
#[wasm_bindgen]
impl EmbedMoleculeResult {
    #[wasm_bindgen(getter,js_name=confId)]
    pub fn conf_id(&self) -> Option<u32> {
        self.0.conf_id
    }
    pub fn ok(&self) -> bool {
        self.0.ok()
    }
}
#[wasm_bindgen]
impl EmbedMultipleConfsResult {
    #[wasm_bindgen(getter,js_name=confIds)]
    pub fn conf_ids(&self) -> Result<Ts<crate::fields::U32Array>, JsValue> {
        Vec::<u32>::get(NativeField::owned(self.0.conf_ids.clone()))
    }
    #[wasm_bindgen(getter,js_name=requestedNumConfs)]
    pub fn requested_num_confs(&self) -> u32 {
        self.0.requested_num_confs
    }
    #[wasm_bindgen(getter,js_name=generatedCount)]
    pub fn generated_count(&self) -> usize {
        self.0.generated_count()
    }
}
crate::errors::enum_error!(EmbedFailureCountError,EmbedFailureCountErrorKind,ck::EmbedFailureCountError,{InvalidCause=>ck::EmbedFailureCountError::InvalidCause,MissingCounter=>ck::EmbedFailureCountError::MissingCounter{..}},{
    #[wasm_bindgen(getter)]pub fn index(&self)->Option<usize>{match *self.inner{ck::EmbedFailureCountError::MissingCounter{index,..}=>Some(index),_=>None}}
    #[wasm_bindgen(getter)]pub fn available(&self)->Option<usize>{match *self.inner{ck::EmbedFailureCountError::MissingCounter{available,..}=>Some(available),_=>None}}
});
