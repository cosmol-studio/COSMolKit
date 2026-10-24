//! Ordered native queries and their synchronous progress callback.
#[cfg(feature = "cap-fingerprints")]
use crate::morgan_reusable::{
    MorganAtomInvariantsGenerator, MorganAtomInvariantsGeneratorInput,
    MorganBondInvariantsGenerator, MorganBondInvariantsGeneratorInput,
};
use crate::{
    batch::{BatchValidationError, MoleculeBatch},
    callback,
    fields::Field,
    input,
    storage::NativeField,
};
#[cfg(feature = "cap-fingerprints")]
use crate::{morgan_configuration::*, pair_torsion::*, path_fingerprints::*, pattern::*};
use cosmolkit as ck;
use std::sync::Arc;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "(() => void)")]
    pub type ProgressCallback;
    #[wasm_bindgen(typescript_type = "ReadonlyArray<string | undefined>")]
    pub type BatchTexts;
    #[wasm_bindgen(
        typescript_type = "ReadonlyArray<ReadonlyArray<ReadonlyArray<number>> | undefined>"
    )]
    pub type BatchBoundsMatrices;
}
impl Field for Arc<dyn Fn() + Send + Sync> {
    type Input = ProgressCallback;
    type Output = ProgressCallback;
    fn input(value: Self::Input) -> Result<Self, JsValue> {
        Self::decode(value.as_ref())
    }
    fn decode(value: &JsValue) -> Result<Self, JsValue> {
        let handle = callback::Handle::new(value)?;
        let retained = Arc::clone(&handle);
        let native: Self = Arc::new(move || {
            if let Err(error) = retained.call(|| Ok(js_sys::Array::new())) {
                callback::record(error);
            }
        });
        handle.bind(&native);
        Ok(native)
    }
    fn get(field: NativeField<Self>) -> Result<Self::Output, JsValue> {
        field.read(callback::original)?.map(JsCast::unchecked_into)
    }
}
ck::__ck_configuration_BatchQueryParams!(crate::project_parameters, ck::BatchQueryParams, {}, {
    BatchQueryParams
});

fn query<R>(
    params: Option<&BatchQueryParamsInput>,
    n_jobs: Option<Ts<input::IndexInput>>,
    progress_bar: Option<Ts<input::BoolInput>>,
    call: impl FnOnce(&ck::BatchQueryParams) -> Result<R, ck::BatchValidationError>,
) -> Result<R, JsValue> {
    if params.is_some() && (n_jobs.is_some() || progress_bar.is_some()) {
        return Err(js_sys::TypeError::new(
            "query Params and explicit scheduling options are mutually exclusive",
        )
        .into());
    }
    let n_jobs = n_jobs.map(input::read::<usize>).transpose()?;
    let progress_bar = progress_bar.map(input::read::<bool>).transpose()?;
    callback::run(|| {
        if params.is_some() {
            BatchQueryParams::with_input(params, |p| {
                call(p).map_err(BatchValidationError::from_rust)
            })
        } else {
            call(&ck::BatchQueryParams {
                n_jobs,
                progress_bar,
                progress_callback: None,
            })
            .map_err(BatchValidationError::from_rust)
        }
    })
}

// Each declaration supplies only the canonical parameter/result association.
macro_rules! fingerprint {
    ($method:ident,$js:ident,$native:ident,$params:ident,$input:ident,$output:ty,$convert:expr) => {
        paste::paste! {
        #[wasm_bindgen]
        impl MoleculeBatch {
            #[wasm_bindgen(js_name=$js)]
            pub fn $method(
                &self,
                options: Option<$input>,
                runtime: Option<BatchQueryParamsInput>,
                n_jobs: Option<Ts<input::IndexInput>>,
                progress_bar: Option<Ts<input::BoolInput>>,
            ) -> Result<$output, JsValue> {
                $params::with_input(options.as_ref(), |p| {
                    let values = query(runtime.as_ref(), n_jobs, progress_bar, |r| {
                        self.0.$native(p, r)
                    })?;
                    ($convert)(values)
                })
            }
            #[wasm_bindgen(js_name=[<$js WithParams>])]
            pub fn [<$method _with_params>](&self,options:$input,runtime:BatchQueryParamsInput)->Result<$output,JsValue>{
                self.$method(Some(options),Some(runtime),None,None)
            }
        }
        }
    };
}
#[cfg(feature = "cap-fingerprints")]
fingerprint!(
    fingerprint_morgan,
    fingerprintMorgan,
    fingerprint_morgan_with_params,
    MorganFingerprintParams,
    MorganFingerprintParamsInput,
    crate::fingerprint_batch::FingerprintBatch,
    |v| Ok(crate::fingerprint_batch::FingerprintBatch(v))
);
#[cfg(feature = "cap-fingerprints")]
fingerprint!(
    fingerprint_atom_pair,
    fingerprintAtomPair,
    fingerprint_atom_pair_with_params,
    AtomPairFingerprintParams,
    AtomPairFingerprintParamsInput,
    crate::fingerprint_batch::FingerprintBatch,
    |v| Ok(crate::fingerprint_batch::FingerprintBatch(v))
);
#[cfg(feature = "cap-fingerprints")]
fingerprint!(
    fingerprint_topological_torsion,
    fingerprintTopologicalTorsion,
    fingerprint_topological_torsion_with_params,
    TopologicalTorsionFingerprintParams,
    TopologicalTorsionFingerprintParamsInput,
    crate::fingerprint_batch::FingerprintBatch,
    |v| Ok(crate::fingerprint_batch::FingerprintBatch(v))
);
#[cfg(feature = "cap-fingerprints")]
fingerprint!(
    fingerprint_layered,
    fingerprintLayered,
    fingerprint_layered_with_params,
    LayeredFingerprintParams,
    LayeredFingerprintParamsInput,
    crate::fingerprint_batch::FingerprintBatch,
    |v| Ok(crate::fingerprint_batch::FingerprintBatch(v))
);
#[cfg(feature = "cap-fingerprints")]
fingerprint!(
    fingerprint_pattern,
    fingerprintPattern,
    fingerprint_pattern_with_params,
    PatternFingerprintParams,
    PatternFingerprintParamsInput,
    crate::fingerprint_batch::FingerprintBatch,
    |v| Ok(crate::fingerprint_batch::FingerprintBatch(v))
);
#[cfg(feature = "cap-fingerprints")]
fingerprint!(
    fingerprint_atom_pair_count,
    fingerprintAtomPairCount,
    fingerprint_atom_pair_count_with_params,
    AtomPairFingerprintParams,
    AtomPairFingerprintParamsInput,
    crate::sparse_fingerprints::SparseCountFingerprint32Array,
    |v| Ok(crate::sparse_fingerprints::SparseCountFingerprint32Array::from_native(v))
);
#[cfg(feature = "cap-fingerprints")]
fingerprint!(
    fingerprint_atom_pair_sparse_count,
    fingerprintAtomPairSparseCount,
    fingerprint_atom_pair_sparse_count_with_params,
    AtomPairFingerprintParams,
    AtomPairFingerprintParamsInput,
    crate::sparse_fingerprints::SparseCountFingerprintArray,
    |v| Ok(crate::sparse_fingerprints::SparseCountFingerprintArray::from_native(v))
);
#[cfg(feature = "cap-fingerprints")]
fingerprint!(
    fingerprint_atom_pair_sparse,
    fingerprintAtomPairSparse,
    fingerprint_atom_pair_sparse_with_params,
    AtomPairFingerprintParams,
    AtomPairFingerprintParamsInput,
    crate::sparse_fingerprints::SparseBitFingerprintArray,
    |v| Ok(crate::sparse_fingerprints::SparseBitFingerprintArray::from_native(v))
);

#[wasm_bindgen]
impl MoleculeBatch {
    #[cfg(feature = "cap-fingerprints")]
    #[wasm_bindgen(js_name=fingerprintMorganWithGeneratorParams)]
    pub fn fingerprint_morgan_with_generator_params(
        &self,
        options: MorganGeneratorParamsInput,
        atom: Option<MorganAtomInvariantsGeneratorInput>,
        bond: Option<MorganBondInvariantsGeneratorInput>,
        call: Option<MorganCallParamsInput>,
        runtime: Option<BatchQueryParamsInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<crate::fingerprint_batch::FingerprintBatch, JsValue> {
        MorganGeneratorParams::with_input(Some(&options), |p| {
            MorganAtomInvariantsGenerator::with_optional_ref(atom.as_ref(), |a| {
                MorganBondInvariantsGenerator::with_optional_ref(bond.as_ref(), |b| {
                    MorganCallParams::with_input(call.as_ref(), |c| {
                        query(runtime.as_ref(), n_jobs, progress_bar, |r| {
                            self.0
                                .fingerprint_morgan_with_generator_params(p, a, b, c, r)
                        })
                        .map(crate::fingerprint_batch::FingerprintBatch)
                    })
                })
            })
        })
    }
    #[cfg(feature = "cap-fingerprints")]
    #[wasm_bindgen(js_name=fingerprintMorganWithOutputWithGeneratorParams)]
    pub fn fingerprint_morgan_with_output_with_generator_params(
        &self,
        options: MorganGeneratorParamsInput,
        atom: Option<MorganAtomInvariantsGeneratorInput>,
        bond: Option<MorganBondInvariantsGeneratorInput>,
        call: Option<MorganCallParamsInput>,
        additional_output: Option<Ts<input::BoolInput>>,
        runtime: Option<BatchQueryParamsInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<crate::batch_output::BatchFingerprintOutputs, JsValue> {
        let collect = additional_output
            .map(input::read::<bool>)
            .transpose()?
            .unwrap_or(true);
        MorganGeneratorParams::with_input(Some(&options), |p| {
            MorganAtomInvariantsGenerator::with_optional_ref(atom.as_ref(), |a| {
                MorganBondInvariantsGenerator::with_optional_ref(bond.as_ref(), |b| {
                    MorganCallParams::with_input(call.as_ref(), |c| {
                        query(runtime.as_ref(), n_jobs, progress_bar, |r| {
                            self.0.fingerprint_morgan_with_output_with_generator_params(
                                p, a, b, c, collect, r,
                            )
                        })
                        .map(crate::batch_output::BatchFingerprintOutputs::from_native)
                    })
                })
            })
        })
    }
    #[cfg(feature = "cap-fingerprints")]
    #[wasm_bindgen(js_name=fingerprintMorganWithOutputWithParams)]
    pub fn fingerprint_morgan_with_output_with_params(
        &self,
        options: MorganFingerprintParamsInput,
        additional_output: Ts<input::BoolInput>,
        runtime: BatchQueryParamsInput,
    ) -> Result<crate::batch_output::BatchFingerprintOutputs, JsValue> {
        self.fingerprint_morgan_with_output(
            Some(options),
            Some(additional_output),
            Some(runtime),
            None,
            None,
        )
    }
    #[cfg(feature = "cap-fingerprints")]
    #[wasm_bindgen(js_name=fingerprintAtomPairWithOutputWithParams)]
    pub fn fingerprint_atom_pair_with_output_with_params(
        &self,
        options: AtomPairFingerprintParamsInput,
        additional_output: Ts<input::BoolInput>,
        runtime: BatchQueryParamsInput,
    ) -> Result<crate::batch_output::BatchFingerprintOutputs, JsValue> {
        self.fingerprint_atom_pair_with_output(
            Some(options),
            Some(additional_output),
            Some(runtime),
            None,
            None,
        )
    }
    #[cfg(feature = "cap-fingerprints")]
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutputWithParams)]
    pub fn fingerprint_layered_with_output_with_params(
        &self,
        options: LayeredFingerprintParamsInput,
        runtime: BatchQueryParamsInput,
    ) -> Result<crate::path_fingerprints::LayeredFingerprintResults, JsValue> {
        self.fingerprint_layered_with_output(Some(options), Some(runtime), None, None)
    }
    #[wasm_bindgen(js_name=toSmilesWithParams)]
    pub fn to_smiles_with_params(
        &self,
        options: crate::SmilesWriteParamsInput,
        runtime: BatchQueryParamsInput,
    ) -> Result<BatchTexts, JsValue> {
        self.to_smiles(Some(options), Some(runtime), None, None)
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=toSvgWithParams)]
    pub fn to_svg_with_params(
        &self,
        width: Ts<input::U32Input>,
        height: Ts<input::U32Input>,
        runtime: BatchQueryParamsInput,
    ) -> Result<BatchTexts, JsValue> {
        self.to_svg(width, height, Some(runtime), None, None)
    }
    #[cfg(feature = "cap-conformer")]
    #[wasm_bindgen(js_name=dgBoundsMatrix)]
    pub fn dg_bounds_matrix(
        &self,
        runtime: Option<BatchQueryParamsInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<BatchBoundsMatrices, JsValue> {
        let values = query(runtime.as_ref(), n_jobs, progress_bar, |r| {
            self.0.dg_bounds_matrix_with_params(r)
        })?;
        let array = values
            .into_iter()
            .map(|v| {
                v.map(crate::conformer_calls::matrix_rows)
                    .map(JsValue::from)
                    .unwrap_or(JsValue::UNDEFINED)
            })
            .collect::<js_sys::Array>();
        Ok(js_sys::Object::freeze(&array).unchecked_into())
    }
    #[cfg(feature = "cap-conformer")]
    #[wasm_bindgen(js_name=dgBoundsMatrixWithParams)]
    pub fn dg_bounds_matrix_with_params(
        &self,
        runtime: BatchQueryParamsInput,
    ) -> Result<BatchBoundsMatrices, JsValue> {
        self.dg_bounds_matrix(Some(runtime), None, None)
    }
    #[cfg(feature = "cap-fingerprints")]
    #[wasm_bindgen(js_name=fingerprintMorganWithOutput)]
    pub fn fingerprint_morgan_with_output(
        &self,
        options: Option<MorganFingerprintParamsInput>,
        additional_output: Option<Ts<input::BoolInput>>,
        runtime: Option<BatchQueryParamsInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<crate::batch_output::BatchFingerprintOutputs, JsValue> {
        let collect = additional_output
            .map(input::read::<bool>)
            .transpose()?
            .unwrap_or(true);
        MorganFingerprintParams::with_input(options.as_ref(), |p| {
            query(runtime.as_ref(), n_jobs, progress_bar, |r| {
                self.0
                    .fingerprint_morgan_with_output_with_params(p, collect, r)
            })
            .map(crate::batch_output::BatchFingerprintOutputs::from_native)
        })
    }
    #[cfg(feature = "cap-fingerprints")]
    #[wasm_bindgen(js_name=fingerprintAtomPairWithOutput)]
    pub fn fingerprint_atom_pair_with_output(
        &self,
        options: Option<AtomPairFingerprintParamsInput>,
        additional_output: Option<Ts<input::BoolInput>>,
        runtime: Option<BatchQueryParamsInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<crate::batch_output::BatchFingerprintOutputs, JsValue> {
        let collect = additional_output
            .map(input::read::<bool>)
            .transpose()?
            .unwrap_or(true);
        AtomPairFingerprintParams::with_input(options.as_ref(), |p| {
            query(runtime.as_ref(), n_jobs, progress_bar, |r| {
                self.0
                    .fingerprint_atom_pair_with_output_with_params(p, collect, r)
            })
            .map(crate::batch_output::BatchFingerprintOutputs::from_native)
        })
    }
    #[cfg(feature = "cap-fingerprints")]
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutput)]
    pub fn fingerprint_layered_with_output(
        &self,
        options: Option<LayeredFingerprintParamsInput>,
        runtime: Option<BatchQueryParamsInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<LayeredFingerprintResults, JsValue> {
        LayeredFingerprintParams::with_input(options.as_ref(), |p| {
            query(runtime.as_ref(), n_jobs, progress_bar, |r| {
                self.0.fingerprint_layered_with_output_with_params(p, r)
            })
            .map(LayeredFingerprintResults::from_native)
        })
    }
    #[wasm_bindgen(js_name=toSmiles)]
    pub fn to_smiles(
        &self,
        options: Option<crate::SmilesWriteParamsInput>,
        runtime: Option<BatchQueryParamsInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<BatchTexts, JsValue> {
        crate::SmilesWriteParams::with_input(options.as_ref(), |p| {
            let values = query(runtime.as_ref(), n_jobs, progress_bar, |r| {
                self.0.to_smiles_with_params(p, r)
            })?;
            let array = values
                .into_iter()
                .map(|v| {
                    v.map(|v| crate::property_errors::text(&v).map(JsValue::from))
                        .unwrap_or(Ok(JsValue::UNDEFINED))
                })
                .collect::<Result<js_sys::Array, JsValue>>()?;
            Ok(js_sys::Object::freeze(&array).unchecked_into())
        })
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=toSvg)]
    pub fn to_svg(
        &self,
        width: Ts<input::U32Input>,
        height: Ts<input::U32Input>,
        runtime: Option<BatchQueryParamsInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<BatchTexts, JsValue> {
        let width = input::read::<u32>(width)?;
        let height = input::read::<u32>(height)?;
        let values = query(runtime.as_ref(), n_jobs, progress_bar, |r| {
            self.0.to_svg_with_params(width, height, r)
        })?;
        let array = values
            .into_iter()
            .map(|v| v.map(JsValue::from).unwrap_or(JsValue::UNDEFINED))
            .collect::<js_sys::Array>();
        Ok(js_sys::Object::freeze(&array).unchecked_into())
    }
}
