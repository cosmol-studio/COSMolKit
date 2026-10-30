//! Ordered native batch ownership and transformations.
use crate::{
    Molecule, SmilesParseParams, SmilesParseParamsInput,
    chemistry::*,
    coordinates::*,
    errors::{error_projection, exception},
    fields::{BoolArray, Field, IndexArray, TextArray},
    input,
    sanitize::*,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use std::error::Error;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
crate::enums::enum_value!(BatchErrorMode,"BatchErrorMode"{Strict="strict",KeepErrors="keep_errors"});
ck::__ck_configuration_BatchParams!(crate::project_parameters, ck::BatchParams, {}, {
    BatchParams
});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "Molecule | BatchError")]
    pub type BatchRecord;
    #[wasm_bindgen(typescript_type = "ReadonlyArray<BatchError>")]
    pub type BatchErrors;
    #[wasm_bindgen(typescript_type = "ReadonlyArray<Molecule | BatchError>")]
    pub type BatchRecords;
}
error_projection!(
    BatchError,
    ck::BatchError,
    {
        #[wasm_bindgen(getter)]
        pub fn index(&self) -> usize {
            self.inner.index
        }
        #[wasm_bindgen(getter)]
        pub fn operation(&self) -> String {
            self.inner.operation.into()
        }
        #[wasm_bindgen(getter)]
        pub fn message(&self) -> String {
            self.inner.message.clone()
        }
    },
    source = BatchCause,
    cause = |owner: &ReadonlyValue<ck::BatchError>| owner.source().map(|_| {
        let node = owner.child(|e| e.source().expect("immutable native batch source"));
        if node.is::<ck::SmilesError>() {
            return crate::smiles_errors::SmilesError::from_owner(
                node.child(|e| e.downcast_ref().expect("SMILES source")),
            )
            .unchecked_into();
        }
        if node.is::<ck::OperationError>() {
            return crate::operation_errors::OperationError::from_owner(
                node.child(|e| e.downcast_ref().expect("operation source")),
            )
            .unchecked_into();
        }
        if node.is::<ck::SdfReadError>() {
            return crate::sdf_errors::SdfReadError::from_owner(
                node.child(|e| e.downcast_ref().expect("SDF read source")),
            )
            .unchecked_into();
        }
        if node.is::<ck::MolPostError>() {
            return crate::mol_post_errors::MolPostError::from_owner(
                node.child(|e| e.downcast_ref().expect("molecule finalization source")),
            )
            .unchecked_into();
        }
        macro_rules! project {
            ($native:ty,$binding:path) => {
                if node.is::<$native>() {
                    return <$binding>::from_owner(
                        node.child(|e| e.downcast_ref().expect("checked batch source type")),
                    )
                    .unchecked_into();
                }
            };
        }
        project!(ck::SmilesWriteError, crate::smiles_errors::SmilesWriteError);
        #[cfg(feature = "cap-fingerprints")]
        project!(ck::MorganReadError, crate::morgan_errors::MorganReadError);
        #[cfg(feature = "cap-fingerprints")]
        project!(
            ck::AtomPairReadError,
            crate::pair_torsion::AtomPairReadError
        );
        #[cfg(feature = "cap-fingerprints")]
        project!(ck::AtomPairError, crate::morgan_errors::AtomPairError);
        #[cfg(feature = "cap-fingerprints")]
        project!(
            ck::TopologicalTorsionReadError,
            crate::pair_torsion::TopologicalTorsionReadError
        );
        #[cfg(feature = "cap-fingerprints")]
        project!(
            ck::LayeredFingerprintError,
            crate::path_fingerprints::LayeredFingerprintError
        );
        #[cfg(feature = "cap-fingerprints")]
        project!(
            ck::PatternFingerprintError,
            crate::pattern::PatternFingerprintError
        );
        #[cfg(feature = "cap-depict")]
        project!(ck::DrawingError, crate::drawing::DrawingError);
        #[cfg(feature = "cap-conformer")]
        project!(
            ck::ConformerRunError,
            crate::conformer_errors::ConformerRunError
        );
        // Scheduler thread-pool failures have external-library sources.
        crate::external_errors::standard_source(node).unchecked_into()
    })
);
fn errors(rows: &[ck::BatchError]) -> BatchErrors {
    let out = rows
        .iter()
        .cloned()
        .map(BatchError::from_rust)
        .collect::<js_sys::Array>();
    js_sys::Object::freeze(&out);
    out.unchecked_into()
}
error_projection!(BatchValidationError, ck::BatchValidationError, {
    #[wasm_bindgen(getter)]
    pub fn errors(&self) -> usize {
        self.inner.errors
    }
    #[wasm_bindgen(getter,js_name=recordErrors)]
    pub fn record_errors(&self) -> BatchErrors {
        errors(&self.inner.record_errors)
    }
});
#[wasm_bindgen]
pub struct MoleculeBatch(pub(crate) ck::MoleculeBatch);
#[wasm_bindgen(js_name = molsFromSmilesList)]
pub fn mols_from_smiles_list(
    smiles: Ts<TextArray>,
    parse: Option<SmilesParseParamsInput>,
    params: Option<BatchParamsInput>,
) -> Result<MoleculeBatch, JsValue> {
    MoleculeBatch::from_smiles_list(smiles, parse, params)
}
#[wasm_bindgen(js_name = molsFromRecords)]
pub fn mols_from_records(
    records: BatchRecords,
    errors: BatchErrorModeInput,
) -> Result<MoleculeBatch, JsValue> {
    MoleculeBatch::from_records(records, errors)
}
#[wasm_bindgen(js_name=molsFromSdfRecords)]
pub fn mols_from_sdf_records(
    text: Ts<input::TextInput>,
    read: Option<crate::molecular_io::SdfReadParamsInput>,
    errors: Option<BatchErrorModeInput>,
    n_jobs: Option<Ts<input::IndexInput>>,
    progress_bar: Option<Ts<input::BoolInput>>,
) -> Result<MoleculeBatch, JsValue> {
    MoleculeBatch::from_sdf_records(text, read, errors, n_jobs, progress_bar)
}
#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name=fromSdfRecords)]
    pub fn from_sdf_records(
        text: Ts<input::TextInput>,
        read: Option<crate::molecular_io::SdfReadParamsInput>,
        errors: Option<BatchErrorModeInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<Self, JsValue> {
        let text = input::read::<String>(text)?;
        let errors = errors
            .map(ck::BatchErrorMode::input)
            .transpose()?
            .unwrap_or_default();
        let n_jobs = n_jobs.map(input::read::<usize>).transpose()?;
        let progress_bar = progress_bar
            .map(input::read::<bool>)
            .transpose()?
            .unwrap_or(false);
        crate::molecular_io::SdfReadParams::with_input(read.as_ref(), |read| {
            ck::MoleculeBatch::from_sdf_records_with_params(
                &text,
                read,
                errors,
                n_jobs,
                progress_bar,
            )
            .map(Self)
            .map_err(BatchValidationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=fromSdfRecordsWithParams)]
    pub fn from_sdf_records_with_params(
        text: Ts<input::TextInput>,
        read: crate::molecular_io::SdfReadParamsInput,
        errors: Option<BatchErrorModeInput>,
        n_jobs: Option<Ts<input::IndexInput>>,
        progress_bar: Option<Ts<input::BoolInput>>,
    ) -> Result<Self, JsValue> {
        Self::from_sdf_records(text, Some(read), errors, n_jobs, progress_bar)
    }
    #[wasm_bindgen(js_name=fromSmilesListWithParams)]
    pub fn from_smiles_list_with_params(
        smiles: Ts<TextArray>,
        parse: SmilesParseParamsInput,
        params: BatchParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_smiles_list(smiles, Some(parse), Some(params))
    }
    #[wasm_bindgen(js_name=fromRecords)]
    pub fn from_records(
        records: BatchRecords,
        errors: BatchErrorModeInput,
    ) -> Result<Self, JsValue> {
        if !js_sys::Array::is_array(records.as_ref()) {
            return Err(js_sys::TypeError::new("expected an array of batch records").into());
        }
        let rows = records
            .unchecked_ref::<js_sys::Array>()
            .iter()
            .map(|row| {
                if let Ok(molecule) = ck::Molecule::decode(&row) {
                    Ok(ck::BatchRecord::Molecule(molecule))
                } else {
                    BatchError::with_native(&row, |error| ck::BatchRecord::Error(error.clone()))
                }
            })
            .collect::<Result<Vec<_>, JsValue>>()?;
        ck::MoleculeBatch::from_records(rows, ck::BatchErrorMode::input(errors)?)
            .map(Self)
            .map_err(BatchValidationError::from_rust)
    }
    #[wasm_bindgen(js_name=toList)]
    pub fn to_list(&self) -> crate::MoleculeArray {
        crate::MoleculeArray::from_native(self.0.to_list())
    }
    #[wasm_bindgen(getter)]
    pub fn records(&self) -> BatchRecords {
        let rows = self
            .0
            .records()
            .iter()
            .map(|row| match row {
                ck::BatchRecord::Molecule(molecule) => JsValue::from(Molecule(molecule.clone())),
                ck::BatchRecord::Error(error) => BatchError::from_rust(error.clone()),
            })
            .collect::<js_sys::Array>();
        js_sys::Object::freeze(&rows).unchecked_into()
    }
    #[wasm_bindgen(js_name=fromSmilesList)]
    pub fn from_smiles_list(
        smiles: Ts<TextArray>,
        parse: Option<SmilesParseParamsInput>,
        params: Option<BatchParamsInput>,
    ) -> Result<Self, JsValue> {
        let smiles = Vec::<String>::input(smiles)?;
        SmilesParseParams::with_input(parse.as_ref(), |parse| {
            BatchParams::with_input(params.as_ref(), |p| {
                ck::MoleculeBatch::from_smiles_list_with_params(&smiles, parse, p)
                    .map(Self)
                    .map_err(BatchValidationError::from_rust)
            })
        })
    }
    #[wasm_bindgen]
    pub fn copy(&self) -> Self {
        Self(self.0.clone())
    }
    #[wasm_bindgen(getter)]
    pub fn len(&self) -> usize {
        self.0.len()
    }
    #[wasm_bindgen(getter,js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    #[wasm_bindgen(getter,js_name=errorMode)]
    pub fn error_mode(&self) -> BatchErrorMode {
        self.0.error_mode().into()
    }
    #[wasm_bindgen]
    pub fn get(&self, index: Ts<input::IndexInput>) -> Result<BatchRecord, JsValue> {
        let index = input::read::<usize>(index)?;
        match self.0.get(index) {
            Some(ck::BatchRecord::Molecule(v)) => {
                Ok(JsValue::from(Molecule(v.clone())).unchecked_into())
            }
            Some(ck::BatchRecord::Error(v)) => {
                Ok(BatchError::from_rust(v.clone()).unchecked_into())
            }
            None => Err(js_sys::RangeError::new("batch index is out of range").into()),
        }
    }
    #[wasm_bindgen]
    pub fn select(&self, indices: Ts<IndexArray>) -> Result<Self, JsValue> {
        self.0
            .select(&Vec::<usize>::input(indices)?)
            .map(Self)
            .ok_or_else(|| js_sys::RangeError::new("batch index is out of range").into())
    }
    #[wasm_bindgen(js_name=validMask)]
    pub fn valid_mask(&self) -> Result<Ts<BoolArray>, JsValue> {
        Vec::<bool>::get(crate::storage::NativeField::owned(self.0.valid_mask()))
    }
    #[wasm_bindgen(js_name=invalidMask)]
    pub fn invalid_mask(&self) -> Result<Ts<BoolArray>, JsValue> {
        Vec::<bool>::get(crate::storage::NativeField::owned(self.0.invalid_mask()))
    }
    #[wasm_bindgen(getter,js_name=validCount)]
    pub fn valid_count(&self) -> usize {
        self.0.valid_count()
    }
    #[wasm_bindgen(getter,js_name=invalidCount)]
    pub fn invalid_count(&self) -> usize {
        self.0.invalid_count()
    }
    #[wasm_bindgen]
    pub fn errors(&self) -> BatchErrors {
        errors(&self.0.errors())
    }
    #[wasm_bindgen(js_name=withValidRecords)]
    pub fn with_valid_records(&self) -> Self {
        Self(self.0.with_valid_records())
    }
    #[wasm_bindgen(js_name=withParallelJobs)]
    pub fn with_parallel_jobs(&self, jobs: Option<Ts<input::IndexInput>>) -> Result<Self, JsValue> {
        let jobs = jobs.map(input::read::<usize>).transpose()?;
        self.0
            .clone()
            .with_parallel_jobs(jobs)
            .map(Self)
            .map_err(BatchValidationError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=parallelJobs)]
    pub fn parallel_jobs(&self) -> Option<usize> {
        self.0.parallel_jobs()
    }
    #[wasm_bindgen(js_name=withProgressBar)]
    pub fn with_progress_bar(
        &self,
        enabled: Option<Ts<input::BoolInput>>,
    ) -> Result<Self, JsValue> {
        Ok(Self(self.0.clone().with_progress_bar(
            enabled.map(input::read::<bool>).transpose()?,
        )))
    }
    #[wasm_bindgen(getter,js_name=progressBar)]
    pub fn progress_bar(&self) -> Option<bool> {
        self.0.progress_bar()
    }
    #[wasm_bindgen(js_name=sanitize)]
    pub fn sanitize(
        &self,
        options: Option<SanitizeParamsInput>,
        params: Option<BatchParamsInput>,
    ) -> Result<Self, JsValue> {
        SanitizeParams::with_input(options.as_ref(), |o| {
            BatchParams::with_input(params.as_ref(), |p| {
                self.0
                    .sanitize_with_params(o, p)
                    .map(Self)
                    .map_err(BatchValidationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=withHydrogens)]
    pub fn with_hydrogens(
        &self,
        options: Option<AddHsParamsInput>,
        params: Option<BatchParamsInput>,
    ) -> Result<Self, JsValue> {
        AddHsParams::with_input(options.as_ref(), |o| {
            BatchParams::with_input(params.as_ref(), |p| {
                self.0
                    .with_hydrogens_with_params(o, p)
                    .map(Self)
                    .map_err(BatchValidationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=withoutHydrogens)]
    pub fn without_hydrogens(
        &self,
        options: Option<RemoveHsParamsInput>,
        params: Option<BatchParamsInput>,
    ) -> Result<Self, JsValue> {
        RemoveHsParams::with_input(options.as_ref(), |o| {
            BatchParams::with_input(params.as_ref(), |p| {
                self.0
                    .without_hydrogens_with_params(o, p)
                    .map(Self)
                    .map_err(BatchValidationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=withKekulizedBonds)]
    pub fn with_kekulized_bonds(
        &self,
        options: Option<KekulizeParamsInput>,
        params: Option<BatchParamsInput>,
    ) -> Result<Self, JsValue> {
        KekulizeParams::with_input(options.as_ref(), |o| {
            BatchParams::with_input(params.as_ref(), |p| {
                self.0
                    .with_kekulized_bonds_with_params(o, p)
                    .map(Self)
                    .map_err(BatchValidationError::from_rust)
            })
        })
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=with2dCoordinates)]
    pub fn with_2d_coordinates(
        &self,
        options: Option<Coordinate2DParamsInput>,
        params: Option<BatchParamsInput>,
    ) -> Result<Self, JsValue> {
        Coordinate2DParams::with_input(options.as_ref(), |o| {
            BatchParams::with_input(params.as_ref(), |p| {
                self.0
                    .with_2d_coordinates_with_params(o, p)
                    .map(Self)
                    .map_err(BatchValidationError::from_rust)
            })
        })
    }
}

cosmolkit_tsify_macros::type_union! {
    BatchCause;
    "SmilesError | OperationError | SdfReadError | MolPostError | SmilesWriteError | IoError | ExternalError";
    #[cfg(feature = "cap-fingerprints")]
    "MorganReadError | AtomPairReadError | AtomPairError | TopologicalTorsionReadError | LayeredFingerprintError | PatternFingerprintError";
    #[cfg(feature = "cap-depict")]
    "DrawingError";
    #[cfg(feature = "cap-conformer")]
    "ConformerRunError";
}
#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name=sanitizeWithParams)]
    pub fn sanitize_with_params(
        &self,
        options: SanitizeParamsInput,
        params: BatchParamsInput,
    ) -> Result<Self, JsValue> {
        self.sanitize(Some(options), Some(params))
    }
    #[wasm_bindgen(js_name=withHydrogensWithParams)]
    pub fn with_hydrogens_with_params(
        &self,
        options: AddHsParamsInput,
        params: BatchParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_hydrogens(Some(options), Some(params))
    }
    #[wasm_bindgen(js_name=withoutHydrogensWithParams)]
    pub fn without_hydrogens_with_params(
        &self,
        options: RemoveHsParamsInput,
        params: BatchParamsInput,
    ) -> Result<Self, JsValue> {
        self.without_hydrogens(Some(options), Some(params))
    }
    #[wasm_bindgen(js_name=withKekulizedBondsWithParams)]
    pub fn with_kekulized_bonds_with_params(
        &self,
        options: KekulizeParamsInput,
        params: BatchParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_kekulized_bonds(Some(options), Some(params))
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=with2dCoordinatesWithParams)]
    pub fn with_2d_coordinates_with_params(
        &self,
        options: Coordinate2DParamsInput,
        params: BatchParamsInput,
    ) -> Result<Self, JsValue> {
        self.with_2d_coordinates(Some(options), Some(params))
    }
}
