//! AtomPair and torsion fingerprint calls share the native output boundary.
use crate::{
    Molecule,
    errors::{enum_error, error_projection, exception},
    fingerprint_output::{FingerprintAdditionalOutput, FingerprintAdditionalOutputInput},
    morgan_errors::child,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub struct AtomPairsParameters;
#[wasm_bindgen]
impl AtomPairsParameters {
    pub fn version() -> String {
        ck::AtomPairsParameters::version().into()
    }
    #[wasm_bindgen(js_name=atomTypes)]
    pub fn atom_types() -> Vec<u32> {
        ck::AtomPairsParameters::atom_types()
    }
}
macro_rules! encoding_value{($($method:ident,$js:ident);+$(;)?)=>{#[wasm_bindgen]impl AtomPairsParameters{$(
    #[wasm_bindgen(js_name=$js)]pub fn $method()->u32{ck::AtomPairsParameters::$method()}
)+}};}
encoding_value! {
    num_type_bits,numTypeBits;num_pi_bits,numPiBits;num_branch_bits,numBranchBits;
    num_chiral_bits,numChiralBits;code_size,codeSize;num_path_bits,numPathBits;
    max_path_length,maxPathLength;num_atom_pair_fingerprint_bits,numAtomPairFingerprintBits;
}
error_projection!(AllocationError, std::collections::TryReserveError, {});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "FingerprintPreparationError | AtomPairError")]
    #[derive(Clone)]
    pub type PairReadCause;
    #[wasm_bindgen(
        typescript_type = "FingerprintError | AtomPairError | PathError | MorganError | FingerprintJsonError | ThreadCountError | IoError | AllocationError"
    )]
    #[derive(Clone)]
    pub type TorsionCause;
    #[wasm_bindgen(typescript_type = "FingerprintPreparationError | TopologicalTorsionError")]
    #[derive(Clone)]
    pub type TorsionReadCause;
}
fn pair_read(owner: &ReadonlyValue<ck::AtomPairReadError>) -> Option<PairReadCause> {
    Some(match owner.as_ref() {
        ck::AtomPairReadError::Preparation(_) => child!(
            owner,
            ck::AtomPairReadError::Preparation,
            crate::morgan_errors::FingerprintPreparationError
        ),
        ck::AtomPairReadError::Generator(_) => child!(
            owner,
            ck::AtomPairReadError::Generator,
            crate::morgan_errors::AtomPairError
        ),
    })
}
enum_error!(AtomPairReadError,AtomPairReadErrorKind,ck::AtomPairReadError,{
    Preparation=>ck::AtomPairReadError::Preparation(_),Generator=>ck::AtomPairReadError::Generator(_)
},{#[wasm_bindgen(getter)] pub fn source(&self)->Option<PairReadCause>{pair_read(&self.inner)}},source=PairReadCause,cause=pair_read);
fn torsion(owner: &ReadonlyValue<ck::TopologicalTorsionError>) -> Option<TorsionCause> {
    use ck::TopologicalTorsionError as E;
    Some(match owner.as_ref() {
        E::Fingerprint(_) => child!(
            owner,
            E::Fingerprint,
            crate::fingerprint_errors::FingerprintError
        ),
        E::AtomInvariants(_) => child!(
            owner,
            E::AtomInvariants,
            crate::morgan_errors::AtomPairError
        ),
        E::Path(_) => child!(owner, E::Path, crate::path_errors::PathError),
        E::Preparation(_) => child!(owner, E::Preparation, crate::morgan_errors::MorganError),
        E::Json(_) => child!(owner, E::Json, crate::morgan_errors::FingerprintJsonError),
        E::ThreadCount(_) => child!(
            owner,
            E::ThreadCount,
            crate::worker_errors::ThreadCountError
        ),
        E::ThreadSpawn(_) => child!(owner, E::ThreadSpawn, crate::external_errors::IoError),
        E::OutputAllocation(_) => child!(owner, E::OutputAllocation, AllocationError),
        E::StatePoisoned | E::WorkerPanic | E::WorkerProtocol => return None,
    })
}
enum_error!(TopologicalTorsionError,TopologicalTorsionErrorKind,ck::TopologicalTorsionError,{
    Fingerprint=>ck::TopologicalTorsionError::Fingerprint(_),AtomInvariants=>ck::TopologicalTorsionError::AtomInvariants(_),Path=>ck::TopologicalTorsionError::Path(_),Preparation=>ck::TopologicalTorsionError::Preparation(_),Json=>ck::TopologicalTorsionError::Json(_),StatePoisoned=>ck::TopologicalTorsionError::StatePoisoned,ThreadCount=>ck::TopologicalTorsionError::ThreadCount(_),ThreadSpawn=>ck::TopologicalTorsionError::ThreadSpawn(_),WorkerPanic=>ck::TopologicalTorsionError::WorkerPanic,WorkerProtocol=>ck::TopologicalTorsionError::WorkerProtocol,OutputAllocation=>ck::TopologicalTorsionError::OutputAllocation(_)
},{#[wasm_bindgen(getter)] pub fn source(&self)->Option<TorsionCause>{torsion(&self.inner)}},source=TorsionCause,cause=torsion);
fn torsion_read(
    owner: &ReadonlyValue<ck::TopologicalTorsionReadError>,
) -> Option<TorsionReadCause> {
    Some(match owner.as_ref() {
        ck::TopologicalTorsionReadError::Preparation(_) => child!(
            owner,
            ck::TopologicalTorsionReadError::Preparation,
            crate::morgan_errors::FingerprintPreparationError
        ),
        ck::TopologicalTorsionReadError::Generator(_) => child!(
            owner,
            ck::TopologicalTorsionReadError::Generator,
            TopologicalTorsionError
        ),
    })
}
enum_error!(TopologicalTorsionReadError,TopologicalTorsionReadErrorKind,ck::TopologicalTorsionReadError,{
    Preparation=>ck::TopologicalTorsionReadError::Preparation(_),Generator=>ck::TopologicalTorsionReadError::Generator(_)
},{#[wasm_bindgen(getter)] pub fn source(&self)->Option<TorsionReadCause>{torsion_read(&self.inner)}},source=TorsionReadCause,cause=torsion_read);
ck::__ck_configuration_AtomPairGeneratorParams!(
    crate::project_parameters,
    ck::AtomPairGeneratorParams,
    {},
    { AtomPairGeneratorParams }
);
ck::__ck_configuration_AtomPairFingerprintParams!(
    crate::project_parameters,
    ck::AtomPairFingerprintParams,
    {},
    { AtomPairFingerprintParams }
);
ck::__ck_configuration_TopologicalTorsionGeneratorParams!(
    crate::project_parameters,
    ck::TopologicalTorsionGeneratorParams,
    {},
    { TopologicalTorsionGeneratorParams }
);
ck::__ck_configuration_TopologicalTorsionFingerprintParams!(
    crate::project_parameters,
    ck::TopologicalTorsionFingerprintParams,
    {},
    { TopologicalTorsionFingerprintParams }
);
macro_rules! call {
    ($method:ident,$js:ident,$native:ident,$params:ident,$input:ident,$error:ident,$output:path) => {
        #[wasm_bindgen]
        impl Molecule {
            #[wasm_bindgen(js_name=$js)]
            pub fn $method(
                &self,
                params: Option<$input>,
                additional_output: Option<FingerprintAdditionalOutputInput>,
            ) -> Result<$output, JsValue> {
                $params::with_input(params.as_ref(), |params| {
                    FingerprintAdditionalOutput::with_optional(
                        additional_output.as_ref(),
                        |output| {
                            self.0
                                .$native(params, output)
                                .map($output)
                                .map_err($error::from_rust)
                        },
                    )
                })
            }
        }
        paste::paste! {
            #[wasm_bindgen]
            impl Molecule {
                #[wasm_bindgen(js_name=[<$js WithParams>])]
                pub fn $native(&self, params: $input, additional_output: Option<FingerprintAdditionalOutputInput>) -> Result<$output, JsValue> {
                    self.$method(Some(params), additional_output)
                }
            }
        }
    };
}
call!(
    fingerprint_atom_pair,
    fingerprintAtomPair,
    fingerprint_atom_pair_with_params,
    AtomPairFingerprintParams,
    AtomPairFingerprintParamsInput,
    AtomPairReadError,
    crate::fingerprint_values::Fingerprint
);
call!(
    fingerprint_atom_pair_sparse,
    fingerprintAtomPairSparse,
    fingerprint_atom_pair_sparse_with_params,
    AtomPairFingerprintParams,
    AtomPairFingerprintParamsInput,
    AtomPairReadError,
    crate::sparse_fingerprints::SparseBitFingerprint
);
call!(
    fingerprint_atom_pair_count,
    fingerprintAtomPairCount,
    fingerprint_atom_pair_count_with_params,
    AtomPairFingerprintParams,
    AtomPairFingerprintParamsInput,
    AtomPairReadError,
    crate::sparse_fingerprints::SparseCountFingerprint32
);
call!(
    fingerprint_atom_pair_sparse_count,
    fingerprintAtomPairSparseCount,
    fingerprint_atom_pair_sparse_count_with_params,
    AtomPairFingerprintParams,
    AtomPairFingerprintParamsInput,
    AtomPairReadError,
    crate::sparse_fingerprints::SparseCountFingerprint
);
call!(
    fingerprint_topological_torsion,
    fingerprintTopologicalTorsion,
    fingerprint_topological_torsion_with_params,
    TopologicalTorsionFingerprintParams,
    TopologicalTorsionFingerprintParamsInput,
    TopologicalTorsionReadError,
    crate::fingerprint_values::Fingerprint
);
call!(
    fingerprint_topological_torsion_sparse,
    fingerprintTopologicalTorsionSparse,
    fingerprint_topological_torsion_sparse_with_params,
    TopologicalTorsionFingerprintParams,
    TopologicalTorsionFingerprintParamsInput,
    TopologicalTorsionReadError,
    crate::sparse_fingerprints::SparseBitFingerprint
);
call!(
    fingerprint_topological_torsion_count,
    fingerprintTopologicalTorsionCount,
    fingerprint_topological_torsion_count_with_params,
    TopologicalTorsionFingerprintParams,
    TopologicalTorsionFingerprintParamsInput,
    TopologicalTorsionReadError,
    crate::sparse_fingerprints::SparseCountFingerprint32
);
call!(
    fingerprint_topological_torsion_sparse_count,
    fingerprintTopologicalTorsionSparseCount,
    fingerprint_topological_torsion_sparse_count_with_params,
    TopologicalTorsionFingerprintParams,
    TopologicalTorsionFingerprintParamsInput,
    TopologicalTorsionReadError,
    crate::sparse_fingerprints::SparseCountFingerprint
);
#[wasm_bindgen]
impl AtomPairGeneratorParams {
    #[wasm_bindgen(js_name = infoString)]
    pub fn info_string(&self) -> Result<String, JsValue> {
        self.inner.read(ck::AtomPairGeneratorParams::info_string)
    }
    #[wasm_bindgen(js_name = toJson)]
    pub fn to_json(&self) -> Result<String, JsValue> {
        self.inner.read(ck::AtomPairGeneratorParams::to_json)
    }
    #[wasm_bindgen(js_name = withJson)]
    pub fn with_json(&self, json: tsify::Ts<crate::input::TextInput>) -> Result<Self, JsValue> {
        let json = crate::input::read::<String>(json)?;
        let native = self
            .inner
            .read(|value| value.with_json(&json))?
            .map_err(crate::morgan_errors::FingerprintJsonError::from_rust)?;
        Ok(Self {
            inner: crate::storage::NativeField::owned(native),
        })
    }
}

#[wasm_bindgen]
impl TopologicalTorsionGeneratorParams {
    #[wasm_bindgen(js_name = infoString)]
    pub fn info_string(&self) -> Result<String, JsValue> {
        self.inner
            .read(ck::TopologicalTorsionGeneratorParams::info_string)
    }
    #[wasm_bindgen(js_name = toJson)]
    pub fn to_json(&self) -> Result<String, JsValue> {
        self.inner
            .read(ck::TopologicalTorsionGeneratorParams::to_json)
    }
    #[wasm_bindgen(js_name = withJson)]
    pub fn with_json(&self, json: tsify::Ts<crate::input::TextInput>) -> Result<Self, JsValue> {
        let json = crate::input::read::<String>(json)?;
        let native = self
            .inner
            .read(|value| value.with_json(&json))?
            .map_err(crate::morgan_errors::FingerprintJsonError::from_rust)?;
        Ok(Self {
            inner: crate::storage::NativeField::owned(native),
        })
    }
}
