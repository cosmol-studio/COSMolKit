//! Avalon projection retains the native conversion and engine failures.
use crate::{
    Molecule,
    errors::{enum_error, exception},
    fingerprint_values::Fingerprint,
    input,
    morgan_errors::child,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
enum_error!(AvalonEngineError,AvalonEngineErrorKind,ck::AvalonEngineError,{
    InvalidArguments=>ck::AvalonEngineError::InvalidArguments{..},AvalonConversion=>ck::AvalonEngineError::AvalonConversion{..}
},{#[wasm_bindgen(getter)]pub fn reason(&self)->String{match self.inner.as_ref(){ck::AvalonEngineError::InvalidArguments{reason}=>(*reason).into(),ck::AvalonEngineError::AvalonConversion{reason}=>reason.clone()}}});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "MolecularIoError | AvalonEngineError")]
    #[derive(Clone)]
    pub type AvalonCause;
}
fn cause(owner: &ReadonlyValue<ck::AvalonFingerprintError>) -> Option<AvalonCause> {
    Some(match owner.as_ref() {
        ck::AvalonFingerprintError::Input(_) => child!(
            owner,
            ck::AvalonFingerprintError::Input,
            crate::molecular_io_errors::MolecularIoError
        ),
        ck::AvalonFingerprintError::Engine(_) => {
            child!(owner, ck::AvalonFingerprintError::Engine, AvalonEngineError)
        }
    })
}
enum_error!(AvalonFingerprintError,AvalonFingerprintErrorKind,ck::AvalonFingerprintError,{
    Input=>ck::AvalonFingerprintError::Input(_),Engine=>ck::AvalonFingerprintError::Engine(_)
},{#[wasm_bindgen(getter)]pub fn source(&self)->Option<AvalonCause>{cause(&self.inner)}},source=AvalonCause,cause=cause);
#[wasm_bindgen]
pub struct AvalonFingerprintFlags(pub(crate) ck::AvalonFingerprintFlags);
crate::fields::value_field!(
    ck::AvalonFingerprintFlags,
    AvalonFingerprintFlags,
    AvalonFingerprintFlagsInput,
    "AvalonFingerprintFlags"
);
#[wasm_bindgen]
impl AvalonFingerprintFlags {
    #[wasm_bindgen(js_name=fromBits)]
    pub fn from_bits(bits: tsify::Ts<input::U32Input>) -> Result<Option<Self>, JsValue> {
        Ok(ck::AvalonFingerprintFlags::from_bits(input::read(bits)?).map(Self))
    }
    #[wasm_bindgen(js_name=fromBitsRetain)]
    pub fn from_bits_retain(bits: tsify::Ts<input::U32Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::AvalonFingerprintFlags::from_bits_retain(
            input::read(bits)?,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn bits(&self) -> u32 {
        self.0.bits()
    }
    pub fn contains(&self, other: AvalonFingerprintFlagsInput) -> Result<bool, JsValue> {
        Self::with_ref(&other, |o| Ok(self.0.contains(*o)))
    }
    pub fn union(&self, other: AvalonFingerprintFlagsInput) -> Result<Self, JsValue> {
        Self::with_ref(&other, |o| Ok(Self(self.0 | *o)))
    }
}
macro_rules! flags{($($n:ident),*)=>{#[wasm_bindgen]impl AvalonFingerprintFlags{$(#[wasm_bindgen(getter,js_name=$n)]pub fn $n()->Self{Self(ck::AvalonFingerprintFlags::$n)})*}};}
#[allow(non_snake_case)]
const _: () = {
    flags!(
        RING_PATTERN,
        RING_PATH,
        ATOM_SYMBOL_PATH,
        ATOM_CLASS_PATH,
        ATOM_COUNT,
        AUGMENTED_ATOM,
        HCOUNT_PATH,
        HCOUNT_CLASS_PATH,
        HCOUNT_PAIR,
        BOND_PATH,
        AUGMENTED_BOND,
        RING_SIZE_COUNTS,
        DEGREE_PATH,
        CLASS_SPIDERS,
        FEATURE_PAIRS,
        SCAFFOLD_IDS,
        SCAFFOLD_COLORS,
        SCAFFOLD_LINKS,
        SHORTCUT_LABELS,
        ALL_FEATURES,
        NON_SSS_BITS,
        SIMILARITY
    );
};
ck::__ck_configuration_AvalonFingerprintParams!(
    crate::project_parameters,
    ck::AvalonFingerprintParams,
    {},
    { AvalonFingerprintParams }
);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintAvalon)]
    pub fn fingerprint_avalon(
        &self,
        params: Option<AvalonFingerprintParamsInput>,
    ) -> Result<Fingerprint, JsValue> {
        AvalonFingerprintParams::with_input(params.as_ref(), |p| {
            self.0
                .fingerprint_avalon_with_params(p)
                .map(Fingerprint)
                .map_err(AvalonFingerprintError::from_rust)
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = fingerprintAvalonWithParams)]
    pub fn fingerprint_avalon_with_params(
        &self,
        params: AvalonFingerprintParamsInput,
    ) -> Result<Fingerprint, JsValue> {
        self.fingerprint_avalon(Some(params))
    }
}
