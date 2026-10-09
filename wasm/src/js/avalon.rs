//! Avalon configuration/result/error transport; no fingerprint chemistry.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use crate::fingerprint_values::Fingerprint;
use crate::host_values::{bool_value, u32_value};
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct AvalonFingerprintFlags { inner: ck::AvalonFingerprintFlags }
#[wasm_bindgen]
impl AvalonFingerprintFlags {
    #[wasm_bindgen(js_name=fromBitsRetain)]
    pub fn from_bits_retain(#[wasm_bindgen(unchecked_param_type="number")] bits: JsValue) -> Result<Self, JsValue> {
        Ok(Self { inner: ck::AvalonFingerprintFlags::from_bits_retain(u32_value(&bits, "bits")?) })
    }
    pub fn bits(&self) -> u32 { self.inner.bits() }
    #[wasm_bindgen(js_name=toString)]
    pub fn repr(&self) -> String { format!("AvalonFingerprintFlags({:#x})", self.inner.bits()) }
}
#[wasm_bindgen]
pub struct AvalonFingerprintParams { inner: ck::AvalonFingerprintParams }
#[cosmolkit_wasm::javascript_options("AvalonFingerprintOptions")]
#[wasm_bindgen]
impl AvalonFingerprintParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type="number")] n_bits: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="boolean")] is_query: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="number")] bit_flags: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::AvalonFingerprintParams::default();
        if !n_bits.is_undefined() { inner.n_bits = u32_value(&n_bits, "nBits")?; }
        if !is_query.is_undefined() { inner.is_query = bool_value(&is_query, "isQuery")?; }
        if !bit_flags.is_undefined() { inner.bit_flags = ck::AvalonFingerprintFlags::from_bits_retain(u32_value(&bit_flags, "bitFlags")?); }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=nBits)]
    pub fn n_bits(&self) -> u32 { self.inner.n_bits }
    #[wasm_bindgen(setter,js_name=nBits)]
    pub fn set_n_bits(&mut self, #[wasm_bindgen(unchecked_param_type="number")] value: JsValue) -> Result<(), JsValue> {
        self.inner.n_bits = u32_value(&value, "nBits")?; Ok(())
    }
    #[wasm_bindgen(getter,js_name=isQuery)]
    pub fn is_query(&self) -> bool { self.inner.is_query }
    #[wasm_bindgen(setter,js_name=isQuery)]
    pub fn set_is_query(&mut self, #[wasm_bindgen(unchecked_param_type="boolean")] value: JsValue) -> Result<(), JsValue> {
        self.inner.is_query = bool_value(&value, "isQuery")?; Ok(())
    }
    #[wasm_bindgen(getter,js_name=bitFlags)]
    pub fn bit_flags(&self) -> u32 { self.inner.bit_flags.bits() }
    #[wasm_bindgen(setter,js_name=bitFlags)]
    pub fn set_bit_flags(&mut self, #[wasm_bindgen(unchecked_param_type="number")] value: JsValue) -> Result<(), JsValue> {
        self.inner.bit_flags = ck::AvalonFingerprintFlags::from_bits_retain(u32_value(&value, "bitFlags")?); Ok(())
    }
    #[wasm_bindgen(js_name=toString)]
    pub fn repr(&self) -> String {
        format!("AvalonFingerprintParams(nBits={}, isQuery={}, bitFlags={:#x})", self.inner.n_bits, self.inner.is_query, self.inner.bit_flags.bits())
    }
}

#[wasm_bindgen(inline_js="export function visitAvalonParams(value,visit){visit(value);}")]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitAvalonParams)]
    fn visit(value: &JsValue, callback: &mut dyn FnMut(&AvalonFingerprintParams)) -> Result<(), JsValue>;
}
impl AvalonFingerprintParams {
    fn from_configuration(value: &JsValue) -> Result<Self, JsValue> {
        let mut inner = None;
        if visit(value, &mut |params: &Self| inner = Some(params.inner.clone())).is_ok() {
            return inner.map(|inner| Self { inner }).ok_or_else(|| js_sys::TypeError::new("invalid AvalonFingerprintParams").into());
        }
        Self::from_js_options(value)
    }
}
#[wasm_bindgen]
pub struct AvalonEngineError { kind: String, reason: String }
#[wasm_bindgen]
impl AvalonEngineError {
    #[wasm_bindgen(getter)] pub fn domain(&self) -> String { "Avalon".into() }
    #[wasm_bindgen(getter)] pub fn kind(&self) -> String { self.kind.clone() }
    #[wasm_bindgen(getter)] pub fn reason(&self) -> String { self.reason.clone() }
}
#[wasm_bindgen]
pub struct AvalonFingerprintError { kind: String, message: String, cause: JsValue }
#[wasm_bindgen]
impl AvalonFingerprintError {
    #[wasm_bindgen(getter)] pub fn domain(&self) -> String { "Fingerprint".into() }
    #[wasm_bindgen(getter)] pub fn kind(&self) -> String { self.kind.clone() }
    #[wasm_bindgen(getter)] pub fn message(&self) -> String { self.message.clone() }
    #[wasm_bindgen(getter,unchecked_return_type="Error | null")] pub fn cause(&self) -> JsValue { self.cause.clone() }
}
fn error(source: ck::AvalonFingerprintError) -> JsValue {
    let convert = || -> Result<JsValue, JsValue> {
        let (kind, cause) = match &source {
            ck::AvalonFingerprintError::Input(e) => ("Input", source_error(e)?),
            ck::AvalonFingerprintError::Engine(e) => {
                let (kind, reason) = match e {
                    ck::AvalonEngineError::InvalidArguments { reason } => ("InvalidArguments", *reason),
                    ck::AvalonEngineError::AvalonConversion { reason } => ("Conversion", reason.as_str()),
                };
                let cause: JsValue = js_sys::Error::new(&e.to_string()).into();
                set(&cause, "name", "AvalonEngineError".into())?;
                set(&cause, "domain", "Avalon".into())?;
                set(&cause, "kind", kind.into())?;
                set(&cause, "reason", reason.into())?;
                set(&cause, "detail", AvalonEngineError { kind: kind.into(), reason: reason.into() }.into())?;
                (kind, cause)
            }
        };
        let result: JsValue = js_sys::Error::new(&source.to_string()).into();
        set(&result, "name", "AvalonFingerprintError".into())?;
        set(&result, "domain", "Fingerprint".into())?;
        set(&result, "kind", kind.into())?;
        set(&result, "cause", cause.clone())?;
        set(&result, "detail", AvalonFingerprintError { kind: kind.into(), message: source.to_string(), cause }.into())?;
        Ok(result)
    };
    convert().unwrap_or_else(|error| error)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintAvalon)]
    pub fn fingerprint_avalon(&self,
        #[wasm_bindgen(unchecked_optional_param_type="AvalonFingerprintParams | AvalonFingerprintOptions")] params: JsValue,
    ) -> Result<Fingerprint, JsValue> {
        let configured = (!params.is_undefined()).then(|| AvalonFingerprintParams::from_configuration(&params)).transpose()?;
        configured.as_ref().map_or_else(|| self.inner.fingerprint_avalon(),
            |params| self.inner.fingerprint_avalon_with_params(&params.inner))
            .map(|inner| Fingerprint { inner }).map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintAvalonWithParams)]
    pub fn fingerprint_avalon_with_params(&self, params: &AvalonFingerprintParams) -> Result<Fingerprint, JsValue> {
        self.inner.fingerprint_avalon_with_params(&params.inner).map(|inner| Fingerprint { inner }).map_err(error)
    }
}
