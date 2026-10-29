//! Complete structural read error categories, source phases and causal contexts.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub enum BioPdbReadStage {
    Stream = 0,
    Record = 1,
    Finalization = 2,
}
#[wasm_bindgen]
pub enum BioMmcifReadStage {
    CifDocument = 0,
    CoordinateBlock = 1,
    CrystalCell = 2,
    Refinement = 3,
    Tls = 4,
    Experimental = 5,
    Reflections = 6,
    Software = 7,
    Ncs = 8,
    FractionalTransform = 9,
    Origx = 10,
    AnisotropicU = 11,
    AtomSites = 12,
    EntitySequence = 13,
    Helices = 14,
    Sheets = 15,
    Connections = 16,
    CisPeptides = 17,
    ModifiedResidues = 18,
    Assemblies = 19,
    SiftsUnp = 20,
    CcdRestoration = 21,
    Materialization = 22,
    StructureValidation = 23,
}
#[wasm_bindgen]
pub struct BioReadError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl BioReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
#[wasm_bindgen]
pub struct BioPdbReadError {
    kind: String,
    message: String,
    cause: JsValue,
    stage: u32,
    line: Option<i32>,
    tag: Option<[u8; 4]>,
}
#[wasm_bindgen]
impl BioPdbReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(unchecked_return_type = "BioPdbReadStage")]
    pub fn stage(&self) -> u32 {
        self.stage
    }
    #[wasm_bindgen(js_name=lineNumber,unchecked_return_type="number | null")]
    pub fn line_number(&self) -> JsValue {
        self.line.map_or(JsValue::NULL, JsValue::from)
    }
    #[wasm_bindgen(js_name=recordTag,unchecked_return_type="Uint8Array | null")]
    pub fn record_tag(&self) -> JsValue {
        self.tag.map_or(JsValue::NULL, |v| {
            js_sys::Uint8Array::from(v.as_slice()).into()
        })
    }
}
#[wasm_bindgen]
pub struct BioMmcifReadError {
    kind: String,
    message: String,
    cause: JsValue,
    stage: u32,
}
#[wasm_bindgen]
impl BioMmcifReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(unchecked_return_type = "BioMmcifReadStage")]
    pub fn stage(&self) -> u32 {
        self.stage
    }
}
fn annotated(
    name: &str,
    kind: &str,
    message: &str,
    cause: &JsValue,
    detail: JsValue,
) -> Result<JsValue, JsValue> {
    let e = js_sys::Error::new(message);
    e.set_name(name);
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    set(&e, "detail", detail)?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    Ok(e)
}
pub(crate) fn pdb_error(source: &ck::BioPdbReadError) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let detail = BioPdbReadError {
        kind: "Pdb".into(),
        message: source.to_string(),
        cause: cause.clone(),
        stage: source.stage() as u32,
        line: source.line_number(),
        tag: source.record_tag(),
    };
    let line = detail.line_number();
    let tag = detail.record_tag();
    let error = annotated(
        "BioPdbReadError",
        "Pdb",
        &source.to_string(),
        &cause,
        detail.into(),
    )?;
    set(&error, "stage", (source.stage() as u32).into())?;
    set(&error, "lineNumber", line)?;
    set(&error, "recordTag", tag)?;
    Ok(error)
}
pub(crate) fn mmcif_error(source: &ck::BioMmcifReadError) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let detail = BioMmcifReadError {
        kind: "Mmcif".into(),
        message: source.to_string(),
        cause: cause.clone(),
        stage: source.stage() as u32,
    };
    let error = annotated(
        "BioMmcifReadError",
        "Mmcif",
        &source.to_string(),
        &cause,
        detail.into(),
    )?;
    set(&error, "stage", (source.stage() as u32).into())?;
    Ok(error)
}
pub(crate) fn read_error(source: &ck::BioReadError) -> Result<JsValue, JsValue> {
    use ck::BioReadError as E;
    let kind = match source {
        E::Io { .. } => "Io",
        E::Utf8 { .. } => "Utf8",
        E::Pdb(_) => "Pdb",
        E::Cif(_) => "Cif",
        E::Mmcif(_) => "Mmcif",
        E::Mmjson(_) => "Mmjson",
        E::ChemComp(_) => "ChemComp",
        E::WrongFormat { .. } => "WrongFormat",
        E::UnknownFileFormat(_) => "UnknownFileFormat",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let detail = BioReadError {
        kind: kind.into(),
        message: source.to_string(),
        cause: cause.clone(),
    };
    let error = annotated(
        "BioReadError",
        kind,
        &source.to_string(),
        &cause,
        detail.into(),
    )?;
    match source {
        E::Io { path, .. } | E::Utf8 { path, .. } | E::UnknownFileFormat(path) => {
            set(&error, "path", path.to_string_lossy().into_owned().into())?
        }
        E::WrongFormat {
            source_name,
            format,
        } => {
            set(&error, "sourceName", source_name.clone().into())?;
            set(&error, "format", (*format as u32).into())?;
        }
        _ => (),
    }
    Ok(error)
}
