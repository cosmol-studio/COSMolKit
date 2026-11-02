//! Binary value transfer, with no pointer exchange or binding-owned codec.
use crate::errors::enum_error;
use crate::{Molecule, input};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "StereoGroupError")]
    #[derive(Clone)]
    pub type StereoGroupCause;
    #[wasm_bindgen(typescript_type = "PickleError")]
    #[derive(Clone)]
    pub type CodecCause;
}

enum_error!(PickleError, PickleErrorKind, ck::PickleError, {
    StereoGroup => ck::PickleError::StereoGroup(_),
    UnexpectedEof => ck::PickleError::UnexpectedEof,
    UnsupportedVersion => ck::PickleError::UnsupportedVersion(_),
    UnsupportedArchiveVersion => ck::PickleError::UnsupportedArchiveVersion { .. },
    UnsupportedSectionVersion => ck::PickleError::UnsupportedSectionVersion { .. },
    InvalidArchive => ck::PickleError::InvalidArchive(_),
    MissingRequiredSection => ck::PickleError::MissingRequiredSection(_),
    DuplicateSection => ck::PickleError::DuplicateSection(_),
    UnknownRequiredSection => ck::PickleError::UnknownRequiredSection(_),
    DataLengthMismatch => ck::PickleError::DataLengthMismatch { .. },
    InvalidEnumValue => ck::PickleError::InvalidEnumValue { .. },
    InvalidMolecule => ck::PickleError::InvalidMolecule(_),
    TooManyAtoms => ck::PickleError::TooManyAtoms(_),
    TooManyBonds => ck::PickleError::TooManyBonds(_),
    StringTooLong => ck::PickleError::StringTooLong(_),
}, {
    #[wasm_bindgen(getter)]
    pub fn version(&self) -> Option<u16> { match self.inner.as_ref() {
        ck::PickleError::UnsupportedVersion(value) => Some(u16::from(*value)),
        ck::PickleError::UnsupportedSectionVersion { version, .. } => Some(*version),
        _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn section(&self) -> Option<u16> { match self.inner.as_ref() {
        ck::PickleError::UnsupportedSectionVersion { section, .. }
        | ck::PickleError::MissingRequiredSection(section)
        | ck::PickleError::DuplicateSection(section)
        | ck::PickleError::UnknownRequiredSection(section) => Some(*section), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn major(&self) -> Option<u16> { match self.inner.as_ref() {
        ck::PickleError::UnsupportedArchiveVersion { major, .. } => Some(*major), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn minor(&self) -> Option<u16> { match self.inner.as_ref() {
        ck::PickleError::UnsupportedArchiveVersion { minor, .. } => Some(*minor), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn expected(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::PickleError::DataLengthMismatch { expected, .. } => Some(*expected), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn actual(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::PickleError::DataLengthMismatch { actual, .. } => Some(*actual), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn value(&self) -> Option<u8> { match self.inner.as_ref() {
        ck::PickleError::InvalidEnumValue { value, .. } => Some(*value), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = typeName)]
    pub fn type_name(&self) -> Option<String> { match self.inner.as_ref() {
        ck::PickleError::InvalidEnumValue { type_name, .. } => Some((*type_name).into()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> Option<String> { match self.inner.as_ref() {
        ck::PickleError::InvalidArchive(message) | ck::PickleError::InvalidMolecule(message) => Some(message.clone()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn length(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::PickleError::TooManyAtoms(length) | ck::PickleError::TooManyBonds(length)
        | ck::PickleError::StringTooLong(length) => Some(*length), _ => None,
    } }
}, source = StereoGroupCause, cause = |inner: &crate::storage::ReadonlyValue<ck::PickleError>| match inner.as_ref() {
    ck::PickleError::StereoGroup(error) => Some(crate::model_errors::StereoGroupError::from_rust(*error).unchecked_into()),
    _ => None,
});

enum_error!(MoleculeTransferError, MoleculeTransferErrorKind, ck::MoleculeTransferError, {
    VersionMismatch => ck::MoleculeTransferError::VersionMismatch { .. },
    Codec => ck::MoleculeTransferError::Codec(_),
}, {
    #[wasm_bindgen(getter)]
    pub fn expected(&self) -> Option<String> { match self.inner.as_ref() {
        ck::MoleculeTransferError::VersionMismatch { expected, .. } => Some(expected.clone()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn actual(&self) -> Option<String> { match self.inner.as_ref() {
        ck::MoleculeTransferError::VersionMismatch { actual, .. } => Some(actual.clone()), _ => None,
    } }
}, source = CodecCause, cause = |inner: &crate::storage::ReadonlyValue<ck::MoleculeTransferError>| match inner.as_ref() {
    ck::MoleculeTransferError::Codec(error) => Some(PickleError::from_rust(error.clone()).unchecked_into()),
    ck::MoleculeTransferError::VersionMismatch { .. } => None,
});

#[wasm_bindgen]
impl Molecule {
    /// @internal
    /// Transfer this value to another CK module of exactly the same version.
    /// Includes topology, metadata and all coordinates, but not derived caches.
    #[wasm_bindgen(js_name = toTransferBinary)]
    pub fn to_transfer_binary(&self) -> Result<js_sys::Uint8Array, JsValue> {
        self.0.to_transfer_binary()
            .map(|bytes| js_sys::Uint8Array::from(bytes.as_slice()))
            .map_err(MoleculeTransferError::from_rust)
    }

    /// @internal
    /// Restore a value owned by this module; reject every different CK version.
    #[wasm_bindgen(js_name = fromTransferBinary)]
    pub fn from_transfer_binary(data: Ts<input::BytesInput>) -> Result<Molecule, JsValue> {
        let data = input::read::<Vec<u8>>(data)?;
        ck::Molecule::from_transfer_binary(&data)
            .map(Molecule)
            .map_err(MoleculeTransferError::from_rust)
    }
}
