//! Exact BIO source identifiers with native constructors and immutable reads.
use crate::{fields::Field, input};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct AtomName(pub(crate) ck::AtomName);
crate::fields::value_field!(ck::AtomName, AtomName, AtomNameInput, "AtomName");
#[wasm_bindgen]
pub struct ResidueName(pub(crate) ck::ResidueName);
crate::fields::value_field!(
    ck::ResidueName,
    ResidueName,
    ResidueNameInput,
    "ResidueName"
);
#[wasm_bindgen]
pub struct PdbChainId(pub(crate) ck::PdbChainId);
crate::fields::value_field!(ck::PdbChainId, PdbChainId, PdbChainIdInput, "PdbChainId");
#[wasm_bindgen]
pub struct PdbAtomSerial(pub(crate) ck::PdbAtomSerial);
crate::fields::value_field!(
    ck::PdbAtomSerial,
    PdbAtomSerial,
    PdbAtomSerialInput,
    "PdbAtomSerial"
);
#[wasm_bindgen]
pub struct PdbSeqId(pub(crate) ck::PdbSeqId);
crate::fields::value_field!(ck::PdbSeqId, PdbSeqId, PdbSeqIdInput, "PdbSeqId");
#[wasm_bindgen]
pub struct AltLocLabel(pub(crate) ck::AltLocLabel);
crate::fields::value_field!(
    ck::AltLocLabel,
    AltLocLabel,
    AltLocLabelInput,
    "AltLocLabel"
);
#[wasm_bindgen]
pub struct AtomSourceIds(pub(crate) ck::AtomSourceIds);
crate::fields::value_field!(
    ck::AtomSourceIds,
    AtomSourceIds,
    AtomSourceIdsInput,
    "AtomSourceIds"
);
#[wasm_bindgen]
pub struct ChainSourceIds(pub(crate) ck::ChainSourceIds);
crate::fields::value_field!(
    ck::ChainSourceIds,
    ChainSourceIds,
    ChainSourceIdsInput,
    "ChainSourceIds"
);
#[wasm_bindgen]
pub struct ResidueSourceIds(pub(crate) ck::ResidueSourceIds);
crate::fields::value_field!(
    ck::ResidueSourceIds,
    ResidueSourceIds,
    ResidueSourceIdsInput,
    "ResidueSourceIds"
);
#[wasm_bindgen]
impl AtomName {
    #[wasm_bindgen(js_name=fromAscii)]
    pub fn from_ascii(bytes: Ts<input::BytesInput>) -> Result<Option<Self>, JsValue> {
        Ok(ck::AtomName::from_ascii(&input::read::<Vec<u8>>(bytes)?).map(Self))
    }
    #[wasm_bindgen(js_name=asStr)]
    pub fn as_str(&self) -> String {
        self.0.as_str().into()
    }
    #[wasm_bindgen(js_name=asBytes)]
    pub fn as_bytes(&self) -> js_sys::Uint8Array {
        js_sys::Uint8Array::from(self.0.as_bytes())
    }
}
#[wasm_bindgen]
impl ResidueName {
    #[wasm_bindgen(js_name=fromAscii)]
    pub fn from_ascii(bytes: Ts<input::BytesInput>) -> Result<Option<Self>, JsValue> {
        Ok(ck::ResidueName::from_ascii(&input::read::<Vec<u8>>(bytes)?).map(Self))
    }
    #[wasm_bindgen(js_name=asStr)]
    pub fn as_str(&self) -> String {
        self.0.as_str().into()
    }
    #[wasm_bindgen(js_name=asBytes)]
    pub fn as_bytes(&self) -> js_sys::Uint8Array {
        js_sys::Uint8Array::from(self.0.as_bytes())
    }
}
#[wasm_bindgen]
impl PdbChainId {
    #[wasm_bindgen(js_name=fromAscii)]
    pub fn from_ascii(bytes: Ts<input::BytesInput>) -> Result<Option<Self>, JsValue> {
        Ok(ck::PdbChainId::from_ascii(&input::read::<Vec<u8>>(bytes)?).map(Self))
    }
    #[wasm_bindgen(js_name=asStr)]
    pub fn as_str(&self) -> String {
        self.0.as_str().into()
    }
    #[wasm_bindgen(js_name=asBytes)]
    pub fn as_bytes(&self) -> js_sys::Uint8Array {
        js_sys::Uint8Array::from(self.0.as_bytes())
    }
}
#[wasm_bindgen]
impl PdbAtomSerial {
    #[wasm_bindgen(constructor)]
    pub fn new(value: Ts<input::I32Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::PdbAtomSerial::new(input::read::<i32>(value)?)))
    }
    #[wasm_bindgen(getter)]
    pub fn value(&self) -> i32 {
        self.0.value()
    }
}
#[wasm_bindgen]
impl AltLocLabel {
    #[wasm_bindgen(constructor)]
    pub fn new(value: Ts<input::U8Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::AltLocLabel::new(input::read::<u8>(value)?)))
    }
    #[wasm_bindgen(getter)]
    pub fn value(&self) -> u8 {
        self.0.value()
    }
}
#[wasm_bindgen]
impl PdbSeqId {
    #[wasm_bindgen(constructor)]
    pub fn new(
        seq_num: Ts<input::I32Input>,
        ins_code: Option<Ts<input::U8Input>>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::PdbSeqId::new(
            input::read::<i32>(seq_num)?,
            ins_code.map(input::read::<u8>).transpose()?,
        )))
    }
    #[wasm_bindgen(getter,js_name=seqNum)]
    pub fn seq_num(&self) -> i32 {
        self.0.seq_num()
    }
    #[wasm_bindgen(getter,js_name=insCode)]
    pub fn ins_code(&self) -> Option<u8> {
        self.0.ins_code()
    }
}
#[wasm_bindgen]
impl AtomSourceIds {
    #[wasm_bindgen(constructor)]
    pub fn new(serial: Option<PdbAtomSerialInput>) -> Result<Self, JsValue> {
        PdbAtomSerial::with_optional_ref(serial.as_ref(), |v| {
            Ok(Self(ck::AtomSourceIds::new(v.copied())))
        })
    }
    #[wasm_bindgen(getter)]
    pub fn serial(&self) -> Option<PdbAtomSerial> {
        self.0.serial().map(PdbAtomSerial)
    }
}
#[wasm_bindgen]
impl ChainSourceIds {
    #[wasm_bindgen(constructor)]
    pub fn new(
        auth_chain_id: Option<PdbChainIdInput>,
        label_asym_id: Option<Ts<input::TextInput>>,
    ) -> Result<Self, JsValue> {
        let label = label_asym_id.map(input::read::<String>).transpose()?;
        PdbChainId::with_optional_ref(auth_chain_id.as_ref(), |v| {
            Ok(Self(ck::ChainSourceIds::new(v.copied(), label)))
        })
    }
    #[wasm_bindgen(getter,js_name=authChainId)]
    pub fn auth_chain_id(&self) -> Option<PdbChainId> {
        self.0.auth_chain_id().map(PdbChainId)
    }
    #[wasm_bindgen(getter,js_name=labelAsymId)]
    pub fn label_asym_id(&self) -> Option<String> {
        self.0.label_asym_id().map(str::to_owned)
    }
}
#[derive(serde::Serialize, serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct SegmentId(#[serde(deserialize_with = "input::exact_array")] [u8; 4]);
#[wasm_bindgen]
impl ResidueSourceIds {
    #[wasm_bindgen(js_name=new)]
    pub fn new(
        seq_id: Option<PdbSeqIdInput>,
        label_seq_id: Option<Ts<input::I32Input>>,
        segment_id: Option<Ts<SegmentId>>,
        subchain_id: Option<Ts<input::TextInput>>,
        label_entity_id: Option<Ts<input::TextInput>>,
    ) -> Result<Option<Self>, JsValue> {
        let label_seq_id = label_seq_id.map(input::read::<i32>).transpose()?;
        let segment_id = segment_id.map(input::decode).transpose()?.map(|v| v.0);
        let subchain_id = subchain_id.map(input::read::<String>).transpose()?;
        let label_entity_id = label_entity_id.map(input::read::<String>).transpose()?;
        PdbSeqId::with_optional_ref(seq_id.as_ref(), |v| {
            Ok(ck::ResidueSourceIds::new(
                v.copied(),
                label_seq_id,
                segment_id,
                subchain_id,
                label_entity_id,
            )
            .map(Self))
        })
    }
    #[wasm_bindgen(getter,js_name=seqId)]
    pub fn seq_id(&self) -> Option<PdbSeqId> {
        self.0.seq_id().map(PdbSeqId)
    }
    #[wasm_bindgen(getter,js_name=labelSeqId)]
    pub fn label_seq_id(&self) -> Option<i32> {
        self.0.label_seq_id()
    }
    #[wasm_bindgen(getter,js_name=segmentId)]
    pub fn segment_id(&self) -> Result<Option<Ts<SegmentId>>, JsValue> {
        self.0
            .segment_id()
            .map(|v| {
                SegmentId(*v)
                    .into_ts()
                    .map_err(JsError::from)
                    .map_err(Into::into)
            })
            .transpose()
    }
    #[wasm_bindgen(getter,js_name=subchainId)]
    pub fn subchain_id(&self) -> Option<String> {
        self.0.subchain_id().map(str::to_owned)
    }
    #[wasm_bindgen(getter,js_name=labelEntityId)]
    pub fn label_entity_id(&self) -> Option<String> {
        self.0.label_entity_id().map(str::to_owned)
    }
}
