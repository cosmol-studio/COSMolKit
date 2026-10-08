//! Immutable BIO hierarchy projections; snapshots retain their public owner.
use crate::{
    bio_readers::BioStructure,
    host_values::{i32_value, integer, sequence, type_error, u32_value},
};
use cosmolkit_wasm::rust as ck;
use std::rc::Rc;
use wasm_bindgen::prelude::*;
fn nullable<T: Into<JsValue>>(value: Option<T>) -> JsValue {
    value.map_or(JsValue::NULL, Into::into)
}
fn bytes(value: &JsValue) -> Result<Vec<u8>, JsValue> {
    sequence(value, "bytes")?
        .iter()
        .map(|v| integer(&v, "byte", 0., 255.).map(|n| n as u8))
        .collect()
}

#[wasm_bindgen]
pub enum BioCalcFlag {
    NotSet = 0,
    NoHydrogen = 1,
    Determined = 2,
    Calculated = 3,
    Dummy = 4,
}
#[wasm_bindgen]
pub enum EntityKind {
    Unknown = 0,
    Polymer = 1,
    NonPolymer = 2,
    Branched = 3,
    Water = 4,
}
#[wasm_bindgen]
pub enum PolymerKind {
    Unknown = 0,
    PeptideL = 1,
    PeptideD = 2,
    Dna = 3,
    Rna = 4,
    DnaRnaHybrid = 5,
    SaccharideD = 6,
    SaccharideL = 7,
    Pna = 8,
    CyclicPseudoPeptide = 9,
    Other = 10,
}
#[wasm_bindgen]
pub enum ResidueKind {
    AminoAcid = 0,
    Dna = 1,
    Rna = 2,
    Saccharide = 3,
    Water = 4,
    Buffer = 5,
    Ligand = 6,
    Unknown = 7,
}
#[wasm_bindgen]
pub enum ChainKind {
    Protein = 0,
    Dna = 1,
    Rna = 2,
    ProteinDnaComplex = 3,
    ProteinRnaComplex = 4,
    LigandOnly = 5,
    WaterOnly = 6,
    Mixed = 7,
    Unknown = 8,
}
#[wasm_bindgen]
pub struct BioAtomId {
    inner: ck::BioAtomId,
}
#[wasm_bindgen]
impl BioAtomId {
    #[wasm_bindgen(constructor)]
    pub fn new(value: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::BioAtomId::new(u32_value(&value, "value")?),
        })
    }
    pub fn value(&self) -> u32 {
        self.inner.value()
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
}
#[wasm_bindgen]
pub struct BioResidueId {
    inner: ck::BioResidueId,
}
#[wasm_bindgen]
impl BioResidueId {
    #[wasm_bindgen(constructor)]
    pub fn new(value: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::BioResidueId::new(u32_value(&value, "value")?),
        })
    }
    pub fn value(&self) -> u32 {
        self.inner.value()
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
}
#[wasm_bindgen]
pub struct BioChainId {
    inner: ck::BioChainId,
}
#[wasm_bindgen]
impl BioChainId {
    #[wasm_bindgen(constructor)]
    pub fn new(value: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::BioChainId::new(u32_value(&value, "value")?),
        })
    }
    pub fn value(&self) -> u32 {
        self.inner.value()
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
}
#[wasm_bindgen]
pub struct BioEntityId {
    inner: ck::BioEntityId,
}
#[wasm_bindgen]
impl BioEntityId {
    #[wasm_bindgen(constructor)]
    pub fn new(value: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::BioEntityId::new(u32_value(&value, "value")?),
        })
    }
    pub fn value(&self) -> u32 {
        self.inner.value()
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
}
#[wasm_bindgen]
pub struct BioModelId {
    inner: ck::BioModelId,
}
#[wasm_bindgen]
impl BioModelId {
    #[wasm_bindgen(constructor)]
    pub fn new(value: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::BioModelId::new(u32_value(&value, "value")?),
        })
    }
    pub fn value(&self) -> u32 {
        self.inner.value()
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
}
#[wasm_bindgen]
pub struct BioAssemblyId {
    inner: ck::BioAssemblyId,
}
#[wasm_bindgen]
impl BioAssemblyId {
    #[wasm_bindgen(constructor)]
    pub fn new(value: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::BioAssemblyId::new(u32_value(&value, "value")?),
        })
    }
    pub fn value(&self) -> u32 {
        self.inner.value()
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
}
#[wasm_bindgen]
pub struct BioAltLocGroupId {
    inner: ck::BioAltLocGroupId,
}
#[wasm_bindgen]
impl BioAltLocGroupId {
    #[wasm_bindgen(constructor)]
    pub fn new(value: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::BioAltLocGroupId::new(u32_value(&value, "value")?),
        })
    }
    pub fn value(&self) -> u32 {
        self.inner.value()
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
}
#[wasm_bindgen]
pub struct AtomName {
    pub(crate) inner: ck::AtomName,
}
#[wasm_bindgen]
impl AtomName {
    #[wasm_bindgen(js_name=fromAscii,unchecked_return_type="AtomName | null")]
    pub fn from_ascii(
        #[wasm_bindgen(unchecked_param_type = "Uint8Array | number[]")] value: JsValue,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: ck::AtomName::from_ascii(bytes.as_bytes()).map(|inner| Self { inner })
        Ok(nullable(
            ck::AtomName::from_ascii(&bytes(&value)?).map(|inner| Self { inner }),
        ))
    }
    #[wasm_bindgen(js_name=asBytes)]
    pub fn as_bytes(&self) -> Vec<u8> {
        self.inner.as_bytes().to_vec()
    }
    #[wasm_bindgen(js_name=asStr)]
    pub fn as_str(&self) -> String {
        self.inner.as_str().into()
    }
}
#[wasm_bindgen]
pub struct ResidueName {
    pub(crate) inner: ck::ResidueName,
}
#[wasm_bindgen]
impl ResidueName {
    #[wasm_bindgen(js_name=fromAscii,unchecked_return_type="ResidueName | null")]
    pub fn from_ascii(
        #[wasm_bindgen(unchecked_param_type = "Uint8Array | number[]")] value: JsValue,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: ck::ResidueName::from_ascii(bytes.as_bytes()).map(|inner| Self { inner })
        Ok(nullable(
            ck::ResidueName::from_ascii(&bytes(&value)?).map(|inner| Self { inner }),
        ))
    }
    #[wasm_bindgen(js_name=asBytes)]
    pub fn as_bytes(&self) -> Vec<u8> {
        self.inner.as_bytes().to_vec()
    }
    #[wasm_bindgen(js_name=asStr)]
    pub fn as_str(&self) -> String {
        self.inner.as_str().into()
    }
}
#[wasm_bindgen]
pub struct PdbChainId {
    pub(crate) inner: ck::PdbChainId,
}
#[wasm_bindgen]
impl PdbChainId {
    #[wasm_bindgen(js_name=fromAscii,unchecked_return_type="PdbChainId | null")]
    pub fn from_ascii(
        #[wasm_bindgen(unchecked_param_type = "Uint8Array | number[]")] value: JsValue,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: ck::PdbChainId::from_ascii(bytes.as_bytes()).map(|inner| Self { inner })
        Ok(nullable(
            ck::PdbChainId::from_ascii(&bytes(&value)?).map(|inner| Self { inner }),
        ))
    }
    #[wasm_bindgen(js_name=asBytes)]
    pub fn as_bytes(&self) -> Vec<u8> {
        self.inner.as_bytes().to_vec()
    }
    #[wasm_bindgen(js_name=asStr)]
    pub fn as_str(&self) -> String {
        self.inner.as_str().into()
    }
}
#[wasm_bindgen]
pub struct PdbAtomSerial {
    inner: ck::PdbAtomSerial,
}
#[wasm_bindgen]
impl PdbAtomSerial {
    #[wasm_bindgen(constructor)]
    pub fn new(value: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::PdbAtomSerial::new(i32_value(&value, "serial")?),
        })
    }
    pub fn value(&self) -> i32 {
        self.inner.value()
    }
}
#[wasm_bindgen]
pub struct PdbSeqId {
    pub(crate) inner: ck::PdbSeqId,
}
#[wasm_bindgen]
impl PdbSeqId {
    #[wasm_bindgen(js_name=seqNum)]
    pub fn seq_num(&self) -> i32 {
        self.inner.seq_num()
    }
    #[wasm_bindgen(js_name=insCode,unchecked_return_type="string | null")]
    pub fn ins_code(&self) -> JsValue {
        nullable(self.inner.ins_code().map(|x| char::from(x).to_string()))
    }
}
#[wasm_bindgen]
pub struct AltLocLabel {
    pub(crate) inner: ck::AltLocLabel,
}
#[wasm_bindgen]
impl AltLocLabel {
    #[wasm_bindgen(constructor)]
    pub fn new(value: JsValue) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::AltLocLabel::new(integer(&value, "altloc", 0., 255.)? as u8),
        })
    }
    pub fn value(&self) -> u8 {
        self.inner.value()
    }
}
#[wasm_bindgen]
pub struct AltLocRequest {
    pub(crate) inner: ck::AltLocRequest,
}
#[wasm_bindgen]
impl AltLocRequest {
    #[wasm_bindgen(getter,js_name=Any)]
    pub fn any() -> Self {
        Self {
            inner: ck::AltLocRequest::Any,
        }
    }
    #[wasm_bindgen(js_name=Exact)]
    pub fn exact(
        #[wasm_bindgen(unchecked_param_type = "AltLocLabel | null")] altloc: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::AltLocRequest::Exact(optional_label(&altloc)?),
        })
    }
}
#[wasm_bindgen]
pub struct BioRowSpan {
    inner: ck::BioRowSpan<ck::BioAtomId>,
}
#[wasm_bindgen]
impl BioRowSpan {
    #[wasm_bindgen(constructor)]
    pub fn new(start: JsValue, len: JsValue) -> Result<Self, JsValue> {
        ck::BioRowSpan::new(u32_value(&start, "start")?, u32_value(&len, "len")?)
            .map(|inner| Self { inner })
            .map_err(|e| crate::bio_hierarchy_errors::structure_error(&e).unwrap_or_else(|e| e))
    }
    pub fn start(&self) -> u32 {
        self.inner.start()
    }
    pub fn len(&self) -> u32 {
        self.inner.len()
    }
    pub fn end(&self) -> u32 {
        self.inner.end()
    }
    #[wasm_bindgen(js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

#[wasm_bindgen]
pub struct AtomSourceIds {
    inner: ck::AtomSourceIds,
}
#[wasm_bindgen]
impl AtomSourceIds {
    #[wasm_bindgen(js_name=serial,unchecked_return_type="number | null")]
    pub fn serial(&self) -> JsValue {
        nullable(self.inner.serial().map(|v| v.value()))
    }
}
#[wasm_bindgen]
pub struct ChainSourceIds {
    pub(crate) inner: ck::ChainSourceIds,
}
#[wasm_bindgen]
impl ChainSourceIds {
    #[wasm_bindgen(js_name=authChainId,unchecked_return_type="string | null")]
    pub fn auth_chain_id(&self) -> JsValue {
        nullable(self.inner.auth_chain_id().map(|v| v.as_str().to_owned()))
    }
    #[wasm_bindgen(js_name=labelAsymId,unchecked_return_type="string | null")]
    pub fn label_asym_id(&self) -> JsValue {
        nullable(self.inner.label_asym_id().map(str::to_owned))
    }
}
#[wasm_bindgen]
pub struct EntitySourceIds {
    inner: ck::EntitySourceIds,
}
#[wasm_bindgen]
impl EntitySourceIds {
    #[wasm_bindgen(js_name=sourceEntityId,unchecked_return_type="string")]
    pub fn source_entity_id(&self) -> JsValue {
        self.inner.source_entity_id().into()
    }
}
#[wasm_bindgen]
pub struct ResidueSourceIds {
    inner: ck::ResidueSourceIds,
}
#[wasm_bindgen]
impl ResidueSourceIds {
    #[wasm_bindgen(js_name=seqId,unchecked_return_type="PdbSeqId | null")]
    pub fn seq_id(&self) -> JsValue {
        nullable(self.inner.seq_id().map(|inner| PdbSeqId { inner }))
    }
    #[wasm_bindgen(js_name=labelSeqId,unchecked_return_type="number | null")]
    pub fn label_seq_id(&self) -> JsValue {
        nullable(self.inner.label_seq_id())
    }
    #[wasm_bindgen(js_name=segmentId,unchecked_return_type="Uint8Array | null")]
    pub fn segment_id(&self) -> JsValue {
        nullable(
            self.inner
                .segment_id()
                .map(|v| js_sys::Uint8Array::from(v.as_slice())),
        )
    }
    #[wasm_bindgen(js_name=subchainId,unchecked_return_type="string | null")]
    pub fn subchain_id(&self) -> JsValue {
        nullable(self.inner.subchain_id().map(str::to_owned))
    }
    #[wasm_bindgen(js_name=labelEntityId,unchecked_return_type="string | null")]
    pub fn label_entity_id(&self) -> JsValue {
        nullable(self.inner.label_entity_id().map(str::to_owned))
    }
}
#[wasm_bindgen]
pub struct BioModelRow {
    pub(crate) inner: Rc<ck::BioStructure>,
    pub(crate) index: usize,
}
impl BioModelRow {
    fn row(&self) -> &ck::BioModelRow {
        &self.inner.models()[self.index]
    }
}
#[wasm_bindgen]
impl BioModelRow {
    pub fn id(&self) -> usize {
        self.index
    }
    #[wasm_bindgen(js_name=sourceModelNumber,unchecked_return_type="number | null")]
    pub fn source_model_number(&self) -> Result<JsValue, JsValue> {
        Ok(nullable(self.row().source_model_number()))
    }
    #[wasm_bindgen(unchecked_return_type = "BioChainRow[]")]
    pub fn chains(&self) -> js_sys::Array {
        // COSMolKit❗✔️: inner: self.inner.clone(), index,
        // Each row retains an O(1) shared owner reference; output allocation follows Python.
        let span = self.row().chain_span();
        (span.start() as usize..span.end() as usize)
            .map(|index| -> JsValue {
                BioChainRow {
                    inner: Rc::clone(&self.inner),
                    index,
                }
                .into()
            })
            .collect()
    }
    pub fn len(&self) -> u32 {
        self.row().chain_span().len()
    }
}
#[wasm_bindgen]
pub struct BioChainRow {
    pub(crate) inner: Rc<ck::BioStructure>,
    pub(crate) index: usize,
}
impl BioChainRow {
    fn row(&self) -> &ck::BioChainRow {
        &self.inner.chains()[self.index]
    }
}
#[wasm_bindgen]
impl BioChainRow {
    pub fn id(&self) -> usize {
        self.index
    }
    #[wasm_bindgen(js_name=modelId,unchecked_return_type="number")]
    pub fn model_id(&self) -> Result<JsValue, JsValue> {
        Ok(self.row().model_id().value().into())
    }
    #[wasm_bindgen(js_name=entityId,unchecked_return_type="number | null")]
    pub fn entity_id(&self) -> Result<JsValue, JsValue> {
        Ok(nullable(self.row().entity_id().map(|x| x.value())))
    }
    #[wasm_bindgen(js_name=kind,unchecked_return_type="ChainKind")]
    pub fn kind(&self) -> Result<JsValue, JsValue> {
        Ok((self.row().kind() as u32).into())
    }
    #[wasm_bindgen(js_name=source,unchecked_return_type="ChainSourceIds")]
    pub fn source(&self) -> Result<JsValue, JsValue> {
        Ok(ChainSourceIds {
            inner: self.row().source().clone(),
        }
        .into())
    }
    #[wasm_bindgen(unchecked_return_type = "BioResidueRow[]")]
    pub fn residues(&self) -> js_sys::Array {
        // COSMolKit❗✔️: inner: self.inner.clone(), index,
        // Each row retains an O(1) shared owner reference; output allocation follows Python.
        let span = self.row().residue_span();
        (span.start() as usize..span.end() as usize)
            .map(|index| -> JsValue {
                BioResidueRow {
                    inner: Rc::clone(&self.inner),
                    index,
                }
                .into()
            })
            .collect()
    }
    pub fn len(&self) -> u32 {
        self.row().residue_span().len()
    }
    #[wasm_bindgen(unchecked_return_type = "BioAtomRow[]")]
    pub fn atoms(&self) -> js_sys::Array {
        let s = self.row().residue_span();
        self.inner.residues()[s.start() as usize..s.end() as usize]
            .iter()
            .flat_map(|r| r.atom_span().start() as usize..r.atom_span().end() as usize)
            .map(|index| -> JsValue {
                BioAtomRow {
                    inner: Rc::clone(&self.inner),
                    index,
                }
                .into()
            })
            .collect()
    }
}
#[wasm_bindgen]
pub struct BioResidueRow {
    pub(crate) inner: Rc<ck::BioStructure>,
    pub(crate) index: usize,
}
impl BioResidueRow {
    fn row(&self) -> &ck::BioResidueRow {
        &self.inner.residues()[self.index]
    }
}
#[wasm_bindgen]
impl BioResidueRow {
    pub fn id(&self) -> usize {
        self.index
    }
    #[wasm_bindgen(js_name=chainId,unchecked_return_type="number")]
    pub fn chain_id(&self) -> Result<JsValue, JsValue> {
        Ok(self.row().chain_id().value().into())
    }
    #[wasm_bindgen(js_name=name,unchecked_return_type="string")]
    pub fn name(&self) -> Result<JsValue, JsValue> {
        Ok(self.row().name().as_str().into())
    }
    #[wasm_bindgen(js_name=kind,unchecked_return_type="ResidueKind")]
    pub fn kind(&self) -> Result<JsValue, JsValue> {
        Ok((self.row().kind() as u32).into())
    }
    #[wasm_bindgen(js_name=entityKind,unchecked_return_type="EntityKind")]
    pub fn entity_kind(&self) -> Result<JsValue, JsValue> {
        Ok((self.row().entity_kind() as u32).into())
    }
    #[wasm_bindgen(js_name=info,unchecked_return_type="ResidueInfo")]
    pub fn info(&self) -> Result<JsValue, JsValue> {
        Ok(crate::bio_residue::ResidueInfo {
            inner: ck::find_residue_info(self.row().name().as_str()),
        }
        .into())
    }
    #[wasm_bindgen(js_name=code,unchecked_return_type="ResidueCode")]
    pub fn code(&self) -> Result<JsValue, JsValue> {
        Ok((ck::residue_code(self.row().name().as_str()).as_u16() as u32).into())
    }
    #[wasm_bindgen(js_name=source,unchecked_return_type="ResidueSourceIds")]
    pub fn source(&self) -> Result<JsValue, JsValue> {
        Ok(ResidueSourceIds {
            inner: self.row().source().clone(),
        }
        .into())
    }
    #[wasm_bindgen(unchecked_return_type = "BioAtomRow[]")]
    pub fn atoms(&self) -> js_sys::Array {
        // COSMolKit❗✔️: inner: self.inner.clone(), index,
        // Each row retains an O(1) shared owner reference; output allocation follows Python.
        let span = self.row().atom_span();
        (span.start() as usize..span.end() as usize)
            .map(|index| -> JsValue {
                BioAtomRow {
                    inner: Rc::clone(&self.inner),
                    index,
                }
                .into()
            })
            .collect()
    }
    pub fn len(&self) -> u32 {
        self.row().atom_span().len()
    }
}
#[wasm_bindgen]
pub struct BioAtomRow {
    pub(crate) inner: Rc<ck::BioStructure>,
    pub(crate) index: usize,
}
impl BioAtomRow {
    fn row(&self) -> &ck::BioAtomRow {
        &self.inner.atoms()[self.index]
    }
}
#[wasm_bindgen]
impl BioAtomRow {
    pub fn id(&self) -> usize {
        self.index
    }
    #[wasm_bindgen(js_name=residueId,unchecked_return_type="number")]
    pub fn residue_id(&self) -> Result<JsValue, JsValue> {
        Ok(self.row().residue_id().value().into())
    }
    #[wasm_bindgen(js_name=name,unchecked_return_type="string")]
    pub fn name(&self) -> Result<JsValue, JsValue> {
        Ok(self.row().name().as_str().trim().into())
    }
    #[wasm_bindgen(js_name=element,unchecked_return_type="Element")]
    pub fn element(&self) -> Result<JsValue, JsValue> {
        Ok(crate::element_boundary::Element::from_atomic_number(
            self.row().element().atomic_number() as f64,
        )?)
    }
    #[wasm_bindgen(js_name=elementSymbol,unchecked_return_type="string")]
    pub fn element_symbol(&self) -> Result<JsValue, JsValue> {
        Ok(self.row().element().symbol().into())
    }
    #[wasm_bindgen(js_name=position,unchecked_return_type="number[] | null")]
    pub fn position(&self) -> Result<JsValue, JsValue> {
        Ok(nullable(
            self.inner
                .atom_position(ck::BioAtomId::new(self.index as u32))
                .map(|p| p.into_iter().map(JsValue::from).collect::<js_sys::Array>()),
        ))
    }
    #[wasm_bindgen(js_name=altloc,unchecked_return_type="string | null")]
    pub fn altloc(&self) -> Result<JsValue, JsValue> {
        Ok(nullable(
            self.row()
                .altloc()
                .map(|v| char::from(v.value()).to_string()),
        ))
    }
    #[wasm_bindgen(js_name=source,unchecked_return_type="AtomSourceIds")]
    pub fn source(&self) -> Result<JsValue, JsValue> {
        Ok(AtomSourceIds {
            inner: *self.row().source(),
        }
        .into())
    }
    #[wasm_bindgen(js_name=occupancy)]
    pub fn occupancy(&self) -> f64 {
        self.row().occupancy()
    }
    #[wasm_bindgen(js_name=bIso)]
    pub fn b_iso(&self) -> f64 {
        self.row().b_iso()
    }
    #[wasm_bindgen(js_name=formalCharge)]
    pub fn formal_charge(&self) -> i8 {
        self.row().formal_charge()
    }
}
#[wasm_bindgen]
pub struct BioEntityRow {
    pub(crate) inner: Rc<ck::BioStructure>,
    pub(crate) index: usize,
}
impl BioEntityRow {
    fn row(&self) -> &ck::BioEntityRow {
        &self.inner.entities()[self.index]
    }
}
#[wasm_bindgen]
impl BioEntityRow {
    pub fn id(&self) -> usize {
        self.index
    }
    #[wasm_bindgen(js_name=source,unchecked_return_type="EntitySourceIds")]
    pub fn source(&self) -> Result<JsValue, JsValue> {
        Ok(EntitySourceIds {
            inner: self.row().source().clone(),
        }
        .into())
    }
    #[wasm_bindgen(js_name=kind,unchecked_return_type="EntityKind")]
    pub fn kind(&self) -> Result<JsValue, JsValue> {
        Ok((self.row().kind() as u32).into())
    }
    #[wasm_bindgen(js_name=polymerKind,unchecked_return_type="PolymerKind")]
    pub fn polymer_kind(&self) -> Result<JsValue, JsValue> {
        Ok((self.row().polymer_kind() as u32).into())
    }
    #[wasm_bindgen(js_name=fullSequence,unchecked_return_type="string[]")]
    pub fn full_sequence(&self) -> Result<JsValue, JsValue> {
        Ok(self
            .row()
            .full_sequence()
            .iter()
            .map(|v| JsValue::from(v.as_str()))
            .collect::<js_sys::Array>()
            .into())
    }
    #[wasm_bindgen(js_name=subchains,unchecked_return_type="string[]")]
    pub fn subchains(&self) -> Result<JsValue, JsValue> {
        Ok(self
            .row()
            .subchains()
            .iter()
            .map(|v| JsValue::from(v.as_str()))
            .collect::<js_sys::Array>()
            .into())
    }
    pub fn len(&self) -> usize {
        self.row().full_sequence().len()
    }
}
#[wasm_bindgen]
pub struct BioCoordinateBlock {
    inner: Rc<ck::BioStructure>,
}
#[wasm_bindgen]
impl BioCoordinateBlock {
    #[wasm_bindgen(getter, unchecked_return_type = "number[][]")]
    pub fn positions(&self) -> js_sys::Array {
        self.inner
            .coordinates()
            .positions()
            .iter()
            .map(|p| -> JsValue {
                p.iter()
                    .copied()
                    .map(JsValue::from)
                    .collect::<js_sys::Array>()
                    .into()
            })
            .collect()
    }
    pub fn len(&self) -> usize {
        self.inner.coordinates().len()
    }
    #[wasm_bindgen(js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.inner.coordinates().is_empty()
    }
}

impl BioCoordinateBlock {
    pub(crate) fn from_owner(inner: Rc<ck::BioStructure>) -> Self {
        Self { inner }
    }
}

#[wasm_bindgen(
    inline_js = "export function visitBioAltloc(v,f){try{f(v);}catch(cause){throw new TypeError('invalid AltLocLabel',{cause});}} export function visitBioElement(v,f){try{f(v);}catch(cause){throw new TypeError('invalid Element',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitBioAltloc)]
    fn visit_label(v: &JsValue, f: &mut dyn FnMut(&AltLocLabel)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitBioElement)]
    fn visit_element(
        v: &JsValue,
        f: &mut dyn FnMut(&crate::element_boundary::Element),
    ) -> Result<(), JsValue>;
}
pub(crate) fn optional_label(v: &JsValue) -> Result<Option<ck::AltLocLabel>, JsValue> {
    let mut out = None;
    if !v.is_null() && !v.is_undefined() {
        visit_label(v, &mut |v: &AltLocLabel| out = Some(v.inner))?;
    }
    Ok(out)
}
pub(crate) fn optional_element(v: &JsValue) -> Result<Option<ck::Element>, JsValue> {
    let mut out = None;
    if !v.is_null() && !v.is_undefined() {
        visit_element(v, &mut |v: &crate::element_boundary::Element| {
            out = ck::Element::from_atomic_number(v.atomic_number())
        })?;
    }
    Ok(out)
}
