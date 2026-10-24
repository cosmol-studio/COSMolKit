//! Detached native BIO rows; hierarchy relationships stay in the facade.
use crate::{
    bio_hierarchy::{ChainKind, ChainKindInput, ResidueKind},
    bio_read_errors::BioStructureError,
    bio_source_values::*,
    fields::Field,
    input,
    residue_values::{ResidueInfoKind, ResidueInfoKindInput},
};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
crate::enums::enum_value!(BioCalcFlag,"BioCalcFlag"{NotSet="not_set",NoHydrogen="no_hydrogen",Determined="determined",Calculated="calculated",Dummy="dummy"});
crate::enums::enum_value!(EntityKind,"EntityKind"{Unknown="unknown",Polymer="polymer",NonPolymer="non_polymer",Branched="branched",Water="water"});
crate::enums::enum_value!(PolymerKind,"PolymerKind"{Unknown="unknown",PeptideL="peptide_l",PeptideD="peptide_d",Dna="dna",Rna="rna",DnaRnaHybrid="dna_rna_hybrid",SaccharideD="saccharide_d",SaccharideL="saccharide_l",Pna="pna",CyclicPseudoPeptide="cyclic_pseudo_peptide",Other="other"});
#[derive(serde::Serialize, serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct Anisou(#[serde(deserialize_with = "input::exact_array")] pub [f64; 6]);
#[wasm_bindgen]
pub struct BioRowSpan(ck::BioRowSpan<ck::BioAtomId>);
crate::fields::value_field!(
    ck::BioRowSpan<ck::BioAtomId>,
    BioRowSpan,
    BioRowSpanInput,
    "BioRowSpan"
);
impl BioRowSpan {
    pub(crate) fn from_native<I>(span: ck::BioRowSpan<I>) -> Result<Self, JsValue> {
        ck::BioRowSpan::new(span.start(), span.len())
            .map(Self)
            .map_err(BioStructureError::from_rust)
    }
}
#[wasm_bindgen]
impl BioRowSpan {
    #[wasm_bindgen(constructor)]
    pub fn new(start: Ts<input::U32Input>, len: Ts<input::U32Input>) -> Result<Self, JsValue> {
        ck::BioRowSpan::new(input::read::<u32>(start)?, input::read::<u32>(len)?)
            .map(Self)
            .map_err(BioStructureError::from_rust)
    }
    #[wasm_bindgen(js_name=fromUsize)]
    pub fn from_usize(
        start: Ts<input::IndexInput>,
        len: Ts<input::IndexInput>,
    ) -> Result<Self, JsValue> {
        ck::BioRowSpan::from_usize(input::read::<usize>(start)?, input::read::<usize>(len)?)
            .map(Self)
            .map_err(BioStructureError::from_rust)
    }
    #[wasm_bindgen(getter)]
    pub fn start(&self) -> u32 {
        self.0.start()
    }
    #[wasm_bindgen(getter)]
    pub fn len(&self) -> u32 {
        self.0.len()
    }
    #[wasm_bindgen(getter)]
    pub fn end(&self) -> u32 {
        self.0.end()
    }
    #[wasm_bindgen(getter,js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}
#[wasm_bindgen]
pub struct BioSiftsUnpResidue(pub(crate) ck::BioSiftsUnpResidue);
crate::fields::value_field!(
    ck::BioSiftsUnpResidue,
    BioSiftsUnpResidue,
    BioSiftsUnpResidueInput,
    "BioSiftsUnpResidue"
);
#[wasm_bindgen]
impl BioSiftsUnpResidue {
    #[wasm_bindgen(constructor)]
    pub fn new(
        residue: Option<Ts<input::U8Input>>,
        accession_index: Option<Ts<input::U8Input>>,
        number: Option<Ts<input::U16Input>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioSiftsUnpResidue::default();
        Ok(Self(ck::BioSiftsUnpResidue::new(
            residue.map(input::read::<u8>).transpose()?,
            accession_index
                .map(input::read::<u8>)
                .transpose()?
                .unwrap_or(d.accession_index()),
            number
                .map(input::read::<u16>)
                .transpose()?
                .unwrap_or(d.number()),
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn residue(&self) -> Option<u8> {
        self.0.residue()
    }
    #[wasm_bindgen(getter,js_name=accessionIndex)]
    pub fn accession_index(&self) -> u8 {
        self.0.accession_index()
    }
    #[wasm_bindgen(getter)]
    pub fn number(&self) -> u16 {
        self.0.number()
    }
}
#[wasm_bindgen]
pub struct BioModelRow(pub(crate) ck::BioModelRow);
crate::fields::value_field!(
    ck::BioModelRow,
    BioModelRow,
    BioModelRowInput,
    "BioModelRow"
);
#[wasm_bindgen]
impl BioModelRow {
    #[wasm_bindgen(constructor)]
    pub fn new(
        chain_span: BioRowSpanInput,
        source_model_number: Option<Ts<input::I32Input>>,
    ) -> Result<Self, JsValue> {
        BioRowSpan::with_ref(&chain_span, |s| {
            Ok(Self(ck::BioModelRow::new(
                ck::BioRowSpan::new(s.start(), s.len()).map_err(BioStructureError::from_rust)?,
                source_model_number.map(input::read::<i32>).transpose()?,
            )))
        })
    }
    #[wasm_bindgen(getter,js_name=chainSpan)]
    pub fn chain_span(&self) -> Result<BioRowSpan, JsValue> {
        BioRowSpan::from_native(self.0.chain_span())
    }
    #[wasm_bindgen(getter,js_name=sourceModelNumber)]
    pub fn source_model_number(&self) -> Option<i32> {
        self.0.source_model_number()
    }
}
#[wasm_bindgen]
pub struct BioChainRow(pub(crate) ck::BioChainRow);
crate::fields::value_field!(
    ck::BioChainRow,
    BioChainRow,
    BioChainRowInput,
    "BioChainRow"
);
#[wasm_bindgen]
impl BioChainRow {
    #[wasm_bindgen(constructor)]
    pub fn new(
        model_id: Ts<input::U32Input>,
        entity_id: Option<Ts<input::U32Input>>,
        residue_span: BioRowSpanInput,
        kind: ChainKindInput,
        source: ChainSourceIdsInput,
    ) -> Result<Self, JsValue> {
        BioRowSpan::with_ref(&residue_span, |s| {
            Ok(Self(ck::BioChainRow::new(
                ck::BioModelId::new(input::read::<u32>(model_id)?),
                entity_id
                    .map(input::read::<u32>)
                    .transpose()?
                    .map(ck::BioEntityId::new),
                ck::BioRowSpan::new(s.start(), s.len()).map_err(BioStructureError::from_rust)?,
                <ck::ChainKind as Field>::input(kind)?,
                <ck::ChainSourceIds as Field>::input(source)?,
            )))
        })
    }
    #[wasm_bindgen(getter,js_name=modelId)]
    pub fn model_id(&self) -> u32 {
        self.0.model_id().value()
    }
    #[wasm_bindgen(getter,js_name=entityId)]
    pub fn entity_id(&self) -> Option<u32> {
        self.0.entity_id().map(ck::BioEntityId::value)
    }
    #[wasm_bindgen(getter,js_name=residueSpan)]
    pub fn residue_span(&self) -> Result<BioRowSpan, JsValue> {
        BioRowSpan::from_native(self.0.residue_span())
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> ChainKind {
        self.0.kind().into()
    }
    #[wasm_bindgen(getter)]
    pub fn source(&self) -> ChainSourceIds {
        ChainSourceIds(self.0.source().clone())
    }
}
#[wasm_bindgen]
pub struct BioResidueRow(pub(crate) ck::BioResidueRow);
crate::fields::value_field!(
    ck::BioResidueRow,
    BioResidueRow,
    BioResidueRowInput,
    "BioResidueRow"
);
#[wasm_bindgen]
impl BioResidueRow {
    #[wasm_bindgen(constructor)]
    pub fn new(
        chain_id: Ts<input::U32Input>,
        atom_span: BioRowSpanInput,
        name: Ts<input::TextInput>,
        residue_info_kind: ResidueInfoKindInput,
        entity_kind: EntityKindInput,
        entity_id: Option<Ts<input::U32Input>>,
        het_flag: Option<Ts<input::U8Input>>,
        source: ResidueSourceIdsInput,
        sifts_unp: BioSiftsUnpResidueInput,
    ) -> Result<Self, JsValue> {
        BioRowSpan::with_ref(&atom_span, |s| {
            Ok(Self(ck::BioResidueRow::new(
                ck::BioChainId::new(input::read::<u32>(chain_id)?),
                *s,
                ck::ResidueName::from_ascii(input::read::<String>(name)?.as_bytes()).ok_or_else(
                    || JsValue::from(js_sys::RangeError::new("name is not a valid ResidueName")),
                )?,
                <ck::ResidueInfoKind as Field>::input(residue_info_kind)?,
                <ck::EntityKind as Field>::input(entity_kind)?,
                entity_id
                    .map(input::read::<u32>)
                    .transpose()?
                    .map(ck::BioEntityId::new),
                het_flag.map(input::read::<u8>).transpose()?,
                <ck::ResidueSourceIds as Field>::input(source)?,
                <ck::BioSiftsUnpResidue as Field>::input(sifts_unp)?,
            )))
        })
    }
    #[wasm_bindgen(getter,js_name=chainId)]
    pub fn chain_id(&self) -> u32 {
        self.0.chain_id().value()
    }
    #[wasm_bindgen(getter,js_name=atomSpan)]
    pub fn atom_span(&self) -> Result<BioRowSpan, JsValue> {
        BioRowSpan::from_native(self.0.atom_span())
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> ResidueName {
        ResidueName(self.0.name())
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> ResidueKind {
        self.0.kind().into()
    }
    #[wasm_bindgen(getter,js_name=residueInfoKind)]
    pub fn residue_info_kind(&self) -> ResidueInfoKind {
        self.0.residue_info_kind().into()
    }
    #[wasm_bindgen(getter,js_name=entityKind)]
    pub fn entity_kind(&self) -> EntityKind {
        self.0.entity_kind().into()
    }
    #[wasm_bindgen(getter,js_name=entityId)]
    pub fn entity_id(&self) -> Option<u32> {
        self.0.entity_id().map(ck::BioEntityId::value)
    }
    #[wasm_bindgen(getter,js_name=hetFlag)]
    pub fn het_flag(&self) -> Option<u8> {
        self.0.het_flag()
    }
    #[wasm_bindgen(getter)]
    pub fn source(&self) -> ResidueSourceIds {
        ResidueSourceIds(self.0.source().clone())
    }
    #[wasm_bindgen(getter,js_name=siftsUnp)]
    pub fn sifts_unp(&self) -> BioSiftsUnpResidue {
        BioSiftsUnpResidue(self.0.sifts_unp())
    }
}
#[wasm_bindgen]
pub struct BioAtomRow(pub(crate) ck::BioAtomRow);
crate::fields::value_field!(ck::BioAtomRow, BioAtomRow, BioAtomRowInput, "BioAtomRow");
#[wasm_bindgen]
impl BioAtomRow {
    #[wasm_bindgen(constructor)]
    pub fn new(
        residue_id: Ts<input::U32Input>,
        name: AtomNameInput,
        element: crate::atoms::ElementInput,
        isotope_mass_number: Option<Ts<input::U16Input>>,
        altloc: Option<AltLocLabelInput>,
        formal_charge: Ts<input::I8Input>,
        calc_flag: BioCalcFlagInput,
        occupancy: Ts<input::F64Input>,
        b_iso: Ts<input::F64Input>,
        anisou: Ts<Anisou>,
        tls_group_id: Ts<input::I16Input>,
        fraction: Ts<input::F64Input>,
        source: AtomSourceIdsInput,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::BioAtomRow::new(
            ck::BioResidueId::new(input::read::<u32>(residue_id)?),
            <ck::AtomName as Field>::input(name)?,
            <ck::Element as Field>::input(element)?,
            isotope_mass_number.map(input::read::<u16>).transpose()?,
            altloc.map(<ck::AltLocLabel as Field>::input).transpose()?,
            input::read::<i8>(formal_charge)?,
            <ck::BioCalcFlag as Field>::input(calc_flag)?,
            input::read::<f64>(occupancy)?,
            input::read::<f64>(b_iso)?,
            input::decode(anisou)?.0,
            input::read::<i16>(tls_group_id)?,
            input::read::<f64>(fraction)?,
            <ck::AtomSourceIds as Field>::input(source)?,
        )))
    }
    #[wasm_bindgen(getter,js_name=residueId)]
    pub fn residue_id(&self) -> u32 {
        self.0.residue_id().value()
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> AtomName {
        AtomName(self.0.name())
    }
    #[wasm_bindgen(getter)]
    pub fn element(&self) -> crate::atoms::Element {
        crate::atoms::Element(self.0.element())
    }
    #[wasm_bindgen(getter,js_name=isotopeMassNumber)]
    pub fn isotope_mass_number(&self) -> Option<u16> {
        self.0.isotope_mass_number()
    }
    #[wasm_bindgen(getter)]
    pub fn altloc(&self) -> Option<AltLocLabel> {
        self.0.altloc().map(AltLocLabel)
    }
    #[wasm_bindgen(getter,js_name=formalCharge)]
    pub fn formal_charge(&self) -> i8 {
        self.0.formal_charge()
    }
    #[wasm_bindgen(getter,js_name=calcFlag)]
    pub fn calc_flag(&self) -> BioCalcFlag {
        self.0.calc_flag().into()
    }
    #[wasm_bindgen(getter)]
    pub fn occupancy(&self) -> f64 {
        self.0.occupancy()
    }
    #[wasm_bindgen(getter,js_name=bIso)]
    pub fn b_iso(&self) -> f64 {
        self.0.b_iso()
    }
    #[wasm_bindgen(getter)]
    pub fn anisou(&self) -> Result<Ts<Anisou>, JsValue> {
        Anisou(*self.0.anisou())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=tlsGroupId)]
    pub fn tls_group_id(&self) -> i16 {
        self.0.tls_group_id()
    }
    #[wasm_bindgen(getter)]
    pub fn fraction(&self) -> f64 {
        self.0.fraction()
    }
    #[wasm_bindgen(getter)]
    pub fn source(&self) -> AtomSourceIds {
        AtomSourceIds(self.0.source().clone())
    }
}
