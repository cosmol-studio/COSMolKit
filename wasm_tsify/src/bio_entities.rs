//! Native entity annotations and source sequences.
use crate::{
    bio_hierarchy::BioStructure,
    bio_rows::{EntityKind, EntityKindInput, PolymerKind, PolymerKindInput},
    bio_source_values::{PdbSeqId, PdbSeqIdInput},
    fields::{Field, TextArray},
    input,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub struct EntitySourceIds(ck::EntitySourceIds);
crate::fields::value_field!(
    ck::EntitySourceIds,
    EntitySourceIds,
    EntitySourceIdsInput,
    "EntitySourceIds"
);
#[wasm_bindgen]
impl EntitySourceIds {
    #[wasm_bindgen(constructor)]
    pub fn new(source_entity_id: Ts<input::TextInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::EntitySourceIds::new(input::read::<String>(
            source_entity_id,
        )?)))
    }
    #[wasm_bindgen(getter,js_name=sourceEntityId)]
    pub fn source_entity_id(&self) -> String {
        self.0.source_entity_id().into()
    }
}
#[wasm_bindgen]
pub struct BioEntityDbRef(ck::BioEntityDbRef);
crate::fields::value_field!(
    ck::BioEntityDbRef,
    BioEntityDbRef,
    BioEntityDbRefInput,
    "BioEntityDbRef"
);
crate::fields::value_array!(
    ck::BioEntityDbRef,
    BioEntityDbRefs,
    "ReadonlyArray<BioEntityDbRef>"
);
#[wasm_bindgen]
impl BioEntityDbRef {
    #[wasm_bindgen(constructor)]
    pub fn new(
        db_name: Option<Ts<input::TextInput>>,
        accession_code: Option<Ts<input::TextInput>>,
        id_code: Option<Ts<input::TextInput>>,
        isoform: Option<Ts<input::TextInput>>,
        seq_begin: Option<PdbSeqIdInput>,
        seq_end: Option<PdbSeqIdInput>,
        db_begin: Option<PdbSeqIdInput>,
        db_end: Option<PdbSeqIdInput>,
        label_seq_begin: Option<Ts<input::I32Input>>,
        label_seq_end: Option<Ts<input::I32Input>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioEntityDbRef::default();
        Ok(Self(ck::BioEntityDbRef {
            db_name: db_name
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.db_name),
            accession_code: accession_code
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.accession_code),
            id_code: id_code
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.id_code),
            isoform: isoform
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.isoform),
            seq_begin: seq_begin
                .map(<ck::PdbSeqId as Field>::input)
                .transpose()?
                .unwrap_or(d.seq_begin),
            seq_end: seq_end
                .map(<ck::PdbSeqId as Field>::input)
                .transpose()?
                .unwrap_or(d.seq_end),
            db_begin: db_begin
                .map(<ck::PdbSeqId as Field>::input)
                .transpose()?
                .unwrap_or(d.db_begin),
            db_end: db_end
                .map(<ck::PdbSeqId as Field>::input)
                .transpose()?
                .unwrap_or(d.db_end),
            label_seq_begin: label_seq_begin.map(input::read::<i32>).transpose()?,
            label_seq_end: label_seq_end.map(input::read::<i32>).transpose()?,
        }))
    }
    #[wasm_bindgen(getter,js_name=dbName)]
    pub fn db_name(&self) -> String {
        self.0.db_name.clone()
    }
    #[wasm_bindgen(getter,js_name=accessionCode)]
    pub fn accession_code(&self) -> String {
        self.0.accession_code.clone()
    }
    #[wasm_bindgen(getter,js_name=idCode)]
    pub fn id_code(&self) -> String {
        self.0.id_code.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn isoform(&self) -> String {
        self.0.isoform.clone()
    }
    #[wasm_bindgen(getter,js_name=seqBegin)]
    pub fn seq_begin(&self) -> PdbSeqId {
        PdbSeqId(self.0.seq_begin)
    }
    #[wasm_bindgen(getter,js_name=seqEnd)]
    pub fn seq_end(&self) -> PdbSeqId {
        PdbSeqId(self.0.seq_end)
    }
    #[wasm_bindgen(getter,js_name=dbBegin)]
    pub fn db_begin(&self) -> PdbSeqId {
        PdbSeqId(self.0.db_begin)
    }
    #[wasm_bindgen(getter,js_name=dbEnd)]
    pub fn db_end(&self) -> PdbSeqId {
        PdbSeqId(self.0.db_end)
    }
    #[wasm_bindgen(getter,js_name=labelSeqBegin)]
    pub fn label_seq_begin(&self) -> Option<i32> {
        self.0.label_seq_begin
    }
    #[wasm_bindgen(getter,js_name=labelSeqEnd)]
    pub fn label_seq_end(&self) -> Option<i32> {
        self.0.label_seq_end
    }
}
#[wasm_bindgen]
pub struct BioEntityRow(ck::BioEntityRow);
crate::fields::value_field!(
    ck::BioEntityRow,
    BioEntityRow,
    BioEntityRowInput,
    "BioEntityRow"
);
crate::fields::value_array!(
    ck::BioEntityRow,
    BioEntityRows,
    "ReadonlyArray<BioEntityRow>"
);
#[wasm_bindgen]
impl BioEntityRow {
    #[wasm_bindgen(constructor)]
    pub fn new(
        kind: EntityKindInput,
        polymer_kind: PolymerKindInput,
        reflects_microhetero: Ts<input::BoolInput>,
        full_sequence: Ts<TextArray>,
        dbrefs: BioEntityDbRefs,
        sifts_unp_accessions: Ts<TextArray>,
        subchains: Ts<TextArray>,
        source: EntitySourceIdsInput,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::BioEntityRow::new(
            <ck::EntityKind as Field>::input(kind)?,
            <ck::PolymerKind as Field>::input(polymer_kind)?,
            input::read::<bool>(reflects_microhetero)?,
            <Vec<String> as Field>::input(full_sequence)?,
            <Vec<ck::BioEntityDbRef> as Field>::input(dbrefs)?,
            <Vec<String> as Field>::input(sifts_unp_accessions)?,
            <Vec<String> as Field>::input(subchains)?,
            <ck::EntitySourceIds as Field>::input(source)?,
        )))
    }
    #[wasm_bindgen(js_name=firstMon)]
    pub fn first_mon(mon_list: Ts<input::TextInput>) -> Result<String, JsValue> {
        Ok(ck::BioEntityRow::first_mon(&input::read::<String>(mon_list)?).into())
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> EntityKind {
        self.0.kind().into()
    }
    #[wasm_bindgen(getter,js_name=polymerKind)]
    pub fn polymer_kind(&self) -> PolymerKind {
        self.0.polymer_kind().into()
    }
    #[wasm_bindgen(getter,js_name=reflectsMicrohetero)]
    pub fn reflects_microhetero(&self) -> bool {
        self.0.reflects_microhetero()
    }
    #[wasm_bindgen(getter,js_name=fullSequence)]
    pub fn full_sequence(&self) -> Result<Ts<TextArray>, JsValue> {
        <Vec<String> as Field>::get(NativeField::owned(self.0.full_sequence().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn dbrefs(&self) -> Result<BioEntityDbRefs, JsValue> {
        <Vec<ck::BioEntityDbRef> as Field>::get(NativeField::owned(self.0.dbrefs().to_vec()))
    }
    #[wasm_bindgen(getter,js_name=siftsUnpAccessions)]
    pub fn sifts_unp_accessions(&self) -> Result<Ts<TextArray>, JsValue> {
        <Vec<String> as Field>::get(NativeField::owned(self.0.sifts_unp_accessions().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn subchains(&self) -> Result<Ts<TextArray>, JsValue> {
        <Vec<String> as Field>::get(NativeField::owned(self.0.subchains().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn source(&self) -> EntitySourceIds {
        EntitySourceIds(self.0.source().clone())
    }
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "readonly [number, BioEntityRow]")]
    pub type EntityMatch;
}
fn entity_match(value: Option<(ck::BioEntityId, &ck::BioEntityRow)>) -> Option<EntityMatch> {
    value.map(|(id, row)| {
        let array = js_sys::Array::new();
        array.push(&JsValue::from(id.value()));
        array.push(&BioEntityRow(row.clone()).into());
        js_sys::Object::freeze(&array);
        JsCast::unchecked_into(array)
    })
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(getter,js_name=numEntities)]
    pub fn num_entities(&self) -> usize {
        self.0.num_entities()
    }
    #[wasm_bindgen(getter)]
    pub fn entities(&self) -> Result<BioEntityRows, JsValue> {
        <Vec<ck::BioEntityRow> as Field>::get(NativeField::owned(self.0.entities().to_vec()))
    }
    #[wasm_bindgen(js_name=findEntity)]
    pub fn find_entity(
        &self,
        source_id: Ts<input::TextInput>,
    ) -> Result<Option<EntityMatch>, JsValue> {
        Ok(entity_match(
            self.0.find_entity(&input::read::<String>(source_id)?),
        ))
    }
    #[wasm_bindgen(js_name=findEntityOfSubchain)]
    pub fn find_entity_of_subchain(
        &self,
        subchain: Ts<input::TextInput>,
    ) -> Result<Option<EntityMatch>, JsValue> {
        Ok(entity_match(self.0.find_entity_of_subchain(
            &input::read::<String>(subchain)?,
        )))
    }
}
