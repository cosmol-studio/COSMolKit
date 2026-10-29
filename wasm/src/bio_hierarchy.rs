//! Structural reads keep the native COW hierarchy owner.
use crate::{fields::Field, group_values::Point3, input, storage::NativeField};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
ck::__ck_configuration_BioPdbReadParams!(crate::project_parameters, ck::BioPdbReadParams, {}, {
    BioPdbReadParams
});
crate::enums::enum_value!(BioCoordinateFormat,"BioCoordinateFormat"{Unknown="unknown",Detect="detect",Pdb="pdb",Mmcif="mmcif",Mmjson="mmjson",ChemComp="chem_comp"});
crate::enums::enum_value!(ChainKind,"ChainKind"{Protein="protein",Dna="dna",Rna="rna",ProteinDnaComplex="protein_dna_complex",ProteinRnaComplex="protein_rna_complex",LigandOnly="ligand_only",WaterOnly="water_only",Mixed="mixed",Unknown="unknown"});
crate::enums::enum_value!(ResidueKind,"ResidueKind"{AminoAcid="amino_acid",Dna="dna",Rna="rna",Saccharide="saccharide",Water="water",Buffer="buffer",Ligand="ligand",Unknown="unknown"});
#[wasm_bindgen]
pub struct BioStructure(pub(crate) ck::BioStructure);
crate::fields::value_field!(
    ck::BioStructure,
    BioStructure,
    BioStructureInput,
    "BioStructure"
);
crate::fields::value_array!(
    ck::BioStructure,
    BioStructureArray,
    "ReadonlyArray<BioStructure>"
);
#[wasm_bindgen]
pub struct BioModel(pub(crate) ck::BioModel);
crate::fields::value_field!(ck::BioModel, BioModel, BioModelInput, "BioModel");
crate::fields::value_array!(ck::BioModel, BioModelArray, "ReadonlyArray<BioModel>");
#[wasm_bindgen]
pub struct BioChain(pub(crate) ck::BioChain);
crate::fields::value_field!(ck::BioChain, BioChain, BioChainInput, "BioChain");
crate::fields::value_array!(ck::BioChain, BioChainArray, "ReadonlyArray<BioChain>");
#[wasm_bindgen]
pub struct BioResidue(pub(crate) ck::BioResidue);
crate::fields::value_field!(ck::BioResidue, BioResidue, BioResidueInput, "BioResidue");
crate::fields::value_array!(ck::BioResidue, BioResidueArray, "ReadonlyArray<BioResidue>");
#[wasm_bindgen]
pub struct BioAtom(pub(crate) ck::BioAtom);
crate::fields::value_field!(ck::BioAtom, BioAtom, BioAtomInput, "BioAtom");
crate::fields::value_array!(ck::BioAtom, BioAtomArray, "ReadonlyArray<BioAtom>");
#[wasm_bindgen(js_name=bioFromPdb)]
pub fn bio_from_pdb(
    text: Ts<input::TextInput>,
    params: Option<BioPdbReadParamsInput>,
) -> Result<BioStructure, JsValue> {
    BioStructure::from_pdb(text, params)
}
#[wasm_bindgen(js_name=bioFromMmcif)]
pub fn bio_from_mmcif(text: Ts<input::TextInput>) -> Result<BioStructure, JsValue> {
    BioStructure::from_mmcif(text)
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=fromPdbWithParams)]
    pub fn from_pdb_with_params(
        text: Ts<input::TextInput>,
        params: BioPdbReadParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_pdb(text, Some(params))
    }
    #[wasm_bindgen(getter)]
    pub fn coordinates(&self) -> crate::bio_geometry::BioCoordinateBlock {
        crate::bio_geometry::BioCoordinateBlock(self.0.coordinates().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn crystal(&self) -> Option<crate::bio_geometry::BioCrystalInfo> {
        self.0
            .crystal()
            .cloned()
            .map(crate::bio_geometry::BioCrystalInfo)
    }
    #[wasm_bindgen(getter,js_name=ncsOperators)]
    pub fn ncs_operators(&self) -> Result<crate::bio_geometry::BioNcsOperators, JsValue> {
        <Vec<ck::BioNcsOperator> as Field>::get(NativeField::owned(self.0.ncs_operators().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn assemblies(&self) -> Result<crate::bio_assemblies::BioAssemblies, JsValue> {
        <Vec<ck::BioAssembly> as Field>::get(NativeField::owned(self.0.assemblies().to_vec()))
    }
    #[wasm_bindgen(getter,js_name=hasOrigx)]
    pub fn has_origx(&self) -> bool {
        self.0.has_origx()
    }
    #[wasm_bindgen(getter)]
    pub fn origx(&self) -> crate::bio_geometry::BioTransform {
        crate::bio_geometry::BioTransform(*self.0.origx())
    }
    #[wasm_bindgen(getter)]
    pub fn resolution(&self) -> f64 {
        self.0.resolution()
    }
    #[wasm_bindgen(getter,js_name=terStatus)]
    pub fn ter_status(&self) -> u8 {
        self.0.ter_status()
    }
    #[wasm_bindgen(js_name=fromPdb)]
    pub fn from_pdb(
        text: Ts<input::TextInput>,
        params: Option<BioPdbReadParamsInput>,
    ) -> Result<Self, JsValue> {
        let text = input::read::<String>(text)?;
        BioPdbReadParams::with_input(params.as_ref(), |p| {
            ck::BioStructure::from_pdb_with_params(&text, p)
                .map(Self)
                .map_err(crate::pdb_read_errors::BioPdbReadError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=fromMmcif)]
    pub fn from_mmcif(text: Ts<input::TextInput>) -> Result<Self, JsValue> {
        ck::BioStructure::from_mmcif(&input::read::<String>(text)?)
            .map(Self)
            .map_err(crate::bio_read_errors::BioMmcifReadError::from_rust)
    }
    pub fn copy(&self) -> Self {
        Self(self.0.clone())
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name().into()
    }
    #[wasm_bindgen(getter,js_name=inputFormat)]
    pub fn input_format(&self) -> BioCoordinateFormat {
        self.0.input_format().into()
    }
    #[wasm_bindgen(js_name=atomPosition)]
    pub fn atom_position(&self, id: Ts<input::U32Input>) -> Result<Option<Ts<Point3>>, JsValue> {
        self.0
            .atom_position(ck::BioAtomId::new(input::read::<u32>(id)?))
            .map(|v| {
                Point3(v)
                    .into_ts()
                    .map_err(JsError::from)
                    .map_err(Into::into)
            })
            .transpose()
    }
    #[wasm_bindgen(getter)]
    pub fn models(&self) -> Result<BioModelArray, JsValue> {
        <Vec<ck::BioModel> as Field>::get(NativeField::owned(self.0.models()))
    }
    pub fn model(&self, id: Ts<input::U32Input>) -> Result<Option<BioModel>, JsValue> {
        Ok(self
            .0
            .model(ck::BioModelId::new(input::read::<u32>(id)?))
            .map(BioModel))
    }
    #[wasm_bindgen(getter,js_name=numModels)]
    pub fn num_models(&self) -> usize {
        self.0.num_models()
    }
    #[wasm_bindgen(getter)]
    pub fn chains(&self) -> Result<BioChainArray, JsValue> {
        <Vec<ck::BioChain> as Field>::get(NativeField::owned(self.0.chains()))
    }
    pub fn chain(&self, id: Ts<input::U32Input>) -> Result<Option<BioChain>, JsValue> {
        Ok(self
            .0
            .chain(ck::BioChainId::new(input::read::<u32>(id)?))
            .map(BioChain))
    }
    #[wasm_bindgen(getter,js_name=numChains)]
    pub fn num_chains(&self) -> usize {
        self.0.num_chains()
    }
    #[wasm_bindgen(getter)]
    pub fn residues(&self) -> Result<BioResidueArray, JsValue> {
        <Vec<ck::BioResidue> as Field>::get(NativeField::owned(self.0.residues()))
    }
    pub fn residue(&self, id: Ts<input::U32Input>) -> Result<Option<BioResidue>, JsValue> {
        Ok(self
            .0
            .residue(ck::BioResidueId::new(input::read::<u32>(id)?))
            .map(BioResidue))
    }
    #[wasm_bindgen(getter,js_name=numResidues)]
    pub fn num_residues(&self) -> usize {
        self.0.num_residues()
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Result<BioAtomArray, JsValue> {
        <Vec<ck::BioAtom> as Field>::get(NativeField::owned(self.0.atoms()))
    }
    pub fn atom(&self, id: Ts<input::U32Input>) -> Result<Option<BioAtom>, JsValue> {
        Ok(self
            .0
            .atom(ck::BioAtomId::new(input::read::<u32>(id)?))
            .map(BioAtom))
    }
    #[wasm_bindgen(getter,js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }
}
#[wasm_bindgen]
impl BioModel {
    pub fn value(&self) -> crate::bio_rows::BioModelRow {
        crate::bio_rows::BioModelRow(*self.0.value())
    }
    #[wasm_bindgen(getter,js_name=chainSpan)]
    pub fn chain_span(&self) -> Result<crate::bio_rows::BioRowSpan, JsValue> {
        crate::bio_rows::BioRowSpan::from_native(self.0.chain_span())
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.0.id().value()
    }
    #[wasm_bindgen(getter,js_name=numChains)]
    pub fn num_chains(&self) -> usize {
        self.0.num_chains()
    }
    #[wasm_bindgen(getter)]
    pub fn chains(&self) -> Result<BioChainArray, JsValue> {
        <Vec<ck::BioChain> as Field>::get(NativeField::owned(self.0.chains()))
    }
}
#[wasm_bindgen]
impl BioChain {
    pub fn value(&self) -> crate::bio_rows::BioChainRow {
        crate::bio_rows::BioChainRow(self.0.value().clone())
    }
    #[wasm_bindgen(getter,js_name=residueSpan)]
    pub fn residue_span(&self) -> Result<crate::bio_rows::BioRowSpan, JsValue> {
        crate::bio_rows::BioRowSpan::from_native(self.0.residue_span())
    }
    #[wasm_bindgen(getter)]
    pub fn source(&self) -> crate::bio_source_values::ChainSourceIds {
        crate::bio_source_values::ChainSourceIds(self.0.source().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.0.id().value()
    }
    #[wasm_bindgen(getter,js_name=numResidues)]
    pub fn num_residues(&self) -> usize {
        self.0.num_residues()
    }
    #[wasm_bindgen(getter)]
    pub fn residues(&self) -> Result<BioResidueArray, JsValue> {
        <Vec<ck::BioResidue> as Field>::get(NativeField::owned(self.0.residues()))
    }
}
#[wasm_bindgen]
impl BioResidue {
    pub fn value(&self) -> crate::bio_rows::BioResidueRow {
        crate::bio_rows::BioResidueRow(self.0.value().clone())
    }
    #[wasm_bindgen(getter,js_name=atomSpan)]
    pub fn atom_span(&self) -> Result<crate::bio_rows::BioRowSpan, JsValue> {
        crate::bio_rows::BioRowSpan::from_native(self.0.atom_span())
    }
    #[wasm_bindgen(getter,js_name=residueInfoKind)]
    pub fn residue_info_kind(&self) -> crate::residue_values::ResidueInfoKind {
        self.0.residue_info_kind().into()
    }
    #[wasm_bindgen(getter,js_name=entityKind)]
    pub fn entity_kind(&self) -> crate::bio_rows::EntityKind {
        self.0.entity_kind().into()
    }
    #[wasm_bindgen(getter,js_name=siftsUnp)]
    pub fn sifts_unp(&self) -> crate::bio_rows::BioSiftsUnpResidue {
        crate::bio_rows::BioSiftsUnpResidue(self.0.sifts_unp())
    }
    #[wasm_bindgen(getter)]
    pub fn source(&self) -> crate::bio_source_values::ResidueSourceIds {
        crate::bio_source_values::ResidueSourceIds(self.0.source().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.0.id().value()
    }
    #[wasm_bindgen(getter,js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Result<BioAtomArray, JsValue> {
        <Vec<ck::BioAtom> as Field>::get(NativeField::owned(self.0.atoms()))
    }
}
#[wasm_bindgen]
impl BioModel {
    #[wasm_bindgen(getter,js_name=sourceModelNumber)]
    pub fn source_model_number(&self) -> Option<i32> {
        self.0.source_model_number()
    }
}
#[wasm_bindgen]
impl BioChain {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> ChainKind {
        self.0.kind().into()
    }
    #[wasm_bindgen(getter,js_name=modelId)]
    pub fn model_id(&self) -> u32 {
        self.0.model_id().value()
    }
    #[wasm_bindgen(getter,js_name=entityId)]
    pub fn entity_id(&self) -> Option<u32> {
        self.0.entity_id().map(|v| v.value())
    }
}
#[wasm_bindgen]
impl BioResidue {
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name().as_str().into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> ResidueKind {
        self.0.kind().into()
    }
    #[wasm_bindgen(getter,js_name=chainId)]
    pub fn chain_id(&self) -> u32 {
        self.0.chain_id().value()
    }
    #[wasm_bindgen(getter,js_name=entityId)]
    pub fn entity_id(&self) -> Option<u32> {
        self.0.entity_id().map(|v| v.value())
    }
    #[wasm_bindgen(getter,js_name=hetFlag)]
    pub fn het_flag(&self) -> Option<u8> {
        self.0.het_flag()
    }
}
#[wasm_bindgen]
impl BioAtom {
    pub fn value(&self) -> crate::bio_rows::BioAtomRow {
        crate::bio_rows::BioAtomRow(self.0.value().clone())
    }
    #[wasm_bindgen(getter,js_name=calculationFlag)]
    pub fn calculation_flag(&self) -> crate::bio_rows::BioCalcFlag {
        self.0.calculation_flag().into()
    }
    #[wasm_bindgen(getter)]
    pub fn anisou(&self) -> Result<Ts<crate::bio_rows::Anisou>, JsValue> {
        crate::bio_rows::Anisou(*self.0.anisou())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> crate::bio_source_values::AtomName {
        crate::bio_source_values::AtomName(self.0.name())
    }
    #[wasm_bindgen(getter)]
    pub fn source(&self) -> crate::bio_source_values::AtomSourceIds {
        crate::bio_source_values::AtomSourceIds(*self.0.source())
    }
    #[wasm_bindgen(getter)]
    pub fn altloc(&self) -> Option<crate::bio_source_values::AltLocLabel> {
        self.0.altloc().map(crate::bio_source_values::AltLocLabel)
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.0.id().value()
    }
    #[wasm_bindgen(getter,js_name=residueId)]
    pub fn residue_id(&self) -> u32 {
        self.0.residue_id().value()
    }
    #[wasm_bindgen(getter)]
    pub fn element(&self) -> crate::atoms::Element {
        crate::atoms::Element(self.0.element())
    }
    #[wasm_bindgen(getter)]
    pub fn position(&self) -> Result<Ts<Point3>, JsValue> {
        Point3(self.0.position())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=atomicNumber)]
    pub fn atomic_number(&self) -> u8 {
        self.0.atomic_number()
    }
    #[wasm_bindgen(getter,js_name=isotopeMassNumber)]
    pub fn isotope_mass_number(&self) -> Option<u16> {
        self.0.isotope_mass_number()
    }
    #[wasm_bindgen(getter,js_name=formalCharge)]
    pub fn formal_charge(&self) -> i8 {
        self.0.formal_charge()
    }
    #[wasm_bindgen(getter)]
    pub fn occupancy(&self) -> f64 {
        self.0.occupancy()
    }
    #[wasm_bindgen(getter,js_name=bIso)]
    pub fn b_iso(&self) -> f64 {
        self.0.b_iso()
    }
    #[wasm_bindgen(getter,js_name=tlsGroupId)]
    pub fn tls_group_id(&self) -> i16 {
        self.0.tls_group_id()
    }
    #[wasm_bindgen(getter)]
    pub fn fraction(&self) -> f64 {
        self.0.fraction()
    }
}
