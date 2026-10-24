//! Protein projections retain native owners; mutations use native COW.
use crate::{
    Molecule,
    bio_hierarchy::{
        BioCoordinateFormat, BioPdbReadParams, BioPdbReadParamsInput, BioStructure, ChainKind,
        ResidueKind,
    },
    bio_molecule::{BioMoleculeParams, BioMoleculeParamsInput},
    bio_operations::{
        BioOperationError, BioSelection, BioSelectionInput, BioSelectionMatchError,
        ProteinProjectionError,
    },
    bio_text::{BioReadParams, BioReadParamsInput, ProteinReadError},
    fields::Field,
    group_values::Point3,
    input,
    storage::NativeField,
};
use cosmolkit as ck;
use std::sync::Arc;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct Protein(Arc<ck::Protein>);
#[wasm_bindgen(js_name=proteinFromPdb)]
pub fn protein_from_pdb(
    text: Ts<input::TextInput>,
    params: Option<BioPdbReadParamsInput>,
) -> Result<Protein, JsValue> {
    Protein::from_pdb(text, params)
}
#[wasm_bindgen(js_name=proteinFromMmcif)]
pub fn protein_from_mmcif(text: Ts<input::TextInput>) -> Result<Protein, JsValue> {
    Protein::from_mmcif(text)
}
crate::fields::value_field!(Arc<ck::Protein>, Protein, ProteinInput, "Protein");
#[wasm_bindgen]
pub struct ProteinSelectionSummary(ck::ProteinSelectionSummary);
#[wasm_bindgen]
impl ProteinSelectionSummary {
    #[wasm_bindgen(getter)]
    pub fn chains(&self) -> usize {
        self.0.chains
    }
    #[wasm_bindgen(getter)]
    pub fn residues(&self) -> usize {
        self.0.residues
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> usize {
        self.0.atoms
    }
}
#[wasm_bindgen]
impl BioStructure {
    pub fn protein(&self) -> Result<Protein, JsValue> {
        self.0
            .protein()
            .map(|v| Protein(Arc::new(v)))
            .map_err(ProteinProjectionError::from_rust)
    }
}
#[wasm_bindgen]
impl Protein {
    #[wasm_bindgen(js_name=toMoleculeWithParams)]
    pub fn to_molecule_with_params(
        &self,
        params: BioMoleculeParamsInput,
    ) -> Result<Molecule, JsValue> {
        self.to_molecule(Some(params))
    }
    #[wasm_bindgen(js_name=fromTextWithParams)]
    pub fn from_text_with_params(
        text: Ts<input::TextInput>,
        params: BioReadParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_text(text, Some(params))
    }
    #[wasm_bindgen(js_name=fromPdbWithParams)]
    pub fn from_pdb_with_params(
        text: Ts<input::TextInput>,
        params: BioPdbReadParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_pdb(text, Some(params))
    }
    #[wasm_bindgen(js_name=fromPdb)]
    pub fn from_pdb(
        text: Ts<input::TextInput>,
        params: Option<BioPdbReadParamsInput>,
    ) -> Result<Self, JsValue> {
        let text = input::read::<String>(text)?;
        BioPdbReadParams::with_input(params.as_ref(), |p| {
            ck::Protein::from_pdb_with_params(&text, p)
                .map(|v| Self(Arc::new(v)))
                .map_err(ProteinReadError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=fromText)]
    pub fn from_text(
        text: Ts<input::TextInput>,
        params: Option<BioReadParamsInput>,
    ) -> Result<Self, JsValue> {
        let text = input::read::<String>(text)?;
        BioReadParams::with_input(params.as_ref(), |p| {
            ck::Protein::from_text_with_params(&text, p)
                .map(|v| Self(Arc::new(v)))
                .map_err(ProteinReadError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=fromMmcif)]
    pub fn from_mmcif(text: Ts<input::TextInput>) -> Result<Self, JsValue> {
        ck::Protein::from_mmcif(&input::read::<String>(text)?)
            .map(|v| Self(Arc::new(v)))
            .map_err(ProteinReadError::from_rust)
    }
    pub fn copy(&self) -> Self {
        Self(self.0.clone())
    }
    #[wasm_bindgen(js_name=toMolecule)]
    pub fn to_molecule(&self, params: Option<BioMoleculeParamsInput>) -> Result<Molecule, JsValue> {
        BioMoleculeParams::with_input(params.as_ref(), |p| {
            self.0
                .to_molecule_with_params(p)
                .map(Molecule)
                .map_err(crate::bio_conversion_errors::BioMoleculeError::from_rust)
        })
    }
    #[wasm_bindgen(getter,js_name=inputFormat)]
    pub fn input_format(&self) -> BioCoordinateFormat {
        self.0.input_format().into()
    }
    #[wasm_bindgen(js_name=asBioStructure)]
    pub fn as_bio_structure(&self) -> BioStructure {
        BioStructure(self.0.as_bio_structure().clone())
    }
    #[wasm_bindgen(js_name=intoBioStructure)]
    pub fn into_bio_structure(&self) -> BioStructure {
        BioStructure(self.0.as_ref().clone().into_bio_structure())
    }
    #[wasm_bindgen(getter,js_name=numModels)]
    pub fn num_models(&self) -> usize {
        self.0.num_models()
    }
    #[wasm_bindgen(getter,js_name=numChains)]
    pub fn num_chains(&self) -> usize {
        self.0.num_chains()
    }
    #[wasm_bindgen(getter,js_name=numResidues)]
    pub fn num_residues(&self) -> usize {
        self.0.num_residues()
    }
    #[wasm_bindgen(getter,js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }
    #[wasm_bindgen(js_name=selectionSummary)]
    pub fn selection_summary(&self) -> ProteinSelectionSummary {
        ProteinSelectionSummary(self.0.selection_summary())
    }
    #[wasm_bindgen(js_name=withTranslatedCoordinates)]
    pub fn with_translated_coordinates(&self, offset: Ts<Point3>) -> Result<Self, JsValue> {
        self.0
            .with_translated_coordinates(input::decode(offset)?.0)
            .map(|v| Self(Arc::new(v)))
            .map_err(BioOperationError::from_rust)
    }
    #[wasm_bindgen(js_name=translate_)]
    pub fn translate_(&mut self, offset: Ts<Point3>) -> Result<(), JsValue> {
        Arc::make_mut(&mut self.0)
            .translate_(input::decode(offset)?.0)
            .map_err(BioOperationError::from_rust)
    }
    #[wasm_bindgen(js_name=withSelection)]
    pub fn with_selection(&self, selection: BioSelectionInput) -> Result<Self, JsValue> {
        BioSelection::with_ref(&selection, |s| {
            self.0
                .with_selection(s)
                .map(|v| Self(Arc::new(v)))
                .map_err(BioOperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=retainSelection_)]
    pub fn retain_selection_(&mut self, selection: BioSelectionInput) -> Result<(), JsValue> {
        BioSelection::with_ref(&selection, |s| {
            Arc::make_mut(&mut self.0)
                .retain_selection_(s)
                .map_err(BioOperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=selectedAtomIds)]
    pub fn selected_atom_ids(&self, selection: BioSelectionInput) -> Result<Vec<u32>, JsValue> {
        BioSelection::with_ref(&selection, |s| {
            self.0
                .selected_atom_ids(s)
                .map(|v| v.into_iter().map(ck::BioAtomId::value).collect())
                .map_err(BioSelectionMatchError::from_rust)
        })
    }
    #[wasm_bindgen(getter)]
    pub fn chains(&self) -> Result<ProteinChains, JsValue> {
        <Vec<ChainOwner> as Field>::get(NativeField::owned(
            self.0
                .chains()
                .map(|v| ChainOwner(self.0.clone(), v.id().value() as usize))
                .collect(),
        ))
    }
    #[wasm_bindgen(getter)]
    pub fn residues(&self) -> Result<ProteinResidues, JsValue> {
        <Vec<ResidueOwner> as Field>::get(NativeField::owned(
            self.0
                .residues()
                .map(|v| ResidueOwner(self.0.clone(), v.id().value() as usize))
                .collect(),
        ))
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Result<ProteinAtoms, JsValue> {
        <Vec<AtomOwner> as Field>::get(NativeField::owned(
            self.0
                .atoms()
                .map(|v| AtomOwner(self.0.clone(), v.id().value() as usize))
                .collect(),
        ))
    }
    pub fn chain(&self, index: Ts<input::IndexInput>) -> Result<Option<ProteinChainRef>, JsValue> {
        Ok(self
            .0
            .chain(input::read::<usize>(index)?)
            .map(|v| ProteinChainRef(ChainOwner(self.0.clone(), v.id().value() as usize))))
    }
}

// The native iterators implement constant-time nth. Indices originate only
// from the retained owner's canonical iterators; no lifetime extension.
macro_rules! retained {
    ($class:ident,$owner:ident,$input:ident,$array:ident,$iter:ident,$ts:literal,$ats:literal) => {
        #[derive(Clone)]
        pub(crate) struct $owner(Arc<ck::Protein>, usize);
        #[wasm_bindgen]
        pub struct $class($owner);
        crate::fields::value_field!($owner, $class, $input, $ts);
        crate::fields::value_array!($owner, $array, $ats);
        impl $class {
            fn view(&self) -> ck::$class<'_> {
                self.0
                    .0
                    .$iter()
                    .nth(self.0.1)
                    .expect("private index belongs to retained native protein")
            }
        }
    };
}
retained!(
    ProteinChainRef,
    ChainOwner,
    ProteinChainInput,
    ProteinChains,
    chains,
    "ProteinChainRef",
    "ReadonlyArray<ProteinChainRef>"
);
retained!(
    ProteinResidueRef,
    ResidueOwner,
    ProteinResidueInput,
    ProteinResidues,
    residues,
    "ProteinResidueRef",
    "ReadonlyArray<ProteinResidueRef>"
);
retained!(
    ProteinAtomRef,
    AtomOwner,
    ProteinAtomInput,
    ProteinAtoms,
    atoms,
    "ProteinAtomRef",
    "ReadonlyArray<ProteinAtomRef>"
);
#[wasm_bindgen]
impl ProteinChainRef {
    #[wasm_bindgen(getter)]
    pub fn row(&self) -> crate::bio_rows::BioChainRow {
        crate::bio_rows::BioChainRow(self.view().row().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.view().id().value()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> ChainKind {
        self.view().kind().into()
    }
    #[wasm_bindgen(getter)]
    pub fn source(&self) -> crate::bio_source_values::ChainSourceIds {
        crate::bio_source_values::ChainSourceIds(self.view().source().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn residues(&self) -> Result<ProteinResidues, JsValue> {
        <Vec<ResidueOwner> as Field>::get(NativeField::owned(
            self.view()
                .residues()
                .map(|v| ResidueOwner(self.0.0.clone(), v.id().value() as usize))
                .collect(),
        ))
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Result<ProteinAtoms, JsValue> {
        <Vec<AtomOwner> as Field>::get(NativeField::owned(
            self.view()
                .atoms()
                .map(|v| AtomOwner(self.0.0.clone(), v.id().value() as usize))
                .collect(),
        ))
    }
}
#[wasm_bindgen]
impl ProteinResidueRef {
    #[wasm_bindgen(getter)]
    pub fn row(&self) -> crate::bio_rows::BioResidueRow {
        crate::bio_rows::BioResidueRow(self.view().row().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn info(&self) -> crate::residue_values::ResidueInfo {
        crate::residue_values::ResidueInfo(self.view().info())
    }
    #[wasm_bindgen(getter)]
    pub fn code(&self) -> crate::residue_values::ResidueCode {
        self.view().code().into()
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.view().id().value()
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.view().name().as_str().into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> ResidueKind {
        self.view().kind().into()
    }
    #[wasm_bindgen(getter,js_name=oneLetterCode)]
    pub fn one_letter_code(&self) -> char {
        self.view().one_letter_code()
    }
    #[wasm_bindgen(getter,js_name=fastaCode)]
    pub fn fasta_code(&self) -> char {
        self.view().fasta_code()
    }
    #[wasm_bindgen(getter,js_name=isStandard)]
    pub fn is_standard(&self) -> bool {
        self.view().is_standard()
    }
    #[wasm_bindgen(getter)]
    pub fn chain(&self) -> ProteinChainRef {
        ProteinChainRef(ChainOwner(
            self.0.0.clone(),
            self.view().chain().id().value() as usize,
        ))
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> Result<ProteinAtoms, JsValue> {
        <Vec<AtomOwner> as Field>::get(NativeField::owned(
            self.view()
                .atoms()
                .map(|v| AtomOwner(self.0.0.clone(), v.id().value() as usize))
                .collect(),
        ))
    }
}
#[wasm_bindgen]
impl ProteinAtomRef {
    #[wasm_bindgen(getter)]
    pub fn row(&self) -> crate::bio_rows::BioAtomRow {
        crate::bio_rows::BioAtomRow(self.view().row().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> u32 {
        self.view().id().value()
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> crate::bio_source_values::AtomName {
        crate::bio_source_values::AtomName(self.view().name())
    }
    #[wasm_bindgen(getter)]
    pub fn element(&self) -> crate::atoms::Element {
        crate::atoms::Element(self.view().element())
    }
    #[wasm_bindgen(getter)]
    pub fn altloc(&self) -> Option<crate::bio_source_values::AltLocLabel> {
        self.view()
            .altloc()
            .map(crate::bio_source_values::AltLocLabel)
    }
    #[wasm_bindgen(getter)]
    pub fn residue(&self) -> ProteinResidueRef {
        ProteinResidueRef(ResidueOwner(
            self.0.0.clone(),
            self.view().residue().id().value() as usize,
        ))
    }
    #[wasm_bindgen(getter)]
    pub fn position(&self) -> Result<Ts<Point3>, JsValue> {
        Point3(self.view().position())
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
}
