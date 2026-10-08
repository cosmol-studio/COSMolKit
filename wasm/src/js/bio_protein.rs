//! Protein projection, COW operations and source-ordered borrowed row snapshots.
use crate::{
    alignment_values::{set, source_error},
    bio_hierarchy::{BioAtomRow, BioChainRow, BioResidueRow, ChainSourceIds},
    bio_readers::{BioPdbReadParams, BioReadParams, BioStructure, format_value},
    bio_residue::ResidueInfo,
    bio_selection::{BioSelection, point},
    host_values::usize_value,
};
use cosmolkit_wasm::rust as ck;
use std::{error::Error as RustError, path::Path, rc::Rc};
use wasm_bindgen::prelude::*;
#[wasm_bindgen(
    inline_js = "export function attachProteinIterator(v) { Object.defineProperty(v, Symbol.iterator, {value: function() { return this; }}); return v; }"
)]
extern "C" {
    #[wasm_bindgen(js_name=attachProteinIterator)]
    fn attach_iterator(value: JsValue) -> JsValue;
}
#[wasm_bindgen(typescript_custom_section)]
const ITERATORS: &str = r#"
export interface ProteinChainIter extends IterableIterator<ProteinChainRef> {}
export interface ProteinResidueIter extends IterableIterator<ProteinResidueRef> {}
export interface ProteinAtomIter extends IterableIterator<ProteinAtomRef> {}
"#;
#[wasm_bindgen]
pub struct Protein {
    pub(crate) inner: Rc<ck::Protein>,
}
#[wasm_bindgen]
impl BioStructure {
    pub fn protein(&self) -> Result<Protein, JsValue> {
        // COSMolKit❗✔️: self.inner.protein()
        self.inner
            .protein()
            .map(|inner| Protein {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_selection::protein_projection_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
pub struct ProteinSelectionSummary {
    inner: ck::ProteinSelectionSummary,
}
#[wasm_bindgen]
impl ProteinSelectionSummary {
    #[wasm_bindgen(getter)]
    pub fn chains(&self) -> usize {
        self.inner.chains
    }
    #[wasm_bindgen(getter)]
    pub fn residues(&self) -> usize {
        self.inner.residues
    }
    #[wasm_bindgen(getter)]
    pub fn atoms(&self) -> usize {
        self.inner.atoms
    }
}
#[wasm_bindgen]
impl Protein {
    #[wasm_bindgen(js_name=fromPdb)]
    pub fn from_pdb(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::Protein::from_pdb(text)
        ck::Protein::from_pdb(text)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromPdbWithParams)]
    pub fn from_pdb_with_params(text: &str, params: &BioPdbReadParams) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::Protein::from_pdb_with_params(text,&params.inner)
        ck::Protein::from_pdb_with_params(text, &params.inner)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromMmcif)]
    pub fn from_mmcif(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::Protein::from_mmcif(text)
        ck::Protein::from_mmcif(text)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromText)]
    pub fn from_text(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::Protein::from_text(text)
        ck::Protein::from_text(text)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromTextWithParams)]
    pub fn from_text_with_params(text: &str, params: &BioReadParams) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::Protein::from_text_with_params(text,&params.inner)
        ck::Protein::from_text_with_params(text, &params.inner)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=read)]
    pub fn read(path: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::Protein::read(Path::new(path))
        ck::Protein::read(Path::new(path))
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=readWithFormat)]
    pub fn read_with_format(
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "BioCoordinateFormat")] format: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::Protein::read_with_format(Path::new(path),format_value(&format)?)
        ck::Protein::read_with_format(Path::new(path), format_value(&format)?)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=inputFormat,unchecked_return_type="BioCoordinateFormat")]
    pub fn input_format(&self) -> u32 {
        self.inner.input_format() as u32
    }
    #[wasm_bindgen(js_name=asBioStructure)]
    pub fn as_bio_structure(&self) -> BioStructure {
        // COSMolKit❗✔️: inner: Arc::new(self.inner.as_bio_structure().clone()),
        // Canonical block-sharing clone: never materialize all hierarchy or coordinate data.
        BioStructure {
            inner: Rc::new(self.inner.as_bio_structure().clone()),
        }
    }
    #[wasm_bindgen(js_name=intoBioStructure)]
    pub fn into_bio_structure(&self) -> BioStructure {
        // COSMolKit❗✔️: inner: Arc::new(self.inner.as_ref().clone().into_bio_structure()),
        BioStructure {
            inner: Rc::new(self.inner.as_ref().clone().into_bio_structure()),
        }
    }
    #[wasm_bindgen(js_name=selectionSummary)]
    pub fn selection_summary(&self) -> ProteinSelectionSummary {
        ProteinSelectionSummary {
            inner: self.inner.selection_summary(),
        }
    }
    pub fn len(&self) -> usize {
        self.inner.num_chains()
    }
    #[wasm_bindgen(unchecked_return_type = "ProteinChainRef | null")]
    pub fn chain(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<JsValue, JsValue> {
        let index = usize_value(&index, "index")?;
        Ok(self.inner.chain(index).map_or(JsValue::NULL, |v| {
            ProteinChainRef {
                inner: self.inner.clone(),
                index: v.id().index(),
            }
            .into()
        }))
    }
    #[wasm_bindgen(js_name=selectedAtomIds)]
    pub fn selected_atom_ids(&self, selection: &BioSelection) -> Result<Vec<u32>, JsValue> {
        self.inner
            .selected_atom_ids(&selection.inner)
            .map(|v| v.into_iter().map(|id| id.value()).collect())
            .map_err(|e| crate::bio_selection::match_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withSelection)]
    pub fn with_selection(&self, selection: &BioSelection) -> Result<Self, JsValue> {
        self.inner
            .with_selection(&selection.inner)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_selection::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=retainSelection_)]
    pub fn retain_selection_(&mut self, selection: &BioSelection) -> Result<(), JsValue> {
        // COSMolKit❗✔️: Arc::make_mut(&mut self.inner)
        // Source's generated operation owns all validation, COW and atomic commit.
        Rc::make_mut(&mut self.inner)
            .retain_selection_(&selection.inner)
            .map_err(|e| crate::bio_selection::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withTranslatedCoordinates)]
    pub fn with_translated_coordinates(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array")] offset: JsValue,
    ) -> Result<Self, JsValue> {
        self.inner
            .with_translated_coordinates(point(&offset)?)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_selection::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=translate_)]
    pub fn translate_(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array")] offset: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .translate_(offset)
        let offset = point(&offset)?;
        Rc::make_mut(&mut self.inner)
            .translate_(offset)
            .map_err(|e| crate::bio_selection::operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=numModels)]
    pub fn num_models(&self) -> usize {
        self.inner.num_models()
    }
    #[wasm_bindgen(js_name=numChains)]
    pub fn num_chains(&self) -> usize {
        self.inner.num_chains()
    }
    #[wasm_bindgen(js_name=numResidues)]
    pub fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }
    #[wasm_bindgen(js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    #[wasm_bindgen(unchecked_return_type = "ProteinChainIter")]
    pub fn chains(&self) -> JsValue {
        attach_iterator(
            ProteinChainIter {
                inner: self.inner.clone(),
                scope: Scope::Whole,
                cursor: 0,
            }
            .into(),
        )
    }
    #[wasm_bindgen(unchecked_return_type = "ProteinResidueIter")]
    pub fn residues(&self) -> JsValue {
        attach_iterator(
            ProteinResidueIter {
                inner: self.inner.clone(),
                scope: Scope::Whole,
                cursor: 0,
            }
            .into(),
        )
    }
    #[wasm_bindgen(unchecked_return_type = "ProteinAtomIter")]
    pub fn atoms(&self) -> JsValue {
        attach_iterator(
            ProteinAtomIter {
                inner: self.inner.clone(),
                scope: Scope::Whole,
                cursor: 0,
            }
            .into(),
        )
    }
}
#[derive(Clone, Copy)]
enum Scope {
    Whole,
    Chain(usize),
    Residue(usize),
}
#[wasm_bindgen]
pub struct ProteinChainRef {
    inner: Rc<ck::Protein>,
    index: usize,
}
impl ProteinChainRef {
    fn view(&self) -> ck::ProteinChainRef<'_> {
        self.inner
            .chain(self.index)
            .expect("validated protein chain view")
    }
}
#[wasm_bindgen]
impl ProteinChainRef {
    pub fn id(&self) -> usize {
        self.index
    }
    pub fn row(&self) -> BioChainRow {
        // COSMolKit❗✔️: inner: Arc::new(self.inner.as_bio_structure().clone()),
        BioChainRow {
            inner: Rc::new(self.inner.as_bio_structure().clone()),
            index: self.view().id().index(),
        }
    }
    #[wasm_bindgen(unchecked_return_type = "ChainKind")]
    pub fn kind(&self) -> u32 {
        self.view().kind() as u32
    }
    pub fn source(&self) -> ChainSourceIds {
        ChainSourceIds {
            inner: self.view().source().clone(),
        }
    }
    pub fn len(&self) -> usize {
        self.view().residues().len()
    }
    #[wasm_bindgen(unchecked_return_type = "ProteinResidueIter")]
    pub fn residues(&self) -> JsValue {
        attach_iterator(
            ProteinResidueIter {
                inner: self.inner.clone(),
                scope: Scope::Chain(self.index),
                cursor: 0,
            }
            .into(),
        )
    }
    #[wasm_bindgen(unchecked_return_type = "ProteinAtomIter")]
    pub fn atoms(&self) -> JsValue {
        attach_iterator(
            ProteinAtomIter {
                inner: self.inner.clone(),
                scope: Scope::Chain(self.index),
                cursor: 0,
            }
            .into(),
        )
    }
}
#[wasm_bindgen]
pub struct ProteinChainIter {
    inner: Rc<ck::Protein>,
    scope: Scope,
    cursor: usize,
}
impl ProteinChainIter {
    fn view(&self) -> ck::ProteinChainIter<'_> {
        match self.scope {
            Scope::Whole => self.inner.chains(),
            _ => unreachable!("chain iterator only has whole-structure scope"),
        }
    }
}
#[wasm_bindgen]
impl ProteinChainIter {
    pub fn len(&self) -> usize {
        self.view().len().saturating_sub(self.cursor)
    }
    #[wasm_bindgen(unchecked_return_type = "IteratorResult<ProteinChainRef, undefined>")]
    pub fn next(&mut self) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.cursor.nth(n).map(|index| Self::Item {
        // Delegates O(1) skip/next to the public iterator; no row collection or chemistry logic.
        let next = self.view().nth(self.cursor).map(|v| v.id().index());
        let out: JsValue = js_sys::Object::new().into();
        set(&out, "done", next.is_none().into())?;
        let value = if let Some(index) = next {
            self.cursor += 1;
            ProteinChainRef {
                inner: self.inner.clone(),
                index,
            }
            .into()
        } else {
            JsValue::UNDEFINED
        };
        set(&out, "value", value)?;
        Ok(out)
    }
}
#[wasm_bindgen]
pub struct ProteinResidueRef {
    inner: Rc<ck::Protein>,
    index: usize,
}
impl ProteinResidueRef {
    fn view(&self) -> ck::ProteinResidueRef<'_> {
        self.inner
            .residues()
            .nth(self.index)
            .expect("validated protein residue view")
    }
}
#[wasm_bindgen]
impl ProteinResidueRef {
    pub fn id(&self) -> usize {
        self.index
    }
    pub fn row(&self) -> BioResidueRow {
        // COSMolKit❗✔️: inner: Arc::new(self.inner.as_bio_structure().clone()),
        BioResidueRow {
            inner: Rc::new(self.inner.as_bio_structure().clone()),
            index: self.view().id().index(),
        }
    }
    pub fn chain(&self) -> ProteinChainRef {
        ProteinChainRef {
            inner: self.inner.clone(),
            index: self.view().chain().id().index(),
        }
    }
    pub fn name(&self) -> String {
        self.view().name().as_str().to_owned()
    }
    #[wasm_bindgen(unchecked_return_type = "ResidueKind")]
    pub fn kind(&self) -> u32 {
        self.view().kind() as u32
    }
    pub fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: self.view().info(),
        }
    }
    #[wasm_bindgen(unchecked_return_type = "ResidueCode")]
    pub fn code(&self) -> u16 {
        self.view().code().as_u16()
    }
    #[wasm_bindgen(js_name=oneLetterCode)]
    pub fn one_letter_code(&self) -> String {
        self.view().one_letter_code().to_string()
    }
    #[wasm_bindgen(js_name=fastaCode)]
    pub fn fasta_code(&self) -> String {
        self.view().fasta_code().to_string()
    }
    #[wasm_bindgen(js_name=canonicalOneLetterCode,unchecked_return_type="string | null")]
    pub fn canonical_one_letter_code(&self) -> JsValue {
        self.view()
            .info()
            .canonical_one_letter_code()
            .map_or(JsValue::NULL, |v| v.to_string().into())
    }
    #[wasm_bindgen(js_name=parentStandardCode,unchecked_return_type="ResidueCode | null")]
    pub fn parent_standard_code(&self) -> JsValue {
        self.view()
            .info()
            .parent_standard_code()
            .map_or(JsValue::NULL, |v| v.as_u16().into())
    }
    #[wasm_bindgen(js_name=isModifiedAminoAcid)]
    pub fn is_modified_amino_acid(&self) -> bool {
        self.view().info().is_modified_amino_acid()
    }
    #[wasm_bindgen(js_name=isStandard)]
    pub fn is_standard(&self) -> bool {
        self.view().is_standard()
    }
    pub fn len(&self) -> usize {
        self.view().atoms().len()
    }
    #[wasm_bindgen(unchecked_return_type = "ProteinAtomIter")]
    pub fn atoms(&self) -> JsValue {
        attach_iterator(
            ProteinAtomIter {
                inner: self.inner.clone(),
                scope: Scope::Residue(self.index),
                cursor: 0,
            }
            .into(),
        )
    }
}
#[wasm_bindgen]
pub struct ProteinResidueIter {
    inner: Rc<ck::Protein>,
    scope: Scope,
    cursor: usize,
}
impl ProteinResidueIter {
    fn view(&self) -> ck::ProteinResidueIter<'_> {
        match self.scope {
            Scope::Whole => self.inner.residues(),
            Scope::Chain(index) => self.inner.chain(index).expect("validated chain").residues(),
            Scope::Residue(_) => unreachable!("residue iterator cannot contain a residue scope"),
        }
    }
}
#[wasm_bindgen]
impl ProteinResidueIter {
    pub fn len(&self) -> usize {
        self.view().len().saturating_sub(self.cursor)
    }
    #[wasm_bindgen(unchecked_return_type = "IteratorResult<ProteinResidueRef, undefined>")]
    pub fn next(&mut self) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.cursor.nth(n).map(|index| Self::Item {
        // Delegates O(1) skip/next to the public iterator; no row collection or chemistry logic.
        let next = self.view().nth(self.cursor).map(|v| v.id().index());
        let out: JsValue = js_sys::Object::new().into();
        set(&out, "done", next.is_none().into())?;
        let value = if let Some(index) = next {
            self.cursor += 1;
            ProteinResidueRef {
                inner: self.inner.clone(),
                index,
            }
            .into()
        } else {
            JsValue::UNDEFINED
        };
        set(&out, "value", value)?;
        Ok(out)
    }
}
#[wasm_bindgen]
pub struct ProteinAtomRef {
    inner: Rc<ck::Protein>,
    index: usize,
}
impl ProteinAtomRef {
    fn view(&self) -> ck::ProteinAtomRef<'_> {
        self.inner
            .atoms()
            .nth(self.index)
            .expect("validated protein atom view")
    }
}
#[wasm_bindgen]
impl ProteinAtomRef {
    pub fn id(&self) -> usize {
        self.index
    }
    pub fn row(&self) -> BioAtomRow {
        // COSMolKit❗✔️: inner: Arc::new(self.inner.as_bio_structure().clone()),
        BioAtomRow {
            inner: Rc::new(self.inner.as_bio_structure().clone()),
            index: self.view().id().index(),
        }
    }
    pub fn residue(&self) -> ProteinResidueRef {
        ProteinResidueRef {
            inner: self.inner.clone(),
            index: self.view().residue().id().index(),
        }
    }
    pub fn name(&self) -> String {
        self.view().name().as_str().trim().to_owned()
    }
    #[wasm_bindgen(unchecked_return_type = "Element")]
    pub fn element(&self) -> Result<JsValue, JsValue> {
        crate::element_boundary::Element::from_atomic_number(
            self.view().element().atomic_number() as f64
        )
    }
    #[wasm_bindgen(js_name=elementSymbol)]
    pub fn element_symbol(&self) -> String {
        self.view().element().symbol().to_owned()
    }
    #[wasm_bindgen(js_name=atomicNumber)]
    pub fn atomic_number(&self) -> u8 {
        self.view().element().atomic_number()
    }
    #[wasm_bindgen(unchecked_return_type = "number[]")]
    pub fn position(&self) -> JsValue {
        self.view()
            .position()
            .into_iter()
            .map(JsValue::from)
            .collect::<js_sys::Array>()
            .into()
    }
    #[wasm_bindgen(unchecked_return_type = "string | null")]
    pub fn altloc(&self) -> JsValue {
        self.view()
            .altloc()
            .map_or(JsValue::NULL, |v| char::from(v.value()).to_string().into())
    }
}
#[wasm_bindgen]
pub struct ProteinAtomIter {
    inner: Rc<ck::Protein>,
    scope: Scope,
    cursor: usize,
}
impl ProteinAtomIter {
    fn view(&self) -> ck::ProteinAtomIter<'_> {
        match self.scope {
            Scope::Whole => self.inner.atoms(),
            Scope::Chain(index) => self.inner.chain(index).expect("validated chain").atoms(),
            Scope::Residue(index) => self
                .inner
                .residues()
                .nth(index)
                .expect("validated residue")
                .atoms(),
        }
    }
}
#[wasm_bindgen]
impl ProteinAtomIter {
    pub fn len(&self) -> usize {
        self.view().len().saturating_sub(self.cursor)
    }
    #[wasm_bindgen(unchecked_return_type = "IteratorResult<ProteinAtomRef, undefined>")]
    pub fn next(&mut self) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.cursor.nth(n).map(|index| Self::Item {
        // Delegates O(1) skip/next to the public iterator; no row collection or chemistry logic.
        let next = self.view().nth(self.cursor).map(|v| v.id().index());
        let out: JsValue = js_sys::Object::new().into();
        set(&out, "done", next.is_none().into())?;
        let value = if let Some(index) = next {
            self.cursor += 1;
            ProteinAtomRef {
                inner: self.inner.clone(),
                index,
            }
            .into()
        } else {
            JsValue::UNDEFINED
        };
        set(&out, "value", value)?;
        Ok(out)
    }
}
#[wasm_bindgen]
pub struct ProteinReadError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ProteinReadError {
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
pub(crate) fn read_error(source: &ck::ProteinReadError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: let cause = match source {
    // Canonical error.source selects the same four typed source variants.
    let kind = match source {
        ck::ProteinReadError::Structure(_) => "Structure",
        ck::ProteinReadError::Pdb(_) => "Pdb",
        ck::ProteinReadError::Mmcif(_) => "Mmcif",
        ck::ProteinReadError::Projection(_) => "Projection",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("ProteinReadError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    set(
        &e,
        "detail",
        ProteinReadError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
