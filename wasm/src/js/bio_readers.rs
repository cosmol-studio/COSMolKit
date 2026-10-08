//! Read-only structural constructors and complete read options through cosmolkit.
use crate::host_values::{bool_value, i32_value, type_error, u32_value};
use cosmolkit_wasm::rust as ck;
use std::{path::Path, rc::Rc};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum BioCoordinateFormat {
    Unknown = 0,
    Detect = 1,
    Pdb = 2,
    Mmcif = 3,
    Mmjson = 4,
    ChemComp = 5,
}
pub(crate) fn format_value(value: &JsValue) -> Result<ck::BioCoordinateFormat, JsValue> {
    use ck::BioCoordinateFormat as F;
    match u32_value(value, "format")? {
        0 => Ok(F::Unknown),
        1 => Ok(F::Detect),
        2 => Ok(F::Pdb),
        3 => Ok(F::Mmcif),
        4 => Ok(F::Mmjson),
        5 => Ok(F::ChemComp),
        _ => Err(js_sys::RangeError::new("invalid BioCoordinateFormat code").into()),
    }
}
#[wasm_bindgen]
pub struct BioReadParams {
    pub(crate) inner: ck::BioReadParams,
}
#[wasm_bindgen]
impl BioReadParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "BioCoordinateFormat")] format: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string")] source_name: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::BioReadParams::default();
        if !format.is_undefined() {
            inner.format = format_value(&format)?;
        }
        if !source_name.is_undefined() {
            inner.source_name = source_name
                .as_string()
                .ok_or_else(|| type_error("sourceName"))?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, unchecked_return_type = "BioCoordinateFormat")]
    pub fn format(&self) -> u32 {
        self.inner.format as u32
    }
    #[wasm_bindgen(getter,js_name=sourceName)]
    pub fn source_name(&self) -> String {
        self.inner.source_name.clone()
    }
}
#[wasm_bindgen]
pub struct BioPdbReadParams {
    pub(crate) inner: ck::BioPdbReadParams,
}
#[wasm_bindgen]
impl BioPdbReadParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_line_length: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] check_non_ascii: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ignore_ter: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] split_chain_on_ter: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] skip_remarks: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::BioPdbReadParams::default();
        if !max_line_length.is_undefined() {
            inner.max_line_length = i32_value(&max_line_length, "maxLineLength")?;
        }
        if !check_non_ascii.is_undefined() {
            inner.check_non_ascii = bool_value(&check_non_ascii, "checkNonAscii")?;
        }
        if !ignore_ter.is_undefined() {
            inner.ignore_ter = bool_value(&ignore_ter, "ignoreTer")?;
        }
        if !split_chain_on_ter.is_undefined() {
            inner.split_chain_on_ter = bool_value(&split_chain_on_ter, "splitChainOnTer")?;
        }
        if !skip_remarks.is_undefined() {
            inner.skip_remarks = bool_value(&skip_remarks, "skipRemarks")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=maxLineLength)]
    pub fn max_line_length(&self) -> i32 {
        self.inner.max_line_length
    }
    #[wasm_bindgen(getter,js_name=checkNonAscii)]
    pub fn check_non_ascii(&self) -> bool {
        self.inner.check_non_ascii
    }
    #[wasm_bindgen(getter,js_name=ignoreTer)]
    pub fn ignore_ter(&self) -> bool {
        self.inner.ignore_ter
    }
    #[wasm_bindgen(getter,js_name=splitChainOnTer)]
    pub fn split_chain_on_ter(&self) -> bool {
        self.inner.split_chain_on_ter
    }
    #[wasm_bindgen(getter,js_name=skipRemarks)]
    pub fn skip_remarks(&self) -> bool {
        self.inner.skip_remarks
    }
}
#[wasm_bindgen]
pub struct BioStructure {
    pub(crate) inner: Rc<ck::BioStructure>,
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=fromPdb)]
    pub fn from_pdb(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_bio_binding.rs::BioStructure constructors:
        //             .map(|x| Self { inner: Arc::new(x) })
        // Shared immutable hierarchy snapshots; one public owner call, no parser logic.
        ck::BioStructure::from_pdb(text)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_read_errors::pdb_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromPdbWithParams)]
    pub fn from_pdb_with_params(text: &str, params: &BioPdbReadParams) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_bio_binding.rs::BioStructure constructors:
        //             .map(|x| Self { inner: Arc::new(x) })
        // Shared immutable hierarchy snapshots; one public owner call, no parser logic.
        ck::BioStructure::from_pdb_with_params(text, &params.inner)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_read_errors::pdb_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromMmcif)]
    pub fn from_mmcif(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_bio_binding.rs::BioStructure constructors:
        //             .map(|x| Self { inner: Arc::new(x) })
        // Shared immutable hierarchy snapshots; one public owner call, no parser logic.
        ck::BioStructure::from_mmcif(text)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_read_errors::mmcif_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromText)]
    pub fn from_text(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_bio_binding.rs::BioStructure constructors:
        //             .map(|x| Self { inner: Arc::new(x) })
        // Shared immutable hierarchy snapshots; one public owner call, no parser logic.
        ck::BioStructure::from_text(text)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_read_errors::read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromTextWithParams)]
    pub fn from_text_with_params(text: &str, params: &BioReadParams) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_bio_binding.rs::BioStructure constructors:
        //             .map(|x| Self { inner: Arc::new(x) })
        // Shared immutable hierarchy snapshots; one public owner call, no parser logic.
        ck::BioStructure::from_text_with_params(text, &params.inner)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_read_errors::read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=read)]
    pub fn read(path: &str) -> Result<Self, JsValue> {
        ck::BioStructure::read(Path::new(path))
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_read_errors::read_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=readWithFormat)]
    pub fn read_with_format(
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "BioCoordinateFormat")] format: JsValue,
    ) -> Result<Self, JsValue> {
        ck::BioStructure::read_with_format(Path::new(path), format_value(&format)?)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| crate::bio_read_errors::read_error(&e).unwrap_or_else(|e| e))
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
    #[wasm_bindgen(js_name=numEntities)]
    pub fn num_entities(&self) -> usize {
        self.inner.num_entities()
    }
    pub fn name(&self) -> String {
        self.inner.name().into()
    }
    #[wasm_bindgen(js_name=inputFormat,unchecked_return_type="BioCoordinateFormat")]
    pub fn input_format(&self) -> u32 {
        self.inner.input_format() as u32
    }
}

use crate::bio_hierarchy::*;
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(unchecked_return_type = "BioModelRow[]")]
    pub fn models(&self) -> js_sys::Array {
        // COSMolKit❗✔️: (0..self.inner.num_models()).map(|index| BioModelRow { inner: self.inner.clone(), index, }).collect()
        (0..self.inner.num_models())
            .map(|index| -> JsValue {
                BioModelRow {
                    inner: Rc::clone(&self.inner),
                    index,
                }
                .into()
            })
            .collect()
    }
    #[wasm_bindgen(unchecked_return_type = "BioChainRow[]")]
    pub fn chains(&self) -> js_sys::Array {
        // COSMolKit❗✔️: (0..self.inner.num_chains()).map(|index| BioChainRow { inner: self.inner.clone(), index, }).collect()
        (0..self.inner.num_chains())
            .map(|index| -> JsValue {
                BioChainRow {
                    inner: Rc::clone(&self.inner),
                    index,
                }
                .into()
            })
            .collect()
    }
    #[wasm_bindgen(unchecked_return_type = "BioResidueRow[]")]
    pub fn residues(&self) -> js_sys::Array {
        // COSMolKit❗✔️: (0..self.inner.num_residues()).map(|index| BioResidueRow { inner: self.inner.clone(), index, }).collect()
        (0..self.inner.num_residues())
            .map(|index| -> JsValue {
                BioResidueRow {
                    inner: Rc::clone(&self.inner),
                    index,
                }
                .into()
            })
            .collect()
    }
    #[wasm_bindgen(unchecked_return_type = "BioAtomRow[]")]
    pub fn atoms(&self) -> js_sys::Array {
        // COSMolKit❗✔️: (0..self.inner.num_atoms()).map(|index| BioAtomRow { inner: self.inner.clone(), index, }).collect()
        (0..self.inner.num_atoms())
            .map(|index| -> JsValue {
                BioAtomRow {
                    inner: Rc::clone(&self.inner),
                    index,
                }
                .into()
            })
            .collect()
    }
    #[wasm_bindgen(unchecked_return_type = "BioEntityRow[]")]
    pub fn entities(&self) -> js_sys::Array {
        // COSMolKit❗✔️: (0..self.inner.num_entities()).map(|index| BioEntityRow { inner: self.inner.clone(), index, }).collect()
        (0..self.inner.num_entities())
            .map(|index| -> JsValue {
                BioEntityRow {
                    inner: Rc::clone(&self.inner),
                    index,
                }
                .into()
            })
            .collect()
    }
    pub fn coordinates(&self) -> BioCoordinateBlock {
        BioCoordinateBlock::from_owner(Rc::clone(&self.inner))
    }
    #[wasm_bindgen(js_name=atomPosition,unchecked_return_type="number[] | null")]
    pub fn atom_position(&self, atom: JsValue) -> Result<JsValue, JsValue> {
        Ok(self
            .inner
            .atom_position(ck::BioAtomId::new(u32_value(&atom, "atom")?))
            .map_or(JsValue::NULL, |p| {
                p.into_iter()
                    .map(JsValue::from)
                    .collect::<js_sys::Array>()
                    .into()
            }))
    }
    #[wasm_bindgen(js_name=residueAtoms,unchecked_return_type="BioAtomRow[] | null")]
    pub fn residue_atoms(&self, residue: JsValue) -> Result<JsValue, JsValue> {
        let id = u32_value(&residue, "residue")?;
        Ok(self
            .inner
            .residue_atoms(ck::BioResidueId::new(id))
            .map_or(JsValue::NULL, |rows| {
                let start = self.inner.residues()[id as usize].atom_span().start() as usize;
                (start..start + rows.len())
                    .map(|index| -> JsValue {
                        BioAtomRow {
                            inner: Rc::clone(&self.inner),
                            index,
                        }
                        .into()
                    })
                    .collect::<js_sys::Array>()
                    .into()
            }))
    }
    #[wasm_bindgen(js_name=findEntity,unchecked_return_type="[number,BioEntityRow] | null")]
    pub fn find_entity(&self, source_id: &str) -> JsValue {
        self.inner
            .find_entity(source_id)
            .map_or(JsValue::NULL, |(id, _)| {
                js_sys::Array::of2(
                    &JsValue::from(id.value()),
                    &BioEntityRow {
                        inner: Rc::clone(&self.inner),
                        index: id.index(),
                    }
                    .into(),
                )
                .into()
            })
    }
    #[wasm_bindgen(js_name=findEntityOfSubchain,unchecked_return_type="[number,BioEntityRow] | null")]
    pub fn find_entity_of_subchain(&self, subchain: &str) -> JsValue {
        self.inner
            .find_entity_of_subchain(subchain)
            .map_or(JsValue::NULL, |(id, _)| {
                js_sys::Array::of2(
                    &JsValue::from(id.value()),
                    &BioEntityRow {
                        inner: Rc::clone(&self.inner),
                        index: id.index(),
                    }
                    .into(),
                )
                .into()
            })
    }
    #[wasm_bindgen(js_name=findAtom,unchecked_return_type="[number,BioAtomRow] | null")]
    pub fn find_atom(
        &self,
        residue_id: JsValue,
        name: &AtomName,
        request: &AltLocRequest,
        #[wasm_bindgen(unchecked_param_type = "Element | null")] element: JsValue,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.find_atom(ck::BioResidueId::new(residue_id),name.inner,request.inner,element.map(|value|value.inner))
        let element = optional_element(&element)?;
        Ok(self
            .inner
            .find_atom(
                ck::BioResidueId::new(u32_value(&residue_id, "residueId")?),
                name.inner,
                request.inner,
                element,
            )
            .map_or(JsValue::NULL, |(id, _)| {
                js_sys::Array::of2(
                    &JsValue::from(id.value()),
                    &BioAtomRow {
                        inner: Rc::clone(&self.inner),
                        index: id.index(),
                    }
                    .into(),
                )
                .into()
            }))
    }
    #[wasm_bindgen(js_name=atomByAltloc,unchecked_return_type="[number,BioAtomRow]")]
    pub fn atom_by_altloc(
        &self,
        residue_id: JsValue,
        name: &AtomName,
        #[wasm_bindgen(unchecked_param_type = "AltLocLabel | null")] altloc: JsValue,
    ) -> Result<JsValue, JsValue> {
        self.inner
            .atom_by_altloc(
                ck::BioResidueId::new(u32_value(&residue_id, "residueId")?),
                name.inner,
                optional_label(&altloc)?,
            )
            .map(|(id, _)| {
                js_sys::Array::of2(
                    &JsValue::from(id.value()),
                    &BioAtomRow {
                        inner: Rc::clone(&self.inner),
                        index: id.index(),
                    }
                    .into(),
                )
                .into()
            })
            .map_err(|e| crate::bio_hierarchy_errors::structure_error(&e).unwrap_or_else(|e| e))
    }
}
