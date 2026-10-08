//! Readonly detached BIO metadata; every getter projects the public owner's value.
use crate::{
    bio_hierarchy::{BioEntityRow, BioResidueRow, PdbSeqId},
    bio_readers::BioStructure,
};
use cosmolkit_wasm::rust as ck;
use std::rc::Rc;
use wasm_bindgen::prelude::*;
fn nullable<T: Into<JsValue>>(v: Option<T>) -> JsValue {
    v.map_or(JsValue::NULL, Into::into)
}
fn row(v: &[f64]) -> js_sys::Array {
    v.iter().copied().map(JsValue::from).collect()
}
fn matrix(v: &[[f64; 3]; 3]) -> js_sys::Array {
    v.iter().map(|v| -> JsValue { row(v).into() }).collect()
}

#[wasm_bindgen]
pub enum BioConnectionKind {
    Covale = 0,
    Disulf = 1,
    Hydrog = 2,
    MetalC = 3,
    Unknown = 4,
}
#[wasm_bindgen]
pub enum BioAsu {
    Same = 0,
    Different = 1,
    Any = 2,
}
#[wasm_bindgen]
pub enum BioHelixClass {
    UnknownHelix = 0,
    RAlpha = 1,
    ROmega = 2,
    RPi = 3,
    RGamma = 4,
    R310 = 5,
    LAlpha = 6,
    LOmega = 7,
    LGamma = 8,
    Helix27 = 9,
    HelixPolyProlineNone = 10,
}
#[wasm_bindgen]
pub enum BioSoftwareClassification {
    DataCollection = 0,
    DataExtraction = 1,
    DataProcessing = 2,
    DataReduction = 3,
    DataScaling = 4,
    ModelBuilding = 5,
    Phasing = 6,
    Refinement = 7,
    Unspecified = 8,
}
#[wasm_bindgen]
pub enum BioAssemblySpecialKind {
    NotApplicable = 0,
    CompleteIcosahedral = 1,
    RepresentativeHelical = 2,
    CompletePoint = 3,
}
#[wasm_bindgen]
pub struct BioConnection {
    pub(crate) inner: ck::BioConnection,
}
#[wasm_bindgen]
impl BioConnection {
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        // COSMolKit❗✔️: &self.inner.name
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.name.to_owned()
    }
    #[wasm_bindgen(getter,js_name=linkId)]
    pub fn link_id(&self) -> String {
        // COSMolKit❗✔️: &self.inner.link_id
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.link_id.to_owned()
    }
    #[wasm_bindgen(getter,js_name=kind,unchecked_return_type="BioConnectionKind")]
    pub fn kind(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.kind.into()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.kind as u32
    }
    #[wasm_bindgen(getter,js_name=asu,unchecked_return_type="BioAsu")]
    pub fn asu(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.asu.into()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.asu as u32
    }
    #[wasm_bindgen(getter,js_name=partner1)]
    pub fn partner1(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.partner1.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.partner1.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=partner2)]
    pub fn partner2(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.partner2.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.partner2.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=reportedDistance)]
    pub fn reported_distance(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.reported_distance
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.reported_distance
    }
    #[wasm_bindgen(getter,js_name=reportedSym,unchecked_return_type="number[]")]
    pub fn reported_sym(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.reported_sym
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .reported_sym
            .iter()
            .map(|v| JsValue::from(*v))
            .collect()
    }
}
#[wasm_bindgen]
pub struct BioCisPep {
    pub(crate) inner: ck::BioCisPep,
}
#[wasm_bindgen]
impl BioCisPep {
    #[wasm_bindgen(getter,js_name=partnerC)]
    pub fn partner_c(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.partner_c.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.partner_c.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=partnerN)]
    pub fn partner_n(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.partner_n.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.partner_n.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=modelNum)]
    pub fn model_num(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.model_num
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.model_num
    }
    #[wasm_bindgen(getter,js_name=onlyAltloc)]
    pub fn only_altloc(&self) -> u8 {
        // COSMolKit❗✔️: self.inner.only_altloc
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.only_altloc
    }
    #[wasm_bindgen(getter,js_name=reportedAngle)]
    pub fn reported_angle(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.reported_angle
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.reported_angle
    }
}
#[wasm_bindgen]
pub struct BioModRes {
    pub(crate) inner: ck::BioModRes,
}
#[wasm_bindgen]
impl BioModRes {
    #[wasm_bindgen(getter,js_name=chainName)]
    pub fn chain_name(&self) -> String {
        // COSMolKit❗✔️: self.inner.chain_name.as_str().to_owned()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.chain_name.as_str().to_owned()
    }
    #[wasm_bindgen(getter,js_name=resId)]
    pub fn res_id(&self) -> ResidueAddress {
        // COSMolKit❗✔️: inner: self.inner.res_id.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        ResidueAddress {
            inner: self.inner.res_id.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=parentCompId)]
    pub fn parent_comp_id(&self) -> String {
        // COSMolKit❗✔️: &self.inner.parent_comp_id
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.parent_comp_id.to_owned()
    }
    #[wasm_bindgen(getter,js_name=modId)]
    pub fn mod_id(&self) -> String {
        // COSMolKit❗✔️: &self.inner.mod_id
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.mod_id.to_owned()
    }
    #[wasm_bindgen(getter,js_name=details)]
    pub fn details(&self) -> String {
        // COSMolKit❗✔️: &self.inner.details
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.details.to_owned()
    }
}
#[wasm_bindgen]
pub struct BioHelix {
    pub(crate) inner: ck::BioHelix,
}
#[wasm_bindgen]
impl BioHelix {
    #[wasm_bindgen(getter,js_name=start)]
    pub fn start(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.start.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.start.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=end)]
    pub fn end(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.end.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.end.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=pdbHelixClass,unchecked_return_type="BioHelixClass")]
    pub fn pdb_helix_class(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.pdb_helix_class.into()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.pdb_helix_class as u32
    }
    #[wasm_bindgen(getter,js_name=length)]
    pub fn length(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.length
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.length
    }
}
#[wasm_bindgen]
pub struct BioSheet {
    pub(crate) inner: ck::BioSheet,
}
#[wasm_bindgen]
impl BioSheet {
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        // COSMolKit❗✔️: &self.inner.name
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.name.to_owned()
    }
    #[wasm_bindgen(getter,js_name=strands,unchecked_return_type="BioStrand[]")]
    pub fn strands(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .strands
            .iter()
            .map(|v| -> JsValue { BioStrand { inner: v.clone() }.into() })
            .collect()
    }
}
#[wasm_bindgen]
pub struct BioMetadata {
    pub(crate) inner: ck::BioMetadata,
}
#[wasm_bindgen]
impl BioMetadata {
    #[wasm_bindgen(getter,js_name=authors,unchecked_return_type="string[]")]
    pub fn authors(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.authors.clone()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .authors
            .iter()
            .map(|v| -> JsValue { v.as_str().into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=experiments,unchecked_return_type="BioExperimentInfo[]")]
    pub fn experiments(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .experiments
            .iter()
            .map(|v| -> JsValue { BioExperimentInfo { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=crystals,unchecked_return_type="BioExperimentalCrystalInfo[]")]
    pub fn crystals(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .crystals
            .iter()
            .map(|v| -> JsValue { BioExperimentalCrystalInfo { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=refinement,unchecked_return_type="BioRefinementInfo[]")]
    pub fn refinement(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .refinement
            .iter()
            .map(|v| -> JsValue { BioRefinementInfo { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=software,unchecked_return_type="BioSoftwareItem[]")]
    pub fn software(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .software
            .iter()
            .map(|v| -> JsValue { BioSoftwareItem { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=solvedBy)]
    pub fn solved_by(&self) -> String {
        // COSMolKit❗✔️: &self.inner.solved_by
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.solved_by.to_owned()
    }
    #[wasm_bindgen(getter,js_name=startingModel)]
    pub fn starting_model(&self) -> String {
        // COSMolKit❗✔️: &self.inner.starting_model
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.starting_model.to_owned()
    }
    #[wasm_bindgen(getter,js_name=remark_300Detail)]
    pub fn remark_300_detail(&self) -> String {
        // COSMolKit❗✔️: &self.inner.remark_300_detail
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.remark_300_detail.to_owned()
    }
}
#[wasm_bindgen]
pub struct BioStructureSourceState {
    pub(crate) inner: ck::BioStructureSourceState,
}
#[wasm_bindgen]
impl BioStructureSourceState {
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        // COSMolKit❗✔️: &self.inner.name
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.name.to_owned()
    }
    #[wasm_bindgen(getter,js_name=resolution)]
    pub fn resolution(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.resolution
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.resolution
    }
    #[wasm_bindgen(getter,js_name=conectMap,unchecked_return_type="Map<number,number[]>")]
    pub fn conect_map(&self) -> js_sys::Map {
        // COSMolKit❗✔️: self.inner.conect_map.clone()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        {
            let out = js_sys::Map::new();
            for (k, v) in &self.inner.conect_map {
                out.set(
                    &JsValue::from(*k),
                    &v.iter()
                        .copied()
                        .map(JsValue::from)
                        .collect::<js_sys::Array>()
                        .into(),
                );
            }
            out
        }
    }
    #[wasm_bindgen(getter,js_name=hasDFraction)]
    pub fn has_d_fraction(&self) -> bool {
        // COSMolKit❗✔️: self.inner.has_d_fraction
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.has_d_fraction
    }
    #[wasm_bindgen(getter,js_name=nonAsciiLine)]
    pub fn non_ascii_line(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.non_ascii_line
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.non_ascii_line
    }
    #[wasm_bindgen(getter,js_name=terStatus)]
    pub fn ter_status(&self) -> u8 {
        // COSMolKit❗✔️: self.inner.ter_status
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.ter_status
    }
    #[wasm_bindgen(getter,js_name=hasOrigx)]
    pub fn has_origx(&self) -> bool {
        // COSMolKit❗✔️: self.inner.has_origx
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.has_origx
    }
    #[wasm_bindgen(getter,js_name=origx)]
    pub fn origx(&self) -> BioTransform {
        // COSMolKit❗✔️: inner: self.inner.origx,
        // Only boundary value copies; source fields, order and sentinels are preserved.
        BioTransform {
            inner: self.inner.origx.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=info,unchecked_return_type="Map<string,string>")]
    pub fn info(&self) -> js_sys::Map {
        // COSMolKit❗✔️: self.inner.info.clone()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        {
            let out = js_sys::Map::new();
            for (k, v) in &self.inner.info {
                out.set(&JsValue::from(k.as_str()), &v.as_str().into());
            }
            out
        }
    }
    #[wasm_bindgen(getter,js_name=rawRemarks,unchecked_return_type="string[]")]
    pub fn raw_remarks(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.raw_remarks.clone()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .raw_remarks
            .iter()
            .map(|v| -> JsValue { v.as_str().into() })
            .collect()
    }
}
#[wasm_bindgen]
pub struct BioNcsOperator {
    pub(crate) inner: ck::BioNcsOperator,
}
#[wasm_bindgen]
impl BioNcsOperator {
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> String {
        // COSMolKit❗✔️: &self.inner.id
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.id.to_owned()
    }
    #[wasm_bindgen(getter,js_name=given)]
    pub fn given(&self) -> bool {
        // COSMolKit❗✔️: self.inner.given
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.given
    }
    #[wasm_bindgen(getter,js_name=transform)]
    pub fn transform(&self) -> BioTransform {
        // COSMolKit❗✔️: inner: self.inner.transform,
        // Only boundary value copies; source fields, order and sentinels are preserved.
        BioTransform {
            inner: self.inner.transform.clone(),
        }
    }
}
#[wasm_bindgen]
pub struct BioAssembly {
    pub(crate) inner: ck::BioAssembly,
}
#[wasm_bindgen]
impl BioAssembly {
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        // COSMolKit❗✔️: &self.inner.name
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.name.to_owned()
    }
    #[wasm_bindgen(getter,js_name=authorDetermined)]
    pub fn author_determined(&self) -> bool {
        // COSMolKit❗✔️: self.inner.author_determined
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.author_determined
    }
    #[wasm_bindgen(getter,js_name=softwareDetermined)]
    pub fn software_determined(&self) -> bool {
        // COSMolKit❗✔️: self.inner.software_determined
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.software_determined
    }
    #[wasm_bindgen(getter,js_name=specialKind,unchecked_return_type="BioAssemblySpecialKind")]
    pub fn special_kind(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.special_kind.into()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.special_kind as u32
    }
    #[wasm_bindgen(getter,js_name=oligomericCount)]
    pub fn oligomeric_count(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.oligomeric_count
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.oligomeric_count
    }
    #[wasm_bindgen(getter,js_name=oligomericDetails)]
    pub fn oligomeric_details(&self) -> String {
        // COSMolKit❗✔️: &self.inner.oligomeric_details
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.oligomeric_details.to_owned()
    }
    #[wasm_bindgen(getter,js_name=softwareName)]
    pub fn software_name(&self) -> String {
        // COSMolKit❗✔️: &self.inner.software_name
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.software_name.to_owned()
    }
    #[wasm_bindgen(getter,js_name=buriedSurfaceArea)]
    pub fn buried_surface_area(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.buried_surface_area
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.buried_surface_area
    }
    #[wasm_bindgen(getter,js_name=surfaceArea)]
    pub fn surface_area(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.surface_area
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.surface_area
    }
    #[wasm_bindgen(getter,js_name=solventFreeEnergyChange)]
    pub fn solvent_free_energy_change(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.solvent_free_energy_change
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.solvent_free_energy_change
    }
    #[wasm_bindgen(getter,js_name=generators,unchecked_return_type="BioAssemblyGenerator[]")]
    pub fn generators(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .generators
            .iter()
            .map(|v| -> JsValue { BioAssemblyGenerator { inner: v.clone() }.into() })
            .collect()
    }
}
#[wasm_bindgen]
pub struct BioAssemblyGenerator {
    pub(crate) inner: ck::BioAssemblyGenerator,
}
#[wasm_bindgen]
impl BioAssemblyGenerator {
    #[wasm_bindgen(getter,js_name=chains,unchecked_return_type="string[]")]
    pub fn chains(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.chains.clone()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .chains
            .iter()
            .map(|v| -> JsValue { v.as_str().into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=subchains,unchecked_return_type="string[]")]
    pub fn subchains(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.subchains.clone()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .subchains
            .iter()
            .map(|v| -> JsValue { v.as_str().into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=operators,unchecked_return_type="BioAssemblyOperator[]")]
    pub fn operators(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .operators
            .iter()
            .map(|v| -> JsValue { BioAssemblyOperator { inner: v.clone() }.into() })
            .collect()
    }
}
#[wasm_bindgen]
pub struct BioAssemblyOperator {
    pub(crate) inner: ck::BioAssemblyOperator,
}
#[wasm_bindgen]
impl BioAssemblyOperator {
    #[wasm_bindgen(getter,js_name=name,unchecked_return_type="string | null")]
    pub fn name(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.name.clone()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        nullable(self.inner.name.clone())
    }
    #[wasm_bindgen(getter,js_name=operatorType,unchecked_return_type="string | null")]
    pub fn operator_type(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.operator_type.clone()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        nullable(self.inner.operator_type.clone())
    }
    #[wasm_bindgen(getter,js_name=transform)]
    pub fn transform(&self) -> BioTransform {
        // COSMolKit❗✔️: inner: self.inner.transform,
        // Only boundary value copies; source fields, order and sentinels are preserved.
        BioTransform {
            inner: self.inner.transform.clone(),
        }
    }
}
#[wasm_bindgen]
pub struct BioStrand {
    pub(crate) inner: ck::BioStrand,
}
#[wasm_bindgen]
impl BioStrand {
    #[wasm_bindgen(getter,js_name=start)]
    pub fn start(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.start.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.start.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=end)]
    pub fn end(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.end.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.end.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=hbondAtom2)]
    pub fn hbond_atom2(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.hbond_atom2.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.hbond_atom2.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=hbondAtom1)]
    pub fn hbond_atom1(&self) -> AtomAddress {
        // COSMolKit❗✔️: inner: self.inner.hbond_atom1.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        AtomAddress {
            inner: self.inner.hbond_atom1.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=sense)]
    pub fn sense(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.sense
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.sense
    }
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        // COSMolKit❗✔️: &self.inner.name
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.name.to_owned()
    }
}
#[wasm_bindgen]
pub struct BioSoftwareItem {
    pub(crate) inner: ck::BioSoftwareItem,
}
#[wasm_bindgen]
impl BioSoftwareItem {
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        // COSMolKit❗✔️: &self.inner.name
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.name.to_owned()
    }
    #[wasm_bindgen(getter,js_name=version)]
    pub fn version(&self) -> String {
        // COSMolKit❗✔️: &self.inner.version
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.version.to_owned()
    }
    #[wasm_bindgen(getter,js_name=date)]
    pub fn date(&self) -> String {
        // COSMolKit❗✔️: &self.inner.date
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.date.to_owned()
    }
    #[wasm_bindgen(getter,js_name=description)]
    pub fn description(&self) -> String {
        // COSMolKit❗✔️: &self.inner.description
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.description.to_owned()
    }
    #[wasm_bindgen(getter,js_name=contactAuthor)]
    pub fn contact_author(&self) -> String {
        // COSMolKit❗✔️: &self.inner.contact_author
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.contact_author.to_owned()
    }
    #[wasm_bindgen(getter,js_name=contactAuthorEmail)]
    pub fn contact_author_email(&self) -> String {
        // COSMolKit❗✔️: &self.inner.contact_author_email
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.contact_author_email.to_owned()
    }
    #[wasm_bindgen(getter,js_name=classification,unchecked_return_type="BioSoftwareClassification")]
    pub fn classification(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.classification.into()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.classification as u32
    }
}
#[wasm_bindgen]
pub struct BioReflectionsInfo {
    pub(crate) inner: ck::BioReflectionsInfo,
}
#[wasm_bindgen]
impl BioReflectionsInfo {
    #[wasm_bindgen(getter,js_name=resolutionHigh)]
    pub fn resolution_high(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.resolution_high
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.resolution_high
    }
    #[wasm_bindgen(getter,js_name=resolutionLow)]
    pub fn resolution_low(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.resolution_low
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.resolution_low
    }
    #[wasm_bindgen(getter,js_name=completeness)]
    pub fn completeness(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.completeness
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.completeness
    }
    #[wasm_bindgen(getter,js_name=redundancy)]
    pub fn redundancy(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.redundancy
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.redundancy
    }
    #[wasm_bindgen(getter,js_name=rMerge)]
    pub fn r_merge(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.r_merge
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.r_merge
    }
    #[wasm_bindgen(getter,js_name=rSym)]
    pub fn r_sym(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.r_sym
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.r_sym
    }
    #[wasm_bindgen(getter,js_name=meanIOverSigma)]
    pub fn mean_i_over_sigma(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.mean_i_over_sigma
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.mean_i_over_sigma
    }
}
#[wasm_bindgen]
pub struct BioBasicRefinementInfo {
    pub(crate) inner: ck::BioBasicRefinementInfo,
}
#[wasm_bindgen]
impl BioBasicRefinementInfo {
    #[wasm_bindgen(getter,js_name=resolutionHigh)]
    pub fn resolution_high(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.resolution_high
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.resolution_high
    }
    #[wasm_bindgen(getter,js_name=resolutionLow)]
    pub fn resolution_low(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.resolution_low
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.resolution_low
    }
    #[wasm_bindgen(getter,js_name=completeness)]
    pub fn completeness(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.completeness
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.completeness
    }
    #[wasm_bindgen(getter,js_name=reflectionCount)]
    pub fn reflection_count(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.reflection_count
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.reflection_count
    }
    #[wasm_bindgen(getter,js_name=workSetCount)]
    pub fn work_set_count(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.work_set_count
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.work_set_count
    }
    #[wasm_bindgen(getter,js_name=rfreeSetCount)]
    pub fn rfree_set_count(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.rfree_set_count
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.rfree_set_count
    }
    #[wasm_bindgen(getter,js_name=rAll)]
    pub fn r_all(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.r_all
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.r_all
    }
    #[wasm_bindgen(getter,js_name=rWork)]
    pub fn r_work(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.r_work
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.r_work
    }
    #[wasm_bindgen(getter,js_name=rFree)]
    pub fn r_free(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.r_free
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.r_free
    }
    #[wasm_bindgen(getter,js_name=ccFoFcWork)]
    pub fn cc_fo_fc_work(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.cc_fo_fc_work
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.cc_fo_fc_work
    }
    #[wasm_bindgen(getter,js_name=ccFoFcFree)]
    pub fn cc_fo_fc_free(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.cc_fo_fc_free
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.cc_fo_fc_free
    }
    #[wasm_bindgen(getter,js_name=fscWork)]
    pub fn fsc_work(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.fsc_work
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.fsc_work
    }
    #[wasm_bindgen(getter,js_name=fscFree)]
    pub fn fsc_free(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.fsc_free
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.fsc_free
    }
    #[wasm_bindgen(getter,js_name=ccIntensityWork)]
    pub fn cc_intensity_work(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.cc_intensity_work
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.cc_intensity_work
    }
    #[wasm_bindgen(getter,js_name=ccIntensityFree)]
    pub fn cc_intensity_free(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.cc_intensity_free
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.cc_intensity_free
    }
}
#[wasm_bindgen]
pub struct BioRefinementRestraint {
    pub(crate) inner: ck::BioRefinementRestraint,
}
#[wasm_bindgen]
impl BioRefinementRestraint {
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        // COSMolKit❗✔️: &self.inner.name
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.name.to_owned()
    }
    #[wasm_bindgen(getter,js_name=count)]
    pub fn count(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.count
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.count
    }
    #[wasm_bindgen(getter,js_name=weight)]
    pub fn weight(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.weight
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.weight
    }
    #[wasm_bindgen(getter,js_name=function)]
    pub fn function(&self) -> String {
        // COSMolKit❗✔️: &self.inner.function
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.function.to_owned()
    }
    #[wasm_bindgen(getter,js_name=devIdeal)]
    pub fn dev_ideal(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.dev_ideal
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.dev_ideal
    }
}
#[wasm_bindgen]
pub struct BioExperimentInfo {
    pub(crate) inner: ck::BioExperimentInfo,
}
#[wasm_bindgen]
impl BioExperimentInfo {
    #[wasm_bindgen(getter,js_name=method)]
    pub fn method(&self) -> String {
        // COSMolKit❗✔️: &self.inner.method
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.method.to_owned()
    }
    #[wasm_bindgen(getter,js_name=numberOfCrystals)]
    pub fn number_of_crystals(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.number_of_crystals
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.number_of_crystals
    }
    #[wasm_bindgen(getter,js_name=uniqueReflections)]
    pub fn unique_reflections(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.unique_reflections
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.unique_reflections
    }
    #[wasm_bindgen(getter,js_name=reflections)]
    pub fn reflections(&self) -> BioReflectionsInfo {
        // COSMolKit❗✔️: inner: self.inner.reflections.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        BioReflectionsInfo {
            inner: self.inner.reflections.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=bWilson)]
    pub fn b_wilson(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.b_wilson
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.b_wilson
    }
    #[wasm_bindgen(getter,js_name=shells,unchecked_return_type="BioReflectionsInfo[]")]
    pub fn shells(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .shells
            .iter()
            .map(|v| -> JsValue { BioReflectionsInfo { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=diffractionIds,unchecked_return_type="string[]")]
    pub fn diffraction_ids(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.diffraction_ids.clone()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .diffraction_ids
            .iter()
            .map(|v| -> JsValue { v.as_str().into() })
            .collect()
    }
}
#[wasm_bindgen]
pub struct BioDiffractionInfo {
    pub(crate) inner: ck::BioDiffractionInfo,
}
#[wasm_bindgen]
impl BioDiffractionInfo {
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> String {
        // COSMolKit❗✔️: &self.inner.id
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.id.to_owned()
    }
    #[wasm_bindgen(getter,js_name=temperature)]
    pub fn temperature(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.temperature
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.temperature
    }
    #[wasm_bindgen(getter,js_name=source)]
    pub fn source(&self) -> String {
        // COSMolKit❗✔️: &self.inner.source
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.source.to_owned()
    }
    #[wasm_bindgen(getter,js_name=sourceType)]
    pub fn source_type(&self) -> String {
        // COSMolKit❗✔️: &self.inner.source_type
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.source_type.to_owned()
    }
    #[wasm_bindgen(getter,js_name=synchrotron)]
    pub fn synchrotron(&self) -> String {
        // COSMolKit❗✔️: &self.inner.synchrotron
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.synchrotron.to_owned()
    }
    #[wasm_bindgen(getter,js_name=beamline)]
    pub fn beamline(&self) -> String {
        // COSMolKit❗✔️: &self.inner.beamline
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.beamline.to_owned()
    }
    #[wasm_bindgen(getter,js_name=wavelengths)]
    pub fn wavelengths(&self) -> String {
        // COSMolKit❗✔️: &self.inner.wavelengths
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.wavelengths.to_owned()
    }
    #[wasm_bindgen(getter,js_name=scatteringType)]
    pub fn scattering_type(&self) -> String {
        // COSMolKit❗✔️: &self.inner.scattering_type
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.scattering_type.to_owned()
    }
    #[wasm_bindgen(getter,js_name=monoOrLaue)]
    pub fn mono_or_laue(&self) -> u8 {
        // COSMolKit❗✔️: self.inner.mono_or_laue
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.mono_or_laue
    }
    #[wasm_bindgen(getter,js_name=monochromator)]
    pub fn monochromator(&self) -> String {
        // COSMolKit❗✔️: &self.inner.monochromator
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.monochromator.to_owned()
    }
    #[wasm_bindgen(getter,js_name=collectionDate)]
    pub fn collection_date(&self) -> String {
        // COSMolKit❗✔️: &self.inner.collection_date
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.collection_date.to_owned()
    }
    #[wasm_bindgen(getter,js_name=optics)]
    pub fn optics(&self) -> String {
        // COSMolKit❗✔️: &self.inner.optics
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.optics.to_owned()
    }
    #[wasm_bindgen(getter,js_name=detector)]
    pub fn detector(&self) -> String {
        // COSMolKit❗✔️: &self.inner.detector
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.detector.to_owned()
    }
    #[wasm_bindgen(getter,js_name=detectorMake)]
    pub fn detector_make(&self) -> String {
        // COSMolKit❗✔️: &self.inner.detector_make
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.detector_make.to_owned()
    }
}
#[wasm_bindgen]
pub struct BioExperimentalCrystalInfo {
    pub(crate) inner: ck::BioExperimentalCrystalInfo,
}
#[wasm_bindgen]
impl BioExperimentalCrystalInfo {
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> String {
        // COSMolKit❗✔️: &self.inner.id
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.id.to_owned()
    }
    #[wasm_bindgen(getter,js_name=description)]
    pub fn description(&self) -> String {
        // COSMolKit❗✔️: &self.inner.description
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.description.to_owned()
    }
    #[wasm_bindgen(getter,js_name=ph)]
    pub fn ph(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.ph
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.ph
    }
    #[wasm_bindgen(getter,js_name=phRange)]
    pub fn ph_range(&self) -> String {
        // COSMolKit❗✔️: &self.inner.ph_range
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.ph_range.to_owned()
    }
    #[wasm_bindgen(getter,js_name=diffractions,unchecked_return_type="BioDiffractionInfo[]")]
    pub fn diffractions(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .diffractions
            .iter()
            .map(|v| -> JsValue { BioDiffractionInfo { inner: v.clone() }.into() })
            .collect()
    }
}
#[wasm_bindgen]
pub struct BioTlsSelection {
    pub(crate) inner: ck::BioTlsSelection,
}
#[wasm_bindgen]
impl BioTlsSelection {
    #[wasm_bindgen(getter,js_name=chain)]
    pub fn chain(&self) -> String {
        // COSMolKit❗✔️: self.inner.chain.as_str().to_owned()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.chain.as_str().to_owned()
    }
    #[wasm_bindgen(getter,js_name=resBegin)]
    pub fn res_begin(&self) -> PdbSeqId {
        // COSMolKit❗✔️: inner: self.inner.res_begin,
        // Only boundary value copies; source fields, order and sentinels are preserved.
        PdbSeqId {
            inner: self.inner.res_begin.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=resEnd)]
    pub fn res_end(&self) -> PdbSeqId {
        // COSMolKit❗✔️: inner: self.inner.res_end,
        // Only boundary value copies; source fields, order and sentinels are preserved.
        PdbSeqId {
            inner: self.inner.res_end.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=details)]
    pub fn details(&self) -> String {
        // COSMolKit❗✔️: &self.inner.details
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.details.to_owned()
    }
}
#[wasm_bindgen]
pub struct BioTlsGroup {
    pub(crate) inner: ck::BioTlsGroup,
}
#[wasm_bindgen]
impl BioTlsGroup {
    #[wasm_bindgen(getter,js_name=numId)]
    pub fn num_id(&self) -> i16 {
        // COSMolKit❗✔️: self.inner.num_id
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.num_id
    }
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> String {
        // COSMolKit❗✔️: &self.inner.id
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.id.to_owned()
    }
    #[wasm_bindgen(getter,js_name=selections,unchecked_return_type="BioTlsSelection[]")]
    pub fn selections(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .selections
            .iter()
            .map(|v| -> JsValue { BioTlsSelection { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=origin,unchecked_return_type="number[]")]
    pub fn origin(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.origin
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .origin
            .iter()
            .map(|v| JsValue::from(*v))
            .collect()
    }
    #[wasm_bindgen(getter,js_name=t,unchecked_return_type="number[]")]
    pub fn t(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.t
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.t.iter().map(|v| JsValue::from(*v)).collect()
    }
    #[wasm_bindgen(getter,js_name=l,unchecked_return_type="number[]")]
    pub fn l(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.l
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.l.iter().map(|v| JsValue::from(*v)).collect()
    }
    #[wasm_bindgen(getter,js_name=s,unchecked_return_type="number[][]")]
    pub fn s(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.s
        // Only boundary value copies; source fields, order and sentinels are preserved.
        matrix(&self.inner.s)
    }
}
#[wasm_bindgen]
pub struct BioRefinementInfo {
    pub(crate) inner: ck::BioRefinementInfo,
}
#[wasm_bindgen]
impl BioRefinementInfo {
    #[wasm_bindgen(getter,js_name=basic)]
    pub fn basic(&self) -> BioBasicRefinementInfo {
        // COSMolKit❗✔️: inner: self.inner.basic.clone(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        BioBasicRefinementInfo {
            inner: self.inner.basic.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> String {
        // COSMolKit❗✔️: &self.inner.id
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.id.to_owned()
    }
    #[wasm_bindgen(getter,js_name=crossValidationMethod)]
    pub fn cross_validation_method(&self) -> String {
        // COSMolKit❗✔️: &self.inner.cross_validation_method
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.cross_validation_method.to_owned()
    }
    #[wasm_bindgen(getter,js_name=rfreeSelectionMethod)]
    pub fn rfree_selection_method(&self) -> String {
        // COSMolKit❗✔️: &self.inner.rfree_selection_method
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.rfree_selection_method.to_owned()
    }
    #[wasm_bindgen(getter,js_name=binCount)]
    pub fn bin_count(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.bin_count
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.bin_count
    }
    #[wasm_bindgen(getter,js_name=bins,unchecked_return_type="BioBasicRefinementInfo[]")]
    pub fn bins(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .bins
            .iter()
            .map(|v| -> JsValue { BioBasicRefinementInfo { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=meanB)]
    pub fn mean_b(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.mean_b
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.mean_b
    }
    #[wasm_bindgen(getter,js_name=anisoB,unchecked_return_type="number[]")]
    pub fn aniso_b(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner.aniso_b
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .aniso_b
            .iter()
            .map(|v| JsValue::from(*v))
            .collect()
    }
    #[wasm_bindgen(getter,js_name=luzzatiError)]
    pub fn luzzati_error(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.luzzati_error
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.luzzati_error
    }
    #[wasm_bindgen(getter,js_name=dpiBlowR)]
    pub fn dpi_blow_r(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.dpi_blow_r
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.dpi_blow_r
    }
    #[wasm_bindgen(getter,js_name=dpiBlowRfree)]
    pub fn dpi_blow_rfree(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.dpi_blow_rfree
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.dpi_blow_rfree
    }
    #[wasm_bindgen(getter,js_name=dpiCruickshankR)]
    pub fn dpi_cruickshank_r(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.dpi_cruickshank_r
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.dpi_cruickshank_r
    }
    #[wasm_bindgen(getter,js_name=dpiCruickshankRfree)]
    pub fn dpi_cruickshank_rfree(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.dpi_cruickshank_rfree
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.dpi_cruickshank_rfree
    }
    #[wasm_bindgen(getter,js_name=restrStats,unchecked_return_type="BioRefinementRestraint[]")]
    pub fn restr_stats(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .restr_stats
            .iter()
            .map(|v| -> JsValue { BioRefinementRestraint { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=tlsGroups,unchecked_return_type="BioTlsGroup[]")]
    pub fn tls_groups(&self) -> js_sys::Array {
        // COSMolKit❗✔️: self.inner
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner
            .tls_groups
            .iter()
            .map(|v| -> JsValue { BioTlsGroup { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(getter,js_name=remarks)]
    pub fn remarks(&self) -> String {
        // COSMolKit❗✔️: &self.inner.remarks
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.remarks.to_owned()
    }
}
#[wasm_bindgen]
pub struct AtomAddress {
    pub(crate) inner: ck::AtomAddress,
}
#[wasm_bindgen]
impl AtomAddress {
    #[wasm_bindgen(getter,js_name=chainName)]
    pub fn chain_name(&self) -> String {
        // COSMolKit❗✔️: self.inner.chain_name().as_str().to_owned()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.chain_name().as_str().to_owned()
    }
    #[wasm_bindgen(getter,js_name=residue)]
    pub fn residue(&self) -> ResidueAddress {
        // COSMolKit❗✔️: inner: self.inner.residue(),
        // Only boundary value copies; source fields, order and sentinels are preserved.
        ResidueAddress {
            inner: self.inner.residue().clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=logicalAtomName)]
    pub fn logical_atom_name(&self) -> String {
        // COSMolKit❗✔️: self.inner.logical_atom_name()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.logical_atom_name().to_owned()
    }
    #[wasm_bindgen(getter,js_name=altloc)]
    pub fn altloc(&self) -> u8 {
        // COSMolKit❗✔️: self.inner.altloc()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.altloc()
    }
}
#[wasm_bindgen]
pub struct ResidueAddress {
    pub(crate) inner: ck::ResidueAddress,
}
#[wasm_bindgen]
impl ResidueAddress {
    #[wasm_bindgen(getter,js_name=sequenceNumber,unchecked_return_type="number | null")]
    pub fn sequence_number(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.sequence_number()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        nullable(self.inner.sequence_number())
    }
    #[wasm_bindgen(getter,js_name=insertionCode,unchecked_return_type="number | null")]
    pub fn insertion_code(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.insertion_code()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        nullable(self.inner.insertion_code())
    }
    #[wasm_bindgen(getter,js_name=segment)]
    pub fn segment(&self) -> String {
        // COSMolKit❗✔️: self.inner.segment()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.segment().to_owned()
    }
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        // COSMolKit❗✔️: self.inner.name().as_str().to_owned()
        // Only boundary value copies; source fields, order and sentinels are preserved.
        self.inner.name().as_str().to_owned()
    }
}
#[wasm_bindgen]
pub struct BioTransform {
    pub(crate) inner: ck::BioTransform,
}
#[wasm_bindgen]
impl BioTransform {
    #[wasm_bindgen(getter, unchecked_return_type = "number[][]")]
    pub fn matrix(&self) -> js_sys::Array {
        // COSMolKit❗✔️: *self.inner.matrix()
        matrix(self.inner.matrix())
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number[]")]
    pub fn translation(&self) -> js_sys::Array {
        // COSMolKit❗✔️: *self.inner.translation()
        row(self.inner.translation())
    }
}
#[wasm_bindgen]
pub struct BioCrystalCell {
    inner: ck::BioCrystalCell,
}
#[wasm_bindgen]
impl BioCrystalCell {
    #[wasm_bindgen(getter)]
    pub fn a(&self) -> f64 {
        self.inner.a
    }
    #[wasm_bindgen(getter)]
    pub fn b(&self) -> f64 {
        self.inner.b
    }
    #[wasm_bindgen(getter)]
    pub fn c(&self) -> f64 {
        self.inner.c
    }
    #[wasm_bindgen(getter)]
    pub fn alpha(&self) -> f64 {
        self.inner.alpha
    }
    #[wasm_bindgen(getter)]
    pub fn beta(&self) -> f64 {
        self.inner.beta
    }
    #[wasm_bindgen(getter)]
    pub fn gamma(&self) -> f64 {
        self.inner.gamma
    }
}
#[wasm_bindgen]
pub struct BioCrystalInfo {
    inner: Rc<ck::BioStructure>,
}
impl BioCrystalInfo {
    fn value(&self) -> &ck::BioCrystalInfo {
        self.inner
            .crystal()
            .expect("crystal view retains an existing immutable crystal")
    }
}
#[wasm_bindgen]
impl BioCrystalInfo {
    #[wasm_bindgen(js_name=spaceGroupNumber,unchecked_return_type="number | null")]
    pub fn space_group_number(&self) -> JsValue {
        // COSMolKit❗✔️: .space_group_number()
        nullable(self.value().space_group_number())
    }
    #[wasm_bindgen(getter)]
    pub fn cell(&self) -> BioCrystalCell {
        BioCrystalCell {
            inner: self.value().cell(),
        }
    }

    #[wasm_bindgen(getter,js_name=spaceGroupHm,unchecked_return_type="string | null")]
    pub fn space_group_hm(&self) -> JsValue {
        nullable(self.value().space_group_hm().map(str::to_owned))
    }
    #[wasm_bindgen(getter,js_name=zPdb,unchecked_return_type="string | null")]
    pub fn z_pdb(&self) -> JsValue {
        nullable(self.value().z_pdb().map(str::to_owned))
    }
    #[wasm_bindgen(getter)]
    pub fn orthogonal(&self) -> BioTransform {
        BioTransform {
            inner: *self.value().orthogonal(),
        }
    }
    #[wasm_bindgen(getter)]
    pub fn fractional(&self) -> BioTransform {
        BioTransform {
            inner: *self.value().fractional(),
        }
    }
    #[wasm_bindgen(getter,js_name=reciprocalLengths,unchecked_return_type="number[]")]
    pub fn reciprocal_lengths(&self) -> js_sys::Array {
        row(self.value().reciprocal_lengths())
    }
    #[wasm_bindgen(getter,js_name=reciprocalCosines,unchecked_return_type="number[]")]
    pub fn reciprocal_cosines(&self) -> js_sys::Array {
        row(self.value().reciprocal_cosines())
    }
    #[wasm_bindgen(getter,js_name=volume)]
    pub fn volume(&self) -> f64 {
        self.value().volume()
    }
    #[wasm_bindgen(getter,js_name=explicitMatrices)]
    pub fn explicit_matrices(&self) -> bool {
        self.value().explicit_matrices()
    }
    #[wasm_bindgen(getter,js_name=csCount)]
    pub fn cs_count(&self) -> i16 {
        self.value().cs_count()
    }
    #[wasm_bindgen(js_name=isCrystal)]
    pub fn is_crystal(&self) -> bool {
        self.value().is_crystal()
    }
    #[wasm_bindgen(getter,js_name=symmetryImages,unchecked_return_type="BioTransform[]")]
    pub fn symmetry_images(&self) -> js_sys::Array {
        self.value()
            .symmetry_images()
            .iter()
            .map(|v| -> JsValue { BioTransform { inner: *v }.into() })
            .collect()
    }
}
#[wasm_bindgen]
pub struct BioEntityDbRef {
    inner: ck::BioEntityDbRef,
}
#[wasm_bindgen]
impl BioEntityDbRef {
    #[wasm_bindgen(getter,js_name=dbName)]
    pub fn db_name(&self) -> String {
        self.inner.db_name.clone()
    }
    #[wasm_bindgen(getter,js_name=accessionCode)]
    pub fn accession_code(&self) -> String {
        self.inner.accession_code.clone()
    }
    #[wasm_bindgen(getter,js_name=idCode)]
    pub fn id_code(&self) -> String {
        self.inner.id_code.clone()
    }
    #[wasm_bindgen(getter,js_name=isoform)]
    pub fn isoform(&self) -> String {
        self.inner.isoform.clone()
    }
    #[wasm_bindgen(getter,js_name=seqBegin)]
    pub fn seq_begin(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.seq_begin,
        }
    }
    #[wasm_bindgen(getter,js_name=seqEnd)]
    pub fn seq_end(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.seq_end,
        }
    }
    #[wasm_bindgen(getter,js_name=dbBegin)]
    pub fn db_begin(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.db_begin,
        }
    }
    #[wasm_bindgen(getter,js_name=dbEnd)]
    pub fn db_end(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.db_end,
        }
    }
    #[wasm_bindgen(getter,js_name=labelSeqBegin,unchecked_return_type="number | null")]
    pub fn label_seq_begin(&self) -> JsValue {
        nullable(self.inner.label_seq_begin)
    }
    #[wasm_bindgen(getter,js_name=labelSeqEnd,unchecked_return_type="number | null")]
    pub fn label_seq_end(&self) -> JsValue {
        nullable(self.inner.label_seq_end)
    }
}
#[wasm_bindgen]
pub struct BioSiftsUnpResidue {
    inner: ck::BioSiftsUnpResidue,
}
#[wasm_bindgen]
impl BioSiftsUnpResidue {
    #[wasm_bindgen(unchecked_return_type = "number | null")]
    pub fn residue(&self) -> JsValue {
        nullable(self.inner.residue())
    }
    #[wasm_bindgen(js_name=accessionIndex)]
    pub fn accession_index(&self) -> u8 {
        self.inner.accession_index()
    }
    pub fn number(&self) -> u16 {
        self.inner.number()
    }
}
#[wasm_bindgen]
impl BioResidueRow {
    #[wasm_bindgen(js_name=siftsUnp)]
    pub fn sifts_unp(&self) -> BioSiftsUnpResidue {
        BioSiftsUnpResidue {
            inner: self.inner.residues()[self.index].sifts_unp(),
        }
    }
}
#[wasm_bindgen]
impl BioEntityRow {
    #[wasm_bindgen(unchecked_return_type = "BioEntityDbRef[]")]
    pub fn dbrefs(&self) -> js_sys::Array {
        self.inner.entities()[self.index]
            .dbrefs()
            .iter()
            .map(|v| -> JsValue { BioEntityDbRef { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(js_name=siftsUnpAccessions,unchecked_return_type="string[]")]
    pub fn sifts_unp_accessions(&self) -> js_sys::Array {
        self.inner.entities()[self.index]
            .sifts_unp_accessions()
            .iter()
            .map(|v| JsValue::from(v.as_str()))
            .collect()
    }
    #[wasm_bindgen(js_name=reflectsMicrohetero)]
    pub fn reflects_microhetero(&self) -> bool {
        self.inner.entities()[self.index].reflects_microhetero()
    }
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=connections,unchecked_return_type="BioConnection[]")]
    pub fn connections(&self) -> js_sys::Array {
        // COSMolKit❗✔️: .map(|inner| crate::canonical_bio_metadata::BioConnection { inner })
        self.inner
            .connections()
            .iter()
            .map(|v| -> JsValue { BioConnection { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(js_name=cispeps,unchecked_return_type="BioCisPep[]")]
    pub fn cispeps(&self) -> js_sys::Array {
        // COSMolKit❗✔️: .map(|inner| crate::canonical_bio_metadata::BioCisPep { inner })
        self.inner
            .cispeps()
            .iter()
            .map(|v| -> JsValue { BioCisPep { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(js_name=modResidues,unchecked_return_type="BioModRes[]")]
    pub fn mod_residues(&self) -> js_sys::Array {
        // COSMolKit❗✔️: .map(|inner| crate::canonical_bio_metadata::BioModRes { inner })
        self.inner
            .mod_residues()
            .iter()
            .map(|v| -> JsValue { BioModRes { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(js_name=helices,unchecked_return_type="BioHelix[]")]
    pub fn helices(&self) -> js_sys::Array {
        // COSMolKit❗✔️: .map(|inner| crate::canonical_bio_metadata::BioHelix { inner })
        self.inner
            .helices()
            .iter()
            .map(|v| -> JsValue { BioHelix { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(js_name=sheets,unchecked_return_type="BioSheet[]")]
    pub fn sheets(&self) -> js_sys::Array {
        // COSMolKit❗✔️: .map(|inner| crate::canonical_bio_metadata::BioSheet { inner })
        self.inner
            .sheets()
            .iter()
            .map(|v| -> JsValue { BioSheet { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(js_name=ncsOperators,unchecked_return_type="BioNcsOperator[]")]
    pub fn ncs_operators(&self) -> js_sys::Array {
        // COSMolKit❗✔️: .map(|inner| crate::canonical_bio_metadata::BioNcsOperator { inner })
        self.inner
            .ncs_operators()
            .iter()
            .map(|v| -> JsValue { BioNcsOperator { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(js_name=assemblies,unchecked_return_type="BioAssembly[]")]
    pub fn assemblies(&self) -> js_sys::Array {
        // COSMolKit❗✔️: .map(|inner| crate::canonical_bio_metadata::BioAssembly { inner })
        self.inner
            .assemblies()
            .iter()
            .map(|v| -> JsValue { BioAssembly { inner: v.clone() }.into() })
            .collect()
    }
    #[wasm_bindgen(js_name=metadata)]
    pub fn metadata(&self) -> BioMetadata {
        BioMetadata {
            inner: self.inner.metadata().clone(),
        }
    }
    #[wasm_bindgen(js_name=sourceState)]
    pub fn source_state(&self) -> BioStructureSourceState {
        BioStructureSourceState {
            inner: self.inner.source_state().clone(),
        }
    }
    #[wasm_bindgen(unchecked_return_type = "BioCrystalInfo | null")]
    pub fn crystal(&self) -> JsValue {
        nullable(self.inner.crystal().map(|_| BioCrystalInfo {
            inner: Rc::clone(&self.inner),
        }))
    }
    #[wasm_bindgen(js_name=hasOrigx)]
    pub fn has_origx(&self) -> bool {
        self.inner.has_origx()
    }
    pub fn origx(&self) -> BioTransform {
        BioTransform {
            inner: *self.inner.origx(),
        }
    }
    #[wasm_bindgen(js_name=ncsOperIdentityId,unchecked_return_type="string | null")]
    pub fn ncs_oper_identity_id(&self) -> JsValue {
        nullable(self.inner.ncs_oper_identity_id().map(str::to_owned))
    }
    pub fn resolution(&self) -> f64 {
        self.inner.resolution()
    }
    #[wasm_bindgen(js_name=terStatus)]
    pub fn ter_status(&self) -> u8 {
        self.inner.ter_status()
    }
}
