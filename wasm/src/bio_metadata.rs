//! Source metadata records. No interpretation or inferred defaults in bindings.
#[wasm_bindgen]
pub struct BioSoftwareItem(ck::BioSoftwareItem);
crate::fields::value_field!(
    ck::BioSoftwareItem,
    BioSoftwareItem,
    BioSoftwareItemInput,
    "BioSoftwareItem"
);
crate::fields::value_array!(
    ck::BioSoftwareItem,
    BioSoftwareItemArray,
    "ReadonlyArray<BioSoftwareItem>"
);
#[wasm_bindgen]
impl BioSoftwareItem {
    #[wasm_bindgen(constructor)]
    pub fn new(
        name: Option<Ts<input::TextInput>>,
        version: Option<Ts<input::TextInput>>,
        date: Option<Ts<input::TextInput>>,
        description: Option<Ts<input::TextInput>>,
        contact_author: Option<Ts<input::TextInput>>,
        contact_author_email: Option<Ts<input::TextInput>>,
        classification: Option<BioSoftwareClassificationInput>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioSoftwareItem::default();
        Ok(Self(ck::BioSoftwareItem {
            name: name
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.name),
            version: version
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.version),
            date: date
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.date),
            description: description
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.description),
            contact_author: contact_author
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.contact_author),
            contact_author_email: contact_author_email
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.contact_author_email),
            classification: classification
                .map(<ck::BioSoftwareClassification as Field>::input)
                .transpose()?
                .unwrap_or(d.classification),
        }))
    }
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        self.0.name.clone()
    }
    #[wasm_bindgen(getter,js_name=version)]
    pub fn version(&self) -> String {
        self.0.version.clone()
    }
    #[wasm_bindgen(getter,js_name=date)]
    pub fn date(&self) -> String {
        self.0.date.clone()
    }
    #[wasm_bindgen(getter,js_name=description)]
    pub fn description(&self) -> String {
        self.0.description.clone()
    }
    #[wasm_bindgen(getter,js_name=contactAuthor)]
    pub fn contact_author(&self) -> String {
        self.0.contact_author.clone()
    }
    #[wasm_bindgen(getter,js_name=contactAuthorEmail)]
    pub fn contact_author_email(&self) -> String {
        self.0.contact_author_email.clone()
    }
    #[wasm_bindgen(getter,js_name=classification)]
    pub fn classification(&self) -> BioSoftwareClassification {
        self.0.classification.into()
    }
}
#[wasm_bindgen]
pub struct BioExperimentInfo(ck::BioExperimentInfo);
crate::fields::value_field!(
    ck::BioExperimentInfo,
    BioExperimentInfo,
    BioExperimentInfoInput,
    "BioExperimentInfo"
);
crate::fields::value_array!(
    ck::BioExperimentInfo,
    BioExperimentInfoArray,
    "ReadonlyArray<BioExperimentInfo>"
);
#[wasm_bindgen]
impl BioExperimentInfo {
    #[wasm_bindgen(constructor)]
    pub fn new(
        method: Option<Ts<input::TextInput>>,
        number_of_crystals: Option<Ts<input::I32Input>>,
        unique_reflections: Option<Ts<input::I32Input>>,
        b_wilson: Option<Ts<input::F64Input>>,
        diffraction_ids: Option<Ts<TextArray>>,
        reflections: Option<BioReflectionsInfoInput>,
        shells: Option<BioReflectionsInfoArray>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioExperimentInfo::default();
        Ok(Self(ck::BioExperimentInfo {
            method: method
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.method),
            number_of_crystals: number_of_crystals
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.number_of_crystals),
            unique_reflections: unique_reflections
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.unique_reflections),
            b_wilson: b_wilson
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.b_wilson),
            diffraction_ids: diffraction_ids
                .map(<Vec<String> as Field>::input)
                .transpose()?
                .unwrap_or(d.diffraction_ids),
            reflections: reflections
                .map(<ck::BioReflectionsInfo as Field>::input)
                .transpose()?
                .unwrap_or(d.reflections),
            shells: shells
                .map(<Vec<ck::BioReflectionsInfo> as Field>::input)
                .transpose()?
                .unwrap_or(d.shells),
        }))
    }
    #[wasm_bindgen(getter,js_name=method)]
    pub fn method(&self) -> String {
        self.0.method.clone()
    }
    #[wasm_bindgen(getter,js_name=numberOfCrystals)]
    pub fn number_of_crystals(&self) -> i32 {
        self.0.number_of_crystals
    }
    #[wasm_bindgen(getter,js_name=uniqueReflections)]
    pub fn unique_reflections(&self) -> i32 {
        self.0.unique_reflections
    }
    #[wasm_bindgen(getter,js_name=bWilson)]
    pub fn b_wilson(&self) -> f64 {
        self.0.b_wilson
    }
    #[wasm_bindgen(getter,js_name=diffractionIds)]
    pub fn diffraction_ids(&self) -> Result<Ts<TextArray>, JsValue> {
        <Vec<String> as Field>::get(NativeField::owned(self.0.diffraction_ids.clone()))
    }
    #[wasm_bindgen(getter,js_name=reflections)]
    pub fn reflections(&self) -> BioReflectionsInfo {
        BioReflectionsInfo(self.0.reflections.clone())
    }
    #[wasm_bindgen(getter,js_name=shells)]
    pub fn shells(&self) -> Result<BioReflectionsInfoArray, JsValue> {
        <Vec<ck::BioReflectionsInfo> as Field>::get(NativeField::owned(self.0.shells.clone()))
    }
}
#[wasm_bindgen]
pub struct BioExperimentalCrystalInfo(ck::BioExperimentalCrystalInfo);
crate::fields::value_field!(
    ck::BioExperimentalCrystalInfo,
    BioExperimentalCrystalInfo,
    BioExperimentalCrystalInfoInput,
    "BioExperimentalCrystalInfo"
);
crate::fields::value_array!(
    ck::BioExperimentalCrystalInfo,
    BioExperimentalCrystalInfoArray,
    "ReadonlyArray<BioExperimentalCrystalInfo>"
);
#[wasm_bindgen]
impl BioExperimentalCrystalInfo {
    #[wasm_bindgen(constructor)]
    pub fn new(
        id: Option<Ts<input::TextInput>>,
        description: Option<Ts<input::TextInput>>,
        ph: Option<Ts<input::F64Input>>,
        ph_range: Option<Ts<input::TextInput>>,
        diffractions: Option<BioDiffractionInfoArray>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioExperimentalCrystalInfo::default();
        Ok(Self(ck::BioExperimentalCrystalInfo {
            id: id.map(input::read::<String>).transpose()?.unwrap_or(d.id),
            description: description
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.description),
            ph: ph.map(input::read::<f64>).transpose()?.unwrap_or(d.ph),
            ph_range: ph_range
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.ph_range),
            diffractions: diffractions
                .map(<Vec<ck::BioDiffractionInfo> as Field>::input)
                .transpose()?
                .unwrap_or(d.diffractions),
        }))
    }
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> String {
        self.0.id.clone()
    }
    #[wasm_bindgen(getter,js_name=description)]
    pub fn description(&self) -> String {
        self.0.description.clone()
    }
    #[wasm_bindgen(getter,js_name=ph)]
    pub fn ph(&self) -> f64 {
        self.0.ph
    }
    #[wasm_bindgen(getter,js_name=phRange)]
    pub fn ph_range(&self) -> String {
        self.0.ph_range.clone()
    }
    #[wasm_bindgen(getter,js_name=diffractions)]
    pub fn diffractions(&self) -> Result<BioDiffractionInfoArray, JsValue> {
        <Vec<ck::BioDiffractionInfo> as Field>::get(NativeField::owned(self.0.diffractions.clone()))
    }
}
#[wasm_bindgen]
pub struct BioTlsGroup(ck::BioTlsGroup);
crate::fields::value_field!(
    ck::BioTlsGroup,
    BioTlsGroup,
    BioTlsGroupInput,
    "BioTlsGroup"
);
crate::fields::value_array!(
    ck::BioTlsGroup,
    BioTlsGroupArray,
    "ReadonlyArray<BioTlsGroup>"
);
#[wasm_bindgen]
impl BioTlsGroup {
    #[wasm_bindgen(constructor)]
    pub fn new(
        num_id: Option<Ts<input::I16Input>>,
        id: Option<Ts<input::TextInput>>,
        origin: Option<Ts<Point3>>,
        t: Option<Ts<Anisou>>,
        l: Option<Ts<Anisou>>,
        s: Option<Ts<Matrix3>>,
        selections: Option<BioTlsSelectionArray>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioTlsGroup::default();
        Ok(Self(ck::BioTlsGroup {
            num_id: num_id
                .map(input::read::<i16>)
                .transpose()?
                .unwrap_or(d.num_id),
            id: id.map(input::read::<String>).transpose()?.unwrap_or(d.id),
            origin: origin
                .map(input::decode)
                .transpose()?
                .map(|v| v.0)
                .unwrap_or(d.origin),
            t: t.map(input::decode)
                .transpose()?
                .map(|v| v.0)
                .unwrap_or(d.t),
            l: l.map(input::decode)
                .transpose()?
                .map(|v| v.0)
                .unwrap_or(d.l),
            s: s.map(input::decode)
                .transpose()?
                .map(|v| v.0.map(|p| p.0))
                .unwrap_or(d.s),
            selections: selections
                .map(<Vec<ck::BioTlsSelection> as Field>::input)
                .transpose()?
                .unwrap_or(d.selections),
        }))
    }
    #[wasm_bindgen(getter,js_name=numId)]
    pub fn num_id(&self) -> i16 {
        self.0.num_id
    }
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> String {
        self.0.id.clone()
    }
    #[wasm_bindgen(getter,js_name=origin)]
    pub fn origin(&self) -> Result<Ts<Point3>, JsValue> {
        Point3(self.0.origin)
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=t)]
    pub fn t(&self) -> Result<Ts<Anisou>, JsValue> {
        Anisou(self.0.t)
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=l)]
    pub fn l(&self) -> Result<Ts<Anisou>, JsValue> {
        Anisou(self.0.l)
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=s)]
    pub fn s(&self) -> Result<Ts<Matrix3>, JsValue> {
        Matrix3(self.0.s.map(Point3))
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=selections)]
    pub fn selections(&self) -> Result<BioTlsSelectionArray, JsValue> {
        <Vec<ck::BioTlsSelection> as Field>::get(NativeField::owned(self.0.selections.clone()))
    }
}
#[wasm_bindgen]
pub struct BioRefinementInfo(ck::BioRefinementInfo);
crate::fields::value_field!(
    ck::BioRefinementInfo,
    BioRefinementInfo,
    BioRefinementInfoInput,
    "BioRefinementInfo"
);
crate::fields::value_array!(
    ck::BioRefinementInfo,
    BioRefinementInfoArray,
    "ReadonlyArray<BioRefinementInfo>"
);
#[wasm_bindgen]
impl BioRefinementInfo {
    #[wasm_bindgen(constructor)]
    pub fn new(
        id: Option<Ts<input::TextInput>>,
        cross_validation_method: Option<Ts<input::TextInput>>,
        rfree_selection_method: Option<Ts<input::TextInput>>,
        bin_count: Option<Ts<input::I32Input>>,
        mean_b: Option<Ts<input::F64Input>>,
        aniso_b: Option<Ts<Anisou>>,
        luzzati_error: Option<Ts<input::F64Input>>,
        dpi_blow_r: Option<Ts<input::F64Input>>,
        dpi_blow_rfree: Option<Ts<input::F64Input>>,
        dpi_cruickshank_r: Option<Ts<input::F64Input>>,
        dpi_cruickshank_rfree: Option<Ts<input::F64Input>>,
        remarks: Option<Ts<input::TextInput>>,
        basic: Option<BioBasicRefinementInfoInput>,
        bins: Option<BioBasicRefinementInfoArray>,
        restr_stats: Option<BioRefinementRestraintArray>,
        tls_groups: Option<BioTlsGroupArray>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioRefinementInfo::default();
        Ok(Self(ck::BioRefinementInfo {
            id: id.map(input::read::<String>).transpose()?.unwrap_or(d.id),
            cross_validation_method: cross_validation_method
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.cross_validation_method),
            rfree_selection_method: rfree_selection_method
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.rfree_selection_method),
            bin_count: bin_count
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.bin_count),
            mean_b: mean_b
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.mean_b),
            aniso_b: aniso_b
                .map(input::decode)
                .transpose()?
                .map(|v| v.0)
                .unwrap_or(d.aniso_b),
            luzzati_error: luzzati_error
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.luzzati_error),
            dpi_blow_r: dpi_blow_r
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.dpi_blow_r),
            dpi_blow_rfree: dpi_blow_rfree
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.dpi_blow_rfree),
            dpi_cruickshank_r: dpi_cruickshank_r
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.dpi_cruickshank_r),
            dpi_cruickshank_rfree: dpi_cruickshank_rfree
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.dpi_cruickshank_rfree),
            remarks: remarks
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.remarks),
            basic: basic
                .map(<ck::BioBasicRefinementInfo as Field>::input)
                .transpose()?
                .unwrap_or(d.basic),
            bins: bins
                .map(<Vec<ck::BioBasicRefinementInfo> as Field>::input)
                .transpose()?
                .unwrap_or(d.bins),
            restr_stats: restr_stats
                .map(<Vec<ck::BioRefinementRestraint> as Field>::input)
                .transpose()?
                .unwrap_or(d.restr_stats),
            tls_groups: tls_groups
                .map(<Vec<ck::BioTlsGroup> as Field>::input)
                .transpose()?
                .unwrap_or(d.tls_groups),
        }))
    }
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> String {
        self.0.id.clone()
    }
    #[wasm_bindgen(getter,js_name=crossValidationMethod)]
    pub fn cross_validation_method(&self) -> String {
        self.0.cross_validation_method.clone()
    }
    #[wasm_bindgen(getter,js_name=rfreeSelectionMethod)]
    pub fn rfree_selection_method(&self) -> String {
        self.0.rfree_selection_method.clone()
    }
    #[wasm_bindgen(getter,js_name=binCount)]
    pub fn bin_count(&self) -> i32 {
        self.0.bin_count
    }
    #[wasm_bindgen(getter,js_name=meanB)]
    pub fn mean_b(&self) -> f64 {
        self.0.mean_b
    }
    #[wasm_bindgen(getter,js_name=anisoB)]
    pub fn aniso_b(&self) -> Result<Ts<Anisou>, JsValue> {
        Anisou(self.0.aniso_b)
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
    #[wasm_bindgen(getter,js_name=luzzatiError)]
    pub fn luzzati_error(&self) -> f64 {
        self.0.luzzati_error
    }
    #[wasm_bindgen(getter,js_name=dpiBlowR)]
    pub fn dpi_blow_r(&self) -> f64 {
        self.0.dpi_blow_r
    }
    #[wasm_bindgen(getter,js_name=dpiBlowRfree)]
    pub fn dpi_blow_rfree(&self) -> f64 {
        self.0.dpi_blow_rfree
    }
    #[wasm_bindgen(getter,js_name=dpiCruickshankR)]
    pub fn dpi_cruickshank_r(&self) -> f64 {
        self.0.dpi_cruickshank_r
    }
    #[wasm_bindgen(getter,js_name=dpiCruickshankRfree)]
    pub fn dpi_cruickshank_rfree(&self) -> f64 {
        self.0.dpi_cruickshank_rfree
    }
    #[wasm_bindgen(getter,js_name=remarks)]
    pub fn remarks(&self) -> String {
        self.0.remarks.clone()
    }
    #[wasm_bindgen(getter,js_name=basic)]
    pub fn basic(&self) -> BioBasicRefinementInfo {
        BioBasicRefinementInfo(self.0.basic.clone())
    }
    #[wasm_bindgen(getter,js_name=bins)]
    pub fn bins(&self) -> Result<BioBasicRefinementInfoArray, JsValue> {
        <Vec<ck::BioBasicRefinementInfo> as Field>::get(NativeField::owned(self.0.bins.clone()))
    }
    #[wasm_bindgen(getter,js_name=restrStats)]
    pub fn restr_stats(&self) -> Result<BioRefinementRestraintArray, JsValue> {
        <Vec<ck::BioRefinementRestraint> as Field>::get(NativeField::owned(
            self.0.restr_stats.clone(),
        ))
    }
    #[wasm_bindgen(getter,js_name=tlsGroups)]
    pub fn tls_groups(&self) -> Result<BioTlsGroupArray, JsValue> {
        <Vec<ck::BioTlsGroup> as Field>::get(NativeField::owned(self.0.tls_groups.clone()))
    }
}
#[wasm_bindgen]
pub struct BioMetadata(ck::BioMetadata);
crate::fields::value_field!(
    ck::BioMetadata,
    BioMetadata,
    BioMetadataInput,
    "BioMetadata"
);
crate::fields::value_array!(
    ck::BioMetadata,
    BioMetadataArray,
    "ReadonlyArray<BioMetadata>"
);
#[wasm_bindgen]
impl BioMetadata {
    #[wasm_bindgen(constructor)]
    pub fn new(
        authors: Option<Ts<TextArray>>,
        solved_by: Option<Ts<input::TextInput>>,
        starting_model: Option<Ts<input::TextInput>>,
        remark_300_detail: Option<Ts<input::TextInput>>,
        experiments: Option<BioExperimentInfoArray>,
        crystals: Option<BioExperimentalCrystalInfoArray>,
        refinement: Option<BioRefinementInfoArray>,
        software: Option<BioSoftwareItemArray>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioMetadata::default();
        Ok(Self(ck::BioMetadata {
            authors: authors
                .map(<Vec<String> as Field>::input)
                .transpose()?
                .unwrap_or(d.authors),
            solved_by: solved_by
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.solved_by),
            starting_model: starting_model
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.starting_model),
            remark_300_detail: remark_300_detail
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.remark_300_detail),
            experiments: experiments
                .map(<Vec<ck::BioExperimentInfo> as Field>::input)
                .transpose()?
                .unwrap_or(d.experiments),
            crystals: crystals
                .map(<Vec<ck::BioExperimentalCrystalInfo> as Field>::input)
                .transpose()?
                .unwrap_or(d.crystals),
            refinement: refinement
                .map(<Vec<ck::BioRefinementInfo> as Field>::input)
                .transpose()?
                .unwrap_or(d.refinement),
            software: software
                .map(<Vec<ck::BioSoftwareItem> as Field>::input)
                .transpose()?
                .unwrap_or(d.software),
        }))
    }
    #[wasm_bindgen(getter,js_name=authors)]
    pub fn authors(&self) -> Result<Ts<TextArray>, JsValue> {
        <Vec<String> as Field>::get(NativeField::owned(self.0.authors.clone()))
    }
    #[wasm_bindgen(getter,js_name=solvedBy)]
    pub fn solved_by(&self) -> String {
        self.0.solved_by.clone()
    }
    #[wasm_bindgen(getter,js_name=startingModel)]
    pub fn starting_model(&self) -> String {
        self.0.starting_model.clone()
    }
    #[wasm_bindgen(getter,js_name=remark_300Detail)]
    pub fn remark_300_detail(&self) -> String {
        self.0.remark_300_detail.clone()
    }
    #[wasm_bindgen(getter,js_name=experiments)]
    pub fn experiments(&self) -> Result<BioExperimentInfoArray, JsValue> {
        <Vec<ck::BioExperimentInfo> as Field>::get(NativeField::owned(self.0.experiments.clone()))
    }
    #[wasm_bindgen(getter,js_name=crystals)]
    pub fn crystals(&self) -> Result<BioExperimentalCrystalInfoArray, JsValue> {
        <Vec<ck::BioExperimentalCrystalInfo> as Field>::get(NativeField::owned(
            self.0.crystals.clone(),
        ))
    }
    #[wasm_bindgen(getter,js_name=refinement)]
    pub fn refinement(&self) -> Result<BioRefinementInfoArray, JsValue> {
        <Vec<ck::BioRefinementInfo> as Field>::get(NativeField::owned(self.0.refinement.clone()))
    }
    #[wasm_bindgen(getter,js_name=software)]
    pub fn software(&self) -> Result<BioSoftwareItemArray, JsValue> {
        <Vec<ck::BioSoftwareItem> as Field>::get(NativeField::owned(self.0.software.clone()))
    }
}

#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(getter)]
    pub fn metadata(&self) -> BioMetadata {
        BioMetadata(self.0.metadata().clone())
    }
}

use crate::{
    bio_geometry::Matrix3,
    bio_hierarchy::BioStructure,
    bio_rows::Anisou,
    bio_scalar_metadata::*,
    bio_source_values::{PdbSeqId, PdbSeqIdInput},
    fields::{Field, TextArray},
    group_values::Point3,
    input,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
crate::enums::enum_value!(BioSoftwareClassification,"BioSoftwareClassification"{DataCollection="data_collection",DataExtraction="data_extraction",DataProcessing="data_processing",DataReduction="data_reduction",DataScaling="data_scaling",ModelBuilding="model_building",Phasing="phasing",Refinement="refinement",Unspecified="unspecified"});
#[wasm_bindgen]
pub struct BioTlsSelection(ck::BioTlsSelection);
crate::fields::value_field!(
    ck::BioTlsSelection,
    BioTlsSelection,
    BioTlsSelectionInput,
    "BioTlsSelection"
);
crate::fields::value_array!(
    ck::BioTlsSelection,
    BioTlsSelectionArray,
    "ReadonlyArray<BioTlsSelection>"
);
#[wasm_bindgen]
impl BioTlsSelection {
    #[wasm_bindgen(constructor)]
    pub fn new(
        chain: Option<Ts<input::TextInput>>,
        res_begin: Option<PdbSeqIdInput>,
        res_end: Option<PdbSeqIdInput>,
        details: Option<Ts<input::TextInput>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioTlsSelection::default();
        let chain = chain
            .map(input::read::<String>)
            .transpose()?
            .map(|v| {
                ck::PdbChainId::from_ascii(v.as_bytes()).ok_or_else(|| {
                    JsValue::from(js_sys::RangeError::new("chain is not a valid PdbChainId"))
                })
            })
            .transpose()?
            .unwrap_or(d.chain);
        Ok(Self(ck::BioTlsSelection {
            chain,
            res_begin: res_begin
                .map(<ck::PdbSeqId as Field>::input)
                .transpose()?
                .unwrap_or(d.res_begin),
            res_end: res_end
                .map(<ck::PdbSeqId as Field>::input)
                .transpose()?
                .unwrap_or(d.res_end),
            details: details
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.details),
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn chain(&self) -> String {
        self.0.chain.as_str().into()
    }
    #[wasm_bindgen(getter,js_name=resBegin)]
    pub fn res_begin(&self) -> PdbSeqId {
        PdbSeqId(self.0.res_begin)
    }
    #[wasm_bindgen(getter,js_name=resEnd)]
    pub fn res_end(&self) -> PdbSeqId {
        PdbSeqId(self.0.res_end)
    }
    #[wasm_bindgen(getter)]
    pub fn details(&self) -> String {
        self.0.details.clone()
    }
}
