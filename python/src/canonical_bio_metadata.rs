//! Complete readonly projections of the facade's detached BIO metadata.
//! Scalars retain their bits, bytes and missing-value sentinels. Collections
//! retain source order and duplicates; no parsing or chemistry occurs here.

use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
use std::collections::BTreeMap;

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum BioConnectionKind {
    Covale,
    Disulf,
    Hydrog,
    MetalC,
    Unknown,
}
impl From<ck::BioConnectionKind> for BioConnectionKind {
    fn from(value: ck::BioConnectionKind) -> Self {
        match value {
            ck::BioConnectionKind::Covale => Self::Covale,
            ck::BioConnectionKind::Disulf => Self::Disulf,
            ck::BioConnectionKind::Hydrog => Self::Hydrog,
            ck::BioConnectionKind::MetalC => Self::MetalC,
            ck::BioConnectionKind::Unknown => Self::Unknown,
        }
    }
}

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum BioAsu {
    Same,
    Different,
    Any,
}
impl From<ck::BioAsu> for BioAsu {
    fn from(value: ck::BioAsu) -> Self {
        match value {
            ck::BioAsu::Same => Self::Same,
            ck::BioAsu::Different => Self::Different,
            ck::BioAsu::Any => Self::Any,
        }
    }
}

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum BioHelixClass {
    UnknownHelix,
    RAlpha,
    ROmega,
    RPi,
    RGamma,
    R310,
    LAlpha,
    LOmega,
    LGamma,
    Helix27,
    HelixPolyProlineNone,
}
impl From<ck::BioHelixClass> for BioHelixClass {
    fn from(value: ck::BioHelixClass) -> Self {
        match value {
            ck::BioHelixClass::UnknownHelix => Self::UnknownHelix,
            ck::BioHelixClass::RAlpha => Self::RAlpha,
            ck::BioHelixClass::ROmega => Self::ROmega,
            ck::BioHelixClass::RPi => Self::RPi,
            ck::BioHelixClass::RGamma => Self::RGamma,
            ck::BioHelixClass::R310 => Self::R310,
            ck::BioHelixClass::LAlpha => Self::LAlpha,
            ck::BioHelixClass::LOmega => Self::LOmega,
            ck::BioHelixClass::LGamma => Self::LGamma,
            ck::BioHelixClass::Helix27 => Self::Helix27,
            ck::BioHelixClass::HelixPolyProlineNone => Self::HelixPolyProlineNone,
        }
    }
}

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum BioSoftwareClassification {
    DataCollection,
    DataExtraction,
    DataProcessing,
    DataReduction,
    DataScaling,
    ModelBuilding,
    Phasing,
    Refinement,
    Unspecified,
}
impl From<ck::BioSoftwareClassification> for BioSoftwareClassification {
    fn from(value: ck::BioSoftwareClassification) -> Self {
        match value {
            ck::BioSoftwareClassification::DataCollection => Self::DataCollection,
            ck::BioSoftwareClassification::DataExtraction => Self::DataExtraction,
            ck::BioSoftwareClassification::DataProcessing => Self::DataProcessing,
            ck::BioSoftwareClassification::DataReduction => Self::DataReduction,
            ck::BioSoftwareClassification::DataScaling => Self::DataScaling,
            ck::BioSoftwareClassification::ModelBuilding => Self::ModelBuilding,
            ck::BioSoftwareClassification::Phasing => Self::Phasing,
            ck::BioSoftwareClassification::Refinement => Self::Refinement,
            ck::BioSoftwareClassification::Unspecified => Self::Unspecified,
        }
    }
}

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum BioAssemblySpecialKind {
    NotApplicable,
    CompleteIcosahedral,
    RepresentativeHelical,
    CompletePoint,
}
impl From<ck::BioAssemblySpecialKind> for BioAssemblySpecialKind {
    fn from(value: ck::BioAssemblySpecialKind) -> Self {
        match value {
            ck::BioAssemblySpecialKind::NotApplicable => Self::NotApplicable,
            ck::BioAssemblySpecialKind::CompleteIcosahedral => Self::CompleteIcosahedral,
            ck::BioAssemblySpecialKind::RepresentativeHelical => Self::RepresentativeHelical,
            ck::BioAssemblySpecialKind::CompletePoint => Self::CompletePoint,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioConnection {
    pub(crate) inner: ck::BioConnection,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioConnection {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn link_id(&self) -> &str {
        &self.inner.link_id
    }
    #[getter]
    fn kind(&self) -> BioConnectionKind {
        self.inner.kind.into()
    }
    #[getter]
    fn asu(&self) -> BioAsu {
        self.inner.asu.into()
    }
    #[getter]
    fn partner1(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner1.clone(),
        }
    }
    #[getter]
    fn partner2(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner2.clone(),
        }
    }
    #[getter]
    fn reported_distance(&self) -> f64 {
        self.inner.reported_distance
    }
    #[getter]
    fn reported_sym(&self) -> [i16; 4] {
        self.inner.reported_sym
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioCisPep {
    pub(crate) inner: ck::BioCisPep,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioCisPep {
    #[getter]
    fn partner_c(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner_c.clone(),
        }
    }
    #[getter]
    fn partner_n(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner_n.clone(),
        }
    }
    #[getter]
    fn model_num(&self) -> i32 {
        self.inner.model_num
    }
    #[getter]
    fn only_altloc(&self) -> u8 {
        self.inner.only_altloc
    }
    #[getter]
    fn reported_angle(&self) -> f64 {
        self.inner.reported_angle
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioModRes {
    pub(crate) inner: ck::BioModRes,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioModRes {
    #[getter]
    fn chain_name(&self) -> String {
        self.inner.chain_name.as_str().to_owned()
    }
    #[getter]
    fn res_id(&self) -> ResidueAddress {
        ResidueAddress {
            inner: self.inner.res_id.clone(),
        }
    }
    #[getter]
    fn parent_comp_id(&self) -> &str {
        &self.inner.parent_comp_id
    }
    #[getter]
    fn mod_id(&self) -> &str {
        &self.inner.mod_id
    }
    #[getter]
    fn details(&self) -> &str {
        &self.inner.details
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioHelix {
    pub(crate) inner: ck::BioHelix,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioHelix {
    #[getter]
    fn start(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.start.clone(),
        }
    }
    #[getter]
    fn end(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.end.clone(),
        }
    }
    #[getter]
    fn pdb_helix_class(&self) -> BioHelixClass {
        self.inner.pdb_helix_class.into()
    }
    #[getter]
    fn length(&self) -> i32 {
        self.inner.length
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioSheet {
    pub(crate) inner: ck::BioSheet,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioSheet {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn strands(&self) -> Vec<BioStrand> {
        self.inner
            .strands
            .iter()
            .cloned()
            .map(|inner| BioStrand { inner })
            .collect()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioMetadata {
    pub(crate) inner: ck::BioMetadata,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioMetadata {
    #[getter]
    fn authors(&self) -> Vec<String> {
        self.inner.authors.clone()
    }
    #[getter]
    fn experiments(&self) -> Vec<BioExperimentInfo> {
        self.inner
            .experiments
            .iter()
            .cloned()
            .map(|inner| BioExperimentInfo { inner })
            .collect()
    }
    #[getter]
    fn crystals(&self) -> Vec<BioExperimentalCrystalInfo> {
        self.inner
            .crystals
            .iter()
            .cloned()
            .map(|inner| BioExperimentalCrystalInfo { inner })
            .collect()
    }
    #[getter]
    fn refinement(&self) -> Vec<BioRefinementInfo> {
        self.inner
            .refinement
            .iter()
            .cloned()
            .map(|inner| BioRefinementInfo { inner })
            .collect()
    }
    #[getter]
    fn software(&self) -> Vec<BioSoftwareItem> {
        self.inner
            .software
            .iter()
            .cloned()
            .map(|inner| BioSoftwareItem { inner })
            .collect()
    }
    #[getter]
    fn solved_by(&self) -> &str {
        &self.inner.solved_by
    }
    #[getter]
    fn starting_model(&self) -> &str {
        &self.inner.starting_model
    }
    #[getter]
    fn remark_300_detail(&self) -> &str {
        &self.inner.remark_300_detail
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioStructureSourceState {
    pub(crate) inner: ck::BioStructureSourceState,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioStructureSourceState {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn resolution(&self) -> f64 {
        self.inner.resolution
    }
    #[getter]
    fn conect_map(&self) -> BTreeMap<i32, Vec<i32>> {
        self.inner.conect_map.clone()
    }
    #[getter]
    fn has_d_fraction(&self) -> bool {
        self.inner.has_d_fraction
    }
    #[getter]
    fn non_ascii_line(&self) -> i32 {
        self.inner.non_ascii_line
    }
    #[getter]
    fn ter_status(&self) -> u8 {
        self.inner.ter_status
    }
    #[getter]
    fn has_origx(&self) -> bool {
        self.inner.has_origx
    }
    #[getter]
    fn origx(&self) -> crate::canonical_bio_binding::BioTransform {
        crate::canonical_bio_binding::BioTransform {
            inner: self.inner.origx,
        }
    }
    #[getter]
    fn info(&self) -> BTreeMap<String, String> {
        self.inner.info.clone()
    }
    #[getter]
    fn raw_remarks(&self) -> Vec<String> {
        self.inner.raw_remarks.clone()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioNcsOperator {
    pub(crate) inner: ck::BioNcsOperator,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioNcsOperator {
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    #[getter]
    fn given(&self) -> bool {
        self.inner.given
    }
    #[getter]
    fn transform(&self) -> crate::canonical_bio_binding::BioTransform {
        crate::canonical_bio_binding::BioTransform {
            inner: self.inner.transform,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioAssembly {
    pub(crate) inner: ck::BioAssembly,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioAssembly {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn author_determined(&self) -> bool {
        self.inner.author_determined
    }
    #[getter]
    fn software_determined(&self) -> bool {
        self.inner.software_determined
    }
    #[getter]
    fn special_kind(&self) -> BioAssemblySpecialKind {
        self.inner.special_kind.into()
    }
    #[getter]
    fn oligomeric_count(&self) -> i32 {
        self.inner.oligomeric_count
    }
    #[getter]
    fn oligomeric_details(&self) -> &str {
        &self.inner.oligomeric_details
    }
    #[getter]
    fn software_name(&self) -> &str {
        &self.inner.software_name
    }
    #[getter]
    fn buried_surface_area(&self) -> f64 {
        self.inner.buried_surface_area
    }
    #[getter]
    fn surface_area(&self) -> f64 {
        self.inner.surface_area
    }
    #[getter]
    fn solvent_free_energy_change(&self) -> f64 {
        self.inner.solvent_free_energy_change
    }
    #[getter]
    fn generators(&self) -> Vec<BioAssemblyGenerator> {
        self.inner
            .generators
            .iter()
            .cloned()
            .map(|inner| BioAssemblyGenerator { inner })
            .collect()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioAssemblyGenerator {
    pub(crate) inner: ck::BioAssemblyGenerator,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioAssemblyGenerator {
    #[getter]
    fn chains(&self) -> Vec<String> {
        self.inner.chains.clone()
    }
    #[getter]
    fn subchains(&self) -> Vec<String> {
        self.inner.subchains.clone()
    }
    #[getter]
    fn operators(&self) -> Vec<BioAssemblyOperator> {
        self.inner
            .operators
            .iter()
            .cloned()
            .map(|inner| BioAssemblyOperator { inner })
            .collect()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioAssemblyOperator {
    pub(crate) inner: ck::BioAssemblyOperator,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioAssemblyOperator {
    #[getter]
    fn name(&self) -> Option<String> {
        self.inner.name.clone()
    }
    #[getter]
    fn operator_type(&self) -> Option<String> {
        self.inner.operator_type.clone()
    }
    #[getter]
    fn transform(&self) -> crate::canonical_bio_binding::BioTransform {
        crate::canonical_bio_binding::BioTransform {
            inner: self.inner.transform,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioStrand {
    pub(crate) inner: ck::BioStrand,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioStrand {
    #[getter]
    fn start(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.start.clone(),
        }
    }
    #[getter]
    fn end(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.end.clone(),
        }
    }
    #[getter]
    fn hbond_atom2(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.hbond_atom2.clone(),
        }
    }
    #[getter]
    fn hbond_atom1(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.hbond_atom1.clone(),
        }
    }
    #[getter]
    fn sense(&self) -> i32 {
        self.inner.sense
    }
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioSoftwareItem {
    pub(crate) inner: ck::BioSoftwareItem,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioSoftwareItem {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn version(&self) -> &str {
        &self.inner.version
    }
    #[getter]
    fn date(&self) -> &str {
        &self.inner.date
    }
    #[getter]
    fn description(&self) -> &str {
        &self.inner.description
    }
    #[getter]
    fn contact_author(&self) -> &str {
        &self.inner.contact_author
    }
    #[getter]
    fn contact_author_email(&self) -> &str {
        &self.inner.contact_author_email
    }
    #[getter]
    fn classification(&self) -> BioSoftwareClassification {
        self.inner.classification.into()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioReflectionsInfo {
    pub(crate) inner: ck::BioReflectionsInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioReflectionsInfo {
    #[getter]
    fn resolution_high(&self) -> f64 {
        self.inner.resolution_high
    }
    #[getter]
    fn resolution_low(&self) -> f64 {
        self.inner.resolution_low
    }
    #[getter]
    fn completeness(&self) -> f64 {
        self.inner.completeness
    }
    #[getter]
    fn redundancy(&self) -> f64 {
        self.inner.redundancy
    }
    #[getter]
    fn r_merge(&self) -> f64 {
        self.inner.r_merge
    }
    #[getter]
    fn r_sym(&self) -> f64 {
        self.inner.r_sym
    }
    #[getter]
    fn mean_i_over_sigma(&self) -> f64 {
        self.inner.mean_i_over_sigma
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioBasicRefinementInfo {
    pub(crate) inner: ck::BioBasicRefinementInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioBasicRefinementInfo {
    #[getter]
    fn resolution_high(&self) -> f64 {
        self.inner.resolution_high
    }
    #[getter]
    fn resolution_low(&self) -> f64 {
        self.inner.resolution_low
    }
    #[getter]
    fn completeness(&self) -> f64 {
        self.inner.completeness
    }
    #[getter]
    fn reflection_count(&self) -> i32 {
        self.inner.reflection_count
    }
    #[getter]
    fn work_set_count(&self) -> i32 {
        self.inner.work_set_count
    }
    #[getter]
    fn rfree_set_count(&self) -> i32 {
        self.inner.rfree_set_count
    }
    #[getter]
    fn r_all(&self) -> f64 {
        self.inner.r_all
    }
    #[getter]
    fn r_work(&self) -> f64 {
        self.inner.r_work
    }
    #[getter]
    fn r_free(&self) -> f64 {
        self.inner.r_free
    }
    #[getter]
    fn cc_fo_fc_work(&self) -> f64 {
        self.inner.cc_fo_fc_work
    }
    #[getter]
    fn cc_fo_fc_free(&self) -> f64 {
        self.inner.cc_fo_fc_free
    }
    #[getter]
    fn fsc_work(&self) -> f64 {
        self.inner.fsc_work
    }
    #[getter]
    fn fsc_free(&self) -> f64 {
        self.inner.fsc_free
    }
    #[getter]
    fn cc_intensity_work(&self) -> f64 {
        self.inner.cc_intensity_work
    }
    #[getter]
    fn cc_intensity_free(&self) -> f64 {
        self.inner.cc_intensity_free
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioRefinementRestraint {
    pub(crate) inner: ck::BioRefinementRestraint,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioRefinementRestraint {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn count(&self) -> i32 {
        self.inner.count
    }
    #[getter]
    fn weight(&self) -> f64 {
        self.inner.weight
    }
    #[getter]
    fn function(&self) -> &str {
        &self.inner.function
    }
    #[getter]
    fn dev_ideal(&self) -> f64 {
        self.inner.dev_ideal
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioExperimentInfo {
    pub(crate) inner: ck::BioExperimentInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioExperimentInfo {
    #[getter]
    fn method(&self) -> &str {
        &self.inner.method
    }
    #[getter]
    fn number_of_crystals(&self) -> i32 {
        self.inner.number_of_crystals
    }
    #[getter]
    fn unique_reflections(&self) -> i32 {
        self.inner.unique_reflections
    }
    #[getter]
    fn reflections(&self) -> BioReflectionsInfo {
        BioReflectionsInfo {
            inner: self.inner.reflections.clone(),
        }
    }
    #[getter]
    fn b_wilson(&self) -> f64 {
        self.inner.b_wilson
    }
    #[getter]
    fn shells(&self) -> Vec<BioReflectionsInfo> {
        self.inner
            .shells
            .iter()
            .cloned()
            .map(|inner| BioReflectionsInfo { inner })
            .collect()
    }
    #[getter]
    fn diffraction_ids(&self) -> Vec<String> {
        self.inner.diffraction_ids.clone()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioDiffractionInfo {
    pub(crate) inner: ck::BioDiffractionInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioDiffractionInfo {
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    #[getter]
    fn temperature(&self) -> f64 {
        self.inner.temperature
    }
    #[getter]
    fn source(&self) -> &str {
        &self.inner.source
    }
    #[getter]
    fn source_type(&self) -> &str {
        &self.inner.source_type
    }
    #[getter]
    fn synchrotron(&self) -> &str {
        &self.inner.synchrotron
    }
    #[getter]
    fn beamline(&self) -> &str {
        &self.inner.beamline
    }
    #[getter]
    fn wavelengths(&self) -> &str {
        &self.inner.wavelengths
    }
    #[getter]
    fn scattering_type(&self) -> &str {
        &self.inner.scattering_type
    }
    #[getter]
    fn mono_or_laue(&self) -> u8 {
        self.inner.mono_or_laue
    }
    #[getter]
    fn monochromator(&self) -> &str {
        &self.inner.monochromator
    }
    #[getter]
    fn collection_date(&self) -> &str {
        &self.inner.collection_date
    }
    #[getter]
    fn optics(&self) -> &str {
        &self.inner.optics
    }
    #[getter]
    fn detector(&self) -> &str {
        &self.inner.detector
    }
    #[getter]
    fn detector_make(&self) -> &str {
        &self.inner.detector_make
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioExperimentalCrystalInfo {
    pub(crate) inner: ck::BioExperimentalCrystalInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioExperimentalCrystalInfo {
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    #[getter]
    fn description(&self) -> &str {
        &self.inner.description
    }
    #[getter]
    fn ph(&self) -> f64 {
        self.inner.ph
    }
    #[getter]
    fn ph_range(&self) -> &str {
        &self.inner.ph_range
    }
    #[getter]
    fn diffractions(&self) -> Vec<BioDiffractionInfo> {
        self.inner
            .diffractions
            .iter()
            .cloned()
            .map(|inner| BioDiffractionInfo { inner })
            .collect()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioTlsSelection {
    pub(crate) inner: ck::BioTlsSelection,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioTlsSelection {
    #[getter]
    fn chain(&self) -> String {
        self.inner.chain.as_str().to_owned()
    }
    #[getter]
    fn res_begin(&self) -> crate::canonical_bio_binding::PdbSeqId {
        crate::canonical_bio_binding::PdbSeqId {
            inner: self.inner.res_begin,
        }
    }
    #[getter]
    fn res_end(&self) -> crate::canonical_bio_binding::PdbSeqId {
        crate::canonical_bio_binding::PdbSeqId {
            inner: self.inner.res_end,
        }
    }
    #[getter]
    fn details(&self) -> &str {
        &self.inner.details
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioTlsGroup {
    pub(crate) inner: ck::BioTlsGroup,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioTlsGroup {
    #[getter]
    fn num_id(&self) -> i16 {
        self.inner.num_id
    }
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    #[getter]
    fn selections(&self) -> Vec<BioTlsSelection> {
        self.inner
            .selections
            .iter()
            .cloned()
            .map(|inner| BioTlsSelection { inner })
            .collect()
    }
    #[getter]
    fn origin(&self) -> [f64; 3] {
        self.inner.origin
    }
    #[getter]
    fn t(&self) -> [f64; 6] {
        self.inner.t
    }
    #[getter]
    fn l(&self) -> [f64; 6] {
        self.inner.l
    }
    #[getter]
    fn s(&self) -> [[f64; 3]; 3] {
        self.inner.s
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioRefinementInfo {
    pub(crate) inner: ck::BioRefinementInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioRefinementInfo {
    #[getter]
    fn basic(&self) -> BioBasicRefinementInfo {
        BioBasicRefinementInfo {
            inner: self.inner.basic.clone(),
        }
    }
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    #[getter]
    fn cross_validation_method(&self) -> &str {
        &self.inner.cross_validation_method
    }
    #[getter]
    fn rfree_selection_method(&self) -> &str {
        &self.inner.rfree_selection_method
    }
    #[getter]
    fn bin_count(&self) -> i32 {
        self.inner.bin_count
    }
    #[getter]
    fn bins(&self) -> Vec<BioBasicRefinementInfo> {
        self.inner
            .bins
            .iter()
            .cloned()
            .map(|inner| BioBasicRefinementInfo { inner })
            .collect()
    }
    #[getter]
    fn mean_b(&self) -> f64 {
        self.inner.mean_b
    }
    #[getter]
    fn aniso_b(&self) -> [f64; 6] {
        self.inner.aniso_b
    }
    #[getter]
    fn luzzati_error(&self) -> f64 {
        self.inner.luzzati_error
    }
    #[getter]
    fn dpi_blow_r(&self) -> f64 {
        self.inner.dpi_blow_r
    }
    #[getter]
    fn dpi_blow_rfree(&self) -> f64 {
        self.inner.dpi_blow_rfree
    }
    #[getter]
    fn dpi_cruickshank_r(&self) -> f64 {
        self.inner.dpi_cruickshank_r
    }
    #[getter]
    fn dpi_cruickshank_rfree(&self) -> f64 {
        self.inner.dpi_cruickshank_rfree
    }
    #[getter]
    fn restr_stats(&self) -> Vec<BioRefinementRestraint> {
        self.inner
            .restr_stats
            .iter()
            .cloned()
            .map(|inner| BioRefinementRestraint { inner })
            .collect()
    }
    #[getter]
    fn tls_groups(&self) -> Vec<BioTlsGroup> {
        self.inner
            .tls_groups
            .iter()
            .cloned()
            .map(|inner| BioTlsGroup { inner })
            .collect()
    }
    #[getter]
    fn remarks(&self) -> &str {
        &self.inner.remarks
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct AtomAddress {
    pub(crate) inner: ck::AtomAddress,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomAddress {
    #[getter]
    fn chain_name(&self) -> String {
        self.inner.chain_name().as_str().to_owned()
    }
    #[getter]
    fn residue(&self) -> ResidueAddress {
        ResidueAddress {
            inner: self.inner.residue(),
        }
    }
    #[getter]
    fn logical_atom_name(&self) -> &str {
        self.inner.logical_atom_name()
    }
    #[getter]
    fn altloc(&self) -> u8 {
        self.inner.altloc()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct ResidueAddress {
    pub(crate) inner: ck::ResidueAddress,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ResidueAddress {
    #[getter]
    fn sequence_number(&self) -> Option<i32> {
        self.inner.sequence_number()
    }
    #[getter]
    fn insertion_code(&self) -> Option<u8> {
        self.inner.insertion_code()
    }
    #[getter]
    fn segment(&self) -> &str {
        self.inner.segment()
    }
    #[getter]
    fn name(&self) -> String {
        self.inner.name().as_str().to_owned()
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    // String-valued enums preserve the existing Python string comparisons
    // while exposing the facade's closed, typed vocabulary.
    let enum_factory = module.py().import("enum")?.getattr("Enum")?;
    for (name, variants) in [
        (
            "BioCalcFlag",
            &["NotSet", "NoHydrogen", "Determined", "Calculated", "Dummy"][..],
        ),
        (
            "EntityKind",
            &["Unknown", "Polymer", "NonPolymer", "Branched", "Water"][..],
        ),
        (
            "PolymerKind",
            &[
                "Unknown",
                "PeptideL",
                "PeptideD",
                "Dna",
                "Rna",
                "DnaRnaHybrid",
                "SaccharideD",
                "SaccharideL",
                "Pna",
                "CyclicPseudoPeptide",
                "Other",
            ][..],
        ),
        (
            "ResidueKind",
            &[
                "AminoAcid",
                "Dna",
                "Rna",
                "Saccharide",
                "Water",
                "Buffer",
                "Ligand",
                "Unknown",
            ][..],
        ),
        (
            "ChainKind",
            &[
                "Protein",
                "Dna",
                "Rna",
                "ProteinDnaComplex",
                "ProteinRnaComplex",
                "LigandOnly",
                "WaterOnly",
                "Mixed",
                "Unknown",
            ][..],
        ),
    ] {
        let members = pyo3::types::PyDict::new(module.py());
        for variant in variants {
            members.set_item(*variant, *variant)?;
        }
        let options = pyo3::types::PyDict::new(module.py());
        options.set_item("type", module.py().get_type::<pyo3::types::PyString>())?;
        options.set_item("module", "cosmolkit")?;
        let kind = enum_factory.call((name, members), Some(&options))?;
        kind.setattr(
            "__str__",
            module
                .py()
                .get_type::<pyo3::types::PyString>()
                .getattr("__str__")?,
        )?;
        module.add(name, kind)?;
    }
    module.add_class::<BioConnection>()?;
    module.add_class::<BioCisPep>()?;
    module.add_class::<BioModRes>()?;
    module.add_class::<BioHelix>()?;
    module.add_class::<BioSheet>()?;
    module.add_class::<BioMetadata>()?;
    module.add_class::<BioStructureSourceState>()?;
    module.add_class::<BioNcsOperator>()?;
    module.add_class::<BioAssembly>()?;
    module.add_class::<BioAssemblyGenerator>()?;
    module.add_class::<BioAssemblyOperator>()?;
    module.add_class::<BioStrand>()?;
    module.add_class::<BioSoftwareItem>()?;
    module.add_class::<BioReflectionsInfo>()?;
    module.add_class::<BioBasicRefinementInfo>()?;
    module.add_class::<BioRefinementRestraint>()?;
    module.add_class::<BioExperimentInfo>()?;
    module.add_class::<BioDiffractionInfo>()?;
    module.add_class::<BioExperimentalCrystalInfo>()?;
    module.add_class::<BioTlsSelection>()?;
    module.add_class::<BioTlsGroup>()?;
    module.add_class::<BioRefinementInfo>()?;
    module.add_class::<BioConnectionKind>()?;
    module.add_class::<BioAsu>()?;
    module.add_class::<BioHelixClass>()?;
    module.add_class::<BioSoftwareClassification>()?;
    module.add_class::<BioAssemblySpecialKind>()?;
    module.add_class::<AtomAddress>()?;
    module.add_class::<ResidueAddress>()?;
    Ok(())
}
