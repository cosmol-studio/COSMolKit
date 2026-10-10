//! Complete readonly projections of the facade's detached BIO metadata.
//! Scalars retain their bits, bytes and missing-value sentinels. Collections
//! retain source order and duplicates; no parsing or chemistry occurs here.

use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
use std::collections::BTreeMap;

/// Relationship kind retained from Gemmi's ``_struct_conn`` model.
///
/// Declared values: ``Covale``, ``Disulf``, ``Hydrog``, ``MetalC``, ``Unknown``.
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

/// Unit-cell ASU relation for a structural connection.
///
/// Declared values: ``Same``, ``Different``, ``Any``.
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

/// PDB helix class values used by Gemmi's ``Helix`` metadata.
///
/// Declared values: ``UnknownHelix``, ``RAlpha``, ``ROmega``, ``RPi``, ``RGamma``, ``R310``, ``LAlpha``, ``LOmega``, ``LGamma``, ``Helix27``, ``HelixPolyProlineNone``.
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

/// Classification of the software entry recorded in structural metadata.
///
/// The explicit representation and declaration order preserve Gemmi's
/// ``SoftwareItem::Classification`` integer values.
///
/// Declared values: ``DataCollection``, ``DataExtraction``, ``DataProcessing``, ``DataReduction``, ``DataScaling``, ``ModelBuilding``, ``Phasing``, ``Refinement``, ``Unspecified``.
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

/// Special assembly classification retained from structural metadata.
///
/// Declared values: ``NotApplicable``, ``CompleteIcosahedral``, ``RepresentativeHelical``, ``CompletePoint``.
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

/// Typed ``_struct_conn`` relationship with source atom addresses.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioConnection {
    pub(crate) inner: ck::BioConnection,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioConnection {
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    /// Connection link identifier retained from structural metadata.
    #[getter]
    fn link_id(&self) -> &str {
        &self.inner.link_id
    }
    /// Return the structural connection type: covalent, disulfide, hydrogen, metal or unknown.
    #[getter]
    fn kind(&self) -> BioConnectionKind {
        self.inner.kind.into()
    }
    /// Whether the structural connection is within one asymmetric unit or between symmetry-related units.
    #[getter]
    fn asu(&self) -> BioAsu {
        self.inner.asu.into()
    }
    /// Structural address of the first connection partner.
    #[getter]
    fn partner1(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner1.clone(),
        }
    }
    /// Structural address of the second connection partner.
    #[getter]
    fn partner2(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner2.clone(),
        }
    }
    /// Connection distance in angstroms as reported by the source structure.
    #[getter]
    fn reported_distance(&self) -> f64 {
        self.inner.reported_distance
    }
    /// Source symmetry annotation associated with the connection.
    #[getter]
    fn reported_sym(&self) -> [i16; 4] {
        self.inner.reported_sym
    }
}

/// Cis-peptide relationship retained from Gemmi's PDBx/mmCIF metadata.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioCisPep {
    pub(crate) inner: ck::BioCisPep,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioCisPep {
    /// Structural address of the peptide carbonyl-carbon partner.
    #[getter]
    fn partner_c(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner_c.clone(),
        }
    }
    /// Structural address of the peptide nitrogen partner.
    #[getter]
    fn partner_n(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.partner_n.clone(),
        }
    }
    /// Source model number associated with this annotation.
    #[getter]
    fn model_num(&self) -> i32 {
        self.inner.model_num
    }
    /// Alternate-location byte restricting this annotation, or the source unspecified sentinel.
    #[getter]
    fn only_altloc(&self) -> u8 {
        self.inner.only_altloc
    }
    /// Cis-peptide angle in degrees reported by the source structure.
    #[getter]
    fn reported_angle(&self) -> f64 {
        self.inner.reported_angle
    }
}

/// Modified-residue metadata retained from Gemmi's ``ModRes`` value.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioModRes {
    pub(crate) inner: ck::BioModRes,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioModRes {
    /// Chain name recorded by the source structural address.
    #[getter]
    fn chain_name(&self) -> String {
        self.inner.chain_name.as_str().to_owned()
    }
    /// Structural residue address of the modified residue.
    #[getter]
    fn res_id(&self) -> ResidueAddress {
        ResidueAddress {
            inner: self.inner.res_id.clone(),
        }
    }
    /// Chemical component identifier of the unmodified parent residue.
    #[getter]
    fn parent_comp_id(&self) -> &str {
        &self.inner.parent_comp_id
    }
    /// Identifier of the residue modification.
    #[getter]
    fn mod_id(&self) -> &str {
        &self.inner.mod_id
    }
    /// Source-provided annotation detail text.
    #[getter]
    fn details(&self) -> &str {
        &self.inner.details
    }
}

/// Helix endpoints and reported classification from PDB/mmCIF metadata.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioHelix {
    pub(crate) inner: ck::BioHelix,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioHelix {
    /// Return the structural atom address marking the helix beginning.
    #[getter]
    fn start(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.start.clone(),
        }
    }
    /// Return the structural atom address marking the helix end.
    #[getter]
    fn end(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.end.clone(),
        }
    }
    /// PDB helix-class designation.
    #[getter]
    fn pdb_helix_class(&self) -> BioHelixClass {
        self.inner.pdb_helix_class.into()
    }
    /// Return the source-reported helix length.
    #[getter]
    fn length(&self) -> i32 {
        self.inner.length
    }
}

/// A sheet containing strands in source order.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioSheet {
    pub(crate) inner: ck::BioSheet,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioSheet {
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    /// Sheet strand annotations in source order.
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

/// Source-ordered metadata attached to a biological structure.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioMetadata {
    pub(crate) inner: ck::BioMetadata,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioMetadata {
    /// Structure authors in stored order.
    #[getter]
    fn authors(&self) -> Vec<String> {
        self.inner.authors.clone()
    }
    /// Experimental method and reflection metadata records.
    #[getter]
    fn experiments(&self) -> Vec<BioExperimentInfo> {
        self.inner
            .experiments
            .iter()
            .cloned()
            .map(|inner| BioExperimentInfo { inner })
            .collect()
    }
    /// Experimental crystal and diffraction records.
    #[getter]
    fn crystals(&self) -> Vec<BioExperimentalCrystalInfo> {
        self.inner
            .crystals
            .iter()
            .cloned()
            .map(|inner| BioExperimentalCrystalInfo { inner })
            .collect()
    }
    /// Structure refinement metadata records.
    #[getter]
    fn refinement(&self) -> Vec<BioRefinementInfo> {
        self.inner
            .refinement
            .iter()
            .cloned()
            .map(|inner| BioRefinementInfo { inner })
            .collect()
    }
    /// Software metadata records, or the writer option selecting that category.
    #[getter]
    fn software(&self) -> Vec<BioSoftwareItem> {
        self.inner
            .software
            .iter()
            .cloned()
            .map(|inner| BioSoftwareItem { inner })
            .collect()
    }
    /// Method/software reported for solving the structure.
    #[getter]
    fn solved_by(&self) -> &str {
        &self.inner.solved_by
    }
    /// Starting model description reported in structural metadata.
    #[getter]
    fn starting_model(&self) -> &str {
        &self.inner.starting_model
    }
    /// Biological assembly detail retained from PDB REMARK 300.
    #[getter]
    fn remark_300_detail(&self) -> &str {
        &self.inner.remark_300_detail
    }
}

/// Source-defined structure state that is not part of the atom hierarchy.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioStructureSourceState {
    pub(crate) inner: ck::BioStructureSourceState,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioStructureSourceState {
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    /// Experimental resolution in angstroms as recorded in the structure.
    #[getter]
    fn resolution(&self) -> f64 {
        self.inner.resolution
    }
    /// PDB serial-to-serial connectivity; duplicate targets encode repeated
    /// connectivity entries and remain distinct in each sorted target list.
    #[getter]
    fn conect_map(&self) -> BTreeMap<i32, Vec<i32>> {
        self.inner.conect_map.clone()
    }
    /// Whether the source includes deuterium-fraction annotations.
    #[getter]
    fn has_d_fraction(&self) -> bool {
        self.inner.has_d_fraction
    }
    /// Source line associated with non-ASCII input detection, when recorded.
    #[getter]
    fn non_ascii_line(&self) -> i32 {
        self.inner.non_ascii_line
    }
    /// Source TER-record handling state.
    #[getter]
    fn ter_status(&self) -> u8 {
        self.inner.ter_status
    }
    /// Whether the structure contains an original-coordinate transformation.
    #[getter]
    fn has_origx(&self) -> bool {
        self.inner.has_origx
    }
    /// Original-coordinate transformation stored in the structure.
    #[getter]
    fn origx(&self) -> crate::canonical_bio_binding::BioTransform {
        crate::canonical_bio_binding::BioTransform {
            inner: self.inner.origx,
        }
    }
    /// Return the associated residue dictionary information, including its classification and sequence codes.
    #[getter]
    fn info(&self) -> BTreeMap<String, String> {
        self.inner.info.clone()
    }
    /// Original PDB REMARK records, in source order.
    #[getter]
    fn raw_remarks(&self) -> Vec<String> {
        self.inner.raw_remarks.clone()
    }
}

/// Noncrystallographic-symmetry operator with identity, supplied/generated flag and 3D transform.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioNcsOperator {
    pub(crate) inner: ck::BioNcsOperator,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioNcsOperator {
    /// Identifier of this non-crystallographic symmetry operator in the structural metadata.
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    /// Whether the NCS operator was explicitly supplied by the source structure.
    #[getter]
    fn given(&self) -> bool {
        self.inner.given
    }
    /// Stored spatial transformation associated with this result or structural operator.
    #[getter]
    fn transform(&self) -> crate::canonical_bio_binding::BioTransform {
        crate::canonical_bio_binding::BioTransform {
            inner: self.inner.transform,
        }
    }
}

/// Biological assembly definition with provenance, oligomeric classification and transformation generators.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioAssembly {
    pub(crate) inner: ck::BioAssembly,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioAssembly {
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    /// Whether the biological assembly was identified by the structure authors.
    #[getter]
    fn author_determined(&self) -> bool {
        self.inner.author_determined
    }
    /// Whether the biological assembly was predicted by software.
    #[getter]
    fn software_determined(&self) -> bool {
        self.inner.software_determined
    }
    /// Special biological-assembly classification.
    #[getter]
    fn special_kind(&self) -> BioAssemblySpecialKind {
        self.inner.special_kind.into()
    }
    /// Reported number of subunits in the biological assembly.
    #[getter]
    fn oligomeric_count(&self) -> i32 {
        self.inner.oligomeric_count
    }
    /// Source description of the assembly oligomeric state.
    #[getter]
    fn oligomeric_details(&self) -> &str {
        &self.inner.oligomeric_details
    }
    /// Name of the software that determined the assembly.
    #[getter]
    fn software_name(&self) -> &str {
        &self.inner.software_name
    }
    /// Reported buried surface area of the assembly.
    #[getter]
    fn buried_surface_area(&self) -> f64 {
        self.inner.buried_surface_area
    }
    /// Reported assembly surface area.
    #[getter]
    fn surface_area(&self) -> f64 {
        self.inner.surface_area
    }
    /// Reported solvent free-energy change for assembly formation.
    #[getter]
    fn solvent_free_energy_change(&self) -> f64 {
        self.inner.solvent_free_energy_change
    }
    /// Assembly generators selecting subchains and symmetry operators.
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

/// Biological assembly generator selecting chains/subchains and operator combinations.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioAssemblyGenerator {
    pub(crate) inner: ck::BioAssemblyGenerator,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioAssemblyGenerator {
    /// Return chain rows/references in stored hierarchy order.
    #[getter]
    fn chains(&self) -> Vec<String> {
        self.inner.chains.clone()
    }
    /// Source subchain identifiers associated with the entity.
    #[getter]
    fn subchains(&self) -> Vec<String> {
        self.inner.subchains.clone()
    }
    /// Spatial operators used by this assembly generator.
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

/// Named biological-assembly spatial operator.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioAssemblyOperator {
    pub(crate) inner: ck::BioAssemblyOperator,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioAssemblyOperator {
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> Option<String> {
        self.inner.name.clone()
    }
    /// Source classification of the assembly operator.
    #[getter]
    fn operator_type(&self) -> Option<String> {
        self.inner.operator_type.clone()
    }
    /// Stored spatial transformation associated with this result or structural operator.
    #[getter]
    fn transform(&self) -> crate::canonical_bio_binding::BioTransform {
        crate::canonical_bio_binding::BioTransform {
            inner: self.inner.transform,
        }
    }
}

/// One ordered strand and its source endpoint/hydrogen-bond addresses.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioStrand {
    pub(crate) inner: ck::BioStrand,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioStrand {
    /// Return the structural atom address marking the strand beginning.
    #[getter]
    fn start(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.start.clone(),
        }
    }
    /// Return the structural atom address marking the strand end.
    #[getter]
    fn end(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.end.clone(),
        }
    }
    /// Structural address of the second sheet hydrogen-bond atom.
    #[getter]
    fn hbond_atom2(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.hbond_atom2.clone(),
        }
    }
    /// Structural address of the first sheet hydrogen-bond atom.
    #[getter]
    fn hbond_atom1(&self) -> AtomAddress {
        AtomAddress {
            inner: self.inner.hbond_atom1.clone(),
        }
    }
    /// Strand orientation relative to the previous strand as encoded by the source.
    #[getter]
    fn sense(&self) -> i32 {
        self.inner.sense
    }
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
}

/// One source-ordered software record from structural metadata.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioSoftwareItem {
    pub(crate) inner: ck::BioSoftwareItem,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioSoftwareItem {
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    /// Return the version identifier recorded by this value or implementation.
    #[getter]
    fn version(&self) -> &str {
        &self.inner.version
    }
    /// Date recorded in the source metadata.
    #[getter]
    fn date(&self) -> &str {
        &self.inner.date
    }
    /// Source-provided description of this metadata record.
    #[getter]
    fn description(&self) -> &str {
        &self.inner.description
    }
    /// Contact author recorded for this software entry.
    #[getter]
    fn contact_author(&self) -> &str {
        &self.inner.contact_author
    }
    /// Contact-author email recorded for this software entry.
    #[getter]
    fn contact_author_email(&self) -> &str {
        &self.inner.contact_author_email
    }
    /// Role of the software in collection, processing, model building or refinement.
    #[getter]
    fn classification(&self) -> BioSoftwareClassification {
        self.inner.classification.into()
    }
}

/// Reflection statistics retained from diffraction metadata.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioReflectionsInfo {
    pub(crate) inner: ck::BioReflectionsInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioReflectionsInfo {
    /// High-resolution limit in angstroms reported for this dataset/bin.
    #[getter]
    fn resolution_high(&self) -> f64 {
        self.inner.resolution_high
    }
    /// Low-resolution limit in angstroms reported for this dataset/bin.
    #[getter]
    fn resolution_low(&self) -> f64 {
        self.inner.resolution_low
    }
    /// Reflection completeness reported by the source; source units/sentinels are preserved.
    #[getter]
    fn completeness(&self) -> f64 {
        self.inner.completeness
    }
    /// Reported average reflection multiplicity.
    #[getter]
    fn redundancy(&self) -> f64 {
        self.inner.redundancy
    }
    /// Reported merging R statistic.
    #[getter]
    fn r_merge(&self) -> f64 {
        self.inner.r_merge
    }
    /// Reported symmetry-related R statistic.
    #[getter]
    fn r_sym(&self) -> f64 {
        self.inner.r_sym
    }
    /// Reported mean intensity divided by its uncertainty.
    #[getter]
    fn mean_i_over_sigma(&self) -> f64 {
        self.inner.mean_i_over_sigma
    }
}

/// Refinement statistics shared by overall refinement and per-bin values.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioBasicRefinementInfo {
    pub(crate) inner: ck::BioBasicRefinementInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioBasicRefinementInfo {
    /// High-resolution limit in angstroms reported for this dataset/bin.
    #[getter]
    fn resolution_high(&self) -> f64 {
        self.inner.resolution_high
    }
    /// Low-resolution limit in angstroms reported for this dataset/bin.
    #[getter]
    fn resolution_low(&self) -> f64 {
        self.inner.resolution_low
    }
    /// Reflection completeness reported by the source; source units/sentinels are preserved.
    #[getter]
    fn completeness(&self) -> f64 {
        self.inner.completeness
    }
    /// Number of reflections used in this refinement record.
    #[getter]
    fn reflection_count(&self) -> i32 {
        self.inner.reflection_count
    }
    /// Number of reflections in the refinement working set.
    #[getter]
    fn work_set_count(&self) -> i32 {
        self.inner.work_set_count
    }
    /// Number of reflections in the R-free validation set.
    #[getter]
    fn rfree_set_count(&self) -> i32 {
        self.inner.rfree_set_count
    }
    /// Reported R factor over all reflections.
    #[getter]
    fn r_all(&self) -> f64 {
        self.inner.r_all
    }
    /// Reported working-set R factor.
    #[getter]
    fn r_work(&self) -> f64 {
        self.inner.r_work
    }
    /// Reported free-set R factor.
    #[getter]
    fn r_free(&self) -> f64 {
        self.inner.r_free
    }
    /// Observed/calculated structure-factor correlation for the working set.
    #[getter]
    fn cc_fo_fc_work(&self) -> f64 {
        self.inner.cc_fo_fc_work
    }
    /// Observed/calculated structure-factor correlation for the free set.
    #[getter]
    fn cc_fo_fc_free(&self) -> f64 {
        self.inner.cc_fo_fc_free
    }
    /// Reported Fourier-shell correlation for the working set.
    #[getter]
    fn fsc_work(&self) -> f64 {
        self.inner.fsc_work
    }
    /// Reported Fourier-shell correlation for the free set.
    #[getter]
    fn fsc_free(&self) -> f64 {
        self.inner.fsc_free
    }
    /// Reported intensity correlation for the working set.
    #[getter]
    fn cc_intensity_work(&self) -> f64 {
        self.inner.cc_intensity_work
    }
    /// Reported intensity correlation for the free set.
    #[getter]
    fn cc_intensity_free(&self) -> f64 {
        self.inner.cc_intensity_free
    }
}

/// One named restraint statistic in refinement metadata.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioRefinementRestraint {
    pub(crate) inner: ck::BioRefinementRestraint,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioRefinementRestraint {
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    /// Number of restraints represented by this refinement statistic.
    #[getter]
    fn count(&self) -> i32 {
        self.inner.count
    }
    /// Residue molecular weight from the residue dictionary.
    #[getter]
    fn weight(&self) -> f64 {
        self.inner.weight
    }
    /// Name of the restraint function or failing chemical function.
    #[getter]
    fn function(&self) -> &str {
        &self.inner.function
    }
    /// Reported deviation from ideal restraint geometry.
    #[getter]
    fn dev_ideal(&self) -> f64 {
        self.inner.dev_ideal
    }
}

/// Experimental method and its associated reflection statistics.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioExperimentInfo {
    pub(crate) inner: ck::BioExperimentInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioExperimentInfo {
    /// Experimental method, or canonical operation name in operation metadata.
    #[getter]
    fn method(&self) -> &str {
        &self.inner.method
    }
    /// Number of crystals contributing to the experiment.
    #[getter]
    fn number_of_crystals(&self) -> i32 {
        self.inner.number_of_crystals
    }
    /// Reported number of unique reflections.
    #[getter]
    fn unique_reflections(&self) -> i32 {
        self.inner.unique_reflections
    }
    /// Overall reflection statistics associated with the experiment.
    #[getter]
    fn reflections(&self) -> BioReflectionsInfo {
        BioReflectionsInfo {
            inner: self.inner.reflections.clone(),
        }
    }
    /// Reported Wilson B factor.
    #[getter]
    fn b_wilson(&self) -> f64 {
        self.inner.b_wilson
    }
    /// Resolution-shell reflection statistics in stored order.
    #[getter]
    fn shells(&self) -> Vec<BioReflectionsInfo> {
        self.inner
            .shells
            .iter()
            .cloned()
            .map(|inner| BioReflectionsInfo { inner })
            .collect()
    }
    /// Source identifiers of diffraction datasets associated with the experiment.
    #[getter]
    fn diffraction_ids(&self) -> Vec<String> {
        self.inner.diffraction_ids.clone()
    }
}

/// Source-ordered diffraction metadata for an experiment.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioDiffractionInfo {
    pub(crate) inner: ck::BioDiffractionInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioDiffractionInfo {
    /// Identifier of this diffraction dataset in the structural metadata.
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    /// Diffraction collection temperature in kelvin, as recorded.
    #[getter]
    fn temperature(&self) -> f64 {
        self.inner.temperature
    }
    /// Source-format identifiers retained separately from local row identifiers.
    #[getter]
    fn source(&self) -> &str {
        &self.inner.source
    }
    /// Classification of the diffraction radiation source.
    #[getter]
    fn source_type(&self) -> &str {
        &self.inner.source_type
    }
    /// Synchrotron facility name retained from source metadata.
    #[getter]
    fn synchrotron(&self) -> &str {
        &self.inner.synchrotron
    }
    /// Beamline used for diffraction collection.
    #[getter]
    fn beamline(&self) -> &str {
        &self.inner.beamline
    }
    /// Radiation wavelengths in angstroms reported for the diffraction dataset.
    #[getter]
    fn wavelengths(&self) -> &str {
        &self.inner.wavelengths
    }
    /// Scattering/radiation type used by the diffraction experiment.
    #[getter]
    fn scattering_type(&self) -> &str {
        &self.inner.scattering_type
    }
    /// The source's single-byte monochromatic/Laue code; NUL means unset.
    #[getter]
    fn mono_or_laue(&self) -> u8 {
        self.inner.mono_or_laue
    }
    /// Monochromator description recorded for the experiment.
    #[getter]
    fn monochromator(&self) -> &str {
        &self.inner.monochromator
    }
    /// Date of diffraction data collection.
    #[getter]
    fn collection_date(&self) -> &str {
        &self.inner.collection_date
    }
    /// Source description of diffraction beam optics.
    #[getter]
    fn optics(&self) -> &str {
        &self.inner.optics
    }
    /// Detector type recorded for the diffraction experiment.
    #[getter]
    fn detector(&self) -> &str {
        &self.inner.detector
    }
    /// Detector manufacturer/model recorded for the experiment.
    #[getter]
    fn detector_make(&self) -> &str {
        &self.inner.detector_make
    }
}

/// Experimental crystal metadata, distinct from lattice/coordinate crystal state.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioExperimentalCrystalInfo {
    pub(crate) inner: ck::BioExperimentalCrystalInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioExperimentalCrystalInfo {
    /// Identifier of this experimental crystal in the structural metadata.
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    /// Source-provided description of this metadata record.
    #[getter]
    fn description(&self) -> &str {
        &self.inner.description
    }
    /// Crystallization pH reported in the source metadata.
    #[getter]
    fn ph(&self) -> f64 {
        self.inner.ph
    }
    /// Reported crystallization pH range.
    #[getter]
    fn ph_range(&self) -> &str {
        &self.inner.ph_range
    }
    /// Diffraction dataset metadata associated with this crystal.
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

/// Residue range selected for TLS refinement, with chain name, sequence endpoints and insertion codes.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioTlsSelection {
    pub(crate) inner: ck::BioTlsSelection,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioTlsSelection {
    /// Return the chain reference for a local chain identifier.
    #[getter]
    fn chain(&self) -> String {
        self.inner.chain.as_str().to_owned()
    }
    /// Beginning residue address of the TLS selection.
    #[getter]
    fn res_begin(&self) -> crate::canonical_bio_binding::PdbSeqId {
        crate::canonical_bio_binding::PdbSeqId {
            inner: self.inner.res_begin,
        }
    }
    /// Ending residue address of the TLS selection.
    #[getter]
    fn res_end(&self) -> crate::canonical_bio_binding::PdbSeqId {
        crate::canonical_bio_binding::PdbSeqId {
            inner: self.inner.res_end,
        }
    }
    /// Source-provided annotation detail text.
    #[getter]
    fn details(&self) -> &str {
        &self.inner.details
    }
}

/// Tensor, origin, and ordered selections for one TLS group.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioTlsGroup {
    pub(crate) inner: ck::BioTlsGroup,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioTlsGroup {
    /// Source's small numeric ID; ``-1`` means no numeric ID was assigned.
    #[getter]
    fn num_id(&self) -> i16 {
        self.inner.num_id
    }
    /// Identifier of this TLS refinement group in the structural metadata.
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    /// Residue selections included in this TLS group.
    #[getter]
    fn selections(&self) -> Vec<BioTlsSelection> {
        self.inner
            .selections
            .iter()
            .cloned()
            .map(|inner| BioTlsSelection { inner })
            .collect()
    }
    /// Cartesian origin of the TLS group in angstroms.
    #[getter]
    fn origin(&self) -> [f64; 3] {
        self.inner.origin
    }
    /// Symmetric components in source order: u11, u22, u33, u12, u13, u23.
    #[getter]
    fn t(&self) -> [f64; 6] {
        self.inner.t
    }
    /// Symmetric components in source order: u11, u22, u33, u12, u13, u23.
    #[getter]
    fn l(&self) -> [f64; 6] {
        self.inner.l
    }
    /// Full row-major 3x3 S matrix.
    #[getter]
    fn s(&self) -> [[f64; 3]; 3] {
        self.inner.s
    }
}

/// Overall refinement metadata, including the source-defined per-bin and
/// restraint records.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct BioRefinementInfo {
    pub(crate) inner: ck::BioRefinementInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BioRefinementInfo {
    /// The inherited overall refinement statistics.
    #[getter]
    fn basic(&self) -> BioBasicRefinementInfo {
        BioBasicRefinementInfo {
            inner: self.inner.basic.clone(),
        }
    }
    /// Identifier of this refinement record in the structural metadata.
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    /// Cross-validation method reported for refinement.
    #[getter]
    fn cross_validation_method(&self) -> &str {
        &self.inner.cross_validation_method
    }
    /// Method used to choose the R-free reflection set.
    #[getter]
    fn rfree_selection_method(&self) -> &str {
        &self.inner.rfree_selection_method
    }
    /// Number of refinement resolution bins.
    #[getter]
    fn bin_count(&self) -> i32 {
        self.inner.bin_count
    }
    /// Per-resolution-bin refinement statistics in stored order.
    #[getter]
    fn bins(&self) -> Vec<BioBasicRefinementInfo> {
        self.inner
            .bins
            .iter()
            .cloned()
            .map(|inner| BioBasicRefinementInfo { inner })
            .collect()
    }
    /// Reported mean atomic B factor in square angstroms.
    #[getter]
    fn mean_b(&self) -> f64 {
        self.inner.mean_b
    }
    /// Symmetric components in Gemmi's source order: u11, u22, u33, u12, u13, u23.
    #[getter]
    fn aniso_b(&self) -> [f64; 6] {
        self.inner.aniso_b
    }
    /// Reported Luzzati coordinate error estimate.
    #[getter]
    fn luzzati_error(&self) -> f64 {
        self.inner.luzzati_error
    }
    /// Blow diffraction precision index calculated from the R factor.
    #[getter]
    fn dpi_blow_r(&self) -> f64 {
        self.inner.dpi_blow_r
    }
    /// Blow diffraction precision index calculated from R-free.
    #[getter]
    fn dpi_blow_rfree(&self) -> f64 {
        self.inner.dpi_blow_rfree
    }
    /// Cruickshank diffraction precision index calculated from the R factor.
    #[getter]
    fn dpi_cruickshank_r(&self) -> f64 {
        self.inner.dpi_cruickshank_r
    }
    /// Cruickshank diffraction precision index calculated from R-free.
    #[getter]
    fn dpi_cruickshank_rfree(&self) -> f64 {
        self.inner.dpi_cruickshank_rfree
    }
    /// Refinement restraint statistics in stored order.
    #[getter]
    fn restr_stats(&self) -> Vec<BioRefinementRestraint> {
        self.inner
            .restr_stats
            .iter()
            .cloned()
            .map(|inner| BioRefinementRestraint { inner })
            .collect()
    }
    /// TLS refinement groups in stored order.
    #[getter]
    fn tls_groups(&self) -> Vec<BioTlsGroup> {
        self.inner
            .tls_groups
            .iter()
            .cloned()
            .map(|inner| BioTlsGroup { inner })
            .collect()
    }
    /// Source refinement remarks.
    #[getter]
    fn remarks(&self) -> &str {
        &self.inner.remarks
    }
}

/// Logical structural atom address identified by chain, residue, atom name and alternate-location label.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct AtomAddress {
    pub(crate) inner: ck::AtomAddress,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomAddress {
    /// Chain name recorded by the source structural address.
    #[getter]
    fn chain_name(&self) -> String {
        self.inner.chain_name().as_str().to_owned()
    }
    /// Return the parent residue reference.
    #[getter]
    fn residue(&self) -> ResidueAddress {
        ResidueAddress {
            inner: self.inner.residue(),
        }
    }
    /// Atom name used for structural-address matching, excluding storage padding.
    #[getter]
    fn logical_atom_name(&self) -> &str {
        self.inner.logical_atom_name()
    }
    /// Alternate-location label represented as an integer character code.
    #[getter]
    fn altloc(&self) -> u8 {
        self.inner.altloc()
    }
}

/// Structural residue address carrying sequence number, insertion code, segment and residue name; not a local row index.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub(crate) struct ResidueAddress {
    pub(crate) inner: ck::ResidueAddress,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ResidueAddress {
    /// Source residue sequence number, not a local row identifier.
    #[getter]
    fn sequence_number(&self) -> Option<i32> {
        self.inner.sequence_number()
    }
    /// Source residue insertion code.
    #[getter]
    fn insertion_code(&self) -> Option<u8> {
        self.inner.insertion_code()
    }
    /// Source segment identifier, when present.
    #[getter]
    fn segment(&self) -> &str {
        self.inner.segment()
    }
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> String {
        self.inner.name().as_str().to_owned()
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    // String-valued enums preserve the existing Python string comparisons
    // while exposing the facade's closed, typed vocabulary.
    let enum_factory = module.py().import("enum")?.getattr("Enum")?;
    for (name, variants, description) in [
        (
            "BioCalcFlag",
            &["NotSet", "NoHydrogen", "Determined", "Calculated", "Dummy"][..],
            "Origin/classification of a structural atom site, including determined, calculated and dummy positions.",
        ),
        (
            "EntityKind",
            &["Unknown", "Polymer", "NonPolymer", "Branched", "Water"][..],
            "Biological entity classification: polymer, non-polymer, branched, water or unknown.",
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
            "Polymer chemistry classification, distinguishing peptide handedness, DNA, RNA, saccharides and other polymer types.",
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
            "Structural residue classification, including amino acids, nucleotides, saccharides, water, buffers and ligands.",
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
            "Structural chain classification based on its residue composition, including protein, nucleic acids, ligands, water and mixed chains.",
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
        kind.setattr("__doc__", description)?;
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
