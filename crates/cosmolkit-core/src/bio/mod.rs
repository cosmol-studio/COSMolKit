//! Biomolecular structure primitives.
//!
//! `BioStructure` is a flat-row, hierarchy-indexed representation for proteins,
//! DNA, RNA, and complexes. It is NOT a giant `Molecule`; it is a hierarchy +
//! coordinate + assembly object. See `dev/bio_structure_operation_contract_design.md`.
//!
//! Gemmi marker convention is defined in `dev/source_reproduction_protocol.md`.

use std::collections::HashMap;
use std::marker::PhantomData;

pub mod invariants;
pub mod ops;
pub mod protein;
pub mod resinfo;

// ---------------------------------------------------------------------------
// Stable row IDs
// ---------------------------------------------------------------------------

macro_rules! row_id {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
            #[must_use]
            pub const fn new(index: u32) -> Self {
                Self(index)
            }

            #[must_use]
            pub const fn index(self) -> u32 {
                self.0
            }
        }
    };
}

row_id!(AtomId);
row_id!(ResidueId);
row_id!(ChainId);
row_id!(EntityId);
row_id!(ModelId);
row_id!(BondId);
row_id!(AssemblyId);
row_id!(AltLocGroupId);

// ---------------------------------------------------------------------------
// RowSpan — contiguous child range within a parent block
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RowSpan<T> {
    pub start: u32,
    pub len: u32,
    _marker: PhantomData<T>,
}

impl<T> RowSpan<T> {
    #[must_use]
    pub const fn new(start: u32, len: u32) -> Self {
        Self {
            start,
            len,
            _marker: PhantomData,
        }
    }

    #[must_use]
    pub const fn end(self) -> u32 {
        self.start + self.len
    }

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.len == 0
    }
}

// ---------------------------------------------------------------------------
// Polymer classification enums
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResidueKind {
    AminoAcid,
    DNA,
    RNA,
    Saccharide,
    Water,
    Ligand,
    Ion,
    Cofactor,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PolymerKind {
    Peptide,
    DNA,
    RNA,
    PeptideLike,
    NucleicAcidHybrid,
    Saccharide,
    NonPolymer,
    Water,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EntityKind {
    Polymer,
    NonPolymer,
    Branched,
    Water,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChainKind {
    Protein,
    DNA,
    RNA,
    ProteinDNAComplex,
    ProteinRNAComplex,
    LigandOnly,
    WaterOnly,
    Mixed,
    Unknown,
}

// ---------------------------------------------------------------------------
// Source identifier types (PDB/mmCIF provenance, NOT row ids)
// ---------------------------------------------------------------------------

/// PDB atom serial number (source provenance only).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PdbAtomSerial(pub i32);

/// PDB/mmCIF chain identifier string (up to 4 chars).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PdbChainId(pub [u8; 4], pub u8);

impl PdbChainId {
    #[must_use]
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.0[..self.1 as usize]).unwrap_or("")
    }
}

/// PDB residue sequence number + insertion code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PdbSeqId {
    pub seq_num: i32,
    pub ins_code: Option<u8>,
}

// ---------------------------------------------------------------------------
// Atom / residue / chain / model name types
// ---------------------------------------------------------------------------

/// Up to 4-char atom name (e.g. " CA ", " N  ").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomName(pub [u8; 4]);

/// Up to 3-char residue name (e.g. "ALA", "GLY").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ResidueName(pub [u8; 4], pub u8);

impl ResidueName {
    #[must_use]
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.0[..self.1 as usize]).unwrap_or("")
    }
}

#[must_use]
pub fn classify_residue_name(name: ResidueName) -> ResidueKind {
    let info = resinfo::find_tabulated_residue(name.as_str());
    if info.is_amino_acid() {
        ResidueKind::AminoAcid
    } else if info.is_water() {
        ResidueKind::Water
    } else {
        ResidueKind::Unknown
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn residue_name(value: &str) -> ResidueName {
        let mut bytes = [0; 4];
        bytes[..value.len()].copy_from_slice(value.as_bytes());
        ResidueName(bytes, value.len() as u8)
    }

    #[test]
    fn classifies_complete_gemmi_amino_acid_vocabulary() {
        let amino_acids = resinfo::RESIDUE_INFO_TABLE
            .iter()
            .filter(|info| info.is_amino_acid())
            .collect::<Vec<_>>();
        assert_eq!(amino_acids.len(), 128);
        for info in amino_acids {
            assert_eq!(
                classify_residue_name(residue_name(info.name)),
                ResidueKind::AminoAcid
            );
        }
    }

    #[test]
    fn classifies_gemmi_water_names_without_guessing_other_residues() {
        for name in ["HOH", "DOD", "WAT", "H2O"] {
            assert_eq!(
                classify_residue_name(residue_name(name)),
                ResidueKind::Water
            );
        }
        assert_eq!(
            classify_residue_name(residue_name("XYZ")),
            ResidueKind::Unknown
        );
    }
}

/// Single-char altloc label (e.g. b'A', b'B').
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AltLocLabel(pub u8);

// ---------------------------------------------------------------------------
// Source identifier bundles (per row)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomSourceIds {
    pub serial: Option<PdbAtomSerial>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BioCalcFlag {
    #[default]
    NotSet,
    NoHydrogen,
    Determined,
    Calculated,
    Dummy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ResidueSourceIds {
    pub seq_id: Option<PdbSeqId>,
    pub label_seq_id: Option<i32>,
    pub segment_id: Option<[u8; 4]>,
    pub subchain_id: Option<PdbChainId>,
    pub label_entity_id: Option<EntityId>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChainSourceIds {
    pub auth_chain_id: Option<PdbChainId>,
    pub label_asym_id: Option<PdbChainId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EntitySourceIds {
    pub source_entity_id: String,
}

// ---------------------------------------------------------------------------
// Flat row types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct AtomRow {
    pub residue_id: ResidueId,
    pub name: AtomName,
    pub element: crate::Element,
    pub altloc: Option<AltLocLabel>,
    pub occupancy: Option<f32>,
    pub b_iso: Option<f32>,
    pub formal_charge: Option<i8>,
    pub anisou: Option<[f32; 6]>,
    pub calc_flag: BioCalcFlag,
    pub tls_group_id: Option<i16>,
    pub fraction: Option<f32>,
    pub source: AtomSourceIds,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ResidueRow {
    pub chain_id: ChainId,
    pub atom_span: RowSpan<AtomId>,
    pub name: ResidueName,
    pub kind: ResidueKind,
    pub entity_kind: EntityKind,
    pub het_flag: Option<char>,
    pub source: ResidueSourceIds,
    pub sifts_unp: Option<BioSiftsUnpResidue>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChainRow {
    pub model_id: ModelId,
    pub entity_id: Option<EntityId>,
    pub residue_span: RowSpan<ResidueId>,
    pub kind: ChainKind,
    pub source: ChainSourceIds,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EntityRow {
    pub kind: EntityKind,
    pub polymer_kind: PolymerKind,
    pub reflects_microhetero: bool,
    pub sequence: Vec<String>,
    pub dbrefs: Vec<BioEntityDbRef>,
    pub sifts_unp_acc: Vec<String>,
    pub subchains: Vec<PdbChainId>,
    pub source: EntitySourceIds,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BioEntityDbRef {
    pub db_name: String,
    pub accession_code: String,
    pub id_code: String,
    pub isoform: String,
    pub seq_begin: PdbSeqId,
    pub seq_end: PdbSeqId,
    pub db_begin: PdbSeqId,
    pub db_end: PdbSeqId,
    pub label_seq_begin: Option<i32>,
    pub label_seq_end: Option<i32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BioSiftsUnpResidue {
    pub res: Option<char>,
    pub acc_index: u8,
    pub num: u16,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BioModRes {
    pub chain_name: String,
    pub res_id: PdbSeqId,
    pub parent_comp_id: String,
    pub mod_id: String,
    pub details: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ModelRow {
    pub chain_span: RowSpan<ChainId>,
    pub source_model_number: Option<i32>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioMetadata {
    pub entry_id: Option<String>,
    pub title: Option<String>,
    pub pdbx_keywords: Option<String>,
    pub keywords: Option<String>,
    pub experimental_method: Option<String>,
    pub received_initial_deposition_date: Option<String>,
    pub authors: Vec<String>,
    pub software: Vec<BioSoftwareItem>,
    pub refinement: Vec<BioRefinementInfo>,
    pub experiments: Vec<BioExperimentInfo>,
    pub experiment_crystals: Vec<BioExperimentCrystalInfo>,
    pub solved_by: Option<String>,
    pub starting_model: Option<String>,
    pub remark_300_detail: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BioSoftwareClassification {
    DataCollection,
    DataExtraction,
    DataProcessing,
    DataReduction,
    DataScaling,
    ModelBuilding,
    Phasing,
    Refinement,
    #[default]
    Unspecified,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioSoftwareItem {
    pub name: String,
    pub version: String,
    pub date: String,
    pub description: String,
    pub contact_author: String,
    pub contact_author_email: String,
    pub classification: BioSoftwareClassification,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioRefinementBin {
    pub resolution_high: Option<f64>,
    pub resolution_low: Option<f64>,
    pub completeness: Option<f64>,
    pub reflection_count: Option<i32>,
    pub work_set_count: Option<i32>,
    pub rfree_set_count: Option<i32>,
    pub r_all: Option<f64>,
    pub r_work: Option<f64>,
    pub r_free: Option<f64>,
    pub cc_fo_fc_work: Option<f64>,
    pub cc_fo_fc_free: Option<f64>,
    pub fsc_work: Option<f64>,
    pub fsc_free: Option<f64>,
    pub cc_intensity_work: Option<f64>,
    pub cc_intensity_free: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioRefinementRestraint {
    pub name: String,
    pub count: Option<i32>,
    pub weight: Option<f64>,
    pub function: String,
    pub dev_ideal: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioTlsSelection {
    pub chain: String,
    pub res_begin: String,
    pub res_end: String,
    pub details: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioTlsGroup {
    pub num_id: Option<i16>,
    pub id: String,
    pub selections: Vec<BioTlsSelection>,
    pub origin: [f64; 3],
    pub t: [[f64; 3]; 3],
    pub l: [[f64; 3]; 3],
    pub s: [[f64; 3]; 3],
}

impl Default for BioTlsGroup {
    fn default() -> Self {
        Self {
            num_id: None,
            id: String::new(),
            selections: Vec::new(),
            origin: [f64::NAN; 3],
            t: [[f64::NAN; 3]; 3],
            l: [[f64::NAN; 3]; 3],
            s: [[f64::NAN; 3]; 3],
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioAnisotropicB {
    pub u11: f64,
    pub u22: f64,
    pub u33: f64,
    pub u12: f64,
    pub u13: f64,
    pub u23: f64,
}

impl Default for BioAnisotropicB {
    fn default() -> Self {
        Self {
            u11: f64::NAN,
            u22: f64::NAN,
            u33: f64::NAN,
            u12: f64::NAN,
            u13: f64::NAN,
            u23: f64::NAN,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioRefinementInfo {
    pub id: String,
    pub resolution_high: Option<f64>,
    pub resolution_low: Option<f64>,
    pub completeness: Option<f64>,
    pub reflection_count: Option<i32>,
    pub work_set_count: Option<i32>,
    pub rfree_set_count: Option<i32>,
    pub r_all: Option<f64>,
    pub r_work: Option<f64>,
    pub r_free: Option<f64>,
    pub cross_validation_method: String,
    pub rfree_selection_method: String,
    pub bin_count: Option<i32>,
    pub bins: Vec<BioRefinementBin>,
    pub mean_b: Option<f64>,
    pub aniso_b: BioAnisotropicB,
    pub luzzati_error: Option<f64>,
    pub dpi_blow_r: Option<f64>,
    pub dpi_blow_rfree: Option<f64>,
    pub dpi_cruickshank_r: Option<f64>,
    pub dpi_cruickshank_rfree: Option<f64>,
    pub cc_fo_fc_work: Option<f64>,
    pub cc_fo_fc_free: Option<f64>,
    pub fsc_work: Option<f64>,
    pub fsc_free: Option<f64>,
    pub cc_intensity_work: Option<f64>,
    pub cc_intensity_free: Option<f64>,
    pub restr_stats: Vec<BioRefinementRestraint>,
    pub tls_groups: Vec<BioTlsGroup>,
    pub remarks: String,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioReflectionsInfo {
    pub resolution_high: Option<f64>,
    pub resolution_low: Option<f64>,
    pub completeness: Option<f64>,
    pub redundancy: Option<f64>,
    pub r_merge: Option<f64>,
    pub r_sym: Option<f64>,
    pub mean_i_over_sigma: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioExperimentInfo {
    pub method: String,
    pub number_of_crystals: Option<i32>,
    pub unique_reflections: Option<i32>,
    pub diffraction_ids: Vec<String>,
    pub reflections: BioReflectionsInfo,
    pub b_wilson: Option<f64>,
    pub shells: Vec<BioReflectionsInfo>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioDiffractionInfo {
    pub id: String,
    pub collection_date: String,
    pub temperature: Option<f64>,
    pub source: String,
    pub source_type: String,
    pub synchrotron: String,
    pub beamline: String,
    pub wavelengths: String,
    pub scattering_type: String,
    pub monochromator: String,
    pub optics: String,
    pub detector: String,
    pub detector_make: String,
    pub mono_or_laue: Option<char>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioExperimentCrystalInfo {
    pub id: String,
    pub description: String,
    pub ph: Option<f64>,
    pub ph_range: String,
    pub diffractions: Vec<BioDiffractionInfo>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BioAsu {
    #[default]
    Any,
    Same,
    Different,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BioAtomAddress {
    pub chain_name: String,
    pub seq_id: Option<PdbSeqId>,
    pub atom_name: String,
    pub altloc: Option<AltLocLabel>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BioConnectionType {
    #[default]
    Disulf,
    Covale,
    MetalC,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioConnection {
    pub name: String,
    pub type_: BioConnectionType,
    pub partner1: BioAtomAddress,
    pub partner2: BioAtomAddress,
    pub asu: BioAsu,
    pub reported_sym: [i16; 4],
    pub reported_distance: Option<f32>,
    pub link_id: String,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioCisPep {
    pub partner_c: BioAtomAddress,
    pub partner_n: BioAtomAddress,
    pub model_num: i32,
    pub only_altloc: Option<AltLocLabel>,
    pub reported_angle: Option<f32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BioHelixClass {
    #[default]
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

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioHelix {
    pub start: BioAtomAddress,
    pub end: BioAtomAddress,
    pub helix_class: BioHelixClass,
    pub length: i32,
}

impl BioHelix {
    pub fn set_helix_class_as_int(&mut self, n: i32) {
        self.helix_class = match n {
            1 => BioHelixClass::RAlpha,
            2 => BioHelixClass::ROmega,
            3 => BioHelixClass::RPi,
            4 => BioHelixClass::RGamma,
            5 => BioHelixClass::R310,
            6 => BioHelixClass::LAlpha,
            7 => BioHelixClass::LOmega,
            8 => BioHelixClass::LGamma,
            9 => BioHelixClass::Helix27,
            10 => BioHelixClass::HelixPolyProlineNone,
            _ => BioHelixClass::UnknownHelix,
        };
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioSheetStrand {
    pub start: BioAtomAddress,
    pub end: BioAtomAddress,
    pub hbond_atom2: BioAtomAddress,
    pub hbond_atom1: BioAtomAddress,
    pub sense: i32,
    pub name: String,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioSheet {
    pub name: String,
    pub strands: Vec<BioSheetStrand>,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BioTransform {
    pub mat: [[f32; 3]; 3],
    pub vec: [f32; 3],
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioNcsOperator {
    pub id: String,
    pub given: bool,
    pub transform: BioTransform,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioAssemblyOperator {
    pub name: String,
    pub type_: String,
    pub transform: BioTransform,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioAssemblyGenerator {
    pub chains: Vec<String>,
    pub subchains: Vec<String>,
    pub operators: Vec<BioAssemblyOperator>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BioAssemblySpecialKind {
    #[default]
    NA,
    CompleteIcosahedral,
    RepresentativeHelical,
    CompletePoint,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioAssembly {
    pub name: String,
    pub author_determined: bool,
    pub software_determined: bool,
    pub special_kind: BioAssemblySpecialKind,
    pub oligomeric_count: i32,
    pub oligomeric_details: String,
    pub software_name: String,
    pub absa: Option<f64>,
    pub ssa: Option<f64>,
    pub more: Option<f64>,
    pub generators: Vec<BioAssemblyGenerator>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrystalCell {
    pub a: f32,
    pub b: f32,
    pub c: f32,
    pub alpha: f32,
    pub beta: f32,
    pub gamma: f32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CrystalInfo {
    pub cell: CrystalCell,
    pub spacegroup_hm: Option<String>,
    pub z_pdb: Option<String>,
    pub scale: Option<BioTransform>,
    pub frac: BioTransform,
    pub orth: BioTransform,
    pub explicit_matrices: bool,
    pub cs_count: i16,
    pub cell_images: Vec<BioTransform>,
}

// ---------------------------------------------------------------------------
// Coordinate block
// ---------------------------------------------------------------------------

/// 3D coordinates for all atoms, indexed by AtomId.
/// Invariant: `len() == atoms.len()` in the owning BioStructure.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CoordinateBlock {
    pub(crate) positions: Vec<[f32; 3]>,
}

impl CoordinateBlock {
    #[must_use]
    pub fn positions(&self) -> &[[f32; 3]] {
        &self.positions
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BioCoorFormat {
    Pdb,
    Mmcif,
    Mmjson,
    ChemComp,
    #[default]
    Unknown,
    Detect,
}

// ---------------------------------------------------------------------------
// Top-level BioStructure
// ---------------------------------------------------------------------------

/// Flat-row biomolecular structure.
///
/// Hierarchy: Model → Chain → Residue → Atom (all stored as flat Vecs).
/// Row ids are snapshot-local indices; PDB/mmCIF source identifiers are stored
/// in the `source` fields of each row.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioStructure {
    pub(crate) name: String,
    pub(crate) input_format: BioCoorFormat,
    pub(crate) atoms: Vec<AtomRow>,
    pub(crate) residues: Vec<ResidueRow>,
    pub(crate) chains: Vec<ChainRow>,
    pub(crate) entities: Vec<EntityRow>,
    pub(crate) models: Vec<ModelRow>,
    pub(crate) coordinates: CoordinateBlock,
    pub(crate) metadata: BioMetadata,
    pub(crate) crystal: Option<CrystalInfo>,
    pub(crate) resolution: Option<f64>,
    pub(crate) non_ascii_line: Option<usize>,
    pub(crate) raw_remarks: Vec<String>,
    pub(crate) ter_status: char,
    pub(crate) has_d_fraction: bool,
    pub(crate) mod_residues: Vec<BioModRes>,
    pub(crate) shortened_ccd_codes: Vec<(String, String)>,
    pub(crate) conect_map: HashMap<i32, Vec<i32>>,
    pub(crate) deferred_conn_records: Vec<String>,
    pub(crate) connections: Vec<BioConnection>,
    pub(crate) cispeps: Vec<BioCisPep>,
    pub(crate) helices: Vec<BioHelix>,
    pub(crate) sheets: Vec<BioSheet>,
    pub(crate) remark_290_operators: Vec<String>,
    pub(crate) assemblies: Vec<BioAssembly>,
    pub(crate) has_origx: bool,
    pub(crate) origx: BioTransform,
    pub(crate) ncs_operators: Vec<BioNcsOperator>,
    pub(crate) ncs_oper_identity_id: Option<String>,
}

impl BioStructure {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Reads a Gemmi-aligned PDB structural record stream into a `BioStructure`.
    pub fn from_pdb_str(text: &str) -> Result<Self, crate::io::bio::BioReadError> {
        Self::from_pdb_str_with_params(text, crate::io::bio::BioPdbReadParams::default())
    }

    /// Reads a Gemmi-aligned PDB structural record stream with explicit PDB reader parameters.
    pub fn from_pdb_str_with_params(
        text: &str,
        params: crate::io::bio::BioPdbReadParams,
    ) -> Result<Self, crate::io::bio::BioReadError> {
        crate::io::bio::read_pdb_bio_structure_from_str_with_params(text, params)
    }

    /// Reads a Gemmi-aligned mmCIF structural document into a `BioStructure`.
    pub fn from_mmcif_str(text: &str, path: &str) -> Result<Self, crate::io::bio::BioReadError> {
        Self::from_str_with_format(text, path, BioCoorFormat::Mmcif)
    }

    /// Reads a Gemmi-aligned structure by detecting the input format from text.
    pub fn from_structure_str(
        text: &str,
        path: &str,
    ) -> Result<Self, crate::io::bio::BioReadError> {
        Self::from_str_with_format(text, path, BioCoorFormat::Detect)
    }

    /// Reads a Gemmi-aligned structure using the requested coordinate format.
    pub fn from_str_with_format(
        text: &str,
        path: &str,
        format: BioCoorFormat,
    ) -> Result<Self, crate::io::bio::BioReadError> {
        crate::io::bio::read_structure_from_memory(text, path, format)
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    #[must_use]
    pub fn num_residues(&self) -> usize {
        self.residues.len()
    }

    #[must_use]
    pub fn num_chains(&self) -> usize {
        self.chains.len()
    }

    #[must_use]
    pub fn num_models(&self) -> usize {
        self.models.len()
    }

    #[must_use]
    pub fn num_entities(&self) -> usize {
        self.entities.len()
    }

    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub fn input_format(&self) -> BioCoorFormat {
        self.input_format
    }

    #[must_use]
    pub fn atoms(&self) -> &[AtomRow] {
        &self.atoms
    }

    #[must_use]
    pub fn residues(&self) -> &[ResidueRow] {
        &self.residues
    }

    #[must_use]
    pub fn chains(&self) -> &[ChainRow] {
        &self.chains
    }

    #[must_use]
    pub fn models(&self) -> &[ModelRow] {
        &self.models
    }

    #[must_use]
    pub fn entities(&self) -> &[EntityRow] {
        &self.entities
    }

    #[must_use]
    pub fn metadata(&self) -> &BioMetadata {
        &self.metadata
    }

    #[must_use]
    pub fn crystal(&self) -> Option<&CrystalInfo> {
        self.crystal.as_ref()
    }

    #[must_use]
    pub fn resolution(&self) -> Option<f64> {
        self.resolution
    }

    #[must_use]
    pub fn ter_status(&self) -> char {
        self.ter_status
    }

    #[must_use]
    pub fn connections(&self) -> &[BioConnection] {
        &self.connections
    }

    #[must_use]
    pub fn cispeps(&self) -> &[BioCisPep] {
        &self.cispeps
    }

    #[must_use]
    pub fn mod_residues(&self) -> &[BioModRes] {
        &self.mod_residues
    }

    #[must_use]
    pub fn assemblies(&self) -> &[BioAssembly] {
        &self.assemblies
    }

    #[must_use]
    pub fn has_origx(&self) -> bool {
        self.has_origx
    }

    #[must_use]
    pub fn origx(&self) -> &BioTransform {
        &self.origx
    }

    #[must_use]
    pub fn ncs_operators(&self) -> &[BioNcsOperator] {
        &self.ncs_operators
    }

    #[must_use]
    pub fn ncs_oper_identity_id(&self) -> Option<&str> {
        self.ncs_oper_identity_id.as_deref()
    }

    #[must_use]
    pub fn coordinates(&self) -> &CoordinateBlock {
        &self.coordinates
    }

    #[must_use]
    pub fn atom_position(&self, atom: AtomId) -> Option<[f32; 3]> {
        self.coordinates
            .positions
            .get(atom.index() as usize)
            .copied()
    }

    #[must_use]
    pub fn residue_atoms(&self, residue: ResidueId) -> Option<&[AtomRow]> {
        let row = self.residues.get(residue.index() as usize)?;
        let start = row.atom_span.start as usize;
        let end = row.atom_span.end() as usize;
        self.atoms.get(start..end)
    }
}
