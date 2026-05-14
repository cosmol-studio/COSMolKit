//! Biomolecular structure primitives.
//!
//! `BioStructure` is a flat-row, hierarchy-indexed representation for proteins,
//! DNA, RNA, and complexes. It is NOT a giant `Molecule`; it is a hierarchy +
//! coordinate + assembly object. See `dev/BioStructureOperationContractDesign.md`.
//!
//! Gemmi marker convention is defined in `dev/source_reproduction_protocol.md`.

use std::marker::PhantomData;

pub mod invariants;
pub mod ops;

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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

macro_rules! gemmi_residue_classification {
    (
        amino_acids: [$($amino_acid:literal),+ $(,)?],
        waters: [$($water:literal),+ $(,)?] $(,)?
    ) => {
        #[cfg(test)]
        const GEMMI_AMINO_ACID_RESIDUE_NAMES: &[&str] = &[$($amino_acid),+];
        #[cfg(test)]
        const GEMMI_WATER_RESIDUE_NAMES: &[&str] = &[$($water),+];

        #[must_use]
        pub fn classify_residue_name(name: ResidueName) -> ResidueKind {
            match name.as_str() {
                $($amino_acid)|+ => ResidueKind::AminoAcid,
                $($water)|+ => ResidueKind::Water,
                _ => ResidueKind::Unknown,
            }
        }
    };
}

// BEGIN GEMMI CPP TABLE gemmi::get_residue_info RI::AA/RI::AAD/RI::PAA/RI::MAA
// Gemmi✔️✔️: ResidueInfo::is_amino_acid returns true for AA, AAD, PAA, and MAA.
// Gemmi✔️✔️: The names below are the complete unique residue-name set from
// Gemmi✔️✔️: `third_party/gemmi/src/resinfo.cpp` entries with these four kinds.
// Gemmi❌❌: one-letter code, linking type, hydrogen count, molecular weight,
// Gemmi❌❌: and standard-vs-modified residue semantics are not modeled here.
// END GEMMI CPP TABLE gemmi::get_residue_info RI::AA/RI::AAD/RI::PAA/RI::MAA

// BEGIN GEMMI CPP FUNCTION gemmi::Residue::is_water
// Gemmi✔️✔️: return id == ialpha4_id("HOH") || id == ialpha4_id("DOD") ||
// Gemmi✔️✔️:        id == ialpha4_id("WAT") || id == ialpha4_id("H2O");
// END GEMMI CPP FUNCTION gemmi::Residue::is_water
gemmi_residue_classification! {
    amino_acids: [
        "0AF", "0TD", "3FG", "ABA", "AGM", "AIB", "ALA", "ALC", "ALY", "ARG", "ASN", "ASP",
        "ASX", "B3E", "BFD", "BMT", "CAF", "CAS", "CGU", "CIR", "CME", "CR2", "CR8", "CRF",
        "CRO", "CRQ", "CSD", "CSH", "CSO", "CSS", "CSX", "CXM", "CYS", "DAB", "DAL", "DAR",
        "DAS", "DCY", "DGL", "DGN", "DHA", "DHI", "DIL", "DLE", "DLY", "DPN", "DPR", "DSG",
        "DSN", "DTH", "DTR", "DTY", "DVA", "FGA", "FME", "FVA", "GHP", "GL3", "GLN", "GLU",
        "GLX", "GLY", "GYS", "HIC", "HIS", "HYP", "IAS", "ILE", "KCX", "KPI", "LEU", "LLP",
        "LYS", "M3L", "MAA", "MDO", "MEA", "MED", "MEN", "MEQ", "MET", "MHO", "MHS", "MK8",
        "MLE", "MLU", "MLY", "MLZ", "MSE", "MVA", "NEP", "NLE", "NRQ", "OAS", "OCS", "OMY",
        "OMZ", "ORN", "PCA", "PHD", "PHE", "PHI", "PHL", "PRO", "PTR", "PYL", "SAC", "SAR",
        "SCH", "SCY", "SEC", "SEP", "SER", "SMC", "SME", "SNC", "SNN", "THR", "TOX", "TPO",
        "TPQ", "TRP", "TRQ", "TYR", "TYS", "UNK", "VAL", "YCM",
    ],
    waters: ["HOH", "DOD", "WAT", "H2O"],
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
        assert_eq!(GEMMI_AMINO_ACID_RESIDUE_NAMES.len(), 128);
        for name in GEMMI_AMINO_ACID_RESIDUE_NAMES {
            assert_eq!(
                classify_residue_name(residue_name(name)),
                ResidueKind::AminoAcid
            );
        }
    }

    #[test]
    fn classifies_gemmi_water_names_without_guessing_other_residues() {
        assert_eq!(GEMMI_WATER_RESIDUE_NAMES.len(), 4);
        for name in GEMMI_WATER_RESIDUE_NAMES {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ResidueSourceIds {
    pub seq_id: Option<PdbSeqId>,
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
    pub source: AtomSourceIds,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ResidueRow {
    pub chain_id: ChainId,
    pub atom_span: RowSpan<AtomId>,
    pub name: ResidueName,
    pub kind: ResidueKind,
    pub source: ResidueSourceIds,
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
    pub sequence: Vec<String>,
    pub subchains: Vec<PdbChainId>,
    pub source: EntitySourceIds,
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
    pub(crate) atoms: Vec<AtomRow>,
    pub(crate) residues: Vec<ResidueRow>,
    pub(crate) chains: Vec<ChainRow>,
    pub(crate) entities: Vec<EntityRow>,
    pub(crate) models: Vec<ModelRow>,
    pub(crate) coordinates: CoordinateBlock,
    pub(crate) metadata: BioMetadata,
    pub(crate) crystal: Option<CrystalInfo>,
}

impl BioStructure {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
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
