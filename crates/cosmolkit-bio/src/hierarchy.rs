//! Validated detached macromolecular hierarchy values.

use std::{fmt, marker::PhantomData};

use cosmolkit_types::Element;

use crate::{
    AltLocLabel, AtomName, AtomSourceIds, ChainSourceIds, EntitySourceIds, PdbSeqId,
    ResidueInfoKind, ResidueName, ResidueSourceIds,
};

macro_rules! row_id {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
            #[must_use]
            pub const fn new(value: u32) -> Self {
                Self(value)
            }

            #[must_use]
            pub const fn index(self) -> usize {
                self.0 as usize
            }

            #[must_use]
            pub const fn value(self) -> u32 {
                self.0
            }
        }
    };
}

row_id!(BioAtomId);
row_id!(BioResidueId);
row_id!(BioChainId);
row_id!(BioEntityId);
row_id!(BioModelId);
row_id!(BioAssemblyId);
row_id!(BioAltLocGroupId);

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct BioRowSpan<I> {
    start: u32,
    len: u32,
    marker: PhantomData<fn() -> I>,
}

impl<I> Copy for BioRowSpan<I> {}

impl<I> Clone for BioRowSpan<I> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<I> BioRowSpan<I> {
    pub fn new(start: u32, len: u32) -> Result<Self, BioStructureError> {
        start
            .checked_add(len)
            .ok_or(BioStructureError::RowSpanOverflow { start, len })?;
        Ok(Self {
            start,
            len,
            marker: PhantomData,
        })
    }

    pub fn from_usize(start: usize, len: usize) -> Result<Self, BioStructureError> {
        let start = u32::try_from(start)
            .map_err(|_| BioStructureError::RowIndexTooLarge { value: start })?;
        let len =
            u32::try_from(len).map_err(|_| BioStructureError::RowIndexTooLarge { value: len })?;
        Self::new(start, len)
    }

    #[must_use]
    pub const fn start(self) -> u32 {
        self.start
    }

    #[must_use]
    pub const fn len(self) -> u32 {
        self.len
    }

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.len == 0
    }

    #[must_use]
    pub const fn end(self) -> u32 {
        // Construction proves this addition cannot overflow.
        self.start + self.len
    }

    pub fn slice<'a, T>(self, rows: &'a [T]) -> Result<&'a [T], BioStructureError> {
        rows.get(self.start as usize..self.end() as usize).ok_or(
            BioStructureError::RowSpanOutOfBounds {
                start: self.start,
                len: self.len,
                table_len: rows.len(),
            },
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BioCoordinateFormat {
    #[default]
    Unknown,
    Detect,
    Pdb,
    Mmcif,
    Mmjson,
    ChemComp,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BioCalcFlag {
    #[default]
    NotSet,
    NoHydrogen,
    Determined,
    Calculated,
    Dummy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum EntityKind {
    #[default]
    Unknown,
    Polymer,
    NonPolymer,
    Branched,
    Water,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum PolymerKind {
    #[default]
    Unknown,
    PeptideL,
    PeptideD,
    Dna,
    Rna,
    DnaRnaHybrid,
    SaccharideD,
    SaccharideL,
    Pna,
    CyclicPseudoPeptide,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ResidueKind {
    AminoAcid,
    Dna,
    Rna,
    Saccharide,
    Water,
    Buffer,
    Ligand,
    #[default]
    Unknown,
}

impl From<ResidueInfoKind> for ResidueKind {
    fn from(value: ResidueInfoKind) -> Self {
        match value {
            ResidueInfoKind::Aa
            | ResidueInfoKind::Aad
            | ResidueInfoKind::Paa
            | ResidueInfoKind::Maa => Self::AminoAcid,
            ResidueInfoKind::Dna => Self::Dna,
            ResidueInfoKind::Rna => Self::Rna,
            ResidueInfoKind::Hoh => Self::Water,
            ResidueInfoKind::Buf => Self::Buffer,
            ResidueInfoKind::Pyr | ResidueInfoKind::Ket => Self::Saccharide,
            ResidueInfoKind::Els => Self::Ligand,
            ResidueInfoKind::Unknown => Self::Unknown,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ChainKind {
    Protein,
    Dna,
    Rna,
    ProteinDnaComplex,
    ProteinRnaComplex,
    LigandOnly,
    WaterOnly,
    Mixed,
    #[default]
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AltLocRequest {
    Any,
    Exact(Option<AltLocLabel>),
}

#[must_use]
pub fn is_same_conformer(first: Option<AltLocLabel>, second: Option<AltLocLabel>) -> bool {
    // Gemmi✔️✔️: inline bool is_same_conformer(char altloc1, char altloc2) {
    // Gemmi✔️✔️:   return altloc1 == '\0' || altloc2 == '\0' || altloc1 == altloc2;
    // Gemmi✔️✔️: }
    // Behavior review: None is exactly the source NUL state and equality is byte-exact.
    // Complexity review: both implementations use a constant number of comparisons.
    first.is_none() || second.is_none() || first == second
}

#[must_use]
pub fn altloc_matches(stored: Option<AltLocLabel>, request: AltLocRequest) -> bool {
    // Gemmi✔️✔️: bool altloc_matches(char request) const {
    // Gemmi✔️✔️:   return request == '*' || altloc == '\0' || altloc == request;
    // Gemmi✔️✔️: }
    // Behavior review: wildcard is represented by the request enum, never as stored data.
    // Complexity review: both implementations use constant-time comparisons.
    match request {
        AltLocRequest::Any => true,
        AltLocRequest::Exact(requested) => stored.is_none() || stored == requested,
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioAtomRow {
    residue_id: BioResidueId,
    name: AtomName,
    element: Element,
    altloc: Option<AltLocLabel>,
    formal_charge: i8,
    calc_flag: BioCalcFlag,
    occupancy: f64,
    b_iso: f64,
    anisou: [f64; 6],
    tls_group_id: i16,
    fraction: f64,
    source: AtomSourceIds,
}

impl BioAtomRow {
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        residue_id: BioResidueId,
        name: AtomName,
        element: Element,
        altloc: Option<AltLocLabel>,
        formal_charge: i8,
        calc_flag: BioCalcFlag,
        occupancy: f64,
        b_iso: f64,
        anisou: [f64; 6],
        tls_group_id: i16,
        fraction: f64,
        source: AtomSourceIds,
    ) -> Self {
        // Gemmi✔️✔️: char altloc = '\0'; // 0 if not set
        // Gemmi✔️✔️: signed char charge = 0;  // [-8, +8]
        // Gemmi✔️✔️: CalcFlag calc_flag = CalcFlag::NotSet;  // mmCIF _atom_site.calc_flag
        // Gemmi✔️✔️: short tls_group_id = -1;
        // Gemmi✔️✔️: float fraction = 0.f;
        // Gemmi✔️✔️: float occ = 1.0f;
        // Gemmi✔️✔️: float b_iso = 20.0f; // arbitrary default value
        // Gemmi✔️✔️: SMat33<float> aniso = {0, 0, 0, 0, 0, 0};
        // Behavior review: IO supplies source values; no value-layer normalization is performed.
        // Complexity review: construction is constant-time with no allocation.
        Self {
            residue_id,
            name,
            element,
            altloc,
            formal_charge,
            calc_flag,
            occupancy,
            b_iso,
            anisou,
            tls_group_id,
            fraction,
            source,
        }
    }

    #[must_use]
    pub const fn residue_id(&self) -> BioResidueId {
        self.residue_id
    }
    #[must_use]
    pub const fn name(&self) -> AtomName {
        self.name
    }
    #[must_use]
    pub const fn element(&self) -> Element {
        self.element
    }
    #[must_use]
    pub const fn altloc(&self) -> Option<AltLocLabel> {
        self.altloc
    }
    #[must_use]
    pub const fn formal_charge(&self) -> i8 {
        self.formal_charge
    }
    #[must_use]
    pub const fn calc_flag(&self) -> BioCalcFlag {
        self.calc_flag
    }
    #[must_use]
    pub const fn occupancy(&self) -> f64 {
        self.occupancy
    }
    #[must_use]
    pub const fn b_iso(&self) -> f64 {
        self.b_iso
    }
    #[must_use]
    pub const fn anisou(&self) -> &[f64; 6] {
        &self.anisou
    }
    #[must_use]
    pub const fn tls_group_id(&self) -> i16 {
        self.tls_group_id
    }
    #[must_use]
    pub const fn fraction(&self) -> f64 {
        self.fraction
    }
    #[must_use]
    pub const fn source(&self) -> &AtomSourceIds {
        &self.source
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct BioSiftsUnpResidue {
    residue: Option<u8>,
    accession_index: u8,
    number: u16,
}

impl BioSiftsUnpResidue {
    #[must_use]
    pub const fn new(residue: Option<u8>, accession_index: u8, number: u16) -> Self {
        // Gemmi✔️✔️: char res = '\0';             // _pdbx_sifts_xref_db.unp_res
        // Gemmi✔️✔️: std::uint8_t acc_index = 0;  // index of Entity::sifts_unp_acc
        // Gemmi✔️✔️: std::uint16_t num = 0;       // _pdbx_sifts_xref_db.unp_num
        // Behavior review: None represents the source NUL sentinel without losing byte values.
        // Complexity review: construction is constant-time.
        Self {
            residue,
            accession_index,
            number,
        }
    }

    #[must_use]
    pub const fn residue(self) -> Option<u8> {
        self.residue
    }
    #[must_use]
    pub const fn accession_index(self) -> u8 {
        self.accession_index
    }
    #[must_use]
    pub const fn number(self) -> u16 {
        self.number
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioResidueRow {
    chain_id: BioChainId,
    atom_span: BioRowSpan<BioAtomId>,
    name: ResidueName,
    residue_info_kind: ResidueInfoKind,
    kind: ResidueKind,
    entity_kind: EntityKind,
    entity_id: Option<BioEntityId>,
    het_flag: Option<u8>,
    source: ResidueSourceIds,
    sifts_unp: BioSiftsUnpResidue,
}

impl BioResidueRow {
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        chain_id: BioChainId,
        atom_span: BioRowSpan<BioAtomId>,
        name: ResidueName,
        residue_info_kind: ResidueInfoKind,
        entity_kind: EntityKind,
        entity_id: Option<BioEntityId>,
        het_flag: Option<u8>,
        source: ResidueSourceIds,
        sifts_unp: BioSiftsUnpResidue,
    ) -> Self {
        Self {
            chain_id,
            atom_span,
            name,
            residue_info_kind,
            kind: residue_info_kind.into(),
            entity_kind,
            entity_id,
            het_flag,
            source,
            sifts_unp,
        }
    }

    #[must_use]
    pub const fn chain_id(&self) -> BioChainId {
        self.chain_id
    }
    #[must_use]
    pub const fn atom_span(&self) -> BioRowSpan<BioAtomId> {
        self.atom_span
    }
    #[must_use]
    pub const fn name(&self) -> ResidueName {
        self.name
    }
    #[must_use]
    pub const fn residue_info_kind(&self) -> ResidueInfoKind {
        self.residue_info_kind
    }
    #[must_use]
    pub const fn kind(&self) -> ResidueKind {
        self.kind
    }
    #[must_use]
    pub const fn entity_kind(&self) -> EntityKind {
        self.entity_kind
    }
    #[must_use]
    pub const fn entity_id(&self) -> Option<BioEntityId> {
        self.entity_id
    }
    #[must_use]
    pub const fn het_flag(&self) -> Option<u8> {
        self.het_flag
    }
    #[must_use]
    pub const fn source(&self) -> &ResidueSourceIds {
        &self.source
    }
    #[must_use]
    pub const fn sifts_unp(&self) -> BioSiftsUnpResidue {
        self.sifts_unp
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioChainRow {
    model_id: BioModelId,
    entity_id: Option<BioEntityId>,
    residue_span: BioRowSpan<BioResidueId>,
    kind: ChainKind,
    source: ChainSourceIds,
}

impl BioChainRow {
    #[must_use]
    pub fn new(
        model_id: BioModelId,
        entity_id: Option<BioEntityId>,
        residue_span: BioRowSpan<BioResidueId>,
        kind: ChainKind,
        source: ChainSourceIds,
    ) -> Self {
        Self {
            model_id,
            entity_id,
            residue_span,
            kind,
            source,
        }
    }

    #[must_use]
    pub const fn model_id(&self) -> BioModelId {
        self.model_id
    }
    #[must_use]
    pub const fn entity_id(&self) -> Option<BioEntityId> {
        self.entity_id
    }
    #[must_use]
    pub const fn residue_span(&self) -> BioRowSpan<BioResidueId> {
        self.residue_span
    }
    #[must_use]
    pub const fn kind(&self) -> ChainKind {
        self.kind
    }
    #[must_use]
    pub const fn source(&self) -> &ChainSourceIds {
        &self.source
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BioModelRow {
    chain_span: BioRowSpan<BioChainId>,
    source_model_number: Option<i32>,
}

impl BioModelRow {
    #[must_use]
    pub const fn new(chain_span: BioRowSpan<BioChainId>, source_model_number: Option<i32>) -> Self {
        Self {
            chain_span,
            source_model_number,
        }
    }

    #[must_use]
    pub const fn chain_span(self) -> BioRowSpan<BioChainId> {
        self.chain_span
    }
    #[must_use]
    pub const fn source_model_number(self) -> Option<i32> {
        self.source_model_number
    }
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BioEntityRow {
    kind: EntityKind,
    polymer_kind: PolymerKind,
    reflects_microhetero: bool,
    full_sequence: Vec<String>,
    dbrefs: Vec<BioEntityDbRef>,
    sifts_unp_accessions: Vec<String>,
    subchains: Vec<String>,
    source: EntitySourceIds,
}

impl BioEntityRow {
    #[must_use]
    pub fn new(
        kind: EntityKind,
        polymer_kind: PolymerKind,
        reflects_microhetero: bool,
        full_sequence: Vec<String>,
        dbrefs: Vec<BioEntityDbRef>,
        sifts_unp_accessions: Vec<String>,
        subchains: Vec<String>,
        source: EntitySourceIds,
    ) -> Self {
        // Gemmi✔️✔️: std::vector<std::string> subchains;
        // Gemmi✔️✔️: bool reflects_microhetero = false;
        // Gemmi✔️✔️: std::vector<DbRef> dbrefs;
        // Gemmi✔️✔️: std::vector<std::string> sifts_unp_acc;
        // Gemmi✔️✔️: std::vector<std::string> full_sequence;
        // Behavior review: source order, duplicates, and comma-bearing sequence rows are retained.
        // Complexity review: each input vector is moved without rescanning or cloning.
        Self {
            kind,
            polymer_kind,
            reflects_microhetero,
            full_sequence,
            dbrefs,
            sifts_unp_accessions,
            subchains,
            source,
        }
    }

    #[must_use]
    pub fn first_mon(mon_list: &str) -> &str {
        // Gemmi✔️✔️: static std::string first_mon(const std::string& mon_list) {
        // Gemmi✔️✔️:   return mon_list.substr(0, mon_list.find(','));
        // Gemmi✔️✔️: }
        // Behavior review: split_once reproduces substr-to-first-comma, including empty prefixes.
        // Complexity review: both implementations scan once to the first comma.
        mon_list
            .split_once(',')
            .map_or(mon_list, |(first, _)| first)
    }

    #[must_use]
    pub const fn kind(&self) -> EntityKind {
        self.kind
    }
    #[must_use]
    pub const fn polymer_kind(&self) -> PolymerKind {
        self.polymer_kind
    }
    #[must_use]
    pub const fn reflects_microhetero(&self) -> bool {
        self.reflects_microhetero
    }
    #[must_use]
    pub fn full_sequence(&self) -> &[String] {
        &self.full_sequence
    }
    #[must_use]
    pub fn dbrefs(&self) -> &[BioEntityDbRef] {
        &self.dbrefs
    }
    #[must_use]
    pub fn sifts_unp_accessions(&self) -> &[String] {
        &self.sifts_unp_accessions
    }
    #[must_use]
    pub fn subchains(&self) -> &[String] {
        &self.subchains
    }
    #[must_use]
    pub const fn source(&self) -> &EntitySourceIds {
        &self.source
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioCoordinateBlock {
    positions: Vec<[f64; 3]>,
}

impl BioCoordinateBlock {
    #[must_use]
    pub fn new(positions: Vec<[f64; 3]>) -> Self {
        Self { positions }
    }
    #[must_use]
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }
    #[must_use]
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    #[must_use]
    pub fn into_positions(self) -> Vec<[f64; 3]> {
        self.positions
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BioTransform {
    matrix: [[f64; 3]; 3],
    translation: [f64; 3],
}

impl Default for BioTransform {
    fn default() -> Self {
        Self::identity()
    }
}

impl BioTransform {
    #[must_use]
    pub const fn new(matrix: [[f64; 3]; 3], translation: [f64; 3]) -> Self {
        Self {
            matrix,
            translation,
        }
    }

    #[must_use]
    pub const fn identity() -> Self {
        Self::new(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [0.0, 0.0, 0.0],
        )
    }

    #[must_use]
    pub const fn matrix(&self) -> &[[f64; 3]; 3] {
        &self.matrix
    }
    #[must_use]
    pub const fn translation(&self) -> &[f64; 3] {
        &self.translation
    }

    #[must_use]
    pub fn apply(&self, point: [f64; 3]) -> [f64; 3] {
        // Gemmi✔️✔️: Vec3 apply(const Vec3& x) const { return mat.multiply(x) + vec; }
        // Behavior review: row-major matrix-vector multiplication precedes translation.
        // Complexity review: both implementations perform a fixed 3x3 multiply and add.
        let multiplied = multiply_matrix_vector(self.matrix, point);
        [
            multiplied[0] + self.translation[0],
            multiplied[1] + self.translation[1],
            multiplied[2] + self.translation[2],
        ]
    }

    #[must_use]
    pub fn combine(&self, other: &Self) -> Self {
        // Gemmi✔️✔️: Transform combine(const Transform& b) const {
        // Gemmi✔️✔️:   return {mat.multiply(b.mat), vec + mat.multiply(b.vec)};
        // Gemmi✔️✔️: }
        // Behavior review: result application is self(other(point)), in source order.
        // Complexity review: fixed 3x3 multiplication and vector arithmetic are equivalent.
        let translated = multiply_matrix_vector(self.matrix, other.translation);
        Self {
            matrix: multiply_matrices(self.matrix, other.matrix),
            translation: [
                self.translation[0] + translated[0],
                self.translation[1] + translated[1],
                self.translation[2] + translated[2],
            ],
        }
    }
}

fn multiply_matrix_vector(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ]
}

fn multiply_matrices(first: [[f64; 3]; 3], second: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            for inner in 0..3 {
                result[row][column] += first[row][inner] * second[inner][column];
            }
        }
    }
    result
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BioCrystalCell {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
}

impl Default for BioCrystalCell {
    fn default() -> Self {
        // Gemmi✔️✔️: double a = 1.0, b = 1.0, c = 1.0;
        // Gemmi✔️✔️: double alpha = 90.0, beta = 90.0, gamma = 90.0;
        // Behavior review: all six source defaults are retained exactly as f64 values.
        // Complexity review: construction is constant-time.
        Self {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioCrystalInfo {
    cell: BioCrystalCell,
    space_group_hm: Option<String>,
    z_pdb: Option<String>,
    orthogonal: BioTransform,
    fractional: BioTransform,
    volume: f64,
    reciprocal_lengths: [f64; 3],
    reciprocal_cosines: [f64; 3],
    explicit_matrices: bool,
    cs_count: i16,
    symmetry_images: Vec<BioTransform>,
}

impl BioCrystalInfo {
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        cell: BioCrystalCell,
        space_group_hm: Option<String>,
        z_pdb: Option<String>,
        orthogonal: BioTransform,
        fractional: BioTransform,
        explicit_matrices: bool,
        cs_count: i16,
        symmetry_images: Vec<BioTransform>,
    ) -> Self {
        Self {
            cell,
            space_group_hm,
            z_pdb,
            orthogonal,
            fractional,
            volume: 1.0,
            reciprocal_lengths: [1.0; 3],
            reciprocal_cosines: [0.0; 3],
            explicit_matrices,
            cs_count,
            symmetry_images,
        }
    }

    pub fn calculate_properties(&mut self) -> Result<(), BioStructureError> {
        // Gemmi✔️✔️: // ensure exact values for right angles
        // Gemmi✔️✔️: double cos_alpha = alpha == 90. ? 0. : std::cos(rad(alpha));
        // Gemmi✔️✔️: double cos_beta  = beta  == 90. ? 0. : std::cos(rad(beta));
        // Gemmi✔️✔️: double cos_gamma = gamma == 90. ? 0. : std::cos(rad(gamma));
        // Gemmi✔️✔️: double sin_alpha = alpha == 90. ? 1. : std::sin(rad(alpha));
        // Gemmi✔️✔️: double sin_beta  = beta  == 90. ? 1. : std::sin(rad(beta));
        // Gemmi✔️✔️: double sin_gamma = gamma == 90. ? 1. : std::sin(rad(gamma));
        // Gemmi✔️✔️: if (sin_alpha == 0 || sin_beta == 0 || sin_gamma == 0)
        // Gemmi✔️✔️:   fail("Impossible angle - N*180deg.");
        // Gemmi✔️✔️: volume = a * b * c * std::sqrt(1 - cos_alpha * cos_alpha
        // Gemmi✔️✔️:                                - cos_beta * cos_beta - cos_gamma * cos_gamma
        // Gemmi✔️✔️:                                + 2 * cos_alpha * cos_beta * cos_gamma);
        // Gemmi✔️✔️: constexpr double rad(double angle) { return pi() / 180.0 * angle; }
        let radians = |degrees: f64| std::f64::consts::PI / 180.0 * degrees;
        let cosine = |angle: f64| {
            if angle == 90.0 {
                0.0
            } else {
                radians(angle).cos()
            }
        };
        let sine = |angle: f64| {
            if angle == 90.0 {
                1.0
            } else {
                radians(angle).sin()
            }
        };
        let cos_alpha = cosine(self.cell.alpha);
        let cos_beta = cosine(self.cell.beta);
        let cos_gamma = cosine(self.cell.gamma);
        let sin_alpha = sine(self.cell.alpha);
        let sin_beta = sine(self.cell.beta);
        let sin_gamma = sine(self.cell.gamma);
        if sin_alpha == 0.0 || sin_beta == 0.0 || sin_gamma == 0.0 {
            return Err(BioStructureError::ImpossibleCrystalAngle);
        }
        let cos_alphar_sin_beta = (cos_beta * cos_gamma - cos_alpha) / sin_gamma;
        let volume = self.cell.a
            * self.cell.b
            * self.cell.c
            * (1.0 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma
                + 2.0 * cos_alpha * cos_beta * cos_gamma)
                .sqrt();
        let reciprocal_lengths = [
            self.cell.b * self.cell.c * sin_alpha / volume,
            self.cell.a * self.cell.c * sin_beta / volume,
            self.cell.a * self.cell.b * sin_gamma / volume,
        ];
        let cos_alphar = cos_alphar_sin_beta / sin_beta;
        let reciprocal_cosines = [
            cos_alphar,
            (cos_alpha * cos_gamma - cos_beta) / (sin_alpha * sin_gamma),
            (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta),
        ];
        self.volume = volume;
        self.reciprocal_lengths = reciprocal_lengths;
        self.reciprocal_cosines = reciprocal_cosines;

        // Gemmi✔️✔️: if (explicit_matrices)
        // Gemmi✔️✔️:   return;
        if self.explicit_matrices {
            return Ok(());
        }

        // Gemmi✔️✔️: double sin_alphar = std::sqrt(1.0 - cos_alphar * cos_alphar);
        // Gemmi✔️✔️: orth.mat = {a,  b * cos_gamma,  c * cos_beta,
        // Gemmi✔️✔️:             0., b * sin_gamma, -c * cos_alphar_sin_beta,
        // Gemmi✔️✔️:             0., 0.           ,  c * sin_beta * sin_alphar};
        // Gemmi✔️✔️: orth.vec = {0., 0., 0.};
        // Gemmi✔️✔️: double o12 = -cos_gamma / (sin_gamma * a);
        // Gemmi✔️✔️: double o13 = -(cos_gamma * cos_alphar_sin_beta + cos_beta * sin_gamma)
        // Gemmi✔️✔️:               / (sin_alphar * sin_beta * sin_gamma * a);
        // Gemmi✔️✔️: double o23 = cos_alphar / (sin_alphar * sin_gamma * b);
        // Gemmi✔️✔️: frac.mat = {1 / a,  o12,                 o13,
        // Gemmi✔️✔️:             0.,     1 / orth.mat[1][1],  o23,
        // Gemmi✔️✔️:             0.,     0.,                  1 / orth.mat[2][2]};
        // Gemmi✔️✔️: frac.vec = {0., 0., 0.};
        // Behavior review: all exact-angle, explicit-matrix, and multiplication branches match.
        // Complexity review: both implementations perform fixed-size scalar arithmetic only.
        let sin_alphar = (1.0 - cos_alphar * cos_alphar).sqrt();
        let orthogonal = BioTransform::new(
            [
                [self.cell.a, self.cell.b * cos_gamma, self.cell.c * cos_beta],
                [
                    0.0,
                    self.cell.b * sin_gamma,
                    -self.cell.c * cos_alphar_sin_beta,
                ],
                [0.0, 0.0, self.cell.c * sin_beta * sin_alphar],
            ],
            [0.0; 3],
        );
        let o12 = -cos_gamma / (sin_gamma * self.cell.a);
        let o13 = -(cos_gamma * cos_alphar_sin_beta + cos_beta * sin_gamma)
            / (sin_alphar * sin_beta * sin_gamma * self.cell.a);
        let o23 = cos_alphar / (sin_alphar * sin_gamma * self.cell.b);
        self.fractional = BioTransform::new(
            [
                [1.0 / self.cell.a, o12, o13],
                [0.0, 1.0 / orthogonal.matrix[1][1], o23],
                [0.0, 0.0, 1.0 / orthogonal.matrix[2][2]],
            ],
            [0.0; 3],
        );
        self.orthogonal = orthogonal;
        Ok(())
    }

    #[must_use]
    pub fn is_crystal(&self) -> bool {
        // Gemmi✔️✔️: bool is_crystal() const { return a != 1.0 && frac.mat[0][0] != 1.0; }
        // Behavior review: exact comparisons reproduce the source predicate without tolerance.
        // Complexity review: both implementations perform two comparisons.
        self.cell.a != 1.0 && self.fractional.matrix[0][0] != 1.0
    }

    #[must_use]
    pub const fn cell(&self) -> BioCrystalCell {
        self.cell
    }
    #[must_use]
    pub const fn orthogonal(&self) -> &BioTransform {
        &self.orthogonal
    }
    #[must_use]
    pub const fn fractional(&self) -> &BioTransform {
        &self.fractional
    }
    #[must_use]
    pub const fn volume(&self) -> f64 {
        self.volume
    }
    #[must_use]
    pub const fn reciprocal_lengths(&self) -> &[f64; 3] {
        &self.reciprocal_lengths
    }
    #[must_use]
    pub const fn reciprocal_cosines(&self) -> &[f64; 3] {
        &self.reciprocal_cosines
    }
    #[must_use]
    pub const fn explicit_matrices(&self) -> bool {
        self.explicit_matrices
    }
    #[must_use]
    pub const fn cs_count(&self) -> i16 {
        self.cs_count
    }
    #[must_use]
    pub fn symmetry_images(&self) -> &[BioTransform] {
        &self.symmetry_images
    }
    #[must_use]
    pub fn space_group_hm(&self) -> Option<&str> {
        self.space_group_hm.as_deref()
    }
    #[must_use]
    pub fn z_pdb(&self) -> Option<&str> {
        self.z_pdb.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioNcsOperator {
    pub id: String,
    pub given: bool,
    pub transform: BioTransform,
}

impl BioNcsOperator {
    #[must_use]
    pub fn new(id: String, given: bool, transform: BioTransform) -> Self {
        // Gemmi✔️✔️: struct NcsOp {
        // Gemmi✔️✔️:   std::string id;
        // Gemmi✔️✔️:   bool given;
        // Gemmi✔️✔️:   Transform tr;
        // Gemmi✔️✔️: };
        // Behavior review: all source fields are retained without transform narrowing.
        // Complexity review: the identifier is moved and fixed-size values are copied.
        Self {
            id,
            given,
            transform,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioAssemblyOperator {
    pub name: Option<String>,
    pub operator_type: Option<String>,
    pub transform: BioTransform,
}

impl BioAssemblyOperator {
    #[must_use]
    pub fn new(
        name: Option<String>,
        operator_type: Option<String>,
        transform: BioTransform,
    ) -> Self {
        // Gemmi✔️✔️: struct Operator {
        // Gemmi✔️✔️:   std::string name; // optional
        // Gemmi✔️✔️:   std::string type; // optional (from mmCIF only)
        // Gemmi✔️✔️:   Transform transform;
        // Gemmi✔️✔️: };
        // Behavior review: absence is distinct from a present empty source string.
        // Complexity review: strings are moved and the fixed-size transform is copied.
        Self {
            name,
            operator_type,
            transform,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioAssemblyGenerator {
    pub chains: Vec<String>,
    pub subchains: Vec<String>,
    pub operators: Vec<BioAssemblyOperator>,
}

impl BioAssemblyGenerator {
    #[must_use]
    pub fn new(
        chains: Vec<String>,
        subchains: Vec<String>,
        operators: Vec<BioAssemblyOperator>,
    ) -> Self {
        // Gemmi✔️✔️: struct Gen {
        // Gemmi✔️✔️:   std::vector<std::string> chains;
        // Gemmi✔️✔️:   std::vector<std::string> subchains;
        // Gemmi✔️✔️:   std::vector<Operator> operators;
        // Gemmi✔️✔️: };
        // Behavior review: source order, duplicates, and empty vectors are retained.
        // Complexity review: all vectors are moved without element cloning.
        Self {
            chains,
            subchains,
            operators,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BioAssemblySpecialKind {
    #[default]
    NotApplicable,
    CompleteIcosahedral,
    RepresentativeHelical,
    CompletePoint,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioAssembly {
    pub name: String,
    pub author_determined: bool,
    pub software_determined: bool,
    pub special_kind: BioAssemblySpecialKind,
    pub oligomeric_count: i32,
    pub oligomeric_details: String,
    pub software_name: String,
    pub buried_surface_area: f64,
    pub surface_area: f64,
    pub solvent_free_energy_change: f64,
    pub generators: Vec<BioAssemblyGenerator>,
}

impl BioAssembly {
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        name: String,
        author_determined: bool,
        software_determined: bool,
        special_kind: BioAssemblySpecialKind,
        oligomeric_count: i32,
        oligomeric_details: String,
        software_name: String,
        buried_surface_area: f64,
        surface_area: f64,
        solvent_free_energy_change: f64,
        generators: Vec<BioAssemblyGenerator>,
    ) -> Self {
        // Gemmi✔️✔️: std::string name;
        // Gemmi✔️✔️: bool author_determined = false;
        // Gemmi✔️✔️: bool software_determined = false;
        // Gemmi✔️✔️: SpecialKind special_kind = SpecialKind::NA;
        // Gemmi✔️✔️: int oligomeric_count = 0;
        // Gemmi✔️✔️: std::string oligomeric_details;
        // Gemmi✔️✔️: std::string software_name;
        // Gemmi✔️✔️: double absa = NAN;
        // Gemmi✔️✔️: double ssa = NAN;
        // Gemmi✔️✔️: double more = NAN;
        // Gemmi✔️✔️: std::vector<Gen> generators;
        // Behavior review: source ordering, duplicates, and NaN payload values remain observable.
        // Complexity review: owned strings/vectors move once and scalar state is copied.
        Self {
            name,
            author_determined,
            software_determined,
            special_kind,
            oligomeric_count,
            oligomeric_details,
            software_name,
            buried_surface_area,
            surface_area,
            solvent_free_energy_change,
            generators,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioStructureParts {
    pub input_format: BioCoordinateFormat,
    pub models: Vec<BioModelRow>,
    pub chains: Vec<BioChainRow>,
    pub residues: Vec<BioResidueRow>,
    pub atoms: Vec<BioAtomRow>,
    pub entities: Vec<BioEntityRow>,
    pub coordinates: BioCoordinateBlock,
    pub crystal: Option<BioCrystalInfo>,
    pub ncs_operators: Vec<BioNcsOperator>,
    pub assemblies: Vec<BioAssembly>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioStructure {
    input_format: BioCoordinateFormat,
    models: Vec<BioModelRow>,
    chains: Vec<BioChainRow>,
    residues: Vec<BioResidueRow>,
    atoms: Vec<BioAtomRow>,
    entities: Vec<BioEntityRow>,
    coordinates: BioCoordinateBlock,
    crystal: Option<BioCrystalInfo>,
    ncs_operators: Vec<BioNcsOperator>,
    assemblies: Vec<BioAssembly>,
}

impl BioStructure {
    pub fn from_parts(parts: BioStructureParts) -> Result<Self, BioStructureError> {
        Self::validate_parts(&parts)?;
        Ok(Self {
            input_format: parts.input_format,
            models: parts.models,
            chains: parts.chains,
            residues: parts.residues,
            atoms: parts.atoms,
            entities: parts.entities,
            coordinates: parts.coordinates,
            crystal: parts.crystal,
            ncs_operators: parts.ncs_operators,
            assemblies: parts.assemblies,
        })
    }

    pub fn validate_parts(parts: &BioStructureParts) -> Result<(), BioStructureError> {
        validate_structure(BioStructureView::from(parts))
    }

    pub fn validate(&self) -> Result<(), BioStructureError> {
        validate_structure(BioStructureView::from(self))
    }

    #[must_use]
    pub fn into_parts(self) -> BioStructureParts {
        BioStructureParts {
            input_format: self.input_format,
            models: self.models,
            chains: self.chains,
            residues: self.residues,
            atoms: self.atoms,
            entities: self.entities,
            coordinates: self.coordinates,
            crystal: self.crystal,
            ncs_operators: self.ncs_operators,
            assemblies: self.assemblies,
        }
    }

    #[must_use]
    pub const fn input_format(&self) -> BioCoordinateFormat {
        self.input_format
    }
    #[must_use]
    pub fn models(&self) -> &[BioModelRow] {
        &self.models
    }
    #[must_use]
    pub fn chains(&self) -> &[BioChainRow] {
        &self.chains
    }
    #[must_use]
    pub fn residues(&self) -> &[BioResidueRow] {
        &self.residues
    }
    #[must_use]
    pub fn atoms(&self) -> &[BioAtomRow] {
        &self.atoms
    }
    #[must_use]
    pub fn entities(&self) -> &[BioEntityRow] {
        &self.entities
    }
    #[must_use]
    pub const fn coordinates(&self) -> &BioCoordinateBlock {
        &self.coordinates
    }
    #[must_use]
    pub fn crystal(&self) -> Option<&BioCrystalInfo> {
        self.crystal.as_ref()
    }
    #[must_use]
    pub fn ncs_operators(&self) -> &[BioNcsOperator] {
        &self.ncs_operators
    }
    #[must_use]
    pub fn assemblies(&self) -> &[BioAssembly] {
        &self.assemblies
    }

    #[must_use]
    pub fn find_entity(&self, source_id: &str) -> Option<(BioEntityId, &BioEntityRow)> {
        self.entities
            .iter()
            .enumerate()
            .find_map(|(index, entity)| {
                (entity.source().source_entity_id() == source_id)
                    .then(|| (BioEntityId::new(index as u32), entity))
            })
    }

    #[must_use]
    pub fn find_entity_of_subchain(&self, subchain: &str) -> Option<(BioEntityId, &BioEntityRow)> {
        // Gemmi✔️✔️: if (subchain.empty())
        // Gemmi✔️✔️:   return nullptr;
        // Gemmi✔️✔️: for (Entity& ent : entities)
        // Gemmi✔️✔️:   if (in_vector(ent.subchains, subchain))
        // Gemmi✔️✔️:     return &ent;
        // Gemmi✔️✔️: return nullptr;
        // Behavior review: empty input is rejected and the first exact source-order match wins.
        // Complexity review: both implementations scan entities and their subchain vectors linearly.
        if subchain.is_empty() {
            return None;
        }
        self.entities
            .iter()
            .enumerate()
            .find_map(|(index, entity)| {
                entity
                    .subchains()
                    .iter()
                    .any(|candidate| candidate == subchain)
                    .then(|| (BioEntityId::new(index as u32), entity))
            })
    }

    #[must_use]
    pub fn find_atom(
        &self,
        residue_id: BioResidueId,
        name: AtomName,
        request: AltLocRequest,
        element: Option<Element>,
    ) -> Option<(BioAtomId, &BioAtomRow)> {
        let residue = self.residues.get(residue_id.index())?;
        let atoms = residue.atom_span().slice(&self.atoms).ok()?;
        atoms.iter().enumerate().find_map(|(offset, atom)| {
            (atom.name() == name
                && altloc_matches(atom.altloc(), request)
                && element.is_none_or(|expected| atom.element() == expected))
            .then(|| {
                (
                    BioAtomId::new(residue.atom_span().start() + offset as u32),
                    atom,
                )
            })
        })
    }

    pub fn atom_by_altloc(
        &self,
        residue_id: BioResidueId,
        name: AtomName,
        altloc: Option<AltLocLabel>,
    ) -> Result<(BioAtomId, &BioAtomRow), BioStructureError> {
        let residue = self.residues.get(residue_id.index()).ok_or(
            BioStructureError::RowReferenceOutOfBounds {
                table: "residues",
                index: residue_id.value(),
                table_len: self.residues.len(),
            },
        )?;
        for (offset, atom) in residue.atom_span().slice(&self.atoms)?.iter().enumerate() {
            if atom.name() == name && atom.altloc() == altloc {
                return Ok((
                    BioAtomId::new(residue.atom_span().start() + offset as u32),
                    atom,
                ));
            }
        }
        Err(BioStructureError::AtomNotFound)
    }
}

#[derive(Clone, Copy)]
struct BioStructureView<'a> {
    models: &'a [BioModelRow],
    chains: &'a [BioChainRow],
    residues: &'a [BioResidueRow],
    atoms: &'a [BioAtomRow],
    entities: &'a [BioEntityRow],
    coordinates: &'a BioCoordinateBlock,
    assemblies: &'a [BioAssembly],
}

impl<'a> From<&'a BioStructureParts> for BioStructureView<'a> {
    fn from(parts: &'a BioStructureParts) -> Self {
        Self {
            models: &parts.models,
            chains: &parts.chains,
            residues: &parts.residues,
            atoms: &parts.atoms,
            entities: &parts.entities,
            coordinates: &parts.coordinates,
            assemblies: &parts.assemblies,
        }
    }
}

impl<'a> From<&'a BioStructure> for BioStructureView<'a> {
    fn from(structure: &'a BioStructure) -> Self {
        Self {
            models: &structure.models,
            chains: &structure.chains,
            residues: &structure.residues,
            atoms: &structure.atoms,
            entities: &structure.entities,
            coordinates: &structure.coordinates,
            assemblies: &structure.assemblies,
        }
    }
}

fn validate_structure(view: BioStructureView<'_>) -> Result<(), BioStructureError> {
    validate_table_len("models", view.models.len())?;
    validate_table_len("chains", view.chains.len())?;
    validate_table_len("residues", view.residues.len())?;
    validate_table_len("atoms", view.atoms.len())?;
    validate_table_len("entities", view.entities.len())?;
    validate_table_len("assemblies", view.assemblies.len())?;
    validate_hierarchy(view)
}

fn validate_table_len(table: &'static str, len: usize) -> Result<(), BioStructureError> {
    if u32::try_from(len).is_err() {
        return Err(BioStructureError::TableTooLarge { table, len });
    }
    Ok(())
}

fn validate_hierarchy(view: BioStructureView<'_>) -> Result<(), BioStructureError> {
    let mut chain_cursor = 0_u32;
    for (model_index, model) in view.models.iter().enumerate() {
        validate_span(
            "models->chains",
            model.chain_span(),
            chain_cursor,
            view.chains.len(),
        )?;
        for chain_index in model.chain_span().start()..model.chain_span().end() {
            if view.chains[chain_index as usize].model_id() != BioModelId::new(model_index as u32) {
                return Err(BioStructureError::ParentMismatch {
                    table: "chains",
                    index: chain_index,
                });
            }
        }
        chain_cursor = model.chain_span().end();
    }
    require_full_coverage("models->chains", chain_cursor, view.chains.len())?;

    let mut residue_cursor = 0_u32;
    for (chain_index, chain) in view.chains.iter().enumerate() {
        validate_span(
            "chains->residues",
            chain.residue_span(),
            residue_cursor,
            view.residues.len(),
        )?;
        validate_entity_reference(
            view.entities,
            chain.entity_id(),
            chain.source().label_asym_id(),
        )?;
        for residue_index in chain.residue_span().start()..chain.residue_span().end() {
            if view.residues[residue_index as usize].chain_id()
                != BioChainId::new(chain_index as u32)
            {
                return Err(BioStructureError::ParentMismatch {
                    table: "residues",
                    index: residue_index,
                });
            }
        }
        residue_cursor = chain.residue_span().end();
    }
    require_full_coverage("chains->residues", residue_cursor, view.residues.len())?;

    let mut atom_cursor = 0_u32;
    for (residue_index, residue) in view.residues.iter().enumerate() {
        validate_span(
            "residues->atoms",
            residue.atom_span(),
            atom_cursor,
            view.atoms.len(),
        )?;
        validate_entity_reference(
            view.entities,
            residue.entity_id(),
            residue.source().subchain_id(),
        )?;
        for atom_index in residue.atom_span().start()..residue.atom_span().end() {
            if view.atoms[atom_index as usize].residue_id()
                != BioResidueId::new(residue_index as u32)
            {
                return Err(BioStructureError::ParentMismatch {
                    table: "atoms",
                    index: atom_index,
                });
            }
        }
        atom_cursor = residue.atom_span().end();
    }
    require_full_coverage("residues->atoms", atom_cursor, view.atoms.len())?;
    if view.coordinates.len() != view.atoms.len() {
        return Err(BioStructureError::CoordinateCountMismatch {
            atom_count: view.atoms.len(),
            coordinate_count: view.coordinates.len(),
        });
    }
    validate_assembly_references(view)?;
    Ok(())
}

fn validate_assembly_references(view: BioStructureView<'_>) -> Result<(), BioStructureError> {
    for (assembly_index, assembly) in view.assemblies.iter().enumerate() {
        for generator in &assembly.generators {
            for requested in &generator.chains {
                let present = view.chains.iter().any(|chain| {
                    chain
                        .source()
                        .auth_chain_id()
                        .is_some_and(|id| id.as_str() == requested)
                });
                if !present {
                    return Err(BioStructureError::AssemblyReferenceMissing {
                        assembly: BioAssemblyId::new(assembly_index as u32),
                        kind: "chain",
                        value: requested.clone(),
                    });
                }
            }
            for requested in &generator.subchains {
                let chain_present = view
                    .chains
                    .iter()
                    .any(|chain| chain.source().label_asym_id() == Some(requested.as_str()));
                let residue_present = view
                    .residues
                    .iter()
                    .any(|residue| residue.source().subchain_id() == Some(requested.as_str()));
                if !chain_present && !residue_present {
                    return Err(BioStructureError::AssemblyReferenceMissing {
                        assembly: BioAssemblyId::new(assembly_index as u32),
                        kind: "subchain",
                        value: requested.clone(),
                    });
                }
            }
        }
    }
    Ok(())
}

fn validate_span<I>(
    table: &'static str,
    span: BioRowSpan<I>,
    expected_start: u32,
    child_len: usize,
) -> Result<(), BioStructureError> {
    if span.start() != expected_start {
        return Err(BioStructureError::NonContiguousSpan {
            table,
            expected_start,
            actual_start: span.start(),
        });
    }
    if span.end() as usize > child_len {
        return Err(BioStructureError::RowSpanOutOfBounds {
            start: span.start(),
            len: span.len(),
            table_len: child_len,
        });
    }
    Ok(())
}

fn require_full_coverage(
    table: &'static str,
    covered: u32,
    child_len: usize,
) -> Result<(), BioStructureError> {
    if covered as usize != child_len {
        return Err(BioStructureError::IncompleteCoverage {
            table,
            covered,
            table_len: child_len,
        });
    }
    Ok(())
}

fn validate_entity_reference(
    entities: &[BioEntityRow],
    entity_id: Option<BioEntityId>,
    subchain: Option<&str>,
) -> Result<(), BioStructureError> {
    let Some(entity_id) = entity_id else {
        return Ok(());
    };
    let entity =
        entities
            .get(entity_id.index())
            .ok_or(BioStructureError::RowReferenceOutOfBounds {
                table: "entities",
                index: entity_id.value(),
                table_len: entities.len(),
            })?;
    if let Some(subchain) = subchain
        && !subchain.is_empty()
        && !entity
            .subchains()
            .iter()
            .any(|candidate| candidate == subchain)
    {
        return Err(BioStructureError::EntitySubchainMismatch {
            entity_id,
            subchain: subchain.to_owned(),
        });
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BioStructureError {
    RowIndexTooLarge {
        value: usize,
    },
    RowSpanOverflow {
        start: u32,
        len: u32,
    },
    RowSpanOutOfBounds {
        start: u32,
        len: u32,
        table_len: usize,
    },
    TableTooLarge {
        table: &'static str,
        len: usize,
    },
    NonContiguousSpan {
        table: &'static str,
        expected_start: u32,
        actual_start: u32,
    },
    IncompleteCoverage {
        table: &'static str,
        covered: u32,
        table_len: usize,
    },
    ParentMismatch {
        table: &'static str,
        index: u32,
    },
    RowReferenceOutOfBounds {
        table: &'static str,
        index: u32,
        table_len: usize,
    },
    CoordinateCountMismatch {
        atom_count: usize,
        coordinate_count: usize,
    },
    EntitySubchainMismatch {
        entity_id: BioEntityId,
        subchain: String,
    },
    AssemblyReferenceMissing {
        assembly: BioAssemblyId,
        kind: &'static str,
        value: String,
    },
    ImpossibleCrystalAngle,
    AtomNotFound,
}

impl fmt::Display for BioStructureError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for BioStructureError {}
