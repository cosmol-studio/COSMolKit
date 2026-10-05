//! Validated detached macromolecular hierarchy values.

use std::{fmt, marker::PhantomData};

use cosmolkit_types::Element;

use crate::{
    AltLocLabel, AtomName, AtomSourceIds, BioCisPep, BioConnection, BioHelix, BioMetadata,
    BioModRes, BioSheet, BioStructureSourceState, ChainSourceIds, EntitySourceIds, PdbSeqId,
    ResidueInfoKind, ResidueName, ResidueSourceIds,
};

mod lattice;
mod spacegroup;
pub use lattice::{BioNearestImage, find_nearest_image};
pub use spacegroup::setup_cell_images;

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
    pdb_coordinate_text: Option<[u8; 24]>,
    residue_id: BioResidueId,
    name: AtomName,
    element: Element,
    isotope_mass_number: Option<u16>,
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
    /// Narrow IO provenance: exact PDB coordinate columns, independent of
    /// parsed coordinates and retained when a detached row is copied.
    #[doc(hidden)]
    pub fn with_pdb_coordinate_text(mut self, text: [u8; 24]) -> Self {
        self.pdb_coordinate_text = Some(text);
        self
    }
    #[doc(hidden)]
    pub const fn pdb_coordinate_text(&self) -> Option<&[u8; 24]> {
        self.pdb_coordinate_text.as_ref()
    }
    /// `isotope_mass_number` is independent of `element` and `name`. `None`
    /// means that no isotope mass is specified; it does not mean explicit
    /// protium. This value constructor does not impose physical isotope rules.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        residue_id: BioResidueId,
        name: AtomName,
        element: Element,
        isotope_mass_number: Option<u16>,
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
            pdb_coordinate_text: None,
            residue_id,
            name,
            element,
            isotope_mass_number,
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
    /// Returns the represented isotope mass number, if one was specified.
    /// `None` means unspecified, not explicit mass number 1.
    #[must_use]
    pub const fn isotope_mass_number(&self) -> Option<u16> {
        self.isotope_mass_number
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
    /// Mutable detached coordinate rows; live public objects do not expose this.
    pub fn positions_mut(&mut self) -> &mut [[f64; 3]] {
        &mut self.positions
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
    /// Source affine-transform approximate equality, including its asymmetric NaN treatment.
    #[must_use]
    pub fn approx(&self, other: &Self, epsilon: f64) -> bool {
        // Gemmi❗✔️:   bool approx(const Transform& o, double epsilon) const {
        // Gemmi❗✔️:     return mat.approx(o.mat, epsilon) && vec.approx(o.vec, epsilon);
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   double trace() const { return a[0][0] + a[1][1] + a[2][2]; }
        // Gemmi❗✔️:
        // Gemmi❗✔️:   bool approx(const Mat33& other, double epsilon) const {
        // Gemmi❗✔️:     for (int i = 0; i < 3; ++i)
        // Gemmi❗✔️:       for (int j = 0; j < 3; ++j)
        // Gemmi❗✔️:         if (std::fabs(a[i][j] - other.a[i][j]) > epsilon)
        // Gemmi❗✔️:           return false;
        // Gemmi❗✔️:     return true;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   bool approx(const Vec3_& o, Real epsilon) const {
        // Gemmi❗✔️:     return std::fabs(x - o.x) <= epsilon &&
        // Gemmi❗✔️:            std::fabs(y - o.y) <= epsilon &&
        // Gemmi❗✔️:            std::fabs(z - o.z) <= epsilon;
        // Gemmi❗✔️:   }
        // Cost: fixed 3x3 and three-vector comparisons, without allocation.
        for i in 0..3 {
            for j in 0..3 {
                if (self.matrix[i][j] - other.matrix[i][j]).abs() > epsilon {
                    return false;
                }
            }
        }
        (0..3).all(|i| (self.translation[i] - other.translation[i]).abs() <= epsilon)
    }

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
    // Gemmi❗✔️: Vec3 multiply(const Vec3& p) const {
    // Gemmi❗✔️:   return {a[0][0] * p.x + a[0][1] * p.y + a[0][2] * p.z,
    // Gemmi❗✔️:           a[1][0] * p.x + a[1][1] * p.y + a[1][2] * p.z,
    // Gemmi❗✔️:           a[2][0] * p.x + a[2][1] * p.y + a[2][2] * p.z};
    // Gemmi❗✔️: }
    // Behavior review: row-major dot products retain Gemmi's per-row multiply/add order.
    // Complexity review: both implementations perform three fixed-length dot products.
    [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ]
}

fn matrix_determinant(matrix: [[f64; 3]; 3]) -> f64 {
    // Gemmi❗✔️: double determinant() const {
    // Gemmi❗✔️:   return a[0][0] * (a[1][1]*a[2][2] - a[2][1]*a[1][2]) +
    // Gemmi❗✔️:          a[0][1] * (a[1][2]*a[2][0] - a[2][2]*a[1][0]) +
    // Gemmi❗✔️:          a[0][2] * (a[1][0]*a[2][1] - a[2][0]*a[1][1]);
    // Gemmi❗✔️: }
    // Behavior review: preserve the source cofactor and addition order; parity awaits the
    // dedicated nonsymmetric and singular matrix regressions.
    // Complexity review: both implementations use a fixed number of scalar operations.
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2])
        + matrix[0][1] * (matrix[1][2] * matrix[2][0] - matrix[2][2] * matrix[1][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1])
}

fn inverse_bio_matrix(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    // Gemmi❗✔️: Mat33 inverse() const {
    // Gemmi❗✔️:   Mat33 inv;
    // Gemmi❗✔️:   double inv_det = 1.0 / determinant();
    // Gemmi❗✔️:   inv[0][0] = inv_det * (a[1][1] * a[2][2] - a[2][1] * a[1][2]);
    // Gemmi❗✔️:   inv[0][1] = inv_det * (a[0][2] * a[2][1] - a[0][1] * a[2][2]);
    // Gemmi❗✔️:   inv[0][2] = inv_det * (a[0][1] * a[1][2] - a[0][2] * a[1][1]);
    // Gemmi❗✔️:   inv[1][0] = inv_det * (a[1][2] * a[2][0] - a[1][0] * a[2][2]);
    // Gemmi❗✔️:   inv[1][1] = inv_det * (a[0][0] * a[2][2] - a[0][2] * a[2][0]);
    // Gemmi❗✔️:   inv[1][2] = inv_det * (a[1][0] * a[0][2] - a[0][0] * a[1][2]);
    // Gemmi❗✔️:   inv[2][0] = inv_det * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
    // Gemmi❗✔️:   inv[2][1] = inv_det * (a[2][0] * a[0][1] - a[0][0] * a[2][1]);
    // Gemmi❗✔️:   inv[2][2] = inv_det * (a[0][0] * a[1][1] - a[1][0] * a[0][1]);
    // Gemmi❗✔️:   return inv;
    // Gemmi❗✔️: }
    // Behavior review: all nine outputs follow Gemmi's cofactor expressions and the
    // singular case intentionally uses IEEE division without a Rust-side guard.
    // Complexity review: both implementations are fixed-size, allocation-free 3x3 work.
    let inverse_determinant = 1.0 / matrix_determinant(matrix);
    [
        [
            inverse_determinant * (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2]),
            inverse_determinant * (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]),
            inverse_determinant * (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]),
        ],
        [
            inverse_determinant * (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]),
            inverse_determinant * (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]),
            inverse_determinant * (matrix[1][0] * matrix[0][2] - matrix[0][0] * matrix[1][2]),
        ],
        [
            inverse_determinant * (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]),
            inverse_determinant * (matrix[2][0] * matrix[0][1] - matrix[0][0] * matrix[2][1]),
            inverse_determinant * (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]),
        ],
    ]
}

fn inverse_bio_transform(transform: BioTransform) -> BioTransform {
    // Gemmi❗✔️: Transform inverse() const {
    // Gemmi❗✔️:   Mat33 minv = mat.inverse();
    // Gemmi❗✔️:   return {minv, minv.multiply(vec).negated()};
    // Gemmi❗✔️: }
    // Gemmi❗✔️: Vec3_ negated() const { return {-x, -y, -z}; }
    // Behavior review: the inverse translation is the negated matrix-vector product,
    // in the same operation order; regression coverage is scheduled separately.
    // Complexity review: matrix inversion and translation each have fixed 3D cost.
    let matrix = inverse_bio_matrix(transform.matrix);
    let translated = multiply_matrix_vector(matrix, transform.translation);
    BioTransform::new(matrix, [-translated[0], -translated[1], -translated[2]])
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
    /// International Tables number selected by Gemmi's structure space-group lookup.
    #[must_use]
    pub fn space_group_number(&self) -> Option<i32> {
        // Gemmi❗✔️: if (const SpaceGroup* sg = st.find_spacegroup())
        // Gemmi❗✔️:   span.set_pair("_symmetry.Int_Tables_number", std::to_string(sg->number));
        // Behavior: expose only the selected scalar; lookup stays in the BIO table owner.
        // Cost: no additional table scan or allocation after lookup.
        spacegroup::structure_space_group_number(self)
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

pub fn set_crystal_cell(
    crystal: &mut BioCrystalInfo,
    cell: BioCrystalCell,
) -> Result<(), BioStructureError> {
    // Gemmi❗✔️: void set(double a_, double b_, double c_,
    // Gemmi❗✔️:          double alpha_, double beta_, double gamma_) {
    // Gemmi❗✔️:   if (gamma_ == 0.0)  // ignore empty/partial CRYST1 (example: 3iyp)
    // Gemmi❗✔️:     return;
    // Gemmi❗✔️:   a = a_;
    // Gemmi❗✔️:   b = b_;
    // Gemmi❗✔️:   c = c_;
    // Gemmi❗✔️:   alpha = alpha_;
    // Gemmi❗✔️:   beta = beta_;
    // Gemmi❗✔️:   gamma = gamma_;
    // Gemmi❗✔️:   calculate_properties();
    // Gemmi❗✔️: }
    // Behavior review: the gamma-zero sentinel is a no-op; otherwise the complete
    // six-value cell is installed before calculation, so source-ordered failures
    // retain the new cell and prior derived fields. Existing calculation preserves
    // explicit matrices while refreshing derived values.
    // Complexity review: both paths perform constant-size field and scalar work.
    if cell.gamma == 0.0 {
        return Ok(());
    }
    crystal.cell = cell;
    crystal.calculate_properties()
}

pub fn set_crystal_fractional_transform(crystal: &mut BioCrystalInfo, transform: BioTransform) {
    // Gemmi❗✔️: void set_matrices_from_fract(const Transform& f) {
    // Gemmi❗✔️:   // mmCIF _atom_sites.fract_transf_* and PDB SCALEn records usually contain
    // Gemmi❗✔️:   // fewer significant digits than the unit cell parameters, and sometimes are
    // Gemmi❗✔️:   // just wrong. Use them only if we seem to have non-standard crystal frame.
    // Gemmi❗✔️:   if (f.mat.approx(frac.mat, 1e-4) && f.vec.approx(frac.vec, 1e-6))
    // Gemmi❗✔️:     return;
    // Gemmi❗✔️:   // The SCALE record is sometimes incorrect. Here we only catch cases
    // Gemmi❗✔️:   // when CRYST1 is set as for non-crystal and SCALE is very suspicious.
    // Gemmi❗✔️:   if (frac.mat[0][0] == 1.0 && (f.mat[0][0] == 0.0 || f.mat[0][0] > 1.0))
    // Gemmi❗✔️:     return;
    // Gemmi❗✔️:   frac = f;
    // Gemmi❗✔️:   orth = f.inverse();
    // Gemmi❗✔️:   explicit_matrices = true;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool approx(const Mat33& other, double epsilon) const {
    // Gemmi❗✔️:   for (int i = 0; i < 3; ++i)
    // Gemmi❗✔️:     for (int j = 0; j < 3; ++j)
    // Gemmi❗✔️:       if (std::fabs(a[i][j] - other.a[i][j]) > epsilon)
    // Gemmi❗✔️:         return false;
    // Gemmi❗✔️:   return true;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool approx(const Vec3_& o, Real epsilon) const {
    // Gemmi❗✔️:   return std::fabs(x - o.x) <= epsilon &&
    // Gemmi❗✔️:          std::fabs(y - o.y) <= epsilon &&
    // Gemmi❗✔️:          std::fabs(z - o.z) <= epsilon;
    // Gemmi❗✔️: }
    // Behavior review: matrix comparisons reject only differences strictly greater
    // than 1e-4 (so NaN differences pass that predicate), while vector comparisons
    // require each difference <= 1e-6 (so NaN fails). The source suspicious-matrix
    // guard and assignment -> inverse -> explicit flag order are retained.
    // Complexity review: both implementations use nine matrix comparisons, at most
    // three vector comparisons, and fixed-size inversion; no allocation or rescans.
    let mut matrices_approximate = true;
    for row in 0..3 {
        for column in 0..3 {
            if (transform.matrix[row][column] - crystal.fractional.matrix[row][column]).abs()
                > 1.0e-4
            {
                matrices_approximate = false;
                break;
            }
        }
        if !matrices_approximate {
            break;
        }
    }
    if matrices_approximate
        && (transform.translation[0] - crystal.fractional.translation[0]).abs() <= 1.0e-6
        && (transform.translation[1] - crystal.fractional.translation[1]).abs() <= 1.0e-6
        && (transform.translation[2] - crystal.fractional.translation[2]).abs() <= 1.0e-6
    {
        return;
    }
    if crystal.fractional.matrix[0][0] == 1.0
        && (transform.matrix[0][0] == 0.0 || transform.matrix[0][0] > 1.0)
    {
        return;
    }
    crystal.fractional = transform;
    crystal.orthogonal = inverse_bio_transform(transform);
    crystal.explicit_matrices = true;
}

/// Replace the stored PDB Hermann–Mauguin space-group text.
///
/// The PDB reader applies Gemmi's `len > 56` field-presence condition and
/// `read_string` extraction before calling this BIO-owned state transition.
pub fn set_crystal_space_group_hm(crystal: &mut BioCrystalInfo, value: String) {
    // Gemmi❗✔️: st.spacegroup_hm = read_string(line+55, 11);
    // Behavior review: for a source-present field, replace the prior value even
    // when the already-trimmed source text is empty; Some("") retains that
    // assignment distinctly from a field that was not present.
    // Complexity review: this is one owned String replacement with no scan.
    crystal.space_group_hm = Some(value);
}

/// Replace PDB Z metadata only when the decoded source field is nonempty.
pub fn set_crystal_z_pdb_if_nonempty(crystal: &mut BioCrystalInfo, value: String) {
    // Gemmi❗✔️:         if (!z.empty())
    // Gemmi❗✔️:           st.info["_cell.Z_PDB"] = z;
    // Behavior review: an empty decoded field leaves the previous metadata
    // untouched; a nonempty field replaces it. The PDB caller owns the source
    // length guard and fixed-width read_string operation.
    // Complexity review: one emptiness check and, only when nonempty, one owned
    // String replacement; no collection traversal or allocation on empty input.
    if !value.is_empty() {
        crystal.z_pdb = Some(value);
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
    pub connections: Vec<BioConnection>,
    pub cispeps: Vec<BioCisPep>,
    pub mod_residues: Vec<BioModRes>,
    pub helices: Vec<BioHelix>,
    pub sheets: Vec<BioSheet>,
    pub metadata: BioMetadata,
    pub source_state: BioStructureSourceState,
    pub coordinates: BioCoordinateBlock,
    pub crystal: Option<BioCrystalInfo>,
    pub ncs_operators: Vec<BioNcsOperator>,
    pub assemblies: Vec<BioAssembly>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BioStructureData {
    pub input_format: BioCoordinateFormat,
    pub models: std::sync::Arc<Vec<BioModelRow>>,
    pub chains: std::sync::Arc<Vec<BioChainRow>>,
    pub residues: std::sync::Arc<Vec<BioResidueRow>>,
    pub atoms: std::sync::Arc<Vec<BioAtomRow>>,
    pub entities: std::sync::Arc<Vec<BioEntityRow>>,
    pub connections: std::sync::Arc<Vec<BioConnection>>,
    pub cispeps: std::sync::Arc<Vec<BioCisPep>>,
    pub mod_residues: std::sync::Arc<Vec<BioModRes>>,
    pub helices: std::sync::Arc<Vec<BioHelix>>,
    pub sheets: std::sync::Arc<Vec<BioSheet>>,
    pub metadata: std::sync::Arc<BioMetadata>,
    pub source_state: std::sync::Arc<BioStructureSourceState>,
    pub coordinates: std::sync::Arc<BioCoordinateBlock>,
    pub crystal: std::sync::Arc<Option<BioCrystalInfo>>,
    pub ncs_operators: std::sync::Arc<Vec<BioNcsOperator>>,
    pub assemblies: std::sync::Arc<Vec<BioAssembly>>,
}

/// Borrowed detached input for the selection copier; no live object or mutation authority.
///
/// Arc references preserve unchanged metadata without cloning the source tables.
#[derive(Clone, Copy)]
pub struct BioStructureCopySource<'a> {
    pub input_format: BioCoordinateFormat,
    pub models: &'a std::sync::Arc<Vec<BioModelRow>>,
    pub chains: &'a std::sync::Arc<Vec<BioChainRow>>,
    pub residues: &'a std::sync::Arc<Vec<BioResidueRow>>,
    pub atoms: &'a std::sync::Arc<Vec<BioAtomRow>>,
    pub entities: &'a std::sync::Arc<Vec<BioEntityRow>>,
    pub connections: &'a std::sync::Arc<Vec<BioConnection>>,
    pub cispeps: &'a std::sync::Arc<Vec<BioCisPep>>,
    pub mod_residues: &'a std::sync::Arc<Vec<BioModRes>>,
    pub helices: &'a std::sync::Arc<Vec<BioHelix>>,
    pub sheets: &'a std::sync::Arc<Vec<BioSheet>>,
    pub metadata: &'a std::sync::Arc<BioMetadata>,
    pub source_state: &'a std::sync::Arc<BioStructureSourceState>,
    pub coordinates: &'a std::sync::Arc<BioCoordinateBlock>,
    pub crystal: &'a std::sync::Arc<Option<BioCrystalInfo>>,
    pub ncs_operators: &'a std::sync::Arc<Vec<BioNcsOperator>>,
    pub assemblies: &'a std::sync::Arc<Vec<BioAssembly>>,
}

impl<'a> From<&'a BioStructureData> for BioStructureCopySource<'a> {
    fn from(data: &'a BioStructureData) -> Self {
        Self {
            input_format: data.input_format,
            models: &data.models,
            chains: &data.chains,
            residues: &data.residues,
            atoms: &data.atoms,
            entities: &data.entities,
            connections: &data.connections,
            cispeps: &data.cispeps,
            mod_residues: &data.mod_residues,
            helices: &data.helices,
            sheets: &data.sheets,
            metadata: &data.metadata,
            source_state: &data.source_state,
            coordinates: &data.coordinates,
            crystal: &data.crystal,
            ncs_operators: &data.ncs_operators,
            assemblies: &data.assemblies,
        }
    }
}

impl BioStructureCopySource<'_> {
    pub(crate) fn validate(&self) -> Result<(), BioStructureError> {
        validate_structure(BioStructureView {
            models: self.models,
            chains: self.chains,
            residues: self.residues,
            atoms: self.atoms,
            entities: self.entities,
            coordinates: self.coordinates,
            assemblies: self.assemblies,
        })
    }
}

impl BioStructureData {
    pub fn from_parts(parts: BioStructureParts) -> Result<Self, BioStructureError> {
        Self::validate_parts(&parts)?;
        // Gemmi✔️🔝:     st.connections = connections;
        // Gemmi✔️🔝:     st.cispeps = cispeps;
        // Gemmi✔️🔝:     st.mod_residues = mod_residues;
        // Gemmi✔️🔝:     st.helices = helices;
        // Gemmi✔️🔝:     st.sheets = sheets;
        // Gemmi✔️🔝:     st.meta = meta;
        // Gemmi✔️🔝:     st.input_format = input_format;
        // Gemmi✔️🔝:     st.has_origx = has_origx;
        // Gemmi✔️🔝:     st.origx = origx;
        // Gemmi✔️🔝:     st.info = info;
        // Gemmi✔️🔝:     st.raw_remarks = raw_remarks;
        // Gemmi✔️🔝:     st.resolution = resolution;
        // Behavior review: after full existing structure validation, each
        // owned source relationship/metadata value moves intact into the
        // validated BioStructureData. Its source-address references are not
        // reinterpreted as BIO row ids.
        // Complexity review: moving the vectors and aggregate values is O(1)
        // per field and avoids deep element copies performed by Gemmi's
        // `empty_copy`; validation retains its existing independent cost.
        Ok(Self {
            input_format: parts.input_format,
            models: std::sync::Arc::new(parts.models),
            chains: std::sync::Arc::new(parts.chains),
            residues: std::sync::Arc::new(parts.residues),
            atoms: std::sync::Arc::new(parts.atoms),
            entities: std::sync::Arc::new(parts.entities),
            connections: std::sync::Arc::new(parts.connections),
            cispeps: std::sync::Arc::new(parts.cispeps),
            mod_residues: std::sync::Arc::new(parts.mod_residues),
            helices: std::sync::Arc::new(parts.helices),
            sheets: std::sync::Arc::new(parts.sheets),
            metadata: std::sync::Arc::new(parts.metadata),
            source_state: std::sync::Arc::new(parts.source_state),
            coordinates: std::sync::Arc::new(parts.coordinates),
            crystal: std::sync::Arc::new(parts.crystal),
            ncs_operators: std::sync::Arc::new(parts.ncs_operators),
            assemblies: std::sync::Arc::new(parts.assemblies),
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
            models: std::sync::Arc::unwrap_or_clone(self.models),
            chains: std::sync::Arc::unwrap_or_clone(self.chains),
            residues: std::sync::Arc::unwrap_or_clone(self.residues),
            atoms: std::sync::Arc::unwrap_or_clone(self.atoms),
            entities: std::sync::Arc::unwrap_or_clone(self.entities),
            connections: std::sync::Arc::unwrap_or_clone(self.connections),
            cispeps: std::sync::Arc::unwrap_or_clone(self.cispeps),
            mod_residues: std::sync::Arc::unwrap_or_clone(self.mod_residues),
            helices: std::sync::Arc::unwrap_or_clone(self.helices),
            sheets: std::sync::Arc::unwrap_or_clone(self.sheets),
            metadata: std::sync::Arc::unwrap_or_clone(self.metadata),
            source_state: std::sync::Arc::unwrap_or_clone(self.source_state),
            coordinates: std::sync::Arc::unwrap_or_clone(self.coordinates),
            crystal: std::sync::Arc::unwrap_or_clone(self.crystal),
            ncs_operators: std::sync::Arc::unwrap_or_clone(self.ncs_operators),
            assemblies: std::sync::Arc::unwrap_or_clone(self.assemblies),
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
    pub fn atom_position(&self, atom: BioAtomId) -> Option<[f64; 3]> {
        self.coordinates.positions().get(atom.index()).copied()
    }
    #[must_use]
    pub fn residue_atoms(&self, residue: BioResidueId) -> Option<&[BioAtomRow]> {
        let row = self.residues.get(residue.index())?;
        let span = row.atom_span();
        self.atoms.get(span.start() as usize..span.end() as usize)
    }
    #[must_use]
    pub fn entities(&self) -> &[BioEntityRow] {
        &self.entities
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
    pub fn helices(&self) -> &[BioHelix] {
        &self.helices
    }
    #[must_use]
    pub fn sheets(&self) -> &[BioSheet] {
        &self.sheets
    }
    #[must_use]
    pub fn metadata(&self) -> &BioMetadata {
        &self.metadata
    }
    #[must_use]
    pub fn source_state(&self) -> &BioStructureSourceState {
        &self.source_state
    }
    /// Borrow the last identity NCS operation ID retained by the mmCIF reader.
    #[must_use]
    pub fn ncs_oper_identity_id(&self) -> Option<&str> {
        // Gemmi✔️✔️:             st.info["_struct_ncs_oper.id"] = op.str(12);
        // Behavior: identity IDs reside in the existing source map; the
        // reader overwrites this one key for each later identity row.
        // Complexity: one ordered-map lookup, no allocation or operator scan.
        self.source_state
            .info
            .get("_struct_ncs_oper.id")
            .map(String::as_str)
    }
    #[must_use]
    pub fn coordinates(&self) -> &BioCoordinateBlock {
        &self.coordinates
    }
    #[must_use]
    pub fn crystal(&self) -> Option<&BioCrystalInfo> {
        self.crystal.as_ref().as_ref()
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

    /// Resolve a source address in one model, preserving the first residue match.
    pub fn find_cra(
        &self,
        model_id: BioModelId,
        address: &crate::AtomAddress,
        ignore_segment: bool,
    ) -> Result<Option<(BioChainId, BioResidueId, Option<BioAtomId>)>, BioStructureError> {
        // Gemmi❗✔️:   CRA find_cra(const AtomAddress& address, bool ignore_segment=false) {
        // Gemmi❗✔️:     for (Chain& chain : chains)
        // Gemmi❗✔️:       if (chain.name == address.chain_name) {
        // Gemmi❗✔️:         for (Residue& res : chain.residues)
        // Gemmi❗✔️:           if (address.res_id.matches_noseg(res) &&
        // Gemmi❗✔️:               (ignore_segment || address.res_id.segment == res.segment)) {
        // Gemmi❗✔️:             Atom *at = nullptr;
        // Gemmi❗✔️:             if (!address.atom_name.empty())
        // Gemmi❗✔️:               at = res.find_atom(address.atom_name, address.altloc);
        // Gemmi❗✔️:             return {&chain, &res, at};
        // Gemmi❗✔️:           }
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     return {nullptr, nullptr, nullptr};
        // Gemmi❗✔️:   }
        // Behavior: delegates residue matching and altloc matching to the canonical BIO owners.
        // Cost: one ordered chain/residue/atom traversal, no cloned rows or temporary atom names.
        let model = self.models.get(model_id.index()).ok_or(
            BioStructureError::RowReferenceOutOfBounds {
                table: "models",
                index: model_id.value(),
                table_len: self.models.len(),
            },
        )?;
        for (chain_offset, chain) in model.chain_span().slice(&self.chains)?.iter().enumerate() {
            if !chain.source().auth_chain_id().map_or_else(
                || address.chain_name().as_str().is_empty(),
                |id| id == address.chain_name(),
            ) {
                continue;
            }
            for (residue_offset, residue) in chain
                .residue_span()
                .slice(&self.residues)?
                .iter()
                .enumerate()
            {
                let seq = residue.source().seq_id();
                let segment = residue.source().segment_id().map_or(&[][..], |bytes| {
                    let end = bytes
                        .iter()
                        .rposition(|b| *b != 0 && *b != b' ')
                        .map_or(0, |i| i + 1);
                    &bytes[..end]
                });
                let candidate = crate::ResidueAddress::new(
                    seq.map(|s| s.seq_num()),
                    seq.and_then(|s| s.ins_code()),
                    segment,
                    residue.name(),
                )
                .expect("validated BIO source identifiers are bounded ASCII");
                if !address.residue().matches_without_segment(&candidate)
                    || (!ignore_segment && address.residue().segment() != candidate.segment())
                {
                    continue;
                }
                let chain_id = BioChainId::new(model.chain_span().start() + chain_offset as u32);
                let residue_id =
                    BioResidueId::new(chain.residue_span().start() + residue_offset as u32);
                let atom_id = if address.logical_atom_name().is_empty() {
                    None
                } else {
                    residue
                        .atom_span()
                        .slice(&self.atoms)?
                        .iter()
                        .enumerate()
                        .find_map(|(offset, atom)| {
                            let requested = if address.altloc() == b'*' {
                                AltLocRequest::Any
                            } else {
                                AltLocRequest::Exact(
                                    (address.altloc() != 0)
                                        .then(|| crate::AltLocLabel::new(address.altloc())),
                                )
                            };
                            (atom_name_logical_view(&atom.name(), self.input_format)
                                == address.logical_atom_name()
                                && altloc_matches(atom.altloc(), requested))
                            .then(|| BioAtomId::new(residue.atom_span().start() + offset as u32))
                        })
                };
                return Ok(Some((chain_id, residue_id, atom_id)));
            }
        }
        Ok(None)
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

impl<'a> From<&'a BioStructureData> for BioStructureView<'a> {
    fn from(structure: &'a BioStructureData) -> Self {
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
    // Assembly chain/subchain names are SOURCE METADATA, not live CK row
    // references (approved contract, bio_architecture.md "Selection copying
    // follows pinned Gemmi"): construction must not reject, prune or
    // rewrite them merely because rows are absent — selection copying can
    // legitimately drop every named row while retaining the assigned
    // assembly. Local row IDs, parents/spans, entity-row references and
    // coordinate alignment remain validated above; a consumer needing
    // actual assembly targets must resolve the names itself and handle
    // absence explicitly (no such consumer exists yet).
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
    EmptyResidueSpan {
        operation: &'static str,
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

/// Shared private `read_string` logical view over raw stored ASCII bytes
/// (ROW-NAME Step 2): reproduces the pinned PDB reader helper in exact
/// source order — is_space left trim (bytes 9-13 and 32), LF/CR/NUL
/// termination of the remaining field, then is_space right trim. Borrowed,
/// allocation-free.
pub(crate) fn pdb_read_string_view(stored: &str) -> &str {
    // Gemmi✔️✔️: inline bool is_space(char c) {
    // Gemmi✔️✔️:   static const std::uint8_t table[256] = { // 1 for 9-13 and 32
    // Gemmi✔️✔️:     0,0,0,0,0,0,0,0, 0,1,1,1,1,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi✔️✔️:     1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, ... remaining rows all zero ... };
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: std::string read_string(const char* p, int field_length) {
    // Gemmi✔️✔️:   // left trim
    // Gemmi✔️✔️:   while (field_length != 0 && is_space(*p)) {
    // Gemmi✔️✔️:     ++p;
    // Gemmi✔️✔️:     --field_length;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   // EOL/EOF ends the string
    // Gemmi✔️✔️:   for (int i = 0; i < field_length; ++i)
    // Gemmi✔️✔️:     if (p[i] == '\n' || p[i] == '\r' || p[i] == '\0') {
    // Gemmi✔️✔️:       field_length = i;
    // Gemmi✔️✔️:       break;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   // right trim
    // Gemmi✔️✔️:   while (field_length != 0 && is_space(p[field_length-1]))
    // Gemmi✔️✔️:     --field_length;
    // Gemmi✔️✔️:   return std::string(p, field_length);
    // Gemmi✔️✔️: }
    // Behavior: exact source order. The left trim consumes bytes 9-13 and
    // 32; the termination scan runs over the REMAINING field (p advanced
    // by the left trim) and cuts at the first NUL/LF/CR, discarding
    // everything after it; the right trim then removes trailing is_space
    // bytes of the cut field. Interior TAB/VT/FF (9/11/12) are is_space
    // but NOT terminators, so an interior TAB is preserved while an
    // interior LF/CR or NUL truncates. The stored name is
    // constructor-guaranteed ASCII (<= 4 bytes), so the byte view is valid
    // UTF-8; the debug_assert documents that invariant rather than
    // substituting for it.
    // Complexity: at most two linear passes over <= 4 bytes, borrowed
    // output, no allocation.
    fn is_space(b: u8) -> bool {
        matches!(b, 9..=13 | 32)
    }
    let bytes = stored.as_bytes();
    // left trim
    let mut start = 0;
    while start != bytes.len() && is_space(bytes[start]) {
        start += 1;
    }
    // EOL/EOF ends the string (scan the remaining field only)
    let mut end = bytes.len();
    for (i, &b) in bytes[start..].iter().enumerate() {
        if b == b'\n' || b == b'\r' || b == 0 {
            end = start + i;
            break;
        }
    }
    // right trim
    while end != start && is_space(bytes[end - 1]) {
        end -= 1;
    }
    debug_assert!(stored.is_ascii());
    std::str::from_utf8(&bytes[start..end]).unwrap_or("")
}

/// Exact PDB four-column atom-name logical view (BIO-ROWS R01).
///
/// Borrowed helper owning the name conversion that Gemmi performs at READ
/// time via `read_string`: only PDB-provenance documents carry raw
/// four-column names; every CIF-family format stores the decoded logical
/// name that must be preserved byte-for-byte. The IO writer delegates
/// here instead of maintaining a duplicate.
#[must_use]
pub fn atom_name_logical_view(
    name: &crate::source_ids::AtomName,
    input_format: BioCoordinateFormat,
) -> &str {
    // Gemmi✔️✔️: atom.name = read_string(line+12, 4);          // PDB reader
    // Gemmi✔️✔️: a.atom_name = row.str(kLabelAtomId+i);         // mmCIF reader
    // Behavior: Gemmi's stored name for PDB input is the read_string
    // view of the raw four columns (is_space bytes 9-13/32 trimmed at the
    // edges, first NUL/LF/CR terminating), while CIF-family readers store
    // the decoded logical name verbatim. CK's PDB reader keeps the RAW
    // four columns, so this view applies the exact read_string sequence
    // (shared private `pdb_read_string_view`, which carries the verbatim
    // anchors) only for Pdb provenance; every other format (and Unknown,
    // which never has PDB column padding) preserves the stored bytes
    // verbatim. No length/shape inference: input_format is the
    // authoritative per-document discriminator.
    // Complexity: O(1) over at most four bytes, borrowed output — no
    // allocation.
    if input_format == BioCoordinateFormat::Pdb {
        pdb_read_string_view(name.as_str())
    } else {
        name.as_str()
    }
}

/// Exact residue-name logical view (BIO-ROWS R02): the same provenance
/// rule as [`atom_name_logical_view`], applied to the stored residue
/// name. CK's PDB reader already stores the trimmed residue name, so the
/// Pdb branch is the identity there, but the rule is applied explicitly;
/// CIF-family logical names (e.g. quoted `' A '`) preserve their spaces
/// byte-for-byte.
#[must_use]
pub fn residue_name_logical_view(
    name: &crate::source_ids::ResidueName,
    input_format: BioCoordinateFormat,
) -> &str {
    // Gemmi✔️✔️: return {read_seq_id(seq_id), {}, read_string(name, 3)};  // PDB reader
    // Gemmi✔️✔️: vv.emplace_back(cif::quote(res.name));                   // writer
    // Behavior: identical provenance discriminator and the same exact
    // read_string sequence (shared private `pdb_read_string_view`);
    // borrowed output, O(1) over at most four bytes.
    if input_format == BioCoordinateFormat::Pdb {
        pdb_read_string_view(name.as_str())
    } else {
        name.as_str()
    }
}

#[cfg(test)]
mod bio_rows_namefix_tests {
    use super::{BioCoordinateFormat, atom_name_logical_view, residue_name_logical_view};
    use crate::source_ids::{AtomName, ResidueName};

    fn aname(bytes: &[u8]) -> AtomName {
        AtomName::from_ascii(bytes).unwrap()
    }
    fn rname(bytes: &[u8]) -> ResidueName {
        ResidueName::from_ascii(bytes).unwrap()
    }

    #[test]
    fn bio_rows_namefix_all_whitespace_bytes_at_each_edge() {
        // is_space table (atox.hpp) marks bytes 9-13 and 32; the left and
        // right trim loops consume exactly those at the field edges.
        let pdb = BioCoordinateFormat::Pdb;
        for ws in [9u8, 10, 11, 12, 13, 32] {
            let w = [ws];
            // leading: trim consumes it, no terminator, no right trim
            let leading = [&w[..], b"CA"].concat();
            assert_eq!(
                atom_name_logical_view(&aname(&leading), pdb),
                "CA",
                "leading byte {ws}"
            );
            // trailing: bytes 9/11/12 fall to the right trim; 10/13 hit the
            // termination scan first — the pinned source yields "CA" either
            // way (field cut at the byte, nothing after it).
            let trailing = [b"CA".as_slice(), &w].concat();
            assert_eq!(
                atom_name_logical_view(&aname(&trailing), pdb),
                "CA",
                "trailing byte {ws}"
            );
            assert_eq!(
                residue_name_logical_view(&rname(&leading), pdb),
                "CA",
                "residue leading byte {ws}"
            );
            assert_eq!(
                residue_name_logical_view(&rname(&trailing), pdb),
                "CA",
                "residue trailing byte {ws}"
            );
        }
        // A field made only of whitespace bytes trims to empty.
        assert_eq!(atom_name_logical_view(&aname(b"\t\n\x0b\x0c"), pdb), "");
        assert_eq!(residue_name_logical_view(&rname(b" \t\r "), pdb), "");
        // Full-width name with no whitespace survives intact.
        assert_eq!(atom_name_logical_view(&aname(b"ABCD"), pdb), "ABCD");
    }

    #[test]
    fn bio_rows_namefix_interior_terminators_and_tab_family() {
        let pdb = BioCoordinateFormat::Pdb;
        // The EOL/EOF scan cuts at the FIRST NUL/LF/CR in the remaining
        // field; everything after it is discarded.
        assert_eq!(atom_name_logical_view(&aname(b"AB\nC"), pdb), "AB");
        assert_eq!(atom_name_logical_view(&aname(b"A\rBC"), pdb), "A");
        assert_eq!(atom_name_logical_view(&aname(b"A\0BC"), pdb), "A");
        assert_eq!(residue_name_logical_view(&rname(b"AL\nA"), pdb), "AL");
        // NUL as the first byte of the remaining field yields the empty
        // string (left trim does not consume NUL; is_space(0) == 0).
        assert_eq!(atom_name_logical_view(&aname(b"\0ABC"), pdb), "");
        // Interior TAB/VT/FF are is_space but NOT terminators: interior
        // whitespace survives because the trims only touch the edges.
        assert_eq!(atom_name_logical_view(&aname(b"A\tB"), pdb), "A\tB");
        assert_eq!(atom_name_logical_view(&aname(b"A\x0bB"), pdb), "A\x0bB");
        assert_eq!(residue_name_logical_view(&rname(b"A\x0cB"), pdb), "A\x0cB");
        // A leading LF is whitespace for the LEFT trim (is_space) and is
        // consumed there — termination only applies to interior bytes.
        assert_eq!(atom_name_logical_view(&aname(b"\n CA"), pdb), "CA");
        // Leading trim advancing past whitespace keeps a later terminator
        // interior: "\tA\nC" -> field "A\nC" -> cut at LF -> "A".
        assert_eq!(atom_name_logical_view(&aname(b"\tA\nC"), pdb), "A");
    }

    #[test]
    fn bio_rows_namefix_non_pdb_preserves_these_bytes_verbatim() {
        // CIF-family/Unknown provenance never applies read_string: the very
        // same stored bytes are returned untouched.
        for format in [
            BioCoordinateFormat::Unknown,
            BioCoordinateFormat::Detect,
            BioCoordinateFormat::Mmcif,
            BioCoordinateFormat::Mmjson,
            BioCoordinateFormat::ChemComp,
        ] {
            for stored in [&b"\tCA"[..], b"CA\n", b"AB\0C", b"A\tB", b"\n CA", b" CA "] {
                assert_eq!(
                    atom_name_logical_view(&aname(stored), format),
                    std::str::from_utf8(stored).unwrap(),
                    "{format:?} atom {stored:?}"
                );
                assert_eq!(
                    residue_name_logical_view(&rname(stored), format),
                    std::str::from_utf8(stored).unwrap(),
                    "{format:?} residue {stored:?}"
                );
            }
        }
    }
}

#[cfg(test)]
mod bio_rows_r02_tests {
    use super::{BioCoordinateFormat, residue_name_logical_view};
    use crate::source_ids::ResidueName;

    fn name(bytes: &[u8]) -> ResidueName {
        ResidueName::from_ascii(bytes).unwrap()
    }

    #[test]
    fn bio_rows_r02_pdb_padding_versus_cif_spaces() {
        let pdb = BioCoordinateFormat::Pdb;
        // PDB three-column provenance: read_string trims outer spaces.
        assert_eq!(residue_name_logical_view(&name(b"ALA "), pdb), "ALA");
        assert_eq!(residue_name_logical_view(&name(b" ALA"), pdb), "ALA");
        assert_eq!(residue_name_logical_view(&name(b"   "), pdb), "");
        // CIF-family logical names keep literal spaces byte-for-byte.
        for format in [
            BioCoordinateFormat::Unknown,
            BioCoordinateFormat::Detect,
            BioCoordinateFormat::Mmcif,
            BioCoordinateFormat::Mmjson,
            BioCoordinateFormat::ChemComp,
        ] {
            assert_eq!(
                residue_name_logical_view(&name(b" ALA"), format),
                " ALA",
                "{format:?}"
            );
            assert_eq!(
                residue_name_logical_view(&name(b"   "), format),
                "   ",
                "{format:?}"
            );
        }
    }

    #[test]
    fn bio_rows_r02_blank_and_maximal_stored_names() {
        let pdb = BioCoordinateFormat::Pdb;
        // Blank stored name stays blank in both provenances.
        assert_eq!(residue_name_logical_view(&name(b""), pdb), "");
        assert_eq!(
            residue_name_logical_view(&name(b""), BioCoordinateFormat::Mmcif),
            ""
        );
        // Maximal four-byte stored names keep every byte except the PDB
        // outer-space trim; embedded spaces are preserved in both.
        assert_eq!(residue_name_logical_view(&name(b"MSE "), pdb), "MSE");
        assert_eq!(
            residue_name_logical_view(&name(b"AB D"), BioCoordinateFormat::Mmcif),
            "AB D"
        );
        assert_eq!(residue_name_logical_view(&name(b"A B "), pdb), "A B");
        // An all-space CIF name is a literal logical name, not blanked.
        assert_eq!(
            residue_name_logical_view(&name(b" "), BioCoordinateFormat::Mmcif),
            " "
        );
    }
}

#[cfg(test)]
mod bio_rows_r01_tests {
    use super::{BioCoordinateFormat, atom_name_logical_view};
    use crate::source_ids::AtomName;

    fn name(bytes: &[u8]) -> AtomName {
        AtomName::from_ascii(bytes).unwrap()
    }

    #[test]
    fn bio_rows_r01_pdb_padding_rules() {
        // read_string trims leading and trailing spaces of the raw four
        // columns (pdb.cpp read_string), independently derived from the
        // pinned left-trim/EOL/right-trim sequence.
        let pdb = BioCoordinateFormat::Pdb;
        assert_eq!(atom_name_logical_view(&name(b" CA "), pdb), "CA");
        assert_eq!(atom_name_logical_view(&name(b"CA  "), pdb), "CA");
        assert_eq!(atom_name_logical_view(&name(b"  CA"), pdb), "CA");
        assert_eq!(atom_name_logical_view(&name(b"    "), pdb), "");
        assert_eq!(atom_name_logical_view(&name(b"CA"), pdb), "CA");
        // Embedded spaces are NOT trimmed: only the outer run.
        assert_eq!(atom_name_logical_view(&name(b"C A "), pdb), "C A");
        assert_eq!(atom_name_logical_view(&name(b" C A"), pdb), "C A");
        // A name that is entirely spaces between letters stays verbatim.
        assert_eq!(atom_name_logical_view(&name(b"HB1 "), pdb), "HB1");
    }

    #[test]
    fn bio_rows_r01_non_pdb_preserves_bytes() {
        // Every CIF-family format and Unknown/Detect preserve the stored
        // logical name byte-for-byte — no trims at all.
        for format in [
            BioCoordinateFormat::Unknown,
            BioCoordinateFormat::Detect,
            BioCoordinateFormat::Mmcif,
            BioCoordinateFormat::Mmjson,
            BioCoordinateFormat::ChemComp,
        ] {
            assert_eq!(
                atom_name_logical_view(&name(b" CA "), format),
                " CA ",
                "{format:?}"
            );
            assert_eq!(
                atom_name_logical_view(&name(b"  "), format),
                "  ",
                "{format:?}"
            );
            assert_eq!(
                atom_name_logical_view(&name(b"C A"), format),
                "C A",
                "{format:?}"
            );
        }
    }
}

#[cfg(test)]
mod bio_legacy_n01_tests {
    use super::{
        BioCoordinateBlock, BioCoordinateFormat, BioMetadata, BioStructureData, BioStructureParts,
        BioStructureSourceState,
    };

    #[test]
    fn bio_legacy_n01_identity_id_borrows_existing_source_map() {
        let mut parts = BioStructureParts {
            input_format: BioCoordinateFormat::Mmcif,
            models: vec![],
            chains: vec![],
            residues: vec![],
            atoms: vec![],
            entities: vec![],
            connections: vec![],
            cispeps: vec![],
            mod_residues: vec![],
            helices: vec![],
            sheets: vec![],
            metadata: BioMetadata::default(),
            source_state: BioStructureSourceState::default(),
            coordinates: BioCoordinateBlock::default(),
            crystal: None,
            ncs_operators: vec![],
            assemblies: vec![],
        };
        assert_eq!(
            BioStructureData::from_parts(parts.clone())
                .unwrap()
                .ncs_oper_identity_id(),
            None
        );
        for value in ["", "identity-last"] {
            parts
                .source_state
                .info
                .insert("_struct_ncs_oper.id".into(), value.into());
            parts
                .source_state
                .info
                .insert("unrelated".into(), "kept".into());
            let data = BioStructureData::from_parts(parts.clone()).unwrap();
            assert_eq!(data.ncs_oper_identity_id(), Some(value));
            assert_eq!(data.source_state().info["_struct_ncs_oper.id"], value);
            assert_eq!(data.source_state().info["unrelated"], "kept");
            assert_eq!(
                data.ncs_oper_identity_id().unwrap().as_ptr(),
                data.source_state().info["_struct_ncs_oper.id"].as_ptr()
            );
        }
    }
}

#[cfg(test)]
mod crystal_transition_tests {
    use super::{
        BioCrystalCell, BioCrystalInfo, BioStructureError, BioTransform, inverse_bio_matrix,
        inverse_bio_transform, set_crystal_cell, set_crystal_fractional_transform,
        set_crystal_space_group_hm, set_crystal_z_pdb_if_nonempty,
    };

    fn crystal_with_matrices(
        fractional: BioTransform,
        orthogonal: BioTransform,
        explicit_matrices: bool,
    ) -> BioCrystalInfo {
        BioCrystalInfo::new(
            BioCrystalCell::default(),
            Some("P 1".to_owned()),
            Some("1".to_owned()),
            orthogonal,
            fractional,
            explicit_matrices,
            0,
            Vec::new(),
        )
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1.0e-15,
            "actual {actual:.17e} differs from expected {expected:.17e}"
        );
    }

    #[test]
    fn crystal_transition_inverse_nonsymmetric_matches_gemmi_cofactors() {
        let matrix = [[2.0, 1.0, 0.0], [0.0, 3.0, 1.0], [1.0, 0.0, 4.0]];
        let expected = [
            [12.0 / 25.0, -4.0 / 25.0, 1.0 / 25.0],
            [1.0 / 25.0, 8.0 / 25.0, -2.0 / 25.0],
            [-3.0 / 25.0, 1.0 / 25.0, 6.0 / 25.0],
        ];

        let actual = inverse_bio_matrix(matrix);
        for row in 0..3 {
            for column in 0..3 {
                assert_close(actual[row][column], expected[row][column]);
            }
        }
    }

    #[test]
    fn crystal_transition_inverse_translation_uses_inverse_matrix() {
        let transform = BioTransform::new(
            [[2.0, 1.0, 0.0], [0.0, 3.0, 1.0], [1.0, 0.0, 4.0]],
            [1.0, 2.0, -3.0],
        );

        let actual = inverse_bio_transform(transform);
        for (value, expected) in
            actual
                .translation()
                .iter()
                .zip([-1.0 / 25.0, -23.0 / 25.0, 19.0 / 25.0])
        {
            assert_close(*value, expected);
        }
    }

    #[test]
    fn crystal_transition_identity_inverse_is_identity() {
        let actual = inverse_bio_transform(BioTransform::identity());
        assert_eq!(actual, BioTransform::identity());
        assert_eq!(
            inverse_bio_matrix([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
    }

    #[test]
    fn crystal_transition_singular_inverse_preserves_ieee_non_finite_values() {
        let actual = inverse_bio_matrix([[0.0; 3]; 3]);
        assert!(actual.into_iter().flatten().all(f64::is_nan));
        let transform = inverse_bio_transform(BioTransform::new([[0.0; 3]; 3], [1.0, -2.0, 3.0]));
        assert!(
            transform
                .matrix()
                .iter()
                .flatten()
                .all(|value| value.is_nan())
        );
        assert!(transform.translation().iter().all(|value| value.is_nan()));
    }

    #[test]
    fn crystal_transition_cell_gamma_zero_is_an_exact_noop() {
        let transform = BioTransform::new(
            [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            [5.0, 6.0, 7.0],
        );
        let mut crystal = BioCrystalInfo::new(
            BioCrystalCell {
                a: 4.0,
                b: 5.0,
                c: 6.0,
                alpha: 80.0,
                beta: 95.0,
                gamma: 100.0,
            },
            Some("P 21 21 21".to_owned()),
            Some("8".to_owned()),
            transform,
            transform,
            true,
            3,
            vec![BioTransform::identity()],
        );
        crystal.calculate_properties().unwrap();
        let before = crystal.clone();

        assert_eq!(
            set_crystal_cell(
                &mut crystal,
                BioCrystalCell {
                    a: 0.0,
                    b: 0.0,
                    c: 0.0,
                    alpha: 0.0,
                    beta: 0.0,
                    gamma: -0.0,
                }
            ),
            Ok(())
        );
        assert_eq!(crystal, before);
    }

    #[test]
    fn crystal_transition_cell_replacement_refreshes_orthogonal_properties() {
        let mut crystal = BioCrystalInfo::new(
            BioCrystalCell::default(),
            None,
            None,
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            0,
            Vec::new(),
        );
        crystal.calculate_properties().unwrap();

        let replacement = BioCrystalCell {
            a: 4.0,
            b: 5.0,
            c: 6.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        };
        set_crystal_cell(&mut crystal, replacement).unwrap();

        assert_eq!(crystal.cell(), replacement);
        assert_eq!(crystal.volume(), 120.0);
        assert_eq!(*crystal.reciprocal_lengths(), [0.25, 0.2, 1.0 / 6.0]);
        assert_eq!(*crystal.reciprocal_cosines(), [0.0; 3]);
        assert_eq!(
            *crystal.orthogonal().matrix(),
            [[4.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 6.0]]
        );
    }

    #[test]
    fn crystal_transition_cell_error_keeps_source_ordered_partial_state() {
        let mut crystal = BioCrystalInfo::new(
            BioCrystalCell::default(),
            Some("P 1".to_owned()),
            Some("1".to_owned()),
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            2,
            Vec::new(),
        );
        crystal.calculate_properties().unwrap();
        let old_volume = crystal.volume();
        let old_lengths = *crystal.reciprocal_lengths();
        let old_cosines = *crystal.reciprocal_cosines();
        let old_orthogonal = *crystal.orthogonal();
        let old_fractional = *crystal.fractional();
        let invalid = BioCrystalCell {
            a: 2.0,
            b: 3.0,
            c: 4.0,
            alpha: 0.0,
            beta: 91.0,
            gamma: 92.0,
        };

        assert_eq!(
            set_crystal_cell(&mut crystal, invalid),
            Err(BioStructureError::ImpossibleCrystalAngle)
        );
        assert_eq!(crystal.cell(), invalid);
        assert_eq!(crystal.volume(), old_volume);
        assert_eq!(*crystal.reciprocal_lengths(), old_lengths);
        assert_eq!(*crystal.reciprocal_cosines(), old_cosines);
        assert_eq!(*crystal.orthogonal(), old_orthogonal);
        assert_eq!(*crystal.fractional(), old_fractional);
        assert_eq!(crystal.space_group_hm(), Some("P 1"));
        assert_eq!(crystal.z_pdb(), Some("1"));
        assert_eq!(crystal.cs_count(), 2);
    }

    #[test]
    fn crystal_transition_cell_refreshes_scalars_without_replacing_explicit_matrices() {
        let explicit = BioTransform::new(
            [[7.0, 1.0, 2.0], [3.0, 8.0, 4.0], [5.0, 6.0, 9.0]],
            [10.0, 11.0, 12.0],
        );
        let mut crystal = BioCrystalInfo::new(
            BioCrystalCell {
                a: 2.0,
                b: 3.0,
                c: 4.0,
                alpha: 80.0,
                beta: 90.0,
                gamma: 100.0,
            },
            Some("P 1".to_owned()),
            Some("4".to_owned()),
            explicit,
            explicit,
            true,
            1,
            Vec::new(),
        );
        crystal.calculate_properties().unwrap();

        set_crystal_cell(
            &mut crystal,
            BioCrystalCell {
                a: 5.0,
                b: 6.0,
                c: 7.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            },
        )
        .unwrap();

        assert!(crystal.explicit_matrices());
        assert_eq!(*crystal.orthogonal(), explicit);
        assert_eq!(*crystal.fractional(), explicit);
        assert_eq!(crystal.volume(), 210.0);
        assert_eq!(*crystal.reciprocal_lengths(), [0.2, 1.0 / 6.0, 1.0 / 7.0]);
        assert_eq!(*crystal.reciprocal_cosines(), [0.0; 3]);
        assert_eq!(crystal.space_group_hm(), Some("P 1"));
        assert_eq!(crystal.z_pdb(), Some("4"));
        assert_eq!(crystal.cs_count(), 1);
    }

    #[test]
    fn crystal_transition_fractional_thresholds_are_inclusive() {
        let orthogonal = BioTransform::new(
            [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            [5.0, 6.0, 7.0],
        );
        let mut crystal = crystal_with_matrices(BioTransform::identity(), orthogonal, false);
        let before = crystal.clone();
        let mut just_inside = BioTransform::identity();
        just_inside.matrix[0][1] = 0.999e-4;
        just_inside.translation[0] = 0.999e-6;
        set_crystal_fractional_transform(&mut crystal, just_inside);
        assert_eq!(crystal, before);

        let mut at_threshold = BioTransform::identity();
        at_threshold.matrix[0][1] = 1.0e-4;
        at_threshold.translation[0] = 1.0e-6;
        set_crystal_fractional_transform(&mut crystal, at_threshold);
        assert_eq!(crystal, before);
    }

    #[test]
    fn crystal_transition_fractional_differences_above_each_threshold_are_installed() {
        let mut matrix_above = BioTransform::identity();
        matrix_above.matrix[0][1] = 1.0001e-4;
        let mut crystal =
            crystal_with_matrices(BioTransform::identity(), BioTransform::identity(), false);
        set_crystal_fractional_transform(&mut crystal, matrix_above);
        assert_eq!(*crystal.fractional(), matrix_above);
        assert!(crystal.explicit_matrices());

        let translation_above =
            BioTransform::new(*BioTransform::identity().matrix(), [1.0001e-6, 0.0, 0.0]);
        let mut crystal =
            crystal_with_matrices(BioTransform::identity(), BioTransform::identity(), false);
        set_crystal_fractional_transform(&mut crystal, translation_above);
        assert_eq!(*crystal.fractional(), translation_above);
        assert_eq!(crystal.orthogonal().translation(), &[-1.0001e-6, 0.0, 0.0]);
        assert!(crystal.explicit_matrices());
    }

    #[test]
    fn crystal_transition_fractional_matrix_nan_difference_is_source_approximate() {
        let orthogonal = BioTransform::new(
            [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            [5.0, 6.0, 7.0],
        );
        let mut crystal = crystal_with_matrices(BioTransform::identity(), orthogonal, false);
        let before = crystal.clone();
        let transform = BioTransform::new(
            [[1.0, f64::NAN, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [0.0; 3],
        );

        set_crystal_fractional_transform(&mut crystal, transform);

        assert_eq!(crystal, before);
    }

    #[test]
    fn crystal_transition_fractional_vector_nan_is_not_source_approximate() {
        let mut crystal =
            crystal_with_matrices(BioTransform::identity(), BioTransform::identity(), false);
        let transform = BioTransform::new(*BioTransform::identity().matrix(), [f64::NAN, 0.0, 0.0]);

        set_crystal_fractional_transform(&mut crystal, transform);

        assert!(crystal.fractional().translation()[0].is_nan());
        assert!(crystal.orthogonal().translation()[0].is_nan());
        assert!(crystal.explicit_matrices());
    }

    #[test]
    fn crystal_transition_fractional_suspicious_first_matrix_values_are_ignored() {
        for first_value in [0.0, 1.25] {
            let orthogonal = BioTransform::new(
                [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
                [5.0, 6.0, 7.0],
            );
            let mut crystal = crystal_with_matrices(BioTransform::identity(), orthogonal, false);
            let before = crystal.clone();
            let transform = BioTransform::new(
                [[first_value, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.0; 3],
            );

            set_crystal_fractional_transform(&mut crystal, transform);

            assert_eq!(crystal, before);
        }
    }

    #[test]
    fn crystal_transition_fractional_replacement_inverts_translation_and_sets_explicit_flag() {
        let initial_fractional = BioTransform::new(
            [[0.5, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [0.0; 3],
        );
        let mut crystal =
            crystal_with_matrices(initial_fractional, BioTransform::identity(), false);
        let replacement = BioTransform::new(
            [[2.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 5.0]],
            [2.0, 4.0, 10.0],
        );

        set_crystal_fractional_transform(&mut crystal, replacement);

        assert_eq!(*crystal.fractional(), replacement);
        assert_eq!(
            *crystal.orthogonal(),
            BioTransform::new(
                [[0.5, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.2]],
                [-1.0, -1.0, -2.0]
            )
        );
        assert!(crystal.explicit_matrices());
        assert_eq!(crystal.cell(), BioCrystalCell::default());
        assert_eq!(crystal.space_group_hm(), Some("P 1"));
        assert_eq!(crystal.z_pdb(), Some("1"));
    }

    #[test]
    fn crystal_transition_space_group_replaces_only_its_metadata() {
        let orthogonal = BioTransform::new(
            [[7.0, 1.0, 2.0], [3.0, 8.0, 4.0], [5.0, 6.0, 9.0]],
            [10.0, 11.0, 12.0],
        );
        let fractional = BioTransform::new(
            [[0.5, 0.1, 0.2], [0.3, 0.4, 0.6], [0.7, 0.8, 0.9]],
            [-1.0, -2.0, -3.0],
        );
        let mut crystal = BioCrystalInfo::new(
            BioCrystalCell {
                a: 6.0,
                b: 7.0,
                c: 8.0,
                alpha: 90.0,
                beta: 100.0,
                gamma: 110.0,
            },
            Some("old-group".to_owned()),
            Some("17".to_owned()),
            orthogonal,
            fractional,
            true,
            9,
            vec![BioTransform::identity(), orthogonal],
        );
        crystal.calculate_properties().unwrap();

        let before = crystal.clone();
        set_crystal_space_group_hm(&mut crystal, "new-group".to_owned());
        let mut expected = before.clone();
        expected.space_group_hm = Some("new-group".to_owned());
        assert_eq!(crystal, expected);

        let before_empty_replacement = crystal.clone();
        set_crystal_space_group_hm(&mut crystal, String::new());
        let mut expected_empty = before_empty_replacement;
        expected_empty.space_group_hm = Some(String::new());
        assert_eq!(crystal, expected_empty);
    }

    #[test]
    fn crystal_transition_z_replaces_nonempty_and_retains_empty() {
        let mut crystal = BioCrystalInfo::new(
            BioCrystalCell {
                a: 5.0,
                b: 6.0,
                c: 7.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            },
            Some("P 21 21 21".to_owned()),
            Some("8".to_owned()),
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            4,
            vec![BioTransform::identity()],
        );
        crystal.calculate_properties().unwrap();

        let before = crystal.clone();
        set_crystal_z_pdb_if_nonempty(&mut crystal, "12".to_owned());
        let mut expected = before;
        expected.z_pdb = Some("12".to_owned());
        assert_eq!(crystal, expected);

        let before_empty = crystal.clone();
        set_crystal_z_pdb_if_nonempty(&mut crystal, String::new());
        assert_eq!(crystal, before_empty);
    }

    #[test]
    fn crystal_transition_empty_z_keeps_absent_metadata_absent() {
        let mut crystal = BioCrystalInfo::new(
            BioCrystalCell::default(),
            None,
            None,
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            0,
            Vec::new(),
        );
        crystal.calculate_properties().unwrap();
        let before = crystal.clone();

        set_crystal_z_pdb_if_nonempty(&mut crystal, String::new());

        assert_eq!(crystal, before);
        assert_eq!(crystal.z_pdb(), None);
    }
}
