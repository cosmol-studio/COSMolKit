//! Source-address values used to relate BioStructure records.

use crate::{PdbChainId, ResidueName};

/// A residue's source address, distinct from its stable [`crate::BioResidueId`].
///
/// Sequence-number absence and insertion code are independent, matching
/// Gemmi's `SeqId`; segment and residue-name bytes retain their source spelling
/// within COSMolKit's approved four-byte ASCII boundary.
#[derive(Debug, Clone, Copy)]
pub struct ResidueAddress {
    sequence_number: Option<i32>,
    insertion_code: Option<u8>,
    segment: [u8; 4],
    segment_len: u8,
    name: ResidueName,
}

impl ResidueAddress {
    /// Construct an address, returning `None` when the segment exceeds the
    /// approved four-byte ASCII representation.
    #[must_use]
    pub fn new(
        sequence_number: Option<i32>,
        insertion_code: Option<u8>,
        segment: &[u8],
        name: ResidueName,
    ) -> Option<Self> {
        // Gemmi❗✔️: struct ResidueId {
        // Gemmi❗✔️:   SeqId seqid;
        // Gemmi❗✔️:   std::string segment; // segid - up to 4 characters in the PDB file
        // Gemmi❗✔️:   std::string name;
        // Gemmi❗✔️: };
        // Behavior review: the source owns unrestricted strings; the approved
        // model boundary retains exact empty-through-four-byte ASCII values.
        // Complexity review: one bounded validation and fixed-array copy.
        if segment.len() > 4 || !segment.is_ascii() {
            return None;
        }
        let mut stored_segment = [0; 4];
        stored_segment[..segment.len()].copy_from_slice(segment);

        // Gemmi❗✔️: using OptionalNum = OptionalInt<INT_MIN>;
        // Gemmi❗✔️: OptionalNum num;   // sequence number
        // Gemmi❗✔️: char icode = ' ';  // insertion code
        // Behavior review: `INT_MIN` is Gemmi's unset sentinel and space is
        // its unset insertion code, so equivalent Rust inputs are canonicalized.
        // Complexity review: constant-time scalar normalization.
        Some(Self {
            sequence_number: sequence_number.filter(|number| *number != i32::MIN),
            insertion_code: insertion_code.filter(|code| *code != b' '),
            segment: stored_segment,
            segment_len: segment.len() as u8,
            name,
        })
    }

    #[must_use]
    pub const fn sequence_number(self) -> Option<i32> {
        self.sequence_number
    }

    #[must_use]
    pub const fn insertion_code(self) -> Option<u8> {
        self.insertion_code
    }

    #[must_use]
    pub fn segment(&self) -> &str {
        // Construction accepts ASCII only, so the stored slice is valid UTF-8.
        std::str::from_utf8(&self.segment[..self.segment_len as usize])
            .expect("ResidueAddress segment invariant")
    }

    #[must_use]
    pub const fn name(self) -> ResidueName {
        self.name
    }

    /// Match sequence, insertion code, segment, and residue name using
    /// Gemmi's `ResidueId::matches` semantics.
    #[must_use]
    pub fn matches(&self, other: &Self) -> bool {
        // Gemmi✔️✔️: bool matches(const ResidueId& o) const {
        // Gemmi✔️✔️:   return seqid == o.seqid && segment == o.segment && name == o.name;
        // Gemmi✔️✔️: }
        // Behavior review: equality retains the source's case-insensitive
        // insertion-code comparison and exact segment/name comparison.
        // Complexity review: all fields are bounded to four bytes; no scans
        // or allocations grow with unbounded input.
        source_seq_id_matches(
            self.sequence_number,
            self.insertion_code,
            other.sequence_number,
            other.insertion_code,
        ) && self.segment() == other.segment()
            && self.name == other.name
    }

    /// Match sequence, insertion code, and residue name, ignoring segment.
    #[must_use]
    pub fn matches_without_segment(&self, other: &Self) -> bool {
        // Gemmi✔️✔️: bool matches_noseg(const ResidueId& o) const {
        // Gemmi✔️✔️:   return seqid == o.seqid && name == o.name;
        // Gemmi✔️✔️: }
        // Behavior review: only the source's segment comparison is omitted.
        // Complexity review: constant bounded scalar and name comparisons.
        source_seq_id_matches(
            self.sequence_number,
            self.insertion_code,
            other.sequence_number,
            other.insertion_code,
        ) && self.name == other.name
    }
}

impl PartialEq for ResidueAddress {
    fn eq(&self, other: &Self) -> bool {
        // Gemmi✔️✔️: bool operator==(const ResidueId& o) const { return matches(o); }
        // Behavior review: `matches` is the unique implementation of source
        // address equality. Complexity review: one bounded comparison pass.
        self.matches(other)
    }
}

impl Eq for ResidueAddress {}

/// A logical atom address, distinct from the fixed four-column [`crate::AtomName`].
#[derive(Debug, Clone)]
pub struct AtomAddress {
    chain_name: PdbChainId,
    residue: ResidueAddress,
    logical_atom_name: String,
    altloc: u8,
}

impl AtomAddress {
    /// Construct a logical address. `None` selects Gemmi's default NUL altloc.
    #[must_use]
    pub fn new(
        chain_name: PdbChainId,
        residue: ResidueAddress,
        logical_atom_name: impl Into<String>,
        altloc: Option<u8>,
    ) -> Self {
        // Gemmi✔️✔️: AtomAddress(const std::string& ch, const ResidueId& resid,
        // Gemmi✔️✔️:               const std::string& atom, char alt='\0')
        // Gemmi✔️✔️:   : chain_name(ch), res_id(resid), atom_name(atom), altloc(alt) {}
        // Behavior review: logical atom-name text is retained exactly; no
        // four-column padding is inferred. `None` maps only to source altloc NUL.
        // Complexity review: the source copies its strings; Rust performs one
        // owned atom-name conversion and keeps the bounded value fields inline.
        Self {
            chain_name,
            residue,
            logical_atom_name: logical_atom_name.into(),
            altloc: altloc.unwrap_or(0),
        }
    }

    #[must_use]
    pub const fn chain_name(&self) -> PdbChainId {
        self.chain_name
    }

    #[must_use]
    pub const fn residue(&self) -> ResidueAddress {
        self.residue
    }

    #[must_use]
    pub fn logical_atom_name(&self) -> &str {
        &self.logical_atom_name
    }

    #[must_use]
    pub const fn altloc(&self) -> u8 {
        self.altloc
    }
}

impl PartialEq for AtomAddress {
    fn eq(&self, other: &Self) -> bool {
        // Gemmi✔️✔️: bool operator==(const AtomAddress& o) const {
        // Gemmi✔️✔️:   return chain_name == o.chain_name && res_id.matches(o.res_id) &&
        // Gemmi✔️✔️:          atom_name == o.atom_name && altloc == o.altloc;
        // Gemmi✔️✔️: }
        // Behavior review: chain, logical atom name, and altloc compare exactly;
        // residue equality delegates to the source-shaped ResidueId matcher.
        // Complexity review: bounded chain/residue comparisons plus string
        // comparisons linear in the logical atom-name length, as for std::string.
        self.chain_name == other.chain_name
            && self.residue.matches(&other.residue)
            && self.logical_atom_name == other.logical_atom_name
            && self.altloc == other.altloc
    }
}

impl Eq for AtomAddress {}

impl Default for AtomAddress {
    fn default() -> Self {
        // Gemmi✔️✔️: AtomAddress() = default;
        // Gemmi✔️✔️: ResidueId res_id;
        // Gemmi✔️✔️: std::string atom_name;
        // Gemmi✔️✔️: char altloc = '\0';
        // Behavior review: the default chain/name strings are empty; its
        // default residue has no sequence, blank insertion code, empty segment
        // and name. The source NUL altloc is represented as byte zero.
        // Complexity review: constant-size value construction; empty Strings
        // allocate no heap storage.
        let empty_chain = PdbChainId::from_ascii(b"")
            .expect("empty chain id is within the approved representation");
        let empty_residue_name = ResidueName::from_ascii(b"")
            .expect("empty residue name is within the approved representation");
        let empty_residue = ResidueAddress::new(None, None, b"", empty_residue_name)
            .expect("empty residue address is within the approved representation");
        Self {
            chain_name: empty_chain,
            residue: empty_residue,
            logical_atom_name: String::new(),
            altloc: 0,
        }
    }
}

/// Relationship kind retained from Gemmi's `_struct_conn` model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum BioConnectionKind {
    Covale = 0,
    Disulf,
    Hydrog,
    MetalC,
    Unknown,
}

impl Default for BioConnectionKind {
    fn default() -> Self {
        // Gemmi✔️✔️: enum Type : unsigned char { Covale=0, Disulf, Hydrog, MetalC, Unknown };
        // Gemmi✔️✔️: Type type = Unknown;
        // Behavior review: the source ordinal and default sentinel are retained.
        // Complexity review: constant-time scalar selection.
        Self::Unknown
    }
}

/// Unit-cell ASU relation for a structural connection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum BioAsu {
    Same = 0,
    Different,
    Any,
}

impl Default for BioAsu {
    fn default() -> Self {
        // Gemmi✔️✔️: enum class Asu : unsigned char { Same, Different, Any };
        // Gemmi✔️✔️: Asu asu = Asu::Any;
        // Behavior review: variants retain source ordinals and Any is the default.
        // Complexity review: constant-time scalar selection.
        Self::Any
    }
}

/// Typed `_struct_conn` relationship with source atom addresses.
#[derive(Debug, Clone, PartialEq)]
pub struct BioConnection {
    pub name: String,
    pub link_id: String,
    pub kind: BioConnectionKind,
    pub asu: BioAsu,
    pub partner1: AtomAddress,
    pub partner2: AtomAddress,
    pub reported_distance: f64,
    pub reported_sym: [i16; 4],
}

impl Default for BioConnection {
    fn default() -> Self {
        // Gemmi✔️✔️: std::string name;
        // Gemmi✔️✔️: std::string link_id;  // _struct_conn.ccp4_link_id (== _chem_link.id)
        // Gemmi✔️✔️: Type type = Unknown;
        // Gemmi✔️✔️: Asu asu = Asu::Any;
        // Gemmi✔️✔️: AtomAddress partner1, partner2;
        // Gemmi✔️✔️: double reported_distance = 0.0;
        // Gemmi✔️✔️: short reported_sym[4] = {};  // don't rely on it, for internal use only
        // Behavior review: Rust preserves both empty addresses, source default
        // kind/ASU, zero distance and all four zero symmetry components.
        // Complexity review: fixed-size field assembly; Strings and addresses
        // remain empty inline/zero-allocation values.
        Self {
            name: String::new(),
            link_id: String::new(),
            kind: BioConnectionKind::default(),
            asu: BioAsu::default(),
            partner1: AtomAddress::default(),
            partner2: AtomAddress::default(),
            reported_distance: 0.0,
            reported_sym: [0; 4],
        }
    }
}

/// Cis-peptide relationship retained from Gemmi's PDBx/mmCIF metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct BioCisPep {
    pub partner_c: AtomAddress,
    pub partner_n: AtomAddress,
    pub model_num: i32,
    pub only_altloc: u8,
    pub reported_angle: f64,
}

impl Default for BioCisPep {
    fn default() -> Self {
        // Gemmi✔️✔️: struct CisPep {
        // Gemmi✔️✔️:   AtomAddress partner_c, partner_n;
        // Gemmi✔️✔️:   int model_num = 0;
        // Gemmi✔️✔️:   // mmCIF has (unused by the PDB) tag _struct_mon_prot_cis.label_alt_id
        // Gemmi✔️✔️:   // that enables defining CIS link per conformation.
        // Gemmi✔️✔️:   char only_altloc = '\0';
        // Gemmi✔️✔️:   double reported_angle = NAN;
        // Gemmi✔️✔️: };
        // Behavior review: both addresses default independently, the model and
        // altloc retain their source zero values, and the angle is NaN.
        // Complexity review: fixed-field construction with two zero-allocation
        // empty address values.
        Self {
            partner_c: AtomAddress::default(),
            partner_n: AtomAddress::default(),
            model_num: 0,
            only_altloc: 0,
            reported_angle: f64::NAN,
        }
    }
}

/// Modified-residue metadata retained from Gemmi's `ModRes` value.
#[derive(Debug, Clone, PartialEq)]
pub struct BioModRes {
    pub chain_name: PdbChainId,
    pub res_id: ResidueAddress,
    pub parent_comp_id: String,
    pub mod_id: String,
    pub details: String,
}

impl Default for BioModRes {
    fn default() -> Self {
        // Gemmi❗✔️: struct ModRes {
        // Gemmi❗✔️:   std::string chain_name;
        // Gemmi❗✔️:   ResidueId res_id;
        // Gemmi❗✔️:   std::string parent_comp_id;
        // Gemmi❗✔️:   std::string mod_id;  // non-standard extension used in Refmac
        // Gemmi❗✔️:   std::string details;
        // Gemmi❗✔️: };
        // Behavior review: source strings default empty and `ResidueId`
        // default state is represented by absent sequence/insertion, empty
        // segment, and empty name. Chain and residue-name widths intentionally
        // remain within the approved BIO value boundary.
        // Complexity review: empty String construction is allocation-free;
        // bounded identifiers and the residue address are fixed-size values.
        let empty_chain = PdbChainId::from_ascii(b"")
            .expect("empty chain id is within the approved representation");
        let empty_residue_name = ResidueName::from_ascii(b"")
            .expect("empty residue name is within the approved representation");
        let empty_residue = ResidueAddress::new(None, None, b"", empty_residue_name)
            .expect("empty residue address is within the approved representation");
        Self {
            chain_name: empty_chain,
            res_id: empty_residue,
            parent_comp_id: String::new(),
            mod_id: String::new(),
            details: String::new(),
        }
    }
}

fn source_seq_id_matches(
    left_number: Option<i32>,
    left_insertion_code: Option<u8>,
    right_number: Option<i32>,
    right_insertion_code: Option<u8>,
) -> bool {
    // Gemmi✔️✔️: bool operator==(const OptionalInt& o) const { return value == o.value; }
    // Gemmi✔️✔️: bool operator==(const SeqId& o) const {
    // Gemmi✔️✔️:   return num == o.num && ((icode ^ o.icode) & ~0x20) == 0;
    // Gemmi✔️✔️: }
    // Behavior review: constructor normalization maps Gemmi's INT_MIN/space
    // sentinels to `None`; signed-char promotion preserves the pinned GNU
    // profile's bit-0x20 comparison for every insertion byte.
    // Complexity review: fixed scalar comparisons only.
    let left_char = i32::from(i8::from_ne_bytes([left_insertion_code.unwrap_or(b' ')]));
    let right_char = i32::from(i8::from_ne_bytes([right_insertion_code.unwrap_or(b' ')]));
    left_number == right_number && (left_char ^ right_char) & !0x20 == 0
}
