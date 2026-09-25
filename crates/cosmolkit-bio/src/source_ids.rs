//! Immutable source identifiers retained from PDB/mmCIF input.

/// PDB atom serial number. This is source provenance, not a row identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PdbAtomSerial(i32);

impl PdbAtomSerial {
    #[must_use]
    pub const fn new(value: i32) -> Self {
        // Gemmi✔️✔️: int serial = 0;
        Self(value)
    }

    #[must_use]
    pub const fn value(self) -> i32 {
        self.0
    }
}

/// PDB/mmCIF chain identifier within COSMolKit's declared four-byte boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PdbChainId {
    bytes: [u8; 4],
    len: u8,
}

impl PdbChainId {
    /// Constructs an identifier without trimming or case conversion.
    #[must_use]
    pub fn from_ascii(bytes: &[u8]) -> Option<Self> {
        // Gemmi❗✔️: struct Chain {
        // Gemmi❗✔️:   static const char* what() { return "Chain"; }
        // Gemmi❗✔️:   std::string name;
        //
        // Behavior review: Gemmi retains an unrestricted string. COSMolKit's
        // approved representational boundary is zero through four ASCII bytes;
        // the IO owner reports wider/non-ASCII source identifiers structurally.
        // Complexity review: validation and the bounded copy are O(n), n <= 4,
        // with no heap allocation.
        if bytes.len() > 4 || !bytes.is_ascii() {
            return None;
        }
        let mut stored = [0; 4];
        stored[..bytes.len()].copy_from_slice(bytes);
        Some(Self {
            bytes: stored,
            len: bytes.len() as u8,
        })
    }

    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.bytes.split_at(self.len as usize).0
    }

    #[must_use]
    pub fn as_str(&self) -> &str {
        // Construction accepts ASCII only, so this conversion is infallible.
        std::str::from_utf8(self.as_bytes()).expect("PdbChainId invariant")
    }
}

/// PDB author residue sequence number and optional insertion code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct PdbSeqId {
    seq_num: i32,
    ins_code: Option<u8>,
}

impl PdbSeqId {
    #[must_use]
    pub const fn new(seq_num: i32, ins_code: Option<u8>) -> Self {
        // Gemmi✔️✔️: struct SeqId {
        // Gemmi✔️✔️:   using OptionalNum = OptionalInt<INT_MIN>;
        // Gemmi✔️✔️:
        // Gemmi✔️✔️:   OptionalNum num;   // sequence number
        // Gemmi✔️✔️:   char icode = ' ';  // insertion code
        Self { seq_num, ins_code }
    }

    #[must_use]
    pub const fn seq_num(self) -> i32 {
        self.seq_num
    }

    #[must_use]
    pub const fn ins_code(self) -> Option<u8> {
        self.ins_code
    }
}

/// Exact four-column PDB atom name, including leading and trailing spaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AtomName([u8; 4]);

impl AtomName {
    #[must_use]
    pub const fn from_ascii(bytes: [u8; 4]) -> Option<Self> {
        // Gemmi❗✔️: struct Atom {
        // Gemmi❗✔️:   static const char* what() { return "Atom"; }
        // Gemmi❗✔️:   std::string name;
        //
        // Behavior review: the four-byte PDB spelling is retained exactly.
        // Wider mmCIF names are rejected by the IO representational boundary,
        // rather than silently truncated here.
        if bytes.is_ascii() {
            Some(Self(bytes))
        } else {
            None
        }
    }

    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 4] {
        &self.0
    }

    #[must_use]
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.0).expect("AtomName invariant")
    }
}

/// Residue name within COSMolKit's declared four-byte source boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ResidueName {
    bytes: [u8; 4],
    len: u8,
}

impl ResidueName {
    #[must_use]
    pub fn from_ascii(bytes: &[u8]) -> Option<Self> {
        // Gemmi❗✔️: struct ResidueId {
        // Gemmi❗✔️:   SeqId seqid;
        // Gemmi❗✔️:   std::string segment; // segid - up to 4 characters in the PDB file
        // Gemmi❗✔️:   std::string name;
        //
        // Behavior review: exact bytes are preserved within the approved
        // zero-through-four-byte boundary. Gemmi itself stores a string.
        if bytes.len() > 4 || !bytes.is_ascii() {
            return None;
        }
        let mut stored = [0; 4];
        stored[..bytes.len()].copy_from_slice(bytes);
        Some(Self {
            bytes: stored,
            len: bytes.len() as u8,
        })
    }

    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.bytes.split_at(self.len as usize).0
    }

    #[must_use]
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(self.as_bytes()).expect("ResidueName invariant")
    }
}

/// Present alternate-location byte after IO has interpreted missing markers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AltLocLabel(u8);

impl AltLocLabel {
    #[must_use]
    pub const fn new(value: u8) -> Self {
        // Gemmi✔️✔️: char altloc = '\0'; // 0 if not set
        Self(value)
    }

    #[must_use]
    pub const fn value(self) -> u8 {
        self.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct AtomSourceIds {
    serial: Option<PdbAtomSerial>,
}

impl AtomSourceIds {
    #[must_use]
    pub const fn new(serial: Option<PdbAtomSerial>) -> Self {
        Self { serial }
    }

    #[must_use]
    pub const fn serial(self) -> Option<PdbAtomSerial> {
        self.serial
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct ResidueSourceIds {
    seq_id: Option<PdbSeqId>,
    label_seq_id: Option<i32>,
    segment_id: Option<[u8; 4]>,
    subchain_id: Option<String>,
    label_entity_id: Option<String>,
}

impl ResidueSourceIds {
    #[must_use]
    pub fn new(
        seq_id: Option<PdbSeqId>,
        label_seq_id: Option<i32>,
        segment_id: Option<[u8; 4]>,
        subchain_id: Option<String>,
        label_entity_id: Option<String>,
    ) -> Option<Self> {
        // Gemmi✔️✔️: std::string subchain;   // mmCIF _atom_site.label_asym_id
        // Gemmi✔️✔️: std::string entity_id;  // mmCIF _atom_site.label_entity_id
        // Gemmi✔️✔️: OptionalNum label_seq;  // mmCIF _atom_site.label_seq_id
        //
        // Source entity spelling is retained as source text, not confused
        // with the hierarchy unit's future stable EntityId row identity.
        if segment_id.is_some_and(|segment| !segment.is_ascii()) {
            return None;
        }
        Some(Self {
            seq_id,
            label_seq_id,
            segment_id,
            subchain_id,
            label_entity_id,
        })
    }

    #[must_use]
    pub const fn seq_id(&self) -> Option<PdbSeqId> {
        self.seq_id
    }

    #[must_use]
    pub const fn label_seq_id(&self) -> Option<i32> {
        self.label_seq_id
    }

    #[must_use]
    pub const fn segment_id(&self) -> Option<&[u8; 4]> {
        self.segment_id.as_ref()
    }

    #[must_use]
    pub fn subchain_id(&self) -> Option<&str> {
        self.subchain_id.as_deref()
    }

    #[must_use]
    pub fn label_entity_id(&self) -> Option<&str> {
        self.label_entity_id.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct ChainSourceIds {
    auth_chain_id: Option<PdbChainId>,
    label_asym_id: Option<String>,
}

impl ChainSourceIds {
    #[must_use]
    pub fn new(auth_chain_id: Option<PdbChainId>, label_asym_id: Option<String>) -> Self {
        Self {
            auth_chain_id,
            label_asym_id,
        }
    }

    #[must_use]
    pub const fn auth_chain_id(&self) -> Option<PdbChainId> {
        self.auth_chain_id
    }

    #[must_use]
    pub fn label_asym_id(&self) -> Option<&str> {
        self.label_asym_id.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct EntitySourceIds {
    source_entity_id: String,
}

impl EntitySourceIds {
    #[must_use]
    pub fn new(source_entity_id: String) -> Self {
        Self { source_entity_id }
    }

    #[must_use]
    pub fn source_entity_id(&self) -> &str {
        &self.source_entity_id
    }
}
