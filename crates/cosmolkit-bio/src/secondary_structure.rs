//! Source-addressed secondary-structure values retained from Gemmi.

use crate::AtomAddress;

/// PDB helix class values used by Gemmi's `Helix` metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
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

/// Helix endpoints and reported classification from PDB/mmCIF metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BioHelix {
    pub start: AtomAddress,
    pub end: AtomAddress,
    pub pdb_helix_class: BioHelixClass,
    pub length: i32,
}

impl BioHelix {
    /// Apply Gemmi's integer helix-class setter; values outside 1 through 10
    /// leave the previous class unchanged.
    pub fn set_helix_class_as_int(&mut self, n: i32) {
        // Gemmi✔️✔️: enum HelixClass {
        // Gemmi✔️✔️:   UnknownHelix, RAlpha, ROmega, RPi, RGamma, R310,
        // Gemmi✔️✔️:   LAlpha, LOmega, LGamma, Helix27, HelixPolyProlineNone
        // Gemmi✔️✔️: };
        // Gemmi✔️✔️: void set_helix_class_as_int(int n) {
        // Gemmi✔️✔️:   if (n >= 1 && n <= 10)
        // Gemmi✔️✔️:     pdb_helix_class = static_cast<HelixClass>(n);
        // Gemmi✔️✔️: }
        // Behavior review: matching the source's inclusive guard preserves
        // prior state for every other i32; the enum discriminants mirror the
        // source declaration order used by its cast.
        // Complexity review: one scalar range check and a bounded branch,
        // with no allocation, scan, or data-dependent lookup.
        self.pdb_helix_class = match n {
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
            _ => return,
        };
    }
}

impl Default for BioHelix {
    fn default() -> Self {
        // Gemmi✔️✔️: AtomAddress start, end;
        // Gemmi✔️✔️: HelixClass pdb_helix_class = UnknownHelix;
        // Gemmi✔️✔️: int length = -1;
        // Behavior review: default atom addresses carry Gemmi's unresolved
        // empty/NUL state; class and length match the source member initializers.
        // Complexity review: fixed-field construction; default empty strings
        // require no heap allocation.
        Self {
            start: AtomAddress::default(),
            end: AtomAddress::default(),
            pdb_helix_class: BioHelixClass::UnknownHelix,
            length: -1,
        }
    }
}

/// One ordered strand and its source endpoint/hydrogen-bond addresses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BioStrand {
    pub start: AtomAddress,
    pub end: AtomAddress,
    pub hbond_atom2: AtomAddress,
    pub hbond_atom1: AtomAddress,
    pub sense: i32,
    pub name: String,
}

impl BioStrand {
    /// Construct a strand with an explicit, source-defined sense value.
    /// Gemmi leaves this field uninitialized in a default-constructed strand,
    /// so COSMolKit intentionally has no implicit `Default` for this value.
    #[must_use]
    pub fn new(
        start: AtomAddress,
        end: AtomAddress,
        hbond_atom2: AtomAddress,
        hbond_atom1: AtomAddress,
        sense: i32,
        name: String,
    ) -> Self {
        // Gemmi❗🔝: struct Strand {
        // Gemmi❗🔝:   AtomAddress start, end;
        // Gemmi❗🔝:   AtomAddress hbond_atom2, hbond_atom1;
        // Gemmi❗🔝:   int sense;  // 0 = first strand, 1 = parallel, -1 = anti-parallel.
        // Gemmi❗🔝:   std::string name; // optional, _struct_sheet_range.id if from mmCIF
        // Gemmi❗🔝: };
        // Behavior review: all four source addresses and the full signed sense
        // are retained; the absent C++ initializer for `sense` is not replaced
        // by a guessed Rust default. Callers must supply a defined value.
        // Complexity review: owned address/string fields move into the value,
        // avoiding the C++ parser's subsequent per-string assignments/copies.
        Self {
            start,
            end,
            hbond_atom2,
            hbond_atom1,
            sense,
            name,
        }
    }
}

/// A sheet containing strands in source order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BioSheet {
    pub name: String,
    pub strands: Vec<BioStrand>,
}

impl BioSheet {
    /// Create an empty sheet with its source identifier.
    #[must_use]
    pub fn new(sheet_id: &str) -> Self {
        // Gemmi✔️✔️: explicit Sheet(const std::string& sheet_id) noexcept : name(sheet_id) {}
        // Behavior review: the supplied sheet id is copied exactly and the
        // strand sequence begins empty, matching the source constructor.
        // Complexity review: copying the identifier is linear in its byte
        // length, as with the source's std::string copy construction.
        Self {
            name: String::from(sheet_id),
            strands: Vec::new(),
        }
    }
}

impl Default for BioSheet {
    fn default() -> Self {
        // Gemmi✔️✔️: std::string name;
        // Gemmi✔️✔️: std::vector<Strand> strands;
        // Gemmi✔️✔️: Sheet() = default;
        // Behavior review: both owned collections start empty, and no strand
        // (whose source `sense` has no initializer) is fabricated.
        // Complexity review: empty Rust String/Vec construction performs no
        // heap allocation.
        Self {
            name: String::new(),
            strands: Vec::new(),
        }
    }
}
