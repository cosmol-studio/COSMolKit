//! Pinned MACCS keys over detached inputs and the single SEARCH matcher.
use crate::{Fingerprint, FingerprintError};
use cosmolkit_core::{PathError, RingFindingError, RingInfo, ValenceAssignment, ValenceError};
use cosmolkit_model::{CoordinateBlock, TopologyBlock, TopologyValidationError};
use cosmolkit_search::{
    CompiledQuery, QueryCompileError, QueryMatchContext, QueryMatchContextError, SearchTarget,
    SmartsParseError, SubstructMatchError, SubstructMatchParams,
    build_prepared_query_match_context, parse_smarts,
    try_get_substruct_atom_matches_with_compiled_query_and_context,
};
use std::{error::Error, fmt, sync::OnceLock};

/// Exact source fixed-width public projection. Other widths remain source-defined rejected options.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaccsFingerprintParams {
    pub n_bits: usize,
}
impl Default for MaccsFingerprintParams {
    fn default() -> Self {
        // COSMolKit✔️✔️: impl Default for MaccsFingerprintParams {
        // COSMolKit✔️✔️:     fn default() -> Self {
        // COSMolKit✔️✔️:         Self {
        // COSMolKit✔️✔️:             n_bits: COSMOLKIT_MACCS_PUBLIC_BITS,
        // COSMolKit✔️✔️:         }
        // COSMolKit✔️✔️:     }
        // COSMolKit✔️✔️: }
        Self {
            n_bits: COSMOLKIT_MACCS_PUBLIC_BITS,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum MaccsFingerprintError {
    UnsupportedOption {
        option: &'static str,
        reason: &'static str,
    },
    MissingPattern {
        bit: usize,
    },
    Topology(TopologyValidationError),
    Rings(RingFindingError),
    Valence(ValenceError),
    Paths(PathError),
    Smarts(SmartsParseError),
    QueryCompile(QueryCompileError),
    QueryContext(QueryMatchContextError),
    Match(SubstructMatchError),
    Value(FingerprintError),
}
impl fmt::Display for MaccsFingerprintError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedOption { option, reason } => {
                write!(f, "unsupported fingerprint option {option}: {reason}")
            }
            Self::MissingPattern { bit } => write!(f, "missing static MACCS pattern {bit}"),
            Self::Topology(e) => e.fmt(f),
            Self::Rings(e) => e.fmt(f),
            Self::Valence(e) => e.fmt(f),
            Self::Paths(e) => e.fmt(f),
            Self::Smarts(e) => e.fmt(f),
            Self::QueryCompile(e) => e.fmt(f),
            Self::QueryContext(e) => e.fmt(f),
            Self::Match(e) => e.fmt(f),
            Self::Value(e) => e.fmt(f),
        }
    }
}
impl Error for MaccsFingerprintError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::UnsupportedOption { .. } | Self::MissingPattern { .. } => None,
            Self::Topology(e) => Some(e),
            Self::Rings(e) => Some(e),
            Self::Valence(e) => Some(e),
            Self::Paths(e) => Some(e),
            Self::Smarts(e) => Some(e),
            Self::QueryCompile(e) => Some(e),
            Self::QueryContext(e) => Some(e),
            Self::Match(e) => Some(e),
            Self::Value(e) => Some(e),
        }
    }
}
macro_rules! maccs_error_from {
    ($source:ty,$variant:ident) => {
        impl From<$source> for MaccsFingerprintError {
            fn from(e: $source) -> Self {
                Self::$variant(e)
            }
        }
    };
}
maccs_error_from!(TopologyValidationError, Topology);
maccs_error_from!(RingFindingError, Rings);
maccs_error_from!(ValenceError, Valence);
maccs_error_from!(PathError, Paths);
maccs_error_from!(SmartsParseError, Smarts);
maccs_error_from!(QueryCompileError, QueryCompile);
maccs_error_from!(QueryMatchContextError, QueryContext);
maccs_error_from!(SubstructMatchError, Match);
maccs_error_from!(FingerprintError, Value);

const RDKIT_MACCS_RAW_BITS: usize = 167;
const COSMOLKIT_MACCS_PUBLIC_BITS: usize = 166;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MaccsPatternSpec {
    bit: usize,
    smarts: &'static str,
}
pub(crate) const RDKIT_MACCS_PATTERNS: &[MaccsPatternSpec] = &[
    MaccsPatternSpec {
        bit: 8,
        smarts: "[!#6!#1]1~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 11,
        smarts: "*1~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 13,
        smarts: "[#8]~[#7](~[#6])~[#6]",
    },
    MaccsPatternSpec {
        bit: 14,
        smarts: "[#16]-[#16]",
    },
    MaccsPatternSpec {
        bit: 15,
        smarts: "[#8]~[#6](~[#8])~[#8]",
    },
    MaccsPatternSpec {
        bit: 16,
        smarts: "[!#6!#1]1~*~*~1",
    },
    MaccsPatternSpec {
        bit: 17,
        smarts: "[#6]#[#6]",
    },
    MaccsPatternSpec {
        bit: 19,
        smarts: "*1~*~*~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 20,
        smarts: "[#14]",
    },
    MaccsPatternSpec {
        bit: 21,
        smarts: "[#6]=[#6](~[!#6!#1])~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 22,
        smarts: "*1~*~*~1",
    },
    MaccsPatternSpec {
        bit: 23,
        smarts: "[#7]~[#6](~[#8])~[#8]",
    },
    MaccsPatternSpec {
        bit: 24,
        smarts: "[#7]-[#8]",
    },
    MaccsPatternSpec {
        bit: 25,
        smarts: "[#7]~[#6](~[#7])~[#7]",
    },
    MaccsPatternSpec {
        bit: 26,
        smarts: "[#6]=@[#6](@*)@*",
    },
    MaccsPatternSpec {
        bit: 28,
        smarts: "[!#6!#1]~[CH2]~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 30,
        smarts: "[#6]~[!#6!#1](~[#6])(~[#6])~*",
    },
    MaccsPatternSpec {
        bit: 31,
        smarts: "[!#6!#1]~[F,Cl,Br,I]",
    },
    MaccsPatternSpec {
        bit: 32,
        smarts: "[#6]~[#16]~[#7]",
    },
    MaccsPatternSpec {
        bit: 33,
        smarts: "[#7]~[#16]",
    },
    MaccsPatternSpec {
        bit: 34,
        smarts: "[CH2]=*",
    },
    MaccsPatternSpec {
        bit: 36,
        smarts: "[#16R]",
    },
    MaccsPatternSpec {
        bit: 37,
        smarts: "[#7]~[#6](~[#8])~[#7]",
    },
    MaccsPatternSpec {
        bit: 38,
        smarts: "[#7]~[#6](~[#6])~[#7]",
    },
    MaccsPatternSpec {
        bit: 39,
        smarts: "[#8]~[#16](~[#8])~[#8]",
    },
    MaccsPatternSpec {
        bit: 40,
        smarts: "[#16]-[#8]",
    },
    MaccsPatternSpec {
        bit: 41,
        smarts: "[#6]#[#7]",
    },
    MaccsPatternSpec {
        bit: 43,
        smarts: "[!#6!#1!H0]~*~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 44,
        smarts: "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]",
    },
    MaccsPatternSpec {
        bit: 45,
        smarts: "[#6]=[#6]~[#7]",
    },
    MaccsPatternSpec {
        bit: 47,
        smarts: "[#16]~*~[#7]",
    },
    MaccsPatternSpec {
        bit: 48,
        smarts: "[#8]~[!#6!#1](~[#8])~[#8]",
    },
    MaccsPatternSpec {
        bit: 49,
        smarts: "[!+0]",
    },
    MaccsPatternSpec {
        bit: 50,
        smarts: "[#6]=[#6](~[#6])~[#6]",
    },
    MaccsPatternSpec {
        bit: 51,
        smarts: "[#6]~[#16]~[#8]",
    },
    MaccsPatternSpec {
        bit: 52,
        smarts: "[#7]~[#7]",
    },
    MaccsPatternSpec {
        bit: 53,
        smarts: "[!#6!#1!H0]~*~*~*~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 54,
        smarts: "[!#6!#1!H0]~*~*~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 55,
        smarts: "[#8]~[#16]~[#8]",
    },
    MaccsPatternSpec {
        bit: 56,
        smarts: "[#8]~[#7](~[#8])~[#6]",
    },
    MaccsPatternSpec {
        bit: 57,
        smarts: "[#8R]",
    },
    MaccsPatternSpec {
        bit: 58,
        smarts: "[!#6!#1]~[#16]~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 59,
        smarts: "[#16]!:*:*",
    },
    MaccsPatternSpec {
        bit: 60,
        smarts: "[#16]=[#8]",
    },
    MaccsPatternSpec {
        bit: 61,
        smarts: "*~[#16](~*)~*",
    },
    MaccsPatternSpec {
        bit: 62,
        smarts: "*@*!@*@*",
    },
    MaccsPatternSpec {
        bit: 63,
        smarts: "[#7]=[#8]",
    },
    MaccsPatternSpec {
        bit: 64,
        smarts: "*@*!@[#16]",
    },
    MaccsPatternSpec {
        bit: 65,
        smarts: "c:n",
    },
    MaccsPatternSpec {
        bit: 66,
        smarts: "[#6]~[#6](~[#6])(~[#6])~*",
    },
    MaccsPatternSpec {
        bit: 67,
        smarts: "[!#6!#1]~[#16]",
    },
    MaccsPatternSpec {
        bit: 68,
        smarts: "[!#6!#1!H0]~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 69,
        smarts: "[!#6!#1]~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 70,
        smarts: "[!#6!#1]~[#7]~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 71,
        smarts: "[#7]~[#8]",
    },
    MaccsPatternSpec {
        bit: 72,
        smarts: "[#8]~*~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 73,
        smarts: "[#16]=*",
    },
    MaccsPatternSpec {
        bit: 74,
        smarts: "[CH3]~*~[CH3]",
    },
    MaccsPatternSpec {
        bit: 75,
        smarts: "*!@[#7]@*",
    },
    MaccsPatternSpec {
        bit: 76,
        smarts: "[#6]=[#6](~*)~*",
    },
    MaccsPatternSpec {
        bit: 77,
        smarts: "[#7]~*~[#7]",
    },
    MaccsPatternSpec {
        bit: 78,
        smarts: "[#6]=[#7]",
    },
    MaccsPatternSpec {
        bit: 79,
        smarts: "[#7]~*~*~[#7]",
    },
    MaccsPatternSpec {
        bit: 80,
        smarts: "[#7]~*~*~*~[#7]",
    },
    MaccsPatternSpec {
        bit: 81,
        smarts: "[#16]~*(~*)~*",
    },
    MaccsPatternSpec {
        bit: 82,
        smarts: "*~[CH2]~[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 83,
        smarts: "[!#6!#1]1~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 84,
        smarts: "[NH2]",
    },
    MaccsPatternSpec {
        bit: 85,
        smarts: "[#6]~[#7](~[#6])~[#6]",
    },
    MaccsPatternSpec {
        bit: 86,
        smarts: "[C;H2,H3][!#6!#1][C;H2,H3]",
    },
    MaccsPatternSpec {
        bit: 87,
        smarts: "[F,Cl,Br,I]!@*@*",
    },
    MaccsPatternSpec {
        bit: 89,
        smarts: "[#8]~*~*~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 90,
        smarts: "[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 91,
        smarts: "[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 92,
        smarts: "[#8]~[#6](~[#7])~[#6]",
    },
    MaccsPatternSpec {
        bit: 93,
        smarts: "[!#6!#1]~[CH3]",
    },
    MaccsPatternSpec {
        bit: 94,
        smarts: "[!#6!#1]~[#7]",
    },
    MaccsPatternSpec {
        bit: 95,
        smarts: "[#7]~*~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 96,
        smarts: "*1~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 97,
        smarts: "[#7]~*~*~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 98,
        smarts: "[!#6!#1]1~*~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 99,
        smarts: "[#6]=[#6]",
    },
    MaccsPatternSpec {
        bit: 100,
        smarts: "*~[CH2]~[#7]",
    },
    MaccsPatternSpec {
        bit: 101,
        smarts: "[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]",
    },
    MaccsPatternSpec {
        bit: 102,
        smarts: "[!#6!#1]~[#8]",
    },
    MaccsPatternSpec {
        bit: 104,
        smarts: "[!#6!#1!H0]~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 105,
        smarts: "*@*(@*)@*",
    },
    MaccsPatternSpec {
        bit: 106,
        smarts: "[!#6!#1]~*(~[!#6!#1])~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 107,
        smarts: "[F,Cl,Br,I]~*(~*)~*",
    },
    MaccsPatternSpec {
        bit: 108,
        smarts: "[CH3]~*~*~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 109,
        smarts: "*~[CH2]~[#8]",
    },
    MaccsPatternSpec {
        bit: 110,
        smarts: "[#7]~[#6]~[#8]",
    },
    MaccsPatternSpec {
        bit: 111,
        smarts: "[#7]~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 112,
        smarts: "*~*(~*)(~*)~*",
    },
    MaccsPatternSpec {
        bit: 113,
        smarts: "[#8]!:*:*",
    },
    MaccsPatternSpec {
        bit: 114,
        smarts: "[CH3]~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 115,
        smarts: "[CH3]~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 116,
        smarts: "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]",
    },
    MaccsPatternSpec {
        bit: 117,
        smarts: "[#7]~*~[#8]",
    },
    MaccsPatternSpec {
        bit: 118,
        smarts: "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]",
    },
    MaccsPatternSpec {
        bit: 119,
        smarts: "[#7]=*",
    },
    MaccsPatternSpec {
        bit: 120,
        smarts: "[!#6R]",
    },
    MaccsPatternSpec {
        bit: 121,
        smarts: "[#7R]",
    },
    MaccsPatternSpec {
        bit: 122,
        smarts: "*~[#7](~*)~*",
    },
    MaccsPatternSpec {
        bit: 123,
        smarts: "[#8]~[#6]~[#8]",
    },
    MaccsPatternSpec {
        bit: 124,
        smarts: "[!#6!#1]~[!#6!#1]",
    },
    MaccsPatternSpec {
        bit: 126,
        smarts: "*!@[#8]!@*",
    },
    MaccsPatternSpec {
        bit: 127,
        smarts: "*@*!@[#8]",
    },
    MaccsPatternSpec {
        bit: 128,
        smarts: "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 129,
        smarts: "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 131,
        smarts: "[!#6!#1!H0]",
    },
    MaccsPatternSpec {
        bit: 132,
        smarts: "[#8]~*~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 133,
        smarts: "*@*!@[#7]",
    },
    MaccsPatternSpec {
        bit: 135,
        smarts: "[#7]!:*:*",
    },
    MaccsPatternSpec {
        bit: 136,
        smarts: "[#8]=*",
    },
    MaccsPatternSpec {
        bit: 137,
        smarts: "[!C!cR]",
    },
    MaccsPatternSpec {
        bit: 138,
        smarts: "[!#6!#1]~[CH2]~*",
    },
    MaccsPatternSpec {
        bit: 139,
        smarts: "[O!H0]",
    },
    MaccsPatternSpec {
        bit: 140,
        smarts: "[#8]",
    },
    MaccsPatternSpec {
        bit: 141,
        smarts: "[CH3]",
    },
    MaccsPatternSpec {
        bit: 142,
        smarts: "[#7]",
    },
    MaccsPatternSpec {
        bit: 144,
        smarts: "*!:*:*!:*",
    },
    MaccsPatternSpec {
        bit: 145,
        smarts: "*1~*~*~*~*~*~1",
    },
    MaccsPatternSpec {
        bit: 147,
        smarts: "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]",
    },
    MaccsPatternSpec {
        bit: 148,
        smarts: "*~[!#6!#1](~*)~*",
    },
    MaccsPatternSpec {
        bit: 149,
        smarts: "[C;H3,H4]",
    },
    MaccsPatternSpec {
        bit: 150,
        smarts: "*!@*@*!@*",
    },
    MaccsPatternSpec {
        bit: 151,
        smarts: "[#7!H0]",
    },
    MaccsPatternSpec {
        bit: 152,
        smarts: "[#8]~[#6](~[#6])~[#6]",
    },
    MaccsPatternSpec {
        bit: 154,
        smarts: "[#6]=[#8]",
    },
    MaccsPatternSpec {
        bit: 155,
        smarts: "*!@[CH2]!@*",
    },
    MaccsPatternSpec {
        bit: 156,
        smarts: "[#7]~*(~*)~*",
    },
    MaccsPatternSpec {
        bit: 157,
        smarts: "[#6]-[#8]",
    },
    MaccsPatternSpec {
        bit: 158,
        smarts: "[#6]-[#7]",
    },
    MaccsPatternSpec {
        bit: 162,
        smarts: "a",
    },
    MaccsPatternSpec {
        bit: 165,
        smarts: "[R]",
    },
];

fn rdkit_maccs_pattern(bit: usize) -> Option<&'static MaccsPatternSpec> {
    RDKIT_MACCS_PATTERNS.iter().find(|p| p.bit == bit)
}
fn rdkit_maccs_public_index(bit: usize) -> Option<usize> {
    (1..RDKIT_MACCS_RAW_BITS).contains(&bit).then(|| bit - 1)
}

fn cached_maccs_queries() -> Result<&'static [(usize, CompiledQuery)], MaccsFingerprintError> {
    // RDKit❗✔️: struct Patterns {
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_8 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]1~*~*~*~1"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_11 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*1~*~*~*~1"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_13 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#8]~[#7](~[#6])~[#6]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_14 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]-[#16]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_15 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#8]~[#6](~[#8])~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_16 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]1~*~*~1"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_17 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]#[#6]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_19 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*1~*~*~*~*~*~*~1"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_20 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#14]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_21 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#6]=[#6](~[!#6!#1])~[!#6!#1]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_22 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*1~*~*~1"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_23 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#7]~[#6](~[#8])~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_24 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]-[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_25 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#7]~[#6](~[#7])~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_26 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=@[#6](@*)@*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_28 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1]~[CH2]~[!#6!#1]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_30 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#6]~[!#6!#1](~[#6])(~[#6])~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_31 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]~[F,Cl,Br,I]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_32 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#16]~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_33 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[#16]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_34 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[CH2]=*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_36 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16R]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_37 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#7]~[#6](~[#8])~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_38 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#7]~[#6](~[#6])~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_39 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#8]~[#16](~[#8])~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_40 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]-[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_41 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]#[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_43 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1!H0]~*~[!#6!#1!H0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_44 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
    // RDKit❗✔️:           "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_45 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=[#6]~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_47 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]~*~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_48 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#8]~[!#6!#1](~[#8])~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_49 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!+0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_50 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#6]=[#6](~[#6])~[#6]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_51 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]~[#16]~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_52 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_53 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1!H0]~*~*~*~[!#6!#1!H0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_54 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1!H0]~*~*~[!#6!#1!H0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_55 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~[#16]~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_56 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#8]~[#7](~[#8])~[#6]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_57 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8R]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_58 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1]~[#16]~[!#6!#1]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_59 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]!:*:*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_60 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]=[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_61 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~[#16](~*)~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_62 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*@*!@*@*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_63 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]=[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_64 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*@*!@[#16]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_65 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("c:n"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_66 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#6]~[#6](~[#6])(~[#6])~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_67 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]~[#16]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_68 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1!H0]~[!#6!#1!H0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_69 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]~[!#6!#1!H0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_70 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1]~[#7]~[!#6!#1]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_71 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_72 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~*~*~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_73 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]=*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_74 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[CH3]~*~[CH3]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_75 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*!@[#7]@*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_76 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=[#6](~*)~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_77 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~*~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_78 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_79 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~*~*~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_80 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~*~*~*~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_81 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#16]~*(~*)~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_82 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~[CH2]~[!#6!#1!H0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_83 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]1~*~*~*~*~1"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_84 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[NH2]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_85 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#6]~[#7](~[#6])~[#6]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_86 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[C;H2,H3][!#6!#1][C;H2,H3]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_87 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[F,Cl,Br,I]!@*@*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_89 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~*~*~*~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_90 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@["
    // RDKit❗✔️:                          "CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_91 = std::unique_ptr<
    // RDKit❗✔️:       RDKit::ROMol>(RDKit::SmartsToMol(
    // RDKit❗✔️:       "[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#"
    // RDKit❗✔️:       "6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_92 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#8]~[#6](~[#7])~[#6]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_93 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]~[CH3]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_94 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_95 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~*~*~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_96 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*1~*~*~*~*~1"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_97 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~*~*~*~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_98 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1]1~*~*~*~*~*~1"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_99 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=[#6]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_100 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~[CH2]~[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_101 = std::unique_ptr<
    // RDKit❗✔️:       RDKit::ROMol>(RDKit::SmartsToMol(
    // RDKit❗✔️:       "[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@["
    // RDKit❗✔️:       "R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@["
    // RDKit❗✔️:       "R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]"
    // RDKit❗✔️:       "@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@"
    // RDKit❗✔️:       "1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_102 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_104 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1!H0]~*~[CH2]~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_105 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*@*(@*)@*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_106 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[!#6!#1]~*(~[!#6!#1])~[!#6!#1]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_107 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[F,Cl,Br,I]~*(~*)~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_108 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[CH3]~*~*~*~[CH2]~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_109 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~[CH2]~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_110 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~[#6]~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_111 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~*~[CH2]~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_112 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~*(~*)(~*)~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_113 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]!:*:*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_114 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[CH3]~[CH2]~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_115 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[CH3]~*~[CH2]~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_116 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_117 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~*~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_118 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_119 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]=*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_120 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6R]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_121 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7R]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_122 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~[#7](~*)~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_123 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~[#6]~[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_124 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]~[!#6!#1]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_126 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*!@[#8]!@*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_127 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*@*!@[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_128 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol(
    // RDKit❗✔️:           "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~["
    // RDKit❗✔️:           "CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_129 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1)"
    // RDKit❗✔️:                          ",$(*~[CH2]~[R]1@[R]@[CH2R]1)]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_131 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1!H0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_132 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]~*~[CH2]~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_133 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*@*!@[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_135 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]!:*:*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_136 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]=*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_137 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!C!cR]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_138 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[!#6!#1]~[CH2]~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_139 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[O!H0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_140 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_141 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[CH3]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_142 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_144 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*!:*:*!:*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_145 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*1~*~*~*~*~*~1"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_147 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_148 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*~[!#6!#1](~*)~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_149 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[C;H3,H4]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_150 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*!@*@*!@*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_151 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7!H0]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_152 = std::unique_ptr<RDKit::ROMol>(
    // RDKit❗✔️:       RDKit::SmartsToMol("[#8]~[#6](~[#6])~[#6]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_154 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]=[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_155 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("*!@[CH2]!@*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_156 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#7]~*(~*)~*"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_157 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]-[#8]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_158 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[#6]-[#7]"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_162 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("a"));
    // RDKit❗✔️:   std::unique_ptr<RDKit::ROMol> bit_165 =
    // RDKit❗✔️:       std::unique_ptr<RDKit::ROMol>(RDKit::SmartsToMol("[R]"));
    // RDKit❗✔️: };
    // Complexity: source caches the 136 parsed patterns; Rust caches those same patterns
    // plus query-side adjacency/order once. No target data or match state enters this cache.
    static CACHE: OnceLock<Result<Vec<(usize, CompiledQuery)>, MaccsFingerprintError>> =
        OnceLock::new();
    match CACHE.get_or_init(|| {
        RDKIT_MACCS_PATTERNS
            .iter()
            .map(|p| {
                let graph = parse_smarts(p.smarts, &Default::default())?;
                Ok((p.bit, CompiledQuery::compile(graph)?))
            })
            .collect()
    }) {
        Ok(v) => Ok(v.as_slice()),
        Err(e) => Err(e.clone()),
    }
}

fn maccs_match_count(
    bit: usize,
    queries: &[(usize, CompiledQuery)],
    target: &SearchTarget<'_>,
    context: &QueryMatchContext<'_>,
    max_matches: usize,
) -> Result<usize, MaccsFingerprintError> {
    // RDKit❗❌:   count = RDKit::SubstructMatch(mol, *pats.bit_145, matches, true, true);
    // Cost: each call scans up to 136 cached queries instead of the source's O(1)
    // direct pattern-member access. SEARCH matching still borrows the cached query/context.
    // Exact source maxMatches 1000 for counting and 1 for first-match queries; uniquify true.
    // Source MATCHER remains SEARCH-owned; cached query/context borrowed, no chemistry copy.
    let query = &queries
        .iter()
        .find(|(k, _)| *k == bit)
        .ok_or(MaccsFingerprintError::MissingPattern { bit })?
        .1;
    let params = SubstructMatchParams {
        max_matches,
        ..Default::default()
    };
    Ok(
        try_get_substruct_atom_matches_with_compiled_query_and_context(
            target, query, &params, context,
        )?
        .len(),
    )
}
/// Raw source 167 bits, bit zero unused. Private runtime caches are never mutated.
#[doc(hidden)]
pub fn maccs_fingerprint_raw(
    topology: &TopologyBlock,
    rings: Option<&RingInfo>,
    valence: Option<&ValenceAssignment>,
) -> Result<Fingerprint, MaccsFingerprintError> {
    // RDKit❗❌: void GenerateFP(const RDKit::ROMol &mol, ExplicitBitVect &fp) {
    // RDKit❗❌:   if (!gpats.get()) {
    // RDKit❗❌:     gpats = std::unique_ptr<Patterns>(new Patterns());
    // RDKit❗❌:   }
    // RDKit❗❌:   const Patterns &pats = *(gpats.get());
    // RDKit❗❌:   PRECONDITION(fp.size() == 167, "bad fingerprint");
    // RDKit❗❌:   fp.clearBits();
    // RDKit❗❌:
    // RDKit❗❌:   if (!mol.getNumAtoms()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<RDKit::MatchVectType> matches;
    // RDKit❗❌:   RDKit::RWMol::ConstAtomIterator atom;
    // RDKit❗❌:   RDKit::MatchVectType match;
    // RDKit❗❌:   unsigned int count;
    // RDKit❗❌:
    // RDKit❗❌:   for (atom = mol.beginAtoms(); atom != mol.endAtoms(); ++atom) {
    // RDKit❗❌:     switch ((*atom)->getAtomicNum()) {
    // RDKit❗❌:       case 3:
    // RDKit❗❌:       case 11:
    // RDKit❗❌:       case 19:
    // RDKit❗❌:       case 37:
    // RDKit❗❌:       case 55:
    // RDKit❗❌:       case 87:
    // RDKit❗❌:         fp.setBit(35);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 4:
    // RDKit❗❌:       case 12:
    // RDKit❗❌:       case 20:
    // RDKit❗❌:       case 38:
    // RDKit❗❌:       case 56:
    // RDKit❗❌:       case 88:
    // RDKit❗❌:         fp.setBit(10);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 5:
    // RDKit❗❌:       case 13:
    // RDKit❗❌:       case 31:
    // RDKit❗❌:       case 49:
    // RDKit❗❌:       case 81:
    // RDKit❗❌:         fp.setBit(18);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 9:
    // RDKit❗❌:         fp.setBit(42);
    // RDKit❗❌:         fp.setBit(134);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 15:
    // RDKit❗❌:         fp.setBit(29);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 16:
    // RDKit❗❌:         fp.setBit(88);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 17:
    // RDKit❗❌:         fp.setBit(103);
    // RDKit❗❌:         fp.setBit(134);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 21:
    // RDKit❗❌:       case 22:
    // RDKit❗❌:       case 39:
    // RDKit❗❌:       case 40:
    // RDKit❗❌:       case 72:
    // RDKit❗❌:         fp.setBit(5);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 23:
    // RDKit❗❌:       case 24:
    // RDKit❗❌:       case 25:
    // RDKit❗❌:       case 41:
    // RDKit❗❌:       case 42:
    // RDKit❗❌:       case 43:
    // RDKit❗❌:       case 73:
    // RDKit❗❌:       case 74:
    // RDKit❗❌:       case 75:
    // RDKit❗❌:         fp.setBit(7);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 26:
    // RDKit❗❌:       case 27:
    // RDKit❗❌:       case 28:
    // RDKit❗❌:       case 44:
    // RDKit❗❌:       case 45:
    // RDKit❗❌:       case 46:
    // RDKit❗❌:       case 76:
    // RDKit❗❌:       case 77:
    // RDKit❗❌:       case 78:
    // RDKit❗❌:         fp.setBit(9);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 29:
    // RDKit❗❌:       case 30:
    // RDKit❗❌:       case 47:
    // RDKit❗❌:       case 48:
    // RDKit❗❌:       case 79:
    // RDKit❗❌:       case 80:
    // RDKit❗❌:         fp.setBit(12);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 32:
    // RDKit❗❌:       case 33:
    // RDKit❗❌:       case 34:
    // RDKit❗❌:       case 50:
    // RDKit❗❌:       case 51:
    // RDKit❗❌:       case 52:
    // RDKit❗❌:       case 82:
    // RDKit❗❌:       case 83:
    // RDKit❗❌:       case 84:
    // RDKit❗❌:         fp.setBit(3);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 35:
    // RDKit❗❌:         fp.setBit(46);
    // RDKit❗❌:         fp.setBit(134);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 53:
    // RDKit❗❌:         fp.setBit(27);
    // RDKit❗❌:         fp.setBit(134);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 57:
    // RDKit❗❌:       case 58:
    // RDKit❗❌:       case 59:
    // RDKit❗❌:       case 60:
    // RDKit❗❌:       case 61:
    // RDKit❗❌:       case 62:
    // RDKit❗❌:       case 63:
    // RDKit❗❌:       case 64:
    // RDKit❗❌:       case 65:
    // RDKit❗❌:       case 66:
    // RDKit❗❌:       case 67:
    // RDKit❗❌:       case 68:
    // RDKit❗❌:       case 69:
    // RDKit❗❌:       case 70:
    // RDKit❗❌:       case 71:
    // RDKit❗❌:         fp.setBit(6);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 89:
    // RDKit❗❌:       case 90:
    // RDKit❗❌:       case 91:
    // RDKit❗❌:       case 92:
    // RDKit❗❌:       case 93:
    // RDKit❗❌:       case 94:
    // RDKit❗❌:       case 95:
    // RDKit❗❌:       case 96:
    // RDKit❗❌:       case 97:
    // RDKit❗❌:       case 98:
    // RDKit❗❌:       case 99:
    // RDKit❗❌:       case 100:
    // RDKit❗❌:       case 101:
    // RDKit❗❌:       case 102:
    // RDKit❗❌:       case 103:
    // RDKit❗❌:         fp.setBit(4);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case 104:
    // RDKit❗❌:         fp.setBit(2);
    // RDKit❗❌:         break;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_8, match, true)) {
    // RDKit❗❌:     fp.setBit(8);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_11, match, true)) {
    // RDKit❗❌:     fp.setBit(11);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_13, match, true)) {
    // RDKit❗❌:     fp.setBit(13);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_14, match, true)) {
    // RDKit❗❌:     fp.setBit(14);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_15, match, true)) {
    // RDKit❗❌:     fp.setBit(15);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_16, match, true)) {
    // RDKit❗❌:     fp.setBit(16);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_17, match, true)) {
    // RDKit❗❌:     fp.setBit(17);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_19, match, true)) {
    // RDKit❗❌:     fp.setBit(19);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_20, match, true)) {
    // RDKit❗❌:     fp.setBit(20);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_21, match, true)) {
    // RDKit❗❌:     fp.setBit(21);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_22, match, true)) {
    // RDKit❗❌:     fp.setBit(22);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_23, match, true)) {
    // RDKit❗❌:     fp.setBit(23);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_24, match, true)) {
    // RDKit❗❌:     fp.setBit(24);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_25, match, true)) {
    // RDKit❗❌:     fp.setBit(25);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_26, match, true)) {
    // RDKit❗❌:     fp.setBit(26);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_28, match, true)) {
    // RDKit❗❌:     fp.setBit(28);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_30, match, true)) {
    // RDKit❗❌:     fp.setBit(30);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_31, match, true)) {
    // RDKit❗❌:     fp.setBit(31);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_32, match, true)) {
    // RDKit❗❌:     fp.setBit(32);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_33, match, true)) {
    // RDKit❗❌:     fp.setBit(33);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_34, match, true)) {
    // RDKit❗❌:     fp.setBit(34);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_36, match, true)) {
    // RDKit❗❌:     fp.setBit(36);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_37, match, true)) {
    // RDKit❗❌:     fp.setBit(37);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_38, match, true)) {
    // RDKit❗❌:     fp.setBit(38);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_39, match, true)) {
    // RDKit❗❌:     fp.setBit(39);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_40, match, true)) {
    // RDKit❗❌:     fp.setBit(40);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_41, match, true)) {
    // RDKit❗❌:     fp.setBit(41);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_43, match, true)) {
    // RDKit❗❌:     fp.setBit(43);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_44, match, true)) {
    // RDKit❗❌:     fp.setBit(44);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_45, match, true)) {
    // RDKit❗❌:     fp.setBit(45);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_47, match, true)) {
    // RDKit❗❌:     fp.setBit(47);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_48, match, true)) {
    // RDKit❗❌:     fp.setBit(48);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_49, match, true)) {
    // RDKit❗❌:     fp.setBit(49);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_50, match, true)) {
    // RDKit❗❌:     fp.setBit(50);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_51, match, true)) {
    // RDKit❗❌:     fp.setBit(51);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_52, match, true)) {
    // RDKit❗❌:     fp.setBit(52);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_53, match, true)) {
    // RDKit❗❌:     fp.setBit(53);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_54, match, true)) {
    // RDKit❗❌:     fp.setBit(54);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_55, match, true)) {
    // RDKit❗❌:     fp.setBit(55);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_56, match, true)) {
    // RDKit❗❌:     fp.setBit(56);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_57, match, true)) {
    // RDKit❗❌:     fp.setBit(57);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_58, match, true)) {
    // RDKit❗❌:     fp.setBit(58);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_59, match, true)) {
    // RDKit❗❌:     fp.setBit(59);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_60, match, true)) {
    // RDKit❗❌:     fp.setBit(60);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_61, match, true)) {
    // RDKit❗❌:     fp.setBit(61);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_62, match, true)) {
    // RDKit❗❌:     fp.setBit(62);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_63, match, true)) {
    // RDKit❗❌:     fp.setBit(63);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_64, match, true)) {
    // RDKit❗❌:     fp.setBit(64);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_65, match, true)) {
    // RDKit❗❌:     fp.setBit(65);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_66, match, true)) {
    // RDKit❗❌:     fp.setBit(66);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_67, match, true)) {
    // RDKit❗❌:     fp.setBit(67);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_68, match, true)) {
    // RDKit❗❌:     fp.setBit(68);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_69, match, true)) {
    // RDKit❗❌:     fp.setBit(69);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_70, match, true)) {
    // RDKit❗❌:     fp.setBit(70);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_71, match, true)) {
    // RDKit❗❌:     fp.setBit(71);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_72, match, true)) {
    // RDKit❗❌:     fp.setBit(72);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_73, match, true)) {
    // RDKit❗❌:     fp.setBit(73);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_74, match, true)) {
    // RDKit❗❌:     fp.setBit(74);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_75, match, true)) {
    // RDKit❗❌:     fp.setBit(75);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_76, match, true)) {
    // RDKit❗❌:     fp.setBit(76);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_77, match, true)) {
    // RDKit❗❌:     fp.setBit(77);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_78, match, true)) {
    // RDKit❗❌:     fp.setBit(78);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_79, match, true)) {
    // RDKit❗❌:     fp.setBit(79);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_80, match, true)) {
    // RDKit❗❌:     fp.setBit(80);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_81, match, true)) {
    // RDKit❗❌:     fp.setBit(81);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_82, match, true)) {
    // RDKit❗❌:     fp.setBit(82);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_83, match, true)) {
    // RDKit❗❌:     fp.setBit(83);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_84, match, true)) {
    // RDKit❗❌:     fp.setBit(84);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_85, match, true)) {
    // RDKit❗❌:     fp.setBit(85);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_86, match, true)) {
    // RDKit❗❌:     fp.setBit(86);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_87, match, true)) {
    // RDKit❗❌:     fp.setBit(87);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_89, match, true)) {
    // RDKit❗❌:     fp.setBit(89);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_90, match, true)) {
    // RDKit❗❌:     fp.setBit(90);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_91, match, true)) {
    // RDKit❗❌:     fp.setBit(91);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_92, match, true)) {
    // RDKit❗❌:     fp.setBit(92);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_93, match, true)) {
    // RDKit❗❌:     fp.setBit(93);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_94, match, true)) {
    // RDKit❗❌:     fp.setBit(94);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_95, match, true)) {
    // RDKit❗❌:     fp.setBit(95);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_96, match, true)) {
    // RDKit❗❌:     fp.setBit(96);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_97, match, true)) {
    // RDKit❗❌:     fp.setBit(97);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_98, match, true)) {
    // RDKit❗❌:     fp.setBit(98);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_99, match, true)) {
    // RDKit❗❌:     fp.setBit(99);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_100, match, true)) {
    // RDKit❗❌:     fp.setBit(100);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_101, match, true)) {
    // RDKit❗❌:     fp.setBit(101);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_102, match, true)) {
    // RDKit❗❌:     fp.setBit(102);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_104, match, true)) {
    // RDKit❗❌:     fp.setBit(104);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_105, match, true)) {
    // RDKit❗❌:     fp.setBit(105);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_106, match, true)) {
    // RDKit❗❌:     fp.setBit(106);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_107, match, true)) {
    // RDKit❗❌:     fp.setBit(107);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_108, match, true)) {
    // RDKit❗❌:     fp.setBit(108);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_109, match, true)) {
    // RDKit❗❌:     fp.setBit(109);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_110, match, true)) {
    // RDKit❗❌:     fp.setBit(110);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_111, match, true)) {
    // RDKit❗❌:     fp.setBit(111);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_112, match, true)) {
    // RDKit❗❌:     fp.setBit(112);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_113, match, true)) {
    // RDKit❗❌:     fp.setBit(113);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_114, match, true)) {
    // RDKit❗❌:     fp.setBit(114);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_115, match, true)) {
    // RDKit❗❌:     fp.setBit(115);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_116, match, true)) {
    // RDKit❗❌:     fp.setBit(116);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_117, match, true)) {
    // RDKit❗❌:     fp.setBit(117);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_118, matches, true, true) > 1) {
    // RDKit❗❌:     fp.setBit(118);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_119, match, true)) {
    // RDKit❗❌:     fp.setBit(119);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_120, matches, true, true) > 1) {
    // RDKit❗❌:     fp.setBit(120);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_121, match, true)) {
    // RDKit❗❌:     fp.setBit(121);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_122, match, true)) {
    // RDKit❗❌:     fp.setBit(122);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_123, match, true)) {
    // RDKit❗❌:     fp.setBit(123);
    // RDKit❗❌:   }
    // RDKit❗❌:   count = RDKit::SubstructMatch(mol, *pats.bit_124, matches, true, true);
    // RDKit❗❌:   if (count > 0) {
    // RDKit❗❌:     fp.setBit(124);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (count > 1) {
    // RDKit❗❌:     fp.setBit(130);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_126, match, true)) {
    // RDKit❗❌:     fp.setBit(126);
    // RDKit❗❌:   }
    // RDKit❗❌:   count = RDKit::SubstructMatch(mol, *pats.bit_127, matches, true, true);
    // RDKit❗❌:   if (count > 1) {
    // RDKit❗❌:     fp.setBit(127);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (count > 0) {
    // RDKit❗❌:     fp.setBit(143);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_128, match, true)) {
    // RDKit❗❌:     fp.setBit(128);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_129, match, true)) {
    // RDKit❗❌:     fp.setBit(129);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_131, matches, true, true) > 1) {
    // RDKit❗❌:     fp.setBit(131);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_132, match, true)) {
    // RDKit❗❌:     fp.setBit(132);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_133, match, true)) {
    // RDKit❗❌:     fp.setBit(133);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_135, match, true)) {
    // RDKit❗❌:     fp.setBit(135);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_136, matches, true, true) > 1) {
    // RDKit❗❌:     fp.setBit(136);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_137, match, true)) {
    // RDKit❗❌:     fp.setBit(137);
    // RDKit❗❌:   }
    // RDKit❗❌:   count = RDKit::SubstructMatch(mol, *pats.bit_138, matches, true, true);
    // RDKit❗❌:   if (count > 1) {
    // RDKit❗❌:     fp.setBit(138);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (count > 0) {
    // RDKit❗❌:     fp.setBit(153);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_139, match, true)) {
    // RDKit❗❌:     fp.setBit(139);
    // RDKit❗❌:   }
    // RDKit❗❌:   count = RDKit::SubstructMatch(mol, *pats.bit_140, matches, true, true);
    // RDKit❗❌:   if (count > 3) {
    // RDKit❗❌:     fp.setBit(140);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (count > 2) {
    // RDKit❗❌:     fp.setBit(146);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (count > 1) {
    // RDKit❗❌:     fp.setBit(159);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (count > 0) {
    // RDKit❗❌:     fp.setBit(164);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_141, matches, true, true) > 2) {
    // RDKit❗❌:     fp.setBit(141);
    // RDKit❗❌:   }
    // RDKit❗❌:   count = RDKit::SubstructMatch(mol, *pats.bit_142, matches, true, true);
    // RDKit❗❌:   if (count > 1) {
    // RDKit❗❌:     fp.setBit(142);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (count > 0) {
    // RDKit❗❌:     fp.setBit(161);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_144, match, true)) {
    // RDKit❗❌:     fp.setBit(144);
    // RDKit❗❌:   }
    // RDKit❗❌:   count = RDKit::SubstructMatch(mol, *pats.bit_145, matches, true, true);
    // RDKit❗❌:   if (count > 1) {
    // RDKit❗❌:     fp.setBit(145);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (count > 0) {
    // RDKit❗❌:     fp.setBit(163);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_147, match, true)) {
    // RDKit❗❌:     fp.setBit(147);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_148, match, true)) {
    // RDKit❗❌:     fp.setBit(148);
    // RDKit❗❌:   }
    // RDKit❗❌:   count = RDKit::SubstructMatch(mol, *pats.bit_149, matches, true, true);
    // RDKit❗❌:   if (count > 1) {
    // RDKit❗❌:     fp.setBit(149);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (count > 0) {
    // RDKit❗❌:     fp.setBit(160);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_150, match, true)) {
    // RDKit❗❌:     fp.setBit(150);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_151, match, true)) {
    // RDKit❗❌:     fp.setBit(151);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_152, match, true)) {
    // RDKit❗❌:     fp.setBit(152);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_154, match, true)) {
    // RDKit❗❌:     fp.setBit(154);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_155, match, true)) {
    // RDKit❗❌:     fp.setBit(155);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_156, match, true)) {
    // RDKit❗❌:     fp.setBit(156);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_157, match, true)) {
    // RDKit❗❌:     fp.setBit(157);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_158, match, true)) {
    // RDKit❗❌:     fp.setBit(158);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_162, match, true)) {
    // RDKit❗❌:     fp.setBit(162);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (RDKit::SubstructMatch(mol, *pats.bit_165, match, true)) {
    // RDKit❗❌:     fp.setBit(165);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   /* BIT 125 */
    // RDKit❗❌:   RDKit::RingInfo *info = mol.getRingInfo();
    // RDKit❗❌:   unsigned int ringcount = info->numRings();
    // RDKit❗❌:   unsigned int nArom = 0;
    // RDKit❗❌:   for (unsigned int i = 0; i < ringcount; i++) {
    // RDKit❗❌:     bool isArom = true;
    // RDKit❗❌:     const std::vector<int> *ring = &info->bondRings()[i];
    // RDKit❗❌:     std::vector<int>::const_iterator iter;
    // RDKit❗❌:     for (iter = ring->begin(); iter != ring->end(); ++iter) {
    // RDKit❗❌:       if (!mol.getBondWithIdx(*iter)->getIsAromatic()) {
    // RDKit❗❌:         isArom = false;
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (isArom) {
    // RDKit❗❌:       if (nArom) {
    // RDKit❗❌:         fp.setBit(125);
    // RDKit❗❌:         break;
    // RDKit❗❌:       } else {
    // RDKit❗❌:         nArom++;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   /* BIT 166 */
    // RDKit❗❌:   std::vector<int> mapping;
    // RDKit❗❌:   if (RDKit::MolOps::getMolFrags(mol, mapping) > 1) {
    // RDKit❗❌:     fp.setBit(166);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Behavior: source direct element cases and all 136 exact patterns, count thresholds,
    // aromatic-ring and fragment keys, including source thresholds for 145/163
    // and 149/160 from pinned MACCS.cpp 857..876. Fixture bytes are unchanged.
    // Cost: O(1) bit writes use the single Fingerprint implementation, matching source.
    // Warm ring/valence contexts borrow data without clones. Cold input adds one
    // symmetrized SSSR/valence pass, explicitly worse than source initialized caches.
    // CORE component enumeration additionally owns component lists besides source mapping.
    // SEARCH performs all VF2 matching; no local matcher or duplicate graph algorithm.
    topology.validate()?;
    let pattern_matchers = cached_maccs_queries()?;
    let mut fp = Fingerprint::new(167);
    if topology.atoms.is_empty() {
        return Ok(fp);
    }
    let owned_rings;
    let ring_info = match rings {
        Some(v) => v,
        None => {
            owned_rings = cosmolkit_core::symmetrize_sssr_with_options_from_parts(
                topology.atoms.len(),
                &topology.bonds,
                &topology.adjacency,
                false,
                false,
            )?;
            &owned_rings
        }
    };
    let owned_valence;
    let valence = match valence {
        Some(v) => v,
        None => {
            owned_valence = cosmolkit_core::assign_valence_with_options_for_topology(
                topology,
                cosmolkit_core::ValenceModel::RdkitLike,
                false,
            )?;
            &owned_valence
        }
    };
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        topology,
        &coordinates,
        &topology.stereo_groups,
        Some(ring_info),
        Some(valence),
    );
    let context = build_prepared_query_match_context(topology, ring_info, valence)?;
    for atom in topology.atoms.iter() {
        match atom.atomic_number() {
            3 | 11 | 19 | 37 | 55 | 87 => {
                fp.set_bit(u32::try_from(35).expect("fixed MACCS bit"))?;
            }
            4 | 12 | 20 | 38 | 56 | 88 => {
                fp.set_bit(u32::try_from(10).expect("fixed MACCS bit"))?;
            }
            5 | 13 | 31 | 49 | 81 => {
                fp.set_bit(u32::try_from(18).expect("fixed MACCS bit"))?;
            }
            9 => {
                fp.set_bit(u32::try_from(42).expect("fixed MACCS bit"))?;
                fp.set_bit(u32::try_from(134).expect("fixed MACCS bit"))?;
            }
            15 => {
                fp.set_bit(u32::try_from(29).expect("fixed MACCS bit"))?;
            }
            16 => {
                fp.set_bit(u32::try_from(88).expect("fixed MACCS bit"))?;
            }
            17 => {
                fp.set_bit(u32::try_from(103).expect("fixed MACCS bit"))?;
                fp.set_bit(u32::try_from(134).expect("fixed MACCS bit"))?;
            }
            21 | 22 | 39 | 40 | 72 => {
                fp.set_bit(u32::try_from(5).expect("fixed MACCS bit"))?;
            }
            23 | 24 | 25 | 41 | 42 | 43 | 73 | 74 | 75 => {
                fp.set_bit(u32::try_from(7).expect("fixed MACCS bit"))?;
            }
            26 | 27 | 28 | 44 | 45 | 46 | 76 | 77 | 78 => {
                fp.set_bit(u32::try_from(9).expect("fixed MACCS bit"))?;
            }
            29 | 30 | 47 | 48 | 79 | 80 => {
                fp.set_bit(u32::try_from(12).expect("fixed MACCS bit"))?;
            }
            32 | 33 | 34 | 50 | 51 | 52 | 82 | 83 | 84 => {
                fp.set_bit(u32::try_from(3).expect("fixed MACCS bit"))?;
            }
            35 => {
                fp.set_bit(u32::try_from(46).expect("fixed MACCS bit"))?;
                fp.set_bit(u32::try_from(134).expect("fixed MACCS bit"))?;
            }
            53 => {
                fp.set_bit(u32::try_from(27).expect("fixed MACCS bit"))?;
                fp.set_bit(u32::try_from(134).expect("fixed MACCS bit"))?;
            }
            57..=71 => {
                fp.set_bit(u32::try_from(6).expect("fixed MACCS bit"))?;
            }
            89..=103 => {
                fp.set_bit(u32::try_from(4).expect("fixed MACCS bit"))?;
            }
            104 => {
                fp.set_bit(u32::try_from(2).expect("fixed MACCS bit"))?;
            }
            _ => {}
        }
    }

    for &(bit, ref matcher) in pattern_matchers
        .iter()
        .filter(|&&(bit, _)| bit <= 120 && bit != 118 && bit != 120)
    {
        if maccs_match_count(bit, pattern_matchers, &target, &context, 1)? > 0 {
            fp.set_bit(u32::try_from(bit).expect("fixed MACCS bit"))?;
        }
    }

    for bit in [118usize, 120] {
        if maccs_match_count(bit, pattern_matchers, &target, &context, 1000)? > 1 {
            fp.set_bit(u32::try_from(bit).expect("fixed MACCS bit"))?;
        }
    }
    let maccs_match_count = |bit: usize| -> Result<usize, MaccsFingerprintError> {
        maccs_match_count(bit, pattern_matchers, &target, &context, 1000)
    };
    let has_maccs_match = |bit: usize| -> Result<bool, MaccsFingerprintError> {
        Ok(self::maccs_match_count(bit, pattern_matchers, &target, &context, 1)? > 0)
    };
    for bit in [
        121usize, 122, 123, 126, 128, 129, 132, 133, 135, 137, 139, 144, 147, 148, 150, 151, 152,
        154, 155, 156, 157, 158, 162, 165,
    ] {
        if has_maccs_match(bit)? {
            fp.set_bit(u32::try_from(bit).expect("fixed MACCS bit"))?;
        }
    }

    let count = maccs_match_count(124)?;
    if count > 0 {
        fp.set_bit(u32::try_from(124).expect("fixed MACCS bit"))?;
    }
    if count > 1 {
        fp.set_bit(u32::try_from(130).expect("fixed MACCS bit"))?;
    }

    let count = maccs_match_count(127)?;
    if count > 1 {
        fp.set_bit(u32::try_from(127).expect("fixed MACCS bit"))?;
    }
    if count > 0 {
        fp.set_bit(u32::try_from(143).expect("fixed MACCS bit"))?;
    }

    if maccs_match_count(131)? > 1 {
        fp.set_bit(u32::try_from(131).expect("fixed MACCS bit"))?;
    }

    if maccs_match_count(136)? > 1 {
        fp.set_bit(u32::try_from(136).expect("fixed MACCS bit"))?;
    }

    let count = maccs_match_count(138)?;
    if count > 1 {
        fp.set_bit(u32::try_from(138).expect("fixed MACCS bit"))?;
    }
    if count > 0 {
        fp.set_bit(u32::try_from(153).expect("fixed MACCS bit"))?;
    }

    let count = maccs_match_count(140)?;
    if count > 3 {
        fp.set_bit(u32::try_from(140).expect("fixed MACCS bit"))?;
    }
    if count > 2 {
        fp.set_bit(u32::try_from(146).expect("fixed MACCS bit"))?;
    }
    if count > 1 {
        fp.set_bit(u32::try_from(159).expect("fixed MACCS bit"))?;
    }
    if count > 0 {
        fp.set_bit(u32::try_from(164).expect("fixed MACCS bit"))?;
    }

    if maccs_match_count(141)? > 2 {
        fp.set_bit(u32::try_from(141).expect("fixed MACCS bit"))?;
    }

    let count = maccs_match_count(142)?;
    if count > 1 {
        fp.set_bit(u32::try_from(142).expect("fixed MACCS bit"))?;
    }
    if count > 0 {
        fp.set_bit(u32::try_from(161).expect("fixed MACCS bit"))?;
    }

    let count = maccs_match_count(145)?;
    if count > 1 {
        fp.set_bit(u32::try_from(145).expect("fixed MACCS bit"))?;
    }
    if count > 0 {
        fp.set_bit(u32::try_from(163).expect("fixed MACCS bit"))?;
    }

    let count = maccs_match_count(149)?;
    if count > 1 {
        fp.set_bit(u32::try_from(149).expect("fixed MACCS bit"))?;
    }
    if count > 0 {
        fp.set_bit(u32::try_from(160).expect("fixed MACCS bit"))?;
    }

    let mut aromatic_ring_count = 0usize;
    for bond_ring in ring_info.bond_rings() {
        if bond_ring
            .iter()
            .all(|b| topology.bonds[b.index()].is_aromatic())
        {
            if aromatic_ring_count > 0 {
                fp.set_bit(u32::try_from(125).expect("fixed MACCS bit"))?;
                break;
            }
            aromatic_ring_count += 1;
        }
    }

    if cosmolkit_core::connected_components(topology)?
        .components
        .len()
        > 1
    {
        fp.set_bit(u32::try_from(166).expect("fixed MACCS bit"))?;
    }

    Ok(fp)
}

pub fn maccs_fingerprint(
    topology: &TopologyBlock,
    rings: Option<&RingInfo>,
    valence: Option<&ValenceAssignment>,
    params: &MaccsFingerprintParams,
) -> Result<Fingerprint, MaccsFingerprintError> {
    // RDKit❗❌: ExplicitBitVect *getFingerprintAsBitVect(const ROMol &mol) {
    // RDKit❗❌:   std::unique_ptr<ExplicitBitVect> fp(new ExplicitBitVect(167));
    // RDKit❗❌:   GenerateFP(mol, *fp);
    // RDKit❗❌:   return fp.release();
    // RDKit❗❌: }
    // COSMolKit✔️✔️:     if n_bits != COSMOLKIT_MACCS_PUBLIC_BITS {
    // COSMolKit✔️✔️:         return Err(FingerprintError::UnsupportedOption {
    // COSMolKit✔️✔️:             option: "MaccsFingerprintParams.n_bits",
    // COSMolKit✔️✔️:             reason: "RDKit MACCS exposes a fixed 167-bit raw vector with bit 0 unused; COSMolKit only exposes the exact 166-bit public projection",
    // COSMolKit✔️✔️:         });
    // COSMolKit✔️✔️:     }
    // Exact old public projection subtracts one from raw bits 1..166, with width166.
    if params.n_bits != COSMOLKIT_MACCS_PUBLIC_BITS {
        return Err(MaccsFingerprintError::UnsupportedOption {
            option: "MaccsFingerprintParams.n_bits",
            reason: "RDKit MACCS exposes a fixed 167-bit raw vector with bit 0 unused; COSMolKit only exposes the exact 166-bit public projection",
        });
    }
    if topology.atoms.is_empty() {
        return Ok(Fingerprint::new(166));
    }
    let raw = maccs_fingerprint_raw(topology, rings, valence)?;
    Ok(Fingerprint::from_on_bits(
        166,
        raw.on_bits().into_iter().map(|bit| bit - 1),
    )?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, Element};
    struct Fixture {
        topology: TopologyBlock,
    }
    impl Fixture {
        fn new() -> Self {
            Self {
                topology: TopologyBlock::default(),
            }
        }
        fn from_smiles(s: &str) -> Result<Self, cosmolkit_smiles::SmilesParseError> {
            Ok(Self {
                topology: cosmolkit_smiles::parse_smiles(s, &Default::default())?.topology,
            })
        }
    }
    fn maccs_fingerprint(
        input: &Fixture,
        p: &MaccsFingerprintParams,
    ) -> Result<Fingerprint, MaccsFingerprintError> {
        super::maccs_fingerprint(&input.topology, None, None, p)
    }
    fn maccs_get_fingerprint_as_bit_vect(
        input: &Fixture,
    ) -> Result<Fingerprint, MaccsFingerprintError> {
        super::maccs_fingerprint_raw(&input.topology, None, None)
    }
    fn rdkit_maccs_patterns_oracle() -> &'static [(usize, &'static str)] {
        &[
            (8usize, "[!#6!#1]1~*~*~*~1"),
            (11usize, "*1~*~*~*~1"),
            (13usize, "[#8]~[#7](~[#6])~[#6]"),
            (14usize, "[#16]-[#16]"),
            (15usize, "[#8]~[#6](~[#8])~[#8]"),
            (16usize, "[!#6!#1]1~*~*~1"),
            (17usize, "[#6]#[#6]"),
            (19usize, "*1~*~*~*~*~*~*~1"),
            (20usize, "[#14]"),
            (21usize, "[#6]=[#6](~[!#6!#1])~[!#6!#1]"),
            (22usize, "*1~*~*~1"),
            (23usize, "[#7]~[#6](~[#8])~[#8]"),
            (24usize, "[#7]-[#8]"),
            (25usize, "[#7]~[#6](~[#7])~[#7]"),
            (26usize, "[#6]=@[#6](@*)@*"),
            (28usize, "[!#6!#1]~[CH2]~[!#6!#1]"),
            (30usize, "[#6]~[!#6!#1](~[#6])(~[#6])~*"),
            (31usize, "[!#6!#1]~[F,Cl,Br,I]"),
            (32usize, "[#6]~[#16]~[#7]"),
            (33usize, "[#7]~[#16]"),
            (34usize, "[CH2]=*"),
            (36usize, "[#16R]"),
            (37usize, "[#7]~[#6](~[#8])~[#7]"),
            (38usize, "[#7]~[#6](~[#6])~[#7]"),
            (39usize, "[#8]~[#16](~[#8])~[#8]"),
            (40usize, "[#16]-[#8]"),
            (41usize, "[#6]#[#7]"),
            (43usize, "[!#6!#1!H0]~*~[!#6!#1!H0]"),
            (
                44usize,
                "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]",
            ),
            (45usize, "[#6]=[#6]~[#7]"),
            (47usize, "[#16]~*~[#7]"),
            (48usize, "[#8]~[!#6!#1](~[#8])~[#8]"),
            (49usize, "[!+0]"),
            (50usize, "[#6]=[#6](~[#6])~[#6]"),
            (51usize, "[#6]~[#16]~[#8]"),
            (52usize, "[#7]~[#7]"),
            (53usize, "[!#6!#1!H0]~*~*~*~[!#6!#1!H0]"),
            (54usize, "[!#6!#1!H0]~*~*~[!#6!#1!H0]"),
            (55usize, "[#8]~[#16]~[#8]"),
            (56usize, "[#8]~[#7](~[#8])~[#6]"),
            (57usize, "[#8R]"),
            (58usize, "[!#6!#1]~[#16]~[!#6!#1]"),
            (59usize, "[#16]!:*:*"),
            (60usize, "[#16]=[#8]"),
            (61usize, "*~[#16](~*)~*"),
            (62usize, "*@*!@*@*"),
            (63usize, "[#7]=[#8]"),
            (64usize, "*@*!@[#16]"),
            (65usize, "c:n"),
            (66usize, "[#6]~[#6](~[#6])(~[#6])~*"),
            (67usize, "[!#6!#1]~[#16]"),
            (68usize, "[!#6!#1!H0]~[!#6!#1!H0]"),
            (69usize, "[!#6!#1]~[!#6!#1!H0]"),
            (70usize, "[!#6!#1]~[#7]~[!#6!#1]"),
            (71usize, "[#7]~[#8]"),
            (72usize, "[#8]~*~*~[#8]"),
            (73usize, "[#16]=*"),
            (74usize, "[CH3]~*~[CH3]"),
            (75usize, "*!@[#7]@*"),
            (76usize, "[#6]=[#6](~*)~*"),
            (77usize, "[#7]~*~[#7]"),
            (78usize, "[#6]=[#7]"),
            (79usize, "[#7]~*~*~[#7]"),
            (80usize, "[#7]~*~*~*~[#7]"),
            (81usize, "[#16]~*(~*)~*"),
            (82usize, "*~[CH2]~[!#6!#1!H0]"),
            (83usize, "[!#6!#1]1~*~*~*~*~1"),
            (84usize, "[NH2]"),
            (85usize, "[#6]~[#7](~[#6])~[#6]"),
            (86usize, "[C;H2,H3][!#6!#1][C;H2,H3]"),
            (87usize, "[F,Cl,Br,I]!@*@*"),
            (89usize, "[#8]~*~*~*~[#8]"),
            (
                90usize,
                "[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]",
            ),
            (
                91usize,
                "[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]",
            ),
            (92usize, "[#8]~[#6](~[#7])~[#6]"),
            (93usize, "[!#6!#1]~[CH3]"),
            (94usize, "[!#6!#1]~[#7]"),
            (95usize, "[#7]~*~*~[#8]"),
            (96usize, "*1~*~*~*~*~1"),
            (97usize, "[#7]~*~*~*~[#8]"),
            (98usize, "[!#6!#1]1~*~*~*~*~*~1"),
            (99usize, "[#6]=[#6]"),
            (100usize, "*~[CH2]~[#7]"),
            (
                101usize,
                "[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]",
            ),
            (102usize, "[!#6!#1]~[#8]"),
            (104usize, "[!#6!#1!H0]~*~[CH2]~*"),
            (105usize, "*@*(@*)@*"),
            (106usize, "[!#6!#1]~*(~[!#6!#1])~[!#6!#1]"),
            (107usize, "[F,Cl,Br,I]~*(~*)~*"),
            (108usize, "[CH3]~*~*~*~[CH2]~*"),
            (109usize, "*~[CH2]~[#8]"),
            (110usize, "[#7]~[#6]~[#8]"),
            (111usize, "[#7]~*~[CH2]~*"),
            (112usize, "*~*(~*)(~*)~*"),
            (113usize, "[#8]!:*:*"),
            (114usize, "[CH3]~[CH2]~*"),
            (115usize, "[CH3]~*~[CH2]~*"),
            (116usize, "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]"),
            (117usize, "[#7]~*~[#8]"),
            (118usize, "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]"),
            (119usize, "[#7]=*"),
            (120usize, "[!#6R]"),
            (121usize, "[#7R]"),
            (122usize, "*~[#7](~*)~*"),
            (123usize, "[#8]~[#6]~[#8]"),
            (124usize, "[!#6!#1]~[!#6!#1]"),
            (126usize, "*!@[#8]!@*"),
            (127usize, "*@*!@[#8]"),
            (
                128usize,
                "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]",
            ),
            (
                129usize,
                "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[CH2R]1)]",
            ),
            (131usize, "[!#6!#1!H0]"),
            (132usize, "[#8]~*~[CH2]~*"),
            (133usize, "*@*!@[#7]"),
            (135usize, "[#7]!:*:*"),
            (136usize, "[#8]=*"),
            (137usize, "[!C!cR]"),
            (138usize, "[!#6!#1]~[CH2]~*"),
            (139usize, "[O!H0]"),
            (140usize, "[#8]"),
            (141usize, "[CH3]"),
            (142usize, "[#7]"),
            (144usize, "*!:*:*!:*"),
            (145usize, "*1~*~*~*~*~*~1"),
            (147usize, "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]"),
            (148usize, "*~[!#6!#1](~*)~*"),
            (149usize, "[C;H3,H4]"),
            (150usize, "*!@*@*!@*"),
            (151usize, "[#7!H0]"),
            (152usize, "[#8]~[#6](~[#6])~[#6]"),
            (154usize, "[#6]=[#8]"),
            (155usize, "*!@[CH2]!@*"),
            (156usize, "[#7]~*(~*)~*"),
            (157usize, "[#6]-[#8]"),
            (158usize, "[#6]-[#7]"),
            (162usize, "a"),
            (165usize, "[R]"),
        ]
    }

    #[test]
    fn maccs_pattern_table_matches_rdkit_source_patterns() {
        assert_eq!(RDKIT_MACCS_RAW_BITS, 167);
        assert_eq!(COSMOLKIT_MACCS_PUBLIC_BITS, 166);
        assert_eq!(rdkit_maccs_public_index(0), None);
        assert_eq!(rdkit_maccs_public_index(1), Some(0));
        assert_eq!(rdkit_maccs_public_index(166), Some(165));
        assert_eq!(rdkit_maccs_public_index(167), None);

        let expected = rdkit_maccs_patterns_oracle();
        assert_eq!(RDKIT_MACCS_PATTERNS.len(), expected.len());
        assert_eq!(expected.len(), 136);

        for (actual, &(expected_bit, expected_smarts)) in
            RDKIT_MACCS_PATTERNS.iter().zip(expected.iter())
        {
            assert_eq!(actual.bit, expected_bit);
            assert_eq!(actual.smarts, expected_smarts);

            let looked_up =
                rdkit_maccs_pattern(expected_bit).expect("RDKit MACCS pattern bit is present");
            assert_eq!(looked_up.bit, expected_bit);
            assert_eq!(looked_up.smarts, expected_smarts);
        }

        for bit in 0..=RDKIT_MACCS_RAW_BITS {
            let expected_entry = expected
                .iter()
                .find(|&&(expected_bit, _)| expected_bit == bit);
            assert_eq!(
                rdkit_maccs_pattern(bit).map(|pattern| (pattern.bit, pattern.smarts)),
                expected_entry.copied(),
                "pattern lookup mismatch for RDKit MACCS bit {bit}"
            );
        }
    }

    fn single_atom_molecule(atomic_number: u8) -> Fixture {
        let element = Element::from_atomic_number(atomic_number).unwrap();
        let atoms = vec![cosmolkit_model::Atom::from_spec(
            cosmolkit_model::AtomId::new(0),
            AtomSpec::new(element),
        )];
        Fixture {
            topology: TopologyBlock {
                adjacency: cosmolkit_model::AdjacencyList::from_topology(atoms.len(), &[]),
                atoms,
                ..Default::default()
            },
        }
    }

    fn assert_maccs_keys_001_040(smiles: &str, expected_public_bits: &[usize]) {
        let mol = Fixture::from_smiles(smiles).unwrap_or_else(|err| {
            panic!("MACCS keys 001-040 fixture {smiles} should parse: {err}")
        });
        assert_maccs_keys_001_040_for_mol(smiles, &mol, expected_public_bits);
    }

    fn assert_maccs_keys_001_040_for_mol(
        label: &str,
        mol: &Fixture,
        expected_public_bits: &[usize],
    ) {
        let params = MaccsFingerprintParams::default();
        let actual: Vec<usize> = maccs_fingerprint(mol, &params)
            .expect("MACCS fingerprint")
            .on_bits()
            .into_iter()
            .map(|bit| bit as usize)
            .collect::<Vec<_>>()
            .into_iter()
            .filter(|&bit| bit < 40)
            .collect();
        assert_eq!(
            actual, expected_public_bits,
            "RDKit MACCS keys 001-040 public projection mismatch for {label}"
        );
    }

    fn assert_maccs_keys_041_080(smiles: &str, expected_public_bits: &[usize]) {
        let mol = Fixture::from_smiles(smiles).unwrap_or_else(|err| {
            panic!("MACCS keys 041-080 fixture {smiles} should parse: {err}")
        });
        let params = MaccsFingerprintParams::default();
        let actual: Vec<usize> = maccs_fingerprint(&mol, &params)
            .expect("MACCS fingerprint")
            .on_bits()
            .into_iter()
            .map(|bit| bit as usize)
            .collect::<Vec<_>>()
            .into_iter()
            .filter(|&bit| (40..80).contains(&bit))
            .collect();
        assert_eq!(
            actual, expected_public_bits,
            "RDKit MACCS keys 041-080 public projection mismatch for {smiles}"
        );
    }

    fn assert_maccs_keys_081_120(smiles: &str, expected_public_bits: &[usize]) {
        let mol = Fixture::from_smiles(smiles).unwrap_or_else(|err| {
            panic!("MACCS keys 081-120 fixture {smiles} should parse: {err}")
        });
        let params = MaccsFingerprintParams::default();
        let actual: Vec<usize> = maccs_fingerprint(&mol, &params)
            .expect("MACCS fingerprint")
            .on_bits()
            .into_iter()
            .map(|bit| bit as usize)
            .collect::<Vec<_>>()
            .into_iter()
            .filter(|&bit| (80..120).contains(&bit))
            .collect();
        assert_eq!(
            actual, expected_public_bits,
            "RDKit MACCS keys 081-120 public projection mismatch for {smiles}"
        );
    }

    fn assert_maccs_keys_121_166(smiles: &str, expected_public_bits: &[usize]) {
        let mol = Fixture::from_smiles(smiles).unwrap_or_else(|err| {
            panic!("MACCS keys 121-166 fixture {smiles} should parse: {err}")
        });
        let params = MaccsFingerprintParams::default();
        let actual: Vec<usize> = maccs_fingerprint(&mol, &params)
            .expect("MACCS fingerprint")
            .on_bits()
            .into_iter()
            .map(|bit| bit as usize)
            .collect::<Vec<_>>()
            .into_iter()
            .filter(|&bit| (120..166).contains(&bit))
            .collect();
        assert_eq!(
            actual, expected_public_bits,
            "RDKit MACCS keys 121-166 public projection mismatch for {smiles}"
        );
    }

    fn assert_maccs_full_vector_for_mol(
        label: &str,
        mol: &Fixture,
        expected_raw_bits: &[usize],
        expected_public_bits: &[usize],
    ) {
        let raw = maccs_get_fingerprint_as_bit_vect(mol).expect("raw MACCS fingerprint");
        assert_eq!(raw.n_bits() as usize, RDKIT_MACCS_RAW_BITS);
        assert!(
            !raw.on_bits()
                .into_iter()
                .map(|bit| bit as usize)
                .collect::<Vec<_>>()
                .contains(&0),
            "RDKit MACCS raw bit 0 must stay unused for {label}"
        );
        assert_eq!(
            raw.on_bits()
                .into_iter()
                .map(|bit| bit as usize)
                .collect::<Vec<_>>(),
            expected_raw_bits,
            "RDKit MACCS raw 167-bit vector mismatch for {label}"
        );

        let params = MaccsFingerprintParams::default();
        let public = maccs_fingerprint(mol, &params).expect("MACCS fingerprint");
        assert_eq!(public.n_bits() as usize, COSMOLKIT_MACCS_PUBLIC_BITS);
        assert_eq!(
            public
                .on_bits()
                .into_iter()
                .map(|bit| bit as usize)
                .collect::<Vec<_>>(),
            expected_public_bits,
            "COSMolKit MACCS public 166-bit projection mismatch for {label}"
        );

        let projected_from_raw: Vec<usize> = raw
            .on_bits()
            .into_iter()
            .map(|bit| bit as usize)
            .collect::<Vec<_>>()
            .into_iter()
            .filter_map(rdkit_maccs_public_index)
            .collect();
        assert_eq!(
            projected_from_raw, expected_public_bits,
            "RDKit raw-to-public MACCS projection mismatch for {label}"
        );
    }

    fn assert_maccs_full_vector(
        label: &str,
        smiles: &str,
        expected_raw_bits: &[usize],
        expected_public_bits: &[usize],
    ) {
        let mol = Fixture::from_smiles(smiles)
            .unwrap_or_else(|err| panic!("MACCS full-vector fixture {smiles} should parse: {err}"));
        assert_maccs_full_vector_for_mol(label, &mol, expected_raw_bits, expected_public_bits);
    }

    #[test]
    fn maccs_keys_001_040_direct_element_keys_match_rdkit() {
        let fixtures: &[(&str, u8, &[usize])] = &[
            ("key_002_rf", 104, &[1]),
            ("key_003_ge", 32, &[2]),
            ("key_004_ac", 89, &[3]),
            ("key_005_sc", 21, &[4]),
            ("key_006_la", 57, &[5]),
            ("key_007_v", 23, &[6]),
            ("key_009_fe", 26, &[8]),
            ("key_010_be", 4, &[9]),
            ("key_012_cu", 29, &[11]),
            ("key_018_boron", 5, &[17]),
            ("key_020_silicon", 14, &[19]),
            ("key_027_iodine", 53, &[26]),
            ("key_029_phosphorus", 15, &[28]),
            ("key_035_lithium", 3, &[34]),
        ];
        for &(label, atomic_number, expected_public_bits) in fixtures {
            let mol = single_atom_molecule(atomic_number);
            assert_maccs_keys_001_040_for_mol(label, &mol, expected_public_bits);
        }

        assert_maccs_keys_001_040("C", &[]);
    }

    #[test]
    fn maccs_fingerprint_full_bit_vectors_match_rdkit_raw_and_public_projection() {
        let empty = Fixture::new();
        assert_maccs_full_vector_for_mol("empty", &empty, &[], &[]);

        let fixtures: &[(&str, &str, &[usize], &[usize])] = &[
            ("methane", "C", &[160], &[159]),
            ("fluorine_atom", "F", &[42, 134], &[41, 133]),
            ("sulfur_atom", "S", &[88], &[87]),
            ("chlorine_atom", "Cl", &[103, 134], &[102, 133]),
            (
                "salt_fragments",
                "CCO.Cl",
                &[
                    82, 103, 109, 114, 131, 134, 139, 153, 155, 157, 160, 164, 166,
                ],
                &[
                    81, 102, 108, 113, 130, 133, 138, 152, 154, 156, 159, 163, 165,
                ],
            ),
            ("benzene", "c1ccccc1", &[162, 163, 165], &[161, 162, 164]),
            (
                "biphenyl",
                "c1ccccc1c2ccccc2",
                &[62, 125, 145, 162, 163, 165],
                &[61, 124, 144, 161, 162, 164],
            ),
            (
                "pyridine",
                "c1ncccc1",
                &[65, 98, 121, 137, 161, 162, 163, 165],
                &[64, 97, 120, 136, 160, 161, 162, 164],
            ),
            (
                "morpholine",
                "O1CCNCC1",
                &[
                    57, 82, 86, 91, 95, 98, 100, 104, 109, 111, 118, 120, 121, 128, 129, 132, 137,
                    138, 147, 151, 153, 157, 158, 161, 163, 164, 165,
                ],
                &[
                    56, 81, 85, 90, 94, 97, 99, 103, 108, 110, 117, 119, 120, 127, 128, 131, 136,
                    137, 146, 150, 152, 156, 157, 160, 162, 163, 164,
                ],
            ),
            ("ammonium", "[NH4+]", &[49, 151, 161], &[48, 150, 160]),
            (
                "acetate",
                "CC(=O)[O-]",
                &[49, 123, 154, 157, 159, 160, 164],
                &[48, 122, 153, 156, 158, 159, 163],
            ),
            ("isotopic_methane", "[13CH4]", &[160], &[159]),
            (
                "nitro",
                "N=O",
                &[63, 69, 71, 94, 102, 119, 124, 151, 161, 164],
                &[62, 68, 70, 93, 101, 118, 123, 150, 160, 163],
            ),
            (
                "cyclopropanol",
                "C1CC1O",
                &[22, 90, 104, 127, 132, 139, 143, 147, 152, 157, 164, 165],
                &[21, 89, 103, 126, 131, 138, 142, 146, 151, 156, 163, 164],
            ),
            (
                "fragment_methanes",
                "C.C",
                &[149, 160, 166],
                &[148, 159, 165],
            ),
            (
                "all_key_low_mix",
                "NCCO",
                &[
                    54, 82, 84, 95, 100, 104, 109, 111, 118, 131, 132, 138, 139, 147, 151, 153,
                    155, 157, 158, 161, 164,
                ],
                &[
                    53, 81, 83, 94, 99, 103, 108, 110, 117, 130, 131, 137, 138, 146, 150, 152, 154,
                    156, 157, 160, 163,
                ],
            ),
            (
                "all_key_high_mix",
                "OCOCOCO",
                &[
                    28, 82, 86, 89, 90, 109, 123, 126, 128, 131, 138, 139, 140, 146, 153, 155, 157,
                    159, 164,
                ],
                &[
                    27, 81, 85, 88, 89, 108, 122, 125, 127, 130, 137, 138, 139, 145, 152, 154, 156,
                    158, 163,
                ],
            ),
        ];

        for &(label, smiles, expected_raw_bits, expected_public_bits) in fixtures {
            assert_maccs_full_vector(label, smiles, expected_raw_bits, expected_public_bits);
        }

        let err = maccs_fingerprint(
            &Fixture::from_smiles("NCCO").unwrap(),
            &MaccsFingerprintParams { n_bits: 64 },
        )
        .unwrap_err();
        assert!(matches!(
            err,
            MaccsFingerprintError::UnsupportedOption {
                option: "MaccsFingerprintParams.n_bits",
                ..
            }
        ));
    }

    #[test]
    fn maccs_keys_001_040_pattern_keys_match_rdkit() {
        let fixtures: &[(&str, &str, &[usize])] = &[
            ("key_008_hetero_four_ring", "O1CCC1", &[7, 10]),
            ("key_011_four_ring", "C1CCC1", &[10]),
            ("key_013_o_n_c_c", "ON(C)C", &[12, 23]),
            ("key_014_disulfide", "CSSC", &[13]),
            ("key_015_o_c_o_o", "O=C(O)O", &[14]),
            ("key_016_hetero_three_ring", "O1CC1", &[15, 21]),
            ("key_017_alkyne", "C#C", &[16]),
            ("key_019_seven_ring", "C1CCCCCC1", &[18]),
            ("key_021_alkene_dihetero", "C=C(O)O", &[20, 33]),
            ("key_022_three_ring", "C1CC1", &[21]),
            ("key_023_n_c_o_o", "NC(=O)O", &[22]),
            ("key_024_n_o", "ON(C)C", &[12, 23]),
            ("key_025_n_c_n_n", "NC(N)N", &[24]),
            ("key_026_cyclic_alkene", "C1=C2CCCC2C1", &[10, 18, 25]),
            ("key_028_hetero_ch2_hetero", "OCO", &[27]),
            ("key_030_c_hetero_c_c_any", "C[S](C)(C)C", &[29]),
            ("key_031_hetero_halogen", "N[Pt](Cl)(Cl)N", &[8, 30]),
            ("key_032_c_s_n", "CSN", &[31, 32]),
            ("key_033_n_s", "CSN", &[31, 32]),
            ("key_034_ch2_double", "C=C", &[33]),
            ("key_036_s_ring", "S1CC1", &[15, 21, 35]),
            ("key_037_n_c_o_n", "NC(=O)N", &[36]),
            ("key_038_n_c_c_n", "NC(C)N", &[37]),
            ("key_039_o_s_o_o", "COS(=O)(=O)O", &[38, 39]),
            ("key_040_s_o", "CSO", &[39]),
        ];
        for &(label, smiles, expected_public_bits) in fixtures {
            assert_maccs_keys_001_040(smiles, expected_public_bits);
            assert!(
                expected_public_bits.iter().any(|&bit| bit + 1 <= 40),
                "{label} should exercise at least one RDKit raw key in 1..=40"
            );
        }
    }

    #[test]
    fn maccs_keys_041_080_match_rdkit() {
        let fixtures: &[(&str, &str, &[usize])] = &[
            ("key_041_c_n_triple", "C#N", &[40]),
            ("key_042_fluorine", "F", &[41]),
            ("key_043_hetero_bridge_h", "OCO", &[42]),
            ("key_044_exotic_element", "[SeH2]", &[43]),
            ("key_045_c_c_n", "C=CN", &[44]),
            ("key_046_bromine", "Br", &[45]),
            ("key_047_s_x_n", "SCN", &[42, 46]),
            (
                "key_048_o_hetero_o_o",
                "COS(=O)(=O)O",
                &[47, 54, 57, 59, 60, 66, 68, 72],
            ),
            ("key_049_charged", "[NH4+]", &[48]),
            ("key_050_substituted_alkene", "CC(C)=C", &[49, 73, 75]),
            ("key_051_c_s_o", "CSO", &[50, 66, 68]),
            ("key_052_n_n", "NNO", &[42, 51, 67, 68, 69, 70]),
            ("key_053_hetero_bridge_3", "NCCCO", &[52]),
            ("key_054_hetero_bridge_2", "NCCN", &[53, 78]),
            (
                "key_055_o_s_o",
                "CS(=O)(=O)C",
                &[50, 54, 57, 59, 60, 66, 72, 73],
            ),
            ("key_056_o_n_o_c", "ON(O)C", &[42, 55, 68, 69, 70]),
            ("key_057_o_ring", "O1CC1", &[56]),
            (
                "key_058_hetero_s_hetero",
                "CS(=O)(=O)C",
                &[50, 54, 57, 59, 60, 66, 72, 73],
            ),
            ("key_059_s_aromatic_chain", "Sc1ccccc1", &[58, 63]),
            ("key_060_s_o_double", "CS(=O)C", &[50, 59, 60, 66, 72, 73]),
            ("key_061_s_three_neighbors", "C[S](C)(C)C", &[60, 73]),
            ("key_062_ring_nonring_ring", "C1CC1C1CC1", &[61]),
            ("key_063_n_o_double", "N=O", &[62, 68, 70]),
            ("key_064_ring_to_s", "Sc1ccccc1", &[58, 63]),
            ("key_065_aromatic_c_n", "c1ncccc1", &[64]),
            ("key_066_quaternary_carbon", "CC(C)(C)C", &[65, 73]),
            ("key_067_hetero_s", "CSSC", &[66]),
            (
                "key_068_hetero_h_hetero_h",
                "NNO",
                &[42, 51, 67, 68, 69, 70],
            ),
            ("key_069_hetero_hetero_h", "CSN", &[66, 68]),
            ("key_070_hetero_n_hetero", "NNO", &[42, 51, 67, 68, 69, 70]),
            ("key_071_n_o", "N=O", &[62, 68, 70]),
            ("key_072_o_x_x_o", "OCCO", &[53, 71]),
            ("key_073_s_double", "CS(=O)C", &[50, 59, 60, 66, 72, 73]),
            ("key_074_methyl_bridge", "CCC", &[73]),
            ("key_075_exocyclic_n_ring", "CN1CC1", &[74]),
            ("key_076_substituted_alkene", "CC(C)=C", &[49, 73, 75]),
            ("key_077_n_x_n", "NC(=O)N", &[42, 76]),
            ("key_078_c_n_double", "C=N", &[77]),
            ("key_079_n_x_x_n", "NCCN", &[53, 78]),
            ("key_080_n_x_x_x_n", "NCCCN", &[52, 79]),
        ];
        for &(label, smiles, expected_public_bits) in fixtures {
            assert_maccs_keys_041_080(smiles, expected_public_bits);
            assert!(
                expected_public_bits
                    .iter()
                    .any(|&bit| (40..80).contains(&bit)),
                "{label} should exercise at least one RDKit raw key in 41..=80"
            );
        }
    }

    #[test]
    fn maccs_keys_081_120_match_rdkit() {
        let fixtures: &[(&str, &str, &[usize])] = &[
            (
                "key_081_s_three_neighbors",
                "C[S](C)(C)C",
                &[85, 87, 92, 111],
            ),
            ("key_082_ch2_hetero_h", "CCN", &[81, 83, 99, 113]),
            (
                "key_083_hetero_five_ring",
                "O1CCCC1",
                &[82, 85, 95, 108, 117],
            ),
            ("key_084_nh2", "N", &[]),
            ("key_085_c_n_c_c", "CN(C)C", &[84, 85, 92]),
            ("key_086_c_h2_h3_hetero_c", "CN(C)C", &[84, 85, 92]),
            ("key_087_halogen_ring_chain", "FC1CC1", &[86, 106]),
            ("key_088_sulfur", "S", &[87]),
            ("key_089_o_bridge_3", "OCCCO", &[81, 88, 89, 103, 108, 117]),
            (
                "key_090_hetero_ch2_bridge",
                "NCCO",
                &[81, 83, 94, 99, 103, 108, 110, 117],
            ),
            (
                "key_091_hetero_ch2_bridge_4",
                "NCCCO",
                &[81, 83, 89, 96, 99, 103, 108, 110, 117],
            ),
            ("key_092_o_c_n_c", "OC(N)C", &[83, 91, 109, 116]),
            ("key_093_hetero_methyl", "CN", &[83, 92]),
            ("key_094_hetero_n", "CN", &[83, 92]),
            (
                "key_095_n_bridge_o",
                "NCCO",
                &[81, 83, 94, 99, 103, 108, 110, 117],
            ),
            ("key_096_five_ring", "C1CCCC1", &[95, 117]),
            (
                "key_097_n_bridge3_o",
                "NCCCO",
                &[81, 83, 89, 96, 99, 103, 108, 110, 117],
            ),
            ("key_098_hetero_six_ring", "O1CCCCC1", &[85, 97, 108, 117]),
            ("key_099_alkene", "C=C", &[98]),
            ("key_100_ch2_n", "CCN", &[81, 83, 99, 113]),
            ("key_101_large_ring", "C1CCCCCCC1", &[100, 117]),
            ("key_102_hetero_o", "CO", &[92]),
            ("key_103_chlorine", "Cl", &[102]),
            ("key_104_hetero_ch2_chain", "NCC", &[81, 83, 99, 113]),
            ("key_105_ring_branch_ring", "C1CC(C1)C1CC1", &[117]),
            (
                "key_106_hetero_three_neighbors",
                "N(O)(O)O",
                &[93, 101, 105],
            ),
            (
                "key_107_halogen_three_neighbors",
                "FC(F)(F)F",
                &[105, 106, 111],
            ),
            ("key_108_methyl_chain_ch2", "CCCC", &[113, 114, 117]),
            ("key_109_ch2_o", "CCO", &[81, 108, 113]),
            ("key_110_n_c_o", "NCO", &[81, 83, 99, 108, 109, 116]),
            ("key_111_n_ch2_chain", "NCC", &[81, 83, 99, 113]),
            ("key_112_quaternary_any", "CC(C)(C)C", &[111]),
            ("key_113_o_aromatic_chain", "Oc1ccccc1", &[112]),
            ("key_114_methyl_ch2_any", "CCC", &[113]),
            ("key_115_methyl_any_ch2_any", "CC(C)C", &[]),
            ("key_116_methyl_ch2_bridge", "CCCC", &[113, 114, 117]),
            ("key_117_n_x_o", "NCO", &[81, 83, 99, 108, 109, 116]),
            ("key_118_two_ch2_paths", "CCCC", &[113, 114, 117]),
            ("key_119_n_double_any", "N=O", &[93, 101, 118]),
            (
                "key_120_two_noncarbon_ring_atoms",
                "N1CCO1",
                &[81, 89, 93, 94, 99, 101, 103, 108, 110, 117, 119],
            ),
        ];
        for &(label, smiles, expected_public_bits) in fixtures {
            assert_maccs_keys_081_120(smiles, expected_public_bits);
            if !expected_public_bits.is_empty() {
                assert!(
                    expected_public_bits
                        .iter()
                        .any(|&bit| (80..120).contains(&bit)),
                    "{label} should exercise at least one RDKit raw key in 81..=120"
                );
            }
        }
    }

    #[test]
    fn maccs_keys_121_166_match_rdkit() {
        let fixtures: &[(&str, &str, &[usize])] = &[
            (
                "key_121_n_ring",
                "N1CCCC1",
                &[120, 128, 136, 137, 146, 150, 152, 157, 160, 164],
            ),
            (
                "key_122_n_three_neighbors",
                "CN(C)C",
                &[121, 140, 147, 148, 157, 159, 160],
            ),
            ("key_123_o_c_o", "COC", &[125, 148, 156, 159, 163]),
            ("key_124_two_hetero", "NO", &[123, 130, 138, 150, 160, 163]),
            (
                "key_125_two_aromatic_rings",
                "c1ccccc1c2ccccc2",
                &[124, 144, 161, 162, 164],
            ),
            ("key_126_o_nonring", "COC", &[125, 148, 156, 159, 163]),
            (
                "key_127_ring_nonring_o",
                "C1CC1O",
                &[126, 131, 138, 142, 146, 151, 156, 163, 164],
            ),
            (
                "key_128_ch2_bridge_5",
                "CCCCCCC",
                &[127, 128, 146, 148, 154, 159],
            ),
            ("key_129_ch2_bridge_4", "CCCCCC", &[128, 146, 148, 154, 159]),
            (
                "key_130_two_hetero_pairs",
                "NOON",
                &[123, 125, 129, 130, 141, 150, 158, 160, 163],
            ),
            (
                "key_131_two_hetero_h",
                "NNO",
                &[123, 129, 130, 138, 141, 150, 160, 163],
            ),
            (
                "key_132_o_ch2_chain",
                "OCCC",
                &[131, 138, 146, 152, 154, 156, 159, 163],
            ),
            (
                "key_133_ring_nonring_n",
                "C1CC1N",
                &[132, 146, 150, 155, 157, 160, 164],
            ),
            ("key_134_halogen", "F", &[133]),
            (
                "key_135_n_aromatic_chain",
                "Nc1ccccc1",
                &[132, 134, 150, 155, 157, 160, 161, 162, 164],
            ),
            ("key_136_two_o_double", "O=CC=O", &[135, 153, 158, 163]),
            (
                "key_137_noncarbon_ring",
                "O1CC1",
                &[136, 146, 152, 156, 163, 164],
            ),
            (
                "key_138_two_hetero_ch2",
                "NCCO",
                &[130, 131, 137, 138, 146, 150, 152, 154, 156, 157, 160, 163],
            ),
            ("key_139_o_no_h", "COC", &[125, 148, 156, 159, 163]),
            (
                "key_140_four_oxygens",
                "OCOCOCO",
                &[
                    122, 125, 127, 130, 137, 138, 139, 145, 152, 154, 156, 158, 163,
                ],
            ),
            ("key_141_three_methyl", "CC(C)(C)C", &[140, 148, 159]),
            ("key_142_two_nitrogens", "NN", &[123, 130, 141, 150, 160]),
            (
                "key_143_ring_nonring_o_once",
                "C1CC1O",
                &[126, 131, 138, 142, 146, 151, 156, 163, 164],
            ),
            ("key_144_aromatic_chain_three", "c1ccccc1", &[161, 162, 164]),
            (
                "key_145_two_six_rings",
                "C1CCCCC1C1CCCCC1",
                &[127, 128, 144, 146, 162, 164],
            ),
            (
                "key_146_three_oxygens",
                "OCOCO",
                &[122, 125, 130, 137, 138, 145, 152, 154, 156, 158, 163],
            ),
            ("key_147_two_ch2", "CCCC", &[146, 148, 154, 159]),
            (
                "key_148_hetero_three_neighbors",
                "N(C)(C)C",
                &[121, 140, 147, 148, 157, 159, 160],
            ),
            ("key_149_two_methyl", "CCC", &[148, 154, 159]),
            (
                "key_150_nonring_ring_path",
                "C1CC1CC1CC1",
                &[127, 128, 146, 154, 164],
            ),
            (
                "key_151_n_no_h",
                "CN(C)C",
                &[121, 140, 147, 148, 157, 159, 160],
            ),
            ("key_152_o_c_c_c", "OC(C)C", &[138, 148, 151, 156, 159, 163]),
            (
                "key_153_hetero_ch2_once",
                "CCN",
                &[150, 152, 154, 157, 159, 160],
            ),
            ("key_154_carbonyl", "C=O", &[153, 163]),
            ("key_155_ch2_nonring", "CCC", &[148, 154, 159]),
            (
                "key_156_n_three_neighbors",
                "N(C)(C)C",
                &[121, 140, 147, 148, 157, 159, 160],
            ),
            ("key_157_c_o_single", "CO", &[138, 156, 159, 163]),
            ("key_158_c_n_single", "CN", &[150, 157, 159, 160]),
            (
                "key_159_two_oxygens",
                "OCO",
                &[122, 130, 138, 152, 154, 156, 158, 163],
            ),
            ("key_160_methyl_once", "CC", &[148, 159]),
            ("key_161_n_once", "CN", &[150, 157, 159, 160]),
            ("key_162_aromatic", "c1ccccc1", &[161, 162, 164]),
            (
                "key_163_six_ring_once",
                "C1CCCCC1",
                &[127, 128, 146, 162, 164],
            ),
            ("key_164_o_once", "CO", &[138, 156, 159, 163]),
            ("key_165_ring", "C1CC1", &[146, 164]),
            ("key_166_fragments", "C.C", &[148, 159, 165]),
        ];
        for &(label, smiles, expected_public_bits) in fixtures {
            assert_maccs_keys_121_166(smiles, expected_public_bits);
            assert!(
                expected_public_bits
                    .iter()
                    .any(|&bit| (120..166).contains(&bit)),
                "{label} should exercise at least one RDKit raw key in 121..=166"
            );
        }
    }
}
