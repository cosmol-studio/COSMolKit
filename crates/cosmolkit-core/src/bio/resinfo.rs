//! Gemmi-derived residue information table.
//!
//! Source: `third_party/gemmi/src/resinfo.cpp` and `third_party/gemmi/include/gemmi/resinfo.hpp`.

#![allow(non_camel_case_types)]

use std::{convert::Infallible, fmt, str::FromStr};

use serde::{Deserialize, Deserializer, Serialize, Serializer};
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u16)]
pub enum ResidueCode {
    ALA = 0,
    ARG = 1,
    ASN = 2,
    ABA = 3,
    ASP = 4,
    ASX = 5,
    CYS = 6,
    CSH = 7,
    GLN = 8,
    GLU = 9,
    GLX = 10,
    GLY = 11,
    HIS = 12,
    ILE = 13,
    LEU = 14,
    LYS = 15,
    MET = 16,
    MSE = 17,
    ORN = 18,
    PHE = 19,
    PRO = 20,
    SER = 21,
    THR = 22,
    TRP = 23,
    TYR = 24,
    UNK = 25,
    VAL = 26,
    SEC = 27,
    PYL = 28,
    SEP = 29,
    TPO = 30,
    PCA = 31,
    CSO = 32,
    PTR = 33,
    KCX = 34,
    CSD = 35,
    LLP = 36,
    CME = 37,
    MLY = 38,
    DAL = 39,
    TYS = 40,
    OCS = 41,
    M3L = 42,
    FME = 43,
    ALY = 44,
    HYP = 45,
    CAS = 46,
    CRO = 47,
    CSX = 48,
    DPR = 49,
    DGL = 50,
    DVA = 51,
    CSS = 52,
    DPN = 53,
    DSN = 54,
    DLE = 55,
    HIC = 56,
    NLE = 57,
    MVA = 58,
    MLZ = 59,
    CR2 = 60,
    SAR = 61,
    DAR = 62,
    DLY = 63,
    YCM = 64,
    NRQ = 65,
    CGU = 66,
    R0TD = 67,
    MLE = 68,
    DAS = 69,
    DTR = 70,
    CXM = 71,
    TPQ = 72,
    DCY = 73,
    DSG = 74,
    DTY = 75,
    DHI = 76,
    MEN = 77,
    DTH = 78,
    SAC = 79,
    DGN = 80,
    AIB = 81,
    SMC = 82,
    IAS = 83,
    CIR = 84,
    BMT = 85,
    DIL = 86,
    FGA = 87,
    PHI = 88,
    CRQ = 89,
    SME = 90,
    GHP = 91,
    MHO = 92,
    NEP = 93,
    TRQ = 94,
    TOX = 95,
    ALC = 96,
    R3FG = 97,
    SCH = 98,
    MDO = 99,
    MAA = 100,
    GYS = 101,
    MK8 = 102,
    CR8 = 103,
    KPI = 104,
    SCY = 105,
    DHA = 106,
    OMY = 107,
    CAF = 108,
    R0AF = 109,
    SNN = 110,
    MHS = 111,
    MLU = 112,
    SNC = 113,
    PHD = 114,
    B3E = 115,
    MEA = 116,
    MED = 117,
    OAS = 118,
    GL3 = 119,
    FVA = 120,
    PHL = 121,
    CRF = 122,
    OMZ = 123,
    BFD = 124,
    MEQ = 125,
    DAB = 126,
    AGM = 127,
    PSU = 128,
    R5MU = 129,
    R7MG = 130,
    OMG = 131,
    UR3 = 132,
    OMC = 133,
    R2MG = 134,
    H2U = 135,
    R4SU = 136,
    OMU = 137,
    R4OC = 138,
    MA6 = 139,
    M2G = 140,
    R1MA = 141,
    R6MZ = 142,
    CCC = 143,
    R2MA = 144,
    R1MG = 145,
    R5BU = 146,
    MIA = 147,
    DOC = 148,
    R8OG = 149,
    R5CM = 150,
    R3DR = 151,
    BRU = 152,
    CBR = 153,
    HOH = 154,
    DOD = 155,
    HEM = 156,
    SO4 = 157,
    GOL = 158,
    EDO = 159,
    NAG = 160,
    PO4 = 161,
    ACT = 162,
    PEG = 163,
    MAN = 164,
    FAD = 165,
    BMA = 166,
    ADP = 167,
    DMS = 168,
    ACE = 169,
    NH2 = 170,
    MPD = 171,
    MES = 172,
    NAD = 173,
    NAP = 174,
    TRS = 175,
    ATP = 176,
    PG4 = 177,
    GDP = 178,
    FUC = 179,
    FMT = 180,
    GAL = 181,
    PGE = 182,
    FMN = 183,
    PLP = 184,
    EPE = 185,
    SF4 = 186,
    BME = 187,
    CIT = 188,
    BE7 = 189,
    MRD = 190,
    MHA = 191,
    BU3 = 192,
    PGO = 193,
    BU2 = 194,
    PDO = 195,
    BU1 = 196,
    PG6 = 197,
    R1BO = 198,
    PE7 = 199,
    PG5 = 200,
    TFP = 201,
    DHD = 202,
    PEU = 203,
    TAU = 204,
    SBT = 205,
    SAL = 206,
    IOH = 207,
    IPA = 208,
    PIG = 209,
    B3P = 210,
    BTB = 211,
    NHE = 212,
    C8E = 213,
    OTE = 214,
    PE4 = 215,
    XPE = 216,
    PE8 = 217,
    P33 = 218,
    N8E = 219,
    R2OS = 220,
    R1PS = 221,
    CPS = 222,
    DMX = 223,
    MPO = 224,
    GCD = 225,
    DXG = 226,
    CM5 = 227,
    ACA = 228,
    ACN = 229,
    CCN = 230,
    GLC = 231,
    DR6 = 232,
    NH4 = 233,
    AZI = 234,
    BNG = 235,
    BOG = 236,
    BGC = 237,
    BCN = 238,
    BRO = 239,
    CAC = 240,
    CBX = 241,
    ACY = 242,
    CBM = 243,
    CLO = 244,
    R3CO = 245,
    NCO = 246,
    CU1 = 247,
    CYN = 248,
    MA4 = 249,
    TAR = 250,
    GLO = 251,
    MTL = 252,
    SOR = 253,
    DMU = 254,
    DDQ = 255,
    DMF = 256,
    DIO = 257,
    DOX = 258,
    R12P = 259,
    SDS = 260,
    LMT = 261,
    EOH = 262,
    EEE = 263,
    EGL = 264,
    FLO = 265,
    TRT = 266,
    FCY = 267,
    FRU = 268,
    GBL = 269,
    GPX = 270,
    HTO = 271,
    HTG = 272,
    B7G = 273,
    C10 = 274,
    R16D = 275,
    HEZ = 276,
    IOD = 277,
    IDO = 278,
    ICI = 279,
    ICT = 280,
    TLA = 281,
    LAT = 282,
    LBT = 283,
    LDA = 284,
    MN3 = 285,
    MRY = 286,
    MOH = 287,
    BEQ = 288,
    C15 = 289,
    MG8 = 290,
    POL = 291,
    NO3 = 292,
    JEF = 293,
    P4C = 294,
    CE1 = 295,
    DIA = 296,
    CXE = 297,
    IPH = 298,
    PIN = 299,
    R15P = 300,
    CRY = 301,
    PGR = 302,
    PGQ = 303,
    SPD = 304,
    SPK = 305,
    SPM = 306,
    SUC = 307,
    TBU = 308,
    TMA = 309,
    TEP = 310,
    SCN = 311,
    TRE = 312,
    ETF = 313,
    R144 = 314,
    UMQ = 315,
    URE = 316,
    YT3 = 317,
    ZN2 = 318,
    FE2 = 319,
    R3NI = 320,
    SIA = 321,
    XYP = 322,
    A2G = 323,
    GLA = 324,
    NDG = 325,
    NGA = 326,
    A = 327,
    C = 328,
    G = 329,
    I = 330,
    U = 331,
    N = 332,
    F = 333,
    K = 334,
    DA = 335,
    DC = 336,
    DG = 337,
    DI = 338,
    DT = 339,
    DU = 340,
    DN = 341,
    AG = 342,
    AL = 343,
    BA = 344,
    BR = 345,
    CA = 346,
    CD = 347,
    CL = 348,
    CM = 349,
    CN = 350,
    CO = 351,
    CS = 352,
    CU = 353,
    FE = 354,
    HG = 355,
    LI = 356,
    MG = 357,
    MN = 358,
    NA = 359,
    NI = 360,
    NO = 361,
    PB = 362,
    RB = 363,
    SR = 364,
    Y1 = 365,
    ZN = 366,
    UNKNOWN = 367,
}

impl ResidueCode {
    /// Return the Gemmi table index used by this build.
    ///
    /// This number is an implementation detail, not a stable serialization
    /// format. Use [`Display`](fmt::Display), [`FromStr`], or serde when a
    /// persistent representation is required.
    #[must_use]
    pub const fn as_u16(self) -> u16 {
        self as u16
    }

    /// Return the source-derived metadata for this residue code.
    #[must_use]
    pub const fn info(self) -> ResidueInfo {
        RESIDUE_INFO_TABLE[self as usize]
    }

    /// Return the stable residue name used for display and serialization.
    ///
    /// [`ResidueCode::UNKNOWN`] is the absence of a tabulated residue rather
    /// than the real tabulated residue `UNK`, so it uses the explicit protocol
    /// token `"UNKNOWN"`.
    #[must_use]
    pub const fn name(self) -> &'static str {
        if matches!(self, Self::UNKNOWN) {
            "UNKNOWN"
        } else {
            self.info().name
        }
    }
}

impl fmt::Display for ResidueCode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

/// Error returned when a name is not present in the source residue table.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error("unknown residue code name '{input}'")]
pub struct ResidueCodeParseError {
    input: String,
}

impl ResidueCodeParseError {
    #[must_use]
    pub fn input(&self) -> &str {
        &self.input
    }
}

impl FromStr for ResidueCode {
    type Err = ResidueCodeParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        if input.eq_ignore_ascii_case("UNKNOWN") {
            return Ok(Self::UNKNOWN);
        }
        let code = residue_code_from_name(input);
        if matches!(code, Self::UNKNOWN) {
            Err(ResidueCodeParseError {
                input: input.to_string(),
            })
        } else {
            Ok(code)
        }
    }
}

impl Serialize for ResidueCode {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.name())
    }
}

impl<'de> Deserialize<'de> for ResidueCode {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let name = String::deserialize(deserializer)?;
        name.parse().map_err(serde::de::Error::custom)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ResidueInfoKind {
    Unknown = 0,
    Aa = 1,
    Aad = 2,
    Paa = 3,
    Maa = 4,
    Rna = 5,
    Dna = 6,
    Buf = 7,
    Hoh = 8,
    Pyr = 9,
    Ket = 10,
    Els = 11,
}

impl ResidueInfoKind {
    /// Return the stable Gemmi residue-kind name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Unknown => "UNKNOWN",
            Self::Aa => "AA",
            Self::Aad => "AAD",
            Self::Paa => "PAA",
            Self::Maa => "MAA",
            Self::Rna => "RNA",
            Self::Dna => "DNA",
            Self::Buf => "BUF",
            Self::Hoh => "HOH",
            Self::Pyr => "PYR",
            Self::Ket => "KET",
            Self::Els => "ELS",
        }
    }
}

impl fmt::Display for ResidueInfoKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl Serialize for ResidueInfoKind {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.name())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
#[non_exhaustive]
pub struct ResidueInfo {
    pub code: ResidueCode,
    pub name: &'static str,
    pub kind: ResidueInfoKind,
    pub linking_type: u8,
    pub one_letter_code: char,
    pub hydrogen_count: u8,
    pub weight: f32,
}

/// A residue name together with its table classification.
///
/// The original name remains authoritative so unrecognized component names
/// survive roundtrips. The classified code is private and is always derived
/// by the source-aligned table lookup, preventing contradictory identities.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResidueIdentity {
    name: String,
    code: ResidueCode,
}

impl ResidueIdentity {
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        let name = name.into();
        let code = residue_code_from_name(&name);
        Self { name, code }
    }

    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub const fn code(&self) -> ResidueCode {
        self.code
    }

    #[must_use]
    pub const fn info(&self) -> ResidueInfo {
        self.code.info()
    }

    #[must_use]
    pub const fn is_tabulated(&self) -> bool {
        !matches!(self.code, ResidueCode::UNKNOWN)
    }
}

impl From<&str> for ResidueIdentity {
    fn from(name: &str) -> Self {
        Self::new(name)
    }
}

impl From<String> for ResidueIdentity {
    fn from(name: String) -> Self {
        Self::new(name)
    }
}

impl FromStr for ResidueIdentity {
    type Err = Infallible;

    fn from_str(name: &str) -> Result<Self, Self::Err> {
        Ok(Self::new(name))
    }
}

impl fmt::Display for ResidueIdentity {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.name)
    }
}

impl Serialize for ResidueIdentity {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&self.name)
    }
}

impl<'de> Deserialize<'de> for ResidueIdentity {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        String::deserialize(deserializer).map(Self::new)
    }
}

impl ResidueInfo {
    #[must_use]
    pub const fn found(self) -> bool {
        // Gemmi✔️✔️: bool found() const { return kind != ResidueKind::UNKNOWN; }
        !matches!(self.kind, ResidueInfoKind::Unknown)
    }
    #[must_use]
    pub const fn is_water(self) -> bool {
        // Gemmi✔️✔️: bool is_water() const { return kind == ResidueKind::HOH; }
        matches!(self.kind, ResidueInfoKind::Hoh)
    }
    #[must_use]
    pub const fn is_dna(self) -> bool {
        // Gemmi✔️✔️: bool is_dna() const { return kind == ResidueKind::DNA; }
        matches!(self.kind, ResidueInfoKind::Dna)
    }
    #[must_use]
    pub const fn is_rna(self) -> bool {
        // Gemmi✔️✔️: bool is_rna() const { return kind == ResidueKind::RNA; }
        matches!(self.kind, ResidueInfoKind::Rna)
    }
    #[must_use]
    pub const fn is_nucleic_acid(self) -> bool {
        // Gemmi✔️✔️: bool is_nucleic_acid() const { return is_dna() || is_rna(); }
        self.is_dna() || self.is_rna()
    }
    #[must_use]
    pub const fn is_amino_acid(self) -> bool {
        // Gemmi✔️✔️: return kind == ResidueKind::AA || kind == ResidueKind::AAD ||
        // Gemmi✔️✔️:        kind == ResidueKind::PAA || kind == ResidueKind::MAA;
        matches!(
            self.kind,
            ResidueInfoKind::Aa
                | ResidueInfoKind::Aad
                | ResidueInfoKind::Paa
                | ResidueInfoKind::Maa
        )
    }
    #[must_use]
    pub const fn is_buffer_or_water(self) -> bool {
        // Gemmi✔️✔️: return kind == ResidueKind::HOH || kind == ResidueKind::BUF;
        matches!(self.kind, ResidueInfoKind::Hoh | ResidueInfoKind::Buf)
    }
    #[must_use]
    pub const fn is_standard(self) -> bool {
        // Gemmi✔️✔️: bool is_standard() const { return (one_letter_code & 0x20) == 0; }
        (self.one_letter_code as u32 & 0x20) == 0
    }
    #[must_use]
    pub const fn fasta_code(self) -> char {
        // Gemmi✔️✔️: char fasta_code() const { return is_standard() ? one_letter_code : 'X'; }
        if self.is_standard() {
            self.one_letter_code
        } else {
            'X'
        }
    }

    /// Return the canonical one-letter amino-acid code when Gemmi provides
    /// enough source metadata to identify one.
    ///
    /// Gemmi stores standard codes in uppercase and the corresponding code
    /// for modified amino acids in lowercase. Blank or unmappable entries do
    /// not acquire an inferred parent.
    #[must_use]
    pub fn canonical_one_letter_code(self) -> Option<char> {
        if !self.is_amino_acid() {
            return None;
        }
        let code = self.one_letter_code.to_ascii_uppercase();
        expand_one_letter(code, ResidueInfoKind::Aa).map(|_| code)
    }

    /// Return the standard amino-acid code represented by Gemmi's one-letter
    /// metadata.
    ///
    /// Standard residues return themselves; modified residues return their
    /// table-defined parent. Entries without a mappable one-letter code return
    /// `None` rather than using a name-based guess.
    #[must_use]
    pub fn parent_standard_code(self) -> Option<ResidueCode> {
        let code = self.canonical_one_letter_code()?;
        let name = expand_one_letter(code, ResidueInfoKind::Aa)?;
        name.parse().ok()
    }

    /// Return whether Gemmi classifies this as a non-standard amino acid.
    #[must_use]
    pub const fn is_modified_amino_acid(self) -> bool {
        self.is_amino_acid() && !self.is_standard()
    }
    #[must_use]
    pub const fn is_peptide_linking(self) -> bool {
        // Gemmi✔️✔️: bool is_peptide_linking() const { return (linking_type & 1); }
        (self.linking_type & 1) != 0
    }
    #[must_use]
    pub const fn is_na_linking(self) -> bool {
        // Gemmi✔️✔️: bool is_na_linking() const { return (linking_type & 2); }
        (self.linking_type & 2) != 0
    }
}

pub const UNKNOWN_TABULATED_RESIDUE_IDX: usize = 367;

// BEGIN GEMMI CPP TABLE gemmi::get_residue_info
// Gemmi✔️✔️: static ResidueInfo array[368] = {
pub const RESIDUE_INFO_TABLE: [ResidueInfo; 368] = [
    ResidueInfo {
        code: ResidueCode::ALA,
        name: "ALA",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'A',
        hydrogen_count: 7,
        weight: 89.0932f32,
    },
    ResidueInfo {
        code: ResidueCode::ARG,
        name: "ARG",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'R',
        hydrogen_count: 15,
        weight: 175.209f32,
    },
    ResidueInfo {
        code: ResidueCode::ASN,
        name: "ASN",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'N',
        hydrogen_count: 8,
        weight: 132.118f32,
    },
    ResidueInfo {
        code: ResidueCode::ABA,
        name: "ABA",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'a',
        hydrogen_count: 9,
        weight: 103.120f32,
    },
    ResidueInfo {
        code: ResidueCode::ASP,
        name: "ASP",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'D',
        hydrogen_count: 7,
        weight: 133.103f32,
    },
    ResidueInfo {
        code: ResidueCode::ASX,
        name: "ASX",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'B',
        hydrogen_count: 6,
        weight: 100.096f32,
    },
    ResidueInfo {
        code: ResidueCode::CYS,
        name: "CYS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'C',
        hydrogen_count: 7,
        weight: 121.158f32,
    },
    ResidueInfo {
        code: ResidueCode::CSH,
        name: "CSH",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 's',
        hydrogen_count: 17,
        weight: 283.284f32,
    },
    ResidueInfo {
        code: ResidueCode::GLN,
        name: "GLN",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'Q',
        hydrogen_count: 10,
        weight: 146.144f32,
    },
    ResidueInfo {
        code: ResidueCode::GLU,
        name: "GLU",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'E',
        hydrogen_count: 9,
        weight: 147.129f32,
    },
    ResidueInfo {
        code: ResidueCode::GLX,
        name: "GLX",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'Z',
        hydrogen_count: 8,
        weight: 114.123f32,
    },
    ResidueInfo {
        code: ResidueCode::GLY,
        name: "GLY",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'G',
        hydrogen_count: 5,
        weight: 75.0666f32,
    },
    ResidueInfo {
        code: ResidueCode::HIS,
        name: "HIS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'H',
        hydrogen_count: 10,
        weight: 156.162f32,
    },
    ResidueInfo {
        code: ResidueCode::ILE,
        name: "ILE",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'I',
        hydrogen_count: 13,
        weight: 131.173f32,
    },
    ResidueInfo {
        code: ResidueCode::LEU,
        name: "LEU",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'L',
        hydrogen_count: 13,
        weight: 131.173f32,
    },
    ResidueInfo {
        code: ResidueCode::LYS,
        name: "LYS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'K',
        hydrogen_count: 15,
        weight: 147.196f32,
    },
    ResidueInfo {
        code: ResidueCode::MET,
        name: "MET",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'M',
        hydrogen_count: 11,
        weight: 149.211f32,
    },
    ResidueInfo {
        code: ResidueCode::MSE,
        name: "MSE",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'm',
        hydrogen_count: 11,
        weight: 196.106f32,
    },
    ResidueInfo {
        code: ResidueCode::ORN,
        name: "ORN",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'a',
        hydrogen_count: 12,
        weight: 132.161f32,
    },
    ResidueInfo {
        code: ResidueCode::PHE,
        name: "PHE",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'F',
        hydrogen_count: 11,
        weight: 165.189f32,
    },
    ResidueInfo {
        code: ResidueCode::PRO,
        name: "PRO",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'P',
        hydrogen_count: 9,
        weight: 115.130f32,
    },
    ResidueInfo {
        code: ResidueCode::SER,
        name: "SER",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'S',
        hydrogen_count: 7,
        weight: 105.093f32,
    },
    ResidueInfo {
        code: ResidueCode::THR,
        name: "THR",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'T',
        hydrogen_count: 9,
        weight: 119.119f32,
    },
    ResidueInfo {
        code: ResidueCode::TRP,
        name: "TRP",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'W',
        hydrogen_count: 12,
        weight: 204.225f32,
    },
    ResidueInfo {
        code: ResidueCode::TYR,
        name: "TYR",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'Y',
        hydrogen_count: 11,
        weight: 181.189f32,
    },
    ResidueInfo {
        code: ResidueCode::UNK,
        name: "UNK",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'X',
        hydrogen_count: 9,
        weight: 103.120f32,
    },
    ResidueInfo {
        code: ResidueCode::VAL,
        name: "VAL",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'V',
        hydrogen_count: 11,
        weight: 117.146f32,
    },
    ResidueInfo {
        code: ResidueCode::SEC,
        name: "SEC",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'U',
        hydrogen_count: 7,
        weight: 168.053f32,
    },
    ResidueInfo {
        code: ResidueCode::PYL,
        name: "PYL",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'O',
        hydrogen_count: 21,
        weight: 255.313f32,
    },
    ResidueInfo {
        code: ResidueCode::SEP,
        name: "SEP",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 's',
        hydrogen_count: 8,
        weight: 185.072f32,
    },
    ResidueInfo {
        code: ResidueCode::TPO,
        name: "TPO",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 't',
        hydrogen_count: 10,
        weight: 199.099f32,
    },
    ResidueInfo {
        code: ResidueCode::PCA,
        name: "PCA",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'e',
        hydrogen_count: 7,
        weight: 129.114f32,
    },
    ResidueInfo {
        code: ResidueCode::CSO,
        name: "CSO",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 7,
        weight: 137.158f32,
    },
    ResidueInfo {
        code: ResidueCode::PTR,
        name: "PTR",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'y',
        hydrogen_count: 12,
        weight: 261.168f32,
    },
    ResidueInfo {
        code: ResidueCode::KCX,
        name: "KCX",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'k',
        hydrogen_count: 14,
        weight: 190.197f32,
    },
    ResidueInfo {
        code: ResidueCode::CSD,
        name: "CSD",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 7,
        weight: 153.157f32,
    },
    ResidueInfo {
        code: ResidueCode::LLP,
        name: "LLP",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'k',
        hydrogen_count: 22,
        weight: 375.314f32,
    },
    ResidueInfo {
        code: ResidueCode::CME,
        name: "CME",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 11,
        weight: 197.276f32,
    },
    ResidueInfo {
        code: ResidueCode::MLY,
        name: "MLY",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'k',
        hydrogen_count: 18,
        weight: 174.241f32,
    },
    ResidueInfo {
        code: ResidueCode::DAL,
        name: "DAL",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'a',
        hydrogen_count: 7,
        weight: 89.0932f32,
    },
    ResidueInfo {
        code: ResidueCode::TYS,
        name: "TYS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'y',
        hydrogen_count: 11,
        weight: 261.252f32,
    },
    ResidueInfo {
        code: ResidueCode::OCS,
        name: "OCS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 7,
        weight: 169.156f32,
    },
    ResidueInfo {
        code: ResidueCode::M3L,
        name: "M3L",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'k',
        hydrogen_count: 21,
        weight: 189.275f32,
    },
    ResidueInfo {
        code: ResidueCode::FME,
        name: "FME",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'm',
        hydrogen_count: 11,
        weight: 177.221f32,
    },
    ResidueInfo {
        code: ResidueCode::ALY,
        name: "ALY",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'k',
        hydrogen_count: 16,
        weight: 188.224f32,
    },
    ResidueInfo {
        code: ResidueCode::HYP,
        name: "HYP",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'p',
        hydrogen_count: 9,
        weight: 131.130f32,
    },
    ResidueInfo {
        code: ResidueCode::CAS,
        name: "CAS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 12,
        weight: 225.141f32,
    },
    ResidueInfo {
        code: ResidueCode::CRO,
        name: "CRO",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 't',
        hydrogen_count: 17,
        weight: 319.313f32,
    },
    ResidueInfo {
        code: ResidueCode::CSX,
        name: "CSX",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 7,
        weight: 137.158f32,
    },
    ResidueInfo {
        code: ResidueCode::DPR,
        name: "DPR",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'p',
        hydrogen_count: 9,
        weight: 115.130f32,
    },
    ResidueInfo {
        code: ResidueCode::DGL,
        name: "DGL",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'e',
        hydrogen_count: 9,
        weight: 147.129f32,
    },
    ResidueInfo {
        code: ResidueCode::DVA,
        name: "DVA",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'v',
        hydrogen_count: 11,
        weight: 117.146f32,
    },
    ResidueInfo {
        code: ResidueCode::CSS,
        name: "CSS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 7,
        weight: 153.223f32,
    },
    ResidueInfo {
        code: ResidueCode::DPN,
        name: "DPN",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'f',
        hydrogen_count: 11,
        weight: 165.189f32,
    },
    ResidueInfo {
        code: ResidueCode::DSN,
        name: "DSN",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 's',
        hydrogen_count: 7,
        weight: 105.093f32,
    },
    ResidueInfo {
        code: ResidueCode::DLE,
        name: "DLE",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'l',
        hydrogen_count: 13,
        weight: 131.173f32,
    },
    ResidueInfo {
        code: ResidueCode::HIC,
        name: "HIC",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'h',
        hydrogen_count: 11,
        weight: 169.181f32,
    },
    ResidueInfo {
        code: ResidueCode::NLE,
        name: "NLE",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'l',
        hydrogen_count: 13,
        weight: 131.173f32,
    },
    ResidueInfo {
        code: ResidueCode::MVA,
        name: "MVA",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'v',
        hydrogen_count: 13,
        weight: 131.173f32,
    },
    ResidueInfo {
        code: ResidueCode::MLZ,
        name: "MLZ",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'k',
        hydrogen_count: 16,
        weight: 160.214f32,
    },
    ResidueInfo {
        code: ResidueCode::CR2,
        name: "CR2",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'g',
        hydrogen_count: 13,
        weight: 275.260f32,
    },
    ResidueInfo {
        code: ResidueCode::SAR,
        name: "SAR",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'g',
        hydrogen_count: 7,
        weight: 89.0932f32,
    },
    ResidueInfo {
        code: ResidueCode::DAR,
        name: "DAR",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'r',
        hydrogen_count: 15,
        weight: 175.209f32,
    },
    ResidueInfo {
        code: ResidueCode::DLY,
        name: "DLY",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'k',
        hydrogen_count: 14,
        weight: 146.188f32,
    },
    ResidueInfo {
        code: ResidueCode::YCM,
        name: "YCM",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 10,
        weight: 178.209f32,
    },
    ResidueInfo {
        code: ResidueCode::NRQ,
        name: "NRQ",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'm',
        hydrogen_count: 17,
        weight: 347.389f32,
    },
    ResidueInfo {
        code: ResidueCode::CGU,
        name: "CGU",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'e',
        hydrogen_count: 9,
        weight: 191.139f32,
    },
    ResidueInfo {
        code: ResidueCode::R0TD,
        name: "0TD",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'd',
        hydrogen_count: 9,
        weight: 179.194f32,
    },
    ResidueInfo {
        code: ResidueCode::MLE,
        name: "MLE",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'l',
        hydrogen_count: 15,
        weight: 145.200f32,
    },
    ResidueInfo {
        code: ResidueCode::DAS,
        name: "DAS",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'd',
        hydrogen_count: 7,
        weight: 133.103f32,
    },
    ResidueInfo {
        code: ResidueCode::DTR,
        name: "DTR",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'w',
        hydrogen_count: 12,
        weight: 204.225f32,
    },
    ResidueInfo {
        code: ResidueCode::CXM,
        name: "CXM",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'm',
        hydrogen_count: 11,
        weight: 193.221f32,
    },
    ResidueInfo {
        code: ResidueCode::TPQ,
        name: "TPQ",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'y',
        hydrogen_count: 9,
        weight: 211.171f32,
    },
    ResidueInfo {
        code: ResidueCode::DCY,
        name: "DCY",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 7,
        weight: 121.158f32,
    },
    ResidueInfo {
        code: ResidueCode::DSG,
        name: "DSG",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'n',
        hydrogen_count: 8,
        weight: 132.118f32,
    },
    ResidueInfo {
        code: ResidueCode::DTY,
        name: "DTY",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'y',
        hydrogen_count: 11,
        weight: 181.189f32,
    },
    ResidueInfo {
        code: ResidueCode::DHI,
        name: "DHI",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'h',
        hydrogen_count: 10,
        weight: 156.162f32,
    },
    ResidueInfo {
        code: ResidueCode::MEN,
        name: "MEN",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'n',
        hydrogen_count: 10,
        weight: 146.144f32,
    },
    ResidueInfo {
        code: ResidueCode::DTH,
        name: "DTH",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 't',
        hydrogen_count: 9,
        weight: 119.119f32,
    },
    ResidueInfo {
        code: ResidueCode::SAC,
        name: "SAC",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 's',
        hydrogen_count: 9,
        weight: 147.129f32,
    },
    ResidueInfo {
        code: ResidueCode::DGN,
        name: "DGN",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'q',
        hydrogen_count: 10,
        weight: 146.144f32,
    },
    ResidueInfo {
        code: ResidueCode::AIB,
        name: "AIB",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'a',
        hydrogen_count: 9,
        weight: 103.120f32,
    },
    ResidueInfo {
        code: ResidueCode::SMC,
        name: "SMC",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 9,
        weight: 135.185f32,
    },
    ResidueInfo {
        code: ResidueCode::IAS,
        name: "IAS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'd',
        hydrogen_count: 7,
        weight: 133.103f32,
    },
    ResidueInfo {
        code: ResidueCode::CIR,
        name: "CIR",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'r',
        hydrogen_count: 13,
        weight: 175.186f32,
    },
    ResidueInfo {
        code: ResidueCode::BMT,
        name: "BMT",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 't',
        hydrogen_count: 19,
        weight: 201.263f32,
    },
    ResidueInfo {
        code: ResidueCode::DIL,
        name: "DIL",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'i',
        hydrogen_count: 13,
        weight: 131.173f32,
    },
    ResidueInfo {
        code: ResidueCode::FGA,
        name: "FGA",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'e',
        hydrogen_count: 9,
        weight: 147.129f32,
    },
    ResidueInfo {
        code: ResidueCode::PHI,
        name: "PHI",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'f',
        hydrogen_count: 10,
        weight: 291.086f32,
    },
    ResidueInfo {
        code: ResidueCode::CRQ,
        name: "CRQ",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'q',
        hydrogen_count: 16,
        weight: 344.322f32,
    },
    ResidueInfo {
        code: ResidueCode::SME,
        name: "SME",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'm',
        hydrogen_count: 11,
        weight: 165.211f32,
    },
    ResidueInfo {
        code: ResidueCode::GHP,
        name: "GHP",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'g',
        hydrogen_count: 9,
        weight: 167.162f32,
    },
    ResidueInfo {
        code: ResidueCode::MHO,
        name: "MHO",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'm',
        hydrogen_count: 11,
        weight: 165.211f32,
    },
    ResidueInfo {
        code: ResidueCode::NEP,
        name: "NEP",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'h',
        hydrogen_count: 10,
        weight: 235.134f32,
    },
    ResidueInfo {
        code: ResidueCode::TRQ,
        name: "TRQ",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'w',
        hydrogen_count: 10,
        weight: 234.208f32,
    },
    ResidueInfo {
        code: ResidueCode::TOX,
        name: "TOX",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'w',
        hydrogen_count: 12,
        weight: 236.224f32,
    },
    ResidueInfo {
        code: ResidueCode::ALC,
        name: "ALC",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'a',
        hydrogen_count: 17,
        weight: 171.237f32,
    },
    ResidueInfo {
        code: ResidueCode::R3FG,
        name: "3FG",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: ' ',
        hydrogen_count: 9,
        weight: 183.161f32,
    },
    ResidueInfo {
        code: ResidueCode::SCH,
        name: "SCH",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 9,
        weight: 167.250f32,
    },
    ResidueInfo {
        code: ResidueCode::MDO,
        name: "MDO",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'a',
        hydrogen_count: 11,
        weight: 197.191f32,
    },
    ResidueInfo {
        code: ResidueCode::MAA,
        name: "MAA",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'a',
        hydrogen_count: 9,
        weight: 103.120f32,
    },
    ResidueInfo {
        code: ResidueCode::GYS,
        name: "GYS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 's',
        hydrogen_count: 15,
        weight: 305.286f32,
    },
    ResidueInfo {
        code: ResidueCode::MK8,
        name: "MK8",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'l',
        hydrogen_count: 15,
        weight: 145.200f32,
    },
    ResidueInfo {
        code: ResidueCode::CR8,
        name: "CR8",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'h',
        hydrogen_count: 16,
        weight: 354.340f32,
    },
    ResidueInfo {
        code: ResidueCode::KPI,
        name: "KPI",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'k',
        hydrogen_count: 16,
        weight: 216.234f32,
    },
    ResidueInfo {
        code: ResidueCode::SCY,
        name: "SCY",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 9,
        weight: 163.195f32,
    },
    ResidueInfo {
        code: ResidueCode::DHA,
        name: "DHA",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 's',
        hydrogen_count: 5,
        weight: 87.0773f32,
    },
    ResidueInfo {
        code: ResidueCode::OMY,
        name: "OMY",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'y',
        hydrogen_count: 10,
        weight: 231.633f32,
    },
    ResidueInfo {
        code: ResidueCode::CAF,
        name: "CAF",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 12,
        weight: 241.140f32,
    },
    ResidueInfo {
        code: ResidueCode::R0AF,
        name: "0AF",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'w',
        hydrogen_count: 12,
        weight: 220.225f32,
    },
    ResidueInfo {
        code: ResidueCode::SNN,
        name: "SNN",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'n',
        hydrogen_count: 6,
        weight: 114.103f32,
    },
    ResidueInfo {
        code: ResidueCode::MHS,
        name: "MHS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'h',
        hydrogen_count: 11,
        weight: 169.181f32,
    },
    ResidueInfo {
        code: ResidueCode::MLU,
        name: "MLU",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: ' ',
        hydrogen_count: 15,
        weight: 145.200f32,
    },
    ResidueInfo {
        code: ResidueCode::SNC,
        name: "SNC",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'c',
        hydrogen_count: 6,
        weight: 150.156f32,
    },
    ResidueInfo {
        code: ResidueCode::PHD,
        name: "PHD",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'd',
        hydrogen_count: 8,
        weight: 213.083f32,
    },
    ResidueInfo {
        code: ResidueCode::B3E,
        name: "B3E",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'e',
        hydrogen_count: 11,
        weight: 161.156f32,
    },
    ResidueInfo {
        code: ResidueCode::MEA,
        name: "MEA",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'f',
        hydrogen_count: 13,
        weight: 179.216f32,
    },
    ResidueInfo {
        code: ResidueCode::MED,
        name: "MED",
        kind: ResidueInfoKind::Aad,
        linking_type: 1,
        one_letter_code: 'm',
        hydrogen_count: 11,
        weight: 149.211f32,
    },
    ResidueInfo {
        code: ResidueCode::OAS,
        name: "OAS",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 's',
        hydrogen_count: 9,
        weight: 147.129f32,
    },
    ResidueInfo {
        code: ResidueCode::GL3,
        name: "GL3",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'g',
        hydrogen_count: 5,
        weight: 91.1322f32,
    },
    ResidueInfo {
        code: ResidueCode::FVA,
        name: "FVA",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'v',
        hydrogen_count: 11,
        weight: 145.156f32,
    },
    ResidueInfo {
        code: ResidueCode::PHL,
        name: "PHL",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'f',
        hydrogen_count: 13,
        weight: 151.206f32,
    },
    ResidueInfo {
        code: ResidueCode::CRF,
        name: "CRF",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 't',
        hydrogen_count: 18,
        weight: 342.349f32,
    },
    ResidueInfo {
        code: ResidueCode::OMZ,
        name: "OMZ",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 231.633f32,
    },
    ResidueInfo {
        code: ResidueCode::BFD,
        name: "BFD",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'd',
        hydrogen_count: 6,
        weight: 198.102f32,
    },
    ResidueInfo {
        code: ResidueCode::MEQ,
        name: "MEQ",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'q',
        hydrogen_count: 12,
        weight: 160.171f32,
    },
    ResidueInfo {
        code: ResidueCode::DAB,
        name: "DAB",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'a',
        hydrogen_count: 10,
        weight: 118.134f32,
    },
    ResidueInfo {
        code: ResidueCode::AGM,
        name: "AGM",
        kind: ResidueInfoKind::Aa,
        linking_type: 1,
        one_letter_code: 'r',
        hydrogen_count: 17,
        weight: 189.235f32,
    },
    ResidueInfo {
        code: ResidueCode::PSU,
        name: "PSU",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'u',
        hydrogen_count: 13,
        weight: 324.181f32,
    },
    ResidueInfo {
        code: ResidueCode::R5MU,
        name: "5MU",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'u',
        hydrogen_count: 15,
        weight: 338.208f32,
    },
    ResidueInfo {
        code: ResidueCode::R7MG,
        name: "7MG",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'g',
        hydrogen_count: 18,
        weight: 379.263f32,
    },
    ResidueInfo {
        code: ResidueCode::OMG,
        name: "OMG",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'g',
        hydrogen_count: 16,
        weight: 377.247f32,
    },
    ResidueInfo {
        code: ResidueCode::UR3,
        name: "UR3",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'u',
        hydrogen_count: 15,
        weight: 338.208f32,
    },
    ResidueInfo {
        code: ResidueCode::OMC,
        name: "OMC",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'c',
        hydrogen_count: 16,
        weight: 337.223f32,
    },
    ResidueInfo {
        code: ResidueCode::R2MG,
        name: "2MG",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'g',
        hydrogen_count: 16,
        weight: 377.247f32,
    },
    ResidueInfo {
        code: ResidueCode::H2U,
        name: "H2U",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'u',
        hydrogen_count: 15,
        weight: 326.197f32,
    },
    ResidueInfo {
        code: ResidueCode::R4SU,
        name: "4SU",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'u',
        hydrogen_count: 13,
        weight: 340.247f32,
    },
    ResidueInfo {
        code: ResidueCode::OMU,
        name: "OMU",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'u',
        hydrogen_count: 15,
        weight: 338.208f32,
    },
    ResidueInfo {
        code: ResidueCode::R4OC,
        name: "4OC",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'c',
        hydrogen_count: 18,
        weight: 351.250f32,
    },
    ResidueInfo {
        code: ResidueCode::MA6,
        name: "MA6",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'a',
        hydrogen_count: 18,
        weight: 375.274f32,
    },
    ResidueInfo {
        code: ResidueCode::M2G,
        name: "M2G",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'g',
        hydrogen_count: 18,
        weight: 391.274f32,
    },
    ResidueInfo {
        code: ResidueCode::R1MA,
        name: "1MA",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'a',
        hydrogen_count: 16,
        weight: 361.248f32,
    },
    ResidueInfo {
        code: ResidueCode::R6MZ,
        name: "6MZ",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'a',
        hydrogen_count: 16,
        weight: 361.248f32,
    },
    ResidueInfo {
        code: ResidueCode::CCC,
        name: "CCC",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'c',
        hydrogen_count: 13,
        weight: 385.161f32,
    },
    ResidueInfo {
        code: ResidueCode::R2MA,
        name: "2MA",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'a',
        hydrogen_count: 16,
        weight: 361.248f32,
    },
    ResidueInfo {
        code: ResidueCode::R1MG,
        name: "1MG",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'g',
        hydrogen_count: 16,
        weight: 377.247f32,
    },
    ResidueInfo {
        code: ResidueCode::R5BU,
        name: "5BU",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'u',
        hydrogen_count: 12,
        weight: 403.077f32,
    },
    ResidueInfo {
        code: ResidueCode::MIA,
        name: "MIA",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'a',
        hydrogen_count: 24,
        weight: 461.430f32,
    },
    ResidueInfo {
        code: ResidueCode::DOC,
        name: "DOC",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'c',
        hydrogen_count: 14,
        weight: 291.198f32,
    },
    ResidueInfo {
        code: ResidueCode::R8OG,
        name: "8OG",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'g',
        hydrogen_count: 14,
        weight: 363.221f32,
    },
    ResidueInfo {
        code: ResidueCode::R5CM,
        name: "5CM",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'c',
        hydrogen_count: 16,
        weight: 321.224f32,
    },
    ResidueInfo {
        code: ResidueCode::R3DR,
        name: "3DR",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: ' ',
        hydrogen_count: 11,
        weight: 198.111f32,
    },
    ResidueInfo {
        code: ResidueCode::BRU,
        name: "BRU",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'u',
        hydrogen_count: 12,
        weight: 387.078f32,
    },
    ResidueInfo {
        code: ResidueCode::CBR,
        name: "CBR",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'c',
        hydrogen_count: 13,
        weight: 386.093f32,
    },
    ResidueInfo {
        code: ResidueCode::HOH,
        name: "HOH",
        kind: ResidueInfoKind::Hoh,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 2,
        weight: 18.0153f32,
    },
    ResidueInfo {
        code: ResidueCode::DOD,
        name: "DOD",
        kind: ResidueInfoKind::Hoh,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 2,
        weight: 20.0276f32,
    },
    ResidueInfo {
        code: ResidueCode::HEM,
        name: "HEM",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 32,
        weight: 616.487f32,
    },
    ResidueInfo {
        code: ResidueCode::SO4,
        name: "SO4",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 96.0626f32,
    },
    ResidueInfo {
        code: ResidueCode::GOL,
        name: "GOL",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 92.0938f32,
    },
    ResidueInfo {
        code: ResidueCode::EDO,
        name: "EDO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 62.0678f32,
    },
    ResidueInfo {
        code: ResidueCode::NAG,
        name: "NAG",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 15,
        weight: 221.208f32,
    },
    ResidueInfo {
        code: ResidueCode::PO4,
        name: "PO4",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 94.9714f32,
    },
    ResidueInfo {
        code: ResidueCode::ACT,
        name: "ACT",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 3,
        weight: 59.0440f32,
    },
    ResidueInfo {
        code: ResidueCode::PEG,
        name: "PEG",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 106.120f32,
    },
    ResidueInfo {
        code: ResidueCode::MAN,
        name: "MAN",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 180.156f32,
    },
    ResidueInfo {
        code: ResidueCode::FAD,
        name: "FAD",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 33,
        weight: 785.550f32,
    },
    ResidueInfo {
        code: ResidueCode::BMA,
        name: "BMA",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 180.156f32,
    },
    ResidueInfo {
        code: ResidueCode::ADP,
        name: "ADP",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 15,
        weight: 427.201f32,
    },
    ResidueInfo {
        code: ResidueCode::DMS,
        name: "DMS",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 78.1334f32,
    },
    ResidueInfo {
        code: ResidueCode::ACE,
        name: "ACE",
        kind: ResidueInfoKind::Els,
        linking_type: 1,
        one_letter_code: ' ',
        hydrogen_count: 4,
        weight: 44.0526f32,
    },
    ResidueInfo {
        code: ResidueCode::NH2,
        name: "NH2",
        kind: ResidueInfoKind::Els,
        linking_type: 1,
        one_letter_code: ' ',
        hydrogen_count: 2,
        weight: 16.0226f32,
    },
    ResidueInfo {
        code: ResidueCode::MPD,
        name: "MPD",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 14,
        weight: 118.174f32,
    },
    ResidueInfo {
        code: ResidueCode::MES,
        name: "MES",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 13,
        weight: 195.237f32,
    },
    ResidueInfo {
        code: ResidueCode::NAD,
        name: "NAD",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 27,
        weight: 663.425f32,
    },
    ResidueInfo {
        code: ResidueCode::NAP,
        name: "NAP",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 28,
        weight: 743.405f32,
    },
    ResidueInfo {
        code: ResidueCode::TRS,
        name: "TRS",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 122.143f32,
    },
    ResidueInfo {
        code: ResidueCode::ATP,
        name: "ATP",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 16,
        weight: 507.181f32,
    },
    ResidueInfo {
        code: ResidueCode::PG4,
        name: "PG4",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 18,
        weight: 194.226f32,
    },
    ResidueInfo {
        code: ResidueCode::GDP,
        name: "GDP",
        kind: ResidueInfoKind::Els,
        linking_type: 2,
        one_letter_code: 'g',
        hydrogen_count: 15,
        weight: 443.201f32,
    },
    ResidueInfo {
        code: ResidueCode::FUC,
        name: "FUC",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 164.156f32,
    },
    ResidueInfo {
        code: ResidueCode::FMT,
        name: "FMT",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 2,
        weight: 46.0254f32,
    },
    ResidueInfo {
        code: ResidueCode::GAL,
        name: "GAL",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 180.156f32,
    },
    ResidueInfo {
        code: ResidueCode::PGE,
        name: "PGE",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 14,
        weight: 150.173f32,
    },
    ResidueInfo {
        code: ResidueCode::FMN,
        name: "FMN",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 21,
        weight: 456.344f32,
    },
    ResidueInfo {
        code: ResidueCode::PLP,
        name: "PLP",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 247.142f32,
    },
    ResidueInfo {
        code: ResidueCode::EPE,
        name: "EPE",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 18,
        weight: 238.305f32,
    },
    ResidueInfo {
        code: ResidueCode::SF4,
        name: "SF4",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 351.640f32,
    },
    ResidueInfo {
        code: ResidueCode::BME,
        name: "BME",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 78.1334f32,
    },
    ResidueInfo {
        code: ResidueCode::CIT,
        name: "CIT",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 192.124f32,
    },
    ResidueInfo {
        code: ResidueCode::BE7,
        name: "BE7",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 5,
        weight: 357.156f32,
    },
    ResidueInfo {
        code: ResidueCode::MRD,
        name: "MRD",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 14,
        weight: 118.174f32,
    },
    ResidueInfo {
        code: ResidueCode::MHA,
        name: "MHA",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 190.154f32,
    },
    ResidueInfo {
        code: ResidueCode::BU3,
        name: "BU3",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 90.1210f32,
    },
    ResidueInfo {
        code: ResidueCode::PGO,
        name: "PGO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 76.0944f32,
    },
    ResidueInfo {
        code: ResidueCode::BU2,
        name: "BU2",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 90.1210f32,
    },
    ResidueInfo {
        code: ResidueCode::PDO,
        name: "PDO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 76.0944f32,
    },
    ResidueInfo {
        code: ResidueCode::BU1,
        name: "BU1",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 90.1210f32,
    },
    ResidueInfo {
        code: ResidueCode::PG6,
        name: "PG6",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 26,
        weight: 266.331f32,
    },
    ResidueInfo {
        code: ResidueCode::R1BO,
        name: "1BO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 74.1216f32,
    },
    ResidueInfo {
        code: ResidueCode::PE7,
        name: "PE7",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 30,
        weight: 342.449f32,
    },
    ResidueInfo {
        code: ResidueCode::PG5,
        name: "PG5",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 18,
        weight: 178.226f32,
    },
    ResidueInfo {
        code: ResidueCode::TFP,
        name: "TFP",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 24,
        weight: 407.496f32,
    },
    ResidueInfo {
        code: ResidueCode::DHD,
        name: "DHD",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 4,
        weight: 160.082f32,
    },
    ResidueInfo {
        code: ResidueCode::PEU,
        name: "PEU",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 112,
        weight: 1221.46f32,
    },
    ResidueInfo {
        code: ResidueCode::TAU,
        name: "TAU",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 7,
        weight: 125.147f32,
    },
    ResidueInfo {
        code: ResidueCode::SBT,
        name: "SBT",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 74.1216f32,
    },
    ResidueInfo {
        code: ResidueCode::SAL,
        name: "SAL",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 138.121f32,
    },
    ResidueInfo {
        code: ResidueCode::IOH,
        name: "IOH",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 60.0950f32,
    },
    ResidueInfo {
        code: ResidueCode::IPA,
        name: "IPA",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 60.0950f32,
    },
    ResidueInfo {
        code: ResidueCode::PIG,
        name: "PIG",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 14,
        weight: 150.173f32,
    },
    ResidueInfo {
        code: ResidueCode::B3P,
        name: "B3P",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 26,
        weight: 282.334f32,
    },
    ResidueInfo {
        code: ResidueCode::BTB,
        name: "BTB",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 19,
        weight: 209.240f32,
    },
    ResidueInfo {
        code: ResidueCode::NHE,
        name: "NHE",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 17,
        weight: 207.290f32,
    },
    ResidueInfo {
        code: ResidueCode::C8E,
        name: "C8E",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 34,
        weight: 306.438f32,
    },
    ResidueInfo {
        code: ResidueCode::OTE,
        name: "OTE",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 34,
        weight: 306.438f32,
    },
    ResidueInfo {
        code: ResidueCode::PE4,
        name: "PE4",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 34,
        weight: 354.436f32,
    },
    ResidueInfo {
        code: ResidueCode::XPE,
        name: "XPE",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 42,
        weight: 458.541f32,
    },
    ResidueInfo {
        code: ResidueCode::PE8,
        name: "PE8",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 34,
        weight: 370.436f32,
    },
    ResidueInfo {
        code: ResidueCode::P33,
        name: "P33",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 30,
        weight: 326.383f32,
    },
    ResidueInfo {
        code: ResidueCode::N8E,
        name: "N8E",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 38,
        weight: 350.491f32,
    },
    ResidueInfo {
        code: ResidueCode::R2OS,
        name: "2OS",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 36,
        weight: 468.493f32,
    },
    ResidueInfo {
        code: ResidueCode::R1PS,
        name: "1PS",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 11,
        weight: 201.243f32,
    },
    ResidueInfo {
        code: ResidueCode::CPS,
        name: "CPS",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 58,
        weight: 614.877f32,
    },
    ResidueInfo {
        code: ResidueCode::DMX,
        name: "DMX",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 19,
        weight: 257.349f32,
    },
    ResidueInfo {
        code: ResidueCode::MPO,
        name: "MPO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 15,
        weight: 209.263f32,
    },
    ResidueInfo {
        code: ResidueCode::GCD,
        name: "GCD",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 176.124f32,
    },
    ResidueInfo {
        code: ResidueCode::DXG,
        name: "DXG",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 192.124f32,
    },
    ResidueInfo {
        code: ResidueCode::CM5,
        name: "CM5",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 42,
        weight: 494.573f32,
    },
    ResidueInfo {
        code: ResidueCode::ACA,
        name: "ACA",
        kind: ResidueInfoKind::Buf,
        linking_type: 1,
        one_letter_code: ' ',
        hydrogen_count: 13,
        weight: 131.173f32,
    },
    ResidueInfo {
        code: ResidueCode::ACN,
        name: "ACN",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 58.0791f32,
    },
    ResidueInfo {
        code: ResidueCode::CCN,
        name: "CCN",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 3,
        weight: 41.0519f32,
    },
    ResidueInfo {
        code: ResidueCode::GLC,
        name: "GLC",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 180.156f32,
    },
    ResidueInfo {
        code: ResidueCode::DR6,
        name: "DR6",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 142,
        weight: 1527.90f32,
    },
    ResidueInfo {
        code: ResidueCode::NH4,
        name: "NH4",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 4,
        weight: 18.0385f32,
    },
    ResidueInfo {
        code: ResidueCode::AZI,
        name: "AZI",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 42.0201f32,
    },
    ResidueInfo {
        code: ResidueCode::BNG,
        name: "BNG",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 30,
        weight: 306.395f32,
    },
    ResidueInfo {
        code: ResidueCode::BOG,
        name: "BOG",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 28,
        weight: 292.369f32,
    },
    ResidueInfo {
        code: ResidueCode::BGC,
        name: "BGC",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 180.156f32,
    },
    ResidueInfo {
        code: ResidueCode::BCN,
        name: "BCN",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 13,
        weight: 163.172f32,
    },
    ResidueInfo {
        code: ResidueCode::BRO,
        name: "BRO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 79.9040f32,
    },
    ResidueInfo {
        code: ResidueCode::CAC,
        name: "CAC",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 136.989f32,
    },
    ResidueInfo {
        code: ResidueCode::CBX,
        name: "CBX",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 2,
        weight: 46.0254f32,
    },
    ResidueInfo {
        code: ResidueCode::ACY,
        name: "ACY",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 4,
        weight: 60.0520f32,
    },
    ResidueInfo {
        code: ResidueCode::CBM,
        name: "CBM",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 4,
        weight: 60.0520f32,
    },
    ResidueInfo {
        code: ResidueCode::CLO,
        name: "CLO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 35.4530f32,
    },
    ResidueInfo {
        code: ResidueCode::R3CO,
        name: "3CO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 58.9332f32,
    },
    ResidueInfo {
        code: ResidueCode::NCO,
        name: "NCO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 18,
        weight: 161.116f32,
    },
    ResidueInfo {
        code: ResidueCode::CU1,
        name: "CU1",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 63.5460f32,
    },
    ResidueInfo {
        code: ResidueCode::CYN,
        name: "CYN",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 26.0174f32,
    },
    ResidueInfo {
        code: ResidueCode::MA4,
        name: "MA4",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 44,
        weight: 508.600f32,
    },
    ResidueInfo {
        code: ResidueCode::TAR,
        name: "TAR",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 150.087f32,
    },
    ResidueInfo {
        code: ResidueCode::GLO,
        name: "GLO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 180.156f32,
    },
    ResidueInfo {
        code: ResidueCode::MTL,
        name: "MTL",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 14,
        weight: 182.172f32,
    },
    ResidueInfo {
        code: ResidueCode::SOR,
        name: "SOR",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 14,
        weight: 182.172f32,
    },
    ResidueInfo {
        code: ResidueCode::DMU,
        name: "DMU",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 42,
        weight: 482.562f32,
    },
    ResidueInfo {
        code: ResidueCode::DDQ,
        name: "DDQ",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 27,
        weight: 201.349f32,
    },
    ResidueInfo {
        code: ResidueCode::DMF,
        name: "DMF",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 7,
        weight: 73.0938f32,
    },
    ResidueInfo {
        code: ResidueCode::DIO,
        name: "DIO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 88.1051f32,
    },
    ResidueInfo {
        code: ResidueCode::DOX,
        name: "DOX",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 88.1051f32,
    },
    ResidueInfo {
        code: ResidueCode::R12P,
        name: "12P",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 50,
        weight: 546.646f32,
    },
    ResidueInfo {
        code: ResidueCode::SDS,
        name: "SDS",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 26,
        weight: 266.397f32,
    },
    ResidueInfo {
        code: ResidueCode::LMT,
        name: "LMT",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 46,
        weight: 510.615f32,
    },
    ResidueInfo {
        code: ResidueCode::EOH,
        name: "EOH",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 46.0684f32,
    },
    ResidueInfo {
        code: ResidueCode::EEE,
        name: "EEE",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 88.1051f32,
    },
    ResidueInfo {
        code: ResidueCode::EGL,
        name: "EGL",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 62.0678f32,
    },
    ResidueInfo {
        code: ResidueCode::FLO,
        name: "FLO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 18.9984f32,
    },
    ResidueInfo {
        code: ResidueCode::TRT,
        name: "TRT",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 36,
        weight: 352.508f32,
    },
    ResidueInfo {
        code: ResidueCode::FCY,
        name: "FCY",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 7,
        weight: 121.158f32,
    },
    ResidueInfo {
        code: ResidueCode::FRU,
        name: "FRU",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 180.156f32,
    },
    ResidueInfo {
        code: ResidueCode::GBL,
        name: "GBL",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 86.0892f32,
    },
    ResidueInfo {
        code: ResidueCode::GPX,
        name: "GPX",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 14,
        weight: 505.165f32,
    },
    ResidueInfo {
        code: ResidueCode::HTO,
        name: "HTO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 16,
        weight: 148.200f32,
    },
    ResidueInfo {
        code: ResidueCode::HTG,
        name: "HTG",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 26,
        weight: 294.408f32,
    },
    ResidueInfo {
        code: ResidueCode::B7G,
        name: "B7G",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 26,
        weight: 278.342f32,
    },
    ResidueInfo {
        code: ResidueCode::C10,
        name: "C10",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 46,
        weight: 422.596f32,
    },
    ResidueInfo {
        code: ResidueCode::R16D,
        name: "16D",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 16,
        weight: 116.205f32,
    },
    ResidueInfo {
        code: ResidueCode::HEZ,
        name: "HEZ",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 14,
        weight: 118.174f32,
    },
    ResidueInfo {
        code: ResidueCode::IOD,
        name: "IOD",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 126.904f32,
    },
    ResidueInfo {
        code: ResidueCode::IDO,
        name: "IDO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 126.904f32,
    },
    ResidueInfo {
        code: ResidueCode::ICI,
        name: "ICI",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 192.124f32,
    },
    ResidueInfo {
        code: ResidueCode::ICT,
        name: "ICT",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 192.124f32,
    },
    ResidueInfo {
        code: ResidueCode::TLA,
        name: "TLA",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 150.087f32,
    },
    ResidueInfo {
        code: ResidueCode::LAT,
        name: "LAT",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 22,
        weight: 342.296f32,
    },
    ResidueInfo {
        code: ResidueCode::LBT,
        name: "LBT",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 22,
        weight: 342.296f32,
    },
    ResidueInfo {
        code: ResidueCode::LDA,
        name: "LDA",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 31,
        weight: 229.402f32,
    },
    ResidueInfo {
        code: ResidueCode::MN3,
        name: "MN3",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 54.9380f32,
    },
    ResidueInfo {
        code: ResidueCode::MRY,
        name: "MRY",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 122.120f32,
    },
    ResidueInfo {
        code: ResidueCode::MOH,
        name: "MOH",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 4,
        weight: 32.0419f32,
    },
    ResidueInfo {
        code: ResidueCode::BEQ,
        name: "BEQ",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 38,
        weight: 342.517f32,
    },
    ResidueInfo {
        code: ResidueCode::C15,
        name: "C15",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 38,
        weight: 336.554f32,
    },
    ResidueInfo {
        code: ResidueCode::MG8,
        name: "MG8",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 31,
        weight: 321.410f32,
    },
    ResidueInfo {
        code: ResidueCode::POL,
        name: "POL",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 60.0950f32,
    },
    ResidueInfo {
        code: ResidueCode::NO3,
        name: "NO3",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 62.0049f32,
    },
    ResidueInfo {
        code: ResidueCode::JEF,
        name: "JEF",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 63,
        weight: 597.822f32,
    },
    ResidueInfo {
        code: ResidueCode::P4C,
        name: "P4C",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 28,
        weight: 324.367f32,
    },
    ResidueInfo {
        code: ResidueCode::CE1,
        name: "CE1",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 58,
        weight: 538.755f32,
    },
    ResidueInfo {
        code: ResidueCode::DIA,
        name: "DIA",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 20,
        weight: 144.258f32,
    },
    ResidueInfo {
        code: ResidueCode::CXE,
        name: "CXE",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 42,
        weight: 378.544f32,
    },
    ResidueInfo {
        code: ResidueCode::IPH,
        name: "IPH",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 6,
        weight: 94.1112f32,
    },
    ResidueInfo {
        code: ResidueCode::PIN,
        name: "PIN",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 18,
        weight: 302.368f32,
    },
    ResidueInfo {
        code: ResidueCode::R15P,
        name: "15P",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 140,
        weight: 1529.83f32,
    },
    ResidueInfo {
        code: ResidueCode::CRY,
        name: "CRY",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 92.0938f32,
    },
    ResidueInfo {
        code: ResidueCode::PGR,
        name: "PGR",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 76.0944f32,
    },
    ResidueInfo {
        code: ResidueCode::PGQ,
        name: "PGQ",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 76.0944f32,
    },
    ResidueInfo {
        code: ResidueCode::SPD,
        name: "SPD",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 19,
        weight: 145.246f32,
    },
    ResidueInfo {
        code: ResidueCode::SPK,
        name: "SPK",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 30,
        weight: 206.372f32,
    },
    ResidueInfo {
        code: ResidueCode::SPM,
        name: "SPM",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 26,
        weight: 202.340f32,
    },
    ResidueInfo {
        code: ResidueCode::SUC,
        name: "SUC",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 22,
        weight: 342.296f32,
    },
    ResidueInfo {
        code: ResidueCode::TBU,
        name: "TBU",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 10,
        weight: 74.1216f32,
    },
    ResidueInfo {
        code: ResidueCode::TMA,
        name: "TMA",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 74.1448f32,
    },
    ResidueInfo {
        code: ResidueCode::TEP,
        name: "TEP",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 8,
        weight: 180.164f32,
    },
    ResidueInfo {
        code: ResidueCode::SCN,
        name: "SCN",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 58.0824f32,
    },
    ResidueInfo {
        code: ResidueCode::TRE,
        name: "TRE",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 22,
        weight: 342.296f32,
    },
    ResidueInfo {
        code: ResidueCode::ETF,
        name: "ETF",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 3,
        weight: 100.040f32,
    },
    ResidueInfo {
        code: ResidueCode::R144,
        name: "144",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 12,
        weight: 122.143f32,
    },
    ResidueInfo {
        code: ResidueCode::UMQ,
        name: "UMQ",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 44,
        weight: 496.589f32,
    },
    ResidueInfo {
        code: ResidueCode::URE,
        name: "URE",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 4,
        weight: 60.0553f32,
    },
    ResidueInfo {
        code: ResidueCode::YT3,
        name: "YT3",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 88.9059f32,
    },
    ResidueInfo {
        code: ResidueCode::ZN2,
        name: "ZN2",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 65.3800f32,
    },
    ResidueInfo {
        code: ResidueCode::FE2,
        name: "FE2",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 55.8450f32,
    },
    ResidueInfo {
        code: ResidueCode::R3NI,
        name: "3NI",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 58.6934f32,
    },
    ResidueInfo {
        code: ResidueCode::SIA,
        name: "SIA",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 0.0f32,
    },
    ResidueInfo {
        code: ResidueCode::XYP,
        name: "XYP",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 0.0f32,
    },
    ResidueInfo {
        code: ResidueCode::A2G,
        name: "A2G",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 0.0f32,
    },
    ResidueInfo {
        code: ResidueCode::GLA,
        name: "GLA",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 0.0f32,
    },
    ResidueInfo {
        code: ResidueCode::NDG,
        name: "NDG",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 0.0f32,
    },
    ResidueInfo {
        code: ResidueCode::NGA,
        name: "NGA",
        kind: ResidueInfoKind::Pyr,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 0.0f32,
    },
    ResidueInfo {
        code: ResidueCode::A,
        name: "A",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'A',
        hydrogen_count: 14,
        weight: 347.221f32,
    },
    ResidueInfo {
        code: ResidueCode::C,
        name: "C",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'C',
        hydrogen_count: 14,
        weight: 323.197f32,
    },
    ResidueInfo {
        code: ResidueCode::G,
        name: "G",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'G',
        hydrogen_count: 14,
        weight: 363.221f32,
    },
    ResidueInfo {
        code: ResidueCode::I,
        name: "I",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'I',
        hydrogen_count: 13,
        weight: 348.206f32,
    },
    ResidueInfo {
        code: ResidueCode::U,
        name: "U",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'U',
        hydrogen_count: 13,
        weight: 324.181f32,
    },
    ResidueInfo {
        code: ResidueCode::N,
        name: "N",
        kind: ResidueInfoKind::Rna,
        linking_type: 2,
        one_letter_code: 'N',
        hydrogen_count: 11,
        weight: 214.11f32,
    },
    ResidueInfo {
        code: ResidueCode::F,
        name: "F",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 18.9984f32,
    },
    ResidueInfo {
        code: ResidueCode::K,
        name: "K",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 39.0983f32,
    },
    ResidueInfo {
        code: ResidueCode::DA,
        name: "DA",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'A',
        hydrogen_count: 14,
        weight: 331.222f32,
    },
    ResidueInfo {
        code: ResidueCode::DC,
        name: "DC",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'C',
        hydrogen_count: 14,
        weight: 307.197f32,
    },
    ResidueInfo {
        code: ResidueCode::DG,
        name: "DG",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'G',
        hydrogen_count: 14,
        weight: 347.221f32,
    },
    ResidueInfo {
        code: ResidueCode::DI,
        name: "DI",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'I',
        hydrogen_count: 13,
        weight: 332.207f32,
    },
    ResidueInfo {
        code: ResidueCode::DT,
        name: "DT",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'T',
        hydrogen_count: 15,
        weight: 322.208f32,
    },
    ResidueInfo {
        code: ResidueCode::DU,
        name: "DU",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'U',
        hydrogen_count: 13,
        weight: 308.182f32,
    },
    ResidueInfo {
        code: ResidueCode::DN,
        name: "DN",
        kind: ResidueInfoKind::Dna,
        linking_type: 2,
        one_letter_code: 'N',
        hydrogen_count: 14,
        weight: 198.111f32,
    },
    ResidueInfo {
        code: ResidueCode::AG,
        name: "AG",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 107.868f32,
    },
    ResidueInfo {
        code: ResidueCode::AL,
        name: "AL",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 26.9815f32,
    },
    ResidueInfo {
        code: ResidueCode::BA,
        name: "BA",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 137.327f32,
    },
    ResidueInfo {
        code: ResidueCode::BR,
        name: "BR",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 79.9040f32,
    },
    ResidueInfo {
        code: ResidueCode::CA,
        name: "CA",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 40.0780f32,
    },
    ResidueInfo {
        code: ResidueCode::CD,
        name: "CD",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 112.411f32,
    },
    ResidueInfo {
        code: ResidueCode::CL,
        name: "CL",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 35.4530f32,
    },
    ResidueInfo {
        code: ResidueCode::CM,
        name: "CM",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 4,
        weight: 60.0520f32,
    },
    ResidueInfo {
        code: ResidueCode::CN,
        name: "CN",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 27.0253f32,
    },
    ResidueInfo {
        code: ResidueCode::CO,
        name: "CO",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 58.9332f32,
    },
    ResidueInfo {
        code: ResidueCode::CS,
        name: "CS",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 132.905f32,
    },
    ResidueInfo {
        code: ResidueCode::CU,
        name: "CU",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 63.5460f32,
    },
    ResidueInfo {
        code: ResidueCode::FE,
        name: "FE",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 55.8450f32,
    },
    ResidueInfo {
        code: ResidueCode::HG,
        name: "HG",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 200.590f32,
    },
    ResidueInfo {
        code: ResidueCode::LI,
        name: "LI",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 6.94100f32,
    },
    ResidueInfo {
        code: ResidueCode::MG,
        name: "MG",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 24.3050f32,
    },
    ResidueInfo {
        code: ResidueCode::MN,
        name: "MN",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 54.9380f32,
    },
    ResidueInfo {
        code: ResidueCode::NA,
        name: "NA",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 22.9898f32,
    },
    ResidueInfo {
        code: ResidueCode::NI,
        name: "NI",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 58.6934f32,
    },
    ResidueInfo {
        code: ResidueCode::NO,
        name: "NO",
        kind: ResidueInfoKind::Els,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 30.0061f32,
    },
    ResidueInfo {
        code: ResidueCode::PB,
        name: "PB",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 207.200f32,
    },
    ResidueInfo {
        code: ResidueCode::RB,
        name: "RB",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 85.4678f32,
    },
    ResidueInfo {
        code: ResidueCode::SR,
        name: "SR",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 87.6200f32,
    },
    ResidueInfo {
        code: ResidueCode::Y1,
        name: "Y1",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 88.9059f32,
    },
    ResidueInfo {
        code: ResidueCode::ZN,
        name: "ZN",
        kind: ResidueInfoKind::Buf,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 65.3800f32,
    },
    ResidueInfo {
        code: ResidueCode::UNKNOWN,
        name: "",
        kind: ResidueInfoKind::Unknown,
        linking_type: 0,
        one_letter_code: ' ',
        hydrogen_count: 0,
        weight: 0.0f32,
    },
];
// Gemmi✔️✔️: };
// END GEMMI CPP TABLE gemmi::get_residue_info

#[must_use]
pub fn get_residue_info(idx: usize) -> ResidueInfo {
    // Gemmi✔️✔️: ResidueInfo& get_residue_info(size_t idx) {
    // Gemmi✔️✔️:   static ResidueInfo array[368] = {
    // Gemmi✔️✔️:   return array[idx];
    RESIDUE_INFO_TABLE[idx]
}

#[must_use]
pub fn get_residue_info_checked(idx: usize) -> Option<ResidueInfo> {
    RESIDUE_INFO_TABLE.get(idx).copied()
}

#[must_use]
pub fn find_tabulated_residue_idx(name: &str) -> usize {
    // BEGIN GEMMI CPP FUNCTION gemmi::find_tabulated_residue_idx
    // Gemmi✔️✔️: if (name.size() == 3) {
    // Gemmi✔️✔️: #define ID(s) (((s)[0] << 16 | (s)[1] << 8 | (s)[2]) & ~0x202020)
    // Gemmi✔️✔️:     switch (ID(name.c_str())) {
    // Gemmi✔️✔️:       case ID("TRY"):
    // Gemmi✔️✔️:       case ID("TRP"): return 23;
    // Gemmi✔️✔️:       case ID("WAT"):
    // Gemmi✔️✔️:       case ID("H2O"):
    // Gemmi✔️✔️:       case ID("HOH"): return 154;
    // Gemmi✔️✔️: }} else if (name.size() == 1) {
    // Gemmi✔️✔️:     switch (name[0]& ~0x20) {
    // Gemmi✔️✔️:   } else if (name.size() == 2) {
    // Gemmi✔️✔️:     if (name[0] == 'D' || name[0] == '+')
    // Gemmi✔️✔️:       switch (name[1]) {
    // Gemmi✔️✔️: #define ID(s) ((s)[0] << 8 | (s)[1])
    // Gemmi✔️✔️:     switch (ID(name.c_str())) {
    // Gemmi✔️✔️:     return 367;
    match name.len() {
        3 => match name.to_ascii_uppercase().as_str() {
            "ALA" => 0,
            "ARG" => 1,
            "ASN" => 2,
            "ABA" => 3,
            "ASP" => 4,
            "ASX" => 5,
            "CYS" => 6,
            "CSH" => 7,
            "GLN" => 8,
            "GLU" => 9,
            "GLX" => 10,
            "GLY" => 11,
            "HIS" => 12,
            "ILE" => 13,
            "LEU" => 14,
            "LYS" => 15,
            "MET" => 16,
            "MSE" => 17,
            "ORN" => 18,
            "PHE" => 19,
            "PRO" => 20,
            "SER" => 21,
            "THR" => 22,
            "TRY" => 23,
            "TRP" => 23,
            "TYR" => 24,
            "UNK" => 25,
            "VAL" => 26,
            "SEC" => 27,
            "PYL" => 28,
            "SEP" => 29,
            "TPO" => 30,
            "PCA" => 31,
            "CSO" => 32,
            "PTR" => 33,
            "KCX" => 34,
            "CSD" => 35,
            "LLP" => 36,
            "CME" => 37,
            "MLY" => 38,
            "DAL" => 39,
            "TYS" => 40,
            "OCS" => 41,
            "M3L" => 42,
            "FME" => 43,
            "ALY" => 44,
            "HYP" => 45,
            "CAS" => 46,
            "CRO" => 47,
            "CSX" => 48,
            "DPR" => 49,
            "DGL" => 50,
            "DVA" => 51,
            "CSS" => 52,
            "DPN" => 53,
            "DSN" => 54,
            "DLE" => 55,
            "HIC" => 56,
            "NLE" => 57,
            "MVA" => 58,
            "MLZ" => 59,
            "CR2" => 60,
            "SAR" => 61,
            "DAR" => 62,
            "DLY" => 63,
            "YCM" => 64,
            "NRQ" => 65,
            "CGU" => 66,
            "0TD" => 67,
            "MLE" => 68,
            "DAS" => 69,
            "DTR" => 70,
            "CXM" => 71,
            "TPQ" => 72,
            "DCY" => 73,
            "DSG" => 74,
            "DTY" => 75,
            "DHI" => 76,
            "MEN" => 77,
            "DTH" => 78,
            "SAC" => 79,
            "DGN" => 80,
            "AIB" => 81,
            "SMC" => 82,
            "IAS" => 83,
            "CIR" => 84,
            "BMT" => 85,
            "DIL" => 86,
            "FGA" => 87,
            "PHI" => 88,
            "CRQ" => 89,
            "SME" => 90,
            "GHP" => 91,
            "MHO" => 92,
            "NEP" => 93,
            "TRQ" => 94,
            "TOX" => 95,
            "ALC" => 96,
            "3FG" => 97,
            "SCH" => 98,
            "MDO" => 99,
            "MAA" => 100,
            "GYS" => 101,
            "MK8" => 102,
            "CR8" => 103,
            "KPI" => 104,
            "SCY" => 105,
            "DHA" => 106,
            "OMY" => 107,
            "CAF" => 108,
            "0AF" => 109,
            "SNN" => 110,
            "MHS" => 111,
            "MLU" => 112,
            "SNC" => 113,
            "PHD" => 114,
            "B3E" => 115,
            "MEA" => 116,
            "MED" => 117,
            "OAS" => 118,
            "GL3" => 119,
            "FVA" => 120,
            "PHL" => 121,
            "CRF" => 122,
            "OMZ" => 123,
            "BFD" => 124,
            "MEQ" => 125,
            "DAB" => 126,
            "AGM" => 127,
            "PSU" => 128,
            "5MU" => 129,
            "7MG" => 130,
            "OMG" => 131,
            "UR3" => 132,
            "OMC" => 133,
            "2MG" => 134,
            "H2U" => 135,
            "4SU" => 136,
            "OMU" => 137,
            "4OC" => 138,
            "MA6" => 139,
            "M2G" => 140,
            "1MA" => 141,
            "6MZ" => 142,
            "CCC" => 143,
            "2MA" => 144,
            "1MG" => 145,
            "5BU" => 146,
            "MIA" => 147,
            "DOC" => 148,
            "8OG" => 149,
            "5CM" => 150,
            "3DR" => 151,
            "BRU" => 152,
            "CBR" => 153,
            "WAT" => 154,
            "H2O" => 154,
            "HOH" => 154,
            "DOD" => 155,
            "HEM" => 156,
            "SO4" => 157,
            "GOL" => 158,
            "EDO" => 159,
            "NAG" => 160,
            "PO4" => 161,
            "ACT" => 162,
            "PEG" => 163,
            "MAN" => 164,
            "FAD" => 165,
            "BMA" => 166,
            "ADP" => 167,
            "DMS" => 168,
            "ACE" => 169,
            "NH2" => 170,
            "MPD" => 171,
            "MES" => 172,
            "NAD" => 173,
            "NAP" => 174,
            "TRS" => 175,
            "ATP" => 176,
            "PG4" => 177,
            "GDP" => 178,
            "FUC" => 179,
            "FMT" => 180,
            "GAL" => 181,
            "PGE" => 182,
            "FMN" => 183,
            "PLP" => 184,
            "EPE" => 185,
            "SF4" => 186,
            "BME" => 187,
            "CIT" => 188,
            "BE7" => 189,
            "MRD" => 190,
            "MHA" => 191,
            "BU3" => 192,
            "PGO" => 193,
            "BU2" => 194,
            "PDO" => 195,
            "BU1" => 196,
            "PG6" => 197,
            "1BO" => 198,
            "PE7" => 199,
            "PG5" => 200,
            "TFP" => 201,
            "DHD" => 202,
            "PEU" => 203,
            "TAU" => 204,
            "SBT" => 205,
            "SAL" => 206,
            "IOH" => 207,
            "IPA" => 208,
            "PIG" => 209,
            "B3P" => 210,
            "BTB" => 211,
            "NHE" => 212,
            "C8E" => 213,
            "OTE" => 214,
            "PE4" => 215,
            "XPE" => 216,
            "PE8" => 217,
            "P33" => 218,
            "N8E" => 219,
            "2OS" => 220,
            "1PS" => 221,
            "CPS" => 222,
            "DMX" => 223,
            "MPO" => 224,
            "GCD" => 225,
            "DXG" => 226,
            "CM5" => 227,
            "ACA" => 228,
            "ACN" => 229,
            "CCN" => 230,
            "GLC" => 231,
            "DR6" => 232,
            "NH4" => 233,
            "AZI" => 234,
            "BNG" => 235,
            "BOG" => 236,
            "BGC" => 237,
            "BCN" => 238,
            "BRO" => 239,
            "CAC" => 240,
            "CBX" => 241,
            "ACY" => 242,
            "CBM" => 243,
            "CLO" => 244,
            "3CO" => 245,
            "NCO" => 246,
            "CU1" => 247,
            "CYN" => 248,
            "MA4" => 249,
            "TAR" => 250,
            "GLO" => 251,
            "MTL" => 252,
            "SOR" => 253,
            "DMU" => 254,
            "DDQ" => 255,
            "DMF" => 256,
            "DIO" => 257,
            "DOX" => 258,
            "12P" => 259,
            "SDS" => 260,
            "LMT" => 261,
            "EOH" => 262,
            "EEE" => 263,
            "EGL" => 264,
            "FLO" => 265,
            "TRT" => 266,
            "FCY" => 267,
            "FRU" => 268,
            "GBL" => 269,
            "GPX" => 270,
            "HTO" => 271,
            "HTG" => 272,
            "B7G" => 273,
            "C10" => 274,
            "16D" => 275,
            "HEZ" => 276,
            "IOD" => 277,
            "IDO" => 278,
            "ICI" => 279,
            "ICT" => 280,
            "TLA" => 281,
            "LAT" => 282,
            "LBT" => 283,
            "LDA" => 284,
            "MN3" => 285,
            "MRY" => 286,
            "MOH" => 287,
            "BEQ" => 288,
            "C15" => 289,
            "MG8" => 290,
            "POL" => 291,
            "NO3" => 292,
            "JEF" => 293,
            "P4C" => 294,
            "CE1" => 295,
            "DIA" => 296,
            "CXE" => 297,
            "IPH" => 298,
            "PIN" => 299,
            "15P" => 300,
            "CRY" => 301,
            "PGR" => 302,
            "PGQ" => 303,
            "SPD" => 304,
            "SPK" => 305,
            "SPM" => 306,
            "SUC" => 307,
            "TBU" => 308,
            "TMA" => 309,
            "TEP" => 310,
            "SCN" => 311,
            "TRE" => 312,
            "ETF" => 313,
            "144" => 314,
            "UMQ" => 315,
            "URE" => 316,
            "YT3" => 317,
            "ZN2" => 318,
            "FE2" => 319,
            "3NI" => 320,
            "SIA" => 321,
            "XYP" => 322,
            "A2G" => 323,
            "GLA" => 324,
            "NDG" => 325,
            "NGA" => 326,
            _ => UNKNOWN_TABULATED_RESIDUE_IDX,
        },
        1 => match name.to_ascii_uppercase().as_str() {
            "A" => 327,
            "C" => 328,
            "G" => 329,
            "I" => 330,
            "U" => 331,
            "N" => 332,
            "F" => 333,
            "K" => 334,
            _ => UNKNOWN_TABULATED_RESIDUE_IDX,
        },
        2 => {
            let bytes = name.as_bytes();
            if bytes[0] == b'D' || bytes[0] == b'+' {
                match bytes[1] as char {
                    'A' => return 335,
                    'C' => return 336,
                    'G' => return 337,
                    'I' => return 338,
                    'T' => return 339,
                    'U' => return 340,
                    'N' => return 341,
                    _ => {}
                }
            }
            match name {
                "AG" => 342,
                "AL" => 343,
                "BA" => 344,
                "BR" => 345,
                "CA" => 346,
                "CD" => 347,
                "CL" => 348,
                "CM" => 349,
                "CN" => 350,
                "CO" => 351,
                "CS" => 352,
                "CU" => 353,
                "FE" => 354,
                "HG" => 355,
                "LI" => 356,
                "MG" => 357,
                "MN" => 358,
                "NA" => 359,
                "NI" => 360,
                "NO" => 361,
                "PB" => 362,
                "RB" => 363,
                "SR" => 364,
                "Y1" => 365,
                "ZN" => 366,
                _ => UNKNOWN_TABULATED_RESIDUE_IDX,
            }
        }
        _ => UNKNOWN_TABULATED_RESIDUE_IDX,
    }
    // END GEMMI CPP FUNCTION gemmi::find_tabulated_residue_idx
}

#[must_use]
pub fn find_tabulated_residue(name: &str) -> ResidueInfo {
    // Gemmi✔️✔️: ResidueInfo& find_tabulated_residue(const std::string& name);
    get_residue_info(find_tabulated_residue_idx(name))
}

#[must_use]
pub fn residue_code_from_name(name: &str) -> ResidueCode {
    find_tabulated_residue(name).code
}

#[must_use]
pub fn expand_one_letter(code: char, kind: ResidueInfoKind) -> Option<&'static str> {
    // BEGIN GEMMI CPP FUNCTION gemmi::expand_one_letter
    // Gemmi✔️✔️: static const char* names =
    // Gemmi✔️✔️:   "ALA\0ASX\0CYS\0ASP\0GLU\0PHE\0GLY\0HIS\0ILE\0\0   LYS\0LEU\0MET\0"
    // Gemmi✔️✔️:   "ASN\0PYL\0PRO\0GLN\0ARG\0SER\0THR\0SEC\0VAL\0TRP\0UNK\0TYR\0GLX\0"
    // Gemmi✔️✔️:   "DA\0 \0\0  DC\0 \0\0  \0\0  \0\0  DG\0 \0\0  DI\0 \0\0  \0\0  \0\0  \0\0  "
    // Gemmi✔️✔️:   "DN\0 \0\0  \0\0  \0\0  \0\0  \0\0  DT\0 DU\0 \0\0  \0\0  \0\0  \0\0  \0\0  ";
    // Gemmi✔️✔️: c &= ~0x20;
    // Gemmi✔️✔️: if (c >= 'A' && c <= 'Z') {
    // Gemmi✔️✔️:   ret = &names[4 * (c - 'A')];
    // Gemmi✔️✔️:   if (kind == ResidueKind::AA) {
    // Gemmi✔️✔️:   } else if (kind == ResidueKind::DNA) {
    // Gemmi✔️✔️:     ret += 4 * 26;
    // Gemmi✔️✔️:   } else if (kind == ResidueKind::RNA && c != 'T') {
    // Gemmi✔️✔️:     ret += 4 * 26 + 1;
    // Gemmi✔️✔️:   } else {
    // Gemmi✔️✔️:     ret = nullptr;
    // Gemmi✔️✔️: return (ret && *ret) ? ret : nullptr;
    let c = code.to_ascii_uppercase();
    match kind {
        ResidueInfoKind::Aa => match c {
            'A' => Some("ALA"),
            'B' => Some("ASX"),
            'C' => Some("CYS"),
            'D' => Some("ASP"),
            'E' => Some("GLU"),
            'F' => Some("PHE"),
            'G' => Some("GLY"),
            'H' => Some("HIS"),
            'I' => Some("ILE"),
            'K' => Some("LYS"),
            'L' => Some("LEU"),
            'M' => Some("MET"),
            'N' => Some("ASN"),
            'O' => Some("PYL"),
            'P' => Some("PRO"),
            'Q' => Some("GLN"),
            'R' => Some("ARG"),
            'S' => Some("SER"),
            'T' => Some("THR"),
            'U' => Some("SEC"),
            'V' => Some("VAL"),
            'W' => Some("TRP"),
            'X' => Some("UNK"),
            'Y' => Some("TYR"),
            'Z' => Some("GLX"),
            _ => None,
        },
        ResidueInfoKind::Dna => match c {
            'A' => Some("DA"),
            'C' => Some("DC"),
            'G' => Some("DG"),
            'I' => Some("DI"),
            'N' => Some("DN"),
            'T' => Some("DT"),
            'U' => Some("DU"),
            _ => None,
        },
        ResidueInfoKind::Rna if c != 'T' => match c {
            'A' => Some("A"),
            'C' => Some("C"),
            'G' => Some("G"),
            'I' => Some("I"),
            'N' => Some("N"),
            'U' => Some("U"),
            _ => None,
        },
        _ => None,
    }
    // END GEMMI CPP FUNCTION gemmi::expand_one_letter
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum ResidueInfoSequenceError {
    #[error("unmatched '(' in sequence")]
    UnmatchedParenthesis,
    #[error("unexpected letter in {kind} sequence: {letter} ({code})")]
    UnexpectedLetter {
        kind: &'static str,
        letter: char,
        code: u8,
    },
}

fn residue_sequence_kind_str(kind: ResidueInfoKind) -> &'static str {
    match kind {
        // Gemmi✔️✔️: case ResidueKind::AA: return "peptide";
        ResidueInfoKind::Aa => "peptide",
        // Gemmi✔️✔️: case ResidueKind::RNA: return "RNA";
        ResidueInfoKind::Rna => "RNA",
        // Gemmi✔️✔️: case ResidueKind::DNA: return "DNA";
        ResidueInfoKind::Dna => "DNA",
        // Gemmi✔️✔️: default: return "unknown";
        _ => "unknown",
    }
}

fn is_gemmi_space(byte: u8) -> bool {
    // Gemmi✔️✔️: inline bool is_space(char c) {
    // Gemmi✔️✔️:   static const std::uint8_t table[256] = { // 1 for 9-13 and 32
    // Gemmi✔️✔️:   return table[(std::uint8_t)c] != 0;
    matches!(byte, b'\t' | b'\n' | 0x0b | 0x0c | b'\r' | b' ')
}

#[must_use]
pub fn expand_protein_one_letter(code: char) -> Option<&'static str> {
    // Gemmi✔️✔️: inline const char* expand_protein_one_letter(char c) {
    // Gemmi✔️✔️:   return expand_one_letter(c, ResidueKind::AA);
    expand_one_letter(code, ResidueInfoKind::Aa)
}

pub fn expand_one_letter_sequence(
    seq: &str,
    kind: ResidueInfoKind,
) -> Result<Vec<String>, ResidueInfoSequenceError> {
    // BEGIN GEMMI CPP FUNCTION gemmi::expand_one_letter_sequence
    // Gemmi✔️✔️: std::vector<std::string> expand_one_letter_sequence(const std::string& seq,
    // Gemmi✔️✔️:                                                     ResidueKind kind) {
    // Gemmi✔️✔️:  std::vector<std::string> r;
    // Gemmi✔️✔️:   r.reserve(seq.size());
    // Gemmi✔️✔️:   auto kind_str = [&]() {
    // Gemmi✔️✔️:     switch (kind) {
    // Gemmi✔️✔️:       case ResidueKind::AA: return "peptide";
    // Gemmi✔️✔️:       case ResidueKind::RNA: return "RNA";
    // Gemmi✔️✔️:       case ResidueKind::DNA: return "DNA";
    // Gemmi✔️✔️:       default: return "unknown";
    // Gemmi✔️✔️:   for (size_t i = 0; i != seq.size(); ++i) {
    // Gemmi✔️✔️:     char c = seq[i];
    // Gemmi✔️✔️:     if (is_space(c))
    // Gemmi✔️✔️:       continue;
    // Gemmi✔️✔️:     if (c == '(') { // special case, e.g. (MSE)
    // Gemmi✔️✔️:       size_t start = i + 1;
    // Gemmi✔️✔️:       i = seq.find(')', start);
    // Gemmi✔️✔️:       if (i == std::string::npos)
    // Gemmi✔️✔️:         gemmi::fail("unmatched '(' in sequence");
    // Gemmi✔️✔️:       r.emplace_back(seq, start, i - start);
    // Gemmi✔️✔️:     } else {
    // Gemmi✔️✔️:       const char* str = gemmi::expand_one_letter(c, kind);
    // Gemmi✔️✔️:       if (str == nullptr)
    // Gemmi✔️✔️:         gemmi::fail("unexpected letter in ", kind_str(), " sequence: ", c,
    // Gemmi✔️✔️:              " (", std::to_string(int(c)), ')');
    // Gemmi✔️✔️:       r.emplace_back(str);
    // Gemmi✔️✔️:   return r;
    let bytes = seq.as_bytes();
    let mut residues = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i != bytes.len() {
        let c = bytes[i];
        if is_gemmi_space(c) {
            i += 1;
            continue;
        }
        if c == b'(' {
            let start = i + 1;
            let Some(offset) = bytes[start..].iter().position(|&byte| byte == b')') else {
                return Err(ResidueInfoSequenceError::UnmatchedParenthesis);
            };
            let end = start + offset;
            residues.push(seq[start..end].to_string());
            i = end + 1;
        } else if let Some(name) = expand_one_letter(c as char, kind) {
            residues.push(name.to_string());
            i += 1;
        } else {
            return Err(ResidueInfoSequenceError::UnexpectedLetter {
                kind: residue_sequence_kind_str(kind),
                letter: c as char,
                code: c,
            });
        }
    }
    Ok(residues)
    // END GEMMI CPP FUNCTION gemmi::expand_one_letter_sequence
}

pub fn expand_protein_one_letter_string(
    seq: &str,
) -> Result<Vec<String>, ResidueInfoSequenceError> {
    // Gemmi✔️✔️: std::vector<std::string> expand_protein_one_letter_string(const std::string& s);
    expand_one_letter_sequence(seq, ResidueInfoKind::Aa)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gemmi_residue_info_table_shape_and_unknown_row_match_source() {
        assert_eq!(RESIDUE_INFO_TABLE.len(), 368);
        assert_eq!(UNKNOWN_TABULATED_RESIDUE_IDX, 367);
        let unknown = get_residue_info(UNKNOWN_TABULATED_RESIDUE_IDX);
        assert_eq!(unknown.code, ResidueCode::UNKNOWN);
        assert_eq!(unknown.name, "");
        assert_eq!(unknown.kind, ResidueInfoKind::Unknown);
        assert!(!unknown.found());
        assert_eq!(get_residue_info_checked(999), None);
    }

    #[test]
    fn every_residue_code_maps_to_its_own_table_row_and_stable_name() {
        for (index, info) in RESIDUE_INFO_TABLE.iter().copied().enumerate() {
            assert_eq!(usize::from(info.code.as_u16()), index);
            assert_eq!(info.code.info(), info);
            if matches!(info.code, ResidueCode::UNKNOWN) {
                assert_eq!(info.code.name(), "UNKNOWN");
                assert_eq!(info.code.to_string(), "UNKNOWN");
                assert_eq!("UNKNOWN".parse(), Ok(ResidueCode::UNKNOWN));
            } else {
                assert_eq!(info.code.name(), info.name);
                assert_eq!(info.code.to_string(), info.name);
                assert_eq!(info.name.parse(), Ok(info.code));
            }
        }
    }

    #[test]
    fn residue_code_parsing_keeps_source_aliases_but_rejects_unknown_names() {
        assert_eq!("TRY".parse(), Ok(ResidueCode::TRP));
        assert_eq!("wat".parse(), Ok(ResidueCode::HOH));
        assert_eq!("H2O".parse(), Ok(ResidueCode::HOH));
        assert_eq!("+A".parse(), Ok(ResidueCode::DA));
        assert_eq!("unknown".parse(), Ok(ResidueCode::UNKNOWN));

        let error = "not-tabulated"
            .parse::<ResidueCode>()
            .expect_err("an unrecognized residue name must not become UNKNOWN implicitly");
        assert_eq!(error.input(), "not-tabulated");
    }

    #[test]
    fn residue_code_and_info_serde_use_stable_source_names() {
        assert_eq!(
            serde_json::to_string(&ResidueCode::MSE).unwrap(),
            r#""MSE""#
        );
        assert_eq!(
            serde_json::to_string(&ResidueCode::R3FG).unwrap(),
            r#""3FG""#
        );
        assert_eq!(
            serde_json::to_string(&ResidueCode::UNKNOWN).unwrap(),
            r#""UNKNOWN""#
        );
        assert_eq!(
            serde_json::from_str::<ResidueCode>(r#""MSE""#).unwrap(),
            ResidueCode::MSE
        );
        assert_eq!(
            serde_json::from_str::<ResidueCode>(r#""3FG""#).unwrap(),
            ResidueCode::R3FG
        );
        assert!(serde_json::from_str::<ResidueCode>(r#""XYZ""#).is_err());

        let serialized = serde_json::to_value(ResidueCode::MSE.info()).unwrap();
        assert_eq!(serialized["code"], "MSE");
        assert_eq!(serialized["name"], "MSE");
        assert_eq!(serialized["kind"], "AA");
        assert_eq!(serialized["one_letter_code"], "m");
        assert_eq!(serialized["hydrogen_count"], 11);
    }

    #[test]
    fn modified_amino_acid_metadata_maps_only_through_source_one_letter_codes() {
        for (modified, one_letter, parent) in [
            (ResidueCode::MSE, 'M', ResidueCode::MET),
            (ResidueCode::HYP, 'P', ResidueCode::PRO),
            (ResidueCode::SEP, 'S', ResidueCode::SER),
        ] {
            let info = modified.info();
            assert!(info.is_modified_amino_acid());
            assert_eq!(info.canonical_one_letter_code(), Some(one_letter));
            assert_eq!(info.parent_standard_code(), Some(parent));
        }

        let proline = ResidueCode::PRO.info();
        assert!(!proline.is_modified_amino_acid());
        assert_eq!(proline.canonical_one_letter_code(), Some('P'));
        assert_eq!(proline.parent_standard_code(), Some(ResidueCode::PRO));

        let unmappable_amino_acid = ResidueCode::R3FG.info();
        assert!(unmappable_amino_acid.is_modified_amino_acid());
        assert_eq!(unmappable_amino_acid.canonical_one_letter_code(), None);
        assert_eq!(unmappable_amino_acid.parent_standard_code(), None);

        let water = ResidueCode::HOH.info();
        assert!(!water.is_modified_amino_acid());
        assert_eq!(water.canonical_one_letter_code(), None);
        assert_eq!(water.parent_standard_code(), None);
    }

    #[test]
    fn residue_identity_preserves_raw_known_aliases_and_unknown_names() {
        let alias = ResidueIdentity::new("wat");
        assert_eq!(alias.name(), "wat");
        assert_eq!(alias.code(), ResidueCode::HOH);
        assert_eq!(alias.info(), ResidueCode::HOH.info());
        assert!(alias.is_tabulated());

        let unknown = ResidueIdentity::new("CUSTOM_COMPONENT");
        assert_eq!(unknown.name(), "CUSTOM_COMPONENT");
        assert_eq!(unknown.code(), ResidueCode::UNKNOWN);
        assert!(!unknown.is_tabulated());
        assert_eq!(unknown.to_string(), "CUSTOM_COMPONENT");
        assert_eq!(
            serde_json::to_string(&unknown).unwrap(),
            r#""CUSTOM_COMPONENT""#
        );
        assert_eq!(
            serde_json::from_str::<ResidueIdentity>(r#""CUSTOM_COMPONENT""#).unwrap(),
            unknown
        );
    }

    #[test]
    fn gemmi_residue_info_selected_rows_preserve_source_fields() {
        let gln = find_tabulated_residue("GLN");
        assert_eq!(gln.code, ResidueCode::GLN);
        assert_eq!(gln.kind, ResidueInfoKind::Aa);
        assert_eq!(gln.linking_type, 1);
        assert_eq!(gln.one_letter_code, 'Q');
        assert_eq!(gln.hydrogen_count, 10);
        assert_eq!(gln.weight, 146.144f32);
        assert!(gln.found());
        assert!(gln.is_amino_acid());
        assert!(gln.is_peptide_linking());
        assert!(!gln.is_na_linking());
        assert!(gln.is_standard());
        assert_eq!(gln.fasta_code(), 'Q');

        let mse = find_tabulated_residue("MSE");
        assert_eq!(mse.code, ResidueCode::MSE);
        assert_eq!(mse.one_letter_code, 'm');
        assert!(!mse.is_standard());
        assert_eq!(mse.fasta_code(), 'X');

        let hoh = find_tabulated_residue("HOH");
        assert_eq!(hoh.code, ResidueCode::HOH);
        assert!(hoh.is_water());
        assert!(hoh.is_buffer_or_water());
        assert_eq!(hoh.one_letter_code, ' ');
        assert_eq!(hoh.fasta_code(), 'X');

        let da = find_tabulated_residue("DA");
        assert_eq!(da.code, ResidueCode::DA);
        assert!(da.is_dna());
        assert!(da.is_nucleic_acid());
        assert!(da.is_na_linking());
    }

    #[test]
    fn gemmi_find_tabulated_residue_alias_and_case_behavior_matches_source() {
        assert_eq!(find_tabulated_residue_idx("TRY"), 23);
        assert_eq!(find_tabulated_residue_idx("trp"), 23);
        assert_eq!(find_tabulated_residue_idx("WAT"), 154);
        assert_eq!(find_tabulated_residue_idx("h2o"), 154);

        assert_eq!(find_tabulated_residue_idx("a"), 327);
        assert_eq!(find_tabulated_residue_idx("+A"), 335);
        assert_eq!(find_tabulated_residue_idx("DA"), 335);

        assert_eq!(
            find_tabulated_residue_idx("ag"),
            UNKNOWN_TABULATED_RESIDUE_IDX
        );
        assert_eq!(
            find_tabulated_residue_idx("da"),
            UNKNOWN_TABULATED_RESIDUE_IDX
        );
        assert_eq!(
            find_tabulated_residue_idx("XYZ"),
            UNKNOWN_TABULATED_RESIDUE_IDX
        );
    }

    #[test]
    fn gemmi_expand_one_letter_matches_source_offsets() {
        assert_eq!(expand_one_letter('q', ResidueInfoKind::Aa), Some("GLN"));
        assert_eq!(expand_one_letter('J', ResidueInfoKind::Aa), None);
        assert_eq!(expand_one_letter('t', ResidueInfoKind::Dna), Some("DT"));
        assert_eq!(expand_one_letter('u', ResidueInfoKind::Dna), Some("DU"));
        assert_eq!(expand_one_letter('u', ResidueInfoKind::Rna), Some("U"));
        assert_eq!(expand_one_letter('t', ResidueInfoKind::Rna), None);
        assert_eq!(expand_one_letter('A', ResidueInfoKind::Buf), None);
    }

    #[test]
    fn gemmi_expand_one_letter_sequence_matches_source_behavior() {
        assert_eq!(
            expand_one_letter_sequence("ACD (MSE)\nU", ResidueInfoKind::Aa).unwrap(),
            ["ALA", "CYS", "ASP", "MSE", "SEC"]
        );
        assert_eq!(
            expand_one_letter_sequence("ACGTU", ResidueInfoKind::Dna).unwrap(),
            ["DA", "DC", "DG", "DT", "DU"]
        );
        assert_eq!(
            expand_one_letter_sequence("ACGU", ResidueInfoKind::Rna).unwrap(),
            ["A", "C", "G", "U"]
        );
        assert_eq!(
            expand_one_letter_sequence("T", ResidueInfoKind::Rna).unwrap_err(),
            ResidueInfoSequenceError::UnexpectedLetter {
                kind: "RNA",
                letter: 'T',
                code: b'T'
            }
        );
        assert_eq!(
            expand_one_letter_sequence("(MSE", ResidueInfoKind::Aa).unwrap_err(),
            ResidueInfoSequenceError::UnmatchedParenthesis
        );
        assert_eq!(
            expand_protein_one_letter_string("m").unwrap(),
            ["MET".to_string()]
        );
    }
}
