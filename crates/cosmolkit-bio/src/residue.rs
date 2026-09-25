//! Gemmi-derived residue information table.
//!
//! Source: `third_party/gemmi/src/resinfo.cpp` and `third_party/gemmi/include/gemmi/resinfo.hpp`.

#![allow(non_camel_case_types)]

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
    #[must_use]
    pub const fn as_u16(self) -> u16 {
        self as u16
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResidueInfo {
    pub code: ResidueCode,
    pub name: &'static str,
    pub kind: ResidueInfoKind,
    pub linking_type: u8,
    pub one_letter_code: char,
    pub hydrogen_count: u8,
    pub weight: f32,
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

pub const UNKNOWN_TABULATED_RESIDUE_INDEX: usize = 367;

// BEGIN GEMMI CPP TABLE gemmi::residue_info
// Gemmi✔️✔️: static ResidueInfo array[368] = {
const RESIDUE_INFO_TABLE: [ResidueInfo; 368] = [
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
// END GEMMI CPP TABLE gemmi::residue_info

#[must_use]
pub fn residue_info(idx: usize) -> ResidueInfo {
    // Gemmi✔️✔️: ResidueInfo& get_residue_info(size_t idx) {
    // Gemmi✔️✔️:   static ResidueInfo array[368] = {
    // Gemmi✔️✔️:     // hydrogen_count needs to be verified
    // Gemmi✔️✔️:     {"ALA", RI::AA,  1, 'A',   7, 89.0932f },
    // Gemmi✔️✔️:     {"ARG", RI::AA,  1, 'R',  15, 175.209f },
    // Gemmi✔️✔️:     {"ASN", RI::AA,  1, 'N',   8, 132.118f },
    // Gemmi✔️✔️:     {"ABA", RI::AA,  1, 'a',   9, 103.120f },
    // Gemmi✔️✔️:     {"ASP", RI::AA,  1, 'D',   7, 133.103f },
    // Gemmi✔️✔️:     {"ASX", RI::AA,  1, 'B',   6, 100.096f },
    // Gemmi✔️✔️:     {"CYS", RI::AA,  1, 'C',   7, 121.158f },  // also BUF
    // Gemmi✔️✔️:     {"CSH", RI::AA,  1, 's',  17, 283.284f },
    // Gemmi✔️✔️:     {"GLN", RI::AA,  1, 'Q',  10, 146.144f },
    // Gemmi✔️✔️:     {"GLU", RI::AA,  1, 'E',   9, 147.129f },
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"GLX", RI::AA,  1, 'Z',   8, 114.123f },
    // Gemmi✔️✔️:     {"GLY", RI::AA,  1, 'G',   5, 75.0666f },  // also BUF
    // Gemmi✔️✔️:     {"HIS", RI::AA,  1, 'H',  10, 156.162f },
    // Gemmi✔️✔️:     {"ILE", RI::AA,  1, 'I',  13, 131.173f },
    // Gemmi✔️✔️:     {"LEU", RI::AA,  1, 'L',  13, 131.173f },
    // Gemmi✔️✔️:     {"LYS", RI::AA,  1, 'K',  15, 147.196f },
    // Gemmi✔️✔️:     {"MET", RI::AA,  1, 'M',  11, 149.211f },
    // Gemmi✔️✔️:     {"MSE", RI::AA,  1, 'm',  11, 196.106f },
    // Gemmi✔️✔️:     {"ORN", RI::AA,  1, 'a',  12, 132.161f },
    // Gemmi✔️✔️:     {"PHE", RI::AA,  1, 'F',  11, 165.189f },
    // Gemmi✔️✔️:     //20
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"PRO", RI::AA,  1, 'P',   9, 115.130f },
    // Gemmi✔️✔️:     {"SER", RI::AA,  1, 'S',   7, 105.093f },
    // Gemmi✔️✔️:     {"THR", RI::AA,  1, 'T',   9, 119.119f },
    // Gemmi✔️✔️:     {"TRP", RI::AA,  1, 'W',  12, 204.225f },
    // Gemmi✔️✔️:     {"TYR", RI::AA,  1, 'Y',  11, 181.189f },
    // Gemmi✔️✔️:     {"UNK", RI::AA,  1, 'X',   9, 103.120f },
    // Gemmi✔️✔️:     {"VAL", RI::AA,  1, 'V',  11, 117.146f },
    // Gemmi✔️✔️:     {"SEC", RI::AA,  1, 'U',   7, 168.053f },
    // Gemmi✔️✔️:     {"PYL", RI::AA,  1, 'O',  21, 255.313f },
    // Gemmi✔️✔️:     {"SEP", RI::AA,  1, 's',   8, 185.072f },
    // Gemmi✔️✔️:     // 30
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"TPO", RI::AA,  1, 't',  10, 199.099f },
    // Gemmi✔️✔️:     {"PCA", RI::AA,  1, 'e',   7, 129.114f },
    // Gemmi✔️✔️:     {"CSO", RI::AA,  1, 'c',   7, 137.158f },
    // Gemmi✔️✔️:     {"PTR", RI::AA,  1, 'y',  12, 261.168f },
    // Gemmi✔️✔️:     {"KCX", RI::AA,  1, 'k',  14, 190.197f },
    // Gemmi✔️✔️:     {"CSD", RI::AA,  1, 'c',   7, 153.157f },
    // Gemmi✔️✔️:     {"LLP", RI::AA,  1, 'k',  22, 375.314f },
    // Gemmi✔️✔️:     {"CME", RI::AA,  1, 'c',  11, 197.276f },
    // Gemmi✔️✔️:     {"MLY", RI::AA,  1, 'k',  18, 174.241f },
    // Gemmi✔️✔️:     {"DAL", RI::AAD, 1, 'a',   7, 89.0932f },
    // Gemmi✔️✔️:     {"TYS", RI::AA,  1, 'y',  11, 261.252f },
    // Gemmi✔️✔️:     {"OCS", RI::AA,  1, 'c',   7, 169.156f },
    // Gemmi✔️✔️:     // 40
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"M3L", RI::AA,  1, 'k',  21, 189.275f },
    // Gemmi✔️✔️:     {"FME", RI::AA,  1, 'm',  11, 177.221f },
    // Gemmi✔️✔️:     {"ALY", RI::AA,  1, 'k',  16, 188.224f },
    // Gemmi✔️✔️:     {"HYP", RI::AA,  1, 'p',   9, 131.130f },
    // Gemmi✔️✔️:     {"CAS", RI::AA,  1, 'c',  12, 225.141f },
    // Gemmi✔️✔️:     {"CRO", RI::AA,  1, 't',  17, 319.313f },
    // Gemmi✔️✔️:     {"CSX", RI::AA,  1, 'c',   7, 137.158f },
    // Gemmi✔️✔️:     {"DPR", RI::AAD, 1, 'p',   9, 115.130f },  // also BUF
    // Gemmi✔️✔️:     {"DGL", RI::AAD, 1, 'e',   9, 147.129f },
    // Gemmi✔️✔️:     {"DVA", RI::AAD, 1, 'v',  11, 117.146f },
    // Gemmi✔️✔️:     {"CSS", RI::AA,  1, 'c',   7, 153.223f },
    // Gemmi✔️✔️:     {"DPN", RI::AAD, 1, 'f',  11, 165.189f },
    // Gemmi✔️✔️:     {"DSN", RI::AAD, 1, 's',   7, 105.093f },
    // Gemmi✔️✔️:     // 50
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"DLE", RI::AAD, 1, 'l',  13, 131.173f },
    // Gemmi✔️✔️:     {"HIC", RI::AA,  1, 'h',  11, 169.181f },
    // Gemmi✔️✔️:     {"NLE", RI::AA,  1, 'l',  13, 131.173f },
    // Gemmi✔️✔️:     {"MVA", RI::AA,  1, 'v',  13, 131.173f },
    // Gemmi✔️✔️:     {"MLZ", RI::AA,  1, 'k',  16, 160.214f },
    // Gemmi✔️✔️:     {"CR2", RI::AA,  1, 'g',  13, 275.260f },
    // Gemmi✔️✔️:     {"SAR", RI::AA,  1, 'g',   7, 89.0932f },
    // Gemmi✔️✔️:     {"DAR", RI::AAD, 1, 'r',  15, 175.209f },
    // Gemmi✔️✔️:     {"DLY", RI::AAD, 1, 'k',  14, 146.188f },
    // Gemmi✔️✔️:     {"YCM", RI::AA,  1, 'c',  10, 178.209f },
    // Gemmi✔️✔️:     // 60
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"NRQ", RI::AA,  1, 'm',  17, 347.389f },
    // Gemmi✔️✔️:     {"CGU", RI::AA,  1, 'e',   9, 191.139f },
    // Gemmi✔️✔️:     {"0TD", RI::AA,  1, 'd',   9, 179.194f },
    // Gemmi✔️✔️:     {"MLE", RI::AA,  1, 'l',  15, 145.200f },
    // Gemmi✔️✔️:     {"DAS", RI::AAD, 1, 'd',   7, 133.103f },
    // Gemmi✔️✔️:     {"DTR", RI::AAD, 1, 'w',  12, 204.225f },
    // Gemmi✔️✔️:     {"CXM", RI::AA,  1, 'm',  11, 193.221f },
    // Gemmi✔️✔️:     {"TPQ", RI::AA,  1, 'y',   9, 211.171f },
    // Gemmi✔️✔️:     {"DCY", RI::AAD, 1, 'c',   7, 121.158f },
    // Gemmi✔️✔️:     {"DSG", RI::AAD, 1, 'n',   8, 132.118f },
    // Gemmi✔️✔️:     {"DTY", RI::AAD, 1, 'y',  11, 181.189f },
    // Gemmi✔️✔️:     // 70
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"DHI", RI::AAD, 1, 'h',  10, 156.162f },
    // Gemmi✔️✔️:     {"MEN", RI::AA,  1, 'n',  10, 146.144f },
    // Gemmi✔️✔️:     {"DTH", RI::AAD, 1, 't',   9, 119.119f },
    // Gemmi✔️✔️:     {"SAC", RI::AA,  1, 's',   9, 147.129f },
    // Gemmi✔️✔️:     {"DGN", RI::AAD, 1, 'q',  10, 146.144f },
    // Gemmi✔️✔️:     {"AIB", RI::AA,  1, 'a',   9, 103.120f },
    // Gemmi✔️✔️:     {"SMC", RI::AA,  1, 'c',   9, 135.185f },
    // Gemmi✔️✔️:     {"IAS", RI::AA,  1, 'd',   7, 133.103f },
    // Gemmi✔️✔️:     {"CIR", RI::AA,  1, 'r',  13, 175.186f },
    // Gemmi✔️✔️:     {"BMT", RI::AA,  1, 't',  19, 201.263f },
    // Gemmi✔️✔️:     {"DIL", RI::AAD, 1, 'i',  13, 131.173f },
    // Gemmi✔️✔️:     {"FGA", RI::AA,  1, 'e',   9, 147.129f },
    // Gemmi✔️✔️:     //80
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"PHI", RI::AA,  1, 'f',  10, 291.086f },
    // Gemmi✔️✔️:     {"CRQ", RI::AA,  1, 'q',  16, 344.322f },
    // Gemmi✔️✔️:     {"SME", RI::AA,  1, 'm',  11, 165.211f },
    // Gemmi✔️✔️:     {"GHP", RI::AA,  1, 'g',   9, 167.162f },  // d-peptide in CCD
    // Gemmi✔️✔️:     {"MHO", RI::AA,  1, 'm',  11, 165.211f },
    // Gemmi✔️✔️:     {"NEP", RI::AA,  1, 'h',  10, 235.134f },
    // Gemmi✔️✔️:     {"TRQ", RI::AA,  1, 'w',  10, 234.208f },
    // Gemmi✔️✔️:     {"TOX", RI::AA,  1, 'w',  12, 236.224f },
    // Gemmi✔️✔️:     {"ALC", RI::AA,  1, 'a',  17, 171.237f },
    // Gemmi✔️✔️:     {"3FG", RI::AA,  1, ' ',   9, 183.161f },
    // Gemmi✔️✔️:     {"SCH", RI::AA,  1, 'c',   9, 167.250f },
    // Gemmi✔️✔️:     {"MDO", RI::AA,  1, 'a',  11, 197.191f },
    // Gemmi✔️✔️:     {"MAA", RI::AA,  1, 'a',   9, 103.120f },
    // Gemmi✔️✔️:     //90
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"GYS", RI::AA,  1, 's',  15, 305.286f },
    // Gemmi✔️✔️:     {"MK8", RI::AA,  1, 'l',  15, 145.200f },
    // Gemmi✔️✔️:     {"CR8", RI::AA,  1, 'h',  16, 354.340f },
    // Gemmi✔️✔️:     {"KPI", RI::AA,  1, 'k',  16, 216.234f },
    // Gemmi✔️✔️:     {"SCY", RI::AA,  1, 'c',   9, 163.195f },
    // Gemmi✔️✔️:     {"DHA", RI::AA,  1, 's',   5, 87.0773f },
    // Gemmi✔️✔️:     {"OMY", RI::AA,  1, 'y',  10, 231.633f },
    // Gemmi✔️✔️:     {"CAF", RI::AA,  1, 'c',  12, 241.140f },
    // Gemmi✔️✔️:     {"0AF", RI::AA,  1, 'w',  12, 220.225f },
    // Gemmi✔️✔️:     {"SNN", RI::AA,  1, 'n',   6, 114.103f },
    // Gemmi✔️✔️:     // 100
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"MHS", RI::AA,  1, 'h',  11, 169.181f },
    // Gemmi✔️✔️:     {"MLU", RI::AAD, 1, ' ',  15, 145.200f },
    // Gemmi✔️✔️:     {"SNC", RI::AA,  1, 'c',   6, 150.156f },
    // Gemmi✔️✔️:     {"PHD", RI::AA,  1, 'd',   8, 213.083f },
    // Gemmi✔️✔️:     {"B3E", RI::AA,  1, 'e',  11, 161.156f },
    // Gemmi✔️✔️:     {"MEA", RI::AA,  1, 'f',  13, 179.216f },
    // Gemmi✔️✔️:     {"MED", RI::AAD, 1, 'm',  11, 149.211f },
    // Gemmi✔️✔️:     {"OAS", RI::AA,  1, 's',   9, 147.129f },
    // Gemmi✔️✔️:     {"GL3", RI::AA,  1, 'g',   5, 91.1322f },
    // Gemmi✔️✔️:     {"FVA", RI::AA,  1, 'v',  11, 145.156f },
    // Gemmi✔️✔️:     // 110
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"PHL", RI::AA,  1, 'f',  13, 151.206f },
    // Gemmi✔️✔️:     {"CRF", RI::AA,  1, 't',  18, 342.349f },
    // Gemmi✔️✔️:     {"OMZ", RI::AA,  1, ' ',  10, 231.633f },  // d-peptide in CCD
    // Gemmi✔️✔️:     {"BFD", RI::AA,  1, 'd',   6, 198.102f },
    // Gemmi✔️✔️:     {"MEQ", RI::AA,  1, 'q',  12, 160.171f },
    // Gemmi✔️✔️:     {"DAB", RI::AA,  1, 'a',  10, 118.134f },
    // Gemmi✔️✔️:     {"AGM", RI::AA,  1, 'r',  17, 189.235f },
    // Gemmi✔️✔️:     {"PSU", RI::RNA, 2, 'u',  13, 324.181f },
    // Gemmi✔️✔️:     {"5MU", RI::RNA, 2, 'u',  15, 338.208f },
    // Gemmi✔️✔️:     {"7MG", RI::RNA, 2, 'g',  18, 379.263f },
    // Gemmi✔️✔️:     // 120
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"OMG", RI::RNA, 2, 'g',  16, 377.247f },
    // Gemmi✔️✔️:     {"UR3", RI::RNA, 2, 'u',  15, 338.208f },
    // Gemmi✔️✔️:     {"OMC", RI::RNA, 2, 'c',  16, 337.223f },
    // Gemmi✔️✔️:     {"2MG", RI::RNA, 2, 'g',  16, 377.247f },
    // Gemmi✔️✔️:     {"H2U", RI::RNA, 2, 'u',  15, 326.197f },
    // Gemmi✔️✔️:     {"4SU", RI::RNA, 2, 'u',  13, 340.247f },
    // Gemmi✔️✔️:     {"OMU", RI::RNA, 2, 'u',  15, 338.208f },
    // Gemmi✔️✔️:     {"4OC", RI::RNA, 2, 'c',  18, 351.250f },
    // Gemmi✔️✔️:     {"MA6", RI::RNA, 2, 'a',  18, 375.274f },
    // Gemmi✔️✔️:     {"M2G", RI::RNA, 2, 'g',  18, 391.274f },
    // Gemmi✔️✔️:     // 130
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"1MA", RI::RNA, 2, 'a',  16, 361.248f },
    // Gemmi✔️✔️:     {"6MZ", RI::RNA, 2, 'a',  16, 361.248f },
    // Gemmi✔️✔️:     {"CCC", RI::RNA, 2, 'c',  13, 385.161f },
    // Gemmi✔️✔️:     {"2MA", RI::RNA, 2, 'a',  16, 361.248f },
    // Gemmi✔️✔️:     {"1MG", RI::RNA, 2, 'g',  16, 377.247f },
    // Gemmi✔️✔️:     {"5BU", RI::RNA, 2, 'u',  12, 403.077f },
    // Gemmi✔️✔️:     {"MIA", RI::RNA, 2, 'a',  24, 461.430f },
    // Gemmi✔️✔️:     {"DOC", RI::DNA, 2, 'c',  14, 291.198f },
    // Gemmi✔️✔️:     {"8OG", RI::DNA, 2, 'g',  14, 363.221f },
    // Gemmi✔️✔️:     {"5CM", RI::DNA, 2, 'c',  16, 321.224f },
    // Gemmi✔️✔️:     // 140
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"3DR", RI::DNA, 2, ' ',  11, 198.111f },
    // Gemmi✔️✔️:     {"BRU", RI::DNA, 2, 'u',  12, 387.078f },
    // Gemmi✔️✔️:     {"CBR", RI::DNA, 2, 'c',  13, 386.093f },
    // Gemmi✔️✔️:     {"HOH", RI::HOH, 0, ' ',   2, 18.0153f },
    // Gemmi✔️✔️:     {"DOD", RI::HOH, 0, ' ',   2, 20.0276f },
    // Gemmi✔️✔️:     {"HEM", RI::ELS, 0, ' ',  32, 616.487f },
    // Gemmi✔️✔️:     {"SO4", RI::BUF, 0, ' ',   0, 96.0626f },
    // Gemmi✔️✔️:     {"GOL", RI::BUF, 0, ' ',   8, 92.0938f },
    // Gemmi✔️✔️:     {"EDO", RI::BUF, 0, ' ',   6, 62.0678f },
    // Gemmi✔️✔️:     {"NAG", RI::PYR, 0, ' ',  15, 221.208f },
    // Gemmi✔️✔️:     // 150
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"PO4", RI::ELS, 0, ' ',   0, 94.9714f },
    // Gemmi✔️✔️:     {"ACT", RI::BUF, 0, ' ',   3, 59.0440f },
    // Gemmi✔️✔️:     {"PEG", RI::ELS, 0, ' ',  10, 106.120f },
    // Gemmi✔️✔️:     {"MAN", RI::PYR, 0, ' ',  12, 180.156f },  // also BUF
    // Gemmi✔️✔️:     {"FAD", RI::ELS, 0, ' ',  33, 785.550f },
    // Gemmi✔️✔️:     {"BMA", RI::PYR, 0, ' ',  12, 180.156f },  // also BUF
    // Gemmi✔️✔️:     {"ADP", RI::ELS, 0, ' ',  15, 427.201f },
    // Gemmi✔️✔️:     {"DMS", RI::BUF, 0, ' ',   6, 78.1334f },
    // Gemmi✔️✔️:     {"ACE", RI::ELS, 1, ' ',   4, 44.0526f },
    // Gemmi✔️✔️:     {"NH2", RI::ELS, 1, ' ',   2, 16.0226f },  // ?
    // Gemmi✔️✔️:     {"MPD", RI::BUF, 0, ' ',  14, 118.174f },
    // Gemmi✔️✔️:     {"MES", RI::ELS, 0, ' ',  13, 195.237f },
    // Gemmi✔️✔️:     // 160
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"NAD", RI::ELS, 0, ' ',  27, 663.425f },
    // Gemmi✔️✔️:     {"NAP", RI::ELS, 0, ' ',  28, 743.405f },
    // Gemmi✔️✔️:     {"TRS", RI::BUF, 0, ' ',  12, 122.143f },
    // Gemmi✔️✔️:     {"ATP", RI::ELS, 0, ' ',  16, 507.181f },
    // Gemmi✔️✔️:     {"PG4", RI::ELS, 0, ' ',  18, 194.226f },
    // Gemmi✔️✔️:     {"GDP", RI::ELS, 2, 'g',  15, 443.201f },  // RNA in CCD
    // Gemmi✔️✔️:     {"FUC", RI::PYR, 0, ' ',  12, 164.156f },
    // Gemmi✔️✔️:     {"FMT", RI::BUF, 0, ' ',   2, 46.0254f },
    // Gemmi✔️✔️:     {"GAL", RI::PYR, 0, ' ',  12, 180.156f },
    // Gemmi✔️✔️:     {"PGE", RI::BUF, 0, ' ',  14, 150.173f },
    // Gemmi✔️✔️:     {"FMN", RI::ELS, 0, ' ',  21, 456.344f },
    // Gemmi✔️✔️:     // 170
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"PLP", RI::ELS, 0, ' ',  10, 247.142f },
    // Gemmi✔️✔️:     {"EPE", RI::ELS, 0, ' ',  18, 238.305f },
    // Gemmi✔️✔️:     {"SF4", RI::ELS, 0, ' ',   0, 351.640f },
    // Gemmi✔️✔️:     {"BME", RI::ELS, 0, ' ',   6, 78.1334f },
    // Gemmi✔️✔️:     {"CIT", RI::BUF, 0, ' ',   8, 192.124f },
    // Gemmi✔️✔️:     {"BE7", RI::BUF, 0, ' ',   5, 357.156f },
    // Gemmi✔️✔️:     {"MRD", RI::BUF, 0, ' ',  14, 118.174f },
    // Gemmi✔️✔️:     {"MHA", RI::BUF, 0, ' ',  10, 190.154f },
    // Gemmi✔️✔️:     {"BU3", RI::BUF, 0, ' ',  10, 90.1210f },
    // Gemmi✔️✔️:     {"PGO", RI::BUF, 0, ' ',   8, 76.0944f },
    // Gemmi✔️✔️:     {"BU2", RI::BUF, 0, ' ',  10, 90.1210f },
    // Gemmi✔️✔️:     // 180
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"PDO", RI::BUF, 0, ' ',   8, 76.0944f },
    // Gemmi✔️✔️:     {"BU1", RI::BUF, 0, ' ',  10, 90.1210f },
    // Gemmi✔️✔️:     {"PG6", RI::BUF, 0, ' ',  26, 266.331f },
    // Gemmi✔️✔️:     {"1BO", RI::BUF, 0, ' ',  10, 74.1216f },
    // Gemmi✔️✔️:     {"PE7", RI::BUF, 0, ' ',  30, 342.449f },
    // Gemmi✔️✔️:     {"PG5", RI::BUF, 0, ' ',  18, 178.226f },
    // Gemmi✔️✔️:     {"TFP", RI::BUF, 0, ' ',  24, 407.496f },
    // Gemmi✔️✔️:     {"DHD", RI::BUF, 0, ' ',   4, 160.082f },
    // Gemmi✔️✔️:     {"PEU", RI::BUF, 0, ' ', 112, 1221.46f },
    // Gemmi✔️✔️:     {"TAU", RI::BUF, 0, ' ',   7, 125.147f },
    // Gemmi✔️✔️:     {"SBT", RI::BUF, 0, ' ',  10, 74.1216f },
    // Gemmi✔️✔️:     // 180
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"SAL", RI::BUF, 0, ' ',   6, 138.121f },
    // Gemmi✔️✔️:     {"IOH", RI::BUF, 0, ' ',   8, 60.0950f },
    // Gemmi✔️✔️:     {"IPA", RI::BUF, 0, ' ',   8, 60.0950f },
    // Gemmi✔️✔️:     {"PIG", RI::BUF, 0, ' ',  14, 150.173f },
    // Gemmi✔️✔️:     {"B3P", RI::BUF, 0, ' ',  26, 282.334f },
    // Gemmi✔️✔️:     {"BTB", RI::BUF, 0, ' ',  19, 209.240f },
    // Gemmi✔️✔️:     {"NHE", RI::BUF, 0, ' ',  17, 207.290f },
    // Gemmi✔️✔️:     {"C8E", RI::BUF, 0, ' ',  34, 306.438f },
    // Gemmi✔️✔️:     {"OTE", RI::BUF, 0, ' ',  34, 306.438f },
    // Gemmi✔️✔️:     {"PE4", RI::BUF, 0, ' ',  34, 354.436f },
    // Gemmi✔️✔️:     {"XPE", RI::BUF, 0, ' ',  42, 458.541f },
    // Gemmi✔️✔️:     // 200
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"PE8", RI::BUF, 0, ' ',  34, 370.436f },
    // Gemmi✔️✔️:     {"P33", RI::BUF, 0, ' ',  30, 326.383f },
    // Gemmi✔️✔️:     {"N8E", RI::BUF, 0, ' ',  38, 350.491f },
    // Gemmi✔️✔️:     {"2OS", RI::BUF, 0, ' ',  36, 468.493f },
    // Gemmi✔️✔️:     {"1PS", RI::BUF, 0, ' ',  11, 201.243f },
    // Gemmi✔️✔️:     {"CPS", RI::BUF, 0, ' ',  58, 614.877f },
    // Gemmi✔️✔️:     {"DMX", RI::BUF, 0, ' ',  19, 257.349f },
    // Gemmi✔️✔️:     {"MPO", RI::BUF, 0, ' ',  15, 209.263f },
    // Gemmi✔️✔️:     {"GCD", RI::PYR, 0, ' ',   8, 176.124f },
    // Gemmi✔️✔️:     {"DXG", RI::BUF, 0, ' ',   8, 192.124f },
    // Gemmi✔️✔️:     {"CM5", RI::BUF, 0, ' ',  42, 494.573f },
    // Gemmi✔️✔️:     // 210
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"ACA", RI::BUF, 1, ' ',  13, 131.173f }, // peptide linking
    // Gemmi✔️✔️:     {"ACN", RI::BUF, 0, ' ',   6, 58.0791f },
    // Gemmi✔️✔️:     {"CCN", RI::BUF, 0, ' ',   3, 41.0519f },
    // Gemmi✔️✔️:     {"GLC", RI::PYR, 0, ' ',  12, 180.156f },
    // Gemmi✔️✔️:     {"DR6", RI::BUF, 0, ' ', 142, 1527.90f },
    // Gemmi✔️✔️:     {"NH4", RI::BUF, 0, ' ',   4, 18.0385f },
    // Gemmi✔️✔️:     {"AZI", RI::BUF, 0, ' ',   0, 42.0201f },
    // Gemmi✔️✔️:     {"BNG", RI::PYR, 0, ' ',  30, 306.395f },
    // Gemmi✔️✔️:     {"BOG", RI::PYR, 0, ' ',  28, 292.369f },
    // Gemmi✔️✔️:     {"BGC", RI::PYR, 0, ' ',  12, 180.156f },
    // Gemmi✔️✔️:     {"BCN", RI::BUF, 0, ' ',  13, 163.172f },
    // Gemmi✔️✔️:     // 220
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"BRO", RI::BUF, 0, ' ',   0, 79.9040f },
    // Gemmi✔️✔️:     {"CAC", RI::BUF, 0, ' ',   6, 136.989f },
    // Gemmi✔️✔️:     {"CBX", RI::BUF, 0, ' ',   2, 46.0254f },
    // Gemmi✔️✔️:     {"ACY", RI::BUF, 0, ' ',   4, 60.0520f },
    // Gemmi✔️✔️:     {"CBM", RI::BUF, 0, ' ',   4, 60.0520f },
    // Gemmi✔️✔️:     {"CLO", RI::BUF, 0, ' ',   0, 35.4530f },
    // Gemmi✔️✔️:     {"3CO", RI::BUF, 0, ' ',   0, 58.9332f },
    // Gemmi✔️✔️:     {"NCO", RI::BUF, 0, ' ',  18, 161.116f },
    // Gemmi✔️✔️:     {"CU1", RI::BUF, 0, ' ',   0, 63.5460f },
    // Gemmi✔️✔️:     {"CYN", RI::BUF, 0, ' ',   0, 26.0174f },
    // Gemmi✔️✔️:     {"MA4", RI::BUF, 0, ' ',  44, 508.600f },
    // Gemmi✔️✔️:     // 230
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"TAR", RI::BUF, 0, ' ',   6, 150.087f },
    // Gemmi✔️✔️:     {"GLO", RI::BUF, 0, ' ',  12, 180.156f },  // d-saccharide
    // Gemmi✔️✔️:     {"MTL", RI::BUF, 0, ' ',  14, 182.172f },
    // Gemmi✔️✔️:     {"SOR", RI::BUF, 0, ' ',  14, 182.172f },
    // Gemmi✔️✔️:     {"DMU", RI::BUF, 0, ' ',  42, 482.562f },  // d-saccharide
    // Gemmi✔️✔️:     {"DDQ", RI::BUF, 0, ' ',  27, 201.349f },
    // Gemmi✔️✔️:     {"DMF", RI::BUF, 0, ' ',   7, 73.0938f },
    // Gemmi✔️✔️:     {"DIO", RI::BUF, 0, ' ',   8, 88.1051f },
    // Gemmi✔️✔️:     {"DOX", RI::BUF, 0, ' ',   8, 88.1051f },
    // Gemmi✔️✔️:     {"12P", RI::BUF, 0, ' ',  50, 546.646f },
    // Gemmi✔️✔️:     {"SDS", RI::BUF, 0, ' ',  26, 266.397f },
    // Gemmi✔️✔️:     // 240
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"LMT", RI::BUF, 0, ' ',  46, 510.615f },  // d-saccharide
    // Gemmi✔️✔️:     {"EOH", RI::BUF, 0, ' ',   6, 46.0684f },
    // Gemmi✔️✔️:     {"EEE", RI::BUF, 0, ' ',   8, 88.1051f },
    // Gemmi✔️✔️:     {"EGL", RI::BUF, 0, ' ',   6, 62.0678f },
    // Gemmi✔️✔️:     {"FLO", RI::BUF, 0, ' ',   0, 18.9984f },
    // Gemmi✔️✔️:     {"TRT", RI::BUF, 0, ' ',  36, 352.508f },
    // Gemmi✔️✔️:     {"FCY", RI::BUF, 0, ' ',   7, 121.158f },
    // Gemmi✔️✔️:     {"FRU", RI::BUF, 0, ' ',  12, 180.156f },  // saccharide
    // Gemmi✔️✔️:     {"GBL", RI::BUF, 0, ' ',   6, 86.0892f },
    // Gemmi✔️✔️:     {"GPX", RI::BUF, 0, ' ',  14, 505.165f },
    // Gemmi✔️✔️:     {"HTO", RI::BUF, 0, ' ',  16, 148.200f },
    // Gemmi✔️✔️:     // 250
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"HTG", RI::BUF, 0, ' ',  26, 294.408f },
    // Gemmi✔️✔️:     {"B7G", RI::BUF, 0, ' ',  26, 278.342f },
    // Gemmi✔️✔️:     {"C10", RI::BUF, 0, ' ',  46, 422.596f },
    // Gemmi✔️✔️:     {"16D", RI::BUF, 0, ' ',  16, 116.205f },
    // Gemmi✔️✔️:     {"HEZ", RI::BUF, 0, ' ',  14, 118.174f },
    // Gemmi✔️✔️:     {"IOD", RI::BUF, 0, ' ',   0, 126.904f },
    // Gemmi✔️✔️:     {"IDO", RI::BUF, 0, ' ',   0, 126.904f },
    // Gemmi✔️✔️:     {"ICI", RI::BUF, 0, ' ',   8, 192.124f },
    // Gemmi✔️✔️:     {"ICT", RI::BUF, 0, ' ',   8, 192.124f },
    // Gemmi✔️✔️:     {"TLA", RI::BUF, 0, ' ',   6, 150.087f },
    // Gemmi✔️✔️:     {"LAT", RI::BUF, 0, ' ',  22, 342.296f },  // saccharide
    // Gemmi✔️✔️:     // 260
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"LBT", RI::BUF, 0, ' ',  22, 342.296f },  // saccharide
    // Gemmi✔️✔️:     {"LDA", RI::BUF, 0, ' ',  31, 229.402f },
    // Gemmi✔️✔️:     {"MN3", RI::BUF, 0, ' ',   0, 54.9380f },
    // Gemmi✔️✔️:     {"MRY", RI::BUF, 0, ' ',  10, 122.120f },
    // Gemmi✔️✔️:     {"MOH", RI::BUF, 0, ' ',   4, 32.0419f },
    // Gemmi✔️✔️:     {"BEQ", RI::BUF, 0, ' ',  38, 342.517f },
    // Gemmi✔️✔️:     {"C15", RI::BUF, 0, ' ',  38, 336.554f },
    // Gemmi✔️✔️:     {"MG8", RI::BUF, 0, ' ',  31, 321.410f },
    // Gemmi✔️✔️:     {"POL", RI::BUF, 0, ' ',   8, 60.0950f },
    // Gemmi✔️✔️:     {"NO3", RI::BUF, 0, ' ',   0, 62.0049f },
    // Gemmi✔️✔️:     {"JEF", RI::BUF, 0, ' ',  63, 597.822f },
    // Gemmi✔️✔️:     // 270
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"P4C", RI::BUF, 0, ' ',  28, 324.367f },
    // Gemmi✔️✔️:     {"CE1", RI::BUF, 0, ' ',  58, 538.755f },
    // Gemmi✔️✔️:     {"DIA", RI::BUF, 0, ' ',  20, 144.258f },
    // Gemmi✔️✔️:     {"CXE", RI::BUF, 0, ' ',  42, 378.544f },
    // Gemmi✔️✔️:     {"IPH", RI::BUF, 0, ' ',   6, 94.1112f },
    // Gemmi✔️✔️:     {"PIN", RI::BUF, 0, ' ',  18, 302.368f },
    // Gemmi✔️✔️:     {"15P", RI::BUF, 0, ' ', 140, 1529.83f },
    // Gemmi✔️✔️:     {"CRY", RI::BUF, 0, ' ',   8, 92.0938f },
    // Gemmi✔️✔️:     {"PGR", RI::BUF, 0, ' ',   8, 76.0944f },
    // Gemmi✔️✔️:     {"PGQ", RI::BUF, 0, ' ',   8, 76.0944f },
    // Gemmi✔️✔️:     {"SPD", RI::BUF, 0, ' ',  19, 145.246f },
    // Gemmi✔️✔️:     // 270
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"SPK", RI::BUF, 0, ' ',  30, 206.372f },
    // Gemmi✔️✔️:     {"SPM", RI::BUF, 0, ' ',  26, 202.340f },
    // Gemmi✔️✔️:     {"SUC", RI::PYR, 0, ' ',  22, 342.296f },
    // Gemmi✔️✔️:     {"TBU", RI::BUF, 0, ' ',  10, 74.1216f },
    // Gemmi✔️✔️:     {"TMA", RI::BUF, 0, ' ',  12, 74.1448f },
    // Gemmi✔️✔️:     {"TEP", RI::BUF, 0, ' ',   8, 180.164f },
    // Gemmi✔️✔️:     {"SCN", RI::BUF, 0, ' ',   0, 58.0824f },
    // Gemmi✔️✔️:     {"TRE", RI::PYR, 0, ' ',  22, 342.296f },
    // Gemmi✔️✔️:     {"ETF", RI::BUF, 0, ' ',   3, 100.040f },
    // Gemmi✔️✔️:     {"144", RI::BUF, 0, ' ',  12, 122.143f },
    // Gemmi✔️✔️:     {"UMQ", RI::BUF, 0, ' ',  44, 496.589f },
    // Gemmi✔️✔️:     // 280
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"URE", RI::BUF, 0, ' ',   4, 60.0553f },
    // Gemmi✔️✔️:     {"YT3", RI::BUF, 0, ' ',   0, 88.9059f },
    // Gemmi✔️✔️:     {"ZN2", RI::BUF, 0, ' ',   0, 65.3800f },
    // Gemmi✔️✔️:     {"FE2", RI::BUF, 0, ' ',   0, 55.8450f },
    // Gemmi✔️✔️:     {"3NI", RI::BUF, 0, ' ',   0, 58.6934f },
    // Gemmi✔️✔️:     {"SIA", RI::PYR, 0, ' ',   0, 0.0f },
    // Gemmi✔️✔️:     {"XYP", RI::PYR, 0, ' ',   0, 0.0f },
    // Gemmi✔️✔️:     {"A2G", RI::PYR, 0, ' ',   0, 0.0f },
    // Gemmi✔️✔️:     {"GLA", RI::PYR, 0, ' ',   0, 0.0f },
    // Gemmi✔️✔️:     {"NDG", RI::PYR, 0, ' ',   0, 0.0f },
    // Gemmi✔️✔️:     // 290
    // Gemmi✔️✔️:     {"NGA", RI::PYR, 0, ' ',   0, 0.0f },
    // Gemmi✔️✔️:     {"A",   RI::RNA, 2, 'A',  14, 347.221f },
    // Gemmi✔️✔️:     {"C",   RI::RNA, 2, 'C',  14, 323.197f },
    // Gemmi✔️✔️:     {"G",   RI::RNA, 2, 'G',  14, 363.221f },
    // Gemmi✔️✔️:     {"I",   RI::RNA, 2, 'I',  13, 348.206f },
    // Gemmi✔️✔️:     {"U",   RI::RNA, 2, 'U',  13, 324.181f },
    // Gemmi✔️✔️:     {"N",   RI::RNA, 2, 'N',  11, 214.11f },
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"F",   RI::BUF, 0, ' ',   0, 18.9984f },
    // Gemmi✔️✔️:     {"K",   RI::BUF, 0, ' ',   0, 39.0983f },
    // Gemmi✔️✔️:     {"DA", RI::DNA, 2, 'A',  14, 331.222f },
    // Gemmi✔️✔️:     // 300
    // Gemmi✔️✔️:     {"DC", RI::DNA, 2, 'C',  14, 307.197f },
    // Gemmi✔️✔️:     {"DG", RI::DNA, 2, 'G',  14, 347.221f },
    // Gemmi✔️✔️:     {"DI", RI::DNA, 2, 'I',  13, 332.207f },
    // Gemmi✔️✔️:     {"DT", RI::DNA, 2, 'T',  15, 322.208f },
    // Gemmi✔️✔️:     {"DU", RI::DNA, 2, 'U',  13, 308.182f },
    // Gemmi✔️✔️:     {"DN", RI::DNA, 2, 'N',  14, 198.111f },  // unknown DNA
    // Gemmi✔️✔️:     {"AG",  RI::BUF, 0, ' ',   0, 107.868f },
    // Gemmi✔️✔️:     {"AL",  RI::BUF, 0, ' ',   0, 26.9815f },
    // Gemmi✔️✔️:     // 308
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"BA",  RI::BUF, 0, ' ',   0, 137.327f },
    // Gemmi✔️✔️:     {"BR",  RI::BUF, 0, ' ',   0, 79.9040f },
    // Gemmi✔️✔️:     {"CA",  RI::BUF, 0, ' ',   0, 40.0780f },
    // Gemmi✔️✔️:     {"CD",  RI::BUF, 0, ' ',   0, 112.411f },
    // Gemmi✔️✔️:     {"CL",  RI::BUF, 0, ' ',   0, 35.4530f },
    // Gemmi✔️✔️:     {"CM",  RI::BUF, 0, ' ',   4, 60.0520f },
    // Gemmi✔️✔️:     {"CN",  RI::BUF, 0, ' ',   0, 27.0253f },
    // Gemmi✔️✔️:     {"CO",  RI::BUF, 0, ' ',   0, 58.9332f },
    // Gemmi✔️✔️:     {"CS",  RI::BUF, 0, ' ',   0, 132.905f },
    // Gemmi✔️✔️:     {"CU",  RI::BUF, 0, ' ',   0, 63.5460f },
    // Gemmi✔️✔️:     {"FE",  RI::BUF, 0, ' ',   0, 55.8450f },
    // Gemmi✔️✔️:     // 318
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"HG",  RI::BUF, 0, ' ',   0, 200.590f },
    // Gemmi✔️✔️:     {"LI",  RI::BUF, 0, ' ',   0, 6.94100f },
    // Gemmi✔️✔️:     {"MG",  RI::BUF, 0, ' ',   0, 24.3050f },
    // Gemmi✔️✔️:     {"MN",  RI::BUF, 0, ' ',   0, 54.9380f },
    // Gemmi✔️✔️:     {"NA",  RI::BUF, 0, ' ',   0, 22.9898f },
    // Gemmi✔️✔️:     {"NI",  RI::BUF, 0, ' ',   0, 58.6934f },
    // Gemmi✔️✔️:     {"NO",  RI::ELS, 0, ' ',   0, 30.0061f },
    // Gemmi✔️✔️:     {"PB",  RI::BUF, 0, ' ',   0, 207.200f },
    // Gemmi✔️✔️:     {"RB",  RI::BUF, 0, ' ',   0, 85.4678f },
    // Gemmi✔️✔️:     {"SR",  RI::BUF, 0, ' ',   0, 87.6200f },
    // Gemmi✔️✔️:     {"Y1",  RI::BUF, 0, ' ',   0, 88.9059f },
    // Gemmi✔️✔️:     // 328
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:     {"ZN",  RI::BUF, 0, ' ',   0, 65.3800f },
    // Gemmi✔️✔️:     {"",    RI::UNKNOWN, 0, ' ', 0, 0.0f }
    // Gemmi✔️✔️:   };
    // Gemmi✔️✔️:   return array[idx];
    // Gemmi✔️✔️: }
    // Behavior review: the frozen 368-row source table is represented field-for-field;
    // the valid-index source precondition is retained by direct indexing.
    // Complexity review: one static-table index is O(1), with no allocation or scan.
    RESIDUE_INFO_TABLE[idx]
}

#[must_use]
pub fn residue_info_checked(idx: usize) -> Option<ResidueInfo> {
    RESIDUE_INFO_TABLE.get(idx).copied()
}

#[must_use]
pub fn find_residue_info_index(name: &str) -> usize {
    // Gemmi✔️❌:  size_t find_tabulated_residue_idx(const std::string& name) {
    // Gemmi✔️❌:   if (name.size() == 3) {
    // Gemmi✔️❌:
    // Gemmi✔️❌: #define ID(s) (((s)[0] << 16 | (s)[1] << 8 | (s)[2]) & ~0x202020)
    // Gemmi✔️❌:     //printf(">>> %x %x %x\n", ID(name.c_str()), ID("ALA"), ID("GLX"));
    // Gemmi✔️❌:     switch (ID(name.c_str())) {
    // Gemmi✔️❌:       case ID("ALA"): return 0;
    // Gemmi✔️❌:       case ID("ARG"): return 1;
    // Gemmi✔️❌:       case ID("ASN"): return 2;
    // Gemmi✔️❌:       case ID("ABA"): return 3;
    // Gemmi✔️❌:       case ID("ASP"): return 4;
    // Gemmi✔️❌:       case ID("ASX"): return 5;
    // Gemmi✔️❌:       case ID("CYS"): return 6;
    // Gemmi✔️❌:       case ID("CSH"): return 7;
    // Gemmi✔️❌:       case ID("GLN"): return 8;
    // Gemmi✔️❌:       case ID("GLU"): return 9;
    // Gemmi✔️❌:       case ID("GLX"): return 10;
    // Gemmi✔️❌:       case ID("GLY"): return 11;
    // Gemmi✔️❌:       case ID("HIS"): return 12;
    // Gemmi✔️❌:       case ID("ILE"): return 13;
    // Gemmi✔️❌:       case ID("LEU"): return 14;
    // Gemmi✔️❌:       case ID("LYS"): return 15;
    // Gemmi✔️❌:       case ID("MET"): return 16;
    // Gemmi✔️❌:       case ID("MSE"): return 17;
    // Gemmi✔️❌:       case ID("ORN"): return 18;
    // Gemmi✔️❌:       case ID("PHE"): return 19;
    // Gemmi✔️❌:       case ID("PRO"): return 20;
    // Gemmi✔️❌:       case ID("SER"): return 21;
    // Gemmi✔️❌:       case ID("THR"): return 22;
    // Gemmi✔️❌:       case ID("TRY"):
    // Gemmi✔️❌:       case ID("TRP"): return 23;
    // Gemmi✔️❌:       case ID("TYR"): return 24;
    // Gemmi✔️❌:       case ID("UNK"): return 25;
    // Gemmi✔️❌:       case ID("VAL"): return 26;
    // Gemmi✔️❌:       case ID("SEC"): return 27;
    // Gemmi✔️❌:       case ID("PYL"): return 28;
    // Gemmi✔️❌:       case ID("SEP"): return 29;
    // Gemmi✔️❌:       case ID("TPO"): return 30;
    // Gemmi✔️❌:       case ID("PCA"): return 31;
    // Gemmi✔️❌:       case ID("CSO"): return 32;
    // Gemmi✔️❌:       case ID("PTR"): return 33;
    // Gemmi✔️❌:       case ID("KCX"): return 34;
    // Gemmi✔️❌:       case ID("CSD"): return 35;
    // Gemmi✔️❌:       case ID("LLP"): return 36;
    // Gemmi✔️❌:       case ID("CME"): return 37;
    // Gemmi✔️❌:       case ID("MLY"): return 38;
    // Gemmi✔️❌:       case ID("DAL"): return 39;
    // Gemmi✔️❌:       case ID("TYS"): return 40;
    // Gemmi✔️❌:       case ID("OCS"): return 41;
    // Gemmi✔️❌:       case ID("M3L"): return 42;
    // Gemmi✔️❌:       case ID("FME"): return 43;
    // Gemmi✔️❌:       case ID("ALY"): return 44;
    // Gemmi✔️❌:       case ID("HYP"): return 45;
    // Gemmi✔️❌:       case ID("CAS"): return 46;
    // Gemmi✔️❌:       case ID("CRO"): return 47;
    // Gemmi✔️❌:       case ID("CSX"): return 48;
    // Gemmi✔️❌:       case ID("DPR"): return 49;
    // Gemmi✔️❌:       case ID("DGL"): return 50;
    // Gemmi✔️❌:       case ID("DVA"): return 51;
    // Gemmi✔️❌:       case ID("CSS"): return 52;
    // Gemmi✔️❌:       case ID("DPN"): return 53;
    // Gemmi✔️❌:       case ID("DSN"): return 54;
    // Gemmi✔️❌:       case ID("DLE"): return 55;
    // Gemmi✔️❌:       case ID("HIC"): return 56;
    // Gemmi✔️❌:       case ID("NLE"): return 57;
    // Gemmi✔️❌:       case ID("MVA"): return 58;
    // Gemmi✔️❌:       case ID("MLZ"): return 59;
    // Gemmi✔️❌:       case ID("CR2"): return 60;
    // Gemmi✔️❌:       case ID("SAR"): return 61;
    // Gemmi✔️❌:       case ID("DAR"): return 62;
    // Gemmi✔️❌:       case ID("DLY"): return 63;
    // Gemmi✔️❌:       case ID("YCM"): return 64;
    // Gemmi✔️❌:       case ID("NRQ"): return 65;
    // Gemmi✔️❌:       case ID("CGU"): return 66;
    // Gemmi✔️❌:       case ID("0TD"): return 67;
    // Gemmi✔️❌:       case ID("MLE"): return 68;
    // Gemmi✔️❌:       case ID("DAS"): return 69;
    // Gemmi✔️❌:       case ID("DTR"): return 70;
    // Gemmi✔️❌:       case ID("CXM"): return 71;
    // Gemmi✔️❌:       case ID("TPQ"): return 72;
    // Gemmi✔️❌:       case ID("DCY"): return 73;
    // Gemmi✔️❌:       case ID("DSG"): return 74;
    // Gemmi✔️❌:       case ID("DTY"): return 75;
    // Gemmi✔️❌:       case ID("DHI"): return 76;
    // Gemmi✔️❌:       case ID("MEN"): return 77;
    // Gemmi✔️❌:       case ID("DTH"): return 78;
    // Gemmi✔️❌:       case ID("SAC"): return 79;
    // Gemmi✔️❌:       case ID("DGN"): return 80;
    // Gemmi✔️❌:       case ID("AIB"): return 81;
    // Gemmi✔️❌:       case ID("SMC"): return 82;
    // Gemmi✔️❌:       case ID("IAS"): return 83;
    // Gemmi✔️❌:       case ID("CIR"): return 84;
    // Gemmi✔️❌:       case ID("BMT"): return 85;
    // Gemmi✔️❌:       case ID("DIL"): return 86;
    // Gemmi✔️❌:       case ID("FGA"): return 87;
    // Gemmi✔️❌:       case ID("PHI"): return 88;
    // Gemmi✔️❌:       case ID("CRQ"): return 89;
    // Gemmi✔️❌:       case ID("SME"): return 90;
    // Gemmi✔️❌:       case ID("GHP"): return 91;
    // Gemmi✔️❌:       case ID("MHO"): return 92;
    // Gemmi✔️❌:       case ID("NEP"): return 93;
    // Gemmi✔️❌:       case ID("TRQ"): return 94;
    // Gemmi✔️❌:       case ID("TOX"): return 95;
    // Gemmi✔️❌:       case ID("ALC"): return 96;
    // Gemmi✔️❌:       case ID("3FG"): return 97;
    // Gemmi✔️❌:       case ID("SCH"): return 98;
    // Gemmi✔️❌:       case ID("MDO"): return 99;
    // Gemmi✔️❌:       case ID("MAA"): return 100;
    // Gemmi✔️❌:       case ID("GYS"): return 101;
    // Gemmi✔️❌:       case ID("MK8"): return 102;
    // Gemmi✔️❌:       case ID("CR8"): return 103;
    // Gemmi✔️❌:       case ID("KPI"): return 104;
    // Gemmi✔️❌:       case ID("SCY"): return 105;
    // Gemmi✔️❌:       case ID("DHA"): return 106;
    // Gemmi✔️❌:       case ID("OMY"): return 107;
    // Gemmi✔️❌:       case ID("CAF"): return 108;
    // Gemmi✔️❌:       case ID("0AF"): return 109;
    // Gemmi✔️❌:       case ID("SNN"): return 110;
    // Gemmi✔️❌:       case ID("MHS"): return 111;
    // Gemmi✔️❌:       case ID("MLU"): return 112;
    // Gemmi✔️❌:       case ID("SNC"): return 113;
    // Gemmi✔️❌:       case ID("PHD"): return 114;
    // Gemmi✔️❌:       case ID("B3E"): return 115;
    // Gemmi✔️❌:       case ID("MEA"): return 116;
    // Gemmi✔️❌:       case ID("MED"): return 117;
    // Gemmi✔️❌:       case ID("OAS"): return 118;
    // Gemmi✔️❌:       case ID("GL3"): return 119;
    // Gemmi✔️❌:       case ID("FVA"): return 120;
    // Gemmi✔️❌:       case ID("PHL"): return 121;
    // Gemmi✔️❌:       case ID("CRF"): return 122;
    // Gemmi✔️❌:       case ID("OMZ"): return 123;
    // Gemmi✔️❌:       case ID("BFD"): return 124;
    // Gemmi✔️❌:       case ID("MEQ"): return 125;
    // Gemmi✔️❌:       case ID("DAB"): return 126;
    // Gemmi✔️❌:       case ID("AGM"): return 127;
    // Gemmi✔️❌:       case ID("PSU"): return 128;
    // Gemmi✔️❌:       case ID("5MU"): return 129;
    // Gemmi✔️❌:       case ID("7MG"): return 130;
    // Gemmi✔️❌:       case ID("OMG"): return 131;
    // Gemmi✔️❌:       case ID("UR3"): return 132;
    // Gemmi✔️❌:       case ID("OMC"): return 133;
    // Gemmi✔️❌:       case ID("2MG"): return 134;
    // Gemmi✔️❌:       case ID("H2U"): return 135;
    // Gemmi✔️❌:       case ID("4SU"): return 136;
    // Gemmi✔️❌:       case ID("OMU"): return 137;
    // Gemmi✔️❌:       case ID("4OC"): return 138;
    // Gemmi✔️❌:       case ID("MA6"): return 139;
    // Gemmi✔️❌:       case ID("M2G"): return 140;
    // Gemmi✔️❌:       case ID("1MA"): return 141;
    // Gemmi✔️❌:       case ID("6MZ"): return 142;
    // Gemmi✔️❌:       case ID("CCC"): return 143;
    // Gemmi✔️❌:       case ID("2MA"): return 144;
    // Gemmi✔️❌:       case ID("1MG"): return 145;
    // Gemmi✔️❌:       case ID("5BU"): return 146;
    // Gemmi✔️❌:       case ID("MIA"): return 147;
    // Gemmi✔️❌:       case ID("DOC"): return 148;
    // Gemmi✔️❌:       case ID("8OG"): return 149;
    // Gemmi✔️❌:       case ID("5CM"): return 150;
    // Gemmi✔️❌:       case ID("3DR"): return 151;
    // Gemmi✔️❌:       case ID("BRU"): return 152;
    // Gemmi✔️❌:       case ID("CBR"): return 153;
    // Gemmi✔️❌:       case ID("WAT"):
    // Gemmi✔️❌:       case ID("H2O"):
    // Gemmi✔️❌:       case ID("HOH"): return 154;
    // Gemmi✔️❌:       case ID("DOD"): return 155;
    // Gemmi✔️❌:       case ID("HEM"): return 156;
    // Gemmi✔️❌:       case ID("SO4"): return 157;
    // Gemmi✔️❌:       case ID("GOL"): return 158;
    // Gemmi✔️❌:       case ID("EDO"): return 159;
    // Gemmi✔️❌:       case ID("NAG"): return 160;
    // Gemmi✔️❌:       case ID("PO4"): return 161;
    // Gemmi✔️❌:       case ID("ACT"): return 162;
    // Gemmi✔️❌:       case ID("PEG"): return 163;
    // Gemmi✔️❌:       case ID("MAN"): return 164;
    // Gemmi✔️❌:       case ID("FAD"): return 165;
    // Gemmi✔️❌:       case ID("BMA"): return 166;
    // Gemmi✔️❌:       case ID("ADP"): return 167;
    // Gemmi✔️❌:       case ID("DMS"): return 168;
    // Gemmi✔️❌:       case ID("ACE"): return 169;
    // Gemmi✔️❌:       case ID("NH2"): return 170;
    // Gemmi✔️❌:       case ID("MPD"): return 171;
    // Gemmi✔️❌:       case ID("MES"): return 172;
    // Gemmi✔️❌:       case ID("NAD"): return 173;
    // Gemmi✔️❌:       case ID("NAP"): return 174;
    // Gemmi✔️❌:       case ID("TRS"): return 175;
    // Gemmi✔️❌:       case ID("ATP"): return 176;
    // Gemmi✔️❌:       case ID("PG4"): return 177;
    // Gemmi✔️❌:       case ID("GDP"): return 178;
    // Gemmi✔️❌:       case ID("FUC"): return 179;
    // Gemmi✔️❌:       case ID("FMT"): return 180;
    // Gemmi✔️❌:       case ID("GAL"): return 181;
    // Gemmi✔️❌:       case ID("PGE"): return 182;
    // Gemmi✔️❌:       case ID("FMN"): return 183;
    // Gemmi✔️❌:       case ID("PLP"): return 184;
    // Gemmi✔️❌:       case ID("EPE"): return 185;
    // Gemmi✔️❌:       case ID("SF4"): return 186;
    // Gemmi✔️❌:       case ID("BME"): return 187;
    // Gemmi✔️❌:       case ID("CIT"): return 188;
    // Gemmi✔️❌:       case ID("BE7"): return 189;
    // Gemmi✔️❌:       case ID("MRD"): return 190;
    // Gemmi✔️❌:       case ID("MHA"): return 191;
    // Gemmi✔️❌:       case ID("BU3"): return 192;
    // Gemmi✔️❌:       case ID("PGO"): return 193;
    // Gemmi✔️❌:       case ID("BU2"): return 194;
    // Gemmi✔️❌:       case ID("PDO"): return 195;
    // Gemmi✔️❌:       case ID("BU1"): return 196;
    // Gemmi✔️❌:       case ID("PG6"): return 197;
    // Gemmi✔️❌:       case ID("1BO"): return 198;
    // Gemmi✔️❌:       case ID("PE7"): return 199;
    // Gemmi✔️❌:       case ID("PG5"): return 200;
    // Gemmi✔️❌:       case ID("TFP"): return 201;
    // Gemmi✔️❌:       case ID("DHD"): return 202;
    // Gemmi✔️❌:       case ID("PEU"): return 203;
    // Gemmi✔️❌:       case ID("TAU"): return 204;
    // Gemmi✔️❌:       case ID("SBT"): return 205;
    // Gemmi✔️❌:       case ID("SAL"): return 206;
    // Gemmi✔️❌:       case ID("IOH"): return 207;
    // Gemmi✔️❌:       case ID("IPA"): return 208;
    // Gemmi✔️❌:       case ID("PIG"): return 209;
    // Gemmi✔️❌:       case ID("B3P"): return 210;
    // Gemmi✔️❌:       case ID("BTB"): return 211;
    // Gemmi✔️❌:       case ID("NHE"): return 212;
    // Gemmi✔️❌:       case ID("C8E"): return 213;
    // Gemmi✔️❌:       case ID("OTE"): return 214;
    // Gemmi✔️❌:       case ID("PE4"): return 215;
    // Gemmi✔️❌:       case ID("XPE"): return 216;
    // Gemmi✔️❌:       case ID("PE8"): return 217;
    // Gemmi✔️❌:       case ID("P33"): return 218;
    // Gemmi✔️❌:       case ID("N8E"): return 219;
    // Gemmi✔️❌:       case ID("2OS"): return 220;
    // Gemmi✔️❌:       case ID("1PS"): return 221;
    // Gemmi✔️❌:       case ID("CPS"): return 222;
    // Gemmi✔️❌:       case ID("DMX"): return 223;
    // Gemmi✔️❌:       case ID("MPO"): return 224;
    // Gemmi✔️❌:       case ID("GCD"): return 225;
    // Gemmi✔️❌:       case ID("DXG"): return 226;
    // Gemmi✔️❌:       case ID("CM5"): return 227;
    // Gemmi✔️❌:       case ID("ACA"): return 228;
    // Gemmi✔️❌:       case ID("ACN"): return 229;
    // Gemmi✔️❌:       case ID("CCN"): return 230;
    // Gemmi✔️❌:       case ID("GLC"): return 231;
    // Gemmi✔️❌:       case ID("DR6"): return 232;
    // Gemmi✔️❌:       case ID("NH4"): return 233;
    // Gemmi✔️❌:       case ID("AZI"): return 234;
    // Gemmi✔️❌:       case ID("BNG"): return 235;
    // Gemmi✔️❌:       case ID("BOG"): return 236;
    // Gemmi✔️❌:       case ID("BGC"): return 237;
    // Gemmi✔️❌:       case ID("BCN"): return 238;
    // Gemmi✔️❌:       case ID("BRO"): return 239;
    // Gemmi✔️❌:       case ID("CAC"): return 240;
    // Gemmi✔️❌:       case ID("CBX"): return 241;
    // Gemmi✔️❌:       case ID("ACY"): return 242;
    // Gemmi✔️❌:       case ID("CBM"): return 243;
    // Gemmi✔️❌:       case ID("CLO"): return 244;
    // Gemmi✔️❌:       case ID("3CO"): return 245;
    // Gemmi✔️❌:       case ID("NCO"): return 246;
    // Gemmi✔️❌:       case ID("CU1"): return 247;
    // Gemmi✔️❌:       case ID("CYN"): return 248;
    // Gemmi✔️❌:       case ID("MA4"): return 249;
    // Gemmi✔️❌:       case ID("TAR"): return 250;
    // Gemmi✔️❌:       case ID("GLO"): return 251;
    // Gemmi✔️❌:       case ID("MTL"): return 252;
    // Gemmi✔️❌:       case ID("SOR"): return 253;
    // Gemmi✔️❌:       case ID("DMU"): return 254;
    // Gemmi✔️❌:       case ID("DDQ"): return 255;
    // Gemmi✔️❌:       case ID("DMF"): return 256;
    // Gemmi✔️❌:       case ID("DIO"): return 257;
    // Gemmi✔️❌:       case ID("DOX"): return 258;
    // Gemmi✔️❌:       case ID("12P"): return 259;
    // Gemmi✔️❌:       case ID("SDS"): return 260;
    // Gemmi✔️❌:       case ID("LMT"): return 261;
    // Gemmi✔️❌:       case ID("EOH"): return 262;
    // Gemmi✔️❌:       case ID("EEE"): return 263;
    // Gemmi✔️❌:       case ID("EGL"): return 264;
    // Gemmi✔️❌:       case ID("FLO"): return 265;
    // Gemmi✔️❌:       case ID("TRT"): return 266;
    // Gemmi✔️❌:       case ID("FCY"): return 267;
    // Gemmi✔️❌:       case ID("FRU"): return 268;
    // Gemmi✔️❌:       case ID("GBL"): return 269;
    // Gemmi✔️❌:       case ID("GPX"): return 270;
    // Gemmi✔️❌:       case ID("HTO"): return 271;
    // Gemmi✔️❌:       case ID("HTG"): return 272;
    // Gemmi✔️❌:       case ID("B7G"): return 273;
    // Gemmi✔️❌:       case ID("C10"): return 274;
    // Gemmi✔️❌:       case ID("16D"): return 275;
    // Gemmi✔️❌:       case ID("HEZ"): return 276;
    // Gemmi✔️❌:       case ID("IOD"): return 277;
    // Gemmi✔️❌:       case ID("IDO"): return 278;
    // Gemmi✔️❌:       case ID("ICI"): return 279;
    // Gemmi✔️❌:       case ID("ICT"): return 280;
    // Gemmi✔️❌:       case ID("TLA"): return 281;
    // Gemmi✔️❌:       case ID("LAT"): return 282;
    // Gemmi✔️❌:       case ID("LBT"): return 283;
    // Gemmi✔️❌:       case ID("LDA"): return 284;
    // Gemmi✔️❌:       case ID("MN3"): return 285;
    // Gemmi✔️❌:       case ID("MRY"): return 286;
    // Gemmi✔️❌:       case ID("MOH"): return 287;
    // Gemmi✔️❌:       case ID("BEQ"): return 288;
    // Gemmi✔️❌:       case ID("C15"): return 289;
    // Gemmi✔️❌:       case ID("MG8"): return 290;
    // Gemmi✔️❌:       case ID("POL"): return 291;
    // Gemmi✔️❌:       case ID("NO3"): return 292;
    // Gemmi✔️❌:       case ID("JEF"): return 293;
    // Gemmi✔️❌:       case ID("P4C"): return 294;
    // Gemmi✔️❌:       case ID("CE1"): return 295;
    // Gemmi✔️❌:       case ID("DIA"): return 296;
    // Gemmi✔️❌:       case ID("CXE"): return 297;
    // Gemmi✔️❌:       case ID("IPH"): return 298;
    // Gemmi✔️❌:       case ID("PIN"): return 299;
    // Gemmi✔️❌:       case ID("15P"): return 300;
    // Gemmi✔️❌:       case ID("CRY"): return 301;
    // Gemmi✔️❌:       case ID("PGR"): return 302;
    // Gemmi✔️❌:       case ID("PGQ"): return 303;
    // Gemmi✔️❌:       case ID("SPD"): return 304;
    // Gemmi✔️❌:       case ID("SPK"): return 305;
    // Gemmi✔️❌:       case ID("SPM"): return 306;
    // Gemmi✔️❌:       case ID("SUC"): return 307;
    // Gemmi✔️❌:       case ID("TBU"): return 308;
    // Gemmi✔️❌:       case ID("TMA"): return 309;
    // Gemmi✔️❌:       case ID("TEP"): return 310;
    // Gemmi✔️❌:       case ID("SCN"): return 311;
    // Gemmi✔️❌:       case ID("TRE"): return 312;
    // Gemmi✔️❌:       case ID("ETF"): return 313;
    // Gemmi✔️❌:       case ID("144"): return 314;
    // Gemmi✔️❌:       case ID("UMQ"): return 315;
    // Gemmi✔️❌:       case ID("URE"): return 316;
    // Gemmi✔️❌:       case ID("YT3"): return 317;
    // Gemmi✔️❌:       case ID("ZN2"): return 318;
    // Gemmi✔️❌:       case ID("FE2"): return 319;
    // Gemmi✔️❌:       case ID("3NI"): return 320;
    // Gemmi✔️❌:       case ID("SIA"): return 321;
    // Gemmi✔️❌:       case ID("XYP"): return 322;
    // Gemmi✔️❌:       case ID("A2G"): return 323;
    // Gemmi✔️❌:       case ID("GLA"): return 324;
    // Gemmi✔️❌:       case ID("NDG"): return 325;
    // Gemmi✔️❌:       case ID("NGA"): return 326;
    // Gemmi✔️❌: #undef ID
    // Gemmi✔️❌:     }} else if (name.size() == 1) {
    // Gemmi✔️❌:     switch (name[0]& ~0x20) {
    // Gemmi✔️❌:       case 'A': return 327;
    // Gemmi✔️❌:       case 'C': return 328;
    // Gemmi✔️❌:       case 'G': return 329;
    // Gemmi✔️❌:       case 'I': return 330;
    // Gemmi✔️❌:       case 'U': return 331;
    // Gemmi✔️❌:       case 'N': return 332;
    // Gemmi✔️❌:       case 'F': return 333;
    // Gemmi✔️❌:       case 'K': return 334;
    // Gemmi✔️❌:     }
    // Gemmi✔️❌:   } else if (name.size() == 2) {
    // Gemmi✔️❌:     if (name[0] == 'D' || name[0] == '+')
    // Gemmi✔️❌:       switch (name[1]) {
    // Gemmi✔️❌:         case 'A': return 335;
    // Gemmi✔️❌:         case 'C': return 336;
    // Gemmi✔️❌:         case 'G': return 337;
    // Gemmi✔️❌:         case 'I': return 338;
    // Gemmi✔️❌:         case 'T': return 339;
    // Gemmi✔️❌:         case 'U': return 340;
    // Gemmi✔️❌:         case 'N': return 341;
    // Gemmi✔️❌:       }
    // Gemmi✔️❌: #define ID(s) ((s)[0] << 8 | (s)[1])
    // Gemmi✔️❌:     switch (ID(name.c_str())) {
    // Gemmi✔️❌:         case ID("AG"): return 342;
    // Gemmi✔️❌:         case ID("AL"): return 343;
    // Gemmi✔️❌:         case ID("BA"): return 344;
    // Gemmi✔️❌:         case ID("BR"): return 345;
    // Gemmi✔️❌:         case ID("CA"): return 346;
    // Gemmi✔️❌:         case ID("CD"): return 347;
    // Gemmi✔️❌:         case ID("CL"): return 348;
    // Gemmi✔️❌:         case ID("CM"): return 349;
    // Gemmi✔️❌:         case ID("CN"): return 350;
    // Gemmi✔️❌:         case ID("CO"): return 351;
    // Gemmi✔️❌:         case ID("CS"): return 352;
    // Gemmi✔️❌:         case ID("CU"): return 353;
    // Gemmi✔️❌:         case ID("FE"): return 354;
    // Gemmi✔️❌:         case ID("HG"): return 355;
    // Gemmi✔️❌:         case ID("LI"): return 356;
    // Gemmi✔️❌:         case ID("MG"): return 357;
    // Gemmi✔️❌:         case ID("MN"): return 358;
    // Gemmi✔️❌:         case ID("NA"): return 359;
    // Gemmi✔️❌:         case ID("NI"): return 360;
    // Gemmi✔️❌:         case ID("NO"): return 361;
    // Gemmi✔️❌:         case ID("PB"): return 362;
    // Gemmi✔️❌:         case ID("RB"): return 363;
    // Gemmi✔️❌:         case ID("SR"): return 364;
    // Gemmi✔️❌:         case ID("Y1"): return 365;
    // Gemmi✔️❌:         case ID("ZN"): return 366;
    // Gemmi✔️❌:       }
    // Gemmi✔️❌: #undef ID
    // Gemmi✔️❌:     }
    // Gemmi✔️❌:     return 367;
    // Gemmi✔️❌:  }
    // Behavior review: every source length branch, alias, byte-case rule and fallback
    // is retained below for the modeled ASCII residue-name input space.
    // Complexity review: the three-byte branch currently allocates one tiny uppercase
    // String whereas the source masks bytes in place; asymptotics remain O(1), but this
    // is a real avoidable allocation and therefore not marked performance-equivalent.
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
            _ => UNKNOWN_TABULATED_RESIDUE_INDEX,
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
            _ => UNKNOWN_TABULATED_RESIDUE_INDEX,
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
                _ => UNKNOWN_TABULATED_RESIDUE_INDEX,
            }
        }
        _ => UNKNOWN_TABULATED_RESIDUE_INDEX,
    }
    // END GEMMI CPP FUNCTION gemmi::find_residue_info_index
}

#[must_use]
pub fn find_residue_info(name: &str) -> ResidueInfo {
    // Gemmi✔️✔️: ResidueInfo& find_tabulated_residue(const std::string& name) {
    // Gemmi✔️✔️:   size_t idx = find_tabulated_residue_idx(name);
    // Gemmi✔️✔️:   return get_residue_info(idx);
    // Gemmi✔️✔️: }
    // Behavior and complexity review: one source lookup followed by one O(1) table index.
    residue_info(find_residue_info_index(name))
}

#[must_use]
pub fn residue_code(name: &str) -> ResidueCode {
    find_residue_info(name).code
}

#[must_use]
pub fn expand_one_letter(code: char, kind: ResidueInfoKind) -> Option<&'static str> {
    // Gemmi✔️✔️: inline const char* expand_one_letter(char c, ResidueKind kind) {
    // Gemmi✔️✔️:   static const char* names =
    // Gemmi✔️✔️:     // amino-acids (all letters but J are used)
    // Gemmi✔️✔️:     "ALA\0ASX\0CYS\0ASP\0GLU\0PHE\0GLY\0HIS\0ILE\0\0   LYS\0LEU\0MET\0"  // A-M
    // Gemmi✔️✔️:     "ASN\0PYL\0PRO\0GLN\0ARG\0SER\0THR\0SEC\0VAL\0TRP\0UNK\0TYR\0GLX\0"  // N-Z
    // Gemmi✔️✔️:     // DNA
    // Gemmi✔️✔️:     "DA\0 \0\0  DC\0 \0\0  \0\0  \0\0  DG\0 \0\0  DI\0 \0\0  \0\0  \0\0  \0\0  "   // A-M
    // Gemmi✔️✔️:     "DN\0 \0\0  \0\0  \0\0  \0\0  \0\0  DT\0 DU\0 \0\0  \0\0  \0\0  \0\0  \0\0  "; // N-Z
    // Gemmi✔️✔️:   c &= ~0x20;
    // Gemmi✔️✔️:   const char* ret = nullptr;
    // Gemmi✔️✔️:   if (c >= 'A' && c <= 'Z') {
    // Gemmi✔️✔️:     ret = &names[4 * (c - 'A')];
    // Gemmi✔️✔️:     if (kind == ResidueKind::AA) {
    // Gemmi✔️✔️:       // ret is already set
    // Gemmi✔️✔️:     } else if (kind == ResidueKind::DNA) {
    // Gemmi✔️✔️:       ret += 4 * 26;
    // Gemmi✔️✔️:     } else if (kind == ResidueKind::RNA && c != 'T') {
    // Gemmi✔️✔️:       ret += 4 * 26 + 1;
    // Gemmi✔️✔️:     } else {
    // Gemmi✔️✔️:       ret = nullptr;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return (ret && *ret) ? ret : nullptr;
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️:
    // Gemmi✔️✔️: /// kind can be AA, RNA or DNA
    // Gemmi✔️✔️: GEMMI_DLL std::vector<std::string> expand_one_letter_sequence(const std::string& seq,
    // Gemmi✔️✔️:                                                               ResidueKind kind);
    // Gemmi✔️✔️:
    // Behavior review: the same AA/DNA/RNA offsets are projected as explicit matches;
    // non-ASCII Rust scalar values are outside the source char domain and return None.
    // Complexity review: both forms use constant-time bounded dispatch and no allocation.
    let c = code.to_ascii_uppercase();
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResidueSequenceError {
    UnmatchedParenthesis,
    UnexpectedLetter {
        kind: &'static str,
        letter: char,
        source_code: i16,
    },
}

impl std::fmt::Display for ResidueSequenceError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnmatchedParenthesis => formatter.write_str("unmatched '(' in sequence"),
            Self::UnexpectedLetter {
                kind,
                letter,
                source_code,
            } => write!(
                formatter,
                "unexpected letter in {kind} sequence: {letter} ({source_code})"
            ),
        }
    }
}

impl std::error::Error for ResidueSequenceError {}

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
    // Gemmi✔️✔️: // equivalent of std::isspace for C locale (no handling of EOF)
    // Gemmi✔️✔️: inline bool is_space(char c) {
    // Gemmi✔️✔️:   static const std::uint8_t table[256] = { // 1 for 9-13 and 32
    // Gemmi✔️✔️:     0,0,0,0,0,0,0,0, 0,1,1,1,1,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi✔️✔️:     1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi✔️✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi✔️✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi✔️✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi✔️✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi✔️✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi✔️✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0
    // Gemmi✔️✔️:   };
    // Gemmi✔️✔️:   return table[(std::uint8_t)c] != 0;
    // Gemmi✔️✔️: }
    // Behavior and complexity review: the six C-locale whitespace bytes are tested in
    // constant time without locale state or allocation.
    matches!(byte, b'\t' | b'\n' | 0x0b | 0x0c | b'\r' | b' ')
}

pub fn expand_one_letter_sequence(
    seq: &str,
    kind: ResidueInfoKind,
) -> Result<Vec<String>, ResidueSequenceError> {
    // Gemmi✔️✔️: std::vector<std::string> expand_one_letter_sequence(const std::string& seq,
    // Gemmi✔️✔️:                                                     ResidueKind kind) {
    // Gemmi✔️✔️:
    // Gemmi✔️✔️:  std::vector<std::string> r;
    // Gemmi✔️✔️:   r.reserve(seq.size());
    // Gemmi✔️✔️:   auto kind_str = [&]() {
    // Gemmi✔️✔️:     switch (kind) {
    // Gemmi✔️✔️:       case ResidueKind::AA: return "peptide";
    // Gemmi✔️✔️:       case ResidueKind::RNA: return "RNA";
    // Gemmi✔️✔️:       case ResidueKind::DNA: return "DNA";
    // Gemmi✔️✔️:       default: return "unknown";
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   };
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
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return r;
    // Gemmi✔️✔️: }
    // Behavior review: byte iteration, C-locale whitespace, parenthesized tokens,
    // first-error ordering and source error payloads are retained for UTF-8 Rust input.
    // Complexity review: one linear scan and one output allocation match the source shape.
    let bytes = seq.as_bytes();
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
                return Err(ResidueSequenceError::UnmatchedParenthesis);
            };
            let end = start + offset;
            residues.push(seq[start..end].to_string());
            i = end + 1;
        } else if let Some(name) = expand_one_letter(c as char, kind) {
            residues.push(name.to_string());
            i += 1;
        } else {
            return Err(ResidueSequenceError::UnexpectedLetter {
                kind: residue_sequence_kind_str(kind),
                letter: c as char,
                // The pinned Linux C++ reference uses signed `char`; preserve
                // its integer promotion while retaining the original byte in
                // `letter` through the one-byte Unicode scalar projection.
                source_code: i16::from(c as i8),
            });
        }
    }
    Ok(residues)
    // END GEMMI CPP FUNCTION gemmi::expand_one_letter_sequence
}

pub fn expand_protein_one_letter_string(seq: &str) -> Result<Vec<String>, ResidueSequenceError> {
    // Gemmi✔️✔️: std::vector<std::string> expand_protein_one_letter_string(const std::string& s);
    expand_one_letter_sequence(seq, ResidueInfoKind::Aa)
}
