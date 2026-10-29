//! Thin residue vocabulary projections through the canonical public BIO owner.
use crate::alignment_values::set;
use crate::host_values::usize_value;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[allow(non_camel_case_types)]
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
#[wasm_bindgen]
pub struct ResidueInfoKind {
    pub(crate) inner: ck::ResidueInfoKind,
}
#[wasm_bindgen]
impl ResidueInfoKind {
    pub fn name(&self) -> String {
        self.inner.name().into()
    }
    #[wasm_bindgen(getter)]
    pub fn value(&self) -> u8 {
        self.inner as u8
    }
    #[wasm_bindgen(getter,js_name=UNKNOWN)]
    pub fn unknown() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Unknown,
        }
    }
    #[wasm_bindgen(getter,js_name=AA)]
    pub fn aa() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Aa,
        }
    }
    #[wasm_bindgen(getter,js_name=AAD)]
    pub fn aad() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Aad,
        }
    }
    #[wasm_bindgen(getter,js_name=PAA)]
    pub fn paa() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Paa,
        }
    }
    #[wasm_bindgen(getter,js_name=MAA)]
    pub fn maa() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Maa,
        }
    }
    #[wasm_bindgen(getter,js_name=RNA)]
    pub fn rna() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Rna,
        }
    }
    #[wasm_bindgen(getter,js_name=DNA)]
    pub fn dna() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Dna,
        }
    }
    #[wasm_bindgen(getter,js_name=BUF)]
    pub fn buf() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Buf,
        }
    }
    #[wasm_bindgen(getter,js_name=HOH)]
    pub fn hoh() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Hoh,
        }
    }
    #[wasm_bindgen(getter,js_name=PYR)]
    pub fn pyr() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Pyr,
        }
    }
    #[wasm_bindgen(getter,js_name=KET)]
    pub fn ket() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Ket,
        }
    }
    #[wasm_bindgen(getter,js_name=ELS)]
    pub fn els() -> Self {
        Self {
            inner: ck::ResidueInfoKind::Els,
        }
    }
}
#[wasm_bindgen]
pub struct ResidueInfo {
    pub(crate) inner: ck::ResidueInfo,
}
#[wasm_bindgen]
impl ResidueInfo {
    #[wasm_bindgen(unchecked_return_type = "ResidueCode")]
    pub fn code(&self) -> u16 {
        self.inner.code.as_u16()
    }
    pub fn name(&self) -> String {
        self.inner.name.into()
    }
    pub fn kind(&self) -> ResidueInfoKind {
        ResidueInfoKind {
            inner: self.inner.kind,
        }
    }
    #[wasm_bindgen(js_name=kindName)]
    pub fn kind_name(&self) -> String {
        self.inner.kind.name().into()
    }
    #[wasm_bindgen(js_name=linkingType)]
    pub fn linking_type(&self) -> u8 {
        self.inner.linking_type
    }
    #[wasm_bindgen(js_name=oneLetterCode)]
    pub fn one_letter_code(&self) -> String {
        self.inner.one_letter_code.to_string()
    }
    #[wasm_bindgen(js_name=hydrogenCount)]
    pub fn hydrogen_count(&self) -> u8 {
        self.inner.hydrogen_count
    }
    pub fn weight(&self) -> f32 {
        self.inner.weight
    }
    #[wasm_bindgen(js_name=found)]
    pub fn found(&self) -> bool {
        self.inner.found()
    }
    #[wasm_bindgen(js_name=isWater)]
    pub fn is_water(&self) -> bool {
        self.inner.is_water()
    }
    #[wasm_bindgen(js_name=isDna)]
    pub fn is_dna(&self) -> bool {
        self.inner.is_dna()
    }
    #[wasm_bindgen(js_name=isRna)]
    pub fn is_rna(&self) -> bool {
        self.inner.is_rna()
    }
    #[wasm_bindgen(js_name=isNucleicAcid)]
    pub fn is_nucleic_acid(&self) -> bool {
        self.inner.is_nucleic_acid()
    }
    #[wasm_bindgen(js_name=isAminoAcid)]
    pub fn is_amino_acid(&self) -> bool {
        self.inner.is_amino_acid()
    }
    #[wasm_bindgen(js_name=isBufferOrWater)]
    pub fn is_buffer_or_water(&self) -> bool {
        self.inner.is_buffer_or_water()
    }
    #[wasm_bindgen(js_name=isStandard)]
    pub fn is_standard(&self) -> bool {
        self.inner.is_standard()
    }
    #[wasm_bindgen(js_name=isModifiedAminoAcid)]
    pub fn is_modified_amino_acid(&self) -> bool {
        self.inner.is_modified_amino_acid()
    }
    #[wasm_bindgen(js_name=isPeptideLinking)]
    pub fn is_peptide_linking(&self) -> bool {
        self.inner.is_peptide_linking()
    }
    #[wasm_bindgen(js_name=isNaLinking)]
    pub fn is_na_linking(&self) -> bool {
        self.inner.is_na_linking()
    }
    #[wasm_bindgen(js_name=fastaCode)]
    pub fn fasta_code(&self) -> String {
        self.inner.fasta_code().to_string()
    }
    #[wasm_bindgen(js_name=canonicalOneLetterCode,unchecked_return_type="string | null")]
    pub fn canonical_one_letter_code(&self) -> JsValue {
        self.inner
            .canonical_one_letter_code()
            .map_or(JsValue::NULL, |c| c.to_string().into())
    }
    #[wasm_bindgen(js_name=parentStandardCode,unchecked_return_type="ResidueCode | null")]
    pub fn parent_standard_code(&self) -> JsValue {
        self.inner
            .parent_standard_code()
            .map_or(JsValue::NULL, |c| c.as_u16().into())
    }
}
#[wasm_bindgen]
pub struct ResidueIdentity {
    inner: ck::ResidueIdentity,
}
#[wasm_bindgen]
impl ResidueIdentity {
    #[wasm_bindgen(constructor)]
    pub fn construct(name: String) -> Self {
        Self::new(name)
    }
    pub fn new(name: String) -> Self {
        // COSMolKit❗✔️: canonical_bio_residue.rs::ResidueIdentity::new:
        //             inner: ck::ResidueIdentity::new(name),
        // Transport one owned string; all lookup rules remain in the sole BIO owner.
        Self {
            inner: ck::ResidueIdentity::new(name),
        }
    }
    pub fn name(&self) -> String {
        self.inner.name().into()
    }
    #[wasm_bindgen(unchecked_return_type = "ResidueCode")]
    pub fn code(&self) -> u16 {
        self.inner.code().as_u16()
    }
    pub fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: self.inner.info(),
        }
    }
    #[wasm_bindgen(js_name=isTabulated)]
    pub fn is_tabulated(&self) -> bool {
        self.inner.is_tabulated()
    }
}
#[wasm_bindgen(js_name=findResidueInfo)]
pub fn find_residue_info(name: &str) -> ResidueInfo {
    ResidueInfo {
        inner: ck::find_residue_info(name),
    }
}
#[wasm_bindgen(js_name=findResidueInfoIndex)]
pub fn find_residue_info_index(name: &str) -> usize {
    ck::find_residue_info_index(name)
}
#[wasm_bindgen(js_name=residueInfo)]
pub fn residue_info(
    #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
) -> Result<ResidueInfo, JsValue> {
    let index = usize_value(&index, "index")?;
    ck::residue_info_checked(index)
        .map(|inner| ResidueInfo { inner })
        .ok_or_else(|| {
            js_sys::RangeError::new(&format!("residue info index {index} out of range")).into()
        })
}
#[wasm_bindgen(js_name=residueInfoChecked,unchecked_return_type="ResidueInfo | null")]
pub fn residue_info_checked(
    #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
) -> Result<JsValue, JsValue> {
    Ok(ck::residue_info_checked(usize_value(&index, "index")?)
        .map_or(JsValue::NULL, |inner| ResidueInfo { inner }.into()))
}
#[wasm_bindgen(js_name=residueCode,unchecked_return_type="ResidueCode")]
pub fn residue_code(name: &str) -> u16 {
    ck::residue_code(name).as_u16()
}
#[wasm_bindgen(js_name=expandOneLetter,unchecked_return_type="string | null")]
pub fn expand_one_letter(code: &str, kind: &ResidueInfoKind) -> Result<JsValue, JsValue> {
    let mut chars = code.chars();
    let c = chars
        .next()
        .filter(|_| chars.next().is_none())
        .ok_or_else(|| js_sys::RangeError::new("code must contain exactly one character"))?;
    Ok(ck::expand_one_letter(c, kind.inner).map_or(JsValue::NULL, JsValue::from_str))
}
#[wasm_bindgen(js_name=expandOneLetterSequence,unchecked_return_type="string[]")]
pub fn expand_one_letter_sequence(
    seq: &str,
    kind: &ResidueInfoKind,
) -> Result<js_sys::Array, JsValue> {
    ck::expand_one_letter_sequence(seq, kind.inner)
        .map(|values| values.into_iter().map(JsValue::from).collect())
        .map_err(|e| sequence_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen]
pub struct ResidueSequenceError {
    kind: String,
    message: String,
}
#[wasm_bindgen]
impl ResidueSequenceError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
}
pub(crate) fn sequence_error(source: &ck::ResidueSequenceError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::ResidueSequenceError::UnmatchedParenthesis => "UnmatchedParenthesis",
        ck::ResidueSequenceError::UnexpectedLetter { .. } => "UnexpectedLetter",
    };
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("ResidueSequenceError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    set(
        &e,
        "detail",
        ResidueSequenceError {
            kind: kind.into(),
            message: source.to_string(),
        }
        .into(),
    )?;
    if let ck::ResidueSequenceError::UnexpectedLetter {
        kind,
        letter,
        source_code,
    } = source
    {
        set(&e, "sequenceKind", (*kind).into())?;
        set(&e, "letter", letter.to_string().into())?;
        set(&e, "sourceCode", (*source_code).into())?;
    }
    Ok(e)
}
#[wasm_bindgen]
pub struct ResidueCodeParseError {
    inner: ck::ResidueCodeParseError,
}
#[wasm_bindgen]
impl ResidueCodeParseError {
    pub fn input(&self) -> String {
        self.inner.input().into()
    }
}
