#![allow(non_camel_case_types)]
use crate::chemical_values::{EnumValue, enum_value};
use crate::error_framework::error_projection;
use cosmolkit as ck;
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use std::error::Error;

error_projection!(ResidueSequenceError, ResidueSequenceErrorKind, ck::ResidueSequenceError, "bio", {
    UnmatchedParenthesis,
    UnexpectedLetter { kind: &'static str, letter: char, source_code: i16 },
}, cause = |_inner: &ck::ResidueSequenceError| None::<PyErr>);

// Exhaustive native language projection of the canonical Rust discriminants,
// not a residue recognition/table implementation.
enum_value!(ResidueCode, {
    ALA, ARG, ASN, ABA, ASP, ASX, CYS, CSH, GLN,
    GLU, GLX, GLY, HIS, ILE, LEU, LYS, MET, MSE,
    ORN, PHE, PRO, SER, THR, TRP, TYR, UNK, VAL,
    SEC, PYL, SEP, TPO, PCA, CSO, PTR, KCX, CSD,
    LLP, CME, MLY, DAL, TYS, OCS, M3L, FME, ALY,
    HYP, CAS, CRO, CSX, DPR, DGL, DVA, CSS, DPN,
    DSN, DLE, HIC, NLE, MVA, MLZ, CR2, SAR, DAR,
    DLY, YCM, NRQ, CGU, R0TD, MLE, DAS, DTR, CXM,
    TPQ, DCY, DSG, DTY, DHI, MEN, DTH, SAC, DGN,
    AIB, SMC, IAS, CIR, BMT, DIL, FGA, PHI, CRQ,
    SME, GHP, MHO, NEP, TRQ, TOX, ALC, R3FG, SCH,
    MDO, MAA, GYS, MK8, CR8, KPI, SCY, DHA, OMY,
    CAF, R0AF, SNN, MHS, MLU, SNC, PHD, B3E, MEA,
    MED, OAS, GL3, FVA, PHL, CRF, OMZ, BFD, MEQ,
    DAB, AGM, PSU, R5MU, R7MG, OMG, UR3, OMC, R2MG,
    H2U, R4SU, OMU, R4OC, MA6, M2G, R1MA, R6MZ, CCC,
    R2MA, R1MG, R5BU, MIA, DOC, R8OG, R5CM, R3DR, BRU,
    CBR, HOH, DOD, HEM, SO4, GOL, EDO, NAG, PO4,
    ACT, PEG, MAN, FAD, BMA, ADP, DMS, ACE, NH2,
    MPD, MES, NAD, NAP, TRS, ATP, PG4, GDP, FUC,
    FMT, GAL, PGE, FMN, PLP, EPE, SF4, BME, CIT,
    BE7, MRD, MHA, BU3, PGO, BU2, PDO, BU1, PG6,
    R1BO, PE7, PG5, TFP, DHD, PEU, TAU, SBT, SAL,
    IOH, IPA, PIG, B3P, BTB, NHE, C8E, OTE, PE4,
    XPE, PE8, P33, N8E, R2OS, R1PS, CPS, DMX, MPO,
    GCD, DXG, CM5, ACA, ACN, CCN, GLC, DR6, NH4,
    AZI, BNG, BOG, BGC, BCN, BRO, CAC, CBX, ACY,
    CBM, CLO, R3CO, NCO, CU1, CYN, MA4, TAR, GLO,
    MTL, SOR, DMU, DDQ, DMF, DIO, DOX, R12P, SDS,
    LMT, EOH, EEE, EGL, FLO, TRT, FCY, FRU, GBL,
    GPX, HTO, HTG, B7G, C10, R16D, HEZ, IOD, IDO,
    ICI, ICT, TLA, LAT, LBT, LDA, MN3, MRY, MOH,
    BEQ, C15, MG8, POL, NO3, JEF, P4C, CE1, DIA,
    CXE, IPH, PIN, R15P, CRY, PGR, PGQ, SPD, SPK,
    SPM, SUC, TBU, TMA, TEP, SCN, TRE, ETF, R144,
    UMQ, URE, YT3, ZN2, FE2, R3NI, SIA, XYP, A2G,
    GLA, NDG, NGA, A, C, G, I, U, N,
    F, K, DA, DC, DG, DI, DT, DU, DN,
    AG, AL, BA, BR, CA, CD, CL, CM, CN,
    CO, CS, CU, FE, HG, LI, MG, MN, NA,
    NI, NO, PB, RB, SR, Y1, ZN, UNKNOWN,
}, constructor = custom, methods = {
    #[new]
    fn new(value: ResidueCodeInput) -> PyResult<Self> {
        match value {
            ResidueCodeInput::Name(name) => name
                .parse::<ck::ResidueCode>()
                .map(Into::into)
                .map_err(ResidueCodeParseError::from_rust),
            ResidueCodeInput::Code(code) => ck::residue_info_checked(usize::from(code))
                .map(|info| info.code.into())
                .ok_or_else(|| PyValueError::new_err(format!("invalid ResidueCode: {code}"))),
        }
    }
    /// Native numeric residue-table code.
    fn as_u16(&self) -> u16 {
        ck::ResidueCode::from(*self).as_u16()
    }
});
enum_value!(ResidueInfoKind, {
    UNKNOWN => Unknown, AA => Aa, AAD => Aad, PAA => Paa, MAA => Maa,
    RNA => Rna, DNA => Dna, BUF => Buf, HOH => Hoh, PYR => Pyr,
    KET => Ket, ELS => Els,
});

#[derive(FromPyObject, IntoPyObject)]
pub enum ResidueInfoKindInput {
    Enum(ResidueInfoKind),
    Name(String),
    Code(u8),
}
impl ResidueInfoKindInput {
    pub(crate) fn into_rust(self) -> PyResult<ck::ResidueInfoKind> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(name) => ResidueInfoKind::from_name(&name).map(Into::into),
            Self::Code(code) => match code {
                0 => Ok(ck::ResidueInfoKind::Unknown),
                1 => Ok(ck::ResidueInfoKind::Aa),
                2 => Ok(ck::ResidueInfoKind::Aad),
                3 => Ok(ck::ResidueInfoKind::Paa),
                4 => Ok(ck::ResidueInfoKind::Maa),
                5 => Ok(ck::ResidueInfoKind::Rna),
                6 => Ok(ck::ResidueInfoKind::Dna),
                7 => Ok(ck::ResidueInfoKind::Buf),
                8 => Ok(ck::ResidueInfoKind::Hoh),
                9 => Ok(ck::ResidueInfoKind::Pyr),
                10 => Ok(ck::ResidueInfoKind::Ket),
                11 => Ok(ck::ResidueInfoKind::Els),
                _ => Err(PyValueError::new_err(format!(
                    "invalid ResidueInfoKind code: {code}"
                ))),
            },
        }
    }
}

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct ResidueCodeParseError {
    inner: ck::ResidueCodeParseError,
}
}
crate::error_framework::error_projection! {methods impl ResidueCodeParseError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    /// Original residue-code text rejected by the parser.
    fn input(&self) -> &str {
        self.inner.input()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl ResidueCodeParseError {
    fn from_rust(inner: ck::ResidueCodeParseError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}

#[derive(FromPyObject, IntoPyObject)]
pub enum ResidueCodeInput {
    Name(String),
    Code(u16),
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Readonly residue-table entry containing chemical classification, sequence codes and tabulated composition information.
    pub struct ResidueInfo {
        pub(crate) inner: ck::ResidueInfo,
    }
}
#[pymethods]
impl ResidueInfo {
    #[getter]
    /// Canonical residue code associated with this table entry or identity.
    fn code(&self) -> ResidueCode {
        self.inner.code.into()
    }
    #[getter]
    /// Stored residue name.
    fn name(&self) -> &'static str {
        self.inner.name
    }
    #[getter]
    /// Chemical category of the residue-table entry.
    fn kind(&self) -> ResidueInfoKind {
        self.inner.kind.into()
    }
    #[getter]
    /// Native linking-type code for this residue.
    fn linking_type(&self) -> u8 {
        self.inner.linking_type
    }
    #[getter]
    /// One-letter code stored in the native residue table.
    fn one_letter_code(&self) -> char {
        self.inner.one_letter_code
    }
    #[getter]
    /// Native FASTA code, using X for a nonstandard residue.
    fn fasta_code(&self) -> char {
        self.inner.fasta_code()
    }
    #[getter]
    /// Hydrogen count recorded in the residue table.
    fn hydrogen_count(&self) -> u8 {
        self.inner.hydrogen_count
    }
    #[getter]
    /// Residue molecular weight recorded in the native table.
    fn weight(&self) -> f32 {
        self.inner.weight
    }
    #[getter]
    /// Canonical parent one-letter code, or None when no canonical code is known.
    fn canonical_one_letter_code(&self) -> Option<char> {
        self.inner.canonical_one_letter_code()
    }
    #[getter]
    /// Standard parent residue code for this entry, or None when unavailable.
    fn parent_standard_code(&self) -> Option<ResidueCode> {
        self.inner.parent_standard_code().map(Into::into)
    }
    #[getter]
    /// Whether the native table classifies this as a modified amino acid.
    fn is_modified_amino_acid(&self) -> bool {
        self.inner.is_modified_amino_acid()
    }
    #[getter]
    /// Whether the native residue table classifies this entry as an amino acid.
    fn is_amino_acid(&self) -> bool {
        self.inner.is_amino_acid()
    }
    #[getter]
    /// Whether this is a recognized entry rather than the unknown-residue sentinel.
    fn found(&self) -> bool {
        self.inner.found()
    }
    #[getter]
    /// Whether this entry is classified as water.
    fn is_water(&self) -> bool {
        self.inner.is_water()
    }
    #[getter]
    /// Whether this entry is classified as a DNA residue.
    fn is_dna(&self) -> bool {
        self.inner.is_dna()
    }
    #[getter]
    /// Whether this entry is classified as an RNA residue.
    fn is_rna(&self) -> bool {
        self.inner.is_rna()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Residue name preserved together with its canonical table lookup, including unrecognized names.
    pub struct ResidueIdentity {
        pub(crate) inner: ck::ResidueIdentity,
    }
}
#[pymethods]
impl ResidueIdentity {
    #[new]
    fn construct(name: String) -> Self {
        Self::new(name)
    }
    #[staticmethod]
    /// Create a residue identity preserving the supplied name and its native table lookup.
    fn new(name: String) -> Self {
        Self {
            inner: ck::ResidueIdentity::new(name),
        }
    }
    #[getter]
    /// Stored residue name.
    fn name(&self) -> &str {
        self.inner.name()
    }
    #[getter]
    /// Canonical residue code associated with this table entry or identity.
    fn code(&self) -> ResidueCode {
        self.inner.code().into()
    }
    #[getter]
    /// Native table information for this residue identity.
    fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: self.inner.info(),
        }
    }
    #[getter]
    /// Whether the identity matches a known native residue-table entry.
    fn is_tabulated(&self) -> bool {
        self.inner.is_tabulated()
    }
}

#[pyfunction]
/// Return residue-table information at the given index; invalid indices raise IndexError.
pub fn residue_info(index: usize) -> PyResult<ResidueInfo> {
    // Rust's raw table index has a valid-index precondition. Its checked public
    // counterpart validates that boundary without a copied length or table.
    ck::residue_info_checked(index)
        .map(|inner| ResidueInfo { inner })
        .ok_or_else(|| PyIndexError::new_err("residue information index out of range"))
}
#[pyfunction]
/// Return residue-table information at the given index, or None when out of bounds.
pub fn residue_info_checked(index: usize) -> Option<ResidueInfo> {
    ck::residue_info_checked(index).map(|inner| ResidueInfo { inner })
}
#[pyfunction]
/// Look up a residue name in the canonical table, preserving the native unknown-name sentinel.
pub fn find_residue_info_index(name: &str) -> usize {
    ck::find_residue_info_index(name)
}
#[pyfunction]
/// Look up canonical residue information by name, returning the native unknown-residue entry when unrecognized.
pub fn find_residue_info(name: &str) -> ResidueInfo {
    ResidueInfo {
        inner: ck::find_residue_info(name),
    }
}
#[pyfunction]
/// Convert a residue name to the canonical ResidueCode, preserving the native unknown-name value.
pub fn residue_code(name: &str) -> ResidueCode {
    ck::residue_code(name).into()
}
#[pyfunction]
/// Expand a one-letter residue code in the selected residue category; return None when the native table has no expansion.
pub fn expand_one_letter(code: char, kind: ResidueInfoKindInput) -> PyResult<Option<&'static str>> {
    Ok(ck::expand_one_letter(code, kind.into_rust()?))
}

#[pyfunction]
/// Expand one-letter sequence codes using the selected residue category and native sequence expansion rules.
pub fn expand_one_letter_sequence(
    sequence: &str,
    kind: ResidueInfoKindInput,
) -> PyResult<Vec<String>> {
    ck::expand_one_letter_sequence(sequence, kind.into_rust()?)
        .map_err(ResidueSequenceError::from_rust)
}
