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
});
enum_value!(ResidueInfoKind, {
    UNKNOWN => Unknown, AA => Aa, AAD => Aad, PAA => Paa, MAA => Maa,
    RNA => Rna, DNA => Dna, BUF => Buf, HOH => Hoh, PYR => Pyr,
    KET => Ket, ELS => Els,
});
#[pymethods]
impl ResidueInfoKind {
    #[getter]
    fn name(&self) -> &'static str {
        ck::ResidueInfoKind::from(*self).name()
    }
}

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

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="fixed_projection._native", extends=PyValueError)]
pub struct ResidueCodeParseError {
    inner: ck::ResidueCodeParseError,
}
#[pymethods]
impl ResidueCodeParseError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    fn input(&self) -> &str {
        self.inner.input()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl ResidueCodeParseError {
    fn from_rust(inner: ck::ResidueCodeParseError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

#[derive(FromPyObject, IntoPyObject)]
pub enum ResidueCodeInput {
    Name(String),
    Code(u16),
}
#[pymethods]
impl ResidueCode {
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
    fn as_u16(&self) -> u16 {
        ck::ResidueCode::from(*self).as_u16()
    }
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct ResidueInfo {
    pub(crate) inner: ck::ResidueInfo,
}
#[pymethods]
impl ResidueInfo {
    #[getter]
    fn code(&self) -> ResidueCode {
        self.inner.code.into()
    }
    #[getter]
    fn name(&self) -> &'static str {
        self.inner.name
    }
    #[getter]
    fn kind(&self) -> ResidueInfoKind {
        self.inner.kind.into()
    }
    #[getter]
    fn linking_type(&self) -> u8 {
        self.inner.linking_type
    }
    #[getter]
    fn one_letter_code(&self) -> char {
        self.inner.one_letter_code
    }
    #[getter]
    fn hydrogen_count(&self) -> u8 {
        self.inner.hydrogen_count
    }
    #[getter]
    fn weight(&self) -> f32 {
        self.inner.weight
    }
    fn canonical_one_letter_code(&self) -> Option<char> {
        self.inner.canonical_one_letter_code()
    }
    fn parent_standard_code(&self) -> Option<ResidueCode> {
        self.inner.parent_standard_code().map(Into::into)
    }
    fn is_modified_amino_acid(&self) -> bool {
        self.inner.is_modified_amino_acid()
    }
    fn found(&self) -> bool {
        self.inner.found()
    }
    fn is_water(&self) -> bool {
        self.inner.is_water()
    }
    fn is_dna(&self) -> bool {
        self.inner.is_dna()
    }
    fn is_rna(&self) -> bool {
        self.inner.is_rna()
    }
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct ResidueIdentity {
    pub(crate) inner: ck::ResidueIdentity,
}
#[pymethods]
impl ResidueIdentity {
    #[new]
    fn construct(name: String) -> Self {
        Self::new(name)
    }
    #[staticmethod]
    fn new(name: String) -> Self {
        Self {
            inner: ck::ResidueIdentity::new(name),
        }
    }
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn code(&self) -> ResidueCode {
        self.inner.code().into()
    }
    fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: self.inner.info(),
        }
    }
    fn is_tabulated(&self) -> bool {
        self.inner.is_tabulated()
    }
}

#[pyfunction]
pub fn residue_info(index: usize) -> PyResult<ResidueInfo> {
    // Rust's raw table index has a valid-index precondition. Its checked public
    // counterpart validates that boundary without a copied length or table.
    ck::residue_info_checked(index)
        .map(|inner| ResidueInfo { inner })
        .ok_or_else(|| PyIndexError::new_err("residue information index out of range"))
}
#[pyfunction]
pub fn residue_info_checked(index: usize) -> Option<ResidueInfo> {
    ck::residue_info_checked(index).map(|inner| ResidueInfo { inner })
}
#[pyfunction]
pub fn find_residue_info_index(name: &str) -> usize {
    ck::find_residue_info_index(name)
}
#[pyfunction]
pub fn find_residue_info(name: &str) -> ResidueInfo {
    ResidueInfo {
        inner: ck::find_residue_info(name),
    }
}
#[pyfunction]
pub fn residue_code(name: &str) -> ResidueCode {
    ck::residue_code(name).into()
}
#[pyfunction]
pub fn expand_one_letter(code: char, kind: ResidueInfoKindInput) -> PyResult<Option<&'static str>> {
    Ok(ck::expand_one_letter(code, kind.into_rust()?))
}

#[pyfunction]
pub fn expand_one_letter_sequence(
    sequence: &str,
    kind: ResidueInfoKindInput,
) -> PyResult<Vec<String>> {
    ck::expand_one_letter_sequence(sequence, kind.into_rust()?)
        .map_err(ResidueSequenceError::from_rust)
}
