//! RDKit tautomer-enumeration source reproduction.
//!
//! This module is the sole owner of tautomer-specific transforms, catalogs,
//! scoring, and enumeration state. SMARTS parsing and matching remain owned by
//! `crate::search`; this module stores and reuses their compiled query value.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{self, BufRead, BufReader, Write},
    path::{Path, PathBuf},
    string::FromUtf8Error,
    sync::{Arc, OnceLock},
};

use super::tautomer_transforms::{
    CURRENT_TAUTOMER_TRANSFORM_DEFINITIONS, V1_TAUTOMER_TRANSFORM_DEFINITIONS,
};
use crate::{
    AtomId, BondDirection, BondId, BondOrder, BondStereo, ChiralTag, Hybridization, Molecule,
    QueryGraph, RingFindingError, SmartsParseError, SmartsParseParams, StereoError, parse_smarts,
    read_parts::MoleculeReadParts,
};

pub(super) type TautomerTransformDefinition<'a> = (&'a str, &'a str, &'a str, &'a str);

#[derive(Debug, thiserror::Error)]
pub enum TautomerCatalogError {
    #[error("Bad input file {path}")]
    BadInputFile {
        path: PathBuf,
        #[source]
        source: io::Error,
    },
    #[error("Bad stream contents.")]
    BadStreamContents(#[source] io::Error),
    #[error("tautomer catalog input is not UTF-8")]
    InvalidUtf8(#[from] FromUtf8Error),
    #[error(transparent)]
    Transform(#[from] TautomerTransformError),
    #[error("tautomer transform index {index} is outside catalog size {len}")]
    TransformIndexOutOfRange { index: usize, len: usize },
    #[error("tautomer catalog deserialization is under construction in the source library")]
    DeserializationUnderConstruction,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TautomerScoreError {
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TautomerEnumerationError {
    #[error("tautomer index {index} is outside result size {len}")]
    IndexOutOfRange { index: usize, len: usize },
    #[error("candidate {canonical_smiles} has no materialized tautomer molecule")]
    MissingCandidateMolecule { canonical_smiles: String },
}

/// Structured failures produced while executing a tautomer-enumeration run.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TautomerRunError {
    #[error(transparent)]
    Smiles(#[from] crate::SmilesWriteError),
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Kekulize(#[from] crate::KekulizeError),
    #[error(transparent)]
    Substructure(#[from] crate::SubstructMatchError),
    #[error(transparent)]
    Stereo(#[from] StereoError),
    #[error(transparent)]
    MoleculeInvariant(#[from] crate::InvariantError),
    #[error(transparent)]
    Enumeration(#[from] TautomerEnumerationError),
    #[error("no tautomer candidate can be selected canonically")]
    NoCanonicalTautomer,
    #[error("tautomer atom table has {actual} rows, but the source has {expected}")]
    AtomCountMismatch { expected: usize, actual: usize },
    #[error("tautomer bond table has {actual} rows, but the source has {expected}")]
    BondCountMismatch { expected: usize, actual: usize },
    #[error("tautomer match has {actual} atom rows, expected {expected}")]
    AtomMappingCount { expected: usize, actual: usize },
    #[error("tautomer match has {actual} bond rows, expected {expected}")]
    BondMappingCount { expected: usize, actual: usize },
    #[error("tautomer match atom index {atom} is outside candidate atom count {atom_count}")]
    AtomOutOfRange { atom: usize, atom_count: usize },
    #[error("tautomer match bond index {bond} is outside candidate bond count {bond_count}")]
    BondOutOfRange { bond: usize, bond_count: usize },
    #[error("tautomer hydrogen count {count} on atom {atom} exceeds the modeled atom field")]
    HydrogenCountOutOfRange { atom: AtomId, count: u32 },
    #[error("tautomer formal charge {charge} on atom {atom} exceeds the modeled atom field")]
    FormalChargeOutOfRange { atom: AtomId, charge: i32 },
    #[error("tautomer partial sanitization failed: {0}")]
    Sanitize(Box<crate::OperationError>),
    #[error("tautomer candidate {canonical_smiles} has no kekulized branch")]
    MissingKekulizedBranch { canonical_smiles: String },
    #[error("tautomer candidate {canonical_smiles} has no materialized branch")]
    MissingTautomerBranch { canonical_smiles: String },
    #[error("tautomer expansion attempted to replace existing canonical key {canonical_smiles}")]
    DuplicateProductKey { canonical_smiles: String },
}

#[derive(Debug, thiserror::Error)]
pub enum TautomerCanonicalizationError {
    #[error("tautomer enumeration failed during canonicalization: {source}")]
    Enumeration {
        #[source]
        source: Box<crate::OperationError>,
    },
    #[error(transparent)]
    Selection(#[from] TautomerRunError),
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
enum TautomerStereoTransitionError {
    #[error("tautomer atom table has {actual} rows, but the source has {expected}")]
    AtomCountMismatch { expected: usize, actual: usize },
    #[error("tautomer bond table has {actual} rows, but the source has {expected}")]
    BondCountMismatch { expected: usize, actual: usize },
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Stereo(#[from] StereoError),
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub(crate) enum TautomerInitializationError {
    #[error(transparent)]
    Smiles(#[from] crate::SmilesWriteError),
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Kekulize(#[from] crate::KekulizeError),
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub(crate) enum TautomerTransformApplicationError {
    #[error("tautomer match has {actual} atom rows, expected {expected}")]
    AtomMappingCount { expected: usize, actual: usize },
    #[error("tautomer match has {actual} bond rows, expected {expected}")]
    BondMappingCount { expected: usize, actual: usize },
    #[error("tautomer match atom index {atom} is outside candidate atom count {atom_count}")]
    AtomOutOfRange { atom: usize, atom_count: usize },
    #[error("tautomer match bond index {bond} is outside candidate bond count {bond_count}")]
    BondOutOfRange { bond: usize, bond_count: usize },
    #[error("tautomer hydrogen count {count} on atom {atom} exceeds the modeled atom field")]
    HydrogenCountOutOfRange { atom: AtomId, count: u32 },
    #[error("tautomer formal charge {charge} on atom {atom} exceeds the modeled atom field")]
    FormalChargeOutOfRange { atom: AtomId, charge: i32 },
    #[error(transparent)]
    MoleculeInvariant(#[from] crate::InvariantError),
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error("tautomer partial sanitization failed: {0}")]
    Sanitize(crate::OperationError),
    #[error(transparent)]
    Stereo(#[from] TautomerStereoTransitionError),
    #[error(transparent)]
    Smiles(#[from] crate::SmilesWriteError),
    #[error(transparent)]
    Kekulize(#[from] crate::KekulizeError),
}

impl From<TautomerInitializationError> for TautomerRunError {
    fn from(error: TautomerInitializationError) -> Self {
        match error {
            TautomerInitializationError::Smiles(source) => Self::Smiles(source),
            TautomerInitializationError::Valence(source) => Self::Valence(source),
            TautomerInitializationError::RingFinding(source) => Self::RingFinding(source),
            TautomerInitializationError::Kekulize(source) => Self::Kekulize(source),
        }
    }
}

impl From<TautomerScoreError> for TautomerRunError {
    fn from(error: TautomerScoreError) -> Self {
        match error {
            TautomerScoreError::RingFinding(source) => Self::RingFinding(source),
            TautomerScoreError::Valence(source) => Self::Valence(source),
        }
    }
}

impl From<TautomerStereoTransitionError> for TautomerRunError {
    fn from(error: TautomerStereoTransitionError) -> Self {
        match error {
            TautomerStereoTransitionError::AtomCountMismatch { expected, actual } => {
                Self::AtomCountMismatch { expected, actual }
            }
            TautomerStereoTransitionError::BondCountMismatch { expected, actual } => {
                Self::BondCountMismatch { expected, actual }
            }
            TautomerStereoTransitionError::RingFinding(source) => Self::RingFinding(source),
            TautomerStereoTransitionError::Stereo(source) => Self::Stereo(source),
            TautomerStereoTransitionError::Valence(source) => Self::Valence(source),
        }
    }
}

impl From<TautomerTransformApplicationError> for TautomerRunError {
    fn from(error: TautomerTransformApplicationError) -> Self {
        match error {
            TautomerTransformApplicationError::AtomMappingCount { expected, actual } => {
                Self::AtomMappingCount { expected, actual }
            }
            TautomerTransformApplicationError::BondMappingCount { expected, actual } => {
                Self::BondMappingCount { expected, actual }
            }
            TautomerTransformApplicationError::AtomOutOfRange { atom, atom_count } => {
                Self::AtomOutOfRange { atom, atom_count }
            }
            TautomerTransformApplicationError::BondOutOfRange { bond, bond_count } => {
                Self::BondOutOfRange { bond, bond_count }
            }
            TautomerTransformApplicationError::HydrogenCountOutOfRange { atom, count } => {
                Self::HydrogenCountOutOfRange { atom, count }
            }
            TautomerTransformApplicationError::FormalChargeOutOfRange { atom, charge } => {
                Self::FormalChargeOutOfRange { atom, charge }
            }
            TautomerTransformApplicationError::MoleculeInvariant(source) => {
                Self::MoleculeInvariant(source)
            }
            TautomerTransformApplicationError::Valence(source) => Self::Valence(source),
            TautomerTransformApplicationError::RingFinding(source) => Self::RingFinding(source),
            TautomerTransformApplicationError::Sanitize(source) => Self::Sanitize(Box::new(source)),
            TautomerTransformApplicationError::Stereo(source) => source.into(),
            TautomerTransformApplicationError::Smiles(source) => Self::Smiles(source),
            TautomerTransformApplicationError::Kekulize(source) => Self::Kekulize(source),
        }
    }
}

impl From<TautomerExpansionError<TautomerRunError>> for TautomerRunError {
    fn from(error: TautomerExpansionError<TautomerRunError>) -> Self {
        match error {
            TautomerExpansionError::MissingKekulizedBranch { canonical_smiles } => {
                Self::MissingKekulizedBranch { canonical_smiles }
            }
            TautomerExpansionError::DuplicateProductKey { canonical_smiles } => {
                Self::DuplicateProductKey { canonical_smiles }
            }
            TautomerExpansionError::Backend(source) => source,
        }
    }
}

impl From<TautomerPruningError<TautomerRunError>> for TautomerRunError {
    fn from(error: TautomerPruningError<TautomerRunError>) -> Self {
        match error {
            TautomerPruningError::MissingTautomerBranch { canonical_smiles } => {
                Self::MissingTautomerBranch { canonical_smiles }
            }
            TautomerPruningError::Backend(source) => source,
        }
    }
}

/// A structurally invalid tautomer transform.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TautomerTransformError {
    #[error("a tautomer transform requires donor and acceptor atoms; query has {actual} atoms")]
    MissingDonorOrAcceptor { actual: usize },
    #[error("tautomer transform has {actual} bond edits, but its query has {expected} bonds")]
    BondEditCount { expected: usize, actual: usize },
    #[error("tautomer transform has {actual} charge edits, but its query has {expected} atoms")]
    ChargeEditCount { expected: usize, actual: usize },
    #[error("Charge symbol not recognised.")]
    ChargeSymbolNotRecognised,
    #[error("cannot parse tautomer SMARTS: {smarts}")]
    CannotParseSmarts {
        smarts: String,
        #[source]
        source: SmartsParseError,
    },
}

fn string_to_bond_types(bond_str: &str) -> Vec<BondOrder> {
    // RDKit✔️🔝: std::vector<Bond::BondType> stringToBondType(std::string bond_str) {
    // RDKit✔️🔝:   std::vector<Bond::BondType> bonds;
    // RDKit✔️🔝:   for (const auto &c : bond_str) {
    // RDKit✔️🔝:     switch (c) {
    // RDKit✔️🔝:       case '-':
    // RDKit✔️🔝:         bonds.push_back(Bond::SINGLE);
    // RDKit✔️🔝:         break;
    // RDKit✔️🔝:       case '=':
    // RDKit✔️🔝:         bonds.push_back(Bond::DOUBLE);
    // RDKit✔️🔝:         break;
    // RDKit✔️🔝:       case '#':
    // RDKit✔️🔝:         bonds.push_back(Bond::TRIPLE);
    // RDKit✔️🔝:         break;
    // RDKit✔️🔝:       case ':':
    // RDKit✔️🔝:         bonds.push_back(Bond::AROMATIC);
    // RDKit✔️🔝:         break;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return bonds;
    // RDKit✔️🔝: }
    // Borrowing avoids the source function's by-value string copy without
    // changing its single-pass scan, source ordering, or ignored-byte rules.
    let mut bonds = Vec::new();
    for byte in bond_str.bytes() {
        match byte {
            b'-' => bonds.push(BondOrder::Single),
            b'=' => bonds.push(BondOrder::Double),
            b'#' => bonds.push(BondOrder::Triple),
            b':' => bonds.push(BondOrder::Aromatic),
            _ => {}
        }
    }
    bonds
}

fn string_to_charges(charge_str: &str) -> Result<Vec<i32>, TautomerTransformError> {
    // RDKit✔️🔝: std::vector<int> stringToCharge(std::string charge_str) {
    // RDKit✔️🔝:   std::vector<int> charges;
    // RDKit✔️🔝:   for (const auto &c : charge_str) {
    // RDKit✔️🔝:     switch (c) {
    // RDKit✔️🔝:       case '+':
    // RDKit✔️🔝:         charges.push_back(1);
    // RDKit✔️🔝:         break;
    // RDKit✔️🔝:       case '0':
    // RDKit✔️🔝:         charges.push_back(0);
    // RDKit✔️🔝:         break;
    // RDKit✔️🔝:       case '-':
    // RDKit✔️🔝:         charges.push_back(-1);
    // RDKit✔️🔝:         break;
    // RDKit✔️🔝:       default:
    // RDKit✔️🔝:         throw ValueErrorException("Charge symbol not recognised.");
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return charges;
    // RDKit✔️🔝: }
    // Borrowing avoids the source function's by-value string copy. Byte-wise
    // iteration deliberately retains C++ `std::string` behavior for non-ASCII
    // input: the first unrecognized byte returns the same error category.
    let mut charges = Vec::new();
    for byte in charge_str.bytes() {
        match byte {
            b'+' => charges.push(1),
            b'0' => charges.push(0),
            b'-' => charges.push(-1),
            _ => return Err(TautomerTransformError::ChargeSymbolNotRecognised),
        }
    }
    Ok(charges)
}

fn transform_from_fields(
    name: &str,
    smarts: &str,
    bond_str: &str,
    charge_str: &str,
) -> Result<TautomerTransform, TautomerTransformError> {
    // RDKit✔️🔝: std::unique_ptr<MolStandardize::TautomerTransform> getTautomer(
    // RDKit✔️🔝:     const std::string &name, const std::string &smarts,
    // RDKit✔️🔝:     const std::string &bond_str, const std::string &charge_str) {
    // RDKit✔️🔝:   std::vector<Bond::BondType> bond_types =
    // RDKit✔️🔝:       MolStandardize::stringToBondType(bond_str);
    // RDKit✔️🔝:   std::vector<int> charges = MolStandardize::stringToCharge(charge_str);
    // RDKit✔️🔝:
    // RDKit✔️🔝:   ROMol *tautomer = SmartsToMol(smarts);
    // RDKit✔️🔝:   if (!tautomer) {
    // RDKit✔️🔝:     throw ValueErrorException("cannot parse tautomer SMARTS: " + smarts);
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   tautomer->setProp(common_properties::_Name, name);
    // RDKit✔️🔝:   return std::make_unique<MolStandardize::TautomerTransform>(
    // RDKit✔️🔝:       tautomer, bond_types, charges);
    // RDKit✔️🔝: }
    // Borrowed fields avoid four source-side string copies. The parser,
    // token conversion order, compiled query ownership, and error order remain
    // unchanged.
    let bond_types = string_to_bond_types(bond_str);
    let charges = string_to_charges(charge_str)?;
    let query = parse_smarts(smarts, &SmartsParseParams::default()).map_err(|source| {
        TautomerTransformError::CannotParseSmarts {
            smarts: smarts.to_owned(),
            source,
        }
    })?;
    TautomerTransform::new(name, query, bond_types, charges)
}

fn transform_from_line(line: &str) -> Result<Option<TautomerTransform>, TautomerTransformError> {
    // RDKit✔️✔️: std::unique_ptr<MolStandardize::TautomerTransform> getTautomer(
    // RDKit✔️✔️:     const std::string &tmpStr) {
    // RDKit✔️✔️:   if (tmpStr.length() == 0 || tmpStr.substr(0, 2) == "//") {
    // RDKit✔️✔️:     // empty or comment line
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    if line.is_empty() || line.starts_with("//") {
        return Ok(None);
    }

    // RDKit✔️✔️:   boost::char_separator<char> tabSep("\t");
    // RDKit✔️✔️:   tokenizer tokens(tmpStr, tabSep);
    // RDKit✔️✔️:   std::vector<std::string> result(tokens.begin(), tokens.end());
    // `char_separator` drops delimiters and empty tokens by default.
    let fields: Vec<&str> = line.split('\t').filter(|field| !field.is_empty()).collect();

    // RDKit✔️✔️:   // tautomer information to collect from each line
    // RDKit✔️✔️:   std::string name;
    // RDKit✔️✔️:   std::string smarts;
    // RDKit✔️✔️:   std::string bond_str;
    // RDKit✔️✔️:   std::string charge_str;
    let (mut name, mut smarts, mut bond_str, mut charge_str) = ("", "", "", "");

    // RDKit✔️✔️:   // line must have at least two tab separated values
    // RDKit✔️✔️:   if (result.size() < 2) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << "Invalid line: " << tmpStr << std::endl;
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    if fields.len() < 2 {
        return Ok(None);
    }

    // RDKit✔️✔️:   // line only has name and smarts
    // RDKit✔️✔️:   if (result.size() == 2) {
    // RDKit✔️✔️:     name = result[0];
    // RDKit✔️✔️:     smarts = result[1];
    // RDKit✔️✔️:   }
    if fields.len() == 2 {
        name = fields[0];
        smarts = fields[1];
    }
    // RDKit✔️✔️:   // line has name, smarts, bonds
    // RDKit✔️✔️:   if (result.size() == 3) {
    // RDKit✔️✔️:     name = result[0];
    // RDKit✔️✔️:     smarts = result[1];
    // RDKit✔️✔️:     bond_str = result[2];
    // RDKit✔️✔️:   }
    if fields.len() == 3 {
        name = fields[0];
        smarts = fields[1];
        bond_str = fields[2];
    }
    // RDKit✔️✔️:   // line has name, smarts, bonds, charges
    // RDKit✔️✔️:   if (result.size() == 4) {
    // RDKit✔️✔️:     name = result[0];
    // RDKit✔️✔️:     smarts = result[1];
    // RDKit✔️✔️:     bond_str = result[2];
    // RDKit✔️✔️:     charge_str = result[3];
    // RDKit✔️✔️:   }
    if fields.len() == 4 {
        name = fields[0];
        smarts = fields[1];
        bond_str = fields[2];
        charge_str = fields[3];
    }

    // RDKit✔️✔️:   boost::erase_all(smarts, " ");
    // RDKit✔️✔️:   boost::erase_all(name, " ");
    // RDKit✔️✔️:   boost::erase_all(bond_str, " ");
    // RDKit✔️✔️:   boost::erase_all(charge_str, " ");
    let name = name.replace(' ', "");
    let smarts = smarts.replace(' ', "");
    let bond_str = bond_str.replace(' ', "");
    let charge_str = charge_str.replace(' ', "");

    // RDKit✔️✔️:   return getTautomer(name, smarts, bond_str, charge_str);
    // RDKit✔️✔️: }
    transform_from_fields(&name, &smarts, &bond_str, &charge_str).map(Some)
}

fn read_source_line<R: BufRead>(reader: &mut R) -> io::Result<(Option<Vec<u8>>, bool)> {
    const SOURCE_PAYLOAD_LIMIT: usize = 511;
    let mut line = Vec::with_capacity(SOURCE_PAYLOAD_LIMIT);

    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            return Ok(((!line.is_empty()).then_some(line), false));
        }

        let remaining = SOURCE_PAYLOAD_LIMIT - line.len();
        let visible = &available[..available.len().min(remaining)];
        if let Some(newline) = visible.iter().position(|byte| *byte == b'\n') {
            line.extend_from_slice(&visible[..newline]);
            reader.consume(newline + 1);
            return Ok((Some(line), false));
        }

        let consumed = visible.len();
        line.extend_from_slice(visible);
        reader.consume(consumed);
        if line.len() == SOURCE_PAYLOAD_LIMIT {
            return Ok((Some(line), true));
        }
    }
}

fn read_transforms<R: BufRead>(
    reader: &mut R,
    n_to_read: i32,
) -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    // RDKit✔️✔️: std::vector<TautomerTransform> readTautomers(std::istream &inStream,
    // RDKit✔️✔️:                                              int nToRead) {
    // RDKit✔️✔️:   if (inStream.bad()) {
    // RDKit✔️✔️:     throw BadFileException("Bad stream contents.");
    // RDKit✔️✔️:   }
    // Rust reports a `BufRead` failure at the first operation that observes it,
    // because the trait has no independent pre-read `badbit` query.
    let mut transforms = if n_to_read > 0 {
        // RDKit✔️✔️:   std::vector<TautomerTransform> tautomers;
        // RDKit✔️✔️:   if (nToRead > 0) {
        // RDKit✔️✔️:     tautomers.reserve(nToRead);
        // RDKit✔️✔️:   }
        Vec::with_capacity(n_to_read as usize)
    } else {
        Vec::new()
    };

    // RDKit✔️✔️:   const int MAX_LINE_LEN = 512;
    // RDKit✔️✔️:   char inLine[MAX_LINE_LEN];
    // RDKit✔️✔️:   std::string tmpstr;
    // RDKit✔️✔️:   int nRead = 0;
    // RDKit✔️✔️:   while (!inStream.eof() && !inStream.fail() &&
    // RDKit✔️✔️:          (nToRead < 0 || nRead < nToRead)) {
    while n_to_read < 0 || transforms.len() < n_to_read as usize {
        // RDKit✔️✔️:     inStream.getline(inLine, MAX_LINE_LEN, '\n');
        // RDKit✔️✔️:     tmpstr = inLine;
        let (line, source_failbit) =
            read_source_line(reader).map_err(TautomerCatalogError::BadStreamContents)?;
        let Some(line) = line else {
            break;
        };
        let line = String::from_utf8(line)?;

        // RDKit✔️✔️:     // parse the tautomer on this line (if there is one)
        // RDKit✔️✔️:     auto transform = getTautomer(tmpstr);
        // RDKit✔️✔️:     if (transform) {
        // RDKit✔️✔️:       tautomers.emplace_back(*transform);
        // RDKit✔️✔️:       nRead++;
        // RDKit✔️✔️:     }
        if let Some(transform) = transform_from_line(&line)? {
            transforms.push(transform);
        }
        // RDKit✔️✔️:   }
        if source_failbit {
            break;
        }
    }

    // RDKit✔️✔️:   return tautomers;
    // RDKit✔️✔️: }
    Ok(transforms)
}

fn read_transforms_from_file(
    file_name: impl AsRef<Path>,
) -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    // RDKit✔️✔️: std::vector<TautomerTransform> readTautomers(std::string fileName) {
    // RDKit✔️✔️:   std::ifstream inStream(fileName.c_str());
    // RDKit✔️✔️:   if ((!inStream) || (inStream.bad())) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Bad input file " << fileName;
    // RDKit✔️✔️:     throw BadFileException(errout.str());
    // RDKit✔️✔️:   }
    let path = file_name.as_ref();
    let file = File::open(path).map_err(|source| TautomerCatalogError::BadInputFile {
        path: path.to_owned(),
        source,
    })?;
    // RDKit✔️✔️:   std::vector<TautomerTransform> tautomers = readTautomers(inStream);
    // RDKit✔️✔️:   return tautomers;
    // RDKit✔️✔️: }
    read_transforms(&mut BufReader::new(file), -1)
}

fn read_transforms_from_definitions(
    data: &[TautomerTransformDefinition<'_>],
) -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    // RDKit✔️✔️: std::vector<TautomerTransform> readTautomers(
    // RDKit✔️✔️:     const TautomerTransformDefs &data) {
    // RDKit✔️✔️:   std::vector<TautomerTransform> tautomers;
    let mut transforms = Vec::new();
    // RDKit✔️✔️:   for (const auto &tpl : data) {
    // RDKit✔️✔️:     auto transform = getTautomer(std::get<0>(tpl), std::get<1>(tpl),
    // RDKit✔️✔️:                                  std::get<2>(tpl), std::get<3>(tpl));
    // RDKit✔️✔️:     if (transform) {
    // RDKit✔️✔️:       tautomers.emplace_back(*transform);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for &(name, smarts, bonds, charges) in data {
        transforms.push(transform_from_fields(name, smarts, bonds, charges)?);
    }
    // RDKit✔️✔️:   return tautomers;
    // RDKit✔️✔️: }
    Ok(transforms)
}

fn current_builtin_transforms() -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    read_transforms_from_definitions(CURRENT_TAUTOMER_TRANSFORM_DEFINITIONS)
}

fn v1_builtin_transforms() -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    read_transforms_from_definitions(V1_TAUTOMER_TRANSFORM_DEFINITIONS)
}

/// One compiled, source-ordered tautomer transformation.
///
/// The first and last query atoms are respectively the hydrogen donor and
/// acceptor. An empty bond-edit vector selects RDKit's alternating
/// single/double update; an empty charge-edit vector leaves formal charges
/// unchanged.
#[derive(Debug, PartialEq)]
pub struct TautomerTransform {
    query: QueryGraph,
    bond_types: Vec<BondOrder>,
    charges: Vec<i32>,
}

/// An ordered collection of compiled tautomer transforms.
#[derive(Debug, PartialEq)]
pub struct TautomerCatalog {
    transforms: Vec<TautomerTransform>,
}

/// One named SMARTS contribution to tautomer canonicalization.
#[derive(Debug, Clone)]
pub struct TautomerScoreTerm {
    name: String,
    smarts: String,
    score: i32,
    matcher: Option<QueryGraph>,
}

/// Source-defined components of a tautomer canonicalization score.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TautomerScore {
    ring: i32,
    substructure: i32,
    hetero_hydrogen: i32,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TautomerCandidate<M = Molecule> {
    pub(crate) tautomer: Option<M>,
    pub(crate) kekulized: Option<M>,
    pub(crate) num_modified_atoms: usize,
    pub(crate) num_modified_bonds: usize,
    pub(crate) done: bool,
}

pub(crate) type SmilesTautomerMap<M = Molecule> = BTreeMap<String, TautomerCandidate<M>>;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TautomerInitializationPlan {
    pub(crate) canonical_smiles: String,
    pub(crate) valence_update: Option<crate::ValenceAssignment>,
    pub(crate) rings_update: Option<crate::RingInfo>,
    pub(crate) target_derived_cache: crate::molecule::DerivedCacheBlock,
    pub(crate) kekulize_assignment: crate::kekulize::KekulizeAssignment,
    pub(crate) atom_count: usize,
    pub(crate) bond_count: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TautomerProductPlan {
    pub(crate) topology: crate::molecule::TopologyBlock,
    pub(crate) properties: crate::MoleculeProperties,
    pub(crate) derived_cache: crate::molecule::DerivedCacheBlock,
    pub(crate) canonical_smiles: String,
    pub(crate) kekulize_assignment: crate::kekulize::KekulizeAssignment,
    pub(crate) modified_atoms: BTreeSet<AtomId>,
    pub(crate) modified_bonds: BTreeSet<BondId>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TautomerStereoUpdatePlan {
    pub(crate) topology: crate::molecule::TopologyBlock,
    pub(crate) properties: crate::MoleculeProperties,
    pub(crate) derived_cache: crate::molecule::DerivedCacheBlock,
    pub(crate) changed: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum TautomerTransformAttempt {
    RecoverableKekulizeFailure {
        modified_atoms: BTreeSet<AtomId>,
        modified_bonds: BTreeSet<BondId>,
    },
    Duplicate {
        canonical_smiles: String,
        modified_atoms: BTreeSet<AtomId>,
        modified_bonds: BTreeSet<BondId>,
    },
    Product(TautomerProductPlan),
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TautomerExpandedProduct<M> {
    pub(crate) tautomer: M,
    pub(crate) kekulized: M,
    pub(crate) canonical_smiles: String,
    pub(crate) modified_atoms: BTreeSet<AtomId>,
    pub(crate) modified_bonds: BTreeSet<BondId>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum TautomerExpansionAttempt<M> {
    RecoverableKekulizeFailure {
        modified_atoms: BTreeSet<AtomId>,
        modified_bonds: BTreeSet<BondId>,
    },
    Duplicate {
        canonical_smiles: String,
        modified_atoms: BTreeSet<AtomId>,
        modified_bonds: BTreeSet<BondId>,
    },
    Product(TautomerExpandedProduct<M>),
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TautomerExpansionState<M> {
    pub(crate) candidates: SmilesTautomerMap<M>,
    pub(crate) modified_atoms: BTreeSet<AtomId>,
    pub(crate) modified_bonds: BTreeSet<BondId>,
    pub(crate) status: TautomerEnumerationStatus,
    pub(crate) num_transforms: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct TautomerExpansionPass {
    pub(crate) bailed_out: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct TautomerPruningPass {
    pub(crate) completed: bool,
    pub(crate) bailed_out: bool,
}

#[derive(Debug, thiserror::Error)]
pub(crate) enum TautomerExpansionError<E> {
    #[error("tautomer candidate {canonical_smiles} has no kekulized branch")]
    MissingKekulizedBranch { canonical_smiles: String },
    #[error("tautomer expansion attempted to replace existing canonical key {canonical_smiles}")]
    DuplicateProductKey { canonical_smiles: String },
    #[error(transparent)]
    Backend(E),
}

#[derive(Debug, thiserror::Error)]
pub(crate) enum TautomerPruningError<E> {
    #[error("tautomer candidate {canonical_smiles} has no materialized branch")]
    MissingTautomerBranch { canonical_smiles: String },
    #[error(transparent)]
    Backend(E),
}

impl TautomerInitializationPlan {
    pub(crate) fn into_candidate_map<M>(self, tautomer: M, kekulized: M) -> SmilesTautomerMap<M> {
        BTreeMap::from([(
            self.canonical_smiles,
            TautomerCandidate {
                tautomer: Some(tautomer),
                kekulized: Some(kekulized),
                num_modified_atoms: 0,
                num_modified_bonds: 0,
                done: false,
            },
        )])
    }
}

/// Completion state of a tautomer-enumeration run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TautomerEnumerationStatus {
    #[default]
    Completed,
    MaxTautomersReached,
    MaxTransformsReached,
    Canceled,
}

/// Source-compatible limits and stereochemistry policies for enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TautomerOptions {
    max_tautomers: u32,
    max_transforms: u32,
    remove_sp3_stereo: bool,
    remove_bond_stereo: bool,
    remove_isotopic_hydrogens: bool,
    reassign_stereo: bool,
}

impl Default for TautomerOptions {
    fn default() -> Self {
        Self {
            max_tautomers: 1000,
            max_transforms: 1000,
            remove_sp3_stereo: true,
            remove_bond_stereo: true,
            remove_isotopic_hydrogens: true,
            reassign_stereo: true,
        }
    }
}

impl TautomerOptions {
    #[must_use]
    pub const fn max_tautomers(self) -> u32 {
        self.max_tautomers
    }

    pub fn set_max_tautomers(&mut self, value: u32) {
        self.max_tautomers = value;
    }

    #[must_use]
    pub const fn with_max_tautomers(mut self, value: u32) -> Self {
        self.max_tautomers = value;
        self
    }

    #[must_use]
    pub const fn max_transforms(self) -> u32 {
        self.max_transforms
    }

    pub fn set_max_transforms(&mut self, value: u32) {
        self.max_transforms = value;
    }

    #[must_use]
    pub const fn with_max_transforms(mut self, value: u32) -> Self {
        self.max_transforms = value;
        self
    }

    #[must_use]
    pub const fn remove_sp3_stereo(self) -> bool {
        self.remove_sp3_stereo
    }

    pub fn set_remove_sp3_stereo(&mut self, value: bool) {
        self.remove_sp3_stereo = value;
    }

    #[must_use]
    pub const fn with_remove_sp3_stereo(mut self, value: bool) -> Self {
        self.remove_sp3_stereo = value;
        self
    }

    #[must_use]
    pub const fn remove_bond_stereo(self) -> bool {
        self.remove_bond_stereo
    }

    pub fn set_remove_bond_stereo(&mut self, value: bool) {
        self.remove_bond_stereo = value;
    }

    #[must_use]
    pub const fn with_remove_bond_stereo(mut self, value: bool) -> Self {
        self.remove_bond_stereo = value;
        self
    }

    #[must_use]
    pub const fn remove_isotopic_hydrogens(self) -> bool {
        self.remove_isotopic_hydrogens
    }

    pub fn set_remove_isotopic_hydrogens(&mut self, value: bool) {
        self.remove_isotopic_hydrogens = value;
    }

    #[must_use]
    pub const fn with_remove_isotopic_hydrogens(mut self, value: bool) -> Self {
        self.remove_isotopic_hydrogens = value;
        self
    }

    #[must_use]
    pub const fn reassign_stereo(self) -> bool {
        self.reassign_stereo
    }

    pub fn set_reassign_stereo(&mut self, value: bool) {
        self.reassign_stereo = value;
    }

    #[must_use]
    pub const fn with_reassign_stereo(mut self, value: bool) -> Self {
        self.reassign_stereo = value;
        self
    }
}

/// Borrowed cancellation hook evaluated at the source callback point.
pub trait TautomerEnumerationCallback: Send + Sync {
    fn should_continue(&self, molecule: &Molecule, result: &TautomerEnumeration) -> bool;
}

/// Configured tautomer enumerator backed by one shared immutable catalog.
pub struct TautomerEnumerator<'callback> {
    catalog: Arc<TautomerCatalog>,
    options: TautomerOptions,
    callback: Option<&'callback dyn TautomerEnumerationCallback>,
}

impl std::fmt::Debug for TautomerEnumerator<'_> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("TautomerEnumerator")
            .field("catalog", &self.catalog)
            .field("options", &self.options)
            .field("has_callback", &self.callback.is_some())
            .finish()
    }
}

impl<'callback> Clone for TautomerEnumerator<'callback> {
    fn clone(&self) -> Self {
        // RDKit✔️✔️: TautomerEnumerator(const TautomerEnumerator &other)
        // RDKit✔️✔️:     : dp_catalog(other.dp_catalog),
        // RDKit✔️✔️:       d_callback(other.d_callback),
        // RDKit✔️✔️:       d_maxTautomers(other.d_maxTautomers),
        // RDKit✔️✔️:       d_maxTransforms(other.d_maxTransforms),
        // RDKit✔️✔️:       d_removeSp3Stereo(other.d_removeSp3Stereo),
        // RDKit✔️✔️:       d_removeBondStereo(other.d_removeBondStereo),
        // RDKit✔️✔️:       d_removeIsotopicHs(other.d_removeIsotopicHs),
        // RDKit✔️✔️:       d_reassignStereo(other.d_reassignStereo) {}
        Self {
            catalog: Arc::clone(&self.catalog),
            options: self.options,
            callback: self.callback,
        }
    }

    fn clone_from(&mut self, source: &Self) {
        // RDKit✔️✔️: TautomerEnumerator &operator=(const TautomerEnumerator &other) {
        // RDKit✔️✔️:   if (this == &other) {
        // RDKit✔️✔️:     return *this;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   dp_catalog = other.dp_catalog;
        // RDKit✔️✔️:   d_callback = other.d_callback;
        // RDKit✔️✔️:   d_maxTautomers = other.d_maxTautomers;
        // RDKit✔️✔️:   d_maxTransforms = other.d_maxTransforms;
        // RDKit✔️✔️:   d_removeSp3Stereo = other.d_removeSp3Stereo;
        // RDKit✔️✔️:   d_removeBondStereo = other.d_removeBondStereo;
        // RDKit✔️✔️:   d_removeIsotopicHs = other.d_removeIsotopicHs;
        // RDKit✔️✔️:   d_reassignStereo = other.d_reassignStereo;
        // RDKit✔️✔️:   return *this;
        // RDKit✔️✔️: }
        Arc::clone_from(&mut self.catalog, &source.catalog);
        self.options = source.options;
        self.callback = source.callback;
    }
}

impl Default for TautomerEnumerator<'static> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'callback> TautomerEnumerator<'callback> {
    #[must_use]
    pub fn new() -> Self {
        Self::from_options(TautomerOptions::default())
    }

    #[must_use]
    pub fn from_options(options: TautomerOptions) -> Self {
        // RDKit✔️✔️: inline TautomerEnumerator *tautomerEnumeratorFromParams(
        // RDKit✔️✔️:     const CleanupParameters &params) {
        // RDKit✔️✔️:   return new TautomerEnumerator(params);
        // RDKit✔️✔️: }
        Self {
            catalog: Arc::new(TautomerCatalog::default()),
            options,
            callback: None,
        }
    }

    #[must_use]
    pub fn from_catalog(catalog: TautomerCatalog) -> Self {
        // RDKit✔️✔️: TautomerEnumerator(TautomerCatalog *tautCat)
        // RDKit✔️✔️:     : dp_catalog(tautCat),
        // RDKit✔️✔️:       d_maxTautomers(1000),
        // RDKit✔️✔️:       d_maxTransforms(1000),
        // RDKit✔️✔️:       d_removeSp3Stereo(true),
        // RDKit✔️✔️:       d_removeBondStereo(true),
        // RDKit✔️✔️:       d_removeIsotopicHs(true),
        // RDKit✔️✔️:       d_reassignStereo(true) {}
        Self {
            catalog: Arc::new(catalog),
            options: TautomerOptions::default(),
            callback: None,
        }
    }

    #[must_use]
    pub fn from_catalog_and_options(catalog: TautomerCatalog, options: TautomerOptions) -> Self {
        Self {
            catalog: Arc::new(catalog),
            options,
            callback: None,
        }
    }

    pub fn v1() -> Result<Self, TautomerCatalogError> {
        // RDKit✔️✔️: inline TautomerEnumerator *getV1TautomerEnumerator() {
        // RDKit✔️✔️:   TautomerCatalogParams tparms(
        // RDKit✔️✔️:       MolStandardize::defaults::defaultTautomerTransformsv1);
        // RDKit✔️✔️:   return new TautomerEnumerator(new TautomerCatalog(&tparms));
        // RDKit✔️✔️: }
        Ok(Self::from_catalog(TautomerCatalog::v1()?))
    }

    #[must_use]
    pub fn catalog(&self) -> &TautomerCatalog {
        &self.catalog
    }

    #[must_use]
    pub const fn options(&self) -> TautomerOptions {
        self.options
    }

    pub fn set_max_tautomers(&mut self, value: u32) {
        // RDKit✔️✔️: void setMaxTautomers(unsigned int maxTautomers) {
        // RDKit✔️✔️:   d_maxTautomers = maxTautomers;
        // RDKit✔️✔️: }
        self.options.set_max_tautomers(value);
    }

    #[must_use]
    pub const fn max_tautomers(&self) -> u32 {
        // RDKit✔️✔️: unsigned int getMaxTautomers() { return d_maxTautomers; }
        self.options.max_tautomers()
    }

    pub fn set_max_transforms(&mut self, value: u32) {
        // RDKit✔️✔️: void setMaxTransforms(unsigned int maxTransforms) {
        // RDKit✔️✔️:   d_maxTransforms = maxTransforms;
        // RDKit✔️✔️: }
        self.options.set_max_transforms(value);
    }

    #[must_use]
    pub const fn max_transforms(&self) -> u32 {
        // RDKit✔️✔️: unsigned int getMaxTransforms() { return d_maxTransforms; }
        self.options.max_transforms()
    }

    pub fn set_remove_sp3_stereo(&mut self, value: bool) {
        // RDKit✔️✔️: void setRemoveSp3Stereo(bool removeSp3Stereo) {
        // RDKit✔️✔️:   d_removeSp3Stereo = removeSp3Stereo;
        // RDKit✔️✔️: }
        self.options.set_remove_sp3_stereo(value);
    }

    #[must_use]
    pub const fn remove_sp3_stereo(&self) -> bool {
        // RDKit✔️✔️: bool getRemoveSp3Stereo() { return d_removeSp3Stereo; }
        self.options.remove_sp3_stereo()
    }

    pub fn set_remove_bond_stereo(&mut self, value: bool) {
        // RDKit✔️✔️: void setRemoveBondStereo(bool removeBondStereo) {
        // RDKit✔️✔️:   d_removeBondStereo = removeBondStereo;
        // RDKit✔️✔️: }
        self.options.set_remove_bond_stereo(value);
    }

    #[must_use]
    pub const fn remove_bond_stereo(&self) -> bool {
        // RDKit✔️✔️: bool getRemoveBondStereo() { return d_removeBondStereo; }
        self.options.remove_bond_stereo()
    }

    pub fn set_remove_isotopic_hydrogens(&mut self, value: bool) {
        // RDKit✔️✔️: void setRemoveIsotopicHs(bool removeIsotopicHs) {
        // RDKit✔️✔️:   d_removeIsotopicHs = removeIsotopicHs;
        // RDKit✔️✔️: }
        self.options.set_remove_isotopic_hydrogens(value);
    }

    #[must_use]
    pub const fn remove_isotopic_hydrogens(&self) -> bool {
        // RDKit✔️✔️: bool getRemoveIsotopicHs() { return d_removeIsotopicHs; }
        self.options.remove_isotopic_hydrogens()
    }

    pub fn set_reassign_stereo(&mut self, value: bool) {
        // RDKit✔️✔️: void setReassignStereo(bool reassignStereo) {
        // RDKit✔️✔️:   d_reassignStereo = reassignStereo;
        // RDKit✔️✔️: }
        self.options.set_reassign_stereo(value);
    }

    #[must_use]
    pub const fn reassign_stereo(&self) -> bool {
        // RDKit✔️✔️: bool getReassignStereo() { return d_reassignStereo; }
        self.options.reassign_stereo()
    }

    pub fn set_callback(&mut self, callback: Option<&'callback dyn TautomerEnumerationCallback>) {
        // RDKit✔️🔝: void setCallback(TautomerEnumeratorCallback *callback) {
        // RDKit✔️🔝:   d_callback.reset(callback);
        // RDKit✔️🔝: }
        // Borrowing makes lifetime and ownership explicit and eliminates the
        // raw-pointer ownership transfer without changing callback identity.
        self.callback = callback;
    }

    #[must_use]
    pub const fn callback(&self) -> Option<&'callback dyn TautomerEnumerationCallback> {
        // RDKit✔️✔️: TautomerEnumeratorCallback *getCallback() const { return d_callback.get(); }
        self.callback
    }

    pub fn enumerate(
        &self,
        molecule: &Molecule,
    ) -> Result<TautomerEnumeration, crate::OperationError> {
        // RDKit✔️❌: TautomerEnumeratorResult TautomerEnumerator::enumerate(
        // RDKit✔️❌:     const ROMol &mol) const;
        molecule.enumerate_tautomers_with_options(self)
    }

    #[allow(dead_code)]
    fn enumerate_deprecated_compat(
        &self,
        molecule: &Molecule,
        modified_atoms: Option<&mut BTreeSet<AtomId>>,
        modified_bonds: Option<&mut BTreeSet<BondId>>,
    ) -> Result<Vec<Molecule>, crate::OperationError> {
        // RDKit✔️✔️: std::vector<ROMOL_SPTR> TautomerEnumerator::enumerate(
        // RDKit✔️✔️:     const ROMol &mol, boost::dynamic_bitset<> *modifiedAtoms,
        // RDKit✔️✔️:     boost::dynamic_bitset<> *modifiedBonds) const {
        // RDKit✔️✔️:   TautomerEnumeratorResult tresult = enumerate(mol);
        let result = self.enumerate(molecule)?;
        // RDKit✔️✔️:   if (modifiedAtoms) {
        // RDKit✔️✔️:     *modifiedAtoms = tresult.modifiedAtoms();
        // RDKit✔️✔️:   }
        if let Some(modified_atoms) = modified_atoms {
            modified_atoms.clone_from(result.modified_atoms());
        }
        // RDKit✔️✔️:   if (modifiedBonds) {
        // RDKit✔️✔️:     *modifiedBonds = tresult.modifiedBonds();
        // RDKit✔️✔️:   }
        if let Some(modified_bonds) = modified_bonds {
            modified_bonds.clone_from(result.modified_bonds());
        }
        // RDKit✔️✔️:   return tresult.tautomers();
        // RDKit✔️✔️: }
        Ok(result.molecules())
    }

    pub fn pick_canonical(
        &self,
        result: &TautomerEnumeration,
    ) -> Result<Molecule, TautomerRunError> {
        self.pick_canonical_with(result, |molecule| Ok(score_tautomer(molecule)?.total()))
    }

    pub fn pick_canonical_with(
        &self,
        result: &TautomerEnumeration,
        score_func: impl FnMut(&Molecule) -> Result<i32, TautomerRunError>,
    ) -> Result<Molecule, TautomerRunError> {
        // RDKit✔️✔️: ROMol *TautomerEnumerator::pickCanonical(
        // RDKit✔️✔️:     const TautomerEnumeratorResult &tautRes,
        // RDKit✔️✔️:     boost::function<int(const ROMol &mol)> scoreFunc) const {
        // RDKit✔️✔️:   ROMOL_SPTR bestMol;
        // RDKit✔️✔️:   if (tautRes.d_tautomers.size() == 1) {
        // RDKit✔️✔️:     bestMol = tautRes.d_tautomers.begin()->second.tautomer;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     // Calculate score for each tautomer
        // RDKit✔️✔️:     int bestScore = std::numeric_limits<int>::min();
        // RDKit✔️✔️:     std::string bestSmiles = "";
        // RDKit✔️✔️:     for (const auto &t : tautRes.d_tautomers) {
        // RDKit✔️✔️:       auto score = scoreFunc(*t.second.tautomer);
        // RDKit✔️✔️: #ifdef VERBOSE_ENUMERATION
        // RDKit✔️✔️:       std::cerr << "  " << t.first << " " << score << std::endl;
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       if (score > bestScore) {
        // RDKit✔️✔️:         bestScore = score;
        // RDKit✔️✔️:         bestSmiles = t.first;
        // RDKit✔️✔️:         bestMol = t.second.tautomer;
        // RDKit✔️✔️:       } else if (score == bestScore) {
        // RDKit✔️✔️:         if (t.first < bestSmiles) {
        // RDKit✔️✔️:           bestSmiles = t.first;
        // RDKit✔️✔️:           bestMol = t.second.tautomer;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        pick_canonical_candidates(
            result
                .iter_with_smiles()
                .map(|(smiles, molecule)| CanonicalCandidate {
                    molecule,
                    retained_smiles: Some(smiles),
                }),
            score_func,
        )
    }

    pub fn pick_canonical_from_iterable<'a, I>(
        &self,
        tautomers: I,
    ) -> Result<Molecule, TautomerRunError>
    where
        I: IntoIterator<Item = &'a Molecule>,
        I::IntoIter: ExactSizeIterator,
    {
        self.pick_canonical_from_iterable_with(tautomers, |molecule| {
            Ok(score_tautomer(molecule)?.total())
        })
    }

    pub fn pick_canonical_from_iterable_with<'a, I>(
        &self,
        tautomers: I,
        score_func: impl FnMut(&Molecule) -> Result<i32, TautomerRunError>,
    ) -> Result<Molecule, TautomerRunError>
    where
        I: IntoIterator<Item = &'a Molecule>,
        I::IntoIter: ExactSizeIterator,
    {
        // RDKit✔️✔️: template <class Iterable,
        // RDKit✔️✔️:           typename std::enable_if<
        // RDKit✔️✔️:               !std::is_same<Iterable, TautomerEnumeratorResult>::value,
        // RDKit✔️✔️:               int>::type = 0>
        // RDKit✔️✔️: ROMol *pickCanonical(const Iterable &tautomers,
        // RDKit✔️✔️:                      boost::function<int(const ROMol &mol)> scoreFunc =
        // RDKit✔️✔️:                          TautomerScoringFunctions::scoreTautomer) const {
        // RDKit✔️✔️:   ROMOL_SPTR bestMol;
        // RDKit✔️✔️:   if (tautomers.size() == 1) {
        // RDKit✔️✔️:     bestMol = *tautomers.begin();
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     // Calculate score for each tautomer
        // RDKit✔️✔️:     int bestScore = std::numeric_limits<int>::min();
        // RDKit✔️✔️:     std::string bestSmiles = "";
        // RDKit✔️✔️:     for (const auto &t : tautomers) {
        // RDKit✔️✔️:       auto score = scoreFunc(*t);
        // RDKit✔️✔️: #ifdef VERBOSE_ENUMERATION
        // RDKit✔️✔️:       std::cerr << "  " << MolToSmiles(*t) << " " << score << std::endl;
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:       if (score > bestScore) {
        // RDKit✔️✔️:         bestScore = score;
        // RDKit✔️✔️:         bestSmiles = MolToSmiles(*t);
        // RDKit✔️✔️:         bestMol = t;
        // RDKit✔️✔️:       } else if (score == bestScore) {
        // RDKit✔️✔️:         auto smiles = MolToSmiles(*t);
        // RDKit✔️✔️:         if (smiles < bestSmiles) {
        // RDKit✔️✔️:           bestSmiles = smiles;
        // RDKit✔️✔️:           bestMol = t;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        pick_canonical_candidates(
            tautomers.into_iter().map(|molecule| CanonicalCandidate {
                molecule,
                retained_smiles: None,
            }),
            score_func,
        )
    }

    pub fn canonicalize(
        &self,
        molecule: &Molecule,
    ) -> Result<Molecule, TautomerCanonicalizationError> {
        self.canonicalize_with(molecule, |candidate| Ok(score_tautomer(candidate)?.total()))
    }

    pub fn canonicalize_with(
        &self,
        molecule: &Molecule,
        score_func: impl FnMut(&Molecule) -> Result<i32, TautomerRunError>,
    ) -> Result<Molecule, TautomerCanonicalizationError> {
        // RDKit✔️✔️: ROMol *TautomerEnumerator::canonicalize(
        // RDKit✔️✔️:     const ROMol &mol, boost::function<int(const ROMol &mol)> scoreFunc) const {
        // RDKit✔️✔️:   auto thisCopy = TautomerEnumerator(*this);
        let mut enumerator = self.clone();
        // RDKit✔️✔️:   thisCopy.setReassignStereo(false);
        enumerator.set_reassign_stereo(false);
        // RDKit✔️✔️:   auto res = thisCopy.enumerate(mol);
        let result = enumerator.enumerate(molecule).map_err(|source| {
            TautomerCanonicalizationError::Enumeration {
                source: Box::new(source),
            }
        })?;
        // RDKit✔️✔️:   if (res.empty()) {
        if result.is_empty() {
            // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
            // RDKit✔️✔️:         << "no tautomers found, returning input molecule" << std::endl;
            // RDKit✔️✔️:     return new ROMol(mol);
            return Ok(molecule.clone());
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   return pickCanonical(res, scoreFunc);
        // RDKit✔️✔️: }
        Ok(self.pick_canonical_with(&result, score_func)?)
    }

    #[cfg(test)]
    fn canonicalize_in_place_compat_with(
        &self,
        molecule: &mut Molecule,
        score_func: impl FnMut(&Molecule) -> Result<i32, TautomerRunError>,
    ) -> Result<(), TautomerCanonicalizationError> {
        // RDKit✔️✔️: void TautomerEnumerator::canonicalizeInPlace(
        // RDKit✔️✔️:     RWMol &mol, boost::function<int(const ROMol &mol)> scoreFunc) const {
        // RDKit✔️✔️:   auto thisCopy = TautomerEnumerator(*this);
        let mut enumerator = self.clone();
        // RDKit✔️✔️:   thisCopy.setReassignStereo(false);
        enumerator.set_reassign_stereo(false);
        // RDKit✔️✔️:   auto res = thisCopy.enumerate(mol);
        let result = enumerator.enumerate(molecule).map_err(|source| {
            TautomerCanonicalizationError::Enumeration {
                source: Box::new(source),
            }
        })?;
        // RDKit✔️✔️:   if (res.empty()) {
        if result.is_empty() {
            // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
            // RDKit✔️✔️:         << "no tautomers found, molecule unchanged" << std::endl;
            // RDKit✔️✔️:     return;
            return Ok(());
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   std::unique_ptr<ROMol> tmp{pickCanonical(res, scoreFunc)};
        let canonical = self.pick_canonical_with(&result, score_func)?;

        // RDKit✔️✔️:   TEST_ASSERT(tmp->getNumAtoms() == mol.getNumAtoms());
        assert_eq!(canonical.num_atoms(), molecule.num_atoms());
        // RDKit✔️✔️:   TEST_ASSERT(tmp->getNumBonds() == mol.getNumBonds());
        assert_eq!(canonical.num_bonds(), molecule.num_bonds());
        // RDKit✔️✔️:   // now copy the info from the canonical tautomer over to the input molecule
        // RDKit✔️✔️:   for (const auto tmpAtom : tmp->atoms()) {
        for (index, canonical_atom) in canonical.atoms().iter().enumerate() {
            // RDKit✔️✔️:     auto atom = mol.getAtomWithIdx(tmpAtom->getIdx());
            let atom = &mut molecule.topology_block_mut().atoms[index];
            // RDKit✔️✔️:     TEST_ASSERT(tmpAtom->getAtomicNum() == atom->getAtomicNum());
            assert_eq!(canonical_atom.atomic_number(), atom.atomic_number());
            // RDKit✔️✔️:     atom->setFormalCharge(tmpAtom->getFormalCharge());
            atom.set_formal_charge(canonical_atom.formal_charge());
            // RDKit✔️✔️:     atom->setNoImplicit(tmpAtom->getNoImplicit());
            atom.set_no_implicit(canonical_atom.no_implicit());
            // RDKit✔️✔️:     atom->setIsAromatic(tmpAtom->getIsAromatic());
            atom.set_aromatic(canonical_atom.is_aromatic());
            // RDKit✔️✔️:     atom->setNumExplicitHs(tmpAtom->getNumExplicitHs());
            atom.set_explicit_hydrogens(canonical_atom.explicit_hydrogens());
            // RDKit✔️✔️:     atom->setNumRadicalElectrons(tmpAtom->getNumRadicalElectrons());
            atom.set_radical_electrons(canonical_atom.radical_electrons());
            // RDKit✔️✔️:     atom->setChiralTag(tmpAtom->getChiralTag());
            atom.set_chiral_tag(canonical_atom.chiral_tag());
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   for (const auto tmpBond : tmp->bonds()) {
        for (index, canonical_bond) in canonical.bonds().iter().enumerate() {
            // RDKit✔️✔️:     auto bond = mol.getBondWithIdx(tmpBond->getIdx());
            let bond = &mut molecule.topology_block_mut().bonds[index];
            // RDKit✔️✔️:     TEST_ASSERT(tmpBond->getBeginAtomIdx() == bond->getBeginAtomIdx());
            assert_eq!(canonical_bond.begin(), bond.begin());
            // RDKit✔️✔️:     TEST_ASSERT(tmpBond->getEndAtomIdx() == bond->getEndAtomIdx());
            assert_eq!(canonical_bond.end(), bond.end());
            // RDKit✔️✔️:     bond->setBondType(tmpBond->getBondType());
            bond.set_order(canonical_bond.order());
            // RDKit✔️✔️:     bond->setBondDir(tmpBond->getBondDir());
            bond.set_direction(canonical_bond.direction());
            // RDKit✔️✔️:     bond->setIsAromatic(tmpBond->getIsAromatic());
            bond.set_aromatic(canonical_bond.is_aromatic());
            // RDKit✔️✔️:     bond->setIsConjugated(tmpBond->getIsConjugated());
            bond.set_conjugated(canonical_bond.is_conjugated());
            // RDKit✔️✔️:     if (tmpBond->getStereoAtoms().size() == 2) {
            if let Some(stereo_atoms) = canonical_bond.stereo_atoms() {
                // RDKit✔️✔️:       bond->setStereoAtoms(tmpBond->getStereoAtoms()[0],
                // RDKit✔️✔️:                            tmpBond->getStereoAtoms()[1]);
                bond.set_stereo_atoms(Some(stereo_atoms));
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     bond->setStereo(tmpBond->getStereo());
            bond.set_stereo(canonical_bond.stereo());
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   mol.updatePropertyCache(false);
        let valence = crate::valence::assign_valence_with_options(
            molecule,
            crate::ValenceModel::RdkitLike,
            false,
        )
        .map_err(TautomerRunError::from)?;
        molecule.derived_cache_mut().valence = Some(valence);
        // RDKit✔️✔️: }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct CanonicalCandidate<'a> {
    molecule: &'a Molecule,
    retained_smiles: Option<&'a str>,
}

fn pick_canonical_candidates<'a, I>(
    mut candidates: I,
    mut score_func: impl FnMut(&Molecule) -> Result<i32, TautomerRunError>,
) -> Result<Molecule, TautomerRunError>
where
    I: ExactSizeIterator<Item = CanonicalCandidate<'a>>,
{
    let selected = if candidates.len() == 1 {
        candidates
            .next()
            .map(|candidate| candidate.molecule)
            .ok_or(TautomerRunError::NoCanonicalTautomer)?
    } else {
        let mut best_score = i32::MIN;
        let mut best_smiles = String::new();
        let mut best_molecule = None;
        for candidate in candidates {
            let score = score_func(candidate.molecule)?;
            if score > best_score {
                best_score = score;
                best_smiles = canonical_candidate_smiles(candidate)?;
                best_molecule = Some(candidate.molecule);
            } else if score == best_score {
                let smiles = canonical_candidate_smiles(candidate)?;
                if smiles < best_smiles {
                    best_smiles = smiles;
                    best_molecule = Some(candidate.molecule);
                }
            }
        }
        best_molecule.ok_or(TautomerRunError::NoCanonicalTautomer)?
    };

    // RDKit✔️✔️:   ROMol *res = new ROMol(*bestMol);
    // RDKit✔️✔️:   static const bool cleanIt = true;
    // RDKit✔️✔️:   static const bool force = true;
    // RDKit✔️✔️:   MolOps::assignStereochemistry(*res, cleanIt, force);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // The shared stereo implementation always executes here, reproducing the
    // source's `force=true` call shape for the tautomer canonicalization
    // boundary. This does not widen the claim to unrelated coordinate-driven
    // `assignStereochemistry()` entrypoints.
    let mut result = selected.clone();
    crate::smiles::assign_stereochemistry_cleanup_subset(&mut result, true)?;
    Ok(result)
}

fn canonical_candidate_smiles(
    candidate: CanonicalCandidate<'_>,
) -> Result<String, TautomerRunError> {
    match candidate.retained_smiles {
        Some(smiles) => Ok(smiles.to_owned()),
        None => {
            Ok(MoleculeReadParts::from_molecule(candidate.molecule).canonical_isomeric_smiles()?)
        }
    }
}

pub(crate) fn plan_tautomer_initialization(
    read: MoleculeReadParts<'_>,
) -> Result<TautomerInitializationPlan, TautomerInitializationError> {
    // RDKit✔️✔️:   TautomerEnumeratorResult res;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const std::vector<TautomerTransform> &transforms =
    // RDKit✔️✔️:       tautparams->getTransforms();
    // The result and transform references are owned by later state-machine
    // stages. This helper reproduces the complete input/candidate block only.

    // RDKit✔️✔️:   // Enumerate all possible tautomers and return them as a vector.
    // RDKit✔️✔️:   // smi is the input molecule SMILES
    // RDKit✔️✔️:   std::string smi = MolToSmiles(mol, true);
    let canonical_smiles = read.canonical_isomeric_smiles()?;
    plan_tautomer_initialization_with_canonical_smiles(read, canonical_smiles)
}

fn plan_tautomer_initialization_with_canonical_smiles(
    read: MoleculeReadParts<'_>,
    canonical_smiles: String,
) -> Result<TautomerInitializationPlan, TautomerInitializationError> {
    // RDKit✔️✔️:   // taut is a copy of the input molecule
    // RDKit✔️✔️:   ROMOL_SPTR taut(new ROMol(mol));
    // Copy-on-write branch materialization is deliberately deferred to
    // `MultiMoleculeOpParts`; this plan contains only source-derived updates.

    // RDKit✔️✔️:   // do whatever sanitization bits are required
    // RDKit✔️✔️:   if (taut->needsUpdatePropertyCache()) {
    // RDKit✔️✔️:     taut->updatePropertyCache(false);
    // RDKit✔️✔️:   }
    let valence_update = if crate::valence::needs_update_property_cache(read) {
        Some(read.assign_valence_with_options(crate::ValenceModel::RdkitLike, false)?)
    } else {
        None
    };

    // RDKit✔️✔️:   if (!taut->getRingInfo()->isSymmSssr()) {
    // RDKit✔️✔️:     MolOps::symmetrizeSSSR(*taut);
    // RDKit✔️✔️:   }
    let rings_update = match read.derived_cache().rings.as_ref() {
        Some(rings) if rings.is_symm_sssr() => None,
        _ => Some(read.symmetrize_sssr()?),
    };

    let effective_valence = valence_update
        .as_ref()
        .or(read.derived_cache().valence.as_ref());
    let effective_rings = rings_update
        .as_ref()
        .or(read.derived_cache().rings.as_ref())
        .expect("SymmSSSR is either retained or computed above");

    // RDKit✔️✔️:   // Create a kekulized form of the molecule to match the SMARTS against.
    // RDKit✔️✔️:   // canonical=true is required so that tautomer deduplication is independent
    // RDKit✔️✔️:   // of atom ordering in the molecule.
    // RDKit✔️✔️:   RWMOL_SPTR kekulized(new RWMol(*taut));
    // RDKit✔️✔️:   MolOps::Kekulize(*kekulized, false, true);
    let kekulize_assignment = read.kekulize_assignment_with_valence(
        Some(effective_rings),
        effective_valence,
        false,
        true,
        100,
    )?;

    // RDKit✔️✔️:   res.d_tautomers = {{smi, Tautomer(taut, kekulized, 0, 0)}};
    // RDKit✔️✔️:   res.d_modifiedAtoms.resize(mol.getNumAtoms());
    // RDKit✔️✔️:   res.d_modifiedBonds.resize(mol.getNumBonds());
    // `into_candidate_map()` performs the same one-entry ordered-map
    // insertion once the operation layer supplies its validated branch
    // handles. Typed modified sets are sparse, so their valid source index
    // bounds are retained explicitly instead of allocating false-filled bits.
    let mut target_derived_cache = read.derived_cache().clone();
    if let Some(valence) = &valence_update {
        target_derived_cache.valence = Some(valence.clone());
    }
    if let Some(rings) = &rings_update {
        target_derived_cache.rings = Some(rings.clone());
    }

    Ok(TautomerInitializationPlan {
        canonical_smiles,
        valence_update,
        rings_update,
        target_derived_cache,
        kekulize_assignment,
        atom_count: read.num_atoms(),
        bond_count: read.num_bonds(),
    })
}

fn molecule_snapshot_from_read_parts(
    read: MoleculeReadParts<'_>,
) -> Result<Molecule, crate::InvariantError> {
    Molecule::from_operation_blocks(
        read.topology().clone(),
        read.coordinates().clone(),
        read.properties().clone(),
        read.derived_cache().clone(),
        read.capabilities(),
    )
}

pub(crate) fn find_tautomer_transform_matches(
    candidate_read: MoleculeReadParts<'_>,
    transform: &TautomerTransform,
) -> Result<Vec<crate::SubstructMatchResult>, TautomerRunError> {
    let candidate = molecule_snapshot_from_read_parts(candidate_read)?;
    crate::try_get_substruct_matches_with_params(
        &candidate,
        transform.query(),
        &crate::SubstructMatchParams::default(),
    )
    .map_err(TautomerRunError::Substructure)
}

pub(crate) fn evaluate_tautomer_callback(
    callback: &dyn TautomerEnumerationCallback,
    source_read: MoleculeReadParts<'_>,
    candidate_reads: Vec<(String, MoleculeReadParts<'_>)>,
    status: TautomerEnumerationStatus,
    modified_atoms: &BTreeSet<AtomId>,
    modified_bonds: &BTreeSet<BondId>,
) -> Result<bool, TautomerRunError> {
    let source = molecule_snapshot_from_read_parts(source_read)?;
    let entries = candidate_reads
        .into_iter()
        .map(|(canonical_smiles, read)| {
            molecule_snapshot_from_read_parts(read).map(|molecule| (canonical_smiles, molecule))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let result = TautomerEnumeration::from_ordered_entries(
        entries,
        status,
        modified_atoms.clone(),
        modified_bonds.clone(),
    );
    Ok(callback.should_continue(&source, &result))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn apply_tautomer_transform_match(
    source_read: MoleculeReadParts<'_>,
    candidate_read: MoleculeReadParts<'_>,
    transform: &TautomerTransform,
    match_result: &crate::SubstructMatchResult,
    current_modified_atoms: &BTreeSet<AtomId>,
    current_modified_bonds: &BTreeSet<BondId>,
    existing_smiles: &BTreeSet<String>,
    options: TautomerOptions,
) -> Result<TautomerTransformAttempt, TautomerTransformApplicationError> {
    // RDKit✔️❌:           // Create a copy of in the input molecule so we can modify it
    // RDKit✔️❌:           // Use kekule form so bonds are explicitly single/double instead of
    // RDKit✔️❌:           // aromatic
    // RDKit✔️❌:           RWMOL_SPTR product(new RWMol(*kmol));
    // The operation-local scratch value is built only inside this narrowed
    // helper and never crosses back into an operation body. This currently
    // deep-clones block state before the nested shared sanitize operation, so
    // the complexity axis remains explicitly below the source shared-pointer
    // clone path even though the behavior axis is source-validated.
    let source = Molecule::from_operation_blocks(
        source_read.topology().clone(),
        source_read.coordinates().clone(),
        source_read.properties().clone(),
        source_read.derived_cache().clone(),
        source_read.capabilities(),
    )?;
    let mut topology = candidate_read.topology().clone();
    let candidate = Molecule::from_operation_blocks(
        candidate_read.topology().clone(),
        candidate_read.coordinates().clone(),
        candidate_read.properties().clone(),
        candidate_read.derived_cache().clone(),
        candidate_read.capabilities(),
    )?;

    if match_result.atom_mapping.len() != transform.query().num_atoms() {
        return Err(TautomerTransformApplicationError::AtomMappingCount {
            expected: transform.query().num_atoms(),
            actual: match_result.atom_mapping.len(),
        });
    }
    if match_result.bond_mapping.len() != transform.query().num_bonds() {
        return Err(TautomerTransformApplicationError::BondMappingCount {
            expected: transform.query().num_bonds(),
            actual: match_result.bond_mapping.len(),
        });
    }
    for &atom in &match_result.atom_mapping {
        if atom >= topology.atoms.len() {
            return Err(TautomerTransformApplicationError::AtomOutOfRange {
                atom,
                atom_count: topology.atoms.len(),
            });
        }
    }
    for &bond in &match_result.bond_mapping {
        if bond >= topology.bonds.len() {
            return Err(TautomerTransformApplicationError::BondOutOfRange {
                bond,
                bond_count: topology.bonds.len(),
            });
        }
    }

    // RDKit✔️❌:           // Remove a hydrogen from the first matched atom and add one to the
    // RDKit✔️❌:           // last
    // RDKit✔️❌:           int firstIdx = match.front().second;
    // RDKit✔️❌:           int lastIdx = match.back().second;
    // RDKit✔️❌:           Atom *first = product->getAtomWithIdx(firstIdx);
    // RDKit✔️❌:           Atom *last = product->getAtomWithIdx(lastIdx);
    // RDKit✔️❌:           res.d_modifiedAtoms.set(firstIdx);
    // RDKit✔️❌:           res.d_modifiedAtoms.set(lastIdx);
    let first = AtomId::new(match_result.atom_mapping[0]);
    let last = AtomId::new(match_result.atom_mapping[match_result.atom_mapping.len() - 1]);
    let mut modified_atoms = current_modified_atoms.clone();
    let mut modified_bonds = current_modified_bonds.clone();
    modified_atoms.insert(first);
    modified_atoms.insert(last);

    // RDKit✔️❌:           first->setNumExplicitHs(
    // RDKit✔️❌:               std::max(0, static_cast<int>(first->getTotalNumHs()) - 1));
    // RDKit✔️❌:           last->setNumExplicitHs(last->getTotalNumHs() + 1);
    let first_hydrogens = crate::valence::total_num_hydrogens(&candidate, first, false)?;
    let last_hydrogens = crate::valence::total_num_hydrogens(&candidate, last, false)?;
    let first_explicit = u8::try_from(first_hydrogens.saturating_sub(1)).map_err(|_| {
        TautomerTransformApplicationError::HydrogenCountOutOfRange {
            atom: first,
            count: first_hydrogens.saturating_sub(1),
        }
    })?;
    let last_total = last_hydrogens.checked_add(1).ok_or(
        TautomerTransformApplicationError::HydrogenCountOutOfRange {
            atom: last,
            count: u32::MAX,
        },
    )?;
    let last_explicit = u8::try_from(last_total).map_err(|_| {
        TautomerTransformApplicationError::HydrogenCountOutOfRange {
            atom: last,
            count: last_total,
        }
    })?;
    topology.atoms[first.index()].set_explicit_hydrogens(first_explicit);
    topology.atoms[last.index()].set_explicit_hydrogens(last_explicit);

    // RDKit✔️❌:           // Remove any implicit hydrogens from the first and last atoms
    // RDKit✔️❌:           // now we have set the count explicitly
    // RDKit✔️❌:           first->setNoImplicit(true);
    // RDKit✔️❌:           last->setNoImplicit(true);
    topology.atoms[first.index()].set_no_implicit(true);
    topology.atoms[last.index()].set_no_implicit(true);

    // RDKit✔️❌:           // Adjust bond orders
    // RDKit✔️❌:           unsigned int bi = 0;
    // RDKit✔️❌:           for (size_t i = 0; i < transform.Mol->getNumBonds(); ++i) {
    // RDKit✔️❌:             const auto tbond = transform.Mol->getBondWithIdx(i);
    // RDKit✔️❌:             Bond *bond = product->getBondBetweenAtoms(
    // RDKit✔️❌:                 match[tbond->getBeginAtomIdx()].second,
    // RDKit✔️❌:                 match[tbond->getEndAtomIdx()].second);
    // RDKit✔️❌:             ASSERT_INVARIANT(bond, "required bond not found");
    for (query_bond, &target_bond) in match_result.bond_mapping.iter().enumerate() {
        let bond_id = BondId::new(target_bond);
        let bond = &mut topology.bonds[target_bond];
        // RDKit✔️❌:             // check if bonds is specified in tautomer.in file
        // RDKit✔️❌:             if (!transform.BondTypes.empty()) {
        // RDKit✔️❌:               bond->setBondType(transform.BondTypes[bi]);
        // RDKit✔️❌:               ++bi;
        // RDKit✔️❌:             } else {
        if let Some(&bond_type) = transform.bond_types().get(query_bond) {
            bond.set_order(bond_type);
        } else {
            // RDKit✔️❌:               Bond::BondType bondtype = bond->getBondType();
            // RDKit✔️❌:               if (bondtype == Bond::SINGLE) {
            // RDKit✔️❌:                 bond->setBondType(Bond::DOUBLE);
            // RDKit✔️❌:               }
            // RDKit✔️❌:               if (bondtype == Bond::DOUBLE) {
            // RDKit✔️❌:                 bond->setBondType(Bond::SINGLE);
            // RDKit✔️❌:               }
            let bond_type = bond.order();
            if bond_type == BondOrder::Single {
                bond.set_order(BondOrder::Double);
            }
            if bond_type == BondOrder::Double {
                bond.set_order(BondOrder::Single);
            }
        }
        // RDKit✔️❌:             }
        // RDKit✔️❌:             res.d_modifiedBonds.set(bond->getIdx());
        modified_bonds.insert(bond_id);
        // RDKit✔️❌:           }
    }

    // RDKit✔️❌:           // TODO adjust charges
    // RDKit✔️❌:           if (!transform.Charges.empty()) {
    // RDKit✔️❌:             unsigned int ci = 0;
    // RDKit✔️❌:             for (const auto &pair : match) {
    // RDKit✔️❌:               Atom *atom = product->getAtomWithIdx(pair.second);
    // RDKit✔️❌:               atom->setFormalCharge(atom->getFormalCharge() +
    // RDKit✔️❌:                                     transform.Charges[ci++]);
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    for (&target_atom, &delta) in match_result.atom_mapping.iter().zip(transform.charges()) {
        let atom_id = AtomId::new(target_atom);
        let charge = i32::from(topology.atoms[target_atom].formal_charge()) + delta;
        let charge = i8::try_from(charge).map_err(|_| {
            TautomerTransformApplicationError::FormalChargeOutOfRange {
                atom: atom_id,
                charge,
            }
        })?;
        topology.atoms[target_atom].set_formal_charge(charge);
    }

    let product = Molecule::from_operation_blocks(
        topology,
        candidate_read.coordinates().clone(),
        candidate_read.properties().clone(),
        candidate_read.derived_cache().clone(),
        candidate_read.capabilities(),
    )?;
    let sanitize_ops = crate::SanitizeOps::KEKULIZE
        | crate::SanitizeOps::SET_AROMATICITY
        | crate::SanitizeOps::SET_CONJUGATION
        | crate::SanitizeOps::SET_HYBRIDIZATION
        | crate::SanitizeOps::ADJUST_HYDROGENS;

    // RDKit✔️❌:           unsigned int failedOp;
    // RDKit✔️❌:           try {
    // RDKit✔️❌:             MolOps::sanitizeMol(*product, failedOp,
    // RDKit✔️❌:                                 MolOps::SANITIZE_KEKULIZE |
    // RDKit✔️❌:                                     MolOps::SANITIZE_SETAROMATICITY |
    // RDKit✔️❌:                                     MolOps::SANITIZE_SETCONJUGATION |
    // RDKit✔️❌:                                     MolOps::SANITIZE_SETHYBRIDIZATION |
    // RDKit✔️❌:                                     MolOps::SANITIZE_ADJUSTHS);
    // RDKit✔️❌:           } catch (const KekulizeException &) {
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    let mut product = match product.sanitize_with_ops(sanitize_ops) {
        Ok(product) => product,
        Err(error) => {
            if matches!(
                &error,
                crate::OperationError::Sanitize { source, .. }
                    if source.step == crate::SanitizeStep::Kekulize
            ) {
                return Ok(TautomerTransformAttempt::RecoverableKekulizeFailure {
                    modified_atoms,
                    modified_bonds,
                });
            }
            return Err(TautomerTransformApplicationError::Sanitize(error));
        }
    };
    // RDKit✔️❌:           setTautomerStereoAndIsoHs(mol, *product, res);
    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut product,
        &modified_atoms,
        &modified_bonds,
        options,
    )?;
    // RDKit✔️❌:           tsmiles = MolToSmiles(*product, true);
    let canonical_smiles =
        MoleculeReadParts::from_molecule(&product).canonical_isomeric_smiles()?;

    // RDKit✔️❌:           if (res.d_tautomers.find(tsmiles) != res.d_tautomers.end()) {
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    if existing_smiles.contains(&canonical_smiles) {
        return Ok(TautomerTransformAttempt::Duplicate {
            canonical_smiles,
            modified_atoms,
            modified_bonds,
        });
    }

    // RDKit✔️❌:           // in addition to the above transformations, sanitization may modify
    // RDKit✔️❌:           // bonds, e.g. Cc1nc2ccccc2[nH]1
    // RDKit✔️❌:           for (size_t i = 0; i < mol.getNumBonds(); i++) {
    // RDKit✔️❌:             auto molBondType = mol.getBondWithIdx(i)->getBondType();
    // RDKit✔️❌:             auto tautBondType = product->getBondWithIdx(i)->getBondType();
    // RDKit✔️❌:             if (molBondType != tautBondType && !res.d_modifiedBonds.test(i)) {
    // RDKit✔️❌:               res.d_modifiedBonds.set(i);
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    for (index, (source_bond, product_bond)) in
        source.bonds().iter().zip(product.bonds()).enumerate()
    {
        if source_bond.order() != product_bond.order() {
            modified_bonds.insert(BondId::new(index));
        }
    }

    // RDKit✔️❌:           RWMOL_SPTR kekulized_product(new RWMol(*product));
    // RDKit✔️❌:           // canonical=true for order-independent tautomer deduplication
    // RDKit✔️❌:           MolOps::Kekulize(*kekulized_product, false, true);
    let product_read = MoleculeReadParts::from_molecule(&product);
    let owned_rings;
    let rings = match product_read.derived_cache().rings.as_ref() {
        Some(rings) => rings,
        None => {
            owned_rings = product_read.symmetrize_sssr()?;
            &owned_rings
        }
    };
    let kekulize_assignment = product_read.kekulize_assignment_with_valence(
        Some(rings),
        product_read.derived_cache().valence.as_ref(),
        false,
        true,
        100,
    )?;

    // RDKit✔️❌:           res.d_tautomers[tsmiles] = Tautomer(
    // RDKit✔️❌:               std::move(product), std::move(kekulized_product),
    // RDKit✔️❌:               res.d_modifiedAtoms.count(), res.d_modifiedBonds.count());
    Ok(TautomerTransformAttempt::Product(TautomerProductPlan {
        topology: product.topology_block().clone(),
        properties: product.properties().clone(),
        derived_cache: product.derived_cache().clone(),
        canonical_smiles,
        kekulize_assignment,
        modified_atoms,
        modified_bonds,
    }))
}

pub(crate) fn expand_tautomer_candidates_in_source_order<M: Clone, E>(
    state: &mut TautomerExpansionState<M>,
    transforms: &[TautomerTransform],
    options: TautomerOptions,
    mut find_matches: impl FnMut(&M, &TautomerTransform) -> Result<Vec<crate::SubstructMatchResult>, E>,
    mut apply_match: impl FnMut(
        &M,
        &TautomerTransform,
        &crate::SubstructMatchResult,
        &BTreeSet<AtomId>,
        &BTreeSet<BondId>,
        &BTreeSet<String>,
    ) -> Result<TautomerExpansionAttempt<M>, E>,
    mut callback: impl FnMut(&TautomerExpansionState<M>) -> Result<bool, E>,
) -> Result<TautomerExpansionPass, TautomerExpansionError<E>> {
    // RDKit✔️✔️:   bool completed = false;
    // RDKit✔️✔️:   bool bailOut = false;
    // RDKit✔️✔️:   unsigned int nTransforms = 0;
    // `completed` belongs to the post-expansion pruning/rekeying stage. This
    // helper owns one exact ordered expansion pass and retains the other two
    // source variables as its return flag and cumulative state counter.
    let mut bail_out = false;

    // RDKit✔️✔️:   while (!completed && !bailOut) {
    // RDKit✔️✔️:     // std::map automatically sorts res.d_tautomers into alphabetical order
    // RDKit✔️✔️:     // (SMILES)
    // RDKit✔️✔️:     for (auto &smilesTautomerPair : res.d_tautomers) {
    // A BTreeMap range cursor, rather than a snapshot of keys, reproduces
    // std::map iterator behavior when a transform inserts during traversal:
    // later keys are visible in this pass and earlier keys wait for the next.
    let mut previous_key: Option<String> = None;
    loop {
        let next_key = match previous_key.as_deref() {
            Some(key) => state
                .candidates
                .range::<str, _>((std::ops::Bound::Excluded(key), std::ops::Bound::Unbounded))
                .next()
                .map(|(key, _)| key.clone()),
            None => state.candidates.keys().next().cloned(),
        };
        let Some(current_key) = next_key else {
            break;
        };
        previous_key = Some(current_key.clone());

        // RDKit✔️✔️:       if (smilesTautomerPair.second.d_done) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        let (done, kekulized) = {
            let candidate = &state.candidates[&current_key];
            (candidate.done, candidate.kekulized.clone())
        };
        if done {
            continue;
        }
        let kekulized =
            kekulized.ok_or_else(|| TautomerExpansionError::MissingKekulizedBranch {
                canonical_smiles: current_key.clone(),
            })?;

        // RDKit✔️✔️:       // tautomer not yet done
        // RDKit✔️✔️:       for (const auto &transform : transforms) {
        for transform in transforms {
            // RDKit✔️✔️:         if (bailOut) {
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         }
            if bail_out {
                break;
            }

            // RDKit✔️✔️:         // kmol is the kekulized version of the tautomer
            // RDKit✔️✔️:         const auto &kmol = smilesTautomerPair.second.kekulized;
            // RDKit✔️✔️:         std::vector<MatchVectType> matches;
            // RDKit✔️✔️:         unsigned int matched =
            // RDKit✔️✔️:             SubstructMatch(*kmol, *(transform.Mol), matches);
            let matches =
                find_matches(&kekulized, transform).map_err(TautomerExpansionError::Backend)?;

            // RDKit✔️✔️:         if (!matched) {
            // RDKit✔️✔️:           continue;
            // RDKit✔️✔️:         }
            if matches.is_empty() {
                continue;
            }

            // RDKit✔️✔️:         ++nTransforms;
            state.num_transforms = state.num_transforms.wrapping_add(1);

            // RDKit✔️✔️:         // loop over transform matches
            // RDKit✔️✔️:         for (const auto &match : matches) {
            for matched in &matches {
                // RDKit✔️✔️:           if (nTransforms >= d_maxTransforms) {
                // RDKit✔️✔️:             res.d_status =
                // RDKit✔️✔️:                 TautomerEnumeratorStatus::MaxTransformsReached;
                // RDKit✔️✔️:             bailOut = true;
                // RDKit✔️✔️:           } else if (res.d_tautomers.size() >= d_maxTautomers) {
                // RDKit✔️✔️:             res.d_status =
                // RDKit✔️✔️:                 TautomerEnumeratorStatus::MaxTautomersReached;
                // RDKit✔️✔️:             bailOut = true;
                // RDKit✔️✔️:           } else if (d_callback.get() &&
                // RDKit✔️✔️:                      !(*d_callback)(mol, res)) {
                // RDKit✔️✔️:             res.d_status = TautomerEnumeratorStatus::Canceled;
                // RDKit✔️✔️:             bailOut = true;
                // RDKit✔️✔️:           }
                if state.num_transforms >= options.max_transforms() {
                    state.status = TautomerEnumerationStatus::MaxTransformsReached;
                    bail_out = true;
                } else if state.candidates.len() >= options.max_tautomers() as usize {
                    state.status = TautomerEnumerationStatus::MaxTautomersReached;
                    bail_out = true;
                } else if !callback(state).map_err(TautomerExpansionError::Backend)? {
                    state.status = TautomerEnumerationStatus::Canceled;
                    bail_out = true;
                }

                // RDKit✔️✔️:           if (bailOut) {
                // RDKit✔️✔️:             break;
                // RDKit✔️✔️:           }
                if bail_out {
                    break;
                }

                let existing_smiles = state.candidates.keys().cloned().collect::<BTreeSet<_>>();
                let attempt = apply_match(
                    &kekulized,
                    transform,
                    matched,
                    &state.modified_atoms,
                    &state.modified_bonds,
                    &existing_smiles,
                )
                .map_err(TautomerExpansionError::Backend)?;

                match attempt {
                    TautomerExpansionAttempt::RecoverableKekulizeFailure {
                        modified_atoms,
                        modified_bonds,
                    } => {
                        // RDKit sets the endpoint and directly edited bond bits
                        // before either source `continue` branch.
                        state.modified_atoms = modified_atoms;
                        state.modified_bonds = modified_bonds;
                    }
                    TautomerExpansionAttempt::Duplicate {
                        canonical_smiles: _,
                        modified_atoms,
                        modified_bonds,
                    } => {
                        // RDKit sets the endpoint and directly edited bond bits
                        // before either source `continue` branch.
                        state.modified_atoms = modified_atoms;
                        state.modified_bonds = modified_bonds;
                    }
                    TautomerExpansionAttempt::Product(product) => {
                        state.modified_atoms = product.modified_atoms;
                        state.modified_bonds = product.modified_bonds;
                        let canonical_smiles = product.canonical_smiles;

                        // RDKit✔️✔️:           res.d_tautomers[tsmiles] = Tautomer(
                        // RDKit✔️✔️:               std::move(product),
                        // RDKit✔️✔️:               std::move(kekulized_product),
                        // RDKit✔️✔️:               res.d_modifiedAtoms.count(),
                        // RDKit✔️✔️:               res.d_modifiedBonds.count());
                        if state.candidates.contains_key(&canonical_smiles) {
                            return Err(TautomerExpansionError::DuplicateProductKey {
                                canonical_smiles,
                            });
                        }
                        state.candidates.insert(
                            canonical_smiles,
                            TautomerCandidate {
                                tautomer: Some(product.tautomer),
                                kekulized: Some(product.kekulized),
                                num_modified_atoms: state.modified_atoms.len(),
                                num_modified_bonds: state.modified_bonds.len(),
                                done: false,
                            },
                        );
                    }
                }
            }
        }

        // RDKit✔️✔️:       smilesTautomerPair.second.d_done = true;
        // RDKit✔️✔️:     }
        state
            .candidates
            .get_mut(&current_key)
            .expect("the current ordered-map entry cannot be removed during expansion")
            .mark_done();
    }

    Ok(TautomerExpansionPass {
        bailed_out: bail_out,
    })
}

pub(crate) fn prune_and_rekey_tautomer_candidates_in_source_order<M, E>(
    state: &mut TautomerExpansionState<M>,
    options: TautomerOptions,
    mut bail_out: bool,
    mut set_stereo_and_isotopic_hydrogens: impl FnMut(
        &mut M,
        &BTreeSet<AtomId>,
        &BTreeSet<BondId>,
    ) -> Result<bool, E>,
    mut canonical_isomeric_smiles: impl FnMut(&M) -> Result<String, E>,
) -> Result<TautomerPruningPass, TautomerPruningError<E>> {
    // RDKit✔️✔️:     completed = true;
    // RDKit✔️✔️:     size_t maxNumModifiedAtoms = res.d_modifiedAtoms.count();
    // RDKit✔️✔️:     size_t maxNumModifiedBonds = res.d_modifiedBonds.count();
    let mut completed = true;
    let max_num_modified_atoms = state.modified_atoms.len();
    let max_num_modified_bonds = state.modified_bonds.len();

    // RDKit✔️✔️:     for (auto it = res.d_tautomers.begin(); it != res.d_tautomers.end();) {
    let mut current_key = state.candidates.keys().next().cloned();
    while let Some(key) = current_key {
        let (done, num_modified_atoms, num_modified_bonds) = {
            let candidate = &state.candidates[&key];
            (
                candidate.done,
                candidate.num_modified_atoms,
                candidate.num_modified_bonds,
            )
        };

        // RDKit✔️✔️:       auto &taut = it->second;
        // RDKit✔️✔️:       if (!taut.d_done) {
        // RDKit✔️✔️:         completed = false;
        // RDKit✔️✔️:       }
        if !done {
            completed = false;
        }

        // RDKit✔️✔️:       if ((taut.d_numModifiedAtoms < maxNumModifiedAtoms ||
        // RDKit✔️✔️:            taut.d_numModifiedBonds < maxNumModifiedBonds) &&
        // RDKit✔️✔️:           setTautomerStereoAndIsoHs(mol, *taut.tautomer, res)) {
        let needs_stereo_update = num_modified_atoms < max_num_modified_atoms
            || num_modified_bonds < max_num_modified_bonds;
        let stereo_changed = if needs_stereo_update {
            let candidate = state
                .candidates
                .get_mut(&key)
                .expect("the current ordered-map entry exists while pruning");
            let tautomer = candidate.tautomer.as_mut().ok_or_else(|| {
                TautomerPruningError::MissingTautomerBranch {
                    canonical_smiles: key.clone(),
                }
            })?;
            set_stereo_and_isotopic_hydrogens(
                tautomer,
                &state.modified_atoms,
                &state.modified_bonds,
            )
            .map_err(TautomerPruningError::Backend)?
        } else {
            false
        };

        if stereo_changed {
            let new_key = {
                let candidate = &state.candidates[&key];
                let tautomer = candidate.tautomer.as_ref().ok_or_else(|| {
                    TautomerPruningError::MissingTautomerBranch {
                        canonical_smiles: key.clone(),
                    }
                })?;
                canonical_isomeric_smiles(tautomer).map_err(TautomerPruningError::Backend)?
            };
            // RDKit✔️✔️:         Tautomer tautStored = std::move(taut);
            // RDKit✔️✔️:         it = res.d_tautomers.erase(it);
            let mut candidate = state
                .candidates
                .remove(&key)
                .expect("the current ordered-map entry exists while rekeying");
            let next_after_erased = state
                .candidates
                .range::<str, _>((
                    std::ops::Bound::Excluded(key.as_str()),
                    std::ops::Bound::Unbounded,
                ))
                .next()
                .map(|(next_key, _)| next_key.clone());

            // RDKit✔️✔️:         tautStored.d_numModifiedAtoms = maxNumModifiedAtoms;
            // RDKit✔️✔️:         tautStored.d_numModifiedBonds = maxNumModifiedBonds;
            candidate.update_modified_counts(max_num_modified_atoms, max_num_modified_bonds);

            // RDKit✔️✔️:         auto insertRes = res.d_tautomers.insert(std::make_pair(
            // RDKit✔️✔️:             MolToSmiles(*tautStored.tautomer), std::move(tautStored)));
            // RDKit✔️✔️:         if (insertRes.second) {
            // RDKit✔️✔️:           it = insertRes.first;
            // RDKit✔️✔️:         }
            if state.candidates.contains_key(&new_key) {
                current_key = next_after_erased;
            } else {
                state.candidates.insert(new_key.clone(), candidate);
                current_key = Some(new_key);
            }
        } else {
            // RDKit✔️✔️:       } else {
            // RDKit✔️✔️:         ++it;
            // RDKit✔️✔️:       }
            current_key = state
                .candidates
                .range::<str, _>((
                    std::ops::Bound::Excluded(key.as_str()),
                    std::ops::Bound::Unbounded,
                ))
                .next()
                .map(|(next_key, _)| next_key.clone());
        }
        // RDKit✔️✔️:     }
    }

    // RDKit✔️✔️:     if (bailOut && res.d_tautomers.size() < d_maxTautomers &&
    // RDKit✔️✔️:         res.d_status == TautomerEnumeratorStatus::MaxTautomersReached) {
    // RDKit✔️✔️:       res.d_status = TautomerEnumeratorStatus::Completed;
    // RDKit✔️✔️:       bailOut = false;
    // RDKit✔️✔️:     }
    if bail_out
        && state.candidates.len() < options.max_tautomers() as usize
        && state.status == TautomerEnumerationStatus::MaxTautomersReached
    {
        state.status = TautomerEnumerationStatus::Completed;
        bail_out = false;
    }

    Ok(TautomerPruningPass {
        completed,
        bailed_out: bail_out,
    })
}

pub(crate) fn materialize_tautomer_candidates_in_source_order<M>(
    candidates: SmilesTautomerMap<M>,
) -> Result<Vec<(String, M)>, TautomerEnumerationError> {
    // RDKit✔️✔️:   res.fillTautomersItVec();
    // RDKit✔️✔️:   void fillTautomersItVec() {
    // RDKit✔️✔️:     for (auto it = d_tautomers.begin(); it != d_tautomers.end(); ++it) {
    // RDKit✔️✔️:       d_tautomersItVec.push_back(it);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut entries = Vec::with_capacity(candidates.len());
    for (canonical_smiles, candidate) in candidates {
        let tautomer = candidate.tautomer.ok_or_else(|| {
            TautomerEnumerationError::MissingCandidateMolecule {
                canonical_smiles: canonical_smiles.clone(),
            }
        })?;
        entries.push((canonical_smiles, tautomer));
    }
    Ok(entries)
}

fn is_stereo_beyond_any(stereo: BondStereo) -> bool {
    matches!(
        stereo,
        BondStereo::Z
            | BondStereo::E
            | BondStereo::Cis
            | BondStereo::Trans
            | BondStereo::AtropCw
            | BondStereo::AtropCcw
    )
}

fn tautomer_bond_is_ring(molecule: &Molecule, bond_id: BondId) -> Result<bool, RingFindingError> {
    let bond = &molecule.bonds()[bond_id.index()];
    if bond.order() != BondOrder::Double {
        return Ok(false);
    }
    let owned_rings;
    let rings = match molecule.derived_cache().rings.as_ref() {
        Some(rings) if rings.is_find_fast_or_better() => rings,
        _ => {
            owned_rings = crate::rings::fast_find_rings_from_parts(
                molecule.num_atoms(),
                molecule.bonds(),
                &molecule.topology_block().adjacency,
            )?;
            &owned_rings
        }
    };
    Ok(rings.num_bond_rings(bond_id) != 0
        || (rings.num_atom_rings(bond.begin()) != 0 && rings.num_atom_rings(bond.end()) != 0))
}

fn set_tautomer_stereo_and_isotopic_hydrogens(
    source: &Molecule,
    tautomer: &mut Molecule,
    modified_atoms: &BTreeSet<AtomId>,
    modified_bonds: &BTreeSet<BondId>,
    options: TautomerOptions,
) -> Result<bool, TautomerStereoTransitionError> {
    // RDKit✔️❌: bool TautomerEnumerator::setTautomerStereoAndIsoHs(
    // RDKit✔️❌:     const ROMol &mol, ROMol &taut, const TautomerEnumeratorResult &res) const {
    // RDKit✔️❌:   bool modified = false;
    // The source transition is reproduced below and validated through the
    // complete tautomer stereo matrix. BTreeSet membership is O(log n), versus
    // the source dynamic_bitset's O(1), so the complexity axis is intentionally
    // not claimed equivalent.
    if tautomer.num_atoms() != source.num_atoms() {
        return Err(TautomerStereoTransitionError::AtomCountMismatch {
            expected: source.num_atoms(),
            actual: tautomer.num_atoms(),
        });
    }
    if tautomer.num_bonds() != source.num_bonds() {
        return Err(TautomerStereoTransitionError::BondCountMismatch {
            expected: source.num_bonds(),
            actual: tautomer.num_bonds(),
        });
    }
    let mut modified = false;

    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     auto atomIdx = atom->getIdx();
    // RDKit✔️❌:     if (!res.d_modifiedAtoms.test(atomIdx)) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    for atom_id in source.atoms().iter().map(crate::Atom::id) {
        if !modified_atoms.contains(&atom_id) {
            continue;
        }
        let source_atom = &source.atoms()[atom_id.index()];
        let clear_isotopic_hydrogens = {
            let tautomer_atom = &tautomer.atoms()[atom_id.index()];
            !tautomer_atom.tracked_isotopic_hydrogens().is_empty()
                && (options.remove_isotopic_hydrogens()
                    || crate::valence::total_num_hydrogens(tautomer, atom_id, false)? == 0)
        };
        let tautomer_atom = &mut tautomer.topology_block_mut().atoms[atom_id.index()];

        // RDKit✔️❌:     auto tautAtom = taut.getAtomWithIdx(atomIdx);
        // RDKit✔️❌:     // clear chiral tag on sp2 atoms (also sp3 if d_removeSp3Stereo is true)
        // RDKit✔️❌:     if (tautAtom->getHybridization() == Atom::SP2 || d_removeSp3Stereo) {
        if tautomer_atom.hybridization() == Hybridization::Sp2 || options.remove_sp3_stereo() {
            // RDKit✔️❌:       modified |= (tautAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
            // RDKit✔️❌:       tautAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
            modified |= tautomer_atom.chiral_tag() != ChiralTag::Unspecified;
            tautomer_atom.set_chiral_tag(ChiralTag::Unspecified);
            // RDKit✔️❌:       if (tautAtom->hasProp(common_properties::_CIPCode)) {
            // RDKit✔️❌:         tautAtom->clearProp(common_properties::_CIPCode);
            // RDKit✔️❌:       }
            tautomer_atom.clear_prop("_CIPCode");
        } else {
            // RDKit✔️❌:     } else {
            // RDKit✔️❌:       modified |= (tautAtom->getChiralTag() != atom->getChiralTag());
            // RDKit✔️❌:       tautAtom->setChiralTag(atom->getChiralTag());
            modified |= tautomer_atom.chiral_tag() != source_atom.chiral_tag();
            tautomer_atom.set_chiral_tag(source_atom.chiral_tag());
            // RDKit✔️❌:       if (atom->hasProp(common_properties::_CIPCode)) {
            // RDKit✔️❌:         tautAtom->setProp(
            // RDKit✔️❌:             common_properties::_CIPCode,
            // RDKit✔️❌:             atom->getProp<std::string>(common_properties::_CIPCode));
            // RDKit✔️❌:       }
            if let Some(cip_code) = source_atom.prop("_CIPCode") {
                tautomer_atom.set_prop("_CIPCode", cip_code);
            }
        }
        // RDKit✔️❌:     // remove isotopic Hs if present (and if d_removeIsotopicHs is true)
        // RDKit✔️❌:     if (tautAtom->hasProp(common_properties::_isotopicHs) &&
        // RDKit✔️❌:         (d_removeIsotopicHs || !tautAtom->getTotalNumHs())) {
        // RDKit✔️❌:       tautAtom->clearProp(common_properties::_isotopicHs);
        // RDKit✔️❌:     }
        if clear_isotopic_hydrogens {
            tautomer_atom.set_tracked_isotopic_hydrogens(Vec::new());
        }
        // RDKit✔️❌:   }
    }

    // RDKit✔️❌:   // remove stereochemistry on bonds that are part of a tautomeric path
    // RDKit✔️❌:   for (auto bond : mol.bonds()) {
    // RDKit✔️❌:     auto bondIdx = bond->getIdx();
    // RDKit✔️❌:     if (!res.d_modifiedBonds.test(bondIdx)) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    for bond_id in source.bonds().iter().map(crate::Bond::id) {
        if !modified_bonds.contains(&bond_id) {
            continue;
        }
        let source_bond = &source.bonds()[bond_id.index()];
        // RDKit✔️❌:     std::vector<unsigned int> bondsToClearDirs;
        let mut bonds_to_clear_directions = Vec::new();
        // RDKit✔️❌:     if (bond->getBondType() == Bond::DOUBLE &&
        // RDKit✔️❌:         bond->getStereo() > Bond::STEREOANY) {
        if source_bond.order() == BondOrder::Double && is_stereo_beyond_any(source_bond.stereo()) {
            // RDKit✔️❌:       for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
            // RDKit✔️❌:         for (const auto &nbri :
            // RDKit✔️❌:              boost::make_iterator_range(mol.getAtomBonds(atom))) {
            for atom_id in [source_bond.begin(), source_bond.end()] {
                for neighbor in source
                    .topology_block()
                    .adjacency
                    .neighbors_of(atom_id.index())
                {
                    let adjacent_bond = &source.bonds()[neighbor.bond.index()];
                    // RDKit✔️❌:           const auto &obnd = mol[nbri];
                    // RDKit✔️❌:           if (obnd->getBondDir() == Bond::ENDDOWNRIGHT ||
                    // RDKit✔️❌:               obnd->getBondDir() == Bond::ENDUPRIGHT) {
                    // RDKit✔️❌:             bondsToClearDirs.push_back(obnd->getIdx());
                    // RDKit✔️❌:           }
                    if matches!(
                        adjacent_bond.direction(),
                        BondDirection::EndDownRight | BondDirection::EndUpRight
                    ) {
                        bonds_to_clear_directions.push(adjacent_bond.id());
                    }
                    // RDKit✔️❌:         }
                }
                // RDKit✔️❌:       }
            }
        }

        let tautomer_bond = &tautomer.bonds()[bond_id.index()];
        let remove_stereo =
            tautomer_bond.order() != BondOrder::Double || options.remove_bond_stereo();
        let target_stereo = if remove_stereo {
            let is_ring_bond = tautomer_bond_is_ring(tautomer, bond_id)?;
            if tautomer_bond.order() == BondOrder::Double && !is_ring_bond {
                BondStereo::Any
            } else {
                BondStereo::None
            }
        } else {
            source_bond.stereo()
        };
        let target_stereo_atoms = if remove_stereo {
            None
        } else {
            source_bond.stereo_atoms()
        };
        modified |= tautomer_bond.stereo() != target_stereo
            || (!remove_stereo
                && tautomer_bond.stereo_atoms().is_some() != source_bond.stereo_atoms().is_some());

        // RDKit✔️❌:     auto tautBond = taut.getBondWithIdx(bondIdx);
        // RDKit✔️❌:     if (tautBond->getBondType() != Bond::DOUBLE || d_removeBondStereo) {
        // RDKit✔️❌:       ...
        // RDKit✔️❌:       tautBond->setStereo(targetStereo);
        // RDKit✔️❌:       tautBond->getStereoAtoms().clear();
        // RDKit✔️❌:     } else {
        // RDKit✔️❌:       const INT_VECT &sa = bond->getStereoAtoms();
        // RDKit✔️❌:       if (sa.size() == 2) {
        // RDKit✔️❌:         tautBond->setStereoAtoms(sa.front(), sa.back());
        // RDKit✔️❌:       }
        // RDKit✔️❌:       tautBond->setStereo(bond->getStereo());
        // RDKit✔️❌:     }
        let tautomer_bond = &mut tautomer.topology_block_mut().bonds[bond_id.index()];
        tautomer_bond.set_stereo_atoms(target_stereo_atoms);
        tautomer_bond.set_stereo(target_stereo);
        for adjacent_bond_id in bonds_to_clear_directions {
            // RDKit✔️❌:       for (auto bi : bondsToClearDirs) {
            // RDKit✔️❌:         taut.getBondWithIdx(bi)->setBondDir(
            // RDKit✔️❌:             mol.getBondWithIdx(bi)->getBondDir());
            // RDKit✔️❌:       }
            let direction = if remove_stereo {
                BondDirection::None
            } else {
                source.bonds()[adjacent_bond_id.index()].direction()
            };
            tautomer.topology_block_mut().bonds[adjacent_bond_id.index()].set_direction(direction);
        }
        // RDKit✔️❌:   }
    }

    // RDKit✔️❌:   if (d_reassignStereo) {
    if options.reassign_stereo() {
        // RDKit✔️❌:     static const bool cleanIt = true;
        // RDKit✔️❌:     static const bool force = true;
        // RDKit✔️❌:     MolOps::assignStereochemistry(taut, cleanIt, force);
        crate::smiles::assign_stereochemistry_cleanup_subset(tautomer, true)?;

        // RDKit✔️❌:     if (d_removeBondStereo) {
        if options.remove_bond_stereo() {
            for bond_id in modified_bonds.iter().copied() {
                let bond = &tautomer.bonds()[bond_id.index()];
                let target_stereo = if bond.order() != BondOrder::Double {
                    BondStereo::None
                } else if tautomer_bond_is_ring(tautomer, bond_id)? {
                    BondStereo::None
                } else {
                    BondStereo::Any
                };
                // RDKit✔️❌:       for (auto bond : taut.bonds()) {
                // RDKit✔️❌:         const auto bondIdx = bond->getIdx();
                // RDKit✔️❌:         if (!res.d_modifiedBonds.test(bondIdx)) {
                // RDKit✔️❌:           continue;
                // RDKit✔️❌:         }
                // RDKit✔️❌:         if (bond->getBondType() != Bond::DOUBLE) {
                // RDKit✔️❌:           bond->setStereo(Bond::STEREONONE);
                // RDKit✔️❌:           bond->getStereoAtoms().clear();
                // RDKit✔️❌:           continue;
                // RDKit✔️❌:         }
                // RDKit✔️❌:         bond->setStereo(isRingBond ? Bond::STEREONONE : Bond::STEREOANY);
                // RDKit✔️❌:         bond->getStereoAtoms().clear();
                let bond = &mut tautomer.topology_block_mut().bonds[bond_id.index()];
                bond.set_stereo_atoms(None);
                bond.set_stereo(target_stereo);
                // RDKit✔️❌:       }
            }
            // RDKit✔️❌:     }
        }
    } else {
        // RDKit✔️❌:   } else {
        // RDKit✔️❌:     taut.setProp(common_properties::_StereochemDone, 1);
        tautomer.properties_mut().set_prop("_StereochemDone", "1");
        // RDKit✔️❌:   }
    }
    // RDKit✔️❌:   return modified;
    // RDKit✔️❌: }
    Ok(modified)
}

pub(crate) fn plan_tautomer_stereo_update(
    source_read: MoleculeReadParts<'_>,
    tautomer_read: MoleculeReadParts<'_>,
    modified_atoms: &BTreeSet<AtomId>,
    modified_bonds: &BTreeSet<BondId>,
    options: TautomerOptions,
) -> Result<TautomerStereoUpdatePlan, TautomerRunError> {
    let source = molecule_snapshot_from_read_parts(source_read)?;
    let mut tautomer = molecule_snapshot_from_read_parts(tautomer_read)?;
    let changed = set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut tautomer,
        modified_atoms,
        modified_bonds,
        options,
    )
    .map_err(TautomerRunError::from)?;
    Ok(TautomerStereoUpdatePlan {
        topology: tautomer.topology_block().clone(),
        properties: tautomer.properties().clone(),
        derived_cache: tautomer.derived_cache().clone(),
        changed,
    })
}

/// Ordered molecules and metadata produced by one tautomer-enumeration run.
#[derive(Debug, PartialEq)]
pub struct TautomerEnumeration {
    entries: Vec<(String, Molecule)>,
    status: TautomerEnumerationStatus,
    modified_atoms: BTreeSet<AtomId>,
    modified_bonds: BTreeSet<BondId>,
}

impl Default for TautomerEnumeration {
    fn default() -> Self {
        // RDKit✔️✔️: TautomerEnumeratorResult() : d_status(TautomerEnumeratorStatus::Completed) {}
        Self {
            entries: Vec::new(),
            status: TautomerEnumerationStatus::Completed,
            modified_atoms: BTreeSet::new(),
            modified_bonds: BTreeSet::new(),
        }
    }
}

impl Clone for TautomerEnumeration {
    fn clone(&self) -> Self {
        // RDKit✔️🔝: TautomerEnumeratorResult(const TautomerEnumeratorResult &other)
        // RDKit✔️🔝:     : d_tautomers(other.d_tautomers),
        // RDKit✔️🔝:       d_status(other.d_status),
        // RDKit✔️🔝:       d_modifiedAtoms(other.d_modifiedAtoms),
        // RDKit✔️🔝:       d_modifiedBonds(other.d_modifiedBonds) {
        // RDKit✔️🔝:   fillTautomersItVec();
        // RDKit✔️🔝: }
        // The sorted vector is both the owning representation and the random-
        // access index, so cloning needs no second iterator-vector rebuild.
        Self {
            entries: self.entries.clone(),
            status: self.status,
            modified_atoms: self.modified_atoms.clone(),
            modified_bonds: self.modified_bonds.clone(),
        }
    }
}

impl TautomerEnumeration {
    pub(crate) fn from_ordered_entries(
        entries: Vec<(String, Molecule)>,
        status: TautomerEnumerationStatus,
        modified_atoms: BTreeSet<AtomId>,
        modified_bonds: BTreeSet<BondId>,
    ) -> Self {
        Self {
            entries,
            status,
            modified_atoms,
            modified_bonds,
        }
    }

    fn from_candidates(
        candidates: SmilesTautomerMap,
        status: TautomerEnumerationStatus,
        modified_atoms: BTreeSet<AtomId>,
        modified_bonds: BTreeSet<BondId>,
    ) -> Result<Self, TautomerEnumerationError> {
        let entries = materialize_tautomer_candidates_in_source_order(candidates)?;
        Ok(Self {
            entries,
            status,
            modified_atoms,
            modified_bonds,
        })
    }

    #[must_use]
    pub fn len(&self) -> usize {
        // RDKit✔️✔️: size_t size() const { return d_tautomers.size(); }
        self.entries.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        // RDKit✔️✔️: bool empty() const { return d_tautomers.empty(); }
        self.entries.is_empty()
    }

    pub fn at(&self, index: usize) -> Result<&Molecule, TautomerEnumerationError> {
        // RDKit✔️✔️: const ROMOL_SPTR &at(size_t pos) const {
        // RDKit✔️✔️:   PRECONDITION(pos < d_tautomers.size(), "index out of bounds");
        // RDKit✔️✔️:   return d_tautomersItVec.at(pos)->second.tautomer;
        // RDKit✔️✔️: }
        self.entries.get(index).map(|(_, molecule)| molecule).ok_or(
            TautomerEnumerationError::IndexOutOfRange {
                index,
                len: self.entries.len(),
            },
        )
    }

    #[must_use]
    pub fn get(&self, index: usize) -> Option<&Molecule> {
        self.entries.get(index).map(|(_, molecule)| molecule)
    }

    pub fn iter(&self) -> impl DoubleEndedIterator<Item = &Molecule> + ExactSizeIterator {
        // RDKit✔️✔️: const const_iterator begin() const {
        // RDKit✔️✔️:   return const_iterator(d_tautomers.begin());
        // RDKit✔️✔️: }
        // RDKit✔️✔️: const const_iterator end() const { return const_iterator(d_tautomers.end()); }
        self.entries.iter().map(|(_, molecule)| molecule)
    }

    pub fn iter_with_smiles(
        &self,
    ) -> impl DoubleEndedIterator<Item = (&str, &Molecule)> + ExactSizeIterator {
        self.entries
            .iter()
            .map(|(smiles, molecule)| (smiles.as_str(), molecule))
    }

    #[must_use]
    pub fn molecules(&self) -> Vec<Molecule> {
        // RDKit✔️🔝: std::vector<ROMOL_SPTR> tautomers() const {
        // RDKit✔️🔝:   std::vector<ROMOL_SPTR> tautomerVec;
        // RDKit✔️🔝:   tautomerVec.reserve(d_tautomers.size());
        // RDKit✔️🔝:   std::transform(
        // RDKit✔️🔝:       d_tautomers.begin(), d_tautomers.end(), std::back_inserter(tautomerVec),
        // RDKit✔️🔝:       [](const SmilesTautomerPair &t) { return t.second.tautomer; });
        // RDKit✔️🔝:   return tautomerVec;
        // RDKit✔️🔝: }
        // Molecule cloning is COW-equivalent to copying RDKit shared pointers.
        self.iter().cloned().collect()
    }

    #[must_use]
    pub fn canonical_smiles(&self) -> Vec<String> {
        // RDKit✔️✔️: std::vector<std::string> smiles() const {
        // RDKit✔️✔️:   std::vector<std::string> smilesVec;
        // RDKit✔️✔️:   smilesVec.reserve(d_tautomers.size());
        // RDKit✔️✔️:   std::transform(d_tautomers.begin(), d_tautomers.end(),
        // RDKit✔️✔️:                  std::back_inserter(smilesVec),
        // RDKit✔️✔️:                  [](const SmilesTautomerPair &t) { return t.first; });
        // RDKit✔️✔️:   return smilesVec;
        // RDKit✔️✔️: }
        self.entries
            .iter()
            .map(|(smiles, _)| smiles.clone())
            .collect()
    }

    #[must_use]
    pub const fn modified_atoms(&self) -> &BTreeSet<AtomId> {
        // RDKit✔️✔️: const boost::dynamic_bitset<> &modifiedAtoms() const {
        // RDKit✔️✔️:   return d_modifiedAtoms;
        // RDKit✔️✔️: }
        &self.modified_atoms
    }

    #[must_use]
    pub const fn modified_bonds(&self) -> &BTreeSet<BondId> {
        // RDKit✔️✔️: const boost::dynamic_bitset<> &modifiedBonds() const {
        // RDKit✔️✔️:   return d_modifiedBonds;
        // RDKit✔️✔️: }
        &self.modified_bonds
    }

    #[must_use]
    pub const fn status(&self) -> TautomerEnumerationStatus {
        // RDKit✔️✔️: TautomerEnumeratorStatus status() const { return d_status; }
        self.status
    }
}

impl std::ops::Index<usize> for TautomerEnumeration {
    type Output = Molecule;

    fn index(&self, index: usize) -> &Self::Output {
        // RDKit✔️✔️: const ROMOL_SPTR &operator[](size_t pos) const { return at(pos); }
        self.at(index).expect("tautomer result index out of bounds")
    }
}

impl<'a> IntoIterator for &'a TautomerEnumeration {
    type Item = &'a Molecule;
    type IntoIter = std::iter::Map<
        std::slice::Iter<'a, (String, Molecule)>,
        fn(&(String, Molecule)) -> &Molecule,
    >;

    fn into_iter(self) -> Self::IntoIter {
        fn molecule(entry: &(String, Molecule)) -> &Molecule {
            &entry.1
        }
        self.entries.iter().map(molecule)
    }
}

impl TautomerCandidate {
    fn empty() -> Self {
        // RDKit✔️✔️: Tautomer() : d_numModifiedAtoms(0), d_numModifiedBonds(0), d_done(false) {}
        Self {
            tautomer: None,
            kekulized: None,
            num_modified_atoms: 0,
            num_modified_bonds: 0,
            done: false,
        }
    }

    fn new(
        tautomer: Molecule,
        kekulized: Molecule,
        num_modified_atoms: usize,
        num_modified_bonds: usize,
    ) -> Self {
        // RDKit✔️✔️: Tautomer(ROMOL_SPTR t, ROMOL_SPTR k, size_t a = 0, size_t b = 0)
        // RDKit✔️✔️:     : tautomer(std::move(t)),
        // RDKit✔️✔️:       kekulized(std::move(k)),
        // RDKit✔️✔️:       d_numModifiedAtoms(a),
        // RDKit✔️✔️:       d_numModifiedBonds(b),
        // RDKit✔️✔️:       d_done(false) {}
        Self {
            tautomer: Some(tautomer),
            kekulized: Some(kekulized),
            num_modified_atoms,
            num_modified_bonds,
            done: false,
        }
    }
}

impl<M> TautomerCandidate<M> {
    fn mark_done(&mut self) {
        // RDKit✔️✔️: smilesTautomerPair.second.d_done = true;
        self.done = true;
    }

    fn update_modified_counts(&mut self, num_modified_atoms: usize, num_modified_bonds: usize) {
        // RDKit✔️✔️: tautStored.d_numModifiedAtoms = maxNumModifiedAtoms;
        // RDKit✔️✔️: tautStored.d_numModifiedBonds = maxNumModifiedBonds;
        self.num_modified_atoms = num_modified_atoms;
        self.num_modified_bonds = num_modified_bonds;
    }
}

impl TautomerScore {
    #[must_use]
    pub const fn ring(self) -> i32 {
        self.ring
    }

    #[must_use]
    pub const fn substructure(self) -> i32 {
        self.substructure
    }

    #[must_use]
    pub const fn hetero_hydrogen(self) -> i32 {
        self.hetero_hydrogen
    }

    #[must_use]
    pub const fn total(self) -> i32 {
        self.ring
            .wrapping_add(self.substructure)
            .wrapping_add(self.hetero_hydrogen)
    }
}

impl TautomerScoreTerm {
    #[must_use]
    pub fn new(name: impl Into<String>, smarts: impl Into<String>, score: i32) -> Self {
        // RDKit✔️✔️: SubstructTerm::SubstructTerm(std::string aname, std::string asmarts, int ascore)
        // RDKit✔️✔️:     : name(std::move(aname)), smarts(std::move(asmarts)), score(ascore) {
        // RDKit✔️✔️:   std::unique_ptr<ROMol> pattern(SmartsToMol(smarts));
        // RDKit✔️✔️:   if (pattern) {
        // RDKit✔️✔️:     matcher = std::move(*pattern);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        let name = name.into();
        let smarts = smarts.into();
        let matcher = parse_smarts(&smarts, &SmartsParseParams::default()).ok();
        Self {
            name,
            smarts,
            score,
            matcher,
        }
    }

    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub fn smarts(&self) -> &str {
        &self.smarts
    }

    #[must_use]
    pub const fn score(&self) -> i32 {
        self.score
    }

    pub(crate) const fn matcher(&self) -> Option<&QueryGraph> {
        self.matcher.as_ref()
    }
}

impl PartialEq for TautomerScoreTerm {
    fn eq(&self, other: &Self) -> bool {
        // RDKit✔️✔️: bool operator==(const SubstructTerm &rhs) const {
        // RDKit✔️✔️:   return name == rhs.name && smarts == rhs.smarts && score == rhs.score;
        // RDKit✔️✔️: }
        self.name == other.name && self.smarts == other.smarts && self.score == other.score
    }
}

impl Eq for TautomerScoreTerm {}

#[must_use]
pub fn default_tautomer_score_terms() -> &'static [TautomerScoreTerm] {
    // RDKit✔️✔️: const std::vector<SubstructTerm> &getDefaultTautomerScoreSubstructs() {
    // RDKit✔️✔️:   static std::vector<SubstructTerm> substructureTerms{
    static TERMS: OnceLock<Vec<TautomerScoreTerm>> = OnceLock::new();
    TERMS.get_or_init(|| {
        vec![
            // RDKit✔️✔️:       {"benzoquinone", "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
            // RDKit✔️✔️:        25},
            TautomerScoreTerm::new(
                "benzoquinone",
                "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
                25,
            ),
            // RDKit✔️✔️:       {"oxim", "[#6]=[N][OH]", 4},
            TautomerScoreTerm::new("oxim", "[#6]=[N][OH]", 4),
            // RDKit✔️✔️:       {"C=O", "[#6]=,:[#8]", 2},
            TautomerScoreTerm::new("C=O", "[#6]=,:[#8]", 2),
            // RDKit✔️✔️:       {"N=O", "[#7]=,:[#8]", 2},
            TautomerScoreTerm::new("N=O", "[#7]=,:[#8]", 2),
            // RDKit✔️✔️:       {"P=O", "[#15]=,:[#8]", 2},
            TautomerScoreTerm::new("P=O", "[#15]=,:[#8]", 2),
            // RDKit✔️✔️:       {"C=hetero", "[C]=[!#1;!#6]", 1},
            TautomerScoreTerm::new("C=hetero", "[C]=[!#1;!#6]", 1),
            // RDKit✔️✔️:       {"C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2},
            TautomerScoreTerm::new("C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2),
            // RDKit✔️✔️:       {"aromatic C = exocyclic N", "[c]=!@[N]", -1},
            TautomerScoreTerm::new("aromatic C = exocyclic N", "[c]=!@[N]", -1),
            // RDKit✔️✔️:       {"methyl", "[CX4H3]", 1},
            TautomerScoreTerm::new("methyl", "[CX4H3]", 1),
            // RDKit✔️✔️:       {"guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1},
            TautomerScoreTerm::new("guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1),
            // RDKit✔️✔️:       {"guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2},
            TautomerScoreTerm::new("guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2),
            // RDKit✔️✔️:       {"aci-nitro", "[#6]=[N+]([O-])[OH]", -4}};
            TautomerScoreTerm::new("aci-nitro", "[#6]=[N+]([O-])[OH]", -4),
        ]
    })
    // RDKit✔️✔️:   return substructureTerms;
    // RDKit✔️✔️: }
}

/// Score aromatic rings using RDKit's tautomer-canonicalization weights.
pub fn score_tautomer_rings(molecule: &Molecule) -> Result<i32, RingFindingError> {
    // RDKit✔️✔️: int scoreRings(const ROMol &mol) {
    // RDKit✔️✔️:   int score = 0;
    let mut score = 0_i32;

    // RDKit✔️✔️:   auto ringInfo = mol.getRingInfo();
    // RDKit✔️✔️:   std::unique_ptr<ROMol> cp;
    // RDKit✔️✔️:   if (!ringInfo->isSymmSssr()) {
    // RDKit✔️✔️:     cp.reset(new ROMol(mol));
    // RDKit✔️✔️:     MolOps::symmetrizeSSSR(*cp);
    // RDKit✔️✔️:     ringInfo = cp->getRingInfo();
    // RDKit✔️✔️:   }
    let computed_rings;
    let ring_info = match molecule.derived_cache().rings.as_ref() {
        Some(rings) if rings.is_symm_sssr() => rings,
        _ => {
            computed_rings = crate::symmetrize_sssr(molecule)?;
            &computed_rings
        }
    };

    // RDKit✔️✔️:   boost::dynamic_bitset<> isArom(mol.getNumBonds());
    // RDKit✔️✔️:   boost::dynamic_bitset<> bothCarbon(mol.getNumBonds());
    let mut is_aromatic = vec![false; molecule.num_bonds()];
    let mut both_carbon = vec![false; molecule.num_bonds()];
    // RDKit✔️✔️:   for (const auto &bnd : mol.bonds()) {
    for bond in molecule.bonds() {
        // RDKit✔️✔️:     if (bnd->getIsAromatic()) {
        if bond.is_aromatic() {
            // RDKit✔️✔️:       isArom.set(bnd->getIdx());
            is_aromatic[bond.id().index()] = true;
            // RDKit✔️✔️:       if (bnd->getBeginAtom()->getAtomicNum() == 6 &&
            // RDKit✔️✔️:           bnd->getEndAtom()->getAtomicNum() == 6) {
            if molecule.atoms()[bond.begin().index()].atomic_number() == 6
                && molecule.atoms()[bond.end().index()].atomic_number() == 6
            {
                // RDKit✔️✔️:         bothCarbon.set(bnd->getIdx());
                both_carbon[bond.id().index()] = true;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   for (const auto &bring : ringInfo->bondRings()) {
    for bond_ring in ring_info.bond_rings() {
        // RDKit✔️✔️:     bool allC = true;
        // RDKit✔️✔️:     bool allAromatic = true;
        let mut all_carbon = true;
        let mut all_aromatic = true;
        // RDKit✔️✔️:     for (const auto bidx : bring) {
        for bond_id in bond_ring {
            // RDKit✔️✔️:       if (!isArom[bidx]) {
            if !is_aromatic[bond_id.index()] {
                // RDKit✔️✔️:         allAromatic = false;
                all_aromatic = false;
                // RDKit✔️✔️:         break;
                break;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       if (!bothCarbon[bidx]) {
            if !both_carbon[bond_id.index()] {
                // RDKit✔️✔️:         allC = false;
                all_carbon = false;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (allAromatic) {
        if all_aromatic {
            // RDKit✔️✔️:       score += 100;
            score += 100;
            // RDKit✔️✔️:       if (allC) {
            if all_carbon {
                // RDKit✔️✔️:         score += 150;
                score += 150;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return score;
    Ok(score)
    // RDKit✔️✔️: };
}

/// Score every match of the supplied source-ordered tautomer SMARTS terms.
#[must_use]
pub fn score_tautomer_substructures(molecule: &Molecule, terms: &[TautomerScoreTerm]) -> i32 {
    // RDKit✔️✔️: int scoreSubstructs(const ROMol &mol,
    // RDKit✔️✔️:                     const std::vector<SubstructTerm> &substructureTerms) {
    // RDKit✔️✔️:   int score = 0;
    let mut score = 0_i32;
    // RDKit✔️✔️:   for (const auto &term : substructureTerms) {
    for term in terms {
        // RDKit✔️✔️:     if (!term.matcher.getNumAtoms()) {
        let Some(matcher) = term.matcher() else {
            // RDKit✔️✔️:       BOOST_LOG(rdErrorLog) << " matcher for term " << term.name
            // RDKit✔️✔️:                             << " is invalid, ignoring it." << std::endl;
            // RDKit✔️✔️:       continue;
            continue;
            // RDKit✔️✔️:     }
        };
        // RDKit✔️✔️:     SubstructMatchParameters params;
        // RDKit✔️✔️:     const auto matches = SubstructMatch(mol, term.matcher, params);
        let matches = crate::get_substruct_matches(molecule, matcher);
        // RDKit✔️✔️:     // if (!matches.empty()) {
        // RDKit✔️✔️:     //   std::cerr << " " << matches.size() << " matches to " << term.name
        // RDKit✔️✔️:     //             << std::endl;
        // RDKit✔️✔️:     // }
        // RDKit✔️✔️:     score += static_cast<int>(matches.size()) * term.score;
        // The pinned 32-bit source ABI narrows `size_t` to `int` here; Rust's
        // cast and wrapping arithmetic retain that release-build machine
        // behavior instead of introducing a debug-only overflow branch.
        let match_count = matches.len() as i32;
        score = score.wrapping_add(match_count.wrapping_mul(term.score));
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return score;
    score
    // RDKit✔️✔️: }
}

/// Penalize hydrogens attached to phosphorus, sulfur, selenium, and tellurium.
pub fn score_tautomer_hetero_hydrogens(molecule: &Molecule) -> Result<i32, TautomerScoreError> {
    // RDKit✔️✔️: int scoreHeteroHs(const ROMol &mol) {
    // RDKit✔️✔️:   int score = 0;
    let mut score = 0_i32;
    // RDKit✔️✔️:   for (const auto &at : mol.atoms()) {
    for atom in molecule.atoms() {
        // RDKit✔️✔️:     int anum = at->getAtomicNum();
        let atomic_number = atom.atomic_number();
        // RDKit✔️✔️:     if (anum == 15 || anum == 16 || anum == 34 || anum == 52) {
        if matches!(atomic_number, 15 | 16 | 34 | 52) {
            // RDKit✔️✔️:       score -= at->getTotalNumHs();
            let total_hydrogens = crate::valence::total_num_hydrogens(molecule, atom.id(), false)?;
            score = score.wrapping_sub(total_hydrogens as i32);
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return score;
    Ok(score)
    // RDKit✔️✔️: }
}

/// Compute all source tautomer score components and their signed aggregate.
pub fn score_tautomer_with_terms(
    molecule: &Molecule,
    terms: &[TautomerScoreTerm],
) -> Result<TautomerScore, TautomerScoreError> {
    // RDKit✔️✔️: inline int scoreTautomer(const ROMol &mol) {
    // RDKit✔️✔️:   return scoreRings(mol) + scoreSubstructs(mol) + scoreHeteroHs(mol);
    // RDKit✔️✔️: }
    Ok(TautomerScore {
        ring: score_tautomer_rings(molecule)?,
        substructure: score_tautomer_substructures(molecule, terms),
        hetero_hydrogen: score_tautomer_hetero_hydrogens(molecule)?,
    })
}

/// Compute the source tautomer score using the pinned default SMARTS terms.
pub fn score_tautomer(molecule: &Molecule) -> Result<TautomerScore, TautomerScoreError> {
    score_tautomer_with_terms(molecule, default_tautomer_score_terms())
}

impl TautomerCatalog {
    /// Construct the current pinned built-in catalog.
    pub fn current() -> Result<Self, TautomerCatalogError> {
        // RDKit✔️✔️: TautomerCatalogParams::TautomerCatalogParams(const std::string &tautomerFile) {
        // RDKit✔️✔️:   d_transforms.clear();
        // RDKit✔️✔️:   if (tautomerFile.empty()) {
        // RDKit✔️✔️:     d_transforms = readTautomers(defaults::defaultTautomerTransforms);
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     d_transforms = readTautomers(tautomerFile);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        Ok(Self {
            transforms: current_builtin_transforms()?,
        })
    }

    /// Construct the pinned pre-2021.09 catalog.
    pub fn v1() -> Result<Self, TautomerCatalogError> {
        Ok(Self {
            transforms: v1_builtin_transforms()?,
        })
    }

    /// Construct a catalog from the source four-column data representation.
    pub fn from_data(data: &[(&str, &str, &str, &str)]) -> Result<Self, TautomerCatalogError> {
        // RDKit✔️✔️: TautomerCatalogParams::TautomerCatalogParams(
        // RDKit✔️✔️:     const TautomerTransformDefs &data) {
        // RDKit✔️✔️:   d_transforms.clear();
        // RDKit✔️✔️:   d_transforms = readTautomers(data);
        // RDKit✔️✔️: }
        Ok(Self {
            transforms: read_transforms_from_definitions(data)?,
        })
    }

    /// Construct a catalog from an RDKit-format transform file.
    pub fn from_file(file_name: impl AsRef<Path>) -> Result<Self, TautomerCatalogError> {
        if file_name.as_ref().as_os_str().is_empty() {
            return Self::current();
        }
        Ok(Self {
            transforms: read_transforms_from_file(file_name)?,
        })
    }

    #[must_use]
    pub fn transforms(&self) -> &[TautomerTransform] {
        // RDKit✔️✔️: const std::vector<TautomerTransform> &TautomerCatalogParams::getTransforms()
        // RDKit✔️✔️:     const {
        // RDKit✔️✔️:   return d_transforms;
        // RDKit✔️✔️: }
        &self.transforms
    }

    pub fn transform(&self, index: usize) -> Result<TautomerTransform, TautomerCatalogError> {
        // RDKit✔️🔝: const TautomerTransform TautomerCatalogParams::getTransform(
        // RDKit✔️🔝:     unsigned int fid) const {
        // RDKit✔️🔝:   URANGE_CHECK(fid, d_transforms.size());
        // RDKit✔️🔝:   return d_transforms[fid];  //.get();
        // RDKit✔️🔝: }
        // Cloning retains RDKit's returned-value independence while COW shares
        // immutable query blocks until a caller actually changes its copy.
        self.transforms
            .get(index)
            .cloned()
            .ok_or(TautomerCatalogError::TransformIndexOutOfRange {
                index,
                len: self.transforms.len(),
            })
    }

    pub fn write_to(&self, stream: &mut impl Write) -> io::Result<()> {
        // RDKit✔️✔️: void TautomerCatalogParams::toStream(std::ostream &ss) const {
        // RDKit✔️✔️:   ss << d_transforms.size() << "\n";
        // RDKit✔️✔️: }
        writeln!(stream, "{}", self.transforms.len())
    }

    #[must_use]
    pub fn serialize(&self) -> String {
        // RDKit✔️✔️: std::string TautomerCatalogParams::Serialize() const {
        // RDKit✔️✔️:   std::stringstream ss;
        // RDKit✔️✔️:   toStream(ss);
        // RDKit✔️✔️:   return ss.str();
        // RDKit✔️✔️: }
        format!("{}\n", self.transforms.len())
    }

    pub fn deserialize(_serialized: &str) -> Result<Self, TautomerCatalogError> {
        Err(TautomerCatalogError::DeserializationUnderConstruction)
    }
}

impl Default for TautomerCatalog {
    fn default() -> Self {
        Self::current().expect("pinned built-in tautomer catalog must compile")
    }
}

impl Clone for TautomerCatalog {
    fn clone(&self) -> Self {
        // RDKit✔️🔝: TautomerCatalogParams::TautomerCatalogParams(
        // RDKit✔️🔝:     const TautomerCatalogParams &other) {
        // RDKit✔️🔝:   d_typeStr = other.d_typeStr;
        // RDKit✔️🔝:   d_transforms.clear();
        // RDKit✔️🔝:
        // RDKit✔️🔝:   const std::vector<TautomerTransform> &transforms = other.getTransforms();
        // RDKit✔️🔝:   for (const auto &transform : transforms) {
        // RDKit✔️🔝:     d_transforms.push_back(transform);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // One `Vec::clone` preallocates exactly once; each transform retains
        // source value semantics and COW-shares its immutable compiled query.
        Self {
            transforms: self.transforms.clone(),
        }
    }

    fn clone_from(&mut self, source: &Self) {
        self.transforms.clone_from(&source.transforms);
    }
}

impl TautomerTransform {
    /// Construct a transform from an already compiled SMARTS query.
    pub fn new(
        name: impl Into<String>,
        query: QueryGraph,
        bond_types: Vec<BondOrder>,
        charges: Vec<i32>,
    ) -> Result<Self, TautomerTransformError> {
        // RDKit✔️✔️: TautomerTransform(ROMol *mol, std::vector<Bond::BondType> bondtypes,
        // RDKit✔️✔️:                   std::vector<int> charges)
        // RDKit✔️✔️:     : Mol(mol),
        // RDKit✔️✔️:       BondTypes(std::move(bondtypes)),
        // RDKit✔️✔️:       Charges(std::move(charges)) {}
        //
        // Rust owns the compiled query directly and validates the sizes that
        // RDKit later indexes without checking. This changes only invalid
        // construction from undefined/out-of-range behavior into a structured
        // error; every source-reachable transform has the validated shape.
        if query.num_atoms() < 2 {
            return Err(TautomerTransformError::MissingDonorOrAcceptor {
                actual: query.num_atoms(),
            });
        }
        if !bond_types.is_empty() && bond_types.len() != query.num_bonds() {
            return Err(TautomerTransformError::BondEditCount {
                expected: query.num_bonds(),
                actual: bond_types.len(),
            });
        }
        if !charges.is_empty() && charges.len() != query.num_atoms() {
            return Err(TautomerTransformError::ChargeEditCount {
                expected: query.num_atoms(),
                actual: charges.len(),
            });
        }

        Ok(Self {
            query: query.with_name(name),
            bond_types,
            charges,
        })
    }

    #[must_use]
    pub fn name(&self) -> &str {
        self.query.prop("_Name").unwrap_or("")
    }

    #[must_use]
    pub const fn query(&self) -> &QueryGraph {
        &self.query
    }

    #[must_use]
    pub fn bond_types(&self) -> &[BondOrder] {
        &self.bond_types
    }

    #[must_use]
    pub fn charges(&self) -> &[i32] {
        &self.charges
    }
}

impl Clone for TautomerTransform {
    fn clone(&self) -> Self {
        // RDKit✔️🔝: TautomerTransform(const TautomerTransform &other)
        // RDKit✔️🔝:     : BondTypes(other.BondTypes), Charges(other.Charges) {
        // RDKit✔️🔝:   Mol = new ROMol(*other.Mol);
        // RDKit✔️🔝: }
        //
        // `Molecule::clone()` uses immutable Arc-backed blocks. It preserves
        // the source's independent value semantics while avoiding RDKit's
        // eager deep copy; any later value-style change detaches the block.
        Self {
            query: self.query.clone(),
            bond_types: self.bond_types.clone(),
            charges: self.charges.clone(),
        }
    }

    fn clone_from(&mut self, source: &Self) {
        // RDKit✔️🔝: TautomerTransform &operator=(const TautomerTransform &other) {
        // RDKit✔️🔝:   if (this != &other) {
        // RDKit✔️🔝:     delete Mol;
        // RDKit✔️🔝:     Mol = new ROMol(*other.Mol);
        // RDKit✔️🔝:     BondTypes = other.BondTypes;
        // RDKit✔️🔝:     Charges = other.Charges;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:   return *this;
        // RDKit✔️🔝: }
        //
        // Rust borrowing prevents an ordinary self-assignment alias. Reusing
        // each allocation through `clone_from` is at least as efficient as the
        // source delete-and-deep-copy path and leaves identical value state.
        self.query.clone_from(&source.query);
        self.bond_types.clone_from(&source.bond_types);
        self.charges.clone_from(&source.charges);
    }
}

#[cfg(test)]
mod tests {
    use std::io::{self, BufRead, Cursor, Read};

    use super::*;
    use crate::{SmartsParseParams, parse_smarts};

    fn query(smarts: &str) -> QueryGraph {
        parse_smarts(smarts, &SmartsParseParams::default()).expect("compile transform SMARTS")
    }

    #[test]
    fn catalog_tokens_map_every_source_bond_symbol_in_order() {
        assert_eq!(
            string_to_bond_types("-=#:"),
            vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Triple,
                BondOrder::Aromatic,
            ]
        );
    }

    #[test]
    fn catalog_tokens_map_every_source_charge_symbol_in_order() {
        assert_eq!(
            string_to_charges("+0-").expect("valid charges"),
            &[1, 0, -1]
        );
    }

    #[test]
    fn catalog_tokens_accept_empty_bond_and_charge_fields() {
        assert!(string_to_bond_types("").is_empty());
        assert!(string_to_charges("").expect("empty charges").is_empty());
    }

    #[test]
    fn catalog_tokens_bond_parser_ignores_whitespace_and_unknown_bytes() {
        assert_eq!(
            string_to_bond_types(" - x\t=\n?#:µ"),
            vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Triple,
                BondOrder::Aromatic,
            ]
        );
    }

    #[test]
    fn catalog_tokens_charge_parser_rejects_whitespace() {
        assert_eq!(
            string_to_charges("+ -"),
            Err(TautomerTransformError::ChargeSymbolNotRecognised)
        );
    }

    #[test]
    fn catalog_tokens_charge_parser_rejects_invalid_ascii_and_non_ascii_bytes() {
        for invalid in ["1", "9", "x", "µ"] {
            assert_eq!(
                string_to_charges(invalid),
                Err(TautomerTransformError::ChargeSymbolNotRecognised),
                "invalid token {invalid:?}"
            );
        }
    }

    #[test]
    fn catalog_tokens_charge_parser_does_not_accept_signed_integer_spellings() {
        assert_eq!(
            string_to_charges("+1"),
            Err(TautomerTransformError::ChargeSymbolNotRecognised)
        );
        assert_eq!(
            string_to_charges("-1"),
            Err(TautomerTransformError::ChargeSymbolNotRecognised)
        );
    }

    #[test]
    fn catalog_tokens_charge_parser_returns_at_first_invalid_byte() {
        assert_eq!(
            string_to_charges("+0x-"),
            Err(TautomerTransformError::ChargeSymbolNotRecognised)
        );
    }

    #[test]
    fn transform_lines_skip_empty_comment_and_too_short_lines() {
        for line in ["", "// comment", " ", "only-a-name"] {
            assert!(
                transform_from_line(line)
                    .expect("skipped line is not an error")
                    .is_none(),
                "line {line:?}"
            );
        }
    }

    #[test]
    fn transform_lines_parse_two_columns_with_source_defaults() {
        let transform = transform_from_line("Simple name\t[O]-[C]")
            .expect("valid line")
            .expect("transform");

        assert_eq!(transform.name(), "Simplename");
        assert_eq!(transform.query().num_atoms(), 2);
        assert!(transform.bond_types().is_empty());
        assert!(transform.charges().is_empty());
    }

    #[test]
    fn transform_lines_parse_three_columns_and_remove_ascii_spaces() {
        let transform = transform_from_line("Bond edit\t[O] - [C] = [C]\t= -")
            .expect("valid line")
            .expect("transform");

        assert_eq!(transform.name(), "Bondedit");
        assert_eq!(
            transform.bond_types(),
            &[BondOrder::Double, BondOrder::Single]
        );
        assert!(transform.charges().is_empty());
    }

    #[test]
    fn transform_lines_parse_four_columns_in_source_order() {
        let transform = transform_from_line("Charged\t[N]-[C]=[O]\t= -\t+ 0 -")
            .expect("valid line")
            .expect("transform");

        assert_eq!(
            transform.bond_types(),
            &[BondOrder::Double, BondOrder::Single]
        );
        assert_eq!(transform.charges(), &[1, 0, -1]);
    }

    #[test]
    fn transform_lines_collapse_omitted_optional_tab_columns_like_boost_tokenizer() {
        let transform = transform_from_line("Collapsed\t[O]-[C]\t\t+")
            .expect("source treats the fourth token as the third nonempty token")
            .expect("transform");

        assert!(transform.bond_types().is_empty());
        assert!(transform.charges().is_empty());
    }

    #[test]
    fn transform_lines_do_not_treat_indented_slashes_as_a_comment() {
        let error = transform_from_line(" // name\tinvalid smarts")
            .expect_err("only columns beginning with two slashes are comments");

        assert!(matches!(
            error,
            TautomerTransformError::CannotParseSmarts { .. }
        ));
    }

    #[test]
    fn transform_lines_report_malformed_smarts_with_the_source_text() {
        let error =
            transform_from_line("broken\t[C").expect_err("malformed SMARTS must remain visible");

        assert!(matches!(
            error,
            TautomerTransformError::CannotParseSmarts { ref smarts, .. } if smarts == "[C"
        ));
    }

    #[test]
    fn transform_lines_report_invalid_charge_symbols_before_smarts_parsing() {
        let error = transform_from_fields("broken", "[", "", "x")
            .expect_err("source converts charge tokens before parsing SMARTS");

        assert_eq!(error, TautomerTransformError::ChargeSymbolNotRecognised);
    }

    #[test]
    fn transform_lines_report_more_than_four_nonempty_columns_as_invalid_shape() {
        let error = transform_from_line("one\ttwo\tthree\tfour\tfive")
            .expect_err("source leaves all fields defaulted for this invalid line");

        assert_eq!(
            error,
            TautomerTransformError::MissingDonorOrAcceptor { actual: 0 }
        );
    }

    #[test]
    fn catalog_reader_reads_stream_comments_blank_lines_and_early_eof() {
        let data = "// header\n\nfirst\t[O]-[C]\n // not a comment\nsecond\t[N]-[C]";
        let mut reader = Cursor::new(data.as_bytes());
        let transforms = read_transforms(&mut reader, -1).expect("read stream");

        assert_eq!(
            transforms
                .iter()
                .map(TautomerTransform::name)
                .collect::<Vec<_>>(),
            ["first", "second"]
        );
    }

    #[test]
    fn catalog_reader_bounded_stream_counts_only_emitted_transforms() {
        let data = "// header\nfirst\t[O]-[C]\n\n// middle\nsecond\t[N]-[C]\nthird\t[S]-[C]\n";
        let mut reader = Cursor::new(data.as_bytes());
        let transforms = read_transforms(&mut reader, 2).expect("read bounded stream");

        assert_eq!(
            transforms
                .iter()
                .map(TautomerTransform::name)
                .collect::<Vec<_>>(),
            ["first", "second"]
        );
    }

    #[test]
    fn catalog_reader_zero_limit_reads_nothing() {
        let mut reader = Cursor::new(b"first\t[O]-[C]\n".as_slice());
        assert!(
            read_transforms(&mut reader, 0)
                .expect("zero limit")
                .is_empty()
        );
        assert_eq!(reader.position(), 0);
    }

    #[test]
    fn catalog_reader_file_entry_preserves_source_order() {
        let directory = tempfile::tempdir().expect("temporary directory");
        let path = directory.path().join("tautomers.in");
        std::fs::write(&path, "first\t[O]-[C]\nsecond\t[N]-[C]\nthird\t[S]-[C]\n")
            .expect("write fixture");

        let transforms = read_transforms_from_file(&path).expect("read file");
        assert_eq!(
            transforms
                .iter()
                .map(TautomerTransform::name)
                .collect::<Vec<_>>(),
            ["first", "second", "third"]
        );
    }

    #[test]
    fn catalog_reader_embedded_definitions_use_the_same_constructor() {
        let definitions = [
            ("first", "[O]-[C]", "=", "+0"),
            ("second", "[N]-[C]", "-", "0-"),
        ];
        let transforms =
            read_transforms_from_definitions(&definitions).expect("read embedded definitions");

        assert_eq!(transforms[0].name(), "first");
        assert_eq!(transforms[0].bond_types(), &[BondOrder::Double]);
        assert_eq!(transforms[0].charges(), &[1, 0]);
        assert_eq!(transforms[1].name(), "second");
        assert_eq!(transforms[1].bond_types(), &[BondOrder::Single]);
        assert_eq!(transforms[1].charges(), &[0, -1]);
    }

    #[test]
    fn catalog_reader_missing_file_returns_structured_io_error() {
        let directory = tempfile::tempdir().expect("temporary directory");
        let path = directory.path().join("missing.in");
        let error = read_transforms_from_file(&path).expect_err("missing file must fail");

        assert!(matches!(
            error,
            TautomerCatalogError::BadInputFile { path: actual, .. } if actual == path
        ));
    }

    struct BrokenReader;

    impl Read for BrokenReader {
        fn read(&mut self, _buffer: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("broken catalog stream"))
        }
    }

    impl BufRead for BrokenReader {
        fn fill_buf(&mut self) -> io::Result<&[u8]> {
            Err(io::Error::other("broken catalog stream"))
        }

        fn consume(&mut self, _amount: usize) {}
    }

    #[test]
    fn catalog_reader_stream_failure_returns_structured_io_error() {
        let error = read_transforms(&mut BrokenReader, -1).expect_err("bad stream must fail");
        assert!(matches!(error, TautomerCatalogError::BadStreamContents(_)));
    }

    #[test]
    fn catalog_reader_rejects_non_utf8_catalog_data() {
        let mut reader = Cursor::new([0xff, b'\n']);
        let error = read_transforms(&mut reader, -1).expect_err("catalog text must be UTF-8");
        assert!(matches!(error, TautomerCatalogError::InvalidUtf8(_)));
    }

    #[test]
    fn catalog_reader_source_line_limit_sets_fail_state_after_truncated_line() {
        let mut bytes = vec![b'/'; 511];
        bytes.extend_from_slice(b"\nvalid\t[O]-[C]\n");
        let mut reader = Cursor::new(bytes);

        let transforms = read_transforms(&mut reader, -1).expect("read truncated source line");
        assert!(transforms.is_empty());
    }

    fn assert_builtin_catalog(
        actual: &[TautomerTransform],
        expected: &[TautomerTransformDefinition<'_>],
    ) {
        assert_eq!(actual.len(), expected.len());
        for (index, (transform, &(name, smarts, bonds, charges))) in
            actual.iter().zip(expected).enumerate()
        {
            assert_eq!(transform.name(), name, "name at source row {index}");
            let expected_query = query(smarts);
            assert_eq!(
                transform.query().atoms(),
                expected_query.atoms(),
                "compiled query atoms at source row {index}: {name}"
            );
            assert_eq!(
                transform.query().bonds(),
                expected_query.bonds(),
                "compiled query bonds at source row {index}: {name}"
            );
            let expected_bonds: Vec<_> = bonds
                .bytes()
                .filter_map(|byte| match byte {
                    b'-' => Some(BondOrder::Single),
                    b'=' => Some(BondOrder::Double),
                    b'#' => Some(BondOrder::Triple),
                    b':' => Some(BondOrder::Aromatic),
                    _ => None,
                })
                .collect();
            assert_eq!(
                transform.bond_types(),
                expected_bonds,
                "bond edits at source row {index}: {name}"
            );
            let expected_charges: Vec<_> = charges
                .bytes()
                .map(|byte| match byte {
                    b'+' => 1,
                    b'0' => 0,
                    b'-' => -1,
                    _ => panic!("invalid generated charge byte at source row {index}"),
                })
                .collect();
            assert_eq!(
                transform.charges(),
                expected_charges,
                "charge edits at source row {index}: {name}"
            );
        }
    }

    #[test]
    fn builtin_catalogs_current_contains_every_compiled_source_row_in_order() {
        assert_eq!(CURRENT_TAUTOMER_TRANSFORM_DEFINITIONS.len(), 37);
        let transforms = current_builtin_transforms().expect("compile current catalog");
        assert_builtin_catalog(&transforms, CURRENT_TAUTOMER_TRANSFORM_DEFINITIONS);
    }

    #[test]
    fn builtin_catalogs_v1_contains_every_compiled_source_row_in_order() {
        assert_eq!(V1_TAUTOMER_TRANSFORM_DEFINITIONS.len(), 36);
        let transforms = v1_builtin_transforms().expect("compile V1 catalog");
        assert_builtin_catalog(&transforms, V1_TAUTOMER_TRANSFORM_DEFINITIONS);
    }

    #[test]
    fn catalog_object_default_and_empty_file_name_select_current_catalog() {
        let default_catalog = TautomerCatalog::default();
        let empty_file_catalog =
            TautomerCatalog::from_file("").expect("empty file selects default");

        assert_eq!(default_catalog.transforms().len(), 37);
        assert_eq!(empty_file_catalog, default_catalog);
    }

    #[test]
    fn catalog_object_constructs_from_file_data_and_v1_sources() {
        let directory = tempfile::tempdir().expect("temporary directory");
        let path = directory.path().join("custom.in");
        std::fs::write(&path, "file-first\t[O]-[C]\nfile-second\t[N]-[C]\n")
            .expect("write catalog");
        let file_catalog = TautomerCatalog::from_file(&path).expect("file catalog");
        let data_catalog =
            TautomerCatalog::from_data(&[("data", "[S]-[C]", "=", "+0")]).expect("data catalog");
        let v1_catalog = TautomerCatalog::v1().expect("V1 catalog");

        assert_eq!(
            file_catalog
                .transforms()
                .iter()
                .map(TautomerTransform::name)
                .collect::<Vec<_>>(),
            ["file-first", "file-second"]
        );
        assert_eq!(data_catalog.transforms()[0].name(), "data");
        assert_eq!(
            data_catalog.transforms()[0].bond_types(),
            &[BondOrder::Double]
        );
        assert_eq!(data_catalog.transforms()[0].charges(), &[1, 0]);
        assert_eq!(v1_catalog.transforms().len(), 36);
    }

    #[test]
    fn catalog_object_indexing_returns_an_independent_value_and_checks_bounds() {
        let catalog = TautomerCatalog::from_data(&[("one", "[O]-[C]", "-", "")]).expect("catalog");
        let mut transform = catalog.transform(0).expect("first transform");
        transform
            .query
            .atom_mut(0)
            .expect("query atom")
            .atom_mut()
            .set_formal_charge(-1);
        transform.bond_types[0] = BondOrder::Double;

        assert_eq!(
            catalog.transforms()[0].query().atoms()[0]
                .atom()
                .formal_charge(),
            0
        );
        assert_eq!(catalog.transforms()[0].bond_types(), &[BondOrder::Single]);
        assert!(matches!(
            catalog.transform(1),
            Err(TautomerCatalogError::TransformIndexOutOfRange { index: 1, len: 1 })
        ));
    }

    #[test]
    fn catalog_object_clone_and_clone_from_have_independent_value_semantics() {
        let source =
            TautomerCatalog::from_data(&[("source", "[O]-[C]", "-", "")]).expect("source catalog");
        let mut cloned = source.clone();
        cloned.transforms[0].bond_types[0] = BondOrder::Double;
        assert_eq!(source.transforms()[0].bond_types(), &[BondOrder::Single]);

        let mut target =
            TautomerCatalog::from_data(&[("target", "[N]-[C]", "=", "")]).expect("target catalog");
        target.clone_from(&source);
        assert_eq!(target, source);
    }

    #[test]
    fn catalog_object_serialization_is_exactly_count_and_newline() {
        let catalog = TautomerCatalog::from_data(&[
            ("first", "[O]-[C]", "", ""),
            ("second", "[N]-[C]", "", ""),
        ])
        .expect("catalog");
        let mut bytes = Vec::new();

        catalog.write_to(&mut bytes).expect("write catalog count");
        assert_eq!(bytes, b"2\n");
        assert_eq!(catalog.serialize(), "2\n");
    }

    #[test]
    fn catalog_object_deserialization_is_a_structured_upstream_boundary() {
        assert!(matches!(
            TautomerCatalog::deserialize("37\n"),
            Err(TautomerCatalogError::DeserializationUnderConstruction)
        ));
    }

    #[test]
    fn score_terms_match_all_twelve_source_rows_in_order() {
        let expected = [
            (
                "benzoquinone",
                "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
                25,
            ),
            ("oxim", "[#6]=[N][OH]", 4),
            ("C=O", "[#6]=,:[#8]", 2),
            ("N=O", "[#7]=,:[#8]", 2),
            ("P=O", "[#15]=,:[#8]", 2),
            ("C=hetero", "[C]=[!#1;!#6]", 1),
            ("C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2),
            ("aromatic C = exocyclic N", "[c]=!@[N]", -1),
            ("methyl", "[CX4H3]", 1),
            ("guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1),
            ("guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2),
            ("aci-nitro", "[#6]=[N+]([O-])[OH]", -4),
        ];
        let terms = default_tautomer_score_terms();

        assert_eq!(terms.len(), expected.len());
        for (index, (term, &(name, smarts, score))) in terms.iter().zip(&expected).enumerate() {
            assert_eq!(term.name(), name, "name at source row {index}");
            assert_eq!(term.smarts(), smarts, "SMARTS at source row {index}");
            assert_eq!(term.score(), score, "score at source row {index}");
            let expected_matcher = query(smarts);
            assert_eq!(
                term.matcher().expect("built-in matcher").atoms(),
                expected_matcher.atoms(),
                "matcher atoms at source row {index}"
            );
            assert_eq!(
                term.matcher().expect("built-in matcher").bonds(),
                expected_matcher.bonds(),
                "matcher bonds at source row {index}"
            );
        }
    }

    #[test]
    fn score_terms_compile_the_default_table_only_once() {
        let first = default_tautomer_score_terms();
        let second = default_tautomer_score_terms();

        assert!(std::ptr::eq(first, second));
        assert!(std::ptr::eq(
            first[0].matcher().expect("matcher").atoms().as_ptr(),
            second[0].matcher().expect("matcher").atoms().as_ptr()
        ));
    }

    #[test]
    fn score_terms_invalid_smarts_keep_metadata_and_an_empty_matcher() {
        let invalid = TautomerScoreTerm::new("invalid", "[", -7);

        assert_eq!(invalid.name(), "invalid");
        assert_eq!(invalid.smarts(), "[");
        assert_eq!(invalid.score(), -7);
        assert!(invalid.matcher().is_none());
    }

    #[test]
    fn score_terms_equality_uses_only_source_metadata() {
        let compiled = TautomerScoreTerm::new("same", "[#6]=[#8]", 2);
        let mut without_matcher = compiled.clone();
        without_matcher.matcher = None;

        assert_eq!(compiled, without_matcher);
        assert_ne!(
            compiled,
            TautomerScoreTerm::new("different", "[#6]=[#8]", 2)
        );
        assert_ne!(compiled, TautomerScoreTerm::new("same", "[#7]=[#8]", 2));
        assert_ne!(compiled, TautomerScoreTerm::new("same", "[#6]=[#8]", -2));
    }

    #[test]
    fn ring_score_is_zero_without_rings() {
        let molecule = Molecule::from_smiles("CCO").expect("parse acyclic molecule");

        assert_eq!(score_tautomer_rings(&molecule).expect("score rings"), 0);
    }

    #[test]
    fn ring_score_counts_an_aromatic_heteroring_without_the_all_carbon_bonus() {
        let molecule = Molecule::from_smiles("c1ccncc1").expect("parse pyridine");

        assert_eq!(score_tautomer_rings(&molecule).expect("score rings"), 100);
    }

    #[test]
    fn ring_score_adds_both_aromatic_and_all_carbon_weights() {
        let molecule = Molecule::from_smiles("c1ccccc1").expect("parse benzene");

        assert_eq!(score_tautomer_rings(&molecule).expect("score rings"), 250);
    }

    #[test]
    fn ring_score_traverses_every_symmetrized_fused_ring() {
        let molecule = Molecule::from_smiles("c1ccc2ccccc2c1").expect("parse naphthalene");

        assert_eq!(
            crate::symmetrize_sssr(&molecule)
                .expect("find fused rings")
                .bond_rings()
                .len(),
            2
        );
        assert_eq!(score_tautomer_rings(&molecule).expect("score rings"), 500);
    }

    #[test]
    fn ring_score_ignores_nonaromatic_rings() {
        let molecule = Molecule::from_smiles("C1CCCCC1").expect("parse cyclohexane");

        assert_eq!(score_tautomer_rings(&molecule).expect("score rings"), 0);
    }

    #[test]
    fn ring_score_matches_with_or_without_a_precomputed_symm_sssr_cache() {
        let without_cache = Molecule::from_smiles_with_sanitize("c1ccncc1", false)
            .expect("parse uncached pyridine");
        assert!(without_cache.derived_cache().rings.is_none());
        let with_cache = without_cache
            .with_assigned_rings()
            .expect("assign symmetric SSSR cache");
        assert!(
            with_cache
                .derived_cache()
                .rings
                .as_ref()
                .is_some_and(crate::RingInfo::is_symm_sssr)
        );

        assert_eq!(
            score_tautomer_rings(&without_cache).expect("score uncached"),
            100
        );
        assert_eq!(
            score_tautomer_rings(&with_cache).expect("score cached"),
            100
        );
    }

    #[test]
    fn ring_score_does_not_modify_the_source_or_materialize_its_cache() {
        let molecule =
            Molecule::from_smiles_with_sanitize("c1ccccc1", false).expect("parse uncached benzene");
        let topology_before = molecule.topology_block().clone();
        assert!(molecule.derived_cache().rings.is_none());

        assert_eq!(score_tautomer_rings(&molecule).expect("score rings"), 250);
        assert_eq!(molecule.topology_block(), &topology_before);
        assert!(molecule.derived_cache().rings.is_none());
    }

    #[test]
    fn substructure_score_is_zero_when_no_term_matches() {
        let molecule = Molecule::from_smiles("CCC").expect("parse propane");
        let terms = [TautomerScoreTerm::new("oxygen", "[#8]", 7)];

        assert_eq!(score_tautomer_substructures(&molecule, &terms), 0);
    }

    #[test]
    fn substructure_score_counts_one_and_multiple_matches() {
        let ethanol = Molecule::from_smiles("CCO").expect("parse ethanol");
        let oxygen = [TautomerScoreTerm::new("oxygen", "[#8]", 7)];
        let carbon = [TautomerScoreTerm::new("carbon", "[#6]", 2)];

        assert_eq!(score_tautomer_substructures(&ethanol, &oxygen), 7);
        assert_eq!(score_tautomer_substructures(&ethanol, &carbon), 4);
    }

    #[test]
    fn substructure_score_counts_distinct_overlapping_matches() {
        let propane = Molecule::from_smiles("CCC").expect("parse propane");
        let terms = [TautomerScoreTerm::new("carbon bond", "[#6]-[#6]", 3)];

        assert_eq!(score_tautomer_substructures(&propane, &terms), 6);
    }

    #[test]
    fn substructure_score_evaluates_every_builtin_term_with_shared_match_semantics() {
        let molecule = Molecule::from_smiles("CC(=O)NC(=N)N.c1ccccc1")
            .expect("parse representative scoring molecule");

        for term in default_tautomer_score_terms() {
            let matcher = term.matcher().expect("built-in matcher compiles");
            let expected = (crate::get_substruct_matches(&molecule, matcher).len() as i32)
                .wrapping_mul(term.score());
            assert_eq!(
                score_tautomer_substructures(&molecule, std::slice::from_ref(term)),
                expected,
                "built-in term {}",
                term.name()
            );
        }
    }

    #[test]
    fn substructure_score_accumulates_custom_positive_and_negative_terms_in_source_order() {
        let molecule = Molecule::from_smiles("CCO").expect("parse ethanol");
        let terms = [
            TautomerScoreTerm::new("carbon", "[#6]", 3),
            TautomerScoreTerm::new("oxygen", "[#8]", -5),
            TautomerScoreTerm::new("no match", "[#15]", 100),
        ];

        assert_eq!(score_tautomer_substructures(&molecule, &terms), 1);
        assert_eq!(
            terms
                .iter()
                .map(TautomerScoreTerm::name)
                .collect::<Vec<_>>(),
            ["carbon", "oxygen", "no match"]
        );
    }

    #[test]
    fn substructure_score_skips_invalid_matchers_without_dropping_other_terms() {
        let molecule = Molecule::from_smiles("CCO").expect("parse ethanol");
        let terms = [
            TautomerScoreTerm::new("invalid", "[", i32::MAX),
            TautomerScoreTerm::new("oxygen", "[#8]", -5),
        ];

        assert_eq!(score_tautomer_substructures(&molecule, &terms), -5);
    }

    #[test]
    fn substructure_score_retains_source_int_accumulation_boundaries() {
        let methane = Molecule::from_smiles("C").expect("parse methane");
        let terms = [
            TautomerScoreTerm::new("maximum", "[#6]", i32::MAX),
            TautomerScoreTerm::new("one", "[#6]", 1),
        ];

        assert_eq!(score_tautomer_substructures(&methane, &terms), i32::MIN);
    }

    #[test]
    fn aggregate_score_penalizes_implicit_hydrogens_on_each_source_element() {
        for (smiles, expected) in [("P", -3), ("S", -2), ("[SeH2]", -2), ("[TeH3]", -3)] {
            let molecule = Molecule::from_smiles(smiles).expect("parse heteroatom");
            assert_eq!(
                score_tautomer_hetero_hydrogens(&molecule).expect("score hetero hydrogens"),
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn aggregate_score_counts_bracket_explicit_but_not_neighbor_isotopic_hydrogen() {
        let explicit = Molecule::from_smiles("[PH3]").expect("parse explicit phosphorus Hs");
        let isotopic_neighbor =
            Molecule::from_smiles("[2H]P").expect("parse isotopic hydrogen neighbor");

        assert_eq!(score_tautomer_hetero_hydrogens(&explicit).unwrap(), -3);
        assert_eq!(
            score_tautomer_hetero_hydrogens(&isotopic_neighbor).unwrap(),
            -2
        );
    }

    #[test]
    fn aggregate_score_counts_hydrogens_on_charged_penalized_atoms() {
        for (smiles, expected) in [("[PH4+]", -4), ("[SH-]", -1)] {
            let molecule = Molecule::from_smiles(smiles).expect("parse charged heteroatom");
            assert_eq!(
                score_tautomer_hetero_hydrogens(&molecule).expect("score charged atom"),
                expected,
                "{smiles}"
            );
        }
    }

    #[test]
    fn aggregate_score_excludes_hydrogen_bearing_elements_outside_the_source_set() {
        for smiles in ["C", "N", "O", "F", "Cl"] {
            let molecule = Molecule::from_smiles(smiles).expect("parse excluded element");
            assert_eq!(
                score_tautomer_hetero_hydrogens(&molecule).expect("score excluded element"),
                0,
                "{smiles}"
            );
        }
    }

    #[test]
    fn aggregate_score_reports_missing_valence_only_when_a_penalized_atom_needs_it() {
        let phosphorus =
            Molecule::from_smiles_with_sanitize("P", false).expect("parse unsanitized phosphorus");
        let oxygen =
            Molecule::from_smiles_with_sanitize("O", false).expect("parse unsanitized oxygen");

        assert_eq!(
            score_tautomer_hetero_hydrogens(&phosphorus),
            Err(TautomerScoreError::Valence(
                crate::ValenceError::ImplicitValenceCacheNotInitialized {
                    atom: AtomId::new(0)
                }
            ))
        );
        assert_eq!(score_tautomer_hetero_hydrogens(&oxygen).unwrap(), 0);
    }

    #[test]
    fn aggregate_score_exposes_exact_components_and_signed_sum() {
        let benzene = Molecule::from_smiles("c1ccccc1").expect("parse benzene");
        let carbon_terms = [TautomerScoreTerm::new("carbon", "[#6]", 2)];
        let score = score_tautomer_with_terms(&benzene, &carbon_terms).expect("score benzene");

        assert_eq!(score.ring(), 250);
        assert_eq!(score.substructure(), 12);
        assert_eq!(score.hetero_hydrogen(), 0);
        assert_eq!(score.total(), 262);
    }

    #[test]
    fn aggregate_score_default_delegate_uses_the_single_builtin_term_table() {
        let molecule = Molecule::from_smiles("CC(=O)N").expect("parse acetamide");
        let delegated = score_tautomer(&molecule).expect("score with defaults");
        let explicit = score_tautomer_with_terms(&molecule, default_tautomer_score_terms())
            .expect("score with explicit defaults");

        assert_eq!(delegated, explicit);
        assert_eq!(
            delegated.total(),
            delegated
                .ring()
                .wrapping_add(delegated.substructure())
                .wrapping_add(delegated.hetero_hydrogen())
        );
    }

    fn marked_atoms(indices: impl IntoIterator<Item = usize>) -> BTreeSet<AtomId> {
        indices.into_iter().map(AtomId::new).collect()
    }

    fn marked_bonds(indices: impl IntoIterator<Item = usize>) -> BTreeSet<BondId> {
        indices.into_iter().map(BondId::new).collect()
    }

    #[test]
    fn stereo_and_isotopic_hydrogens_sp2_and_remove_sp3_clear_chiral_and_cip_state() {
        let mut source = Molecule::from_smiles("[C@H](F)(Cl)Br").expect("chiral source");
        source.topology_block_mut().atoms[0].set_prop("_CIPCode", "R");
        let source_before = source.clone();

        for (hybridization, remove_sp3_stereo) in
            [(Hybridization::Sp2, false), (Hybridization::Sp3, true)]
        {
            let mut tautomer = source.clone();
            tautomer.topology_block_mut().atoms[0].set_hybridization(hybridization);
            let options = TautomerOptions::default()
                .with_remove_sp3_stereo(remove_sp3_stereo)
                .with_reassign_stereo(false);
            let changed = set_tautomer_stereo_and_isotopic_hydrogens(
                &source,
                &mut tautomer,
                &marked_atoms([0]),
                &BTreeSet::new(),
                options,
            )
            .expect("apply atom stereo transition");

            assert!(changed);
            assert_eq!(tautomer.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
            assert_eq!(tautomer.atoms()[0].prop("_CIPCode"), None);
            assert_eq!(tautomer.prop("_StereochemDone"), Some("1"));
            assert_eq!(source, source_before);
        }
    }

    #[test]
    fn stereo_and_isotopic_hydrogens_sp3_restore_copies_source_tag_and_present_cip_only() {
        let mut source = Molecule::from_smiles("[C@H](F)(Cl)Br").expect("chiral source");
        source.topology_block_mut().atoms[0].set_prop("_CIPCode", "S");
        let mut tautomer = source.clone();
        let atom = &mut tautomer.topology_block_mut().atoms[0];
        atom.set_hybridization(Hybridization::Sp3);
        atom.set_chiral_tag(match source.atoms()[0].chiral_tag() {
            ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
            _ => ChiralTag::TetrahedralCw,
        });
        atom.set_prop("_CIPCode", "R");

        let changed = set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut tautomer,
            &marked_atoms([0]),
            &BTreeSet::new(),
            TautomerOptions::default()
                .with_remove_sp3_stereo(false)
                .with_reassign_stereo(false),
        )
        .expect("restore source stereo");

        assert!(changed);
        assert_eq!(
            tautomer.atoms()[0].chiral_tag(),
            source.atoms()[0].chiral_tag()
        );
        assert_eq!(tautomer.atoms()[0].prop("_CIPCode"), Some("S"));

        source.topology_block_mut().atoms[0].clear_prop("_CIPCode");
        tautomer.topology_block_mut().atoms[0].set_prop("_CIPCode", "R");
        set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut tautomer,
            &marked_atoms([0]),
            &BTreeSet::new(),
            TautomerOptions::default()
                .with_remove_sp3_stereo(false)
                .with_reassign_stereo(false),
        )
        .expect("source-absent CIP branch");
        assert_eq!(tautomer.atoms()[0].prop("_CIPCode"), Some("R"));
    }

    #[test]
    fn stereo_and_isotopic_hydrogens_remove_or_zero_total_h_clear_tracked_isotopes() {
        let mut builder = crate::MoleculeBuilder::new();
        let atom = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(1)
                .with_tracked_isotopic_hydrogens(vec![2, 3]),
        );
        let source = builder.build().expect("isotopic-H source");

        let mut retained = source.clone();
        set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut retained,
            &marked_atoms([atom.index()]),
            &BTreeSet::new(),
            TautomerOptions::default()
                .with_remove_isotopic_hydrogens(false)
                .with_reassign_stereo(false),
        )
        .expect("retain tracked isotopes");
        assert_eq!(retained.atoms()[0].tracked_isotopic_hydrogens(), &[2, 3]);

        let mut removed_by_option = source.clone();
        set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut removed_by_option,
            &marked_atoms([0]),
            &BTreeSet::new(),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect("remove tracked isotopes by option");
        assert!(
            removed_by_option.atoms()[0]
                .tracked_isotopic_hydrogens()
                .is_empty()
        );

        let mut removed_at_zero_h = source.clone();
        removed_at_zero_h.topology_block_mut().atoms[0].set_explicit_hydrogens(0);
        set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut removed_at_zero_h,
            &marked_atoms([0]),
            &BTreeSet::new(),
            TautomerOptions::default()
                .with_remove_isotopic_hydrogens(false)
                .with_reassign_stereo(false),
        )
        .expect("remove tracked isotopes at zero total H");
        assert!(
            removed_at_zero_h.atoms()[0]
                .tracked_isotopic_hydrogens()
                .is_empty()
        );
    }

    #[test]
    fn stereo_and_isotopic_hydrogens_non_ring_double_bond_removal_sets_any_and_clears_dirs() {
        let source = Molecule::from_smiles("F/C=C/Cl").expect("E/Z source");
        let source_before = source.clone();
        let double_bond = source
            .bonds()
            .iter()
            .find(|bond| bond.order() == BondOrder::Double)
            .expect("double bond")
            .id();
        assert!(is_stereo_beyond_any(
            source.bonds()[double_bond.index()].stereo()
        ));
        let directional = source
            .bonds()
            .iter()
            .filter(|bond| {
                matches!(
                    bond.direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                )
            })
            .map(crate::Bond::id)
            .collect::<Vec<_>>();
        assert!(!directional.is_empty());
        let mut tautomer = source.clone();

        set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut tautomer,
            &BTreeSet::new(),
            &marked_bonds([double_bond.index()]),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect("remove non-ring bond stereo");

        assert_eq!(
            tautomer.bonds()[double_bond.index()].stereo(),
            BondStereo::Any
        );
        assert_eq!(tautomer.bonds()[double_bond.index()].stereo_atoms(), None);
        for bond_id in directional {
            assert_eq!(
                tautomer.bonds()[bond_id.index()].direction(),
                BondDirection::None
            );
        }
        assert_eq!(source, source_before);
    }

    #[test]
    fn stereo_and_isotopic_hydrogens_preservation_restores_bond_stereo_atoms_and_dirs() {
        let source = Molecule::from_smiles("F/C=C/Cl").expect("E/Z source");
        let double_bond = source
            .bonds()
            .iter()
            .find(|bond| bond.order() == BondOrder::Double)
            .expect("double bond")
            .id();
        let mut tautomer = source.clone();
        tautomer.topology_block_mut().bonds[double_bond.index()].set_stereo_atoms(None);
        tautomer.topology_block_mut().bonds[double_bond.index()].set_stereo(BondStereo::Any);
        for bond in &mut tautomer.topology_block_mut().bonds {
            if matches!(
                bond.direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            ) {
                bond.set_direction(BondDirection::None);
            }
        }

        let changed = set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut tautomer,
            &BTreeSet::new(),
            &marked_bonds([double_bond.index()]),
            TautomerOptions::default()
                .with_remove_bond_stereo(false)
                .with_reassign_stereo(false),
        )
        .expect("restore source bond stereo");

        assert!(changed);
        assert_eq!(
            tautomer.bonds()[double_bond.index()].stereo(),
            source.bonds()[double_bond.index()].stereo()
        );
        assert_eq!(
            tautomer.bonds()[double_bond.index()].stereo_atoms(),
            source.bonds()[double_bond.index()].stereo_atoms()
        );
        for (actual, expected) in tautomer.bonds().iter().zip(source.bonds()) {
            assert_eq!(actual.direction(), expected.direction());
        }
    }

    #[test]
    fn stereo_and_isotopic_hydrogens_ring_and_non_double_bonds_use_none() {
        let ring_source = Molecule::from_smiles("C1=CCCCC1").expect("ring alkene");
        let ring_bond = ring_source
            .bonds()
            .iter()
            .find(|bond| bond.order() == BondOrder::Double)
            .expect("ring double bond")
            .id();
        let mut ring_tautomer = ring_source.clone();
        ring_tautomer.topology_block_mut().bonds[ring_bond.index()].set_stereo(BondStereo::Any);
        set_tautomer_stereo_and_isotopic_hydrogens(
            &ring_source,
            &mut ring_tautomer,
            &BTreeSet::new(),
            &marked_bonds([ring_bond.index()]),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect("ring fallback");
        assert_eq!(
            ring_tautomer.bonds()[ring_bond.index()].stereo(),
            BondStereo::None
        );

        let source = Molecule::from_smiles("CC").expect("single bond");
        let mut tautomer = source.clone();
        tautomer.topology_block_mut().bonds[0].set_stereo(BondStereo::Any);
        set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut tautomer,
            &BTreeSet::new(),
            &marked_bonds([0]),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect("single bond cleanup");
        assert_eq!(tautomer.bonds()[0].stereo(), BondStereo::None);
        assert_eq!(tautomer.bonds()[0].stereo_atoms(), None);
    }

    #[test]
    fn stereo_and_isotopic_hydrogens_reassignment_reapplies_any_contract() {
        let source = Molecule::from_smiles("F/C=C/Cl").expect("E/Z source");
        let double_bond = source
            .bonds()
            .iter()
            .find(|bond| bond.order() == BondOrder::Double)
            .expect("double bond")
            .id();
        let mut tautomer = source.clone();

        set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut tautomer,
            &BTreeSet::new(),
            &marked_bonds([double_bond.index()]),
            TautomerOptions::default(),
        )
        .expect("reassign and reapply explicit undefined stereo");

        assert_eq!(
            tautomer.bonds()[double_bond.index()].stereo(),
            BondStereo::Any
        );
        assert_eq!(tautomer.bonds()[double_bond.index()].stereo_atoms(), None);
        assert_eq!(tautomer.prop("_StereochemDone"), Some("1"));
        assert!(tautomer.is_prop_computed("_StereochemDone"));
    }

    #[test]
    fn stereo_and_isotopic_hydrogens_executes_every_option_combination_without_source_mutation() {
        let source = Molecule::from_smiles("F/C=C/[C@H](Cl)Br").expect("combined source");
        let source_before = source.clone();
        let double_bond = source
            .bonds()
            .iter()
            .find(|bond| bond.order() == BondOrder::Double)
            .expect("double bond")
            .id();
        let chiral_atom = source
            .atoms()
            .iter()
            .find(|atom| atom.chiral_tag() != ChiralTag::Unspecified)
            .expect("chiral atom")
            .id();

        for bits in 0_u8..16 {
            let options = TautomerOptions::default()
                .with_remove_sp3_stereo(bits & 1 != 0)
                .with_remove_bond_stereo(bits & 2 != 0)
                .with_remove_isotopic_hydrogens(bits & 4 != 0)
                .with_reassign_stereo(bits & 8 != 0);
            let mut tautomer = source.clone();
            set_tautomer_stereo_and_isotopic_hydrogens(
                &source,
                &mut tautomer,
                &marked_atoms([chiral_atom.index()]),
                &marked_bonds([double_bond.index()]),
                options,
            )
            .unwrap_or_else(|error| panic!("option mask {bits:#06b}: {error}"));
            assert_eq!(source, source_before, "option mask {bits:#06b}");
            if options.remove_bond_stereo() {
                assert_eq!(
                    tautomer.bonds()[double_bond.index()].stereo(),
                    BondStereo::Any,
                    "option mask {bits:#06b}"
                );
            }
            if !options.reassign_stereo() {
                assert_eq!(tautomer.prop("_StereochemDone"), Some("1"));
                assert_eq!(
                    tautomer.is_prop_computed("_StereochemDone"),
                    source.is_prop_computed("_StereochemDone")
                );
            }
        }
    }

    fn initialization_plan(molecule: &Molecule) -> TautomerInitializationPlan {
        plan_tautomer_initialization(MoleculeReadParts::from_molecule(molecule))
            .expect("tautomer initialization")
    }

    #[test]
    fn enumeration_initialization_updates_only_missing_valence_and_symm_sssr_caches() {
        let unsanitized = Molecule::from_smiles_with_sanitize("c1ccccc1", false)
            .expect("parse unsanitized benzene");
        assert!(unsanitized.derived_cache().valence.is_none());
        assert!(unsanitized.derived_cache().rings.is_none());
        let uncached = initialization_plan(&unsanitized);
        assert!(uncached.valence_update.is_some());
        assert!(
            uncached
                .rings_update
                .as_ref()
                .is_some_and(crate::RingInfo::is_symm_sssr)
        );

        let cached = Molecule::from_smiles("c1ccccc1").expect("parse sanitized benzene");
        assert!(cached.derived_cache().valence.is_some());
        assert!(
            cached
                .derived_cache()
                .rings
                .as_ref()
                .is_some_and(crate::RingInfo::is_symm_sssr)
        );
        let retained = initialization_plan(&cached);
        assert!(retained.valence_update.is_none());
        assert!(retained.rings_update.is_none());
    }

    #[test]
    fn enumeration_initialization_canonical_kekulizes_aromatic_input_without_clearing_flags() {
        let molecule = Molecule::from_smiles("c1ccccc1").expect("parse benzene");
        let plan = initialization_plan(&molecule);
        assert_eq!(plan.canonical_smiles, "c1ccccc1");

        let mut topology = molecule.topology_block().clone();
        crate::kekulize::apply_kekulize_assignment(&mut topology, &plan.kekulize_assignment);
        assert!(topology.atoms.iter().all(crate::Atom::is_aromatic));
        assert!(topology.bonds.iter().all(crate::Bond::is_aromatic));
        assert_eq!(
            topology
                .bonds
                .iter()
                .filter(|bond| bond.order() == BondOrder::Double)
                .count(),
            3
        );
        assert_eq!(
            topology
                .bonds
                .iter()
                .filter(|bond| bond.order() == BondOrder::Single)
                .count(),
            3
        );
    }

    #[test]
    fn enumeration_initialization_handles_empty_and_disconnected_inputs() {
        let empty = initialization_plan(&Molecule::new());
        assert_eq!(empty.canonical_smiles, "");
        assert_eq!(empty.atom_count, 0);
        assert_eq!(empty.bond_count, 0);

        let disconnected = Molecule::from_smiles("O.CC").expect("parse disconnected molecule");
        let plan = initialization_plan(&disconnected);
        assert_eq!(plan.canonical_smiles, "CC.O");
        assert_eq!(plan.atom_count, 3);
        assert_eq!(plan.bond_count, 1);
    }

    #[test]
    fn enumeration_initialization_is_independent_of_input_atom_order() {
        let first = Molecule::from_smiles("CC(=O)O").expect("parse first atom order");
        let second = Molecule::from_smiles("OC(=O)C").expect("parse second atom order");

        let first_plan = initialization_plan(&first);
        let second_plan = initialization_plan(&second);
        assert_eq!(first_plan.canonical_smiles, "CC(=O)O");
        assert_eq!(second_plan.canonical_smiles, first_plan.canonical_smiles);
    }

    #[test]
    fn enumeration_initialization_inserts_one_unfinished_candidate_in_key_order() {
        let molecule = Molecule::from_smiles("OC(=O)C").expect("parse candidate");
        let plan = initialization_plan(&molecule);
        let candidates = plan.into_candidate_map(17usize, 23usize);

        assert_eq!(
            candidates.keys().map(String::as_str).collect::<Vec<_>>(),
            ["CC(=O)O"]
        );
        let candidate = &candidates["CC(=O)O"];
        assert_eq!(candidate.tautomer, Some(17));
        assert_eq!(candidate.kekulized, Some(23));
        assert_eq!(candidate.num_modified_atoms, 0);
        assert_eq!(candidate.num_modified_bonds, 0);
        assert!(!candidate.done);
    }

    #[test]
    fn enumeration_initialization_reports_invalid_aromatic_input_as_a_structured_error() {
        let molecule = Molecule::from_smiles_with_sanitize("c1cccc1", false)
            .expect("parse deliberately unkekulizable odd aromatic ring");
        let error = plan_tautomer_initialization_with_canonical_smiles(
            MoleculeReadParts::from_molecule(&molecule),
            "c1cccc1".to_owned(),
        )
        .expect_err("initial canonical kekulization must reject the invalid aromatic state");

        assert!(matches!(error, TautomerInitializationError::Kekulize(_)));
    }

    #[test]
    fn enumeration_initialization_never_mutates_the_input() {
        for smiles in ["c1ccccc1", "CC(=O)O", "O.CC"] {
            let molecule = Molecule::from_smiles(smiles).expect("parse immutable input");
            let before = molecule.clone();
            let _ = initialization_plan(&molecule);
            assert_eq!(molecule, before, "input {smiles}");
        }
    }

    fn builtin_transform(name: &str) -> TautomerTransform {
        TautomerCatalog::current()
            .expect("load built-in tautomer catalog")
            .transforms()
            .iter()
            .find(|transform| transform.name() == name)
            .unwrap_or_else(|| panic!("missing built-in tautomer transform {name}"))
            .clone()
    }

    fn prepared_tautomer_candidate(molecule: &Molecule) -> Molecule {
        let plan = initialization_plan(molecule);
        let mut topology = molecule.topology_block().clone();
        crate::kekulize::apply_kekulize_assignment(&mut topology, &plan.kekulize_assignment);
        let mut derived_cache = molecule.derived_cache().clone();
        if let Some(valence) = plan.valence_update {
            derived_cache.valence = Some(valence);
        }
        if let Some(rings) = plan.rings_update {
            derived_cache.rings = Some(rings);
        }
        Molecule::from_operation_blocks(
            topology,
            molecule.coordinate_block().clone(),
            molecule.properties().clone(),
            derived_cache,
            molecule.capabilities_block(),
        )
        .expect("prepared tautomer candidate satisfies molecule invariants")
    }

    fn first_transform_match(
        candidate: &Molecule,
        transform: &TautomerTransform,
    ) -> crate::SubstructMatchResult {
        crate::get_substruct_matches(candidate, transform.query())
            .into_iter()
            .next()
            .unwrap_or_else(|| panic!("{} must match focused fixture", transform.name()))
    }

    fn apply_single_transform(
        source: &Molecule,
        candidate: &Molecule,
        transform: &TautomerTransform,
        matched: &crate::SubstructMatchResult,
        modified_atoms: &BTreeSet<AtomId>,
        modified_bonds: &BTreeSet<BondId>,
        existing_smiles: &BTreeSet<String>,
        options: TautomerOptions,
    ) -> Result<TautomerTransformAttempt, TautomerTransformApplicationError> {
        apply_tautomer_transform_match(
            MoleculeReadParts::from_molecule(source),
            MoleculeReadParts::from_molecule(candidate),
            transform,
            matched,
            modified_atoms,
            modified_bonds,
            existing_smiles,
            options,
        )
    }

    fn molecule_from_product_plan(source: &Molecule, product: &TautomerProductPlan) -> Molecule {
        Molecule::from_operation_blocks(
            product.topology.clone(),
            source.coordinate_block().clone(),
            product.properties.clone(),
            product.derived_cache.clone(),
            source.capabilities_block(),
        )
        .expect("single-transform product satisfies molecule invariants")
    }

    #[test]
    fn single_transform_application_reproduces_every_builtin_bond_edit_shape() {
        let cases = [
            ("CC=O", "1,3 (thio)keto/enol f", "C=CO", ""),
            ("CC=C=O", "keten/ynol f", "CC#CO", "#-"),
            ("CC#CO", "keten/ynol r", "CC=C=O", "=="),
            (
                "NC(N)=S(=O)=O",
                "formamidinesulfinic acid f",
                "N=C(N)S(=O)O",
                "=--",
            ),
            (
                "NC(=N)S(=O)O",
                "formamidinesulfinic acid r",
                "NC(N)=S(=O)=O",
                "==-",
            ),
            ("C#N", "isocyanide f", "[C-]#[NH+]", "#"),
            ("OP(O)O", "phosphonic acid f", "O=[PH](O)O", "="),
            ("[PH](=O)(O)(O)", "phosphonic acid r", "OP(O)O", "-"),
        ];

        for (input, transform_name, expected_smiles, expected_bonds) in cases {
            let source = Molecule::from_smiles(input).expect("parse edit-shape fixture");
            let source_before = source.clone();
            let candidate = prepared_tautomer_candidate(&source);
            let candidate_before = candidate.clone();
            let transform = builtin_transform(transform_name);
            assert_eq!(
                transform.bond_types(),
                string_to_bond_types(expected_bonds),
                "transform {transform_name}"
            );
            let matched = first_transform_match(&candidate, &transform);

            let attempt = apply_single_transform(
                &source,
                &candidate,
                &transform,
                &matched,
                &BTreeSet::new(),
                &BTreeSet::new(),
                &BTreeSet::new(),
                TautomerOptions::default().with_reassign_stereo(false),
            )
            .unwrap_or_else(|error| panic!("transform {transform_name}: {error}"));
            let TautomerTransformAttempt::Product(product) = attempt else {
                panic!("transform {transform_name} did not produce a unique tautomer");
            };
            let molecule = molecule_from_product_plan(&source, &product);

            assert_eq!(
                product.canonical_smiles, expected_smiles,
                "{transform_name}"
            );
            assert_eq!(
                MoleculeReadParts::from_molecule(&molecule)
                    .canonical_isomeric_smiles()
                    .expect("write product canonical SMILES"),
                expected_smiles,
                "{transform_name}"
            );
            assert_eq!(source, source_before, "source changed for {transform_name}");
            assert_eq!(
                candidate, candidate_before,
                "candidate changed for {transform_name}"
            );
        }
    }

    #[test]
    fn single_transform_application_uses_match_endpoint_order_and_unions_modified_sets() {
        let source = Molecule::from_smiles("CC=O").expect("parse acetaldehyde");
        let candidate = prepared_tautomer_candidate(&source);
        let transform = builtin_transform("1,3 (thio)keto/enol f");
        let matched = first_transform_match(&candidate, &transform);
        let donor = AtomId::new(matched.atom_mapping[0]);
        let acceptor = AtomId::new(*matched.atom_mapping.last().expect("acceptor"));
        let retained_atom = AtomId::new(matched.atom_mapping[1]);
        let retained_bond = BondId::new(matched.bond_mapping[0]);

        let attempt = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &matched,
            &BTreeSet::from([retained_atom]),
            &BTreeSet::from([retained_bond]),
            &BTreeSet::new(),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect("apply ordered endpoint transform");
        let TautomerTransformAttempt::Product(product) = attempt else {
            panic!("ordered endpoint transform must produce a product");
        };

        assert_eq!(
            product.topology.atoms[donor.index()].explicit_hydrogens(),
            2
        );
        assert_eq!(
            product.topology.atoms[acceptor.index()].explicit_hydrogens(),
            1
        );
        assert!(product.topology.atoms[donor.index()].no_implicit());
        assert!(product.topology.atoms[acceptor.index()].no_implicit());
        assert_eq!(
            product.modified_atoms,
            BTreeSet::from([donor, retained_atom, acceptor])
        );
        assert!(product.modified_bonds.contains(&retained_bond));
        assert!(
            matched
                .bond_mapping
                .iter()
                .all(|&bond| product.modified_bonds.contains(&BondId::new(bond)))
        );
    }

    #[test]
    fn single_transform_application_handles_explicit_implicit_and_isotopic_hydrogens() {
        for input in ["CC=O", "[CH3]C=O"] {
            let source = Molecule::from_smiles(input).expect("parse hydrogen fixture");
            let candidate = prepared_tautomer_candidate(&source);
            let transform = builtin_transform("1,3 (thio)keto/enol f");
            let matched = first_transform_match(&candidate, &transform);
            let attempt = apply_single_transform(
                &source,
                &candidate,
                &transform,
                &matched,
                &BTreeSet::new(),
                &BTreeSet::new(),
                &BTreeSet::new(),
                TautomerOptions::default().with_reassign_stereo(false),
            )
            .expect("apply implicit/explicit hydrogen transform");
            let TautomerTransformAttempt::Product(product) = attempt else {
                panic!("hydrogen transform must produce a product");
            };
            assert_eq!(product.canonical_smiles, "C=CO");
        }

        let mut source = Molecule::from_smiles("CC=O").expect("parse isotopic-H fixture");
        source.topology_block_mut().atoms[0].set_tracked_isotopic_hydrogens(vec![2, 3]);
        let candidate = prepared_tautomer_candidate(&source);
        let transform = builtin_transform("1,3 (thio)keto/enol f");
        let matched = first_transform_match(&candidate, &transform);
        for (remove, expected) in [(true, &[][..]), (false, &[2, 3][..])] {
            let attempt = apply_single_transform(
                &source,
                &candidate,
                &transform,
                &matched,
                &BTreeSet::new(),
                &BTreeSet::new(),
                &BTreeSet::new(),
                TautomerOptions::default()
                    .with_remove_isotopic_hydrogens(remove)
                    .with_reassign_stereo(false),
            )
            .expect("apply isotopic-H transform");
            let TautomerTransformAttempt::Product(product) = attempt else {
                panic!("isotopic-H transform must produce a product");
            };
            assert_eq!(
                product.topology.atoms[0].tracked_isotopic_hydrogens(),
                expected
            );
        }
    }

    #[test]
    fn single_transform_application_applies_source_ordered_charge_deltas() {
        let source = Molecule::from_smiles("C#N").expect("parse hydrogen cyanide");
        let candidate = prepared_tautomer_candidate(&source);
        let transform = builtin_transform("isocyanide f");
        let matched = first_transform_match(&candidate, &transform);
        assert_eq!(transform.charges(), &[-1, 1]);

        let attempt = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect("apply isocyanide charge transform");
        let TautomerTransformAttempt::Product(product) = attempt else {
            panic!("isocyanide transform must produce a product");
        };
        assert_eq!(
            product.topology.atoms[matched.atom_mapping[0]].formal_charge(),
            -1
        );
        assert_eq!(
            product.topology.atoms[matched.atom_mapping[1]].formal_charge(),
            1
        );
        assert_eq!(product.canonical_smiles, "[C-]#[NH+]");
    }

    #[test]
    fn single_transform_application_records_bonds_changed_only_by_sanitization() {
        let source = Molecule::from_smiles("Cc1nc2ccccc2[nH]1")
            .expect("parse source-comment sanitization fixture");
        let source_before = source.clone();
        let candidate = prepared_tautomer_candidate(&source);
        let catalog = TautomerCatalog::current().expect("load catalog");
        let mut witnessed = None;

        'transforms: for transform in catalog.transforms() {
            for matched in crate::get_substruct_matches(&candidate, transform.query()) {
                let directly_edited = matched
                    .bond_mapping
                    .iter()
                    .copied()
                    .map(BondId::new)
                    .collect::<BTreeSet<_>>();
                let attempt = apply_single_transform(
                    &source,
                    &candidate,
                    transform,
                    &matched,
                    &BTreeSet::new(),
                    &BTreeSet::new(),
                    &BTreeSet::new(),
                    TautomerOptions::default().with_reassign_stereo(false),
                )
                .unwrap_or_else(|error| panic!("{}: {error}", transform.name()));
                if let TautomerTransformAttempt::Product(product) = attempt {
                    let sanitize_only = product
                        .modified_bonds
                        .difference(&directly_edited)
                        .copied()
                        .collect::<BTreeSet<_>>();
                    if !sanitize_only.is_empty() {
                        witnessed = Some((product, sanitize_only));
                        break 'transforms;
                    }
                }
            }
        }

        let (product, sanitize_only) = witnessed
            .expect("source-comment fixture must expose a sanitization-only bond-order change");
        assert!(sanitize_only.iter().all(|bond| {
            source.bonds()[bond.index()].order() != product.topology.bonds[bond.index()].order()
        }));
        let _ = molecule_from_product_plan(&source, &product);
        assert_eq!(source, source_before);
    }

    #[test]
    fn single_transform_application_reports_duplicate_after_stereo_and_before_product_kekulization()
    {
        let source = Molecule::from_smiles("CC=O").expect("parse duplicate fixture");
        let candidate = prepared_tautomer_candidate(&source);
        let transform = builtin_transform("1,3 (thio)keto/enol f");
        let matched = first_transform_match(&candidate, &transform);
        let first = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect("compute duplicate key");
        let TautomerTransformAttempt::Product(first) = first else {
            panic!("first transform must be unique");
        };

        let duplicate = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::from([first.canonical_smiles.clone()]),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect("duplicate is a non-error branch");
        let TautomerTransformAttempt::Duplicate {
            canonical_smiles,
            modified_atoms,
            modified_bonds,
        } = duplicate
        else {
            panic!("existing canonical key must select duplicate branch");
        };
        assert_eq!(canonical_smiles, first.canonical_smiles);
        assert!(!modified_atoms.is_empty());
        assert!(!modified_bonds.is_empty());
    }

    #[test]
    fn single_transform_application_catches_only_the_kekulize_failure_branch() {
        let source = Molecule::from_smiles_with_sanitize("c1cccc1.CC", false)
            .expect("parse odd aromatic ring fixture");
        let candidate = source
            .with_assigned_valence_strict(false)
            .expect("prepare non-strict valence cache");
        let transform = TautomerTransform::new(
            "focused kekulize failure",
            query("[C]-[C]"),
            Vec::new(),
            Vec::new(),
        )
        .expect("construct focused transform");
        let matched = first_transform_match(&candidate, &transform);
        let source_before = source.clone();
        let candidate_before = candidate.clone();

        let attempt = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect("Kekulize failure is the source-defined recoverable branch");
        let TautomerTransformAttempt::RecoverableKekulizeFailure {
            modified_atoms,
            modified_bonds,
        } = attempt
        else {
            panic!("odd aromatic ring must fail in partial-sanitize Kekulize stage");
        };
        assert_eq!(modified_atoms.len(), 2);
        assert_eq!(modified_bonds.len(), 1);
        assert_eq!(source, source_before);
        assert_eq!(candidate, candidate_before);
    }

    #[test]
    fn single_transform_application_propagates_non_kekulize_sanitize_failures() {
        let source = Molecule::from_smiles("CCC").expect("parse propagation fixture");
        let mut candidate = source.clone();
        candidate.topology_block_mut().bonds[1].set_order(BondOrder::Other);
        let transform = TautomerTransform::new(
            "focused sanitize propagation",
            query("[C]-[C]"),
            Vec::new(),
            Vec::new(),
        )
        .expect("construct focused transform");
        let matched = crate::SubstructMatchResult {
            atom_mapping: vec![0, 1],
            bond_mapping: vec![0],
        };
        let source_before = source.clone();
        let candidate_before = candidate.clone();

        let error = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            TautomerOptions::default().with_reassign_stereo(false),
        )
        .expect_err("bad unrelated bond type must propagate from property-cache sanitization");
        assert!(matches!(
            error,
            TautomerTransformApplicationError::Sanitize(crate::OperationError::Sanitize {
                source: crate::SanitizeError {
                    step: crate::SanitizeStep::Properties,
                    ..
                },
                ..
            })
        ));
        assert_eq!(source, source_before);
        assert_eq!(candidate, candidate_before);
    }

    #[test]
    fn single_transform_application_rejects_invalid_mappings_without_mutating_inputs() {
        let source = Molecule::from_smiles("CC=O").expect("parse mapping fixture");
        let candidate = prepared_tautomer_candidate(&source);
        let transform = builtin_transform("1,3 (thio)keto/enol f");
        let source_before = source.clone();
        let candidate_before = candidate.clone();
        let error = apply_single_transform(
            &source,
            &candidate,
            &transform,
            &crate::SubstructMatchResult {
                atom_mapping: vec![0, 1],
                bond_mapping: vec![0, 1],
            },
            &BTreeSet::new(),
            &BTreeSet::new(),
            &BTreeSet::new(),
            TautomerOptions::default(),
        )
        .expect_err("short atom mapping must return a structured error");

        assert!(matches!(
            error,
            TautomerTransformApplicationError::AtomMappingCount {
                expected: 3,
                actual: 2
            }
        ));
        assert_eq!(source, source_before);
        assert_eq!(candidate, candidate_before);
    }

    fn expansion_transform(name: &str) -> TautomerTransform {
        TautomerTransform::new(name, query("[C]-[C]"), Vec::new(), Vec::new())
            .expect("construct expansion-control transform")
    }

    fn expansion_candidate(handle: &str, done: bool) -> TautomerCandidate<String> {
        TautomerCandidate {
            tautomer: Some(format!("{handle}-tautomer")),
            kekulized: Some(handle.to_owned()),
            num_modified_atoms: 0,
            num_modified_bonds: 0,
            done,
        }
    }

    fn expansion_state(entries: &[(&str, &str, bool)]) -> TautomerExpansionState<String> {
        TautomerExpansionState {
            candidates: entries
                .iter()
                .map(|&(key, handle, done)| (key.to_owned(), expansion_candidate(handle, done)))
                .collect(),
            modified_atoms: BTreeSet::new(),
            modified_bonds: BTreeSet::new(),
            status: TautomerEnumerationStatus::Completed,
            num_transforms: 0,
        }
    }

    fn expansion_match(tag: usize) -> crate::SubstructMatchResult {
        crate::SubstructMatchResult {
            atom_mapping: vec![tag, tag],
            bond_mapping: vec![tag],
        }
    }

    fn expansion_product(key: &str, tag: usize) -> TautomerExpansionAttempt<String> {
        TautomerExpansionAttempt::Product(TautomerExpandedProduct {
            tautomer: format!("{key}-tautomer"),
            kekulized: key.to_owned(),
            canonical_smiles: key.to_owned(),
            modified_atoms: BTreeSet::from([AtomId::new(tag)]),
            modified_bonds: BTreeSet::from([BondId::new(tag)]),
        })
    }

    #[test]
    fn enumeration_expansion_visits_later_insertions_in_pass_and_defers_earlier_insertions() {
        let transforms = [expansion_transform("ordered")];
        let mut state = expansion_state(&[("m", "m", false)]);
        let mut visited = Vec::new();

        let first = expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerOptions::default(),
            |handle, _| {
                visited.push(handle.clone());
                Ok::<_, &'static str>(match handle.as_str() {
                    "m" | "z" => vec![expansion_match(0)],
                    _ => Vec::new(),
                })
            },
            |handle, _, _, _, _, _| {
                Ok::<_, &'static str>(match handle.as_str() {
                    "m" => expansion_product("z", 0),
                    "z" => expansion_product("a", 1),
                    _ => unreachable!("only matching handles are applied"),
                })
            },
            |_| Ok(true),
        )
        .expect("first ordered expansion pass");

        assert!(!first.bailed_out);
        assert_eq!(visited, ["m", "z"]);
        assert_eq!(
            state
                .candidates
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["a", "m", "z"]
        );
        assert!(!state.candidates["a"].done);
        assert!(state.candidates["m"].done);
        assert!(state.candidates["z"].done);

        visited.clear();
        let second = expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerOptions::default(),
            |handle, _| {
                visited.push(handle.clone());
                Ok::<_, &'static str>(Vec::new())
            },
            |_, _, _, _, _, _| -> Result<_, &'static str> {
                unreachable!("the deferred key has no match")
            },
            |_| Ok(true),
        )
        .expect("second ordered expansion pass");
        assert!(!second.bailed_out);
        assert_eq!(visited, ["a"]);
        assert!(state.candidates.values().all(|candidate| candidate.done));
    }

    #[test]
    fn enumeration_expansion_counts_one_transform_for_multiple_and_duplicate_matches() {
        let transforms = [expansion_transform("multi-match")];
        let mut state = expansion_state(&[("m", "m", false)]);
        let mut applied = 0;
        let mut callback_sizes = Vec::new();

        expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerOptions::default(),
            |handle, _| {
                Ok::<_, &'static str>(if handle == "m" {
                    vec![expansion_match(0), expansion_match(1)]
                } else {
                    Vec::new()
                })
            },
            |_, _, matched, modified_atoms, modified_bonds, existing| {
                applied += 1;
                if matched.bond_mapping[0] == 0 {
                    Ok::<_, &'static str>(expansion_product("z", 0))
                } else {
                    assert!(existing.contains("z"));
                    Ok(TautomerExpansionAttempt::Duplicate {
                        canonical_smiles: "z".to_owned(),
                        modified_atoms: modified_atoms
                            .union(&BTreeSet::from([AtomId::new(1)]))
                            .copied()
                            .collect(),
                        modified_bonds: modified_bonds
                            .union(&BTreeSet::from([BondId::new(1)]))
                            .copied()
                            .collect(),
                    })
                }
            },
            |view| {
                callback_sizes.push((view.num_transforms, view.candidates.len()));
                Ok(true)
            },
        )
        .expect("expand multiple matches");

        assert_eq!(state.num_transforms, 1);
        assert_eq!(applied, 2);
        assert_eq!(callback_sizes, [(1, 1), (1, 2)]);
        assert_eq!(state.modified_atoms, marked_atoms([0, 1]));
        assert_eq!(state.modified_bonds, marked_bonds([0, 1]));
    }

    #[test]
    fn enumeration_expansion_zero_and_exact_transform_limits_use_source_increment_order() {
        for limit in [0, 1] {
            let transforms = [expansion_transform("limited")];
            let mut state = expansion_state(&[("m", "m", false)]);
            let mut callbacks = 0;
            let mut applications = 0;
            let pass = expand_tautomer_candidates_in_source_order(
                &mut state,
                &transforms,
                TautomerOptions::default().with_max_transforms(limit),
                |_, _| Ok::<_, &'static str>(vec![expansion_match(0)]),
                |_, _, _, _, _, _| {
                    applications += 1;
                    Ok::<_, &'static str>(expansion_product("z", 0))
                },
                |_| {
                    callbacks += 1;
                    Ok(true)
                },
            )
            .expect("apply transform limit");

            assert!(pass.bailed_out, "limit {limit}");
            assert_eq!(state.num_transforms, 1, "limit {limit}");
            assert_eq!(
                state.status,
                TautomerEnumerationStatus::MaxTransformsReached,
                "limit {limit}"
            );
            assert_eq!(callbacks, 0, "limit {limit}");
            assert_eq!(applications, 0, "limit {limit}");
            assert!(state.candidates["m"].done, "limit {limit}");
        }

        let transforms = [expansion_transform("first"), expansion_transform("second")];
        let mut state = expansion_state(&[("m", "m", false)]);
        let mut applied = Vec::new();
        expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerOptions::default().with_max_transforms(2),
            |_, _| Ok::<_, &'static str>(vec![expansion_match(0)]),
            |_, transform, _, _, _, _| {
                applied.push(transform.name().to_owned());
                Ok::<_, &'static str>(TautomerExpansionAttempt::Duplicate {
                    canonical_smiles: "m".to_owned(),
                    modified_atoms: BTreeSet::new(),
                    modified_bonds: BTreeSet::new(),
                })
            },
            |_| Ok(true),
        )
        .expect("exact transform boundary");
        assert_eq!(applied, ["first"]);
        assert_eq!(state.num_transforms, 2);
        assert_eq!(
            state.status,
            TautomerEnumerationStatus::MaxTransformsReached
        );
    }

    #[test]
    fn enumeration_expansion_tautomer_limits_are_checked_before_callback_and_each_match() {
        let transforms = [expansion_transform("two matches")];
        for limit in [0, 1] {
            let mut state = expansion_state(&[("m", "m", false)]);
            let mut callbacks = 0;
            let pass = expand_tautomer_candidates_in_source_order(
                &mut state,
                &transforms,
                TautomerOptions::default().with_max_tautomers(limit),
                |_, _| Ok::<_, &'static str>(vec![expansion_match(0)]),
                |_, _, _, _, _, _| Ok::<_, &'static str>(expansion_product("z", 0)),
                |_| {
                    callbacks += 1;
                    Ok(true)
                },
            )
            .expect("apply immediate tautomer limit");
            assert!(pass.bailed_out, "limit {limit}");
            assert_eq!(
                state.status,
                TautomerEnumerationStatus::MaxTautomersReached,
                "limit {limit}"
            );
            assert_eq!(callbacks, 0, "limit {limit}");
        }

        let mut state = expansion_state(&[("m", "m", false)]);
        let mut callbacks = Vec::new();
        let mut applications = 0;
        expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerOptions::default().with_max_tautomers(2),
            |_, _| Ok::<_, &'static str>(vec![expansion_match(0), expansion_match(1)]),
            |_, _, _, _, _, _| {
                applications += 1;
                Ok::<_, &'static str>(expansion_product("z", 0))
            },
            |view| {
                callbacks.push(view.candidates.len());
                Ok(true)
            },
        )
        .expect("apply per-match tautomer limit");
        assert_eq!(callbacks, [1]);
        assert_eq!(applications, 1);
        assert_eq!(state.status, TautomerEnumerationStatus::MaxTautomersReached);
    }

    #[test]
    fn enumeration_expansion_callback_observes_preapplication_state_and_cancels_deterministically()
    {
        let transforms = [expansion_transform("callback")];
        let run = || {
            let mut state = expansion_state(&[("m", "m", false)]);
            let mut observations = Vec::new();
            let mut applications = 0;
            let pass = expand_tautomer_candidates_in_source_order(
                &mut state,
                &transforms,
                TautomerOptions::default(),
                |_, _| Ok::<_, &'static str>(vec![expansion_match(0)]),
                |_, _, _, _, _, _| {
                    applications += 1;
                    Ok::<_, &'static str>(expansion_product("z", 0))
                },
                |view| {
                    observations.push((
                        view.num_transforms,
                        view.candidates.keys().cloned().collect::<Vec<_>>(),
                        view.status,
                    ));
                    Ok(false)
                },
            )
            .expect("callback cancellation");
            (state, pass, observations, applications)
        };

        let first = run();
        let second = run();
        assert_eq!(first, second);
        assert!(first.1.bailed_out);
        assert_eq!(first.0.status, TautomerEnumerationStatus::Canceled);
        assert_eq!(
            first.2,
            [(
                1,
                vec!["m".to_owned()],
                TautomerEnumerationStatus::Completed
            )]
        );
        assert_eq!(first.3, 0);
        assert!(first.0.candidates["m"].done);
    }

    #[test]
    fn enumeration_expansion_skips_done_candidates_and_marks_later_keys_done_after_bailout() {
        let transforms = [expansion_transform("stop")];
        let mut state = expansion_state(&[("a", "a", true), ("m", "m", false), ("z", "z", false)]);
        let mut matched_handles = Vec::new();
        expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerOptions::default().with_max_transforms(1),
            |handle, _| {
                matched_handles.push(handle.clone());
                Ok::<_, &'static str>(vec![expansion_match(0)])
            },
            |_, _, _, _, _, _| -> Result<_, &'static str> {
                unreachable!("limit is checked before application")
            },
            |_| Ok(true),
        )
        .expect("bail out in ordered traversal");

        assert_eq!(matched_handles, ["m"]);
        assert!(state.candidates["a"].done);
        assert!(state.candidates["m"].done);
        assert!(state.candidates["z"].done);
    }

    #[test]
    fn enumeration_pruning_rekeys_only_when_stereo_changes_and_tracks_completion() {
        let mut state = expansion_state(&[("old", "old", true), ("stable", "stable", false)]);
        state.modified_atoms = marked_atoms([0]);
        let mut stereo_calls = Vec::new();
        let mut key_calls = 0;

        let pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            TautomerOptions::default(),
            false,
            |tautomer, modified_atoms, modified_bonds| {
                stereo_calls.push(tautomer.clone());
                assert_eq!(modified_atoms, &marked_atoms([0]));
                assert!(modified_bonds.is_empty());
                if tautomer == "old-tautomer" {
                    *tautomer = "new".to_owned();
                    Ok::<_, &'static str>(true)
                } else {
                    Ok(false)
                }
            },
            |tautomer| {
                key_calls += 1;
                Ok::<_, &'static str>(tautomer.clone())
            },
        )
        .expect("prune changed and unchanged stereo branches");

        assert!(!pass.completed);
        assert!(!pass.bailed_out);
        assert_eq!(stereo_calls, ["old-tautomer", "stable-tautomer"]);
        assert_eq!(key_calls, 1);
        assert!(!state.candidates.contains_key("old"));
        assert_eq!(state.candidates["new"].num_modified_atoms, 1);
        assert_eq!(state.candidates["new"].num_modified_bonds, 0);
        assert_eq!(state.candidates["stable"].num_modified_atoms, 0);
        assert_eq!(state.candidates["stable"].num_modified_bonds, 0);
    }

    #[test]
    fn enumeration_pruning_reproduces_source_rekey_iterator_order() {
        let mut forward = expansion_state(&[("a", "a", true), ("b", "b", true), ("c", "c", true)]);
        forward.modified_atoms = marked_atoms([0]);
        let mut forward_calls = Vec::new();
        let forward_pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut forward,
            TautomerOptions::default(),
            false,
            |tautomer, _, _| {
                forward_calls.push(tautomer.clone());
                if tautomer == "a-tautomer" {
                    *tautomer = "z".to_owned();
                    Ok::<_, &'static str>(true)
                } else {
                    Ok(false)
                }
            },
            |tautomer| Ok::<_, &'static str>(tautomer.clone()),
        )
        .expect("forward rekey traversal");
        assert!(forward_pass.completed);
        assert_eq!(forward_calls, ["a-tautomer"]);
        assert_eq!(
            forward
                .candidates
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["b", "c", "z"]
        );

        let mut backward = expansion_state(&[("a", "a", true), ("b", "b", true), ("c", "c", true)]);
        backward.modified_atoms = marked_atoms([0]);
        let mut backward_calls = Vec::new();
        let backward_pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut backward,
            TautomerOptions::default(),
            false,
            |tautomer, _, _| {
                backward_calls.push(tautomer.clone());
                if tautomer == "c-tautomer" {
                    *tautomer = "aa".to_owned();
                    Ok::<_, &'static str>(true)
                } else {
                    Ok(false)
                }
            },
            |tautomer| Ok::<_, &'static str>(tautomer.clone()),
        )
        .expect("backward rekey traversal");
        assert!(backward_pass.completed);
        assert_eq!(
            backward_calls,
            ["a-tautomer", "b-tautomer", "c-tautomer", "b-tautomer"]
        );
        assert_eq!(
            backward
                .candidates
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["a", "aa", "b"]
        );
    }

    #[test]
    fn enumeration_pruning_collapses_duplicates_and_corrects_only_tautomer_limit_status() {
        let mut state = expansion_state(&[("a", "a", true), ("b", "b", true)]);
        state.modified_atoms = marked_atoms([0]);
        state.candidates.get_mut("a").unwrap().num_modified_atoms = 1;
        state.status = TautomerEnumerationStatus::MaxTautomersReached;

        let pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            TautomerOptions::default().with_max_tautomers(2),
            true,
            |tautomer, _, _| {
                assert_eq!(tautomer, "b-tautomer");
                *tautomer = "a".to_owned();
                Ok::<_, &'static str>(true)
            },
            |tautomer| Ok::<_, &'static str>(tautomer.clone()),
        )
        .expect("collapse duplicate and correct status");

        assert!(pass.completed);
        assert!(!pass.bailed_out);
        assert_eq!(state.status, TautomerEnumerationStatus::Completed);
        assert_eq!(state.candidates.len(), 1);
        assert_eq!(
            state.candidates["a"].tautomer.as_deref(),
            Some("a-tautomer")
        );

        let mut transform_limited = expansion_state(&[("a", "a", true), ("b", "b", true)]);
        transform_limited.modified_atoms = marked_atoms([0]);
        transform_limited
            .candidates
            .get_mut("a")
            .unwrap()
            .num_modified_atoms = 1;
        transform_limited.status = TautomerEnumerationStatus::MaxTransformsReached;
        let pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut transform_limited,
            TautomerOptions::default().with_max_tautomers(2),
            true,
            |tautomer, _, _| {
                *tautomer = "a".to_owned();
                Ok::<_, &'static str>(true)
            },
            |tautomer| Ok::<_, &'static str>(tautomer.clone()),
        )
        .expect("retain non-tautomer limit status");
        assert!(pass.bailed_out);
        assert_eq!(
            transform_limited.status,
            TautomerEnumerationStatus::MaxTransformsReached
        );
    }

    #[test]
    fn enumeration_pruning_reapplies_after_modified_sets_grow_across_rounds() {
        let mut state = expansion_state(&[("m", "m", true)]);
        state.modified_atoms = marked_atoms([0]);
        state.modified_bonds = marked_bonds([0]);
        {
            let candidate = state.candidates.get_mut("m").unwrap();
            candidate.num_modified_atoms = 1;
            candidate.num_modified_bonds = 1;
        }

        let first = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            TautomerOptions::default(),
            false,
            |_, _, _| -> Result<_, &'static str> {
                unreachable!("equal modified-set snapshots do not reapply stereo")
            },
            |_| -> Result<_, &'static str> {
                unreachable!("an unchanged candidate is not rekeyed")
            },
        )
        .expect("prune with unchanged modified sets");
        assert!(first.completed);

        state.modified_atoms.insert(AtomId::new(1));
        state.modified_bonds.insert(BondId::new(1));
        let mut calls = 0;
        let second = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            TautomerOptions::default(),
            false,
            |tautomer, modified_atoms, modified_bonds| {
                calls += 1;
                assert_eq!(modified_atoms, &marked_atoms([0, 1]));
                assert_eq!(modified_bonds, &marked_bonds([0, 1]));
                tautomer.push_str("-updated");
                Ok::<_, &'static str>(true)
            },
            |_| Ok::<_, &'static str>("m".to_owned()),
        )
        .expect("prune after modified sets grow");

        assert!(second.completed);
        assert_eq!(calls, 1);
        assert_eq!(state.candidates["m"].num_modified_atoms, 2);
        assert_eq!(state.candidates["m"].num_modified_bonds, 2);
        assert_eq!(
            state.candidates["m"].tautomer.as_deref(),
            Some("m-tautomer-updated")
        );
    }

    #[test]
    fn enumeration_pruning_materializes_final_candidates_in_key_order() {
        let state = expansion_state(&[("z", "z", true), ("a", "a", true), ("m", "m", true)]);
        let entries = materialize_tautomer_candidates_in_source_order(state.candidates)
            .expect("materialize ordered candidates");
        assert_eq!(
            entries,
            [
                ("a".to_owned(), "a-tautomer".to_owned()),
                ("m".to_owned(), "m-tautomer".to_owned()),
                ("z".to_owned(), "z-tautomer".to_owned()),
            ]
        );

        let mut missing = expansion_state(&[("missing", "missing", true)]).candidates;
        missing.get_mut("missing").unwrap().tautomer = None;
        assert!(matches!(
            materialize_tautomer_candidates_in_source_order(missing),
            Err(TautomerEnumerationError::MissingCandidateMolecule { canonical_smiles })
                if canonical_smiles == "missing"
        ));
    }

    #[test]
    fn candidate_record_default_constructor_has_source_empty_state() {
        let candidate = TautomerCandidate::empty();

        assert!(candidate.tautomer.is_none());
        assert!(candidate.kekulized.is_none());
        assert_eq!(candidate.num_modified_atoms, 0);
        assert_eq!(candidate.num_modified_bonds, 0);
        assert!(!candidate.done);
    }

    #[test]
    fn candidate_record_explicit_constructor_preserves_values_and_initial_state() {
        let tautomer = Molecule::from_smiles("CCO").expect("parse tautomer");
        let kekulized = Molecule::from_smiles("CCO").expect("parse kekulized tautomer");
        let candidate = TautomerCandidate::new(tautomer.clone(), kekulized.clone(), 3, 2);

        assert_eq!(candidate.tautomer.as_ref(), Some(&tautomer));
        assert_eq!(candidate.kekulized.as_ref(), Some(&kekulized));
        assert_eq!(candidate.num_modified_atoms, 3);
        assert_eq!(candidate.num_modified_bonds, 2);
        assert!(!candidate.done);
    }

    #[test]
    fn candidate_record_clone_and_state_transitions_are_independent() {
        let molecule = Molecule::from_smiles("CC=O").expect("parse candidate");
        let original = TautomerCandidate::new(molecule.clone(), molecule, 1, 1);
        let mut changed = original.clone();

        changed.update_modified_counts(4, 5);
        changed.mark_done();

        assert_eq!(original.num_modified_atoms, 1);
        assert_eq!(original.num_modified_bonds, 1);
        assert!(!original.done);
        assert_eq!(changed.num_modified_atoms, 4);
        assert_eq!(changed.num_modified_bonds, 5);
        assert!(changed.done);
    }

    #[test]
    fn candidate_record_ordered_map_sorts_keys_and_replaces_one_canonical_key() {
        let molecule = Molecule::from_smiles("C").expect("parse candidate");
        let candidate = |modified_atoms| {
            TautomerCandidate::new(molecule.clone(), molecule.clone(), modified_atoms, 0)
        };
        let mut candidates = SmilesTautomerMap::new();
        candidates.insert("z".to_owned(), candidate(1));
        candidates.insert("a".to_owned(), candidate(2));
        candidates.insert("m".to_owned(), candidate(3));
        let replaced = candidates.insert("a".to_owned(), candidate(9));

        assert_eq!(replaced.expect("existing key").num_modified_atoms, 2);
        assert_eq!(
            candidates.keys().map(String::as_str).collect::<Vec<_>>(),
            ["a", "m", "z"]
        );
        assert_eq!(candidates["a"].num_modified_atoms, 9);
    }

    fn enumeration_result_fixture(status: TautomerEnumerationStatus) -> TautomerEnumeration {
        let mut candidates = SmilesTautomerMap::new();
        for (smiles, molecule_smiles) in [("z", "CCC"), ("a", "C"), ("m", "CC")] {
            let molecule = Molecule::from_smiles(molecule_smiles).expect("parse result molecule");
            candidates.insert(
                smiles.to_owned(),
                TautomerCandidate::new(molecule.clone(), molecule, 1, 1),
            );
        }
        TautomerEnumeration::from_candidates(
            candidates,
            status,
            BTreeSet::from([AtomId::new(0), AtomId::new(2)]),
            BTreeSet::from([BondId::new(1)]),
        )
        .expect("assemble result")
    }

    #[test]
    fn enumeration_result_default_is_empty_and_completed() {
        let result = TautomerEnumeration::default();

        assert!(result.is_empty());
        assert_eq!(result.len(), 0);
        assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
        assert!(result.modified_atoms().is_empty());
        assert!(result.modified_bonds().is_empty());
        assert!(result.iter().next().is_none());
    }

    #[test]
    fn enumeration_result_single_entry_supports_all_access_paths() {
        let molecule = Molecule::from_smiles("CCO").expect("parse sole tautomer");
        let mut candidates = SmilesTautomerMap::new();
        candidates.insert(
            "CCO".to_owned(),
            TautomerCandidate::new(molecule.clone(), molecule, 0, 0),
        );
        let result = TautomerEnumeration::from_candidates(
            candidates,
            TautomerEnumerationStatus::Completed,
            BTreeSet::new(),
            BTreeSet::new(),
        )
        .expect("assemble single result");

        assert_eq!(result.len(), 1);
        assert!(!result.is_empty());
        assert_eq!(result.canonical_smiles(), ["CCO"]);
        assert_eq!(result.at(0).unwrap().num_atoms(), 3);
    }

    #[test]
    fn enumeration_result_preserves_canonical_smiles_order_in_every_projection() {
        let result = enumeration_result_fixture(TautomerEnumerationStatus::Completed);

        assert_eq!(result.canonical_smiles(), ["a", "m", "z"]);
        assert_eq!(
            result
                .iter_with_smiles()
                .map(|(smiles, molecule)| (smiles, molecule.num_atoms()))
                .collect::<Vec<_>>(),
            [("a", 1), ("m", 2), ("z", 3)]
        );
        assert_eq!(
            result.iter().map(Molecule::num_atoms).collect::<Vec<_>>(),
            [1, 2, 3]
        );
        assert_eq!(
            result
                .iter()
                .rev()
                .map(Molecule::num_atoms)
                .collect::<Vec<_>>(),
            [3, 2, 1]
        );
        assert_eq!(
            (&result)
                .into_iter()
                .map(Molecule::num_atoms)
                .collect::<Vec<_>>(),
            [1, 2, 3]
        );
    }

    #[test]
    fn enumeration_result_random_access_and_projections_cover_bounds() {
        let result = enumeration_result_fixture(TautomerEnumerationStatus::Completed);

        assert_eq!(result.at(0).unwrap().num_atoms(), 1);
        assert_eq!(result[1].num_atoms(), 2);
        assert_eq!(result.get(2).map(Molecule::num_atoms), Some(3));
        assert!(result.get(3).is_none());
        assert_eq!(
            result.at(3),
            Err(TautomerEnumerationError::IndexOutOfRange { index: 3, len: 3 })
        );
        assert_eq!(
            result
                .molecules()
                .iter()
                .map(Molecule::num_atoms)
                .collect::<Vec<_>>(),
            [1, 2, 3]
        );
    }

    #[test]
    fn enumeration_result_clone_has_independent_entries_and_typed_modified_sets() {
        let original = enumeration_result_fixture(TautomerEnumerationStatus::Completed);
        let mut cloned = original.clone();
        cloned.entries.remove(0);
        cloned.modified_atoms.insert(AtomId::new(9));
        cloned.modified_bonds.clear();

        assert_eq!(original.len(), 3);
        assert_eq!(cloned.len(), 2);
        assert_eq!(
            original.modified_atoms(),
            &BTreeSet::from([AtomId::new(0), AtomId::new(2)])
        );
        assert_eq!(original.modified_bonds(), &BTreeSet::from([BondId::new(1)]));
        assert!(cloned.modified_atoms().contains(&AtomId::new(9)));
        assert!(cloned.modified_bonds().is_empty());
    }

    #[test]
    fn enumeration_result_exposes_every_source_status_without_reinterpretation() {
        for status in [
            TautomerEnumerationStatus::Completed,
            TautomerEnumerationStatus::MaxTautomersReached,
            TautomerEnumerationStatus::MaxTransformsReached,
            TautomerEnumerationStatus::Canceled,
        ] {
            assert_eq!(enumeration_result_fixture(status).status(), status);
        }
    }

    #[test]
    fn deprecated_enumerate_correspondence_matches_every_rich_result_molecule_and_modified_set() {
        let cases = [
            ("C", TautomerOptions::default()),
            ("CC(C)=O", TautomerOptions::default()),
            ("CC(C)=O", TautomerOptions::default().with_max_transforms(1)),
            (
                "OC(C)=C(C)C",
                TautomerOptions::default().with_max_tautomers(2),
            ),
        ];

        for (smiles, options) in cases {
            let molecule = Molecule::from_smiles(smiles).expect("parse correspondence fixture");
            let before = molecule.clone();
            let enumerator = TautomerEnumerator::from_options(options);
            let rich = enumerator.enumerate(&molecule).expect("rich enumeration");
            let mut modified_atoms = BTreeSet::from([AtomId::new(molecule.num_atoms() + 7)]);
            let mut modified_bonds = BTreeSet::from([BondId::new(molecule.num_bonds() + 7)]);

            let molecules = enumerator
                .enumerate_deprecated_compat(
                    &molecule,
                    Some(&mut modified_atoms),
                    Some(&mut modified_bonds),
                )
                .expect("deprecated compatibility enumeration");

            assert_eq!(molecules, rich.molecules(), "molecules for {smiles}");
            assert_eq!(modified_atoms, *rich.modified_atoms(), "atoms for {smiles}");
            assert_eq!(modified_bonds, *rich.modified_bonds(), "bonds for {smiles}");
            assert_eq!(molecule, before, "source for {smiles}");
        }
    }

    #[test]
    fn deprecated_enumerate_correspondence_honors_optional_outputs_and_failure_order() {
        let molecule = Molecule::from_smiles("CC(C)=O").expect("parse optional-output fixture");
        let enumerator = TautomerEnumerator::new();
        let rich = enumerator.enumerate(&molecule).expect("rich enumeration");

        let mut atoms_only = BTreeSet::from([AtomId::new(99)]);
        assert_eq!(
            enumerator
                .enumerate_deprecated_compat(&molecule, Some(&mut atoms_only), None)
                .unwrap(),
            rich.molecules()
        );
        assert_eq!(atoms_only, *rich.modified_atoms());

        let mut bonds_only = BTreeSet::from([BondId::new(99)]);
        assert_eq!(
            enumerator
                .enumerate_deprecated_compat(&molecule, None, Some(&mut bonds_only))
                .unwrap(),
            rich.molecules()
        );
        assert_eq!(bonds_only, *rich.modified_bonds());
        assert_eq!(
            enumerator
                .enumerate_deprecated_compat(&molecule, None, None)
                .unwrap(),
            rich.molecules()
        );

        let invalid = Molecule::from_smiles_with_sanitize("c1cccc1", false)
            .expect("parse invalid aromatic fixture");
        let mut unchanged_atoms = BTreeSet::from([AtomId::new(17)]);
        let mut unchanged_bonds = BTreeSet::from([BondId::new(19)]);
        assert!(
            enumerator
                .enumerate_deprecated_compat(
                    &invalid,
                    Some(&mut unchanged_atoms),
                    Some(&mut unchanged_bonds),
                )
                .is_err()
        );
        assert_eq!(unchanged_atoms, BTreeSet::from([AtomId::new(17)]));
        assert_eq!(unchanged_bonds, BTreeSet::from([BondId::new(19)]));
    }

    fn canonical_selection_smiles(molecule: &Molecule) -> String {
        MoleculeReadParts::from_molecule(molecule)
            .canonical_isomeric_smiles()
            .expect("write canonical-selection molecule")
    }

    #[test]
    fn enumeration_clears_computed_ring_stereo_before_transforming_candidates_like_rdkit() {
        let molecule = Molecule::from_smiles(
            "N[C@H](C(=O)N1CCCC1)[C@H]1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1",
        )
        .expect("parse CHEMBL23979 ring-stereo regression");
        let enumerator = TautomerEnumerator::from_options(
            TautomerOptions::default().with_reassign_stereo(false),
        );

        let result = enumerator
            .enumerate(&molecule)
            .expect("enumerate CHEMBL23979 ring-stereo regression");

        assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
        assert_eq!(
            result.canonical_smiles(),
            [
                "N=C(C1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1)C(O)N1CCCC1",
                "NC(=C(O)N1CCCC1)C1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1",
                "NC(=C1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1)C(O)N1CCCC1",
                "NC(C(=O)N1CCCC1)C1CC[C@H](NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1",
            ]
        );
        assert_eq!(
            result.modified_atoms(),
            &BTreeSet::from([
                AtomId::new(0),
                AtomId::new(1),
                AtomId::new(2),
                AtomId::new(3),
                AtomId::new(9),
            ])
        );
        assert_eq!(
            result.modified_bonds(),
            &BTreeSet::from([
                BondId::new(0),
                BondId::new(1),
                BondId::new(2),
                BondId::new(8),
            ])
        );
        assert_eq!(
            result
                .iter()
                .map(|candidate| score_tautomer(candidate).unwrap().total())
                .collect::<Vec<_>>(),
            [251, 250, 250, 255]
        );
        assert_eq!(
            enumerator
                .pick_canonical(&result)
                .expect("select CHEMBL23979 canonical tautomer")
                .to_smiles(true)
                .expect("write CHEMBL23979 canonical tautomer"),
            "NC(C(=O)N1CCCC1)C1CCC(NS(=O)(=O)c2ccc(OC(F)(F)F)cc2)CC1"
        );
    }

    #[test]
    fn enumeration_rekeys_stale_computed_double_bond_stereo_like_rdkit_chembl12724() {
        let molecule =
            Molecule::from_smiles("COc1ccc(OC)c(/C=N/N=C(\\N)NO)c1.Cc1ccc(S(=O)(=O)O)cc1").unwrap();
        let enumerator = TautomerEnumerator::from_options(
            TautomerOptions::default().with_reassign_stereo(false),
        );
        let result = enumerator.enumerate(&molecule).unwrap();
        assert_eq!(
            result.canonical_smiles(),
            &[
                "COc1ccc(OC)c(/C=N/N=C(N)NO)c1.Cc1ccc(S(=O)(=O)O)cc1",
                "COc1ccc(OC)c(/C=N/NC(=N)NO)c1.Cc1ccc(S(=O)(=O)O)cc1",
                "COc1ccc(OC)c(/C=N/NC(N)=NO)c1.Cc1ccc(S(=O)(=O)O)cc1",
                "COc1ccc(OC)c(/C=N/NC(N)N=O)c1.Cc1ccc(S(=O)(=O)O)cc1",
                "COc1ccc(OC)c(C=NN=C(N)NO)c1.Cc1ccc(S(=O)(=O)O)cc1",
            ]
        );
        assert_eq!(
            result.modified_atoms(),
            &BTreeSet::from([
                AtomId::new(11),
                AtomId::new(12),
                AtomId::new(13),
                AtomId::new(14),
                AtomId::new(15),
            ])
        );
        assert_eq!(
            result.modified_bonds(),
            &BTreeSet::from([
                BondId::new(11),
                BondId::new(12),
                BondId::new(13),
                BondId::new(14),
            ])
        );
        assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
        assert_eq!(
            result
                .iter()
                .map(|candidate| score_tautomer(candidate).unwrap().total())
                .collect::<Vec<_>>(),
            [509, 509, 513, 506, 509]
        );
        assert_eq!(
            enumerator
                .pick_canonical(&result)
                .expect("select CHEMBL12724 canonical tautomer")
                .to_smiles(true)
                .expect("write CHEMBL12724 canonical tautomer"),
            "COc1ccc(OC)c(C=NNC(N)=NO)c1.Cc1ccc(S(=O)(=O)O)cc1"
        );
    }

    #[test]
    fn canonical_selection_rejects_empty_inputs_and_unselectable_minimum_scores() {
        let enumerator = TautomerEnumerator::new();
        let empty = TautomerEnumeration::default();
        assert!(matches!(
            enumerator.pick_canonical(&empty),
            Err(TautomerRunError::NoCanonicalTautomer)
        ));

        let no_molecules: Vec<Molecule> = Vec::new();
        assert!(matches!(
            enumerator.pick_canonical_from_iterable(no_molecules.iter()),
            Err(TautomerRunError::NoCanonicalTautomer)
        ));

        let methane = Molecule::from_smiles("C").expect("parse methane");
        let ethane = Molecule::from_smiles("CC").expect("parse ethane");
        assert!(matches!(
            enumerator.pick_canonical_from_iterable_with([&methane, &ethane], |_| Ok(i32::MIN)),
            Err(TautomerRunError::NoCanonicalTautomer)
        ));
    }

    #[test]
    fn canonical_selection_single_item_skips_scoring_and_returns_a_clean_copy() {
        let enumerator = TautomerEnumerator::new();
        let mut source = Molecule::from_smiles("F[C@H](Cl)Br").expect("parse chiral molecule");
        source.properties_mut().clear_prop("_StereochemDone");
        let before = source.clone();
        let result = TautomerEnumeration::from_ordered_entries(
            vec![("retained-without-rewrite".to_owned(), source.clone())],
            TautomerEnumerationStatus::Completed,
            BTreeSet::new(),
            BTreeSet::new(),
        );

        let mut result_score_calls = 0;
        let mut selected = enumerator
            .pick_canonical_with(&result, |_| {
                result_score_calls += 1;
                Ok(7)
            })
            .expect("select sole result tautomer");
        assert_eq!(result_score_calls, 0);
        assert_eq!(source, before);
        assert_eq!(source.prop("_StereochemDone"), None);
        assert_eq!(selected.prop("_StereochemDone"), Some("1"));
        assert!(selected.is_prop_computed("_StereochemDone"));

        selected.properties_mut().set_prop("selected-only", "yes");
        assert_eq!(source.prop("selected-only"), None);

        let mut iterable_score_calls = 0;
        let iterable_selected = enumerator
            .pick_canonical_from_iterable_with([&source], |_| {
                iterable_score_calls += 1;
                Ok(-11)
            })
            .expect("select sole iterable tautomer");
        assert_eq!(iterable_score_calls, 0);
        assert_eq!(source, before);
        assert_eq!(iterable_selected.prop("_StereochemDone"), Some("1"));
    }

    #[test]
    fn canonical_selection_uses_signed_maximum_and_retained_result_keys() {
        let enumerator = TautomerEnumerator::new();
        let methane = Molecule::from_smiles("C").expect("parse methane");
        let ethane = Molecule::from_smiles("CC").expect("parse ethane");
        let propane = Molecule::from_smiles("CCC").expect("parse propane");
        let result = TautomerEnumeration::from_ordered_entries(
            vec![
                ("z-retained".to_owned(), methane.clone()),
                ("a-retained".to_owned(), ethane.clone()),
                ("m-retained".to_owned(), propane.clone()),
            ],
            TautomerEnumerationStatus::Completed,
            BTreeSet::new(),
            BTreeSet::new(),
        );

        let negative = enumerator
            .pick_canonical_with(&result, |molecule| {
                Ok(match molecule.num_atoms() {
                    1 => -10,
                    2 => -1,
                    _ => -5,
                })
            })
            .expect("select greatest negative score");
        assert_eq!(canonical_selection_smiles(&negative), "CC");

        let positive = enumerator
            .pick_canonical_with(&result, |molecule| Ok(molecule.num_atoms() as i32))
            .expect("select greatest positive score");
        assert_eq!(canonical_selection_smiles(&positive), "CCC");

        let tied = enumerator
            .pick_canonical_with(&result, |_| Ok(4))
            .expect("break retained-key tie");
        assert_eq!(canonical_selection_smiles(&tied), "CC");
    }

    #[test]
    fn canonical_selection_iterable_computes_lexical_ties_and_matches_result_path() {
        let enumerator = TautomerEnumerator::new();
        let propane = Molecule::from_smiles("CCC").expect("parse propane");
        let ethane = Molecule::from_smiles("CC").expect("parse ethane");
        let methane = Molecule::from_smiles("C").expect("parse methane");
        let inputs = [&propane, &ethane, &methane];

        let iterable = enumerator
            .pick_canonical_from_iterable_with(inputs, |_| Ok(9))
            .expect("break computed-SMILES tie");
        assert_eq!(canonical_selection_smiles(&iterable), "C");

        let mut entries = inputs
            .into_iter()
            .map(|molecule| (canonical_selection_smiles(molecule), molecule.clone()))
            .collect::<Vec<_>>();
        entries.sort_by(|left, right| left.0.cmp(&right.0));
        let result = TautomerEnumeration::from_ordered_entries(
            entries,
            TautomerEnumerationStatus::Completed,
            BTreeSet::new(),
            BTreeSet::new(),
        );
        let retained = enumerator
            .pick_canonical_with(&result, |_| Ok(9))
            .expect("break retained-SMILES tie");
        assert_eq!(retained, iterable);
        assert_eq!(propane, Molecule::from_smiles("CCC").unwrap());
        assert_eq!(ethane, Molecule::from_smiles("CC").unwrap());
        assert_eq!(methane, Molecule::from_smiles("C").unwrap());
    }

    #[test]
    fn canonical_selection_default_and_custom_scorers_share_finalization() {
        let molecule = Molecule::from_smiles("CC(C)=O").expect("parse tautomerizable molecule");
        let before = molecule.clone();
        let enumerator = TautomerEnumerator::new();
        let result = enumerator
            .enumerate(&molecule)
            .expect("enumerate tautomers");

        let default_selected = enumerator
            .pick_canonical(&result)
            .expect("select with default score");
        let custom_selected = enumerator
            .pick_canonical_from_iterable_with(result.iter(), |candidate| {
                Ok(score_tautomer(candidate)?.total())
            })
            .expect("select with custom default-equivalent score");

        assert_eq!(default_selected, custom_selected);
        assert_eq!(default_selected.prop("_StereochemDone"), Some("1"));
        assert!(default_selected.is_prop_computed("_StereochemDone"));
        assert_eq!(molecule, before);
    }

    #[test]
    fn canonicalization_and_factories_share_current_and_v1_catalog_paths() {
        let molecule = Molecule::from_smiles("CC(C)=O").expect("parse factory fixture");
        let default = TautomerEnumerator::new();
        let current = TautomerEnumerator::from_catalog(
            TautomerCatalog::current().expect("construct current catalog"),
        );
        let v1 = TautomerEnumerator::v1().expect("construct V1 enumerator");

        assert_eq!(default.catalog().transforms().len(), 37);
        assert_eq!(current.catalog().transforms().len(), 37);
        assert_eq!(v1.catalog().transforms().len(), 36);
        assert_eq!(default.options(), TautomerOptions::default());
        assert_eq!(v1.options(), TautomerOptions::default());

        let default_endpoint = default.canonicalize(&molecule).unwrap();
        let current_endpoint = current.canonicalize(&molecule).unwrap();
        let v1_endpoint = v1.canonicalize(&molecule).unwrap();
        assert_eq!(default_endpoint, current_endpoint);
        assert_eq!(
            canonical_selection_smiles(&default_endpoint),
            canonical_selection_smiles(&v1_endpoint)
        );
    }

    #[test]
    fn canonicalization_and_factories_honor_options_without_mutating_configuration_or_source() {
        let molecule = Molecule::from_smiles("CC(C)=O")
            .expect("parse option fixture")
            .with_prop("source_id", "canonical-options");
        let before = molecule.clone();
        let options = TautomerOptions::default()
            .with_max_transforms(0)
            .with_remove_sp3_stereo(false)
            .with_remove_bond_stereo(false)
            .with_remove_isotopic_hydrogens(false)
            .with_reassign_stereo(true);
        let enumerator = TautomerEnumerator::from_options(options);

        let endpoint = enumerator
            .canonicalize(&molecule)
            .expect("canonicalize transform-limited input");

        assert_eq!(canonical_selection_smiles(&endpoint), "CC(C)=O");
        assert_eq!(endpoint.prop("source_id"), Some("canonical-options"));
        assert_eq!(enumerator.options(), options);
        assert!(enumerator.reassign_stereo());
        assert_eq!(molecule, before);
    }

    #[test]
    fn canonicalization_and_factories_custom_scorer_selects_from_the_sole_enumeration() {
        let molecule = Molecule::from_smiles("CC(C)=O").expect("parse custom-score fixture");
        let before = molecule.clone();
        let enumerator = TautomerEnumerator::new();
        let mut scored = Vec::new();

        let endpoint = enumerator
            .canonicalize_with(&molecule, |candidate| {
                let smiles = canonical_selection_smiles(candidate);
                scored.push(smiles.clone());
                Ok(if smiles == "C=C(C)O" { 100 } else { -100 })
            })
            .expect("canonicalize with custom score");

        scored.sort();
        assert_eq!(scored, ["C=C(C)O", "CC(C)=O"]);
        assert_eq!(canonical_selection_smiles(&endpoint), "C=C(C)O");
        assert_eq!(endpoint.prop("_StereochemDone"), Some("1"));
        assert_eq!(molecule, before);
    }

    #[test]
    fn canonicalization_and_factories_are_independent_of_input_tautomer_and_atom_order() {
        let enumerator = TautomerEnumerator::new();
        let inputs = ["CC(C)=O", "C=C(C)O", "O=C(C)C"];
        let endpoints = inputs
            .into_iter()
            .map(|smiles| {
                let molecule = Molecule::from_smiles(smiles).expect("parse endpoint fixture");
                canonical_selection_smiles(
                    &enumerator
                        .canonicalize(&molecule)
                        .expect("canonicalize endpoint fixture"),
                )
            })
            .collect::<Vec<_>>();

        assert!(endpoints.iter().all(|endpoint| endpoint == &endpoints[0]));
        assert_eq!(endpoints[0], "CC(C)=O");
    }

    #[test]
    fn canonicalization_and_factories_private_in_place_correspondence_preserves_outer_state() {
        let source = Molecule::from_smiles("CC(C)=O")
            .expect("parse in-place fixture")
            .with_2d_coordinates()
            .expect("generate 2D coordinates");
        let coordinates_3d = (0..source.num_atoms())
            .map(|index| [index as f64, index as f64 + 0.5, -(index as f64)])
            .collect::<Vec<_>>();
        let source = source
            .with_only_3d_conformer(coordinates_3d, true)
            .expect("attach 3D coordinates")
            .with_name("in-place-source")
            .with_prop("source_id", "in-place-correspondence")
            .with_sdf_data_field("dataset", "tautomer");
        let original_coordinates_2d = source.conformers_2d().to_vec();
        let original_coordinates_3d = source.conformers_3d().to_vec();
        let original_properties = source.properties().clone();
        let enumerator = TautomerEnumerator::new();
        let expected = enumerator.canonicalize(&source).unwrap();
        let mut in_place = source.clone();

        enumerator
            .canonicalize_in_place_compat_with(&mut in_place, |candidate| {
                Ok(score_tautomer(candidate)?.total())
            })
            .expect("run private in-place correspondence");

        for (actual, expected) in in_place.atoms().iter().zip(expected.atoms()) {
            assert_eq!(actual.atomic_number(), expected.atomic_number());
            assert_eq!(actual.formal_charge(), expected.formal_charge());
            assert_eq!(actual.no_implicit(), expected.no_implicit());
            assert_eq!(actual.is_aromatic(), expected.is_aromatic());
            assert_eq!(actual.explicit_hydrogens(), expected.explicit_hydrogens());
            assert_eq!(actual.radical_electrons(), expected.radical_electrons());
            assert_eq!(actual.chiral_tag(), expected.chiral_tag());
        }
        for (actual, expected) in in_place.bonds().iter().zip(expected.bonds()) {
            assert_eq!(actual.begin(), expected.begin());
            assert_eq!(actual.end(), expected.end());
            assert_eq!(actual.order(), expected.order());
            assert_eq!(actual.direction(), expected.direction());
            assert_eq!(actual.is_aromatic(), expected.is_aromatic());
            assert_eq!(actual.is_conjugated(), expected.is_conjugated());
            assert_eq!(actual.stereo(), expected.stereo());
            if expected.stereo_atoms().is_some() {
                assert_eq!(actual.stereo_atoms(), expected.stereo_atoms());
            }
        }
        assert_eq!(in_place.conformers_2d(), original_coordinates_2d);
        assert_eq!(in_place.conformers_3d(), original_coordinates_3d);
        assert_eq!(in_place.properties(), &original_properties);
        assert!(crate::valence::cached_valence_assignment(&in_place).is_some());
    }

    #[test]
    fn canonicalization_and_factories_do_not_expose_a_duplicate_in_place_api() {
        let source = include_str!("tautomer.rs");
        let forbidden = ["pub fn canonicalize", "_in_place"].concat();
        assert!(!source.contains(&forbidden));
        assert!(source.contains("fn canonicalize_in_place_compat_with"));
        assert!(source.contains("#[cfg(test)]"));
    }

    #[test]
    fn enumeration_result_rejects_an_unmaterialized_candidate() {
        let mut candidates = SmilesTautomerMap::new();
        candidates.insert("C".to_owned(), TautomerCandidate::empty());

        assert_eq!(
            TautomerEnumeration::from_candidates(
                candidates,
                TautomerEnumerationStatus::Completed,
                BTreeSet::new(),
                BTreeSet::new(),
            ),
            Err(TautomerEnumerationError::MissingCandidateMolecule {
                canonical_smiles: "C".to_owned(),
            })
        );
    }

    #[test]
    fn enumerator_configuration_defaults_match_every_source_default() {
        let options = TautomerOptions::default();
        assert_eq!(options.max_tautomers(), 1000);
        assert_eq!(options.max_transforms(), 1000);
        assert!(options.remove_sp3_stereo());
        assert!(options.remove_bond_stereo());
        assert!(options.remove_isotopic_hydrogens());
        assert!(options.reassign_stereo());

        let enumerator = TautomerEnumerator::new();
        assert_eq!(enumerator.options(), options);
        assert_eq!(enumerator.catalog().transforms().len(), 37);
        assert!(enumerator.callback().is_none());
    }

    #[test]
    fn enumeration_matches_pcs_fused_ring_max_transform_boundary() {
        let molecule = Molecule::from_smiles("Cc1nc2c(nc1C)C(=O)C1=C(C2=O)C2C=CC1CC2")
            .expect("parse PCS regression molecule");

        let result = TautomerEnumerator::default()
            .enumerate(&molecule)
            .expect("enumerate PCS regression molecule");

        assert_eq!(result.len(), 272);
        assert_eq!(
            result.status(),
            TautomerEnumerationStatus::MaxTransformsReached
        );
    }

    #[test]
    fn enumerator_configuration_setters_cover_zero_maximum_and_every_boolean() {
        let mut enumerator = TautomerEnumerator::new();
        enumerator.set_max_tautomers(0);
        enumerator.set_max_transforms(u32::MAX);
        enumerator.set_remove_sp3_stereo(false);
        enumerator.set_remove_bond_stereo(false);
        enumerator.set_remove_isotopic_hydrogens(false);
        enumerator.set_reassign_stereo(false);

        assert_eq!(enumerator.max_tautomers(), 0);
        assert_eq!(enumerator.max_transforms(), u32::MAX);
        assert!(!enumerator.remove_sp3_stereo());
        assert!(!enumerator.remove_bond_stereo());
        assert!(!enumerator.remove_isotopic_hydrogens());
        assert!(!enumerator.reassign_stereo());

        let options = TautomerOptions::default()
            .with_max_tautomers(u32::MAX)
            .with_max_transforms(0)
            .with_remove_sp3_stereo(false)
            .with_remove_bond_stereo(false)
            .with_remove_isotopic_hydrogens(false)
            .with_reassign_stereo(false);
        assert_eq!(options.max_tautomers(), u32::MAX);
        assert_eq!(options.max_transforms(), 0);
        assert!(!options.remove_sp3_stereo());
        assert!(!options.remove_bond_stereo());
        assert!(!options.remove_isotopic_hydrogens());
        assert!(!options.reassign_stereo());
    }

    #[test]
    fn enumerator_configuration_clone_shares_catalog_but_not_option_state() {
        let definitions = [("custom", "[O]-[C]", "=", "")];
        let catalog = TautomerCatalog::from_data(&definitions).expect("compile custom catalog");
        let options = TautomerOptions::default().with_max_tautomers(17);
        let original = TautomerEnumerator::from_catalog_and_options(catalog.clone(), options);
        let mut cloned = original.clone();

        assert!(Arc::ptr_eq(&original.catalog, &cloned.catalog));
        assert_eq!(original.catalog(), &catalog);
        cloned.set_max_tautomers(2);
        cloned.set_remove_sp3_stereo(false);
        assert_eq!(original.max_tautomers(), 17);
        assert!(original.remove_sp3_stereo());
        assert_eq!(cloned.max_tautomers(), 2);
        assert!(!cloned.remove_sp3_stereo());
    }

    #[test]
    fn enumerator_configuration_clone_from_replaces_catalog_options_and_callback() {
        struct Continue;
        impl TautomerEnumerationCallback for Continue {
            fn should_continue(&self, _molecule: &Molecule, _result: &TautomerEnumeration) -> bool {
                true
            }
        }

        let callback = Continue;
        let mut source =
            TautomerEnumerator::from_options(TautomerOptions::default().with_max_transforms(9));
        source.set_callback(Some(&callback));
        let mut target = TautomerEnumerator::from_catalog(TautomerCatalog::v1().unwrap());
        target.clone_from(&source);

        assert!(Arc::ptr_eq(&source.catalog, &target.catalog));
        assert_eq!(target.max_transforms(), 9);
        assert!(std::ptr::eq(
            target.callback().expect("copied callback"),
            &callback as &dyn TautomerEnumerationCallback
        ));
    }

    #[test]
    fn enumerator_configuration_callback_is_borrowed_and_replaceable() {
        struct Decision(bool);
        impl TautomerEnumerationCallback for Decision {
            fn should_continue(&self, _molecule: &Molecule, _result: &TautomerEnumeration) -> bool {
                self.0
            }
        }

        let first = Decision(true);
        let second = Decision(false);
        let mut enumerator = TautomerEnumerator::new();
        enumerator.set_callback(Some(&first));
        assert!(std::ptr::eq(
            enumerator.callback().expect("first callback"),
            &first as &dyn TautomerEnumerationCallback
        ));
        enumerator.set_callback(Some(&second));
        assert!(std::ptr::eq(
            enumerator.callback().expect("replacement callback"),
            &second as &dyn TautomerEnumerationCallback
        ));
        enumerator.set_callback(None);
        assert!(enumerator.callback().is_none());
    }

    #[test]
    fn transform_construction_preserves_name_query_and_source_ordered_edits() {
        let transform = TautomerTransform::new(
            "ordered",
            query("[O]-[C]=[C]"),
            vec![BondOrder::Double, BondOrder::Single],
            vec![-1, 0, 1],
        )
        .expect("valid transform");

        assert_eq!(transform.name(), "ordered");
        assert_eq!(transform.query().num_atoms(), 3);
        assert_eq!(transform.query().num_bonds(), 2);
        assert_eq!(
            transform.bond_types(),
            &[BondOrder::Double, BondOrder::Single]
        );
        assert_eq!(transform.charges(), &[-1, 0, 1]);
    }

    #[test]
    fn transform_clone_has_independent_value_semantics() {
        let original = TautomerTransform::new(
            "original",
            query("[O]-[C]=[C]"),
            vec![BondOrder::Double, BondOrder::Single],
            Vec::new(),
        )
        .expect("valid transform");
        let mut cloned = original.clone();

        cloned
            .query
            .atom_mut(0)
            .expect("query atom")
            .atom_mut()
            .set_formal_charge(-1);
        cloned.query = cloned.query.with_name("clone");
        cloned.bond_types[0] = BondOrder::Triple;

        assert_eq!(original.name(), "original");
        assert_eq!(original.query().atoms()[0].atom().formal_charge(), 0);
        assert_eq!(original.bond_types()[0], BondOrder::Double);
        assert_eq!(cloned.name(), "clone");
        assert_eq!(cloned.query().atoms()[0].atom().formal_charge(), -1);
        assert_eq!(cloned.bond_types()[0], BondOrder::Triple);
    }

    #[test]
    fn transform_clone_from_replaces_every_source_owned_field() {
        let source = TautomerTransform::new(
            "source",
            query("[N]-[C]=[O]"),
            vec![BondOrder::Double, BondOrder::Single],
            vec![1, 0, -1],
        )
        .expect("valid source transform");
        let mut target =
            TautomerTransform::new("target", query("[O]-[C]=[C]"), Vec::new(), Vec::new())
                .expect("valid target transform");

        target.clone_from(&source);

        assert_eq!(target, source);
    }

    #[test]
    fn transform_keeps_and_reuses_the_compiled_query_value() {
        let compiled = query("[O]-[C]=[C]").with_prop("compiled-sentinel", "present");
        let shared_atoms = compiled.atoms().as_ptr();

        let transform = TautomerTransform::new("compiled", compiled, Vec::new(), Vec::new())
            .expect("valid transform");

        assert_eq!(transform.query().prop("compiled-sentinel"), Some("present"));
        assert_eq!(transform.query.atoms().as_ptr(), shared_atoms);
    }

    #[test]
    fn transform_accepts_empty_bond_and_charge_edit_vectors() {
        let transform =
            TautomerTransform::new("alternating", query("[O]-[C]=[C]"), Vec::new(), Vec::new())
                .expect("empty edits select source defaults");

        assert!(transform.bond_types().is_empty());
        assert!(transform.charges().is_empty());
    }

    #[test]
    fn transform_rejects_missing_donor_or_acceptor() {
        let error = TautomerTransform::new("invalid", query("[O]"), Vec::new(), Vec::new())
            .expect_err("one query atom cannot provide both endpoints");

        assert_eq!(
            error,
            TautomerTransformError::MissingDonorOrAcceptor { actual: 1 }
        );
    }

    #[test]
    fn transform_rejects_mismatched_bond_edit_count() {
        let error = TautomerTransform::new(
            "invalid",
            query("[O]-[C]=[C]"),
            vec![BondOrder::Double],
            Vec::new(),
        )
        .expect_err("explicit bond edits must align to query bonds");

        assert_eq!(
            error,
            TautomerTransformError::BondEditCount {
                expected: 2,
                actual: 1
            }
        );
    }

    #[test]
    fn transform_rejects_mismatched_charge_edit_count() {
        let error =
            TautomerTransform::new("invalid", query("[O]-[C]=[C]"), Vec::new(), vec![1, -1])
                .expect_err("explicit charge edits must align to query atoms");

        assert_eq!(
            error,
            TautomerTransformError::ChargeEditCount {
                expected: 3,
                actual: 2
            }
        );
    }
}
