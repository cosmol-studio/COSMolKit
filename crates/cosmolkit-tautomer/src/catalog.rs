//! Source-backed tautomer transform catalog over detached query values.
use crate::transforms::{
    CURRENT_TAUTOMER_TRANSFORM_DEFINITIONS, TautomerTransformDefinition,
    V1_TAUTOMER_TRANSFORM_DEFINITIONS,
};
use cosmolkit_model::QueryGraph;
use cosmolkit_search::{
    CompiledQuery, QueryCompileError, SmartsParseError, SmartsParseParams, parse_smarts,
};
use cosmolkit_types::BondOrder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::string::FromUtf8Error;

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
pub enum TautomerTransformError {
    #[error("Charge symbol not recognised.")]
    ChargeSymbolNotRecognised,
    #[error("cannot parse tautomer SMARTS: {smarts}")]
    CannotParseSmarts {
        smarts: String,
        #[source]
        source: SmartsParseError,
    },
    #[error(transparent)]
    QueryCompile(#[from] QueryCompileError),
}

fn string_to_bond_types(bond_str: &str) -> Vec<BondOrder> {
    // RDKit✔️✔️: std::vector<Bond::BondType> stringToBondType(std::string bond_str) {
    // RDKit✔️✔️:   std::vector<Bond::BondType> bonds;
    // RDKit✔️✔️:   for (const auto &c : bond_str) {
    // RDKit✔️✔️:     switch (c) {
    // RDKit✔️✔️:       case '-':
    // RDKit✔️✔️:         bonds.push_back(Bond::SINGLE);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case '=':
    // RDKit✔️✔️:         bonds.push_back(Bond::DOUBLE);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case '#':
    // RDKit✔️✔️:         bonds.push_back(Bond::TRIPLE);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case ':':
    // RDKit✔️✔️:         bonds.push_back(Bond::AROMATIC);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return bonds;
    // RDKit✔️✔️: }
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
    // RDKit✔️✔️: std::vector<int> stringToCharge(std::string charge_str) {
    // RDKit✔️✔️:   std::vector<int> charges;
    // RDKit✔️✔️:   for (const auto &c : charge_str) {
    // RDKit✔️✔️:     switch (c) {
    // RDKit✔️✔️:       case '+':
    // RDKit✔️✔️:         charges.push_back(1);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case '0':
    // RDKit✔️✔️:         charges.push_back(0);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case '-':
    // RDKit✔️✔️:         charges.push_back(-1);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         throw ValueErrorException("Charge symbol not recognised.");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return charges;
    // RDKit✔️✔️: }
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
    // RDKit✔️✔️: std::unique_ptr<MolStandardize::TautomerTransform> getTautomer(
    // RDKit✔️✔️:     const std::string &name, const std::string &smarts,
    // RDKit✔️✔️:     const std::string &bond_str, const std::string &charge_str) {
    // RDKit✔️✔️:   std::vector<Bond::BondType> bond_types =
    // RDKit✔️✔️:       MolStandardize::stringToBondType(bond_str);
    // RDKit✔️✔️:   std::vector<int> charges = MolStandardize::stringToCharge(charge_str);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ROMol *tautomer = SmartsToMol(smarts);
    // RDKit✔️✔️:   if (!tautomer) {
    // RDKit✔️✔️:     throw ValueErrorException("cannot parse tautomer SMARTS: " + smarts);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   tautomer->setProp(common_properties::_Name, name);
    // RDKit✔️✔️:   return std::make_unique<MolStandardize::TautomerTransform>(
    // RDKit✔️✔️:       tautomer, bond_types, charges);
    // RDKit✔️✔️: }
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
    // RDKit✔️✔️:   const int MAX_LINE_LEN = 512;
    // RDKit✔️✔️:   char inLine[MAX_LINE_LEN];
    // RDKit✔️✔️:     inStream.getline(inLine, MAX_LINE_LEN, '\n');
    // C++ getline tests EOF and the delimiter before the capacity limit.
    // ROOT-approved migration corrects the archived helper's 511+newline bug.
    // O(bytes) buffered scan, fixed 511-byte allocation; no unbounded line read.
    const SOURCE_PAYLOAD_LIMIT: usize = 511;
    let mut line = Vec::with_capacity(SOURCE_PAYLOAD_LIMIT);
    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            return Ok(((!line.is_empty()).then_some(line), false));
        }
        if available[0] == b'\n' {
            reader.consume(1);
            return Ok((Some(line), false));
        }
        if line.len() == SOURCE_PAYLOAD_LIMIT {
            return Ok((Some(line), true));
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
        // RDKit✔️✔️:     tmpstr = inLine;
        // std::string(char*) stops at the first NUL, although getline consumes
        // the entire physical line. Ignore only that source-defined suffix.
        let prefix_len = line
            .iter()
            .position(|byte| *byte == 0)
            .unwrap_or(line.len());
        let line = String::from_utf8(line[..prefix_len].to_vec())?;

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

/// One immutable source transform. Catalog construction does not validate
/// donor/acceptor or edit dimensions; application owns reached unsafe boundaries.
#[derive(Debug, PartialEq)]
pub struct TautomerTransform {
    pub(crate) query: CompiledQuery,
    bond_types: Vec<BondOrder>,
    charges: Vec<i32>,
}
#[derive(Debug, PartialEq)]
pub struct TautomerCatalog {
    transforms: Vec<TautomerTransform>,
}
impl TautomerCatalog {
    /// Source zero-argument CatalogParams constructor, distinct from file="".
    pub fn empty() -> Self {
        // RDKit✔️✔️:   TautomerCatalogParams() {
        // RDKit✔️✔️:     d_typeStr = "Tautomer Catalog Parameters";
        // RDKit✔️✔️:     d_transforms.clear();
        // RDKit✔️✔️:   }
        Self {
            transforms: Vec::new(),
        }
    }
    /// Read a bounded source-format stream; negative limits mean unbounded.
    pub fn from_reader(
        reader: &mut impl BufRead,
        n_to_read: i32,
    ) -> Result<Self, TautomerCatalogError> {
        Ok(Self {
            transforms: read_transforms(reader, n_to_read)?,
        })
    }
    /// Source stream deserialization has no implementation.
    pub fn deserialize_from_reader(
        _reader: &mut impl BufRead,
    ) -> Result<Self, TautomerCatalogError> {
        // RDKit✔️✔️: void TautomerCatalogParams::initFromStream(std::istream &) {
        // RDKit✔️✔️:   UNDER_CONSTRUCTION("not implemented");
        // RDKit✔️✔️: }
        Err(TautomerCatalogError::DeserializationUnderConstruction)
    }
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
        // RDKit✔️✔️: const TautomerTransform TautomerCatalogParams::getTransform(
        // RDKit✔️✔️:     unsigned int fid) const {
        // RDKit✔️✔️:   URANGE_CHECK(fid, d_transforms.size());
        // RDKit✔️✔️:   return d_transforms[fid];  //.get();
        // RDKit✔️✔️: }
        // Detached query/plan and edit vectors are independently cloned.
        // Source and Rust are O(query + edits); no Arc/COW speedup claim.
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
        // RDKit✔️✔️: void TautomerCatalogParams::initFromString(const std::string &) {
        // RDKit✔️✔️:   UNDER_CONSTRUCTION("not implemented");
        // RDKit✔️✔️: }
        Err(TautomerCatalogError::DeserializationUnderConstruction)
    }
}

impl Default for TautomerCatalog {
    fn default() -> Self {
        Self::empty()
    }
}

impl Clone for TautomerCatalog {
    fn clone(&self) -> Self {
        // RDKit✔️✔️: TautomerCatalogParams::TautomerCatalogParams(
        // RDKit✔️✔️:     const TautomerCatalogParams &other) {
        // RDKit✔️✔️:   d_typeStr = other.d_typeStr;
        // RDKit✔️✔️:   d_transforms.clear();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const std::vector<TautomerTransform> &transforms = other.getTransforms();
        // RDKit✔️✔️:   for (const auto &transform : transforms) {
        // RDKit✔️✔️:     d_transforms.push_back(transform);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // One Vec allocation plus independent query/plan clones, O(total rows).
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
        name: impl Into<cosmolkit_model::PropertyText>,
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
        // The source constructor stores these values without donor or edit
        // dimension checks. Query compilation is delegated to the existing Q
        // owner; valid empty/single queries and arbitrary edit rows are kept.

        Ok(Self {
            query: CompiledQuery::compile(query.with_name(name))?,
            bond_types,
            charges,
        })
    }

    #[must_use]
    pub fn name(&self) -> &cosmolkit_model::PropertyText {
        // The private constructor always installs this String-tag property,
        // including an empty name; callers borrow the immutable query only.
        // Return the same counted bytes without decoding or absent fallback.
        self.query
            .query()
            .name()
            .expect("transform constructor installs a String-tag name")
            .expect("transform constructor installs its name")
    }

    #[must_use]
    pub fn query(&self) -> &QueryGraph {
        self.query.query()
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
        // RDKit✔️✔️: TautomerTransform(const TautomerTransform &other)
        // RDKit✔️✔️:     : BondTypes(other.BondTypes), Charges(other.Charges) {
        // RDKit✔️✔️:   Mol = new ROMol(*other.Mol);
        // RDKit✔️✔️: }
        //
        // QueryGraph and the compiled VF2 plan own Vec storage. Cloning the
        // extra plan allocates more than the source ROMol copy: RDKit✔️❌.
        Self {
            query: self.query.clone(),
            bond_types: self.bond_types.clone(),
            charges: self.charges.clone(),
        }
    }

    fn clone_from(&mut self, source: &Self) {
        // RDKit✔️✔️: TautomerTransform &operator=(const TautomerTransform &other) {
        // RDKit✔️✔️:   if (this != &other) {
        // RDKit✔️✔️:     delete Mol;
        // RDKit✔️✔️:     Mol = new ROMol(*other.Mol);
        // RDKit✔️✔️:     BondTypes = other.BondTypes;
        // RDKit✔️✔️:     Charges = other.Charges;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return *this;
        // RDKit✔️✔️: }
        //
        // Rust borrowing prevents an ordinary self-assignment alias. Reusing
        // each allocation through `clone_from` is at least as efficient as the
        // source delete-and-deep-copy path and leaves identical value state.
        self.query.clone_from(&source.query);
        self.bond_types.clone_from(&source.bond_types);
        self.charges.clone_from(&source.charges);
    }
}
