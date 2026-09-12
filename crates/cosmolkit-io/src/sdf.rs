//! Detached V2000 molfile/SDF reader.

use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate,
    BondSpec, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, MoleculeProperties,
    QueryAtom, QueryBond, QueryGraph, QueryNode, SdfPropertyList, SdfPropertyListTarget,
    SubstanceGroup, TopologyBlock,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, Element};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SdfReadError {
    #[error("empty molfile block")]
    Empty,
    #[error("invalid V2000 counts line")]
    Counts,
    #[error("invalid {kind} field on line {line}: {value}")]
    Field {
        kind: &'static str,
        line: usize,
        value: String,
    },
    #[error("unsupported molfile format: {0}")]
    Unsupported(&'static str),
    #[error("{0}")]
    Parse(String),
    #[error("invalid detached topology: {0}")]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error("invalid detached coordinates: {0}")]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error("invalid detached query graph: {0}")]
    QueryGraph(#[from] cosmolkit_model::QueryGraphError),
    #[error("invalid atom property: {0}")]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error("invalid bond value: {0}")]
    BondValue(#[from] cosmolkit_model::BondValueError),
    #[error("invalid molecule property: {0}")]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SdfWriteError {
    #[error("V2000 mol block supports at most 999 atoms, bonds, and substance groups")]
    CountLimit,
    #[error("detached topology is invalid: {0}")]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error("detached coordinates are invalid: {0}")]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error("unsupported bond order {0:?}")]
    BondOrder(BondOrder),
    #[error("unsupported detached molfile atom state: {0}")]
    Atom(&'static str),
    #[error("unsupported detached substance-group state: {0}")]
    SubstanceGroup(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    pub data_fields: Vec<(String, String)>,
}

/// Query-bearing MolBlock state which cannot be installed into a concrete
/// `Molecule` without losing query semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct QueryMolBlockRecord {
    pub query: QueryGraph,
    pub substance_groups: Vec<SubstanceGroup>,
    pub properties: MoleculeProperties,
    pub source_coordinate_dim: Option<CoordinateDimension>,
}

/// Detached MolBlock parse result. Concrete and query graphs remain distinct
/// at the type boundary; neither is lowered into the other.
#[derive(Debug, Clone, PartialEq)]
pub enum MolBlockRecord {
    Concrete {
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    },
    Query(QueryMolBlockRecord),
}

/// Syntax-only MolBlock reader policy.
///
/// Sanitization, hydrogen removal, attachment-point expansion, and coordinate
/// coercion require live-molecule operations and intentionally remain outside
/// this detached parameter object.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolBlockReadParams {
    /// Apply RDKit's strict MolBlock syntax checks. Non-strict mode accepts
    /// the explicitly source-backed legacy relaxations only.
    pub strict_parsing: bool,
}

impl Default for MolBlockReadParams {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION MolFileParserParams
        // RDKit✔️✔️:   bool strictParsing = true; /**< if set to false, the parser is more lax about
        // RDKit✔️✔️:                                  correctness of the contents. */
        // END RDKIT CPP FUNCTION
        Self {
            strict_parsing: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfGraphRecord {
    pub mol_block: MolBlockRecord,
    pub data_fields: Vec<(String, String)>,
}

/// Controls SDF data-field handling after a MolBlock has been parsed.
///
/// Chemistry finalization options such as sanitization, hydrogen removal, and
/// attachment-point expansion belong to the live-molecule runtime and are
/// intentionally not represented at this detached syntax boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SdfDataReadParams {
    /// Apply strict MolBlock checks and reject non-header content between SDF
    /// data fields.
    pub strict_parsing: bool,
    /// Apply `atom.*prop.*` and `bond.*prop.*` lists to detached graph items.
    pub process_property_lists: bool,
}

impl Default for SdfDataReadParams {
    fn default() -> Self {
        // RDKit source: FileParsers.h `MolFileParserParams` and
        // ForwardSDMolSupplier.cpp `ForwardSDMolSupplier::init`.
        // RDKit✔️✔️:   bool strictParsing = true; /**< if set to false, the parser is more lax about
        // RDKit✔️✔️:   df_processPropertyLists = true;
        Self {
            strict_parsing: true,
            process_property_lists: true,
        }
    }
}

impl From<SdfDataReadParams> for MolBlockReadParams {
    fn from(params: SdfDataReadParams) -> Self {
        Self {
            strict_parsing: params.strict_parsing,
        }
    }
}

/// Location and title information for one indexed SDF record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfRecordMetadata {
    pub index: usize,
    pub byte_offset: u64,
    pub byte_len: u64,
    pub line_offset: usize,
    pub line_len: usize,
    pub title: Option<String>,
}

/// Exact text and consumed span for one forward SDF record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfRecordText {
    pub text: String,
    pub byte_len: u64,
    pub line_len: usize,
    pub hit_eof: bool,
}

/// Forward-only, query-preserving detached SDF reader.
pub struct SdfGraphReader<R> {
    reader: R,
    params: SdfDataReadParams,
    next_index: usize,
    byte_offset: u64,
    line_offset: usize,
    end: bool,
}

impl<R: BufRead> SdfGraphReader<R> {
    #[must_use]
    pub fn new(reader: R) -> Self {
        Self::with_params(reader, SdfDataReadParams::default())
    }

    #[must_use]
    pub const fn with_params(reader: R, params: SdfDataReadParams) -> Self {
        Self {
            reader,
            params,
            next_index: 0,
            byte_offset: 0,
            line_offset: 0,
            end: false,
        }
    }

    #[must_use]
    pub const fn params(&self) -> SdfDataReadParams {
        self.params
    }

    #[must_use]
    pub const fn is_end(&self) -> bool {
        self.end
    }

    /// Read and consume the next record.
    ///
    /// A parse error consumes that record through its delimiter, so a
    /// subsequent call continues at the following record like RDKit's forward
    /// supplier.
    pub fn next_record(&mut self) -> Result<Option<SdfGraphRecord>, SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ForwardSDMolSupplier::_next
        // RDKit✔️❌:   if (dp_inStream->eof()) {
        // RDKit✔️❌:     df_end = true;
        // RDKit✔️❌:     return res;
        // RDKit✔️❌:   }
        // RDKit✔️❌:   unsigned int line = d_line;
        // RDKit✔️❌:     MolFromMolDataStream(*dp_inStream, line, d_params).swap(res);
        // RDKit✔️❌:     d_line = line;
        // RDKit✔️❌:     if (res) {
        // RDKit✔️❌:       this->readMolProps(*res);
        // RDKit✔️❌:     }
        // This detached reader buffers one full record before parsing it, so
        // recovery and ordering match while memory use is worse than RDKit's
        // stream-positioned MolBlock/data-field path.
        if self.end {
            return Ok(None);
        }
        let Some(raw) = read_sdf_record_text(&mut self.reader)? else {
            self.end = true;
            return Ok(None);
        };
        self.next_index += 1;
        self.byte_offset += raw.byte_len;
        self.line_offset += raw.line_len;
        self.end = raw.hit_eof;
        let record = read_sdf_graph_record_detached_with_params(&raw.text, self.params)?;
        // END RDKIT CPP FUNCTION
        Ok(Some(record))
    }

    #[must_use]
    pub const fn records_consumed(&self) -> usize {
        self.next_index
    }

    #[must_use]
    pub const fn bytes_consumed(&self) -> u64 {
        self.byte_offset
    }

    #[must_use]
    pub const fn lines_consumed(&self) -> usize {
        self.line_offset
    }
}

/// Random-access index over a detached SDF file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfGraphDataset {
    path: PathBuf,
    params: SdfDataReadParams,
    metadata: Vec<SdfRecordMetadata>,
}

impl SdfGraphDataset {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, SdfReadError> {
        Self::open_with_params(path, SdfDataReadParams::default())
    }

    pub fn open_with_params(
        path: impl AsRef<Path>,
        params: SdfDataReadParams,
    ) -> Result<Self, SdfReadError> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|error| SdfReadError::Parse(error.to_string()))?;
        let mut reader = BufReader::new(file);
        let metadata = index_sdf_records(&mut reader)?;
        Ok(Self {
            path: path.to_path_buf(),
            params,
            metadata,
        })
    }

    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    #[must_use]
    pub const fn params(&self) -> SdfDataReadParams {
        self.params
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.metadata.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.metadata.is_empty()
    }

    #[must_use]
    pub fn metadata(&self, index: usize) -> Option<&SdfRecordMetadata> {
        self.metadata.get(index)
    }

    pub fn metadata_iter(&self) -> impl Iterator<Item = &SdfRecordMetadata> {
        self.metadata.iter()
    }

    pub fn record(&self, index: usize) -> Result<SdfGraphRecord, SdfReadError> {
        self.record_with_params(index, self.params)
    }

    pub fn record_with_params(
        &self,
        index: usize,
        params: SdfDataReadParams,
    ) -> Result<SdfGraphRecord, SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION SDMolSupplier::operator[]
        // RDKit✔️✔️:   PRECONDITION(dp_inStream, "no stream");
        // RDKit✔️✔️:   // get the molecule with index idx
        // RDKit✔️✔️:   moveTo(idx);
        // RDKit✔️✔️:   return next();
        let text = self.record_text(index)?;
        let record = read_sdf_graph_record_detached_with_params(&text, params)?;
        // END RDKIT CPP FUNCTION
        Ok(record)
    }

    /// Read the exact indexed record text, including its delimiter line.
    ///
    /// This is the adapter boundary for runtimes that need to apply a
    /// different chemistry-finalization policy after detached record framing.
    pub fn record_text(&self, index: usize) -> Result<String, SdfReadError> {
        let metadata = self.metadata.get(index).ok_or_else(|| {
            SdfReadError::Parse(format!(
                "ERROR: Index error (idx = {index}) :  we do not have enough mol blocks"
            ))
        })?;
        let mut file =
            File::open(&self.path).map_err(|error| SdfReadError::Parse(error.to_string()))?;
        file.seek(SeekFrom::Start(metadata.byte_offset))
            .map_err(|error| SdfReadError::Parse(error.to_string()))?;
        let mut text = String::with_capacity(metadata.byte_len as usize);
        file.take(metadata.byte_len)
            .read_to_string(&mut text)
            .map_err(|error| SdfReadError::Parse(error.to_string()))?;
        Ok(text)
    }
}

fn field<T: std::str::FromStr>(
    value: &str,
    kind: &'static str,
    line: usize,
) -> Result<T, SdfReadError> {
    value.trim().parse().map_err(|_| SdfReadError::Field {
        kind,
        line,
        value: value.to_string(),
    })
}

pub(super) fn rdkit_substr(text: &str, start: usize, len: usize) -> &str {
    if start >= text.len() {
        return "";
    }
    let end = start.saturating_add(len).min(text.len());
    text.get(start..end).unwrap_or("")
}

pub(super) fn parse_rdkit_unsigned(text: &str) -> Result<u32, ()> {
    let checked = text.split_once('\0').map_or(text, |(prefix, _)| prefix);
    if !checked
        .bytes()
        .all(|byte| byte.is_ascii_digit() || byte == b'+' || byte == b' ')
    {
        return Err(());
    }
    let input = checked.trim_start_matches(' ');
    if input.is_empty() {
        return Ok(0);
    }
    let input = input.strip_prefix('+').unwrap_or(input);
    let digit_count = input.bytes().take_while(u8::is_ascii_digit).count();
    if digit_count == 0 {
        return Ok(0);
    }
    Ok(input[..digit_count].parse().unwrap_or(0))
}

pub(super) fn parse_rdkit_int(text: &str) -> Result<i32, ()> {
    let checked = text.split_once('\0').map_or(text, |(prefix, _)| prefix);
    if !checked
        .bytes()
        .all(|byte| byte.is_ascii_digit() || byte == b'+' || byte == b'-' || byte == b' ')
    {
        return Err(());
    }
    let input = checked.trim_start_matches(' ');
    if input.is_empty() {
        return Ok(0);
    }
    let sign_len = usize::from(input.starts_with(['+', '-']));
    let digit_count = input[sign_len..]
        .bytes()
        .take_while(u8::is_ascii_digit)
        .count();
    if digit_count == 0 {
        return Ok(0);
    }
    Ok(input[..sign_len + digit_count].parse().unwrap_or(0))
}

pub(super) fn parse_rdkit_double(text: &str) -> Result<f64, ()> {
    let checked = text.split_once('\0').map_or(text, |(prefix, _)| prefix);
    if !checked
        .bytes()
        .all(|byte| byte.is_ascii_digit() || matches!(byte, b'+' | b'-' | b'.' | b',' | b' '))
    {
        return Err(());
    }
    let input = checked.trim_start_matches(' ');
    let numeric_prefix = input.split_once(',').map_or(input, |(prefix, _)| prefix);
    if numeric_prefix.is_empty() || matches!(numeric_prefix, "+" | "-" | ".") {
        return Ok(0.0);
    }
    numeric_prefix.parse().map_err(|_| ())
}

fn parse_required_int(
    line: &str,
    start: usize,
    len: usize,
    line_number: usize,
) -> Result<i32, SdfReadError> {
    let value = rdkit_substr(line, start, len);
    parse_rdkit_int(value).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{value}' to int on line {line_number}"
        ))
    })
}

fn parse_required_unsigned(
    line: &str,
    start: usize,
    len: usize,
    line_number: usize,
) -> Result<u32, SdfReadError> {
    let value = rdkit_substr(line, start, len);
    parse_rdkit_unsigned(value).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{value}' to unsigned int on line {line_number}"
        ))
    })
}

#[derive(Debug, Clone)]
struct ParsedV2000Atom {
    spec: AtomSpec,
    query: Option<QueryNode<AtomQueryPredicate>>,
    coordinate: [f64; 3],
}

fn update_v2000_atom_spec(
    atom: &mut ParsedV2000Atom,
    update: impl FnOnce(AtomSpec) -> Result<AtomSpec, cosmolkit_model::AtomPropertyError>,
) -> Result<(), SdfReadError> {
    let spec = std::mem::replace(&mut atom.spec, AtomSpec::new(Element::DUMMY));
    atom.spec = update(spec)?;
    Ok(())
}

#[derive(Debug, Clone)]
struct ParsedV2000Bond {
    bond: Bond,
    query: Option<QueryNode<BondQueryPredicate>>,
}

fn strip_sdf_line(line: &str) -> &str {
    line.trim_matches([' ', '\t', '\r', '\n'])
}

fn strip_terminal_cr(line: &str) -> &str {
    line.strip_suffix('\r').unwrap_or(line)
}

fn starts_with_sdf_continuation_space(line: &str) -> bool {
    line.starts_with(' ') || line.starts_with('\t')
}

fn is_sdf_record_delimiter(line: &str) -> bool {
    line.starts_with("$$$$")
}

pub fn read_sdf_record_text<R: BufRead>(
    reader: &mut R,
) -> Result<Option<SdfRecordText>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION MultithreadedSDMolSupplier::extractNextRecord
    // RDKit✔️❌:   std::string currentStr, prevStr;
    // RDKit✔️❌:   record = "";
    // RDKit✔️❌:   lineNum = d_line;
    // RDKit✔️❌:   while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❌:          ((prevStr.find_first_not_of(" \t\r\n") != std::string::npos &&
    // RDKit✔️❌:            prevStr.find("M  END") != 0) ||
    // RDKit✔️❌:           currentStr[0] != '$' || currentStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❌:     prevStr = currentStr;
    // RDKit✔️❌:     std::getline(*dp_inStream, currentStr);
    // RDKit✔️❌:     record += currentStr + "\n";
    // RDKit✔️❌:     ++d_line;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (record.find_first_not_of("\n\r") == std::string::npos) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // Line-oriented buffering preserves record boundaries and recovery, but
    // allocates a second full copy before detached parsing.
    let mut text = String::new();
    let mut byte_len = 0_u64;
    let mut line_len = 0_usize;
    let mut hit_eof = false;
    loop {
        let mut line = String::new();
        let read = reader
            .read_line(&mut line)
            .map_err(|error| SdfReadError::Parse(error.to_string()))?;
        if read == 0 {
            hit_eof = true;
            break;
        }
        byte_len += read as u64;
        line_len += 1;
        let content = strip_terminal_cr(line.strip_suffix('\n').unwrap_or(&line));
        let delimiter = is_sdf_record_delimiter(content);
        text.push_str(&line);
        if !line.ends_with('\n') {
            hit_eof = true;
        }
        if delimiter {
            break;
        }
        if hit_eof {
            break;
        }
    }
    if text
        .chars()
        .all(|character| matches!(character, '\n' | '\r'))
    {
        return Ok(None);
    }
    // END RDKIT CPP FUNCTION
    Ok(Some(SdfRecordText {
        text,
        byte_len,
        line_len,
        hit_eof,
    }))
}

/// Scan SDF record boundaries and return reusable byte/line metadata.
pub fn index_sdf_records<R: BufRead>(
    reader: &mut R,
) -> Result<Vec<SdfRecordMetadata>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION SDMolSupplier::buildIndexTo
    // RDKit✔️❌:       constexpr char dollarSigns[]{"$$$$"};
    // RDKit✔️❌:       auto match = std::search(ptr, bufEnd, dollarSigns, dollarSigns + 4);
    // RDKit✔️❌:       if (match == bufEnd) {
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (*(match - 1) == '\n') {  // ensure $$$$ is at start of line
    // RDKit✔️❌:             d_molpos.push_back(posHold);
    // RDKit scans fixed-size byte chunks; this implementation performs the
    // same O(file-size) boundary scan one line at a time and therefore has
    // higher allocation and dispatch overhead.
    let mut metadata = Vec::new();
    let mut byte_offset = 0_u64;
    let mut line_offset = 0_usize;
    while let Some(raw) = read_sdf_record_text(reader)? {
        let title = raw
            .text
            .lines()
            .next()
            .map(strip_terminal_cr)
            .filter(|title| !title.is_empty())
            .map(ToOwned::to_owned);
        metadata.push(SdfRecordMetadata {
            index: metadata.len(),
            byte_offset,
            byte_len: raw.byte_len,
            line_offset,
            line_len: raw.line_len,
            title,
        });
        byte_offset += raw.byte_len;
        line_offset += raw.line_len;
        if raw.hit_eof {
            break;
        }
    }
    // END RDKIT CPP FUNCTION
    Ok(metadata)
}

fn parse_sdf_data_header(line: &str) -> Option<String> {
    // BEGIN RDKIT CPP FUNCTION ForwardSDMolSupplier::readMolProps header extraction
    // RDKit✔️✔️: tempStr = FileParserUtils::strip(tempStr);
    // RDKit✔️✔️: if (!tempStr.empty()) {
    // RDKit✔️✔️:   if (tempStr.at(0) == '>') {  // data header line: start of a data item
    // RDKit✔️✔️:     tempStr = tempStr.substr(1);     // remove the first ">" sign
    // RDKit✔️✔️:     size_t sl = tempStr.find("<");   // begin datalabel
    // RDKit✔️✔️:     size_t se = tempStr.rfind(">");  // end datalabel
    // RDKit✔️✔️:     if ((sl == std::string::npos) || (se == std::string::npos) ||
    // RDKit✔️✔️:         (se == (sl + 1))) {
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       dlabel = tempStr.substr(sl + 1, se - sl - 1);
    let temp_str = strip_sdf_line(line).strip_prefix('>')?;
    let start = temp_str.find('<')?;
    let end = temp_str.rfind('>')?;
    if end == start + 1 {
        return None;
    }
    if end < start {
        // `size_t` subtraction underflows in the source, so `substr()` keeps
        // the remainder beginning after '<'.
        return Some(temp_str[start + 1..].to_owned());
    }
    Some(temp_str[start + 1..end].to_owned())
    // END RDKIT CPP FUNCTION
}

fn parse_sdf_data_fields(
    lines: &[&str],
    params: SdfDataReadParams,
) -> Result<Vec<(String, String)>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ForwardSDMolSupplier::readMolProps
    // RDKit✔️❌: while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❌:        (tempStr.empty() || tempStr.at(0) != '$' ||
    // RDKit✔️❌:         tempStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❌:   tempStr = FileParserUtils::strip(tempStr);
    // RDKit✔️❌:   if (!tempStr.empty()) {
    // RDKit✔️❌:     if (tempStr.at(0) == '>') {  // data header line: start of a data item
    // RDKit✔️❌:       hasProp = true;
    // RDKit✔️❌:       warningIssued = false;
    // RDKit✔️❌:       tempStr = tempStr.substr(1);     // remove the first ">" sign
    // RDKit✔️❌:       size_t sl = tempStr.find("<");   // begin datalabel
    // RDKit✔️❌:       size_t se = tempStr.rfind(">");  // end datalabel
    // RDKit✔️❌:       if ((sl == std::string::npos) || (se == std::string::npos) ||
    // RDKit✔️❌:           (se == (sl + 1))) {
    // RDKit✔️❌:         d_line++;
    // RDKit✔️❌:         std::getline(*dp_inStream, inl);
    // RDKit✔️❌:         tempStr = inl;
    // RDKit✔️❌:         auto stmp = FileParserUtils::strip(tempStr);
    // RDKit✔️❌:         while (stmp.length() != 0) {
    // RDKit✔️❌:           d_line++;
    // RDKit✔️❌:           std::getline(*dp_inStream, inl);
    // RDKit✔️❌:           tempStr = inl;
    // RDKit✔️❌:           if (dp_inStream->eof()) {
    // RDKit✔️❌:             throw FileParseException("End of data field name not found");
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         dlabel = tempStr.substr(sl + 1, se - sl - 1);
    // RDKit✔️❌:         d_line++;
    // RDKit✔️❌:         std::getline(*dp_inStream, inl);
    // RDKit✔️❌:         tempStr = inl;
    // RDKit✔️❌:         std::string prop = "";
    // RDKit✔️❌:         auto stmp = FileParserUtils::strip(tempStr);
    // RDKit✔️❌:         int nplines = 0;  // number of lines for this property
    // RDKit✔️❌:         while (!stmp.empty() ||
    // RDKit✔️❌:                (!tempStr.empty() &&
    // RDKit✔️❌:                 (tempStr.at(0) == ' ' || tempStr.at(0) == '\t'))) {
    // RDKit✔️❌:           nplines++;
    // RDKit✔️❌:           if (nplines > 1) {
    // RDKit✔️❌:             prop += "\n";
    // RDKit✔️❌:           }
    // RDKit✔️❌:           if (!tempStr.empty()) {
    // RDKit✔️❌:             if (tempStr.back() == '\r') {
    // RDKit✔️❌:               tempStr = tempStr.substr(0, tempStr.size() - 1);
    // RDKit✔️❌:             }
    // RDKit✔️❌:             prop += tempStr;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           d_line++;
    // RDKit✔️❌:           inl.erase();
    // RDKit✔️❌:           std::getline(*dp_inStream, inl);
    // RDKit✔️❌:           tempStr = inl;
    // RDKit✔️❌:           if (tempStr.empty()) {
    // RDKit✔️❌:             stmp = tempStr;
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             stmp = FileParserUtils::strip(tempStr);
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:         mol.setProp(dlabel, prop);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       if (d_params.strictParsing) {
    // RDKit✔️❌:         throw FileParseException("Problems encountered parsing data fields");
    // RDKit✔️❌:       }
    // COSMolKit buffers a record and its line references before parsing, while
    // RDKit consumes a stream in-place. Behavior matches, but the extra record
    // allocation is materially more expensive.
    let mut fields = Vec::new();
    let mut index = 0;
    while index < lines.len() {
        let line = strip_terminal_cr(lines[index]);
        if is_sdf_record_delimiter(line) {
            break;
        }
        let stripped = strip_sdf_line(line);
        if stripped.is_empty() {
            index += 1;
            continue;
        }
        if !stripped.starts_with('>') && params.strict_parsing {
            return Err(SdfReadError::Parse(
                "Problems encountered parsing data fields".to_owned(),
            ));
        }

        if !stripped.starts_with('>') {
            index += 1;
            continue;
        }

        let Some(name) = parse_sdf_data_header(stripped) else {
            index += 1;
            while index < lines.len() && !strip_sdf_line(lines[index]).is_empty() {
                if is_sdf_record_delimiter(strip_terminal_cr(lines[index])) {
                    return Err(SdfReadError::Parse(
                        "End of data field name not found".to_owned(),
                    ));
                }
                index += 1;
            }
            if index >= lines.len() {
                return Err(SdfReadError::Parse(
                    "End of data field name not found".to_owned(),
                ));
            }
            index += 1;
            continue;
        };

        index += 1;
        let mut value = String::new();
        let mut value_line_count = 0;
        while index < lines.len() {
            let line = strip_terminal_cr(lines[index]);
            if is_sdf_record_delimiter(line) {
                break;
            }
            if strip_sdf_line(line).is_empty() && !starts_with_sdf_continuation_space(line) {
                break;
            }
            if value_line_count > 0 {
                value.push('\n');
            }
            value.push_str(line);
            value_line_count += 1;
            index += 1;
        }
        fields.push((name, value));
        index += usize::from(index < lines.len());
    }
    // END RDKIT CPP FUNCTION
    Ok(fields)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SdfPropertyListValueKind {
    String,
    Int,
    Double,
    Bool,
}

fn sdf_property_list_target(
    name: &str,
) -> Option<(SdfPropertyListTarget, &str, SdfPropertyListValueKind)> {
    // BEGIN RDKIT CPP FUNCTION processMolPropertyList
    // RDKit✔️✔️:   auto propSetter = [&](const std::string &propPrefix, auto getter,
    // RDKit✔️✔️:                         size_t nItems) {
    // RDKit✔️✔️:     std::string prefix = propPrefix + "prop.";
    // RDKit✔️✔️:     if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
    // RDKit✔️✔️:       applyMolListProp<std::string>(mol, pn, prefix, missingValueMarker, nItems,
    // RDKit✔️✔️:                                     getter);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       prefix = propPrefix + "iprop.";
    // RDKit✔️✔️:       if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
    // RDKit✔️✔️:         applyMolListProp<int>(mol, pn, prefix, missingValueMarker, nItems,
    // RDKit✔️✔️:                               getter);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         prefix = propPrefix + "dprop.";
    // RDKit✔️✔️:         if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
    // RDKit✔️✔️:           applyMolListProp<double>(mol, pn, prefix, missingValueMarker, nItems,
    // RDKit✔️✔️:                                    getter);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           prefix = propPrefix + "bprop.";
    // RDKit✔️✔️:           if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
    // RDKit✔️✔️:             applyMolListProp<bool>(mol, pn, prefix, missingValueMarker, nItems,
    // RDKit✔️✔️:                                    getter);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   if (pn.find(atomPropPrefix) == 0 && pn.length() > atomPropPrefixLength) {
    // RDKit✔️✔️:     propSetter(
    // RDKit✔️✔️:         atomPropPrefix,
    // RDKit✔️✔️:         [&mol](size_t which) { return mol.getAtomWithIdx(which); },
    // RDKit✔️✔️:         mol.getNumAtoms());
    // RDKit✔️✔️:   } else if (pn.find(bondPropPrefix) == 0 &&
    // RDKit✔️✔️:              pn.length() > bondPropPrefixLength) {
    // RDKit✔️✔️:     propSetter(
    // RDKit✔️✔️:         bondPropPrefix,
    // RDKit✔️✔️:         [&mol](size_t which) { return mol.getBondWithIdx(which); },
    // RDKit✔️✔️:         mol.getNumBonds());
    // RDKit✔️✔️:   }
    const PREFIXES: [(&str, SdfPropertyListTarget, SdfPropertyListValueKind); 8] = [
        (
            "atom.prop.",
            SdfPropertyListTarget::Atom,
            SdfPropertyListValueKind::String,
        ),
        (
            "atom.iprop.",
            SdfPropertyListTarget::Atom,
            SdfPropertyListValueKind::Int,
        ),
        (
            "atom.dprop.",
            SdfPropertyListTarget::Atom,
            SdfPropertyListValueKind::Double,
        ),
        (
            "atom.bprop.",
            SdfPropertyListTarget::Atom,
            SdfPropertyListValueKind::Bool,
        ),
        (
            "bond.prop.",
            SdfPropertyListTarget::Bond,
            SdfPropertyListValueKind::String,
        ),
        (
            "bond.iprop.",
            SdfPropertyListTarget::Bond,
            SdfPropertyListValueKind::Int,
        ),
        (
            "bond.dprop.",
            SdfPropertyListTarget::Bond,
            SdfPropertyListValueKind::Double,
        ),
        (
            "bond.bprop.",
            SdfPropertyListTarget::Bond,
            SdfPropertyListValueKind::Bool,
        ),
    ];
    let result = PREFIXES.iter().find_map(|(prefix, target, kind)| {
        name.strip_prefix(prefix)
            .filter(|property_name| !property_name.is_empty())
            .map(|property_name| (*target, property_name, *kind))
    });
    // END RDKIT CPP FUNCTION
    result
}

fn split_sdf_property_list_tokens(value: &str) -> Vec<&str> {
    let bytes = value.as_bytes();
    let mut tokens = Vec::new();
    let mut start = 0;
    let mut cursor = 0;
    while cursor < bytes.len() {
        if matches!(bytes[cursor], b' ' | b'\t' | b'\n') {
            tokens.push(&value[start..cursor]);
            while cursor < bytes.len() && matches!(bytes[cursor], b' ' | b'\t' | b'\n') {
                cursor += 1;
            }
            start = cursor;
        } else {
            cursor += 1;
        }
    }
    tokens.push(&value[start..]);
    tokens
}

fn parse_sdf_bool_property(value: &str) -> Option<String> {
    match value {
        "1" => Some("true".to_owned()),
        "0" => Some("false".to_owned()),
        _ => None,
    }
}

fn parse_sdf_property_list_values(
    value: &str,
    item_count: usize,
    value_kind: SdfPropertyListValueKind,
) -> Option<Vec<Option<String>>> {
    // BEGIN RDKIT CPP FUNCTION applyMolListProp
    // RDKit✔️✔️: void applyMolListProp(ROMol &mol, const std::string &pn,
    // RDKit✔️✔️:                       const std::string &prefix,
    // RDKit✔️✔️:                       const std::string &missingValueMarker, size_t nItems,
    // RDKit✔️✔️:                       U getter) {
    // RDKit✔️✔️:   std::string itempn = pn.substr(prefix.size());
    // RDKit✔️✔️:   std::string strVect = mol.getProp<std::string>(pn);
    // RDKit✔️✔️:   std::vector<std::string> tokens;
    // RDKit✔️✔️:   boost::split(tokens, strVect, boost::is_any_of(" \t\n"),
    // RDKit✔️✔️:                boost::token_compress_on);
    // RDKit✔️✔️:   std::string mv = missingValueMarker;
    // RDKit✔️✔️:   size_t first_token = 0;
    // RDKit✔️✔️:   if (tokens.size() == nItems + 1 && tokens[0].front() == '[' &&
    // RDKit✔️✔️:       tokens[0].back() == ']') {
    // RDKit✔️✔️:     mv = std::string(tokens[0].begin() + 1, tokens[0].end() - 1);
    // RDKit✔️✔️:     first_token = 1;
    // RDKit✔️✔️:   }
    // RDKit❗✔️:   if (mv.empty()) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "Missing value marker for property " << pn
    // RDKit❗✔️:                             << " is empty." << std::endl;
    // RDKit❗✔️:   }
    // RDKit✔️✔️:   if(tokens.size() - first_token != nItems) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "Property list " << pn << " has incompatible size, "
    // RDKit❗✔️:                             << tokens.size() << " elements found; expecting "
    // RDKit❗✔️:                             << nItems << ". Ignoring it." << std::endl;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (size_t i = first_token; i < tokens.size(); ++i) {
    // RDKit✔️✔️:     if (tokens[i] != mv) {
    // RDKit✔️✔️:       unsigned int itemid = i - first_token;
    // RDKit✔️✔️:       try {
    // RDKit✔️✔️:         T apv = boost::lexical_cast<T>(tokens[i]);
    // RDKit✔️✔️:         getter(itemid)->setProp(itempn, apv);
    // RDKit✔️✔️:       } catch (const boost::bad_lexical_cast &) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "Value " << tokens[i] << " for property " << pn << " of item "
    // RDKit❗✔️:             << itemid << " can not be parsed. Ignoring it." << std::endl;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let tokens = split_sdf_property_list_tokens(value);
    let mut missing_value = "n/a";
    let mut first_token = 0;
    if tokens.len() == item_count + 1
        && tokens[0].starts_with('[')
        && tokens[0].ends_with(']')
        && tokens[0].len() >= 2
    {
        missing_value = &tokens[0][1..tokens[0].len() - 1];
        first_token = 1;
    }
    if tokens.len().saturating_sub(first_token) != item_count {
        return None;
    }
    let values = tokens[first_token..]
        .iter()
        .map(|token| {
            if *token == missing_value {
                return None;
            }
            match value_kind {
                SdfPropertyListValueKind::String => Some((*token).to_owned()),
                SdfPropertyListValueKind::Int => {
                    token.parse::<i32>().ok().map(|_| (*token).to_owned())
                }
                SdfPropertyListValueKind::Double => {
                    token.parse::<f64>().ok().map(|_| (*token).to_owned())
                }
                SdfPropertyListValueKind::Bool => parse_sdf_bool_property(token),
            }
        })
        .collect();
    // END RDKIT CPP FUNCTION
    Some(values)
}

fn apply_sdf_property_lists(record: &mut MolBlockRecord, data_fields: &[(String, String)]) {
    // BEGIN RDKIT CPP FUNCTION processMolPropertyLists
    // RDKit✔️✔️: inline void processMolPropertyLists(
    // RDKit✔️✔️:     ROMol &mol, const std::string &missingValueMarker = "n/a") {
    // RDKit✔️✔️:   for (const auto &pn : mol.getPropList()) {
    // RDKit✔️✔️:     processMolPropertyList(mol, pn, missingValueMarker);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for (field_name, field_value) in data_fields {
        let Some((target, property_name, value_kind)) = sdf_property_list_target(field_name) else {
            continue;
        };
        let item_count = match (&*record, target) {
            (MolBlockRecord::Concrete { topology, .. }, SdfPropertyListTarget::Atom) => {
                topology.atoms.len()
            }
            (MolBlockRecord::Concrete { topology, .. }, SdfPropertyListTarget::Bond) => {
                topology.bonds.len()
            }
            (MolBlockRecord::Query(query, ..), SdfPropertyListTarget::Atom) => {
                query.query.num_atoms()
            }
            (MolBlockRecord::Query(query, ..), SdfPropertyListTarget::Bond) => {
                query.query.num_bonds()
            }
        };
        let Some(values) = parse_sdf_property_list_values(field_value, item_count, value_kind)
        else {
            continue;
        };
        match (&mut *record, target) {
            (MolBlockRecord::Concrete { topology, .. }, SdfPropertyListTarget::Atom) => {
                for (atom, value) in topology.atoms.iter_mut().zip(&values) {
                    if let Some(value) = value {
                        atom.set_prop(property_name, value);
                    }
                }
            }
            (MolBlockRecord::Concrete { topology, .. }, SdfPropertyListTarget::Bond) => {
                for (bond, value) in topology.bonds.iter_mut().zip(&values) {
                    if let Some(value) = value {
                        bond.set_prop(property_name, value);
                    }
                }
            }
            (MolBlockRecord::Query(query, ..), SdfPropertyListTarget::Atom) => {
                for (atom, value) in query.query.atoms_mut().iter_mut().zip(&values) {
                    if let Some(value) = value {
                        atom.atom_mut().set_prop(property_name, value);
                    }
                }
            }
            (MolBlockRecord::Query(query, ..), SdfPropertyListTarget::Bond) => {
                for (bond, value) in query.query.bonds_mut().iter_mut().zip(&values) {
                    if let Some(value) = value {
                        bond.bond_mut().set_prop(property_name, value);
                    }
                }
            }
        }
        let property_list = SdfPropertyList::new(target, property_name, values);
        match record {
            MolBlockRecord::Concrete { properties, .. } => {
                *properties = properties.clone().with_sdf_property_list(property_list);
            }
            MolBlockRecord::Query(query) => {
                query.properties = query
                    .properties
                    .clone()
                    .with_sdf_property_list(property_list);
            }
        }
    }
    // END RDKIT CPP FUNCTION
}

fn sdf_data_lines(block: &str) -> Result<Vec<&str>, SdfReadError> {
    let mut byte_offset = 0;
    for line in block.split_inclusive('\n') {
        let content = strip_terminal_cr(line.strip_suffix('\n').unwrap_or(line));
        byte_offset += line.len();
        if content.starts_with("M  END") {
            return Ok(block[byte_offset..].lines().collect());
        }
    }
    Err(SdfReadError::Parse(
        "mol block terminator 'M  END' not found".to_owned(),
    ))
}

fn sdf_record_slices(block: &str) -> Vec<&str> {
    // RDKit✔️❌: while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❌:        (tempStr.empty() || tempStr.at(0) != '$' ||
    // RDKit✔️❌:         tempStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❌:   std::getline(*dp_inStream, tempStr);
    // The source supplier reads incrementally. This detached convenience API
    // keeps borrowed whole-record slices, adding one linear scan but no record
    // string copies.
    let mut records = Vec::new();
    let mut record_start = 0;
    let mut cursor = 0;
    for line in block.split_inclusive('\n') {
        let content = strip_terminal_cr(line.strip_suffix('\n').unwrap_or(line));
        cursor += line.len();
        if is_sdf_record_delimiter(content) {
            let record = &block[record_start..cursor];
            if !strip_sdf_line(record).is_empty() {
                records.push(record);
            }
            record_start = cursor;
        }
    }
    if record_start < block.len() {
        let record = &block[record_start..];
        if !strip_sdf_line(record).is_empty() {
            records.push(record);
        }
    }
    records
}

fn atomic_number_query(number: u8) -> QueryNode<AtomQueryPredicate> {
    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(number))
}

fn complex_molfile_atom_query(symbol: &str) -> Option<QueryNode<AtomQueryPredicate>> {
    // BEGIN RDKIT CPP FUNCTION convertComplexNameToQuery
    // RDKit✔️✔️: if (symb == "Q") {
    // RDKit✔️✔️:   query->setQuery(makeQAtomQuery());
    // RDKit✔️✔️: } else if (symb == "QH") {
    // RDKit✔️✔️:   query->setQuery(makeQHAtomQuery());
    // RDKit✔️✔️: } else if (symb == "A") {
    // RDKit✔️✔️:   query->setQuery(makeAAtomQuery());
    // RDKit✔️✔️: } else if (symb == "AH") {
    // RDKit✔️✔️:   query->setQuery(makeAHAtomQuery());
    // RDKit✔️✔️: } else if (symb == "X") {
    // RDKit✔️✔️:   query->setQuery(makeXAtomQuery());
    // RDKit✔️✔️: } else if (symb == "XH") {
    // RDKit✔️✔️:   query->setQuery(makeXHAtomQuery());
    // RDKit✔️✔️: } else if (symb == "M") {
    // RDKit✔️✔️:   query->setQuery(makeMAtomQuery());
    // RDKit✔️✔️: } else if (symb == "MH") {
    // RDKit✔️✔️:   query->setQuery(makeMHAtomQuery());
    // RDKit✔️✔️: }
    let query = match symbol {
        "Q" => QueryNode::not(QueryNode::or(vec![
            atomic_number_query(6),
            atomic_number_query(1),
        ])),
        "QH" => QueryNode::not(atomic_number_query(6)),
        "A" => QueryNode::not(atomic_number_query(1)),
        "AH" => QueryNode::predicate(AtomQueryPredicate::Any),
        "X" | "XH" => {
            let mut values = vec![9, 17, 35, 53, 85];
            if symbol == "XH" {
                values.push(1);
            }
            QueryNode::or(values.into_iter().map(atomic_number_query).collect())
        }
        "M" | "MH" => {
            let mut excluded = vec![
                0, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86,
            ];
            if symbol == "M" {
                excluded.push(1);
            }
            QueryNode::not(QueryNode::or(
                excluded.into_iter().map(atomic_number_query).collect(),
            ))
        }
        _ => return None,
    };
    Some(query)
    // END RDKIT CPP FUNCTION
}

fn v2000_element(
    symbol: &str,
    strict_parsing: bool,
) -> Result<
    (
        Element,
        Option<u16>,
        Option<QueryNode<AtomQueryPredicate>>,
        Option<String>,
    ),
    SdfReadError,
> {
    if let Some(query) = complex_molfile_atom_query(symbol) {
        return Ok((Element::DUMMY, None, Some(query), None));
    }
    if matches!(symbol, "*" | "R") {
        return Ok((
            Element::DUMMY,
            None,
            Some(QueryNode::predicate(AtomQueryPredicate::Any)),
            None,
        ));
    }
    if let Some(label) = symbol.strip_prefix('R')
        && !label.is_empty()
        && label.len() <= 2
        && label.bytes().all(|byte| byte.is_ascii_digit())
    {
        let label = label.parse::<u32>().unwrap_or(0);
        return Ok((
            Element::DUMMY,
            Some(label as u16),
            Some(QueryNode::predicate(AtomQueryPredicate::RGroupLabel(label))),
            Some(symbol.to_owned()),
        ));
    }
    if matches!(symbol, "L" | "LP" | "R#") {
        return Ok((
            Element::DUMMY,
            None,
            Some(QueryNode::predicate(AtomQueryPredicate::MolFileAlias(
                symbol.to_owned(),
            ))),
            (symbol.starts_with('R')).then(|| symbol.to_owned()),
        ));
    }
    if symbol == "D" {
        return Ok((Element::H, Some(2), None, None));
    }
    if symbol == "T" {
        return Ok((Element::H, Some(3), None, None));
    }
    if matches!(symbol, "Pol" | "Mod") {
        return Ok((Element::DUMMY, None, None, Some(symbol.to_owned())));
    }
    let normalized;
    let lookup = if symbol.len() == 2 && symbol.as_bytes()[1].is_ascii_uppercase() {
        normalized = format!(
            "{}{}",
            &symbol[..1],
            (symbol.as_bytes()[1] as char).to_ascii_lowercase()
        );
        normalized.as_str()
    } else {
        symbol
    };
    match Element::from_symbol(lookup) {
        Some(element) => Ok((element, None, None, None)),
        None if !strict_parsing && !symbol.is_empty() => {
            // BEGIN RDKIT CPP FUNCTION lookupAtomicNumber
            // RDKit✔️✔️:   } catch (const Invar::Invariant &e) {
            // RDKit✔️✔️:     if (strictParsing || symb.empty()) {
            // RDKit✔️✔️:       throw FileParseException(e.what());
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       res->setAtomicNum(0);
            // RDKit✔️✔️:       res->setProp(common_properties::dummyLabel, symb);
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:   }
            // END RDKIT CPP FUNCTION
            Ok((Element::DUMMY, None, None, Some(symbol.to_owned())))
        }
        None => Err(SdfReadError::Parse(format!("Element '{symbol}' not found"))),
    }
}

fn parse_v2000_atom_line(
    line: &str,
    line_number: usize,
    strict_parsing: bool,
) -> Result<ParsedV2000Atom, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseMolFileAtomLine
    // RDKit✔️✔️: if ((strictParsing && text.size() < 34) || text.size() < 32) {
    // RDKit✔️✔️:   std::ostringstream errout;
    // RDKit✔️✔️:   errout << "Atom line too short: '" << text << "' on line " << line;
    // RDKit✔️✔️:   throw FileParseException(errout.str());
    // RDKit✔️✔️: }
    if (strict_parsing && line.len() < 34) || line.len() < 32 {
        return Err(SdfReadError::Parse(format!(
            "Atom line too short: '{line}' on line {line_number}"
        )));
    }

    // RDKit✔️✔️: pos.x = FileParserUtils::toDouble(text.substr(0, 10));
    // RDKit✔️✔️: pos.y = FileParserUtils::toDouble(text.substr(10, 10));
    // RDKit✔️✔️: pos.z = FileParserUtils::toDouble(text.substr(20, 10));
    let coordinate = [
        parse_rdkit_double(rdkit_substr(line, 0, 10)).map_err(|()| {
            SdfReadError::Parse(format!("Cannot process coordinates on line {line_number}"))
        })?,
        parse_rdkit_double(rdkit_substr(line, 10, 10)).map_err(|()| {
            SdfReadError::Parse(format!("Cannot process coordinates on line {line_number}"))
        })?,
        parse_rdkit_double(rdkit_substr(line, 20, 10)).map_err(|()| {
            SdfReadError::Parse(format!("Cannot process coordinates on line {line_number}"))
        })?,
    ];
    // RDKit✔️✔️: symb = text.substr(31, 3);
    // RDKit✔️✔️: boost::trim(symb);
    let symbol = rdkit_substr(line, 31, 3).trim();
    let (element, shorthand_isotope, mut query, dummy_label) =
        v2000_element(symbol, strict_parsing)?;
    let mut spec = AtomSpec::new(element);
    if query.is_some() {
        spec = spec.with_no_implicit(true);
    }
    if let Some(isotope) = shorthand_isotope {
        spec = spec.with_isotope(isotope);
    }
    if let Some(dummy_label) = dummy_label {
        spec = spec.with_prop("dummyLabel", dummy_label)?;
    }

    // RDKit✔️✔️: massDiff = 0;
    // RDKit✔️✔️: if (text.size() >= 36 && text.substr(34, 2) != " 0") {
    // RDKit✔️✔️:   massDiff = FileParserUtils::toInt(text.substr(34, 2), true);
    // RDKit✔️✔️: }
    let mass_diff = if line.len() >= 36 && rdkit_substr(line, 34, 2) != " 0" {
        parse_required_int(line, 34, 2, line_number)?
    } else {
        0
    };
    if mass_diff != 0 {
        let base = cosmolkit_core::most_common_isotope(element) as i32;
        let isotope = base + mass_diff;
        if isotope >= 0 {
            spec = spec.with_isotope(isotope as u16);
        }
    }

    // RDKit✔️✔️: chg = 0;
    // RDKit✔️✔️: if (text.size() >= 39 && text.substr(36, 3) != "  0") {
    // RDKit✔️✔️:   chg = FileParserUtils::toInt(text.substr(36, 3), true);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (chg != 0) {
    // RDKit✔️✔️:   res->setFormalCharge(4 - chg);
    // RDKit✔️✔️: }
    if line.len() >= 39 && rdkit_substr(line, 36, 3) != "  0" {
        let charge_code = parse_required_int(line, 36, 3, line_number)?;
        if charge_code != 0 {
            spec = spec.with_formal_charge((4 - charge_code) as i8);
        }
    }

    // RDKit✔️✔️: if (text.size() >= 42 && text.substr(39, 3) != "  0") {
    // RDKit✔️✔️:   parity = FileParserUtils::toInt(text.substr(39, 3), true);
    // RDKit✔️✔️:   res->setProp(common_properties::molParity, parity);
    // RDKit✔️✔️: }
    if line.len() >= 42 && rdkit_substr(line, 39, 3) != "  0" {
        let parity = parse_required_int(line, 39, 3, line_number)?;
        spec = spec
            .with_mol_parity(parity)
            .with_prop("molParity", parity.to_string())?;
    }

    // RDKit✔️✔️: hCount = 0;
    // RDKit✔️✔️: if (text.size() >= 45 && text.substr(42, 3) != "  0") {
    // RDKit✔️✔️:   hCount = FileParserUtils::toInt(text.substr(42, 3), true);
    // RDKit✔️✔️: }
    if line.len() >= 45 && rdkit_substr(line, 42, 3) != "  0" {
        let h_count = parse_required_int(line, 42, 3, line_number)?;
        if h_count >= 1 {
            spec = spec.with_no_implicit(true);
            let predicate = if h_count > 1 {
                AtomQueryPredicate::ImplicitHydrogenCountLessEqual((h_count - 1) as u8)
            } else {
                AtomQueryPredicate::ImplicitHydrogenCount(0)
            };
            let predicate = QueryNode::predicate(predicate);
            query = Some(match query {
                Some(existing) => QueryNode::and(vec![existing, predicate]),
                None => predicate,
            });
        }
    }

    for (start, key) in [
        (45, "molStereoCare"),
        (48, "molTotValence"),
        (54, "molRxnRole"),
        (57, "molRxnComponent"),
        (66, "molRxnExactChange"),
    ] {
        if line.len() >= start + 3 && rdkit_substr(line, start, 3) != "  0" {
            let value = parse_required_int(line, start, 3, line_number)?;
            if value != 0 || key == "molStereoCare" {
                spec = spec.with_prop(key, value.to_string())?;
            }
        }
    }
    // RDKit✔️✔️: atomMapNumber = FileParserUtils::toInt(text.substr(60, 3), true);
    // RDKit✔️✔️: res->setProp(common_properties::molAtomMapNumber, atomMapNumber);
    if line.len() >= 63 && rdkit_substr(line, 60, 3) != "  0" {
        let atom_map = parse_required_int(line, 60, 3, line_number)?;
        spec = spec
            .with_atom_map(atom_map as u32)
            .with_prop("molAtomMapNumber", atom_map.to_string())?;
    }
    // RDKit✔️✔️: inversionFlag = FileParserUtils::toInt(text.substr(63, 3), true);
    // RDKit✔️✔️: res->setProp(common_properties::molInversionFlag, inversionFlag);
    if line.len() >= 66 && rdkit_substr(line, 63, 3) != "  0" {
        let inversion = parse_required_int(line, 63, 3, line_number)?;
        spec = spec.with_mol_inversion_flag(inversion);
    }
    // END RDKIT CPP FUNCTION
    Ok(ParsedV2000Atom {
        spec,
        query,
        coordinate,
    })
}

fn parse_v2000_bond_line(
    line: &str,
    line_number: usize,
    atom_count: usize,
    bond_id: usize,
) -> Result<ParsedV2000Bond, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseMolFileBondLine
    // RDKit✔️✔️: if (text.size() < 9) {
    // RDKit✔️✔️:   errout << "Bond line too short: '" << text << "' on line " << line;
    // RDKit✔️✔️:   throw FileParseException(errout.str());
    // RDKit✔️✔️: }
    if line.len() < 9 {
        return Err(SdfReadError::Parse(format!(
            "Bond line too short: '{line}' on line {line_number}"
        )));
    }
    // RDKit✔️✔️: idx1 = FileParserUtils::toUnsigned(text.substr(spos, 3));
    // RDKit✔️✔️: idx2 = FileParserUtils::toUnsigned(text.substr(spos, 3));
    // RDKit✔️✔️: bType = FileParserUtils::toUnsigned(text.substr(spos, 3));
    // RDKit✔️✔️: idx1--;
    // RDKit✔️✔️: idx2--;
    let begin = parse_required_unsigned(line, 0, 3, line_number)? as usize;
    let end = parse_required_unsigned(line, 3, 3, line_number)? as usize;
    if begin == 0 || end == 0 || begin > atom_count || end > atom_count {
        return Err(SdfReadError::Field {
            kind: "bond endpoint",
            line: line_number,
            value: line.to_owned(),
        });
    }
    let bond_type = parse_required_unsigned(line, 6, 3, line_number)?;
    // RDKit✔️✔️: switch (bType) {
    // RDKit✔️✔️:   case 1:
    // RDKit✔️✔️:     type = Bond::SINGLE;
    // RDKit✔️✔️:     res = new Bond;
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   case 2:
    // RDKit✔️✔️:     type = Bond::DOUBLE;
    // RDKit✔️✔️:     res = new Bond;
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   case 3:
    // RDKit✔️✔️:     type = Bond::TRIPLE;
    // RDKit✔️✔️:     res = new Bond;
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   case 4:
    // RDKit✔️✔️:     type = Bond::AROMATIC;
    // RDKit✔️✔️:     res = new Bond;
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   case 9:
    // RDKit✔️✔️:     type = Bond::DATIVE;
    // RDKit✔️✔️:     res = new Bond;
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   case 0:
    // RDKit✔️✔️:     type = Bond::UNSPECIFIED;
    // RDKit✔️✔️:     res = new Bond;
    // RDKit✔️✔️:   default:
    // RDKit✔️✔️:     type = Bond::UNSPECIFIED;
    // RDKit✔️✔️:     // it's a query bond of some type
    // RDKit✔️✔️:     res = new QueryBond;
    // RDKit✔️✔️:     if (bType == 8) {
    // RDKit✔️✔️:       BOND_NULL_QUERY *q;
    // RDKit✔️✔️:       q = makeBondNullQuery();
    // RDKit✔️✔️:       res->setQuery(q);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       BOND_NULL_QUERY *q;
    // RDKit✔️✔️:       q = makeBondNullQuery();
    // RDKit✔️✔️:       res->setQuery(q);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     break;
    let (order, mut query) = match bond_type {
        0 => (BondOrder::Unspecified, None),
        1 => (BondOrder::Single, None),
        2 => (BondOrder::Double, None),
        3 => (BondOrder::Triple, None),
        4 => (BondOrder::Aromatic, None),
        9 => (BondOrder::Dative, None),
        5 => (
            BondOrder::Unspecified,
            Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ]))),
        ),
        6 => (
            BondOrder::Unspecified,
            Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ]))),
        ),
        7 => (
            BondOrder::Unspecified,
            Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Double,
                BondOrder::Aromatic,
            ]))),
        ),
        8 => (
            BondOrder::Unspecified,
            Some(QueryNode::predicate(BondQueryPredicate::Any)),
        ),
        _ => (
            BondOrder::Unspecified,
            Some(QueryNode::predicate(BondQueryPredicate::Any)),
        ),
    };
    let mut spec = BondSpec::new(AtomId::new(begin - 1), AtomId::new(end - 1), order)
        .with_prop("_MolFileBondType", bond_type.to_string())?;
    if order == BondOrder::Aromatic {
        spec = spec.with_aromatic(true);
    }
    if matches!(bond_type, 5..=7) {
        spec = spec.with_prop("_MolFileBondQuery", "1")?;
    }

    // RDKit✔️✔️: if (text.size() >= 12 && text.substr(9, 3) != "  0") {
    // RDKit✔️✔️:   stereo = FileParserUtils::toUnsigned(text.substr(9, 3));
    // RDKit✔️✔️:   switch (stereo) {
    // RDKit✔️✔️:     case 0:
    // RDKit✔️✔️:       res->setBondDir(Bond::NONE);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 1:
    // RDKit✔️✔️:       res->setBondDir(Bond::BEGINWEDGE);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 6:
    // RDKit✔️✔️:       res->setBondDir(Bond::BEGINDASH);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 3:  // "either" double bond
    // RDKit✔️✔️:       res->setBondDir(Bond::EITHERDOUBLE);
    // RDKit✔️✔️:       res->setStereo(Bond::STEREOANY);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 4:  // "either" single bond
    // RDKit✔️✔️:       res->setBondDir(Bond::UNKNOWN);
    // RDKit✔️✔️:       break;
    if line.len() >= 12
        && rdkit_substr(line, 9, 3) != "  0"
        && let Ok(stereo) = parse_rdkit_unsigned(rdkit_substr(line, 9, 3))
    {
        spec = spec.with_prop("_MolFileBondStereo", stereo.to_string())?;
        spec = match stereo {
            0 => spec.with_direction(BondDirection::None),
            1 => spec.with_direction(BondDirection::BeginWedge),
            6 => spec.with_direction(BondDirection::BeginDash),
            3 => spec
                .with_direction(BondDirection::EitherDouble)
                .with_stereo(BondStereo::Any),
            4 => spec.with_direction(BondDirection::Unknown),
            _ => spec,
        };
    }
    if line.len() >= 21
        && rdkit_substr(line, 18, 3) != "  0"
        && let Ok(status) = parse_rdkit_int(rdkit_substr(line, 18, 3))
    {
        spec = spec.with_prop("molReactStatus", status.to_string())?;
    }
    if line.len() >= 18
        && rdkit_substr(line, 15, 3) != "  0"
        && let Ok(topology) = parse_rdkit_int(rdkit_substr(line, 15, 3))
        && topology != 0
    {
        let topology_query = match topology {
            1 => QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            2 => QueryNode::predicate(BondQueryPredicate::IsInRing(false)),
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Unrecognized bond topology specifier: {topology} on line {line_number}"
                )));
            }
        };
        query = Some(match query {
            Some(existing) => QueryNode::and(vec![existing, topology_query]),
            None => topology_query,
        });
    }
    // END RDKIT CPP FUNCTION
    Ok(ParsedV2000Bond {
        bond: Bond::from_spec(BondId::new(bond_id), spec),
        query,
    })
}

fn v2000_atom_mut(
    atoms: &mut [ParsedV2000Atom],
    one_based_id: i32,
    line_number: usize,
) -> Result<&mut ParsedV2000Atom, SdfReadError> {
    let index = one_based_id.checked_sub(1).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Atom index {one_based_id} out of range on line {line_number}"
        ))
    })? as usize;
    atoms.get_mut(index).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Atom index {one_based_id} out of range on line {line_number}"
        ))
    })
}

fn merge_v2000_atom_query(
    existing: Option<QueryNode<AtomQueryPredicate>>,
    next: QueryNode<AtomQueryPredicate>,
) -> QueryNode<AtomQueryPredicate> {
    match existing {
        Some(existing) => QueryNode::and(vec![existing, next]),
        None => next,
    }
}

fn v2000_list_atomic_number(symbol: &str, line_number: usize) -> Result<u8, SdfReadError> {
    let normalized;
    let lookup = if symbol.len() == 2 && symbol.as_bytes()[1].is_ascii_uppercase() {
        normalized = format!(
            "{}{}",
            &symbol[..1],
            (symbol.as_bytes()[1] as char).to_ascii_lowercase()
        );
        normalized.as_str()
    } else {
        symbol
    };
    Element::from_symbol(lookup)
        .map(Element::atomic_number)
        .ok_or_else(|| {
            SdfReadError::Parse(format!(
                "Element '{symbol}' not found in atom list on line {line_number}"
            ))
        })
}

fn non_atomic_and_components(
    query: QueryNode<AtomQueryPredicate>,
) -> Option<Vec<QueryNode<AtomQueryPredicate>>> {
    match query {
        QueryNode::And(children) => {
            let mut result = Vec::new();
            for child in children {
                result.extend(non_atomic_and_components(child)?);
            }
            Some(result)
        }
        QueryNode::Predicate(
            AtomQueryPredicate::AtomicNumber(_)
            | AtomQueryPredicate::AtomicNumberIn(_)
            | AtomQueryPredicate::AtomicNumberNotIn(_),
        ) => Some(Vec::new()),
        QueryNode::Predicate(_) => Some(vec![query]),
        QueryNode::Or(_) | QueryNode::Xor(_) | QueryNode::Not(_) => None,
    }
}

fn replace_v2000_atom_list_query(
    existing: Option<QueryNode<AtomQueryPredicate>>,
    list_query: QueryNode<AtomQueryPredicate>,
) -> QueryNode<AtomQueryPredicate> {
    let Some(mut preserved) = existing.and_then(non_atomic_and_components) else {
        return list_query;
    };
    if preserved.is_empty() {
        list_query
    } else {
        preserved.insert(0, list_query);
        QueryNode::and(preserved)
    }
}

fn parse_v2000_atom_list_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseNewAtomList
    // RDKit❗✔️: if (text.size() < 15) {
    // RDKit❗✔️:   std::ostringstream errout;
    // RDKit❗✔️:   errout << "Atom list line too short: '" << text << "'";
    // RDKit❗✔️:   throw FileParseException(errout.str());
    // RDKit❗✔️: }
    // RDKit❗✔️:     idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(7, 3)) -
    // RDKit❗✔️:           1;
    // RDKit❗✔️: nQueries = FileParserUtils::toInt(text.substr(10, 3));
    // RDKit❗✔️: if (!nQueries) {
    // RDKit❗✔️:   return;
    // RDKit❗✔️: }
    // RDKit❗✔️: if (nQueries < 0) {
    // RDKit❗✔️:   throw FileParseException(errout.str());
    // RDKit❗✔️: }
    // RDKit❗✔️: for (unsigned int i = 0; i < static_cast<unsigned int>(nQueries); i++) {
    // RDKit❗✔️:   unsigned int pos = 16 + i * 4;
    // RDKit❗✔️:   if (text.size() < pos + 4) {
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::string atSymb = text.substr(pos, 4);
    // RDKit❗✔️:   atSymb.erase(atSymb.find(' '), atSymb.size());
    // RDKit❗✔️:   int atNum = PeriodicTable::getTable()->getAtomicNumber(atSymb);
    // RDKit❗✔️:   if (!i) {
    // RDKit❗✔️:     a->setAtomicNum(atNum);
    // RDKit❗✔️:     a->setQuery(makeAtomNumQuery(atNum));
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     a->expandQuery(makeAtomNumQuery(atNum), Queries::COMPOSITE_OR, true);
    // RDKit❗✔️:     a->setAtomicNum(0);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: a->setProp(common_properties::_MolFileAtomQuery, 1);
    // RDKit❗✔️: switch (text[14]) {
    // RDKit❗✔️:   case 'T':
    // RDKit❗✔️:     a->getQuery()->setNegation(true);
    // RDKit❗✔️:     break;
    // RDKit❗✔️:   case 'F':
    // RDKit❗✔️:     a->getQuery()->setNegation(false);
    // RDKit❗✔️:     break;
    // RDKit❗✔️: }
    if line.len() < 15 {
        return Err(SdfReadError::Parse(format!(
            "Atom list line too short: '{line}'"
        )));
    }
    let atom_id = parse_required_unsigned(line, 7, 3, line_number)?;
    let query_count = parse_required_int(line, 10, 3, line_number)?;
    if query_count == 0 {
        return Ok(());
    }
    if query_count < 0 {
        return Err(SdfReadError::Parse(format!(
            "negative length atom list: '{line}' on line {line_number}."
        )));
    }
    let mut atomic_numbers = Vec::with_capacity(query_count as usize);
    for index in 0..query_count as usize {
        let position = 16 + index * 4;
        if line.len() < position + 4 {
            return Err(SdfReadError::Parse(format!(
                "Atom list line too short: '{line}' on line {line_number}"
            )));
        }
        let symbol = rdkit_substr(line, position, 4)
            .split_once(' ')
            .map_or(rdkit_substr(line, position, 4), |(symbol, _)| symbol);
        atomic_numbers.push(v2000_list_atomic_number(symbol, line_number)?);
    }
    let modifier = line.as_bytes()[14];
    let predicate = match modifier {
        b'T' => AtomQueryPredicate::AtomicNumberNotIn(atomic_numbers),
        b'F' => AtomQueryPredicate::AtomicNumberIn(atomic_numbers),
        other => {
            return Err(SdfReadError::Parse(format!(
                "Unrecognized atom-list query modifier: '{}' on line {line_number}",
                other as char
            )));
        }
    };
    let element = match &predicate {
        AtomQueryPredicate::AtomicNumberIn(numbers) if numbers.len() == 1 => {
            Element::from_atomic_number(numbers[0]).unwrap_or(Element::DUMMY)
        }
        _ => Element::DUMMY,
    };
    let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
    update_v2000_atom_spec(atom, |spec| {
        spec.with_element(element)
            .with_no_implicit(true)
            .with_prop("_MolFileAtomQuery", "1")
    })?;
    atom.query = Some(replace_v2000_atom_list_query(
        atom.query.take(),
        QueryNode::predicate(predicate),
    ));
    // END RDKIT CPP FUNCTION
    Ok(())
}

fn parse_v2000_old_atom_list_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseOldAtomList
    // RDKit❗✔️: void ParseOldAtomList(RWMol *mol, const std::string_view &text,
    // RDKit❗✔️:                       unsigned int line) {
    // RDKit❗✔️:   PRECONDITION(mol, "bad mol");
    // RDKit❗✔️:   unsigned int idx;
    // RDKit❗✔️:   try {
    // RDKit❗✔️:     idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(0, 3)) -
    // RDKit❗✔️:           1;
    // RDKit❗✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit❗✔️:     std::ostringstream errout;
    // RDKit❗✔️:     errout << "Cannot convert '" << text.substr(0, 3) << "' to int on line "
    // RDKit❗✔️:            << line;
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   URANGE_CHECK(idx, mol->getNumAtoms());
    // RDKit❗✔️:   QueryAtom a(*(mol->getAtomWithIdx(idx)));
    // RDKit❗✔️:
    // RDKit❗✔️:   auto *q = new ATOM_OR_QUERY;
    // RDKit❗✔️:   q->setDescription("AtomOr");
    // RDKit❗✔️:
    // RDKit❗✔️:   switch (text[4]) {
    // RDKit❗✔️:     case 'T':
    // RDKit❗✔️:       q->setNegation(true);
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 'F':
    // RDKit❗✔️:       q->setNegation(false);
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       delete q;
    // RDKit❗✔️:       std::ostringstream errout;
    // RDKit❗✔️:       errout << "Unrecognized atom-list query modifier: '" << text[4]
    // RDKit❗✔️:              << "' on line " << line;
    // RDKit❗✔️:       throw FileParseException(errout.str());
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   int nQueries;
    // RDKit❗✔️:   try {
    // RDKit❗✔️:     nQueries = FileParserUtils::toInt(text.substr(9, 1));
    // RDKit❗✔️:   } catch (const std::out_of_range &) {
    // RDKit❗✔️:     delete q;
    // RDKit❗✔️:     std::ostringstream errout;
    // RDKit❗✔️:     errout << "Cannot convert position 9 of '" << text << "' to int on line "
    // RDKit❗✔️:            << line;
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit❗✔️:     delete q;
    // RDKit❗✔️:     std::ostringstream errout;
    // RDKit❗✔️:     errout << "Cannot convert '" << text.substr(9, 1) << "' to int on line "
    // RDKit❗✔️:            << line;
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   RANGE_CHECK(0, nQueries, 5);
    // RDKit❗✔️:   for (int i = 0; i < nQueries; i++) {
    // RDKit❗✔️:     int pos = 11 + i * 4;
    // RDKit❗✔️:     int atNum;
    // RDKit❗✔️:     try {
    // RDKit❗✔️:       atNum = FileParserUtils::toInt(text.substr(pos, 3));
    // RDKit❗✔️:     } catch (const std::out_of_range &) {
    // RDKit❗✔️:       delete q;
    // RDKit❗✔️:       std::ostringstream errout;
    // RDKit❗✔️:       errout << "Cannot convert position " << pos << " of '" << text
    // RDKit❗✔️:              << "' to int on line " << line;
    // RDKit❗✔️:       throw FileParseException(errout.str());
    // RDKit❗✔️:     } catch (boost::bad_lexical_cast &) {
    // RDKit❗✔️:       delete q;
    // RDKit❗✔️:       std::ostringstream errout;
    // RDKit❗✔️:       errout << "Cannot convert '" << text.substr(pos, 3) << "' to int on line "
    // RDKit❗✔️:              << line;
    // RDKit❗✔️:       throw FileParseException(errout.str());
    // RDKit❗✔️:     }
    // RDKit❗✔️:     RANGE_CHECK(0, atNum, 200);  // goofy!
    // RDKit❗✔️:     q->addChild(
    // RDKit❗✔️:         QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(atNum)));
    // RDKit❗✔️:     if (!i) {
    // RDKit❗✔️:       a.setAtomicNum(atNum);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   a.setQuery(q);
    // RDKit❗✔️:   a.setProp(common_properties::_MolFileAtomQuery, 1);
    // RDKit❗✔️:
    // RDKit❗✔️:   mol->replaceAtom(idx, &a);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION
    let atom_id = parse_required_unsigned(line, 0, 3, line_number)?;
    let modifier = line
        .as_bytes()
        .get(4)
        .copied()
        .ok_or_else(|| SdfReadError::Parse(format!("Atom list line too short: '{line}'")))?;
    if !matches!(modifier, b'T' | b'F') {
        return Err(SdfReadError::Parse(format!(
            "Unrecognized atom-list query modifier: '{}' on line {line_number}",
            modifier as char
        )));
    }
    let query_count = parse_required_int(line, 9, 1, line_number)?;
    if !(0..=5).contains(&query_count) {
        return Err(SdfReadError::Parse(format!(
            "atom-list query count {query_count} out of range on line {line_number}"
        )));
    }
    let mut atomic_numbers = Vec::with_capacity(query_count as usize);
    for index in 0..query_count as usize {
        let position = 11 + index * 4;
        let atomic_number = parse_required_int(line, position, 3, line_number)?;
        if !(0..=200).contains(&atomic_number) {
            return Err(SdfReadError::Parse(format!(
                "atom-list atomic number {atomic_number} out of range on line {line_number}"
            )));
        }
        atomic_numbers.push(atomic_number as u8);
    }
    let element = atomic_numbers
        .first()
        .copied()
        .map(|atomic_number| {
            Element::from_atomic_number(atomic_number).ok_or(SdfReadError::Unsupported(
                "old-style V2000 atom-list element outside the detached periodic table",
            ))
        })
        .transpose()?;
    let predicate = if modifier == b'T' {
        AtomQueryPredicate::AtomicNumberNotIn(atomic_numbers)
    } else {
        AtomQueryPredicate::AtomicNumberIn(atomic_numbers)
    };
    let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
    update_v2000_atom_spec(atom, |mut spec| {
        if let Some(element) = element {
            spec = spec.with_element(element);
        }
        spec.with_no_implicit(true)
            .with_prop("_MolFileAtomQuery", "1")
    })?;
    atom.query = Some(QueryNode::predicate(predicate));
    Ok(())
}

fn parse_v2000_rgroup_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseRGroupLabels
    // RDKit❗✔️: nLabels = FileParserUtils::toInt(text.substr(6, 3));
    // RDKit❗✔️: for (int i = 0; i < nLabels; i++) {
    // RDKit❗✔️:   int pos = 10 + i * 8;
    // RDKit❗✔️:     atIdx = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit❗✔️:         text.substr(pos, 3));
    // RDKit❗✔️:     rLabel = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit❗✔️:         text.substr(pos + 4, 3));
    // RDKit❗✔️:   atIdx -= 1;
    // RDKit❗✔️:   QueryAtom qatom(*(mol->getAtomWithIdx(atIdx)));
    // RDKit❗✔️:   qatom.setProp(common_properties::_MolFileRLabel, rLabel);
    // RDKit❗✔️:   std::string dLabel = "R" + std::to_string(rLabel);
    // RDKit❗✔️:   qatom.setProp(common_properties::dummyLabel, dLabel);
    // RDKit❗✔️:   if (rLabel > 0 && rLabel < 999) {
    // RDKit❗✔️:     qatom.setIsotope(rLabel);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   qatom.setQuery(makeAtomNullQuery());
    // RDKit❗✔️:   mol->replaceAtom(atIdx, &qatom);
    // RDKit❗✔️: }
    let label_count = parse_required_int(line, 6, 3, line_number)?;
    for index in 0..label_count {
        let position = 10 + index as usize * 8;
        let atom_id = parse_required_unsigned(line, position, 3, line_number)?;
        let label = parse_required_unsigned(line, position + 4, 3, line_number)?;
        let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
        update_v2000_atom_spec(atom, |spec| {
            spec.with_prop("_MolFileRLabel", label.to_string())
                .and_then(|spec| spec.with_prop("dummyLabel", format!("R{label}")))
        })?;
        if (1..999).contains(&label) {
            update_v2000_atom_spec(atom, |spec| Ok(spec.with_isotope(label as u16)))?;
        }
        atom.query = Some(QueryNode::predicate(AtomQueryPredicate::Any));
    }
    // END RDKIT CPP FUNCTION
    Ok(())
}

fn parse_v2000_query_count_entries(
    line: &str,
    line_number: usize,
) -> Result<Vec<(u32, i32)>, SdfReadError> {
    let entry_count = parse_required_unsigned(line, 6, 3, line_number)?;
    let mut position = 9;
    let mut entries = Vec::with_capacity(entry_count as usize);
    for _ in 0..entry_count {
        let atom_id = parse_required_unsigned(line, position, 4, line_number)?;
        position += 4;
        let count = if line.len() >= position + 4 && rdkit_substr(line, position, 4) != "    " {
            parse_required_int(line, position, 4, line_number)?
        } else {
            0
        };
        position += 4;
        entries.push((atom_id, count));
    }
    Ok(entries)
}

fn parse_v2000_substitution_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
    explicit_degrees: &[u8],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseSubstitutionCountLine
    // RDKit❗✔️: nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
    // RDKit❗✔️: unsigned int spos = 9;
    // RDKit❗✔️: for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit❗✔️:   aid = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(spos, 4));
    // RDKit❗✔️:   spos += 4;
    // RDKit❗✔️:   Atom *atom = mol->getAtomWithIdx(aid - 1);
    // RDKit❗✔️:   count = FileParserUtils::toInt(text.substr(spos, 4));
    // RDKit❗✔️:   spos += 4;
    // RDKit❗✔️:   if (count == 0) {
    // RDKit❗✔️:     continue;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ATOM_EQUALS_QUERY *q = makeAtomExplicitDegreeQuery(0);
    // RDKit❗✔️:   switch (count) {
    // RDKit❗✔️:     case -1:
    // RDKit❗✔️:       q->setVal(0);
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case -2:
    // RDKit❗✔️:       q->setVal(atom->getDegree());
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 1:
    // RDKit❗✔️:     case 2:
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:     case 4:
    // RDKit❗✔️:     case 5:
    // RDKit❗✔️:       q->setVal(count);
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 6:
    // RDKit❗✔️:       q->setVal(6);
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       throw FileParseException(errout.str());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   atom->expandQuery(q, Queries::COMPOSITE_AND);
    // RDKit❗✔️: }
    for (atom_id, count) in parse_v2000_query_count_entries(line, line_number)? {
        if count == 0 {
            continue;
        }
        let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
        let degree = match count {
            -1 => 0,
            -2 => explicit_degrees[atom_id as usize - 1],
            1..=6 => count as u8,
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Value {count} is not supported as a degree query. line: {line_number}"
                )));
            }
        };
        update_v2000_atom_spec(atom, |spec| {
            spec.with_prop("molSubstCount", count.to_string())
        })?;
        atom.query = Some(merge_v2000_atom_query(
            atom.query.take(),
            QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(degree)),
        ));
    }
    // END RDKIT CPP FUNCTION
    Ok(())
}

fn parse_v2000_unsaturation_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseUnsaturationLine
    // RDKit❗✔️: if (count == 0) {
    // RDKit❗✔️:   continue;
    // RDKit❗✔️: } else if (count == 1) {
    // RDKit❗✔️:   ATOM_EQUALS_QUERY *q = makeAtomUnsaturatedQuery();
    // RDKit❗✔️:   atom->expandQuery(q, Queries::COMPOSITE_AND);
    // RDKit❗✔️: } else {
    // RDKit❗✔️:   throw FileParseException(errout.str());
    // RDKit❗✔️: }
    for (atom_id, count) in parse_v2000_query_count_entries(line, line_number)? {
        if count == 0 {
            continue;
        }
        if count != 1 {
            return Err(SdfReadError::Parse(format!(
                "Value {count} is not supported as an unsaturation query (only 0 and 1 are allowed). line: {line_number}"
            )));
        }
        let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
        atom.query = Some(merge_v2000_atom_query(
            atom.query.take(),
            QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
        ));
    }
    // END RDKIT CPP FUNCTION
    Ok(())
}

fn parse_v2000_ring_bond_count_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
    needs_query_scan: &mut bool,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseRingBondCountLine
    // RDKit❗✔️: ATOM_EQUALS_QUERY *q = makeAtomRingBondCountQuery(0);
    // RDKit❗✔️: switch (count) {
    // RDKit❗✔️:   case -1:
    // RDKit❗✔️:     q->setVal(0);
    // RDKit❗✔️:     break;
    // RDKit❗✔️:   case -2:
    // RDKit❗✔️:     q->setVal(0xDEADBEEF);
    // RDKit❗✔️:     mol->setProp(common_properties::_NeedsQueryScan, 1);
    // RDKit❗✔️:     break;
    // RDKit❗✔️:   case 1:
    // RDKit❗✔️:   case 2:
    // RDKit❗✔️:   case 3:
    // RDKit❗✔️:     q->setVal(count);
    // RDKit❗✔️:     break;
    // RDKit❗✔️:   case 4:
    // RDKit❗✔️:     q = static_cast<ATOM_EQUALS_QUERY *>(new ATOM_LESSEQUAL_QUERY);
    // RDKit❗✔️:     q->setVal(4);
    // RDKit❗✔️:     q->setDescription("AtomRingBondCount");
    // RDKit❗✔️:     q->setDataFunc(queryAtomRingBondCount);
    // RDKit❗✔️:     break;
    // RDKit❗✔️:   default:
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️: }
    // RDKit❗✔️: atom->expandQuery(q, Queries::COMPOSITE_AND);
    for (atom_id, count) in parse_v2000_query_count_entries(line, line_number)? {
        if count == 0 {
            continue;
        }
        let predicate = match count {
            -1 => AtomQueryPredicate::RingBondCount(0),
            -2 => {
                *needs_query_scan = true;
                AtomQueryPredicate::RingBondCount(0xDEAD_BEEF)
            }
            1..=3 => AtomQueryPredicate::RingBondCount(count as u32),
            4 => AtomQueryPredicate::RingBondCountLessEqual(4),
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Value {count} is not supported as a ring-bond count query. line: {line_number}"
                )));
            }
        };
        let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
        update_v2000_atom_spec(atom, |spec| {
            spec.with_prop("molRingBondCount", count.to_string())
        })?;
        atom.query = Some(merge_v2000_atom_query(
            atom.query.take(),
            QueryNode::predicate(predicate),
        ));
    }
    // END RDKIT CPP FUNCTION
    Ok(())
}

fn parse_v2000_pxa_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParsePXALine
    // RDKit✔️✔️:   unsigned int pos = 7;
    // RDKit✔️✔️:     auto atIdx =
    // RDKit✔️✔️:         FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(pos, 3));
    // RDKit✔️✔️:     pos += 3;
    // RDKit✔️✔️:     mol->getAtomWithIdx(atIdx - 1)->setProp(
    // RDKit✔️✔️:         "_MolFile_PXA", text.substr(pos, text.length() - pos));
    // END RDKIT CPP FUNCTION
    let mut position = 7;
    let atom_id = parse_required_unsigned(line, position, 3, line_number)?;
    position += 3;
    let value = line.get(position..).unwrap_or_default();
    let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
    update_v2000_atom_spec(atom, |spec| spec.with_prop("_MolFile_PXA", value))?;
    Ok(())
}

fn parse_v2000_zch_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseZCHLine
    // RDKit✔️✔️:   unsigned int spos = 9;
    // RDKit✔️✔️:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️✔️:       aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️✔️:           text.substr(spos, 4));
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️✔️:         val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos, 4));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!aid || aid > mol->getNumAtoms()) {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       --aid;
    // RDKit✔️✔️:       Atom *atom = mol->getAtomWithIdx(aid);
    // RDKit✔️✔️:         atom->setFormalCharge(val);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    let count = parse_required_unsigned(line, 6, 3, line_number)?;
    let mut position = 9;
    for _ in 0..count {
        let atom_id = parse_required_unsigned(line, position, 4, line_number)?;
        position += 4;
        let charge = if line.len() >= position + 4 && rdkit_substr(line, position, 4) != "    " {
            parse_required_int(line, position, 4, line_number)?
        } else {
            0
        };
        position += 4;
        let charge = i8::try_from(charge).map_err(|_| {
            SdfReadError::Unsupported("V2000 ZCH charge outside the detached i8 charge model")
        })?;
        let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
        update_v2000_atom_spec(atom, |spec| Ok(spec.with_formal_charge(charge)))?;
    }
    Ok(())
}

fn parse_v2000_hyd_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseHYDLine
    // RDKit✔️✔️:   unsigned int spos = 9;
    // RDKit✔️✔️:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️✔️:       aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️✔️:           text.substr(spos, 4));
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️✔️:         val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos, 4));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!aid || aid > mol->getNumAtoms()) {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       --aid;
    // RDKit✔️✔️:       Atom *atom = mol->getAtomWithIdx(aid);
    // RDKit✔️✔️:         if (val >= 0) {
    // RDKit✔️✔️:           atom->setProp("_ZBO_H", true);
    // RDKit✔️✔️:           atom->setNumExplicitHs(val);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    let count = parse_required_unsigned(line, 6, 3, line_number)?;
    let mut position = 9;
    for _ in 0..count {
        let atom_id = parse_required_unsigned(line, position, 4, line_number)?;
        position += 4;
        let explicit_hydrogens =
            if line.len() >= position + 4 && rdkit_substr(line, position, 4) != "    " {
                parse_required_int(line, position, 4, line_number)?
            } else {
                -1
            };
        position += 4;
        let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
        if explicit_hydrogens >= 0 {
            let explicit_hydrogens = u8::try_from(explicit_hydrogens).map_err(|_| {
                SdfReadError::Unsupported(
                    "V2000 HYD count outside the detached u8 explicit-hydrogen model",
                )
            })?;
            update_v2000_atom_spec(atom, |spec| {
                Ok(spec
                    .with_prop("_ZBO_H", "1")?
                    .with_explicit_hydrogens(explicit_hydrogens))
            })?;
        }
    }
    Ok(())
}

fn parse_v2000_zbo_line(
    line: &str,
    line_number: usize,
    bonds: &mut [ParsedV2000Bond],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseZBOLine
    // RDKit✔️✔️:   unsigned int spos = 9;
    // RDKit✔️✔️:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️✔️:       bid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️✔️:           text.substr(spos, 4));
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️✔️:         order = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️✔️:             text.substr(spos, 4));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!bid || bid > mol->getNumBonds()) {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       --bid;
    // RDKit✔️✔️:       Bond *bnd = mol->getBondWithIdx(bid);
    // RDKit✔️✔️:       if (!bnd) {
    // RDKit✔️✔️:         std::ostringstream errout;
    // RDKit✔️✔️:         errout << "Bond " << bid << " from ZBO specification on line " << line
    // RDKit✔️✔️:                << " not found";
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (order == 0) {
    // RDKit✔️✔️:           bnd->setBondType(Bond::ZERO);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           bnd->setBondType(static_cast<Bond::BondType>(order));
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    let count = parse_required_unsigned(line, 6, 3, line_number)?;
    let mut position = 9;
    for _ in 0..count {
        let bond_id = parse_required_unsigned(line, position, 4, line_number)?;
        position += 4;
        let order = if line.len() >= position + 4 && rdkit_substr(line, position, 4) != "    " {
            parse_required_unsigned(line, position, 4, line_number)?
        } else {
            0
        };
        position += 4;
        let bond = bond_id
            .checked_sub(1)
            .and_then(|index| bonds.get_mut(index as usize))
            .ok_or_else(|| {
                SdfReadError::Parse(format!("Bad ZBO specification on line {line_number}"))
            })?;
        let order = match order {
            0 | 21 => BondOrder::Zero,
            1 => BondOrder::Single,
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Quadruple,
            5 => BondOrder::Quintuple,
            6 => BondOrder::Hextuple,
            7 => BondOrder::OneAndHalf,
            8 => BondOrder::TwoAndHalf,
            9 => BondOrder::ThreeAndHalf,
            10 => BondOrder::FourAndHalf,
            11 => BondOrder::FiveAndHalf,
            12 => BondOrder::Aromatic,
            13 => BondOrder::Ionic,
            14 => BondOrder::Hydrogen,
            15 => BondOrder::ThreeCenter,
            16 => BondOrder::DativeOne,
            17 => BondOrder::Dative,
            18 => BondOrder::DativeLeft,
            19 => BondOrder::DativeRight,
            20 => BondOrder::Other,
            _ => {
                return Err(SdfReadError::Unsupported(
                    "V2000 ZBO bond type outside the detached BondOrder model",
                ));
            }
        };
        bond.bond.set_order(order);
    }
    Ok(())
}

fn parse_v2000_atom_alias(
    line: &str,
    value: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseAtomAlias
    // RDKit✔️✔️:     idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(3, 3)) -
    // RDKit✔️✔️:           1;
    // RDKit✔️✔️:   URANGE_CHECK(idx, mol->getNumAtoms());
    // RDKit✔️✔️:   Atom *at = mol->getAtomWithIdx(idx);
    // RDKit✔️✔️:   at->setProp(common_properties::molFileAlias, nextLine);
    // END RDKIT CPP FUNCTION
    let atom_id = parse_required_unsigned(line, 3, 3, line_number)?;
    let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
    update_v2000_atom_spec(atom, |spec| spec.with_prop("molFileAlias", value))?;
    Ok(())
}

fn parse_v2000_atom_value(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseAtomValue
    // RDKit✔️✔️:     idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(3, 3)) -
    // RDKit✔️✔️:           1;
    // RDKit✔️✔️:   URANGE_CHECK(idx, mol->getNumAtoms());
    // RDKit✔️✔️:   Atom *at = mol->getAtomWithIdx(idx);
    // RDKit✔️✔️:   at->setProp(common_properties::molFileValue,
    // RDKit✔️✔️:               text.substr(7, text.length() - 7));
    // END RDKIT CPP FUNCTION
    let atom_id = parse_required_unsigned(line, 3, 3, line_number)?;
    let value = line.get(7..).unwrap_or_default();
    let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
    update_v2000_atom_spec(atom, |spec| spec.with_prop("molFileValue", value))?;
    Ok(())
}

fn parse_v2000_apo_line(
    line: &str,
    line_number: usize,
    atoms: &mut [ParsedV2000Atom],
    params: MolBlockReadParams,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseAttachPointLine
    // RDKit✔️✔️:   unsigned int spos = 9;
    // RDKit✔️✔️:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️✔️:       aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️✔️:           text.substr(spos, 4));
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️✔️:         val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos, 4));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!aid || aid > mol->getNumAtoms()) {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       --aid;
    // RDKit✔️✔️:       Atom *atom = mol->getAtomWithIdx(aid);
    // RDKit✔️✔️:         if (val < 0 || val > 3) {
    // RDKit✔️✔️:           throw FileParseException(errout.str());
    // RDKit✔️✔️:         } else if (val) {
    // RDKit✔️✔️:           if (val == 3) {
    // RDKit✔️✔️:             val = -1;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           if (atom->hasProp(common_properties::molAttachPoint)) {
    // RDKit❗✔️:             if (strictParsing) {
    // RDKit✔️✔️:               throw FileParseException(errout.str());
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:             }
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             atom->setProp(common_properties::molAttachPoint, val);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    let count = parse_required_unsigned(line, 6, 3, line_number)?;
    let mut position = 9;
    for _ in 0..count {
        let atom_id = parse_required_unsigned(line, position, 4, line_number)?;
        position += 4;
        let mut value = if line.len() >= position + 4 && rdkit_substr(line, position, 4) != "    " {
            parse_required_int(line, position, 4, line_number)?
        } else {
            0
        };
        position += 4;
        if !(0..=3).contains(&value) {
            return Err(SdfReadError::Parse(format!(
                "Value {value} from APO specification on line {line_number} is invalid"
            )));
        }
        if value == 0 {
            let _ = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
            continue;
        }
        if value == 3 {
            value = -1;
        }
        let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
        if atom.spec.prop("molAttachPoint").is_some() {
            if params.strict_parsing {
                return Err(SdfReadError::Parse(format!(
                    "Multiple ATTCHPT values for atom {atom_id} on line {line_number}"
                )));
            }
        } else {
            update_v2000_atom_spec(atom, |spec| {
                spec.with_prop("molAttachPoint", value.to_string())
            })?;
        }
    }
    Ok(())
}

fn parse_v2000_lin_line(
    line: &str,
    line_number: usize,
    atom_count: usize,
) -> Result<String, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseLinkNodeLine
    // RDKit✔️✔️:   std::string propVal = "";
    // RDKit✔️✔️:   unsigned int spos = 9;
    // RDKit✔️✔️:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️✔️:       auto aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️✔️:           text.substr(spos, 4));
    // RDKit✔️✔️:       if (!aid || aid > mol->getNumAtoms()) {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       if (text.size() < spos + 4 || text.substr(spos, 4) == "    ") {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       auto repeatCount = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️✔️:           text.substr(spos, 4));
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       if (repeatCount < 2) {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       unsigned int substB = 0;
    // RDKit✔️✔️:       unsigned int substC = 0;
    // RDKit✔️✔️:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️✔️:         substB = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️✔️:             text.substr(spos, 4));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️✔️:         substC = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️✔️:             text.substr(spos, 4));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       spos += 4;
    // RDKit✔️✔️:       if (!substB || substB > mol->getNumAtoms() ||
    // RDKit✔️✔️:           substC > mol->getNumAtoms()) {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (substC) {
    // RDKit✔️✔️:         formatter = boost::format("1 %1% 2 %2% %3% %2% %4%") % repeatCount %
    // RDKit✔️✔️:                     aid % substB % substC;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         formatter = boost::format("1 %1% 1 %2% %3%") % repeatCount % aid %
    // RDKit✔️✔️:                     substB % substC;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!propVal.empty()) {
    // RDKit✔️✔️:         propVal += "|";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       propVal += formatter.str();
    // RDKit✔️✔️:     mol->setProp(common_properties::molFileLinkNodes, propVal);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    let count = parse_required_unsigned(line, 6, 3, line_number)?;
    let mut value = String::new();
    let mut position = 9;
    for _ in 0..count {
        let atom_id = parse_required_unsigned(line, position, 4, line_number)?;
        if atom_id == 0 || atom_id as usize > atom_count {
            return Err(SdfReadError::Parse(format!(
                "LIN specification has bad atom idx on line {line_number}"
            )));
        }
        position += 4;
        if line.len() < position + 4 || rdkit_substr(line, position, 4) == "    " {
            return Err(SdfReadError::Parse(format!(
                "LIN specification missing repeat count on line {line_number}"
            )));
        }
        let repeat_count = parse_required_unsigned(line, position, 4, line_number)?;
        position += 4;
        if repeat_count < 2 {
            return Err(SdfReadError::Parse(format!(
                "LIN specification: repeat count must be >=2 on line {line_number}"
            )));
        }
        let substituent_b =
            if line.len() >= position + 4 && rdkit_substr(line, position, 4) != "    " {
                parse_required_unsigned(line, position, 4, line_number)?
            } else {
                0
            };
        position += 4;
        let substituent_c =
            if line.len() >= position + 4 && rdkit_substr(line, position, 4) != "    " {
                parse_required_unsigned(line, position, 4, line_number)?
            } else {
                0
            };
        position += 4;
        if substituent_b == 0
            || substituent_b as usize > atom_count
            || substituent_c as usize > atom_count
        {
            return Err(SdfReadError::Parse(format!(
                "LIN specification has bad substituent idx on line {line_number}"
            )));
        }
        if !value.is_empty() {
            value.push('|');
        }
        if substituent_c == 0 {
            value.push_str(&format!("1 {repeat_count} 1 {atom_id} {substituent_b}"));
        } else {
            value.push_str(&format!(
                "1 {repeat_count} 2 {atom_id} {substituent_b} {atom_id} {substituent_c}"
            ));
        }
    }
    Ok(value)
}

fn apply_v2000_property_lines(
    lines: &[&str],
    start: usize,
    atoms: &mut [ParsedV2000Atom],
    bonds: &mut [ParsedV2000Bond],
    params: MolBlockReadParams,
) -> Result<(Vec<SubstanceGroup>, bool, BTreeMap<String, String>), SdfReadError> {
    let mut first_charge_line = true;
    let mut needs_query_scan = false;
    let mut molfile_properties = BTreeMap::new();
    let mut explicit_degrees = vec![0_u8; atoms.len()];
    for bond in bonds.iter() {
        for atom in [bond.bond.begin(), bond.bond.end()] {
            explicit_degrees[atom.index()] = explicit_degrees[atom.index()].saturating_add(1);
        }
    }
    let bond_endpoints = bonds
        .iter()
        .map(|bond| (bond.bond.begin(), bond.bond.end()))
        .collect::<Vec<_>>();
    let mut sgroup_state = crate::sdf_sgroups::V2000SgroupState::new(params.strict_parsing);
    let mut cursor = start;
    while cursor < lines.len() {
        let line_number = cursor + 1;
        let line = strip_terminal_cr(lines[cursor]);
        if line.starts_with("M  END") {
            return Ok((
                sgroup_state.finish(&bond_endpoints)?,
                needs_query_scan,
                molfile_properties,
            ));
        }
        if line.is_empty() {
            // BEGIN RDKIT CPP FUNCTION ParseMolBlockProperties
            // RDKit✔️✔️:   if (!tempStr.size()) {
            // RDKit❗✔️:     if (!strictParsing) {
            // RDKit✔️✔️:       tempStr = getLine(inStream);
            // RDKit✔️✔️:       ++line;
            // RDKit❗✔️:     } else {
            // RDKit✔️✔️:       throw FileParseException(errout.str());
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:   }
            // END RDKIT CPP FUNCTION
            if params.strict_parsing {
                return Err(SdfReadError::Parse(format!(
                    "Problems encountered parsing Mol data, unexpected blank line found at line {line_number}"
                )));
            }
            cursor += 1;
            continue;
        }
        if line.starts_with('A') {
            let value_line_number = line_number + 1;
            let value = lines
                .get(cursor + 1)
                .copied()
                .map(strip_terminal_cr)
                .ok_or_else(|| {
                    SdfReadError::Parse("EOF hit while reading atom alias line".to_owned())
                })?;
            parse_v2000_atom_alias(line, value, value_line_number, atoms)?;
            cursor += 2;
            continue;
        }
        if line.starts_with('G') {
            if lines.get(cursor + 1).is_none() {
                return Err(SdfReadError::Parse(
                    "EOF hit while skipping deprecated group abbreviation".to_owned(),
                ));
            }
            cursor += 2;
            continue;
        }
        if line.starts_with('V') {
            parse_v2000_atom_value(line, line_number, atoms)?;
            cursor += 1;
            continue;
        }
        if line.starts_with("S  SKP") {
            let skip = parse_required_int(line, 6, 3, line_number)?;
            if skip < 0 {
                return Err(SdfReadError::Parse(format!(
                    "negative skip value {skip} on line {line_number}"
                )));
            }
            cursor = cursor.saturating_add(skip as usize + 1);
            continue;
        }
        if !line.starts_with(['M', 'S']) {
            // RDKit❗✔️:   // older mol files can have an atom list block here
            // RDKit❗✔️:   } else {
            // RDKit❗✔️:     if (tempStr[0] != 'M' && tempStr[0] != 'A' && tempStr[0] != 'V' &&
            // RDKit❗✔️:         tempStr[0] != 'G' && tempStr[0] != 'S') {
            // RDKit❗✔️:       ParseOldAtomList(mol, std::string_view(tempStr.c_str()), line);
            // RDKit❗✔️:     }
            // RDKit❗✔️:   }
            if cursor != start {
                return Err(SdfReadError::Unsupported(
                    "old-style V2000 atom-list records are only source-defined at the start of the property block",
                ));
            }
            parse_v2000_old_atom_list_line(line, line_number, atoms)?;
            cursor += 1;
            continue;
        }
        let prefix = rdkit_substr(line, 0, 6);
        match prefix {
            "M  CHG" => {
                // BEGIN RDKIT CPP FUNCTION ParseChargeLine
                // RDKit✔️✔️: if (firstCall) {
                // RDKit✔️✔️:   for (ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms();
                // RDKit✔️✔️:        ++ai) {
                // RDKit✔️✔️:     (*ai)->setFormalCharge(0);
                // RDKit✔️✔️:   }
                // RDKit✔️✔️: }
                if first_charge_line {
                    for atom in atoms.iter_mut() {
                        update_v2000_atom_spec(atom, |spec| Ok(spec.with_formal_charge(0)))?;
                    }
                }
                let count = parse_required_int(line, 6, 3, line_number)?;
                let mut position = 9;
                // RDKit✔️✔️: aid = FileParserUtils::toInt(text.substr(spos, 4));
                // RDKit✔️✔️: chg = FileParserUtils::toInt(text.substr(spos, 4));
                // RDKit✔️✔️: mol->getAtomWithIdx(aid - 1)->setFormalCharge(chg);
                for _ in 0..count {
                    let atom_id = parse_required_int(line, position, 4, line_number)?;
                    position += 4;
                    let charge = parse_required_int(line, position, 4, line_number)?;
                    position += 4;
                    let atom = v2000_atom_mut(atoms, atom_id, line_number)?;
                    update_v2000_atom_spec(atom, |spec| Ok(spec.with_formal_charge(charge as i8)))?;
                }
                first_charge_line = false;
                // END RDKIT CPP FUNCTION
            }
            "M  RAD" => {
                // BEGIN RDKIT CPP FUNCTION ParseRadicalLine
                // RDKit✔️✔️: if (firstCall) {
                // RDKit✔️✔️:   for (ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms();
                // RDKit✔️✔️:        ++ai) {
                // RDKit✔️✔️:     (*ai)->setFormalCharge(0);
                // RDKit✔️✔️:   }
                // RDKit✔️✔️: }
                if first_charge_line {
                    for atom in atoms.iter_mut() {
                        update_v2000_atom_spec(atom, |spec| Ok(spec.with_formal_charge(0)))?;
                    }
                }
                let count = parse_required_int(line, 6, 3, line_number)?;
                let mut position = 9;
                for _ in 0..count {
                    let atom_id = parse_required_int(line, position, 4, line_number)?;
                    position += 4;
                    let radical = parse_required_int(line, position, 4, line_number)?;
                    position += 4;
                    // RDKit✔️✔️: switch (rad) {
                    // RDKit✔️✔️:   case 0:
                    // RDKit✔️✔️:     // This shouldn't be required, but let's make sure.
                    // RDKit✔️✔️:     mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(0);
                    // RDKit✔️✔️:     break;
                    // RDKit✔️✔️:   case 1:
                    // RDKit✔️✔️:     mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(2);
                    // RDKit✔️✔️:     break;
                    // RDKit✔️✔️:   case 2:
                    // RDKit✔️✔️:     mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(1);
                    // RDKit✔️✔️:     break;
                    // RDKit✔️✔️:   case 3:
                    // RDKit✔️✔️:     mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(2);
                    // RDKit✔️✔️:     break;
                    let electrons = match radical {
                        0 => 0,
                        1 | 3 => 2,
                        2 => 1,
                        _ => {
                            return Err(SdfReadError::Parse(format!(
                                "Unrecognized radical value {radical} for atom {} on line {line_number}",
                                atom_id - 1
                            )));
                        }
                    };
                    let atom = v2000_atom_mut(atoms, atom_id, line_number)?;
                    update_v2000_atom_spec(
                        atom,
                        |spec| Ok(spec.with_radical_electrons(electrons)),
                    )?;
                }
                first_charge_line = false;
                // END RDKIT CPP FUNCTION
            }
            "M  ISO" => {
                // BEGIN RDKIT CPP FUNCTION ParseIsotopeLine
                // RDKit✔️✔️: nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
                // RDKit✔️✔️:       aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
                // RDKit✔️✔️:           text.substr(spos, 4));
                // RDKit✔️✔️: int isotope = FileParserUtils::toInt(text.substr(spos, 4));
                // RDKit✔️✔️: if (isotope < 0) {
                // RDKit✔️✔️: } else {
                // RDKit✔️✔️:   atom->setIsotope(isotope);
                // RDKit✔️✔️: }
                let count = parse_required_unsigned(line, 6, 3, line_number)?;
                let mut position = 9;
                for _ in 0..count {
                    let atom_id = parse_required_unsigned(line, position, 4, line_number)?;
                    position += 4;
                    if line.len() >= position + 4 && rdkit_substr(line, position, 4) != "    " {
                        let isotope = parse_required_int(line, position, 4, line_number)?;
                        if isotope >= 0 {
                            let atom = v2000_atom_mut(atoms, atom_id as i32, line_number)?;
                            update_v2000_atom_spec(atom, |spec| {
                                Ok(spec.with_isotope(isotope as u16))
                            })?;
                        }
                    }
                    position += 4;
                }
                // END RDKIT CPP FUNCTION
            }
            "M  STY" | "M  SST" | "M  SLB" | "M  SCN" | "M  SDS" | "M  SAL" | "M  SBL"
            | "M  SPA" | "M  SMT" | "M  SDI" | "M  SBV" | "M  SDT" | "M  SDD" | "M  SCD"
            | "M  SED" | "M  SPL" | "M  SNC" | "M  SAP" | "M  SCL" | "M  SBT" => {
                sgroup_state.parse_line(line, line_number, atoms.len(), &bond_endpoints)?;
            }
            "M  ALS" => parse_v2000_atom_list_line(line, line_number, atoms)?,
            "M  RGP" => parse_v2000_rgroup_line(line, line_number, atoms)?,
            "M  RBC" => {
                parse_v2000_ring_bond_count_line(line, line_number, atoms, &mut needs_query_scan)?
            }
            "M  SUB" => parse_v2000_substitution_line(line, line_number, atoms, &explicit_degrees)?,
            "M  UNS" => parse_v2000_unsaturation_line(line, line_number, atoms)?,
            "M  PXA" => parse_v2000_pxa_line(line, line_number, atoms)?,
            "M  ZBO" => parse_v2000_zbo_line(line, line_number, bonds)?,
            "M  ZCH" => parse_v2000_zch_line(line, line_number, atoms)?,
            "M  HYD" => parse_v2000_hyd_line(line, line_number, atoms)?,
            "M  APO" => parse_v2000_apo_line(line, line_number, atoms, params)?,
            "M  LIN" => {
                let value = parse_v2000_lin_line(line, line_number, atoms.len())?;
                molfile_properties.insert("_MolFileLinkNodes".to_owned(), value);
            }
            "M  MRV" => {
                return Err(SdfReadError::Unsupported(
                    "V2000 Marvin SMARTS requires the detached SMARTS parser contract",
                ));
            }
            _ => {}
        }
        cursor += 1;
    }
    Err(SdfReadError::Parse(
        "Problems encountered parsing Mol data, M  END missing".to_owned(),
    ))
}

fn molblock_ctab_version(counts: &str, params: MolBlockReadParams) -> Result<u16, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION MolFromMolDataStream (CTAB version dispatch)
    // RDKit✔️✔️:   unsigned int ctabVersion = 2000;
    // RDKit✔️✔️:   if (tempStr.size() > 35) {
    // RDKit✔️✔️:     if (tempStr.size() < 39 || tempStr[34] != 'V') {
    // RDKit✔️✔️:       std::ostringstream errout;
    // RDKit✔️✔️:       errout << "CTAB version string invalid at line " << line;
    // RDKit❗✔️:       if (params.strictParsing) {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:       }
    // RDKit✔️✔️:     } else if (tempStr.substr(34, 5) == "V3000") {
    // RDKit✔️✔️:       ctabVersion = 3000;
    // RDKit✔️✔️:     } else if (tempStr.substr(34, 5) != "V2000") {
    // RDKit✔️✔️:       std::ostringstream errout;
    // RDKit✔️✔️:       errout << "Unsupported CTAB version: '" << tempStr.substr(34, 5)
    // RDKit✔️✔️:              << "' at line " << line;
    // RDKit❗✔️:       if (params.strictParsing) {
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    if counts.len() <= 35 {
        return Ok(2000);
    }
    if counts.len() < 39 || counts.as_bytes().get(34) != Some(&b'V') {
        if params.strict_parsing {
            return Err(SdfReadError::Parse(
                "CTAB version string invalid at line 4".to_owned(),
            ));
        }
        return Ok(2000);
    }
    match rdkit_substr(counts, 34, 5) {
        "V3000" => Ok(3000),
        "V2000" => Ok(2000),
        version if params.strict_parsing => Err(SdfReadError::Parse(format!(
            "Unsupported CTAB version: '{version}' at line 4"
        ))),
        _ => Ok(2000),
    }
}

fn read_v2000_record_detached(
    block: &str,
    params: MolBlockReadParams,
) -> Result<MolBlockRecord, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION MolFromMolDataStream / ParseV2000CTAB
    // RDKit✔️❌: // mol name
    // RDKit✔️❌: line++;
    // RDKit✔️❌: tempStr = getLine(inStream);
    // RDKit✔️❌: res->setProp(common_properties::_Name, tempStr);
    // RDKit✔️❌: // info
    // RDKit✔️❌: line++;
    // RDKit✔️❌: tempStr = getLine(inStream);
    // RDKit✔️❌: res->setProp("_MolFileInfo", tempStr);
    // RDKit✔️❌: // comments
    // RDKit✔️❌: line++;
    // RDKit✔️❌: tempStr = getLine(inStream);
    // RDKit✔️❌: res->setProp("_MolFileComments", tempStr);
    // This convenience API buffers all input lines, unlike the source stream
    // parser, so it preserves behavior with additional O(record-size) memory.
    let lines = block.lines().collect::<Vec<_>>();
    if lines.len() < 4 {
        return Err(SdfReadError::Empty);
    }
    let title = strip_terminal_cr(lines[0]);
    let info = strip_terminal_cr(lines[1]);
    let comments = strip_terminal_cr(lines[2]);
    let counts = lines[3];
    // RDKit✔️✔️: if (tempStr.size() < 6) {
    // RDKit✔️✔️:   errout << "Counts line too short: '" << tempStr << "' on line" << line;
    // RDKit✔️✔️:   throw FileParseException(errout.str());
    // RDKit✔️✔️: }
    if counts.len() < 6 {
        return Err(SdfReadError::Counts);
    }
    if molblock_ctab_version(counts, params)? == 3000 {
        return read_v3000_record_detached(block, params);
    }
    // RDKit✔️✔️:     nAtoms = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️✔️:     spos = 3;
    // RDKit✔️✔️:     nBonds = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    let atom_count = parse_required_unsigned(counts, 0, 3, 4)? as usize;
    let bond_count = parse_required_unsigned(counts, 3, 3, 4)? as usize;
    let chiral_flag = if counts.len() >= 15 {
        parse_rdkit_unsigned(rdkit_substr(counts, 12, 3)).unwrap_or(0)
    } else {
        0
    };
    if lines.len() < 4 + atom_count + bond_count {
        return Err(SdfReadError::Counts);
    }
    let mut parsed_atoms = Vec::with_capacity(atom_count);
    let mut coords2 = Vec::with_capacity(atom_count);
    let mut coords3 = Vec::with_capacity(atom_count);
    let mut has_z = false;
    for index in 0..atom_count {
        let line_no = 5 + index;
        let parsed = parse_v2000_atom_line(lines[4 + index], line_no, params.strict_parsing)?;
        has_z |= parsed.coordinate[2] != 0.0;
        coords2.push([parsed.coordinate[0], parsed.coordinate[1]]);
        coords3.push(parsed.coordinate);
        parsed_atoms.push(parsed);
    }
    let mut bonds = Vec::with_capacity(bond_count);
    for index in 0..bond_count {
        let line_no = 5 + atom_count + index;
        bonds.push(parse_v2000_bond_line(
            lines[4 + atom_count + index],
            line_no,
            atom_count,
            index,
        )?);
    }
    let (substance_groups, needs_query_scan, molfile_properties) = apply_v2000_property_lines(
        &lines,
        4 + atom_count + bond_count,
        &mut parsed_atoms,
        &mut bonds,
        params,
    )?;
    let coordinates = CoordinateBlock {
        conformers_2d: if has_z {
            Vec::new()
        } else {
            vec![Conformer2D::new(0, coords2)]
        },
        conformers_3d: if has_z {
            vec![Conformer3D::new(0, coords3, true)]
        } else {
            Vec::new()
        },
        source_coordinate_dim: Some(if has_z {
            CoordinateDimension::ThreeD
        } else {
            CoordinateDimension::TwoD
        }),
    };
    coordinates.validate_for_atom_count(atom_count)?;
    let mut properties = if title.is_empty() {
        MoleculeProperties::default()
    } else {
        MoleculeProperties::default().with_name(title)
    };
    properties = properties
        .with_prop("_MolFileInfo", info)?
        .with_prop("_MolFileComments", comments)?
        .with_prop("_MolFileChiralFlag", chiral_flag.to_string())?;
    if needs_query_scan {
        properties = properties.with_prop("_NeedsQueryScan", "1")?;
    }
    for (key, value) in molfile_properties {
        properties = properties.with_prop(key, value)?;
    }
    let has_query = parsed_atoms.iter().any(|atom| atom.query.is_some())
        || bonds.iter().any(|bond| bond.query.is_some());
    if has_query {
        let query_atoms = parsed_atoms
            .into_iter()
            .enumerate()
            .map(|(index, parsed)| match parsed.query {
                Some(predicate) => QueryAtom::from_parts(
                    Atom::from_spec(AtomId::new(index), parsed.spec),
                    predicate,
                ),
                None => QueryAtom::new(AtomId::new(index), parsed.spec),
            })
            .collect::<Vec<_>>();
        let query_bonds = bonds
            .into_iter()
            .enumerate()
            .map(|(index, parsed)| match parsed.query {
                Some(predicate) => QueryBond::from_parts(parsed.bond, predicate),
                None => {
                    let predicate = if parsed.bond.order() == BondOrder::Unspecified {
                        QueryNode::predicate(BondQueryPredicate::Any)
                    } else {
                        QueryNode::predicate(BondQueryPredicate::Order(parsed.bond.order()))
                    };
                    debug_assert_eq!(parsed.bond.id(), BondId::new(index));
                    QueryBond::from_parts(parsed.bond, predicate)
                }
            })
            .collect::<Vec<_>>();
        let mut query_props = properties.props().clone();
        if let Some(name) = properties.name() {
            query_props.insert("_Name".to_owned(), name.to_owned());
        }
        let query = QueryGraph::from_parts(
            query_atoms,
            query_bonds,
            query_props,
            coordinates.conformers_2d.clone(),
            coordinates.conformers_3d.clone(),
            Vec::new(),
        )?;
        // END RDKIT CPP FUNCTION
        return Ok(MolBlockRecord::Query(QueryMolBlockRecord {
            query,
            substance_groups,
            properties,
            source_coordinate_dim: coordinates.source_coordinate_dim,
        }));
    }
    let atoms = parsed_atoms
        .into_iter()
        .enumerate()
        .map(|(index, parsed)| Atom::from_spec(AtomId::new(index), parsed.spec))
        .collect::<Vec<_>>();
    let concrete_bonds = bonds
        .into_iter()
        .map(|parsed| parsed.bond)
        .collect::<Vec<_>>();
    let topology = TopologyBlock {
        adjacency: AdjacencyList::from_topology(atom_count, &concrete_bonds),
        atoms,
        bonds: concrete_bonds,
        substance_groups,
        stereo_groups: Vec::new(),
    };
    topology.validate()?;
    // END RDKIT CPP FUNCTION
    Ok(MolBlockRecord::Concrete {
        topology,
        coordinates,
        properties,
    })
}

/// Read a MolBlock into either concrete topology or the canonical detached
/// query graph. Query records are never coerced into a concrete molecule.
pub fn read_mol_block_detached(block: &str) -> Result<MolBlockRecord, SdfReadError> {
    read_mol_block_detached_with_params(block, MolBlockReadParams::default())
}

/// Read a MolBlock with explicit detached syntax policy.
pub fn read_mol_block_detached_with_params(
    block: &str,
    params: MolBlockReadParams,
) -> Result<MolBlockRecord, SdfReadError> {
    read_v2000_record_detached(block, params)
}

/// Read a concrete V2000 mol block into detached model values.
pub fn read_v2000_detached(
    block: &str,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), SdfReadError> {
    read_v2000_detached_with_params(block, MolBlockReadParams::default())
}

/// Read a concrete V2000 MolBlock with explicit detached syntax policy.
pub fn read_v2000_detached_with_params(
    block: &str,
    params: MolBlockReadParams,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), SdfReadError> {
    match read_v2000_record_detached(block, params)? {
        MolBlockRecord::Concrete {
            topology,
            coordinates,
            properties,
        } => Ok((topology, coordinates, properties)),
        MolBlockRecord::Query(_) => Err(SdfReadError::Unsupported(
            "query-bearing MolBlock; use read_mol_block_detached",
        )),
    }
}

pub(super) fn get_v3000_line(
    lines: &[&str],
    cursor: &mut usize,
) -> Result<(String, usize), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION getV3000Line
    // RDKit✔️❌: ++line;
    // RDKit✔️❌: auto inl = getLine(inStream);
    // RDKit✔️❌: std::string_view tempStr = inl;
    // RDKit✔️❌: if (tempStr.size() < 7 || tempStr.substr(0, 7) != "M  V30 ") {
    // RDKit✔️❌:   std::ostringstream errout;
    // RDKit✔️❌:   errout << "Line " << line << " does not start with 'M  V30 '" << std::endl;
    // RDKit✔️❌:   throw FileParseException(errout.str());
    // RDKit✔️❌: }
    // RDKit✔️❌: while (tempStr.back() == '-') {
    // RDKit✔️❌:   res += tempStr.substr(7, tempStr.length() - 8);
    // RDKit✔️❌:   ++line;
    // RDKit✔️❌:   inl = getLine(inStream);
    // RDKit✔️❌:   tempStr = inl;
    // RDKit✔️❌:   if (tempStr.size() < 7 || tempStr.substr(0, 7) != "M  V30 ") {
    // RDKit✔️❌:     std::ostringstream errout;
    // RDKit✔️❌:     errout << "Line " << line << " does not start with 'M  V30 '"
    // RDKit✔️❌:            << std::endl;
    // RDKit✔️❌:     throw FileParseException(errout.str());
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // RDKit✔️❌: res += tempStr.substr(7, tempStr.length() - 7);
    // This slice-based convenience API has source-equivalent control flow but
    // buffers the record and allocates the returned logical line.
    let mut result = String::new();
    loop {
        let line_number = *cursor + 1;
        let line = lines.get(*cursor).copied().ok_or_else(|| {
            SdfReadError::Parse(format!(
                "Line {line_number} does not start with 'M  V30 '\n"
            ))
        })?;
        *cursor += 1;
        let line = strip_terminal_cr(line);
        if line.len() < 7 || !line.starts_with("M  V30 ") {
            return Err(SdfReadError::Parse(format!(
                "Line {line_number} does not start with 'M  V30 '\n"
            )));
        }
        if let Some(continued) = line.strip_suffix('-') {
            result.push_str(&continued[7..]);
        } else {
            result.push_str(&line[7..]);
            return Ok((result, line_number));
        }
    }
    // END RDKIT CPP FUNCTION
}

fn tokenize_v3000_line(line: &str) -> Vec<&str> {
    // BEGIN RDKIT CPP FUNCTION tokenizeV3000Line
    // RDKit✔️✔️: tokens.clear();
    // RDKit✔️✔️: bool inQuotes = false;
    // RDKit✔️✔️: unsigned int parenDepth = 0;
    // RDKit✔️✔️: unsigned int start = 0;
    // RDKit✔️✔️: unsigned int pos = 0;
    // RDKit✔️✔️: while (pos < line.size()) {
    // RDKit✔️✔️:   if (line[pos] == ' ' || line[pos] == '\t') {
    // RDKit✔️✔️:     if (start == pos) {
    // RDKit✔️✔️:       ++start;
    // RDKit✔️✔️:       ++pos;
    // RDKit✔️✔️:     } else if (!inQuotes && parenDepth == 0) {
    // RDKit✔️✔️:       tokens.push_back(line.substr(start, pos - start));
    // RDKit✔️✔️:       ++pos;
    // RDKit✔️✔️:       start = pos;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ++pos;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (line[pos] == ')' && parenDepth > 0) {
    // RDKit✔️✔️:     --parenDepth;
    // RDKit✔️✔️:     ++pos;
    // RDKit✔️✔️:   } else if (line[pos] == '(' && !inQuotes) {
    // RDKit✔️✔️:     ++parenDepth;
    // RDKit✔️✔️:     ++pos;
    // RDKit✔️✔️:   } else if (line[pos] == '"' && parenDepth == 0) {
    // RDKit✔️✔️:     if (pos + 1 < line.size() && line[pos + 1] == '"') {
    // RDKit✔️✔️:       pos += 2;
    // RDKit✔️✔️:     } else if (inQuotes) {
    // RDKit✔️✔️:       tokens.push_back(line.substr(start + 1, pos - start - 1));
    // RDKit✔️✔️:       ++pos;
    // RDKit✔️✔️:       start = pos;
    // RDKit✔️✔️:       inQuotes = false;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ++pos;
    // RDKit✔️✔️:       inQuotes = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     ++pos;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (start != pos) {
    // RDKit✔️✔️:   tokens.push_back(line.substr(start, line.size() - start));
    // RDKit✔️✔️: }
    let bytes = line.as_bytes();
    let mut tokens = Vec::new();
    let mut in_quotes = false;
    let mut paren_depth = 0_usize;
    let mut start = 0_usize;
    let mut pos = 0_usize;
    while pos < bytes.len() {
        match bytes[pos] {
            b' ' | b'\t' if start == pos => {
                start += 1;
                pos += 1;
            }
            b' ' | b'\t' if !in_quotes && paren_depth == 0 => {
                tokens.push(&line[start..pos]);
                pos += 1;
                start = pos;
            }
            b' ' | b'\t' => pos += 1,
            b')' if paren_depth > 0 => {
                paren_depth -= 1;
                pos += 1;
            }
            b'(' if !in_quotes => {
                paren_depth += 1;
                pos += 1;
            }
            b'"' if paren_depth == 0 && pos + 1 < bytes.len() && bytes[pos + 1] == b'"' => {
                pos += 2;
            }
            b'"' if paren_depth == 0 && in_quotes => {
                tokens.push(&line[start + 1..pos]);
                pos += 1;
                start = pos;
                in_quotes = false;
            }
            b'"' if paren_depth == 0 => {
                pos += 1;
                in_quotes = true;
            }
            _ => pos += 1,
        }
    }
    if start != pos {
        tokens.push(&line[start..pos]);
    }
    // END RDKIT CPP FUNCTION
    tokens
}

fn split_v3000_assignment(token: &str) -> Option<(String, &str)> {
    let position = token.find('=')?;
    if position != token.rfind('=')? {
        return None;
    }
    Some((
        token[..position].to_ascii_uppercase(),
        &token[position + 1..],
    ))
}

fn parse_v3000_i32(value: &str, kind: &'static str, line: usize) -> Result<i32, SdfReadError> {
    parse_rdkit_int(value).map_err(|()| SdfReadError::Field {
        kind,
        line,
        value: value.to_owned(),
    })
}

fn v3000_atom_symbol(
    symbol: &str,
    line: usize,
) -> Result<
    (
        Element,
        Option<u16>,
        Option<String>,
        Option<QueryNode<AtomQueryPredicate>>,
    ),
    SdfReadError,
> {
    let mut symbol = symbol.trim();
    let mut negate = false;
    if symbol.len() > 3 && symbol[..3].eq_ignore_ascii_case("NOT") {
        negate = true;
        symbol = symbol[3..].trim();
    }
    if symbol.starts_with('[') {
        if !symbol.ends_with(']') {
            return Err(SdfReadError::Parse(format!(
                "Bad atom token '{symbol}' on line: {line}"
            )));
        }
        let mut numbers = Vec::new();
        for item in symbol[1..symbol.len() - 1].split(',') {
            let item = item.trim();
            if item.is_empty() {
                continue;
            }
            let normalized;
            let lookup = if item.len() == 2 && item.as_bytes()[1].is_ascii_uppercase() {
                normalized = format!(
                    "{}{}",
                    &item[..1],
                    (item.as_bytes()[1] as char).to_ascii_lowercase()
                );
                normalized.as_str()
            } else {
                item
            };
            let element = Element::from_symbol(lookup).ok_or_else(|| {
                SdfReadError::Parse(format!("Element '{item}' not found on line {line}"))
            })?;
            numbers.push(element.atomic_number());
        }
        let predicate = if negate {
            AtomQueryPredicate::AtomicNumberNotIn(numbers)
        } else {
            AtomQueryPredicate::AtomicNumberIn(numbers)
        };
        return Ok((
            Element::DUMMY,
            None,
            None,
            Some(QueryNode::predicate(predicate)),
        ));
    }
    if negate {
        return Err(SdfReadError::Parse(format!(
            "NOT tokens only supported for atom lists. line {line}"
        )));
    }
    if let Some(query) = complex_molfile_atom_query(symbol) {
        return Ok((Element::DUMMY, None, None, Some(query)));
    }
    if symbol == "*" {
        return Ok((
            Element::DUMMY,
            None,
            None,
            Some(QueryNode::predicate(AtomQueryPredicate::Any)),
        ));
    }
    if matches!(symbol, "R" | "R#") {
        return Ok((Element::DUMMY, None, Some(symbol.to_owned()), None));
    }
    if let Some(label) = symbol.strip_prefix('R')
        && !label.is_empty()
        && label.len() <= 2
        && label.bytes().all(|byte| byte.is_ascii_digit())
    {
        return Ok((
            Element::DUMMY,
            Some(label.parse::<u16>().unwrap_or(0)),
            Some(symbol.to_owned()),
            None,
        ));
    }
    if symbol == "D" {
        return Ok((Element::H, Some(2), None, None));
    }
    if symbol == "T" {
        return Ok((Element::H, Some(3), None, None));
    }
    if matches!(symbol, "Pol" | "Mod") {
        return Ok((Element::DUMMY, None, Some(symbol.to_owned()), None));
    }
    let normalized;
    let lookup = if symbol.len() == 2 && symbol.as_bytes()[1].is_ascii_uppercase() {
        normalized = format!(
            "{}{}",
            &symbol[..1],
            (symbol.as_bytes()[1] as char).to_ascii_lowercase()
        );
        normalized.as_str()
    } else {
        symbol
    };
    Element::from_symbol(lookup)
        .map(|element| (element, None, None, None))
        .ok_or_else(|| SdfReadError::Parse(format!("Element '{symbol}' not found")))
}

fn parse_v3000_atom_properties(
    mut spec: AtomSpec,
    mut query: Option<QueryNode<AtomQueryPredicate>>,
    tokens: &[&str],
    atom_index: usize,
    line: usize,
) -> Result<(AtomSpec, Option<QueryNode<AtomQueryPredicate>>), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000AtomProps
    // RDKit❗✔️:   while (token != tokens.end()) {
    // RDKit❗✔️:     std::string prop;
    // RDKit❗✔️:     std::string_view val;
    // RDKit❗✔️:     if (!splitAssignToken(*token, prop, val)) {
    // RDKit❗✔️:       errout << "Invalid atom property: '" << *token << "' for atom "
    // RDKit❗✔️:              << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit❗✔️:       throw FileParseException(errout.str());
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (prop == "CHG") {
    // RDKit❗✔️:       auto charge = FileParserUtils::toInt(val);
    // RDKit❗✔️:       if (!atom->hasQuery()) {
    // RDKit❗✔️:         atom->setFormalCharge(charge);
    // RDKit❗✔️:       } else {
    // RDKit❌❌:         atom->expandQuery(makeAtomFormalChargeQuery(charge));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (prop == "RAD") {
    // RDKit❗✔️:       // FIX handle queries here
    // RDKit❗✔️:       switch (FileParserUtils::toInt(val)) {
    // RDKit❗✔️:         case 0:
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case 1:
    // RDKit❗✔️:           atom->setNumRadicalElectrons(2);
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case 2:
    // RDKit❗✔️:           atom->setNumRadicalElectrons(1);
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case 3:
    // RDKit❗✔️:           atom->setNumRadicalElectrons(2);
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         default:
    // RDKit❗✔️:           errout << "Unrecognized RAD value " << val << " for atom "
    // RDKit❗✔️:                  << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit❗✔️:           throw FileParseException(errout.str());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (prop == "MASS") {
    // RDKit❗✔️:       int v;
    // RDKit❗✔️:       double dv;
    // RDKit❗✔️:       try {
    // RDKit❗✔️:         v = FileParserUtils::toInt(val);
    // RDKit❗✔️:       } catch (boost::bad_lexical_cast &) {
    // RDKit❗✔️:         try {
    // RDKit❗✔️:           dv = FileParserUtils::toDouble(val);
    // RDKit❗✔️:           v = static_cast<int>(floor(dv));
    // RDKit❗✔️:         } catch (boost::bad_lexical_cast &) {
    // RDKit❗✔️:           v = -1;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (v < 0) {
    // RDKit❗✔️:         errout << "Bad value for MASS :" << val << " for atom "
    // RDKit❗✔️:                << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit❗✔️:         throw FileParseException(errout.str());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         if (!atom->hasQuery()) {
    // RDKit❗✔️:           atom->setIsotope(v);
    // RDKit❗✔️:         } else {
    // RDKit❌❌:           atom->expandQuery(makeAtomIsotopeQuery(v));
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (prop == "CFG") {
    // RDKit❗✔️:       auto cfg = FileParserUtils::toInt(val);
    // RDKit❗✔️:       switch (cfg) {
    // RDKit❗✔️:         case 0:
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case 1:
    // RDKit❗✔️:         case 2:
    // RDKit❗✔️:         case 3:
    // RDKit❗✔️:           atom->setProp(common_properties::molParity, cfg);
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         default:
    // RDKit❗✔️:           errout << "Unrecognized CFG value : " << val << " for atom "
    // RDKit❗✔️:                  << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit❗✔️:           throw FileParseException(errout.str());
    // RDKit❗✔️:       }
    // RDKit❌❌:     } else if (prop == "HCOUNT") {
    // RDKit❌❌:       if (val != "0") {
    // RDKit❌❌:         atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
    // RDKit❌❌:         atom->expandQuery(makeAtomImplicitHCountQuery(0));
    // RDKit❌❌:       }
    // RDKit❌❌:     } else if (prop == "UNSAT") {
    // RDKit❌❌:       if (val == "1") {
    // RDKit❌❌:         atom->expandQuery(makeAtomUnsaturatedQuery());
    // RDKit❌❌:       }
    // RDKit❌❌:     } else if (prop == "RBCNT") {
    // RDKit❌❌:       if (val != "0") {
    // RDKit❌❌:         atom->expandQuery(makeAtomRingBondCountQuery(rbcount));
    // RDKit❌❌:       }
    // RDKit❗✔️:     } else if (prop == "VAL") {
    // RDKit❗✔️:       if (val != "0") {
    // RDKit❗✔️:         auto totval = FileParserUtils::toInt(val);
    // RDKit❗✔️:         atom->setProp(common_properties::molTotValence, totval);
    // RDKit❗✔️:       }
    // RDKit❌❌:     } else if (prop == "RGROUPS") {
    // RDKit❌❌:       ParseV3000RGroups(mol, atom, val, line);
    // RDKit❗✔️:     } else if (prop == "STBOX") {
    // RDKit❗✔️:       if (val != "0") {
    // RDKit❗✔️:         auto ival = FileParserUtils::toInt(val);
    // RDKit❗✔️:         atom->setProp(common_properties::molStereoCare, ival);
    // RDKit❗✔️:       }
    // RDKit❌❌:     } else if (prop == "SUBST") {
    // RDKit❌❌:       if (val != "0") {
    // RDKit❌❌:         auto ival = FileParserUtils::toInt(val);
    // RDKit❌❌:         atom->setProp(common_properties::molSubstCount, ival);
    // RDKit❌❌:       }
    // RDKit❗✔️:     } else if (prop == "EXACHG") {
    // RDKit❗✔️:       if (val != "0") {
    // RDKit❗✔️:         auto ival = FileParserUtils::toInt(val);
    // RDKit❗✔️:         atom->setProp(common_properties::molRxnExactChange, ival);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (prop == "INVRET") {
    // RDKit❗✔️:       if (val != "0") {
    // RDKit❗✔️:         auto ival = FileParserUtils::toInt(val);
    // RDKit❗✔️:         atom->setProp(common_properties::molInversionFlag, ival);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (prop == "ATTCHPT") {
    // RDKit❗✔️:       if (val != "0") {
    // RDKit❗✔️:         auto ival = FileParserUtils::toInt(val);
    // RDKit❗✔️:         atom->setProp(common_properties::molAttachPoint, ival);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (prop == "ATTCHORD") {
    // RDKit❗✔️:       auto ival = FileParserUtils::toInt(val);
    // RDKit❗✔️:       atom->setProp(common_properties::molAttachOrder, ival);
    // RDKit❗✔️:     } else if (prop == "CLASS") {
    // RDKit❗✔️:       atom->setProp(common_properties::molAtomClass, std::string(val));
    // RDKit❗✔️:     } else if (prop == "SEQID") {
    // RDKit❗✔️:       if (val != "0") {
    // RDKit❗✔️:         auto ival = FileParserUtils::toInt(val);
    // RDKit❗✔️:         atom->setProp(common_properties::molAtomSeqId, ival);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (prop == "SEQNAME") {
    // RDKit❗✔️:       if (val != "") {
    // RDKit❗✔️:         atom->setProp(common_properties::molAtomSeqName, std::string(val));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ++token;
    // RDKit❗✔️:   }
    // Concrete properties are reproduced. Query-generating branches fail
    // explicitly because this detached API returns `TopologyBlock`, not query
    // topology; this is the documented behavioral boundary for this function.
    let mut has_attach_point = false;
    for token in tokens {
        let (property, value) = split_v3000_assignment(token).ok_or_else(|| {
            SdfReadError::Parse(format!(
                "Invalid atom property: '{token}' for atom {} on line {line}",
                atom_index + 1
            ))
        })?;
        match property.as_str() {
            "CHG" => {
                let charge = parse_v3000_i32(value, "V3000 atom charge", line)?;
                let charge = i8::try_from(charge).map_err(|_| {
                    SdfReadError::Unsupported(
                        "V3000 atom charges outside the detached i8 charge model",
                    )
                })?;
                if query.is_some() {
                    let predicate = QueryNode::predicate(AtomQueryPredicate::FormalCharge(charge));
                    query = Some(match query {
                        Some(existing) => QueryNode::and(vec![existing, predicate]),
                        None => predicate,
                    });
                } else {
                    spec = spec.with_formal_charge(charge);
                }
            }
            "RAD" => {
                let radical = parse_v3000_i32(value, "V3000 radical", line)?;
                let electrons = match radical {
                    0 => 0,
                    1 | 3 => 2,
                    2 => 1,
                    _ => {
                        return Err(SdfReadError::Parse(format!(
                            "Unrecognized RAD value {value} for atom {} on line {line}",
                            atom_index + 1
                        )));
                    }
                };
                spec = spec.with_radical_electrons(electrons);
            }
            "MASS" => {
                let isotope = value
                    .parse::<i32>()
                    .or_else(|_| parse_rdkit_double(value).map(|mass| mass.floor() as i32))
                    .map_err(|_| {
                        SdfReadError::Parse(format!(
                            "Bad value for MASS :{value} for atom {} on line {line}",
                            atom_index + 1
                        ))
                    })?;
                let isotope = u16::try_from(isotope).map_err(|_| {
                    SdfReadError::Unsupported(
                        "V3000 isotope mass outside the detached u16 isotope model",
                    )
                })?;
                if query.is_some() {
                    let predicate = QueryNode::predicate(AtomQueryPredicate::Isotope(isotope));
                    query = Some(match query {
                        Some(existing) => QueryNode::and(vec![existing, predicate]),
                        None => predicate,
                    });
                } else {
                    spec = spec.with_isotope(isotope);
                }
            }
            "CFG" => {
                let configuration = parse_v3000_i32(value, "V3000 atom CFG", line)?;
                if !matches!(configuration, 0..=3) {
                    return Err(SdfReadError::Parse(format!(
                        "Unrecognized CFG value : {value} for atom {} on line {line}",
                        atom_index + 1
                    )));
                }
                if configuration != 0 {
                    spec = spec
                        .with_mol_parity(configuration)
                        .with_prop("molParity", configuration.to_string())?;
                }
            }
            "HCOUNT" if value != "0" => {
                let mut count = parse_v3000_i32(value, "V3000 HCOUNT", line)?;
                if count == -1 {
                    count = 0;
                }
                let predicate = if count > 0 {
                    AtomQueryPredicate::ImplicitHydrogenCountLessEqual(count as u8)
                } else {
                    AtomQueryPredicate::ImplicitHydrogenCount(0)
                };
                let predicate = QueryNode::predicate(predicate);
                query = Some(match query {
                    Some(existing) => QueryNode::and(vec![existing, predicate]),
                    None => predicate,
                });
            }
            "UNSAT" if value == "1" => {
                let predicate = QueryNode::predicate(AtomQueryPredicate::IsUnsaturated);
                query = Some(match query {
                    Some(existing) => QueryNode::and(vec![existing, predicate]),
                    None => predicate,
                });
            }
            "RBCNT" if value != "0" => {
                let mut count = parse_v3000_i32(value, "V3000 RBCNT", line)?;
                spec = spec.with_prop("molRingBondCount", count.to_string())?;
                if count == -1 {
                    count = 0;
                }
                let predicate = if count == -2 {
                    AtomQueryPredicate::RingBondCount(0xDEAD_BEEF)
                } else if count > 4 {
                    AtomQueryPredicate::RingBondCountLessEqual(4)
                } else {
                    AtomQueryPredicate::RingBondCount(count as u32)
                };
                let predicate = QueryNode::predicate(predicate);
                query = Some(match query {
                    Some(existing) => QueryNode::and(vec![existing, predicate]),
                    None => predicate,
                });
            }
            "RGROUPS" => {
                let labels = parse_v3000_rgroups(value, line)?;
                for label in labels {
                    let isotope = u16::try_from(label).map_err(|_| {
                        SdfReadError::Unsupported(
                            "V3000 R-group labels outside the detached u16 isotope model",
                        )
                    })?;
                    spec = spec
                        .with_prop("_MolFileRLabel", label.to_string())?
                        .with_prop("dummyLabel", format!("R{label}"))?
                        .with_isotope(isotope);
                    query = Some(QueryNode::predicate(AtomQueryPredicate::Any));
                }
            }
            "VAL" if value != "0" => {
                let parsed = parse_v3000_i32(value, "V3000 total valence", line)?;
                spec = spec.with_prop("molTotValence", parsed.to_string())?;
            }
            "STBOX" if value != "0" => {
                let parsed = parse_v3000_i32(value, "V3000 stereo care", line)?;
                spec = spec.with_prop("molStereoCare", parsed.to_string())?;
            }
            "SUBST" if value != "0" => {
                let parsed = parse_v3000_i32(value, "V3000 substitution count", line)?;
                spec = spec.with_prop("molSubstCount", parsed.to_string())?;
            }
            "EXACHG" if value != "0" => {
                let parsed = parse_v3000_i32(value, "V3000 exact-change flag", line)?;
                spec = spec.with_prop("molRxnExactChange", parsed.to_string())?;
            }
            "INVRET" if value != "0" => {
                let parsed = parse_v3000_i32(value, "V3000 inversion flag", line)?;
                spec = spec.with_mol_inversion_flag(parsed);
            }
            "ATTCHPT" if value != "0" => {
                if has_attach_point {
                    return Err(SdfReadError::Parse(format!(
                        "Multiple ATTCHPT values for atom {} on line {line}",
                        atom_index + 1
                    )));
                }
                let parsed = parse_v3000_i32(value, "V3000 attachment point", line)?;
                spec = spec.with_prop("molAttachPoint", parsed.to_string())?;
                has_attach_point = true;
            }
            "ATTCHORD" if value.starts_with('(') => {
                return Err(SdfReadError::Unsupported(
                    "V3000 template attachment order is not lowered into detached topology",
                ));
            }
            "ATTCHORD" => {
                let parsed = parse_v3000_i32(value, "V3000 attachment order", line)?;
                spec = spec.with_prop("molAttachOrder", parsed.to_string())?;
            }
            "CLASS" => spec = spec.with_prop("molAtomClass", value)?,
            "SEQID" if value != "0" => {
                let parsed = parse_v3000_i32(value, "V3000 sequence id", line)?;
                spec = spec.with_prop("molAtomSeqId", parsed.to_string())?;
            }
            "SEQNAME" if !value.is_empty() => spec = spec.with_prop("molAtomSeqName", value)?,
            _ => {}
        }
    }
    // END RDKIT CPP FUNCTION
    Ok((spec, query))
}

fn parse_v3000_rgroups(text: &str, line: usize) -> Result<Vec<u32>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000RGroups
    // RDKit✔️✔️:   if (text[0] != '(' || text.back() != ')') {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Bad RGROUPS specification '" << text << "' on line " << line
    // RDKit✔️✔️:            << ". Missing parens.";
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<std::string> splitToken;
    // RDKit✔️✔️:   std::string resid = std::string(text.substr(1, text.size() - 2));
    // RDKit✔️✔️:   boost::split(splitToken, resid, boost::is_any_of(std::string(" ")));
    // RDKit✔️✔️:   if (splitToken.size() < 1) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Bad RGROUPS specification '" << text << "' on line " << line
    // RDKit✔️✔️:            << ". Missing values.";
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int nRs;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     nRs = FileParserUtils::stripSpacesAndCast<unsigned int>(splitToken[0]);
    // RDKit✔️✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Cannot convert '" << splitToken[0] << "' to int on line" << line;
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (splitToken.size() < nRs + 1) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Bad RGROUPS specification '" << text << "' on line " << line
    // RDKit✔️✔️:            << ". Not enough values.";
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nRs; ++i) {
    // RDKit✔️✔️:     unsigned int rLabel;
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       rLabel =
    // RDKit✔️✔️:           FileParserUtils::stripSpacesAndCast<unsigned int>(splitToken[i + 1]);
    // RDKit✔️✔️:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:       std::ostringstream errout;
    // RDKit✔️✔️:       errout << "Cannot convert '" << splitToken[i + 1] << "' to int on line"
    // RDKit✔️✔️:              << line;
    // RDKit✔️✔️:       throw FileParseException(errout.str());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
    // RDKit✔️✔️:     atom->setProp(common_properties::_MolFileRLabel, rLabel);
    // RDKit✔️✔️:     std::string dLabel = "R" + std::to_string(rLabel);
    // RDKit✔️✔️:     atom->setProp(common_properties::dummyLabel, dLabel);
    // RDKit✔️✔️:     atom->setIsotope(rLabel);
    // RDKit✔️✔️:     atom->setQuery(makeAtomNullQuery());
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    if !text.starts_with('(') || !text.ends_with(')') {
        return Err(SdfReadError::Parse(format!(
            "Bad RGROUPS specification '{text}' on line {line}. Missing parens."
        )));
    }
    let tokens = text[1..text.len() - 1].split(' ').collect::<Vec<_>>();
    let count = tokens
        .first()
        .and_then(|token| parse_rdkit_unsigned(token).ok())
        .ok_or_else(|| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line{line}",
                tokens.first().copied().unwrap_or_default()
            ))
        })? as usize;
    if tokens.len() < count + 1 {
        return Err(SdfReadError::Parse(format!(
            "Bad RGROUPS specification '{text}' on line {line}. Not enough values."
        )));
    }
    tokens
        .iter()
        .skip(1)
        .take(count)
        .map(|token| {
            parse_rdkit_unsigned(token).map_err(|()| {
                SdfReadError::Parse(format!("Cannot convert '{token}' to int on line{line}"))
            })
        })
        .collect()
}

fn parse_v3000_bond_properties(
    mut spec: BondSpec,
    mut query: Option<QueryNode<BondQueryPredicate>>,
    bond_type: u32,
    tokens: &[&str],
    line: usize,
) -> Result<(BondSpec, Option<QueryNode<BondQueryPredicate>>, bool), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000BondBlock (property loop)
    // RDKit❗✔️:     while (lPos < splitLine.size()) {
    // RDKit❗✔️:       std::string prop;
    // RDKit❗✔️:       std::string_view val;
    // RDKit❗✔️:       if (!splitAssignToken(splitLine[lPos], prop, val)) {
    // RDKit❗✔️:         errout << "bad bond property '" << splitLine[lPos] << "' on line "
    // RDKit❗✔️:                << line;
    // RDKit❗✔️:         throw FileParseException(errout.str());
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (prop == "CFG") {
    // RDKit❗✔️:         unsigned int cfg = 0;
    // RDKit❗✔️:         std::from_chars(val.data(), val.data() + val.size(), cfg);
    // RDKit❗✔️:         switch (cfg) {
    // RDKit❗✔️:           case 0:
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           case 1:
    // RDKit❗✔️:             bond->setBondDir(Bond::BEGINWEDGE);
    // RDKit❗✔️:             chiralityPossible = true;
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           case 2:
    // RDKit❗✔️:             if (bType == 1) {
    // RDKit❗✔️:               bond->setBondDir(Bond::UNKNOWN);
    // RDKit❗✔️:             } else if (bType == 2) {
    // RDKit❗✔️:               bond->setBondDir(Bond::EITHERDOUBLE);
    // RDKit❗✔️:               bond->setStereo(Bond::STEREOANY);
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           case 3:
    // RDKit❗✔️:             bond->setBondDir(Bond::BEGINDASH);
    // RDKit❗✔️:             chiralityPossible = true;
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           default:
    // RDKit❗✔️:             errout << "bad bond CFG " << val << "' on line " << line;
    // RDKit❗✔️:             throw FileParseException(errout.str());
    // RDKit❗✔️:         }
    // RDKit❗✔️:         bond->setProp(common_properties::_MolFileBondCfg, cfg);
    // RDKit❌❌:       } else if (prop == "TOPO") {
    // RDKit❌❌:         if (val != "0") {
    // RDKit❌❌:           BOND_EQUALS_QUERY *q = makeBondIsInRingQuery();
    // RDKit❌❌:           bond->expandQuery(q);
    // RDKit❌❌:         }
    // RDKit❗✔️:       } else if (prop == "RXCTR") {
    // RDKit❗✔️:         int reactStatus = FileParserUtils::toInt(val);
    // RDKit❗✔️:         bond->setProp(common_properties::molReactStatus, reactStatus);
    // RDKit❗✔️:       } else if (prop == "STBOX") {
    // RDKit❗✔️:         bond->setProp(common_properties::molStereoCare, std::string(val));
    // RDKit❗✔️:       } else if (prop == "ENDPTS") {
    // RDKit❗✔️:         bond->setProp(common_properties::_MolFileBondEndPts, std::string(val));
    // RDKit❗✔️:       } else if (prop == "ATTACH") {
    // RDKit❗✔️:         bond->setProp(common_properties::_MolFileBondAttach, std::string(val));
    // RDKit❗✔️:       }
    // RDKit❗✔️:       ++lPos;
    // RDKit❗✔️:     }
    let mut chirality_possible = false;
    for token in tokens {
        let (property, value) = split_v3000_assignment(token).ok_or_else(|| {
            SdfReadError::Parse(format!("bad bond property '{token}' on line {line}"))
        })?;
        match property.as_str() {
            "CFG" => {
                let configuration = parse_rdkit_unsigned(value).unwrap_or(0);
                spec = match configuration {
                    0 => spec,
                    1 => {
                        chirality_possible = true;
                        spec.with_direction(BondDirection::BeginWedge)
                    }
                    2 if bond_type == 1 => spec.with_direction(BondDirection::Unknown),
                    2 if bond_type == 2 => spec
                        .with_direction(BondDirection::EitherDouble)
                        .with_stereo(BondStereo::Any),
                    2 => spec,
                    3 => {
                        chirality_possible = true;
                        spec.with_direction(BondDirection::BeginDash)
                    }
                    _ => {
                        return Err(SdfReadError::Parse(format!(
                            "bad bond CFG {value}' on line {line}"
                        )));
                    }
                };
                spec = spec.with_prop("_MolFileBondCfg", configuration.to_string())?;
            }
            "TOPO" if value != "0" => {
                let predicate = match parse_v3000_i32(value, "V3000 bond topology", line)? {
                    1 => BondQueryPredicate::IsInRing(true),
                    2 => BondQueryPredicate::IsInRing(false),
                    topology => {
                        return Err(SdfReadError::Parse(format!(
                            "Unrecognized bond topology specifier: {topology} on line {line}"
                        )));
                    }
                };
                let predicate = QueryNode::predicate(predicate);
                query = Some(match query {
                    Some(existing) => QueryNode::and(vec![existing, predicate]),
                    None => predicate,
                });
            }
            "RXCTR" => {
                let parsed = parse_v3000_i32(value, "V3000 bond reaction status", line)?;
                spec = spec.with_prop("molReactStatus", parsed.to_string())?;
            }
            "STBOX" => spec = spec.with_prop("molStereoCare", value)?,
            "ENDPTS" => spec = spec.with_prop("_MolFileBondEndPts", value)?,
            "ATTACH" => spec = spec.with_prop("_MolFileBondAttach", value)?,
            _ => {}
        }
    }
    // END RDKIT CPP FUNCTION
    Ok((spec, query, chirality_possible))
}

fn molfile_info_marks_3d(info: &str) -> bool {
    info.len() >= 22 && matches!(rdkit_substr(info, 20, 2), "3d" | "3D")
}

fn v3000_is_3d(info: &str, coordinates: &[[f64; 3]], chirality_possible: bool) -> bool {
    // BEGIN RDKIT CPP FUNCTION calculate3dFlag
    // RDKit✔️✔️: bool nonzeroZ = hasNonZeroZCoords(conf);
    // RDKit✔️✔️: if (!nonzeroZ && marked3d == 1) {
    // RDKit✔️✔️:   if (chiralityPossible) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: } else if (marked3d == 0 && nonzeroZ) {
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return nonzeroZ;
    let nonzero_z = coordinates.iter().any(|coordinate| coordinate[2] != 0.0);
    if !nonzero_z && molfile_info_marks_3d(info) {
        return !chirality_possible;
    }
    nonzero_z
    // END RDKIT CPP FUNCTION
}

fn read_v3000_record_detached(
    block: &str,
    params: MolBlockReadParams,
) -> Result<MolBlockRecord, SdfReadError> {
    let lines = block.lines().collect::<Vec<_>>();
    if lines.len() < 4 {
        return Err(SdfReadError::Empty);
    }
    let title = strip_terminal_cr(lines[0]);
    let info = strip_terminal_cr(lines[1]);
    let comments = strip_terminal_cr(lines[2]);
    let outer_counts = strip_terminal_cr(lines[3]);
    if outer_counts.len() < 6 {
        return Err(SdfReadError::Counts);
    }
    let outer_atom_count = parse_required_unsigned(outer_counts, 0, 3, 4)?;
    let outer_bond_count = parse_required_unsigned(outer_counts, 3, 3, 4)?;
    // BEGIN RDKIT CPP FUNCTION MolFromMolDataStream (V3000 outer counts)
    // RDKit✔️✔️:       if (nAtoms != 0 || nBonds != 0) {
    // RDKit✔️✔️:         std::ostringstream errout;
    // RDKit✔️✔️:         errout << "V3000 mol blocks should have 0s in the initial counts line. "
    // RDKit✔️✔️:                   "(line: "
    // RDKit✔️✔️:                << line << ")";
    // RDKit❗✔️:         if (params.strictParsing) {
    // RDKit✔️✔️:           throw FileParseException(errout.str());
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:         }
    // RDKit✔️✔️:       }
    // END RDKIT CPP FUNCTION
    if params.strict_parsing && (outer_atom_count != 0 || outer_bond_count != 0) {
        return Err(SdfReadError::Parse(
            "V3000 mol blocks should have 0s in the initial counts line. (line: 4)".to_owned(),
        ));
    }
    let mut cursor = 4_usize;

    // BEGIN RDKIT CPP FUNCTION ParseV3000CTAB
    // RDKit❗✔️:   tempStr = getV3000Line(inStream, line);
    // RDKit❗✔️:   boost::to_upper(tempStr);
    // RDKit❗✔️:   if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN CTAB") {
    // RDKit❗✔️:     std::ostringstream errout;
    // RDKit❗✔️:     errout << "BEGIN CTAB line not found on line " << line;
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️:   }
    let (begin_ctab, _) = get_v3000_line(&lines, &mut cursor)?;
    if !begin_ctab.to_ascii_uppercase().starts_with("BEGIN CTAB") {
        return Err(SdfReadError::Parse("BEGIN CTAB line not found".to_owned()));
    }
    let (counts, counts_line) = get_v3000_line(&lines, &mut cursor)?;
    let upper_counts = counts.to_ascii_uppercase();
    let count_text = upper_counts
        .strip_prefix("COUNTS ")
        .ok_or(SdfReadError::Counts)?;
    let count_fields = count_text.split_whitespace().collect::<Vec<_>>();
    if count_fields.len() < 2 {
        return Err(SdfReadError::Counts);
    }
    // RDKit❗✔️: nAtoms = FileParserUtils::toUnsigned(splitLine[0]);
    // RDKit❗✔️: nBonds = FileParserUtils::toUnsigned(splitLine[1]);
    // RDKit❗✔️: if (splitLine.size() > 2) {
    // RDKit❗✔️:   nSgroups = FileParserUtils::toUnsigned(splitLine[2]);
    // RDKit❗✔️: }
    // RDKit❗✔️: if (splitLine.size() > 3) {
    // RDKit❗✔️:   n3DConstraints = FileParserUtils::toUnsigned(splitLine[3]);
    // RDKit❗✔️: }
    // RDKit❗✔️: if (splitLine.size() > 4) {
    // RDKit❗✔️:   chiralFlag = FileParserUtils::toUnsigned(splitLine[4]);
    // RDKit❗✔️: }
    let atom_count = field(count_fields[0], "atom count", counts_line)?;
    let bond_count = field(count_fields[1], "bond count", counts_line)?;
    let sgroup_count = count_fields
        .get(2)
        .map_or(Ok(0), |value| field(value, "SGroup count", counts_line))?;
    let constraint_count = count_fields.get(3).map_or(Ok(0), |value| {
        field(value, "3D constraint count", counts_line)
    })?;
    let chiral_flag: u32 = count_fields
        .get(4)
        .map_or(Ok(0), |value| field(value, "chiral flag", counts_line))?;
    if constraint_count != 0 {
        return Err(SdfReadError::Unsupported(
            "V3000 OBJ3D constraints are not yet lowered into detached coordinates",
        ));
    }

    let mut atoms = Vec::with_capacity(atom_count);
    let mut coordinates_3d = Vec::with_capacity(atom_count);
    let mut atom_by_bookmark = BTreeMap::new();
    if atom_count != 0 {
        let (begin_atom, _) = get_v3000_line(&lines, &mut cursor)?;
        if !begin_atom.starts_with("BEGIN ATOM") {
            return Err(SdfReadError::Parse("BEGIN ATOM line not found".to_owned()));
        }
        for atom_index in 0..atom_count {
            // RDKit❗✔️: tokenizeV3000Line(trimmed, tokens);
            // RDKit❗✔️: unsigned int molIdx = 0;
            // RDKit❗✔️: std::from_chars(token->data(), token->data() + token->size(), molIdx);
            // RDKit❗✔️: pos.x = atof(std::string(*token).c_str());
            // RDKit❗✔️: pos.y = atof(std::string(*token).c_str());
            // RDKit❗✔️: pos.z = atof(std::string(*token).c_str());
            // RDKit❗✔️: int mapNum = atoi(std::string(*token).c_str());
            // RDKit❗✔️: if (mapNum > 0) {
            // RDKit❗✔️:   atom->setProp(common_properties::molAtomMapNumber, mapNum);
            // RDKit❗✔️: }
            let (atom_line, line_number) = get_v3000_line(&lines, &mut cursor)?;
            let tokens = tokenize_v3000_line(atom_line.trim());
            if tokens.len() < 6 {
                return Err(SdfReadError::Parse(format!(
                    "Bad atom line : '{atom_line}' on line {line_number}"
                )));
            }
            let bookmark = tokens[0].parse::<u32>().unwrap_or(0);
            let (element, isotope, dummy_label, query) = v3000_atom_symbol(tokens[1], line_number)?;
            let mut spec = AtomSpec::new(element);
            if query.is_some() {
                spec = spec.with_no_implicit(true);
            }
            if let Some(isotope) = isotope {
                spec = spec.with_isotope(isotope);
            }
            if let Some(label) = dummy_label {
                spec = spec.with_prop("dummyLabel", label)?;
            }
            let point = [
                tokens[2].parse::<f64>().unwrap_or(0.0),
                tokens[3].parse::<f64>().unwrap_or(0.0),
                tokens[4].parse::<f64>().unwrap_or(0.0),
            ];
            let atom_map = tokens[5].parse::<i64>().unwrap_or(0);
            if atom_map > 0 {
                let atom_map = u32::try_from(atom_map).map_err(|_| {
                    SdfReadError::Unsupported(
                        "V3000 atom map outside the detached u32 atom-map model",
                    )
                })?;
                spec = spec
                    .with_atom_map(atom_map)
                    .with_prop("molAtomMapNumber", atom_map.to_string())?;
            }
            let (spec, query) =
                parse_v3000_atom_properties(spec, query, &tokens[6..], atom_index, line_number)?;
            let atom_id = AtomId::new(atom_index);
            if atom_by_bookmark.insert(bookmark, atom_id).is_some() {
                return Err(SdfReadError::Parse(format!(
                    "duplicate V3000 atom index {bookmark} on line {line_number}"
                )));
            }
            atoms.push((Atom::from_spec(atom_id, spec), query));
            coordinates_3d.push(point);
        }
        let (end_atom, line_number) = get_v3000_line(&lines, &mut cursor)?;
        if !end_atom.starts_with("END ATOM") {
            return Err(SdfReadError::Parse(format!(
                "END ATOM line not found on line {line_number}"
            )));
        }
    }

    let mut bonds = Vec::with_capacity(bond_count);
    let mut bond_by_bookmark = BTreeMap::new();
    let mut chirality_possible = false;
    if bond_count != 0 {
        let (begin_bond, _) = get_v3000_line(&lines, &mut cursor)?;
        if !begin_bond.starts_with("BEGIN BOND") {
            return Err(SdfReadError::Parse("BEGIN BOND line not found".to_owned()));
        }
        for bond_index in 0..bond_count {
            let (bond_line, line_number) = get_v3000_line(&lines, &mut cursor)?;
            let tokens = tokenize_v3000_line(bond_line.trim());
            if tokens.len() < 4 {
                return Err(SdfReadError::Parse(format!(
                    "bond line {line_number} is too short"
                )));
            }
            let bookmark = tokens[0].parse::<u32>().unwrap_or(0);
            let bond_type = tokens[1].parse::<u32>().unwrap_or(0);
            // BEGIN RDKIT CPP FUNCTION ParseV3000BondBlock (bond type switch)
            // RDKit❗✔️:       case 1:
            // RDKit❗✔️:         bond = new Bond(Bond::SINGLE);
            // RDKit❗✔️:         break;
            // RDKit❗✔️:       case 2:
            // RDKit❗✔️:         bond = new Bond(Bond::DOUBLE);
            // RDKit❗✔️:         break;
            // RDKit❗✔️:       case 3:
            // RDKit❗✔️:         bond = new Bond(Bond::TRIPLE);
            // RDKit❗✔️:         break;
            // RDKit❗✔️:       case 4:
            // RDKit❗✔️:         bond = new Bond(Bond::AROMATIC);
            // RDKit❗✔️:         bond->setIsAromatic(true);
            // RDKit❗✔️:         break;
            // RDKit❗✔️:       case 9:
            // RDKit❗✔️:         bond = new Bond(Bond::DATIVE);
            // RDKit❗✔️:         break;
            // RDKit❗✔️:       case 10:
            // RDKit❗✔️:         bond = new Bond(Bond::HYDROGEN);
            // RDKit❗✔️:         break;
            // RDKit❗✔️:       case 0:
            // RDKit❗✔️:         bond = new Bond(Bond::UNSPECIFIED);
            let (order, query) = match bond_type {
                0 => (BondOrder::Unspecified, None),
                1 => (BondOrder::Single, None),
                2 => (BondOrder::Double, None),
                3 => (BondOrder::Triple, None),
                4 => (BondOrder::Aromatic, None),
                9 => (BondOrder::Dative, None),
                10 => (BondOrder::Hydrogen, None),
                5 => (
                    BondOrder::Unspecified,
                    Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                        BondOrder::Single,
                        BondOrder::Double,
                    ]))),
                ),
                6 => (
                    BondOrder::Unspecified,
                    Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                        BondOrder::Single,
                        BondOrder::Aromatic,
                    ]))),
                ),
                7 => (
                    BondOrder::Unspecified,
                    Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                        BondOrder::Double,
                        BondOrder::Aromatic,
                    ]))),
                ),
                8 => (
                    BondOrder::Unspecified,
                    Some(QueryNode::predicate(BondQueryPredicate::Any)),
                ),
                _ => (
                    BondOrder::Unspecified,
                    Some(QueryNode::predicate(BondQueryPredicate::Any)),
                ),
            };
            let begin_bookmark = tokens[2].parse::<u32>().unwrap_or(0);
            let end_bookmark = tokens[3].parse::<u32>().unwrap_or(0);
            let begin = atom_by_bookmark
                .get(&begin_bookmark)
                .copied()
                .ok_or_else(|| SdfReadError::Field {
                    kind: "bond begin atom bookmark",
                    line: line_number,
                    value: tokens[2].to_owned(),
                })?;
            let end = atom_by_bookmark
                .get(&end_bookmark)
                .copied()
                .ok_or_else(|| SdfReadError::Field {
                    kind: "bond end atom bookmark",
                    line: line_number,
                    value: tokens[3].to_owned(),
                })?;
            let mut spec = BondSpec::new(begin, end, order)
                .with_prop("_MolFileBondType", bond_type.to_string())?;
            if order == BondOrder::Aromatic {
                spec = spec.with_aromatic(true);
            }
            if matches!(bond_type, 5..=7) {
                spec = spec.with_prop("_MolFileBondQuery", "1")?;
            }
            let (mut spec, query, bond_chirality) =
                parse_v3000_bond_properties(spec, query, bond_type, &tokens[4..], line_number)?;
            chirality_possible |= bond_chirality;
            if spec.prop("molStereoCare").is_none()
                && let (Some(begin_care), Some(end_care)) = (
                    atoms[begin.index()].0.prop("molStereoCare"),
                    atoms[end.index()].0.prop("molStereoCare"),
                )
                && begin_care == end_care
            {
                spec = spec.with_prop("molStereoCare", begin_care)?;
            }
            let bond_id = BondId::new(bond_index);
            if bond_by_bookmark.insert(bookmark, bond_id).is_some() {
                return Err(SdfReadError::Parse(format!(
                    "duplicate V3000 bond index {bookmark} on line {line_number}"
                )));
            }
            bonds.push(ParsedV2000Bond {
                bond: Bond::from_spec(bond_id, spec),
                query,
            });
        }
        let (end_bond, line_number) = get_v3000_line(&lines, &mut cursor)?;
        if !end_bond.starts_with("END BOND") {
            return Err(SdfReadError::Parse(format!(
                "END BOND line not found at line {line_number}"
            )));
        }
    }

    let mut link_nodes = Vec::new();
    let (mut trailing, mut trailing_line) = get_v3000_line(&lines, &mut cursor)?;
    loop {
        let upper = trailing.to_ascii_uppercase();
        if !upper.starts_with("LINKNODE") {
            break;
        }
        if upper.len() > 9 {
            link_nodes.push(upper[9..].to_owned());
        }
        (trailing, trailing_line) = get_v3000_line(&lines, &mut cursor)?;
    }
    // BEGIN RDKIT CPP FUNCTION ParseV3000CTAB (trailing typed blocks)
    // RDKit❗✔️: bool sgroupFound = false;
    // RDKit❗✔️: bool obj3dFound = false;
    // RDKit❗✔️: while (tempStr.length() > 5 && tempStr.substr(0, 5) == "BEGIN") {
    // RDKit❗✔️:   if (tempStr.length() >= 12 && tempStr.substr(0, 12) == "BEGIN SGROUP") {
    // RDKit❗✔️:     if (sgroupFound) { throw FileParseException(errout.str()); }
    // RDKit❗✔️:     else if (!nSgroups && strictParsing) { throw FileParseException(errout.str()); }
    // RDKit❗✔️:     sgroupFound = true;
    // RDKit❗✔️:     tempStr = ParseV3000SGroupsBlock(inStream, line, nSgroups, mol, strictParsing);
    // RDKit❗✔️:   } else if (tempStr.length() >= 16 &&
    // RDKit❗✔️:              tempStr.substr(0, 16) == "BEGIN COLLECTION") {
    // RDKit❗✔️:     tempStr = parseEnhancedStereo(inStream, line, mol, strictParsing);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let mut substance_groups = Vec::new();
    let mut stereo_groups = Vec::new();
    let mut sgroup_found = false;
    let mut collection_found = false;
    loop {
        let trailing_upper = trailing.to_ascii_uppercase();
        if !trailing_upper.starts_with("BEGIN") {
            break;
        }
        if trailing_upper.starts_with("BEGIN SGROUP") {
            if sgroup_found {
                return Err(SdfReadError::Parse(format!(
                    "BEGIN SGROUP found more than once on line {trailing_line}"
                )));
            }
            if sgroup_count == 0 && params.strict_parsing {
                return Err(SdfReadError::Parse(format!(
                    "BEGIN SGROUP found but Sgroups NOT expected on line {trailing_line}"
                )));
            }
            sgroup_found = true;
            let expected_sgroups = if sgroup_count == 0 {
                usize::MAX
            } else {
                sgroup_count
            };
            substance_groups = crate::sdf_sgroups::parse_v3000_sgroup_block(
                &lines,
                &mut cursor,
                expected_sgroups,
                &atom_by_bookmark,
                &bond_by_bookmark,
                params.strict_parsing,
            )?;
        } else if trailing_upper.starts_with("BEGIN COLLECTION") {
            if collection_found {
                return Err(SdfReadError::Parse(format!(
                    "BEGIN COLLECTION found more than once on line {trailing_line}"
                )));
            }
            collection_found = true;
            stereo_groups = crate::sdf_sgroups::parse_v3000_collection_block(
                &lines,
                &mut cursor,
                &atom_by_bookmark,
            )?;
        } else if trailing_upper.starts_with("BEGIN OBJ3D") {
            return Err(SdfReadError::Unsupported(
                "V3000 OBJ3D constraints are not yet lowered into detached coordinates",
            ));
        } else {
            return Err(SdfReadError::Parse(format!(
                "Unrecognized V3000 block '{trailing}' on line {trailing_line}"
            )));
        }
        (trailing, trailing_line) = get_v3000_line(&lines, &mut cursor)?;
    }
    if sgroup_count != 0 && !sgroup_found && params.strict_parsing {
        return Err(SdfReadError::Parse(format!(
            "BEGIN SGROUP line not found on line {trailing_line}"
        )));
    }
    let trailing_upper = trailing.to_ascii_uppercase();
    if !trailing_upper.starts_with("END CTAB") {
        return Err(SdfReadError::Parse(format!(
            "END CTAB line not found on line {trailing_line}"
        )));
    }
    let mend_line = cursor + 1;
    let mend = lines
        .get(cursor)
        .copied()
        .map(strip_terminal_cr)
        .unwrap_or_default();
    if !mend.starts_with("M  END") {
        return Err(SdfReadError::Parse(format!(
            "M  END line not found on line {mend_line}"
        )));
    }
    // END RDKIT CPP FUNCTION

    let is_3d = v3000_is_3d(info, &coordinates_3d, chirality_possible);
    let coordinates = CoordinateBlock {
        conformers_2d: if is_3d {
            Vec::new()
        } else {
            vec![Conformer2D::new(
                0,
                coordinates_3d
                    .iter()
                    .map(|point| [point[0], point[1]])
                    .collect(),
            )]
        },
        conformers_3d: if is_3d {
            vec![Conformer3D::new(0, coordinates_3d, true)]
        } else {
            Vec::new()
        },
        source_coordinate_dim: Some(if is_3d {
            CoordinateDimension::ThreeD
        } else {
            CoordinateDimension::TwoD
        }),
    };
    coordinates.validate_for_atom_count(atom_count)?;
    let mut properties = if title.is_empty() {
        MoleculeProperties::default()
    } else {
        MoleculeProperties::default().with_name(title)
    };
    properties = properties
        .with_prop("_MolFileInfo", info)?
        .with_prop("_MolFileComments", comments)?
        .with_prop("_MolFileChiralFlag", chiral_flag.to_string())?;
    if !link_nodes.is_empty() {
        // COSMolKit keeps the established cross-format key used by detached
        // CXSMILES and the live-molecule runtime for RDKit's `_molLinkNodes`.
        properties = properties.with_prop("_MolFileLinkNodes", link_nodes.join("|"))?;
    }
    let has_query = atoms.iter().any(|(_, query)| query.is_some())
        || bonds.iter().any(|bond| bond.query.is_some());
    if has_query {
        let query_atoms = atoms
            .into_iter()
            .map(|(atom, query)| match query {
                Some(predicate) => QueryAtom::from_parts(atom, predicate),
                None => {
                    let atomic_number = atom.element().atomic_number();
                    QueryAtom::from_parts(
                        atom,
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atomic_number)),
                    )
                }
            })
            .collect::<Vec<_>>();
        let query_bonds = bonds
            .into_iter()
            .map(|parsed| match parsed.query {
                Some(predicate) => QueryBond::from_parts(parsed.bond, predicate),
                None => {
                    let predicate = if parsed.bond.order() == BondOrder::Unspecified {
                        QueryNode::predicate(BondQueryPredicate::Any)
                    } else {
                        QueryNode::predicate(BondQueryPredicate::Order(parsed.bond.order()))
                    };
                    QueryBond::from_parts(parsed.bond, predicate)
                }
            })
            .collect::<Vec<_>>();
        let mut query_props = properties.props().clone();
        if let Some(name) = properties.name() {
            query_props.insert("_Name".to_owned(), name.to_owned());
        }
        let query = QueryGraph::from_parts(
            query_atoms,
            query_bonds,
            query_props,
            coordinates.conformers_2d.clone(),
            coordinates.conformers_3d.clone(),
            stereo_groups,
        )?;
        return Ok(MolBlockRecord::Query(QueryMolBlockRecord {
            query,
            substance_groups,
            properties,
            source_coordinate_dim: coordinates.source_coordinate_dim,
        }));
    }

    let atoms = atoms.into_iter().map(|(atom, _)| atom).collect::<Vec<_>>();
    let bonds = bonds
        .into_iter()
        .map(|parsed| parsed.bond)
        .collect::<Vec<_>>();
    let topology = TopologyBlock {
        adjacency: AdjacencyList::from_topology(atom_count, &bonds),
        atoms,
        bonds,
        substance_groups,
        stereo_groups,
    };
    topology.validate()?;
    Ok(MolBlockRecord::Concrete {
        topology,
        coordinates,
        properties,
    })
}

/// Read a concrete V3000 mol block into detached model values.
pub fn read_v3000_detached(
    block: &str,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), SdfReadError> {
    read_v3000_detached_with_params(block, MolBlockReadParams::default())
}

/// Read a concrete V3000 MolBlock with explicit detached syntax policy.
pub fn read_v3000_detached_with_params(
    block: &str,
    params: MolBlockReadParams,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), SdfReadError> {
    match read_v3000_record_detached(block, params)? {
        MolBlockRecord::Concrete {
            topology,
            coordinates,
            properties,
        } => Ok((topology, coordinates, properties)),
        MolBlockRecord::Query(_) => Err(SdfReadError::Unsupported(
            "query-bearing MolBlock; use read_mol_block_detached",
        )),
    }
}

/// Read the first detached record in an SDF text block.
pub fn read_sdf_graph_record_detached(block: &str) -> Result<SdfGraphRecord, SdfReadError> {
    read_sdf_graph_record_detached_with_params(block, SdfDataReadParams::default())
}

/// Read the first detached record with explicit SDF data-field policy.
pub fn read_sdf_graph_record_detached_with_params(
    block: &str,
    params: SdfDataReadParams,
) -> Result<SdfGraphRecord, SdfReadError> {
    let mol_text = block.split("M  END").next().ok_or(SdfReadError::Empty)?;
    let mut mol_block =
        read_mol_block_detached_with_params(&format!("{mol_text}M  END\n"), params.into())?;
    let data_fields = parse_sdf_data_fields(&sdf_data_lines(block)?, params)?;
    match &mut mol_block {
        MolBlockRecord::Concrete { properties, .. } => {
            for (name, value) in &data_fields {
                *properties = properties
                    .clone()
                    .with_prop(name, value)?
                    .with_sdf_data_field(name, value);
            }
        }
        MolBlockRecord::Query(record) => {
            for (name, value) in &data_fields {
                record.query.set_prop(name, value);
                record.properties = record
                    .properties
                    .clone()
                    .with_prop(name, value)?
                    .with_sdf_data_field(name, value);
            }
        }
    }
    if params.process_property_lists {
        apply_sdf_property_lists(&mut mol_block, &data_fields);
    }
    Ok(SdfGraphRecord {
        mol_block,
        data_fields,
    })
}

/// Read the first concrete detached record in an SDF text block.
pub fn read_sdf_record_detached(block: &str) -> Result<SdfRecord, SdfReadError> {
    read_sdf_record_detached_with_params(block, SdfDataReadParams::default())
}

/// Read the first concrete detached record with explicit SDF data-field policy.
pub fn read_sdf_record_detached_with_params(
    block: &str,
    params: SdfDataReadParams,
) -> Result<SdfRecord, SdfReadError> {
    let record = read_sdf_graph_record_detached_with_params(block, params)?;
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        properties,
    } = record.mol_block
    else {
        return Err(SdfReadError::Unsupported(
            "query-bearing SDF record; use read_sdf_graph_record_detached",
        ));
    };
    Ok(SdfRecord {
        topology,
        coordinates,
        properties,
        data_fields: record.data_fields,
    })
}

/// Read every non-empty record in an SDF stream into detached values.
pub fn read_sdf_records_detached(block: &str) -> Result<Vec<SdfRecord>, SdfReadError> {
    read_sdf_records_detached_with_params(block, SdfDataReadParams::default())
}

/// Read every non-empty record with explicit SDF data-field policy.
pub fn read_sdf_records_detached_with_params(
    block: &str,
    params: SdfDataReadParams,
) -> Result<Vec<SdfRecord>, SdfReadError> {
    let records = sdf_record_slices(block)
        .into_iter()
        .map(|record| read_sdf_record_detached_with_params(record, params))
        .collect::<Result<Vec<_>, _>>()?;
    if records.is_empty() {
        return Err(SdfReadError::Empty);
    }
    Ok(records)
}

fn v2000_writer_atom_symbol(atom: &Atom) -> Result<&str, SdfWriteError> {
    if atom.element() != Element::DUMMY {
        return Ok(atom.element().symbol());
    }
    match atom.prop("dummyLabel") {
        Some(label @ ("Pol" | "Mod")) => Ok(label),
        _ => Err(SdfWriteError::Atom(
            "dummy/query atoms require query-aware V2000 serialization",
        )),
    }
}

fn atom_int_prop(atom: &Atom, names: &[&str]) -> i32 {
    names
        .iter()
        .find_map(|name| atom.prop(name).and_then(|value| value.parse().ok()))
        .unwrap_or(0)
}

fn v2000_writer_atom_line(atom: &Atom, coordinate: [f64; 3]) -> Result<String, SdfWriteError> {
    // BEGIN RDKIT CPP FUNCTION GetMolFileAtomProperties / GetMolFileAtomLine
    // RDKit✔️✔️: totValence = 0;
    // RDKit✔️✔️: atomMapNumber = 0;
    // RDKit❗✔️: parityFlag = 0;
    // RDKit✔️✔️: x = y = z = 0.0;
    // RDKit✔️✔️: snprintf(dest, 128,
    // RDKit✔️✔️:          "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d", x, y,
    // RDKit✔️✔️:          z, symbol.c_str(), massDiff, chg, parityFlag, hCount, stereoCare,
    // RDKit✔️✔️:          totValence, rxnComponentType, rxnComponentNumber, atomMapNumber,
    // RDKit✔️✔️:          inversionFlag, exactChangeFlag);
    // Raw `molParity` is retained when no detached stereochemical
    // post-processing has replaced it; geometric parity generation remains a
    // separate IO gap, so only that source line is marked behavior-partial.
    let symbol = v2000_writer_atom_symbol(atom)?;
    Ok(format!(
        "{:>10.4}{:>10.4}{:>10.4} {symbol:>3}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}  0{:>3}{:>3}{:>3}{:>3}{:>3}",
        coordinate[0],
        coordinate[1],
        coordinate[2],
        0,
        0,
        atom.mol_parity().unwrap_or(0),
        atom_int_prop(atom, &["_MolFileHCount"]),
        atom_int_prop(atom, &["molStereoCare", "_MolFileStereoCare"]),
        atom_int_prop(atom, &["molTotValence"]),
        atom_int_prop(atom, &["molRxnRole"]),
        atom_int_prop(atom, &["molRxnComponent"]),
        atom.atom_map().unwrap_or(0),
        atom.mol_inversion_flag().unwrap_or(0),
        atom_int_prop(atom, &["molRxnExactChange"]),
    ))
    // END RDKIT CPP FUNCTION
}

fn v2000_writer_bond_type(bond: &Bond) -> Result<u32, SdfWriteError> {
    match bond.order() {
        BondOrder::Single => Ok(if bond.is_aromatic() { 4 } else { 1 }),
        BondOrder::Double => Ok(if bond.is_aromatic() { 4 } else { 2 }),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        BondOrder::Dative => Ok(9),
        BondOrder::Unspecified => Ok(0),
        other => Err(SdfWriteError::BondOrder(other)),
    }
}

fn v2000_writer_bond_line(bond: &Bond) -> Result<String, SdfWriteError> {
    // BEGIN RDKIT CPP FUNCTION BondGetMolFileSymbol / GetMolFileBondLine
    // RDKit✔️✔️: ss << std::setw(3) << bond->getBeginAtomIdx() + 1;
    // RDKit✔️✔️: ss << std::setw(3) << bond->getEndAtomIdx() + 1;
    // RDKit✔️✔️: ss << std::setw(3) << symbol;
    // RDKit✔️✔️: ss << " " << std::setw(2) << dirCode;
    let direction = match bond.direction() {
        BondDirection::None => {
            if bond.order() == BondOrder::Double && bond.stereo() == BondStereo::Any {
                3
            } else {
                0
            }
        }
        BondDirection::BeginWedge => 1,
        BondDirection::BeginDash => 6,
        BondDirection::EitherDouble => 3,
        BondDirection::Unknown => 4,
        BondDirection::EndUpRight | BondDirection::EndDownRight => 0,
    };
    Ok(format!(
        "{:>3}{:>3}{:>3} {:>2}  0  0  0",
        bond.begin().index() + 1,
        bond.end().index() + 1,
        v2000_writer_bond_type(bond)?,
        direction,
    ))
    // END RDKIT CPP FUNCTION
}

fn append_v2000_counted_property(output: &mut String, label: &str, entries: &[(usize, i32)]) {
    for chunk in entries.chunks(8) {
        output.push_str(&format!("M  {label}{:>3}", chunk.len()));
        for (index, value) in chunk {
            output.push_str(&format!(" {index:>3} {value:>3}"));
        }
        output.push('\n');
    }
}

fn append_v2000_atom_properties(output: &mut String, topology: &TopologyBlock) {
    // BEGIN RDKIT CPP FUNCTION GetMolFileChargeInfo
    // RDKit✔️✔️: if (atom->getFormalCharge() != 0) {
    // RDKit✔️✔️:   ++nChgs;
    // RDKit✔️✔️:   chgss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
    // RDKit✔️✔️:                atom->getFormalCharge();
    let charges = topology
        .atoms
        .iter()
        .filter(|atom| atom.formal_charge() != 0)
        .map(|atom| (atom.id().index() + 1, i32::from(atom.formal_charge())))
        .collect::<Vec<_>>();
    append_v2000_counted_property(output, "CHG", &charges);

    // RDKit✔️✔️: unsigned int nRadEs = atom->getNumRadicalElectrons();
    // RDKit✔️✔️: if (nRadEs != 0 && atom->getTotalDegree() != 0) {
    // RDKit✔️✔️:   if (nRadEs % 2) {
    // RDKit✔️✔️:     nRadEs = 2;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     nRadEs = 3;
    // RDKit✔️✔️:   }
    let radicals = topology
        .atoms
        .iter()
        .filter(|atom| {
            atom.radical_electrons() != 0
                && !topology
                    .adjacency
                    .neighbors_of(atom.id().index())
                    .is_empty()
        })
        .map(|atom| {
            (
                atom.id().index() + 1,
                if atom.radical_electrons() % 2 == 1 {
                    2
                } else {
                    3
                },
            )
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(output, "RAD", &radicals);

    // RDKit✔️✔️: int isotope = atom->getIsotope();
    // RDKit✔️✔️: if (isotope != 0) {
    // RDKit✔️✔️:   ++nMassDiffs;
    // RDKit✔️✔️:         massdiffss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
    // RDKit✔️✔️:                           isotope;
    let isotopes = topology
        .atoms
        .iter()
        .filter_map(|atom| {
            atom.isotope()
                .map(|isotope| (atom.id().index() + 1, i32::from(isotope)))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(output, "ISO", &isotopes);
    // END RDKIT CPP FUNCTION
}

pub fn write_v2000_detached(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<String, SdfWriteError> {
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    if topology.atoms.len() > 999
        || topology.bonds.len() > 999
        || topology.substance_groups.len() > 999
    {
        return Err(SdfWriteError::CountLimit);
    }
    if !topology.stereo_groups.is_empty() {
        return Err(SdfWriteError::SubstanceGroup(
            "V2000 cannot encode enhanced stereo groups; use V3000".to_owned(),
        ));
    }
    let title = properties.name().unwrap_or_default();
    let coords = coordinates
        .conformers_3d
        .first()
        .map(|c| {
            c.coordinates()
                .iter()
                .map(|p| [p[0], p[1], p[2]])
                .collect::<Vec<_>>()
        })
        .or_else(|| {
            coordinates.conformers_2d.first().map(|c| {
                c.coordinates()
                    .iter()
                    .map(|p| [p[0], p[1], 0.0])
                    .collect::<Vec<_>>()
            })
        });
    let coords = coords.unwrap_or_else(|| vec![[0.0, 0.0, 0.0]; topology.atoms.len()]);
    // BEGIN RDKIT CPP FUNCTION outputMolToMolBlock V2000 header/body
    // RDKit✔️✔️: if (tmol.getPropIfPresent(common_properties::_Name, text)) {
    // RDKit✔️✔️:   res += text;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: res += "\n";
    let info = properties.prop("_MolFileInfo").unwrap_or("  COSMolKit");
    let comments = properties.prop("_MolFileComments").unwrap_or_default();
    let chiral_flag = properties
        .prop("_MolFileChiralFlag")
        .and_then(|value| value.parse::<u32>().ok())
        .unwrap_or(0);
    let mut output = format!(
        "{title}\n{info}\n{comments}\n{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}999 V2000\n",
        topology.atoms.len(),
        topology.bonds.len(),
        0,
        topology.substance_groups.len(),
        chiral_flag,
        0,
        0,
        0,
        0,
        0,
    );
    for (atom, point) in topology.atoms.iter().zip(coords) {
        output.push_str(&v2000_writer_atom_line(atom, point)?);
        output.push('\n');
    }
    for bond in &topology.bonds {
        output.push_str(&v2000_writer_bond_line(bond)?);
        output.push('\n');
    }
    // RDKit✔️✔️: res += GetMolFileChargeInfo(tmol);
    append_v2000_atom_properties(&mut output, topology);
    output.push_str(&crate::sdf_sgroups::write_v2000_sgroups(topology)?);
    // RDKit✔️✔️: res += "M  END\n";
    output.push_str("M  END\n");
    // END RDKIT CPP FUNCTION
    Ok(output)
}

/// Write one detached SDF record, including data fields held by the model.
pub fn write_sdf_record_detached(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<String, SdfWriteError> {
    let mut output = write_v2000_detached(topology, coordinates, properties)?;
    for (name, value) in properties.sdf_data_fields() {
        // BEGIN RDKIT CPP FUNCTION _writePropToStream
        // RDKit✔️✔️: if (name.find("\n") != std::string::npos) {
        // RDKit✔️✔️:   return;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: if (pval.find("\r\n\r\n") != std::string::npos ||
        // RDKit✔️✔️:     pval.find("\n\n") != std::string::npos) {
        // RDKit✔️✔️:   return;
        // RDKit✔️✔️: }
        if name.contains('\n') || value.contains("\r\n\r\n") || value.contains("\n\n") {
            continue;
        }
        // RDKit✔️✔️: (*dp_ostream) << ">  <" << name << ">  ";
        // RDKit✔️✔️: (*dp_ostream) << "\n";
        // RDKit✔️✔️: (*dp_ostream) << pval << "\n";
        // RDKit✔️✔️: (*dp_ostream) << "\n";
        output.push_str(&format!(">  <{name}>  \n{value}\n\n"));
        // END RDKIT CPP FUNCTION
    }
    // RDKit✔️✔️: (*dp_ostream) << "$$$$\n";
    output.push_str("$$$$\n");
    Ok(output)
}

fn v3000_writer_atom_symbol(atom: &Atom) -> Result<&str, SdfWriteError> {
    if atom.element() != Element::DUMMY {
        return Ok(atom.element().symbol());
    }
    match atom.prop("dummyLabel") {
        Some(label @ ("Pol" | "Mod")) => Ok(label),
        _ => Err(SdfWriteError::Atom(
            "dummy/query atoms require query-aware V3000 serialization",
        )),
    }
}

fn append_v3000_atom_int_prop(output: &mut String, atom: &Atom, key: &str, label: &str) {
    if let Some(value) = atom.prop(key)
        && value != "0"
    {
        output.push_str(&format!(" {label}={value}"));
    }
}

fn v3000_writer_atom_line(
    topology: &TopologyBlock,
    atom: &Atom,
    coordinate: [f64; 3],
) -> Result<String, SdfWriteError> {
    // BEGIN RDKIT CPP FUNCTION GetV3000MolFileAtomLine
    // RDKit❗✔️:   ss << "M  V30 " << atom->getIdx() + 1;
    // RDKit❗✔️:   std::string symbol = AtomGetMolFileSymbol(atom, false, queryListAtoms);
    // RDKit❗✔️:   if (!isAtomListQuery(atom) || queryListAtoms[atom->getIdx()]) {
    // RDKit❗✔️:     ss << " " << symbol;
    // RDKit❗✔️:   } else {
    // RDKit❌❌:     INT_VECT vals;
    // RDKit❌❌:     getAtomListQueryVals(atom->getQuery(), vals);
    // RDKit❗✔️:   }
    // RDKit✔️✔️:   ss << std::fixed;
    // RDKit✔️✔️:   ss << std::setprecision(precision);
    // RDKit✔️✔️:   ss << " " << x << " " << y << " " << z;
    // RDKit✔️✔️:   ss << std::setprecision(currentPrecision);
    // RDKit✔️✔️:   ss << std::defaultfloat;
    // RDKit✔️✔️:   ss << " " << atomMapNumber;
    // RDKit❗✔️:   if (parityFlag != 0) {
    // RDKit❗✔️:     ss << " CFG=" << parityFlag;
    // RDKit❗✔️:   }
    // RDKit✔️✔️:   if (chg != 0) {
    // RDKit✔️✔️:     ss << " CHG=" << chg;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (isotope != 0 && !isAtomRGroup(*atom)) {
    // RDKit✔️✔️:     int mass = static_cast<int>(std::round(atom->getMass()));
    // RDKit✔️✔️:     if (!mass) {
    // RDKit✔️✔️:       mass = isotope;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ss << " MASS=" << mass;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int nRadEs = atom->getNumRadicalElectrons();
    // RDKit✔️✔️:   if (nRadEs != 0 && atom->getTotalDegree() != 0) {
    // RDKit✔️✔️:     if (nRadEs % 2) {
    // RDKit✔️✔️:       nRadEs = 2;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       nRadEs = 3;  // we use triplets, not singlets:
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ss << " RAD=" << nRadEs;
    // RDKit✔️✔️:   }
    // RDKit❗✔️:   if (totValence != 0) {
    // RDKit❗✔️:     if (totValence == 15) {
    // RDKit❗✔️:       ss << " VAL=-1";
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       ss << " VAL=" << totValence;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molAttachOrder, iprop) &&
    // RDKit✔️✔️:         iprop) {
    // RDKit✔️✔️:       ss << " ATTCHORD=" << iprop;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molAttachPoint, iprop) &&
    // RDKit✔️✔️:         iprop) {
    // RDKit✔️✔️:       ss << " ATTCHPT=" << iprop;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molAtomSeqId, iprop) &&
    // RDKit✔️✔️:         iprop) {
    // RDKit✔️✔️:       ss << " SEQID=" << iprop;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:       if (atom->getPropIfPresent(common_properties::molAtomSeqName, sprop)) {
    // RDKit✔️✔️:         ss << " SEQNAME=" << sprop;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molRxnExactChange, iprop) &&
    // RDKit✔️✔️:         iprop) {
    // RDKit✔️✔️:       ss << " EXACHG=" << iprop;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molInversionFlag, iprop) &&
    // RDKit✔️✔️:         iprop) {
    // RDKit✔️✔️:       if (iprop == 1 || iprop == 2) {
    // RDKit✔️✔️:         ss << " INVRET=" << iprop;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molStereoCare, iprop) &&
    // RDKit✔️✔️:         iprop) {
    // RDKit✔️✔️:       ss << " STBOX=" << iprop;
    // RDKit✔️✔️:     }
    // RDKit❌❌:     if (atom->getPropIfPresent(common_properties::molSubstCount, iprop) &&
    // RDKit❌❌:         iprop) {
    // RDKit❌❌:       ss << " SUBST=" << iprop;
    // RDKit❌❌:     }
    // RDKit❌❌:     if (atom->getPropIfPresent(common_properties::molRingBondCount, iprop) &&
    // RDKit❌❌:         iprop) {
    // RDKit❌❌:       ss << " RBCNT=" << iprop;
    // RDKit❌❌:     }
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molAtomClass, sprop)) {
    // RDKit✔️✔️:       ss << " CLASS=" << sprop;
    // RDKit✔️✔️:     }
    // The detached writer preserves parsed parity and total-valence fields; it
    // cannot derive writer-time parity or non-default valence without the
    // runtime geometry/valence services, so those source branches remain
    // partial while their persisted representation is round-trippable.
    let symbol = v3000_writer_atom_symbol(atom)?;
    let mut output = format!(
        "M  V30 {} {} {:.6} {:.6} {:.6} {}",
        atom.id().index() + 1,
        symbol,
        coordinate[0],
        coordinate[1],
        coordinate[2],
        atom.atom_map().unwrap_or(0)
    );
    if let Some(parity) = atom.mol_parity()
        && parity != 0
    {
        output.push_str(&format!(" CFG={parity}"));
    }
    if atom.formal_charge() != 0 {
        output.push_str(&format!(" CHG={}", atom.formal_charge()));
    }
    if let Some(isotope) = atom.isotope() {
        output.push_str(&format!(" MASS={isotope}"));
    }
    let radical_electrons = atom.radical_electrons();
    if radical_electrons != 0
        && !topology
            .adjacency
            .neighbors_of(atom.id().index())
            .is_empty()
    {
        output.push_str(&format!(
            " RAD={}",
            if radical_electrons % 2 == 1 { 2 } else { 3 }
        ));
    }
    if let Some(total_valence) = atom.prop("molTotValence")
        && total_valence != "0"
    {
        output.push_str(" VAL=");
        output.push_str(if total_valence == "15" {
            "-1"
        } else {
            total_valence
        });
    }
    append_v3000_atom_int_prop(&mut output, atom, "molAttachOrder", "ATTCHORD");
    append_v3000_atom_int_prop(&mut output, atom, "molAttachPoint", "ATTCHPT");
    append_v3000_atom_int_prop(&mut output, atom, "molAtomSeqId", "SEQID");
    if let Some(value) = atom.prop("molAtomSeqName") {
        output.push_str(&format!(" SEQNAME={value}"));
    }
    append_v3000_atom_int_prop(&mut output, atom, "molRxnExactChange", "EXACHG");
    if let Some(value) = atom.mol_inversion_flag()
        && matches!(value, 1 | 2)
    {
        output.push_str(&format!(" INVRET={value}"));
    }
    append_v3000_atom_int_prop(&mut output, atom, "molStereoCare", "STBOX");
    if atom.prop("molSubstCount").is_some_and(|value| value != "0")
        || atom
            .prop("molRingBondCount")
            .is_some_and(|value| value != "0")
    {
        return Err(SdfWriteError::Atom(
            "substitution/ring-bond-count query atoms require query-aware V3000 serialization",
        ));
    }
    if let Some(value) = atom.prop("molAtomClass") {
        output.push_str(&format!(" CLASS={value}"));
    }
    // END RDKIT CPP FUNCTION
    Ok(output)
}

fn v3000_writer_bond_type(bond: &Bond) -> Result<u32, SdfWriteError> {
    match bond.order() {
        BondOrder::Single | BondOrder::Double if bond.is_aromatic() => Ok(4),
        BondOrder::Single => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        BondOrder::Dative => Ok(9),
        BondOrder::Hydrogen => Ok(10),
        BondOrder::Unspecified => Ok(0),
        other => Err(SdfWriteError::BondOrder(other)),
    }
}

fn v3000_writer_bond_line(bond: &Bond) -> Result<String, SdfWriteError> {
    // BEGIN RDKIT CPP FUNCTION GetV3000MolFileBondLine
    // RDKit✔️✔️: ss << "M  V30 " << bond->getIdx() + 1;
    // RDKit✔️✔️: ss << " " << GetV3000BondCode(bond);
    // RDKit✔️✔️: ss << " " << bond->getBeginAtomIdx() + 1;
    // RDKit✔️✔️: ss << " " << bond->getEndAtomIdx() + 1;
    // RDKit❗✔️:   if (dirCode != 0) {
    // RDKit❗✔️:     ss << " CFG=" << BondStereoCodeV2000ToV3000(dirCode);
    // RDKit❗✔️:   }
    // RDKit✔️✔️:     if (bond->getPropIfPresent(common_properties::molReactStatus, iprop) &&
    // RDKit✔️✔️:         iprop) {
    // RDKit✔️✔️:       ss << " RXCTR=" << iprop;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (bond->getPropIfPresent(common_properties::molStereoCare, sprop) &&
    // RDKit✔️✔️:         sprop != "0") {
    // RDKit✔️✔️:       ss << " STBOX=" << sprop;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, sprop) &&
    // RDKit✔️✔️:         sprop != "0") {
    // RDKit✔️✔️:       ss << " ENDPTS=" << sprop;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (bond->getPropIfPresent(common_properties::_MolFileBondAttach, sprop) &&
    // RDKit✔️✔️:         sprop != "0") {
    // RDKit✔️✔️:       ss << " ATTACH=" << sprop;
    // RDKit✔️✔️:     }
    // Direction state is serialized directly. Runtime-only wedge selection and
    // endpoint reversal remain outside this detached writer.
    let bond_type = v3000_writer_bond_type(bond)?;
    let mut output = format!(
        "M  V30 {} {} {} {}",
        bond.id().index() + 1,
        bond_type,
        bond.begin().index() + 1,
        bond.end().index() + 1
    );
    let configuration = match bond.direction() {
        BondDirection::BeginWedge => Some(1),
        BondDirection::Unknown | BondDirection::EitherDouble => Some(2),
        BondDirection::BeginDash => Some(3),
        _ => None,
    };
    if let Some(configuration) = configuration {
        output.push_str(&format!(" CFG={configuration}"));
    }
    for (key, label) in [
        ("molReactStatus", "RXCTR"),
        ("molStereoCare", "STBOX"),
        ("_MolFileBondEndPts", "ENDPTS"),
        ("_MolFileBondAttach", "ATTACH"),
    ] {
        if let Some(value) = bond.prop(key)
            && value != "0"
        {
            output.push_str(&format!(" {label}={value}"));
        }
    }
    // END RDKIT CPP FUNCTION
    Ok(output)
}

/// Write a concrete V3000 mol block from detached model values.
pub fn write_v3000_detached(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<String, SdfWriteError> {
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    let title = properties.name().unwrap_or_default();
    let selected_3d = coordinates.conformers_3d.first();
    let selected_2d = coordinates.conformers_2d.first();
    let points = selected_3d
        .map(|conformer| {
            conformer
                .coordinates()
                .iter()
                .map(|p| [p[0], p[1], p[2]])
                .collect::<Vec<_>>()
        })
        .or_else(|| {
            selected_2d.map(|conformer| {
                conformer
                    .coordinates()
                    .iter()
                    .map(|point| [point[0], point[1], 0.0])
                    .collect::<Vec<_>>()
            })
        })
        .unwrap_or_else(|| vec![[0.0, 0.0, 0.0]; topology.atoms.len()]);
    let default_info = if selected_3d.is_some() {
        "     RDKit          3D"
    } else if selected_2d.is_some() {
        "     RDKit          2D"
    } else {
        "     RDKit"
    };
    let info = properties.prop("_MolFileInfo").unwrap_or(default_info);
    let comments = properties.prop("_MolFileComments").unwrap_or_default();
    let chiral_flag = properties
        .prop("_MolFileChiralFlag")
        .and_then(|value| value.parse::<u32>().ok())
        .unwrap_or(0);
    // BEGIN RDKIT CPP FUNCTION outputMolToMolBlock / getV3000CTAB
    // RDKit✔️✔️:   if (tmol.getPropIfPresent(common_properties::_Name, text)) {
    // RDKit✔️✔️:     res += text;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res += "\n";
    // RDKit✔️✔️:   if (tmol.getPropIfPresent(common_properties::MolFileInfo, text)) {
    // RDKit✔️✔️:     res += text;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::stringstream ss;
    // RDKit✔️✔️:     ss << "  " << std::setw(8) << "RDKit";
    // RDKit✔️✔️:     ss << std::setw(10) << "";
    // RDKit✔️✔️:     if (conf) {
    // RDKit✔️✔️:       if (conf->is3D()) {
    // RDKit✔️✔️:         ss << "3D";
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         ss << common_properties::TWOD;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += ss.str();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res += "\n";
    // RDKit✔️✔️:   if (tmol.getPropIfPresent(common_properties::MolFileComments, text)) {
    // RDKit✔️✔️:     res += text;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res += "\n";
    // RDKit✔️✔️:   if (isV3000) {
    // RDKit✔️✔️:     // All counts in the V3000 info line should be 0
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << std::setw(3) << 0;
    // RDKit✔️✔️:     ss << "999 V3000\n";
    // RDKit✔️✔️:   std::string res = "M  V30 BEGIN CTAB\n";
    // RDKit✔️✔️:   ss << "M  V30 COUNTS " << nAtoms << " " << nBonds << " " << nSGroups << " "
    // RDKit✔️✔️:      << num3DConstraints << " " << chiralFlag << "\n";
    let mut output = format!(
        "{title}\n{info}\n{comments}\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {} {} {} 0 {}\n",
        topology.atoms.len(),
        topology.bonds.len(),
        topology.substance_groups.len(),
        chiral_flag
    );
    if !topology.atoms.is_empty() {
        // RDKit✔️✔️: res += "M  V30 BEGIN ATOM\n";
        output.push_str("M  V30 BEGIN ATOM\n");
        for (atom, point) in topology.atoms.iter().zip(points) {
            output.push_str(&v3000_writer_atom_line(topology, atom, point)?);
            output.push('\n');
        }
        // RDKit✔️✔️: res += "M  V30 END ATOM\n";
        output.push_str("M  V30 END ATOM\n");
    }
    if !topology.bonds.is_empty() {
        // RDKit✔️✔️: res += "M  V30 BEGIN BOND\n";
        output.push_str("M  V30 BEGIN BOND\n");
        for bond in &topology.bonds {
            output.push_str(&v3000_writer_bond_line(bond)?);
            output.push('\n');
        }
        // RDKit✔️✔️: res += "M  V30 END BOND\n";
        output.push_str("M  V30 END BOND\n");
    }
    if let Some(link_nodes) = properties.prop("_MolFileLinkNodes") {
        for link_node in link_nodes.split('|') {
            output.push_str(&format!("M  V30 LINKNODE {link_node}\n"));
        }
    }
    output.push_str(&crate::sdf_sgroups::write_v3000_typed_blocks(topology));
    // RDKit✔️✔️: res += "M  V30 END CTAB\n";
    // RDKit✔️✔️: res += "M  END\n";
    output.push_str("M  V30 END CTAB\nM  END\n");
    // END RDKIT CPP FUNCTION
    Ok(output)
}

#[cfg(test)]
mod tests {
    use cosmolkit_model::{
        AtomId, AtomQueryPredicate, BondId, BondQueryPredicate, CoordinateBlock,
        MoleculeProperties, QueryNode, SGroupBondRole, SGroupConnection, SdfPropertyListTarget,
        StereoGroupKind, SubstanceGroupKind,
    };
    use cosmolkit_types::{BondDirection, BondStereo};

    use super::{
        MolBlockReadParams, MolBlockRecord, SdfDataReadParams, SdfGraphDataset, SdfGraphReader,
        read_mol_block_detached, read_mol_block_detached_with_params,
        read_sdf_graph_record_detached, read_sdf_record_detached,
        read_sdf_record_detached_with_params, read_sdf_records_detached, read_v2000_detached,
        read_v2000_detached_with_params, read_v3000_detached, read_v3000_detached_with_params,
        write_sdf_record_detached, write_v2000_detached, write_v3000_detached,
    };

    fn v2000_atom_line(
        symbol: &str,
        mass_diff: i32,
        charge: i32,
        parity: i32,
        h_count: i32,
        stereo_care: i32,
        total_valence: i32,
        reaction_role: i32,
        reaction_component: i32,
        atom_map: i32,
        inversion: i32,
        exact_change: i32,
    ) -> String {
        format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3}{mass_diff:>2}{charge:>3}{parity:>3}{h_count:>3}{stereo_care:>3}{total_valence:>3}{:>3}{reaction_role:>3}{reaction_component:>3}{atom_map:>3}{inversion:>3}{exact_change:>3}",
            0.0, 0.0, 0.0, symbol, 0
        )
    }

    #[test]
    fn v2000_detached_reader_and_writer_round_trip_core_graph() {
        let input = "ethanol\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1000    1.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n";
        let (topology, coordinates, properties) = read_v2000_detached(input).expect("read");
        let output = write_v2000_detached(&topology, &coordinates, &properties).expect("write");
        let (roundtrip, _, roundtrip_properties) = read_v2000_detached(&output).expect("roundtrip");
        assert_eq!(roundtrip.atoms.len(), 3);
        assert_eq!(roundtrip.bonds.len(), 2);
        assert_eq!(roundtrip_properties.name(), Some("ethanol"));
    }

    #[test]
    fn molblock_reader_applies_rdkit_strict_version_dispatch() {
        let malformed = concat!(
            "bad-version\n  test\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 X2000\n",
            "M  END\n",
        );
        assert!(
            read_mol_block_detached(malformed)
                .expect_err("strict malformed CTAB version")
                .to_string()
                .contains("CTAB version string invalid")
        );
        let record = read_mol_block_detached_with_params(
            malformed,
            MolBlockReadParams {
                strict_parsing: false,
            },
        )
        .expect("non-strict malformed version falls back to V2000");
        assert!(matches!(record, MolBlockRecord::Concrete { .. }));

        let unsupported = malformed.replace("X2000", "V4000");
        assert!(
            read_mol_block_detached(&unsupported)
                .expect_err("strict unsupported CTAB version")
                .to_string()
                .contains("Unsupported CTAB version: 'V4000'")
        );

        let sdf = format!("{malformed}>  <ID>\nnon-strict\n\n$$$$\n");
        let record = read_sdf_record_detached_with_params(
            &sdf,
            SdfDataReadParams {
                strict_parsing: false,
                ..SdfDataReadParams::default()
            },
        )
        .expect("SDF strictness reaches the embedded MolBlock parser");
        assert_eq!(
            record.data_fields,
            [("ID".to_owned(), "non-strict".to_owned())]
        );
    }

    #[test]
    fn v2000_atom_line_length_and_unknown_symbol_follow_strict_policy() {
        let short_atom = format!("{:>10}{:>10}{:>10} C", "0", "0", "0");
        assert_eq!(short_atom.len(), 32);
        let short_block = format!(
            "short\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{short_atom}\nM  END\n"
        );
        assert!(
            read_v2000_detached(&short_block)
                .expect_err("strict 32-character atom line")
                .to_string()
                .contains("Atom line too short")
        );
        let (short_topology, _, _) = read_v2000_detached_with_params(
            &short_block,
            MolBlockReadParams {
                strict_parsing: false,
            },
        )
        .expect("non-strict 32-character atom line");
        assert_eq!(short_topology.atoms[0].element().symbol(), "C");

        let unknown_atom = v2000_atom_line("Zz", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let unknown_block = format!(
            "unknown\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{unknown_atom}\nM  END\n"
        );
        assert!(
            read_v2000_detached(&unknown_block)
                .expect_err("strict unknown atom symbol")
                .to_string()
                .contains("Element 'Zz' not found")
        );
        let (unknown_topology, _, _) = read_v2000_detached_with_params(
            &unknown_block,
            MolBlockReadParams {
                strict_parsing: false,
            },
        )
        .expect("non-strict unknown symbol becomes a labeled dummy");
        assert_eq!(unknown_topology.atoms[0].element().atomic_number(), 0);
        assert_eq!(unknown_topology.atoms[0].prop("dummyLabel"), Some("Zz"));
    }

    #[test]
    fn v3000_outer_counts_follow_strict_policy() {
        let input = concat!(
            "outer-counts\n  test\n\n",
            "  1  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        assert!(
            read_v3000_detached(input)
                .expect_err("strict nonzero V3000 outer counts")
                .to_string()
                .contains("should have 0s")
        );
        let (topology, _, _) = read_v3000_detached_with_params(
            input,
            MolBlockReadParams {
                strict_parsing: false,
            },
        )
        .expect("non-strict nonzero V3000 outer counts");
        assert_eq!(topology.atoms.len(), 1);
    }

    #[test]
    fn v2000_reader_lowers_pinned_zbo_charge_and_hydrogen_records() {
        let input = include_str!(
            "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/H3BNH3.mol"
        );
        let (topology, _, _) = read_v2000_detached(input).expect("read pinned H3BNH3 fixture");
        assert_eq!(topology.atoms.len(), 2);
        assert_eq!(topology.atoms[0].formal_charge(), 0);
        assert_eq!(topology.atoms[1].formal_charge(), 0);
        assert_eq!(topology.atoms[0].explicit_hydrogens(), 3);
        assert_eq!(topology.atoms[0].prop("_ZBO_H"), Some("1"));
        assert_eq!(topology.bonds[0].order(), cosmolkit_types::BondOrder::Zero);
        assert!(!topology.bonds[0].is_aromatic());
    }

    #[test]
    fn v2000_zbo_nonzero_values_use_rdkit_bond_type_codes() {
        let atom_one = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let atom_two = v2000_atom_line("O", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            concat!(
                "zbo-order\n  test\n\n",
                "  2  1  0  0  0  0  0  0  0  0999 V2000\n",
                "{}\n{}\n",
                "  1  2  1  0  0  0  0\n",
                "M  ZBO  1   1   2\n",
                "M  END\n",
            ),
            atom_one, atom_two,
        );
        let (topology, _, _) = read_v2000_detached(&input).expect("read nonzero ZBO order");
        assert_eq!(
            topology.bonds[0].order(),
            cosmolkit_types::BondOrder::Double
        );

        let aromatic_then_zero = input
            .replace("  1  2  1  0  0  0  0", "  1  2  4  0  0  0  0")
            .replace("M  ZBO  1   1   2", "M  ZBO  1   1   0");
        let (topology, _, _) =
            read_v2000_detached(&aromatic_then_zero).expect("ZBO only changes RDKit bond type");
        assert_eq!(topology.bonds[0].order(), cosmolkit_types::BondOrder::Zero);
        assert!(topology.bonds[0].is_aromatic());
    }

    #[test]
    fn v2000_reader_preserves_alias_value_pxa_apo_link_and_skip_records() {
        let atom_one = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let atom_two = v2000_atom_line("N", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            concat!(
                "legacy-properties\n  test\n\n",
                "  2  1  0  0  0  0  0  0  0  0999 V2000\n",
                "{}\n{}\n",
                "  1  2  1  0  0  0  0\n",
                "A  {:>3}\ncarbon alias\n",
                "V  {:>3} atom value\n",
                "M  PXA {:>3} pxa payload\n",
                "M  APO  2   1   2   2   1\n",
                "M  LIN  1   1   3   2   0\n",
                "S  SKP  1\n",
                "M  CHG  1   1  -1\n",
                "M  END\n",
            ),
            atom_one, atom_two, 1, 2, 1,
        );
        let (topology, _, properties) = read_v2000_detached(&input).expect("read legacy records");
        assert_eq!(topology.atoms[0].prop("molFileAlias"), Some("carbon alias"));
        assert_eq!(topology.atoms[1].prop("molFileValue"), Some("atom value"));
        assert_eq!(topology.atoms[0].prop("_MolFile_PXA"), Some(" pxa payload"));
        assert_eq!(topology.atoms[0].prop("molAttachPoint"), Some("2"));
        assert_eq!(topology.atoms[1].prop("molAttachPoint"), Some("1"));
        assert_eq!(topology.atoms[0].formal_charge(), 0);
        assert_eq!(properties.prop("_MolFileLinkNodes"), Some("1 3 1 1 2"));
    }

    #[test]
    fn v2000_duplicate_attachment_points_follow_strict_policy() {
        let atom = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            concat!(
                "duplicate-apo\n  test\n\n",
                "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
                "{}\n",
                "M  APO  1   1   2\n",
                "M  APO  1   1   1\n",
                "M  END\n",
            ),
            atom,
        );
        assert!(
            read_v2000_detached(&input)
                .expect_err("strict duplicate APO")
                .to_string()
                .contains("Multiple ATTCHPT values")
        );
        let (topology, _, _) = read_v2000_detached_with_params(
            &input,
            MolBlockReadParams {
                strict_parsing: false,
            },
        )
        .expect("non-strict duplicate APO keeps the first value");
        assert_eq!(topology.atoms[0].prop("molAttachPoint"), Some("2"));
    }

    #[test]
    fn v2000_blank_property_lines_follow_strict_policy() {
        let input = concat!(
            "blank-property\n  test\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V2000\n",
            "\n",
            "M  END\n",
        );
        assert!(
            read_v2000_detached(input)
                .expect_err("strict blank property line")
                .to_string()
                .contains("unexpected blank line")
        );
        read_v2000_detached_with_params(
            input,
            MolBlockReadParams {
                strict_parsing: false,
            },
        )
        .expect("non-strict blank property line");
    }

    #[test]
    fn v2000_reader_preserves_fixed_width_atom_and_property_state() {
        let atom_one = v2000_atom_line("C", 1, 3, 1, 0, 1, 4, 2, 3, 12, 1, 1);
        let atom_two = v2000_atom_line("N", 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            "atoms\n  source            2D\ncomment\n  2  1  0  0  1  0  0  0  0  0999 V2000\n{atom_one}\n{atom_two}\n  1  2  2  3  0  0  0\nM  CHG  1   2  -1\nM  RAD  1   1   2\nM  ISO  1   1  13\nM  END\n"
        );
        let (topology, _, properties) = read_v2000_detached(&input).expect("read atom state");
        let carbon = &topology.atoms[0];
        let nitrogen = &topology.atoms[1];
        assert_eq!(carbon.formal_charge(), 0);
        assert_eq!(nitrogen.formal_charge(), -1);
        assert_eq!(carbon.isotope(), Some(13));
        assert_eq!(carbon.radical_electrons(), 1);
        assert_eq!(carbon.mol_parity(), Some(1));
        assert_eq!(carbon.atom_map(), Some(12));
        assert_eq!(carbon.mol_inversion_flag(), Some(1));
        assert_eq!(carbon.prop("molStereoCare"), Some("1"));
        assert_eq!(carbon.prop("molTotValence"), Some("4"));
        assert_eq!(carbon.prop("molRxnRole"), Some("2"));
        assert_eq!(carbon.prop("molRxnComponent"), Some("3"));
        assert_eq!(carbon.prop("molRxnExactChange"), Some("1"));
        assert_eq!(topology.bonds[0].direction(), BondDirection::EitherDouble);
        assert_eq!(topology.bonds[0].stereo(), BondStereo::Any);
        assert_eq!(topology.bonds[0].prop("_MolFileBondType"), Some("2"));
        assert_eq!(topology.bonds[0].prop("_MolFileBondStereo"), Some("3"));
        assert_eq!(
            properties.prop("_MolFileInfo"),
            Some("  source            2D")
        );
        assert_eq!(properties.prop("_MolFileComments"), Some("comment"));
        assert_eq!(properties.prop("_MolFileChiralFlag"), Some("1"));

        let output = write_v2000_detached(&topology, &CoordinateBlock::default(), &properties)
            .expect("write preserved atom state");
        assert!(output.contains("M  CHG  1   2  -1\n"));
        assert!(output.contains("M  RAD  1   1   2\n"));
        assert!(output.contains("M  ISO  1   1  13\n"));
        let (roundtrip, _, roundtrip_properties) =
            read_v2000_detached(&output).expect("roundtrip atom state");
        assert_eq!(roundtrip.atoms[0].radical_electrons(), 1);
        assert_eq!(roundtrip.atoms[0].isotope(), Some(13));
        assert_eq!(roundtrip.atoms[0].mol_parity(), Some(1));
        assert_eq!(roundtrip.atoms[0].atom_map(), Some(12));
        assert_eq!(roundtrip.atoms[1].formal_charge(), -1);
        assert_eq!(roundtrip.bonds[0].direction(), BondDirection::EitherDouble);
        assert_eq!(roundtrip.bonds[0].stereo(), BondStereo::Any);
        assert_eq!(
            roundtrip_properties.prop("_MolFileComments"),
            Some("comment")
        );
    }

    #[test]
    fn v2000_reader_preserves_basic_sgroup_membership() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let oxygen = v2000_atom_line("O", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            "groups\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n{carbon}\n{oxygen}\n  1  2  1  0  0  0  0\nM  STY  1   1 SUP\nM  SLB  1   1   7\nM  SAL   1  2   1   2\nM  SPA   1  1   1\nM  SBL   1  1   1\nM  SMT   1 Me\nM  END\n"
        );
        let (topology, _, _) = read_v2000_detached(&input).expect("read V2000 SGroup");
        let group = &topology.substance_groups[0];
        assert_eq!(group.kind(), &SubstanceGroupKind::Superatom);
        assert_eq!(group.rdkit_sequence_id(), Some(1));
        assert_eq!(group.external_id(), Some(7));
        assert_eq!(group.atoms(), &[AtomId::new(0), AtomId::new(1)]);
        assert_eq!(group.parent_atoms(), &[AtomId::new(0)]);
        assert_eq!(group.bonds(), &[BondId::new(0)]);
        assert_eq!(group.label(), Some("Me"));
    }

    #[test]
    fn v2000_reader_preserves_typed_sgroup_properties() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let oxygen = v2000_atom_line("O", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let sdi = format!(
            "M  SDI   1  4{:>10.4}{:>10.4}{:>10.4}{:>10.4}",
            0.0, 1.0, 2.0, 3.0
        );
        let sbv = format!("M  SBV   1   1{:>10.4}{:>10.4}", 0.5, 0.25);
        let sdt = format!(
            "M  SDT   2 {:<30}{:<2}{:<20}{:<2}{}",
            "FIELD", "T", "INFO", "Q", "OP"
        );
        let input = format!(
            "groups\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n{carbon}\n{oxygen}\n  1  2  1  0  0  0  0\nM  STY  2   2 DAT   1 SUP\nM  SST  1   1 ALT\nM  SCN  1   1 HT\nM  SDS EXP  1   1\nM  SAL   1  1   1\nM  SBL   1  1   1\n{sdi}\n{sbv}\n{sdt}\nM  SDD   2 display spec\nM  SCD   2 first value\nM  SED   2 second value\nM  SPL  1   2   1\nM  SNC  1   2   5\nM  SAP   1  1   1   2 AP\nM  SCL   2 CLASS\nM  SBT  1   2   1\nM  END\n"
        );
        let (topology, _, _) = read_v2000_detached(&input).expect("read typed SGroups");
        let sup = &topology.substance_groups[0];
        assert_eq!(sup.rdkit_sequence_id(), Some(1));
        assert_eq!(sup.subtype(), Some("ALT"));
        assert_eq!(sup.connection(), Some(&SGroupConnection::HeadToTail));
        assert_eq!(sup.expansion_state(), Some("E"));
        assert_eq!(sup.display().unwrap().brackets[0].p1, [0.0, 1.0]);
        assert_eq!(sup.display().unwrap().brackets[0].p2, [2.0, 3.0]);
        assert_eq!(sup.cstates()[0].bond, BondId::new(0));
        assert_eq!(sup.cstates()[0].vector, [0.5, 0.25]);
        assert_eq!(sup.attach_points()[0].atom, AtomId::new(0));
        assert_eq!(sup.attach_points()[0].leaving_atom, Some(AtomId::new(1)));
        assert_eq!(sup.attach_points()[0].label.as_deref(), Some("AP"));

        let dat = &topology.substance_groups[1];
        assert_eq!(dat.rdkit_sequence_id(), Some(2));
        assert_eq!(dat.parent(), Some(sup.id()));
        assert_eq!(dat.component_number(), Some(5));
        assert_eq!(dat.class(), Some("CLASS"));
        assert_eq!(
            dat.bracket_style(),
            Some(&cosmolkit_model::SGroupBracketStyle::Parenthesis)
        );
        let data = dat.data().unwrap();
        assert_eq!(data.field_name.as_deref(), Some("FIELD"));
        assert_eq!(data.field_type.as_deref(), Some("T"));
        assert_eq!(data.field_info.as_deref(), Some("INFO"));
        assert_eq!(data.query_type.as_deref(), Some("Q"));
        assert_eq!(data.query_op.as_deref(), Some("OP"));
        assert_eq!(data.field_display.as_deref(), Some("display spec"));
        assert_eq!(data.values, ["first valuesecond value"]);

        let output = write_v2000_detached(
            &topology,
            &CoordinateBlock::default(),
            &MoleculeProperties::default(),
        )
        .expect("write typed V2000 SGroups");
        assert_eq!(&output.lines().nth(3).unwrap()[9..12], "  2");
        assert!(output.contains("M  STY"));
        assert!(output.contains("M  SDT"));
        let (roundtrip, _, _) =
            read_v2000_detached(&output).expect("roundtrip typed V2000 SGroups");
        assert_eq!(roundtrip.substance_groups.len(), 2);
        assert_eq!(roundtrip.substance_groups[0].subtype(), Some("ALT"));
        assert_eq!(
            roundtrip.substance_groups[1].parent(),
            Some(roundtrip.substance_groups[0].id())
        );
        assert_eq!(
            roundtrip.substance_groups[1].data().unwrap().values,
            ["first valuesecond value"]
        );
    }

    #[test]
    fn v2000_malformed_sgroups_throw_strictly_and_are_dropped_non_strictly() {
        // Pinned RDKit regression: "do not throw but remove malformed V2000
        // SGroups when strictParsing is false".
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let sbv_one = format!("M  SBV   1   1{:>10.4}{:>10.4}", 0.5, 0.25);
        let sbv_two = format!("M  SBV   2   2{:>10.4}{:>10.4}", -0.5, -0.25);
        let input = format!(
            "malformed-groups\n  test            2D\n\n  4  2  0  0  0  0  0  0  0  0999 V2000\n{carbon}\n{carbon}\n{carbon}\n{carbon}\n  1  2  1  0  0  0  0\n  3  4  1  0  0  0  0\nM  STY  2   1 SUP   2 SUP\nM  SAL   1  1   1\nM  SBL   1  1   1\n{sbv_one}\nM  SAL   2  1   3\nM  SBL   2  1   2\n{sbv_two}\nM  END\n"
        );
        let (base, _, _) = read_v2000_detached(&input).expect("valid source-shaped SGroups");
        assert_eq!(base.substance_groups.len(), 2);

        let malformed = [
            input.replace("M  SBL   1  1   1", "M  SBL   3  1   1"),
            input.replace(&sbv_one, "M  SBV   1   1"),
            input.replace("M  SBL   2  1   2", "M  SBL   2  2   2"),
            input.replace("M  SBL   2  1   2", "M  SBL   2  1  99"),
            input.replace("M  SAL   2  1   3", "M  SAL   2  1  99"),
        ];
        for malformed_input in malformed {
            read_v2000_detached(&malformed_input)
                .expect_err("strict parsing must reject the malformed SGroup");
            let (topology, _, _) = read_v2000_detached_with_params(
                &malformed_input,
                MolBlockReadParams {
                    strict_parsing: false,
                },
            )
            .expect("non-strict parsing must retain the molecule");
            assert_eq!(topology.substance_groups.len(), 1);
        }
    }

    #[test]
    fn v2000_non_strict_sap_infers_single_crossing_bond_or_drops_group() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            "malformed-sap\n  test            2D\n\n  4  2  0  0  0  0  0  0  0  0999 V2000\n{carbon}\n{carbon}\n{carbon}\n{carbon}\n  1  2  1  0  0  0  0\n  3  4  1  0  0  0  0\nM  STY  2   1 SUP   2 SUP\nM  SAL   1  1   1\nM  SAP   1  1   1\nM  SAL   2  1   3\nM  SBL   2  1   2\nM  SAP   2  1   3\nM  END\n"
        );
        read_v2000_detached(&input).expect_err("strict parsing rejects a short SAP line");

        let (topology, _, _) = read_v2000_detached_with_params(
            &input,
            MolBlockReadParams {
                strict_parsing: false,
            },
        )
        .expect("non-strict SAP cleanup");
        assert_eq!(topology.substance_groups.len(), 1);
        let retained = &topology.substance_groups[0];
        assert_eq!(retained.rdkit_sequence_id(), Some(2));
        assert_eq!(retained.attach_points().len(), 1);
        assert_eq!(retained.attach_points()[0].atom, AtomId::new(2));
        assert_eq!(
            retained.attach_points()[0].leaving_atom,
            Some(AtomId::new(3))
        );
        assert_eq!(retained.attach_points()[0].label.as_deref(), Some("  "));
    }

    #[test]
    fn v2000_non_strict_sgroup_field_failures_drop_only_the_affected_group() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let prefix = format!(
            "invalid-fields\n  test            2D\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n{carbon}\n{carbon}\nM  STY  2   1 SUP   2 SUP\nM  SAL   1  1   1\nM  SAL   2  1   2\n"
        );
        let invalid_lines = [
            "M  SST  1   1 BAD",
            "M  SLB  1   1 abc",
            "M  SCN  1   1 XX",
            "M  SMT   1",
            "M  SDI   1  3",
            "M  SNC  1   1 257",
            "M  SCL   1",
            "M  SBT  1   1   2",
        ];
        for invalid_line in invalid_lines {
            let input = format!("{prefix}{invalid_line}\nM  END\n");
            read_v2000_detached(&input).expect_err("strict SGroup field validation");
            let (topology, _, _) = read_v2000_detached_with_params(
                &input,
                MolBlockReadParams {
                    strict_parsing: false,
                },
            )
            .expect("non-strict SGroup field recovery");
            assert_eq!(topology.substance_groups.len(), 1, "{invalid_line}");
            assert_eq!(
                topology.substance_groups[0].rdkit_sequence_id(),
                Some(2),
                "{invalid_line}"
            );
        }
    }

    #[test]
    fn v2000_missing_sgroup_references_are_ignored_and_invalid_types_are_lax_only() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let missing_reference = format!(
            "missing-reference\n  test            2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{carbon}\nM  STY  1   2 SUP\nM  SAL   9  1   1\nM  END\n"
        );
        let (topology, _, _) =
            read_v2000_detached(&missing_reference).expect("FindSgIdx is warning-only");
        assert_eq!(topology.substance_groups.len(), 1);
        assert_eq!(topology.substance_groups[0].rdkit_sequence_id(), Some(2));

        let invalid_type = format!(
            "invalid-type\n  test            2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{carbon}\nM  STY  2   1 BAD   2 SUP\nM  END\n"
        );
        read_v2000_detached(&invalid_type).expect_err("strict invalid SGroup type");
        let (topology, _, _) = read_v2000_detached_with_params(
            &invalid_type,
            MolBlockReadParams {
                strict_parsing: false,
            },
        )
        .expect("non-strict invalid SGroup type is ignored");
        assert_eq!(topology.substance_groups.len(), 1);
        assert_eq!(topology.substance_groups[0].rdkit_sequence_id(), Some(2));
    }

    #[test]
    fn v2000_writer_wraps_utf8_sgroup_data_at_character_boundaries() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let ascii_prefix = "a".repeat(68);
        let input = format!(
            "utf8-group\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{carbon}\nM  STY  1   1 DAT\nM  SDT   1 {:<30}{:<2}{:<20}{:<2}\nM  SCD   1 {ascii_prefix}\nM  SED   1 界\nM  END\n",
            "FIELD", "T", "", ""
        );
        let (topology, coordinates, properties) =
            read_v2000_detached(&input).expect("read UTF-8 SGroup data");
        let expected = format!("{ascii_prefix}界");
        assert_eq!(
            topology.substance_groups[0].data().unwrap().values,
            [expected.as_str()]
        );

        let output = write_v2000_detached(&topology, &coordinates, &properties)
            .expect("write UTF-8 SGroup data");
        assert!(output.contains(&format!("M  SCD   1 {ascii_prefix}\n")));
        assert!(output.contains("M  SED   1 界\n"));
        let (roundtrip, _, _) = read_v2000_detached(&output).expect("roundtrip UTF-8 SGroup data");
        assert_eq!(
            roundtrip.substance_groups[0].data().unwrap().values,
            [expected]
        );
    }

    #[test]
    fn v2000_reader_handles_mass_difference_and_hydrogen_shorthand() {
        let carbon = v2000_atom_line("C", 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let deuterium = v2000_atom_line("D", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let tritium = v2000_atom_line("T", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            "isotopes\n  test\n\n  3  0  0  0  0  0  0  0  0  0999 V2000\n{carbon}\n{deuterium}\n{tritium}\nM  END\n"
        );
        let (topology, _, _) = read_v2000_detached(&input).expect("read isotopes");
        assert_eq!(topology.atoms[0].isotope(), Some(13));
        assert_eq!(topology.atoms[1].isotope(), Some(2));
        assert_eq!(topology.atoms[2].isotope(), Some(3));
    }

    #[test]
    fn v2000_reader_fails_closed_for_query_only_records() {
        let query_atom = v2000_atom_line("*", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            "query\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{query_atom}\nM  END\n"
        );
        assert!(matches!(
            read_v2000_detached(&input),
            Err(super::SdfReadError::Unsupported(_))
        ));
        let MolBlockRecord::Query(query_record) =
            read_mol_block_detached(&input).expect("read wildcard query graph")
        else {
            panic!("wildcard atom must produce a query record");
        };
        assert_eq!(query_record.query.num_atoms(), 1);
        assert_eq!(
            query_record.query.atoms()[0].predicate(),
            &QueryNode::predicate(AtomQueryPredicate::Any)
        );

        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            "query bond\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n{carbon}\n{carbon}\n  1  2  5  0  0  0  0\nM  END\n"
        );
        assert!(matches!(
            read_v2000_detached(&input),
            Err(super::SdfReadError::Unsupported(_))
        ));
        let MolBlockRecord::Query(query_record) =
            read_mol_block_detached(&input).expect("read query bond graph")
        else {
            panic!("query bond must produce a query record");
        };
        assert_eq!(
            query_record.query.bonds()[0].predicate(),
            &QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                cosmolkit_types::BondOrder::Single,
                cosmolkit_types::BondOrder::Double,
            ]))
        );

        let input = format!(
            "unknown query bond\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n{carbon}\n{carbon}\n  1  2 42  0  0  0  0\nM  END\n"
        );
        let MolBlockRecord::Query(query_record) =
            read_mol_block_detached(&input).expect("read unknown query bond")
        else {
            panic!("unknown bond type must use an any query");
        };
        assert_eq!(
            query_record.query.bonds()[0].predicate(),
            &QueryNode::predicate(BondQueryPredicate::Any)
        );
    }

    #[test]
    fn v2000_reader_preserves_atom_list_and_query_property_records() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let r_group = v2000_atom_line("R", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let oxygen = v2000_atom_line("O", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            concat!(
                "query properties\n  test\n\n",
                "  3  2  0  0  0  0  0  0  0  0999 V2000\n",
                "{}\n{}\n{}\n",
                "  1  2  1  0  0  0  0\n",
                "  2  3  1  0  0  0  0\n",
                "M  ALS   1  2 F C   N   \n",
                "M  RGP  1   2   7\n",
                "M  RBC  1   1  -2\n",
                "M  SUB  1   2  -2\n",
                "M  UNS  1   1   1\n",
                "M  END\n",
            ),
            carbon, r_group, oxygen
        );
        assert!(matches!(
            read_v2000_detached(&input),
            Err(super::SdfReadError::Unsupported(_))
        ));
        let MolBlockRecord::Query(record) =
            read_mol_block_detached(&input).expect("read V2000 query properties")
        else {
            panic!("query properties must produce a query record");
        };
        assert_eq!(record.query.num_atoms(), 3);
        assert_eq!(record.properties.prop("_NeedsQueryScan"), Some("1"));
        assert_eq!(record.query.prop("_NeedsQueryScan"), Some("1"));
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(vec![6, 7])),
                    QueryNode::predicate(AtomQueryPredicate::RingBondCount(0xDEAD_BEEF)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
            ])
        );
        assert_eq!(record.query.atoms()[1].prop("_MolFileRLabel"), Some("7"));
        assert_eq!(record.query.atoms()[1].prop("dummyLabel"), Some("R7"));
        assert_eq!(record.query.atoms()[1].atom().isotope(), Some(7));
        assert_eq!(
            record.query.atoms()[1].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::Any),
                QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(2)),
            ])
        );
    }

    #[test]
    fn v2000_reader_lowers_old_atom_lists_and_new_lists_override_them() {
        let not_list = include_str!(
            "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/not-list-query.mol"
        );
        let old_only = not_list.replace("M  ALS   4  2 T N   O   \n", "");
        let MolBlockRecord::Query(record) =
            read_mol_block_detached(&old_only).expect("read pinned old not-list record")
        else {
            panic!("old atom list must remain query state");
        };
        assert_eq!(
            record.query.atoms()[3].predicate(),
            &QueryNode::predicate(AtomQueryPredicate::AtomicNumberNotIn(vec![7, 8]))
        );
        assert_eq!(record.query.atoms()[3].atom().element().atomic_number(), 7);

        let conflicting = include_str!(
            "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/conflicting-list-query.mol"
        );
        let MolBlockRecord::Query(record) = read_mol_block_detached(conflicting)
            .expect("read pinned old/new conflicting-list fixture")
        else {
            panic!("conflicting atom lists must remain query state");
        };
        assert_eq!(
            record.query.atoms()[3].predicate(),
            &QueryNode::predicate(AtomQueryPredicate::AtomicNumberNotIn(vec![7, 8]))
        );
    }

    #[test]
    fn v2000_reader_rejects_invalid_unsaturation_query_value() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            "bad query\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{carbon}\nM  UNS  1   1   2\nM  END\n"
        );
        assert!(matches!(
            read_mol_block_detached(&input),
            Err(super::SdfReadError::Parse(message)) if message.contains("unsaturation query")
        ));
    }

    #[test]
    fn sdf_record_preserves_data_fields() {
        let input = "ethanol\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n>  <ID>\nabc\n\n$$$$\n";
        let record = read_sdf_record_detached(input).expect("read SDF record");
        assert_eq!(
            record.data_fields,
            vec![("ID".to_owned(), "abc".to_owned())]
        );
        assert_eq!(record.properties.prop("ID"), Some("abc"));
        let output =
            write_sdf_record_detached(&record.topology, &record.coordinates, &record.properties)
                .expect("write SDF record");
        assert!(output.contains(">  <ID>  \nabc\n"));
        assert!(output.ends_with("$$$$\n"));
    }

    #[test]
    fn sdf_reader_applies_typed_atom_and_bond_property_lists() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let oxygen = v2000_atom_line("O", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            concat!(
                "property lists\n  test\n\n",
                "  2  1  0  0  0  0  0  0  0  0999 V2000\n",
                "{}\n{}\n",
                "  1  2  1  0  0  0  0\n",
                "M  END\n",
                ">  <atom.prop.Label>\nC1 O1\n\n",
                ">  <atom.iprop.Score>\n[?] 7 ?\n\n",
                ">  <atom.dprop.Partial>\n0.25 invalid\n\n",
                ">  <atom.bprop.Active>\n1 0\n\n",
                ">  <bond.prop.Label>\nsingle\n\n",
                "$$$$\n",
            ),
            carbon, oxygen
        );
        let record = read_sdf_record_detached(&input).expect("read property lists");
        assert_eq!(record.topology.atoms[0].prop("Label"), Some("C1"));
        assert_eq!(record.topology.atoms[1].prop("Label"), Some("O1"));
        assert_eq!(record.topology.atoms[0].prop("Score"), Some("7"));
        assert_eq!(record.topology.atoms[1].prop("Score"), None);
        assert_eq!(record.topology.atoms[0].prop("Partial"), Some("0.25"));
        assert_eq!(record.topology.atoms[1].prop("Partial"), None);
        assert_eq!(record.topology.atoms[0].prop("Active"), Some("true"));
        assert_eq!(record.topology.atoms[1].prop("Active"), Some("false"));
        assert_eq!(record.topology.bonds[0].prop("Label"), Some("single"));
        assert_eq!(record.properties.sdf_property_lists().len(), 5);
        assert_eq!(
            record.properties.sdf_property_lists()[0].target(),
            SdfPropertyListTarget::Atom
        );
        assert_eq!(
            record.properties.sdf_property_lists()[1].values(),
            &[Some("7".to_owned()), None]
        );
    }

    #[test]
    fn sdf_data_params_can_preserve_property_lists_without_applying_them() {
        let carbon = v2000_atom_line("C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let input = format!(
            concat!(
                "property list policy\n  test\n\n",
                "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
                "{}\n",
                "M  END\n",
                ">  <atom.iprop.Score>\n7\n\n",
                "$$$$\n",
            ),
            carbon
        );
        let record = read_sdf_record_detached_with_params(
            &input,
            SdfDataReadParams {
                process_property_lists: false,
                ..SdfDataReadParams::default()
            },
        )
        .expect("read property list without applying it");
        assert_eq!(record.properties.prop("atom.iprop.Score"), Some("7"));
        assert_eq!(record.topology.atoms[0].prop("Score"), None);
        assert!(record.properties.sdf_property_lists().is_empty());
    }

    #[test]
    fn query_sdf_reader_applies_property_lists_without_concrete_lowering() {
        let input = concat!(
            "query property lists\n  test\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 * 0 0 0 0\n",
            "M  V30 2 O 1 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 5 1 2\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
            ">  <atom.prop.Label>\nany oxygen\n\n",
            ">  <bond.iprop.Score>\n4\n\n",
            "$$$$\n",
        );
        let record = read_sdf_graph_record_detached(input).expect("read query property lists");
        let MolBlockRecord::Query(query_record) = record.mol_block else {
            panic!("query property list record must remain a query graph");
        };
        assert_eq!(query_record.query.atoms()[0].prop("Label"), Some("any"));
        assert_eq!(query_record.query.atoms()[1].prop("Label"), Some("oxygen"));
        assert_eq!(
            query_record.query.bonds()[0].bond().prop("Score"),
            Some("4")
        );
        assert_eq!(query_record.properties.sdf_property_lists().len(), 2);
    }

    #[test]
    fn sdf_writer_uses_rdkit_field_header_and_skips_illegal_blank_lines() {
        let input = "fields\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n";
        let (topology, coordinates, properties) = read_v2000_detached(input).expect("read");
        let properties = properties
            .with_sdf_data_field("ID", "abc")
            .with_sdf_data_field("bad\nname", "ignored")
            .with_sdf_data_field("BAD_VALUE", "one\n\ntwo");
        let output =
            write_sdf_record_detached(&topology, &coordinates, &properties).expect("write fields");
        assert!(output.contains(">  <ID>  \nabc\n\n"));
        assert!(!output.contains("bad\nname"));
        assert!(!output.contains("BAD_VALUE"));
    }

    #[test]
    fn v3000_detached_reader_and_writer_round_trip_core_graph() {
        let input = "ethanol\n  test\n\n  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.0 0.0 0.0 0\nM  V30 2 O 1.2 0.0 0.0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        let (topology, coordinates, properties) = read_v3000_detached(input).expect("read V3000");
        let output =
            write_v3000_detached(&topology, &coordinates, &properties).expect("write V3000");
        let (roundtrip, _, roundtrip_properties) =
            read_v3000_detached(&output).expect("roundtrip V3000");
        assert_eq!(roundtrip.atoms.len(), 2);
        assert_eq!(roundtrip.bonds.len(), 1);
        assert_eq!(roundtrip_properties.name(), Some("ethanol"));
    }

    #[test]
    fn v3000_reader_preserves_rgroups_as_query_graph_state() {
        let input = concat!(
            "rgroup\n  test\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 R# 0 0 0 0 RGROUPS=(2 7 12) FOO=BAR\n",
            "M  V30 END ATOM\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let MolBlockRecord::Query(record) =
            read_mol_block_detached(input).expect("read V3000 R-group")
        else {
            panic!("R-group must remain query graph state");
        };
        let atom = &record.query.atoms()[0];
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::Any)
        );
        assert_eq!(atom.atom().prop("_MolFileRLabel"), Some("12"));
        assert_eq!(atom.atom().prop("dummyLabel"), Some("R12"));
        assert_eq!(atom.atom().isotope(), Some(12));
    }

    #[test]
    fn v3000_reader_rejects_malformed_rgroup_lists() {
        let block = |rgroups: &str| {
            format!(
                concat!(
                    "bad rgroup\n  test\n\n",
                    "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
                    "M  V30 BEGIN CTAB\n",
                    "M  V30 COUNTS 1 0 0 0 0\n",
                    "M  V30 BEGIN ATOM\n",
                    "M  V30 1 R# 0 0 0 0 RGROUPS={}\n",
                    "M  V30 END ATOM\n",
                    "M  V30 END CTAB\n",
                    "M  END\n",
                ),
                rgroups
            )
        };

        for (value, expected) in [
            ("1", "Missing parens"),
            ("(2 7)", "Not enough values"),
            ("(1 x)", "Cannot convert 'x'"),
        ] {
            let error = read_mol_block_detached(&block(value)).unwrap_err();
            assert!(error.to_string().contains(expected), "{error}");
        }
    }

    #[test]
    fn v3000_reader_preserves_continued_atom_bond_and_bookmark_state() {
        let input = concat!(
            "state\n",
            "     RDKit          2D\n",
            "comment\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 1\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0.0 0.0 0.0 7 CHG=-1 MASS=13 RAD=-\n",
            "M  V30 2 CFG=1 VAL=4 STBOX=1 EXACHG=2 INVRET=1 ATTCHPT=3 -\n",
            "M  V30 ATTCHORD=2 CLASS=AA SEQID=9 SEQNAME=GLY\n",
            "M  V30 20 N 1.2 0.0 -0.0 0 STBOX=1\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 99 1 20 10 CFG=3 RXCTR=4 STBOX=1 -\n",
            "M  V30 ENDPTS=(2 20 10) ATTACH=ANY\n",
            "M  V30 END BOND\n",
            "M  V30 LINKNODE 1 2 2 10 20\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let (topology, coordinates, properties) =
            read_v3000_detached(input).expect("read continued V3000 state");
        let carbon = &topology.atoms[0];
        assert_eq!(carbon.formal_charge(), -1);
        assert_eq!(carbon.isotope(), Some(13));
        assert_eq!(carbon.radical_electrons(), 1);
        assert_eq!(carbon.atom_map(), Some(7));
        assert_eq!(carbon.mol_parity(), Some(1));
        assert_eq!(carbon.mol_inversion_flag(), Some(1));
        assert_eq!(carbon.prop("molTotValence"), Some("4"));
        assert_eq!(carbon.prop("molStereoCare"), Some("1"));
        assert_eq!(carbon.prop("molAttachPoint"), Some("3"));
        assert_eq!(carbon.prop("molAttachOrder"), Some("2"));
        assert_eq!(carbon.prop("molAtomClass"), Some("AA"));
        assert_eq!(carbon.prop("molAtomSeqId"), Some("9"));
        assert_eq!(carbon.prop("molAtomSeqName"), Some("GLY"));
        let bond = &topology.bonds[0];
        assert_eq!(bond.begin().index(), 1);
        assert_eq!(bond.end().index(), 0);
        assert_eq!(bond.direction(), BondDirection::BeginDash);
        assert_eq!(bond.prop("_MolFileBondCfg"), Some("3"));
        assert_eq!(bond.prop("molReactStatus"), Some("4"));
        assert_eq!(bond.prop("_MolFileBondEndPts"), Some("(2 20 10)"));
        assert_eq!(bond.prop("_MolFileBondAttach"), Some("ANY"));
        assert_eq!(properties.prop("_MolFileComments"), Some("comment"));
        assert_eq!(properties.prop("_MolFileChiralFlag"), Some("1"));
        assert_eq!(properties.prop("_MolFileLinkNodes"), Some("1 2 2 10 20"));
        assert_eq!(coordinates.conformers_2d.len(), 1);
        assert!(coordinates.conformers_3d.is_empty());

        let output = write_v3000_detached(&topology, &coordinates, &properties)
            .expect("write preserved V3000 state");
        assert!(output.contains(" CHG=-1"));
        assert!(output.contains(" MASS=13"));
        assert!(output.contains(" RAD=2"));
        assert!(output.contains(" CFG=3 RXCTR=4 STBOX=1 ENDPTS=(2 20 10) ATTACH=ANY"));
        let (roundtrip, _, roundtrip_properties) =
            read_v3000_detached(&output).expect("roundtrip preserved V3000 state");
        assert_eq!(roundtrip.atoms[0].formal_charge(), -1);
        assert_eq!(roundtrip.atoms[0].radical_electrons(), 1);
        assert_eq!(roundtrip.atoms[0].mol_parity(), Some(1));
        assert_eq!(roundtrip.bonds[0].direction(), BondDirection::BeginDash);
        assert_eq!(
            roundtrip_properties.prop("_MolFileLinkNodes"),
            Some("1 2 2 10 20")
        );
    }

    #[test]
    fn v3000_reader_uses_rdkit_coordinate_dimension_rules() {
        let input = concat!(
            "zero-z\n",
            "     RDKit          3D\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 5 C 0.0 0.0 -0.0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let (_, coordinates, _) = read_v3000_detached(input).expect("read marked 3D");
        assert!(coordinates.conformers_2d.is_empty());
        assert_eq!(coordinates.conformers_3d.len(), 1);
        assert!(coordinates.conformers_3d[0].coordinates()[0][2].is_sign_negative());
    }

    #[test]
    fn v3000_reader_preserves_query_atoms_bonds_and_properties() {
        let input = concat!(
            "query\n  test\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 [C,N] 0 0 0 0 HCOUNT=2 RBCNT=2 UNSAT=1\n",
            "M  V30 2 O 1 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 5 1 2 TOPO=1\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        assert!(matches!(
            read_v3000_detached(input),
            Err(super::SdfReadError::Unsupported(_))
        ));

        let MolBlockRecord::Query(record) =
            read_mol_block_detached(input).expect("read V3000 query graph")
        else {
            panic!("query state must produce a query record");
        };
        assert_eq!(record.query.num_atoms(), 2);
        assert_eq!(record.query.num_bonds(), 1);
        assert_eq!(record.query.name(), Some("query"));
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::and(vec![
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(vec![6, 7])),
                        QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(2),),
                    ]),
                    QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
            ])
        );
        assert_eq!(
            record.query.bonds()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                    cosmolkit_types::BondOrder::Single,
                    cosmolkit_types::BondOrder::Double,
                ])),
                QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            ])
        );
        assert_eq!(
            record.source_coordinate_dim,
            Some(cosmolkit_model::CoordinateDimension::TwoD)
        );
    }

    #[test]
    fn query_sdf_record_preserves_data_fields() {
        let input = concat!(
            "query fields\n  test\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 * 0 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 END CTAB\n",
            "M  END\n",
            ">  <ID>\n",
            "query-17\n\n",
            "$$$$\n",
        );
        let record = read_sdf_graph_record_detached(input).expect("read query SDF record");
        assert_eq!(
            record.data_fields,
            vec![("ID".to_owned(), "query-17".to_owned())]
        );
        let MolBlockRecord::Query(query_record) = record.mol_block else {
            panic!("wildcard SDF record must remain a query graph");
        };
        assert_eq!(query_record.query.prop("ID"), Some("query-17"));
        assert_eq!(
            query_record.properties.sdf_data_fields(),
            &[("ID".to_owned(), "query-17".to_owned())]
        );
    }

    #[test]
    fn v3000_reader_preserves_typed_sgroups_and_stereo_collections() {
        let input = concat!(
            "v3000-sgroup\n",
            "  COSMolKit\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 2 0 1\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0.0 0.0 0.0 0\n",
            "M  V30 20 O 1.25 0.0 0.5 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 99 1 10 20\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 2 DAT 0 FIELDNAME=FIELD FIELDTYPE=T FIELDINFO=INFO -\n",
            "M  V30 QUERYTYPE=Q QUERYOP=OP FIELDDISP=\"display spec\" -\n",
            "M  V30 FIELDDATA=\"payload\" PARENT=1 COMPNO=5\n",
            "M  V30 1 SUP 7 ATOMS=(1 10) XBONDS=(1 99) LABEL=Me -\n",
            "M  V30 CONNECT=HT BRKXYZ=(9 0 1 0 2 3 0 0 0 0) -\n",
            "M  V30 CSTATE=(4 99 0.5 0.25 0.0) SAP=(3 10 20 AP)\n",
            "M  V30 END SGROUP\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/STEABS ATOMS=(1 10)\n",
            "M  V30 MDLV30/STEREL2 ATOMS=(1 20)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let (topology, _, _) = read_v3000_detached(input).expect("read typed V3000 state");
        assert_eq!(topology.substance_groups.len(), 2);
        let sup = &topology.substance_groups[0];
        assert_eq!(sup.kind(), &SubstanceGroupKind::Superatom);
        assert_eq!(sup.rdkit_sequence_id(), Some(1));
        assert_eq!(sup.external_id(), Some(7));
        assert_eq!(sup.atoms(), &[AtomId::new(0)]);
        assert_eq!(sup.bonds(), &[BondId::new(0)]);
        assert_eq!(sup.bond_role(BondId::new(0)), SGroupBondRole::Crossing);
        assert_eq!(sup.label(), Some("Me"));
        assert_eq!(sup.connection(), Some(&SGroupConnection::HeadToTail));
        assert_eq!(sup.display().unwrap().brackets[0].p1, [0.0, 1.0]);
        assert_eq!(sup.display().unwrap().brackets[0].p2, [2.0, 3.0]);
        assert_eq!(sup.cstates()[0].bond, BondId::new(0));
        assert_eq!(sup.cstates()[0].vector, [0.5, 0.25]);
        assert_eq!(sup.attach_points()[0].atom, AtomId::new(0));
        assert_eq!(sup.attach_points()[0].leaving_atom, Some(AtomId::new(1)));
        assert_eq!(sup.attach_points()[0].label.as_deref(), Some("AP"));

        let dat = &topology.substance_groups[1];
        assert_eq!(dat.kind(), &SubstanceGroupKind::Data);
        assert_eq!(dat.parent(), Some(sup.id()));
        assert_eq!(dat.component_number(), Some(5));
        let data = dat.data().unwrap();
        assert_eq!(data.field_name.as_deref(), Some("FIELD"));
        assert_eq!(data.field_type.as_deref(), Some("T"));
        assert_eq!(data.field_info.as_deref(), Some("INFO"));
        assert_eq!(data.query_type.as_deref(), Some("Q"));
        assert_eq!(data.query_op.as_deref(), Some("OP"));
        assert_eq!(data.field_display.as_deref(), Some("display spec"));
        assert_eq!(data.values, ["payload"]);

        assert_eq!(topology.stereo_groups.len(), 2);
        assert_eq!(topology.stereo_groups[0].kind(), StereoGroupKind::Absolute);
        assert_eq!(topology.stereo_groups[0].atoms(), &[AtomId::new(0)]);
        assert_eq!(topology.stereo_groups[1].kind(), StereoGroupKind::Or);
        assert_eq!(topology.stereo_groups[1].id(), Some(2));
        assert_eq!(topology.stereo_groups[1].atoms(), &[AtomId::new(1)]);

        let output = write_v3000_detached(
            &topology,
            &CoordinateBlock::default(),
            &MoleculeProperties::default(),
        )
        .expect("write typed V3000 state");
        assert!(output.contains("M  V30 COUNTS 2 1 2 0 0"));
        assert!(output.contains("M  V30 BEGIN SGROUP"));
        assert!(output.contains("M  V30 BEGIN COLLECTION"));
        let (roundtrip, _, _) = read_v3000_detached(&output).expect("roundtrip typed V3000 state");
        assert_eq!(roundtrip.substance_groups.len(), 2);
        assert_eq!(roundtrip.substance_groups[0].label(), Some("Me"));
        assert_eq!(
            roundtrip.substance_groups[1].parent(),
            Some(roundtrip.substance_groups[0].id())
        );
        assert_eq!(
            roundtrip.substance_groups[1].data().unwrap().values,
            ["payload"]
        );
        assert_eq!(roundtrip.stereo_groups, topology.stereo_groups);
    }

    #[test]
    fn v3000_reader_requires_declared_sgroup_block() {
        let input = concat!(
            "missing\n\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 0 0 1 0 0\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        assert!(
            read_v3000_detached(input)
                .expect_err("missing declared SGroup block")
                .to_string()
                .contains("BEGIN SGROUP line not found")
        );
        let (topology, _, _) = read_v3000_detached_with_params(
            input,
            MolBlockReadParams {
                strict_parsing: false,
            },
        )
        .expect("non-strict parsing accepts a missing declared SGroup block");
        assert!(topology.substance_groups.is_empty());
    }

    #[test]
    fn v3000_malformed_sgroups_throw_strictly_and_are_dropped_non_strictly() {
        let input = concat!(
            "malformed-v3000-groups\n",
            "  COSMolKit          2D\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 4 2 2 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0 0 0 0\n",
            "M  V30 2 C 1 0 0 0\n",
            "M  V30 3 C 2 0 0 0\n",
            "M  V30 4 C 3 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 1 2\n",
            "M  V30 2 1 3 4\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 1 SUP 1 ATOMS=(1 1) XBONDS=(1 1) LABEL=A CSTATE=(4 1 0.5 0.25 0)\n",
            "M  V30 2 SUP 2 ATOMS=(1 3) XBONDS=(1 2) LABEL=B CSTATE=(4 2 -0.5 -0.25 0)\n",
            "M  V30 END SGROUP\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let (base, _, _) = read_v3000_detached(input).expect("valid source-shaped SGroups");
        assert_eq!(base.substance_groups.len(), 2);

        let malformed = [
            input.replace("XBONDS=(1 1)", "XBONDS=(2 1)"),
            input.replace("ATOMS=(1 1)", "ATOMS=(2 1)"),
            input.replace("XBONDS=(1 2)", "XBONDS=(1 99)"),
            input.replace("ATOMS=(1 3)", "ATOMS=(1 99)"),
            input.replace("CSTATE=(4 1 0.5 0.25 0)", "CSTATE=(3 1 0.5 0.25)"),
        ];
        for malformed_input in malformed {
            read_v3000_detached(&malformed_input)
                .expect_err("strict parsing must reject the malformed V3000 SGroup");
            let (topology, _, _) = read_v3000_detached_with_params(
                &malformed_input,
                MolBlockReadParams {
                    strict_parsing: false,
                },
            )
            .expect("non-strict V3000 SGroup cleanup");
            assert_eq!(topology.substance_groups.len(), 1);
        }
    }

    #[test]
    fn v3000_sgroup_block_structure_follows_strict_policy() {
        let input = concat!(
            "v3000-group-structure\n",
            "  COSMolKit          2D\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 1 DAT 0 ATOMS=(1 1)\n",
            "M  V30 END SGROUP\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let lax = MolBlockReadParams {
            strict_parsing: false,
        };

        let missing_end = input.replace("M  V30 END SGROUP\n", "");
        read_v3000_detached(&missing_end).expect_err("strict missing END SGROUP");
        assert_eq!(
            read_v3000_detached_with_params(&missing_end, lax)
                .expect("non-strict missing END SGROUP")
                .0
                .substance_groups
                .len(),
            1
        );

        let missing_group = input.replace("COUNTS 1 0 1 0 0", "COUNTS 1 0 2 0 0");
        read_v3000_detached(&missing_group).expect_err("strict missing declared SGroup");
        assert_eq!(
            read_v3000_detached_with_params(&missing_group, lax)
                .expect("non-strict missing declared SGroup")
                .0
                .substance_groups
                .len(),
            1
        );

        let unexpected_group = input.replace("COUNTS 1 0 1 0 0", "COUNTS 1 0 0 0 0");
        read_v3000_detached(&unexpected_group).expect_err("strict unexpected SGroup block");
        assert_eq!(
            read_v3000_detached_with_params(&unexpected_group, lax)
                .expect("non-strict unexpected SGroup block")
                .0
                .substance_groups
                .len(),
            1
        );

        let invalid_type = input.replace("1 DAT 0", "1 BAD 0");
        read_v3000_detached(&invalid_type).expect_err("strict invalid V3000 SGroup type");
        let (topology, _, _) = read_v3000_detached_with_params(&invalid_type, lax)
            .expect("non-strict generic V3000 SGroup type");
        assert_eq!(
            topology.substance_groups[0].kind(),
            &SubstanceGroupKind::Generic("BAD".to_owned())
        );

        let duplicate = input
            .replace("COUNTS 1 0 1 0 0", "COUNTS 1 0 2 0 0")
            .replace(
                "M  V30 1 DAT 0 ATOMS=(1 1)\n",
                "M  V30 1 DAT 0 ATOMS=(1 1)\nM  V30 1 DAT 0 ATOMS=(1 1)\n",
            );
        read_v3000_detached(&duplicate).expect_err("strict duplicate sequence count mismatch");
        assert_eq!(
            read_v3000_detached_with_params(&duplicate, lax)
                .expect("non-strict duplicate sequence keeps first")
                .0
                .substance_groups
                .len(),
            1
        );
    }

    #[test]
    fn sdf_stream_reader_preserves_record_order() {
        let one = "one\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n";
        let two = "two\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n";
        let stream = format!("{one}$$$$\n{two}$$$$\n");
        let records = read_sdf_records_detached(&stream).expect("read SDF stream");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].properties.name(), Some("one"));
        assert_eq!(records[1].properties.name(), Some("two"));
    }

    #[test]
    fn sdf_data_fields_follow_rdkit_header_and_multiline_rules() {
        let input = "fields\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n> <ID> trailing text\r\nabc\r\n\r\n  >  <NOTE>\nalpha\nbeta\n\n> <SPACES>\n  \n\t\n\n$$$$\n";
        let record = read_sdf_record_detached(input).expect("read fields");
        assert_eq!(
            record.data_fields,
            vec![
                ("ID".to_owned(), "abc".to_owned()),
                ("NOTE".to_owned(), "alpha\nbeta".to_owned()),
                ("SPACES".to_owned(), "  \n\t".to_owned()),
            ]
        );
        assert_eq!(record.properties.sdf_data_fields(), record.data_fields);
    }

    #[test]
    fn sdf_data_fields_ignore_invalid_headers_until_blank_like_rdkit() {
        let input = "fields\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n> <>\nignored\n\n> <ID>\nkept\n\n$$$$\n";
        let record = read_sdf_record_detached(input).expect("read fields");
        assert_eq!(
            record.data_fields,
            vec![("ID".to_owned(), "kept".to_owned())]
        );

        let unterminated = "fields\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n> <>\nunterminated\n$$$$\n";
        assert!(
            read_sdf_record_detached(unterminated)
                .expect_err("unterminated invalid header")
                .to_string()
                .contains("End of data field name not found")
        );
    }

    #[test]
    fn sdf_data_fields_reject_spurious_non_header_data_in_strict_mode() {
        let input =
            "fields\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\nspurious\n$$$$\n";
        assert!(
            read_sdf_record_detached(input)
                .expect_err("spurious data")
                .to_string()
                .contains("Problems encountered parsing data fields")
        );
    }

    #[test]
    fn sdf_data_fields_ignore_spurious_content_in_non_strict_mode() {
        let input = "fields\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\nspurious\n> <ID>\nkept\n\n$$$$\n";
        let record = read_sdf_record_detached_with_params(
            input,
            SdfDataReadParams {
                strict_parsing: false,
                ..SdfDataReadParams::default()
            },
        )
        .expect("non-strict data parsing");
        assert_eq!(
            record.data_fields,
            vec![("ID".to_owned(), "kept".to_owned())]
        );
    }

    #[test]
    fn sdf_stream_framing_only_splits_delimiters_at_line_start() {
        let first = "one\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n> <TEXT>\ninside $$$$ value\n\n$$$$ suffix is accepted\n";
        let second = "two\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n$$$$\n";
        let records =
            read_sdf_records_detached(&format!("{first}{second}")).expect("read framed stream");
        assert_eq!(records.len(), 2);
        assert_eq!(
            records[0].data_fields,
            vec![("TEXT".to_owned(), "inside $$$$ value".to_owned())]
        );
        assert_eq!(records[1].properties.name(), Some("two"));
    }

    #[test]
    fn sdf_stream_framing_preserves_blank_title_header_lines() {
        let first = "\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n$$$$\n";
        let second = "named\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n";
        let records = read_sdf_records_detached(&format!("{first}{second}"))
            .expect("read blank-title stream");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].properties.name(), None);
        assert_eq!(records[1].properties.name(), Some("named"));
    }

    #[test]
    fn detached_forward_reader_preserves_query_records_and_recovers_after_errors() {
        let malformed = "bad\n  test\n\ninvalid\n$$$$\n";
        let query = concat!(
            "query\n  test\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 * 0 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 END CTAB\n",
            "M  END\n",
            "$$$$\n",
        );
        let stream = format!("{malformed}{query}");
        let mut reader = SdfGraphReader::new(std::io::Cursor::new(stream.as_bytes()));
        assert!(reader.next_record().is_err());
        let record = reader
            .next_record()
            .expect("second read")
            .expect("query record");
        assert!(matches!(record.mol_block, MolBlockRecord::Query(_)));
        assert_eq!(reader.records_consumed(), 2);
        assert_eq!(reader.bytes_consumed(), stream.len() as u64);
        assert!(reader.next_record().expect("EOF").is_none());
        assert!(reader.is_end());
    }

    #[test]
    fn detached_dataset_indexes_metadata_and_seeks_directly_to_records() {
        let one = "one\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n$$$$\n";
        let two =
            "two\n  test\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n> <ID>\n2\n\n$$$$\n";
        let file = tempfile::NamedTempFile::new().expect("temporary SDF");
        std::fs::write(file.path(), format!("{one}{two}")).expect("write SDF");

        let dataset = SdfGraphDataset::open(file.path()).expect("index SDF");
        assert_eq!(dataset.len(), 2);
        assert_eq!(dataset.metadata(0).unwrap().title.as_deref(), Some("one"));
        assert_eq!(dataset.metadata(1).unwrap().byte_offset, one.len() as u64);
        assert_eq!(
            dataset.metadata(1).unwrap().line_offset,
            one.lines().count()
        );
        let record = dataset.record(1).expect("read indexed record");
        assert_eq!(record.data_fields, [("ID".to_owned(), "2".to_owned())]);
        assert!(
            dataset
                .record(2)
                .unwrap_err()
                .to_string()
                .contains("Index error")
        );
    }
}
