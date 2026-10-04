//! Private structural read dispatch primitives.

use crate::bio_chemcomp::{
    ChemCompModelError, check_chemcomp_block_number, make_structure_from_chemcomp_block,
    make_structure_from_chemcomp_doc,
};
use crate::bio_mmcif::{BioMmcifReadError, populate_mmcif_bio_structure_document};
use crate::bio_pdb::{BioPdbReadError, BioPdbReadParams, read_pdb_bio_structure};
use crate::cif::{
    CifCheckLevel, CifDocument, CifReadError, MmjsonReadError, read_cif_document,
    read_mmjson_insitu,
};
use cosmolkit_bio::{BioCoordinateFormat, BioStructureData};
use std::error::Error;
use std::fmt;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BioReadParams {
    pub format: BioCoordinateFormat,
    pub source_name: String,
}

impl Default for BioReadParams {
    fn default() -> Self {
        Self {
            format: BioCoordinateFormat::Unknown,
            source_name: "<string>".to_owned(),
        }
    }
}

/// A source-preserving error from structural coordinate reading.
#[derive(Debug)]
pub enum BioReadError {
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    Utf8 {
        path: PathBuf,
        source: std::string::FromUtf8Error,
    },
    Pdb(BioPdbReadError),
    Cif(CifReadError),
    Mmcif(BioMmcifReadError),
    Mmjson(Box<dyn Error + Send + Sync>),
    ChemComp(Box<dyn Error + Send + Sync>),
    WrongFormat {
        source_name: String,
        format: BioCoordinateFormat,
    },
    UnknownFileFormat(PathBuf),
}

impl fmt::Display for BioReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::Utf8 { path, source } => write!(
                f,
                "{}: coordinate text is not UTF-8: {source}",
                path.display()
            ),
            Self::Pdb(source) => source.fmt(f),
            Self::Cif(source) => source.fmt(f),
            Self::Mmcif(source) => source.fmt(f),
            Self::Mmjson(source) => source.fmt(f),
            Self::ChemComp(source) => source.fmt(f),
            Self::WrongFormat {
                source_name,
                format,
            } => write!(
                f,
                "wrong format of coordinate file {source_name}: {format:?}"
            ),
            Self::UnknownFileFormat(path) => write!(f, "Unknown format of {}.", path.display()),
        }
    }
}

impl Error for BioReadError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::Utf8 { source, .. } => Some(source),
            Self::Pdb(source) => Some(source),
            Self::Cif(source) => Some(source),
            Self::Mmcif(source) => Some(source),
            Self::Mmjson(source) => Some(source.as_ref()),
            Self::ChemComp(source) => Some(source.as_ref()),
            Self::WrongFormat { .. } | Self::UnknownFileFormat(_) => None,
        }
    }
}

fn document_error(error: BioDocumentError) -> BioReadError {
    match error {
        BioDocumentError::ChemComp(error) => BioReadError::ChemComp(Box::new(error)),
        BioDocumentError::Mmcif(error) => BioReadError::Mmcif(error),
    }
}

fn memory_error(error: BioMemoryReadError) -> BioReadError {
    match error {
        BioMemoryReadError::Pdb(error) => BioReadError::Pdb(error),
        BioMemoryReadError::Cif(error) => BioReadError::Cif(error),
        BioMemoryReadError::Mmjson(error) => BioReadError::Mmjson(Box::new(error)),
        BioMemoryReadError::Document(error) => document_error(error),
        BioMemoryReadError::WrongFormat {
            source_name,
            format,
        } => BioReadError::WrongFormat {
            source_name,
            format,
        },
    }
}

fn file_error(error: BioFileReadError) -> BioReadError {
    match error {
        BioFileReadError::Load(BioFileLoadError::Io { path, source }) => {
            BioReadError::Io { path, source }
        }
        BioFileReadError::Load(BioFileLoadError::Utf8 { path, source }) => {
            BioReadError::Utf8 { path, source }
        }
        BioFileReadError::Memory(error) => memory_error(error),
        BioFileReadError::Pdb(error) => BioReadError::Pdb(error),
        BioFileReadError::Cif(error) => BioReadError::Cif(error),
        BioFileReadError::Mmjson(error) => BioReadError::Mmjson(Box::new(error)),
        BioFileReadError::Document(error) => document_error(error),
        BioFileReadError::ChemComp(error) => BioReadError::ChemComp(Box::new(error)),
        BioFileReadError::UnknownFormat(path) => BioReadError::UnknownFileFormat(path),
    }
}

pub fn read_bio_structure(
    text: &str,
    params: &BioReadParams,
) -> Result<BioStructureData, BioReadError> {
    read_structure_from_memory(text, &params.source_name, params.format).map_err(memory_error)
}

pub fn read_bio_structure_file(
    path: &Path,
    format: BioCoordinateFormat,
) -> Result<BioStructureData, BioReadError> {
    read_structure_file(path, format).map_err(file_error)
}

#[derive(Debug)]
enum BioFileLoadError {
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    Utf8 {
        path: PathBuf,
        source: std::string::FromUtf8Error,
    },
}

impl fmt::Display for BioFileLoadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::Utf8 { path, source } => {
                write!(
                    f,
                    "{}: coordinate text is not UTF-8: {source}",
                    path.display()
                )
            }
        }
    }
}

impl Error for BioFileLoadError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::Utf8 { source, .. } => Some(source),
        }
    }
}

fn read_regular_file(path: &Path) -> Result<String, BioFileLoadError> {
    // Gemmi❗❌: inline CharArray read_file_into_buffer(const std::string& path) {
    // Gemmi❗❌:   fileptr_t f = file_open(path.c_str(), "rb");
    // Gemmi❗❌:   size_t size = file_size(f.get(), path);
    // Gemmi❗❌:   CharArray buffer(size);
    // Gemmi❗❌:   if (std::fread(buffer.data(), size, 1, f.get()) != 1)
    // Gemmi❗❌:     sys_fail(path + ": fread failed");
    // Gemmi❗❌:   return buffer;
    // Gemmi❗❌: }
    // Behavior: regular-file bytes are loaded once with pathful IO errors.
    // The frozen Rust text boundary additionally rejects invalid UTF-8 with
    // its original bytes retained in FromUtf8Error; Gemmi accepts raw bytes.
    // Complexity: fs::read makes one owned buffer, but validating UTF-8 adds
    // one O(n) scan beyond Gemmi's size-seek and fread.
    let bytes = std::fs::read(path).map_err(|source| BioFileLoadError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    String::from_utf8(bytes).map_err(|source| BioFileLoadError::Utf8 {
        path: path.to_path_buf(),
        source,
    })
}

#[derive(Debug)]
enum BioDocumentError {
    ChemComp(ChemCompModelError),
    Mmcif(BioMmcifReadError),
}

#[derive(Debug)]
enum BioMemoryReadError {
    Pdb(BioPdbReadError),
    Cif(CifReadError),
    Mmjson(MmjsonReadError),
    Document(BioDocumentError),
    WrongFormat {
        source_name: String,
        format: BioCoordinateFormat,
    },
}

impl fmt::Display for BioMemoryReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Pdb(error) => error.fmt(f),
            Self::Cif(error) => error.fmt(f),
            Self::Mmjson(error) => error.fmt(f),
            Self::Document(error) => error.fmt(f),
            Self::WrongFormat {
                source_name,
                format,
            } => {
                write!(
                    f,
                    "wrong format of coordinate file {source_name}: {format:?}"
                )
            }
        }
    }
}

impl Error for BioMemoryReadError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Pdb(error) => Some(error),
            Self::Cif(error) => Some(error),
            Self::Mmjson(error) => Some(error),
            Self::Document(error) => Some(error),
            Self::WrongFormat { .. } => None,
        }
    }
}

impl fmt::Display for BioDocumentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ChemComp(error) => error.fmt(f),
            Self::Mmcif(error) => error.fmt(f),
        }
    }
}

impl Error for BioDocumentError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::ChemComp(error) => Some(error),
            Self::Mmcif(error) => Some(error),
        }
    }
}

fn make_structure_from_doc(
    document: &CifDocument,
    possible_chemcomp: bool,
) -> Result<BioStructureData, BioDocumentError> {
    // Gemmi✔️✔️: inline Structure make_structure_from_doc(cif::Document&& doc, bool possible_chemcomp,
    // Gemmi✔️✔️:                                          cif::Document* save_doc=nullptr) {
    // Gemmi✔️✔️:   if (possible_chemcomp) {
    // Gemmi✔️✔️:     // check for special case - refmac dictionary or CCD file
    // Gemmi✔️✔️:     int n = check_chemcomp_block_number(doc);
    // Gemmi✔️✔️:     if (n != -1)
    // Gemmi✔️✔️:       return make_structure_from_chemcomp_block(doc.blocks[n]);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return make_structure(std::move(doc), save_doc);
    // Gemmi✔️✔️: }
    // Behavior: the private document has no save_doc output; all observable
    // structure selection and typed failure paths retain the source order.
    // Complexity: one fixed-cost chemcomp selector, then exactly one owner
    // population; neither path clones or reparses the document.
    if possible_chemcomp {
        if let Some(index) = check_chemcomp_block_number(document) {
            return make_structure_from_chemcomp_block(&document.blocks()[index], 7)
                .map_err(BioDocumentError::ChemComp);
        }
    }
    populate_mmcif_bio_structure_document(document).map_err(BioDocumentError::Mmcif)
}

fn read_structure_from_memory(
    text: &str,
    source_name: &str,
    format: BioCoordinateFormat,
) -> Result<BioStructureData, BioMemoryReadError> {
    // Gemmi✔️✔️: inline Structure read_structure_from_memory(char* data, size_t size,
    // Gemmi✔️✔️:                                             const std::string& path,
    // Gemmi✔️✔️:                                             CoorFormat format=CoorFormat::Unknown,
    // Gemmi❌❌:                                             cif::Document* save_doc=nullptr) {
    // Gemmi❌❌:   if (save_doc)
    // Gemmi❌❌:     save_doc->clear();
    // Gemmi✔️✔️:   if (format == CoorFormat::Unknown || format == CoorFormat::Detect)
    // Gemmi✔️✔️:     format = coor_format_from_content(data, data + size);
    // Gemmi✔️✔️:   if (format == CoorFormat::Pdb)
    // Gemmi✔️✔️:     return read_pdb_from_memory(data, size, path);
    // Gemmi✔️✔️:   if (format == CoorFormat::Mmcif)
    // Gemmi✔️✔️:     return make_structure_from_doc(cif::read_memory(data, size, path.c_str()),
    // Gemmi❌❌:                                    true, save_doc);
    // Gemmi✔️✔️:   if (format == CoorFormat::Mmjson)
    // Gemmi❌❌:     return make_structure(cif::read_mmjson_insitu(data, size, path), save_doc);
    // Gemmi✔️✔️:   fail("wrong format of coordinate file " + path);
    // Gemmi✔️✔️: }
    // Behavior: the frozen Rust text entry has no save_doc carrier. The PDB,
    // CIF and mmJSON branches use their existing owners and retain typed errors;
    // the mmJSON branch deliberately takes ordinary document finalization.
    // Complexity: classification scans once; each selected owner parses once,
    // and the document is borrowed for one population without reserialization.
    let selected = if matches!(
        format,
        BioCoordinateFormat::Unknown | BioCoordinateFormat::Detect
    ) {
        coor_format_from_content(text.as_bytes())
    } else {
        format
    };
    match selected {
        BioCoordinateFormat::Pdb => {
            read_pdb_bio_structure(text, source_name, &BioPdbReadParams::default())
                .map_err(BioMemoryReadError::Pdb)
        }
        BioCoordinateFormat::Mmcif => {
            let document = read_cif_document(text, source_name, CifCheckLevel::Default)
                .map_err(BioMemoryReadError::Cif)?;
            make_structure_from_doc(&document, true).map_err(BioMemoryReadError::Document)
        }
        BioCoordinateFormat::Mmjson => {
            let document = read_mmjson_insitu(text.as_bytes(), source_name)
                .map_err(BioMemoryReadError::Mmjson)?;
            populate_mmcif_bio_structure_document(&document)
                .map_err(|error| BioMemoryReadError::Document(BioDocumentError::Mmcif(error)))
        }
        BioCoordinateFormat::Unknown
        | BioCoordinateFormat::Detect
        | BioCoordinateFormat::ChemComp => Err(BioMemoryReadError::WrongFormat {
            source_name: source_name.to_owned(),
            format: selected,
        }),
    }
}

#[derive(Debug)]
enum BioFileReadError {
    Load(BioFileLoadError),
    Memory(BioMemoryReadError),
    Pdb(BioPdbReadError),
    Cif(CifReadError),
    Mmjson(MmjsonReadError),
    Document(BioDocumentError),
    ChemComp(ChemCompModelError),
    UnknownFormat(PathBuf),
}

impl fmt::Display for BioFileReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Load(error) => error.fmt(f),
            Self::Memory(error) => error.fmt(f),
            Self::Pdb(error) => error.fmt(f),
            Self::Cif(error) => error.fmt(f),
            Self::Mmjson(error) => error.fmt(f),
            Self::Document(error) => error.fmt(f),
            Self::ChemComp(error) => error.fmt(f),
            Self::UnknownFormat(path) => write!(f, "Unknown format of {}.", path.display()),
        }
    }
}

impl Error for BioFileReadError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Load(error) => Some(error),
            Self::Memory(error) => Some(error),
            Self::Pdb(error) => Some(error),
            Self::Cif(error) => Some(error),
            Self::Mmjson(error) => Some(error),
            Self::Document(error) => Some(error),
            Self::ChemComp(error) => Some(error),
            Self::UnknownFormat(_) => None,
        }
    }
}

fn read_structure_file(
    path: &Path,
    format: BioCoordinateFormat,
) -> Result<BioStructureData, BioFileReadError> {
    // Gemmi✔️✔️:   if (format == CoorFormat::Detect) {
    // Gemmi✔️✔️:     CharArray mem = read_into_buffer(input);
    // Gemmi✔️✔️:     return read_structure_from_memory(mem.data(), mem.size(), input.path(), format, save_doc);
    // Gemmi✔️✔️:   }
    // Gemmi❌❌:   if (save_doc)
    // Gemmi❌❌:     save_doc->clear();
    // Gemmi✔️✔️:   if (format == CoorFormat::Unknown)
    // Gemmi✔️✔️:     format = coor_format_from_ext(input.basepath());
    // Gemmi✔️✔️:   switch (format) {
    // Gemmi✔️✔️:     case CoorFormat::Pdb:
    // Gemmi✔️✔️:       return read_pdb(input);
    // Gemmi✔️✔️:     case CoorFormat::Mmcif:
    // Gemmi✔️✔️:       return make_structure(cif::read(input), save_doc);
    // Gemmi✔️✔️:     case CoorFormat::Mmjson: {
    // Gemmi✔️✔️:       Structure st = make_structure(cif::read_mmjson(input), save_doc);
    // Gemmi✔️✔️:       st.input_format = CoorFormat::Mmjson;
    // Gemmi✔️✔️:       return st;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     case CoorFormat::ChemComp:
    // Gemmi✔️✔️:       return make_structure_from_chemcomp_doc(cif::read(input), save_doc);
    // Gemmi✔️✔️:     case CoorFormat::Unknown:
    // Gemmi✔️✔️:     case CoorFormat::Detect:
    // Gemmi✔️✔️:       fail("Unknown format of " +
    // Gemmi✔️✔️:            (input.path().empty() ? "coordinate file" : input.path()) + ".");
    // Gemmi✔️✔️:   }
    // Behavior: unknown file extension errors before opening; Detect loads and
    // routes through memory. Explicit formats never inspect the file contents
    // for a different format. CIF is populated once, without chemcomp probing
    // except in the explicit ChemComp branch; file mmJSON alone overrides format.
    // Complexity: constant-cost selector, one file buffer and one parse; the
    // owned buffer is moved into String, with no cloned hierarchy or reparse.
    let route = select_file_format(path, format);
    if route == BioFileRoute::Format(BioCoordinateFormat::Unknown) {
        return Err(BioFileReadError::UnknownFormat(path.to_path_buf()));
    }
    let text = read_regular_file(path).map_err(BioFileReadError::Load)?;
    let name = path.to_string_lossy();
    match route {
        BioFileRoute::Memory => {
            read_structure_from_memory(&text, &name, BioCoordinateFormat::Detect)
                .map_err(BioFileReadError::Memory)
        }
        BioFileRoute::Format(BioCoordinateFormat::Pdb) => {
            read_pdb_bio_structure(&text, &name, &BioPdbReadParams::default())
                .map_err(BioFileReadError::Pdb)
        }
        BioFileRoute::Format(BioCoordinateFormat::Mmcif | BioCoordinateFormat::ChemComp) => {
            let document = read_cif_document(&text, &name, CifCheckLevel::Default)
                .map_err(BioFileReadError::Cif)?;
            if route == BioFileRoute::Format(BioCoordinateFormat::ChemComp) {
                make_structure_from_chemcomp_doc(&document, 7).map_err(BioFileReadError::ChemComp)
            } else {
                populate_mmcif_bio_structure_document(&document)
                    .map_err(|error| BioFileReadError::Document(BioDocumentError::Mmcif(error)))
            }
        }
        BioFileRoute::Format(BioCoordinateFormat::Mmjson) => {
            let document =
                read_mmjson_insitu(text.as_bytes(), &name).map_err(BioFileReadError::Mmjson)?;
            let mut data = populate_mmcif_bio_structure_document(&document)
                .map_err(|error| BioFileReadError::Document(BioDocumentError::Mmcif(error)))?;
            data.input_format = BioCoordinateFormat::Mmjson;
            Ok(data)
        }
        BioFileRoute::Format(BioCoordinateFormat::Unknown | BioCoordinateFormat::Detect) => {
            Err(BioFileReadError::UnknownFormat(path.to_path_buf()))
        }
    }
}

fn coor_format_from_ext(path: &str) -> BioCoordinateFormat {
    coor_format_from_ext_bytes(path.as_bytes())
}

fn coor_format_from_ext_bytes(path: &[u8]) -> BioCoordinateFormat {
    // Gemmi✔️✔️: inline CoorFormat coor_format_from_ext(const std::string& path) {
    // Gemmi✔️✔️:   if (iends_with(path, ".pdb") || iends_with(path, ".ent"))
    // Gemmi✔️✔️:     return CoorFormat::Pdb;
    // Gemmi✔️✔️:   if (iends_with(path, ".cif") || iends_with(path, ".mmcif"))
    // Gemmi✔️✔️:     return CoorFormat::Mmcif;
    // Gemmi✔️✔️:   if (iends_with(path, ".json"))
    // Gemmi✔️✔️:     return CoorFormat::Mmjson;
    // Gemmi✔️✔️:   return CoorFormat::Unknown;
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: inline bool iends_with(const std::string& str, const std::string& suffix) {
    // Gemmi✔️✔️:   size_t sl = suffix.length();
    // Gemmi✔️✔️:   return str.length() >= sl &&
    // Gemmi✔️✔️:          std::equal(std::begin(suffix), std::end(suffix), str.end() - sl,
    // Gemmi✔️✔️:                     [](char c1, char c2) { return c1 == lower(c2); });
    // Gemmi✔️✔️: }
    // Behavior: only the terminal suffix is compared, with ASCII case folding.
    // Complexity: each fixed suffix checks at most six bytes; no path allocation.
    let ends_with = |suffix: &str| {
        path.get(path.len().saturating_sub(suffix.len())..)
            .is_some_and(|end| end.eq_ignore_ascii_case(suffix.as_bytes()))
    };
    if ends_with(".pdb") || ends_with(".ent") {
        BioCoordinateFormat::Pdb
    } else if ends_with(".cif") || ends_with(".mmcif") {
        BioCoordinateFormat::Mmcif
    } else if ends_with(".json") {
        BioCoordinateFormat::Mmjson
    } else {
        BioCoordinateFormat::Unknown
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BioFileRoute {
    Memory,
    Format(BioCoordinateFormat),
}

fn select_file_format(path: &Path, requested: BioCoordinateFormat) -> BioFileRoute {
    // Gemmi❗✔️:   if (format == CoorFormat::Detect) {
    // Gemmi❌❌:     CharArray mem = read_into_buffer(input);
    // Gemmi❌❌:     return read_structure_from_memory(mem.data(), mem.size(), input.path(), format, save_doc);
    // Gemmi❗✔️:   }
    // Gemmi❌❌:   if (save_doc)
    // Gemmi❌❌:     save_doc->clear();
    // Gemmi✔️✔️:   if (format == CoorFormat::Unknown)
    // Gemmi✔️✔️:     format = coor_format_from_ext(input.basepath());
    // Behavior: Detect remains the distinct memory route. Unknown alone uses
    // the exact terminal ASCII suffix over OS path bytes; unknown extensions
    // remain Unknown for the later source-defined error branch.
    // Complexity: one constant-size suffix check; no path allocation or
    // lossy conversion, and non-UTF-8 parent components do not affect it.
    match requested {
        BioCoordinateFormat::Detect => BioFileRoute::Memory,
        BioCoordinateFormat::Unknown => BioFileRoute::Format(coor_format_from_ext_bytes(
            path.as_os_str().as_encoded_bytes(),
        )),
        format => BioFileRoute::Format(format),
    }
}

fn coor_format_from_content(bytes: &[u8]) -> BioCoordinateFormat {
    // Gemmi✔️✔️: inline CoorFormat coor_format_from_content(const char* buf, const char* end) {
    // Gemmi✔️✔️:   while (buf < end - 8) {
    // Gemmi✔️✔️:     if (std::isspace(*buf)) {
    // Gemmi✔️✔️:       ++buf;
    // Gemmi✔️✔️:     } else if (*buf == '#') {
    // Gemmi✔️✔️:       while (buf < end - 8 && *buf != '\n')
    // Gemmi✔️✔️:         ++buf;
    // Gemmi✔️✔️:     } else if (*buf == '{') {
    // Gemmi✔️✔️:       return CoorFormat::Mmjson;
    // Gemmi✔️✔️:     } else if (ialpha4_id(buf) == ialpha4_id("data") && buf[4] == '_') {
    // Gemmi✔️✔️:       return CoorFormat::Mmcif;
    // Gemmi✔️✔️:     } else {
    // Gemmi✔️✔️:       return CoorFormat::Pdb;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return CoorFormat::Unknown;
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: constexpr int ialpha4_id(const char* s) {
    // Gemmi✔️✔️:   return (s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) & ~0x20202020;
    // Gemmi✔️✔️: }
    // Behavior: preserve the >8 remaining-byte guard and C-locale ASCII
    // whitespace, including vertical tab. A short tail is not classified.
    // Complexity: a single forward scan with constant-size prefix inspection.
    let mut position = 0usize;
    while bytes.len() - position > 8 {
        match bytes[position] {
            b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c => position += 1,
            b'#' => {
                while bytes.len() - position > 8 && bytes[position] != b'\n' {
                    position += 1;
                }
            }
            b'{' => return BioCoordinateFormat::Mmjson,
            _ => {
                if bytes[position..position + 4].eq_ignore_ascii_case(b"data")
                    && bytes[position + 4] == b'_'
                {
                    return BioCoordinateFormat::Mmcif;
                }
                return BioCoordinateFormat::Pdb;
            }
        }
    }
    BioCoordinateFormat::Unknown
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cif::{CifCheckLevel, read_cif_document};
    use std::fs;

    #[test]
    fn bio_read_m01_chemcomp_switch_selects_special_block_and_exact_format() {
        // mmread.hpp checks for the chemcomp block only when requested;
        // mmcif.hpp::make_structure otherwise populates block zero.
        let atom_block = "data_comp_CMP\n_chem_comp.id 'SELECTED'\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n";
        for prefix in [
            "",
            "data_comp_list\n_chem_comp.id WRONG\n",
            "global_\ndata_comp_list\n_chem_comp.id WRONG\n",
        ] {
            let text = format!("{prefix}{atom_block}");
            let document = read_cif_document(&text, "m01.cif", CifCheckLevel::Syntax).unwrap();
            let special = make_structure_from_doc(&document, true).unwrap();
            assert_eq!(special.input_format(), BioCoordinateFormat::ChemComp);
            assert_eq!(special.source_state().name, "'SELECTED'");
            assert_eq!(special.coordinates().positions(), &[[1.0, 2.0, 3.0]]);
            assert_eq!(special.models().len(), 1);

            if prefix.is_empty() {
                let ordinary = make_structure_from_doc(&document, false).unwrap();
                assert_eq!(ordinary.input_format(), BioCoordinateFormat::Mmcif);
                assert_eq!(ordinary.source_state().name, "comp_CMP");
            }
        }
    }

    #[test]
    fn bio_read_m01_nonchemcomp_document_uses_ordinary_owner_and_keeps_errors() {
        let text = "data_coordinates\n_entry.id DEMO\n";
        let document = read_cif_document(text, "m01.cif", CifCheckLevel::Syntax).unwrap();
        for possible in [false, true] {
            let structure = make_structure_from_doc(&document, possible).unwrap();
            assert_eq!(structure.input_format(), BioCoordinateFormat::Mmcif);
            assert_eq!(structure.source_state().name, "coordinates");
        }

        let bad = read_cif_document(
            "data_first\ndata_second\n_atom_site.id 1\n",
            "later.cif",
            CifCheckLevel::Syntax,
        )
        .unwrap();
        for possible in [false, true] {
            let error = make_structure_from_doc(&bad, possible).unwrap_err();
            assert!(matches!(error, BioDocumentError::Mmcif(_)));
            assert!(error.to_string().contains("block #2: later.cif"));
            assert!(error.source().is_some());
        }
    }

    #[test]
    fn bio_read_m02_six_formats_cross_all_memory_content_families() {
        // mmread.hpp selects a reader by explicit format first; Unknown and
        // Detect classify bytes. Memory ChemComp is never a selected branch.
        let inputs = [
            ("pdb", "HEADER    TEST\n", BioCoordinateFormat::Pdb),
            (
                "cif",
                "data_coordinates\n_entry.id DEMO\n",
                BioCoordinateFormat::Mmcif,
            ),
            (
                "json",
                r#"{"data_coordinates":{"entry":{"id":["DEMO"]}}}"#,
                BioCoordinateFormat::Mmjson,
            ),
            (
                "chemcomp",
                "data_comp_CMP\n_chem_comp.id CMP\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n",
                BioCoordinateFormat::Mmcif,
            ),
            ("short", "", BioCoordinateFormat::Unknown),
        ];
        for requested in [
            BioCoordinateFormat::Unknown,
            BioCoordinateFormat::Detect,
            BioCoordinateFormat::Pdb,
            BioCoordinateFormat::Mmcif,
            BioCoordinateFormat::Mmjson,
            BioCoordinateFormat::ChemComp,
        ] {
            for (name, text, detected) in inputs {
                let source_name = format!("{name}.memory");
                let selected = if matches!(
                    requested,
                    BioCoordinateFormat::Unknown | BioCoordinateFormat::Detect
                ) {
                    detected
                } else {
                    requested
                };
                let result = read_structure_from_memory(text, &source_name, requested);
                match selected {
                    BioCoordinateFormat::ChemComp | BioCoordinateFormat::Unknown => {
                        assert!(
                            matches!(result, Err(BioMemoryReadError::WrongFormat { source_name: ref actual, format }) if actual == &source_name && format == selected),
                            "{name} {requested:?}"
                        );
                    }
                    BioCoordinateFormat::Pdb => match result {
                        Ok(structure) => assert_eq!(
                            structure.input_format(),
                            BioCoordinateFormat::Pdb,
                            "{name} {requested:?}"
                        ),
                        Err(BioMemoryReadError::Pdb(error)) => assert!(error.source().is_some()),
                        other => panic!("PDB branch mismatch {name} {requested:?}: {other:?}"),
                    },
                    BioCoordinateFormat::Mmcif => match result {
                        Ok(structure) => assert_eq!(
                            structure.input_format(),
                            if name == "chemcomp" {
                                BioCoordinateFormat::ChemComp
                            } else {
                                BioCoordinateFormat::Mmcif
                            },
                            "{name} {requested:?}"
                        ),
                        Err(BioMemoryReadError::Cif(error)) => {
                            assert_eq!(error.source(), source_name)
                        }
                        Err(BioMemoryReadError::Document(error)) => {
                            assert!(error.source().is_some())
                        }
                        other => panic!("mmCIF branch mismatch {name} {requested:?}: {other:?}"),
                    },
                    BioCoordinateFormat::Mmjson => match result {
                        Ok(structure) => assert_eq!(
                            structure.input_format(),
                            BioCoordinateFormat::Mmcif,
                            "memory mmJSON must retain ordinary make_structure format"
                        ),
                        Err(BioMemoryReadError::Mmjson(error)) => {
                            assert!(error.to_string().contains(&source_name))
                        }
                        Err(BioMemoryReadError::Document(error)) => {
                            assert!(error.source().is_some())
                        }
                        other => panic!("mmJSON branch mismatch {name} {requested:?}: {other:?}"),
                    },
                    BioCoordinateFormat::Detect => unreachable!(),
                }
            }
        }
    }

    #[test]
    fn bio_read_m02_matching_memory_inputs_preserve_source_name_and_format() {
        let pdb = read_structure_from_memory(
            "HEADER    TEST\n",
            "named.pdb",
            BioCoordinateFormat::Unknown,
        )
        .unwrap();
        assert_eq!(pdb.input_format(), BioCoordinateFormat::Pdb);
        let cif = read_structure_from_memory(
            "data_demo\n_entry.id DEMO\n",
            "named.cif",
            BioCoordinateFormat::Unknown,
        )
        .unwrap();
        assert_eq!(cif.input_format(), BioCoordinateFormat::Mmcif);
        assert_eq!(cif.source_state().name, "demo");
        let json = read_structure_from_memory(
            r#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#,
            "named.json",
            BioCoordinateFormat::Unknown,
        )
        .unwrap();
        assert_eq!(json.input_format(), BioCoordinateFormat::Mmcif);
        assert_eq!(json.source_state().name, "demo");
        let chemcomp = read_structure_from_memory("data_comp_CMP\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n", "named.cif", BioCoordinateFormat::Unknown).unwrap();
        assert_eq!(chemcomp.input_format(), BioCoordinateFormat::ChemComp);
        assert_eq!(chemcomp.models().len(), 1);
    }

    #[test]
    fn bio_read_f01_file_loading_retains_exact_text_and_path() {
        let dir = tempfile::tempdir().unwrap();
        let ordinary = dir.path().join("entry.pdb");
        let misleading = dir.path().join("entry.cif");
        let empty = dir.path().join("empty.json");
        fs::write(&ordinary, "HEADER    café\n").unwrap();
        fs::write(&misleading, "HEADER    TEST\n").unwrap();
        fs::write(&empty, b"").unwrap();
        assert_eq!(read_regular_file(&ordinary).unwrap(), "HEADER    café\n");
        assert_eq!(read_regular_file(&misleading).unwrap(), "HEADER    TEST\n");
        assert_eq!(read_regular_file(&empty).unwrap(), "");

        let missing = dir.path().join("missing.pdb");
        let error = read_regular_file(&missing).unwrap_err();
        match error {
            BioFileLoadError::Io { path, source } => {
                assert_eq!(path, missing);
                assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
            }
            other => panic!("wrong file error: {other:?}"),
        }
    }

    #[test]
    fn bio_read_f01_invalid_utf8_is_structured_and_never_replaced() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("invalid.cif");
        let bytes = b"data_demo\n\xff\n";
        fs::write(&path, bytes).unwrap();
        let error = read_regular_file(&path).unwrap_err();
        assert!(error.to_string().contains(path.to_str().unwrap()));
        match error {
            BioFileLoadError::Utf8 {
                path: actual,
                source,
            } => {
                assert_eq!(actual, path);
                assert_eq!(source.as_bytes(), bytes);
                assert_eq!(source.utf8_error().valid_up_to(), 10);
            }
            other => panic!("wrong decoding error: {other:?}"),
        }
    }

    #[test]
    fn bio_read_f02_unknown_is_extension_only_detect_is_content_route() {
        // mmread.hpp::read_structure tests Detect before extension selection;
        // Unknown selects by basepath only, not by the loaded contents.
        let paths = [
            ("entry.pdb", BioCoordinateFormat::Pdb),
            ("entry.ENT", BioCoordinateFormat::Pdb),
            ("entry.MmCiF", BioCoordinateFormat::Mmcif),
            ("entry.CIF", BioCoordinateFormat::Mmcif),
            ("entry.JsOn", BioCoordinateFormat::Mmjson),
            ("entry.xyz", BioCoordinateFormat::Unknown),
            ("entry.cif.gz", BioCoordinateFormat::Unknown),
        ];
        let contents = [
            ("HEADER    TEST\n", BioCoordinateFormat::Pdb),
            ("data_demo\n_entry.id DEMO\n", BioCoordinateFormat::Mmcif),
            (r#"{"data_demo":{}}"#, BioCoordinateFormat::Mmjson),
            ("", BioCoordinateFormat::Unknown),
        ];
        for (path, extension) in paths {
            for (content, detected) in contents {
                assert_eq!(
                    select_file_format(Path::new(path), BioCoordinateFormat::Unknown),
                    BioFileRoute::Format(extension),
                    "Unknown must ignore {detected:?} content in {path}"
                );
                assert_eq!(
                    select_file_format(Path::new(path), BioCoordinateFormat::Detect),
                    BioFileRoute::Memory,
                    "Detect must ignore {extension:?} extension of {path}"
                );
                assert_eq!(coor_format_from_content(content.as_bytes()), detected);
            }
        }
    }

    #[test]
    fn bio_read_f02_explicit_format_is_never_reclassified() {
        for path in ["entry.pdb", "entry.json", "entry.unknown"] {
            for format in [
                BioCoordinateFormat::Pdb,
                BioCoordinateFormat::Mmcif,
                BioCoordinateFormat::Mmjson,
                BioCoordinateFormat::ChemComp,
            ] {
                assert_eq!(
                    select_file_format(Path::new(path), format),
                    BioFileRoute::Format(format)
                );
            }
        }
    }

    #[test]
    fn bio_read_f03_file_switch_crosses_formats_and_content_families() {
        // mmread.hpp::read_structure switches on explicit format; Detect alone
        // enters memory detection. File Mmcif does not recognize chemcomp;
        // file Mmjson overrides the ordinary structure input_format.
        let dir = tempfile::tempdir().unwrap();
        let inputs = [
            ("pdb", "HEADER    TEST\n"),
            ("cif", "data_demo\n_entry.id DEMO\n"),
            ("json", r#"{"data_demo":{}}"#),
            (
                "chemcomp",
                "data_comp_CMP\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n",
            ),
            ("invalid", ""),
        ];
        for (content_name, text) in inputs {
            let path = dir.path().join(format!("{content_name}.xyz"));
            fs::write(&path, text).unwrap();
            for format in [
                BioCoordinateFormat::Unknown,
                BioCoordinateFormat::Detect,
                BioCoordinateFormat::Pdb,
                BioCoordinateFormat::Mmcif,
                BioCoordinateFormat::Mmjson,
                BioCoordinateFormat::ChemComp,
            ] {
                let result = read_structure_file(&path, format);
                if format == BioCoordinateFormat::Unknown {
                    assert!(
                        matches!(result, Err(BioFileReadError::UnknownFormat(_))),
                        "{content_name}"
                    );
                    continue;
                }
                let expected = match (format, content_name) {
                    (BioCoordinateFormat::Detect, "pdb") | (BioCoordinateFormat::Pdb, "pdb") => {
                        Some(BioCoordinateFormat::Pdb)
                    }
                    (BioCoordinateFormat::Detect, "cif") | (BioCoordinateFormat::Mmcif, "cif") => {
                        Some(BioCoordinateFormat::Mmcif)
                    }
                    (BioCoordinateFormat::Detect, "json") => Some(BioCoordinateFormat::Mmcif),
                    (BioCoordinateFormat::Mmjson, "json") => Some(BioCoordinateFormat::Mmjson),
                    (BioCoordinateFormat::Detect, "chemcomp")
                    | (BioCoordinateFormat::ChemComp, "chemcomp") => {
                        Some(BioCoordinateFormat::ChemComp)
                    }
                    (BioCoordinateFormat::Mmcif, "chemcomp") => Some(BioCoordinateFormat::Mmcif),
                    _ => None,
                };
                if let Some(expected) = expected {
                    assert_eq!(
                        result.unwrap().input_format(),
                        expected,
                        "{format:?}/{content_name}"
                    );
                } else if format == BioCoordinateFormat::Pdb {
                    // read_pdb accepts some unrecognized records (and empty
                    // input) as a valid empty structure; it does not enforce
                    // that the extension/content classify as PDB.
                    match result {
                        Ok(data) => assert_eq!(data.input_format(), BioCoordinateFormat::Pdb),
                        Err(BioFileReadError::Pdb(error)) => assert!(error.source().is_some()),
                        other => panic!("wrong PDB branch for {content_name}: {other:?}"),
                    }
                } else {
                    assert!(result.is_err(), "{format:?}/{content_name}");
                }
            }
        }
    }

    #[test]
    fn bio_read_f03_file_errors_preserve_origin_and_extension_precedence() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("missing.xyz");
        assert!(matches!(
            read_structure_file(&missing, BioCoordinateFormat::Unknown),
            Err(BioFileReadError::UnknownFormat(_))
        ));
        let missing = dir.path().join("missing.cif");
        let error = read_structure_file(&missing, BioCoordinateFormat::Unknown).unwrap_err();
        assert!(matches!(
            error,
            BioFileReadError::Load(BioFileLoadError::Io { .. })
        ));
        assert_eq!(
            error
                .source()
                .unwrap()
                .source()
                .unwrap()
                .downcast_ref::<std::io::Error>()
                .unwrap()
                .kind(),
            std::io::ErrorKind::NotFound
        );
        let invalid = dir.path().join("wrong.cif");
        fs::write(&invalid, "not a cif format\n").unwrap();
        let error = read_structure_file(&invalid, BioCoordinateFormat::Unknown).unwrap_err();
        assert!(matches!(error, BioFileReadError::Cif(_)));
        assert!(error.source().is_some());
        let invalid = dir.path().join("wrong.json");
        fs::write(&invalid, "not json").unwrap();
        let error = read_structure_file(&invalid, BioCoordinateFormat::Unknown).unwrap_err();
        assert!(matches!(error, BioFileReadError::Mmjson(_)));
        assert!(error.source().is_some());
    }

    #[test]
    fn bio_read_params_defaults_route_and_preserve_source() {
        let params = BioReadParams::default();
        assert_eq!(params.format, BioCoordinateFormat::Unknown);
        assert_eq!(params.source_name, "<string>");
        for (text, expected) in [
            ("HEADER    TEST\n", BioCoordinateFormat::Pdb),
            ("data_demo\n_entry.id DEMO\n", BioCoordinateFormat::Mmcif),
            (r#"{"data_demo":{}}"#, BioCoordinateFormat::Mmcif),
            (
                "data_comp_CMP\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n",
                BioCoordinateFormat::ChemComp,
            ),
        ] {
            assert_eq!(
                read_bio_structure(text, &params).unwrap().input_format(),
                expected
            );
        }
        let custom = BioReadParams {
            format: BioCoordinateFormat::Mmcif,
            source_name: "named.cif".into(),
        };
        let data = read_bio_structure("data_custom\n_entry.id NAMED\n", &custom).unwrap();
        assert_eq!(data.input_format(), BioCoordinateFormat::Mmcif);
        let error = read_bio_structure("broken cif", &custom).unwrap_err();
        assert!(matches!(error, BioReadError::Cif(_)));
        assert!(error.to_string().contains("named.cif"));
        let explicit = BioReadParams {
            format: BioCoordinateFormat::ChemComp,
            source_name: "special".into(),
        };
        assert!(
            matches!(read_bio_structure("data_demo\n_entry.id DEMO\n", &explicit), Err(BioReadError::WrongFormat { source_name, format: BioCoordinateFormat::ChemComp }) if source_name == "special")
        );
    }

    #[test]
    fn bio_read_params_file_errors_and_distinct_format_routes() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("record.xyz");
        fs::write(&path, "data_demo\n_entry.id DEMO\n").unwrap();
        assert!(matches!(
            read_bio_structure_file(&path, BioCoordinateFormat::Unknown),
            Err(BioReadError::UnknownFileFormat(_))
        ));
        assert_eq!(
            read_bio_structure_file(&path, BioCoordinateFormat::Detect)
                .unwrap()
                .input_format(),
            BioCoordinateFormat::Mmcif
        );
        assert_eq!(
            read_bio_structure_file(&path, BioCoordinateFormat::Mmcif)
                .unwrap()
                .input_format(),
            BioCoordinateFormat::Mmcif
        );
        let missing = dir.path().join("missing.cif");
        let error = read_bio_structure_file(&missing, BioCoordinateFormat::Unknown).unwrap_err();
        assert!(matches!(error, BioReadError::Io { .. }));
        assert_eq!(
            error
                .source()
                .unwrap()
                .downcast_ref::<std::io::Error>()
                .unwrap()
                .kind(),
            std::io::ErrorKind::NotFound
        );
        let bad = dir.path().join("bad.cif");
        fs::write(&bad, b"data_demo\n\xff").unwrap();
        let error = read_bio_structure_file(&bad, BioCoordinateFormat::Unknown).unwrap_err();
        assert!(matches!(error, BioReadError::Utf8 { .. }));
        assert_eq!(
            error
                .source()
                .unwrap()
                .downcast_ref::<std::string::FromUtf8Error>()
                .unwrap()
                .as_bytes(),
            b"data_demo\n\xff"
        );
    }

    #[test]
    fn bio_read_d01_terminal_extensions_follow_gemmi_case_folding() {
        // mmread.hpp::coor_format_from_ext calls util.hpp::iends_with on the
        // complete path; neither directory names nor a trailing .gz are removed.
        for (suffix, expected) in [
            ("pdb", BioCoordinateFormat::Pdb),
            ("ent", BioCoordinateFormat::Pdb),
            ("cif", BioCoordinateFormat::Mmcif),
            ("mmcif", BioCoordinateFormat::Mmcif),
            ("json", BioCoordinateFormat::Mmjson),
        ] {
            for spelling in [
                suffix.to_owned(),
                suffix.to_ascii_uppercase(),
                suffix
                    .chars()
                    .enumerate()
                    .map(|(index, ch)| {
                        if index % 2 == 0 {
                            ch.to_ascii_uppercase()
                        } else {
                            ch
                        }
                    })
                    .collect(),
            ] {
                assert_eq!(
                    coor_format_from_ext(&format!("/data/entry.{spelling}")),
                    expected
                );
            }
        }
        for path in [
            "entry",
            "entry.xyz",
            "/data/entry.pdb/record",
            "/data/entry.cif/record.txt",
            "entry.pdb.gz",
            "entry.CIF.GZ",
            "/path.json/entry",
        ] {
            assert_eq!(
                coor_format_from_ext(path),
                BioCoordinateFormat::Unknown,
                "{path}"
            );
        }
    }

    #[test]
    fn bio_read_d02_remaining_length_and_ascii_space_follow_gemmi() {
        // mmread.hpp tests `buf < end - 8` before inspecting any byte.
        for len in 0..=9 {
            let bytes = vec![b'X'; len];
            assert_eq!(
                coor_format_from_content(&bytes),
                if len == 9 {
                    BioCoordinateFormat::Pdb
                } else {
                    BioCoordinateFormat::Unknown
                },
                "length {len}"
            );
        }
        for whitespace in [b' ', b'\t', b'\n', b'\r', 0x0b, 0x0c] {
            let mut bytes = vec![whitespace];
            bytes.extend_from_slice(b"DATA_demo_more");
            assert_eq!(coor_format_from_content(&bytes), BioCoordinateFormat::Mmcif);
        }
        for first in [0x00, 0x01, 0x08, 0x7f] {
            let mut bytes = vec![first];
            bytes.extend_from_slice(b"DATA_demo_more");
            assert_eq!(coor_format_from_content(&bytes), BioCoordinateFormat::Pdb);
        }
    }

    #[test]
    fn bio_read_d02_comments_case_and_json_preserve_source_scan_order() {
        for data in [
            b"data_demo_more".as_slice(),
            b"DATA_demo_more",
            b"DaTa_demo_more",
        ] {
            assert_eq!(coor_format_from_content(data), BioCoordinateFormat::Mmcif);
        }
        for data in [b"{12345678".as_slice(), b" \n{12345678"] {
            assert_eq!(coor_format_from_content(data), BioCoordinateFormat::Mmjson);
        }
        for data in [b"ATOM      1".as_slice(), b"dAtX_demo_more"] {
            assert_eq!(coor_format_from_content(data), BioCoordinateFormat::Pdb);
        }
        assert_eq!(
            coor_format_from_content(b"# comment\nDaTa_demo_more"),
            BioCoordinateFormat::Mmcif
        );
        assert_eq!(
            coor_format_from_content(b"# comment without newline"),
            BioCoordinateFormat::Unknown
        );
        assert_eq!(
            coor_format_from_content(b"# comment\n{12345678"),
            BioCoordinateFormat::Mmjson
        );
    }
}
