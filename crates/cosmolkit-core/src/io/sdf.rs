#![allow(dead_code)]

// Source-reproduction protocol.
//
// This file participates in the COSMolKit source-level porting discipline.
// The full marker convention, valid combinations, upgrade rules, and examples
// are defined in:
//
//   dev/source_reproduction_protocol.md
//
// Every copied C++ line block must carry a two-symbol status marker as
// described in that document. Do not deviate from the protocol unless the
// human author explicitly approves an exception.
//
// 11. Do not implement from memory or from a summarized call graph when the
//     RDKit helper body is available. First create the source-level alignment
//     frame with `BEGIN RDKIT CPP FUNCTION` / `END RDKIT CPP FUNCTION` comments
//     and `RDKit❌❌` markers, then change markers only as each line is
//     reproduced and reviewed. If the source-level frame is missing for a
//     newly-touched RDKit helper, add it before making semantic changes.
//
// Agents are not allowed to bypass, weaken, reinterpret, or delete this
// protocol while porting RDKit-derived code. If following it appears to
// conflict with a practical implementation need, stop and ask the human author
// before editing.

use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use crate::search::query::convert_complex_name_to_query;
use crate::{
    AtomId, AtomQueryPredicate, AtomSpec, BondDirection, BondId, BondOrder, BondQueryPredicate,
    BondSpec, BondStereo, Conformer2D, Conformer3D, CoordinateDimension, Element,
    MOLBLOCK_IO_FEATURE, Molecule, MoleculeBuilder, QueryNode, SGroupAttachPoint, SGroupBracket,
    SGroupBracketStyle, SGroupCState, SGroupConnection, SdfPropertyList, SdfPropertyListTarget,
    StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
    UnsupportedFeatureError,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SdfCoordinateMode {
    Preserve,
    Require2D,
    Require3D,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CtabVersion {
    V2000,
    V3000,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SdfReadParams {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub strict_parsing: bool,
    pub expand_attachment_points: bool,
    pub process_property_lists: bool,
    pub coordinate_mode: SdfCoordinateMode,
}

impl Default for SdfReadParams {
    fn default() -> Self {
        Self {
            sanitize: true,
            remove_hs: true,
            strict_parsing: true,
            expand_attachment_points: false,
            process_property_lists: true,
            coordinate_mode: SdfCoordinateMode::Preserve,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    pub molecule: Molecule,
    pub data_fields: Vec<(String, String)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfRecordMetadata {
    pub index: usize,
    pub byte_offset: u64,
    pub byte_len: u64,
    pub line_offset: usize,
    pub line_len: usize,
    pub title: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfDataset {
    path: PathBuf,
    params: SdfReadParams,
    metadata: Vec<SdfRecordMetadata>,
}

impl SdfDataset {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, SdfReadError> {
        open_sdf_dataset(path.as_ref(), SdfReadParams::default())
    }

    pub fn open_with_params(
        path: impl AsRef<Path>,
        params: SdfReadParams,
    ) -> Result<Self, SdfReadError> {
        open_sdf_dataset(path.as_ref(), params)
    }

    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    #[must_use]
    pub const fn params(&self) -> SdfReadParams {
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

    pub fn record(&self, index: usize) -> Result<SdfRecord, SdfReadError> {
        read_indexed_sdf_record(self, index, self.params)
    }

    pub fn record_with_params(
        &self,
        index: usize,
        params: SdfReadParams,
    ) -> Result<SdfRecord, SdfReadError> {
        read_indexed_sdf_record(self, index, params)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SdfReadError {
    #[error("SDF reading is not implemented")]
    NotImplemented,
    #[error(transparent)]
    UnsupportedFeature(#[from] UnsupportedFeatureError),
    #[error(transparent)]
    Operation(#[from] crate::OperationError),
    #[error("{0}")]
    Parse(String),
}

pub struct SdfReader<R> {
    reader: R,
    params: SdfReadParams,
    next_index: usize,
    byte_offset: u64,
    line_number: usize,
    end: bool,
    eof_hit_on_read: bool,
}

impl<R: BufRead> SdfReader<R> {
    #[must_use]
    pub fn new(reader: R) -> Self {
        Self::with_params(reader, SdfReadParams::default())
    }

    #[must_use]
    pub fn with_coordinate_mode(reader: R, coordinate_mode: SdfCoordinateMode) -> Self {
        Self::with_params(
            reader,
            SdfReadParams {
                coordinate_mode,
                ..Default::default()
            },
        )
    }

    #[must_use]
    pub fn with_params(reader: R, params: SdfReadParams) -> Self {
        Self {
            reader,
            params,
            next_index: 0,
            byte_offset: 0,
            line_number: 0,
            end: false,
            eof_hit_on_read: false,
        }
    }

    pub fn next_record(&mut self) -> Result<Option<SdfRecord>, SdfReadError> {
        read_next_sdf_record(self)
    }

    #[must_use]
    pub(crate) const fn is_end(&self) -> bool {
        self.end
    }
}

pub fn read_sdf_from_str(s: &str) -> Result<SdfRecord, SdfReadError> {
    read_sdf_from_str_with_params(s, SdfReadParams::default())
}

pub fn read_sdf_from_str_with_coordinate_mode(
    s: &str,
    coordinate_mode: SdfCoordinateMode,
) -> Result<SdfRecord, SdfReadError> {
    read_sdf_from_str_with_params(
        s,
        SdfReadParams {
            coordinate_mode,
            ..Default::default()
        },
    )
}

pub fn read_sdf_from_str_with_params(
    s: &str,
    params: SdfReadParams,
) -> Result<SdfRecord, SdfReadError> {
    parse_sdf_record_text(s, 0, params)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RawSdfRecord {
    index: usize,
    byte_offset: u64,
    byte_len: u64,
    start_line: usize,
    hit_eof: bool,
    text: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SdfDataField {
    name: String,
    value: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MolHeader {
    name: Option<String>,
    info: String,
    comments: String,
    ctab_version: CtabVersion,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CountsLine {
    atom_count: usize,
    bond_count: usize,
    chiral_flag: u32,
    ctab_version: CtabVersion,
    v3000_sgroup_count: usize,
    v3000_obj3d_count: usize,
}

#[derive(Debug, Clone, PartialEq)]
struct ParsedMolBlock {
    molecule: Molecule,
    header: MolHeader,
    chirality_possible: bool,
}

#[derive(Debug, Clone, PartialEq)]
struct V2000AtomLine {
    line_number: usize,
    text: String,
    spec: AtomSpec,
    coord_3d: [f64; 3],
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct V2000BondLine {
    line_number: usize,
    text: String,
    spec: BondSpec,
    molfile_bond_type: u32,
    molfile_stereo: Option<u32>,
}

#[derive(Debug, Default)]
struct V2000PropertyState {
    first_charge_line: bool,
    sgroups: BTreeMap<u32, SubstanceGroup>,
    scd_counter: u32,
    last_data_sgroup: u32,
    current_data_field: String,
    molecule_props: BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq)]
struct V3000AtomLine {
    line_number: usize,
    mol_idx: u32,
    tokens: Vec<String>,
    spec: AtomSpec,
    coord_3d: [f64; 3],
}

#[derive(Debug, Clone, PartialEq)]
struct V3000BondLine {
    line_number: usize,
    mol_idx: u32,
    tokens: Vec<String>,
    begin_mol_idx: u32,
    end_mol_idx: u32,
    spec: BondSpec,
    molfile_bond_type: u32,
}

fn unsupported<T>() -> Result<T, SdfReadError> {
    Err(UnsupportedFeatureError::from_spec(&MOLBLOCK_IO_FEATURE).into())
}

fn unsupported_feature<T>(feature: &'static crate::FeatureSpec) -> Result<T, SdfReadError> {
    Err(UnsupportedFeatureError::from_spec(feature).into())
}

fn molecule_operation_error(err: crate::OperationError) -> SdfReadError {
    match err {
        crate::OperationError::UnsupportedFeature { source, .. } => source.into(),
        other => other.into(),
    }
}

fn molecule_build_error(err: impl std::fmt::Display) -> SdfReadError {
    SdfReadError::Parse(err.to_string())
}

fn strip_sdf_line(line: &str) -> &str {
    line.trim_matches([' ', '\t', '\r', '\n'])
}

fn strip_terminal_cr(line: &str) -> &str {
    line.strip_suffix('\r').unwrap_or(line)
}

fn read_rdkit_line<R: BufRead>(reader: &mut R) -> Result<Option<String>, SdfReadError> {
    let mut line = String::new();
    let read = reader
        .read_line(&mut line)
        .map_err(|err| SdfReadError::Parse(err.to_string()))?;
    if read == 0 {
        return Ok(None);
    }
    if line.ends_with('\n') {
        line.pop();
    }
    if line.ends_with('\r') {
        line.pop();
    }
    Ok(Some(line))
}

fn rdkit_substr(text: &str, start: usize, len: usize) -> &str {
    if start >= text.len() {
        return "";
    }
    let end = start.saturating_add(len).min(text.len());
    &text[start..end]
}

fn parse_rdkit_unsigned(text: &str, accept_spaces: bool) -> Result<u32, ()> {
    for byte in text.bytes().take_while(|byte| *byte != b'\0') {
        if byte.is_ascii_digit() || byte == b'+' || (accept_spaces && byte == b' ') {
            continue;
        }
        return Err(());
    }

    let mut input = text;
    if accept_spaces {
        input = input.trim_start_matches(' ');
        if input.is_empty() {
            return Ok(0);
        }
    }
    if let Some(rest) = input.strip_prefix('+') {
        input = rest;
    }

    let digits_len = input.bytes().take_while(u8::is_ascii_digit).count();
    if digits_len == 0 {
        return Ok(0);
    }
    Ok(input[..digits_len].parse().unwrap_or(0))
}

fn parse_rdkit_int(text: &str, accept_spaces: bool) -> Result<i32, ()> {
    for byte in text.bytes().take_while(|byte| *byte != b'\0') {
        if byte.is_ascii_digit() || byte == b'+' || byte == b'-' || (accept_spaces && byte == b' ')
        {
            continue;
        }
        return Err(());
    }

    let mut input = text;
    if accept_spaces {
        input = input.trim_start_matches(' ');
        if input.is_empty() {
            return Ok(0);
        }
    }
    let sign_len = usize::from(input.starts_with(['+', '-']));
    let digits_len = input[sign_len..]
        .bytes()
        .take_while(u8::is_ascii_digit)
        .count();
    if digits_len == 0 {
        return Ok(0);
    }
    Ok(input[..sign_len + digits_len].parse().unwrap_or(0))
}

fn parse_required_int_field(
    line: &str,
    start: usize,
    len: usize,
    line_number: usize,
) -> Result<i32, SdfReadError> {
    let field = rdkit_substr(line, start, len);
    parse_rdkit_int(field, true).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{field}' to int on line {line_number}"
        ))
    })
}

fn parse_required_unsigned_field(
    line: &str,
    start: usize,
    len: usize,
    line_number: usize,
) -> Result<u32, SdfReadError> {
    let field = rdkit_substr(line, start, len);
    parse_rdkit_unsigned(field, true).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{field}' to int on line {line_number}"
        ))
    })
}

fn atom_line_mut(
    atoms: &mut [V2000AtomLine],
    one_based_atom_id: u32,
    line_number: usize,
) -> Result<&mut V2000AtomLine, SdfReadError> {
    let index = one_based_atom_id.checked_sub(1).ok_or_else(|| {
        SdfReadError::Parse(format!("Atom index 0 out of range on line {line_number}"))
    })? as usize;
    atoms.get_mut(index).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Atom index {one_based_atom_id} out of range on line {line_number}"
        ))
    })
}

fn parse_rdkit_double(text: &str, accept_spaces: bool) -> Result<f64, ()> {
    for byte in text.bytes().take_while(|byte| *byte != b'\0') {
        if byte.is_ascii_digit()
            || byte == b'+'
            || byte == b'-'
            || byte == b'.'
            || byte == b','
            || (accept_spaces && byte == b' ')
        {
            continue;
        }
        return Err(());
    }
    let input = if accept_spaces {
        text.trim_start_matches(' ')
    } else {
        text
    };
    Ok(input.parse().unwrap_or(0.0))
}

macro_rules! mol_symbol_atomic_numbers {
    (
        zero_aliases: [$($zero_alias:literal),+ $(,)?],
        hydrogen_aliases: [$($hydrogen_alias:literal),+ $(,)?],
        elements: [$($symbol:literal => $atomic_number:literal),+ $(,)?] $(,)?
    ) => {
        fn atomic_number_from_mol_symbol(
            symbol: &str,
            strict_parsing: bool,
        ) -> Result<u8, SdfReadError> {
            fn exact_atomic_number(symbol: &str) -> Option<u8> {
                let atomic_number = match symbol {
                $($zero_alias)|+ => 0,
                $($hydrogen_alias)|+ => 1,
                $($symbol => $atomic_number,)+
                    _ => return None,
                };
                Some(atomic_number)
            }

            let normalized_symbol;
            let lookup_symbol = if exact_atomic_number(symbol).is_none()
                && symbol.len() == 2
                && symbol.as_bytes()[1].is_ascii_uppercase()
            {
                // RDKit✔️✔️: std::string tCopy(symb);
                // RDKit✔️✔️: if (symb.size() == 2 && symb[1] >= 'A' && symb[1] <= 'Z') {
                // RDKit✔️✔️:   tCopy[1] = static_cast<char>(tolower(symb[1]));
                // RDKit✔️✔️: }
                normalized_symbol = format!(
                    "{}{}",
                    &symbol[..1],
                    (symbol.as_bytes()[1] as char).to_ascii_lowercase()
                );
                normalized_symbol.as_str()
            } else {
                symbol
            };

            let atomic_number = match exact_atomic_number(lookup_symbol) {
                Some(atomic_number) => atomic_number,
                None if !strict_parsing => 0,
                _ => {
                    return Err(SdfReadError::Parse(format!(
                        "Element '{}' not found",
                        symbol
                    )));
                }
            };
            Ok(atomic_number)
        }
    };
}

mol_symbol_atomic_numbers! {
    zero_aliases: ["*", "A", "Q", "L", "LP", "R", "R#"],
    hydrogen_aliases: ["H", "D", "T"],
    elements: [
        "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9,
        "Ne" => 10, "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P" => 15, "S" => 16,
        "Cl" => 17, "Ar" => 18, "K" => 19, "Ca" => 20, "Sc" => 21, "Ti" => 22, "V" => 23,
        "Cr" => 24, "Mn" => 25, "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30,
        "Ga" => 31, "Ge" => 32, "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36, "Rb" => 37,
        "Sr" => 38, "Y" => 39, "Zr" => 40, "Nb" => 41, "Mo" => 42, "Tc" => 43, "Ru" => 44,
        "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48, "In" => 49, "Sn" => 50, "Sb" => 51,
        "Te" => 52, "I" => 53, "Xe" => 54, "Cs" => 55, "Ba" => 56, "La" => 57, "Ce" => 58,
        "Pr" => 59, "Nd" => 60, "Pm" => 61, "Sm" => 62, "Eu" => 63, "Gd" => 64, "Tb" => 65,
        "Dy" => 66, "Ho" => 67, "Er" => 68, "Tm" => 69, "Yb" => 70, "Lu" => 71, "Hf" => 72,
        "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77, "Pt" => 78, "Au" => 79,
        "Hg" => 80, "Tl" => 81, "Pb" => 82, "Bi" => 83, "Po" => 84, "At" => 85, "Rn" => 86,
        "Fr" => 87, "Ra" => 88, "Ac" => 89, "Th" => 90, "Pa" => 91, "U" => 92, "Np" => 93,
        "Pu" => 94, "Am" => 95, "Cm" => 96, "Bk" => 97, "Cf" => 98, "Es" => 99, "Fm" => 100,
        "Md" => 101, "No" => 102, "Lr" => 103, "Rf" => 104, "Db" => 105, "Sg" => 106,
        "Bh" => 107, "Hs" => 108, "Mt" => 109, "Ds" => 110, "Rg" => 111, "Cn" => 112,
        "Nh" => 113, "Fl" => 114, "Mc" => 115, "Lv" => 116, "Ts" => 117, "Og" => 118,
    ],
}

fn element_from_query_atomic_number(
    atomic_number: u8,
    line_number: usize,
) -> Result<Element, SdfReadError> {
    Element::from_atomic_number(atomic_number).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Unsupported atom-list atomic number {atomic_number} on line {line_number}"
        ))
    })
}

fn most_common_isotope(atomic_number: u8) -> Option<i32> {
    match atomic_number {
        1 => Some(1),
        2 => Some(4),
        3 => Some(7),
        4 => Some(9),
        5 => Some(11),
        6 => Some(12),
        7 => Some(14),
        8 => Some(16),
        9 => Some(19),
        10 => Some(20),
        11 => Some(23),
        12 => Some(24),
        13 => Some(27),
        14 => Some(28),
        15 => Some(31),
        16 => Some(32),
        17 => Some(35),
        18 => Some(40),
        19 => Some(39),
        20 => Some(40),
        21 => Some(45),
        22 => Some(48),
        23 => Some(51),
        24 => Some(52),
        25 => Some(55),
        26 => Some(56),
        27 => Some(59),
        28 => Some(58),
        29 => Some(63),
        30 => Some(64),
        31 => Some(69),
        32 => Some(74),
        33 => Some(75),
        34 => Some(80),
        36 => Some(84),
        35 => Some(79),
        37 => Some(85),
        38 => Some(88),
        39 => Some(89),
        40 => Some(90),
        41 => Some(93),
        42 => Some(98),
        43 => Some(98),
        44 => Some(102),
        45 => Some(103),
        46 => Some(106),
        47 => Some(107),
        48 => Some(114),
        49 => Some(115),
        50 => Some(120),
        51 => Some(121),
        52 => Some(130),
        54 => Some(132),
        53 => Some(127),
        55 => Some(133),
        56 => Some(138),
        57 => Some(139),
        58 => Some(140),
        59 => Some(141),
        60 => Some(142),
        61 => Some(145),
        62 => Some(152),
        63 => Some(153),
        64 => Some(158),
        65 => Some(159),
        66 => Some(164),
        67 => Some(165),
        68 => Some(166),
        69 => Some(169),
        70 => Some(174),
        71 => Some(175),
        72 => Some(180),
        73 => Some(181),
        74 => Some(184),
        75 => Some(187),
        76 => Some(192),
        77 => Some(193),
        78 => Some(195),
        79 => Some(197),
        80 => Some(202),
        81 => Some(205),
        82 => Some(208),
        83 => Some(209),
        84 => Some(209),
        85 => Some(210),
        86 => Some(222),
        87 => Some(223),
        88 => Some(226),
        89 => Some(227),
        90 => Some(232),
        91 => Some(231),
        92 => Some(238),
        93 => Some(237),
        94 => Some(244),
        95 => Some(243),
        96 => Some(247),
        97 => Some(247),
        98 => Some(251),
        99 => Some(252),
        100 => Some(257),
        101 => Some(258),
        102 => Some(259),
        103 => Some(262),
        104 => Some(267),
        105 => Some(268),
        106 => Some(269),
        107 => Some(270),
        108 => Some(269),
        109 => Some(278),
        110 => Some(281),
        111 => Some(282),
        112 => Some(285),
        113 => Some(286),
        114 => Some(289),
        115 => Some(290),
        116 => Some(293),
        117 => Some(294),
        118 => Some(294),
        _ => None,
    }
}

fn starts_with_sdf_continuation_space(line: &str) -> bool {
    line.starts_with(' ') || line.starts_with('\t')
}

fn is_sdf_record_delimiter(line: &str) -> bool {
    line.starts_with("$$$$")
}

fn open_sdf_dataset(path: &Path, params: SdfReadParams) -> Result<SdfDataset, SdfReadError> {
    // BEGIN RDKIT CPP BODY: open_sdf_dataset
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: SDMolSupplier::SDMolSupplier(const std::string &fileName
    // RDKit✔️❌:
    // RDKit✔️❌:   init();
    // RDKit✔️❌:   dp_inStream = openAndCheckStream(fileName);
    // RDKit✔️❌:   df_owner = true;
    // RDKit✔️❌:   d_molpos.push_back(dp_inStream->tellg());
    // RDKit✔️❌:   d_params = params;
    // RDKit✔️❌:   this->checkForEnd();
    // RDKit✔️❌:   if (df_end) {
    // RDKit✔️❌:     // checkForEnd() sets d_len if we're at EOF. undo that (was GitHub issue
    // RDKit✔️❌:     // 19):
    // RDKit✔️❌:     d_len = 0;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   POSTCONDITION(dp_inStream, "bad instream");
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: SDMolSupplier::SDMolSupplier(const std::string &fileName
    // END RDKIT CPP BODY: open_sdf_dataset

    let file = File::open(path).map_err(|err| SdfReadError::Parse(err.to_string()))?;
    let mut reader = BufReader::new(file);
    let metadata = build_sdf_index(&mut reader, params)?;
    Ok(SdfDataset {
        path: path.to_path_buf(),
        params,
        metadata,
    })
}

fn build_sdf_index<R: BufRead>(
    reader: &mut R,
    params: SdfReadParams,
) -> Result<Vec<SdfRecordMetadata>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: build_sdf_index
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: void SDMolSupplier::buildIndexTo
    // RDKit✔️❌:
    // RDKit✔️❌:   dp_inStream->seekg(d_molpos.back());
    // RDKit✔️❌:   d_last = rdcast<int>(d_molpos.size()) - 1;
    // RDKit✔️❌:   const size_t CHUNK_SIZE = 65536;
    // RDKit✔️❌:   const size_t OVERLAP =
    // RDKit✔️❌:       4;  // to catch "$$$$" at chunk boundaries ("...\n$$ <new chunk> $$...")
    // RDKit✔️❌:   std::vector<char> buffer(CHUNK_SIZE + OVERLAP);
    // RDKit✔️❌:   std::fill(buffer.begin(), buffer.begin() + OVERLAP, '\n');  // safe init
    // RDKit✔️❌:   std::streampos currentStreamPos = dp_inStream->tellg();
    // RDKit✔️❌:   bool foundTarget = false;
    // RDKit✔️❌:
    // RDKit✔️❌:   while (dp_inStream->good() && !foundTarget) {
    // RDKit✔️❌:     std::streampos chunkStartPos = currentStreamPos;
    // RDKit✔️❌:     dp_inStream->read(&buffer[OVERLAP], CHUNK_SIZE);
    // RDKit✔️❌:     std::streamsize bytesRead = dp_inStream->gcount();
    // RDKit✔️❌:     if (bytesRead == 0) {
    // RDKit✔️❌:       break;  // EOF
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     std::streampos chunkEndPos = dp_inStream->tellg();
    // RDKit✔️❌:     // check if the stream is "honest" (binary or text mode with 1 byte newlines
    // RDKit✔️❌:     // (like UNIX), meaning read bytes map 1:1 to disk bytes)
    // RDKit✔️❌:     bool isBinaryLike = (bytesRead == (chunkEndPos - chunkStartPos));
    // RDKit✔️❌:     char *bufStart = &buffer[0];
    // RDKit✔️❌:     char *bufEnd = bufStart + OVERLAP + bytesRead;
    // RDKit✔️❌:     char *ptr = bufStart + 1;
    // RDKit✔️❌:     bool needEOL = false;
    // RDKit✔️❌:
    // RDKit✔️❌:     while (true) {
    // RDKit✔️❌:       constexpr char dollarSigns[]{"$$$$"};
    // RDKit✔️❌:       auto match = std::search(ptr, bufEnd, dollarSigns, dollarSigns + 4);
    // RDKit✔️❌:       if (match == bufEnd) {
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (*(match - 1) == '\n') {  // ensure $$$$ is at start of line
    // RDKit✔️❌:         char *nlPos = match + 4;
    // RDKit✔️❌:         if (nlPos == bufEnd) {
    // RDKit✔️❌:           // corner case, $$$$ is EXACTLY at the end of the buffer
    // RDKit✔️❌:           //  we need the next char in the stream to be a "\n", this is resolved
    // RDKit✔️❌:           //  below.
    // RDKit✔️❌:           needEOL = true;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           while (nlPos < bufEnd && *nlPos != '\n') {
    // RDKit✔️❌:             ++nlPos;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           if (nlPos < bufEnd) {
    // RDKit✔️❌:             ++nlPos;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:
    // RDKit✔️❌:         std::streampos posHold;
    // RDKit✔️❌:         if (isBinaryLike &&
    // RDKit✔️❌:             !needEOL) {  // fast path, math checks out, no need to seek
    // RDKit✔️❌:           posHold = chunkStartPos + std::streamoff(nlPos - bufStart - OVERLAP);
    // RDKit✔️❌:         } else {  // slow path, there is byte translation going on, need to seek
    // RDKit✔️❌:                   // and use the std translation magic to find the actual byte
    // RDKit✔️❌:                   // position
    // RDKit✔️❌:           dp_inStream->clear();
    // RDKit✔️❌:           dp_inStream->seekg(
    // RDKit✔️❌:               chunkStartPos);  // rollback to the start of the chunk
    // RDKit✔️❌:           dp_inStream->ignore(
    // RDKit✔️❌:               nlPos - bufStart -
    // RDKit✔️❌:               OVERLAP);  // advance but with the magic translation in effect now
    // RDKit✔️❌:           posHold =
    // RDKit✔️❌:               dp_inStream
    // RDKit✔️❌:                   ->tellg();  // this is the physical position on disk we want
    // RDKit✔️❌:         }
    // RDKit✔️❌:
    // RDKit✔️❌:         bool atTrueEOF =
    // RDKit✔️❌:             (bytesRead < static_cast<std::streamsize>(CHUNK_SIZE)) &&
    // RDKit✔️❌:             (nlPos >= bufEnd);
    // RDKit✔️❌:         if (!atTrueEOF) {
    // RDKit✔️❌:           if (needEOL) {
    // RDKit✔️❌:             char c = dp_inStream->peek();
    // RDKit✔️❌:             if (c == '\n') {
    // RDKit✔️❌:               posHold = posHold + std::streamoff(1);
    // RDKit✔️❌:               needEOL = false;
    // RDKit✔️❌:             }
    // RDKit✔️❌:             this->checkForEnd();
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             this->peekCheckForEnd(nlPos, bufEnd,
    // RDKit✔️❌:                                   posHold);  // the optimized peek version
    // RDKit✔️❌:           }
    // RDKit✔️❌:           if (!this->df_end) {
    // RDKit✔️❌:             d_molpos.push_back(posHold);
    // RDKit✔️❌:             ++d_last;
    // RDKit✔️❌:             if (static_cast<unsigned int>(d_last) ==
    // RDKit✔️❌:                 targetIdx) {  // not really needed but this way we only index as
    // RDKit✔️❌:                               // much as needed
    // RDKit✔️❌:               foundTarget = true;
    // RDKit✔️❌:               break;
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       ptr = match + 4;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (foundTarget) {
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (!isBinaryLike) {  // need to seek to the end of the chunk again to make
    // RDKit✔️❌:                           // sure next read is from the right position
    // RDKit✔️❌:       dp_inStream->clear();
    // RDKit✔️❌:       dp_inStream->seekg(chunkEndPos);
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (bytesRead >= static_cast<std::streamsize>(OVERLAP)) {
    // RDKit✔️❌:       std::memcpy(&buffer[0], bufEnd - OVERLAP, OVERLAP);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     currentStreamPos = chunkEndPos;
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: void SDMolSupplier::buildIndexTo
    // END RDKIT CPP BODY: build_sdf_index

    let mut metadata = Vec::new();
    let mut next_index = 0_usize;
    let mut byte_offset = 0_u64;
    let mut line_number = 0_usize;

    while let Some(raw) =
        extract_next_indexed_sdf_record(reader, next_index, byte_offset, line_number, params)?
    {
        line_number += raw.text.bytes().filter(|byte| *byte == b'\n').count();
        byte_offset += raw.byte_len;
        metadata.push(SdfRecordMetadata {
            index: raw.index,
            byte_offset: raw.byte_offset,
            byte_len: raw.byte_len,
            line_offset: raw.start_line,
            line_len: raw.text.bytes().filter(|byte| *byte == b'\n').count(),
            title: raw
                .text
                .lines()
                .next()
                .map(strip_terminal_cr)
                .filter(|title| !title.is_empty())
                .map(ToOwned::to_owned),
        });
        next_index += 1;
    }

    Ok(metadata)
}

fn read_indexed_sdf_record(
    dataset: &SdfDataset,
    index: usize,
    params: SdfReadParams,
) -> Result<SdfRecord, SdfReadError> {
    // BEGIN RDKIT CPP BODY: read_indexed_sdf_record
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: std::unique_ptr<RWMol> SDMolSupplier::operator[]
    // RDKit✔️❌:
    // RDKit✔️❌:   PRECONDITION(dp_inStream, "no stream");
    // RDKit✔️❌:   // get the molecule with index idx
    // RDKit✔️❌:   moveTo(idx);
    // RDKit✔️❌:   return next();
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: std::unique_ptr<RWMol> SDMolSupplier::operator[]
    // END RDKIT CPP BODY: read_indexed_sdf_record

    // COSMolKit reopens and rescans the file to the requested record instead
    // of seeking to RDKit's cached stream offset, so parameter behavior matches
    // SDMolSupplier while indexed access remains materially slower.
    let file = File::open(&dataset.path).map_err(|err| SdfReadError::Parse(err.to_string()))?;
    let mut reader = BufReader::new(file);
    let mut next_index = 0_usize;
    let mut byte_offset = 0_u64;
    let mut line_number = 0_usize;

    while let Some(raw) =
        extract_next_indexed_sdf_record(&mut reader, next_index, byte_offset, line_number, params)?
    {
        line_number += raw.text.bytes().filter(|byte| *byte == b'\n').count();
        byte_offset += raw.byte_len;
        if raw.index == index {
            return parse_sdf_record_text(&raw.text, raw.start_line, params);
        }
        next_index += 1;
    }

    Err(SdfReadError::Parse(format!(
        "ERROR: Index error (idx = {index}) :  we do not have enough mol blocks"
    )))
}

fn read_next_sdf_record<R: BufRead>(
    reader: &mut SdfReader<R>,
) -> Result<Option<SdfRecord>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: read_next_sdf_record
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: std::unique_ptr<RWMol> ForwardSDMolSupplier::next
    // RDKit✔️❌:
    // RDKit✔️❌:   PRECONDITION(dp_inStream, "no stream");
    // RDKit✔️❌:
    // RDKit✔️❌:   if (dp_inStream->eof()) {
    // RDKit✔️❌:     // FIX: we should probably be throwing an exception here
    // RDKit✔️❌:     df_end = true;
    // RDKit✔️❌:     return nullptr;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   return _next();
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: std::unique_ptr<RWMol> ForwardSDMolSupplier::next
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: std::unique_ptr<RWMol> ForwardSDMolSupplier::_next
    // RDKit✔️❌:
    // RDKit✔️❌:   PRECONDITION(dp_inStream, "no stream");
    // RDKit✔️❌:
    // RDKit✔️❌:   std::string tempStr;
    // RDKit✔️❌:   std::unique_ptr<RWMol> res;
    // RDKit✔️❌:   if (dp_inStream->eof()) {
    // RDKit✔️❌:     df_end = true;
    // RDKit✔️❌:     return res;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   df_eofHitOnRead = false;
    // RDKit✔️❌:   unsigned int line = d_line;
    // RDKit✔️❌:   try {
    // RDKit✔️❌:     MolFromMolDataStream(*dp_inStream, line, d_params).swap(res);
    // RDKit✔️❌:     // there's a special case when trying to read an empty string that
    // RDKit✔️❌:     // we get an empty molecule after only reading a single line without any
    // RDKit✔️❌:     // additional error state.
    // RDKit✔️❌:     if (!res && dp_inStream->eof() && (line - d_line < 2)) {
    // RDKit✔️❌:       df_eofHitOnRead = true;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     d_line = line;
    // RDKit✔️❌:     if (res) {
    // RDKit✔️❌:       this->readMolProps(*res);
    // RDKit✔️❌:     } else if (!dp_inStream->eof()) {
    // RDKit✔️❌:       // FIX: report files missing the $$$$ marker
    // RDKit✔️❌:       std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:       ++d_line;
    // RDKit✔️❌:       while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❌:              (tempStr.at(0) != '$' || tempStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❌:         std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:         ++d_line;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   } catch (FileParseException &fe) {
    // RDKit✔️❌:     if (d_line < static_cast<int>(line)) {
    // RDKit✔️❌:       d_line = line;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // we couldn't read a mol block or the data for the molecule. In this case
    // RDKit✔️❌:     // advance forward in the stream until we hit the next record and then
    // RDKit✔️❌:     // rethrow
    // RDKit✔️❌:     // the exception. This should allow us to read the next molecule.
    // RDKit✔️❌:     BOOST_LOG(rdErrorLog) << "ERROR: " << fe.what() << std::endl;
    // RDKit✔️❌:     BOOST_LOG(rdErrorLog)
    // RDKit✔️❌:         << "ERROR: moving to the beginning of the next molecule\n";
    // RDKit✔️❌:
    // RDKit✔️❌:     // FIX: report files missing the $$$$ marker
    // RDKit✔️❌:     d_line++;
    // RDKit✔️❌:     std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:     while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❌:            (tempStr.empty() || tempStr.at(0) != '$' ||
    // RDKit✔️❌:             tempStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❌:       d_line++;
    // RDKit✔️❌:       std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   } catch (MolSanitizeException &se) {
    // RDKit✔️❌:     if (d_line < static_cast<int>(line)) {
    // RDKit✔️❌:       d_line = line;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     // We couldn't sanitize a molecule we got - write out an error message and
    // RDKit✔️❌:     // move to
    // RDKit✔️❌:     // the beginning of the next molecule
    // RDKit✔️❌:     BOOST_LOG(rdErrorLog)
    // RDKit✔️❌:         << "ERROR: Could not sanitize molecule ending on line " << d_line
    // RDKit✔️❌:         << std::endl;
    // RDKit✔️❌:     BOOST_LOG(rdErrorLog) << "ERROR: " << se.what() << "\n";
    // RDKit✔️❌:
    // RDKit✔️❌:     d_line++;
    // RDKit✔️❌:     std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:     if (dp_inStream->eof()) {
    // RDKit✔️❌:       df_eofHitOnRead = true;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❌:            (tempStr.empty() || tempStr.at(0) != '$' ||
    // RDKit✔️❌:             tempStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❌:       d_line++;
    // RDKit✔️❌:       std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:     }
    // RDKit❌❌:   } catch (...) {
    // RDKit❌❌:     if (dp_inStream->eof()) {
    // RDKit❌❌:       df_eofHitOnRead = true;
    // RDKit❌❌:     }
    // RDKit❌❌:     if (d_line < static_cast<int>(line)) {
    // RDKit❌❌:       d_line = line;
    // RDKit❌❌:     }
    // RDKit❌❌:
    // RDKit❌❌:     BOOST_LOG(rdErrorLog) << "Unexpected error hit on line " << d_line
    // RDKit❌❌:                           << std::endl;
    // RDKit❌❌:     BOOST_LOG(rdErrorLog)
    // RDKit❌❌:         << "ERROR: moving to the beginning of the next molecule\n";
    // RDKit❌❌:     d_line++;
    // RDKit❌❌:     std::getline(*dp_inStream, tempStr);
    // RDKit❌❌:     if (dp_inStream->eof()) {
    // RDKit❌❌:       df_eofHitOnRead = true;
    // RDKit❌❌:     }
    // RDKit❌❌:     while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit❌❌:            (tempStr.empty() || tempStr.at(0) != '$' ||
    // RDKit❌❌:             tempStr.substr(0, 4) != "$$$$")) {
    // RDKit❌❌:       d_line++;
    // RDKit❌❌:       std::getline(*dp_inStream, tempStr);
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit✔️❌:   if (dp_inStream->eof()) {
    // RDKit✔️❌:     // FIX: we should probably be throwing an exception here
    // RDKit✔️❌:     df_end = true;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: std::unique_ptr<RWMol> ForwardSDMolSupplier::_next
    // END RDKIT CPP BODY: read_next_sdf_record

    // COSMolKit buffers the whole record before parsing it, unlike RDKit's
    // stream-positioned MolFromMolDataStream/readMolProps path, so the normal
    // record behavior is reproduced with extra allocation.
    if reader.end {
        return Ok(None);
    }

    reader.eof_hit_on_read = false;
    let Some(raw) = extract_next_forward_sdf_record(
        &mut reader.reader,
        reader.next_index,
        reader.byte_offset,
        reader.line_number,
        reader.params,
    )?
    else {
        reader.end = true;
        reader.eof_hit_on_read = true;
        return Ok(None);
    };

    reader.next_index += 1;
    reader.byte_offset += raw.byte_len;
    reader.line_number += raw.text.bytes().filter(|byte| *byte == b'\n').count();
    if raw.hit_eof {
        reader.end = true;
    }
    parse_forward_sdf_record_text(&raw, reader.params)
}

fn extract_next_forward_sdf_record<R: BufRead>(
    reader: &mut R,
    next_index: usize,
    byte_offset: u64,
    line_number: usize,
    params: SdfReadParams,
) -> Result<Option<RawSdfRecord>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: extract_next_forward_sdf_record
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: std::unique_ptr<RWMol> ForwardSDMolSupplier::_next
    // RDKit✔️❌:
    // RDKit✔️❌:   MolFromMolDataStream(*dp_inStream, line, d_params).swap(res);
    // RDKit✔️❌:   // there's a special case when trying to read an empty string that
    // RDKit✔️❌:   // we get an empty molecule after only reading a single line without any
    // RDKit✔️❌:   // additional error state.
    // RDKit✔️❌:   if (!res && dp_inStream->eof() && (line - d_line < 2)) {
    // RDKit✔️❌:     df_eofHitOnRead = true;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   d_line = line;
    // RDKit✔️❌:   if (res) {
    // RDKit✔️❌:     this->readMolProps(*res);
    // RDKit✔️❌:   } else if (!dp_inStream->eof()) {
    // RDKit✔️❌:     // FIX: report files missing the $$$$ marker
    // RDKit✔️❌:     std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:     ++d_line;
    // RDKit✔️❌:     while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❌:            (tempStr.at(0) != '$' || tempStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❌:       std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:       ++d_line;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: } catch (FileParseException &fe) {
    // RDKit✔️❌:   if (d_line < static_cast<int>(line)) {
    // RDKit✔️❌:     d_line = line;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // we couldn't read a mol block or the data for the molecule. In this case
    // RDKit✔️❌:   // advance forward in the stream until we hit the next record and then
    // RDKit✔️❌:   // rethrow
    // RDKit✔️❌:   // the exception. This should allow us to read the next molecule.
    // RDKit✔️❌:   BOOST_LOG(rdErrorLog) << "ERROR: " << fe.what() << std::endl;
    // RDKit✔️❌:   BOOST_LOG(rdErrorLog)
    // RDKit✔️❌:       << "ERROR: moving to the beginning of the next molecule\n";
    // RDKit✔️❌:
    // RDKit✔️❌:   // FIX: report files missing the $$$$ marker
    // RDKit✔️❌:   d_line++;
    // RDKit✔️❌:   std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:   while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❌:          (tempStr.empty() || tempStr.at(0) != '$' ||
    // RDKit✔️❌:           tempStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❌:     d_line++;
    // RDKit✔️❌:     std::getline(*dp_inStream, tempStr);
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: std::unique_ptr<RWMol> ForwardSDMolSupplier::_next
    // END RDKIT CPP BODY: extract_next_forward_sdf_record

    // ForwardSDMolSupplier recovers parse failures by consuming through the
    // first following "$$$$" line. COSMolKit still buffers that same span into a
    // temporary record before parsing, so behavior is modeled with extra
    // allocation versus RDKit's stream-positioned path.
    let _ = params;
    let mut record = String::new();
    let mut bytes_read = 0_u64;
    let start_line = line_number;
    let mut saw_any_line = false;
    let mut hit_eof = false;

    loop {
        let mut line = String::new();
        let read = reader
            .read_line(&mut line)
            .map_err(|err| SdfReadError::Parse(err.to_string()))?;
        if read == 0 {
            hit_eof = true;
            break;
        }
        saw_any_line = true;
        bytes_read += read as u64;
        if !line.ends_with('\n') {
            hit_eof = true;
        }
        if line.ends_with('\n') {
            line.pop();
        }
        if line.ends_with('\r') {
            line.pop();
        }
        let is_delimiter = is_sdf_record_delimiter(&line);
        record.push_str(&line);
        record.push('\n');
        if is_delimiter {
            break;
        }
    }

    if !saw_any_line || record.find(|c| c != '\n' && c != '\r').is_none() {
        return Ok(None);
    }

    Ok(Some(RawSdfRecord {
        index: next_index,
        byte_offset,
        byte_len: bytes_read,
        start_line,
        hit_eof,
        text: record,
    }))
}

fn extract_next_indexed_sdf_record<R: BufRead>(
    reader: &mut R,
    next_index: usize,
    byte_offset: u64,
    line_number: usize,
    params: SdfReadParams,
) -> Result<Option<RawSdfRecord>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: extract_next_indexed_sdf_record
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: void SDMolSupplier::buildIndexTo
    // RDKit✔️❌:       constexpr char dollarSigns[]{"$$$$"};
    // RDKit✔️❌:       auto match = std::search(ptr, bufEnd, dollarSigns, dollarSigns + 4);
    // RDKit✔️❌:       if (match == bufEnd) {
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (*(match - 1) == '\n') {  // ensure $$$$ is at start of line
    // RDKit✔️❌:         char *nlPos = match + 4;
    // RDKit✔️❌:         if (nlPos == bufEnd) {
    // RDKit✔️❌:           // corner case, $$$$ is EXACTLY at the end of the buffer
    // RDKit✔️❌:           //  we need the next char in the stream to be a "\n", this is resolved
    // RDKit✔️❌:           //  below.
    // RDKit✔️❌:           needEOL = true;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           while (nlPos < bufEnd && *nlPos != '\n') {
    // RDKit✔️❌:             ++nlPos;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           if (nlPos < bufEnd) {
    // RDKit✔️❌:             ++nlPos;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:
    // RDKit✔️❌:         std::streampos posHold;
    // RDKit✔️❌:         if (isBinaryLike &&
    // RDKit✔️❌:             !needEOL) {  // fast path, math checks out, no need to seek
    // RDKit✔️❌:           posHold = chunkStartPos + std::streamoff(nlPos - bufStart - OVERLAP);
    // RDKit✔️❌:         } else {  // slow path, there is byte translation going on, need to seek
    // RDKit✔️❌:                   // and use the std translation magic to find the actual byte
    // RDKit✔️❌:                   // position
    // RDKit✔️❌:           dp_inStream->clear();
    // RDKit✔️❌:           dp_inStream->seekg(
    // RDKit✔️❌:               chunkStartPos);  // rollback to the start of the chunk
    // RDKit✔️❌:           dp_inStream->ignore(
    // RDKit✔️❌:               nlPos - bufStart -
    // RDKit✔️❌:               OVERLAP);  // advance but with the magic translation in effect now
    // RDKit✔️❌:           posHold =
    // RDKit✔️❌:               dp_inStream
    // RDKit✔️❌:                   ->tellg();  // this is the physical position on disk we want
    // RDKit✔️❌:         }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: void SDMolSupplier::buildIndexTo
    // END RDKIT CPP BODY: extract_next_indexed_sdf_record

    // SDMolSupplier indexes the next record by a "$$$$" marker at the start of
    // a line, independent of whether the preceding mol block ended with
    // "M  END". COSMolKit reads line-by-line instead of RDKit's chunk scanner,
    // so behavior matches with extra buffering and worse seek/index cost.
    let _ = params;
    let mut record = String::new();
    let mut bytes_read = 0_u64;
    let start_line = line_number;
    let mut saw_any_line = false;
    let mut hit_eof = false;

    loop {
        let mut line = String::new();
        let read = reader
            .read_line(&mut line)
            .map_err(|err| SdfReadError::Parse(err.to_string()))?;
        if read == 0 {
            hit_eof = true;
            break;
        }
        saw_any_line = true;
        bytes_read += read as u64;
        if !line.ends_with('\n') {
            hit_eof = true;
        }
        if line.ends_with('\n') {
            line.pop();
        }
        if line.ends_with('\r') {
            line.pop();
        }
        let is_delimiter = is_sdf_record_delimiter(&line);
        record.push_str(&line);
        record.push('\n');
        if is_delimiter {
            break;
        }
    }

    if !saw_any_line || record.find(|c| c != '\n' && c != '\r').is_none() {
        return Ok(None);
    }

    Ok(Some(RawSdfRecord {
        index: next_index,
        byte_offset,
        byte_len: bytes_read,
        start_line,
        hit_eof,
        text: record,
    }))
}

fn extract_next_raw_sdf_record<R: BufRead>(
    reader: &mut R,
    next_index: usize,
    byte_offset: u64,
    line_number: usize,
    params: SdfReadParams,
) -> Result<Option<RawSdfRecord>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: extract_next_raw_sdf_record
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MultithreadedSDMolSupplier.cpp :: bool MultithreadedSDMolSupplier::extractNextRecord
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(dp_inStream, "no stream");
    // RDKit✔️❗:   if (dp_inStream->eof()) {
    // RDKit✔️❗:     df_end = true;
    // RDKit✔️❗:     return false;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   std::string currentStr, prevStr;
    // RDKit✔️❗:   record = "";
    // RDKit✔️❗:   lineNum = d_line;
    // RDKit✔️❗:   while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❗:          ((prevStr.find_first_not_of(" \t\r\n") != std::string::npos &&
    // RDKit✔️❗:            prevStr.find("M  END") != 0) ||
    // RDKit✔️❗:           currentStr[0] != '$' || currentStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❗:     prevStr = currentStr;
    // RDKit✔️❗:     std::getline(*dp_inStream, currentStr);
    // RDKit✔️❗:     record += currentStr + "\n";
    // RDKit✔️❗:     ++d_line;
    // RDKit✔️❗:     if (prevStr.find_first_not_of(" \t\r\n") == std::string::npos &&
    // RDKit✔️❗:         currentStr[0] == '$' && currentStr.substr(0, 4) == "$$$$") {
    // RDKit❗✔️:       this->checkForEnd();
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   // ignore trailing new lines
    // RDKit✔️❗:   if (record.find_first_not_of("\n\r") == std::string::npos) {
    // RDKit✔️❗:     return false;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   index = d_currentRecordId;
    // RDKit✔️❗:   ++d_currentRecordId;
    // RDKit✔️❗:   return true;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MultithreadedSDMolSupplier.cpp :: bool MultithreadedSDMolSupplier::extractNextRecord
    // END RDKIT CPP BODY: extract_next_raw_sdf_record

    // RDKit allocates Bond/QueryBond objects here. COSMolKit records the same
    // modeled semantics in a BondSpec and keeps raw molfile codes beside it for
    // later property projection.
    let _ = params;
    let mut record = String::new();
    let mut bytes_read = 0_u64;
    let mut current = String::new();
    let mut previous = String::new();
    let start_line = line_number;
    let mut saw_any_line = false;
    let mut hit_eof = false;

    loop {
        let previous_nonblank = previous.find(|c: char| !" \t\r\n".contains(c)).is_some();
        let previous_is_m_end = previous.starts_with("M  END");
        let current_is_delimiter = current.starts_with("$$$$");
        if (!previous_nonblank || previous_is_m_end) && current_is_delimiter {
            break;
        }

        previous.clone_from(&current);
        current.clear();
        let read = reader
            .read_line(&mut current)
            .map_err(|err| SdfReadError::Parse(err.to_string()))?;
        if read == 0 {
            hit_eof = true;
            break;
        }
        saw_any_line = true;
        bytes_read += read as u64;

        // `std::getline()` removes '\n' but keeps a preceding '\r'. RDKit then
        // appends exactly one '\n' to the accumulated record.
        if !current.ends_with('\n') {
            hit_eof = true;
        }
        if current.ends_with('\n') {
            current.pop();
        }
        record.push_str(&current);
        record.push('\n');

        // RDKit calls checkForEnd() here. COSMolKit's raw-record extraction has
        // no supplier-level end flag yet, so that state update remains ❌❌.
        if previous.find(|c: char| !" \t\r\n".contains(c)).is_none() && current.starts_with("$$$$")
        {
            break;
        }
    }

    if !saw_any_line || record.find(|c| c != '\n' && c != '\r').is_none() {
        return Ok(None);
    }

    Ok(Some(RawSdfRecord {
        index: next_index,
        byte_offset,
        byte_len: bytes_read,
        start_line,
        hit_eof,
        text: record,
    }))
}

fn seek_to_next_sdf_record<R: BufRead>(
    reader: &mut R,
    line_number: usize,
    params: SdfReadParams,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: seek_to_next_sdf_record
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: void SDMolSupplier::moveTo
    // RDKit❗❗:
    // RDKit❗❗:   PRECONDITION(dp_inStream, "no stream");
    // RDKit❗❗:
    // RDKit❗❗:   // dp_inStream->seekg() is called for all idx values
    // RDKit❗❗:   // and earlier calls to next() may have put the stream into a bad state
    // RDKit❗❗:   dp_inStream->clear();
    // RDKit❗❗:
    // RDKit❗❗:   // move until we hit the desired idx
    // RDKit❗❗:   if (idx < d_molpos.size()) {
    // RDKit❗❗:     dp_inStream->seekg(d_molpos[idx]);
    // RDKit❗❗:     d_last = idx;
    // RDKit❗❗:   }
    // RDKit❗❗:   // actually scan with buffering
    // RDKit❗❗:   else {
    // RDKit❗❗:     buildIndexTo(idx);
    // RDKit❗❗:
    // RDKit❗❗:     if (idx < d_molpos.size()) {
    // RDKit❗❗:       dp_inStream->clear();
    // RDKit❗❗:       dp_inStream->seekg(d_molpos[idx]);
    // RDKit❗❗:       d_last = idx;
    // RDKit❗❗:     } else {
    // RDKit❗❗:       /*Unfortunately, the FileParseException is not being catched and thrown on
    // RDKit❗❗:       python directly. Instead, we use this df_end flag workaround to indicate
    // RDKit❗❗:       that we reached the end of file (and signal the error). There's a comment
    // RDKit❗❗:       on MolSupplier.h about problems with Boost exception handling and the full
    // RDKit❗❗:       explanation That's the only reason for the following line*/
    // RDKit❗❗:       df_end = true;
    // RDKit❗❗:
    // RDKit❗❗:       // if we reached end of file without reaching "idx" we have an index error
    // RDKit❗❗:       d_len = rdcast<int>(d_molpos.size());
    // RDKit❗❗:       std::ostringstream errout;
    // RDKit❗❗:       errout << "ERROR: Index error (idx = " << idx << ") : "
    // RDKit❗❗:              << " we do not have enough mol blocks";
    // RDKit❗❗:       throw FileParseException(errout.str());
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp :: void SDMolSupplier::moveTo
    // END RDKIT CPP BODY: seek_to_next_sdf_record

    let _ = params;
    let mut line_number = line_number;
    while let Some(line) = read_rdkit_line(reader)? {
        line_number += 1;
        if is_sdf_record_delimiter(&line) {
            break;
        }
    }
    let _ = line_number;
    Ok(())
}

fn parse_sdf_record_text(
    text: &str,
    start_line: usize,
    params: SdfReadParams,
) -> Result<SdfRecord, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_sdf_record_text
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MultithreadedSDMolSupplier.cpp :: RWMol *MultithreadedSDMolSupplier::processMoleculeRecord
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(dp_inStream, "no stream");
    // RDKit✔️❗:   std::istringstream inStream(record);
    // RDKit✔️❗:   auto res =
    // RDKit✔️❗:       v2::FileParsers::MolFromMolDataStream(inStream, lineNum, d_parseParams);
    // RDKit✔️❗:   if (res) {
    // RDKit✔️❗:     this->readMolProps(*res, inStream);
    // RDKit✔️❗:     return res.release();
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return nullptr;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MultithreadedSDMolSupplier.cpp :: RWMol *MultithreadedSDMolSupplier::processMoleculeRecord
    // END RDKIT CPP BODY: parse_sdf_record_text

    let raw = RawSdfRecord {
        index: 0,
        byte_offset: 0,
        byte_len: text.len() as u64,
        start_line,
        hit_eof: false,
        text: text.to_string(),
    };
    let (mol_block, fields) = split_sdf_record(&raw, params)?;
    let parsed = mol_from_mol_block(mol_block, params)?;
    let mut molecule = parsed.molecule;
    for field in &fields {
        molecule = molecule
            .with_prop(field.name.clone(), field.value.clone())
            .with_sdf_data_field(field.name.clone(), field.value.clone());
    }
    process_sdf_property_lists(&mut molecule, &fields, params)?;
    Ok(SdfRecord {
        molecule,
        data_fields: fields
            .into_iter()
            .map(|field| (field.name, field.value))
            .collect(),
    })
}

fn parse_forward_sdf_record_text(
    raw: &RawSdfRecord,
    params: SdfReadParams,
) -> Result<Option<SdfRecord>, SdfReadError> {
    let (mol_block, data_lines, data_start_line) = split_sdf_record_parts(raw);
    let parsed = match mol_from_mol_block(mol_block, params) {
        Ok(parsed) => parsed,
        Err(SdfReadError::Parse(_))
        | Err(SdfReadError::Operation(crate::OperationError::Sanitize { .. })) => {
            return Ok(None);
        }
        Err(error) => return Err(error),
    };

    let mut molecule = parsed.molecule;
    let fields = match parse_sdf_data_fields(&data_lines, data_start_line, params) {
        Ok(fields) => fields,
        Err(SdfReadError::Parse(_)) => Vec::new(),
        Err(error) => return Err(error),
    };
    for field in &fields {
        molecule = molecule
            .with_prop(field.name.clone(), field.value.clone())
            .with_sdf_data_field(field.name.clone(), field.value.clone());
    }
    process_sdf_property_lists(&mut molecule, &fields, params)?;
    Ok(Some(SdfRecord {
        molecule,
        data_fields: fields
            .into_iter()
            .map(|field| (field.name, field.value))
            .collect(),
    }))
}

fn split_sdf_record(
    raw: &RawSdfRecord,
    params: SdfReadParams,
) -> Result<(&str, Vec<SdfDataField>), SdfReadError> {
    let (mol_block, data_lines, data_start_line) = split_sdf_record_parts(raw);
    let fields = parse_sdf_data_fields(&data_lines, data_start_line, params)?;
    Ok((mol_block, fields))
}

fn split_sdf_record_parts(raw: &RawSdfRecord) -> (&str, Vec<&str>, usize) {
    let mut data_start = raw.text.len();
    let mut data_start_line = raw.start_line;
    let mut offset = 0;
    for (line_idx, line) in raw.text.split_inclusive('\n').enumerate() {
        let content = line.trim_end_matches('\n');
        let content = content.strip_suffix('\r').unwrap_or(content);
        let stripped = strip_sdf_line(content);
        if stripped.starts_with('>') || is_sdf_record_delimiter(stripped) {
            data_start = offset;
            data_start_line = raw.start_line + line_idx;
            break;
        }
        offset += line.len();
    }

    let mol_block = &raw.text[..data_start];
    let data_text = &raw.text[data_start..];
    let data_lines: Vec<&str> = data_text.lines().collect();
    (mol_block, data_lines, data_start_line)
}

fn parse_sdf_data_fields(
    lines: &[&str],
    start_line: usize,
    params: SdfReadParams,
) -> Result<Vec<SdfDataField>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_sdf_data_fields
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: void ForwardSDMolSupplier::readMolProps
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(dp_inStream, "no stream");
    // RDKit✔️❗:   d_line++;
    // RDKit✔️❗:   bool hasProp = false;
    // RDKit✔️❗:   bool warningIssued = false;
    // RDKit✔️❗:   std::string dlabel = "";
    // RDKit✔️❗:   std::string inl;
    // RDKit✔️❗:   std::getline(*dp_inStream, inl);
    // RDKit✔️❗:   std::string_view tempStr = inl;
    // RDKit✔️❗:
    // RDKit✔️❗:   // FIX: report files missing the $$$$ marker
    // RDKit✔️❗:   while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❗:          (tempStr.empty() || tempStr.at(0) != '$' ||
    // RDKit✔️❗:           tempStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❗:     tempStr = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:     if (!tempStr.empty()) {
    // RDKit✔️❗:       if (tempStr.at(0) == '>') {  // data header line: start of a data item
    // RDKit✔️❗:         // ignore all other crap and seek for for a data label enclosed
    // RDKit✔️❗:         // by '<' and '>'
    // RDKit✔️❗:         // FIX: "CTfile.pdf" (page 51) says that the data header line does not
    // RDKit✔️❗:         // have to contain a data label (instead can have something line field
    // RDKit✔️❗:         // id into a MACCS db). But we do not currently know what to do in this
    // RDKit✔️❗:         // situation - so ignore such data items for now
    // RDKit✔️❗:         hasProp = true;
    // RDKit✔️❗:         warningIssued = false;
    // RDKit✔️❗:         tempStr = tempStr.substr(1);     // remove the first ">" sign
    // RDKit✔️❗:         size_t sl = tempStr.find("<");   // begin datalabel
    // RDKit✔️❗:         size_t se = tempStr.rfind(">");  // end datalabel
    // RDKit✔️❗:         if ((sl == std::string::npos) || (se == std::string::npos) ||
    // RDKit✔️❗:             (se == (sl + 1))) {
    // RDKit✔️❗:           // we either do not have a data label or the label is empty
    // RDKit✔️❗:           // no data label ignore until next data item
    // RDKit✔️❗:           // i.e. until we hit a blank line
    // RDKit✔️❗:           d_line++;
    // RDKit✔️❗:           std::getline(*dp_inStream, inl);
    // RDKit✔️❗:           tempStr = inl;
    // RDKit✔️❗:           auto stmp = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:           while (stmp.length() != 0) {
    // RDKit✔️❗:             d_line++;
    // RDKit✔️❗:             std::getline(*dp_inStream, inl);
    // RDKit✔️❗:             tempStr = inl;
    // RDKit✔️❗:             if (dp_inStream->eof()) {
    // RDKit✔️❗:               throw FileParseException("End of data field name not found");
    // RDKit✔️❗:             }
    // RDKit✔️❗:           }
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           dlabel = tempStr.substr(sl + 1, se - sl - 1);
    // RDKit✔️❗:           // we know the label - now read in the relevant properties
    // RDKit✔️❗:           // until we hit a blank line
    // RDKit✔️❗:           d_line++;
    // RDKit✔️❗:           std::getline(*dp_inStream, inl);
    // RDKit✔️❗:           tempStr = inl;
    // RDKit✔️❗:
    // RDKit✔️❗:           std::string prop = "";
    // RDKit✔️❗:           auto stmp = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:           int nplines = 0;  // number of lines for this property
    // RDKit✔️❗:           while (!stmp.empty() ||
    // RDKit✔️❗:                  (!tempStr.empty() &&
    // RDKit✔️❗:                   (tempStr.at(0) == ' ' || tempStr.at(0) == '\t'))) {
    // RDKit✔️❗:             nplines++;
    // RDKit✔️❗:             if (nplines > 1) {
    // RDKit✔️❗:               prop += "\n";
    // RDKit✔️❗:             }
    // RDKit✔️❗:             // take off \r if it's still in the property:
    // RDKit✔️❗:             if (!tempStr.empty()) {
    // RDKit✔️❗:               if (tempStr.back() == '\r') {
    // RDKit✔️❗:                 tempStr = tempStr.substr(0, tempStr.size() - 1);
    // RDKit✔️❗:               }
    // RDKit✔️❗:               prop += tempStr;
    // RDKit✔️❗:             }
    // RDKit✔️❗:             d_line++;
    // RDKit✔️❗:             // erase inl in case the file does not end with a carriage
    // RDKit✔️❗:             // return (we will end up in an infinite loop if we don't do
    // RDKit✔️❗:             // this and we do not check for EOF in this while loop body)
    // RDKit✔️❗:             inl.erase();
    // RDKit✔️❗:             std::getline(*dp_inStream, inl);
    // RDKit✔️❗:             tempStr = inl;
    // RDKit✔️❗:             if (tempStr.empty()) {
    // RDKit✔️❗:               stmp = tempStr;
    // RDKit✔️❗:             } else {
    // RDKit✔️❗:               stmp = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:             }
    // RDKit✔️❗:           }
    // RDKit✔️❗:           mol.setProp(dlabel, prop);
    // RDKit✔️❗:           if (df_processPropertyLists) {
    // RDKit✔️❗:             // apply this as an atom property list if that's appropriate
    // RDKit✔️❗:             FileParserUtils::processMolPropertyList(mol, dlabel);
    // RDKit✔️❗:           }
    // RDKit✔️❗:         }
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         if (d_params.strictParsing) {
    // RDKit✔️❗:           // at this point we should always be at a line starting with '>'
    // RDKit✔️❗:           // following a blank line. If this is not true and df_strictParsing
    // RDKit✔️❗:           // is true, then throw an exception, otherwise truncate the rest of
    // RDKit✔️❗:           // the data field following the blank line until the next '>' or EOF
    // RDKit✔️❗:           // and issue a warning
    // RDKit✔️❗:           // FIX: should we be deleting the molecule (which is probably fine)
    // RDKit✔️❗:           // because we couldn't read the data ???
    // RDKit✔️❗:           throw FileParseException("Problems encountered parsing data fields");
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           if (!warningIssued) {
    // RDKit❗✔️:             if (hasProp) {
    // RDKit❗✔️:               BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                   << "Property <" << dlabel << "> will be truncated after "
    // RDKit❗✔️:                   << "the first blank line" << std::endl;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:                   << "Spurious data before the first property will be "
    // RDKit❗✔️:                      "ignored"
    // RDKit❗✔️:                   << std::endl;
    // RDKit❗✔️:             }
    // RDKit✔️❗:             warningIssued = true;
    // RDKit✔️❗:           }
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:     d_line++;
    // RDKit✔️❗:     std::getline(*dp_inStream, inl);
    // RDKit✔️❗:     tempStr = inl;
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: void ForwardSDMolSupplier::readMolProps
    // END RDKIT CPP BODY: parse_sdf_data_fields

    // RDKit consumes an owning stream and mutates ROMol properties in-place.
    // COSMolKit currently parses an already split record and returns data
    // fields for later molecule-property projection, so stream movement and
    // `mol.setProp()` are marked ✔️❗ instead of ✔️✔️.
    let mut fields = Vec::new();
    let mut idx = 0;
    let mut has_prop = false;
    let mut warning_issued = false;
    let mut current_label = String::new();

    while idx < lines.len() {
        let mut temp_str = strip_sdf_line(lines[idx]);
        if is_sdf_record_delimiter(temp_str) {
            break;
        }

        if !temp_str.is_empty() {
            if temp_str.starts_with('>') {
                has_prop = true;
                warning_issued = false;
                match parse_sdf_data_header(temp_str, start_line + idx)? {
                    Some(label) => {
                        current_label = label;
                        idx += 1;
                        let mut prop = String::new();
                        let mut property_line_count = 0;
                        while idx < lines.len() {
                            let line = lines[idx];
                            if is_sdf_record_delimiter(strip_sdf_line(line)) {
                                break;
                            }
                            let stripped = strip_sdf_line(line);
                            if stripped.is_empty() && !starts_with_sdf_continuation_space(line) {
                                break;
                            }
                            property_line_count += 1;
                            if property_line_count > 1 {
                                prop.push('\n');
                            }
                            prop.push_str(strip_terminal_cr(line));
                            idx += 1;
                        }
                        fields.push(SdfDataField {
                            name: current_label.clone(),
                            value: prop,
                        });
                    }
                    None => {
                        idx += 1;
                        while idx < lines.len() && !strip_sdf_line(lines[idx]).is_empty() {
                            if is_sdf_record_delimiter(strip_sdf_line(lines[idx])) {
                                return Err(SdfReadError::Parse(
                                    "End of data field name not found".to_string(),
                                ));
                            }
                            idx += 1;
                        }
                        if idx >= lines.len() {
                            return Err(SdfReadError::Parse(
                                "End of data field name not found".to_string(),
                            ));
                        }
                    }
                }
            } else if params.strict_parsing {
                return Err(SdfReadError::Parse(
                    "Problems encountered parsing data fields".to_string(),
                ));
            } else if !warning_issued {
                warning_issued = true;
            }
        }

        idx += 1;
        if idx < lines.len() {
            temp_str = strip_sdf_line(lines[idx]);
            let _ = temp_str;
        }
    }

    let _ = (has_prop, warning_issued, current_label);
    Ok(fields)
}

fn parse_sdf_data_header(line: &str, line_number: usize) -> Result<Option<String>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_sdf_data_header
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: void ForwardSDMolSupplier::readMolProps header extraction
    // RDKit✔️❗:     tempStr = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:     if (!tempStr.empty()) {
    // RDKit✔️❗:       if (tempStr.at(0) == '>') {  // data header line: start of a data item
    // RDKit✔️❗:         // ignore all other crap and seek for for a data label enclosed
    // RDKit✔️❗:         // by '<' and '>'
    // RDKit✔️❗:         // FIX: "CTfile.pdf" (page 51) says that the data header line does not
    // RDKit✔️❗:         // have to contain a data label (instead can have something line field
    // RDKit✔️❗:         // id into a MACCS db). But we do not currently know what to do in this
    // RDKit✔️❗:         // situation - so ignore such data items for now
    // RDKit✔️❗:         tempStr = tempStr.substr(1);     // remove the first ">" sign
    // RDKit✔️❗:         size_t sl = tempStr.find("<");   // begin datalabel
    // RDKit✔️❗:         size_t se = tempStr.rfind(">");  // end datalabel
    // RDKit✔️❗:         if ((sl == std::string::npos) || (se == std::string::npos) ||
    // RDKit✔️❗:             (se == (sl + 1))) {
    // RDKit✔️❗:           // we either do not have a data label or the label is empty
    // RDKit✔️❗:           // no data label ignore until next data item
    // RDKit✔️❗:           // i.e. until we hit a blank line
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           dlabel = tempStr.substr(sl + 1, se - sl - 1);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: void ForwardSDMolSupplier::readMolProps header extraction
    // END RDKIT CPP BODY: parse_sdf_data_header

    let _ = line_number;
    let temp_str = line.trim_matches([' ', '\t', '\r', '\n']);
    if temp_str.is_empty() {
        return Ok(None);
    }
    if !temp_str.starts_with('>') {
        return Ok(None);
    }

    let temp_str = &temp_str[1..];
    let Some(start) = temp_str.find('<') else {
        return Ok(None);
    };
    let Some(end) = temp_str.rfind('>') else {
        return Ok(None);
    };
    if end == start + 1 {
        return Ok(None);
    }
    if end < start {
        // C++ uses unsigned size_t arithmetic here. When the last '>' is before
        // '<', `se - sl - 1` underflows and `substr()` returns from `sl + 1`
        // through the end. Rust must spell that behavior out to avoid panic.
        return Ok(Some(temp_str[start + 1..].to_string()));
    }

    Ok(Some(temp_str[start + 1..end].to_string()))
}

fn process_sdf_property_lists(
    molecule: &mut Molecule,
    fields: &[SdfDataField],
    params: SdfReadParams,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_sdf_property_lists
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: void ForwardSDMolSupplier::readMolProps
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(dp_inStream, "no stream");
    // RDKit✔️❗:   d_line++;
    // RDKit✔️❗:   bool hasProp = false;
    // RDKit✔️❗:   bool warningIssued = false;
    // RDKit✔️❗:   std::string dlabel = "";
    // RDKit✔️❗:   std::string inl;
    // RDKit✔️❗:   std::getline(*dp_inStream, inl);
    // RDKit✔️❗:   std::string_view tempStr = inl;
    // RDKit✔️❗:
    // RDKit✔️❗:   // FIX: report files missing the $$$$ marker
    // RDKit✔️❗:   while (!dp_inStream->eof() && !dp_inStream->fail() &&
    // RDKit✔️❗:          (tempStr.empty() || tempStr.at(0) != '$' ||
    // RDKit✔️❗:           tempStr.substr(0, 4) != "$$$$")) {
    // RDKit✔️❗:     tempStr = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:     if (!tempStr.empty()) {
    // RDKit✔️❗:       if (tempStr.at(0) == '>') {  // data header line: start of a data item
    // RDKit✔️❗:         // ignore all other crap and seek for for a data label enclosed
    // RDKit✔️❗:         // by '<' and '>'
    // RDKit✔️❗:         // FIX: "CTfile.pdf" (page 51) says that the data header line does not
    // RDKit✔️❗:         // have to contain a data label (instead can have something line field
    // RDKit✔️❗:         // id into a MACCS db). But we do not currently know what to do in this
    // RDKit✔️❗:         // situation - so ignore such data items for now
    // RDKit✔️❗:         hasProp = true;
    // RDKit✔️❗:         warningIssued = false;
    // RDKit✔️❗:         tempStr = tempStr.substr(1);     // remove the first ">" sign
    // RDKit✔️❗:         size_t sl = tempStr.find("<");   // begin datalabel
    // RDKit✔️❗:         size_t se = tempStr.rfind(">");  // end datalabel
    // RDKit✔️❗:         if ((sl == std::string::npos) || (se == std::string::npos) ||
    // RDKit✔️❗:             (se == (sl + 1))) {
    // RDKit✔️❗:           // we either do not have a data label or the label is empty
    // RDKit✔️❗:           // no data label ignore until next data item
    // RDKit✔️❗:           // i.e. until we hit a blank line
    // RDKit✔️❗:           d_line++;
    // RDKit✔️❗:           std::getline(*dp_inStream, inl);
    // RDKit✔️❗:           tempStr = inl;
    // RDKit✔️❗:           auto stmp = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:           while (stmp.length() != 0) {
    // RDKit✔️❗:             d_line++;
    // RDKit✔️❗:             std::getline(*dp_inStream, inl);
    // RDKit✔️❗:             tempStr = inl;
    // RDKit✔️❗:             if (dp_inStream->eof()) {
    // RDKit✔️❗:               throw FileParseException("End of data field name not found");
    // RDKit✔️❗:             }
    // RDKit✔️❗:           }
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           dlabel = tempStr.substr(sl + 1, se - sl - 1);
    // RDKit✔️❗:           // we know the label - now read in the relevant properties
    // RDKit✔️❗:           // until we hit a blank line
    // RDKit✔️❗:           d_line++;
    // RDKit✔️❗:           std::getline(*dp_inStream, inl);
    // RDKit✔️❗:           tempStr = inl;
    // RDKit✔️❗:
    // RDKit✔️❗:           std::string prop = "";
    // RDKit✔️❗:           auto stmp = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:           int nplines = 0;  // number of lines for this property
    // RDKit✔️❗:           while (!stmp.empty() ||
    // RDKit✔️❗:                  (!tempStr.empty() &&
    // RDKit✔️❗:                   (tempStr.at(0) == ' ' || tempStr.at(0) == '\t'))) {
    // RDKit✔️❗:             nplines++;
    // RDKit✔️❗:             if (nplines > 1) {
    // RDKit✔️❗:               prop += "\n";
    // RDKit✔️❗:             }
    // RDKit✔️❗:             // take off \r if it's still in the property:
    // RDKit✔️❗:             if (!tempStr.empty()) {
    // RDKit✔️❗:               if (tempStr.back() == '\r') {
    // RDKit✔️❗:                 tempStr = tempStr.substr(0, tempStr.size() - 1);
    // RDKit✔️❗:               }
    // RDKit✔️❗:               prop += tempStr;
    // RDKit✔️❗:             }
    // RDKit✔️❗:             d_line++;
    // RDKit✔️❗:             // erase inl in case the file does not end with a carriage
    // RDKit✔️❗:             // return (we will end up in an infinite loop if we don't do
    // RDKit✔️❗:             // this and we do not check for EOF in this while loop body)
    // RDKit✔️❗:             inl.erase();
    // RDKit✔️❗:             std::getline(*dp_inStream, inl);
    // RDKit✔️❗:             tempStr = inl;
    // RDKit✔️❗:             if (tempStr.empty()) {
    // RDKit✔️❗:               stmp = tempStr;
    // RDKit✔️❗:             } else {
    // RDKit✔️❗:               stmp = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:             }
    // RDKit✔️❗:           }
    // RDKit✔️❗:           mol.setProp(dlabel, prop);
    // RDKit✔️❗:           if (df_processPropertyLists) {
    // RDKit✔️❗:             // apply this as an atom property list if that's appropriate
    // RDKit✔️❗:             FileParserUtils::processMolPropertyList(mol, dlabel);
    // RDKit✔️❗:           }
    // RDKit✔️❗:         }
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         if (d_params.strictParsing) {
    // RDKit✔️❗:           // at this point we should always be at a line starting with '>'
    // RDKit✔️❗:           // following a blank line. If this is not true and df_strictParsing
    // RDKit✔️❗:           // is true, then throw an exception, otherwise truncate the rest of
    // RDKit✔️❗:           // the data field following the blank line until the next '>' or EOF
    // RDKit✔️❗:           // and issue a warning
    // RDKit✔️❗:           // FIX: should we be deleting the molecule (which is probably fine)
    // RDKit✔️❗:           // because we couldn't read the data ???
    // RDKit✔️❗:           throw FileParseException("Problems encountered parsing data fields");
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           if (!warningIssued) {
    // RDKit✔️❗:             if (hasProp) {
    // RDKit✔️❗:               BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:                   << "Property <" << dlabel << "> will be truncated after "
    // RDKit✔️❗:                   << "the first blank line" << std::endl;
    // RDKit✔️❗:             } else {
    // RDKit✔️❗:               BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:                   << "Spurious data before the first property will be "
    // RDKit✔️❗:                      "ignored"
    // RDKit✔️❗:                   << std::endl;
    // RDKit✔️❗:             }
    // RDKit✔️❗:             warningIssued = true;
    // RDKit✔️❗:           }
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:     d_line++;
    // RDKit✔️❗:     std::getline(*dp_inStream, inl);
    // RDKit✔️❗:     tempStr = inl;
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp :: void ForwardSDMolSupplier::readMolProps
    // END RDKIT CPP BODY: process_sdf_property_lists

    if params.process_property_lists {
        for field in fields {
            process_sdf_property_list(molecule, field);
        }
    }
    Ok(())
}

fn process_sdf_property_list(molecule: &mut Molecule, field: &SdfDataField) {
    // BEGIN RDKIT CPP BODY: process_sdf_property_list
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/FileParserUtils.h :: inline void processMolPropertyList
    // RDKit✔️❗: inline void processMolPropertyList(
    // RDKit✔️❗:     ROMol &mol, const std::string &pn,
    // RDKit✔️❗:     const std::string &missingValueMarker = "n/a") {
    // RDKit✔️❗:   auto propSetter = [&](const std::string &propPrefix, auto getter,
    // RDKit✔️❗:                         size_t nItems) {
    // RDKit✔️❗:     std::string prefix = propPrefix + "prop.";
    // RDKit✔️❗:     if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
    // RDKit✔️❗:       applyMolListProp<std::string>(mol, pn, prefix, missingValueMarker, nItems,
    // RDKit✔️❗:                                     getter);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       prefix = propPrefix + "iprop.";
    // RDKit✔️❗:       if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
    // RDKit✔️❗:         applyMolListProp<int>(mol, pn, prefix, missingValueMarker, nItems,
    // RDKit✔️❗:                               getter);
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         prefix = propPrefix + "dprop.";
    // RDKit✔️❗:         if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
    // RDKit✔️❗:           applyMolListProp<double>(mol, pn, prefix, missingValueMarker, nItems,
    // RDKit✔️❗:                                    getter);
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           prefix = propPrefix + "bprop.";
    // RDKit✔️❗:           if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
    // RDKit✔️❗:             applyMolListProp<bool>(mol, pn, prefix, missingValueMarker, nItems,
    // RDKit✔️❗:                                    getter);
    // RDKit✔️❗:           }
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:   };
    // RDKit✔️❗:
    // RDKit✔️❗:   if (pn.find(atomPropPrefix) == 0 && pn.length() > atomPropPrefixLength) {
    // RDKit✔️❗:     propSetter(
    // RDKit✔️❗:         atomPropPrefix,
    // RDKit✔️❗:         [&mol](size_t which) { return mol.getAtomWithIdx(which); },
    // RDKit✔️❗:         mol.getNumAtoms());
    // RDKit✔️❗:   } else if (pn.find(bondPropPrefix) == 0 &&
    // RDKit✔️❗:              pn.length() > bondPropPrefixLength) {
    // RDKit✔️❗:     propSetter(
    // RDKit✔️❗:         bondPropPrefix,
    // RDKit✔️❗:         [&mol](size_t which) { return mol.getBondWithIdx(which); },
    // RDKit✔️❗:         mol.getNumBonds());
    // RDKit✔️❗:   }
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/FileParserUtils.h :: inline void processMolPropertyList
    // END RDKIT CPP BODY: process_sdf_property_list

    // COSMolKit stores atom/bond properties as strings, so successful int,
    // double, and bool lexical casts are preserved as string payloads instead
    // of typed RDProps values. Prefix matching and missing-value behavior
    // follow RDKit's processMolPropertyList/applyMolListProp path.
    let Some((target, prop_name, value_kind)) = sdf_property_list_target(&field.name) else {
        return;
    };
    let n_items = match target {
        SdfPropertyListTarget::Atom => molecule.num_atoms(),
        SdfPropertyListTarget::Bond => molecule.num_bonds(),
    };
    let Some(values) = parse_sdf_property_list_values(&field.value, n_items, value_kind) else {
        return;
    };

    match target {
        SdfPropertyListTarget::Atom => {
            let topology = molecule.topology_block_mut();
            for (idx, value) in values.iter().enumerate() {
                if let Some(value) = value {
                    topology.atoms[idx].set_prop(prop_name, value);
                }
            }
        }
        SdfPropertyListTarget::Bond => {
            let topology = molecule.topology_block_mut();
            for (idx, value) in values.iter().enumerate() {
                if let Some(value) = value {
                    topology.bonds[idx].set_prop(prop_name, value);
                }
            }
        }
    }

    let properties = molecule.properties_mut();
    *properties = properties
        .clone()
        .with_sdf_property_list(SdfPropertyList::new(target, prop_name, values));
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
    PREFIXES.iter().find_map(|(prefix, target, kind)| {
        name.strip_prefix(prefix)
            .filter(|prop_name| !prop_name.is_empty())
            .map(|prop_name| (*target, prop_name, *kind))
    })
}

fn parse_sdf_property_list_values(
    value: &str,
    n_items: usize,
    value_kind: SdfPropertyListValueKind,
) -> Option<Vec<Option<String>>> {
    // BEGIN RDKIT CPP BODY: parse_sdf_property_list_values
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/FileParserUtils.h :: void applyMolListProp
    // RDKit✔️❗: void applyMolListProp(ROMol &mol, const std::string &pn,
    // RDKit✔️❗:                       const std::string &prefix,
    // RDKit✔️❗:                       const std::string &missingValueMarker, size_t nItems,
    // RDKit✔️❗:                       U getter) {
    // RDKit✔️❗:   std::string itempn = pn.substr(prefix.size());
    // RDKit✔️❗:   std::string strVect = mol.getProp<std::string>(pn);
    // RDKit✔️❗:   std::vector<std::string> tokens;
    // RDKit✔️❗:   boost::split(tokens, strVect, boost::is_any_of(" \t\n"),
    // RDKit✔️❗:                boost::token_compress_on);
    // RDKit✔️❗:   std::string mv = missingValueMarker;
    // RDKit✔️❗:   size_t first_token = 0;
    // RDKit✔️❗:   if (tokens.size() == nItems + 1 && tokens[0].front() == '[' &&
    // RDKit✔️❗:       tokens[0].back() == ']') {
    // RDKit✔️❗:     mv = std::string(tokens[0].begin() + 1, tokens[0].end() - 1);
    // RDKit✔️❗:     first_token = 1;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (mv.empty()) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "Missing value marker for property " << pn
    // RDKit❗✔️:                             << " is empty." << std::endl;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if(tokens.size() - first_token != nItems) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "Property list " << pn << " has incompatible size, "
    // RDKit❗✔️:                             << tokens.size() << " elements found; expecting "
    // RDKit❗✔️:                             << nItems << ". Ignoring it." << std::endl;
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   for (size_t i = first_token; i < tokens.size(); ++i) {
    // RDKit✔️❗:     if (tokens[i] != mv) {
    // RDKit✔️❗:       unsigned int itemid = i - first_token;
    // RDKit✔️❗:       try {
    // RDKit✔️❗:         T apv = boost::lexical_cast<T>(tokens[i]);
    // RDKit✔️❗:         getter(itemid)->setProp(itempn, apv);
    // RDKit✔️❗:       } catch (const boost::bad_lexical_cast &) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "Value " << tokens[i] << " for property " << pn << " of item "
    // RDKit❗✔️:             << itemid << " can not be parsed. Ignoring it." << std::endl;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/FileParserUtils.h :: void applyMolListProp
    // END RDKIT CPP BODY: parse_sdf_property_list_values

    // RDKit keeps the original molecule-level property and applies parsed
    // item-level values. COSMolKit returns a parallel vector so callers can
    // set string-backed atom/bond properties and record SdfPropertyList state.
    let tokens = split_rdkit_property_list_tokens(value);
    let mut missing_value = "n/a";
    let mut first_token = 0;
    if tokens.len() == n_items + 1
        && tokens[0].starts_with('[')
        && tokens[0].ends_with(']')
        && tokens[0].len() >= 2
    {
        missing_value = &tokens[0][1..tokens[0].len() - 1];
        first_token = 1;
    }
    if tokens.len().saturating_sub(first_token) != n_items {
        return None;
    }
    let mut values = Vec::with_capacity(n_items);
    for token in &tokens[first_token..] {
        if *token == missing_value {
            values.push(None);
            continue;
        }
        let parsed = match value_kind {
            SdfPropertyListValueKind::String => Some((*token).to_string()),
            SdfPropertyListValueKind::Int => {
                token.parse::<i64>().ok().map(|_| (*token).to_string())
            }
            SdfPropertyListValueKind::Double => {
                token.parse::<f64>().ok().map(|_| (*token).to_string())
            }
            SdfPropertyListValueKind::Bool => parse_rdkit_bool_token(token),
        };
        values.push(parsed);
    }
    Some(values)
}

fn split_rdkit_property_list_tokens(value: &str) -> Vec<&str> {
    let bytes = value.as_bytes();
    let mut tokens = Vec::new();
    let mut start = 0;
    let mut idx = 0;
    while idx < bytes.len() {
        if matches!(bytes[idx], b' ' | b'\t' | b'\n') {
            tokens.push(&value[start..idx]);
            while idx < bytes.len() && matches!(bytes[idx], b' ' | b'\t' | b'\n') {
                idx += 1;
            }
            start = idx;
        } else {
            idx += 1;
        }
    }
    tokens.push(&value[start..]);
    tokens
}

fn parse_rdkit_bool_token(token: &str) -> Option<String> {
    match token {
        "1" => Some("true".to_string()),
        "0" => Some("false".to_string()),
        _ => None,
    }
}

fn mol_from_mol_block(block: &str, params: SdfReadParams) -> Result<ParsedMolBlock, SdfReadError> {
    // BEGIN RDKIT CPP BODY: mol_from_mol_block
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::unique_ptr<RWMol> MolFromMolBlock
    // RDKit✔️✔️: std::unique_ptr<RWMol> MolFromMolBlock(const std::string &molBlock,
    // RDKit✔️✔️:                                        const MolFileParserParams &params) {
    // RDKit✔️✔️:   std::istringstream inStream(molBlock);
    // RDKit✔️✔️:   unsigned int line = 0;
    // RDKit✔️✔️:   return MolFromMolDataStream(inStream, line, params);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::unique_ptr<RWMol> MolFromMolBlock
    // END RDKIT CPP BODY: mol_from_mol_block

    let mut reader = std::io::Cursor::new(block.as_bytes());
    let mut line_number = 0;
    mol_from_mol_data_stream(&mut reader, &mut line_number, params)
}

pub(crate) fn read_mol_block_molecule_with_params(
    block: &str,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    Ok(mol_from_mol_block(block, params)?.molecule)
}

pub(crate) fn read_mol_data_stream_molecule_with_params<R: BufRead>(
    reader: &mut R,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    Ok(mol_from_mol_data_stream(reader, line_number, params)?.molecule)
}

fn mol_from_mol_data_stream<R: BufRead>(
    reader: &mut R,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<ParsedMolBlock, SdfReadError> {
    // BEGIN RDKIT CPP BODY: mol_from_mol_data_stream
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::unique_ptr<RWMol> MolFromMolDataStream
    // RDKit✔️❗:
    // RDKit✔️❗:   std::string tempStr;
    // RDKit✔️❗:   bool fileComplete = false;
    // RDKit✔️❗:   bool chiralityPossible = false;
    // RDKit✔️❗:   Utils::LocaleSwitcher ls;
    // RDKit✔️❗:   // mol name
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:   if (inStream.eof()) {
    // RDKit✔️❗:     return nullptr;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   auto res = std::make_unique<RWMol>();
    // RDKit✔️❗:   res->setProp(common_properties::_Name, tempStr);
    // RDKit✔️❗:
    // RDKit✔️❗:   // info
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:   res->setProp("_MolFileInfo", tempStr);
    // RDKit✔️❗:   if (tempStr.length() >= 22) {
    // RDKit✔️❗:     std::string dimLabel = tempStr.substr(20, 2);
    // RDKit✔️❗:     // Unless labelled as 3D we assume 2D
    // RDKit✔️❗:     if (dimLabel == "3d" || dimLabel == "3D") {
    // RDKit✔️❗:       res->setProp(common_properties::_3DConf, 1);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   // comments
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:   res->setProp("_MolFileComments", tempStr);
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nAtoms = 0, nBonds = 0, nLists = 0, chiralFlag = 0, nsText = 0,
    // RDKit✔️❗:                nRxnComponents = 0;
    // RDKit✔️❗:   int nReactants = 0, nProducts = 0, nIntermediates = 0;
    // RDKit✔️❗:   (void)nLists;  // read from the file but unused
    // RDKit✔️❗:   (void)nsText;
    // RDKit✔️❗:   (void)nRxnComponents;
    // RDKit✔️❗:   (void)nReactants;
    // RDKit✔️❗:   (void)nProducts;
    // RDKit✔️❗:   (void)nIntermediates;
    // RDKit✔️❗:   // counts line, this is where we really get started
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:
    // RDKit✔️❗:   if (tempStr.size() < 6) {
    // RDKit✔️❗:     if (res) {
    // RDKit✔️❗:       res = nullptr;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Counts line too short: '" << tempStr << "' on line" << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int spos = 0;
    // RDKit✔️❗:   // this needs to go into a try block because if the lexical_cast throws an
    // RDKit✔️❗:   // exception we want to catch throw a different exception
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nAtoms = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     spos = 3;
    // RDKit✔️❗:     nBonds = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     spos = 6;
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << tempStr.substr(spos, 3)
    // RDKit✔️❗:            << "' to unsigned int on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     spos = 6;
    // RDKit✔️❗:     if (tempStr.size() >= 9) {
    // RDKit✔️❗:       nLists = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 12;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       chiralFlag = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 15;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nsText = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 18;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nRxnComponents =
    // RDKit✔️❗:           FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 21;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nReactants = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 24;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nProducts = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 27;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nIntermediates =
    // RDKit✔️❗:           FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     // some SD files (such as some from NCI) lack all the extra information
    // RDKit✔️❗:     // on the header line, so ignore problems parsing there.
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int ctabVersion = 2000;
    // RDKit✔️❗:   if (tempStr.size() > 35) {
    // RDKit✔️❗:     if (tempStr.size() < 39 || tempStr[34] != 'V') {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "CTAB version string invalid at line " << line;
    // RDKit✔️❗:       if (params.strictParsing) {
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (tempStr.substr(34, 5) == "V3000") {
    // RDKit✔️❗:       ctabVersion = 3000;
    // RDKit✔️❗:     } else if (tempStr.substr(34, 5) != "V2000") {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Unsupported CTAB version: '" << tempStr.substr(34, 5)
    // RDKit✔️❗:              << "' at line " << line;
    // RDKit✔️❗:       if (params.strictParsing) {
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (params.parsingSCSRMol) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "SCSR Mol files is not V3000 at line" << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   res->setProp(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit✔️❗:
    // RDKit✔️❗:   Conformer *conf = nullptr;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     if (ctabVersion == 2000) {
    // RDKit✔️❗:       fileComplete = FileParserUtils::ParseV2000CTAB(
    // RDKit✔️❗:           &inStream, line, res.get(), conf, chiralityPossible, nAtoms, nBonds,
    // RDKit✔️❗:           params.strictParsing);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       if (nAtoms != 0 || nBonds != 0) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "V3000 mol blocks should have 0s in the initial counts line. "
    // RDKit✔️❗:                   "(line: "
    // RDKit✔️❗:                << line << ")";
    // RDKit✔️❗:         if (params.strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       auto expectMEND = true;
    // RDKit✔️❗:       auto expectMacroAtoms = false;
    // RDKit✔️❗:       if (params.parsingSCSRMol) {
    // RDKit✔️❗:         expectMEND = false;
    // RDKit✔️❗:         expectMacroAtoms = true;
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       fileComplete = FileParserUtils::ParseV3000CTAB(
    // RDKit✔️❗:           &inStream, line, res.get(), conf, chiralityPossible, nAtoms, nBonds,
    // RDKit✔️❗:           params.strictParsing, expectMEND, expectMacroAtoms);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   } catch (MolFileUnhandledFeatureException &e) {
    // RDKit✔️❗:     // unhandled mol file feature, show an error
    // RDKit✔️❗:     res.reset();
    // RDKit✔️❗:     delete conf;
    // RDKit✔️❗:     conf = nullptr;
    // RDKit✔️❗:     BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: '" << e.what()
    // RDKit✔️❗:                           << "'. Molecule skipped." << std::endl;
    // RDKit✔️❗:
    // RDKit✔️❗:     if (!inStream.eof()) {
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     ++line;
    // RDKit✔️❗:     while (!inStream.eof() && !inStream.fail() &&
    // RDKit✔️❗:            tempStr.substr(0, 6) != "M  END" && tempStr.substr(0, 4) != "$$$$") {
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:       ++line;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     fileComplete = !inStream.eof() || tempStr.substr(0, 6) == "M  END" ||
    // RDKit✔️❗:                    tempStr.substr(0, 4) == "$$$$";
    // RDKit✔️❗:   } catch (FileParseException &e) {
    // RDKit✔️❗:     // catch our exceptions and throw them back after cleanup
    // RDKit✔️❗:     delete conf;
    // RDKit✔️❗:     conf = nullptr;
    // RDKit✔️❗:     throw e;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (!fileComplete) {
    // RDKit✔️❗:     delete conf;
    // RDKit✔️❗:     conf = nullptr;
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout
    // RDKit✔️❗:         << "Problems encountered parsing Mol data, M  END missing around line "
    // RDKit✔️❗:         << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (res) {
    // RDKit✔️❗:     FileParserUtils::finishMolProcessing(res.get(), chiralityPossible, params);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::unique_ptr<RWMol> MolFromMolDataStream
    // END RDKIT CPP BODY: mol_from_mol_data_stream

    let mut header = parse_mol_header(reader, line_number, params)?;
    *line_number += 1;
    let counts_line = read_rdkit_line(reader)?.ok_or_else(|| {
        SdfReadError::Parse(format!("Counts line too short: '' on line{}", *line_number))
    })?;
    let counts = parse_counts_line(&counts_line, *line_number, params)?;
    header.ctab_version = counts.ctab_version;

    let parsed = match counts.ctab_version {
        CtabVersion::V2000 => parse_v2000_ctab(reader, header, counts, line_number, params),
        CtabVersion::V3000 => {
            if counts.atom_count != 0 || counts.bond_count != 0 {
                let message = format!(
                    "V3000 mol blocks should have 0s in the initial counts line. (line: {})",
                    *line_number
                );
                if params.strict_parsing {
                    return Err(SdfReadError::Parse(message));
                }
            }
            parse_v3000_ctab(reader, header, counts, line_number, params)
        }
    }?;
    let molecule = finish_mol_processing(parsed.molecule, parsed.chirality_possible, params)?;
    let molecule = finalize_parsed_coordinate_storage(molecule)?;
    Ok(ParsedMolBlock {
        molecule,
        header: parsed.header,
        chirality_possible: parsed.chirality_possible,
    })
}

fn parse_mol_header<R: BufRead>(
    reader: &mut R,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<MolHeader, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_mol_header
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::unique_ptr<RWMol> MolFromMolDataStream
    // RDKit✔️❗:
    // RDKit✔️❗:   std::string tempStr;
    // RDKit✔️❗:   bool fileComplete = false;
    // RDKit✔️❗:   bool chiralityPossible = false;
    // RDKit✔️❗:   Utils::LocaleSwitcher ls;
    // RDKit✔️❗:   // mol name
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:   if (inStream.eof()) {
    // RDKit✔️❗:     return nullptr;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   auto res = std::make_unique<RWMol>();
    // RDKit✔️❗:   res->setProp(common_properties::_Name, tempStr);
    // RDKit✔️❗:
    // RDKit✔️❗:   // info
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:   res->setProp("_MolFileInfo", tempStr);
    // RDKit✔️❗:   if (tempStr.length() >= 22) {
    // RDKit✔️❗:     std::string dimLabel = tempStr.substr(20, 2);
    // RDKit✔️❗:     // Unless labelled as 3D we assume 2D
    // RDKit✔️❗:     if (dimLabel == "3d" || dimLabel == "3D") {
    // RDKit✔️❗:       res->setProp(common_properties::_3DConf, 1);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   // comments
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:   res->setProp("_MolFileComments", tempStr);
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nAtoms = 0, nBonds = 0, nLists = 0, chiralFlag = 0, nsText = 0,
    // RDKit✔️❗:                nRxnComponents = 0;
    // RDKit✔️❗:   int nReactants = 0, nProducts = 0, nIntermediates = 0;
    // RDKit✔️❗:   (void)nLists;  // read from the file but unused
    // RDKit✔️❗:   (void)nsText;
    // RDKit✔️❗:   (void)nRxnComponents;
    // RDKit✔️❗:   (void)nReactants;
    // RDKit✔️❗:   (void)nProducts;
    // RDKit✔️❗:   (void)nIntermediates;
    // RDKit✔️❗:   // counts line, this is where we really get started
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:
    // RDKit✔️❗:   if (tempStr.size() < 6) {
    // RDKit✔️❗:     if (res) {
    // RDKit✔️❗:       res = nullptr;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Counts line too short: '" << tempStr << "' on line" << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int spos = 0;
    // RDKit✔️❗:   // this needs to go into a try block because if the lexical_cast throws an
    // RDKit✔️❗:   // exception we want to catch throw a different exception
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nAtoms = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     spos = 3;
    // RDKit✔️❗:     nBonds = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     spos = 6;
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << tempStr.substr(spos, 3)
    // RDKit✔️❗:            << "' to unsigned int on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     spos = 6;
    // RDKit✔️❗:     if (tempStr.size() >= 9) {
    // RDKit✔️❗:       nLists = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 12;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       chiralFlag = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 15;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nsText = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 18;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nRxnComponents =
    // RDKit✔️❗:           FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 21;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nReactants = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 24;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nProducts = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 27;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nIntermediates =
    // RDKit✔️❗:           FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     // some SD files (such as some from NCI) lack all the extra information
    // RDKit✔️❗:     // on the header line, so ignore problems parsing there.
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int ctabVersion = 2000;
    // RDKit✔️❗:   if (tempStr.size() > 35) {
    // RDKit✔️❗:     if (tempStr.size() < 39 || tempStr[34] != 'V') {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "CTAB version string invalid at line " << line;
    // RDKit✔️❗:       if (params.strictParsing) {
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (tempStr.substr(34, 5) == "V3000") {
    // RDKit✔️❗:       ctabVersion = 3000;
    // RDKit✔️❗:     } else if (tempStr.substr(34, 5) != "V2000") {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Unsupported CTAB version: '" << tempStr.substr(34, 5)
    // RDKit✔️❗:              << "' at line " << line;
    // RDKit✔️❗:       if (params.strictParsing) {
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (params.parsingSCSRMol) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "SCSR Mol files is not V3000 at line" << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   res->setProp(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit✔️❗:
    // RDKit✔️❗:   Conformer *conf = nullptr;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     if (ctabVersion == 2000) {
    // RDKit✔️❗:       fileComplete = FileParserUtils::ParseV2000CTAB(
    // RDKit✔️❗:           &inStream, line, res.get(), conf, chiralityPossible, nAtoms, nBonds,
    // RDKit✔️❗:           params.strictParsing);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       if (nAtoms != 0 || nBonds != 0) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "V3000 mol blocks should have 0s in the initial counts line. "
    // RDKit✔️❗:                   "(line: "
    // RDKit✔️❗:                << line << ")";
    // RDKit✔️❗:         if (params.strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       auto expectMEND = true;
    // RDKit✔️❗:       auto expectMacroAtoms = false;
    // RDKit✔️❗:       if (params.parsingSCSRMol) {
    // RDKit✔️❗:         expectMEND = false;
    // RDKit✔️❗:         expectMacroAtoms = true;
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       fileComplete = FileParserUtils::ParseV3000CTAB(
    // RDKit✔️❗:           &inStream, line, res.get(), conf, chiralityPossible, nAtoms, nBonds,
    // RDKit✔️❗:           params.strictParsing, expectMEND, expectMacroAtoms);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   } catch (MolFileUnhandledFeatureException &e) {
    // RDKit✔️❗:     // unhandled mol file feature, show an error
    // RDKit✔️❗:     res.reset();
    // RDKit✔️❗:     delete conf;
    // RDKit✔️❗:     conf = nullptr;
    // RDKit✔️❗:     BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: '" << e.what()
    // RDKit✔️❗:                           << "'. Molecule skipped." << std::endl;
    // RDKit✔️❗:
    // RDKit✔️❗:     if (!inStream.eof()) {
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     ++line;
    // RDKit✔️❗:     while (!inStream.eof() && !inStream.fail() &&
    // RDKit✔️❗:            tempStr.substr(0, 6) != "M  END" && tempStr.substr(0, 4) != "$$$$") {
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:       ++line;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     fileComplete = !inStream.eof() || tempStr.substr(0, 6) == "M  END" ||
    // RDKit✔️❗:                    tempStr.substr(0, 4) == "$$$$";
    // RDKit✔️❗:   } catch (FileParseException &e) {
    // RDKit✔️❗:     // catch our exceptions and throw them back after cleanup
    // RDKit✔️❗:     delete conf;
    // RDKit✔️❗:     conf = nullptr;
    // RDKit✔️❗:     throw e;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (!fileComplete) {
    // RDKit✔️❗:     delete conf;
    // RDKit✔️❗:     conf = nullptr;
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout
    // RDKit✔️❗:         << "Problems encountered parsing Mol data, M  END missing around line "
    // RDKit✔️❗:         << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (res) {
    // RDKit✔️❗:     FileParserUtils::finishMolProcessing(res.get(), chiralityPossible, params);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::unique_ptr<RWMol> MolFromMolDataStream
    // END RDKIT CPP BODY: parse_mol_header

    *line_number += 1;
    let Some(name) = read_rdkit_line(reader)? else {
        return Err(SdfReadError::Parse(
            "EOF hit while reading mol name".to_string(),
        ));
    };
    *line_number += 1;
    let info = read_rdkit_line(reader)?.unwrap_or_default();
    *line_number += 1;
    let comments = read_rdkit_line(reader)?.unwrap_or_default();

    let _ = params;
    Ok(MolHeader {
        name: Some(name),
        info,
        comments,
        ctab_version: CtabVersion::V2000,
    })
}

fn molfile_info_marks_3d(info: &str) -> bool {
    info.len() >= 22 && matches!(rdkit_substr(info, 20, 2), "3d" | "3D")
}

fn calculate_rdkit_3d_flag(marked_3d: bool, coords: &[[f64; 3]], chirality_possible: bool) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: calculate3dFlag
    // RDKit✔️✔️: int marked3d = 0;
    // RDKit✔️✔️: if (mol.getPropIfPresent(common_properties::_3DConf, marked3d)) {
    // RDKit✔️✔️:   mol.clearProp(common_properties::_3DConf);
    // RDKit✔️✔️: }
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
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: calculate3dFlag
    let nonzero_z = coords.iter().any(|coord| coord[2] != 0.0);
    if !nonzero_z && marked_3d {
        if chirality_possible {
            return false;
        }
        return true;
    } else if !marked_3d && nonzero_z {
        return true;
    }
    nonzero_z
}

fn resolve_coordinate_3d_flag(detected_3d: bool, mode: SdfCoordinateMode) -> bool {
    match mode {
        SdfCoordinateMode::Preserve => detected_3d,
        SdfCoordinateMode::Require2D => false,
        SdfCoordinateMode::Require3D => true,
    }
}

fn finalize_parsed_coordinate_storage(mut molecule: Molecule) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV2000CTAB / bool ParseV3000CTAB
    // RDKit✔️❗:   conf->set3D(is3d);
    // RDKit✔️❗:   mol->addConformer(conf, true);
    // RDKit✔️❗:   conf = nullptr;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV2000CTAB / bool ParseV3000CTAB
    // RDKit has one conformer store. Parsing and postprocessing use one
    // temporary Conformer3D so source routines see the same unified object;
    // this boundary then maps that object into exactly one COSMolKit store.
    let mut coordinates = molecule.take_coordinate_block_or_clone();
    if !coordinates.conformers_2d.is_empty() || coordinates.conformers_3d.len() != 1 {
        return Err(SdfReadError::Parse(
            "parsed molfile must contain exactly one source conformer before coordinate storage finalization"
                .to_string(),
        ));
    }

    if coordinates.conformers_3d[0].is_3d() {
        coordinates.source_coordinate_dim = Some(CoordinateDimension::ThreeD);
    } else {
        let source = coordinates
            .conformers_3d
            .pop()
            .expect("length checked above");
        let mut conformer = Conformer2D::new(
            source.id(),
            source
                .coordinates()
                .iter()
                .map(|coordinate| [coordinate[0], coordinate[1]])
                .collect(),
        );
        for (key, value) in source.props() {
            conformer = conformer.with_prop(key.clone(), value.clone());
        }
        coordinates.conformers_2d.push(conformer);
        coordinates.source_coordinate_dim = Some(CoordinateDimension::TwoD);
    }
    molecule.replace_coordinate_block(coordinates);
    Ok(molecule)
}

fn parse_counts_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
) -> Result<CountsLine, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_counts_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::unique_ptr<RWMol> MolFromMolDataStream
    // RDKit✔️❗:
    // RDKit✔️❗:   std::string tempStr;
    // RDKit✔️❗:   bool fileComplete = false;
    // RDKit✔️❗:   bool chiralityPossible = false;
    // RDKit✔️❗:   Utils::LocaleSwitcher ls;
    // RDKit✔️❗:   // mol name
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:   if (inStream.eof()) {
    // RDKit✔️❗:     return nullptr;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   auto res = std::make_unique<RWMol>();
    // RDKit✔️❗:   res->setProp(common_properties::_Name, tempStr);
    // RDKit✔️❗:
    // RDKit✔️❗:   // info
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:   res->setProp("_MolFileInfo", tempStr);
    // RDKit✔️❗:   if (tempStr.length() >= 22) {
    // RDKit✔️❗:     std::string dimLabel = tempStr.substr(20, 2);
    // RDKit✔️❗:     // Unless labelled as 3D we assume 2D
    // RDKit✔️❗:     if (dimLabel == "3d" || dimLabel == "3D") {
    // RDKit✔️❗:       res->setProp(common_properties::_3DConf, 1);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   // comments
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:   res->setProp("_MolFileComments", tempStr);
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nAtoms = 0, nBonds = 0, nLists = 0, chiralFlag = 0, nsText = 0,
    // RDKit✔️❗:                nRxnComponents = 0;
    // RDKit✔️❗:   int nReactants = 0, nProducts = 0, nIntermediates = 0;
    // RDKit✔️❗:   (void)nLists;  // read from the file but unused
    // RDKit✔️❗:   (void)nsText;
    // RDKit✔️❗:   (void)nRxnComponents;
    // RDKit✔️❗:   (void)nReactants;
    // RDKit✔️❗:   (void)nProducts;
    // RDKit✔️❗:   (void)nIntermediates;
    // RDKit✔️❗:   // counts line, this is where we really get started
    // RDKit✔️❗:   line++;
    // RDKit✔️❗:   tempStr = getLine(inStream);
    // RDKit✔️❗:
    // RDKit✔️❗:   if (tempStr.size() < 6) {
    // RDKit✔️❗:     if (res) {
    // RDKit✔️❗:       res = nullptr;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Counts line too short: '" << tempStr << "' on line" << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int spos = 0;
    // RDKit✔️❗:   // this needs to go into a try block because if the lexical_cast throws an
    // RDKit✔️❗:   // exception we want to catch throw a different exception
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nAtoms = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     spos = 3;
    // RDKit✔️❗:     nBonds = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     spos = 6;
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << tempStr.substr(spos, 3)
    // RDKit✔️❗:            << "' to unsigned int on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     spos = 6;
    // RDKit✔️❗:     if (tempStr.size() >= 9) {
    // RDKit✔️❗:       nLists = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 12;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       chiralFlag = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 15;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nsText = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 18;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nRxnComponents =
    // RDKit✔️❗:           FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 21;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nReactants = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 24;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nProducts = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     spos = 27;
    // RDKit✔️❗:     if (tempStr.size() >= spos + 3) {
    // RDKit✔️❗:       nIntermediates =
    // RDKit✔️❗:           FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     // some SD files (such as some from NCI) lack all the extra information
    // RDKit✔️❗:     // on the header line, so ignore problems parsing there.
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int ctabVersion = 2000;
    // RDKit✔️❗:   if (tempStr.size() > 35) {
    // RDKit✔️❗:     if (tempStr.size() < 39 || tempStr[34] != 'V') {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "CTAB version string invalid at line " << line;
    // RDKit✔️❗:       if (params.strictParsing) {
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (tempStr.substr(34, 5) == "V3000") {
    // RDKit✔️❗:       ctabVersion = 3000;
    // RDKit✔️❗:     } else if (tempStr.substr(34, 5) != "V2000") {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Unsupported CTAB version: '" << tempStr.substr(34, 5)
    // RDKit✔️❗:              << "' at line " << line;
    // RDKit✔️❗:       if (params.strictParsing) {
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (params.parsingSCSRMol) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "SCSR Mol files is not V3000 at line" << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   res->setProp(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit✔️❗:
    // RDKit✔️❗:   Conformer *conf = nullptr;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     if (ctabVersion == 2000) {
    // RDKit✔️❗:       fileComplete = FileParserUtils::ParseV2000CTAB(
    // RDKit✔️❗:           &inStream, line, res.get(), conf, chiralityPossible, nAtoms, nBonds,
    // RDKit✔️❗:           params.strictParsing);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       if (nAtoms != 0 || nBonds != 0) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "V3000 mol blocks should have 0s in the initial counts line. "
    // RDKit✔️❗:                   "(line: "
    // RDKit✔️❗:                << line << ")";
    // RDKit✔️❗:         if (params.strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       auto expectMEND = true;
    // RDKit✔️❗:       auto expectMacroAtoms = false;
    // RDKit✔️❗:       if (params.parsingSCSRMol) {
    // RDKit✔️❗:         expectMEND = false;
    // RDKit✔️❗:         expectMacroAtoms = true;
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       fileComplete = FileParserUtils::ParseV3000CTAB(
    // RDKit✔️❗:           &inStream, line, res.get(), conf, chiralityPossible, nAtoms, nBonds,
    // RDKit✔️❗:           params.strictParsing, expectMEND, expectMacroAtoms);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   } catch (MolFileUnhandledFeatureException &e) {
    // RDKit✔️❗:     // unhandled mol file feature, show an error
    // RDKit✔️❗:     res.reset();
    // RDKit✔️❗:     delete conf;
    // RDKit✔️❗:     conf = nullptr;
    // RDKit✔️❗:     BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: '" << e.what()
    // RDKit✔️❗:                           << "'. Molecule skipped." << std::endl;
    // RDKit✔️❗:
    // RDKit✔️❗:     if (!inStream.eof()) {
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     ++line;
    // RDKit✔️❗:     while (!inStream.eof() && !inStream.fail() &&
    // RDKit✔️❗:            tempStr.substr(0, 6) != "M  END" && tempStr.substr(0, 4) != "$$$$") {
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:       ++line;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     fileComplete = !inStream.eof() || tempStr.substr(0, 6) == "M  END" ||
    // RDKit✔️❗:                    tempStr.substr(0, 4) == "$$$$";
    // RDKit✔️❗:   } catch (FileParseException &e) {
    // RDKit✔️❗:     // catch our exceptions and throw them back after cleanup
    // RDKit✔️❗:     delete conf;
    // RDKit✔️❗:     conf = nullptr;
    // RDKit✔️❗:     throw e;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (!fileComplete) {
    // RDKit✔️❗:     delete conf;
    // RDKit✔️❗:     conf = nullptr;
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout
    // RDKit✔️❗:         << "Problems encountered parsing Mol data, M  END missing around line "
    // RDKit✔️❗:         << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (res) {
    // RDKit✔️❗:     FileParserUtils::finishMolProcessing(res.get(), chiralityPossible, params);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::unique_ptr<RWMol> MolFromMolDataStream
    // END RDKIT CPP BODY: parse_counts_line

    if line.len() < 6 {
        return Err(SdfReadError::Parse(format!(
            "Counts line too short: '{line}' on line{line_number}"
        )));
    }

    let mut spos = 0_usize;
    let atom_count = parse_rdkit_unsigned(rdkit_substr(line, spos, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to unsigned int on line {line_number}",
            rdkit_substr(line, spos, 3)
        ))
    })? as usize;
    spos = 3;
    let bond_count = parse_rdkit_unsigned(rdkit_substr(line, spos, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to unsigned int on line {line_number}",
            rdkit_substr(line, spos, 3)
        ))
    })? as usize;

    let _n_lists = if line.len() >= 9 {
        parse_rdkit_unsigned(rdkit_substr(line, 6, 3), true).unwrap_or(0)
    } else {
        0
    };
    let chiral_flag = if line.len() >= 15 {
        parse_rdkit_unsigned(rdkit_substr(line, 12, 3), true).unwrap_or(0)
    } else {
        0
    };
    let _ns_text = if line.len() >= 18 {
        parse_rdkit_unsigned(rdkit_substr(line, 15, 3), true).unwrap_or(0)
    } else {
        0
    };
    let _n_rxn_components = if line.len() >= 21 {
        parse_rdkit_unsigned(rdkit_substr(line, 18, 3), true).unwrap_or(0)
    } else {
        0
    };
    let _n_reactants = if line.len() >= 24 {
        parse_rdkit_unsigned(rdkit_substr(line, 21, 3), true).unwrap_or(0)
    } else {
        0
    };
    let _n_products = if line.len() >= 27 {
        parse_rdkit_unsigned(rdkit_substr(line, 24, 3), true).unwrap_or(0)
    } else {
        0
    };
    let _n_intermediates = if line.len() >= 30 {
        parse_rdkit_unsigned(rdkit_substr(line, 27, 3), true).unwrap_or(0)
    } else {
        0
    };

    let mut ctab_version = CtabVersion::V2000;
    if line.len() > 35 {
        if line.len() < 39 || line.as_bytes().get(34) != Some(&b'V') {
            if params.strict_parsing {
                return Err(SdfReadError::Parse(format!(
                    "CTAB version string invalid at line {line_number}"
                )));
            }
        } else {
            match rdkit_substr(line, 34, 5) {
                "V3000" => ctab_version = CtabVersion::V3000,
                "V2000" => {}
                other => {
                    if params.strict_parsing {
                        return Err(SdfReadError::Parse(format!(
                            "Unsupported CTAB version: '{other}' at line {line_number}"
                        )));
                    }
                }
            }
        }
    }

    Ok(CountsLine {
        atom_count,
        bond_count,
        chiral_flag,
        ctab_version,
        v3000_sgroup_count: 0,
        v3000_obj3d_count: 0,
    })
}

fn parse_v2000_ctab<R: BufRead>(
    reader: &mut R,
    header: MolHeader,
    counts: CountsLine,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<ParsedMolBlock, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_ctab
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV2000CTAB
    // RDKit✔️❗:
    // RDKit✔️❗:   conf = new Conformer(nAtoms);
    // RDKit✔️❗:
    // RDKit✔️❗:   if (nAtoms == 0) {
    // RDKit✔️❗:     conf->set3D(false);
    // RDKit✔️❗:   } else {
    // RDKit✔️❗:     ParseMolBlockAtoms(inStream, line, nAtoms, mol, conf, strictParsing);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   ParseMolBlockBonds(inStream, line, nBonds, mol, chiralityPossible);
    // RDKit❗✔️:
    // RDKit❗✔️:   auto is3d = calculate3dFlag(*mol, *conf, chiralityPossible);
    // RDKit✔️❗:   conf->set3D(is3d);
    // RDKit✔️❗:   mol->addConformer(conf, true);
    // RDKit✔️❗:   conf = nullptr;
    // RDKit❗✔️:
    // RDKit✔️❗:   bool fileComplete =
    // RDKit✔️❗:       ParseMolBlockProperties(inStream, line, mol, strictParsing);
    // RDKit✔️❗:   return fileComplete;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV2000CTAB
    // END RDKIT CPP BODY: parse_v2000_ctab

    // RDKit stores coordinates in an owning Conformer allocated before atom
    // parsing. COSMolKit keeps the parsed coordinates in a Vec and installs one
    // Conformer3D after the topology is built; the observable modeled state is
    // equivalent for atoms, bonds, coordinates, and the current 3D flag subset.
    let mut atom_lines = if counts.atom_count == 0 {
        Vec::new()
    } else {
        parse_v2000_atom_block(reader, counts.atom_count, line_number, params)?
    };
    let mut bond_lines = parse_v2000_bond_block(reader, counts.bond_count, line_number, params)?;
    let chirality_possible = bond_lines.iter().any(|bond| {
        !matches!(
            bond.spec.direction(),
            BondDirection::None | BondDirection::Unknown
        )
    });

    let mut property_state = V2000PropertyState {
        first_charge_line: true,
        sgroups: BTreeMap::new(),
        scd_counter: 0,
        last_data_sgroup: 0,
        current_data_field: String::new(),
        molecule_props: BTreeMap::new(),
    };
    let file_complete = parse_v2000_property_block(
        reader,
        line_number,
        params,
        &mut atom_lines,
        &mut bond_lines,
        &mut property_state,
    )?;
    if !file_complete {
        return Err(SdfReadError::Parse(format!(
            "Problems encountered parsing Mol data, M  END missing around line {}",
            *line_number
        )));
    }

    let coords = atom_lines
        .iter()
        .map(|atom| atom.coord_3d)
        .collect::<Vec<_>>();
    let is_3d = resolve_coordinate_3d_flag(
        calculate_rdkit_3d_flag(
            molfile_info_marks_3d(&header.info),
            &coords,
            chirality_possible,
        ),
        params.coordinate_mode,
    );

    let mut builder = MoleculeBuilder::new();
    if let Some(name) = header.name.as_deref() {
        builder = builder.with_name(name);
    }
    builder = builder
        .with_property("_MolFileChiralFlag", counts.chiral_flag.to_string())
        .with_property("_MolFileInfoLine", header.info.clone())
        .with_property("_MolFileComments", header.comments.clone());
    for (key, value) in &property_state.molecule_props {
        builder = builder.with_property(key.clone(), value.clone());
    }
    for atom in atom_lines {
        builder.add_atom(atom.spec);
    }
    for bond in bond_lines {
        builder.add_bond(bond.spec).map_err(molecule_build_error)?;
    }
    for substance_group in property_state.sgroups.into_values() {
        builder
            .add_substance_group(substance_group)
            .map_err(molecule_build_error)?;
    }
    builder
        .add_conformer(Conformer3D::new(0, coords, is_3d))
        .map_err(molecule_build_error)?;
    let molecule = builder
        .with_topology_trust(crate::TopologyTrust::TrustedGraph)
        .build()
        .map_err(molecule_build_error)?;

    Ok(ParsedMolBlock {
        molecule,
        header,
        chirality_possible,
    })
}

fn parse_v2000_atom_block<R: BufRead>(
    reader: &mut R,
    atom_count: usize,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<Vec<V2000AtomLine>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_atom_block
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseMolBlockAtoms
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:   PRECONDITION(conf, "bad conformer");
    // RDKit✔️❗:   for (unsigned int i = 1; i <= nAtoms; ++i) {
    // RDKit✔️❗:     ++line;
    // RDKit✔️❗:     std::string tempStr = getLine(inStream);
    // RDKit✔️❗:     if (inStream->eof()) {
    // RDKit✔️❗:       throw FileParseException("EOF hit while reading atoms");
    // RDKit✔️❗:     }
    // RDKit✔️❗:     RDGeom::Point3D pos;
    // RDKit✔️❗:     Atom *atom = ParseMolFileAtomLine(tempStr, pos, line, strictParsing);
    // RDKit✔️❗:     unsigned int aid = mol->addAtom(atom, false, true);
    // RDKit✔️❗:     conf->setAtomPos(aid, pos);
    // RDKit✔️❗:     mol->setAtomBookmark(atom, i);
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseMolBlockAtoms
    // END RDKIT CPP BODY: parse_v2000_atom_block

    let mut atoms = Vec::with_capacity(atom_count);
    for _ in 1..=atom_count {
        *line_number += 1;
        let Some(line) = read_rdkit_line(reader)? else {
            return Err(SdfReadError::Parse(
                "EOF hit while reading atoms".to_string(),
            ));
        };
        atoms.push(parse_v2000_atom_line(&line, *line_number, params)?);
    }
    Ok(atoms)
}

fn parse_v2000_atom_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
) -> Result<V2000AtomLine, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_atom_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: Atom *ParseMolFileAtomLine
    // RDKit✔️❗:
    // RDKit✔️❗:   std::string symb;
    // RDKit✔️❗:   int massDiff, chg, hCount;
    // RDKit✔️❗:
    // RDKit✔️❗:   if ((strictParsing && text.size() < 34) || text.size() < 32) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Atom line too short: '" << text << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     pos.x = FileParserUtils::toDouble(text.substr(0, 10));
    // RDKit✔️❗:     pos.y = FileParserUtils::toDouble(text.substr(10, 10));
    // RDKit✔️❗:     pos.z = FileParserUtils::toDouble(text.substr(20, 10));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot process coordinates on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   symb = text.substr(31, 3);
    // RDKit✔️❗:   boost::trim(symb);
    // RDKit✔️❗:
    // RDKit✔️❗:   // REVIEW: should we handle missing fields at the end of the line?
    // RDKit✔️❗:   massDiff = 0;
    // RDKit✔️❗:   if (text.size() >= 36 && text.substr(34, 2) != " 0") {
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       massDiff = FileParserUtils::toInt(text.substr(34, 2), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(34, 2) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   chg = 0;
    // RDKit✔️❗:   if (text.size() >= 39 && text.substr(36, 3) != "  0") {
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       chg = FileParserUtils::toInt(text.substr(36, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(36, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   hCount = 0;
    // RDKit✔️❗:   if (text.size() >= 45 && text.substr(42, 3) != "  0") {
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       hCount = FileParserUtils::toInt(text.substr(42, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(42, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   std::unique_ptr<Atom> res(new Atom);
    // RDKit✔️❗:   bool isComplexQueryName =
    // RDKit✔️❗:       std::find(complexQueries.begin(), complexQueries.end(), symb) !=
    // RDKit✔️❗:       complexQueries.end();
    // RDKit✔️❗:   if (isComplexQueryName || symb == "L" || symb == "*" || symb == "LP" ||
    // RDKit✔️❗:       symb == "R" || symb == "R#" ||
    // RDKit✔️❗:       (symb[0] == 'R' && symb >= "R0" && symb <= "R99")) {
    // RDKit✔️❗:     if (isComplexQueryName || symb == "*" || symb == "R") {
    // RDKit✔️❗:       auto *query = new QueryAtom(0);
    // RDKit✔️❗:       if (symb == "*" || symb == "R") {
    // RDKit✔️❗:         // according to the MDL spec, these match anything
    // RDKit✔️❗:         query->setQuery(makeAtomNullQuery());
    // RDKit✔️❗:       } else if (isComplexQueryName) {
    // RDKit✔️❗:         convertComplexNameToQuery(query, symb);
    // RDKit✔️❗:       }
    // RDKit✔️❗:       res.reset(query);
    // RDKit✔️❗:       // queries have no implicit Hs:
    // RDKit✔️❗:       res->setNoImplicit(true);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       res->setAtomicNum(0);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     if (massDiff == 0 && symb[0] == 'R') {
    // RDKit✔️❗:       if (symb.length() > 1 && symb >= "R0" && symb <= "R99") {
    // RDKit✔️❗:         std::string rlabel = "";
    // RDKit✔️❗:         rlabel = symb.substr(1, symb.length() - 1);
    // RDKit✔️❗:         int rnumber;
    // RDKit✔️❗:         try {
    // RDKit✔️❗:           rnumber = boost::lexical_cast<int>(rlabel);
    // RDKit✔️❗:         } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:           rnumber = -1;
    // RDKit✔️❗:         }
    // RDKit✔️❗:         if (rnumber >= 0) {
    // RDKit✔️❗:           res->setIsotope(rnumber);
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:     if (symb[0] == 'R') {
    // RDKit✔️❗:       // we used to skip R# here because that really should be handled by an
    // RDKit✔️❗:       // RGP spec, but that turned out to not be permissive enough... <sigh>
    // RDKit✔️❗:       setRGPProps(symb, res.get());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   } else if (symb == "D") {  // mol blocks support "D" and "T" as shorthand...
    // RDKit✔️❗:                              // handle that.
    // RDKit✔️❗:     res->setAtomicNum(1);
    // RDKit✔️❗:     res->setIsotope(2);
    // RDKit✔️❗:   } else if (symb == "T") {  // mol blocks support "D" and "T" as shorthand...
    // RDKit✔️❗:                              // handle that.
    // RDKit✔️❗:     res->setAtomicNum(1);
    // RDKit✔️❗:     res->setIsotope(3);
    // RDKit✔️❗:   } else if (symb == "Pol" || symb == "Mod") {
    // RDKit✔️❗:     res->setAtomicNum(0);
    // RDKit✔️❗:     res->setProp(common_properties::dummyLabel, symb);
    // RDKit✔️❗:   } else if (GenericGroups::genericMatchers.find(symb) !=
    // RDKit✔️❗:              GenericGroups::genericMatchers.end()) {
    // RDKit✔️❗:     res.reset(new QueryAtom(0));
    // RDKit✔️❗:     res->setProp(common_properties::atomLabel, std::string(symb));
    // RDKit✔️❗:   } else {
    // RDKit✔️❗:     lookupAtomicNumber(res.get(), symb, strictParsing);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   // res->setPos(pX,pY,pZ);
    // RDKit✔️❗:   if (chg != 0) {
    // RDKit✔️❗:     res->setFormalCharge(4 - chg);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (hCount >= 1) {
    // RDKit✔️❗:     if (!res->hasQuery()) {
    // RDKit✔️❗:       auto qatom = new QueryAtom(*res);
    // RDKit✔️❗:       res.reset(qatom);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     res->setNoImplicit(true);
    // RDKit✔️❗:     if (hCount > 1) {
    // RDKit✔️❗:       ATOM_EQUALS_QUERY *oq = makeAtomImplicitHCountQuery(hCount - 1);
    // RDKit✔️❗:       auto nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(
    // RDKit✔️❗:           hCount - 1, oq->getDataFunc(),
    // RDKit✔️❗:           std::string("less_") + oq->getDescription());
    // RDKit✔️❗:       res->expandQuery(nq);
    // RDKit✔️❗:       delete oq;
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       res->expandQuery(makeAtomImplicitHCountQuery(0));
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (massDiff != 0) {
    // RDKit✔️❗:     int defIso =
    // RDKit✔️❗:         PeriodicTable::getTable()->getMostCommonIsotope(res->getAtomicNum());
    // RDKit✔️❗:     int dIso = defIso + massDiff;
    // RDKit✔️❗:     if (dIso < 0) {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:           << " atom " << res->getIdx()
    // RDKit✔️❗:           << " has a negative isotope offset. line:  " << line << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     res->setIsotope(dIso);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (text.size() >= 42 && text.substr(39, 3) != "  0") {
    // RDKit✔️❗:     int parity = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       parity = FileParserUtils::toInt(text.substr(39, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(39, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     res->setProp(common_properties::molParity, parity);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (text.size() >= 48 && text.substr(45, 3) != "  0") {
    // RDKit✔️❗:     int stereoCare = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       stereoCare = FileParserUtils::toInt(text.substr(45, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(45, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     res->setProp(common_properties::molStereoCare, stereoCare);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (text.size() >= 51 && text.substr(48, 3) != "  0") {
    // RDKit✔️❗:     int totValence = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       totValence = FileParserUtils::toInt(text.substr(48, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(48, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     if (totValence != 0) {
    // RDKit✔️❗:       // only set if it's a non-default value
    // RDKit✔️❗:       res->setProp(common_properties::molTotValence, totValence);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (text.size() >= 57 && text.substr(54, 3) != "  0") {
    // RDKit✔️❗:     int rxnRole = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       rxnRole = FileParserUtils::toInt(text.substr(54, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(54, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     if (rxnRole != 0) {
    // RDKit✔️❗:       // only set if it's a non-default value
    // RDKit✔️❗:       res->setProp(common_properties::molRxnRole, rxnRole);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (text.size() >= 60 && text.substr(57, 3) != "  0") {
    // RDKit✔️❗:     int rxnComponent = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       rxnComponent = FileParserUtils::toInt(text.substr(57, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(57, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     if (rxnComponent != 0) {
    // RDKit✔️❗:       // only set if it's a non-default value
    // RDKit✔️❗:       res->setProp(common_properties::molRxnComponent, rxnComponent);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (text.size() >= 63 && text.substr(60, 3) != "  0") {
    // RDKit✔️❗:     int atomMapNumber = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       atomMapNumber = FileParserUtils::toInt(text.substr(60, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(60, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     res->setProp(common_properties::molAtomMapNumber, atomMapNumber);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (text.size() >= 66 && text.substr(63, 3) != "  0") {
    // RDKit✔️❗:     int inversionFlag = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       inversionFlag = FileParserUtils::toInt(text.substr(63, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(63, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     res->setProp(common_properties::molInversionFlag, inversionFlag);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (text.size() >= 69 && text.substr(66, 3) != "  0") {
    // RDKit✔️❗:     int exactChangeFlag = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       exactChangeFlag = FileParserUtils::toInt(text.substr(66, 3), true);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(66, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     res->setProp(common_properties::molRxnExactChange, exactChangeFlag);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res.release();
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: Atom *ParseMolFileAtomLine
    // END RDKIT CPP BODY: parse_v2000_atom_line

    // RDKit owns a heap Atom/QueryAtom here; COSMolKit builds an AtomSpec plus
    // coordinate payload and lets MoleculeBuilder assign the final AtomId.
    if (params.strict_parsing && line.len() < 34) || line.len() < 32 {
        return Err(SdfReadError::Parse(format!(
            "Atom line too short: '{line}' on line {line_number}"
        )));
    }

    let x = parse_rdkit_double(rdkit_substr(line, 0, 10), true).map_err(|_| {
        SdfReadError::Parse(format!("Cannot process coordinates on line {line_number}"))
    })?;
    let y = parse_rdkit_double(rdkit_substr(line, 10, 10), true).map_err(|_| {
        SdfReadError::Parse(format!("Cannot process coordinates on line {line_number}"))
    })?;
    let z = parse_rdkit_double(rdkit_substr(line, 20, 10), true).map_err(|_| {
        SdfReadError::Parse(format!("Cannot process coordinates on line {line_number}"))
    })?;
    let symbol = rdkit_substr(line, 31, 3).trim();

    let mass_diff = if line.len() >= 36 && rdkit_substr(line, 34, 2) != " 0" {
        parse_rdkit_int(rdkit_substr(line, 34, 2), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, 34, 2)
            ))
        })?
    } else {
        0
    };
    let charge_code = if line.len() >= 39 && rdkit_substr(line, 36, 3) != "  0" {
        parse_rdkit_int(rdkit_substr(line, 36, 3), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, 36, 3)
            ))
        })?
    } else {
        0
    };
    let h_count = if line.len() >= 45 && rdkit_substr(line, 42, 3) != "  0" {
        parse_rdkit_int(rdkit_substr(line, 42, 3), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, 42, 3)
            ))
        })?
    } else {
        0
    };
    let mol_parity = if line.len() >= 42 && rdkit_substr(line, 39, 3) != "  0" {
        Some(
            parse_rdkit_int(rdkit_substr(line, 39, 3), true).map_err(|_| {
                SdfReadError::Parse(format!(
                    "Cannot convert '{}' to int on line {line_number}",
                    rdkit_substr(line, 39, 3)
                ))
            })?,
        )
    } else {
        None
    };
    let mol_inversion_flag = if line.len() >= 66 && rdkit_substr(line, 63, 3) != "  0" {
        Some(
            parse_rdkit_int(rdkit_substr(line, 63, 3), true).map_err(|_| {
                SdfReadError::Parse(format!(
                    "Cannot convert '{}' to int on line {line_number}",
                    rdkit_substr(line, 63, 3)
                ))
            })?,
        )
    } else {
        None
    };

    let complex_query = molfile_complex_atom_query(symbol);
    let atomic_number = if complex_query.is_some() {
        0
    } else {
        atomic_number_from_mol_symbol(symbol, params.strict_parsing)?
    };
    let mut spec = if atomic_number == 0 {
        AtomSpec::new(Element::DUMMY)
    } else {
        AtomSpec::new(Element::from_atomic_number(atomic_number).unwrap())
    };
    let mut query: Option<QueryNode<AtomQueryPredicate>> = None;

    if let Some(complex_query) = complex_query {
        query = Some(complex_query);
        spec = spec.with_no_implicit(true);
    } else if matches!(symbol, "*" | "R") {
        query = Some(QueryNode::predicate(AtomQueryPredicate::Any));
        spec = spec.with_no_implicit(true);
    } else if symbol == "D" {
        spec = spec.with_isotope(2);
    } else if symbol == "T" {
        spec = spec.with_isotope(3);
    } else if let Some(rest) = symbol.strip_prefix('R')
        && !rest.is_empty()
        && rest.bytes().all(|byte| byte.is_ascii_digit())
        && rest.len() <= 2
    {
        spec = spec.with_isotope(rest.parse::<u16>().unwrap_or(0));
        query = Some(QueryNode::predicate(AtomQueryPredicate::RGroupLabel(
            rest.parse::<u32>().unwrap_or(0),
        )));
    } else if matches!(symbol, "A" | "Q" | "L" | "LP" | "R#") {
        query = Some(QueryNode::predicate(AtomQueryPredicate::MolFileAlias(
            symbol.to_string(),
        )));
        spec = spec.with_no_implicit(true);
    }
    if symbol.starts_with('R') {
        spec = spec.with_prop("dummyLabel", symbol.to_string());
    }

    if charge_code != 0 {
        spec = spec.with_formal_charge((4 - charge_code) as i8);
    }

    if let Some(mol_parity) = mol_parity {
        spec = spec
            .with_mol_parity(mol_parity)
            .with_prop("molParity", mol_parity.to_string());
    }

    if line.len() >= 48 && rdkit_substr(line, 45, 3) != "  0" {
        let stereo_care = parse_rdkit_int(rdkit_substr(line, 45, 3), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, 45, 3)
            ))
        })?;
        spec = spec.with_prop("molStereoCare".to_string(), stereo_care.to_string());
    }

    if line.len() >= 51 && rdkit_substr(line, 48, 3) != "  0" {
        let tot_valence = parse_rdkit_int(rdkit_substr(line, 48, 3), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, 48, 3)
            ))
        })?;
        if tot_valence != 0 {
            spec = spec.with_prop("molTotValence".to_string(), tot_valence.to_string());
        }
    }

    if let Some(mol_inversion_flag) = mol_inversion_flag {
        spec = spec.with_mol_inversion_flag(mol_inversion_flag);
    }

    if h_count >= 1 {
        spec = spec.with_no_implicit(true);
        let predicate = if h_count > 1 {
            AtomQueryPredicate::ImplicitHydrogenCountLessEqual((h_count - 1) as u8)
        } else {
            AtomQueryPredicate::ImplicitHydrogenCount(0)
        };
        query = Some(match query {
            Some(existing) => QueryNode::and(vec![existing, QueryNode::predicate(predicate)]),
            None => QueryNode::predicate(predicate),
        });
    }

    if mass_diff != 0 {
        let Some(default_isotope) = most_common_isotope(atomic_number) else {
            return Err(SdfReadError::Parse(format!(
                "Cannot apply mass difference {mass_diff} to unknown atomic number {atomic_number}"
            )));
        };
        let isotope = default_isotope + mass_diff;
        if isotope >= 0 {
            spec = spec.with_isotope(isotope as u16);
        }
    }

    if line.len() >= 63 && rdkit_substr(line, 60, 3) != "  0" {
        let atom_map = parse_rdkit_int(rdkit_substr(line, 60, 3), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, 60, 3)
            ))
        })?;
        spec = spec
            .with_atom_map(atom_map as u32)
            .with_prop("molAtomMapNumber", atom_map.to_string());
    }

    if let Some(query) = query {
        spec = spec.with_query(query);
    }

    Ok(V2000AtomLine {
        line_number,
        text: line.to_string(),
        spec,
        coord_3d: [x, y, z],
    })
}

fn molfile_complex_atom_query(symbol: &str) -> Option<QueryNode<AtomQueryPredicate>> {
    convert_complex_name_to_query(symbol).ok()
}

fn parse_v2000_bond_block<R: BufRead>(
    reader: &mut R,
    bond_count: usize,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<Vec<V2000BondLine>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_bond_block
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseMolBlockBonds
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:   for (unsigned int i = 1; i <= nBonds; ++i) {
    // RDKit✔️❗:     ++line;
    // RDKit✔️❗:     std::string tempStr = getLine(inStream);
    // RDKit✔️❗:     if (inStream->eof()) {
    // RDKit✔️❗:       throw FileParseException("EOF hit while reading bonds");
    // RDKit✔️❗:     }
    // RDKit✔️❗:     Bond *bond = ParseMolFileBondLine(tempStr, line);
    // RDKit✔️❗:     // if we got an aromatic bond set the flag on the bond and the connected
    // RDKit✔️❗:     // atoms
    // RDKit✔️❗:     if (bond->getBondType() == Bond::AROMATIC) {
    // RDKit✔️❗:       bond->setIsAromatic(true);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     // if the bond might have chirality info associated with it, set a flag:
    // RDKit✔️❗:     if (bond->getBondDir() != Bond::NONE &&
    // RDKit✔️❗:         bond->getBondDir() != Bond::UNKNOWN) {
    // RDKit✔️❗:       chiralityPossible = true;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     // v2k has no way to set stereoCare on bonds, so set the property if both
    // RDKit✔️❗:     // the beginning and end atoms have it set:
    // RDKit✔️❗:     int care1 = 0;
    // RDKit✔️❗:     int care2 = 0;
    // RDKit✔️❗:     if (!bond->hasProp(common_properties::molStereoCare) &&
    // RDKit✔️❗:         mol->getAtomWithIdx(bond->getBeginAtomIdx())
    // RDKit✔️❗:             ->getPropIfPresent(common_properties::molStereoCare, care1) &&
    // RDKit✔️❗:         mol->getAtomWithIdx(bond->getEndAtomIdx())
    // RDKit✔️❗:             ->getPropIfPresent(common_properties::molStereoCare, care2)) {
    // RDKit✔️❗:       if (care1 && care2) {
    // RDKit✔️❗:         bond->setProp(common_properties::molStereoCare, 1);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:     mol->addBond(bond, true);
    // RDKit✔️❗:     mol->setBondBookmark(bond, i);
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseMolBlockBonds
    // END RDKIT CPP BODY: parse_v2000_bond_block

    let mut bonds = Vec::with_capacity(bond_count);
    for _ in 1..=bond_count {
        *line_number += 1;
        let Some(line) = read_rdkit_line(reader)? else {
            return Err(SdfReadError::Parse(
                "EOF hit while reading bonds".to_string(),
            ));
        };
        bonds.push(parse_v2000_bond_line(&line, *line_number, params)?);
    }
    Ok(bonds)
}

fn parse_v2000_bond_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
) -> Result<V2000BondLine, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_bond_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: Bond *ParseMolFileBondLine
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int idx1, idx2, bType, stereo;
    // RDKit✔️❗:   int spos = 0;
    // RDKit✔️❗:
    // RDKit✔️❗:   if (text.size() < 9) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Bond line too short: '" << text << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     idx1 = FileParserUtils::toUnsigned(text.substr(spos, 3));
    // RDKit✔️❗:     spos += 3;
    // RDKit✔️❗:     idx2 = FileParserUtils::toUnsigned(text.substr(spos, 3));
    // RDKit✔️❗:     spos += 3;
    // RDKit✔️❗:     bType = FileParserUtils::toUnsigned(text.substr(spos, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(spos, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   // adjust the numbering
    // RDKit✔️❗:   idx1--;
    // RDKit✔️❗:   idx2--;
    // RDKit✔️❗:
    // RDKit✔️❗:   Bond::BondType type;
    // RDKit✔️❗:   Bond *res = nullptr;
    // RDKit✔️❗:   switch (bType) {
    // RDKit✔️❗:     case 1:
    // RDKit✔️❗:       type = Bond::SINGLE;
    // RDKit✔️❗:       res = new Bond;
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     case 2:
    // RDKit✔️❗:       type = Bond::DOUBLE;
    // RDKit✔️❗:       res = new Bond;
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     case 3:
    // RDKit✔️❗:       type = Bond::TRIPLE;
    // RDKit✔️❗:       res = new Bond;
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     case 4:
    // RDKit✔️❗:       type = Bond::AROMATIC;
    // RDKit✔️❗:       res = new Bond;
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     case 9:
    // RDKit✔️❗:       type = Bond::DATIVE;
    // RDKit✔️❗:       res = new Bond;
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     case 0:
    // RDKit✔️❗:       type = Bond::UNSPECIFIED;
    // RDKit✔️❗:       res = new Bond;
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:           << "bond with order 0 found on line " << line
    // RDKit✔️❗:           << ". This is not part of the MDL specification." << std::endl;
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     default:
    // RDKit✔️❗:       type = Bond::UNSPECIFIED;
    // RDKit✔️❗:       // it's a query bond of some type
    // RDKit✔️❗:       res = new QueryBond;
    // RDKit✔️❗:       if (bType == 8) {
    // RDKit✔️❗:         BOND_NULL_QUERY *q;
    // RDKit✔️❗:         q = makeBondNullQuery();
    // RDKit✔️❗:         res->setQuery(q);
    // RDKit✔️❗:       } else if (bType == 5) {
    // RDKit✔️❗:         res->setQuery(makeSingleOrDoubleBondQuery());
    // RDKit✔️❗:         res->setProp(common_properties::_MolFileBondQuery, 1);
    // RDKit✔️❗:       } else if (bType == 6) {
    // RDKit✔️❗:         res->setQuery(makeSingleOrAromaticBondQuery());
    // RDKit✔️❗:         res->setProp(common_properties::_MolFileBondQuery, 1);
    // RDKit✔️❗:       } else if (bType == 7) {
    // RDKit✔️❗:         res->setQuery(makeDoubleOrAromaticBondQuery());
    // RDKit✔️❗:         res->setProp(common_properties::_MolFileBondQuery, 1);
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         BOND_NULL_QUERY *q;
    // RDKit✔️❗:         q = makeBondNullQuery();
    // RDKit✔️❗:         res->setQuery(q);
    // RDKit✔️❗:         BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:             << "unrecognized query bond type, " << bType << ", found on line "
    // RDKit✔️❗:             << line << ". Using an \"any\" query." << std::endl;
    // RDKit✔️❗:       }
    // RDKit✔️❗:       break;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   res->setBeginAtomIdx(idx1);
    // RDKit✔️❗:   res->setEndAtomIdx(idx2);
    // RDKit✔️❗:   res->setBondType(type);
    // RDKit✔️❗:   res->setProp(common_properties::_MolFileBondType, bType);
    // RDKit✔️❗:
    // RDKit✔️❗:   if (text.size() >= 12 && text.substr(9, 3) != "  0") {
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       stereo = FileParserUtils::toUnsigned(text.substr(9, 3));
    // RDKit✔️❗:       switch (stereo) {
    // RDKit✔️❗:         case 0:
    // RDKit✔️❗:           res->setBondDir(Bond::NONE);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 1:
    // RDKit✔️❗:           res->setBondDir(Bond::BEGINWEDGE);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 6:
    // RDKit✔️❗:           res->setBondDir(Bond::BEGINDASH);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 3:  // "either" double bond
    // RDKit✔️❗:           res->setBondDir(Bond::EITHERDOUBLE);
    // RDKit✔️❗:           res->setStereo(Bond::STEREOANY);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 4:  // "either" single bond
    // RDKit✔️❗:           res->setBondDir(Bond::UNKNOWN);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:       }
    // RDKit✔️❗:       res->setProp(common_properties::_MolFileBondStereo, stereo);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       ;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (text.size() >= 18 && text.substr(15, 3) != "  0") {
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       int topology = FileParserUtils::toInt(text.substr(15, 3));
    // RDKit✔️❗:       if (topology) {
    // RDKit✔️❗:         if (!res->hasQuery()) {
    // RDKit✔️❗:           auto *qBond = new QueryBond(*res);
    // RDKit✔️❗:           delete res;
    // RDKit✔️❗:           res = qBond;
    // RDKit✔️❗:         }
    // RDKit✔️❗:         BOND_EQUALS_QUERY *q = makeBondIsInRingQuery();
    // RDKit✔️❗:         switch (topology) {
    // RDKit✔️❗:           case 1:
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           case 2:
    // RDKit✔️❗:             q->setNegation(true);
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           default:
    // RDKit✔️❗:             std::ostringstream errout;
    // RDKit✔️❗:             errout << "Unrecognized bond topology specifier: " << topology
    // RDKit✔️❗:                    << " on line " << line;
    // RDKit✔️❗:             throw FileParseException(errout.str());
    // RDKit✔️❗:         }
    // RDKit✔️❗:         res->expandQuery(q);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       ;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (text.size() >= 21 && text.substr(18, 3) != "  0") {
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       int reactStatus = FileParserUtils::toInt(text.substr(18, 3));
    // RDKit✔️❗:       res->setProp("molReactStatus", reactStatus);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       ;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: Bond *ParseMolFileBondLine
    // END RDKIT CPP BODY: parse_v2000_bond_line

    let _ = params;
    if line.len() < 9 {
        return Err(SdfReadError::Parse(format!(
            "Bond line too short: '{line}' on line {line_number}"
        )));
    }

    let mut spos = 0_usize;
    let idx1 = parse_rdkit_unsigned(rdkit_substr(line, spos, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {line_number}",
            rdkit_substr(line, spos, 3)
        ))
    })?;
    spos += 3;
    let idx2 = parse_rdkit_unsigned(rdkit_substr(line, spos, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {line_number}",
            rdkit_substr(line, spos, 3)
        ))
    })?;
    spos += 3;
    let bond_type = parse_rdkit_unsigned(rdkit_substr(line, spos, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {line_number}",
            rdkit_substr(line, spos, 3)
        ))
    })?;

    let begin = AtomId::new(idx1.saturating_sub(1) as usize);
    let end = AtomId::new(idx2.saturating_sub(1) as usize);
    let (order, mut query) = match bond_type {
        1 => (BondOrder::Single, None),
        2 => (BondOrder::Double, None),
        3 => (BondOrder::Triple, None),
        4 => (BondOrder::Aromatic, None),
        9 => (BondOrder::Dative, None),
        0 => (BondOrder::Unspecified, None),
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
        other => (
            BondOrder::Unspecified,
            Some(QueryNode::predicate(BondQueryPredicate::MolFileQueryCode(
                other,
            ))),
        ),
    };

    let mut spec = BondSpec::new(begin, end, order);
    if matches!(bond_type, 5..=7) {
        spec = spec.with_prop("_MolFileBondQuery", "1");
    }
    if order == BondOrder::Aromatic {
        spec = spec.with_aromatic(true);
    }

    let mut molfile_stereo = None;
    if line.len() >= 12 && rdkit_substr(line, 9, 3) != "  0" {
        if let Ok(stereo) = parse_rdkit_unsigned(rdkit_substr(line, 9, 3), true) {
            molfile_stereo = Some(stereo);
            match stereo {
                0 => spec = spec.with_direction(BondDirection::None),
                1 => spec = spec.with_direction(BondDirection::BeginWedge),
                6 => spec = spec.with_direction(BondDirection::BeginDash),
                3 => {
                    spec = spec
                        .with_direction(BondDirection::EitherDouble)
                        .with_stereo(BondStereo::Any);
                }
                4 => spec = spec.with_direction(BondDirection::Unknown),
                _ => {}
            }
        }
    }

    if line.len() >= 18
        && rdkit_substr(line, 15, 3) != "  0"
        && let Ok(topology) = parse_rdkit_int(rdkit_substr(line, 15, 3), true)
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

    if let Some(query) = query {
        spec = spec.with_query(query);
    }

    Ok(V2000BondLine {
        line_number,
        text: line.to_string(),
        spec,
        molfile_bond_type: bond_type,
        molfile_stereo,
    })
}

fn parse_v2000_property_block<R: BufRead>(
    reader: &mut R,
    line_number: &mut usize,
    params: SdfReadParams,
    atoms: &mut [V2000AtomLine],
    bonds: &mut [V2000BondLine],
    state: &mut V2000PropertyState,
) -> Result<bool, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_property_block
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseMolBlockProperties
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:   // older mol files can have an atom list block here
    // RDKit✔️❗:   std::string tempStr = getLine(inStream);
    // RDKit✔️❗:   ++line;
    // RDKit✔️❗:   // there is apparently some software out there that puts a
    // RDKit✔️❗:   // blank line in mol blocks before the "M  END". If we aren't
    // RDKit✔️❗:   // doing strict parsing, deal with that here.
    // RDKit✔️❗:   if (!tempStr.size()) {
    // RDKit✔️❗:     if (!strictParsing) {
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:       ++line;
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Problems encountered parsing Mol data, unexpected blank line "
    // RDKit✔️❗:                 "found at line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   } else {
    // RDKit✔️❗:     if (tempStr[0] != 'M' && tempStr[0] != 'A' && tempStr[0] != 'V' &&
    // RDKit✔️❗:         tempStr[0] != 'G' && tempStr[0] != 'S') {
    // RDKit✔️❗:       ParseOldAtomList(mol, std::string_view(tempStr.c_str()), line);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   IDX_TO_SGROUP_MAP sGroupMap;
    // RDKit✔️❗:   IDX_TO_STR_VECT_MAP dataFieldsMap;
    // RDKit✔️❗:   bool fileComplete = false;
    // RDKit✔️❗:   bool firstChargeLine = true;
    // RDKit✔️❗:   unsigned int SCDcounter = 0;
    // RDKit✔️❗:   unsigned int lastDataSGroup = 0;
    // RDKit✔️❗:   std::ostringstream currentDataField;
    // RDKit✔️❗:   std::string lineBeg = tempStr.substr(0, 6);
    // RDKit✔️❗:   while (!inStream->eof() && !inStream->fail() && lineBeg != "M  END" &&
    // RDKit✔️❗:          tempStr.substr(0, 4) != "$$$$") {
    // RDKit✔️❗:     if (tempStr[0] == 'A') {
    // RDKit✔️❗:       line++;
    // RDKit✔️❗:       std::string nextLine = getLine(inStream);
    // RDKit✔️❗:       if (lineBeg != "M  END") {
    // RDKit✔️❗:         ParseAtomAlias(mol, tempStr, nextLine, line);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (tempStr[0] == 'G') {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:           << " deprecated group abbreviation ignored on line " << line
    // RDKit✔️❗:           << std::endl;
    // RDKit✔️❗:       // we need to skip the next line, which holds the abbreviation:
    // RDKit✔️❗:       line++;
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:     } else if (tempStr[0] == 'V') {
    // RDKit✔️❗:       ParseAtomValue(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "S  SKP") {
    // RDKit✔️❗:       int nToSkip = FileParserUtils::toInt(tempStr.substr(6, 3));
    // RDKit✔️❗:       if (nToSkip < 0) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "negative skip value " << nToSkip << " on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       for (unsigned int i = 0; i < static_cast<unsigned int>(nToSkip); ++i) {
    // RDKit✔️❗:         ++line;
    // RDKit✔️❗:         tempStr = getLine(inStream);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (lineBeg == "M  ALS") {
    // RDKit✔️❗:       ParseNewAtomList(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  ISO") {
    // RDKit✔️❗:       ParseIsotopeLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  RGP") {
    // RDKit✔️❗:       ParseRGroupLabels(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  RBC") {
    // RDKit✔️❗:       ParseRingBondCountLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  SUB") {
    // RDKit✔️❗:       ParseSubstitutionCountLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  UNS") {
    // RDKit✔️❗:       ParseUnsaturationLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  CHG") {
    // RDKit✔️❗:       ParseChargeLine(mol, tempStr, firstChargeLine, line);
    // RDKit✔️❗:       firstChargeLine = false;
    // RDKit✔️❗:     } else if (lineBeg == "M  RAD") {
    // RDKit✔️❗:       ParseRadicalLine(mol, tempStr, firstChargeLine, line);
    // RDKit✔️❗:       firstChargeLine = false;
    // RDKit✔️❗:     } else if (lineBeg == "M  PXA") {
    // RDKit✔️❗:       ParsePXALine(mol, tempStr, line);
    // RDKit✔️❗:
    // RDKit✔️❗:       /* SGroup parsing start */
    // RDKit✔️❗:     } else if (lineBeg == "M  STY") {
    // RDKit✔️❗:       ParseSGroupV2000STYLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SST") {
    // RDKit✔️❗:       ParseSGroupV2000SSTLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SLB") {
    // RDKit✔️❗:       ParseSGroupV2000SLBLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SCN") {
    // RDKit✔️❗:       ParseSGroupV2000SCNLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SDS") {
    // RDKit✔️❗:       ParseSGroupV2000SDSLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SAL" || lineBeg == "M  SBL" ||
    // RDKit✔️❗:                lineBeg == "M  SPA") {
    // RDKit✔️❗:       ParseSGroupV2000VectorDataLine(sGroupMap, mol, tempStr, line,
    // RDKit✔️❗:                                      strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SMT") {
    // RDKit✔️❗:       ParseSGroupV2000SMTLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SDI") {
    // RDKit✔️❗:       ParseSGroupV2000SDILine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  CRS") {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Unsupported SGroup subtype '" << lineBeg << "' on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } else if (lineBeg == "M  SBV") {
    // RDKit✔️❗:       ParseSGroupV2000SBVLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SDT") {
    // RDKit✔️❗:       ParseSGroupV2000SDTLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SDD") {
    // RDKit✔️❗:       ParseSGroupV2000SDDLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SCD" || lineBeg == "M  SED") {
    // RDKit✔️❗:       ParseSGroupV2000SCDSEDLine(sGroupMap, dataFieldsMap, mol, tempStr, line,
    // RDKit✔️❗:                                  strictParsing, SCDcounter, lastDataSGroup,
    // RDKit✔️❗:                                  currentDataField);
    // RDKit✔️❗:     } else if (lineBeg == "M  SPL") {
    // RDKit✔️❗:       ParseSGroupV2000SPLLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SNC") {
    // RDKit✔️❗:       ParseSGroupV2000SNCLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SAP") {
    // RDKit✔️❗:       ParseSGroupV2000SAPLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SCL") {
    // RDKit✔️❗:       ParseSGroupV2000SCLLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SBT") {
    // RDKit✔️❗:       ParseSGroupV2000SBTLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:
    // RDKit✔️❗:       /* SGroup parsing end */
    // RDKit✔️❗:     } else if (lineBeg == "M  ZBO") {
    // RDKit✔️❗:       ParseZBOLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  ZCH") {
    // RDKit✔️❗:       ParseZCHLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  HYD") {
    // RDKit✔️❗:       ParseHYDLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  MRV") {
    // RDKit✔️❗:       ParseMarvinSmartsLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  APO") {
    // RDKit✔️❗:       ParseAttachPointLine(mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  LIN") {
    // RDKit✔️❗:       ParseLinkNodeLine(mol, tempStr, line);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     line++;
    // RDKit✔️❗:     tempStr = getLine(inStream);
    // RDKit✔️❗:     lineBeg = tempStr.substr(0, 6);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (tempStr[0] == 'M' && tempStr.substr(0, 6) == "M  END") {
    // RDKit✔️❗:     // All went well, make final updates to SGroups, and add them to Mol
    // RDKit✔️❗:     for (auto &sgroup : sGroupMap) {
    // RDKit✔️❗:       if (sgroup.second.getIsValid()) {
    // RDKit✔️❗:         sgroup.second.setProp("DATAFIELDS", dataFieldsMap[sgroup.first]);
    // RDKit✔️❗:         sgroup.second.setIsValid(checkAttachmentPointsAreValid(mol, sgroup));
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (sgroup.second.getIsValid()) {
    // RDKit✔️❗:         addSubstanceGroup(*mol, sgroup.second);
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "SGroup " << sgroup.first << " is invalid";
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:               << errout.str() << " and will be ignored" << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     fileComplete = true;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return fileComplete;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseMolBlockProperties
    // END RDKIT CPP BODY: parse_v2000_property_block

    let Some(mut line) = read_rdkit_line(reader)? else {
        return Err(SdfReadError::Parse(
            "EOF hit while reading Mol data properties".to_string(),
        ));
    };
    *line_number += 1;

    if line.is_empty() {
        if params.strict_parsing {
            return Err(SdfReadError::Parse(format!(
                "Problems encountered parsing Mol data, unexpected blank line found at line {}",
                *line_number
            )));
        }
        let Some(next_line) = read_rdkit_line(reader)? else {
            return Err(SdfReadError::Parse(
                "EOF hit while reading Mol data properties".to_string(),
            ));
        };
        line = next_line;
        *line_number += 1;
    }

    loop {
        let line_beg = rdkit_substr(&line, 0, 6);
        if line.starts_with('M') && line_beg == "M  END" {
            return Ok(true);
        }
        if rdkit_substr(&line, 0, 4) == "$$$$" {
            return Ok(false);
        }

        if line.starts_with('A') {
            *line_number += 1;
            let next_line = read_rdkit_line(reader)?.ok_or_else(|| {
                SdfReadError::Parse("EOF hit while reading atom alias line".to_string())
            })?;
            parse_atom_alias_line(&line, &next_line, *line_number, atoms)?;
        } else if line.starts_with('G') {
            *line_number += 1;
            let _ = read_rdkit_line(reader)?.ok_or_else(|| {
                SdfReadError::Parse(
                    "EOF hit while skipping deprecated group abbreviation".to_string(),
                )
            })?;
        } else if line.starts_with('V') {
            parse_atom_value_line(&line, *line_number, atoms)?;
        } else if line_beg == "S  SKP" {
            skip_sgroup_lines(reader, &line, line_number)?;
            let Some(next_line) = read_rdkit_line(reader)? else {
                return Err(SdfReadError::Parse(
                    "EOF hit while reading Mol data properties".to_string(),
                ));
            };
            line = next_line;
            *line_number += 1;
            continue;
        } else if line.starts_with('M') {
            match line_beg {
                "M  CHG" => {
                    parse_v2000_charge_line(&line, *line_number, atoms, state.first_charge_line)?;
                    state.first_charge_line = false;
                }
                "M  RAD" => {
                    parse_v2000_radical_line(&line, *line_number, atoms, state.first_charge_line)?;
                    state.first_charge_line = false;
                }
                "M  ISO" => parse_v2000_isotope_line(&line, *line_number, atoms)?,
                "M  RGP" => parse_rgroup_labels_line(&line, *line_number, atoms)?,
                "M  RBC" => parse_ring_bond_count_line(&line, *line_number, atoms, state)?,
                "M  SUB" => parse_substitution_count_line(&line, *line_number, atoms)?,
                "M  UNS" => parse_unsaturation_line(&line, *line_number, atoms)?,
                "M  PXA" => parse_pxa_line(&line, *line_number, atoms)?,
                "M  ALS" => parse_v2000_atom_list_line(&line, *line_number, params, atoms)?,
                "M  ZBO" => parse_zbo_line(&line, *line_number, bonds)?,
                "M  ZCH" => parse_zch_line(&line, *line_number, atoms)?,
                "M  HYD" => parse_hyd_line(&line, *line_number, atoms)?,
                "M  MRV" => parse_marvin_smarts_line(&line, *line_number, atoms)?,
                "M  APO" => parse_attach_point_line(&line, *line_number, params, atoms)?,
                "M  LIN" => parse_link_node_line(&line, *line_number, atoms.len(), state)?,
                "M  STY" | "M  SST" | "M  SLB" | "M  SCN" | "M  SDS" | "M  SAL" | "M  SBL"
                | "M  SPA" | "M  SMT" | "M  SDI" | "M  SBV" | "M  SDT" | "M  SDD" | "M  SCD"
                | "M  SED" | "M  SPL" | "M  SNC" | "M  SAP" | "M  SCL" | "M  SBT" => {
                    parse_v2000_sgroup_line(&line, *line_number, params, state)?
                }
                _ => {}
            }
        } else if !line.starts_with(['S']) {
            parse_v2000_atom_list_line(&line, *line_number, params, atoms)?;
        }

        let Some(next_line) = read_rdkit_line(reader)? else {
            return Ok(false);
        };
        line = next_line;
        *line_number += 1;
    }
}

fn parse_v2000_charge_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
    first_call: bool,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_charge_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseChargeLine
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == std::string("M  CHG"), "bad charge line");
    // RDKit✔️❗:
    // RDKit✔️❗:   // if this line is specified all the atom other than those specified
    // RDKit✔️❗:   // here should carry a charge of 0; but we should only do this once:
    // RDKit✔️❗:   if (firstCall) {
    // RDKit✔️❗:     for (ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms();
    // RDKit✔️❗:          ++ai) {
    // RDKit✔️❗:       (*ai)->setFormalCharge(0);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   int ie, nent;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nent = FileParserUtils::toInt(text.substr(6, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   int spos = 9;
    // RDKit✔️❗:   for (ie = 0; ie < nent; ie++) {
    // RDKit✔️❗:     int aid, chg;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       aid = FileParserUtils::toInt(text.substr(spos, 4));
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       chg = FileParserUtils::toInt(text.substr(spos, 4));
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       mol->getAtomWithIdx(aid - 1)->setFormalCharge(chg);
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(spos, 4)
    // RDKit✔️❗:              << "' to int on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseChargeLine
    // END RDKIT CPP BODY: parse_v2000_charge_line

    if first_call {
        for atom in atoms.iter_mut() {
            atom.spec = atom.spec.clone().with_formal_charge(0);
        }
    }

    let nent = parse_required_int_field(line, 6, 3, line_number)?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_required_int_field(line, spos, 4, line_number)?;
        spos += 4;
        let charge = parse_required_int_field(line, spos, 4, line_number)?;
        spos += 4;
        let atom = atom_line_mut(atoms, aid as u32, line_number)?;
        atom.spec = atom.spec.clone().with_formal_charge(charge as i8);
    }
    Ok(())
}

fn parse_v2000_radical_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
    first_call: bool,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_radical_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseRadicalLine
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == std::string("M  RAD"), "bad charge line");
    // RDKit✔️❗:
    // RDKit✔️❗:   // if this line is specified all the atom other than those specified
    // RDKit✔️❗:   // here should carry a charge of 0; but we should only do this once:
    // RDKit✔️❗:   if (firstCall) {
    // RDKit✔️❗:     for (ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms();
    // RDKit✔️❗:          ++ai) {
    // RDKit✔️❗:       (*ai)->setFormalCharge(0);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   int ie, nent;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nent = FileParserUtils::toInt(text.substr(6, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   int spos = 9;
    // RDKit✔️❗:   for (ie = 0; ie < nent; ie++) {
    // RDKit✔️❗:     int aid, rad;
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       aid = FileParserUtils::toInt(text.substr(spos, 4));
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       rad = FileParserUtils::toInt(text.substr(spos, 4));
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:
    // RDKit✔️❗:       switch (rad) {
    // RDKit✔️❗:         case 0:
    // RDKit✔️❗:           // This shouldn't be required, but let's make sure.
    // RDKit✔️❗:           mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(0);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 1:
    // RDKit✔️❗:           mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(2);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 2:
    // RDKit✔️❗:           mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(1);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 3:
    // RDKit✔️❗:           mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(2);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         default:
    // RDKit✔️❗:           errout << "Unrecognized radical value " << rad << " for atom "
    // RDKit✔️❗:                  << aid - 1 << " on line " << line << std::endl;
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(spos, 4)
    // RDKit✔️❗:              << "' to int on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseRadicalLine
    // END RDKIT CPP BODY: parse_v2000_radical_line

    if first_call {
        for atom in atoms.iter_mut() {
            atom.spec = atom.spec.clone().with_formal_charge(0);
        }
    }

    let nent = parse_required_int_field(line, 6, 3, line_number)?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_required_int_field(line, spos, 4, line_number)?;
        spos += 4;
        let radical_code = parse_required_int_field(line, spos, 4, line_number)?;
        spos += 4;
        let radical_electrons = match radical_code {
            0 => 0,
            1 => 2,
            2 => 1,
            3 => 2,
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Unrecognized radical value {radical_code} for atom {} on line {line_number}",
                    aid - 1
                )));
            }
        };
        let atom = atom_line_mut(atoms, aid as u32, line_number)?;
        atom.spec = atom.spec.clone().with_radical_electrons(radical_electrons);
    }
    Ok(())
}

fn parse_v2000_isotope_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_isotope_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseIsotopeLine
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == std::string("M  ISO"), "bad isotope line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nent;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   unsigned int spos = 9;
    // RDKit✔️❗:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️❗:     unsigned int aid;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️❗:           text.substr(spos, 4));
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       Atom *atom = mol->getAtomWithIdx(aid - 1);
    // RDKit✔️❗:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️❗:         int isotope = FileParserUtils::toInt(text.substr(spos, 4));
    // RDKit✔️❗:         if (isotope < 0) {
    // RDKit✔️❗:           BOOST_LOG(rdErrorLog)
    // RDKit✔️❗:               << " atom " << aid
    // RDKit✔️❗:               << " has a negative isotope value. line:  " << line << std::endl;
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           atom->setIsotope(isotope);
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(spos, 4)
    // RDKit✔️❗:              << "' to int on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseIsotopeLine
    // END RDKIT CPP BODY: parse_v2000_isotope_line

    let nent = parse_required_unsigned_field(line, 6, 3, line_number)?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_required_unsigned_field(line, spos, 4, line_number)?;
        spos += 4;
        let atom = atom_line_mut(atoms, aid, line_number)?;
        if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            let isotope = parse_required_int_field(line, spos, 4, line_number)?;
            if isotope >= 0 {
                atom.spec = atom.spec.clone().with_isotope(isotope as u16);
            }
        }
        spos += 4;
    }
    Ok(())
}

fn parse_rgroup_labels_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_rgroup_labels_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseRGroupLabels
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == std::string("M  RGP"),
    // RDKit✔️❗:                "bad R group label line");
    // RDKit✔️❗:
    // RDKit✔️❗:   int nLabels;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nLabels = FileParserUtils::toInt(text.substr(6, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   for (int i = 0; i < nLabels; i++) {
    // RDKit✔️❗:     int pos = 10 + i * 8;
    // RDKit✔️❗:     unsigned int atIdx;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       atIdx = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️❗:           text.substr(pos, 3));
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(pos, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     unsigned int rLabel;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       rLabel = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️❗:           text.substr(pos + 4, 3));
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(pos + 4, 3)
    // RDKit✔️❗:              << "' to int on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     atIdx -= 1;
    // RDKit✔️❗:     if (atIdx > mol->getNumAtoms()) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Attempt to set R group label on nonexistent atom " << atIdx
    // RDKit✔️❗:              << " on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     QueryAtom qatom(*(mol->getAtomWithIdx(atIdx)));
    // RDKit✔️❗:     qatom.setProp(common_properties::_MolFileRLabel, rLabel);
    // RDKit✔️❗:
    // RDKit✔️❗:     // set the dummy label so that this is shown correctly
    // RDKit✔️❗:     // in other pieces of the code :
    // RDKit✔️❗:     // (this was sf.net issue 3316600)
    // RDKit✔️❗:     std::string dLabel = "R" + std::to_string(rLabel);
    // RDKit✔️❗:     qatom.setProp(common_properties::dummyLabel, dLabel);
    // RDKit✔️❗:
    // RDKit✔️❗:     // the CTFile spec (June 2005 version) technically only allows
    // RDKit✔️❗:     // R labels up to 32. Since there are three digits, we'll accept
    // RDKit✔️❗:     // anything: so long as it's positive and less than 1000:
    // RDKit✔️❗:     if (rLabel > 0 && rLabel < 999) {
    // RDKit✔️❗:       qatom.setIsotope(rLabel);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     qatom.setQuery(makeAtomNullQuery());
    // RDKit✔️❗:     mol->replaceAtom(atIdx, &qatom);
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseRGroupLabels
    // END RDKIT CPP BODY: parse_rgroup_labels_line

    let n_labels = parse_required_int_field(line, 6, 3, line_number)?;
    for i in 0..n_labels {
        let pos = 10 + i as usize * 8;
        let atom_index = parse_required_unsigned_field(line, pos, 3, line_number)?;
        let r_label = parse_required_unsigned_field(line, pos + 4, 3, line_number)?;
        let atom = atom_line_mut(atoms, atom_index, line_number)?;
        let mut spec = atom.spec.clone();
        spec = spec
            .with_prop("_MolFileRLabel", r_label.to_string())
            .with_prop("dummyLabel", format!("R{r_label}"))
            .with_query(QueryNode::predicate(AtomQueryPredicate::RGroupLabel(
                r_label,
            )));
        if (1..999).contains(&r_label) {
            spec = spec.with_isotope(r_label as u16);
        }
        atom.spec = spec;
    }
    Ok(())
}

fn parse_substitution_count_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    let nent = parse_required_unsigned_field(line, 6, 3, line_number)?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_required_unsigned_field(line, spos, 4, line_number)?;
        spos += 4;
        let count = if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            parse_required_int_field(line, spos, 4, line_number)?
        } else {
            0
        };
        spos += 4;
        if count == 0 {
            continue;
        }
        let atom = atom_line_mut(atoms, aid, line_number)?;
        atom.spec = atom
            .spec
            .clone()
            .with_prop("molSubstCount", count.to_string());
    }
    Ok(())
}

fn parse_unsaturation_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    let nent = parse_required_unsigned_field(line, 6, 3, line_number)?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_required_unsigned_field(line, spos, 4, line_number)?;
        spos += 4;
        let count = if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            parse_required_int_field(line, spos, 4, line_number)?
        } else {
            0
        };
        spos += 4;
        if count == 0 {
            continue;
        }
        if count != 1 {
            return Err(SdfReadError::Parse(format!(
                "Value {count} is not supported as an unsaturation query (only 0 and 1 are allowed). line: {line_number}"
            )));
        }
        let atom = atom_line_mut(atoms, aid, line_number)?;
        atom.spec = atom.spec.clone().with_query(merge_atom_query(
            atom.spec.query().cloned(),
            QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
        ));
    }
    Ok(())
}

fn parse_ring_bond_count_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    let nent = parse_required_unsigned_field(line, 6, 3, line_number)?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_required_unsigned_field(line, spos, 4, line_number)?;
        spos += 4;
        let count = if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            parse_required_int_field(line, spos, 4, line_number)?
        } else {
            0
        };
        spos += 4;
        if count == 0 {
            continue;
        }
        let atom = atom_line_mut(atoms, aid, line_number)?;
        let predicate = match count {
            -1 => AtomQueryPredicate::RingBondCount(0),
            -2 => {
                state
                    .molecule_props
                    .insert("_NeedsQueryScan".to_string(), "1".to_string());
                AtomQueryPredicate::RingBondCount(crate::search::query::QUERY_SCAN_MAGIC_VALUE)
            }
            1..=3 => AtomQueryPredicate::RingBondCount(count as u32),
            4 => AtomQueryPredicate::RingBondCountLessEqual(4),
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Value {count} is not supported as a ring-bond count query. line: {line_number}"
                )));
            }
        };
        atom.spec = atom
            .spec
            .clone()
            .with_prop("molRingBondCount", count.to_string())
            .with_query(merge_atom_query(
                atom.spec.query().cloned(),
                QueryNode::predicate(predicate),
            ));
    }
    Ok(())
}

fn parse_pxa_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    let mut pos = 7_usize;
    let atom_index = parse_required_unsigned_field(line, pos, 3, line_number)?;
    pos += 3;
    let atom = atom_line_mut(atoms, atom_index, line_number)?;
    atom.spec = atom
        .spec
        .clone()
        .with_prop("_MolFile_PXA", line[pos..].to_string());
    Ok(())
}

fn parse_marvin_smarts_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_marvin_smarts_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseMarvinSmartsLine
    // RDKit✔️✔️:   if (text.substr(0, 10) != "M  MRV SMA") {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int idx;
    // RDKit✔️✔️:   std::string idxTxt = text.substr(atomNumStart, smartsStart - atomNumStart);
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     idx = FileParserUtils::stripSpacesAndCast<unsigned int>(idxTxt) - 1;
    // RDKit✔️✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Cannot convert '" << idxTxt << "' to an atom index on line "
    // RDKit✔️✔️:            << line;
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   URANGE_CHECK(idx, mol->getNumAtoms());
    // RDKit✔️✔️:   std::string sma = text.substr(smartsStart);
    // RDKit✔️✔️:   Atom *at = mol->getAtomWithIdx(idx);
    // RDKit✔️✔️:   at->setProp(common_properties::MRV_SMA, sma);
    // RDKit✔️✔️:   RWMol *m = nullptr;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     m = SmartsToMol(sma);
    // RDKit✔️✔️:   } catch (...) {
    // RDKit✔️✔️:     // Is this ever used?
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (m) {
    // RDKit✔️✔️:     QueryAtom::QUERYATOM_QUERY *query = new RecursiveStructureQuery(m);
    // RDKit✔️✔️:     if (!at->hasQuery()) {
    // RDKit✔️✔️:       QueryAtom qAt(*at);
    // RDKit✔️✔️:       int oidx = at->getIdx();
    // RDKit✔️✔️:       mol->replaceAtom(oidx, &qAt);
    // RDKit✔️✔️:       at = mol->getAtomWithIdx(oidx);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     at->expandQuery(query, Queries::COMPOSITE_AND);
    // RDKit✔️✔️:     at->setProp(common_properties::_MolFileAtomQuery, 1);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Cannot parse smarts: '" << sma << "' on line " << line;
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseMarvinSmartsLine
    // END RDKIT CPP BODY: parse_marvin_smarts_line

    if rdkit_substr(line, 0, 10) != "M  MRV SMA" {
        return Ok(());
    }
    let atom_num_start = 10_usize;
    let smarts_start = 15_usize;
    let idx_text = rdkit_substr(line, atom_num_start, smarts_start - atom_num_start);
    let atom_index = parse_rdkit_unsigned(idx_text, true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{idx_text}' to an atom index on line {line_number}"
        ))
    })?;
    let sma = line.get(smarts_start..).unwrap_or_default();
    let atom = atom_line_mut(atoms, atom_index, line_number)?;
    // Local complexity review: the Marvin SMARTS is compiled exactly once and
    // the resulting canonical query graph is moved into the recursive query.
    let query_molecule = crate::mol_from_smarts(sma, &crate::SmartsParseParams::default())
        .map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot parse smarts: '{sma}' on line {line_number}"
            ))
        })?;
    atom.spec = atom
        .spec
        .clone()
        .with_prop("MRV_SMA", sma)
        .with_prop("_MolFileAtomQuery", "1")
        .with_query(merge_atom_query(
            atom.spec.query().cloned(),
            QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(
                crate::search::query::RecursiveStructureQuery::from_smarts(
                    sma,
                    query_molecule.into(),
                    0,
                ),
            )),
        ));
    Ok(())
}

fn parse_link_node_line(
    line: &str,
    line_number: usize,
    atom_count: usize,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    let nent = parse_required_unsigned_field(line, 6, 3, line_number)?;
    let mut prop_val = String::new();
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_required_unsigned_field(line, spos, 4, line_number)?;
        if aid == 0 || aid as usize > atom_count {
            return Err(SdfReadError::Parse(format!(
                "LIN specification has bad atom idx on line {line_number}"
            )));
        }
        spos += 4;
        if line.len() < spos + 4 || rdkit_substr(line, spos, 4) == "    " {
            return Err(SdfReadError::Parse(format!(
                "LIN specification missing repeat count on line {line_number}"
            )));
        }
        let repeat_count = parse_required_unsigned_field(line, spos, 4, line_number)?;
        spos += 4;
        if repeat_count < 2 {
            return Err(SdfReadError::Parse(format!(
                "LIN specification: repeat count must be >=2 on line {line_number}"
            )));
        }
        let subst_b = if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            parse_required_unsigned_field(line, spos, 4, line_number)?
        } else {
            0
        };
        spos += 4;
        let subst_c = if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            parse_required_unsigned_field(line, spos, 4, line_number)?
        } else {
            0
        };
        spos += 4;
        if subst_b == 0 || subst_b as usize > atom_count || subst_c as usize > atom_count {
            return Err(SdfReadError::Parse(format!(
                "LIN specification has bad substituent idx on line {line_number}"
            )));
        }
        if !prop_val.is_empty() {
            prop_val.push('|');
        }
        if subst_c != 0 {
            prop_val.push_str(&format!(
                "1 {repeat_count} 2 {aid} {subst_b} {aid} {subst_c}"
            ));
        } else {
            prop_val.push_str(&format!("1 {repeat_count} 1 {aid} {subst_b}"));
        }
    }
    state
        .molecule_props
        .insert("_MolFileLinkNodes".to_string(), prop_val);
    Ok(())
}

fn parse_v2000_atom_list_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_atom_list_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseOldAtomList
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   unsigned int idx;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(0, 3)) -
    // RDKit✔️❗:           1;
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(0, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   URANGE_CHECK(idx, mol->getNumAtoms());
    // RDKit✔️❗:   QueryAtom a(*(mol->getAtomWithIdx(idx)));
    // RDKit✔️❗:
    // RDKit✔️❗:   auto *q = new ATOM_OR_QUERY;
    // RDKit✔️❗:   q->setDescription("AtomOr");
    // RDKit✔️❗:
    // RDKit✔️❗:   switch (text[4]) {
    // RDKit✔️❗:     case 'T':
    // RDKit✔️❗:       q->setNegation(true);
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     case 'F':
    // RDKit✔️❗:       q->setNegation(false);
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     default:
    // RDKit✔️❗:       delete q;
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Unrecognized atom-list query modifier: '" << text[4]
    // RDKit✔️❗:              << "' on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   int nQueries;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nQueries = FileParserUtils::toInt(text.substr(9, 1));
    // RDKit✔️❗:   } catch (const std::out_of_range &) {
    // RDKit✔️❗:     delete q;
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert position 9 of '" << text << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     delete q;
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(9, 1) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   RANGE_CHECK(0, nQueries, 5);
    // RDKit✔️❗:   for (int i = 0; i < nQueries; i++) {
    // RDKit✔️❗:     int pos = 11 + i * 4;
    // RDKit✔️❗:     int atNum;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       atNum = FileParserUtils::toInt(text.substr(pos, 3));
    // RDKit✔️❗:     } catch (const std::out_of_range &) {
    // RDKit✔️❗:       delete q;
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert position " << pos << " of '" << text
    // RDKit✔️❗:              << "' to int on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       delete q;
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(pos, 3) << "' to int on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     RANGE_CHECK(0, atNum, 200);  // goofy!
    // RDKit✔️❗:     q->addChild(
    // RDKit✔️❗:         QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(atNum)));
    // RDKit✔️❗:     if (!i) {
    // RDKit✔️❗:       a.setAtomicNum(atNum);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   a.setQuery(q);
    // RDKit✔️❗:   a.setProp(common_properties::_MolFileAtomQuery, 1);
    // RDKit✔️❗:
    // RDKit✔️❗:   mol->replaceAtom(idx, &a);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseOldAtomList
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseNewAtomList
    // RDKit✔️❗:
    // RDKit✔️❗:   if (text.size() < 15) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Atom list line too short: '" << text << "'";
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == std::string("M  ALS"),
    // RDKit✔️❗:                "bad atom list line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int idx;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(7, 3)) -
    // RDKit✔️❗:           1;
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(7, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   URANGE_CHECK(idx, mol->getNumAtoms());
    // RDKit✔️❗:
    // RDKit✔️❗:   int nQueries;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nQueries = FileParserUtils::toInt(text.substr(10, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(10, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (!nQueries) {
    // RDKit✔️❗:     BOOST_LOG(rdWarningLog) << "Empty atom list: '" << text << "' on line "
    // RDKit✔️❗:                             << line << "." << std::endl;
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (nQueries < 0) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "negative length atom list: '" << text << "' on line " << line
    // RDKit✔️❗:            << "." << std::endl;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   QueryAtom *a = nullptr;
    // RDKit✔️❗:   QueryAtom *qaOrig = nullptr;
    // RDKit✔️❗:   QueryAtom::QUERYATOM_QUERY *qOrig = nullptr;
    // RDKit✔️❗:   Atom *aOrig = mol->getAtomWithIdx(idx);
    // RDKit✔️❗:   for (unsigned int i = 0; i < static_cast<unsigned int>(nQueries); i++) {
    // RDKit✔️❗:     unsigned int pos = 16 + i * 4;
    // RDKit✔️❗:     if (text.size() < pos + 4) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Atom list line too short: '" << text << "' on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     std::string atSymb = text.substr(pos, 4);
    // RDKit✔️❗:     atSymb.erase(atSymb.find(' '), atSymb.size());
    // RDKit✔️❗:     int atNum = PeriodicTable::getTable()->getAtomicNumber(atSymb);
    // RDKit✔️❗:     if (!i) {
    // RDKit✔️❗:       if (aOrig->hasQuery()) {
    // RDKit✔️❗:         qaOrig = dynamic_cast<QueryAtom *>(aOrig);
    // RDKit✔️❗:         if (qaOrig) {
    // RDKit✔️❗:           qOrig = qaOrig->getQuery();
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       a = new QueryAtom(*aOrig);
    // RDKit✔️❗:       a->setAtomicNum(atNum);
    // RDKit✔️❗:       if (!qOrig) {
    // RDKit✔️❗:         qOrig = a->getQuery()->copy();
    // RDKit✔️❗:       }
    // RDKit✔️❗:       a->setQuery(makeAtomNumQuery(atNum));
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       a->expandQuery(makeAtomNumQuery(atNum), Queries::COMPOSITE_OR, true);
    // RDKit✔️❗:       // For COMPOSITE_OR query atoms, reset atomic num to 0 such that they are
    // RDKit✔️❗:       // exported as "*" in SMILES
    // RDKit✔️❗:       a->setAtomicNum(0);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   ASSERT_INVARIANT(a, "no atom built");
    // RDKit✔️❗:   if (qOrig) {
    // RDKit✔️❗:     std::vector<const QueryAtom::QUERYATOM_QUERY *> queryVect;
    // RDKit✔️❗:     if (getAndQueries(qOrig, queryVect)) {
    // RDKit✔️❗:       for (const auto &q : queryVect) {
    // RDKit✔️❗:         if (q->getDescription() != "AtomAtomicNum") {
    // RDKit✔️❗:           a->expandQuery(q->copy(), Queries::COMPOSITE_AND, true);
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:     if (!qaOrig) {
    // RDKit✔️❗:       delete qOrig;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   a->setProp(common_properties::_MolFileAtomQuery, 1);
    // RDKit✔️❗:   switch (text[14]) {
    // RDKit✔️❗:     case 'T':
    // RDKit✔️❗:       a->getQuery()->setNegation(true);
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     case 'F':
    // RDKit✔️❗:       a->getQuery()->setNegation(false);
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     default:
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Unrecognized atom-list query modifier: '" << text[14]
    // RDKit✔️❗:              << "' on line " << line;
    // RDKit✔️❗:       delete a;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   mol->replaceAtom(idx, a);
    // RDKit✔️❗:   delete a;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseNewAtomList
    // END RDKIT CPP BODY: parse_v2000_atom_list_line

    if line.starts_with("M  ALS") {
        parse_v2000_new_atom_list_line(line, line_number, params, atoms)
    } else {
        parse_v2000_old_atom_list_line(line, line_number, params, atoms)
    }
}

fn parse_v2000_new_atom_list_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    if line.len() < 15 {
        return Err(SdfReadError::Parse(format!(
            "Atom list line too short: '{line}'"
        )));
    }

    let atom_index = parse_required_unsigned_field(line, 7, 3, line_number)?;
    let n_queries = parse_required_int_field(line, 10, 3, line_number)?;
    if n_queries == 0 {
        return Ok(());
    }
    if n_queries < 0 {
        return Err(SdfReadError::Parse(format!(
            "negative length atom list: '{line}' on line {line_number}."
        )));
    }

    let mut atomic_numbers = Vec::with_capacity(n_queries as usize);
    for i in 0..n_queries as usize {
        let pos = 16 + i * 4;
        if line.len() < pos + 4 {
            return Err(SdfReadError::Parse(format!(
                "Atom list line too short: '{line}' on line {line_number}"
            )));
        }
        let symbol = rdkit_substr(line, pos, 4)
            .split_once(' ')
            .map_or(rdkit_substr(line, pos, 4), |(head, _)| head);
        atomic_numbers.push(atomic_number_from_mol_symbol(
            symbol,
            params.strict_parsing,
        )?);
    }

    let atom = atom_line_mut(atoms, atom_index, line_number)?;
    let query = match line.as_bytes().get(14) {
        Some(b'T') => QueryNode::predicate(AtomQueryPredicate::AtomicNumberNotIn(atomic_numbers)),
        Some(b'F') => QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(atomic_numbers)),
        Some(other) => {
            return Err(SdfReadError::Parse(format!(
                "Unrecognized atom-list query modifier: '{}' on line {line_number}",
                *other as char
            )));
        }
        None => {
            return Err(SdfReadError::Parse(format!(
                "Atom list line too short: '{line}'"
            )));
        }
    };
    let first_atomic_number = match &query {
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberIn(numbers)) if numbers.len() == 1 => {
            numbers[0]
        }
        _ => 0,
    };
    atom.spec = atom
        .spec
        .clone()
        .with_element(element_from_query_atomic_number(
            first_atomic_number,
            line_number,
        )?)
        .with_query(merge_atom_query(atom.spec.query().cloned(), query))
        .with_prop("_MolFileAtomQuery", "1");
    Ok(())
}

fn parse_v2000_old_atom_list_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    let atom_index = parse_required_unsigned_field(line, 0, 3, line_number)?;
    let n_queries = parse_required_int_field(line, 9, 1, line_number)?;
    if !(0..=5).contains(&n_queries) {
        return Err(SdfReadError::Parse(format!(
            "atom-list query count {n_queries} out of range on line {line_number}"
        )));
    }
    let mut atomic_numbers = Vec::with_capacity(n_queries as usize);
    for i in 0..n_queries as usize {
        let pos = 11 + i * 4;
        let atomic_number = parse_required_int_field(line, pos, 3, line_number)?;
        if !(0..=200).contains(&atomic_number) {
            return Err(SdfReadError::Parse(format!(
                "atom-list atomic number {atomic_number} out of range on line {line_number}"
            )));
        }
        let _ = params;
        atomic_numbers.push(atomic_number as u8);
    }

    let atom = atom_line_mut(atoms, atom_index, line_number)?;
    let query = match line.as_bytes().get(4) {
        Some(b'T') => QueryNode::predicate(AtomQueryPredicate::AtomicNumberNotIn(atomic_numbers)),
        Some(b'F') => QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(atomic_numbers)),
        Some(other) => {
            return Err(SdfReadError::Parse(format!(
                "Unrecognized atom-list query modifier: '{}' on line {line_number}",
                *other as char
            )));
        }
        None => {
            return Err(SdfReadError::Parse(format!(
                "Atom list line too short: '{line}'"
            )));
        }
    };
    let first_atomic_number = match &query {
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberIn(numbers))
            if !numbers.is_empty() =>
        {
            numbers[0]
        }
        _ => 0,
    };
    atom.spec = atom
        .spec
        .clone()
        .with_element(element_from_query_atomic_number(
            first_atomic_number,
            line_number,
        )?)
        .with_query(query)
        .with_prop("_MolFileAtomQuery", "1");
    Ok(())
}

fn parse_v2000_sgroup_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseMolBlockProperties
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:   // older mol files can have an atom list block here
    // RDKit✔️❗:   std::string tempStr = getLine(inStream);
    // RDKit✔️❗:   ++line;
    // RDKit✔️❗:   // there is apparently some software out there that puts a
    // RDKit✔️❗:   // blank line in mol blocks before the "M  END". If we aren't
    // RDKit✔️❗:   // doing strict parsing, deal with that here.
    // RDKit✔️❗:   if (!tempStr.size()) {
    // RDKit✔️❗:     if (!strictParsing) {
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:       ++line;
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Problems encountered parsing Mol data, unexpected blank line "
    // RDKit✔️❗:                 "found at line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   } else {
    // RDKit✔️❗:     if (tempStr[0] != 'M' && tempStr[0] != 'A' && tempStr[0] != 'V' &&
    // RDKit✔️❗:         tempStr[0] != 'G' && tempStr[0] != 'S') {
    // RDKit✔️❗:       ParseOldAtomList(mol, std::string_view(tempStr.c_str()), line);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   IDX_TO_SGROUP_MAP sGroupMap;
    // RDKit✔️❗:   IDX_TO_STR_VECT_MAP dataFieldsMap;
    // RDKit✔️❗:   bool fileComplete = false;
    // RDKit✔️❗:   bool firstChargeLine = true;
    // RDKit✔️❗:   unsigned int SCDcounter = 0;
    // RDKit✔️❗:   unsigned int lastDataSGroup = 0;
    // RDKit✔️❗:   std::ostringstream currentDataField;
    // RDKit✔️❗:   std::string lineBeg = tempStr.substr(0, 6);
    // RDKit✔️❗:   while (!inStream->eof() && !inStream->fail() && lineBeg != "M  END" &&
    // RDKit✔️❗:          tempStr.substr(0, 4) != "$$$$") {
    // RDKit✔️❗:     if (tempStr[0] == 'A') {
    // RDKit✔️❗:       line++;
    // RDKit✔️❗:       std::string nextLine = getLine(inStream);
    // RDKit✔️❗:       if (lineBeg != "M  END") {
    // RDKit✔️❗:         ParseAtomAlias(mol, tempStr, nextLine, line);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (tempStr[0] == 'G') {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:           << " deprecated group abbreviation ignored on line " << line
    // RDKit✔️❗:           << std::endl;
    // RDKit✔️❗:       // we need to skip the next line, which holds the abbreviation:
    // RDKit✔️❗:       line++;
    // RDKit✔️❗:       tempStr = getLine(inStream);
    // RDKit✔️❗:     } else if (tempStr[0] == 'V') {
    // RDKit✔️❗:       ParseAtomValue(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "S  SKP") {
    // RDKit✔️❗:       int nToSkip = FileParserUtils::toInt(tempStr.substr(6, 3));
    // RDKit✔️❗:       if (nToSkip < 0) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "negative skip value " << nToSkip << " on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       for (unsigned int i = 0; i < static_cast<unsigned int>(nToSkip); ++i) {
    // RDKit✔️❗:         ++line;
    // RDKit✔️❗:         tempStr = getLine(inStream);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (lineBeg == "M  ALS") {
    // RDKit✔️❗:       ParseNewAtomList(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  ISO") {
    // RDKit✔️❗:       ParseIsotopeLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  RGP") {
    // RDKit✔️❗:       ParseRGroupLabels(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  RBC") {
    // RDKit✔️❗:       ParseRingBondCountLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  SUB") {
    // RDKit✔️❗:       ParseSubstitutionCountLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  UNS") {
    // RDKit✔️❗:       ParseUnsaturationLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  CHG") {
    // RDKit✔️❗:       ParseChargeLine(mol, tempStr, firstChargeLine, line);
    // RDKit✔️❗:       firstChargeLine = false;
    // RDKit✔️❗:     } else if (lineBeg == "M  RAD") {
    // RDKit✔️❗:       ParseRadicalLine(mol, tempStr, firstChargeLine, line);
    // RDKit✔️❗:       firstChargeLine = false;
    // RDKit✔️❗:     } else if (lineBeg == "M  PXA") {
    // RDKit✔️❗:       ParsePXALine(mol, tempStr, line);
    // RDKit✔️❗:
    // RDKit✔️❗:       /* SGroup parsing start */
    // RDKit✔️❗:     } else if (lineBeg == "M  STY") {
    // RDKit✔️❗:       ParseSGroupV2000STYLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SST") {
    // RDKit✔️❗:       ParseSGroupV2000SSTLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SLB") {
    // RDKit✔️❗:       ParseSGroupV2000SLBLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SCN") {
    // RDKit✔️❗:       ParseSGroupV2000SCNLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SDS") {
    // RDKit✔️❗:       ParseSGroupV2000SDSLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SAL" || lineBeg == "M  SBL" ||
    // RDKit✔️❗:                lineBeg == "M  SPA") {
    // RDKit✔️❗:       ParseSGroupV2000VectorDataLine(sGroupMap, mol, tempStr, line,
    // RDKit✔️❗:                                      strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SMT") {
    // RDKit✔️❗:       ParseSGroupV2000SMTLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SDI") {
    // RDKit✔️❗:       ParseSGroupV2000SDILine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  CRS") {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Unsupported SGroup subtype '" << lineBeg << "' on line "
    // RDKit✔️❗:              << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } else if (lineBeg == "M  SBV") {
    // RDKit✔️❗:       ParseSGroupV2000SBVLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SDT") {
    // RDKit✔️❗:       ParseSGroupV2000SDTLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SDD") {
    // RDKit✔️❗:       ParseSGroupV2000SDDLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SCD" || lineBeg == "M  SED") {
    // RDKit✔️❗:       ParseSGroupV2000SCDSEDLine(sGroupMap, dataFieldsMap, mol, tempStr, line,
    // RDKit✔️❗:                                  strictParsing, SCDcounter, lastDataSGroup,
    // RDKit✔️❗:                                  currentDataField);
    // RDKit✔️❗:     } else if (lineBeg == "M  SPL") {
    // RDKit✔️❗:       ParseSGroupV2000SPLLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SNC") {
    // RDKit✔️❗:       ParseSGroupV2000SNCLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SAP") {
    // RDKit✔️❗:       ParseSGroupV2000SAPLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SCL") {
    // RDKit✔️❗:       ParseSGroupV2000SCLLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  SBT") {
    // RDKit✔️❗:       ParseSGroupV2000SBTLine(sGroupMap, mol, tempStr, line, strictParsing);
    // RDKit✔️❗:
    // RDKit✔️❗:       /* SGroup parsing end */
    // RDKit✔️❗:     } else if (lineBeg == "M  ZBO") {
    // RDKit✔️❗:       ParseZBOLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  ZCH") {
    // RDKit✔️❗:       ParseZCHLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  HYD") {
    // RDKit✔️❗:       ParseHYDLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  MRV") {
    // RDKit✔️❗:       ParseMarvinSmartsLine(mol, tempStr, line);
    // RDKit✔️❗:     } else if (lineBeg == "M  APO") {
    // RDKit✔️❗:       ParseAttachPointLine(mol, tempStr, line, strictParsing);
    // RDKit✔️❗:     } else if (lineBeg == "M  LIN") {
    // RDKit✔️❗:       ParseLinkNodeLine(mol, tempStr, line);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     line++;
    // RDKit✔️❗:     tempStr = getLine(inStream);
    // RDKit✔️❗:     lineBeg = tempStr.substr(0, 6);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (tempStr[0] == 'M' && tempStr.substr(0, 6) == "M  END") {
    // RDKit✔️❗:     // All went well, make final updates to SGroups, and add them to Mol
    // RDKit✔️❗:     for (auto &sgroup : sGroupMap) {
    // RDKit✔️❗:       if (sgroup.second.getIsValid()) {
    // RDKit✔️❗:         sgroup.second.setProp("DATAFIELDS", dataFieldsMap[sgroup.first]);
    // RDKit✔️❗:         sgroup.second.setIsValid(checkAttachmentPointsAreValid(mol, sgroup));
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (sgroup.second.getIsValid()) {
    // RDKit✔️❗:         addSubstanceGroup(*mol, sgroup.second);
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "SGroup " << sgroup.first << " is invalid";
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:               << errout.str() << " and will be ignored" << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     fileComplete = true;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return fileComplete;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseMolBlockProperties
    // END RDKIT CPP BODY: parse_v2000_sgroup_line

    match rdkit_substr(line, 0, 6) {
        "M  STY" => parse_v2000_sgroup_sty_line(line, line_number, params, state),
        "M  SST" => parse_v2000_sgroup_sst_line(line, line_number, params, state),
        "M  SLB" => parse_v2000_sgroup_slb_line(line, line_number, params, state),
        "M  SCN" => parse_v2000_sgroup_scn_line(line, line_number, params, state),
        "M  SDS" => parse_v2000_sgroup_sds_line(line, line_number, params, state),
        "M  SAL" | "M  SBL" | "M  SPA" => {
            parse_v2000_sgroup_vector_data_line(line, line_number, params, state)
        }
        "M  SMT" => parse_v2000_sgroup_smt_line(line, line_number, params, state),
        "M  SDI" => parse_v2000_sgroup_sdi_line(line, line_number, params, state),
        "M  SBV" => parse_v2000_sgroup_sbv_line(line, line_number, params, state),
        "M  SDT" => parse_v2000_sgroup_sdt_line(line, line_number, params, state),
        "M  SDD" => parse_v2000_sgroup_sdd_line(line, line_number, params, state),
        "M  SCD" | "M  SED" => parse_v2000_sgroup_scd_sed_line(line, line_number, params, state),
        "M  SPL" => parse_v2000_sgroup_spl_line(line, line_number, params, state),
        "M  SNC" => parse_v2000_sgroup_snc_line(line, line_number, params, state),
        "M  SAP" => parse_v2000_sgroup_sap_line(line, line_number, params, state),
        "M  SCL" => parse_v2000_sgroup_scl_line(line, line_number, params, state),
        "M  SBT" => parse_v2000_sgroup_sbt_line(line, line_number, params, state),
        _ => Err(SdfReadError::Parse(format!(
            "Unknown SGroup line '{}' on line {line_number}",
            rdkit_substr(line, 0, 6)
        ))),
    }
}

fn parse_sgroup_int_field(
    line: &str,
    line_number: usize,
    pos: &mut usize,
    is_field_counter: bool,
) -> Result<u32, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_sgroup_int_field
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: unsigned int ParseSGroupIntField
    // RDKit✔️❗:   ++pos;  // Account for separation space
    // RDKit✔️❗:   unsigned int fieldValue;
    // RDKit✔️❗:   size_t len = 3 - isFieldCounter;  // field counters are smaller
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     fieldValue = FileParserUtils::toInt(text.substr(pos, len));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(pos, len) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   } catch (const std::out_of_range &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "SGroup line too short: '" << text << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   pos += len;
    // RDKit✔️❗:   return fieldValue;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: unsigned int ParseSGroupIntField
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: unsigned int ParseSGroupIntField(bool &ok, ...)
    // RDKit✔️❗:   ok = true;
    // RDKit✔️❗:   unsigned int res = 0;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     res = ParseSGroupIntField(text, line, pos, isFieldCounter);
    // RDKit✔️❗:   } catch (const std::exception &e) {
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw;
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       ok = false;
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << e.what() << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: unsigned int ParseSGroupIntField(bool &ok, ...)
    // END RDKIT CPP BODY: parse_sgroup_int_field

    // Rust parses into i32 first and casts to u32 to preserve the C++ path
    // where FileParserUtils::toInt feeds an unsigned destination.
    *pos += 1;
    let len = 3 - usize::from(is_field_counter);
    if *pos > line.len() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    let field = rdkit_substr(line, *pos, len);
    if field.is_empty() {
        return Err(SdfReadError::Parse(format!(
            "Cannot convert '{field}' to int on line {line_number}"
        )));
    }
    let value = parse_rdkit_int(field, true).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{field}' to int on line {line_number}"
        ))
    })? as u32;
    *pos += len;
    Ok(value)
}

fn parse_sgroup_double_field(
    line: &str,
    line_number: usize,
    pos: &mut usize,
) -> Result<f64, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_sgroup_double_field
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: double ParseSGroupDoubleField
    // RDKit✔️❗:   size_t len = 10;
    // RDKit✔️❗:   double fieldValue;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     fieldValue = FileParserUtils::toDouble(text.substr(pos, len));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(pos, len)
    // RDKit✔️❗:            << "' to double on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   } catch (const std::out_of_range &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "SGroup line too short: '" << text << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   pos += len;
    // RDKit✔️❗:   return fieldValue;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: double ParseSGroupDoubleField
    // END RDKIT CPP BODY: parse_sgroup_double_field

    let len = 10_usize;
    if *pos > line.len() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    let field = rdkit_substr(line, *pos, len);
    if field.is_empty() {
        return Err(SdfReadError::Parse(format!(
            "Cannot convert '{field}' to double on line {line_number}"
        )));
    }
    let value = parse_rdkit_double(field, true).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{field}' to double on line {line_number}"
        ))
    })?;
    *pos += len;
    Ok(value)
}

fn find_sgroup_mut<'a>(
    state: &'a mut V2000PropertyState,
    sg_idx: u32,
    line_number: usize,
) -> Result<&'a mut SubstanceGroup, SdfReadError> {
    // BEGIN RDKIT CPP BODY: find_sgroup_mut
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: SubstanceGroup *FindSgIdx
    // RDKit✔️❗:   auto sgIt = sGroupMap.find(sgIdx);
    // RDKit✔️❗:   if (sgIt == sGroupMap.end()) {
    // RDKit✔️❗:     BOOST_LOG(rdWarningLog) << "SGroup " << sgIdx << " referenced on line "
    // RDKit✔️❗:                             << line << " not found." << std::endl;
    // RDKit✔️❗:     return nullptr;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return &sgIt->second;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: SubstanceGroup *FindSgIdx
    // END RDKIT CPP BODY: find_sgroup_mut

    state.sgroups.get_mut(&sg_idx).ok_or_else(|| {
        // RDKit logs and returns nullptr. COSMolKit parser errors are explicit
        // because silently ignoring a referenced structure would hide data loss.
        SdfReadError::Parse(format!(
            "SGroup {sg_idx} referenced on line {line_number} not found."
        ))
    })
}

fn sgroup_kind_from_rdkit_type(s: &str) -> SubstanceGroupKind {
    match s {
        "DAT" => SubstanceGroupKind::Data,
        "SUP" => SubstanceGroupKind::Superatom,
        "MUL" => SubstanceGroupKind::MultipleGroup,
        "SRU" => SubstanceGroupKind::StructuralRepeatUnit,
        "MON" => SubstanceGroupKind::Monomer,
        "COP" => SubstanceGroupKind::Copolymer,
        "CRO" => SubstanceGroupKind::Crosslink,
        "GRA" => SubstanceGroupKind::Graft,
        "MOD" => SubstanceGroupKind::Modification,
        "MER" => SubstanceGroupKind::Mer,
        "ANY" => SubstanceGroupKind::AnyPolymer,
        "COM" => SubstanceGroupKind::MixtureComponent,
        "MIX" => SubstanceGroupKind::Mixture,
        "FOR" => SubstanceGroupKind::Formulation,
        "GEN" => SubstanceGroupKind::Generic("GEN".to_string()),
        other => SubstanceGroupKind::Generic(other.to_string()),
    }
}

fn is_valid_rdkit_sgroup_type(s: &str) -> bool {
    matches!(
        s,
        "SRU"
            | "MON"
            | "COP"
            | "CRO"
            | "GRA"
            | "MOD"
            | "MER"
            | "ANY"
            | "COM"
            | "MIX"
            | "FOR"
            | "SUP"
            | "MUL"
            | "DAT"
            | "GEN"
    )
}

fn parse_v2000_sgroup_sty_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_sty_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000STYLine
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == "M  STY", "bad STY line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int pos = 6;
    // RDKit✔️❗:   bool ok;
    // RDKit✔️❗:   unsigned int nent =
    // RDKit✔️❗:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
    // RDKit✔️❗:   if (!ok) {
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   for (unsigned int ie = 0; ie < nent; ++ie) {
    // RDKit✔️❗:     if (text.size() < pos + 8) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "SGroup STY line too short: '" << text << "' on line " << line;
    // RDKit✔️❗:       SGroupWarnOrThrow<>(strictParsing, errout.str());
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     unsigned int sequenceId =
    // RDKit✔️❗:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
    // RDKit✔️❗:     if (!ok) {
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     std::string typ = text.substr(pos + 1, 3);
    // RDKit✔️❗:     if (SubstanceGroupChecks::isValidType(typ)) {
    // RDKit✔️❗:       auto sgroup = SubstanceGroup(mol, typ);
    // RDKit✔️❗:       sgroup.setProp<unsigned int>("index", sequenceId);
    // RDKit✔️❗:       sGroupMap.emplace(sequenceId, sgroup);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "S group " << typ << " on line " << line;
    // RDKit✔️❗:       SGroupWarnOrThrow<MolFileUnhandledFeatureException>(strictParsing,
    // RDKit✔️❗:                                                           errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     pos += 4;
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000STYLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_sty_line

    let mut pos = 6_usize;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 8 {
            return Err(SdfReadError::Parse(format!(
                "SGroup STY line too short: '{line}' on line {line_number}"
            )));
        }
        let sequence_id = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let typ = rdkit_substr(line, pos + 1, 3);
        if !is_valid_rdkit_sgroup_type(typ) {
            return Err(SdfReadError::Parse(format!(
                "S group {typ} on line {line_number}"
            )));
        }
        // RDKit creates a SubstanceGroup keyed by the file sequence id.
        // COSMolKit keeps row ids contiguous for molecule invariants and stores
        // the RDKit sequence id as explicit SGroup metadata.
        let mut sgroup = SubstanceGroup::new(
            SubstanceGroupId::new(state.sgroups.len()),
            sgroup_kind_from_rdkit_type(typ),
        );
        sgroup.set_rdkit_sequence_id(sequence_id);
        sgroup.set_prop("TYPE", typ);
        state.sgroups.insert(sequence_id, sgroup);
        pos += 4;
    }
    Ok(())
}

fn parse_v2000_sgroup_sst_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_sst_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SSTLine
    // RDKit✔️❗:   unsigned int pos = 6;
    // RDKit✔️❗:   unsigned int nent = ParseSGroupIntField(..., true);
    // RDKit✔️❗:   for (...) {
    // RDKit✔️❗:     unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:     SubstanceGroup *sgroup = FindSgIdx(...);
    // RDKit✔️❗:     std::string subType = text.substr(++pos, 3);
    // RDKit✔️❗:     if (!SubstanceGroupChecks::isValidSubType(subType)) { ... }
    // RDKit✔️❗:     sgroup->setProp("SUBTYPE", subType);
    // RDKit✔️❗:     pos += 3;
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SSTLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_sst_line

    let mut pos = 6_usize;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 8 {
            return Err(SdfReadError::Parse(format!(
                "SGroup SST line too short: '{line}' on line {line_number}"
            )));
        }
        let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let subtype = rdkit_substr(line, pos + 1, 3);
        if !matches!(subtype, "ALT" | "RAN" | "BLO") {
            return Err(SdfReadError::Parse(format!(
                "Unsupported SGroup subtype '{subtype}' on line {line_number}"
            )));
        }
        find_sgroup_mut(state, sg_idx, line_number)?.set_subtype(subtype);
        pos += 4;
    }
    Ok(())
}

fn parse_v2000_sgroup_vector_data_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_vector_data_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000VectorDataLine
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:
    // RDKit✔️❗:   std::string typ = text.substr(3, 3);
    // RDKit✔️❗:
    // RDKit✔️❗:   void (SubstanceGroup::*sGroupAddIndexedElement)(const int) = nullptr;
    // RDKit✔️❗:
    // RDKit✔️❗:   if (typ == "SAL") {
    // RDKit✔️❗:     sGroupAddIndexedElement = &SubstanceGroup::addAtomWithBookmark;
    // RDKit✔️❗:   } else if (typ == "SBL") {
    // RDKit✔️❗:     sGroupAddIndexedElement = &SubstanceGroup::addBondWithBookmark;
    // RDKit✔️❗:   } else if (typ == "SPA") {
    // RDKit✔️❗:     sGroupAddIndexedElement = &SubstanceGroup::addParentAtomWithBookmark;
    // RDKit✔️❗:   } else {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Unsupported SGroup line '" << typ
    // RDKit✔️❗:            << "' passed to Vector Data parser ";
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int pos = 6;
    // RDKit✔️❗:   bool ok;
    // RDKit✔️❗:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
    // RDKit✔️❗:   if (!ok) {
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    // RDKit✔️❗:   if (!sgroup) {
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   unsigned int nent =
    // RDKit✔️❗:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
    // RDKit✔️❗:   if (!ok) {
    // RDKit✔️❗:     sgroup->setIsValid(false);
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   for (unsigned int i = 0; i < nent; ++i) {
    // RDKit✔️❗:     if (text.size() < pos + 4) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "SGroup line too short: '" << text << "' on line " << line;
    // RDKit✔️❗:       SGroupWarnOrThrow<>(strictParsing, errout.str());
    // RDKit✔️❗:       sgroup->setIsValid(false);
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     unsigned int nbr = ParseSGroupIntField(ok, strictParsing, text, line, pos);
    // RDKit✔️❗:     if (!ok) {
    // RDKit✔️❗:       sgroup->setIsValid(false);
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       (sgroup->*sGroupAddIndexedElement)(nbr);
    // RDKit✔️❗:     } catch (const std::exception &e) {
    // RDKit✔️❗:       SGroupWarnOrThrow<>(strictParsing, e.what());
    // RDKit✔️❗:       sgroup->setIsValid(false);
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000VectorDataLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_vector_data_line

    // RDKit dispatches through a C++ pointer-to-member and invalidates SGroups
    // on bad references. COSMolKit spells out SAL/SBL branches and returns an
    // explicit parse error instead of carrying invalid hidden SGroup state.
    let typ = rdkit_substr(line, 3, 3);
    let mut pos = 6_usize;
    let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    let sgroup = find_sgroup_mut(state, sg_idx, line_number)?;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 4 {
            return Err(SdfReadError::Parse(format!(
                "SGroup line too short: '{line}' on line {line_number}"
            )));
        }
        let nbr = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        match typ {
            "SAL" => sgroup.push_atom(AtomId::new(nbr.checked_sub(1).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup atom index 0 out of range on line {line_number}"
                ))
            })? as usize)),
            "SBL" => sgroup.push_bond(BondId::new(nbr.checked_sub(1).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup bond index 0 out of range on line {line_number}"
                ))
            })? as usize)),
            "SPA" => sgroup.push_parent_atom(AtomId::new(nbr.checked_sub(1).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup parent atom index 0 out of range on line {line_number}"
                ))
            })? as usize)),
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Unexpected SGroup vector data type '{typ}' on line {line_number}"
                )));
            }
        }
    }
    Ok(())
}

fn parse_v2000_sgroup_slb_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_slb_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SLBLine
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == "M  SLB", "bad SLB line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int pos = 6;
    // RDKit✔️❗:   bool ok;
    // RDKit✔️❗:   unsigned int nent =
    // RDKit✔️❗:       ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
    // RDKit✔️❗:   if (!ok) {
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   for (unsigned int ie = 0; ie < nent; ++ie) {
    // RDKit✔️❗:     if (text.size() < pos + 8) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "SGroup SLB line too short: '" << text << "' on line " << line;
    // RDKit✔️❗:       SGroupWarnOrThrow<>(strictParsing, errout.str());
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     unsigned int sgIdx =
    // RDKit✔️❗:         ParseSGroupIntField(ok, strictParsing, text, line, pos);
    // RDKit✔️❗:     if (!ok) {
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    // RDKit✔️❗:     if (!sgroup) {
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     unsigned int id = ParseSGroupIntField(ok, strictParsing, text, line, pos);
    // RDKit✔️❗:     if (!ok) {
    // RDKit✔️❗:       sgroup->setIsValid(false);
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     if (id != 0 && !SubstanceGroupChecks::isSubstanceGroupIdFree(*mol, id)) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "SGroup ID '" << id
    // RDKit✔️❗:              << "' is assigned to more than one SGroup, on line " << line;
    // RDKit✔️❗:       SGroupWarnOrThrow<>(strictParsing, errout.str());
    // RDKit✔️❗:       sgroup->setIsValid(false);
    // RDKit✔️❗:       return;
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     sgroup->setProp<unsigned int>("ID", id);
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SLBLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_slb_line

    let mut pos = 6_usize;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 8 {
            return Err(SdfReadError::Parse(format!(
                "SGroup SLB line too short: '{line}' on line {line_number}"
            )));
        }
        let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let id = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let sgroup = find_sgroup_mut(state, sg_idx, line_number)?;
        sgroup.set_external_id(id);
    }
    Ok(())
}

fn sgroup_connection_from_rdkit(connect: &str) -> Option<SGroupConnection> {
    match connect {
        "HH" => Some(SGroupConnection::HeadToHead),
        "HT" => Some(SGroupConnection::HeadToTail),
        "EU" => Some(SGroupConnection::Either),
        _ => None,
    }
}

fn parse_v2000_sgroup_scn_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_scn_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SCNLine
    // RDKit✔️❗:   unsigned int pos = 6;
    // RDKit✔️❗:   unsigned int nent = ParseSGroupIntField(..., true);
    // RDKit✔️❗:   for (...) {
    // RDKit✔️❗:     if (text.size() < pos + 7) { ... }
    // RDKit✔️❗:     unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:     SubstanceGroup *sgroup = FindSgIdx(...);
    // RDKit✔️❗:     std::string connect = text.substr(++pos, 2);
    // RDKit✔️❗:     if (!SubstanceGroupChecks::isValidConnectType(connect)) { ... }
    // RDKit✔️❗:     sgroup->setProp("CONNECT", connect);
    // RDKit✔️❗:     pos += 3;
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SCNLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_scn_line

    let mut pos = 6_usize;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 7 {
            return Err(SdfReadError::Parse(format!(
                "SGroup SCN line too short: '{line}' on line {line_number}\n needed: {} found: {}",
                pos + 7,
                line.len()
            )));
        }
        let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let connect = rdkit_substr(line, pos + 1, 2);
        let connection = sgroup_connection_from_rdkit(connect).ok_or_else(|| {
            SdfReadError::Parse(format!(
                "Unsupported SGroup connection type '{connect}' on line {line_number}"
            ))
        })?;
        find_sgroup_mut(state, sg_idx, line_number)?.set_connection(connection);
        pos += 3;
    }
    Ok(())
}

fn parse_v2000_sgroup_sds_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_sds_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SDSLine
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 10) == "M  SDS EXP", "bad SDS line");
    // RDKit✔️❗:   unsigned int pos = 10;
    // RDKit✔️❗:   unsigned int nent = ParseSGroupIntField(..., true);
    // RDKit✔️❗:   for (...) {
    // RDKit✔️❗:     unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:     SubstanceGroup *sgroup = FindSgIdx(...);
    // RDKit✔️❗:     sgroup->setProp("ESTATE", "E");
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SDSLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_sds_line

    if !line.starts_with("M  SDS EXP") {
        return Err(SdfReadError::Parse(format!(
            "bad SDS line on line {line_number}"
        )));
    }
    let mut pos = 10_usize;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 4 {
            return Err(SdfReadError::Parse(format!(
                "SGroup SDS line too short: '{line}' on line {line_number}"
            )));
        }
        let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        find_sgroup_mut(state, sg_idx, line_number)?.set_expansion_state("E");
    }
    Ok(())
}

fn parse_v2000_sgroup_smt_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_smt_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SMTLine
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == "M  SMT", "bad SMT line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int pos = 6;
    // RDKit✔️❗:   bool ok;
    // RDKit✔️❗:   unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
    // RDKit✔️❗:   if (!ok) {
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    // RDKit✔️❗:   if (!sgroup) {
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   ++pos;
    // RDKit✔️❗:
    // RDKit✔️❗:   if (pos >= text.length()) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "SGroup line too short: '" << text << "' on line " << line;
    // RDKit✔️❗:     SGroupWarnOrThrow<>(strictParsing, errout.str());
    // RDKit✔️❗:     sgroup->setIsValid(false);
    // RDKit✔️❗:     return;
    // RDKit✔️❗:   }
    // RDKit✔️❗:   std::string label = text.substr(pos, text.length() - pos);
    // RDKit✔️❗:
    // RDKit✔️❗:   if (sgroup->getProp<std::string>("TYPE") ==
    // RDKit✔️❗:       "MUL") {  // Case of multiple groups
    // RDKit✔️❗:     sgroup->setProp("MULT", label);
    // RDKit✔️❗:
    // RDKit✔️❗:   } else {  // Case of abbreviation groups, but we might not have seen a SCL
    // RDKit✔️❗:             // line yet
    // RDKit✔️❗:     sgroup->setProp("LABEL", label);
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SMTLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_smt_line

    let mut pos = 6_usize;
    let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    let sgroup = find_sgroup_mut(state, sg_idx, line_number)?;
    pos += 1;
    if pos >= line.len() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    let label = &line[pos..];
    if sgroup.props().get("TYPE").is_some_and(|typ| typ == "MUL") {
        sgroup.set_prop("MULT", label);
    } else {
        sgroup.set_label(label);
    }
    Ok(())
}

fn parse_v2000_sgroup_sdi_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_sdi_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SDILine
    // RDKit✔️❗:   unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:   unsigned int nCoords = ParseSGroupIntField(..., true);
    // RDKit✔️❗:   if (nCoords != 4) { ... }
    // RDKit✔️❗:   for (unsigned int i = 0; i < 2; ++i) {
    // RDKit✔️❗:     double x = ParseSGroupDoubleField(...);
    // RDKit✔️❗:     double y = ParseSGroupDoubleField(...);
    // RDKit✔️❗:     bracket[i] = RDGeom::Point3D(x, y, z);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   sgroup->addBracket(bracket);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SDILine
    // END RDKIT CPP BODY: parse_v2000_sgroup_sdi_line

    let mut pos = 6_usize;
    let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    let sgroup = find_sgroup_mut(state, sg_idx, line_number)?;
    let n_coords = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    if n_coords != 4 {
        return Err(SdfReadError::Parse(format!(
            "Unexpected number of coordinates for SDI on line {line_number}"
        )));
    }
    let x1 = parse_sgroup_double_field(line, line_number, &mut pos)?;
    let y1 = parse_sgroup_double_field(line, line_number, &mut pos)?;
    let x2 = parse_sgroup_double_field(line, line_number, &mut pos)?;
    let y2 = parse_sgroup_double_field(line, line_number, &mut pos)?;
    sgroup.display_mut().brackets.push(SGroupBracket {
        p1: [x1, y1],
        p2: [x2, y2],
    });
    Ok(())
}

fn parse_v2000_sgroup_sbv_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_sbv_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SBVLine
    // RDKit✔️❗:   unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:   unsigned int bondMark = ParseSGroupIntField(...);
    // RDKit✔️❗:   Bond *bond = mol->getUniqueBondWithBookmark(bondMark);
    // RDKit✔️❗:   if (sgroup->getProp<std::string>("TYPE") == "SUP") {
    // RDKit✔️❗:     vector.x = ParseSGroupDoubleField(...);
    // RDKit✔️❗:     vector.y = ParseSGroupDoubleField(...);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   sgroup->addCState(bond->getIdx(), vector);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SBVLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_sbv_line

    let mut pos = 6_usize;
    let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    let sgroup = find_sgroup_mut(state, sg_idx, line_number)?;
    let bond_mark = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    let bond = BondId::new(bond_mark.checked_sub(1).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup bond index 0 out of range on line {line_number}"
        ))
    })? as usize);
    let vector = if sgroup.kind() == &SubstanceGroupKind::Superatom {
        [
            parse_sgroup_double_field(line, line_number, &mut pos)?,
            parse_sgroup_double_field(line, line_number, &mut pos)?,
        ]
    } else {
        [0.0, 0.0]
    };
    sgroup.push_cstate(SGroupCState { bond, vector });
    Ok(())
}

fn parse_v2000_sgroup_sdt_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_sdt_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SDTLine
    // RDKit✔️❗:   unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:   fieldName = text.substr(++pos, 30); boost::trim_right(fieldName);
    // RDKit✔️❗:   fieldType = text.substr(pos, 2); boost::trim_right(fieldType);
    // RDKit✔️❗:   fieldInfo = text.substr(pos, 20); boost::trim_right(fieldInfo);
    // RDKit✔️❗:   queryType = text.substr(pos, 2); boost::trim_right(queryType);
    // RDKit✔️❗:   queryOp = text.substr(pos, text.length() - pos); boost::trim_right(queryOp);
    // RDKit✔️❗:   sgroup->setProp("FIELDNAME"/"FIELDTYPE"/"FIELDINFO"/"QUERYTYPE"/"QUERYOP", ...);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SDTLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_sdt_line

    let mut pos = 6_usize;
    let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    let data = find_sgroup_mut(state, sg_idx, line_number)?.data_mut();
    pos += 1;
    let field_name = rdkit_substr(line, pos, 30).trim_end().to_string();
    pos += 30;
    let field_type = rdkit_substr(line, pos, 2).trim_end().to_string();
    pos += 2;
    let field_info = rdkit_substr(line, pos, 20).trim_end().to_string();
    pos += 20;
    let query_type = rdkit_substr(line, pos, 2).trim_end().to_string();
    pos += 2;
    let query_op = rdkit_substr(line, pos, line.len().saturating_sub(pos))
        .trim_end()
        .to_string();
    if !field_name.is_empty() {
        data.field_name = Some(field_name);
    }
    if !field_type.is_empty() {
        data.field_type = Some(field_type);
    }
    if !field_info.is_empty() {
        data.field_info = Some(field_info);
    }
    if !query_type.is_empty() {
        data.query_type = Some(query_type);
    }
    if !query_op.is_empty() {
        data.query_op = Some(query_op);
    }
    Ok(())
}

fn parse_v2000_sgroup_sdd_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_sdd_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SDDLine
    // RDKit✔️❗:   unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:   ++pos;
    // RDKit✔️❗:   if (pos < text.length()) { sgroup->setProp("FIELDDISP", text.substr(pos, ...)); }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SDDLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_sdd_line

    let mut pos = 6_usize;
    let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    pos += 1;
    if pos < line.len() {
        find_sgroup_mut(state, sg_idx, line_number)?
            .data_mut()
            .field_display = Some(line[pos..].to_string());
    }
    Ok(())
}

fn parse_v2000_sgroup_scd_sed_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_scd_sed_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SCDSEDLine
    // RDKit✔️❗:   std::string type = text.substr(3, 3);
    // RDKit✔️❗:   unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:   if (lastDataSGroup != 0 && lastDataSGroup != sgIdx) { ... }
    // RDKit✔️❗:   else if (lastDataSGroup == 0 && type == "SCD") { lastDataSGroup = sgIdx; }
    // RDKit✔️❗:   else if (type == "SED") { lastDataSGroup = 0; }
    // RDKit✔️❗:   if (strictParsing && type == "SCD" && counter > 2) { throw; }
    // RDKit✔️❗:   if (pos + 1 < text.length()) {
    // RDKit✔️❗:     currentDataField << text.substr(++pos, 69);
    // RDKit✔️❗:     if (type == "SED") {
    // RDKit✔️❗:       dataFieldsMap[sgIdx].push_back(trimmedData.substr(0, 200));
    // RDKit✔️❗:       currentDataField.str("");
    // RDKit✔️❗:       counter = 0;
    // RDKit✔️❗:     } else { ++counter; }
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SCDSEDLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_scd_sed_line

    let strict = params.strict_parsing;
    let mut pos = 6_usize;
    let typ = rdkit_substr(line, 3, 3);
    let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    if state.last_data_sgroup != 0 && state.last_data_sgroup != sg_idx {
        return Err(SdfReadError::Parse(format!(
            "Found a Data Field not matching the SGroup of the last Data Field at line {line_number}"
        )));
    } else if state.last_data_sgroup == 0 && typ == "SCD" {
        state.last_data_sgroup = sg_idx;
    } else if typ == "SED" {
        state.last_data_sgroup = 0;
    }
    if strict && typ == "SCD" && state.scd_counter > 2 {
        return Err(SdfReadError::Parse(format!(
            "Found too many consecutive SCD lines, (#{} at line {line_number}) for SGroup {sg_idx}",
            state.scd_counter + 1
        )));
    }
    if pos + 1 < line.len() {
        state
            .current_data_field
            .push_str(rdkit_substr(line, pos + 1, 69));
        if typ == "SED" {
            let trimmed = state.current_data_field.trim_end().to_string();
            find_sgroup_mut(state, sg_idx, line_number)?
                .data_mut()
                .values
                .push(rdkit_substr(&trimmed, 0, 200).to_string());
            state.current_data_field.clear();
            state.scd_counter = 0;
        } else {
            state.scd_counter += 1;
        }
    }
    Ok(())
}

fn parse_v2000_sgroup_spl_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_spl_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SPLLine
    // RDKit✔️❗:   unsigned int nent = ParseSGroupIntField(..., true);
    // RDKit✔️❗:   for (...) {
    // RDKit✔️❗:     unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:     unsigned int parentIdx = ParseSGroupIntField(text, line, pos);
    // RDKit✔️❗:     sgroup->setProp<unsigned int>("PARENT", parentIdx);
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SPLLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_spl_line

    let mut pos = 6_usize;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 8 {
            return Err(SdfReadError::Parse(format!(
                "SGroup SPL line too short: '{line}' on line {line_number}"
            )));
        }
        let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let parent_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let parent = state
            .sgroups
            .get(&parent_idx)
            .map(SubstanceGroup::id)
            .ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup {parent_idx} referenced on line {line_number} not found."
                ))
            })?;
        find_sgroup_mut(state, sg_idx, line_number)?.set_parent(parent);
    }
    Ok(())
}

fn parse_v2000_sgroup_snc_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_snc_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SNCLine
    // RDKit✔️❗:   unsigned int compno = ParseSGroupIntField(...);
    // RDKit✔️❗:   if (compno > 256u) { ... }
    // RDKit✔️❗:   sgroup->setProp<unsigned int>("COMPNO", compno);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SNCLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_snc_line

    let mut pos = 6_usize;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 8 {
            return Err(SdfReadError::Parse(format!(
                "SGroup SNC line too short: '{line}' on line {line_number}"
            )));
        }
        let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let compno = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        if compno > 256 {
            return Err(SdfReadError::Parse(format!(
                "SGroup SNC value over 256: '{compno}' on line {line_number}"
            )));
        }
        find_sgroup_mut(state, sg_idx, line_number)?.set_component_number(compno);
    }
    Ok(())
}

fn parse_v2000_sgroup_sap_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_sap_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SAPLine
    // RDKit✔️❗:   unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:   unsigned int nent = ParseSGroupIntField(..., true);
    // RDKit✔️❗:   unsigned int aIdxMark = ParseSGroupIntField(...);
    // RDKit✔️❗:   unsigned int lvIdxMark = ParseSGroupIntField(...);
    // RDKit✔️❗:   if (text.size() >= pos + 3) { id = text.substr(pos + 1, 2); pos += 3; }
    // RDKit✔️❗:   sgroup->addAttachPoint(aIdx, lvIdx, id);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SAPLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_sap_line

    let mut pos = 6_usize;
    let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    let sgroup = find_sgroup_mut(state, sg_idx, line_number)?;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 11 {
            return Err(SdfReadError::Parse(format!(
                "SGroup SAP line too short: '{line}' on line {line_number}"
            )));
        }
        let atom_mark = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let leaving_mark = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let label = if line.len() >= pos + 3 {
            let id = rdkit_substr(line, pos + 1, 2).to_string();
            pos += 3;
            Some(id)
        } else {
            None
        };
        sgroup.push_attach_point(SGroupAttachPoint {
            atom: AtomId::new(atom_mark.checked_sub(1).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup attach atom index 0 out of range on line {line_number}"
                ))
            })? as usize),
            leaving_atom: if leaving_mark == 0 {
                None
            } else {
                Some(AtomId::new((leaving_mark - 1) as usize))
            },
            label,
            order: None,
        });
    }
    Ok(())
}

fn parse_v2000_sgroup_scl_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_scl_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SCLLine
    // RDKit✔️❗:   unsigned int sgIdx = ParseSGroupIntField(...);
    // RDKit✔️❗:   if (pos + 1 >= text.length()) { ... }
    // RDKit✔️❗:   ++pos;
    // RDKit✔️❗:   sgroup->setProp("CLASS", text.substr(pos, text.length() - pos));
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SCLLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_scl_line

    let mut pos = 6_usize;
    let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
    if pos + 1 >= line.len() {
        return Err(SdfReadError::Parse(format!(
            "SGroup SCL line too short: '{line}' on line {line_number}"
        )));
    }
    pos += 1;
    find_sgroup_mut(state, sg_idx, line_number)?.set_class(&line[pos..]);
    Ok(())
}

fn parse_v2000_sgroup_sbt_line(
    line: &str,
    line_number: usize,
    _params: SdfReadParams,
    state: &mut V2000PropertyState,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v2000_sgroup_sbt_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SBTLine
    // RDKit✔️❗:   unsigned int bracketType = ParseSGroupIntField(...);
    // RDKit✔️❗:   if (bracketType == 0) { sgroup->setProp("BRKTYP", "BRACKET"); }
    // RDKit✔️❗:   else if (bracketType == 1) { sgroup->setProp("BRKTYP", "PAREN"); }
    // RDKit✔️❗:   else { Invalid SBT value ... }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: void ParseSGroupV2000SBTLine
    // END RDKIT CPP BODY: parse_v2000_sgroup_sbt_line

    let mut pos = 6_usize;
    let nent = parse_sgroup_int_field(line, line_number, &mut pos, true)?;
    for _ in 0..nent {
        if line.len() < pos + 8 {
            return Err(SdfReadError::Parse(format!(
                "SGroup SBT line too short: '{line}' on line {line_number}"
            )));
        }
        let sg_idx = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let bracket_type = parse_sgroup_int_field(line, line_number, &mut pos, false)?;
        let style = match bracket_type {
            0 => SGroupBracketStyle::Bracket,
            1 => SGroupBracketStyle::Parenthesis,
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Invalid SBT value '{bracket_type}' on line {line_number}"
                )));
            }
        };
        find_sgroup_mut(state, sg_idx, line_number)?.set_bracket_style(style);
    }
    Ok(())
}

fn parse_attach_point_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_attach_point_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseAttachPointLine
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == std::string("M  APO"), "bad APO line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nent;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   unsigned int spos = 9;
    // RDKit✔️❗:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️❗:     unsigned int aid = 0;
    // RDKit✔️❗:     int val = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️❗:           text.substr(spos, 4));
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️❗:         val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos, 4));
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (!aid || aid > mol->getNumAtoms()) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Bad APO specification on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       --aid;
    // RDKit✔️❗:       Atom *atom = mol->getAtomWithIdx(aid);
    // RDKit✔️❗:       if (!atom) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Atom " << aid << " from APO specification on line " << line
    // RDKit✔️❗:                << " not found";
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         if (val < 0 || val > 3) {
    // RDKit✔️❗:           std::ostringstream errout;
    // RDKit✔️❗:           errout << "Value " << val << " from APO specification on line "
    // RDKit✔️❗:                  << line << " is invalid";
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else if (val) {
    // RDKit✔️❗:           if (val == 3) {
    // RDKit✔️❗:             // this is -1 in v3k mol blocks, so use that:
    // RDKit✔️❗:             val = -1;
    // RDKit✔️❗:           }
    // RDKit✔️❗:           if (atom->hasProp(common_properties::molAttachPoint)) {
    // RDKit✔️❗:             std::ostringstream errout;
    // RDKit✔️❗:             errout << "Multiple ATTCHPT values for atom " << atom->getIdx() + 1
    // RDKit✔️❗:                    << " on line " << line;
    // RDKit✔️❗:             if (strictParsing) {
    // RDKit✔️❗:               throw FileParseException(errout.str());
    // RDKit✔️❗:             } else {
    // RDKit✔️❗:               BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:             }
    // RDKit✔️❗:           } else {
    // RDKit✔️❗:             atom->setProp(common_properties::molAttachPoint, val);
    // RDKit✔️❗:           }
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(spos, 4)
    // RDKit✔️❗:              << "' to int on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseAttachPointLine
    // END RDKIT CPP BODY: parse_attach_point_line

    let nent = parse_rdkit_unsigned(rdkit_substr(line, 6, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {line_number}",
            rdkit_substr(line, 6, 3)
        ))
    })?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_rdkit_unsigned(rdkit_substr(line, spos, 4), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, spos, 4)
            ))
        })?;
        spos += 4;
        let mut val = 0_i32;
        if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            val = parse_rdkit_int(rdkit_substr(line, spos, 4), true).map_err(|_| {
                SdfReadError::Parse(format!(
                    "Cannot convert '{}' to int on line {line_number}",
                    rdkit_substr(line, spos, 4)
                ))
            })?;
        }
        if aid == 0 || aid as usize > atoms.len() {
            return Err(SdfReadError::Parse(format!(
                "Bad APO specification on line {line_number}"
            )));
        }
        spos += 4;
        if !(0..=3).contains(&val) {
            return Err(SdfReadError::Parse(format!(
                "Value {val} from APO specification on line {line_number} is invalid"
            )));
        }
        if val != 0 {
            if val == 3 {
                val = -1;
            }
            let atom = atom_line_mut(atoms, aid, line_number)?;
            if atom.spec.prop("molAttachPoint").is_some() {
                if params.strict_parsing {
                    return Err(SdfReadError::Parse(format!(
                        "Multiple ATTCHPT values for atom {} on line {line_number}",
                        aid
                    )));
                }
            } else {
                atom.spec = atom
                    .spec
                    .clone()
                    .with_prop("molAttachPoint", val.to_string());
            }
        }
    }
    Ok(())
}

fn parse_zch_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_zch_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseZCHLine
    // RDKit✔️❗:   // part of Alex Clark's ZBO proposal
    // RDKit✔️❗:   // from JCIM 51:3149-57 (2011)
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == std::string("M  ZCH"), "bad ZCH line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nent;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   unsigned int spos = 9;
    // RDKit✔️❗:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️❗:     unsigned int aid = 0;
    // RDKit✔️❗:     int val = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️❗:           text.substr(spos, 4));
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️❗:         val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos, 4));
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (!aid || aid > mol->getNumAtoms()) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Bad ZCH specification on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       --aid;
    // RDKit✔️❗:       Atom *atom = mol->getAtomWithIdx(aid);
    // RDKit✔️❗:       if (!atom) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Atom " << aid << " from ZCH specification on line " << line
    // RDKit✔️❗:                << " not found";
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         atom->setFormalCharge(val);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(spos, 4)
    // RDKit✔️❗:              << "' to int on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseZCHLine
    // END RDKIT CPP BODY: parse_zch_line

    let nent = parse_rdkit_unsigned(rdkit_substr(line, 6, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {line_number}",
            rdkit_substr(line, 6, 3)
        ))
    })?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_rdkit_unsigned(rdkit_substr(line, spos, 4), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, spos, 4)
            ))
        })?;
        spos += 4;
        let mut val = 0_i32;
        if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            val = parse_rdkit_int(rdkit_substr(line, spos, 4), true).map_err(|_| {
                SdfReadError::Parse(format!(
                    "Cannot convert '{}' to int on line {line_number}",
                    rdkit_substr(line, spos, 4)
                ))
            })?;
        }
        if aid == 0 || aid as usize > atoms.len() {
            return Err(SdfReadError::Parse(format!(
                "Bad ZCH specification on line {line_number}"
            )));
        }
        spos += 4;
        let atom = atom_line_mut(atoms, aid, line_number)?;
        atom.spec = atom.spec.clone().with_formal_charge(val as i8);
    }
    Ok(())
}

fn parse_hyd_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_hyd_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseHYDLine
    // RDKit✔️❗:   // part of Alex Clark's ZBO proposal
    // RDKit✔️❗:   // from JCIM 51:3149-57 (2011)
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == std::string("M  HYD"), "bad HYD line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nent;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   unsigned int spos = 9;
    // RDKit✔️❗:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️❗:     unsigned int aid = 0;
    // RDKit✔️❗:     int val = -1;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️❗:           text.substr(spos, 4));
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️❗:         val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos, 4));
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (!aid || aid > mol->getNumAtoms()) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Bad HYD specification on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       --aid;
    // RDKit✔️❗:       Atom *atom = mol->getAtomWithIdx(aid);
    // RDKit✔️❗:       if (!atom) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Atom " << aid << " from HYD specification on line " << line
    // RDKit✔️❗:                << " not found";
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         if (val >= 0) {
    // RDKit✔️❗:           atom->setProp("_ZBO_H", true);
    // RDKit✔️❗:           atom->setNumExplicitHs(val);
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(spos, 4)
    // RDKit✔️❗:              << "' to int on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseHYDLine
    // END RDKIT CPP BODY: parse_hyd_line

    let nent = parse_rdkit_unsigned(rdkit_substr(line, 6, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {line_number}",
            rdkit_substr(line, 6, 3)
        ))
    })?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let aid = parse_rdkit_unsigned(rdkit_substr(line, spos, 4), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, spos, 4)
            ))
        })?;
        spos += 4;
        let mut val = -1_i32;
        if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            val = parse_rdkit_int(rdkit_substr(line, spos, 4), true).map_err(|_| {
                SdfReadError::Parse(format!(
                    "Cannot convert '{}' to int on line {line_number}",
                    rdkit_substr(line, spos, 4)
                ))
            })?;
        }
        if aid == 0 || aid as usize > atoms.len() {
            return Err(SdfReadError::Parse(format!(
                "Bad HYD specification on line {line_number}"
            )));
        }
        spos += 4;
        if val >= 0 {
            let atom = atom_line_mut(atoms, aid, line_number)?;
            atom.spec = atom
                .spec
                .clone()
                .with_prop("_ZBO_H", "1")
                .with_explicit_hydrogens(val as u8);
        }
    }
    Ok(())
}

fn parse_zbo_line(
    line: &str,
    line_number: usize,
    bonds: &mut [V2000BondLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_zbo_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseZBOLine
    // RDKit✔️❗:   // part of Alex Clark's ZBO proposal
    // RDKit✔️❗:   // from JCIM 51:3149-57 (2011)
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 6) == std::string("M  ZBO"), "bad ZBO line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nent;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   unsigned int spos = 9;
    // RDKit✔️❗:   for (unsigned int ie = 0; ie < nent; ie++) {
    // RDKit✔️❗:     unsigned int bid = 0;
    // RDKit✔️❗:     unsigned int order = 0;
    // RDKit✔️❗:     try {
    // RDKit✔️❗:       bid = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️❗:           text.substr(spos, 4));
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
    // RDKit✔️❗:         order = FileParserUtils::stripSpacesAndCast<unsigned int>(
    // RDKit✔️❗:             text.substr(spos, 4));
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (!bid || bid > mol->getNumBonds()) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Bad ZBO specification on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       spos += 4;
    // RDKit✔️❗:       --bid;
    // RDKit✔️❗:       Bond *bond = mol->getBondWithIdx(bid);
    // RDKit✔️❗:       if (!bond) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Bond " << bid << " from ZBO specification on line " << line
    // RDKit✔️❗:                << " not found";
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         bond->setBondType(Bond::BondType::ZERO);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Cannot convert '" << text.substr(spos, 4)
    // RDKit✔️❗:              << "' to int on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseZBOLine
    // END RDKIT CPP BODY: parse_zbo_line

    let nent = parse_rdkit_unsigned(rdkit_substr(line, 6, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {line_number}",
            rdkit_substr(line, 6, 3)
        ))
    })?;
    let mut spos = 9_usize;
    for _ in 0..nent {
        let bid = parse_rdkit_unsigned(rdkit_substr(line, spos, 4), true).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                rdkit_substr(line, spos, 4)
            ))
        })?;
        spos += 4;
        if line.len() >= spos + 4 && rdkit_substr(line, spos, 4) != "    " {
            let _ = parse_rdkit_unsigned(rdkit_substr(line, spos, 4), true).map_err(|_| {
                SdfReadError::Parse(format!(
                    "Cannot convert '{}' to int on line {line_number}",
                    rdkit_substr(line, spos, 4)
                ))
            })?;
        }
        if bid == 0 || bid as usize > bonds.len() {
            return Err(SdfReadError::Parse(format!(
                "Bad ZBO specification on line {line_number}"
            )));
        }
        spos += 4;
        let bond = &mut bonds[bid as usize - 1];
        bond.spec = bond
            .spec
            .clone()
            .with_order(BondOrder::Zero)
            .with_aromatic(false)
            .with_prop("_ZBO", "1");
    }
    Ok(())
}

fn parse_atom_alias_line(
    line: &str,
    next_line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_atom_alias_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseAtomAlias
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 2) == std::string("A "), "bad atom alias line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int idx;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(3, 3)) -
    // RDKit✔️❗:           1;
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(3, 3) << "' to int on line "
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   URANGE_CHECK(idx, mol->getNumAtoms());
    // RDKit✔️❗:   Atom *at = mol->getAtomWithIdx(idx);
    // RDKit✔️❗:   at->setProp(common_properties::molFileAlias, nextLine);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseAtomAlias
    // END RDKIT CPP BODY: parse_atom_alias_line

    let idx = parse_rdkit_unsigned(rdkit_substr(line, 3, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {line_number}",
            rdkit_substr(line, 3, 3)
        ))
    })?;
    if idx == 0 || idx as usize > atoms.len() {
        return Err(SdfReadError::Parse(format!(
            "Atom alias index {} out of range on line {line_number}",
            idx.saturating_sub(1)
        )));
    }
    let atom = atom_line_mut(atoms, idx, line_number)?;
    atom.spec = atom
        .spec
        .clone()
        .with_prop("molFileAlias", next_line.to_string());
    Ok(())
}

fn parse_atom_value_line(
    line: &str,
    line_number: usize,
    atoms: &mut [V2000AtomLine],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_atom_value_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseAtomValue
    // RDKit✔️❗:   PRECONDITION(mol, "bad mol");
    // RDKit✔️❗:   PRECONDITION(text.substr(0, 2) == std::string("V "), "bad atom value line");
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int idx;
    // RDKit✔️❗:   try {
    // RDKit✔️❗:     idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(3, 3)) -
    // RDKit✔️❗:           1;
    // RDKit✔️❗:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Cannot convert '" << text.substr(3, 3) << "' to int on line"
    // RDKit✔️❗:            << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   URANGE_CHECK(idx, mol->getNumAtoms());
    // RDKit✔️❗:   Atom *at = mol->getAtomWithIdx(idx);
    // RDKit✔️❗:   at->setProp(common_properties::molFileValue,
    // RDKit✔️❗:               text.substr(7, text.length() - 7));
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseAtomValue
    // END RDKIT CPP BODY: parse_atom_value_line

    let idx = parse_rdkit_unsigned(rdkit_substr(line, 3, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line{line_number}",
            rdkit_substr(line, 3, 3)
        ))
    })?;
    if idx == 0 || idx as usize > atoms.len() {
        return Err(SdfReadError::Parse(format!(
            "Atom value index {} out of range on line {line_number}",
            idx.saturating_sub(1)
        )));
    }
    let atom = atom_line_mut(atoms, idx, line_number)?;
    atom.spec = atom.spec.clone().with_prop(
        "molFileValue",
        line.get(7..).unwrap_or_default().to_string(),
    );
    Ok(())
}

fn skip_sgroup_lines<R: BufRead>(
    reader: &mut R,
    line: &str,
    line_number: &mut usize,
) -> Result<(), SdfReadError> {
    let n_to_skip = parse_rdkit_int(rdkit_substr(line, 6, 3), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {}",
            rdkit_substr(line, 6, 3),
            *line_number
        ))
    })?;
    if n_to_skip < 0 {
        return Err(SdfReadError::Parse(format!(
            "negative skip value {} on line {}",
            n_to_skip, *line_number
        )));
    }
    for _ in 0..n_to_skip as usize {
        *line_number += 1;
        let _ = read_rdkit_line(reader)?.ok_or_else(|| {
            SdfReadError::Parse("EOF hit while skipping S  SKP payload".to_string())
        })?;
    }
    Ok(())
}

fn parse_v3000_ctab<R: BufRead>(
    reader: &mut R,
    header: MolHeader,
    counts: CountsLine,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<ParsedMolBlock, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_ctab
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV3000CTAB
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:
    // RDKit✔️❗:   std::string tempStr;
    // RDKit✔️❗:   std::vector<std::string> splitLine;
    // RDKit✔️❗:
    // RDKit✔️❗:   bool fileComplete = false;
    // RDKit✔️❗:
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN CTAB") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN CTAB line not found on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   if (tempStr.size() < 8 || tempStr.substr(0, 7) != "COUNTS ") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Bad counts line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   std::string trimmed =
    // RDKit✔️❗:       boost::trim_copy(tempStr.substr(7, tempStr.length() - 7));
    // RDKit✔️❗:   boost::split(splitLine, trimmed, boost::is_any_of(" \t"),
    // RDKit✔️❗:                boost::token_compress_on);
    // RDKit✔️❗:   if (splitLine.size() < 2) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Bad counts line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   nAtoms = FileParserUtils::toUnsigned(splitLine[0]);
    // RDKit✔️❗:   nBonds = FileParserUtils::toUnsigned(splitLine[1]);
    // RDKit✔️❗:   conf = new Conformer(nAtoms);
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nSgroups = 0, n3DConstraints = 0, chiralFlag = 0;
    // RDKit✔️❗:
    // RDKit✔️❗:   if (splitLine.size() > 2) {
    // RDKit✔️❗:     nSgroups = FileParserUtils::toUnsigned(splitLine[2]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (splitLine.size() > 3) {
    // RDKit✔️❗:     n3DConstraints = FileParserUtils::toUnsigned(splitLine[3]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (splitLine.size() > 4) {
    // RDKit✔️❗:     chiralFlag = FileParserUtils::toUnsigned(splitLine[4]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   mol->setProp(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit✔️❗:
    // RDKit✔️❗:   if (nAtoms) {
    // RDKit✔️❗:     ParseV3000AtomBlock(inStream, line, nAtoms, mol, conf, strictParsing,
    // RDKit✔️❗:                         expectMacroAtoms);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (nBonds) {
    // RDKit✔️❗:     ParseV3000BondBlock(inStream, line, nBonds, mol, chiralityPossible);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   // do link nodes:
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   while (tempStr.length() > 8 && tempStr.substr(0, 8) == "LINKNODE") {
    // RDKit✔️❗:     boost::to_upper(tempStr);
    // RDKit✔️❗:     // if the line has nothing on it we just ignore it
    // RDKit✔️❗:     if (tempStr.size() > 9) {
    // RDKit✔️❗:       std::string existing = "";
    // RDKit✔️❗:       if (mol->getPropIfPresent(common_properties::molFileLinkNodes,
    // RDKit✔️❗:                                 existing)) {
    // RDKit✔️❗:         existing += "|";
    // RDKit✔️❗:       }
    // RDKit✔️❗:       existing += tempStr.substr(9);  // skip the "LINKNODE "
    // RDKit✔️❗:       mol->setProp(common_properties::molFileLinkNodes, existing);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   bool sgroupFound = false;
    // RDKit✔️❗:   bool obj3dFound = false;
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   while (tempStr.length() > 5 && tempStr.substr(0, 5) == "BEGIN") {
    // RDKit✔️❗:     if (tempStr.length() >= 12 && tempStr.substr(0, 12) == "BEGIN SGROUP") {
    // RDKit✔️❗:       if (sgroupFound) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN SGROUP found more than once on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:
    // RDKit✔️❗:       } else if (!nSgroups) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN SGROUP  found but Sgroups NOT expected on line "
    // RDKit✔️❗:                << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:           // Prepare to read a lot of sgroups
    // RDKit✔️❗:           nSgroups = std::numeric_limits<unsigned int>::max();
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       sgroupFound = true;
    // RDKit✔️❗:       tempStr =
    // RDKit✔️❗:           ParseV3000SGroupsBlock(inStream, line, nSgroups, mol, strictParsing);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:       if (tempStr.length() < 10 || tempStr.substr(0, 10) != "END SGROUP") {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "END SGROUP line not found on line " << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:         boost::to_upper(tempStr);
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:     } else if (tempStr.length() >= 15 &&
    // RDKit✔️❗:                tempStr.substr(6, 10) == "COLLECTION") {
    // RDKit✔️❗:       tempStr = parseEnhancedStereo(inStream, line, mol, strictParsing);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:     } else if (tempStr.length() >= 11 &&
    // RDKit✔️❗:                tempStr.substr(0, 11) == "BEGIN OBJ3D") {
    // RDKit✔️❗:       if (obj3dFound) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN OBJ3D found more than once on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (!n3DConstraints) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN OBJ3D found but 3n3DConstraints NOT expected on line "
    // RDKit✔️❗:                << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:           << "3D constraint information in mol block ignored at line " << line
    // RDKit✔️❗:           << std::endl;
    // RDKit✔️❗:       obj3dFound = true;
    // RDKit✔️❗:       for (unsigned int i = 0; i < n3DConstraints; ++i) {
    // RDKit✔️❗:         tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       }
    // RDKit✔️❗:       tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:       if (tempStr.length() < 9 || tempStr.substr(0, 9) != "END OBJ3D") {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "END OBJ3D line not found on line " << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       // skip blocks we don't know how to read
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << "skipping block at line " << line << ": '"
    // RDKit✔️❗:                               << tempStr << "'" << std::endl;
    // RDKit✔️❗:       while (tempStr.length() < 3 || tempStr.substr(0, 3) != "END") {
    // RDKit✔️❗:         tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       }
    // RDKit✔️❗:       tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (nSgroups && !sgroupFound) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN SGROUP line not found on line " << line;
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (n3DConstraints && !obj3dFound) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN OBJ3D line not found on line " << line;
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END CTAB") {
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException("END CTAB line not found");
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << "END CTAB line not found." << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (expectMEND) {
    // RDKit✔️❗:     tempStr = getLine(inStream);
    // RDKit✔️❗:     ++line;
    // RDKit✔️❗:     if (tempStr[0] == 'M' && tempStr.substr(0, 6) == "M  END") {
    // RDKit✔️❗:       fileComplete = true;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   } else {
    // RDKit✔️❗:     fileComplete = true;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   auto is3d = calculate3dFlag(*mol, *conf, chiralityPossible);
    // RDKit✔️❗:   conf->set3D(is3d);
    // RDKit✔️❗:   mol->addConformer(conf, true);
    // RDKit✔️❗:   conf = nullptr;
    // RDKit✔️❗:
    // RDKit✔️❗:   return fileComplete;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV3000CTAB
    // END RDKIT CPP BODY: parse_v3000_ctab

    let _ = counts;
    let begin_ctab = get_v3000_line(reader, line_number, params)?;
    if begin_ctab.to_ascii_uppercase().len() < 10
        || &begin_ctab.to_ascii_uppercase()[..10] != "BEGIN CTAB"
    {
        return Err(SdfReadError::Parse(format!(
            "BEGIN CTAB line not found on line {}",
            *line_number
        )));
    }

    let counts_line = get_v3000_line(reader, line_number, params)?;
    let v3000_counts = parse_v3000_counts_line(&counts_line, *line_number, params)?;
    let atom_lines = if v3000_counts.atom_count == 0 {
        Vec::new()
    } else {
        parse_v3000_atom_block(reader, v3000_counts.atom_count, line_number, params)?
    };
    let bond_lines = if v3000_counts.bond_count == 0 {
        Vec::new()
    } else {
        parse_v3000_bond_block(reader, v3000_counts.bond_count, line_number, params)?
    };
    let chirality_possible = bond_lines.iter().any(|bond| {
        !matches!(
            bond.spec.direction(),
            BondDirection::None | BondDirection::Unknown
        )
    });

    let mut atom_id_by_mol_idx = BTreeMap::new();
    let mut bond_id_by_mol_idx = BTreeMap::new();
    let mut builder = MoleculeBuilder::new();
    if let Some(name) = header.name.as_ref() {
        builder = builder.with_name(name.clone());
    }
    builder = builder
        .with_property("_MolFileChiralFlag", v3000_counts.chiral_flag.to_string())
        .with_property("_MolFileInfoLine", header.info.clone())
        .with_property("_MolFileComments", header.comments.clone());

    for atom_line in &atom_lines {
        let id = builder.add_atom(atom_line.spec.clone());
        atom_id_by_mol_idx.insert(atom_line.mol_idx, id);
    }
    if atom_lines
        .iter()
        .filter_map(|atom_line| atom_line.spec.query())
        .any(atom_query_needs_scan)
    {
        builder = builder.with_property("_NeedsQueryScan", "1");
    }
    let coords = atom_lines
        .iter()
        .map(|line| line.coord_3d)
        .collect::<Vec<_>>();

    for bond_line in &bond_lines {
        let begin = *atom_id_by_mol_idx
            .get(&bond_line.begin_mol_idx)
            .ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "Bond endpoint {} out of range on line {}",
                    bond_line.begin_mol_idx, bond_line.line_number
                ))
            })?;
        let end = *atom_id_by_mol_idx
            .get(&bond_line.end_mol_idx)
            .ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "Bond endpoint {} out of range on line {}",
                    bond_line.end_mol_idx, bond_line.line_number
                ))
            })?;
        let mut spec = BondSpec::new(begin, end, bond_line.spec.order())
            .with_aromatic(bond_line.spec.is_aromatic())
            .with_conjugated(bond_line.spec.is_conjugated())
            .with_direction(bond_line.spec.direction())
            .with_stereo(bond_line.spec.stereo());
        if let Some(query) = bond_line.spec.query().cloned() {
            spec = spec.with_query(query);
        }
        if let Some(stereo_atoms) = bond_line.spec.stereo_atoms() {
            spec = spec.with_stereo_atoms(stereo_atoms[0], stereo_atoms[1]);
        }
        for (key, value) in bond_line.spec.props() {
            spec = spec.with_prop(key.clone(), value.clone());
        }
        let bond_id = builder.add_bond(spec).map_err(molecule_build_error)?;
        bond_id_by_mol_idx.insert(bond_line.mol_idx, bond_id);
    }

    let mut substance_groups = Vec::new();
    let mut stereo_groups = Vec::new();
    let mut sgroup_found = false;
    let mut collection_found = false;
    let mut obj3d_found = false;
    let mut link_nodes = String::new();

    let mut next_line = get_v3000_line(reader, line_number, params)?;
    let mut next_upper = next_line.to_ascii_uppercase();
    while next_upper.len() > 8 && next_upper.starts_with("LINKNODE") {
        if next_line.len() > 9 {
            if !link_nodes.is_empty() {
                link_nodes.push('|');
            }
            link_nodes.push_str(&next_line[9..]);
        }
        next_line = get_v3000_line(reader, line_number, params)?;
        next_upper = next_line.to_ascii_uppercase();
    }
    if !link_nodes.is_empty() {
        builder = builder.with_property("_MolFileLinkNodes", link_nodes);
    }

    while next_upper.len() > 5 && next_upper.starts_with("BEGIN") {
        if next_upper.len() >= 12 && next_upper.starts_with("BEGIN SGROUP") {
            if sgroup_found {
                return Err(SdfReadError::Parse(format!(
                    "BEGIN SGROUP found more than once on line {}",
                    *line_number
                )));
            }
            if v3000_counts.v3000_sgroup_count == 0 && params.strict_parsing {
                return Err(SdfReadError::Parse(format!(
                    "BEGIN SGROUP  found but Sgroups NOT expected on line {}",
                    *line_number
                )));
            }
            sgroup_found = true;
            substance_groups = parse_v3000_sgroup_block(
                reader,
                v3000_counts.v3000_sgroup_count,
                line_number,
                params,
                &atom_id_by_mol_idx,
                &bond_id_by_mol_idx,
            )?;
            next_line = get_v3000_line(reader, line_number, params)?;
            next_upper = next_line.to_ascii_uppercase();
        } else if next_upper.len() >= 16 && next_upper.starts_with("BEGIN COLLECTION") {
            if collection_found {
                return Err(SdfReadError::Parse(format!(
                    "BEGIN COLLECTION found more than once on line {}",
                    *line_number
                )));
            }
            collection_found = true;
            stereo_groups =
                parse_v3000_collection_block(reader, line_number, params, &atom_id_by_mol_idx)?;
            next_line = get_v3000_line(reader, line_number, params)?;
            next_upper = next_line.to_ascii_uppercase();
        } else if next_upper.len() >= 11 && next_upper.starts_with("BEGIN OBJ3D") {
            if obj3d_found {
                return Err(SdfReadError::Parse(format!(
                    "BEGIN OBJ3D found more than once on line {}",
                    *line_number
                )));
            }
            if v3000_counts.v3000_obj3d_count == 0 && params.strict_parsing {
                return Err(SdfReadError::Parse(format!(
                    "BEGIN OBJ3D found but 3n3DConstraints NOT expected on line {}",
                    *line_number
                )));
            }
            obj3d_found = true;
            parse_v3000_obj3d_block(reader, v3000_counts.v3000_obj3d_count, line_number, params)?;
            next_line = get_v3000_line(reader, line_number, params)?;
            next_upper = next_line.to_ascii_uppercase();
        } else {
            while next_upper.len() < 3 || !next_upper.starts_with("END") {
                next_line = get_v3000_line(reader, line_number, params)?;
                next_upper = next_line.to_ascii_uppercase();
            }
            next_line = get_v3000_line(reader, line_number, params)?;
            next_upper = next_line.to_ascii_uppercase();
        }
    }
    if v3000_counts.v3000_sgroup_count != 0 {
        if !sgroup_found && params.strict_parsing {
            return Err(SdfReadError::Parse(format!(
                "BEGIN SGROUP line not found on line {}",
                *line_number
            )));
        }
    }
    if v3000_counts.v3000_obj3d_count != 0 && !obj3d_found && params.strict_parsing {
        return Err(SdfReadError::Parse(format!(
            "BEGIN OBJ3D line not found on line {}",
            *line_number
        )));
    }
    if (next_upper.len() < 8 || !next_upper.starts_with("END CTAB")) && params.strict_parsing {
        return Err(SdfReadError::Parse("END CTAB line not found".to_string()));
    }

    let m_end = read_rdkit_line(reader)?.unwrap_or_default();
    *line_number += 1;
    if !m_end.starts_with("M  END") {
        return Err(SdfReadError::Parse(
            "Problems encountered parsing Mol data, M  END missing around line ".to_string()
                + &line_number.to_string(),
        ));
    }

    let is_3d = resolve_coordinate_3d_flag(
        calculate_rdkit_3d_flag(
            molfile_info_marks_3d(&header.info),
            &coords,
            chirality_possible,
        ),
        params.coordinate_mode,
    );
    builder
        .add_conformer(Conformer3D::new(0, coords, is_3d))
        .map_err(molecule_build_error)?;
    for substance_group in substance_groups {
        builder
            .add_substance_group(substance_group)
            .map_err(molecule_build_error)?;
    }
    for stereo_group in stereo_groups {
        builder
            .add_stereo_group(stereo_group)
            .map_err(molecule_build_error)?;
    }
    let molecule = builder
        .with_topology_trust(crate::TopologyTrust::TrustedGraph)
        .build()
        .map_err(molecule_build_error)?;
    Ok(ParsedMolBlock {
        molecule,
        header,
        chirality_possible,
    })
}

fn parse_v3000_counts_line(
    line: &str,
    line_number: usize,
    params: SdfReadParams,
) -> Result<CountsLine, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_counts_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV3000CTAB
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:
    // RDKit✔️❗:   std::string tempStr;
    // RDKit✔️❗:   std::vector<std::string> splitLine;
    // RDKit✔️❗:
    // RDKit✔️❗:   bool fileComplete = false;
    // RDKit✔️❗:
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN CTAB") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN CTAB line not found on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   if (tempStr.size() < 8 || tempStr.substr(0, 7) != "COUNTS ") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Bad counts line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   std::string trimmed =
    // RDKit✔️❗:       boost::trim_copy(tempStr.substr(7, tempStr.length() - 7));
    // RDKit✔️❗:   boost::split(splitLine, trimmed, boost::is_any_of(" \t"),
    // RDKit✔️❗:                boost::token_compress_on);
    // RDKit✔️❗:   if (splitLine.size() < 2) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Bad counts line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   nAtoms = FileParserUtils::toUnsigned(splitLine[0]);
    // RDKit✔️❗:   nBonds = FileParserUtils::toUnsigned(splitLine[1]);
    // RDKit✔️❗:   conf = new Conformer(nAtoms);
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nSgroups = 0, n3DConstraints = 0, chiralFlag = 0;
    // RDKit✔️❗:
    // RDKit✔️❗:   if (splitLine.size() > 2) {
    // RDKit✔️❗:     nSgroups = FileParserUtils::toUnsigned(splitLine[2]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (splitLine.size() > 3) {
    // RDKit✔️❗:     n3DConstraints = FileParserUtils::toUnsigned(splitLine[3]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (splitLine.size() > 4) {
    // RDKit✔️❗:     chiralFlag = FileParserUtils::toUnsigned(splitLine[4]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   mol->setProp(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit✔️❗:
    // RDKit✔️❗:   if (nAtoms) {
    // RDKit✔️❗:     ParseV3000AtomBlock(inStream, line, nAtoms, mol, conf, strictParsing,
    // RDKit✔️❗:                         expectMacroAtoms);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (nBonds) {
    // RDKit✔️❗:     ParseV3000BondBlock(inStream, line, nBonds, mol, chiralityPossible);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   // do link nodes:
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   while (tempStr.length() > 8 && tempStr.substr(0, 8) == "LINKNODE") {
    // RDKit✔️❗:     boost::to_upper(tempStr);
    // RDKit✔️❗:     // if the line has nothing on it we just ignore it
    // RDKit✔️❗:     if (tempStr.size() > 9) {
    // RDKit✔️❗:       std::string existing = "";
    // RDKit✔️❗:       if (mol->getPropIfPresent(common_properties::molFileLinkNodes,
    // RDKit✔️❗:                                 existing)) {
    // RDKit✔️❗:         existing += "|";
    // RDKit✔️❗:       }
    // RDKit✔️❗:       existing += tempStr.substr(9);  // skip the "LINKNODE "
    // RDKit✔️❗:       mol->setProp(common_properties::molFileLinkNodes, existing);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   bool sgroupFound = false;
    // RDKit✔️❗:   bool obj3dFound = false;
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   while (tempStr.length() > 5 && tempStr.substr(0, 5) == "BEGIN") {
    // RDKit✔️❗:     if (tempStr.length() >= 12 && tempStr.substr(0, 12) == "BEGIN SGROUP") {
    // RDKit✔️❗:       if (sgroupFound) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN SGROUP found more than once on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:
    // RDKit✔️❗:       } else if (!nSgroups) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN SGROUP  found but Sgroups NOT expected on line "
    // RDKit✔️❗:                << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:           // Prepare to read a lot of sgroups
    // RDKit✔️❗:           nSgroups = std::numeric_limits<unsigned int>::max();
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       sgroupFound = true;
    // RDKit✔️❗:       tempStr =
    // RDKit✔️❗:           ParseV3000SGroupsBlock(inStream, line, nSgroups, mol, strictParsing);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:       if (tempStr.length() < 10 || tempStr.substr(0, 10) != "END SGROUP") {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "END SGROUP line not found on line " << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:         boost::to_upper(tempStr);
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:     } else if (tempStr.length() >= 15 &&
    // RDKit✔️❗:                tempStr.substr(6, 10) == "COLLECTION") {
    // RDKit✔️❗:       tempStr = parseEnhancedStereo(inStream, line, mol, strictParsing);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:     } else if (tempStr.length() >= 11 &&
    // RDKit✔️❗:                tempStr.substr(0, 11) == "BEGIN OBJ3D") {
    // RDKit✔️❗:       if (obj3dFound) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN OBJ3D found more than once on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (!n3DConstraints) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN OBJ3D found but 3n3DConstraints NOT expected on line "
    // RDKit✔️❗:                << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:           << "3D constraint information in mol block ignored at line " << line
    // RDKit✔️❗:           << std::endl;
    // RDKit✔️❗:       obj3dFound = true;
    // RDKit✔️❗:       for (unsigned int i = 0; i < n3DConstraints; ++i) {
    // RDKit✔️❗:         tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       }
    // RDKit✔️❗:       tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:       if (tempStr.length() < 9 || tempStr.substr(0, 9) != "END OBJ3D") {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "END OBJ3D line not found on line " << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       // skip blocks we don't know how to read
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << "skipping block at line " << line << ": '"
    // RDKit✔️❗:                               << tempStr << "'" << std::endl;
    // RDKit✔️❗:       while (tempStr.length() < 3 || tempStr.substr(0, 3) != "END") {
    // RDKit✔️❗:         tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       }
    // RDKit✔️❗:       tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (nSgroups && !sgroupFound) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN SGROUP line not found on line " << line;
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (n3DConstraints && !obj3dFound) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN OBJ3D line not found on line " << line;
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END CTAB") {
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException("END CTAB line not found");
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << "END CTAB line not found." << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (expectMEND) {
    // RDKit✔️❗:     tempStr = getLine(inStream);
    // RDKit✔️❗:     ++line;
    // RDKit✔️❗:     if (tempStr[0] == 'M' && tempStr.substr(0, 6) == "M  END") {
    // RDKit✔️❗:       fileComplete = true;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   } else {
    // RDKit✔️❗:     fileComplete = true;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   auto is3d = calculate3dFlag(*mol, *conf, chiralityPossible);
    // RDKit✔️❗:   conf->set3D(is3d);
    // RDKit✔️❗:   mol->addConformer(conf, true);
    // RDKit✔️❗:   conf = nullptr;
    // RDKit✔️❗:
    // RDKit✔️❗:   return fileComplete;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV3000CTAB
    // END RDKIT CPP BODY: parse_v3000_counts_line

    let _ = params;
    if line.len() < 8 || !line.to_ascii_uppercase().starts_with("COUNTS ") {
        return Err(SdfReadError::Parse(format!(
            "Bad counts line : '{line}' on line {line_number}"
        )));
    }
    let fields = line[7..].split_whitespace().collect::<Vec<_>>();
    if fields.len() < 2 {
        return Err(SdfReadError::Parse(format!(
            "Bad counts line : '{line}' on line {line_number}"
        )));
    }
    let atom_count = parse_rdkit_unsigned(fields[0], false).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to unsigned int on line {line_number}",
            fields[0]
        ))
    })? as usize;
    let bond_count = parse_rdkit_unsigned(fields[1], false).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to unsigned int on line {line_number}",
            fields[1]
        ))
    })? as usize;
    let v3000_sgroup_count = match fields.get(2) {
        Some(field) => parse_rdkit_unsigned(field, false).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{field}' to unsigned int on line {line_number}"
            ))
        })? as usize,
        None => 0,
    };
    let v3000_obj3d_count = match fields.get(3) {
        Some(field) => parse_rdkit_unsigned(field, false).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{field}' to unsigned int on line {line_number}"
            ))
        })? as usize,
        None => 0,
    };
    let chiral_flag = match fields.get(4) {
        Some(field) => parse_rdkit_unsigned(field, false).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{field}' to unsigned int on line {line_number}"
            ))
        })?,
        None => 0,
    };

    Ok(CountsLine {
        atom_count,
        bond_count,
        chiral_flag,
        ctab_version: CtabVersion::V3000,
        v3000_sgroup_count,
        v3000_obj3d_count,
    })
}

fn parse_v3000_atom_block<R: BufRead>(
    reader: &mut R,
    atom_count: usize,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<Vec<V3000AtomLine>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_atom_block
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000AtomBlock
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(nAtoms > 0, "bad atom count");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:   PRECONDITION(conf, "bad conformer");
    // RDKit✔️❗:   std::vector<std::string> splitLine;
    // RDKit✔️❗:
    // RDKit✔️❗:   auto inl = getV3000Line(inStream, line);
    // RDKit✔️❗:   std::string_view tempStr = inl;
    // RDKit✔️❗:   if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN ATOM") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN ATOM line not found on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️❗:     inl = getV3000Line(inStream, line);
    // RDKit✔️❗:     tempStr = inl;
    // RDKit✔️❗:     auto trimmed = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:
    // RDKit✔️❗:     std::vector<std::string_view> tokens;
    // RDKit✔️❗:     std::vector<std::string_view>::iterator token;
    // RDKit✔️❗:
    // RDKit✔️❗:     tokenizeV3000Line(trimmed, tokens);
    // RDKit✔️❗:     token = tokens.begin();
    // RDKit✔️❗:
    // RDKit✔️❗:     if (token == tokens.end()) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Bad atom line : '" << tempStr << "' on line" << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     unsigned int molIdx = 0;
    // RDKit✔️❗:     std::from_chars(token->data(), token->data() + token->size(), molIdx);
    // RDKit✔️❗:
    // RDKit✔️❗:     // start with the symbol:
    // RDKit✔️❗:     ++token;
    // RDKit✔️❗:     if (token == tokens.end()) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Bad atom line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     // before we parse the symbol, we need to know if the atom has a class attr.
    // RDKit✔️❗:     // if it does, it is a macro atom reference, and we do not need to parse the
    // RDKit✔️❗:     // symbol.  (the single letter codes can be the same as element sysmbols or
    // RDKit✔️❗:     // special query names)
    // RDKit✔️❗:
    // RDKit✔️❗:     auto isMacroAtom = false;
    // RDKit✔️❗:     if (expectMacroAtoms) {
    // RDKit✔️❗:       auto lookAheadToken = token + 1;
    // RDKit✔️❗:       while (lookAheadToken != tokens.end()) {
    // RDKit✔️❗:         std::string prop;
    // RDKit✔️❗:         std::string_view val;
    // RDKit✔️❗:         if (splitAssignToken(*lookAheadToken, prop, val) && prop == "CLASS") {
    // RDKit✔️❗:           isMacroAtom = true;
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         }
    // RDKit✔️❗:         ++lookAheadToken;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     Atom *atom = nullptr;
    // RDKit✔️❗:     if (isMacroAtom) {
    // RDKit✔️❗:       atom = new Atom(0);
    // RDKit✔️❗:       atom->setAtomicNum(0);
    // RDKit✔️❗:       std::string tcopy(*token);
    // RDKit✔️❗:       atom->setProp(common_properties::dummyLabel, tcopy);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       atom = ParseV3000AtomSymbol(*token, line, strictParsing);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     // now the position;
    // RDKit✔️❗:     RDGeom::Point3D pos;
    // RDKit✔️❗:     ++token;
    // RDKit✔️❗:     if (token == tokens.end()) {
    // RDKit✔️❗:       delete atom;
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Bad atom line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     pos.x = atof(std::string(*token).c_str());
    // RDKit✔️❗:     ++token;
    // RDKit✔️❗:     if (token == tokens.end()) {
    // RDKit✔️❗:       delete atom;
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Bad atom line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     pos.y = atof(std::string(*token).c_str());
    // RDKit✔️❗:     ++token;
    // RDKit✔️❗:     if (token == tokens.end()) {
    // RDKit✔️❗:       delete atom;
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Bad atom line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     pos.z = atof(std::string(*token).c_str());
    // RDKit✔️❗:     // the map number:
    // RDKit✔️❗:     ++token;
    // RDKit✔️❗:     if (token == tokens.end()) {
    // RDKit✔️❗:       delete atom;
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Bad atom line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     int mapNum = atoi(std::string(*token).c_str());
    // RDKit✔️❗:     if (mapNum > 0) {
    // RDKit✔️❗:       atom->setProp(common_properties::molAtomMapNumber, mapNum);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     ++token;
    // RDKit✔️❗:
    // RDKit✔️❗:     unsigned int aid = mol->addAtom(atom, false, true);
    // RDKit✔️❗:
    // RDKit✔️❗:     // additional properties this may change the atom,
    // RDKit✔️❗:     // so be careful with it:
    // RDKit✔️❗:     ParseV3000AtomProps(mol, atom, token, tokens, line, strictParsing);
    // RDKit✔️❗:
    // RDKit✔️❗:     mol->setAtomBookmark(atom, molIdx);
    // RDKit✔️❗:     conf->setAtomPos(aid, pos);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   inl = getV3000Line(inStream, line);
    // RDKit✔️❗:   tempStr = inl;
    // RDKit✔️❗:   if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END ATOM") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "END ATOM line not found on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000AtomBlock
    // END RDKIT CPP BODY: parse_v3000_atom_block

    let begin = get_v3000_line(reader, line_number, params)?;
    if begin.len() < 10 || !begin.starts_with("BEGIN ATOM") {
        return Err(SdfReadError::Parse(format!(
            "BEGIN ATOM line not found on line {}",
            *line_number
        )));
    }

    let mut atoms = Vec::with_capacity(atom_count);
    for _ in 0..atom_count {
        let temp_str = get_v3000_line(reader, line_number, params)?;
        let trimmed = temp_str.trim().to_string();
        let tokens = tokenize_v3000_line(&trimmed, *line_number)?;
        if tokens.len() < 6 {
            return Err(SdfReadError::Parse(format!(
                "Bad atom line : '{temp_str}' on line {}",
                *line_number
            )));
        }

        let mol_idx = parse_rdkit_unsigned(&tokens[0], false).unwrap_or(0);
        let mut spec = parse_v3000_atom_symbol(&tokens[1], *line_number, params)?;
        let x = parse_rdkit_double(&tokens[2], false).unwrap_or(0.0);
        let y = parse_rdkit_double(&tokens[3], false).unwrap_or(0.0);
        let z = parse_rdkit_double(&tokens[4], false).unwrap_or(0.0);
        let map_num = parse_rdkit_int(&tokens[5], false).unwrap_or(0);
        if map_num > 0 {
            spec = spec
                .with_atom_map(map_num as u32)
                .with_prop("molAtomMapNumber", map_num.to_string());
        }
        spec = parse_v3000_atom_props(&tokens[6..], *line_number, params, spec)?;
        atoms.push(V3000AtomLine {
            line_number: *line_number,
            mol_idx,
            tokens,
            spec,
            coord_3d: [x, y, z],
        });
    }

    let end = get_v3000_line(reader, line_number, params)?;
    if end.len() < 8 || !end.starts_with("END ATOM") {
        return Err(SdfReadError::Parse(format!(
            "END ATOM line not found on line {}",
            *line_number
        )));
    }

    Ok(atoms)
}

fn parse_v3000_atom_symbol(
    token: &str,
    line_number: usize,
    params: SdfReadParams,
) -> Result<AtomSpec, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_atom_symbol
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: Atom *ParseV3000AtomSymbol
    // RDKit✔️❗:
    // RDKit✔️❗:   bool negate = false;
    // RDKit✔️❗:   token = FileParserUtils::strip(token);
    // RDKit✔️❗:   if (token.size() > 3 && (token[0] == 'N' || token[0] == 'n') &&
    // RDKit✔️❗:       (token[1] == 'O' || token[1] == 'o') &&
    // RDKit✔️❗:       (token[2] == 'T' || token[2] == 't')) {
    // RDKit✔️❗:     negate = true;
    // RDKit✔️❗:     token = token.substr(3, token.size() - 3);
    // RDKit✔️❗:     token = FileParserUtils::strip(token);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   std::unique_ptr<Atom> res;
    // RDKit✔️❗:   if (token[0] == '[') {
    // RDKit✔️❗:     // atom list:
    // RDKit✔️❗:     if (token.back() != ']') {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Bad atom token '" << token << "' on line: " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     token = token.substr(1, token.size() - 2);
    // RDKit✔️❗:
    // RDKit✔️❗:     std::vector<std::string> splitToken;
    // RDKit✔️❗:     boost::split(splitToken, token, boost::is_any_of(","));
    // RDKit✔️❗:
    // RDKit✔️❗:     for (std::vector<std::string>::const_iterator stIt = splitToken.begin();
    // RDKit✔️❗:          stIt != splitToken.end(); ++stIt) {
    // RDKit✔️❗:       std::string_view stoken = *stIt;
    // RDKit✔️❗:       std::string atSymb(FileParserUtils::strip(stoken));
    // RDKit✔️❗:       if (atSymb.empty()) {
    // RDKit✔️❗:         continue;
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (atSymb.size() == 2 && atSymb[1] >= 'A' && atSymb[1] <= 'Z') {
    // RDKit✔️❗:         atSymb[1] = static_cast<char>(tolower(atSymb[1]));
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       int atNum = PeriodicTable::getTable()->getAtomicNumber(atSymb);
    // RDKit✔️❗:       if (!res) {
    // RDKit✔️❗:         res.reset(new QueryAtom(atNum));
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         res->expandQuery(makeAtomNumQuery(atNum), Queries::COMPOSITE_OR, true);
    // RDKit✔️❗:       }
    // RDKit✔️❗:       // we want the atomic number of the query itself to always be zero
    // RDKit✔️❗:       // this was Github #8820 and #8823
    // RDKit✔️❗:       res->setAtomicNum(0);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     res->getQuery()->setNegation(negate);
    // RDKit✔️❗:   } else {
    // RDKit✔️❗:     if (negate) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "NOT tokens only supported for atom lists. line " << line;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     // it's a normal CTAB atom symbol:
    // RDKit✔️❗:     // NOTE: "R" and "R0"-"R99" are not in the v3K CTAB spec, but we're going to
    // RDKit✔️❗:     // support them anyway
    // RDKit✔️❗:     bool isComplexQueryName =
    // RDKit✔️❗:         std::find(complexQueries.begin(), complexQueries.end(), token) !=
    // RDKit✔️❗:         complexQueries.end();
    // RDKit✔️❗:     if (isComplexQueryName || token == "R" ||
    // RDKit✔️❗:         (token[0] == 'R' && token >= "R0" && token <= "R99") || token == "R#" ||
    // RDKit✔️❗:         token == "*") {
    // RDKit✔️❗:       if (isComplexQueryName || token == "*") {
    // RDKit✔️❗:         res.reset(new QueryAtom(0));
    // RDKit✔️❗:         if (token == "*") {
    // RDKit✔️❗:           // according to the MDL spec, these match anything
    // RDKit✔️❗:           res->setQuery(makeAtomNullQuery());
    // RDKit✔️❗:         } else if (isComplexQueryName) {
    // RDKit✔️❗:           convertComplexNameToQuery(res.get(), token);
    // RDKit✔️❗:         }
    // RDKit✔️❗:         // queries have no implicit Hs:
    // RDKit✔️❗:         res->setNoImplicit(true);
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         res.reset(new Atom(1));
    // RDKit✔️❗:         res->setAtomicNum(0);
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (token[0] == 'R' && token >= "R0" && token <= "R99") {
    // RDKit✔️❗:         auto rlabel = token.substr(1, token.length() - 1);
    // RDKit✔️❗:         int rnumber;
    // RDKit✔️❗:         try {
    // RDKit✔️❗:           rnumber = boost::lexical_cast<int>(rlabel);
    // RDKit✔️❗:         } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:           rnumber = -1;
    // RDKit✔️❗:         }
    // RDKit✔️❗:         if (rnumber >= 0) {
    // RDKit✔️❗:           res->setIsotope(rnumber);
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (token[0] == 'R') {
    // RDKit✔️❗:         // we used to skip R# here because that really should be handled by an
    // RDKit✔️❗:         // RGP spec, but that turned out to not be permissive enough... <sigh>
    // RDKit✔️❗:         setRGPProps(token, res.get());
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (token == "D") {  // mol blocks support "D" and "T" as
    // RDKit✔️❗:                                 // shorthand... handle that.
    // RDKit✔️❗:       res.reset(new Atom(1));
    // RDKit✔️❗:       res->setIsotope(2);
    // RDKit✔️❗:     } else if (token == "T") {  // mol blocks support "D" and "T" as
    // RDKit✔️❗:                                 // shorthand... handle that.
    // RDKit✔️❗:       res.reset(new Atom(1));
    // RDKit✔️❗:       res->setIsotope(3);
    // RDKit✔️❗:     } else if (token == "Pol" || token == "Mod") {
    // RDKit✔️❗:       res.reset(new Atom(0));
    // RDKit✔️❗:       res->setProp(common_properties::dummyLabel, std::string(token));
    // RDKit✔️❗:     } else if (GenericGroups::genericMatchers.find(std::string(token)) !=
    // RDKit✔️❗:                GenericGroups::genericMatchers.end()) {
    // RDKit✔️❗:       res.reset(new QueryAtom(0));
    // RDKit✔️❗:       res->setProp(common_properties::atomLabel, std::string(token));
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       std::string tcopy(token);
    // RDKit✔️❗:       res.reset(new Atom(0));
    // RDKit✔️❗:       lookupAtomicNumber(res.get(), tcopy, strictParsing);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   POSTCONDITION(res, "no atom built");
    // RDKit✔️❗:   return res.release();
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: Atom *ParseV3000AtomSymbol
    // END RDKIT CPP BODY: parse_v3000_atom_symbol

    let mut token = token.trim();
    let mut negate = false;
    if token.len() > 3 && token[..3].eq_ignore_ascii_case("NOT") {
        negate = true;
        token = token[3..].trim();
    }

    if token.starts_with('[') {
        if !token.ends_with(']') {
            return Err(SdfReadError::Parse(format!(
                "Bad atom token : '{token}' on line {line_number}"
            )));
        }
        let list = &token[1..token.len() - 1];
        let mut atomic_numbers = Vec::new();
        for part in list.split(',') {
            let at_symb = part.trim();
            if at_symb.is_empty() {
                continue;
            }
            let atomic_number = atomic_number_from_mol_symbol(at_symb, params.strict_parsing)?;
            atomic_numbers.push(atomic_number);
        }
        let query = if negate {
            AtomQueryPredicate::AtomicNumberNotIn(atomic_numbers)
        } else {
            AtomQueryPredicate::AtomicNumberIn(atomic_numbers)
        };
        return Ok(AtomSpec::new(Element::DUMMY).with_query(QueryNode::predicate(query)));
    }

    if negate {
        return Err(SdfReadError::Parse(format!(
            "NOT tokens only supported for atom lists. line {line_number}"
        )));
    }

    match token {
        "*" => Ok(AtomSpec::new(Element::DUMMY)
            .with_query(QueryNode::predicate(AtomQueryPredicate::Any))
            .with_no_implicit(true)),
        "D" => Ok(AtomSpec::new(Element::H).with_isotope(2)),
        "T" => Ok(AtomSpec::new(Element::H).with_isotope(3)),
        "R" | "R#" => Ok(AtomSpec::new(Element::DUMMY).with_prop("dummyLabel", token.to_string())),
        _ if token.starts_with('R')
            && token.len() <= 3
            && token[1..].chars().all(|char| char.is_ascii_digit()) =>
        {
            let number = token[1..].parse::<u32>().unwrap_or(0);
            Ok(AtomSpec::new(Element::DUMMY)
                .with_isotope(number as u16)
                .with_prop("dummyLabel", token.to_string()))
        }
        "A" | "Q" | "L" | "LP" => Ok(AtomSpec::new(Element::DUMMY)
            .with_query(QueryNode::predicate(AtomQueryPredicate::MolFileAlias(
                token.to_string(),
            )))
            .with_no_implicit(true)),
        other => {
            let atomic_number = atomic_number_from_mol_symbol(other, params.strict_parsing)?;
            let element = Element::from_atomic_number(atomic_number).unwrap_or(Element::DUMMY);
            Ok(AtomSpec::new(element))
        }
    }
}

fn parse_v3000_atom_props(
    tokens: &[String],
    line_number: usize,
    params: SdfReadParams,
    mut spec: AtomSpec,
) -> Result<AtomSpec, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_atom_props
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000AtomProps
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:   PRECONDITION(atom, "bad atom");
    // RDKit✔️❗:   std::ostringstream errout;
    // RDKit✔️❗:   while (token != tokens.end()) {
    // RDKit✔️❗:     std::string prop;
    // RDKit✔️❗:     std::string_view val;
    // RDKit✔️❗:     if (!splitAssignToken(*token, prop, val)) {
    // RDKit✔️❗:       errout << "Invalid atom property: '" << *token << "' for atom "
    // RDKit✔️❗:              << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     if (prop == "CHG") {
    // RDKit✔️❗:       auto charge = FileParserUtils::toInt(val);
    // RDKit✔️❗:       if (!atom->hasQuery()) {
    // RDKit✔️❗:         atom->setFormalCharge(charge);
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         atom->expandQuery(makeAtomFormalChargeQuery(charge));
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "RAD") {
    // RDKit✔️❗:       // FIX handle queries here
    // RDKit✔️❗:       switch (FileParserUtils::toInt(val)) {
    // RDKit✔️❗:         case 0:
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 1:
    // RDKit✔️❗:           atom->setNumRadicalElectrons(2);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 2:
    // RDKit✔️❗:           atom->setNumRadicalElectrons(1);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 3:
    // RDKit✔️❗:           atom->setNumRadicalElectrons(2);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         default:
    // RDKit✔️❗:           errout << "Unrecognized RAD value " << val << " for atom "
    // RDKit✔️❗:                  << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "MASS") {
    // RDKit✔️❗:       // the documentation for V3000 CTABs says that this should contain the
    // RDKit✔️❗:       // "absolute atomic weight" (whatever that means).
    // RDKit✔️❗:       // Online examples seem to have integer (isotope) values and Marvin
    // RDKit✔️❗:       // won't even read something that has a float. We'll go with the int
    // RDKit✔️❗:       int v;
    // RDKit✔️❗:       double dv;
    // RDKit✔️❗:       try {
    // RDKit✔️❗:         v = FileParserUtils::toInt(val);
    // RDKit✔️❗:       } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:         try {
    // RDKit✔️❗:           dv = FileParserUtils::toDouble(val);
    // RDKit✔️❗:           v = static_cast<int>(floor(dv));
    // RDKit✔️❗:         } catch (boost::bad_lexical_cast &) {
    // RDKit✔️❗:           v = -1;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (v < 0) {
    // RDKit✔️❗:         errout << "Bad value for MASS :" << val << " for atom "
    // RDKit✔️❗:                << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         if (!atom->hasQuery()) {
    // RDKit✔️❗:           atom->setIsotope(v);
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           atom->expandQuery(makeAtomIsotopeQuery(v));
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "CFG") {
    // RDKit✔️❗:       auto cfg = FileParserUtils::toInt(val);
    // RDKit✔️❗:       switch (cfg) {
    // RDKit✔️❗:         case 0:
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         case 1:
    // RDKit✔️❗:         case 2:
    // RDKit✔️❗:         case 3:
    // RDKit✔️❗:           atom->setProp(common_properties::molParity, cfg);
    // RDKit✔️❗:           break;
    // RDKit✔️❗:         default:
    // RDKit✔️❗:           errout << "Unrecognized CFG value : " << val << " for atom "
    // RDKit✔️❗:                  << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "HCOUNT") {
    // RDKit✔️❗:       if (val != "0") {
    // RDKit✔️❗:         auto hcount = FileParserUtils::toInt(val);
    // RDKit✔️❗:         if (!atom->hasQuery()) {
    // RDKit✔️❗:           atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
    // RDKit✔️❗:         }
    // RDKit✔️❗:         if (hcount == -1) {
    // RDKit✔️❗:           hcount = 0;
    // RDKit✔️❗:         }
    // RDKit✔️❗:         if (hcount > 0) {
    // RDKit✔️❗:           ATOM_EQUALS_QUERY *oq = makeAtomImplicitHCountQuery(hcount);
    // RDKit✔️❗:           auto nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(
    // RDKit✔️❗:               hcount, oq->getDataFunc(),
    // RDKit✔️❗:               std::string("less_") + oq->getDescription());
    // RDKit✔️❗:           atom->expandQuery(nq);
    // RDKit✔️❗:           delete oq;
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           atom->expandQuery(makeAtomImplicitHCountQuery(0));
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "UNSAT") {
    // RDKit✔️❗:       if (val == "1") {
    // RDKit✔️❗:         if (!atom->hasQuery()) {
    // RDKit✔️❗:           atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
    // RDKit✔️❗:         }
    // RDKit✔️❗:         atom->expandQuery(makeAtomUnsaturatedQuery());
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "RBCNT") {
    // RDKit✔️❗:       if (val != "0") {
    // RDKit✔️❗:         auto rbcount = FileParserUtils::toInt(val);
    // RDKit✔️❗:         if (!atom->hasQuery()) {
    // RDKit✔️❗:           atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
    // RDKit✔️❗:         }
    // RDKit✔️❗:         atom->setProp(common_properties::molRingBondCount, rbcount);
    // RDKit✔️❗:         if (rbcount == -1) {
    // RDKit✔️❗:           rbcount = 0;
    // RDKit✔️❗:         } else if (rbcount == -2) {
    // RDKit✔️❗:           // Ring bonds can only be counted during post processing
    // RDKit✔️❗:           mol->setProp(common_properties::_NeedsQueryScan, 1);
    // RDKit✔️❗:           rbcount = 0xDEADBEEF;
    // RDKit✔️❗:         } else if (rbcount > 4) {
    // RDKit✔️❗:           rbcount = 4;
    // RDKit✔️❗:         }
    // RDKit✔️❗:         atom->expandQuery(makeAtomRingBondCountQuery(rbcount));
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "VAL") {
    // RDKit✔️❗:       if (val != "0") {
    // RDKit✔️❗:         auto totval = FileParserUtils::toInt(val);
    // RDKit✔️❗:         atom->setProp(common_properties::molTotValence, totval);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "RGROUPS") {
    // RDKit✔️❗:       ParseV3000RGroups(mol, atom, val, line);
    // RDKit✔️❗:       // FIX
    // RDKit✔️❗:     } else if (prop == "STBOX") {
    // RDKit✔️❗:       if (val != "0") {
    // RDKit✔️❗:         auto ival = FileParserUtils::toInt(val);
    // RDKit✔️❗:         atom->setProp(common_properties::molStereoCare, ival);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "SUBST") {
    // RDKit✔️❗:       if (val != "0") {
    // RDKit✔️❗:         auto ival = FileParserUtils::toInt(val);
    // RDKit✔️❗:         atom->setProp(common_properties::molSubstCount, ival);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "EXACHG") {
    // RDKit✔️❗:       if (val != "0") {
    // RDKit✔️❗:         auto ival = FileParserUtils::toInt(val);
    // RDKit✔️❗:         atom->setProp(common_properties::molRxnExactChange, ival);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "INVRET") {
    // RDKit✔️❗:       if (val != "0") {
    // RDKit✔️❗:         auto ival = FileParserUtils::toInt(val);
    // RDKit✔️❗:         atom->setProp(common_properties::molInversionFlag, ival);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "ATTCHPT") {
    // RDKit✔️❗:       if (val != "0") {
    // RDKit✔️❗:         auto ival = FileParserUtils::toInt(val);
    // RDKit✔️❗:         if (atom->hasProp(common_properties::molAttachPoint)) {
    // RDKit✔️❗:           errout << "Multiple ATTCHPT values for atom " << atom->getIdx() + 1
    // RDKit✔️❗:                  << " on line " << line;
    // RDKit✔️❗:           if (strictParsing) {
    // RDKit✔️❗:             throw FileParseException(errout.str());
    // RDKit✔️❗:           } else {
    // RDKit✔️❗:             BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:             errout.str(std::string());
    // RDKit✔️❗:           }
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           atom->setProp(common_properties::molAttachPoint, ival);
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "ATTCHORD") {
    // RDKit✔️❗:       // there are two kinds of ATTCHORD
    // RDKit✔️❗:       // one is for template instances and looks like this: ATTCHORD=(4 1 Al 3
    // RDKit✔️❗:       // Br)
    // RDKit✔️❗:
    // RDKit✔️❗:       if (val.substr(0, 1) == "(") {
    // RDKit✔️❗:         // this is a template instance
    // RDKit✔️❗:
    // RDKit✔️❗:         val = val.substr(1, val.size() - 2);
    // RDKit✔️❗:         std::vector<std::string> splitToken;
    // RDKit✔️❗:         boost::split(splitToken, val, boost::is_any_of(" \t"));
    // RDKit✔️❗:
    // RDKit✔️❗:         unsigned int itemCount = 0;
    // RDKit✔️❗:         if (splitToken.size() > 0) {
    // RDKit✔️❗:           itemCount = FileParserUtils::toInt(splitToken[0]);
    // RDKit✔️❗:         }
    // RDKit✔️❗:
    // RDKit✔️❗:         if (itemCount == 0 || itemCount % 2 != 0 ||
    // RDKit✔️❗:             splitToken.size() != itemCount + 1) {
    // RDKit✔️❗:           errout << "Invalid ATTCHORD value: '" << val << "' for atom "
    // RDKit✔️❗:                  << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         }
    // RDKit✔️❗:         std::vector<std::pair<unsigned int, std::string>> attchOrds;
    // RDKit✔️❗:         for (unsigned int i = 1; i < itemCount; i += 2) {
    // RDKit✔️❗:           unsigned int idx = FileParserUtils::toInt(splitToken[i]);
    // RDKit✔️❗:           // check for uniqueness
    // RDKit✔️❗:           for (const auto &[aidx, lbl] : attchOrds) {
    // RDKit✔️❗:             if (idx == aidx + 1 || splitToken[i + 1] == lbl) {
    // RDKit✔️❗:               errout << "Invalid ATTCHORD value: '" << val << "' for atom "
    // RDKit✔️❗:                      << atom->getIdx() + 1 << " on line " << line << std::endl;
    // RDKit✔️❗:
    // RDKit✔️❗:               throw FileParseException(errout.str());
    // RDKit✔️❗:             }
    // RDKit✔️❗:           }
    // RDKit✔️❗:           attchOrds.emplace_back(idx - 1, splitToken[i + 1]);
    // RDKit✔️❗:         }
    // RDKit✔️❗:         atom->setProp(common_properties::molAttachOrderTemplate, attchOrds);
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         // this is a normal ATTCHORD
    // RDKit✔️❗:         auto ival = FileParserUtils::toInt(val);
    // RDKit✔️❗:         atom->setProp(common_properties::molAttachOrder, ival);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "CLASS") {
    // RDKit✔️❗:       atom->setProp(common_properties::molAtomClass, std::string(val));
    // RDKit✔️❗:     } else if (prop == "SEQID") {
    // RDKit✔️❗:       if (val != "0") {
    // RDKit✔️❗:         auto ival = FileParserUtils::toInt(val);
    // RDKit✔️❗:         atom->setProp(common_properties::molAtomSeqId, ival);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (prop == "SEQNAME") {
    // RDKit✔️❗:       if (val != "") {
    // RDKit✔️❗:         atom->setProp(common_properties::molAtomSeqName, std::string(val));
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:     ++token;
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000AtomProps
    // END RDKIT CPP BODY: parse_v3000_atom_props

    let _ = params;
    for token in tokens {
        let (prop, val) = split_assign_token(token).ok_or_else(|| {
            SdfReadError::Parse(format!(
                "Invalid atom property: '{token}' for atom on line {line_number}"
            ))
        })?;
        match prop {
            "CHG" => {
                let charge = parse_rdkit_int(val, false).map_err(|_| {
                    SdfReadError::Parse(format!(
                        "Cannot convert '{val}' to int on line {line_number}"
                    ))
                })?;
                spec = spec.with_formal_charge(charge as i8);
            }
            "RAD" => {
                let radical_code = parse_rdkit_int(val, false).map_err(|_| {
                    SdfReadError::Parse(format!(
                        "Cannot convert '{val}' to int on line {line_number}"
                    ))
                })?;
                let radical_electrons = match radical_code {
                    0 => 0,
                    1 => 2,
                    2 => 1,
                    3 => 2,
                    _ => {
                        return Err(SdfReadError::Parse(format!(
                            "Unrecognized RAD value {val} for atom on line {line_number}"
                        )));
                    }
                };
                spec = spec.with_radical_electrons(radical_electrons);
            }
            "MASS" => {
                let isotope = parse_rdkit_int(val, false)
                    .or_else(|_| parse_rdkit_double(val, false).map(|dv| dv.floor() as i32))
                    .map_err(|_| {
                        SdfReadError::Parse(format!(
                            "Bad value for MASS :{val} for atom on line {line_number}"
                        ))
                    })?;
                if isotope < 0 {
                    return Err(SdfReadError::Parse(format!(
                        "Bad value for MASS :{val} for atom on line {line_number}"
                    )));
                }
                spec = spec.with_isotope(isotope as u16);
            }
            "CFG" => {
                let cfg = parse_rdkit_int(val, false).map_err(|_| {
                    SdfReadError::Parse(format!(
                        "Cannot convert '{val}' to int on line {line_number}"
                    ))
                })?;
                match cfg {
                    0 => {}
                    1 | 2 | 3 => {
                        spec = spec
                            .with_mol_parity(cfg)
                            .with_prop("molParity", cfg.to_string())
                    }
                    _ => {
                        return Err(SdfReadError::Parse(format!(
                            "Unrecognized CFG value : {val} for atom on line {line_number}"
                        )));
                    }
                }
            }
            "HCOUNT" => {
                if val != "0" {
                    let mut hcount = parse_rdkit_int(val, false).map_err(|_| {
                        SdfReadError::Parse(format!(
                            "Cannot convert '{val}' to int on line {line_number}"
                        ))
                    })?;
                    if hcount == -1 {
                        hcount = 0;
                    }
                    let query = if hcount == 0 {
                        QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(0))
                    } else {
                        QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(
                            hcount as u8,
                        ))
                    };
                    let existing = spec.query().cloned();
                    spec = spec.with_query(merge_atom_query(existing, query));
                }
            }
            "UNSAT" => {
                if val == "1" {
                    let existing = spec.query().cloned();
                    spec = spec.with_query(merge_atom_query(
                        existing,
                        QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
                    ));
                }
            }
            "RBCNT" => {
                if val != "0" {
                    let mut rbcount = parse_rdkit_int(val, false).map_err(|_| {
                        SdfReadError::Parse(format!(
                            "Cannot convert '{val}' to int on line {line_number}"
                        ))
                    })?;
                    spec = spec.with_prop("molRingBondCount", rbcount.to_string());
                    if rbcount == -1 {
                        rbcount = 0;
                    } else if rbcount == -2 {
                        let existing = spec.query().cloned();
                        spec = spec.with_query(merge_atom_query(
                            existing,
                            QueryNode::predicate(AtomQueryPredicate::RingBondCount(
                                crate::search::query::QUERY_SCAN_MAGIC_VALUE,
                            )),
                        ));
                        continue;
                    } else if rbcount > 4 {
                        rbcount = 4;
                    }
                    let predicate = if rbcount == 4 {
                        AtomQueryPredicate::RingBondCountLessEqual(4)
                    } else {
                        AtomQueryPredicate::RingBondCount(rbcount as u32)
                    };
                    let existing = spec.query().cloned();
                    spec = spec
                        .with_query(merge_atom_query(existing, QueryNode::predicate(predicate)));
                }
            }
            "RGROUPS" => {
                let rlabels = parse_v3000_rgroups(val, line_number)?;
                if let Some(rlabel) = rlabels.last().copied() {
                    spec = spec
                        .with_prop("_MolFileRLabel", rlabel.to_string())
                        .with_prop("dummyLabel", format!("R{rlabel}"))
                        .with_isotope(rlabel as u16)
                        .with_query(QueryNode::predicate(AtomQueryPredicate::Any));
                }
            }
            "VAL" if val != "0" => spec = spec.with_prop("molTotValence", val.to_string()),
            "STBOX" if val != "0" => spec = spec.with_prop("molStereoCare", val.to_string()),
            "SUBST" if val != "0" => spec = spec.with_prop("molSubstCount", val.to_string()),
            "EXACHG" if val != "0" => spec = spec.with_prop("molRxnExactChange", val.to_string()),
            "INVRET" if val != "0" => {
                let inversion_flag = parse_rdkit_int(val, false).map_err(|_| {
                    SdfReadError::Parse(format!(
                        "Cannot convert '{val}' to int on line {line_number}"
                    ))
                })?;
                spec = spec.with_mol_inversion_flag(inversion_flag);
            }
            "ATTCHPT" if val != "0" => {
                if spec.prop("molAttachPoint").is_some() {
                    let message = format!("Multiple ATTCHPT values for atom on line {line_number}");
                    if params.strict_parsing {
                        return Err(SdfReadError::Parse(message));
                    }
                } else {
                    spec = spec.with_prop("molAttachPoint", val.to_string());
                }
            }
            "ATTCHORD" => {
                if val.starts_with('(') {
                    spec = spec.with_prop("molAttachOrderTemplate", val.to_string());
                } else {
                    spec = spec.with_prop("molAttachOrder", val.to_string());
                }
            }
            "CLASS" => spec = spec.with_prop("molAtomClass", val.to_string()),
            "SEQID" if val != "0" => spec = spec.with_prop("molAtomSeqId", val.to_string()),
            "SEQNAME" if !val.is_empty() => {
                spec = spec.with_prop("molAtomSeqName", val.to_string())
            }
            _ => {}
        }
    }
    Ok(spec)
}

fn parse_v3000_rgroups(text: &str, line_number: usize) -> Result<Vec<u32>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_rgroups
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000RGroups
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
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000RGroups
    // END RDKIT CPP BODY: parse_v3000_rgroups

    if !text.starts_with('(') || !text.ends_with(')') {
        return Err(SdfReadError::Parse(format!(
            "Bad RGROUPS specification '{text}' on line {line_number}. Missing parens."
        )));
    }
    let resid = &text[1..text.len() - 1];
    let split_token: Vec<_> = resid.split(' ').collect();
    if split_token.is_empty() {
        return Err(SdfReadError::Parse(format!(
            "Bad RGROUPS specification '{text}' on line {line_number}. Missing values."
        )));
    }
    let n_rs = parse_rdkit_unsigned(split_token[0], true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line{}",
            split_token[0], line_number
        ))
    })? as usize;
    if split_token.len() < n_rs + 1 {
        return Err(SdfReadError::Parse(format!(
            "Bad RGROUPS specification '{text}' on line {line_number}. Not enough values."
        )));
    }
    split_token
        .iter()
        .skip(1)
        .take(n_rs)
        .map(|token| {
            parse_rdkit_unsigned(token, true).map_err(|_| {
                SdfReadError::Parse(format!(
                    "Cannot convert '{token}' to int on line{line_number}"
                ))
            })
        })
        .collect()
}

fn split_assign_token(token: &str) -> Option<(&str, &str)> {
    let (prop, val) = token.split_once('=')?;
    if prop.is_empty() {
        return None;
    }
    Some((prop, val))
}

fn v3000_bond_spec_for_type(
    begin: AtomId,
    end: AtomId,
    bond_type: u32,
) -> (BondSpec, Option<QueryNode<BondQueryPredicate>>) {
    match bond_type {
        1 => (BondSpec::new(begin, end, BondOrder::Single), None),
        2 => (BondSpec::new(begin, end, BondOrder::Double), None),
        3 => (BondSpec::new(begin, end, BondOrder::Triple), None),
        4 => (
            BondSpec::new(begin, end, BondOrder::Aromatic).with_aromatic(true),
            None,
        ),
        9 => (BondSpec::new(begin, end, BondOrder::Dative), None),
        10 => (BondSpec::new(begin, end, BondOrder::Hydrogen), None),
        0 => (BondSpec::new(begin, end, BondOrder::Unspecified), None),
        5 => (
            BondSpec::new(begin, end, BondOrder::Unspecified).with_prop("_MolFileBondQuery", "1"),
            Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ]))),
        ),
        6 => (
            BondSpec::new(begin, end, BondOrder::Unspecified).with_prop("_MolFileBondQuery", "1"),
            Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ]))),
        ),
        7 => (
            BondSpec::new(begin, end, BondOrder::Unspecified).with_prop("_MolFileBondQuery", "1"),
            Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Double,
                BondOrder::Aromatic,
            ]))),
        ),
        8 => (
            BondSpec::new(begin, end, BondOrder::Unspecified),
            Some(QueryNode::predicate(BondQueryPredicate::Any)),
        ),
        other => (
            BondSpec::new(begin, end, BondOrder::Unspecified),
            Some(QueryNode::predicate(BondQueryPredicate::MolFileQueryCode(
                other,
            ))),
        ),
    }
}

fn parse_v3000_bond_block<R: BufRead>(
    reader: &mut R,
    bond_count: usize,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<Vec<V3000BondLine>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_bond_block
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000BondBlock
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(nBonds > 0, "bad bond count");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:
    // RDKit✔️❗:   auto inl = getV3000Line(inStream, line);
    // RDKit✔️❗:   std::string_view tempStr = inl;
    // RDKit✔️❗:   if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN BOND") {
    // RDKit✔️❗:     throw FileParseException("BEGIN BOND line not found");
    // RDKit✔️❗:   }
    // RDKit✔️❗:   for (unsigned int i = 0; i < nBonds; ++i) {
    // RDKit✔️❗:     inl = getV3000Line(inStream, line);
    // RDKit✔️❗:     tempStr = inl;
    // RDKit✔️❗:     tempStr = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:     std::vector<std::string_view> splitLine;
    // RDKit✔️❗:     tokenizeV3000Line(tempStr, splitLine);
    // RDKit✔️❗:     if (splitLine.size() < 4) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "bond line " << line << " is too short";
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     Bond *bond;
    // RDKit✔️❗:     unsigned int bondIdx = 0;
    // RDKit✔️❗:     std::from_chars(splitLine[0].data(),
    // RDKit✔️❗:                     splitLine[0].data() + splitLine[0].size(), bondIdx);
    // RDKit✔️❗:     unsigned int bType = 0;
    // RDKit✔️❗:     std::from_chars(splitLine[1].data(),
    // RDKit✔️❗:                     splitLine[1].data() + splitLine[1].size(), bType);
    // RDKit✔️❗:     unsigned int a1Idx = 0;
    // RDKit✔️❗:     std::from_chars(splitLine[2].data(),
    // RDKit✔️❗:                     splitLine[2].data() + splitLine[2].size(), a1Idx);
    // RDKit✔️❗:     unsigned int a2Idx = 0;
    // RDKit✔️❗:     std::from_chars(splitLine[3].data(),
    // RDKit✔️❗:                     splitLine[3].data() + splitLine[3].size(), a2Idx);
    // RDKit✔️❗:
    // RDKit✔️❗:     switch (bType) {
    // RDKit✔️❗:       case 1:
    // RDKit✔️❗:         bond = new Bond(Bond::SINGLE);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 2:
    // RDKit✔️❗:         bond = new Bond(Bond::DOUBLE);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 3:
    // RDKit✔️❗:         bond = new Bond(Bond::TRIPLE);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 4:
    // RDKit✔️❗:         bond = new Bond(Bond::AROMATIC);
    // RDKit✔️❗:         bond->setIsAromatic(true);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 9:
    // RDKit✔️❗:         bond = new Bond(Bond::DATIVE);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 10:
    // RDKit✔️❗:         bond = new Bond(Bond::HYDROGEN);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 0:
    // RDKit✔️❗:         bond = new Bond(Bond::UNSPECIFIED);
    // RDKit✔️❗:         BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:             << "bond with order 0 found on line " << line
    // RDKit✔️❗:             << ". This is not part of the MDL specification." << std::endl;
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       default:
    // RDKit✔️❗:         // it's a query bond of some type
    // RDKit✔️❗:         bond = new QueryBond;
    // RDKit✔️❗:         if (bType == 8) {
    // RDKit✔️❗:           BOND_NULL_QUERY *q;
    // RDKit✔️❗:           q = makeBondNullQuery();
    // RDKit✔️❗:           bond->setQuery(q);
    // RDKit✔️❗:         } else if (bType == 5) {
    // RDKit✔️❗:           bond->setQuery(makeSingleOrDoubleBondQuery());
    // RDKit✔️❗:           bond->setProp(common_properties::_MolFileBondQuery, 1);
    // RDKit✔️❗:         } else if (bType == 6) {
    // RDKit✔️❗:           bond->setQuery(makeSingleOrAromaticBondQuery());
    // RDKit✔️❗:           bond->setProp(common_properties::_MolFileBondQuery, 1);
    // RDKit✔️❗:         } else if (bType == 7) {
    // RDKit✔️❗:           bond->setQuery(makeDoubleOrAromaticBondQuery());
    // RDKit✔️❗:           bond->setProp(common_properties::_MolFileBondQuery, 1);
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOND_NULL_QUERY *q;
    // RDKit✔️❗:           q = makeBondNullQuery();
    // RDKit✔️❗:           bond->setQuery(q);
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:               << "unrecognized query bond type, " << bType << ", found on line "
    // RDKit✔️❗:               << line << ". Using an \"any\" query." << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:         break;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     bond->setProp(common_properties::_MolFileBondType, bType);
    // RDKit✔️❗:
    // RDKit✔️❗:     // additional bond properties:
    // RDKit✔️❗:     unsigned int lPos = 4;
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     while (lPos < splitLine.size()) {
    // RDKit✔️❗:       std::string prop;
    // RDKit✔️❗:       std::string_view val;
    // RDKit✔️❗:       if (!splitAssignToken(splitLine[lPos], prop, val)) {
    // RDKit✔️❗:         errout << "bad bond property '" << splitLine[lPos] << "' on line "
    // RDKit✔️❗:                << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (prop == "CFG") {
    // RDKit✔️❗:         unsigned int cfg = 0;
    // RDKit✔️❗:         std::from_chars(val.data(), val.data() + val.size(), cfg);
    // RDKit✔️❗:         switch (cfg) {
    // RDKit✔️❗:           case 0:
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           case 1:
    // RDKit✔️❗:             bond->setBondDir(Bond::BEGINWEDGE);
    // RDKit✔️❗:             chiralityPossible = true;
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           case 2:
    // RDKit✔️❗:             if (bType == 1) {
    // RDKit✔️❗:               bond->setBondDir(Bond::UNKNOWN);
    // RDKit✔️❗:             } else if (bType == 2) {
    // RDKit✔️❗:               bond->setBondDir(Bond::EITHERDOUBLE);
    // RDKit✔️❗:               bond->setStereo(Bond::STEREOANY);
    // RDKit✔️❗:             }
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           case 3:
    // RDKit✔️❗:             bond->setBondDir(Bond::BEGINDASH);
    // RDKit✔️❗:             chiralityPossible = true;
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           default:
    // RDKit✔️❗:             errout << "bad bond CFG " << val << "' on line " << line;
    // RDKit✔️❗:             throw FileParseException(errout.str());
    // RDKit✔️❗:         }
    // RDKit✔️❗:         bond->setProp(common_properties::_MolFileBondCfg, cfg);
    // RDKit✔️❗:       } else if (prop == "TOPO") {
    // RDKit✔️❗:         if (val != "0") {
    // RDKit✔️❗:           if (!bond->hasQuery()) {
    // RDKit✔️❗:             auto *qBond = new QueryBond(*bond);
    // RDKit✔️❗:             delete bond;
    // RDKit✔️❗:             bond = qBond;
    // RDKit✔️❗:           }
    // RDKit✔️❗:           BOND_EQUALS_QUERY *q = makeBondIsInRingQuery();
    // RDKit✔️❗:           if (val == "1") {
    // RDKit✔️❗:             // nothing
    // RDKit✔️❗:           } else if (val == "2") {
    // RDKit✔️❗:             q->setNegation(true);
    // RDKit✔️❗:           } else {
    // RDKit✔️❗:             errout << "bad bond TOPO " << val << "' on line " << line;
    // RDKit✔️❗:             throw FileParseException(errout.str());
    // RDKit✔️❗:           }
    // RDKit✔️❗:           bond->expandQuery(q);
    // RDKit✔️❗:         }
    // RDKit✔️❗:       } else if (prop == "RXCTR") {
    // RDKit✔️❗:         int reactStatus = FileParserUtils::toInt(val);
    // RDKit✔️❗:         bond->setProp(common_properties::molReactStatus, reactStatus);
    // RDKit✔️❗:       } else if (prop == "STBOX") {
    // RDKit✔️❗:         bond->setProp(common_properties::molStereoCare, std::string(val));
    // RDKit✔️❗:       } else if (prop == "ENDPTS") {
    // RDKit✔️❗:         bond->setProp(common_properties::_MolFileBondEndPts, std::string(val));
    // RDKit✔️❗:       } else if (prop == "ATTACH") {
    // RDKit✔️❗:         bond->setProp(common_properties::_MolFileBondAttach, std::string(val));
    // RDKit✔️❗:       }
    // RDKit✔️❗:       ++lPos;
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     bond->setBeginAtomIdx(mol->getAtomWithBookmark(a1Idx)->getIdx());
    // RDKit✔️❗:     bond->setEndAtomIdx(mol->getAtomWithBookmark(a2Idx)->getIdx());
    // RDKit✔️❗:     mol->addBond(bond, true);
    // RDKit✔️❗:     mol->setBondBookmark(bond, bondIdx);
    // RDKit✔️❗:
    // RDKit✔️❗:     // set the stereoCare property on the bond if it's not set already and
    // RDKit✔️❗:     // both the beginning and end atoms have it set:
    // RDKit✔️❗:     int care1 = 0;
    // RDKit✔️❗:     int care2 = 0;
    // RDKit✔️❗:     if (!bond->hasProp(common_properties::molStereoCare) &&
    // RDKit✔️❗:         mol->getAtomWithIdx(bond->getBeginAtomIdx())
    // RDKit✔️❗:             ->getPropIfPresent(common_properties::molStereoCare, care1) &&
    // RDKit✔️❗:         mol->getAtomWithIdx(bond->getEndAtomIdx())
    // RDKit✔️❗:             ->getPropIfPresent(common_properties::molStereoCare, care2)) {
    // RDKit✔️❗:       if (care1 == care2) {
    // RDKit✔️❗:         bond->setProp(common_properties::molStereoCare, care1);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   inl = getV3000Line(inStream, line);
    // RDKit✔️❗:   tempStr = inl;
    // RDKit✔️❗:   if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END BOND") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "END BOND line not found at line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000BondBlock
    // END RDKIT CPP BODY: parse_v3000_bond_block

    let begin = get_v3000_line(reader, line_number, params)?;
    if begin.len() < 10 || !begin.starts_with("BEGIN BOND") {
        return Err(SdfReadError::Parse("BEGIN BOND line not found".to_string()));
    }

    let mut bonds = Vec::with_capacity(bond_count);
    for _ in 0..bond_count {
        let temp_str = get_v3000_line(reader, line_number, params)?;
        let trimmed = temp_str.trim().to_string();
        let tokens = tokenize_v3000_line(&trimmed, *line_number)?;
        if tokens.len() < 4 {
            return Err(SdfReadError::Parse(format!(
                "bond line {} is too short",
                *line_number
            )));
        }
        let bond_idx = parse_rdkit_unsigned(&tokens[0], false).unwrap_or(0);
        let bond_type = parse_rdkit_unsigned(&tokens[1], false).unwrap_or(0);
        let begin_mol_idx = parse_rdkit_unsigned(&tokens[2], false).unwrap_or(0);
        let end_mol_idx = parse_rdkit_unsigned(&tokens[3], false).unwrap_or(0);
        let (mut spec, query) = v3000_bond_spec_for_type(AtomId::new(0), AtomId::new(0), bond_type);
        if let Some(query) = query {
            spec = spec.with_query(query);
        }
        spec = parse_v3000_bond_props(&tokens[4..], *line_number, params, spec, bond_type)?;
        bonds.push(V3000BondLine {
            line_number: *line_number,
            mol_idx: bond_idx,
            tokens,
            begin_mol_idx,
            end_mol_idx,
            spec,
            molfile_bond_type: bond_type,
        });
    }

    let end = get_v3000_line(reader, line_number, params)?;
    if end.len() < 8 || !end.starts_with("END BOND") {
        return Err(SdfReadError::Parse(format!(
            "END BOND line not found at line {}",
            *line_number
        )));
    }
    Ok(bonds)
}

fn parse_v3000_bond_props(
    tokens: &[String],
    line_number: usize,
    params: SdfReadParams,
    mut spec: BondSpec,
    bond_type: u32,
) -> Result<BondSpec, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_bond_props
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000BondBlock
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(nBonds > 0, "bad bond count");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:
    // RDKit✔️❗:   auto inl = getV3000Line(inStream, line);
    // RDKit✔️❗:   std::string_view tempStr = inl;
    // RDKit✔️❗:   if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN BOND") {
    // RDKit✔️❗:     throw FileParseException("BEGIN BOND line not found");
    // RDKit✔️❗:   }
    // RDKit✔️❗:   for (unsigned int i = 0; i < nBonds; ++i) {
    // RDKit✔️❗:     inl = getV3000Line(inStream, line);
    // RDKit✔️❗:     tempStr = inl;
    // RDKit✔️❗:     tempStr = FileParserUtils::strip(tempStr);
    // RDKit✔️❗:     std::vector<std::string_view> splitLine;
    // RDKit✔️❗:     tokenizeV3000Line(tempStr, splitLine);
    // RDKit✔️❗:     if (splitLine.size() < 4) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "bond line " << line << " is too short";
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:     Bond *bond;
    // RDKit✔️❗:     unsigned int bondIdx = 0;
    // RDKit✔️❗:     std::from_chars(splitLine[0].data(),
    // RDKit✔️❗:                     splitLine[0].data() + splitLine[0].size(), bondIdx);
    // RDKit✔️❗:     unsigned int bType = 0;
    // RDKit✔️❗:     std::from_chars(splitLine[1].data(),
    // RDKit✔️❗:                     splitLine[1].data() + splitLine[1].size(), bType);
    // RDKit✔️❗:     unsigned int a1Idx = 0;
    // RDKit✔️❗:     std::from_chars(splitLine[2].data(),
    // RDKit✔️❗:                     splitLine[2].data() + splitLine[2].size(), a1Idx);
    // RDKit✔️❗:     unsigned int a2Idx = 0;
    // RDKit✔️❗:     std::from_chars(splitLine[3].data(),
    // RDKit✔️❗:                     splitLine[3].data() + splitLine[3].size(), a2Idx);
    // RDKit✔️❗:
    // RDKit✔️❗:     switch (bType) {
    // RDKit✔️❗:       case 1:
    // RDKit✔️❗:         bond = new Bond(Bond::SINGLE);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 2:
    // RDKit✔️❗:         bond = new Bond(Bond::DOUBLE);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 3:
    // RDKit✔️❗:         bond = new Bond(Bond::TRIPLE);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 4:
    // RDKit✔️❗:         bond = new Bond(Bond::AROMATIC);
    // RDKit✔️❗:         bond->setIsAromatic(true);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 9:
    // RDKit✔️❗:         bond = new Bond(Bond::DATIVE);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 10:
    // RDKit✔️❗:         bond = new Bond(Bond::HYDROGEN);
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       case 0:
    // RDKit✔️❗:         bond = new Bond(Bond::UNSPECIFIED);
    // RDKit✔️❗:         BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:             << "bond with order 0 found on line " << line
    // RDKit✔️❗:             << ". This is not part of the MDL specification." << std::endl;
    // RDKit✔️❗:         break;
    // RDKit✔️❗:       default:
    // RDKit✔️❗:         // it's a query bond of some type
    // RDKit✔️❗:         bond = new QueryBond;
    // RDKit✔️❗:         if (bType == 8) {
    // RDKit✔️❗:           BOND_NULL_QUERY *q;
    // RDKit✔️❗:           q = makeBondNullQuery();
    // RDKit✔️❗:           bond->setQuery(q);
    // RDKit✔️❗:         } else if (bType == 5) {
    // RDKit✔️❗:           bond->setQuery(makeSingleOrDoubleBondQuery());
    // RDKit✔️❗:           bond->setProp(common_properties::_MolFileBondQuery, 1);
    // RDKit✔️❗:         } else if (bType == 6) {
    // RDKit✔️❗:           bond->setQuery(makeSingleOrAromaticBondQuery());
    // RDKit✔️❗:           bond->setProp(common_properties::_MolFileBondQuery, 1);
    // RDKit✔️❗:         } else if (bType == 7) {
    // RDKit✔️❗:           bond->setQuery(makeDoubleOrAromaticBondQuery());
    // RDKit✔️❗:           bond->setProp(common_properties::_MolFileBondQuery, 1);
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOND_NULL_QUERY *q;
    // RDKit✔️❗:           q = makeBondNullQuery();
    // RDKit✔️❗:           bond->setQuery(q);
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:               << "unrecognized query bond type, " << bType << ", found on line "
    // RDKit✔️❗:               << line << ". Using an \"any\" query." << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:         break;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     bond->setProp(common_properties::_MolFileBondType, bType);
    // RDKit✔️❗:
    // RDKit✔️❗:     // additional bond properties:
    // RDKit✔️❗:     unsigned int lPos = 4;
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     while (lPos < splitLine.size()) {
    // RDKit✔️❗:       std::string prop;
    // RDKit✔️❗:       std::string_view val;
    // RDKit✔️❗:       if (!splitAssignToken(splitLine[lPos], prop, val)) {
    // RDKit✔️❗:         errout << "bad bond property '" << splitLine[lPos] << "' on line "
    // RDKit✔️❗:                << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (prop == "CFG") {
    // RDKit✔️❗:         unsigned int cfg = 0;
    // RDKit✔️❗:         std::from_chars(val.data(), val.data() + val.size(), cfg);
    // RDKit✔️❗:         switch (cfg) {
    // RDKit✔️❗:           case 0:
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           case 1:
    // RDKit✔️❗:             bond->setBondDir(Bond::BEGINWEDGE);
    // RDKit✔️❗:             chiralityPossible = true;
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           case 2:
    // RDKit✔️❗:             if (bType == 1) {
    // RDKit✔️❗:               bond->setBondDir(Bond::UNKNOWN);
    // RDKit✔️❗:             } else if (bType == 2) {
    // RDKit✔️❗:               bond->setBondDir(Bond::EITHERDOUBLE);
    // RDKit✔️❗:               bond->setStereo(Bond::STEREOANY);
    // RDKit✔️❗:             }
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           case 3:
    // RDKit✔️❗:             bond->setBondDir(Bond::BEGINDASH);
    // RDKit✔️❗:             chiralityPossible = true;
    // RDKit✔️❗:             break;
    // RDKit✔️❗:           default:
    // RDKit✔️❗:             errout << "bad bond CFG " << val << "' on line " << line;
    // RDKit✔️❗:             throw FileParseException(errout.str());
    // RDKit✔️❗:         }
    // RDKit✔️❗:         bond->setProp(common_properties::_MolFileBondCfg, cfg);
    // RDKit✔️❗:       } else if (prop == "TOPO") {
    // RDKit✔️❗:         if (val != "0") {
    // RDKit✔️❗:           if (!bond->hasQuery()) {
    // RDKit✔️❗:             auto *qBond = new QueryBond(*bond);
    // RDKit✔️❗:             delete bond;
    // RDKit✔️❗:             bond = qBond;
    // RDKit✔️❗:           }
    // RDKit✔️❗:           BOND_EQUALS_QUERY *q = makeBondIsInRingQuery();
    // RDKit✔️❗:           if (val == "1") {
    // RDKit✔️❗:             // nothing
    // RDKit✔️❗:           } else if (val == "2") {
    // RDKit✔️❗:             q->setNegation(true);
    // RDKit✔️❗:           } else {
    // RDKit✔️❗:             errout << "bad bond TOPO " << val << "' on line " << line;
    // RDKit✔️❗:             throw FileParseException(errout.str());
    // RDKit✔️❗:           }
    // RDKit✔️❗:           bond->expandQuery(q);
    // RDKit✔️❗:         }
    // RDKit✔️❗:       } else if (prop == "RXCTR") {
    // RDKit✔️❗:         int reactStatus = FileParserUtils::toInt(val);
    // RDKit✔️❗:         bond->setProp(common_properties::molReactStatus, reactStatus);
    // RDKit✔️❗:       } else if (prop == "STBOX") {
    // RDKit✔️❗:         bond->setProp(common_properties::molStereoCare, std::string(val));
    // RDKit✔️❗:       } else if (prop == "ENDPTS") {
    // RDKit✔️❗:         bond->setProp(common_properties::_MolFileBondEndPts, std::string(val));
    // RDKit✔️❗:       } else if (prop == "ATTACH") {
    // RDKit✔️❗:         bond->setProp(common_properties::_MolFileBondAttach, std::string(val));
    // RDKit✔️❗:       }
    // RDKit✔️❗:       ++lPos;
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     bond->setBeginAtomIdx(mol->getAtomWithBookmark(a1Idx)->getIdx());
    // RDKit✔️❗:     bond->setEndAtomIdx(mol->getAtomWithBookmark(a2Idx)->getIdx());
    // RDKit✔️❗:     mol->addBond(bond, true);
    // RDKit✔️❗:     mol->setBondBookmark(bond, bondIdx);
    // RDKit✔️❗:
    // RDKit✔️❗:     // set the stereoCare property on the bond if it's not set already and
    // RDKit✔️❗:     // both the beginning and end atoms have it set:
    // RDKit✔️❗:     int care1 = 0;
    // RDKit✔️❗:     int care2 = 0;
    // RDKit✔️❗:     if (!bond->hasProp(common_properties::molStereoCare) &&
    // RDKit✔️❗:         mol->getAtomWithIdx(bond->getBeginAtomIdx())
    // RDKit✔️❗:             ->getPropIfPresent(common_properties::molStereoCare, care1) &&
    // RDKit✔️❗:         mol->getAtomWithIdx(bond->getEndAtomIdx())
    // RDKit✔️❗:             ->getPropIfPresent(common_properties::molStereoCare, care2)) {
    // RDKit✔️❗:       if (care1 == care2) {
    // RDKit✔️❗:         bond->setProp(common_properties::molStereoCare, care1);
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   inl = getV3000Line(inStream, line);
    // RDKit✔️❗:   tempStr = inl;
    // RDKit✔️❗:   if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END BOND") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "END BOND line not found at line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ParseV3000BondBlock
    // END RDKIT CPP BODY: parse_v3000_bond_props

    let _ = params;
    for token in tokens {
        let (prop, val) = split_assign_token(token).ok_or_else(|| {
            SdfReadError::Parse(format!("bad bond property '{token}' on line {line_number}"))
        })?;
        match prop {
            "CFG" => {
                let cfg = parse_rdkit_unsigned(val, false).unwrap_or(0);
                match cfg {
                    0 => {}
                    1 => spec = spec.with_direction(BondDirection::BeginWedge),
                    2 => {
                        if bond_type == 1 {
                            spec = spec.with_direction(BondDirection::Unknown);
                        } else if bond_type == 2 {
                            spec = spec
                                .with_direction(BondDirection::EitherDouble)
                                .with_stereo(BondStereo::Any);
                        }
                    }
                    3 => spec = spec.with_direction(BondDirection::BeginDash),
                    _ => {
                        return Err(SdfReadError::Parse(format!(
                            "bad bond CFG {val}' on line {line_number}"
                        )));
                    }
                }
            }
            "TOPO" => {
                if val != "0" {
                    let topology_query = match val {
                        "1" => QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
                        "2" => QueryNode::predicate(BondQueryPredicate::IsInRing(false)),
                        _ => {
                            return Err(SdfReadError::Parse(format!(
                                "bad bond TOPO {val}' on line {line_number}"
                            )));
                        }
                    };
                    let query = match spec.query().cloned() {
                        Some(existing) => QueryNode::and(vec![existing, topology_query]),
                        None => topology_query,
                    };
                    spec = spec.with_query(query);
                }
            }
            "RXCTR" => spec = spec.with_prop("molReactStatus", val.to_string()),
            "STBOX" => spec = spec.with_prop("molStereoCare", val.to_string()),
            "ENDPTS" => spec = spec.with_prop("_MolFileBondEndPts", val.to_string()),
            "ATTACH" => spec = spec.with_prop("_MolFileBondAttach", val.to_string()),
            _ => {}
        }
    }
    Ok(spec)
}

fn parse_v3000_sgroup_block<R: BufRead>(
    reader: &mut R,
    sgroup_count: usize,
    line_number: &mut usize,
    params: SdfReadParams,
    atom_id_by_mol_idx: &BTreeMap<u32, AtomId>,
    bond_id_by_mol_idx: &BTreeMap<u32, BondId>,
) -> Result<Vec<SubstanceGroup>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_sgroup_block
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: std::string ParseV3000SGroupsBlock
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "no stream");
    // RDKit✔️❗:   PRECONDITION(mol, "no molecule");
    // RDKit✔️❗:   unsigned int defaultLineNum = 0;
    // RDKit✔️❗:   std::string defaultString;
    // RDKit✔️❗:
    // RDKit✔️❗:   // SGroups may be written in unsorted ID order, according to spec, so we will
    // RDKit✔️❗:   // temporarily store them in a map before adding them to the mol
    // RDKit✔️❗:   IDX_TO_SGROUP_MAP sGroupMap;
    // RDKit✔️❗:
    // RDKit✔️❗:   std::unordered_map<std::string, std::stringstream> defaultLabels;
    // RDKit✔️❗:
    // RDKit✔️❗:   auto tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit✔️❗:
    // RDKit✔️❗:   // Store defaults
    // RDKit✔️❗:   if (tempStr.substr(0, 7) == "DEFAULT" && tempStr.length() > 8) {
    // RDKit✔️❗:     defaultString = tempStr.substr(7);
    // RDKit✔️❗:     defaultLineNum = line;
    // RDKit✔️❗:     boost::trim_right(defaultString);
    // RDKit✔️❗:     tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit✔️❗:     boost::trim_right(tempStr);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   for (unsigned int si = 0; si < nSgroups; ++si) {
    // RDKit✔️❗:     unsigned int sequenceId;
    // RDKit✔️❗:     unsigned int externalId;
    // RDKit✔️❗:     std::string type;
    // RDKit✔️❗:
    // RDKit✔️❗:     std::stringstream lineStream(tempStr);
    // RDKit✔️❗:     lineStream >> sequenceId;
    // RDKit✔️❗:     lineStream >> type;
    // RDKit✔️❗:     lineStream >> externalId;
    // RDKit✔️❗:
    // RDKit✔️❗:     std::set<std::string> parsedLabels;
    // RDKit✔️❗:     if (strictParsing && !SubstanceGroupChecks::isValidType(type)) {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Unsupported SGroup type '" << type << "' on line " << line;
    // RDKit✔️❗:       throw MolFileUnhandledFeatureException(errout.str());
    // RDKit✔️❗:     } else if (!strictParsing &&
    // RDKit✔️❗:                nSgroups == std::numeric_limits<unsigned int>::max() &&
    // RDKit✔️❗:                lineStream.fail()) {
    // RDKit✔️❗:       // something went wrong and we didn't know how many SGroups to expect, and
    // RDKit✔️❗:       // now we have seen something that doesn't look like an SGroup start.
    // RDKit✔️❗:       // So we assume we're done.
    // RDKit✔️❗:       nSgroups = 0;
    // RDKit✔️❗:       break;
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     SubstanceGroup sgroup(mol, type);
    // RDKit✔️❗:     STR_VECT dataFields;
    // RDKit✔️❗:
    // RDKit✔️❗:     sgroup.setProp<unsigned int>("index", sequenceId);
    // RDKit✔️❗:     if (externalId > 0) {
    // RDKit✔️❗:       if (!SubstanceGroupChecks::isSubstanceGroupIdFree(*mol, externalId)) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Existing SGroup ID '" << externalId
    // RDKit✔️❗:                << "' assigned to a second SGroup on line " << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:           sgroup.setIsValid(false);
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       sgroup.setProp<unsigned int>("ID", externalId);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     while (sgroup.getIsValid() && !lineStream.eof() && !lineStream.fail()) {
    // RDKit✔️❗:       char spacer;
    // RDKit✔️❗:       std::string label;
    // RDKit✔️❗:
    // RDKit✔️❗:       lineStream.get(spacer);
    // RDKit✔️❗:       if (lineStream.gcount() == 0) {
    // RDKit✔️❗:         continue;
    // RDKit✔️❗:       } else if (spacer != ' ') {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Found character '" << spacer
    // RDKit✔️❗:                << "' when expecting a separator (space) on line " << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:           sgroup.setIsValid(false);
    // RDKit✔️❗:           continue;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       std::getline(lineStream, label, '=');
    // RDKit✔️❗:       if (label.empty()) {
    // RDKit✔️❗:         continue;
    // RDKit✔️❗:       }
    // RDKit✔️❗:       ParseV3000ParseLabel(label, lineStream, dataFields, line, sgroup,
    // RDKit✔️❗:                            nSgroups, mol, strictParsing);
    // RDKit✔️❗:       parsedLabels.insert(label);
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     // Process defaults
    // RDKit✔️❗:     lineStream.clear();
    // RDKit✔️❗:     lineStream.str(defaultString);
    // RDKit✔️❗:     while (sgroup.getIsValid() && !lineStream.eof() && !lineStream.fail()) {
    // RDKit✔️❗:       char spacer;
    // RDKit✔️❗:       std::string label;
    // RDKit✔️❗:
    // RDKit✔️❗:       lineStream.get(spacer);
    // RDKit✔️❗:       if (lineStream.gcount() == 0) {
    // RDKit✔️❗:         continue;
    // RDKit✔️❗:       } else if (spacer != ' ') {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Found character '" << spacer
    // RDKit✔️❗:                << "' when expecting a separator (space) in DEFAULTS on line "
    // RDKit✔️❗:                << defaultLineNum;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:           sgroup.setIsValid(false);
    // RDKit✔️❗:           continue;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       std::getline(lineStream, label, '=');
    // RDKit✔️❗:       if (label.empty()) {
    // RDKit✔️❗:         continue;
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (std::find(parsedLabels.begin(), parsedLabels.end(), label) ==
    // RDKit✔️❗:           parsedLabels.end()) {
    // RDKit✔️❗:         ParseV3000ParseLabel(label, lineStream, dataFields, defaultLineNum,
    // RDKit✔️❗:                              sgroup, nSgroups, mol, strictParsing);
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         spacer = lineStream.peek();
    // RDKit✔️❗:         if (spacer == ' ') {
    // RDKit✔️❗:           std::ostringstream errout;
    // RDKit✔️❗:           errout << "Found unexpected whitespace at DEFAULT label " << label;
    // RDKit✔️❗:           if (strictParsing) {
    // RDKit✔️❗:             throw FileParseException(errout.str());
    // RDKit✔️❗:           } else {
    // RDKit✔️❗:             BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:             sgroup.setIsValid(false);
    // RDKit✔️❗:             continue;
    // RDKit✔️❗:           }
    // RDKit✔️❗:         } else if (spacer == '(') {
    // RDKit✔️❗:           std::getline(lineStream, label, ')');
    // RDKit✔️❗:           lineStream.get(spacer);
    // RDKit✔️❗:         } else if (spacer == '"') {
    // RDKit✔️❗:           lineStream.get(spacer);
    // RDKit✔️❗:           std::getline(lineStream, label, '"');
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           std::getline(lineStream, label, ' ');
    // RDKit✔️❗:           lineStream.putback(' ');
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:     }
    // RDKit✔️❗:
    // RDKit✔️❗:     sgroup.setProp("DATAFIELDS", dataFields);
    // RDKit✔️❗:     sGroupMap.emplace(sequenceId, sgroup);
    // RDKit✔️❗:
    // RDKit✔️❗:     tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit✔️❗:     boost::trim_right(tempStr);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (sGroupMap.size() != nSgroups) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Found " << sGroupMap.size() << " SGroups when " << nSgroups
    // RDKit✔️❗:            << " were expected." << std::endl;
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   // SGroups successfully parsed, now add them to the molecule
    // RDKit✔️❗:   for (const auto &sg : sGroupMap) {
    // RDKit✔️❗:     if (sg.second.getIsValid()) {
    // RDKit✔️❗:       addSubstanceGroup(*mol, sg.second);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << "SGroup " << sg.first
    // RDKit✔️❗:                               << " is invalid and will be ignored" << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return tempStr;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupParsing.cpp :: std::string ParseV3000SGroupsBlock
    // END RDKIT CPP BODY: parse_v3000_sgroup_block

    let mut sgroups = BTreeMap::<u32, SubstanceGroup>::new();
    let mut parent_by_sequence_id = BTreeMap::<u32, u32>::new();
    let mut temp_str = get_v3000_line(reader, line_number, params)?;
    let mut default_tokens = Vec::new();

    if temp_str.starts_with("DEFAULT") && temp_str.len() > 8 {
        default_tokens = tokenize_v3000_line(temp_str[7..].trim_end(), *line_number)?;
        temp_str = get_v3000_line(reader, line_number, params)?;
    }

    for _ in 0..sgroup_count {
        let trimmed = temp_str.trim_end().to_string();
        let (sequence_token, typ, external_token, label_tokens) =
            split_v3000_sgroup_line(&trimmed, *line_number)?;
        if sequence_token.is_empty() || typ.is_empty() || external_token.is_empty() {
            return Err(SdfReadError::Parse(format!(
                "Found {} SGroups when {sgroup_count} were expected.",
                sgroups.len()
            )));
        }
        let sequence_id = parse_rdkit_unsigned(sequence_token, false).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to unsigned int on line {}",
                sequence_token, *line_number
            ))
        })?;
        if params.strict_parsing && !is_valid_rdkit_sgroup_type(typ) {
            return Err(SdfReadError::Parse(format!(
                "Unsupported SGroup type '{typ}' on line {}",
                *line_number
            )));
        }
        let external_id = parse_rdkit_unsigned(external_token, false).unwrap_or(0);
        let mut sgroup = SubstanceGroup::new(
            SubstanceGroupId::new(sgroups.len()),
            sgroup_kind_from_rdkit_type(typ),
        );
        sgroup.set_rdkit_sequence_id(sequence_id);
        sgroup.set_prop("TYPE", typ);
        if external_id > 0 {
            sgroup.set_external_id(external_id);
        }

        let mut parsed_labels = BTreeMap::<String, ()>::new();
        for token in &label_tokens {
            let Some((label, value)) = split_assign_token(token) else {
                return Err(SdfReadError::Parse(format!(
                    "Found character when expecting a separator (space) on line {}",
                    *line_number
                )));
            };
            parse_v3000_sgroup_label(
                label,
                value,
                *line_number,
                params,
                &mut sgroup,
                &mut parent_by_sequence_id,
                atom_id_by_mol_idx,
                bond_id_by_mol_idx,
            )?;
            parsed_labels.insert(label.to_string(), ());
        }
        let default_label_tokens = merge_v3000_assign_tokens(&default_tokens);
        for token in &default_label_tokens {
            let Some((label, value)) = split_assign_token(token) else {
                return Err(SdfReadError::Parse(format!(
                    "Found character when expecting a separator (space) in DEFAULTS on line {}",
                    *line_number
                )));
            };
            if !parsed_labels.contains_key(label) {
                parse_v3000_sgroup_label(
                    label,
                    value,
                    *line_number,
                    params,
                    &mut sgroup,
                    &mut parent_by_sequence_id,
                    atom_id_by_mol_idx,
                    bond_id_by_mol_idx,
                )?;
            }
        }

        sgroups.insert(sequence_id, sgroup);
        temp_str = get_v3000_line(reader, line_number, params)?;
        if temp_str.to_ascii_uppercase().starts_with("END SGROUP") {
            break;
        }
    }

    if !temp_str.to_ascii_uppercase().starts_with("END SGROUP") {
        return Err(SdfReadError::Parse(format!(
            "END SGROUP line not found on line {}",
            *line_number
        )));
    }
    if sgroups.len() != sgroup_count {
        return Err(SdfReadError::Parse(format!(
            "Found {} SGroups when {sgroup_count} were expected.",
            sgroups.len()
        )));
    }

    // RDKit✔️✔️:   // SGroups successfully parsed, now add them to the molecule
    // RDKit✔️✔️:   for (const auto &sg : sGroupMap) {
    // RDKit✔️✔️:     if (sg.second.getIsValid()) {
    // RDKit✔️✔️:       addSubstanceGroup(*mol, sg.second);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog) << "SGroup " << sg.first
    // RDKit✔️✔️:                               << " is invalid and will be ignored" << std::endl;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    //
    // RDKit iterates the ordered sGroupMap, so COSMolKit's index-backed
    // SubstanceGroupId values must be assigned from that final insertion
    // order, not from the order the V3000 lines appeared in.
    let id_by_sequence_id = sgroups
        .keys()
        .enumerate()
        .map(|(index, sequence_id)| (*sequence_id, SubstanceGroupId::new(index)))
        .collect::<BTreeMap<_, _>>();
    for (sequence_id, id) in &id_by_sequence_id {
        let sgroup = sgroups
            .get_mut(sequence_id)
            .ok_or_else(|| SdfReadError::Parse(format!("SGroup {sequence_id} missing")))?;
        sgroup.set_id(*id);
    }
    for (sequence_id, parent_sequence_id) in parent_by_sequence_id {
        let parent = *id_by_sequence_id.get(&parent_sequence_id).ok_or_else(|| {
            SdfReadError::Parse(format!("Invalid PARENT label found on line {line_number}"))
        })?;
        let sgroup = sgroups
            .get_mut(&sequence_id)
            .ok_or_else(|| SdfReadError::Parse(format!("SGroup {sequence_id} missing")))?;
        sgroup.set_parent(parent);
    }

    Ok(sgroups.into_values().collect())
}

fn parse_v3000_array_values<T>(
    value: &str,
    line_number: usize,
    max_count: Option<usize>,
    parse_value: impl Fn(&str) -> Result<T, SdfReadError>,
) -> Result<Vec<T>, SdfReadError> {
    let trimmed = value.trim();
    if !trimmed.starts_with('(') || !trimmed.ends_with(')') {
        return Err(SdfReadError::Parse(format!(
            "WARNING: first character of V3000 array is not '(' on line {line_number}"
        )));
    }
    let fields = trimmed[1..trimmed.len() - 1]
        .split_whitespace()
        .collect::<Vec<_>>();
    if fields.is_empty() {
        return Ok(Vec::new());
    }
    let count = parse_rdkit_unsigned(fields[0], false).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to unsigned int on line {line_number}",
            fields[0]
        ))
    })? as usize;
    if let Some(max_count) = max_count
        && count > max_count
    {
        return Err(SdfReadError::Parse("invalid count value".to_string()));
    }
    let mut values = Vec::with_capacity(count);
    for field in fields.iter().skip(1).take(count) {
        values.push(parse_value(field)?);
    }
    Ok(values)
}

fn merge_v3000_assign_tokens(tokens: &[String]) -> Vec<String> {
    let mut merged = Vec::new();
    let mut current = String::new();
    let mut in_quoted_value = false;
    for token in tokens {
        if current.is_empty() {
            current.push_str(token);
        } else {
            current.push(' ');
            current.push_str(token);
        }
        let value = current
            .split_once('=')
            .map(|(_, value)| value.trim())
            .unwrap_or_default();
        if value.starts_with('"') {
            let quote_count = value
                .as_bytes()
                .iter()
                .filter(|byte| **byte == b'"')
                .count();
            in_quoted_value = quote_count % 2 == 1;
        }
        if !in_quoted_value {
            merged.push(std::mem::take(&mut current));
        }
    }
    if !current.is_empty() {
        merged.push(current);
    }
    merged
}

fn split_v3000_sgroup_line(
    line: &str,
    line_number: usize,
) -> Result<(&str, &str, &str, Vec<String>), SdfReadError> {
    let mut fields = line.splitn(4, char::is_whitespace);
    let sequence = fields.next().unwrap_or_default();
    let typ = fields.next().unwrap_or_default();
    let external_id = fields.next().unwrap_or_default();
    let rest = fields.next().unwrap_or_default();
    if sequence.is_empty() || typ.is_empty() || external_id.is_empty() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    Ok((
        sequence,
        typ,
        external_id,
        tokenize_v3000_sgroup_labels(rest),
    ))
}

fn tokenize_v3000_sgroup_labels(text: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut start = None;
    let mut paren_depth = 0_usize;
    let mut in_quotes = false;
    let bytes = text.as_bytes();
    let mut pos = 0_usize;
    while pos < bytes.len() {
        match bytes[pos] {
            b' ' | b'\t' if !in_quotes && paren_depth == 0 => {
                if let Some(start_pos) = start.take()
                    && start_pos != pos
                {
                    tokens.push(text[start_pos..pos].to_string());
                }
                pos += 1;
            }
            b'(' if !in_quotes => {
                if start.is_none() {
                    start = Some(pos);
                }
                paren_depth += 1;
                pos += 1;
            }
            b')' if !in_quotes && paren_depth > 0 => {
                paren_depth -= 1;
                pos += 1;
            }
            b'"' => {
                if start.is_none() {
                    start = Some(pos);
                }
                if pos + 1 < bytes.len() && bytes[pos + 1] == b'"' {
                    pos += 2;
                } else {
                    in_quotes = !in_quotes;
                    pos += 1;
                }
            }
            _ => {
                if start.is_none() {
                    start = Some(pos);
                }
                pos += 1;
            }
        }
    }
    if let Some(start_pos) = start
        && start_pos != text.len()
    {
        tokens.push(text[start_pos..].to_string());
    }
    tokens
}

fn parse_v3000_u32_array(
    value: &str,
    line_number: usize,
    max_count: Option<usize>,
) -> Result<Vec<u32>, SdfReadError> {
    parse_v3000_array_values(value, line_number, max_count, |field| {
        parse_rdkit_unsigned(field, false).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{field}' to unsigned int on line {line_number}"
            ))
        })
    })
}

fn parse_v3000_f64_array(
    value: &str,
    line_number: usize,
    max_count: Option<usize>,
) -> Result<Vec<f64>, SdfReadError> {
    parse_v3000_array_values(value, line_number, max_count, |field| {
        parse_rdkit_double(field, false).map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{field}' to double on line {line_number}"
            ))
        })
    })
}

fn parse_v3000_sgroup_label(
    label: &str,
    value: &str,
    line_number: usize,
    params: SdfReadParams,
    sgroup: &mut SubstanceGroup,
    parent_by_sequence_id: &mut BTreeMap<u32, u32>,
    atom_id_by_mol_idx: &BTreeMap<u32, AtomId>,
    bond_id_by_mol_idx: &BTreeMap<u32, BondId>,
) -> Result<(), SdfReadError> {
    match label {
        "ATOMS" => {
            for atom_idx in
                parse_v3000_u32_array(value, line_number, Some(atom_id_by_mol_idx.len()))?
            {
                let atom = *atom_id_by_mol_idx.get(&atom_idx).ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "SGroup atom index {atom_idx} out of range on line {line_number}"
                    ))
                })?;
                sgroup.push_atom(atom);
            }
        }
        "PATOMS" => {
            for atom_idx in
                parse_v3000_u32_array(value, line_number, Some(atom_id_by_mol_idx.len()))?
            {
                let atom = *atom_id_by_mol_idx.get(&atom_idx).ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "SGroup parent atom index {atom_idx} out of range on line {line_number}"
                    ))
                })?;
                sgroup.push_parent_atom(atom);
            }
        }
        "CBONDS" | "XBONDS" => {
            let role = if label == "CBONDS" {
                crate::SGroupBondRole::Contained
            } else {
                crate::SGroupBondRole::Crossing
            };
            for bond_idx in
                parse_v3000_u32_array(value, line_number, Some(bond_id_by_mol_idx.len()))?
            {
                let bond = *bond_id_by_mol_idx.get(&bond_idx).ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "SGroup bond index {bond_idx} out of range on line {line_number}"
                    ))
                })?;
                sgroup.push_bond_with_role(bond, role);
            }
        }
        "BRKXYZ" => {
            let coords = parse_v3000_f64_array(value, line_number, Some(9))?;
            if coords.len() != 9 {
                return Err(SdfReadError::Parse(format!(
                    "Unexpected number of coordinates for BRKXYZ on line {line_number}"
                )));
            }
            sgroup.display_mut().brackets.push(SGroupBracket {
                p1: [coords[0], coords[1]],
                p2: [coords[3], coords[4]],
            });
        }
        "CSTATE" => parse_v3000_cstate_label(value, line_number, sgroup, bond_id_by_mol_idx)?,
        "SAP" => parse_v3000_sap_label(value, line_number, sgroup, atom_id_by_mol_idx)?,
        "PARENT" => {
            let parent_idx = parse_rdkit_unsigned(value, false).map_err(|_| {
                SdfReadError::Parse(format!("Invalid PARENT label found on line {line_number}"))
            })?;
            if let Some(sequence_id) = sgroup.rdkit_sequence_id() {
                parent_by_sequence_id.insert(sequence_id, parent_idx);
            }
        }
        "COMPNO" => {
            let compno = parse_rdkit_unsigned(value, false).map_err(|_| {
                SdfReadError::Parse(format!(
                    "Cannot convert '{value}' to unsigned int on line {line_number}"
                ))
            })?;
            if compno > 256 {
                return Err(SdfReadError::Parse(format!(
                    "SGroup SNC value over 256: '{compno}' on line {line_number}"
                )));
            }
            sgroup.set_component_number(compno);
        }
        "FIELDDATA" => {
            let mut str_value = parse_v3000_string_prop_value(value);
            if params.strict_parsing {
                str_value = rdkit_substr(&str_value, 0, 200).to_string();
            }
            sgroup.data_mut().values.push(str_value);
        }
        "SUBTYPE" => {
            let str_value = parse_v3000_string_prop_value(value);
            if !matches!(str_value.as_str(), "ALT" | "RAN" | "BLO") {
                return Err(SdfReadError::Parse(format!(
                    "Unsupported SGroup subtype '{str_value}' on line {line_number}"
                )));
            }
            sgroup.set_subtype(str_value);
        }
        "CONNECT" => {
            let str_value = parse_v3000_string_prop_value(value);
            let connection = sgroup_connection_from_rdkit(&str_value).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "Unsupported SGroup connection type '{str_value}' on line {line_number}"
                ))
            })?;
            sgroup.set_connection(connection);
        }
        "CLASS" => sgroup.set_class(parse_v3000_string_prop_value(value)),
        "LABEL" => sgroup.set_label(parse_v3000_string_prop_value(value)),
        "FIELDNAME" => sgroup.data_mut().field_name = Some(parse_v3000_string_prop_value(value)),
        "FIELDTYPE" => sgroup.data_mut().field_type = Some(parse_v3000_string_prop_value(value)),
        "FIELDINFO" => sgroup.data_mut().field_info = Some(parse_v3000_string_prop_value(value)),
        "FIELDDISP" => {
            sgroup.data_mut().field_display = Some(parse_v3000_string_prop_value(value));
        }
        "QUERYTYPE" => sgroup.data_mut().query_type = Some(parse_v3000_string_prop_value(value)),
        "QUERYOP" => sgroup.data_mut().query_op = Some(parse_v3000_string_prop_value(value)),
        "ESTATE" => sgroup.set_expansion_state(parse_v3000_string_prop_value(value)),
        "BRKTYP" => {
            let style = match parse_v3000_string_prop_value(value).as_str() {
                "BRACKET" => SGroupBracketStyle::Bracket,
                "PAREN" => SGroupBracketStyle::Parenthesis,
                "" => SGroupBracketStyle::None,
                other => SGroupBracketStyle::Unknown(other.to_string()),
            };
            sgroup.set_bracket_style(style);
        }
        other => sgroup.set_prop(other, parse_v3000_string_prop_value(value)),
    }
    Ok(())
}

fn parse_v3000_cstate_label(
    value: &str,
    line_number: usize,
    sgroup: &mut SubstanceGroup,
    bond_id_by_mol_idx: &BTreeMap<u32, BondId>,
) -> Result<(), SdfReadError> {
    let fields = value
        .trim()
        .trim_start_matches('(')
        .trim_end_matches(')')
        .split_whitespace()
        .collect::<Vec<_>>();
    if fields.len() < 2 {
        return Err(SdfReadError::Parse(format!(
            "Unexpected number of fields for CSTATE field on line {line_number}"
        )));
    }
    let count = parse_rdkit_unsigned(fields[0], false).unwrap_or(0);
    let expected = if sgroup.kind() == &SubstanceGroupKind::Superatom {
        4
    } else {
        1
    };
    if count != expected {
        return Err(SdfReadError::Parse(format!(
            "Unexpected number of fields for CSTATE field on line {line_number}"
        )));
    }
    let bond_mark = parse_rdkit_unsigned(fields[1], false).unwrap_or(0);
    let bond = *bond_id_by_mol_idx.get(&bond_mark).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup bond index {bond_mark} out of range on line {line_number}"
        ))
    })?;
    let vector = if sgroup.kind() == &SubstanceGroupKind::Superatom {
        [
            fields
                .get(2)
                .and_then(|field| parse_rdkit_double(field, false).ok())
                .unwrap_or(0.0),
            fields
                .get(3)
                .and_then(|field| parse_rdkit_double(field, false).ok())
                .unwrap_or(0.0),
        ]
    } else {
        [0.0, 0.0]
    };
    sgroup.push_cstate(SGroupCState { bond, vector });
    Ok(())
}

fn parse_v3000_sap_label(
    value: &str,
    line_number: usize,
    sgroup: &mut SubstanceGroup,
    atom_id_by_mol_idx: &BTreeMap<u32, AtomId>,
) -> Result<(), SdfReadError> {
    let fields = value
        .trim()
        .trim_start_matches('(')
        .trim_end_matches(')')
        .split_whitespace()
        .collect::<Vec<_>>();
    if fields.len() < 4 {
        return Err(SdfReadError::Parse(format!(
            "SGroup SAP line too short on line {line_number}"
        )));
    }
    let atom_mark = parse_rdkit_unsigned(fields[1], false).unwrap_or(0);
    let atom = *atom_id_by_mol_idx.get(&atom_mark).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup attach atom index {atom_mark} out of range on line {line_number}"
        ))
    })?;
    let leaving_atom = if fields[2].eq_ignore_ascii_case("AIDX") {
        Some(atom)
    } else {
        let leaving_mark = parse_rdkit_unsigned(fields[2], false).unwrap_or(0);
        if leaving_mark == 0 {
            None
        } else {
            Some(*atom_id_by_mol_idx.get(&leaving_mark).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup leaving atom index {leaving_mark} out of range on line {line_number}"
                ))
            })?)
        }
    };
    sgroup.push_attach_point(SGroupAttachPoint {
        atom,
        leaving_atom,
        label: Some(fields[3].to_string()),
        order: None,
    });
    Ok(())
}

fn parse_v3000_string_prop_value(value: &str) -> String {
    let trimmed = value.trim_end();
    if trimmed.len() >= 2 && trimmed.starts_with('"') && trimmed.ends_with('"') {
        trimmed[1..trimmed.len() - 1].replace("\"\"", "\"")
    } else if trimmed.len() >= 2 && trimmed.starts_with('\'') && trimmed.ends_with('\'') {
        trimmed[1..trimmed.len() - 1].to_string()
    } else {
        trimmed.trim().to_string()
    }
}

fn parse_v3000_collection_block<R: BufRead>(
    reader: &mut R,
    line_number: &mut usize,
    params: SdfReadParams,
    atom_id_by_mol_idx: &BTreeMap<u32, AtomId>,
) -> Result<Vec<StereoGroup>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_collection_block
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::string parseEnhancedStereo
    // RDKit✔️❗:
    // RDKit✔️❗:   // Lines like (absolute, relative, racemic):
    // RDKit✔️❗:   // M  V30 MDLV30/STEABS ATOMS=(2 2 3)
    // RDKit✔️❗:   // M  V30 MDLV30/STEREL1 ATOMS=(1 12)
    // RDKit✔️❗:   // M  V30 MDLV30/STERAC1 ATOMS=(1 12)
    // RDKit✔️❗:   const regex stereo_label(
    // RDKit✔️❗:       R"regex(MDLV30/STE(...)([0-9]*) +ATOMS=\(([0-9]+) +(.*)\) *)regex");
    // RDKit✔️❗:
    // RDKit✔️❗:   smatch match;
    // RDKit✔️❗:   std::vector<StereoGroup> groups;
    // RDKit✔️❗:
    // RDKit✔️❗:   // Read the collection until the end
    // RDKit✔️❗:   auto tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   unsigned abs_group_seen = 0;
    // RDKit✔️❗:   while (!startsWith(tempStr, "END", 3)) {
    // RDKit✔️❗:     // If this line in the collection is part of a stereo group
    // RDKit✔️❗:     if (regex_match(tempStr, match, stereo_label)) {
    // RDKit✔️❗:       StereoGroupType grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
    // RDKit✔️❗:       unsigned groupid = 0;
    // RDKit✔️❗:
    // RDKit✔️❗:       if (match[1] == "ABS") {
    // RDKit✔️❗:         grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
    // RDKit✔️❗:         // Warn only one per mol about multiple ABS groups
    // RDKit✔️❗:         if (abs_group_seen == 1) {
    // RDKit✔️❗:           std::ostringstream errout;
    // RDKit✔️❗:           errout << "Seen a second ABS stereo group on line " << line
    // RDKit✔️❗:                  << std::endl;
    // RDKit✔️❗:           if (strictParsing) {
    // RDKit✔️❗:             throw FileParseException(errout.str());
    // RDKit✔️❗:           } else {
    // RDKit✔️❗:             BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:           }
    // RDKit✔️❗:         }
    // RDKit✔️❗:         ++abs_group_seen;
    // RDKit✔️❗:       } else if (match[1] == "REL") {
    // RDKit✔️❗:         grouptype = RDKit::StereoGroupType::STEREO_OR;
    // RDKit✔️❗:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit✔️❗:       } else if (match[1] == "RAC") {
    // RDKit✔️❗:         grouptype = RDKit::StereoGroupType::STEREO_AND;
    // RDKit✔️❗:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "Unrecognized stereogroup type : '" << tempStr << "' on line"
    // RDKit✔️❗:                << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:       const unsigned int count = FileParserUtils::toUnsigned(match[3], true);
    // RDKit✔️❗:       std::vector<Atom *> atoms;
    // RDKit✔️❗:       std::stringstream ss(match[4]);
    // RDKit✔️❗:       unsigned int index;
    // RDKit✔️❗:       for (size_t i = 0; i < count; ++i) {
    // RDKit✔️❗:         ss >> index;
    // RDKit✔️❗:         // atoms are 1 indexed in molfiles
    // RDKit✔️❗:         atoms.push_back(mol->getAtomWithIdx(index - 1));
    // RDKit✔️❗:       }
    // RDKit✔️❗:       std::vector<Bond *> newBonds;
    // RDKit✔️❗:       groups.emplace_back(grouptype, std::move(atoms), std::move(newBonds),
    // RDKit✔️❗:                           groupid);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       // skip collection types we don't know how to read. Only one documented
    // RDKit✔️❗:       // is MDLV30/HILITE
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << "Skipping unrecognized collection type at "
    // RDKit✔️❗:                                  "line "
    // RDKit✔️❗:                               << line << ": " << tempStr << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:     tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (!groups.empty()) {
    // RDKit✔️❗:     mol->setStereoGroups(std::move(groups));
    // RDKit✔️❗:   }
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   return tempStr;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::string parseEnhancedStereo
    // END RDKIT CPP BODY: parse_v3000_collection_block

    let mut groups = Vec::new();
    let mut abs_group_seen = 0_u32;
    let mut temp_str = get_v3000_line(reader, line_number, params)?;
    let mut temp_upper = temp_str.to_ascii_uppercase();
    while !temp_upper.starts_with("END") {
        if let Some(group) =
            parse_v3000_stereo_collection_line(&temp_upper, *line_number, atom_id_by_mol_idx)?
        {
            if group.kind() == StereoGroupKind::Absolute {
                if abs_group_seen == 1 && params.strict_parsing {
                    return Err(SdfReadError::Parse(format!(
                        "Seen a second ABS stereo group on line {}\n",
                        *line_number
                    )));
                }
                abs_group_seen += 1;
            }
            groups.push(group);
        }
        temp_str = get_v3000_line(reader, line_number, params)?;
        temp_upper = temp_str.to_ascii_uppercase();
    }
    Ok(groups)
}

fn parse_v3000_stereo_collection_line(
    line: &str,
    line_number: usize,
    atom_id_by_mol_idx: &BTreeMap<u32, AtomId>,
) -> Result<Option<StereoGroup>, SdfReadError> {
    let Some(rest) = line.strip_prefix("MDLV30/STE") else {
        return Ok(None);
    };
    if rest.len() < 3 {
        return Err(SdfReadError::Parse(format!(
            "Unrecognized stereogroup type : '{line}' on line{line_number}"
        )));
    }
    let tag = &rest[..3];
    let after_tag = &rest[3..];
    let group_digits = after_tag
        .chars()
        .take_while(char::is_ascii_digit)
        .collect::<String>();
    let after_digits = &after_tag[group_digits.len()..];
    let Some(atoms_pos) = after_digits.find("ATOMS=") else {
        return Ok(None);
    };
    let atom_value = &after_digits[atoms_pos + 6..];
    let mol_atom_indices =
        parse_v3000_u32_array(atom_value, line_number, Some(atom_id_by_mol_idx.len()))?;
    let mut atoms = Vec::with_capacity(mol_atom_indices.len());
    for atom_idx in mol_atom_indices {
        atoms.push(*atom_id_by_mol_idx.get(&atom_idx).ok_or_else(|| {
            SdfReadError::Parse(format!(
                "Stereo group atom index {atom_idx} out of range on line {line_number}"
            ))
        })?);
    }
    let kind = match tag {
        "ABS" => StereoGroupKind::Absolute,
        "REL" => StereoGroupKind::Or,
        "RAC" => StereoGroupKind::And,
        _ => {
            return Err(SdfReadError::Parse(format!(
                "Unrecognized stereogroup type : '{line}' on line{line_number}"
            )));
        }
    };
    let mut group = StereoGroup::new(kind, atoms, Vec::new());
    if !group_digits.is_empty() {
        group = group.with_id(parse_rdkit_unsigned(&group_digits, true).unwrap_or(0));
    }
    Ok(Some(group))
}

fn parse_v3000_obj3d_block<R: BufRead>(
    reader: &mut R,
    constraint_count: usize,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: parse_v3000_obj3d_block
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV3000CTAB
    // RDKit✔️❗:
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❗:
    // RDKit✔️❗:   std::string tempStr;
    // RDKit✔️❗:   std::vector<std::string> splitLine;
    // RDKit✔️❗:
    // RDKit✔️❗:   bool fileComplete = false;
    // RDKit✔️❗:
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN CTAB") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN CTAB line not found on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   if (tempStr.size() < 8 || tempStr.substr(0, 7) != "COUNTS ") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Bad counts line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   std::string trimmed =
    // RDKit✔️❗:       boost::trim_copy(tempStr.substr(7, tempStr.length() - 7));
    // RDKit✔️❗:   boost::split(splitLine, trimmed, boost::is_any_of(" \t"),
    // RDKit✔️❗:                boost::token_compress_on);
    // RDKit✔️❗:   if (splitLine.size() < 2) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Bad counts line : '" << tempStr << "' on line " << line;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   nAtoms = FileParserUtils::toUnsigned(splitLine[0]);
    // RDKit✔️❗:   nBonds = FileParserUtils::toUnsigned(splitLine[1]);
    // RDKit✔️❗:   conf = new Conformer(nAtoms);
    // RDKit✔️❗:
    // RDKit✔️❗:   unsigned int nSgroups = 0, n3DConstraints = 0, chiralFlag = 0;
    // RDKit✔️❗:
    // RDKit✔️❗:   if (splitLine.size() > 2) {
    // RDKit✔️❗:     nSgroups = FileParserUtils::toUnsigned(splitLine[2]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (splitLine.size() > 3) {
    // RDKit✔️❗:     n3DConstraints = FileParserUtils::toUnsigned(splitLine[3]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (splitLine.size() > 4) {
    // RDKit✔️❗:     chiralFlag = FileParserUtils::toUnsigned(splitLine[4]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   mol->setProp(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit✔️❗:
    // RDKit✔️❗:   if (nAtoms) {
    // RDKit✔️❗:     ParseV3000AtomBlock(inStream, line, nAtoms, mol, conf, strictParsing,
    // RDKit✔️❗:                         expectMacroAtoms);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (nBonds) {
    // RDKit✔️❗:     ParseV3000BondBlock(inStream, line, nBonds, mol, chiralityPossible);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   // do link nodes:
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   while (tempStr.length() > 8 && tempStr.substr(0, 8) == "LINKNODE") {
    // RDKit✔️❗:     boost::to_upper(tempStr);
    // RDKit✔️❗:     // if the line has nothing on it we just ignore it
    // RDKit✔️❗:     if (tempStr.size() > 9) {
    // RDKit✔️❗:       std::string existing = "";
    // RDKit✔️❗:       if (mol->getPropIfPresent(common_properties::molFileLinkNodes,
    // RDKit✔️❗:                                 existing)) {
    // RDKit✔️❗:         existing += "|";
    // RDKit✔️❗:       }
    // RDKit✔️❗:       existing += tempStr.substr(9);  // skip the "LINKNODE "
    // RDKit✔️❗:       mol->setProp(common_properties::molFileLinkNodes, existing);
    // RDKit✔️❗:     }
    // RDKit✔️❗:     tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   bool sgroupFound = false;
    // RDKit✔️❗:   bool obj3dFound = false;
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   while (tempStr.length() > 5 && tempStr.substr(0, 5) == "BEGIN") {
    // RDKit✔️❗:     if (tempStr.length() >= 12 && tempStr.substr(0, 12) == "BEGIN SGROUP") {
    // RDKit✔️❗:       if (sgroupFound) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN SGROUP found more than once on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:
    // RDKit✔️❗:       } else if (!nSgroups) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN SGROUP  found but Sgroups NOT expected on line "
    // RDKit✔️❗:                << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:           // Prepare to read a lot of sgroups
    // RDKit✔️❗:           nSgroups = std::numeric_limits<unsigned int>::max();
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       sgroupFound = true;
    // RDKit✔️❗:       tempStr =
    // RDKit✔️❗:           ParseV3000SGroupsBlock(inStream, line, nSgroups, mol, strictParsing);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:       if (tempStr.length() < 10 || tempStr.substr(0, 10) != "END SGROUP") {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "END SGROUP line not found on line " << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:         boost::to_upper(tempStr);
    // RDKit✔️❗:       }
    // RDKit✔️❗:
    // RDKit✔️❗:     } else if (tempStr.length() >= 15 &&
    // RDKit✔️❗:                tempStr.substr(6, 10) == "COLLECTION") {
    // RDKit✔️❗:       tempStr = parseEnhancedStereo(inStream, line, mol, strictParsing);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:     } else if (tempStr.length() >= 11 &&
    // RDKit✔️❗:                tempStr.substr(0, 11) == "BEGIN OBJ3D") {
    // RDKit✔️❗:       if (obj3dFound) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN OBJ3D found more than once on line " << line;
    // RDKit✔️❗:         throw FileParseException(errout.str());
    // RDKit✔️❗:       }
    // RDKit✔️❗:       if (!n3DConstraints) {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "BEGIN OBJ3D found but 3n3DConstraints NOT expected on line "
    // RDKit✔️❗:                << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❗:           << "3D constraint information in mol block ignored at line " << line
    // RDKit✔️❗:           << std::endl;
    // RDKit✔️❗:       obj3dFound = true;
    // RDKit✔️❗:       for (unsigned int i = 0; i < n3DConstraints; ++i) {
    // RDKit✔️❗:         tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       }
    // RDKit✔️❗:       tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:       if (tempStr.length() < 9 || tempStr.substr(0, 9) != "END OBJ3D") {
    // RDKit✔️❗:         std::ostringstream errout;
    // RDKit✔️❗:         errout << "END OBJ3D line not found on line " << line;
    // RDKit✔️❗:         if (strictParsing) {
    // RDKit✔️❗:           throw FileParseException(errout.str());
    // RDKit✔️❗:         } else {
    // RDKit✔️❗:           BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:         }
    // RDKit✔️❗:       }
    // RDKit✔️❗:       tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       // skip blocks we don't know how to read
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << "skipping block at line " << line << ": '"
    // RDKit✔️❗:                               << tempStr << "'" << std::endl;
    // RDKit✔️❗:       while (tempStr.length() < 3 || tempStr.substr(0, 3) != "END") {
    // RDKit✔️❗:         tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       }
    // RDKit✔️❗:       tempStr = getV3000Line(inStream, line);
    // RDKit✔️❗:       boost::to_upper(tempStr);
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (nSgroups && !sgroupFound) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN SGROUP line not found on line " << line;
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (n3DConstraints && !obj3dFound) {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "BEGIN OBJ3D line not found on line " << line;
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   boost::to_upper(tempStr);
    // RDKit✔️❗:   if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END CTAB") {
    // RDKit✔️❗:     if (strictParsing) {
    // RDKit✔️❗:       throw FileParseException("END CTAB line not found");
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       BOOST_LOG(rdWarningLog) << "END CTAB line not found." << std::endl;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   if (expectMEND) {
    // RDKit✔️❗:     tempStr = getLine(inStream);
    // RDKit✔️❗:     ++line;
    // RDKit✔️❗:     if (tempStr[0] == 'M' && tempStr.substr(0, 6) == "M  END") {
    // RDKit✔️❗:       fileComplete = true;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   } else {
    // RDKit✔️❗:     fileComplete = true;
    // RDKit✔️❗:   }
    // RDKit✔️❗:
    // RDKit✔️❗:   auto is3d = calculate3dFlag(*mol, *conf, chiralityPossible);
    // RDKit✔️❗:   conf->set3D(is3d);
    // RDKit✔️❗:   mol->addConformer(conf, true);
    // RDKit✔️❗:   conf = nullptr;
    // RDKit✔️❗:
    // RDKit✔️❗:   return fileComplete;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: bool ParseV3000CTAB
    // END RDKIT CPP BODY: parse_v3000_obj3d_block

    let _ = params;
    for _ in 0..constraint_count {
        let _ = get_v3000_line(reader, line_number, params)?;
    }
    let end = get_v3000_line(reader, line_number, params)?;
    if (end.to_ascii_uppercase().len() < 9 || !end.to_ascii_uppercase().starts_with("END OBJ3D"))
        && params.strict_parsing
    {
        return Err(SdfReadError::Parse(format!(
            "END OBJ3D line not found on line {}",
            *line_number
        )));
    }
    Ok(())
}

fn tokenize_v3000_line(line: &str, line_number: usize) -> Result<Vec<String>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: tokenize_v3000_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void tokenizeV3000Line
    // RDKit✔️❗:
    // RDKit✔️❗:   tokens.clear();
    // RDKit✔️❗:   bool inQuotes = false;
    // RDKit✔️❗:   unsigned int parenDepth = 0;
    // RDKit✔️❗:   unsigned int start = 0;
    // RDKit✔️❗:   unsigned int pos = 0;
    // RDKit✔️❗:   while (pos < line.size()) {
    // RDKit✔️❗:     if (line[pos] == ' ' || line[pos] == '\t') {
    // RDKit✔️❗:       if (start == pos) {
    // RDKit✔️❗:         ++start;
    // RDKit✔️❗:         ++pos;
    // RDKit✔️❗:       } else if (!inQuotes && parenDepth == 0) {
    // RDKit✔️❗:         tokens.push_back(line.substr(start, pos - start));
    // RDKit✔️❗:         ++pos;
    // RDKit✔️❗:         start = pos;
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         ++pos;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else if (line[pos] == ')' && parenDepth > 0) {
    // RDKit✔️❗:       --parenDepth;
    // RDKit✔️❗:       ++pos;
    // RDKit✔️❗:     } else if (line[pos] == '(' && !inQuotes) {
    // RDKit✔️❗:       ++parenDepth;
    // RDKit✔️❗:       ++pos;
    // RDKit✔️❗:     } else if (line[pos] == '"' && parenDepth == 0) {
    // RDKit✔️❗:       if (pos + 1 < line.size() && line[pos + 1] == '"') {
    // RDKit✔️❗:         pos += 2;
    // RDKit✔️❗:       } else if (inQuotes) {
    // RDKit✔️❗:         // don't push on the quotes themselves
    // RDKit✔️❗:         tokens.push_back(line.substr(start + 1, pos - start - 1));
    // RDKit✔️❗:         ++pos;
    // RDKit✔️❗:         start = pos;
    // RDKit✔️❗:         inQuotes = false;
    // RDKit✔️❗:       } else {
    // RDKit✔️❗:         ++pos;
    // RDKit✔️❗:         inQuotes = true;
    // RDKit✔️❗:       }
    // RDKit✔️❗:     } else {
    // RDKit✔️❗:       ++pos;
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   if (start != pos) {
    // RDKit✔️❗:     tokens.push_back(line.substr(start, line.size() - start));
    // RDKit✔️❗:   }
    // RDKit✔️❗: #if 0
    // RDKit✔️❗:       std::cerr<<"tokens: ";
    // RDKit✔️❗:       std::copy(tokens.begin(),tokens.end(),std::ostream_iterator<std::string>(std::cerr,"|"));
    // RDKit✔️❗:       std::cerr<<std::endl;
    // RDKit✔️❗: #endif
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void tokenizeV3000Line
    // END RDKIT CPP BODY: tokenize_v3000_line

    let _ = line_number;
    let mut tokens = Vec::new();
    let mut in_quotes = false;
    let mut paren_depth = 0_usize;
    let mut start = 0_usize;
    let mut pos = 0_usize;
    let bytes = line.as_bytes();

    while pos < bytes.len() {
        match bytes[pos] {
            b' ' | b'\t' => {
                if start == pos {
                    start += 1;
                    pos += 1;
                } else if !in_quotes && paren_depth == 0 {
                    tokens.push(line[start..pos].to_string());
                    pos += 1;
                    start = pos;
                } else {
                    pos += 1;
                }
            }
            b')' if paren_depth > 0 => {
                paren_depth -= 1;
                pos += 1;
            }
            b'(' if !in_quotes => {
                paren_depth += 1;
                pos += 1;
            }
            b'"' if paren_depth == 0 => {
                if pos + 1 < bytes.len() && bytes[pos + 1] == b'"' {
                    pos += 2;
                } else if in_quotes {
                    tokens.push(line[start + 1..pos].to_string());
                    pos += 1;
                    start = pos;
                    in_quotes = false;
                } else {
                    pos += 1;
                    in_quotes = true;
                }
            }
            _ => {
                pos += 1;
            }
        }
    }

    if start != pos {
        tokens.push(line[start..line.len()].to_string());
    }

    Ok(tokens)
}

fn get_v3000_line<R: BufRead>(
    reader: &mut R,
    line_number: &mut usize,
    params: SdfReadParams,
) -> Result<String, SdfReadError> {
    // BEGIN RDKIT CPP BODY: get_v3000_line
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::string getV3000Line
    // RDKit✔️❗:
    // RDKit✔️❗:   // FIX: technically V3K blocks are case-insensitive. We should really be
    // RDKit✔️❗:   // up-casing everything here.
    // RDKit✔️❗:   PRECONDITION(inStream, "bad stream");
    // RDKit✔️❗:   std::string res;
    // RDKit✔️❗:   ++line;
    // RDKit✔️❗:   auto inl = getLine(inStream);
    // RDKit✔️❗:   std::string_view tempStr = inl;
    // RDKit✔️❗:   if (tempStr.size() < 7 || tempStr.substr(0, 7) != "M  V30 ") {
    // RDKit✔️❗:     std::ostringstream errout;
    // RDKit✔️❗:     errout << "Line " << line << " does not start with 'M  V30 '" << std::endl;
    // RDKit✔️❗:     throw FileParseException(errout.str());
    // RDKit✔️❗:   }
    // RDKit✔️❗:   // FIX: do we need to handle trailing whitespace after a -?
    // RDKit✔️❗:   while (tempStr.back() == '-') {
    // RDKit✔️❗:     // continuation character, append what we read:
    // RDKit✔️❗:     res += tempStr.substr(7, tempStr.length() - 8);
    // RDKit✔️❗:     // and then read another line:
    // RDKit✔️❗:     ++line;
    // RDKit✔️❗:     inl = getLine(inStream);
    // RDKit✔️❗:     tempStr = inl;
    // RDKit✔️❗:     if (tempStr.size() < 7 || tempStr.substr(0, 7) != "M  V30 ") {
    // RDKit✔️❗:       std::ostringstream errout;
    // RDKit✔️❗:       errout << "Line " << line << " does not start with 'M  V30 '"
    // RDKit✔️❗:              << std::endl;
    // RDKit✔️❗:       throw FileParseException(errout.str());
    // RDKit✔️❗:     }
    // RDKit✔️❗:   }
    // RDKit✔️❗:   res += tempStr.substr(7, tempStr.length() - 7);
    // RDKit✔️❗:
    // RDKit✔️❗:   return res;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: std::string getV3000Line
    // END RDKIT CPP BODY: get_v3000_line

    let _ = params;
    let mut res = String::new();
    *line_number += 1;
    let mut temp_str = read_rdkit_line(reader)?.unwrap_or_default();
    if temp_str.len() < 7 || !temp_str.starts_with("M  V30 ") {
        return Err(SdfReadError::Parse(format!(
            "Line {} does not start with 'M  V30 '\n",
            line_number
        )));
    }

    // RDKit stores `tempStr` as a string_view over an owning std::string. Rust
    // keeps an owned String so continuation reads can replace it directly.
    while temp_str.ends_with('-') {
        res.push_str(&temp_str[7..temp_str.len() - 1]);
        *line_number += 1;
        temp_str = read_rdkit_line(reader)?.unwrap_or_default();
        if temp_str.len() < 7 || !temp_str.starts_with("M  V30 ") {
            return Err(SdfReadError::Parse(format!(
                "Line {} does not start with 'M  V30 '\n",
                line_number
            )));
        }
    }
    res.push_str(&temp_str[7..]);

    Ok(res)
}

mod postprocess;

use self::postprocess::*;

#[cfg(test)]
mod tests {
    use crate::{BatchErrorMode, MoleculeBatch};

    use super::*;

    fn write_temp_sdf(input: &str) -> tempfile::NamedTempFile {
        use std::io::Write;

        let mut file = tempfile::NamedTempFile::new().expect("should create temporary SDF file");
        file.write_all(input.as_bytes())
            .expect("should write temporary SDF file");
        file
    }

    fn v2000_atom_line(
        symbol: &str,
        mass: i32,
        charge: i32,
        h_count: i32,
        atom_map: i32,
    ) -> String {
        format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}",
            1.25, -2.5, 0.75, symbol, mass, charge, 0, h_count, 0, 0, 0, 0, 0, atom_map
        )
    }

    fn v2000_atom_line_at(symbol: &str, x: f64, y: f64, z: f64) -> String {
        format!(
            "{x:>10.4}{y:>10.4}{z:>10.4} {symbol:<3}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}",
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        )
    }

    fn v2000_atom_line_with_z(symbol: &str, z: f64) -> String {
        format!(
            "{:>10.4}{:>10.4}{z:>10.4} {:<3}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}",
            1.25, -2.5, symbol, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        )
    }

    fn v2000_bond_line(begin: u32, end: u32, bond_type: u32, stereo: u32, topology: u32) -> String {
        format!(
            "{begin:>3}{end:>3}{bond_type:>3}{stereo:>3}{:>3}{topology:>3}",
            0
        )
    }

    fn single_atom_sdf_record(name: &str, symbol: &str, trailing_newline: bool) -> String {
        let mut record = format!(
            "{name}\n\n\n  1  0  0  0  0  0            999 V2000\n{}\nM  END\n$$$$",
            v2000_atom_line(symbol, 0, 0, 0, 0)
        );
        if trailing_newline {
            record.push('\n');
        }
        record
    }

    fn impossible_valence_sdf_record(name: &str) -> String {
        let mut record = format!("{name}\n\n\n  5  4  0  0  0  0            999 V2000\n");
        for _ in 0..5 {
            record.push_str(&v2000_atom_line("C", 0, 0, 0, 0));
            record.push('\n');
        }
        for end in 2..=5 {
            record.push_str(&v2000_bond_line(1, end, 3, 0, 0));
            record.push('\n');
        }
        record.push_str("M  END\n$$$$\n");
        record
    }

    fn v2000_tetrahedral_wedge_block(stereo: u32, explicit_h: bool) -> String {
        let mut atom_lines = vec![
            v2000_atom_line_at("C", 0.0, 0.0, 0.0),
            v2000_atom_line_at(if explicit_h { "H" } else { "F" }, 1.0, 0.0, 0.0),
            v2000_atom_line_at("Cl", 0.0, 1.0, 0.0),
            v2000_atom_line_at("Br", -1.0, -1.0, 0.0),
        ];
        if explicit_h {
            atom_lines.push(v2000_atom_line_at("F", 0.0, -1.0, 0.0));
        }
        let atom_count = atom_lines.len();
        let bond_count = if explicit_h { 4 } else { 3 };
        let mut block = format!(
            "wedged\n  COSMolKit          2D\ncomment\n{atom_count:>3}{bond_count:>3}  0  0  0  0            999 V2000\n"
        );
        for line in atom_lines {
            block.push_str(&line);
            block.push('\n');
        }
        block.push_str(&v2000_bond_line(1, 2, 1, stereo, 0));
        block.push('\n');
        block.push_str(&v2000_bond_line(1, 3, 1, 0, 0));
        block.push('\n');
        block.push_str(&v2000_bond_line(1, 4, 1, 0, 0));
        block.push('\n');
        if explicit_h {
            block.push_str(&v2000_bond_line(1, 5, 1, 0, 0));
            block.push('\n');
        }
        block.push_str("M  END\n$$$$\n");
        block
    }

    fn v3000_outer_counts_line(atom_count: usize, bond_count: usize) -> String {
        let mut line = format!("{atom_count:>3}{bond_count:>3}  0  0  0            999");
        while line.len() < 34 {
            line.push(' ');
        }
        line.push_str("V3000");
        line
    }

    #[test]
    fn parse_sdf_data_header_extracts_label_like_rdkit() {
        assert_eq!(
            parse_sdf_data_header("> <ID>", 1).unwrap(),
            Some("ID".to_string())
        );
        assert_eq!(
            parse_sdf_data_header(" \t>  <Long Name>  \r\n", 2).unwrap(),
            Some("Long Name".to_string())
        );
        assert_eq!(
            parse_sdf_data_header("> <A> trailing text", 3).unwrap(),
            Some("A".to_string())
        );
    }

    #[test]
    fn parse_sdf_data_header_ignores_missing_or_empty_labels_like_rdkit() {
        assert_eq!(parse_sdf_data_header("", 1).unwrap(), None);
        assert_eq!(parse_sdf_data_header("not a header", 2).unwrap(), None);
        assert_eq!(parse_sdf_data_header(">", 3).unwrap(), None);
        assert_eq!(parse_sdf_data_header("> <>", 4).unwrap(), None);
        assert_eq!(parse_sdf_data_header("> <missing_end", 5).unwrap(), None);
    }

    #[test]
    fn parse_sdf_data_header_matches_rdkit_substr_underflow_case() {
        assert_eq!(
            parse_sdf_data_header(">> <FIELD", 1).unwrap(),
            Some("FIELD".to_string())
        );
    }

    #[test]
    fn parse_sdf_data_fields_reads_single_and_multiline_values_like_rdkit() {
        let lines = [
            "> <ID>", "cmpd-1", "", "> <NOTE>", "alpha", "beta", "", "$$$$",
        ];
        let fields = parse_sdf_data_fields(&lines, 10, SdfReadParams::default()).unwrap();

        assert_eq!(
            fields,
            vec![
                SdfDataField {
                    name: "ID".to_string(),
                    value: "cmpd-1".to_string(),
                },
                SdfDataField {
                    name: "NOTE".to_string(),
                    value: "alpha\nbeta".to_string(),
                },
            ]
        );
    }

    #[test]
    fn parse_sdf_data_fields_keeps_indented_blank_lines_like_rdkit() {
        let lines = ["> <SPACES>", "  ", "\t", "", "$$$$"];
        let fields = parse_sdf_data_fields(&lines, 1, SdfReadParams::default()).unwrap();

        assert_eq!(
            fields,
            vec![SdfDataField {
                name: "SPACES".to_string(),
                value: "  \n\t".to_string(),
            }]
        );
    }

    #[test]
    fn parse_sdf_data_fields_strips_terminal_cr_like_rdkit() {
        let lines = ["> <ID>\r", "value\r", "", "$$$$"];
        let fields = parse_sdf_data_fields(&lines, 1, SdfReadParams::default()).unwrap();

        assert_eq!(
            fields,
            vec![SdfDataField {
                name: "ID".to_string(),
                value: "value".to_string(),
            }]
        );
    }

    #[test]
    fn parse_sdf_data_fields_ignores_empty_or_missing_labels_until_blank_like_rdkit() {
        let lines = ["> <>", "ignored", "", "> <ID>", "kept", "", "$$$$"];
        let fields = parse_sdf_data_fields(&lines, 1, SdfReadParams::default()).unwrap();

        assert_eq!(
            fields,
            vec![SdfDataField {
                name: "ID".to_string(),
                value: "kept".to_string(),
            }]
        );
    }

    #[test]
    fn parse_sdf_data_fields_errors_when_unlabeled_field_reaches_eof_like_rdkit() {
        let lines = ["> <>", "unterminated"];
        let err = parse_sdf_data_fields(&lines, 1, SdfReadParams::default()).unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse("End of data field name not found".to_string())
        );
    }

    #[test]
    fn parse_sdf_data_fields_enforces_strict_spurious_data_like_rdkit() {
        let lines = ["spurious", "$$$$"];
        let err = parse_sdf_data_fields(&lines, 1, SdfReadParams::default()).unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse("Problems encountered parsing data fields".to_string())
        );
    }

    #[test]
    fn parse_sdf_data_fields_ignores_spurious_data_when_not_strict_like_rdkit() {
        let lines = ["spurious", "> <ID>", "kept", "", "$$$$"];
        let fields = parse_sdf_data_fields(
            &lines,
            1,
            SdfReadParams {
                strict_parsing: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(
            fields,
            vec![SdfDataField {
                name: "ID".to_string(),
                value: "kept".to_string(),
            }]
        );
    }

    #[test]
    fn extract_next_raw_sdf_record_reads_through_delimiter_like_rdkit() {
        let input = "mol\nM  END\n$$$$\nnext\nM  END\n$$$$\n";
        let mut reader = std::io::Cursor::new(input.as_bytes());

        let raw = extract_next_raw_sdf_record(&mut reader, 7, 11, 3, SdfReadParams::default())
            .unwrap()
            .unwrap();

        assert_eq!(raw.index, 7);
        assert_eq!(raw.byte_offset, 11);
        assert_eq!(raw.start_line, 3);
        assert_eq!(raw.text, "mol\nM  END\n$$$$\n");
        assert_eq!(raw.byte_len, raw.text.len() as u64);
    }

    #[test]
    fn extract_next_raw_sdf_record_ignores_trailing_newlines_like_rdkit() {
        let mut reader = std::io::Cursor::new("\n\r\n".as_bytes());

        let raw =
            extract_next_raw_sdf_record(&mut reader, 0, 0, 0, SdfReadParams::default()).unwrap();

        assert_eq!(raw, None);
    }

    #[test]
    fn read_next_sdf_record_returns_first_molecule_and_consumes_delimiter_like_rdkit() {
        let first = single_atom_sdf_record("first", "C", true);
        let second = single_atom_sdf_record("second", "N", true);
        let input = format!("{first}{second}");
        let first_len = first.len() as u64;
        let first_line_count = first.bytes().filter(|byte| *byte == b'\n').count();
        let mut reader = SdfReader::new(std::io::Cursor::new(input.as_bytes()));

        let record = read_next_sdf_record(&mut reader).unwrap().unwrap();

        assert_eq!(record.molecule.properties().name(), Some("first"));
        assert_eq!(record.molecule.num_atoms(), 1);
        assert_eq!(reader.next_index, 1);
        assert_eq!(reader.byte_offset, first_len);
        assert_eq!(reader.line_number, first_line_count);
        assert!(!reader.end);
        assert!(!reader.eof_hit_on_read);

        let remaining = reader.reader.fill_buf().unwrap();
        assert!(remaining.starts_with(b"second\n"));
    }

    #[test]
    fn read_next_sdf_record_sets_end_after_last_record_without_trailing_newline_like_rdkit() {
        let input = single_atom_sdf_record("last", "O", false);
        let expected_line_count = input.bytes().filter(|byte| *byte == b'\n').count() + 1;
        let mut reader = SdfReader::new(std::io::Cursor::new(input.as_bytes()));

        let record = read_next_sdf_record(&mut reader).unwrap().unwrap();

        assert_eq!(record.molecule.properties().name(), Some("last"));
        assert_eq!(reader.line_number, expected_line_count);
        assert!(reader.end);
        assert!(!reader.eof_hit_on_read);
        assert!(read_next_sdf_record(&mut reader).unwrap().is_none());
        assert!(reader.end);
    }

    #[test]
    fn read_next_sdf_record_marks_empty_eof_like_rdkit() {
        let mut reader = SdfReader::new(std::io::Cursor::new(Vec::<u8>::new()));

        let record = read_next_sdf_record(&mut reader).unwrap();

        assert!(record.is_none());
        assert_eq!(reader.next_index, 0);
        assert_eq!(reader.byte_offset, 0);
        assert_eq!(reader.line_number, 0);
        assert!(reader.end);
        assert!(reader.eof_hit_on_read);
    }

    #[test]
    fn sdf_supplier_recovery_parse_exception_returns_none_then_next_record_like_rdkit() {
        let malformed = format!(
            "bad\n\n\n  1  0  0  0  0  0            999 V2000\n{}\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0)
        );
        let valid = single_atom_sdf_record("after", "N", true);
        let input = format!("{malformed}{valid}");
        let mut reader = SdfReader::new(std::io::Cursor::new(input.as_bytes()));

        assert!(read_next_sdf_record(&mut reader).unwrap().is_none());
        assert!(!reader.end);

        let record = read_next_sdf_record(&mut reader).unwrap().unwrap();
        assert_eq!(record.molecule.properties().name(), Some("after"));
        assert_eq!(record.molecule.num_atoms(), 1);
    }

    #[test]
    fn sdf_supplier_recovery_sanitize_exception_returns_none_then_next_record_like_rdkit() {
        let malformed = impossible_valence_sdf_record("bad-valence");
        let valid = single_atom_sdf_record("after", "O", true);
        let input = format!("{malformed}{valid}");
        let mut reader = SdfReader::with_params(
            std::io::Cursor::new(input.as_bytes()),
            SdfReadParams {
                sanitize: true,
                remove_hs: false,
                ..Default::default()
            },
        );

        assert!(read_next_sdf_record(&mut reader).unwrap().is_none());
        assert!(!reader.end);

        let record = read_next_sdf_record(&mut reader).unwrap().unwrap();
        assert_eq!(record.molecule.properties().name(), Some("after"));
        assert_eq!(record.molecule.num_atoms(), 1);
    }

    #[test]
    fn sdf_supplier_recovery_malformed_data_field_keeps_molecule_like_rdkit() {
        let input = format!(
            "{}spurious\n\n$$$$\n",
            single_atom_sdf_record("data-error", "C", false).trim_end_matches("$$$$")
        );
        let mut reader = SdfReader::with_params(
            std::io::Cursor::new(input.as_bytes()),
            SdfReadParams {
                strict_parsing: true,
                ..Default::default()
            },
        );

        let record = read_next_sdf_record(&mut reader).unwrap().unwrap();

        assert_eq!(record.molecule.properties().name(), Some("data-error"));
        assert!(record.data_fields.is_empty());
        assert!(!reader.end);
        assert!(read_next_sdf_record(&mut reader).unwrap().is_none());
        assert!(reader.end);
    }

    #[test]
    fn sdf_supplier_recovery_missing_delimiter_at_eof_returns_molecule_like_rdkit() {
        let input = single_atom_sdf_record("missing-delimiter", "C", false)
            .trim_end_matches("$$$$")
            .to_string();
        let mut reader = SdfReader::new(std::io::Cursor::new(input.as_bytes()));

        let record = read_next_sdf_record(&mut reader).unwrap().unwrap();

        assert_eq!(
            record.molecule.properties().name(),
            Some("missing-delimiter")
        );
        assert!(reader.end);
    }

    #[test]
    fn seek_to_next_sdf_record_scans_to_delimiter() {
        let mut reader = std::io::Cursor::new("bad\nrecord\n$$$$\nnext\n".as_bytes());
        seek_to_next_sdf_record(&mut reader, 0, SdfReadParams::default()).unwrap();
        let next = read_rdkit_line(&mut reader).unwrap().unwrap();
        assert_eq!(next, "next");
    }

    #[test]
    fn build_sdf_index_records_offsets_and_lengths_from_raw_records() {
        let input = "mol\nM  END\n$$$$\nnext\nM  END\n$$$$\n";
        let mut reader = std::io::Cursor::new(input.as_bytes());

        let metadata = build_sdf_index(&mut reader, SdfReadParams::default()).unwrap();

        assert_eq!(
            metadata,
            vec![
                SdfRecordMetadata {
                    index: 0,
                    byte_offset: 0,
                    byte_len: "mol\nM  END\n$$$$\n".len() as u64,
                    line_offset: 0,
                    line_len: 3,
                    title: Some("mol".to_string()),
                },
                SdfRecordMetadata {
                    index: 1,
                    byte_offset: "mol\nM  END\n$$$$\n".len() as u64,
                    byte_len: "next\nM  END\n$$$$\n".len() as u64,
                    line_offset: 3,
                    line_len: 3,
                    title: Some("next".to_string()),
                },
            ]
        );
    }

    #[test]
    fn sdf_dataset_random_access_metadata_out_of_range_and_invalid_records_like_rdkit() {
        let valid_first = single_atom_sdf_record("first", "C", true);
        let invalid = format!(
            "bad\n\n\n  1  0  0  0  0  0            999 V2000\n{}\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0)
        );
        let valid_second = single_atom_sdf_record("second", "N", true);
        let input = format!("{valid_first}{invalid}{valid_second}");
        let file = write_temp_sdf(&input);

        let dataset = SdfDataset::open(file.path()).unwrap();

        assert_eq!(dataset.len(), 3);
        assert_eq!(
            dataset.metadata(0),
            Some(&SdfRecordMetadata {
                index: 0,
                byte_offset: 0,
                byte_len: valid_first.len() as u64,
                line_offset: 0,
                line_len: valid_first.bytes().filter(|byte| *byte == b'\n').count(),
                title: Some("first".to_string()),
            })
        );
        assert_eq!(
            dataset.metadata(1),
            Some(&SdfRecordMetadata {
                index: 1,
                byte_offset: valid_first.len() as u64,
                byte_len: invalid.len() as u64,
                line_offset: valid_first.bytes().filter(|byte| *byte == b'\n').count(),
                line_len: invalid.bytes().filter(|byte| *byte == b'\n').count(),
                title: Some("bad".to_string()),
            })
        );

        let second = dataset.record(2).unwrap();
        assert_eq!(second.molecule.properties().name(), Some("second"));
        assert_eq!(second.molecule.num_atoms(), 1);

        let err = dataset.record(1).unwrap_err();
        assert!(matches!(err, SdfReadError::Parse(_)));
        let out_of_range = dataset.record(3).unwrap_err();
        assert_eq!(
            out_of_range,
            SdfReadError::Parse(
                "ERROR: Index error (idx = 3) :  we do not have enough mol blocks".to_string()
            )
        );
    }

    #[test]
    fn sdf_dataset_reuses_open_params_allows_overrides_and_batches_slices_like_rdkit() {
        let explicit_h = format!(
            "with-h\n\n\n  2  1  0  0  0  0            999 V2000\n{}\n{}\n{}\nM  END\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("H", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );
        let valid_second = single_atom_sdf_record("second", "O", true);
        let input = format!("{explicit_h}{valid_second}");
        let file = write_temp_sdf(&input);
        let open_params = SdfReadParams {
            remove_hs: false,
            ..Default::default()
        };

        let dataset = SdfDataset::open_with_params(file.path(), open_params).unwrap();

        let opened_record = dataset.record(0).unwrap();
        assert_eq!(opened_record.molecule.num_atoms(), 2);

        let overridden = dataset
            .record_with_params(
                0,
                SdfReadParams {
                    remove_hs: true,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(overridden.molecule.num_atoms(), 1);

        let batch = MoleculeBatch::read_sdf_dataset_with_params_and_options(
            &dataset,
            open_params,
            BatchErrorMode::KeepErrors,
            Some(1),
            Some(false),
        )
        .unwrap();
        assert_eq!(batch.len(), 2);
        assert_eq!(batch.valid_mask(), vec![true, true]);
    }

    #[test]
    fn tokenize_v3000_line_splits_spaces_tabs_quotes_and_parentheses_like_rdkit() {
        assert_eq!(
            tokenize_v3000_line("1 C 0.0 1.0 2.0 0", 1).unwrap(),
            vec!["1", "C", "0.0", "1.0", "2.0", "0"]
        );
        assert_eq!(
            tokenize_v3000_line("ATOM=(1 2 3) LABEL=\"hello world\" X", 2).unwrap(),
            vec!["ATOM=(1 2 3)", "ABEL=\"hello world", "X"]
        );
        assert_eq!(tokenize_v3000_line("\tA\t\tB", 3).unwrap(), vec!["A", "B"]);
    }

    #[test]
    fn tokenize_v3000_line_keeps_doubled_quotes_inside_quoted_token_like_rdkit() {
        assert_eq!(
            tokenize_v3000_line("LABEL=\"a\"\"b\" NEXT", 1).unwrap(),
            vec!["ABEL=\"a\"\"b", "NEXT"]
        );
    }

    #[test]
    fn get_v3000_line_strips_prefix_and_joins_continuations_like_rdkit() {
        let input = "M  V30 BEGIN ATOM-\nM  V30  1 C 0 0 0 0\n";
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 4;

        let line = get_v3000_line(&mut reader, &mut line_number, SdfReadParams::default()).unwrap();

        assert_eq!(line, "BEGIN ATOM 1 C 0 0 0 0");
        assert_eq!(line_number, 6);
    }

    #[test]
    fn get_v3000_line_rejects_missing_prefix_like_rdkit() {
        let input = "not v3000\n";
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 0;

        let err =
            get_v3000_line(&mut reader, &mut line_number, SdfReadParams::default()).unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse("Line 1 does not start with 'M  V30 '\n".to_string())
        );
    }

    #[test]
    fn parse_counts_line_reads_fixed_width_counts_and_v3000_like_rdkit() {
        let mut line = String::from("  2  1  0  0  1");
        while line.len() < 34 {
            line.push(' ');
        }
        line.push_str("V3000");

        let counts = parse_counts_line(&line, 4, SdfReadParams::default()).unwrap();

        assert_eq!(
            counts,
            CountsLine {
                atom_count: 2,
                bond_count: 1,
                chiral_flag: 1,
                ctab_version: CtabVersion::V3000,
                v3000_sgroup_count: 0,
                v3000_obj3d_count: 0,
            }
        );
    }

    #[test]
    fn parse_counts_line_reports_primary_count_conversion_like_rdkit() {
        let err = parse_counts_line(" xx  1", 7, SdfReadParams::default()).unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse("Cannot convert ' xx' to unsigned int on line 7".to_string())
        );
    }

    #[test]
    fn parse_counts_line_ignores_extra_field_conversion_failures_like_rdkit() {
        let counts = parse_counts_line("  2  1xxx", 4, SdfReadParams::default()).unwrap();

        assert_eq!(counts.atom_count, 2);
        assert_eq!(counts.bond_count, 1);
    }

    #[test]
    fn mol_from_mol_data_stream_reads_header_counts_and_v2000_dispatch_like_rdkit() {
        let input = format!(
            "stream-name\n  COSMolKit          3D\nstream-comment\n  1  0  0  0  1  0            999 V2000\n{}\nM  END\n",
            v2000_atom_line("C", 0, 0, 0, 0)
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 0;

        let parsed = mol_from_mol_data_stream(
            &mut reader,
            &mut line_number,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(parsed.molecule.properties().name(), Some("stream-name"));
        assert_eq!(
            parsed.molecule.prop("_MolFileInfoLine"),
            Some("  COSMolKit          3D")
        );
        assert_eq!(
            parsed.molecule.prop("_MolFileComments"),
            Some("stream-comment")
        );
        assert_eq!(parsed.molecule.prop("_MolFileChiralFlag"), Some("1"));
        assert_eq!(parsed.molecule.num_atoms(), 1);
        assert_eq!(parsed.molecule.num_bonds(), 0);
        assert_eq!(parsed.header.ctab_version, CtabVersion::V2000);
        assert_eq!(line_number, 6);
    }

    #[test]
    fn mol_from_mol_data_stream_reports_missing_counts_line_like_rdkit() {
        let mut reader = std::io::Cursor::new("name\ninfo\ncomments\n".as_bytes());
        let mut line_number = 0;

        let err = mol_from_mol_data_stream(&mut reader, &mut line_number, SdfReadParams::default())
            .unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse("Counts line too short: '' on line4".to_string())
        );
        assert_eq!(line_number, 4);
    }

    #[test]
    fn mol_from_mol_data_stream_applies_strict_counts_version_like_rdkit() {
        let mut bad_version_counts = "  1  0  0  0  0            999".to_string();
        while bad_version_counts.len() < 34 {
            bad_version_counts.push(' ');
        }
        bad_version_counts.push_str("X2000");
        let input = format!(
            "bad-version\ninfo\ncomments\n{bad_version_counts}\n{}\nM  END\n",
            v2000_atom_line("N", 0, 0, 0, 0)
        );

        let mut strict_reader = std::io::Cursor::new(input.as_bytes());
        let mut strict_line = 0;
        let strict_err = mol_from_mol_data_stream(
            &mut strict_reader,
            &mut strict_line,
            SdfReadParams {
                strict_parsing: true,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert_eq!(
            strict_err,
            SdfReadError::Parse("CTAB version string invalid at line 4".to_string())
        );

        let mut nonstrict_reader = std::io::Cursor::new(input.as_bytes());
        let mut nonstrict_line = 0;
        let parsed = mol_from_mol_data_stream(
            &mut nonstrict_reader,
            &mut nonstrict_line,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                strict_parsing: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(parsed.header.ctab_version, CtabVersion::V2000);
        assert_eq!(parsed.molecule.num_atoms(), 1);
    }

    #[test]
    fn mol_from_mol_data_stream_dispatches_v3000_and_checks_initial_counts_like_rdkit() {
        let counts_line = v3000_outer_counts_line(1, 0);
        let input = format!(
            "\
v3000-nonzero
info
comments
{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 1
M  V30 BEGIN ATOM
M  V30 1 O 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 END CTAB
M  END
"
        );

        let mut strict_reader = std::io::Cursor::new(input.as_bytes());
        let mut strict_line = 0;
        let strict_err = mol_from_mol_data_stream(
            &mut strict_reader,
            &mut strict_line,
            SdfReadParams::default(),
        )
        .unwrap_err();
        assert_eq!(
            strict_err,
            SdfReadError::Parse(
                "V3000 mol blocks should have 0s in the initial counts line. (line: 4)".to_string()
            )
        );

        let mut nonstrict_reader = std::io::Cursor::new(input.as_bytes());
        let mut nonstrict_line = 0;
        let parsed = mol_from_mol_data_stream(
            &mut nonstrict_reader,
            &mut nonstrict_line,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                strict_parsing: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(parsed.header.ctab_version, CtabVersion::V3000);
        assert_eq!(parsed.molecule.num_atoms(), 1);
        assert_eq!(parsed.molecule.prop("_MolFileChiralFlag"), Some("1"));
    }

    #[test]
    fn mol_from_mol_data_stream_reports_v3000_counts_conversion_like_rdkit() {
        let counts_line = v3000_outer_counts_line(0, 0);
        let input = format!(
            "\
v3000-bad-counts
info
comments
{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 nope 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 END CTAB
M  END
"
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 0;

        let err = mol_from_mol_data_stream(&mut reader, &mut line_number, SdfReadParams::default())
            .unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse("Cannot convert 'nope' to unsigned int on line 6".to_string())
        );
    }

    #[test]
    fn mol_from_mol_data_stream_reports_v2000_missing_m_end_filecomplete_like_rdkit() {
        let input = format!(
            "missing-mend\ninfo\ncomments\n  1  0  0  0  0  0            999 V2000\n{}\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0)
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 0;

        let err = mol_from_mol_data_stream(&mut reader, &mut line_number, SdfReadParams::default())
            .unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse(
                "Problems encountered parsing Mol data, M  END missing around line 6".to_string()
            )
        );
    }

    #[test]
    fn mol_from_mol_data_stream_reports_v3000_missing_m_end_filecomplete_like_rdkit() {
        let counts_line = v3000_outer_counts_line(0, 0);
        let input = format!(
            "\
missing-v3000-mend
info
comments
{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 END CTAB
"
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 0;

        let err = mol_from_mol_data_stream(&mut reader, &mut line_number, SdfReadParams::default())
            .unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse(
                "Problems encountered parsing Mol data, M  END missing around line 11".to_string()
            )
        );
        assert_eq!(line_number, 11);
    }

    #[test]
    fn mol_from_mol_data_stream_applies_v3000_strict_end_ctab_like_rdkit() {
        let counts_line = v3000_outer_counts_line(0, 0);
        let input = format!(
            "\
v3000-bad-end-ctab
info
comments
{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 NOT END CTAB
M  END
"
        );

        let mut strict_reader = std::io::Cursor::new(input.as_bytes());
        let mut strict_line = 0;
        let strict_err = mol_from_mol_data_stream(
            &mut strict_reader,
            &mut strict_line,
            SdfReadParams::default(),
        )
        .unwrap_err();
        assert_eq!(
            strict_err,
            SdfReadError::Parse("END CTAB line not found".to_string())
        );

        let mut nonstrict_reader = std::io::Cursor::new(input.as_bytes());
        let mut nonstrict_line = 0;
        let parsed = mol_from_mol_data_stream(
            &mut nonstrict_reader,
            &mut nonstrict_line,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                strict_parsing: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(parsed.molecule.num_atoms(), 1);
        assert_eq!(nonstrict_line, 11);
    }

    #[test]
    fn parse_v2000_ctab_builds_atoms_bonds_and_conformer_from_blocks() {
        let input = format!(
            "{}\n{}\n{}\nM  END\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("O", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 4;

        let parsed = parse_v2000_ctab(
            &mut reader,
            MolHeader {
                name: Some("ethanol-fragment".to_string()),
                info: String::new(),
                comments: String::new(),
                ctab_version: CtabVersion::V2000,
            },
            CountsLine {
                atom_count: 2,
                bond_count: 1,
                chiral_flag: 0,
                ctab_version: CtabVersion::V2000,
                v3000_sgroup_count: 0,
                v3000_obj3d_count: 0,
            },
            &mut line_number,
            SdfReadParams::default(),
        )
        .unwrap();

        assert_eq!(parsed.molecule.num_atoms(), 2);
        assert_eq!(parsed.molecule.num_bonds(), 1);
        assert_eq!(
            parsed.molecule.properties().name(),
            Some("ethanol-fragment")
        );
        assert_eq!(parsed.molecule.bonds()[0].order(), BondOrder::Single);
        assert_eq!(parsed.molecule.conformers_3d()[0].coordinates().len(), 2);
        assert_eq!(
            parsed.molecule.conformers_3d()[0].coordinates()[0],
            [1.25, -2.5, 0.75]
        );
        assert!(parsed.molecule.conformers_3d()[0].is_3d());
        assert_eq!(line_number, 8);
    }

    #[test]
    fn read_sdf_from_str_reads_v2000_record_and_data_fields() {
        let input = format!(
            "sample\n  COSMolKit          3D\ncomment\n  2  1  0  0  0  0            999 V2000\n{}\n{}\n{}\nM  END\n> <ID>\ncmpd-1\n\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("O", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.num_atoms(), 2);
        assert_eq!(record.molecule.num_bonds(), 1);
        assert_eq!(record.molecule.properties().name(), Some("sample"));
        assert_eq!(record.molecule.prop("ID"), Some("cmpd-1"));
        assert_eq!(
            record.data_fields,
            vec![("ID".to_string(), "cmpd-1".to_string())]
        );
    }

    #[test]
    fn sdf_data_fields_read_mol_props_named_blank_repeated_and_invalid_headers_like_rdkit() {
        let input = format!(
            "props\n  COSMolKit          2D\ncomment\n  2  1  0  0  0  0            999 V2000\n{}\n{}\n{}\nM  END\n> <ID>\nfirst\n\n> <BLANK>\n\n> <ID>\nsecond\n\n> <>\nignored\n\n> no-label\nignored\n\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("O", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );

        let record = read_sdf_from_str(&input).unwrap();

        assert_eq!(record.molecule.prop("ID"), Some("second"));
        assert_eq!(record.molecule.prop("BLANK"), Some(""));
        assert_eq!(record.molecule.prop(""), None);
        assert_eq!(
            record.data_fields,
            vec![
                ("ID".to_string(), "first".to_string()),
                ("BLANK".to_string(), String::new()),
                ("ID".to_string(), "second".to_string()),
            ]
        );
    }

    #[test]
    fn read_sdf_from_str_applies_atom_and_bond_property_lists_like_rdkit() {
        let input = format!(
            "props\n  COSMolKit          2D\ncomment\n  2  1  0  0  0  0            999 V2000\n{}\n{}\n{}\nM  END\n> <atom.prop.AtomLabel>\nC1 O1\n\n> <atom.iprop.Score>\n[?] 7 ?\n\n> <bond.prop.BondLabel>\nsingle\n\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("O", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );

        let record = read_sdf_from_str(&input).unwrap();

        assert_eq!(record.molecule.atoms()[0].prop("AtomLabel"), Some("C1"));
        assert_eq!(record.molecule.atoms()[1].prop("AtomLabel"), Some("O1"));
        assert_eq!(record.molecule.atoms()[0].prop("Score"), Some("7"));
        assert_eq!(record.molecule.atoms()[1].prop("Score"), None);
        assert_eq!(record.molecule.bonds()[0].prop("BondLabel"), Some("single"));
        assert_eq!(record.molecule.properties().sdf_property_lists().len(), 3);
    }

    #[test]
    fn sdf_data_fields_process_property_lists_atom_bond_and_disable_flag_like_rdkit() {
        let input = format!(
            "props\n  COSMolKit          2D\ncomment\n  3  2  0  0  0  0            999 V2000\n{}\n{}\n{}\n{}\n{}\nM  END\n> <atom.prop.Label>\nC1 O2 N3\n\n> <atom.iprop.Index>\n1 2 3\n\n> <atom.dprop.Partial>\n[?] 0.1 ? 0.3\n\n> <atom.bprop.Active>\n1 0 1\n\n> <bond.prop.BondLabel>\nfirst second\n\n> <bond.iprop.BondIndex>\n5 6\n\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("O", 0, 0, 0, 0),
            v2000_atom_line("N", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0),
            v2000_bond_line(2, 3, 1, 0, 0)
        );

        let record = read_sdf_from_str(&input).unwrap();

        assert_eq!(record.molecule.atoms()[0].prop("Label"), Some("C1"));
        assert_eq!(record.molecule.atoms()[1].prop("Index"), Some("2"));
        assert_eq!(record.molecule.atoms()[1].prop("Partial"), None);
        assert_eq!(record.molecule.atoms()[2].prop("Partial"), Some("0.3"));
        assert_eq!(record.molecule.atoms()[0].prop("Active"), Some("true"));
        assert_eq!(record.molecule.atoms()[1].prop("Active"), Some("false"));
        assert_eq!(record.molecule.bonds()[0].prop("BondLabel"), Some("first"));
        assert_eq!(record.molecule.bonds()[1].prop("BondIndex"), Some("6"));
        assert_eq!(record.molecule.properties().sdf_property_lists().len(), 6);

        let unprocessed = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                process_property_lists: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(
            unprocessed.molecule.prop("atom.prop.Label"),
            Some("C1 O2 N3")
        );
        assert_eq!(unprocessed.molecule.atoms()[0].prop("Label"), None);
        assert!(
            unprocessed
                .molecule
                .properties()
                .sdf_property_lists()
                .is_empty()
        );
    }

    #[test]
    fn read_sdf_property_lists_match_rdkit_bool_and_token_boundaries() {
        let input = format!(
            "props\n  COSMolKit          2D\ncomment\n  2  1  0  0  0  0            999 V2000\n{}\n{}\n{}\nM  END\n> <atom.bprop.Flag>\n1 0\n\n> <atom.bprop.TextFlag>\ntrue false\n\n> <atom.prop.LeadingSpace>\n C1 O1\n\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("O", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );

        let record = read_sdf_from_str(&input).unwrap();

        assert_eq!(record.molecule.atoms()[0].prop("Flag"), Some("true"));
        assert_eq!(record.molecule.atoms()[1].prop("Flag"), Some("false"));
        assert_eq!(record.molecule.atoms()[0].prop("TextFlag"), None);
        assert_eq!(record.molecule.atoms()[1].prop("TextFlag"), None);
        assert_eq!(record.molecule.atoms()[0].prop("LeadingSpace"), None);
        assert_eq!(record.molecule.atoms()[1].prop("LeadingSpace"), None);
        assert_eq!(record.molecule.properties().sdf_property_lists().len(), 2);
    }

    #[test]
    fn sdf_default_removes_simple_explicit_hydrogen_after_parse() {
        let input = format!(
            "methane-fragment\n  COSMolKit          2D\ncomment\n  2  1  0  0  0  0            999 V2000\n{}\n{}\n{}\nM  END\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("H", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );

        let record = read_sdf_from_str(&input).unwrap();

        assert_eq!(record.molecule.num_atoms(), 1);
        assert_eq!(record.molecule.num_bonds(), 0);
        assert_eq!(record.molecule.atoms()[0].atomic_number(), 6);
        assert_eq!(record.molecule.conformers_3d()[0].coordinates().len(), 1);
    }

    #[test]
    fn read_sdf_from_str_preserves_simple_explicit_hydrogen_when_remove_hs_false() {
        let input = format!(
            "methane-fragment\n  COSMolKit          2D\ncomment\n  2  1  0  0  0  0            999 V2000\n{}\n{}\n{}\nM  END\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("H", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.num_atoms(), 2);
        assert_eq!(record.molecule.num_bonds(), 1);
    }

    #[test]
    fn read_sdf_from_str_preserves_hydrogen_only_neighbors_like_rdkit_default() {
        let input = format!(
            "hydrogen\n  COSMolKit          2D\ncomment\n  2  1  0  0  0  0            999 V2000\n{}\n{}\n{}\nM  END\n$$$$\n",
            v2000_atom_line("H", 0, 0, 0, 0),
            v2000_atom_line("H", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.num_atoms(), 2);
        assert_eq!(record.molecule.num_bonds(), 1);
    }

    #[test]
    fn read_sdf_from_str_reads_simple_v3000_ctab() {
        let mut counts_line = "  0  0  0  0  0            999".to_string();
        while counts_line.len() < 34 {
            counts_line.push(' ');
        }
        counts_line.push_str("V3000");
        let input = format!(
            "\
v3000-simple
  COSMolKit

{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0 CHG=1 CFG=3 INVRET=2
M  V30 2 O 1.25 0.0 0.5 0 MASS=18 RAD=2
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2 CFG=2
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.num_atoms(), 2);
        assert_eq!(record.molecule.num_bonds(), 1);
        assert_eq!(record.molecule.atoms()[0].formal_charge(), 1);
        assert_eq!(record.molecule.atoms()[0].mol_parity(), Some(3));
        assert_eq!(record.molecule.atoms()[0].mol_inversion_flag(), Some(2));
        assert_eq!(record.molecule.atoms()[0].prop("_MolFileAtomCfg"), None);
        assert_eq!(record.molecule.atoms()[0].prop("molInversionFlag"), None);
        assert_eq!(record.molecule.atoms()[1].isotope(), Some(18));
        assert_eq!(record.molecule.atoms()[1].radical_electrons(), 1);
        assert_eq!(record.molecule.bonds()[0].order(), BondOrder::Double);
        assert_eq!(
            record.molecule.bonds()[0].direction(),
            BondDirection::EitherDouble
        );
        assert_eq!(record.molecule.bonds()[0].stereo(), BondStereo::Any);
        assert_eq!(
            record.molecule.conformers_3d()[0].coordinates(),
            &[[0.0, 0.0, 0.0], [1.25, 0.0, 0.5]]
        );
    }

    #[test]
    fn read_v3000_unknown_single_bond_sets_crossed_double_bond_controls_like_rdkit() {
        let mut counts_line = "  0  0  0  0  0            999".to_string();
        while counts_line.len() < 34 {
            counts_line.push(' ');
        }
        counts_line.push_str("V3000");
        let input = format!(
            "\
v3000-crossed-double-bond
  COSMolKit

{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.0 0.0 0.0 0
M  V30 2 C 1.0 0.0 0.0 0
M  V30 3 F -1.0 1.0 0.0 0
M  V30 4 Cl 1.0 -1.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 3 1 CFG=2
M  V30 3 1 2 4
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        let double_bond = &record.molecule.bonds()[0];
        assert_eq!(double_bond.stereo(), BondStereo::Any);
        assert_eq!(
            double_bond.stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
    }

    #[test]
    fn read_sdf_from_str_reads_v3000_rgroups_and_ignores_unknown_props_like_rdkit() {
        let mut counts_line = "  0  0  0  0  0            999".to_string();
        while counts_line.len() < 34 {
            counts_line.push(' ');
        }
        counts_line.push_str("V3000");
        let input = format!(
            "\
v3000-rgroup
  COSMolKit

{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 R# 0.0 0.0 0.0 0 RGROUPS=(1 12) FOO=BAR
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        let atom = &record.molecule.atoms()[0];
        assert_eq!(atom.prop("_MolFileRLabel"), Some("12"));
        assert_eq!(atom.prop("dummyLabel"), Some("R12"));
        assert_eq!(atom.isotope(), Some(12));
    }

    #[test]
    fn read_sdf_from_str_reads_v3000_bond_properties() {
        let mut counts_line = "  0  0  0  0  0            999".to_string();
        while counts_line.len() < 34 {
            counts_line.push(' ');
        }
        counts_line.push_str("V3000");
        let input = format!(
            "\
v3000-bond-props
  COSMolKit

{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0 STBOX=1
M  V30 2 O 1.0 0.0 0.0 0 STBOX=1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 RXCTR=4 STBOX=1 ENDPTS=(2 1 2) ATTACH=ALL
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"
        );

        let record = read_sdf_from_str(&input).unwrap();
        let bond = &record.molecule.bonds()[0];
        assert_eq!(bond.prop("molReactStatus"), Some("4"));
        assert_eq!(bond.prop("molStereoCare"), Some("1"));
        assert_eq!(bond.prop("_MolFileBondEndPts"), Some("(2 1 2)"));
        assert_eq!(bond.prop("_MolFileBondAttach"), Some("ALL"));
    }

    #[test]
    fn read_sdf_from_str_completes_v3000_ring_bond_count_scan() {
        let mut counts_line = "  0  0  0  0  0            999".to_string();
        while counts_line.len() < 34 {
            counts_line.push(' ');
        }
        counts_line.push_str("V3000");
        let input = format!(
            "\
v3000-ring-query
  COSMolKit

{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 3 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0 RBCNT=-2
M  V30 2 C 1.0 0.0 0.0 0 HCOUNT=2
M  V30 3 C 0.5 1.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 1
M  V30 END BOND
M  V30 END CTAB
M  END
"
        );

        let record = read_sdf_from_str(&input).unwrap();

        assert_eq!(record.molecule.prop("_NeedsQueryScan"), None);
        assert_eq!(
            record.molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
        );
        assert_eq!(
            record.molecule.atoms()[1].query(),
            Some(&QueryNode::predicate(
                AtomQueryPredicate::ImplicitHydrogenCountLessEqual(2)
            ))
        );
    }

    #[test]
    fn read_sdf_from_str_reads_v3000_sgroups_and_stereo_collections() {
        let mut counts_line = "  0  0  0  0  0            999".to_string();
        while counts_line.len() < 34 {
            counts_line.push(' ');
        }
        counts_line.push_str("V3000");
        let input = format!(
            "\
v3000-sgroup
  COSMolKit

{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 2 0 1
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 O 1.25 0.0 0.5 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 7 ATOMS=(2 1 2) XBONDS=(1 1) LABEL=Me CONNECT=HT BRKXYZ=(9 0 1 0 2 3 0 0 0 0)
M  V30 2 DAT 0 FIELDNAME=FIELD FIELDTYPE=T FIELDINFO=INFO QUERYTYPE=Q QUERYOP=OP FIELDDISP=\"display spec\" FIELDDATA=\"payload\" PARENT=1 COMPNO=5
M  V30 END SGROUP
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 1)
M  V30 MDLV30/STEREL2 ATOMS=(1 2)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$
"
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.substance_groups().len(), 2);
        let sup = &record.molecule.substance_groups()[0];
        assert_eq!(sup.kind(), &SubstanceGroupKind::Superatom);
        assert_eq!(sup.rdkit_sequence_id(), Some(1));
        assert_eq!(sup.external_id(), Some(7));
        assert_eq!(sup.atoms(), &[AtomId::new(0), AtomId::new(1)]);
        assert_eq!(sup.bonds(), &[BondId::new(0)]);
        assert_eq!(
            sup.bond_role(BondId::new(0)),
            crate::SGroupBondRole::Crossing
        );
        assert_eq!(sup.label(), Some("Me"));
        assert_eq!(sup.connection(), Some(&SGroupConnection::HeadToTail));
        assert_eq!(sup.display().unwrap().brackets[0].p1, [0.0, 1.0]);
        assert_eq!(sup.display().unwrap().brackets[0].p2, [2.0, 3.0]);

        let dat = &record.molecule.substance_groups()[1];
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

        assert_eq!(record.molecule.stereo_groups().len(), 2);
        assert_eq!(
            record.molecule.stereo_groups()[0].kind(),
            StereoGroupKind::Absolute
        );
        assert_eq!(
            record.molecule.stereo_groups()[0].atoms(),
            &[AtomId::new(0)]
        );
        assert_eq!(
            record.molecule.stereo_groups()[1].kind(),
            StereoGroupKind::Or
        );
        assert_eq!(record.molecule.stereo_groups()[1].id(), Some(2));
        assert_eq!(
            record.molecule.stereo_groups()[1].atoms(),
            &[AtomId::new(1)]
        );
    }

    #[test]
    fn read_sdf_from_str_rejects_v3000_sgroup_without_sgroup_block_like_rdkit() {
        let mut counts_line = "  0  0  0  0  0            999".to_string();
        while counts_line.len() < 34 {
            counts_line.push(' ');
        }
        counts_line.push_str("V3000");
        let input = format!(
            "\
v3000-sgroup-missing
  COSMolKit

{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 END ATOM
M  V30 END CTAB
M  END
"
        );

        let err = read_sdf_from_str(&input).unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse("BEGIN SGROUP line not found on line 10".to_string())
        );
    }

    #[test]
    fn sdf_reader_streams_multiple_records() {
        let record = format!(
            "mol\n\n\n  1  0  0  0  0  0            999 V2000\n{}\nM  END\n$$$$\n",
            v2000_atom_line("N", 0, 0, 0, 0)
        );
        let input = format!("{record}{record}");
        let mut reader = SdfReader::new(std::io::Cursor::new(input.as_bytes()));

        assert_eq!(
            reader.next_record().unwrap().unwrap().molecule.num_atoms(),
            1
        );
        assert_eq!(
            reader.next_record().unwrap().unwrap().molecule.num_atoms(),
            1
        );
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn parse_v2000_ctab_marks_only_bond_aromatic_for_molfile_type_4_like_rdkit() {
        let input = format!(
            "{}\n{}\n{}\nM  END\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 4, 0, 0)
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 4;

        let parsed = parse_v2000_ctab(
            &mut reader,
            MolHeader {
                name: None,
                info: String::new(),
                comments: String::new(),
                ctab_version: CtabVersion::V2000,
            },
            CountsLine {
                atom_count: 2,
                bond_count: 1,
                chiral_flag: 0,
                ctab_version: CtabVersion::V2000,
                v3000_sgroup_count: 0,
                v3000_obj3d_count: 0,
            },
            &mut line_number,
            SdfReadParams::default(),
        )
        .unwrap();

        assert_eq!(parsed.molecule.bonds()[0].order(), BondOrder::Aromatic);
        assert!(parsed.molecule.bonds()[0].is_aromatic());
        assert!(!parsed.molecule.atoms()[0].is_aromatic());
        assert!(!parsed.molecule.atoms()[1].is_aromatic());
    }

    #[test]
    fn finish_mol_processing_applies_mol_props_queries_and_valence() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_prop("molTotValence", "4")
                .with_prop("molSubstCount", "-2"),
        );
        let a1 = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let processed = finish_mol_processing(
            molecule,
            false,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        let atom = &processed.atoms()[0];
        assert!(atom.no_implicit());
        assert_eq!(atom.explicit_hydrogens(), 3);
        assert_eq!(atom.prop("molTotValence"), None);
        assert_eq!(
            atom.query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(1)))
        );
    }

    #[test]
    fn process_mol_props_finalization_keeps_v2000_query_props_like_rdkit() {
        let input = format!(
            "query-props\n  COSMolKit          2D\ncomment\n  3  3  0  0  1  0            999 V2000\n{}\n{}\n{}\n{}\n{}\n{}\nM  RBC  1   1  -2\nM  UNS  1   2   1\nM  END\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 5, 0, 0),
            v2000_bond_line(2, 3, 1, 0, 0),
            v2000_bond_line(3, 1, 1, 0, 0)
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.prop("_MolFileChiralFlag"), Some("1"));
        assert_eq!(
            record.molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
        );
        assert_eq!(
            record.molecule.atoms()[1].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::IsUnsaturated))
        );
        assert_eq!(
            record.molecule.bonds()[0].query(),
            Some(&QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ])))
        );
        assert_eq!(
            record.molecule.bonds()[0].prop("_MolFileBondQuery"),
            Some("1")
        );
    }

    #[test]
    fn process_sgroups_applies_zbo_zch_and_hyd_extensions() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::N));
        let b0 = builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_bonds(vec![b0])
                    .with_data(crate::SGroupData {
                        field_name: Some("ZBO".to_string()),
                        values: vec!["1".to_string()],
                        ..Default::default()
                    }),
            )
            .unwrap();
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                    .with_atoms(vec![a0, a1])
                    .with_data(crate::SGroupData {
                        field_name: Some("ZCH".to_string()),
                        values: vec!["1;-1".to_string()],
                        ..Default::default()
                    }),
            )
            .unwrap();
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Data)
                    .with_atoms(vec![a0, a1])
                    .with_data(crate::SGroupData {
                        field_name: Some("HYD".to_string()),
                        values: vec!["2;0".to_string()],
                        ..Default::default()
                    }),
            )
            .unwrap();
        let mut molecule = builder.build().unwrap();

        process_sgroups(&mut molecule, SdfReadParams::default()).unwrap();

        assert_eq!(molecule.bonds()[0].order(), BondOrder::Zero);
        assert_eq!(molecule.atoms()[0].formal_charge(), 1);
        assert_eq!(molecule.atoms()[1].formal_charge(), -1);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 2);
        assert_eq!(molecule.atoms()[1].explicit_hydrogens(), 0);
        assert_eq!(molecule.atoms()[0].prop("_ZBO_H"), Some("1"));
        assert!(molecule.substance_groups().is_empty());
    }

    #[test]
    fn smarts_consumer_sdf_smartsq() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_atoms(vec![a0])
                    .with_data(crate::SGroupData {
                        query_type: Some("SMARTSQ".to_string()),
                        query_op: Some("=".to_string()),
                        values: vec!["[#6]".to_string()],
                        ..Default::default()
                    }),
            )
            .unwrap();
        let mut molecule = builder.build().unwrap();

        process_sgroups(&mut molecule, SdfReadParams::default()).unwrap();

        assert!(matches!(
            molecule.atoms()[0].query(),
            Some(QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(_)))
        ));
        assert!(molecule.substance_groups().is_empty());
    }

    #[test]
    fn process_sgroups_ignores_non_equals_smartsq_queryop_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_atoms(vec![a0])
                    .with_data(crate::SGroupData {
                        query_type: Some("SMARTSQ".to_string()),
                        query_op: Some("!=".to_string()),
                        values: vec!["[#6]".to_string()],
                        ..Default::default()
                    }),
            )
            .unwrap();
        let mut molecule = builder.build().unwrap();

        process_sgroups(&mut molecule, SdfReadParams::default()).unwrap();

        assert_eq!(molecule.substance_groups().len(), 0);
        assert_eq!(molecule.atoms()[0].query(), None);
    }

    #[test]
    fn process_mrv_coordinate_bond_zero_index_does_not_rewrite_first_bond_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::O));
        let b0 = builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Unspecified))
            .unwrap();
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_bonds(vec![b0])
                    .with_data(crate::SGroupData {
                        field_name: Some("MRV_COORDINATE_BOND_TYPE".to_string()),
                        values: vec!["0".to_string()],
                        ..Default::default()
                    }),
            )
            .unwrap();
        let mut molecule = builder.build().unwrap();

        process_sgroups(&mut molecule, SdfReadParams::default()).unwrap();

        assert_eq!(molecule.bonds()[0].order(), BondOrder::Unspecified);
    }

    #[test]
    fn expand_attachment_points_materializes_v2000_apo_like_rdkit() {
        let input = format!(
            "attachment\n  COSMolKit          2D\ncomment\n  1  0  0  0  0  0            999 V2000\n{0}\nM  APO  1   1   1\nM  END\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0)
        );

        let result = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                expand_attachment_points: true,
                ..Default::default()
            },
        )
        .unwrap();

        let atom = &result.molecule.atoms()[0];
        assert_eq!(atom.element().atomic_number(), 6); // Carbon
        assert_eq!(atom.prop("molAttachPoint"), None);
        assert_eq!(result.molecule.num_atoms(), 2);
        assert_eq!(result.molecule.num_bonds(), 1);
        let dummy = &result.molecule.atoms()[1];
        assert_eq!(dummy.element().atomic_number(), 0);
        assert_eq!(dummy.prop("_fromAttachPoint"), Some("1"));
        assert_eq!(
            dummy.query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::Any))
        );
        assert_eq!(result.molecule.bonds()[0].begin().index(), 0);
        assert_eq!(result.molecule.bonds()[0].end().index(), 1);
        assert_eq!(result.molecule.bonds()[0].order(), BondOrder::Single);
    }

    #[test]
    fn expand_attachment_points_respects_sdf_read_param_like_rdkit() {
        let input = format!(
            "attachment-off\n  COSMolKit          2D\ncomment\n  1  0  0  0  0  0            999 V2000\n{0}\nM  APO  1   1   1\nM  END\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0)
        );

        let result = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                expand_attachment_points: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(result.molecule.num_atoms(), 1);
        assert_eq!(result.molecule.num_bonds(), 0);
        assert_eq!(result.molecule.atoms()[0].prop("molAttachPoint"), Some("1"));
    }

    #[test]
    fn expand_attachment_points_minus_one_adds_two_query_dummies_like_rdkit() {
        let input = format!(
            "attachment-two\n  COSMolKit          2D\ncomment\n  1  0  0  0  0  0            999 V2000\n{0}\nM  APO  1   1   3\nM  END\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0)
        );

        let result = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                expand_attachment_points: true,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(result.molecule.num_atoms(), 3);
        assert_eq!(result.molecule.num_bonds(), 2);
        assert_eq!(result.molecule.atoms()[0].prop("molAttachPoint"), None);
        assert_eq!(
            result.molecule.atoms()[1].prop("_fromAttachPoint"),
            Some("1")
        );
        assert_eq!(
            result.molecule.atoms()[2].prop("_fromAttachPoint"),
            Some("2")
        );
        assert_eq!(
            result.molecule.atoms()[1].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::Any))
        );
        assert_eq!(
            result.molecule.atoms()[2].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::Any))
        );
    }

    #[test]
    fn expand_attachment_points_invalid_value_is_not_expanded_like_rdkit() {
        let input = "\
attachment-invalid
  COSMolKit          2D
comment
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0 ATTCHPT=4
M  V30 END ATOM
M  V30 END CTAB
M  END
$$$$
";

        let result = read_sdf_from_str_with_params(
            input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                expand_attachment_points: true,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(result.molecule.num_atoms(), 1);
        assert_eq!(result.molecule.num_bonds(), 0);
        assert_eq!(result.molecule.atoms()[0].prop("molAttachPoint"), Some("4"));
    }

    #[test]
    fn expand_attachment_points_creates_2d_coordinates_like_rdkit() {
        let input = "\
attachment-v3000
  COSMolKit          2D
comment
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.7083 3.25 0 0
M  V30 2 O -3.3747 4.02 0 0 ATTCHPT=1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
";

        let result = read_sdf_from_str_with_params(
            input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                expand_attachment_points: true,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(result.molecule.num_atoms(), 3);
        assert!(result.molecule.conformers_3d().is_empty());
        let conformer = result.molecule.conformers_2d().first().unwrap();
        assert_eq!(conformer.coordinates().len(), 3);
        let attachment = conformer.coordinates()[2];
        assert!((attachment[0] - -2.50869).abs() < 1.0e-4);
        assert!((attachment[1] - 4.52002).abs() < 1.0e-4);
        assert_eq!(
            result.molecule.source_coordinate_dim(),
            Some(crate::CoordinateDimension::TwoD)
        );
    }

    #[test]
    fn sdf_sanitize_accepts_aromatic_ring_after_sanitize_port() {
        let input = format!(
            "aromatic\n  COSMolKit          2D\ncomment\n  6  6  0  0  0  0            999 V2000\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\nM  END\n$$$$\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 4, 0, 0),
            v2000_bond_line(2, 3, 4, 0, 0),
            v2000_bond_line(3, 4, 4, 0, 0),
            v2000_bond_line(4, 5, 4, 0, 0),
            v2000_bond_line(5, 6, 4, 0, 0),
            v2000_bond_line(6, 1, 4, 0, 0)
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.num_atoms(), 6);
        assert_eq!(record.molecule.num_bonds(), 6);
        assert!(
            record
                .molecule
                .bonds()
                .iter()
                .all(|bond| bond.order() == BondOrder::Aromatic)
        );
    }

    #[test]
    fn sdf_sanitize_after_remove_hs_aromatizes_kekule_benzene_from_v3000() {
        let input = concat!(
            "\n",
            "     RDKit          2D\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 6 6 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 1.500000 0.000000 0.000000 0\n",
            "M  V30 2 C 0.750000 -1.299038 0.000000 0\n",
            "M  V30 3 C -0.750000 -1.299038 0.000000 0\n",
            "M  V30 4 C -1.500000 0.000000 0.000000 0\n",
            "M  V30 5 C -0.750000 1.299038 0.000000 0\n",
            "M  V30 6 C 0.750000 1.299038 0.000000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 1 2\n",
            "M  V30 2 2 2 3\n",
            "M  V30 3 1 3 4\n",
            "M  V30 4 2 4 5\n",
            "M  V30 5 1 5 6\n",
            "M  V30 6 2 6 1\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
            "$$$$\n",
        );

        let record = read_sdf_from_str_with_params(
            input,
            SdfReadParams {
                sanitize: true,
                remove_hs: true,
                ..Default::default()
            },
        )
        .unwrap();

        assert!(
            record
                .molecule
                .atoms()
                .iter()
                .all(|atom| atom.is_aromatic())
        );
        assert!(
            record
                .molecule
                .bonds()
                .iter()
                .all(|bond| bond.is_aromatic() && bond.order() == BondOrder::Aromatic)
        );
    }

    #[test]
    fn sdf_assigns_stereogenic_double_bond_detection_from_shared_3d_path() {
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let c1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let f = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let cl = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        builder
            .add_bond(BondSpec::new(c0, c1, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c0, f, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, cl, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [-1.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-1.0, 1.0, 0.0],
                    [1.0, -1.0, 0.0],
                ],
                true,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = detect_bond_stereochemistry(
            molecule,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::None);
        assert!(matches!(
            molecule.bonds()[1].direction(),
            BondDirection::EndUpRight | BondDirection::EndDownRight
        ));
        assert!(matches!(
            molecule.bonds()[2].direction(),
            BondDirection::EndUpRight | BondDirection::EndDownRight
        ));
    }

    #[test]
    fn sdf_assigns_wedged_bond_chirality_from_shared_2d_path() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        builder
            .add_bond(
                BondSpec::new(center, fluorine, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, -1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule =
            assign_chiral_types_from_bond_dirs(molecule, SdfReadParams::default()).unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    }

    #[test]
    fn sdf_detect_bond_stereochemistry_uses_first_2d_conformer_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let c1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let f = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let cl = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        builder
            .add_bond(BondSpec::new(c0, c1, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c0, f, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, cl, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [-1.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-1.0, 1.0, 0.0],
                    [1.0, -1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = detect_bond_stereochemistry(
            molecule,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert!(matches!(
            molecule.bonds()[1].direction(),
            BondDirection::EndUpRight | BondDirection::EndDownRight
        ));
        assert!(matches!(
            molecule.bonds()[2].direction(),
            BondDirection::EndUpRight | BondDirection::EndDownRight
        ));
    }

    #[test]
    fn sdf_reader_preserves_2d_coords_in_coords_2d_like_rdkit() {
        let input = concat!(
            "\n",
            "     RDKit          2D\n",
            "\n",
            "  4  3  0  0  0  0  0  0  0  0999 V2000\n",
            "   -1.9796   -0.1365    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "   -0.5994    0.4508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    0.5994   -0.4508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    1.9796    0.1365    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  2  1  0\n",
            "  2  3  2  0\n",
            "  3  4  1  0\n",
            "M  END\n",
            "$$$$\n",
        );

        let record = read_sdf_from_str(input).unwrap();
        let coords_2d = record
            .molecule
            .coordinates_2d()
            .expect("2D coords should be stored");
        assert_eq!(
            coords_2d,
            &[
                [-1.9796, -0.1365],
                [-0.5994, 0.4508],
                [0.5994, -0.4508],
                [1.9796, 0.1365],
            ]
        );
        assert_eq!(record.molecule.conformers_2d().len(), 1);
        assert!(record.molecule.conformers_3d().is_empty());
        assert_eq!(
            record.molecule.source_coordinate_dim(),
            Some(crate::CoordinateDimension::TwoD)
        );
    }

    #[test]
    fn sdf_reader_coordinate_mode_stores_each_source_conformer_once() {
        let input = concat!(
            "zero-z\n",
            "  PubChem          2D\n",
            "\n",
            "  2  1  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  2  1  0\n",
            "M  END\n",
            "$$$$\n",
        );

        for mode in [SdfCoordinateMode::Preserve, SdfCoordinateMode::Require2D] {
            let record = read_sdf_from_str_with_params(
                input,
                SdfReadParams {
                    coordinate_mode: mode,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(record.molecule.conformers_2d().len(), 1);
            assert!(record.molecule.conformers_3d().is_empty());
            assert_eq!(
                record.molecule.source_coordinate_dim(),
                Some(crate::CoordinateDimension::TwoD)
            );
        }

        let record = read_sdf_from_str_with_params(
            input,
            SdfReadParams {
                coordinate_mode: SdfCoordinateMode::Require3D,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(record.molecule.conformers_2d().is_empty());
        assert_eq!(record.molecule.conformers_3d().len(), 1);
        assert!(record.molecule.conformers_3d()[0].is_3d());
        assert_eq!(
            record.molecule.source_coordinate_dim(),
            Some(crate::CoordinateDimension::ThreeD)
        );
    }

    #[test]
    fn sdf_reader_keeps_3d_source_dim_for_zero_z_v3000_like_rdkit_calculate3dflag() {
        let input = concat!(
            "\n",
            "     RDKit          3D\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 3 2 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 1.233224 -0.184393 0.000000 0\n",
            "M  V30 2 C -0.112065 0.430195 -0.000000 0\n",
            "M  V30 3 O -1.121159 -0.245802 -0.000000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 1 2\n",
            "M  V30 2 2 2 3\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
            "$$$$\n",
        );

        let record = read_sdf_from_str(input).unwrap();
        let conformer = record
            .molecule
            .conformers_3d()
            .first()
            .expect("3D source should preserve conformer");
        assert!(conformer.is_3d());
        assert!(conformer.coordinates()[1][2].is_sign_negative());
        assert!(conformer.coordinates()[2][2].is_sign_negative());
        assert_eq!(
            record.molecule.source_coordinate_dim(),
            Some(crate::CoordinateDimension::ThreeD)
        );
    }

    #[test]
    fn sdf_reader_treats_v3000_r_token_as_dummy_label_not_rgp_like_rdkit() {
        let input = concat!(
            "\n",
            "     RDKit          2D\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 R -0.750000 0.000000 0.000000 1 VAL=1\n",
            "M  V30 2 C 0.750000 -0.000000 0.000000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 1 2\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
            "$$$$\n",
        );

        let record = read_sdf_from_str(input).unwrap();
        let atom = &record.molecule.atoms()[0];
        assert_eq!(atom.atomic_number(), 0);
        assert_eq!(atom.prop("dummyLabel"), Some("R"));
        assert_eq!(atom.prop("_MolFileRLabel"), None);
    }

    #[test]
    fn sdf_finish_processing_assigns_double_bond_stereo_from_2d_coords_like_rdkit() {
        let input = concat!(
            "\n",
            "     RDKit          2D\n",
            "\n",
            "  4  3  0  0  0  0  0  0  0  0999 V2000\n",
            "   -1.9796   -0.1365    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "   -0.5994    0.4508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    0.5994   -0.4508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    1.9796    0.1365    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  2  1  0\n",
            "  2  3  2  0\n",
            "  3  4  1  0\n",
            "M  END\n",
            "$$$$\n",
        );

        let record = read_sdf_from_str(input).unwrap();
        assert_ne!(record.molecule.bonds()[0].direction(), BondDirection::None);
        assert_eq!(
            record.molecule.bonds()[0].direction(),
            record.molecule.bonds()[2].direction()
        );
        assert_eq!(record.molecule.bonds()[1].stereo(), BondStereo::E);
        assert_eq!(
            record.molecule.bonds()[1].stereo_atoms(),
            Some([AtomId::new(0), AtomId::new(3)])
        );
    }

    #[test]
    fn molblock_reader_assigns_wedged_atom_chirality_from_2d_bond_dirs_like_rdkit() {
        let input = format!(
            "wedged\n  COSMolKit          2D\ncomment\n  4  3  0  0  0  0            999 V2000\n{}\n{}\n{}\n{}\n{}\n{}\n{}\nM  END\n$$$$\n",
            v2000_atom_line_at("C", 0.0, 0.0, 0.0),
            v2000_atom_line_at("F", 1.0, 0.0, 0.0),
            v2000_atom_line_at("Cl", 0.0, 1.0, 0.0),
            v2000_atom_line_at("Br", -1.0, -1.0, 0.0),
            v2000_bond_line(1, 2, 1, 1, 0),
            v2000_bond_line(1, 3, 1, 0, 0),
            v2000_bond_line(1, 4, 1, 0, 0),
        );

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_ne!(
            record.molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(record.molecule.atoms()[0].explicit_hydrogens(), 1);
        assert_eq!(record.molecule.bonds()[0].direction(), BondDirection::None);
    }

    #[test]
    fn assign_chiral_types_from_bond_dirs_reads_2d_wedge_with_sanitize_false() {
        let input = v2000_tetrahedral_wedge_block(1, false);

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_ne!(
            record.molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(record.molecule.atoms()[0].explicit_hydrogens(), 1);
        assert_eq!(record.molecule.bonds()[0].direction(), BondDirection::None);
    }

    #[test]
    fn assign_chiral_types_from_bond_dirs_reads_2d_dash_with_sanitize_false() {
        let wedged = read_sdf_from_str_with_params(
            &v2000_tetrahedral_wedge_block(1, false),
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        let dashed = read_sdf_from_str_with_params(
            &v2000_tetrahedral_wedge_block(6, false),
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_ne!(
            dashed.molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_ne!(
            wedged.molecule.atoms()[0].chiral_tag(),
            dashed.molecule.atoms()[0].chiral_tag()
        );
        assert_eq!(dashed.molecule.atoms()[0].explicit_hydrogens(), 1);
        assert_eq!(dashed.molecule.bonds()[0].direction(), BondDirection::None);
    }

    #[test]
    fn assign_chiral_types_from_bond_dirs_handles_explicit_wedged_hydrogen() {
        let input = v2000_tetrahedral_wedge_block(1, true);

        let record = read_sdf_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.num_atoms(), 5);
        assert_eq!(record.molecule.atoms()[1].atomic_number(), 1);
        assert_ne!(
            record.molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(record.molecule.atoms()[0].explicit_hydrogens(), 0);
    }

    #[test]
    fn assign_chiral_types_from_bond_dirs_promotes_implicit_hydrogen() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        builder
            .add_bond(
                BondSpec::new(center, fluorine, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, -1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule =
            assign_chiral_types_from_bond_dirs(molecule, SdfReadParams::default()).unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    }

    #[test]
    fn assign_chiral_types_from_bond_dirs_replace_existing_tags_true_clears_unknown() {
        let mut builder = MoleculeBuilder::new();
        let center = builder
            .add_atom(AtomSpec::new(Element::C).with_chiral_tag(crate::ChiralTag::TetrahedralCw));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        builder
            .add_bond(
                BondSpec::new(center, fluorine, BondOrder::Single)
                    .with_direction(BondDirection::Unknown),
            )
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                false,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule =
            assign_chiral_types_from_bond_dirs(molecule, SdfReadParams::default()).unwrap();

        assert_eq!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
    }

    #[test]
    fn sdf_assigns_3d_tetrahedral_chirality_from_shared_3d_path() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        builder
            .add_bond(BondSpec::new(center, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                true,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = assign_chiral_types_from_3d(
            molecule,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
    }

    fn tetrahedral_3d_for_assign_chiral_types_from_3d(center: AtomSpec, is_3d: bool) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(center);
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        builder
            .add_bond(BondSpec::new(center, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                is_3d,
            ))
            .unwrap();
        builder.build().unwrap()
    }

    fn square_planar_3d_for_assign_chiral_types_from_3d() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::P));
        for _ in 0..4 {
            let ligand = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
            builder
                .add_bond(BondSpec::new(center, ligand, BondOrder::Single))
                .unwrap();
        }
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, -1.0, 0.0],
                ],
                true,
            ))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn assign_chiral_types_from_3d_reads_tetrahedral_center_with_sanitize_false() {
        let molecule =
            tetrahedral_3d_for_assign_chiral_types_from_3d(AtomSpec::new(Element::C), true);

        let molecule = assign_chiral_types_from_3d(
            molecule,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(
            molecule.atoms()[0].prop("_NonExplicit3DChirality"),
            Some("1")
        );
    }

    #[test]
    fn assign_chiral_types_from_3d_uses_implicit_hydrogen_like_rdkit() {
        let molecule =
            tetrahedral_3d_for_assign_chiral_types_from_3d(AtomSpec::new(Element::C), true);

        let molecule = assign_chiral_types_from_3d(molecule, SdfReadParams::default()).unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 0);
    }

    #[test]
    fn assign_chiral_types_from_3d_replace_existing_tags_true_reassigns_explicit_tag() {
        let molecule = tetrahedral_3d_for_assign_chiral_types_from_3d(
            AtomSpec::new(Element::C).with_chiral_tag(crate::ChiralTag::TetrahedralCw),
            true,
        );

        let molecule = assign_chiral_types_from_3d(molecule, SdfReadParams::default()).unwrap();

        assert_eq!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::TetrahedralCcw
        );
        assert_eq!(molecule.atoms()[0].prop("_NonExplicit3DChirality"), None);
    }

    #[test]
    fn assign_chiral_types_from_3d_skips_false_3d_conformer_like_rdkit() {
        let molecule = tetrahedral_3d_for_assign_chiral_types_from_3d(
            AtomSpec::new(Element::C).with_chiral_tag(crate::ChiralTag::TetrahedralCw),
            false,
        );

        let molecule = assign_chiral_types_from_3d(molecule, SdfReadParams::default()).unwrap();

        assert_eq!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::TetrahedralCw
        );
    }

    #[test]
    fn assign_chiral_types_from_3d_reads_non_tetrahedral_center_like_rdkit() {
        let molecule = square_planar_3d_for_assign_chiral_types_from_3d();

        let molecule = assign_chiral_types_from_3d(molecule, SdfReadParams::default()).unwrap();

        assert_eq!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::SquarePlanar
        );
        assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(2));
        assert_eq!(
            molecule.atoms()[0].prop("_NonExplicit3DChirality"),
            Some("1")
        );
    }

    fn v3000_atropisomer_block(
        name: &str,
        header_dim: &str,
        atoms: &[(&str, f64, f64, f64)],
        bonds: &[(u32, u32, u32, Option<u32>)],
    ) -> String {
        let mut block = format!(
            "{name}\n  COSMolKit          {header_dim}\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {} {} 0 0 0\nM  V30 BEGIN ATOM\n",
            atoms.len(),
            bonds.len()
        );
        for (idx, (symbol, x, y, z)) in atoms.iter().enumerate() {
            block.push_str(&format!(
                "M  V30 {} {} {:.6} {:.6} {:.6} 0\n",
                idx + 1,
                symbol,
                x,
                y,
                z
            ));
        }
        block.push_str("M  V30 END ATOM\nM  V30 BEGIN BOND\n");
        for (idx, (order, begin, end, cfg)) in bonds.iter().enumerate() {
            block.push_str(&format!("M  V30 {} {} {} {}", idx + 1, order, begin, end));
            if let Some(cfg) = cfg {
                block.push_str(&format!(" CFG={cfg}"));
            }
            block.push('\n');
        }
        block.push_str("M  V30 END BOND\nM  V30 END CTAB\nM  END\n$$$$\n");
        block
    }

    fn minimal_2d_atropisomer_block(axis_cfg: u32) -> String {
        v3000_atropisomer_block(
            "atropisomer-2d",
            "2D",
            &[
                ("Cl", 0.0, -1.0, 0.0),
                ("C", 0.0, 0.0, 0.0),
                ("C", 0.0, 1.0, 0.0),
                ("C", 1.0, 0.0, 0.0),
                ("C", 1.0, -1.0, 0.0),
            ],
            &[
                (1, 2, 1, Some(axis_cfg)),
                (2, 2, 3, None),
                (1, 2, 4, None),
                (2, 4, 5, None),
            ],
        )
    }

    fn minimal_3d_atropisomer_block() -> String {
        v3000_atropisomer_block(
            "atropisomer-3d",
            "3D",
            &[
                ("Cl", 0.0, 0.0, 1.0),
                ("C", 0.0, 0.0, 0.0),
                ("C", 0.0, 1.0, 0.0),
                ("C", 1.0, 0.0, 0.0),
                ("C", 1.0, -1.0, 0.0),
            ],
            &[
                (1, 2, 1, Some(1)),
                (2, 2, 3, None),
                (1, 2, 4, None),
                (2, 4, 5, None),
            ],
        )
    }

    fn inconsistent_2d_atropisomer_block() -> String {
        v3000_atropisomer_block(
            "atropisomer-inconsistent",
            "2D",
            &[
                ("Cl", 0.0, -1.0, 0.0),
                ("C", 0.0, 0.0, 0.0),
                ("C", 0.0, 1.0, 0.0),
                ("C", 1.0, 0.0, 0.0),
                ("C", 1.0, -1.0, 0.0),
            ],
            &[
                (1, 2, 1, Some(1)),
                (2, 2, 3, Some(1)),
                (1, 2, 4, None),
                (2, 4, 5, None),
            ],
        )
    }

    fn non_atropisomer_block() -> String {
        v3000_atropisomer_block(
            "not-atropisomer",
            "2D",
            &[
                ("Cl", 0.0, -1.0, 0.0),
                ("C", 0.0, 0.0, 0.0),
                ("C", 0.0, 1.0, 0.0),
                ("O", 1.0, 0.0, 0.0),
            ],
            &[(1, 2, 1, Some(1)), (2, 2, 3, None), (1, 2, 4, None)],
        )
    }

    #[test]
    fn detect_atropisomer_chirality_assigns_2d_wedge_like_rdkit() {
        let record = read_sdf_from_str_with_params(
            &minimal_2d_atropisomer_block(1),
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.bonds()[2].stereo(), BondStereo::AtropCw);
    }

    #[test]
    fn detect_atropisomer_chirality_flips_2d_hash_like_rdkit() {
        let record = read_sdf_from_str_with_params(
            &minimal_2d_atropisomer_block(3),
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.bonds()[2].stereo(), BondStereo::AtropCcw);
    }

    #[test]
    fn detect_atropisomer_chirality_assigns_3d_geometry_like_rdkit() {
        let record = read_sdf_from_str_with_params(
            &minimal_3d_atropisomer_block(),
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.bonds()[2].stereo(), BondStereo::AtropCw);
    }

    #[test]
    fn detect_atropisomer_chirality_rejects_non_atropisomer_like_rdkit() {
        let record = read_sdf_from_str_with_params(
            &non_atropisomer_block(),
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.bonds()[2].stereo(), BondStereo::None);
    }

    #[test]
    fn detect_atropisomer_chirality_rejects_inconsistent_wedging_like_rdkit() {
        let record = read_sdf_from_str_with_params(
            &inconsistent_2d_atropisomer_block(),
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.bonds()[2].stereo(), BondStereo::None);
    }

    #[test]
    fn detect_atropisomer_chirality_matches_rdkit_fixture_macrocycle() {
        let input = include_str!(
            "../../../../testdata/rdkit_builtin/fixtures/Code/GraphMol/FileParsers/atropisomers/macrocycle-5-meta-Cl-ortho-hash.sdf"
        );
        let record = read_sdf_from_str_with_params(
            input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.bonds()[15].stereo(), BondStereo::AtropCw);
    }

    #[test]
    fn sdf_finish_processing_clears_false_3d_tetrahedral_tag_like_rdkit() {
        let input = concat!(
            "\n",
            "     RDKit          3D\n",
            "\n",
            "  4  3  0  0  0  0  0  0  0  0999 V2000\n",
            "   -1.1465   -0.8248   -0.1082 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    0.0043    0.0355    0.3412 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "    1.3209   -0.5515   -0.1115 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "   -0.1786    1.3409   -0.1215 O   0  0  0  0  0  0  0  0  0  0  0  0\n",
            "  1  2  1  0\n",
            "  2  3  1  0\n",
            "  2  4  1  0\n",
            "M  END\n",
            "$$$$\n",
        );

        let record = read_sdf_from_str(input).unwrap();
        assert_eq!(
            record.molecule.atoms()[1].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
    }

    #[test]
    fn parse_v2000_property_block_ignores_unknown_m_and_s_lines_like_rdkit() {
        let mut reader =
            std::io::Cursor::new("M  XYZ  1 payload\nS  XXX  1 payload\nM  END\n".as_bytes());
        let mut line_number = 8;
        let mut atoms = vec![
            parse_v2000_atom_line(
                &v2000_atom_line("C", 0, 0, 0, 0),
                1,
                SdfReadParams::default(),
            )
            .unwrap(),
        ];
        let mut bonds = Vec::new();
        let mut state = V2000PropertyState::default();

        let file_complete = parse_v2000_property_block(
            &mut reader,
            &mut line_number,
            SdfReadParams::default(),
            &mut atoms,
            &mut bonds,
            &mut state,
        )
        .unwrap();

        assert!(file_complete);
        assert_eq!(line_number, 11);
    }

    #[test]
    fn smarts_consumer_sdf_smarts_line() {
        let input = concat!(
            "M  PXA   1 payload\n",
            "M  RGP  1   2   7\n",
            "M  RBC  1   1  -2\n",
            "M  SUB  1   2   6\n",
            "M  UNS  1   1   1\n",
            "M  MRV SMA   2 [#6]\n",
            "M  LIN  1   1   2   2    \n",
            "M  END\n",
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 8;
        let mut atoms = vec![
            parse_v2000_atom_line(
                &v2000_atom_line("C", 0, 0, 0, 0),
                1,
                SdfReadParams::default(),
            )
            .unwrap(),
            parse_v2000_atom_line(
                &v2000_atom_line("R", 0, 0, 0, 0),
                2,
                SdfReadParams::default(),
            )
            .unwrap(),
        ];
        let mut bonds = Vec::new();
        let mut state = V2000PropertyState::default();

        let file_complete = parse_v2000_property_block(
            &mut reader,
            &mut line_number,
            SdfReadParams::default(),
            &mut atoms,
            &mut bonds,
            &mut state,
        )
        .unwrap();

        assert!(file_complete);
        assert_eq!(atoms[0].spec.prop("_MolFile_PXA"), Some(" payload"));
        assert_eq!(atoms[0].spec.prop("molRingBondCount"), Some("-2"));
        assert_eq!(atoms[1].spec.prop("_MolFileRLabel"), Some("7"));
        assert_eq!(atoms[1].spec.prop("dummyLabel"), Some("R7"));
        assert_eq!(atoms[1].spec.isotope(), Some(7));
        assert_eq!(atoms[1].spec.prop("molSubstCount"), Some("6"));
        assert_eq!(atoms[1].spec.prop("MRV_SMA"), Some("[#6]"));
        assert_eq!(
            state.molecule_props.get("_NeedsQueryScan"),
            Some(&"1".to_string())
        );
        assert_eq!(
            state.molecule_props.get("_MolFileLinkNodes"),
            Some(&"1 2 1 1 2".to_string())
        );
        assert_eq!(line_number, 16);
    }

    #[test]
    fn parse_v2000_property_block_applies_charge_isotope_radical_and_atom_list() {
        let input = concat!(
            "M  CHG  1   1  -1\n",
            "M  ISO  1   2  18\n",
            "M  RAD  1   2   2\n",
            "M  ALS   1  2 F C   N   \n",
            "M  APO  1   1   3\n",
            "M  END\n",
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 8;
        let mut atoms = vec![
            parse_v2000_atom_line(
                &v2000_atom_line("C", 0, 3, 0, 0),
                1,
                SdfReadParams::default(),
            )
            .unwrap(),
            parse_v2000_atom_line(
                &v2000_atom_line("O", 0, 0, 0, 0),
                2,
                SdfReadParams::default(),
            )
            .unwrap(),
        ];
        let mut bonds = Vec::new();
        let mut state = V2000PropertyState {
            first_charge_line: true,
            ..V2000PropertyState::default()
        };

        let file_complete = parse_v2000_property_block(
            &mut reader,
            &mut line_number,
            SdfReadParams::default(),
            &mut atoms,
            &mut bonds,
            &mut state,
        )
        .unwrap();

        assert!(file_complete);
        assert_eq!(atoms[0].spec.formal_charge(), -1);
        assert_eq!(atoms[1].spec.formal_charge(), 0);
        assert_eq!(atoms[1].spec.isotope(), Some(18));
        assert_eq!(atoms[1].spec.radical_electrons(), 1);
        assert_eq!(
            atoms[0].spec.query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(
                vec![6, 7]
            )))
        );
        assert_eq!(atoms[0].spec.prop("_MolFileAtomQuery"), Some("1"));
        assert_eq!(atoms[0].spec.prop("molAttachPoint"), Some("-1"));
        assert_eq!(line_number, 14);
    }

    #[test]
    fn parse_v2000_property_block_applies_zbo_zch_and_hyd_lines() {
        let input = concat!(
            "M  ZCH  1   1  -1\n",
            "M  HYD  1   1   2\n",
            "M  ZBO  1   1   0\n",
            "M  END\n",
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 8;
        let mut atoms = vec![
            parse_v2000_atom_line(
                &v2000_atom_line("C", 0, 0, 0, 0),
                1,
                SdfReadParams::default(),
            )
            .unwrap(),
            parse_v2000_atom_line(
                &v2000_atom_line("O", 0, 0, 0, 0),
                2,
                SdfReadParams::default(),
            )
            .unwrap(),
        ];
        let mut bonds = vec![
            parse_v2000_bond_line(&v2000_bond_line(1, 2, 1, 0, 0), 3, SdfReadParams::default())
                .unwrap(),
        ];
        let mut state = V2000PropertyState::default();

        let file_complete = parse_v2000_property_block(
            &mut reader,
            &mut line_number,
            SdfReadParams::default(),
            &mut atoms,
            &mut bonds,
            &mut state,
        )
        .unwrap();

        assert!(file_complete);
        assert_eq!(atoms[0].spec.formal_charge(), -1);
        assert_eq!(atoms[0].spec.prop("_ZBO_H"), Some("1"));
        assert_eq!(bonds[0].spec.order(), BondOrder::Zero);
    }

    #[test]
    fn parse_v2000_property_block_applies_alias_and_value_lines() {
        let input = concat!("A    1\nAliasLabel\nV    1 payload\nM  END\n");
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 8;
        let mut atoms = vec![
            parse_v2000_atom_line(
                &v2000_atom_line("C", 0, 0, 0, 0),
                1,
                SdfReadParams::default(),
            )
            .unwrap(),
        ];
        let mut bonds = Vec::new();
        let mut state = V2000PropertyState::default();

        let file_complete = parse_v2000_property_block(
            &mut reader,
            &mut line_number,
            SdfReadParams::default(),
            &mut atoms,
            &mut bonds,
            &mut state,
        )
        .unwrap();

        assert!(file_complete);
        assert_eq!(atoms[0].spec.prop("molFileAlias"), Some("AliasLabel"));
        assert_eq!(atoms[0].spec.prop("molFileValue"), Some("payload"));
    }

    #[test]
    fn parse_v2000_property_block_reports_incomplete_at_sdf_delimiter_like_rdkit() {
        let mut reader = std::io::Cursor::new("$$$$\n".as_bytes());
        let mut line_number = 8;
        let mut atoms = vec![
            parse_v2000_atom_line(
                &v2000_atom_line("C", 0, 0, 0, 0),
                1,
                SdfReadParams::default(),
            )
            .unwrap(),
        ];
        let mut bonds = Vec::new();
        let mut state = V2000PropertyState::default();

        let file_complete = parse_v2000_property_block(
            &mut reader,
            &mut line_number,
            SdfReadParams::default(),
            &mut atoms,
            &mut bonds,
            &mut state,
        )
        .unwrap();

        assert!(!file_complete);
        assert_eq!(line_number, 9);
    }

    #[test]
    fn parse_v2000_ctab_preserves_basic_sgroup_membership() {
        let input = format!(
            "{}\n{}\n{}\nM  STY  1   1 SUP\nM  SLB  1   1   7\nM  SAL   1  2   1   2\nM  SPA   1  1   1\nM  SBL   1  1   1\nM  SMT   1 Me\nM  END\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("O", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 4;

        let parsed = parse_v2000_ctab(
            &mut reader,
            MolHeader {
                name: None,
                info: String::new(),
                comments: String::new(),
                ctab_version: CtabVersion::V2000,
            },
            CountsLine {
                atom_count: 2,
                bond_count: 1,
                chiral_flag: 0,
                ctab_version: CtabVersion::V2000,
                v3000_sgroup_count: 0,
                v3000_obj3d_count: 0,
            },
            &mut line_number,
            SdfReadParams::default(),
        )
        .unwrap();

        let sgroup = &parsed.molecule.substance_groups()[0];
        assert_eq!(sgroup.kind(), &SubstanceGroupKind::Superatom);
        assert_eq!(sgroup.rdkit_sequence_id(), Some(1));
        assert_eq!(sgroup.external_id(), Some(7));
        assert_eq!(sgroup.atoms(), &[AtomId::new(0), AtomId::new(1)]);
        assert_eq!(sgroup.parent_atoms(), &[AtomId::new(0)]);
        assert_eq!(sgroup.bonds(), &[BondId::new(0)]);
        assert_eq!(sgroup.label(), Some("Me"));
    }

    #[test]
    fn parse_v2000_ctab_preserves_typed_sgroup_properties() {
        let sdi_line = format!(
            "M  SDI   1  4{:>10.4}{:>10.4}{:>10.4}{:>10.4}",
            0.0, 1.0, 2.0, 3.0
        );
        let sbv_line = format!("M  SBV   1   1{:>10.4}{:>10.4}", 0.5, 0.25);
        let sdt_line = format!(
            "M  SDT   2 {:<30}{:<2}{:<20}{:<2}{}",
            "FIELD", "T", "INFO", "Q", "OP"
        );
        let input = format!(
            "{}\n{}\n{}\nM  STY  2   1 SUP   2 DAT\nM  SST  1   1 ALT\nM  SCN  1   1 HT\nM  SDS EXP  1   1\nM  SAL   1  2   1   2\nM  SBL   1  1   1\n{sdi_line}\n{sbv_line}\n{sdt_line}\nM  SDD   2 display spec\nM  SCD   2 first value\nM  SED   2 second value\nM  SPL  1   2   1\nM  SNC  1   2   5\nM  SAP   1  1   1   2 AP\nM  SCL   2 CLASS\nM  SBT  1   2   1\nM  END\n",
            v2000_atom_line("C", 0, 0, 0, 0),
            v2000_atom_line("O", 0, 0, 0, 0),
            v2000_bond_line(1, 2, 1, 0, 0)
        );
        let mut reader = std::io::Cursor::new(input.as_bytes());
        let mut line_number = 4;

        let parsed = parse_v2000_ctab(
            &mut reader,
            MolHeader {
                name: None,
                info: String::new(),
                comments: String::new(),
                ctab_version: CtabVersion::V2000,
            },
            CountsLine {
                atom_count: 2,
                bond_count: 1,
                chiral_flag: 0,
                ctab_version: CtabVersion::V2000,
                v3000_sgroup_count: 0,
                v3000_obj3d_count: 0,
            },
            &mut line_number,
            SdfReadParams::default(),
        )
        .unwrap();

        let sup = &parsed.molecule.substance_groups()[0];
        assert_eq!(sup.subtype(), Some("ALT"));
        assert_eq!(sup.connection(), Some(&SGroupConnection::HeadToTail));
        assert_eq!(sup.expansion_state(), Some("E"));
        assert_eq!(sup.atoms(), &[AtomId::new(0), AtomId::new(1)]);
        assert_eq!(sup.bonds(), &[BondId::new(0)]);
        assert_eq!(
            sup.bond_role(BondId::new(0)),
            crate::SGroupBondRole::Crossing
        );
        assert_eq!(sup.display().unwrap().brackets[0].p1, [0.0, 1.0]);
        assert_eq!(sup.display().unwrap().brackets[0].p2, [2.0, 3.0]);
        assert_eq!(
            sup.cstates(),
            &[SGroupCState {
                bond: BondId::new(0),
                vector: [0.5, 0.25],
            }]
        );
        assert_eq!(
            sup.attach_points(),
            &[SGroupAttachPoint {
                atom: AtomId::new(0),
                leaving_atom: Some(AtomId::new(1)),
                label: Some("AP".to_string()),
                order: None,
            }]
        );

        let dat = &parsed.molecule.substance_groups()[1];
        assert_eq!(dat.kind(), &SubstanceGroupKind::Data);
        assert_eq!(dat.parent(), Some(sup.id()));
        assert_eq!(dat.component_number(), Some(5));
        assert_eq!(dat.class(), Some("CLASS"));
        assert_eq!(dat.bracket_style(), Some(&SGroupBracketStyle::Parenthesis));
        let data = dat.data().unwrap();
        assert_eq!(data.field_name.as_deref(), Some("FIELD"));
        assert_eq!(data.field_type.as_deref(), Some("T"));
        assert_eq!(data.field_info.as_deref(), Some("INFO"));
        assert_eq!(data.query_type.as_deref(), Some("Q"));
        assert_eq!(data.query_op.as_deref(), Some("OP"));
        assert_eq!(data.field_display.as_deref(), Some("display spec"));
        assert_eq!(data.values, ["first valuesecond value"]);
    }

    #[test]
    fn parse_sgroup_int_field_reports_short_lines_like_rdkit() {
        let mut pos = 6;

        let err = parse_sgroup_int_field("M  STY", 11, &mut pos, true).unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse("SGroup line too short: 'M  STY' on line 11".to_string())
        );
    }

    #[test]
    fn parse_v2000_sgroup_sty_line_rejects_invalid_types_like_rdkit() {
        let mut state = V2000PropertyState::default();

        let err = parse_v2000_sgroup_sty_line(
            "M  STY  1   1 BAD",
            12,
            SdfReadParams::default(),
            &mut state,
        )
        .unwrap_err();

        assert_eq!(
            err,
            SdfReadError::Parse("S group BAD on line 12".to_string())
        );
        assert!(state.sgroups.is_empty());
    }

    #[test]
    fn parse_v2000_atom_line_reads_coordinates_charge_and_atom_map_like_rdkit() {
        let mut line = v2000_atom_line("C", 0, 5, 0, 12);
        line.replace_range(39..42, "  2");
        line.push_str("  1");

        let atom = parse_v2000_atom_line(&line, 9, SdfReadParams::default()).unwrap();

        assert_eq!(atom.line_number, 9);
        assert_eq!(atom.coord_3d, [1.25, -2.5, 0.75]);
        assert_eq!(atom.spec.element().atomic_number(), 6);
        assert_eq!(atom.spec.formal_charge(), -1);
        assert_eq!(atom.spec.atom_map(), Some(12));
        assert_eq!(atom.spec.mol_parity(), Some(2));
        assert_eq!(atom.spec.mol_inversion_flag(), Some(1));
    }

    #[test]
    fn parse_v2000_atom_line_handles_deuterium_and_mass_diff_like_rdkit() {
        let deuterium = parse_v2000_atom_line(
            &v2000_atom_line("D", 0, 0, 0, 0),
            1,
            SdfReadParams::default(),
        )
        .unwrap();
        assert_eq!(deuterium.spec.element().atomic_number(), 1);
        assert_eq!(deuterium.spec.isotope(), Some(2));

        let carbon_13 = parse_v2000_atom_line(
            &v2000_atom_line("C", 1, 0, 0, 0),
            2,
            SdfReadParams::default(),
        )
        .unwrap();
        assert_eq!(carbon_13.spec.isotope(), Some(13));

        let sodium_24 = parse_v2000_atom_line(
            &v2000_atom_line("Na", 1, 0, 0, 0),
            3,
            SdfReadParams::default(),
        )
        .unwrap();
        assert_eq!(sodium_24.spec.isotope(), Some(24));
    }

    #[test]
    fn parse_v2000_atom_line_preserves_query_atom_markers_like_rdkit() {
        let atom = parse_v2000_atom_line(
            &v2000_atom_line("*", 0, 0, 0, 0),
            1,
            SdfReadParams::default(),
        )
        .unwrap();

        assert_eq!(atom.spec.element().atomic_number(), 0);
        assert!(atom.spec.no_implicit());
        assert_eq!(
            atom.spec.query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::Any))
        );
    }

    #[test]
    fn parse_v2000_bond_line_reads_order_stereo_and_query_topology_like_rdkit() {
        let wedge =
            parse_v2000_bond_line(&v2000_bond_line(1, 2, 1, 1, 0), 5, SdfReadParams::default())
                .unwrap();
        assert_eq!(wedge.spec.begin(), AtomId::new(0));
        assert_eq!(wedge.spec.end(), AtomId::new(1));
        assert_eq!(wedge.spec.order(), BondOrder::Single);
        assert_eq!(wedge.spec.direction(), BondDirection::BeginWedge);
        assert_eq!(wedge.molfile_bond_type, 1);
        assert_eq!(wedge.molfile_stereo, Some(1));

        let query =
            parse_v2000_bond_line(&v2000_bond_line(2, 3, 5, 0, 1), 6, SdfReadParams::default())
                .unwrap();
        assert_eq!(
            query.spec.query(),
            Some(&QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                    BondOrder::Single,
                    BondOrder::Double,
                ])),
                QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            ]))
        );
    }
}
