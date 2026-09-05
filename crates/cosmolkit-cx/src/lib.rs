//! Shared CXSMILES extension syntax for COSMolKit parsers.
//!
//! This crate owns representation-independent CX records. SMILES and SMARTS
//! parsers own the lowering of those records into their respective graph
//! states. No molecule, query graph, parser state, or operation-runtime type
//! belongs here.

use std::fmt;

/// A syntax error reported while reading a CX extension block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxParseError {
    /// Byte offset within the supplied CX text, when known.
    pub offset: usize,
    /// Stable source-facing error description.
    pub message: String,
}

impl CxParseError {
    #[must_use]
    pub fn new(offset: usize, message: impl Into<String>) -> Self {
        Self {
            offset,
            message: message.into(),
        }
    }
}

impl fmt::Display for CxParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{} at byte {}", self.message, self.offset)
    }
}

impl std::error::Error for CxParseError {}

/// Parsed CX extension records and the byte position after the closing pipe.
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedCxExtensions {
    records: Vec<CxRecord>,
    consumed: usize,
}

impl ParsedCxExtensions {
    #[must_use]
    pub fn new(records: Vec<CxRecord>, consumed: usize) -> Self {
        Self { records, consumed }
    }

    #[must_use]
    pub fn records(&self) -> &[CxRecord] {
        &self.records
    }

    #[must_use]
    pub fn consumed(&self) -> usize {
        self.consumed
    }

    #[must_use]
    pub fn into_records(self) -> Vec<CxRecord> {
        self.records
    }
}

/// A representation-independent CX extension record.
#[derive(Debug, Clone, PartialEq)]
pub enum CxRecord {
    Coordinates(CxCoordinates),
    AtomLabels(Vec<Option<String>>),
    AtomValues(Vec<Option<String>>),
    AtomProperties(Vec<CxAtomProperty>),
    CoordinateBonds(CxCoordinateBonds),
    ZeroBonds(Vec<usize>),
    EnhancedStereo(CxEnhancedStereo),
    Unsaturation(Vec<usize>),
    RingBonds(Vec<CxRingBond>),
    LinkNodes(Vec<CxLinkNode>),
    DataSGroup(CxDataSGroup),
    SGroupHierarchy(Vec<CxSGroupHierarchy>),
    PolymerSGroup(CxPolymerSGroup),
    Substitution(Vec<CxAtomConstraint>),
    VariableAttachments(Vec<CxVariableAttachment>),
    WedgedBonds(Vec<CxWedgeBond>),
    DoubleBondStereo(CxDoubleBondStereo),
    Radicals(Vec<CxRadical>),
    Unknown(String),
}

/// CX coordinate conformer values. Empty entries represent omitted points.
#[derive(Debug, Clone, PartialEq)]
pub struct CxCoordinates {
    pub conformer: usize,
    pub values: Vec<Option<[f64; 3]>>,
    pub is_3d: bool,
}

/// An atom property assignment from an `atomProp:` record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxAtomProperty {
    pub atom: usize,
    pub name: String,
    pub value: String,
}

/// A CX coordinate or hydrogen bond annotation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxCoordinateBonds {
    pub kind: CxCoordinateBondKind,
    pub bonds: Vec<CxBondReference>,
}

/// Bond kind encoded by the CX `C:` and `H:` records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CxCoordinateBondKind {
    Dative,
    Hydrogen,
}

/// A bond reference uses the source atom and CX/SMILES bond indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxBondReference {
    pub atom: usize,
    pub bond: usize,
}

/// Enhanced stereo group assignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxEnhancedStereo {
    pub kind: CxStereoGroupKind,
    pub group_id: u32,
    pub atoms: Vec<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CxStereoGroupKind {
    Absolute,
    Or,
    And,
}

/// Ring-bond count constraint from `rb:`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxRingBond {
    pub atom: usize,
    pub constraint: CxCountConstraint,
}

/// Substitution count constraint from `s:`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxAtomConstraint {
    pub atom: usize,
    pub constraint: CxCountConstraint,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CxCountConstraint {
    Exact(u32),
    LessEqual(u32),
    QueryScan,
}

/// A CX wedge/dash bond annotation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxWedgeBond {
    pub atom: usize,
    pub bond: usize,
    pub direction: CxWedgeDirection,
    pub configuration: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CxWedgeDirection {
    Unknown,
    BeginWedge,
    BeginDash,
}

/// A cis/trans/unknown double-bond stereo assignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxDoubleBondStereo {
    pub stereo: CxDoubleBondStereoKind,
    pub bonds: Vec<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CxDoubleBondStereoKind {
    Any,
    Cis,
    Trans,
}

/// A radical assignment from a `^n:` section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxRadical {
    pub atom: usize,
    pub electrons: u8,
}

/// One link-node declaration from an `LN:` record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxLinkNode {
    pub atom: usize,
    pub start_repetitions: usize,
    pub end_repetitions: usize,
    /// Explicit outer atoms. When absent, the destination must obtain the two
    /// neighbours of `atom` after validating that its degree is exactly two.
    pub outer_atoms: Option<[usize; 2]>,
}

/// A data substance-group syntax record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxDataSGroup {
    pub atoms: Vec<usize>,
    pub field_name: String,
    pub data: String,
    pub query_op: String,
    pub field_info: String,
    pub field_tag: String,
    pub coordinates: Option<String>,
}

/// One parent-to-children relationship from an `SgH:` record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxSGroupHierarchy {
    pub parent: usize,
    pub children: Vec<usize>,
}

/// A polymer substance-group syntax record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxPolymerSGroup {
    pub type_code: String,
    pub atoms: Vec<usize>,
    pub label: String,
    pub connect: String,
    pub head_crossings: Vec<usize>,
    pub tail_crossings: Vec<usize>,
}

/// One variable-attachment declaration from an `m:` record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CxVariableAttachment {
    pub atom: usize,
    pub endpoints: Vec<usize>,
}

/// Parse one CX extension block without referring to a destination graph.
///
/// RDKit✔️✔️: `parseCXExtensions` validates the outer pipes and dispatches
/// records in source order. Destination-specific atom and bond validity is
/// deliberately deferred to the lowering crate.
pub fn parse_cx_extensions(text: &str) -> Result<ParsedCxExtensions, CxParseError> {
    if text.is_empty() {
        return Ok(ParsedCxExtensions::new(Vec::new(), 0));
    }
    let bytes = text.as_bytes();
    if bytes.first().copied() != Some(b'|') {
        return Err(CxParseError::new(
            0,
            "CXSMILES extension does not start with |",
        ));
    }

    let mut cursor = 1;
    let mut conformer = 0;
    let mut records = Vec::new();
    while cursor < bytes.len() && bytes[cursor] != b'|' {
        if bytes[cursor] == b',' {
            cursor += 1;
            continue;
        }
        let start = cursor;
        let record = if bytes[cursor] == b'(' {
            let mut coordinates = parse_coordinates(text, &mut cursor)?;
            coordinates.conformer = conformer;
            conformer += 1;
            CxRecord::Coordinates(coordinates)
        } else if bytes[cursor] == b'$' {
            parse_labels_or_values(text, &mut cursor)?
        } else if text[cursor..].starts_with("atomProp:") {
            CxRecord::AtomProperties(parse_atom_properties(text, &mut cursor)?)
        } else if matches!(bytes[cursor], b'C' | b'H')
            && cursor + 1 < bytes.len()
            && bytes[cursor + 1] == b':'
        {
            CxRecord::CoordinateBonds(parse_coordinate_bonds(text, &mut cursor)?)
        } else if bytes[cursor] == b'Z' {
            CxRecord::ZeroBonds(parse_index_list(text, &mut cursor, b'Z')?)
        } else if bytes[cursor] == b'^' {
            CxRecord::Radicals(parse_radicals(text, &mut cursor)?)
        } else if matches!(bytes[cursor], b'a' | b'o' | b'&') {
            CxRecord::EnhancedStereo(parse_enhanced_stereo(text, &mut cursor)?)
        } else if text[cursor..].starts_with("rb:") {
            CxRecord::RingBonds(parse_ring_bonds(text, &mut cursor)?)
        } else if bytes[cursor] == b'u' {
            CxRecord::Unsaturation(parse_index_list(text, &mut cursor, b'u')?)
        } else if bytes[cursor] == b's' {
            CxRecord::Substitution(parse_count_constraints(text, &mut cursor, b's')?)
        } else if bytes[cursor] == b'm' {
            CxRecord::VariableAttachments(parse_variable_attachments(text, &mut cursor)?)
        } else if bytes[cursor] == b'w' {
            CxRecord::WedgedBonds(parse_wedge_bonds(text, &mut cursor)?)
        } else if text[cursor..].starts_with("ctu") {
            CxRecord::DoubleBondStereo(parse_double_bond_stereo(
                text,
                &mut cursor,
                CxDoubleBondStereoKind::Any,
            )?)
        } else if bytes[cursor] == b'c' {
            CxRecord::DoubleBondStereo(parse_double_bond_stereo(
                text,
                &mut cursor,
                CxDoubleBondStereoKind::Cis,
            )?)
        } else if bytes[cursor] == b't' {
            CxRecord::DoubleBondStereo(parse_double_bond_stereo(
                text,
                &mut cursor,
                CxDoubleBondStereoKind::Trans,
            )?)
        } else if text[cursor..].starts_with("LN:") {
            CxRecord::LinkNodes(parse_link_nodes(text, &mut cursor)?)
        } else if text[cursor..].starts_with("SgD:") {
            CxRecord::DataSGroup(parse_data_sgroup(text, &mut cursor)?)
        } else if text[cursor..].starts_with("SgH:") {
            CxRecord::SGroupHierarchy(parse_sgroup_hierarchy(text, &mut cursor)?)
        } else if text[cursor..].starts_with("Sg:") {
            CxRecord::PolymerSGroup(parse_polymer_sgroup(text, &mut cursor)?)
        } else {
            let raw = consume_until_record_boundary(text, &mut cursor);
            CxRecord::Unknown(raw)
        };
        if cursor == start {
            return Err(CxParseError::new(cursor, "CX parser made no progress"));
        }
        records.push(record);
    }
    if cursor >= bytes.len() || bytes[cursor] != b'|' {
        return Err(CxParseError::new(
            cursor,
            "failure parsing CXSMILES extensions",
        ));
    }
    Ok(ParsedCxExtensions::new(records, cursor + 1))
}

fn consume_until_record_boundary(text: &str, cursor: &mut usize) -> String {
    let start = *cursor;
    while *cursor < text.len() && text.as_bytes()[*cursor] != b'|' {
        if text.as_bytes()[*cursor] == b',' {
            let next = *cursor + 1;
            if next < text.len() && is_record_start(&text[next..]) {
                break;
            }
        }
        *cursor += 1;
    }
    text[start..*cursor].to_owned()
}

fn is_record_start(text: &str) -> bool {
    text.starts_with('(')
        || text.starts_with('$')
        || text.starts_with("atomProp:")
        || text.starts_with("C:")
        || text.starts_with("H:")
        || text.starts_with("Z:")
        || text.starts_with('^')
        || text.starts_with('a')
        || text.starts_with('o')
        || text.starts_with('&')
        || text.starts_with("rb:")
        || text.starts_with("LN:")
        || text.starts_with("SgD:")
        || text.starts_with("SgH:")
        || text.starts_with("Sg:")
        || text.starts_with("u:")
        || text.starts_with("s:")
        || text.starts_with("m:")
        || text.starts_with("w:")
        || text.starts_with("ctu:")
        || text.starts_with("c:")
        || text.starts_with("t:")
}

fn parse_coordinates(text: &str, cursor: &mut usize) -> Result<CxCoordinates, CxParseError> {
    let start = *cursor;
    *cursor += 1;
    let mut values = Vec::new();
    let mut is_3d = false;
    while *cursor < text.len() && text.as_bytes()[*cursor] != b')' {
        let field_start = *cursor;
        while *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b';' | b')') {
            *cursor += 1;
        }
        let field = &text[field_start..*cursor];
        if field.is_empty() {
            values.push(None);
        } else {
            let mut parts = field.split(',');
            let x = parse_coordinate_component(parts.next(), field_start)?;
            let y = parse_coordinate_component(parts.next(), field_start)?;
            let z = parse_coordinate_component(parts.next(), field_start)?;
            if x.is_none() && y.is_none() {
                return Err(CxParseError::new(field_start, "invalid CX coordinate"));
            }
            if let Some(z) = z {
                is_3d |= z.abs() > 1e-3;
            }
            values.push(Some([
                x.unwrap_or_default(),
                y.unwrap_or_default(),
                z.unwrap_or_default(),
            ]));
        }
        if *cursor < text.len() && text.as_bytes()[*cursor] == b';' {
            *cursor += 1;
        }
    }
    if *cursor >= text.len() || text.as_bytes()[*cursor] != b')' {
        return Err(CxParseError::new(
            start,
            "unterminated CX coordinate record",
        ));
    }
    *cursor += 1;
    Ok(CxCoordinates {
        conformer: 0,
        values,
        is_3d,
    })
}

fn parse_coordinate_component(
    value: Option<&str>,
    offset: usize,
) -> Result<Option<f64>, CxParseError> {
    let Some(value) = value else {
        return Ok(None);
    };
    if value.is_empty() {
        return Ok(None);
    }
    value
        .parse::<f64>()
        .map(Some)
        .map_err(|_| CxParseError::new(offset, "invalid CX coordinate"))
}

fn parse_labels_or_values(text: &str, cursor: &mut usize) -> Result<CxRecord, CxParseError> {
    let value = text[*cursor..].starts_with("$_AV:");
    *cursor += if value { 5 } else { 1 };
    let mut fields = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor] != b'$' {
        let start = *cursor;
        while *cursor < text.len() && !matches!(text.as_bytes()[*cursor], b';' | b'$') {
            *cursor += 1;
        }
        fields.push(if *cursor == start {
            None
        } else {
            Some(text[start..*cursor].to_owned())
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b';' {
            *cursor += 1;
        }
    }
    if *cursor >= text.len() || text.as_bytes()[*cursor] != b'$' {
        return Err(CxParseError::new(*cursor, "unterminated CX atom record"));
    }
    *cursor += 1;
    Ok(if value {
        CxRecord::AtomValues(fields)
    } else {
        CxRecord::AtomLabels(fields)
    })
}

fn parse_atom_properties(
    text: &str,
    cursor: &mut usize,
) -> Result<Vec<CxAtomProperty>, CxParseError> {
    *cursor += "atomProp:".len();
    let start = *cursor;
    let raw = consume_until_record_boundary(text, cursor);
    let mut result = Vec::new();
    for entry in raw.split(':') {
        let mut fields = entry.splitn(3, '.');
        let Some(atom) = fields.next().and_then(|v| v.parse::<usize>().ok()) else {
            return Err(CxParseError::new(start, "invalid CX atomProp atom index"));
        };
        let Some(name) = fields.next() else {
            return Err(CxParseError::new(start, "invalid CX atomProp name"));
        };
        let Some(value) = fields.next() else {
            return Err(CxParseError::new(start, "invalid CX atomProp value"));
        };
        result.push(CxAtomProperty {
            atom,
            name: name.to_owned(),
            value: value.to_owned(),
        });
    }
    Ok(result)
}

fn parse_coordinate_bonds(
    text: &str,
    cursor: &mut usize,
) -> Result<CxCoordinateBonds, CxParseError> {
    let kind = if text.as_bytes()[*cursor] == b'C' {
        CxCoordinateBondKind::Dative
    } else {
        CxCoordinateBondKind::Hydrogen
    };
    *cursor += 2;
    let mut bonds = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b'.')?;
        let bond = read_number(text, cursor)?;
        bonds.push(CxBondReference { atom, bond });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        }
    }
    Ok(CxCoordinateBonds { kind, bonds })
}

fn parse_index_list(
    text: &str,
    cursor: &mut usize,
    marker: u8,
) -> Result<Vec<usize>, CxParseError> {
    expect_byte(text, cursor, marker)?;
    expect_byte(text, cursor, b':')?;
    let mut values = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        values.push(read_number(text, cursor)?);
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(values)
}

fn parse_radicals(text: &str, cursor: &mut usize) -> Result<Vec<CxRadical>, CxParseError> {
    let mut result = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor] == b'^' {
        *cursor += 1;
        let electrons = match text.as_bytes().get(*cursor).copied() {
            Some(b'1') => 1,
            Some(b'2'..=b'4') => 2,
            Some(b'5'..=b'7') => 3,
            _ => return Err(CxParseError::new(*cursor, "invalid CX radical marker")),
        };
        *cursor += 1;
        expect_byte(text, cursor, b':')?;
        loop {
            let atom = read_number(text, cursor)?;
            result.push(CxRadical { atom, electrons });
            if *cursor >= text.len() || text.as_bytes()[*cursor] != b',' {
                break;
            }
            *cursor += 1;
            if !text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
                break;
            }
        }
    }
    Ok(result)
}

fn parse_enhanced_stereo(text: &str, cursor: &mut usize) -> Result<CxEnhancedStereo, CxParseError> {
    let kind = match text.as_bytes().get(*cursor).copied() {
        Some(b'a') => CxStereoGroupKind::Absolute,
        Some(b'o') => CxStereoGroupKind::Or,
        Some(b'&') => CxStereoGroupKind::And,
        _ => return Err(CxParseError::new(*cursor, "invalid CX stereo group")),
    };
    *cursor += 1;
    let group_id = if kind == CxStereoGroupKind::Absolute {
        0
    } else if text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
        read_number(text, cursor)? as u32
    } else {
        0
    };
    expect_byte(text, cursor, b':')?;
    let atoms = parse_number_list(text, cursor);
    Ok(CxEnhancedStereo {
        kind,
        group_id,
        atoms,
    })
}

fn parse_ring_bonds(text: &str, cursor: &mut usize) -> Result<Vec<CxRingBond>, CxParseError> {
    *cursor += 3;
    let mut result = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let constraint = if text.as_bytes().get(*cursor) == Some(&b'*') {
            *cursor += 1;
            CxCountConstraint::QueryScan
        } else {
            let value = read_number(text, cursor)? as u32;
            match value {
                0 | 2 | 3 => CxCountConstraint::Exact(value),
                4 => CxCountConstraint::LessEqual(value),
                _ => {
                    return Err(CxParseError::new(
                        *cursor,
                        "unrecognized CX ring-bond count",
                    ));
                }
            }
        };
        result.push(CxRingBond { atom, constraint });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_count_constraints(
    text: &str,
    cursor: &mut usize,
    marker: u8,
) -> Result<Vec<CxAtomConstraint>, CxParseError> {
    expect_byte(text, cursor, marker)?;
    expect_byte(text, cursor, b':')?;
    let mut result = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let constraint = if text.as_bytes().get(*cursor) == Some(&b'*') {
            *cursor += 1;
            CxCountConstraint::QueryScan
        } else {
            CxCountConstraint::Exact(read_number(text, cursor)? as u32)
        };
        result.push(CxAtomConstraint { atom, constraint });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_link_nodes(text: &str, cursor: &mut usize) -> Result<Vec<CxLinkNode>, CxParseError> {
    // RDKit✔️✔️: `parse_linknodes` reads
    // `atom:min.max[.outer1.outer2]` records separated by commas. Resolution
    // of omitted outer atoms is graph-dependent and remains in the lowerer.
    *cursor += "LN:".len();
    let mut records = Vec::new();
    while text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let start_repetitions = read_number(text, cursor)?;
        expect_byte(text, cursor, b'.')?;
        let end_repetitions = read_number(text, cursor)?;
        let outer_atoms = if text.as_bytes().get(*cursor) == Some(&b'.') {
            *cursor += 1;
            let first = read_number(text, cursor)?;
            expect_byte(text, cursor, b'.')?;
            Some([first, read_number(text, cursor)?])
        } else {
            None
        };
        records.push(CxLinkNode {
            atom,
            start_repetitions,
            end_repetitions,
            outer_atoms,
        });
        if text.as_bytes().get(*cursor) == Some(&b',')
            && text
                .as_bytes()
                .get(*cursor + 1)
                .is_some_and(u8::is_ascii_digit)
        {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(records)
}

fn read_until<'a>(text: &'a str, cursor: &mut usize, delimiters: &[u8]) -> &'a str {
    let start = *cursor;
    while *cursor < text.len() && !delimiters.contains(&text.as_bytes()[*cursor]) {
        *cursor += 1;
    }
    &text[start..*cursor]
}

fn read_colon_field(text: &str, cursor: &mut usize) -> String {
    let value = read_until(text, cursor, b":,|").to_owned();
    if text.as_bytes().get(*cursor) == Some(&b':') {
        *cursor += 1;
    }
    value
}

fn parse_data_sgroup(text: &str, cursor: &mut usize) -> Result<CxDataSGroup, CxParseError> {
    // RDKit✔️✔️: `parse_data_sgroup` reads atom indices followed by FIELDNAME,
    // DATAFIELDS, QUERYOP, FIELDINFO, FIELDTAG, and optional coordinates.
    *cursor += "SgD:".len();
    let atoms = parse_delimited_number_list(text, cursor, b',');
    expect_byte(text, cursor, b':')?;
    let field_name = read_colon_field(text, cursor);
    let data = read_colon_field(text, cursor);
    let query_op = read_colon_field(text, cursor);
    let field_info = read_colon_field(text, cursor);
    let field_tag = read_colon_field(text, cursor);
    let coordinates = if text.as_bytes().get(*cursor) == Some(&b'(') {
        let start = *cursor;
        while *cursor < text.len() && text.as_bytes()[*cursor] != b')' {
            *cursor += 1;
        }
        if *cursor >= text.len() {
            return Err(CxParseError::new(
                start,
                "unterminated CX SGroup coordinates",
            ));
        }
        *cursor += 1;
        Some(text[start..*cursor].to_owned())
    } else {
        None
    };
    Ok(CxDataSGroup {
        atoms,
        field_name,
        data,
        query_op,
        field_info,
        field_tag,
        coordinates,
    })
}

fn parse_sgroup_hierarchy(
    text: &str,
    cursor: &mut usize,
) -> Result<Vec<CxSGroupHierarchy>, CxParseError> {
    // RDKit✔️✔️: `parse_sgroup_hierarchy` reads
    // `parent:child.child` relationships separated by commas.
    *cursor += "SgH:".len();
    let mut result = Vec::new();
    loop {
        let parent = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let children = parse_delimited_number_list(text, cursor, b'.');
        result.push(CxSGroupHierarchy { parent, children });
        if text.as_bytes().get(*cursor) == Some(&b',')
            && text
                .as_bytes()
                .get(*cursor + 1)
                .is_some_and(u8::is_ascii_digit)
        {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_polymer_sgroup(text: &str, cursor: &mut usize) -> Result<CxPolymerSGroup, CxParseError> {
    // RDKit✔️✔️: `parse_polymer_sgroup` reads type, atoms, label, connect,
    // head-crossing bonds, and tail-crossing bonds in this order.
    *cursor += "Sg:".len();
    let type_code = read_until(text, cursor, b":").to_owned();
    expect_byte(text, cursor, b':')?;
    let atoms = parse_delimited_number_list(text, cursor, b',');
    let mut label = String::new();
    let mut connect = String::new();
    let mut head_crossings = Vec::new();
    let mut tail_crossings = Vec::new();
    if text.as_bytes().get(*cursor) == Some(&b':') {
        *cursor += 1;
        label = read_until(text, cursor, b":|").to_owned();
        if text.as_bytes().get(*cursor) == Some(&b':') {
            *cursor += 1;
            connect = read_until(text, cursor, b":|,").to_owned();
            if text.as_bytes().get(*cursor) == Some(&b':') {
                *cursor += 1;
                head_crossings = parse_delimited_number_list(text, cursor, b',');
                if text.as_bytes().get(*cursor) == Some(&b':') {
                    *cursor += 1;
                    tail_crossings = parse_delimited_number_list(text, cursor, b',');
                }
            }
        }
    }
    Ok(CxPolymerSGroup {
        type_code,
        atoms,
        label,
        connect,
        head_crossings,
        tail_crossings,
    })
}

fn parse_variable_attachments(
    text: &str,
    cursor: &mut usize,
) -> Result<Vec<CxVariableAttachment>, CxParseError> {
    // RDKit✔️✔️: `parse_variable_attachments` reads
    // `attachmentAtom:endpoint.endpoint` records separated by commas.
    *cursor += "m:".len();
    let mut result = Vec::new();
    loop {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b':')?;
        let endpoints = parse_delimited_number_list(text, cursor, b'.');
        result.push(CxVariableAttachment { atom, endpoints });
        if text.as_bytes().get(*cursor) == Some(&b',')
            && text
                .as_bytes()
                .get(*cursor + 1)
                .is_some_and(u8::is_ascii_digit)
        {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_delimited_number_list(text: &str, cursor: &mut usize, delimiter: u8) -> Vec<usize> {
    let mut values = Vec::new();
    while text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
        let Ok(value) = read_number(text, cursor) else {
            break;
        };
        values.push(value);
        if text.as_bytes().get(*cursor) == Some(&delimiter)
            && text
                .as_bytes()
                .get(*cursor + 1)
                .is_some_and(u8::is_ascii_digit)
        {
            *cursor += 1;
        } else {
            break;
        }
    }
    values
}

fn parse_wedge_bonds(text: &str, cursor: &mut usize) -> Result<Vec<CxWedgeBond>, CxParseError> {
    expect_byte(text, cursor, b'w')?;
    let (direction, configuration) = match text.as_bytes().get(*cursor).copied() {
        Some(b':') => (CxWedgeDirection::Unknown, 2),
        Some(b'U') => {
            *cursor += 1;
            (CxWedgeDirection::BeginWedge, 1)
        }
        Some(b'D') => {
            *cursor += 1;
            (CxWedgeDirection::BeginDash, 3)
        }
        _ => return Err(CxParseError::new(*cursor, "invalid CX wedge marker")),
    };
    expect_byte(text, cursor, b':')?;
    let mut result = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        let atom = read_number(text, cursor)?;
        expect_byte(text, cursor, b'.')?;
        let bond = read_number(text, cursor)?;
        result.push(CxWedgeBond {
            atom,
            bond,
            direction,
            configuration,
        });
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    Ok(result)
}

fn parse_double_bond_stereo(
    text: &str,
    cursor: &mut usize,
    stereo: CxDoubleBondStereoKind,
) -> Result<CxDoubleBondStereo, CxParseError> {
    while *cursor < text.len() && text.as_bytes()[*cursor] != b':' {
        *cursor += 1;
    }
    expect_byte(text, cursor, b':')?;
    Ok(CxDoubleBondStereo {
        stereo,
        bonds: parse_number_list(text, cursor),
    })
}

fn parse_number_list(text: &str, cursor: &mut usize) -> Vec<usize> {
    let mut values = Vec::new();
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        if let Ok(value) = read_number(text, cursor) {
            values.push(value);
        } else {
            break;
        }
        if *cursor < text.len() && text.as_bytes()[*cursor] == b',' {
            *cursor += 1;
        } else {
            break;
        }
    }
    values
}

fn read_number(text: &str, cursor: &mut usize) -> Result<usize, CxParseError> {
    let start = *cursor;
    while *cursor < text.len() && text.as_bytes()[*cursor].is_ascii_digit() {
        *cursor += 1;
    }
    text[start..*cursor]
        .parse()
        .map_err(|_| CxParseError::new(start, "invalid CX integer"))
}

fn expect_byte(text: &str, cursor: &mut usize, expected: u8) -> Result<(), CxParseError> {
    if text.as_bytes().get(*cursor) == Some(&expected) {
        *cursor += 1;
        Ok(())
    } else {
        Err(CxParseError::new(
            *cursor,
            format!("expected '{}', found CX syntax mismatch", expected as char),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::{CxRecord, CxStereoGroupKind, parse_cx_extensions};

    #[test]
    fn parses_representation_independent_records() {
        let parsed = parse_cx_extensions("|(0,0,;1,0,)$C;O$| name").unwrap();
        assert_eq!(parsed.consumed(), 18);
        assert_eq!(parsed.records().len(), 2);
        assert!(matches!(parsed.records()[0], CxRecord::Coordinates(_)));
        assert!(matches!(parsed.records()[1], CxRecord::AtomLabels(_)));
    }

    #[test]
    fn parses_query_and_stereo_records_without_a_graph() {
        let parsed = parse_cx_extensions("|u:0,rb:1:4,o2:0,1|").unwrap();
        assert!(matches!(parsed.records()[0], CxRecord::Unsaturation(_)));
        assert!(matches!(parsed.records()[1], CxRecord::RingBonds(_)));
        match &parsed.records()[2] {
            CxRecord::EnhancedStereo(stereo) => {
                assert_eq!(stereo.kind, CxStereoGroupKind::Or);
                assert_eq!(stereo.group_id, 2);
            }
            other => panic!("unexpected record: {other:?}"),
        }
    }

    #[test]
    fn rejects_missing_outer_pipe() {
        let error = parse_cx_extensions("u:0|").unwrap_err();
        assert_eq!(error.offset, 0);
    }

    #[test]
    fn parses_structural_records_without_destination_state() {
        let parsed = parse_cx_extensions(
            "|LN:1:1.3.2.6,SgD:2,1:FIELD:info::::,SgH:1:0,Sg:n:0,1:n:ht:2:3,m:2:3.5.4|",
        )
        .unwrap();
        assert!(matches!(parsed.records()[0], CxRecord::LinkNodes(_)));
        assert!(matches!(parsed.records()[1], CxRecord::DataSGroup(_)));
        assert!(matches!(parsed.records()[2], CxRecord::SGroupHierarchy(_)));
        assert!(matches!(parsed.records()[3], CxRecord::PolymerSGroup(_)));
        assert!(matches!(
            parsed.records()[4],
            CxRecord::VariableAttachments(_)
        ));
        match &parsed.records()[1] {
            CxRecord::DataSGroup(group) => {
                assert_eq!(group.atoms, vec![2, 1]);
                assert_eq!(group.field_name, "FIELD");
                assert_eq!(group.data, "info");
            }
            other => panic!("unexpected record: {other:?}"),
        }
    }
}
