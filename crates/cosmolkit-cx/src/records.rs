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

/// Ordered parser progress, including source-committed checkpoints on failure.
///
/// `consumed` is the parser's actual byte iterator position. It is independent
/// of any diagnostic offset carried by `error`.
#[derive(Debug, Clone, PartialEq)]
pub struct CxParseProgress {
    records: Vec<CxRecord>,
    checkpoints: Vec<CxProgressCheckpoint>,
    consumed: usize,
    complete: bool,
    error: Option<CxParseError>,
}

impl CxParseProgress {
    pub(crate) fn from_parts(
        records: Vec<CxRecord>,
        checkpoints: Vec<CxProgressCheckpoint>,
        consumed: usize,
        complete: bool,
        error: Option<CxParseError>,
    ) -> Self {
        Self {
            records,
            checkpoints,
            consumed,
            complete,
            error,
        }
    }

    /// Parsed record state in source order. On failure, the final record may
    /// contain only the prefix that was scanned before the failure.
    #[must_use]
    pub fn records(&self) -> &[CxRecord] {
        &self.records
    }

    /// Source-ordered mutation checkpoints referencing `records()`.
    #[must_use]
    pub fn checkpoints(&self) -> &[CxProgressCheckpoint] {
        &self.checkpoints
    }

    /// Actual byte position of the CX iterator, not an error diagnostic offset.
    #[must_use]
    pub fn consumed(&self) -> usize {
        self.consumed
    }

    /// Whether the parser consumed a complete CX block.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.complete
    }

    /// Syntax failure, if parsing stopped before the closing pipe.
    #[must_use]
    pub fn error(&self) -> Option<&CxParseError> {
        self.error.as_ref()
    }

    pub(crate) fn into_parts(
        self,
    ) -> (
        Vec<CxRecord>,
        Vec<CxProgressCheckpoint>,
        usize,
        bool,
        Option<CxParseError>,
    ) {
        (
            self.records,
            self.checkpoints,
            self.consumed,
            self.complete,
            self.error,
        )
    }
}

/// A source-ordered checkpoint for a destination effect.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxProgressCheckpoint {
    /// Index into [`CxParseProgress::records`].
    pub record_index: usize,
    /// Item index within a record when the source commits one item at a time.
    pub item_index: Option<usize>,
    /// Actual source iterator byte position at this commit.
    pub cursor: usize,
    /// Source commit point represented by this checkpoint.
    pub phase: CxProgressPhase,
}

/// Commit phase reported by a CX progress checkpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CxProgressPhase {
    /// The source helper began and its initial progress record is available.
    Begin,
    /// The progress record received one item; its destination effect may be
    /// deferred until the helper completes.
    Item,
    /// The source helper completed and its final destination effect is ready.
    Complete,
}

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
