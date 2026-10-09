//! COSMolKit native molecule archive, adapted from the complete original
//! properties/mol_pickler.rs, SHA256 a72f13b9b36c6d691b23e2f5eecabb6bff71d166e66ab6f638972e36abcdc6fb.
//! Writers emit archive 2.0 with one complete Müsli molecule block and one
//! derived-state block. Raw1..4 and archive1.0..1.3 remain legacy read formats.
//! This is COS-native serialization, not the RDKit binary wire protocol.

use std::collections::{BTreeMap, BTreeSet};
mod native_state_v2;
use cosmolkit_model::PropertyText;

use serde::{Deserialize, Serialize};

use cosmolkit_model::{
    Atom, AtomId, AtomPdbResidueInfo, AtomSpec, Bond, BondDirection, BondId, BondOrder, BondSpec,
    BondStereo, ChiralTag, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, Element,
    Hybridization, MoleculeProperties, PropertyValue, SGroupAttachPoint, SGroupBondRole,
    SGroupBracket, SGroupBracketStyle, SGroupCState, SGroupConnection, SGroupData, SGroupDisplay,
    SdfPropertyList, SdfPropertyListTarget, StereoGroup, StereoGroupKind, SubstanceGroup,
    SubstanceGroupId, SubstanceGroupKind, TemplateAttachment, TemplateAttachmentOrder,
    TopologyBlock, ordered_atom_properties, ordered_bond_properties,
    replace_atom_template_attachment_order,
};

use cosmolkit_core::{RingFindType, RingInfo, ValenceAssignment, property_value_to_string};

mod archive_v2;

// ──────────────────────────────────────────────
// Format version
// ──────────────────────────────────────────────
const PICKLE_VERSION: u8 = 4;
const ARCHIVE_MAGIC: &[u8; 8] = b"CSMOLPKL";
const ARCHIVE_MAJOR: u16 = 1;
const ARCHIVE_MINOR: u16 = 3;
const SECTION_FLAG_REQUIRED: u8 = 1;
const SECTION_CODEC_RAW: u8 = 0;
const SECTION_CODEC_POSTCARD: u8 = 1;
const SECTION_MANIFEST: u16 = 1;
const SECTION_MOLECULE_STATE: u16 = 2;
const SECTION_DERIVED_STATE: u16 = 3;
const SECTION_CANONICAL_STATE: u16 = 4;
const CANONICAL_STATE_VERSION: u16 = 3;
const MANIFEST_VERSION: u16 = 1;
const MOLECULE_STATE_VERSION: u16 = 1;
const DERIVED_STATE_VERSION: u16 = 1;
const MAX_DERIVED_ROWS: usize = 1_000_000;

// ──────────────────────────────────────────────
// Error type
// ──────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PickleError {
    #[error("{0}")]
    StereoGroup(#[from] cosmolkit_model::StereoGroupError),

    #[error("unexpected end of data while reading pickle")]
    UnexpectedEof,
    #[error("unsupported pickle version: {0}")]
    UnsupportedVersion(u8),
    #[error("unsupported archive version: {major}.{minor}")]
    UnsupportedArchiveVersion { major: u16, minor: u16 },
    #[error("unsupported archive section {section} version: {version}")]
    UnsupportedSectionVersion { section: u16, version: u16 },
    #[error("invalid archive: {0}")]
    InvalidArchive(String),
    #[error("missing required archive section: {0}")]
    MissingRequiredSection(u16),
    #[error("duplicate archive section: {0}")]
    DuplicateSection(u16),
    #[error("unknown required archive section: {0}")]
    UnknownRequiredSection(u16),
    #[error("data length mismatch: expected {expected}, got {actual}")]
    DataLengthMismatch { expected: usize, actual: usize },
    #[error("invalid enum value: {value} for {type_name}")]
    InvalidEnumValue { value: u8, type_name: &'static str },
    #[error("invalid molecule state after unpickling: {0}")]
    InvalidMolecule(String),
    #[error("too many atoms: {0}")]
    TooManyAtoms(usize),
    #[error("too many bonds: {0}")]
    TooManyBonds(usize),
    #[error("string too long: {0}")]
    StringTooLong(usize),
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct ArchiveManifestV1 {
    crate_version: String,
    molecule_state_codec: u8,
    molecule_state_version: u16,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct MoleculeStateV1 {
    encoding: u8,
    encoding_version: u8,
    payload: Vec<u8>,
}

// ──────────────────────────────────────────────
// Binary writer + reader helpers
// ──────────────────────────────────────────────

struct PickleWriter {
    buf: Vec<u8>,
    error: Option<PickleError>,
    legacy_current: Option<LegacyStoreProvenance>,
}

impl PickleWriter {
    fn new() -> Self {
        Self {
            buf: Vec::new(),
            error: None,
            legacy_current: None,
        }
    }

    fn into_inner(self) -> Result<Vec<u8>, PickleError> {
        match self.error {
            Some(error) => Err(error),
            None => Ok(self.buf),
        }
    }

    fn write_u8(&mut self, v: u8) {
        self.buf.push(v);
    }

    fn write_u32(&mut self, v: u32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    fn write_i32(&mut self, v: i32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    fn write_i8(&mut self, v: i8) {
        self.buf.push(v as u8);
    }

    fn write_f64(&mut self, v: f64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    fn write_bool(&mut self, v: bool) {
        self.buf.push(if v { 1 } else { 0 });
    }

    fn write_string(&mut self, value: impl AsRef<[u8]>) {
        let value = value.as_ref();
        // Old raw1..3/canonical1..2 UTF8 remains a checked old wire rule.
        if std::str::from_utf8(value).is_err() {
            self.error.get_or_insert(PickleError::InvalidMolecule(
                "invalid UTF-8 in pickle string".into(),
            ));
            return;
        }
        match u32::try_from(value.len()) {
            Ok(count) if value.len() <= 10_000_000 => {
                self.write_u32(count);
                self.buf.extend_from_slice(value);
            }
            _ => {
                self.error
                    .get_or_insert(PickleError::StringTooLong(value.len()));
            }
        }
    }

    fn write_count(&mut self, count: usize) {
        match u32::try_from(count) {
            Ok(value) => self.write_u32(value),
            Err(_) => {
                self.error
                    .get_or_insert(PickleError::InvalidArchive(format!(
                        "count or index does not fit legacy u32: {count}"
                    )));
            }
        }
    }

    fn write_u64(&mut self, value: u64) {
        self.buf.extend_from_slice(&value.to_le_bytes());
    }

    fn write_typed_props(&mut self, props: &BTreeMap<PropertyText, PropertyValue>) {
        self.write_props(props);
    }

    fn write_option_string(&mut self, s: Option<&PropertyText>) {
        match s {
            Some(val) => {
                self.write_bool(true);
                self.write_string(val);
            }
            None => {
                self.write_bool(false);
            }
        }
    }

    fn write_props<K: AsRef<[u8]>, V: LegacyTextValue>(&mut self, props: &BTreeMap<K, V>) {
        let omit = self.legacy_current.as_ref().is_some_and(|p| !p.reserved);
        let include = |key: &K, value: &V| {
            !(omit && key.as_ref() == b"__computedProps" && value.synthetic_reserved())
        };
        self.write_count(
            props
                .iter()
                .filter(|(key, value)| include(key, value))
                .count(),
        );
        for (key, value) in props {
            if include(key, value) {
                self.write_string(key);
                // Legacy companion verification must reconstruct the original
                // reserved value, even when import normalizes its conflict.
                let original = self.legacy_current.as_ref().and_then(|provenance| {
                    (key.as_ref() == b"__computedProps")
                        .then_some(provenance.collision_reserved.as_ref())
                        .flatten()
                });
                let text = match original {
                    Some(original) => original.legacy_text(),
                    None => value.legacy_text(),
                };
                match text {
                    Ok(v) => self.write_string(v),
                    Err(e) => {
                        self.error.get_or_insert(e);
                    }
                }
            }
        }
    }
    fn write_computed_props(
        &mut self,
        props: Result<Option<&[PropertyText]>, cosmolkit_model::PropertyValueError>,
    ) {
        if let Some(provenance) = self.legacy_current.clone() {
            // Actual old reserved records use their original independent flags
            // for companion verification, including normalized collisions.
            if !provenance.reserved {
                let names = match props {
                    Ok(names) => names.unwrap_or_default(),
                    Err(e) => {
                        self.error.get_or_insert(invalid(e));
                        return;
                    }
                };
                if names.iter().collect::<BTreeSet<_>>()
                    != provenance.computed.iter().collect::<BTreeSet<_>>()
                {
                    self.error
                        .get_or_insert(PickleError::InvalidArchive(format!(
                            "canonical computed membership disagrees with old wire in {}",
                            provenance.context
                        )));
                    return;
                }
            }
            self.write_count(provenance.computed.len());
            for name in provenance.computed {
                self.write_string(name);
            }
            return;
        }
        match props {
            Ok(props) => {
                let sorted: BTreeSet<_> = props.unwrap_or_default().iter().collect();
                self.write_count(sorted.len());
                for key in sorted {
                    self.write_string(key);
                }
            }
            Err(e) => {
                self.error.get_or_insert(invalid(e));
            }
        }
    }
}

#[derive(Clone)]
struct LegacyStoreProvenance {
    reserved: bool,
    collision_reserved: Option<PropertyValue>,
    computed: Vec<PropertyText>,
    context: String,
}
fn migrate_legacy_rows(
    mut rows: Vec<(PropertyText, PropertyValue)>,
    computed: Vec<PropertyText>,
    context: String,
) -> Result<(Vec<(PropertyText, PropertyValue)>, LegacyStoreProvenance), PickleError> {
    let reserved = rows
        .iter_mut()
        .find(|(key, _)| key.as_bytes() == b"__computedProps");
    let has_reserved = reserved.is_some();
    // Explicit user policy: a legacy ordinary reserved value plus independent
    // computed flags imports as an empty computed-name vector, not an error.
    // Keep all other property values. Retain the old value only for validating
    // legacy raw/canonical companion bytes; never persist it in archive 2.0.
    let collision_reserved = reserved
        .filter(|_| !computed.is_empty())
        .map(|(_, value)| std::mem::replace(value, PropertyValue::StringVector(Vec::new())));
    let provenance = LegacyStoreProvenance {
        reserved: has_reserved,
        collision_reserved,
        computed: computed.clone(),
        context,
    };
    if !provenance.reserved && !computed.is_empty() {
        rows.push((
            "__computedProps".into(),
            PropertyValue::StringVector(computed),
        ));
    }
    Ok((rows, provenance))
}
fn legacy_map_rows(
    props: &BTreeMap<String, String>,
    computed: &BTreeSet<String>,
    context: String,
) -> Result<(Vec<(PropertyText, PropertyValue)>, LegacyStoreProvenance), PickleError> {
    validate_computed(props, computed)?;
    migrate_legacy_rows(
        props
            .iter()
            .map(|(k, v)| (k.into(), PropertyValue::from(v)))
            .collect(),
        computed.iter().map(Into::into).collect(),
        context,
    )
}
fn legacy_ordered_rows(
    rows: Vec<(String, PropertyValue, bool)>,
    context: String,
) -> Result<Vec<(PropertyText, PropertyValue)>, PickleError> {
    let computed = rows
        .iter()
        .filter(|(_, _, c)| *c)
        .map(|(k, _, _)| k.into())
        .collect();
    let rows = rows.into_iter().map(|(k, v, _)| (k.into(), v)).collect();
    Ok(migrate_legacy_rows(rows, computed, context)?.0)
}

trait LegacyTextValue {
    fn legacy_text(&self) -> Result<PropertyText, PickleError>;
    fn synthetic_reserved(&self) -> bool {
        false
    }
}
impl LegacyTextValue for PropertyText {
    fn legacy_text(&self) -> Result<PropertyText, PickleError> {
        Ok(self.clone())
    }
}
impl LegacyTextValue for PropertyValue {
    fn legacy_text(&self) -> Result<PropertyText, PickleError> {
        property_value_to_string(self).map_err(invalid)
    }
    fn synthetic_reserved(&self) -> bool {
        matches!(self, PropertyValue::StringVector(_))
    }
}

struct PickleReader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> PickleReader<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    fn ensure(&self, n: usize) -> Result<(), PickleError> {
        if n > self.data.len().saturating_sub(self.pos) {
            Err(PickleError::UnexpectedEof)
        } else {
            Ok(())
        }
    }

    fn read_u8(&mut self) -> Result<u8, PickleError> {
        self.ensure(1)?;
        let v = self.data[self.pos];
        self.pos += 1;
        Ok(v)
    }

    fn read_u64(&mut self) -> Result<u64, PickleError> {
        let bytes = self.read_exact_slice(8)?;
        Ok(u64::from_le_bytes(
            bytes.try_into().expect("checked eight bytes"),
        ))
    }

    fn read_count(&mut self, minimum_bytes: usize) -> Result<usize, PickleError> {
        let count = self.read_u32()? as usize;
        if count > 1_000_000 || count > self.remaining() / minimum_bytes.max(1) {
            return Err(PickleError::InvalidArchive(format!(
                "invalid bounded count: {count}"
            )));
        }
        Ok(count)
    }

    fn read_u32(&mut self) -> Result<u32, PickleError> {
        self.ensure(4)?;
        let bytes: [u8; 4] = self.data[self.pos..self.pos + 4].try_into().unwrap();
        self.pos += 4;
        Ok(u32::from_le_bytes(bytes))
    }

    fn read_i32(&mut self) -> Result<i32, PickleError> {
        self.ensure(4)?;
        let bytes: [u8; 4] = self.data[self.pos..self.pos + 4].try_into().unwrap();
        self.pos += 4;
        Ok(i32::from_le_bytes(bytes))
    }

    fn read_i8(&mut self) -> Result<i8, PickleError> {
        Ok(self.read_u8()? as i8)
    }

    fn read_f64(&mut self) -> Result<f64, PickleError> {
        self.ensure(8)?;
        let bytes: [u8; 8] = self.data[self.pos..self.pos + 8].try_into().unwrap();
        self.pos += 8;
        Ok(f64::from_le_bytes(bytes))
    }

    fn read_bool(&mut self) -> Result<bool, PickleError> {
        match self.read_u8()? {
            0 => Ok(false),
            1 => Ok(true),
            value => Err(PickleError::InvalidEnumValue {
                value,
                type_name: "bool",
            }),
        }
    }

    fn read_string(&mut self) -> Result<String, PickleError> {
        let len = self.read_u32()? as usize;
        if len > 10_000_000 {
            return Err(PickleError::StringTooLong(len));
        }
        self.ensure(len)?;
        let s = std::str::from_utf8(&self.data[self.pos..self.pos + len])
            .map_err(|_| PickleError::InvalidMolecule("invalid UTF-8 in pickle string".into()))?;
        self.pos += len;
        Ok(s.to_string())
    }

    fn read_option_string(&mut self) -> Result<Option<String>, PickleError> {
        if self.read_bool()? {
            Ok(Some(self.read_string()?))
        } else {
            Ok(None)
        }
    }

    fn read_props(&mut self) -> Result<BTreeMap<String, String>, PickleError> {
        let count = self.read_u32()? as usize;
        if count > 1_000_000 {
            return Err(PickleError::StringTooLong(count));
        }
        let mut props = BTreeMap::new();
        for _ in 0..count {
            let key = self.read_string()?;
            let value = self.read_string()?;
            if props.insert(key, value).is_some() {
                return Err(PickleError::InvalidArchive(
                    "duplicate property name".into(),
                ));
            }
        }
        Ok(props)
    }

    fn read_computed_props(&mut self) -> Result<BTreeSet<String>, PickleError> {
        let count = self.read_u32()? as usize;
        if count > 1_000_000 {
            return Err(PickleError::InvalidArchive(format!(
                "unreasonable computed-property count: {count}"
            )));
        }
        let mut props = BTreeSet::new();
        for _ in 0..count {
            if !props.insert(self.read_string()?) {
                return Err(PickleError::InvalidArchive(
                    "duplicate computed property name".into(),
                ));
            }
        }
        Ok(props)
    }

    fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }

    fn read_exact_slice(&mut self, len: usize) -> Result<&'a [u8], PickleError> {
        self.ensure(len)?;
        let slice = &self.data[self.pos..self.pos + len];
        self.pos += len;
        Ok(slice)
    }
}

#[derive(Debug, Clone, Copy)]
struct ArchiveSection<'a> {
    id: u16,
    version: u16,
    flags: u8,
    codec: u8,
    payload: &'a [u8],
}

impl ArchiveSection<'_> {
    fn is_required(self) -> bool {
        self.flags & SECTION_FLAG_REQUIRED != 0
    }
}

fn write_u16_le(buf: &mut Vec<u8>, value: u16) {
    buf.extend_from_slice(&value.to_le_bytes());
}

fn write_u32_le(buf: &mut Vec<u8>, value: u32) {
    buf.extend_from_slice(&value.to_le_bytes());
}

fn read_u16_le(r: &mut PickleReader<'_>) -> Result<u16, PickleError> {
    let lo = r.read_u8()? as u16;
    let hi = r.read_u8()? as u16;
    Ok(lo | (hi << 8))
}

fn read_u32_le(r: &mut PickleReader<'_>) -> Result<u32, PickleError> {
    let bytes = r.read_exact_slice(4)?;
    Ok(u32::from_le_bytes(bytes.try_into().unwrap()))
}

fn checked_u32(value: usize, context: &str) -> Result<u32, PickleError> {
    u32::try_from(value).map_err(|_| {
        PickleError::InvalidArchive(format!(
            "{context} exceeds the archive u32 boundary: {value}"
        ))
    })
}

fn write_atom_id_rows(
    w: &mut PickleWriter,
    rows: &[Vec<AtomId>],
    context: &str,
) -> Result<(), PickleError> {
    w.write_u32(checked_u32(rows.len(), context)?);
    for row in rows {
        w.write_u32(checked_u32(row.len(), context)?);
        for id in row {
            w.write_u32(checked_u32(id.index(), context)?);
        }
    }
    Ok(())
}

fn write_bond_id_rows(
    w: &mut PickleWriter,
    rows: &[Vec<BondId>],
    context: &str,
) -> Result<(), PickleError> {
    w.write_u32(checked_u32(rows.len(), context)?);
    for row in rows {
        w.write_u32(checked_u32(row.len(), context)?);
        for id in row {
            w.write_u32(checked_u32(id.index(), context)?);
        }
    }
    Ok(())
}

fn read_row_count(r: &mut PickleReader<'_>, context: &str) -> Result<usize, PickleError> {
    let count = r.read_count(1)?;
    if count > MAX_DERIVED_ROWS {
        return Err(PickleError::InvalidArchive(format!(
            "{context} row count exceeds the archive limit: {count}"
        )));
    }
    Ok(count)
}

fn read_atom_id_rows(
    r: &mut PickleReader<'_>,
    atom_count: usize,
    context: &str,
) -> Result<Vec<Vec<AtomId>>, PickleError> {
    let row_count = read_row_count(r, context)?;
    let mut rows = Vec::with_capacity(row_count);
    for _ in 0..row_count {
        let count = r.read_count(1)?;
        if count > atom_count {
            return Err(PickleError::InvalidArchive(format!(
                "{context} row length {count} exceeds atom count {atom_count}"
            )));
        }
        let mut row = Vec::with_capacity(count);
        for _ in 0..count {
            let index = r.read_u32()? as usize;
            if index >= atom_count {
                return Err(PickleError::InvalidArchive(format!(
                    "{context} atom index {index} is out of range for {atom_count} atoms"
                )));
            }
            row.push(AtomId::new(index));
        }
        rows.push(row);
    }
    Ok(rows)
}

fn read_bond_id_rows(
    r: &mut PickleReader<'_>,
    bond_count: usize,
    context: &str,
) -> Result<Vec<Vec<BondId>>, PickleError> {
    let row_count = read_row_count(r, context)?;
    let mut rows = Vec::with_capacity(row_count);
    for _ in 0..row_count {
        let count = r.read_count(1)?;
        if count > bond_count {
            return Err(PickleError::InvalidArchive(format!(
                "{context} row length {count} exceeds bond count {bond_count}"
            )));
        }
        let mut row = Vec::with_capacity(count);
        for _ in 0..count {
            let index = r.read_u32()? as usize;
            if index >= bond_count {
                return Err(PickleError::InvalidArchive(format!(
                    "{context} bond index {index} is out of range for {bond_count} bonds"
                )));
            }
            row.push(BondId::new(index));
        }
        rows.push(row);
    }
    Ok(rows)
}

fn write_ring_find_type(w: &mut PickleWriter, find_type: RingFindType) {
    w.write_u8(match find_type {
        RingFindType::OtherOrUnknown => 0,
        RingFindType::Fast => 1,
        RingFindType::Sssr => 2,
        RingFindType::SymmSssr => 3,
    });
}

fn read_ring_find_type(r: &mut PickleReader<'_>) -> Result<RingFindType, PickleError> {
    match r.read_u8()? {
        0 => Ok(RingFindType::OtherOrUnknown),
        1 => Ok(RingFindType::Fast),
        2 => Ok(RingFindType::Sssr),
        3 => Ok(RingFindType::SymmSssr),
        value => Err(PickleError::InvalidEnumValue {
            value,
            type_name: "RingFindType",
        }),
    }
}

fn write_ring_info(w: &mut PickleWriter, ring_info: &RingInfo) -> Result<(), PickleError> {
    w.write_bool(ring_info.is_initialized());
    write_ring_find_type(w, ring_info.persisted_find_type());
    write_atom_id_rows(w, ring_info.atom_rings(), "atom rings")?;
    write_bond_id_rows(w, ring_info.bond_rings(), "bond rings")?;
    write_atom_id_rows(w, ring_info.atom_ring_families(), "atom ring families")?;
    write_bond_id_rows(w, ring_info.bond_ring_families(), "bond ring families")?;
    match ring_info.persisted_relevant_cycle_count() {
        Some(count) => {
            w.write_bool(true);
            w.write_u32(checked_u32(count, "relevant cycle count")?);
        }
        None => w.write_bool(false),
    }
    let fused_rings = ring_info.persisted_fused_rings();
    w.write_u32(checked_u32(fused_rings.len(), "fused-ring matrix")?);
    for row in fused_rings {
        w.write_u32(checked_u32(row.len(), "fused-ring matrix row")?);
        for value in row {
            w.write_bool(*value);
        }
    }
    let num_fused_bonds = ring_info.persisted_num_fused_bonds();
    w.write_u32(checked_u32(
        num_fused_bonds.len(),
        "fused-bond count table",
    )?);
    for count in num_fused_bonds {
        w.write_u32(checked_u32(*count, "fused-bond count")?);
    }
    Ok(())
}

fn read_ring_info(
    r: &mut PickleReader<'_>,
    atom_count: usize,
    bond_count: usize,
) -> Result<RingInfo, PickleError> {
    let initialized = r.read_bool()?;
    let find_type = read_ring_find_type(r)?;
    let atom_rings = read_atom_id_rows(r, atom_count, "atom rings")?;
    let bond_rings = read_bond_id_rows(r, bond_count, "bond rings")?;
    let atom_ring_families = read_atom_id_rows(r, atom_count, "atom ring families")?;
    let bond_ring_families = read_bond_id_rows(r, bond_count, "bond ring families")?;
    let relevant_cycle_count = if r.read_bool()? {
        Some(r.read_u32()? as usize)
    } else {
        None
    };
    let fused_row_count = read_row_count(r, "fused-ring matrix")?;
    let mut fused_rings = Vec::with_capacity(fused_row_count);
    for _ in 0..fused_row_count {
        let count = read_row_count(r, "fused-ring matrix row")?;
        let mut row = Vec::with_capacity(count);
        for _ in 0..count {
            row.push(r.read_bool()?);
        }
        fused_rings.push(row);
    }
    let fused_bond_count = read_row_count(r, "fused-bond count table")?;
    let mut num_fused_bonds = Vec::with_capacity(fused_bond_count);
    for _ in 0..fused_bond_count {
        num_fused_bonds.push(r.read_u32()? as usize);
    }

    RingInfo::from_persisted_components(
        initialized,
        find_type,
        atom_count,
        bond_count,
        atom_rings,
        bond_rings,
        atom_ring_families,
        bond_ring_families,
        relevant_cycle_count,
        fused_rings,
        num_fused_bonds,
    )
    .map_err(|message| PickleError::InvalidArchive(message.to_string()))
}

fn encode_derived_state(mol: &BinaryInput<'_>) -> Result<Vec<u8>, PickleError> {
    let cache = &mol.derived;
    let mut w = PickleWriter::new();

    match &cache.rings {
        Some(rings) => {
            w.write_bool(true);
            write_ring_info(&mut w, rings)?;
        }
        None => w.write_bool(false),
    }
    match &cache.ring_families {
        Some(rings) => {
            w.write_bool(true);
            write_ring_info(&mut w, rings)?;
        }
        None => w.write_bool(false),
    }

    match &cache.valence {
        Some(valence) => {
            if valence.explicit_valence.len() != mol.num_atoms()
                || valence.implicit_hydrogens.len() != mol.num_atoms()
            {
                return Err(PickleError::InvalidMolecule(
                    "valence cache row count does not match atom count".to_string(),
                ));
            }
            w.write_bool(true);
            w.write_u32(checked_u32(
                valence.explicit_valence.len(),
                "explicit valence table",
            )?);
            for value in &valence.explicit_valence {
                w.write_i32(*value);
            }
            w.write_u32(checked_u32(
                valence.implicit_hydrogens.len(),
                "implicit hydrogen table",
            )?);
            for value in &valence.implicit_hydrogens {
                w.write_i32(*value);
            }
        }
        None => w.write_bool(false),
    }
    w.write_bool(cache.aromaticity_valid);
    w.write_bool(cache.stereo_valid);
    w.into_inner()
}

fn decode_derived_state(
    data: &[u8],
    mol: &BinaryRecord,
) -> Result<BinaryDerivedState, PickleError> {
    let mut r = PickleReader::new(data);
    let rings = if r.read_bool()? {
        Some(read_ring_info(&mut r, mol.num_atoms(), mol.num_bonds())?)
    } else {
        None
    };
    let ring_families = if r.read_bool()? {
        Some(read_ring_info(&mut r, mol.num_atoms(), mol.num_bonds())?)
    } else {
        None
    };

    let valence = if r.read_bool()? {
        let explicit_count = r.read_count(1)?;
        if explicit_count != mol.num_atoms() {
            return Err(PickleError::InvalidArchive(format!(
                "explicit valence rows {explicit_count} do not match atom count {}",
                mol.num_atoms()
            )));
        }
        let mut explicit_valence = Vec::with_capacity(explicit_count);
        for _ in 0..explicit_count {
            explicit_valence.push(r.read_i32()?);
        }
        let implicit_count = r.read_count(1)?;
        if implicit_count != mol.num_atoms() {
            return Err(PickleError::InvalidArchive(format!(
                "implicit hydrogen rows {implicit_count} do not match atom count {}",
                mol.num_atoms()
            )));
        }
        let mut implicit_hydrogens = Vec::with_capacity(implicit_count);
        for _ in 0..implicit_count {
            implicit_hydrogens.push(r.read_i32()?);
        }
        Some(ValenceAssignment {
            explicit_valence,
            implicit_hydrogens,
        })
    } else {
        None
    };
    let aromaticity_valid = r.read_bool()?;
    let stereo_valid = r.read_bool()?;
    if r.remaining() != 0 {
        return Err(PickleError::InvalidArchive(format!(
            "trailing bytes after derived state: {}",
            r.remaining()
        )));
    }
    Ok(BinaryDerivedState {
        valid_bits: None,
        rings,
        ring_families,
        valence,
        aromaticity_valid,
        stereo_valid,
    })
}

#[cfg(test)]
fn archive_manifest() -> ArchiveManifestV1 {
    ArchiveManifestV1 {
        crate_version: env!("CARGO_PKG_VERSION").to_string(),
        molecule_state_codec: SECTION_CODEC_POSTCARD,
        molecule_state_version: MOLECULE_STATE_VERSION,
    }
}

#[cfg(test)]
fn encode_manifest() -> Result<Vec<u8>, PickleError> {
    postcard::to_allocvec(&archive_manifest())
        .map_err(|err| PickleError::InvalidArchive(format!("manifest encode failed: {err}")))
}

#[cfg(test)]
fn encode_molecule_state(payload: Vec<u8>) -> Result<Vec<u8>, PickleError> {
    let state = MoleculeStateV1 {
        encoding: SECTION_CODEC_RAW,
        encoding_version: *payload.first().ok_or(PickleError::UnexpectedEof)?,
        payload,
    };
    postcard::to_allocvec(&state)
        .map_err(|err| PickleError::InvalidArchive(format!("molecule state encode failed: {err}")))
}

fn write_archive_section(
    buf: &mut Vec<u8>,
    id: u16,
    version: u16,
    flags: u8,
    codec: u8,
    payload: &[u8],
) -> Result<(), PickleError> {
    if payload.len() > u32::MAX as usize {
        return Err(PickleError::DataLengthMismatch {
            expected: u32::MAX as usize,
            actual: payload.len(),
        });
    }
    write_u16_le(buf, id);
    write_u16_le(buf, version);
    buf.push(flags);
    buf.push(codec);
    write_u32_le(buf, payload.len() as u32);
    buf.extend_from_slice(payload);
    Ok(())
}

#[cfg(test)]
fn encode_sectioned_archive(
    molecule_state: Vec<u8>,
    derived_state: Vec<u8>,
    canonical_state: Vec<u8>,
) -> Result<Vec<u8>, PickleError> {
    let manifest = encode_manifest()?;
    let molecule_state = encode_molecule_state(molecule_state)?;
    let mut buf = Vec::with_capacity(
        ARCHIVE_MAGIC.len()
            + 2
            + 2
            + 2
            + 30
            + manifest.len()
            + molecule_state.len()
            + derived_state.len(),
    );
    buf.extend_from_slice(ARCHIVE_MAGIC);
    write_u16_le(&mut buf, ARCHIVE_MAJOR);
    write_u16_le(&mut buf, ARCHIVE_MINOR);
    write_u16_le(&mut buf, 4);
    write_archive_section(
        &mut buf,
        SECTION_MANIFEST,
        MANIFEST_VERSION,
        0,
        SECTION_CODEC_POSTCARD,
        &manifest,
    )?;
    write_archive_section(
        &mut buf,
        SECTION_MOLECULE_STATE,
        MOLECULE_STATE_VERSION,
        SECTION_FLAG_REQUIRED,
        SECTION_CODEC_POSTCARD,
        &molecule_state,
    )?;
    write_archive_section(
        &mut buf,
        SECTION_DERIVED_STATE,
        DERIVED_STATE_VERSION,
        SECTION_FLAG_REQUIRED,
        SECTION_CODEC_RAW,
        &derived_state,
    )?;
    write_archive_section(
        &mut buf,
        SECTION_CANONICAL_STATE,
        CANONICAL_STATE_VERSION,
        SECTION_FLAG_REQUIRED,
        SECTION_CODEC_RAW,
        &canonical_state,
    )?;
    Ok(buf)
}

fn read_archive_sections<'a>(
    data: &'a [u8],
) -> Result<(u16, Vec<ArchiveSection<'a>>), PickleError> {
    let (major, minor, sections) = read_archive_envelope(data, ARCHIVE_MAGIC)?;
    if major != ARCHIVE_MAJOR || minor > ARCHIVE_MINOR {
        return Err(PickleError::UnsupportedArchiveVersion { major, minor });
    }
    Ok((minor, sections))
}

fn read_archive_envelope<'a>(
    data: &'a [u8],
    expected_magic: &[u8; 8],
) -> Result<(u16, u16, Vec<ArchiveSection<'a>>), PickleError> {
    let mut r = PickleReader::new(data);
    let magic = r.read_exact_slice(expected_magic.len())?;
    if magic != expected_magic {
        return Err(PickleError::InvalidArchive("magic mismatch".to_string()));
    }
    let major = read_u16_le(&mut r)?;
    let minor = read_u16_le(&mut r)?;
    let section_count = read_u16_le(&mut r)? as usize;
    if section_count > 1024 {
        return Err(PickleError::InvalidArchive(format!(
            "unreasonable section count: {section_count}"
        )));
    }
    let mut sections = Vec::with_capacity(section_count);
    for _ in 0..section_count {
        let id = read_u16_le(&mut r)?;
        let version = read_u16_le(&mut r)?;
        let flags = r.read_u8()?;
        if flags & !SECTION_FLAG_REQUIRED != 0 {
            return Err(PickleError::InvalidArchive("unknown section flags".into()));
        }
        let codec = r.read_u8()?;
        let len = read_u32_le(&mut r)? as usize;
        let payload = r.read_exact_slice(len)?;
        sections.push(ArchiveSection {
            id,
            version,
            flags,
            codec,
            payload,
        });
    }
    if r.remaining() != 0 {
        return Err(PickleError::InvalidArchive(format!(
            "trailing bytes after archive sections: {}",
            r.remaining()
        )));
    }
    Ok((major, minor, sections))
}

fn decode_sectioned_archive(data: &[u8]) -> Result<BinaryRecord, PickleError> {
    let (archive_minor, sections) = read_archive_sections(data)?;
    let mut seen = BTreeSet::new();
    let mut canonical_state = None;
    let mut manifest_seen = false;
    let mut molecule_state: Option<&[u8]> = None;
    let mut derived_state: Option<&[u8]> = None;
    for section in sections {
        if !seen.insert(section.id) {
            return Err(PickleError::DuplicateSection(section.id));
        }
        if (section.id == SECTION_MOLECULE_STATE
            || (archive_minor >= 1 && section.id == SECTION_DERIVED_STATE)
            || (archive_minor >= 2 && section.id == SECTION_CANONICAL_STATE))
            && !section.is_required()
        {
            return Err(PickleError::InvalidArchive(
                "required section missing required flag".into(),
            ));
        }
        match section.id {
            SECTION_MANIFEST => {
                if manifest_seen {
                    return Err(PickleError::DuplicateSection(SECTION_MANIFEST));
                }
                manifest_seen = true;
                if section.version != MANIFEST_VERSION {
                    return Err(PickleError::UnsupportedSectionVersion {
                        section: SECTION_MANIFEST,
                        version: section.version,
                    });
                }
                if section.codec != SECTION_CODEC_POSTCARD {
                    return Err(PickleError::InvalidArchive(format!(
                        "manifest uses unsupported codec {}",
                        section.codec
                    )));
                }
                let manifest: ArchiveManifestV1 =
                    postcard_exact(section.payload).map_err(|err| {
                        PickleError::InvalidArchive(format!("manifest decode failed: {err}"))
                    })?;
                if manifest.molecule_state_codec != SECTION_CODEC_POSTCARD {
                    return Err(PickleError::InvalidArchive(format!(
                        "manifest declares unsupported molecule state codec {}",
                        manifest.molecule_state_codec
                    )));
                }
                if manifest.molecule_state_version != MOLECULE_STATE_VERSION {
                    return Err(PickleError::UnsupportedSectionVersion {
                        section: SECTION_MOLECULE_STATE,
                        version: manifest.molecule_state_version,
                    });
                }
            }
            SECTION_MOLECULE_STATE => {
                if molecule_state.is_some() {
                    return Err(PickleError::DuplicateSection(SECTION_MOLECULE_STATE));
                }
                if section.version != MOLECULE_STATE_VERSION {
                    return Err(PickleError::UnsupportedSectionVersion {
                        section: SECTION_MOLECULE_STATE,
                        version: section.version,
                    });
                }
                if section.codec != SECTION_CODEC_POSTCARD {
                    return Err(PickleError::InvalidArchive(format!(
                        "molecule state uses unsupported codec {}",
                        section.codec
                    )));
                }
                molecule_state = Some(section.payload);
            }
            SECTION_DERIVED_STATE => {
                if derived_state.is_some() {
                    return Err(PickleError::DuplicateSection(SECTION_DERIVED_STATE));
                }
                if section.version != DERIVED_STATE_VERSION {
                    return Err(PickleError::UnsupportedSectionVersion {
                        section: SECTION_DERIVED_STATE,
                        version: section.version,
                    });
                }
                if section.codec != SECTION_CODEC_RAW {
                    return Err(PickleError::InvalidArchive(format!(
                        "derived state uses unsupported codec {}",
                        section.codec
                    )));
                }
                derived_state = Some(section.payload);
            }
            SECTION_CANONICAL_STATE => {
                if !((archive_minor == 2 && (1..=2).contains(&section.version))
                    || (archive_minor == 3 && section.version == 3))
                {
                    return Err(PickleError::UnsupportedSectionVersion {
                        section: section.id,
                        version: section.version,
                    });
                }
                if section.codec != SECTION_CODEC_RAW {
                    return Err(PickleError::InvalidArchive("canonical state codec".into()));
                }
                canonical_state = Some((section.version, section.payload));
            }
            unknown => {
                if section.is_required() {
                    return Err(PickleError::UnknownRequiredSection(unknown));
                }
            }
        }
    }
    if !manifest_seen {
        return Err(PickleError::MissingRequiredSection(SECTION_MANIFEST));
    }
    let Some(molecule_state) = molecule_state else {
        return Err(PickleError::MissingRequiredSection(SECTION_MOLECULE_STATE));
    };
    if archive_minor >= 1 && derived_state.is_none() {
        return Err(PickleError::MissingRequiredSection(SECTION_DERIVED_STATE));
    }
    let state: MoleculeStateV1 = postcard_exact(molecule_state).map_err(|err| {
        PickleError::InvalidArchive(format!("molecule state decode failed: {err}"))
    })?;
    if state.encoding != SECTION_CODEC_RAW {
        return Err(PickleError::InvalidArchive(format!(
            "molecule state declares unsupported encoding {}",
            state.encoding
        )));
    }
    if !(1..=PICKLE_VERSION).contains(&state.encoding_version) {
        return Err(PickleError::UnsupportedVersion(state.encoding_version));
    }
    if state.payload.first() != Some(&state.encoding_version) {
        return Err(PickleError::InvalidArchive(
            "payload version disagrees with envelope".into(),
        ));
    }
    if archive_minor == 3 {
        if state.encoding_version != 4 {
            return Err(PickleError::InvalidArchive(
                "archive1.3 requires raw4".into(),
            ));
        }
        let (version, data) =
            canonical_state.ok_or(PickleError::MissingRequiredSection(SECTION_CANONICAL_STATE))?;
        if version != 3 {
            return Err(PickleError::UnsupportedSectionVersion {
                section: 4,
                version,
            });
        }
        let result = native_state_v2::decode(data)?;
        let canonical = native_state_v2::encode_record(&result)?;
        if canonical.as_slice() != data {
            return Err(PickleError::InvalidArchive(
                "canonical3 disagrees with validated native state".into(),
            ));
        }
        let mut raw = Vec::with_capacity(canonical.len() + 1);
        raw.push(state.encoding_version);
        raw.extend_from_slice(&canonical);
        if raw != state.payload {
            return Err(PickleError::InvalidArchive(
                "canonical state disagrees with raw4 companion".into(),
            ));
        }
        let input = BinaryInput {
            topology: &result.topology,
            coordinates: &result.coordinates,
            properties: &result.properties,
            derived: BinaryDerivedView {
                rings: result.derived.rings.as_ref(),
                ring_families: result.derived.ring_families.as_ref(),
                valence: result.derived.valence.as_ref(),
                aromaticity_valid: result.derived.aromaticity_valid,
                stereo_valid: result.derived.stereo_valid,
                valid_bits: result.derived.valid_bits.unwrap_or(0),
            },
        };
        if encode_derived_state(&input)?.as_slice()
            != derived_state.ok_or(PickleError::MissingRequiredSection(SECTION_DERIVED_STATE))?
        {
            return Err(PickleError::InvalidArchive(
                "NativeStateV2 disagrees with derived1 companion".into(),
            ));
        }
        return Ok(result);
    }
    if state.encoding_version > 3 {
        return Err(PickleError::InvalidArchive(
            "old archive cannot contain raw4".into(),
        ));
    }
    let (mut molecule, provenance) = mol_from_legacy_binary_with_provenance(&state.payload)?;
    if let Some(data) = derived_state {
        molecule.derived = decode_derived_state(data, &molecule)?;
    }
    if archive_minor >= 2 {
        let data =
            canonical_state.ok_or(PickleError::MissingRequiredSection(SECTION_CANONICAL_STATE))?;
        decode_canonical_state(data.1, data.0, &mut molecule)?;
        let input = BinaryInput {
            topology: &molecule.topology,
            coordinates: &molecule.coordinates,
            properties: &molecule.properties,
            derived: BinaryDerivedView::default(),
        };
        if mol_to_legacy_binary_version_with_provenance(
            &input,
            state.encoding_version,
            Some(&provenance),
        )? != state.payload
        {
            return Err(PickleError::InvalidArchive(
                "canonical state disagrees with legacy companion".into(),
            ));
        }
    }
    molecule.validate()?;
    Ok(molecule)
}

// ──────────────────────────────────────────────
// Enum serialization helpers
// ──────────────────────────────────────────────

fn write_chiral_tag(w: &mut PickleWriter, tag: ChiralTag) {
    let code: u8 = match tag {
        ChiralTag::Unspecified => 0,
        ChiralTag::TetrahedralCw => 1,
        ChiralTag::TetrahedralCcw => 2,
        ChiralTag::Other => 3,
        ChiralTag::Tetrahedral => 4,
        ChiralTag::Allene => 5,
        ChiralTag::SquarePlanar => 6,
        ChiralTag::TrigonalBipyramidal => 7,
        ChiralTag::Octahedral => 8,
    };
    w.write_u8(code);
}

fn read_chiral_tag(r: &mut PickleReader) -> Result<ChiralTag, PickleError> {
    match r.read_u8()? {
        0 => Ok(ChiralTag::Unspecified),
        1 => Ok(ChiralTag::TetrahedralCw),
        2 => Ok(ChiralTag::TetrahedralCcw),
        3 => Ok(ChiralTag::Other),
        4 => Ok(ChiralTag::Tetrahedral),
        5 => Ok(ChiralTag::Allene),
        6 => Ok(ChiralTag::SquarePlanar),
        7 => Ok(ChiralTag::TrigonalBipyramidal),
        8 => Ok(ChiralTag::Octahedral),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "ChiralTag",
        }),
    }
}

fn write_hybridization(w: &mut PickleWriter, h: Hybridization) {
    let code: u8 = match h {
        Hybridization::Unspecified => 0,
        Hybridization::S => 1,
        Hybridization::Sp => 2,
        Hybridization::Sp2 => 3,
        Hybridization::Sp3 => 4,
        Hybridization::Sp2d => 5,
        Hybridization::Sp3d => 6,
        Hybridization::Sp3d2 => 7,
        Hybridization::Other => 8,
    };
    w.write_u8(code);
}

fn read_hybridization(r: &mut PickleReader) -> Result<Hybridization, PickleError> {
    match r.read_u8()? {
        0 => Ok(Hybridization::Unspecified),
        1 => Ok(Hybridization::S),
        2 => Ok(Hybridization::Sp),
        3 => Ok(Hybridization::Sp2),
        4 => Ok(Hybridization::Sp3),
        5 => Ok(Hybridization::Sp2d),
        6 => Ok(Hybridization::Sp3d),
        7 => Ok(Hybridization::Sp3d2),
        8 => Ok(Hybridization::Other),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "Hybridization",
        }),
    }
}

fn write_bond_order(w: &mut PickleWriter, order: BondOrder) {
    let code: u8 = match order {
        BondOrder::Unspecified => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Quintuple => 5,
        BondOrder::Hextuple => 6,
        BondOrder::OneAndHalf => 7,
        BondOrder::TwoAndHalf => 8,
        BondOrder::ThreeAndHalf => 9,
        BondOrder::FourAndHalf => 10,
        BondOrder::FiveAndHalf => 11,
        BondOrder::Aromatic => 12,
        BondOrder::Ionic => 13,
        BondOrder::Dative => 14,
        BondOrder::DativeOne => 15,
        BondOrder::DativeLeft => 16,
        BondOrder::DativeRight => 17,
        BondOrder::Hydrogen => 18,
        BondOrder::ThreeCenter => 19,
        BondOrder::Other => 20,
        BondOrder::Zero => 21,
    };
    w.write_u8(code);
}

fn read_bond_order(r: &mut PickleReader) -> Result<BondOrder, PickleError> {
    match r.read_u8()? {
        0 => Ok(BondOrder::Unspecified),
        1 => Ok(BondOrder::Single),
        2 => Ok(BondOrder::Double),
        3 => Ok(BondOrder::Triple),
        4 => Ok(BondOrder::Quadruple),
        5 => Ok(BondOrder::Quintuple),
        6 => Ok(BondOrder::Hextuple),
        7 => Ok(BondOrder::OneAndHalf),
        8 => Ok(BondOrder::TwoAndHalf),
        9 => Ok(BondOrder::ThreeAndHalf),
        10 => Ok(BondOrder::FourAndHalf),
        11 => Ok(BondOrder::FiveAndHalf),
        12 => Ok(BondOrder::Aromatic),
        13 => Ok(BondOrder::Ionic),
        14 => Ok(BondOrder::Dative),
        15 => Ok(BondOrder::DativeOne),
        16 => Ok(BondOrder::DativeLeft),
        17 => Ok(BondOrder::DativeRight),
        18 => Ok(BondOrder::Hydrogen),
        19 => Ok(BondOrder::ThreeCenter),
        20 => Ok(BondOrder::Other),
        21 => Ok(BondOrder::Zero),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "BondOrder",
        }),
    }
}

fn write_bond_direction(w: &mut PickleWriter, dir: BondDirection) {
    let code: u8 = match dir {
        BondDirection::None => 0,
        BondDirection::BeginWedge => 1,
        BondDirection::BeginDash => 2,
        BondDirection::EndUpRight => 3,
        BondDirection::EndDownRight => 4,
        BondDirection::EitherDouble => 5,
        BondDirection::Unknown => 6,
    };
    w.write_u8(code);
}

fn read_bond_direction(r: &mut PickleReader) -> Result<BondDirection, PickleError> {
    match r.read_u8()? {
        0 => Ok(BondDirection::None),
        1 => Ok(BondDirection::BeginWedge),
        2 => Ok(BondDirection::BeginDash),
        3 => Ok(BondDirection::EndUpRight),
        4 => Ok(BondDirection::EndDownRight),
        5 => Ok(BondDirection::EitherDouble),
        6 => Ok(BondDirection::Unknown),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "BondDirection",
        }),
    }
}

fn write_bond_stereo(w: &mut PickleWriter, stereo: BondStereo) {
    let code: u8 = match stereo {
        BondStereo::None => 0,
        BondStereo::Any => 1,
        BondStereo::Z => 2,
        BondStereo::E => 3,
        BondStereo::Cis => 4,
        BondStereo::Trans => 5,
        BondStereo::AtropCw => 6,
        BondStereo::AtropCcw => 7,
    };
    w.write_u8(code);
}

fn read_bond_stereo(r: &mut PickleReader) -> Result<BondStereo, PickleError> {
    match r.read_u8()? {
        0 => Ok(BondStereo::None),
        1 => Ok(BondStereo::Any),
        2 => Ok(BondStereo::Z),
        3 => Ok(BondStereo::E),
        4 => Ok(BondStereo::Cis),
        5 => Ok(BondStereo::Trans),
        6 => Ok(BondStereo::AtropCw),
        7 => Ok(BondStereo::AtropCcw),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "BondStereo",
        }),
    }
}

fn write_atom(w: &mut PickleWriter, atom: &Atom, version: u8) {
    // Atomic number
    w.write_u8(atom.atomic_number());

    // Formal charge
    w.write_i8(atom.formal_charge());

    // Isotope
    if let Some(isotope) = atom.isotope() {
        w.write_bool(true);
        w.write_u32(u32::from(isotope));
    } else {
        w.write_bool(false);
    }

    // Chiral tag
    write_chiral_tag(w, atom.chiral_tag());

    // Chiral permutation
    if let Some(perm) = atom.chiral_permutation() {
        w.write_bool(true);
        w.write_u32(perm);
    } else {
        w.write_bool(false);
    }

    // Unknown stereo
    w.write_bool(atom.unknown_stereo());

    // Mol parity
    if let Some(parity) = atom.mol_parity() {
        w.write_bool(true);
        w.write_i32(parity);
    } else {
        w.write_bool(false);
    }

    // Mol inversion flag
    if let Some(inv) = atom.mol_inversion_flag() {
        w.write_bool(true);
        w.write_i32(inv);
    } else {
        w.write_bool(false);
    }

    // Radical electrons
    w.write_u8(atom.radical_electrons());

    // Is aromatic
    w.write_bool(atom.is_aromatic());

    // Hybridization
    write_hybridization(w, atom.hybridization());

    // Atom map
    if let Some(map) = atom.atom_map() {
        w.write_bool(true);
        w.write_u32(map);
    } else {
        w.write_bool(false);
    }

    // No implicit
    w.write_bool(atom.no_implicit());

    // Implicit hydrogen flag
    w.write_bool(atom.implicit_hydrogen());

    // Explicit hydrogens count
    w.write_u8(atom.explicit_hydrogens());

    // Tracked isotopic hydrogens
    let tracked_isotopes = atom.tracked_isotopic_hydrogens();
    w.write_count(tracked_isotopes.len());
    for &iso in tracked_isotopes {
        w.write_u32(u32::from(iso));
    }

    // Concrete atoms do not carry query state; query graphs have their own
    // serialization path and are never encoded as molecule atom flags.
    w.write_bool(false);

    // Properties
    w.write_typed_props(atom.props());
    if version >= 3 {
        w.write_computed_props(atom.computed_prop_names());
    }

    // PDB residue info presence flag
    w.write_bool(atom.pdb_residue_info().is_some());
}

fn read_bond(
    r: &mut PickleReader,
    version: u8,
    id: BondId,
    provenance: &mut Vec<LegacyStoreProvenance>,
) -> Result<Bond, PickleError> {
    let begin_idx = r.read_u32()? as usize;
    let end_idx = r.read_u32()? as usize;
    let order = read_bond_order(r)?;
    let stereo = read_bond_stereo(r)?;
    let direction = read_bond_direction(r)?;
    let is_aromatic = r.read_bool()?;
    let is_conjugated = r.read_bool()?;

    // Stereo atoms
    let stereo_atoms = if r.read_bool()? {
        let sa_begin = AtomId::new(r.read_u32()? as usize);
        let sa_end = AtomId::new(r.read_u32()? as usize);
        Some([sa_begin, sa_end])
    } else {
        None
    };

    let unknown_stereo = r.read_bool()?;

    // Query presence flag from legacy molecule pickles (query graphs are now
    // separate values and are not reconstructed into concrete atoms).
    if r.read_bool()? {
        return Err(PickleError::InvalidMolecule(
            "query state in concrete molecule payload".into(),
        ));
    }

    let props = r.read_props()?;
    let computed_props = if version >= 3 {
        r.read_computed_props()?
    } else {
        BTreeSet::new()
    };

    // Reconstruct via builder-like approach using spec
    // Note: BondSpec is the construction payload; we need id, begin, end
    let spec = BondSpec::new(AtomId::new(begin_idx), AtomId::new(end_idx), order)
        .with_stereo(stereo)
        .with_direction(direction)
        .with_aromatic(is_aromatic)
        .with_conjugated(is_conjugated)
        .with_unknown_stereo(unknown_stereo);

    let spec = if let Some([sa_begin, sa_end]) = stereo_atoms {
        spec.with_stereo_atoms(sa_begin, sa_end)
    } else {
        spec
    };

    let mut spec = spec;
    let (rows, original) = legacy_map_rows(
        &props,
        &computed_props,
        format!("raw{version} bond {}", id.index()),
    )?;
    for (key, value) in rows {
        spec = spec.with_prop(key, value).map_err(invalid)?;
    }
    provenance.push(original);
    // Creating from_spec requires BondId
    Ok(Bond::from_spec(id, spec))
}

fn write_bond(w: &mut PickleWriter, bond: &Bond, version: u8) {
    w.write_count(bond.begin().index());
    w.write_count(bond.end().index());
    write_bond_order(w, bond.order());
    write_bond_stereo(w, bond.stereo());
    write_bond_direction(w, bond.direction());
    w.write_bool(bond.is_aromatic());
    w.write_bool(bond.is_conjugated());

    // Stereo atoms
    if let Some([sa_begin, sa_end]) = bond.stereo_atoms() {
        w.write_bool(true);
        w.write_count(sa_begin.index());
        w.write_count(sa_end.index());
    } else {
        w.write_bool(false);
    }

    w.write_bool(bond.unknown_stereo());

    // Query graphs are separate values, not bond flags on concrete molecules.
    w.write_bool(false);

    // Properties
    w.write_typed_props(bond.props());
    if version >= 3 {
        w.write_computed_props(bond.computed_prop_names());
    }
}

// ──────────────────────────────────────────────
// Substance group serialization
// ──────────────────────────────────────────────

fn write_substance_group_kind(w: &mut PickleWriter, kind: &SubstanceGroupKind) {
    let (code, generic_name): (u8, Option<&PropertyText>) = match kind {
        SubstanceGroupKind::Data => (0, None),
        SubstanceGroupKind::Superatom => (1, None),
        SubstanceGroupKind::MultipleGroup => (2, None),
        SubstanceGroupKind::StructuralRepeatUnit => (3, None),
        SubstanceGroupKind::Monomer => (4, None),
        SubstanceGroupKind::Copolymer => (5, None),
        SubstanceGroupKind::Crosslink => (6, None),
        SubstanceGroupKind::Graft => (7, None),
        SubstanceGroupKind::Modification => (8, None),
        SubstanceGroupKind::Mer => (9, None),
        SubstanceGroupKind::AnyPolymer => (10, None),
        SubstanceGroupKind::MixtureComponent => (11, None),
        SubstanceGroupKind::Mixture => (12, None),
        SubstanceGroupKind::Formulation => (13, None),
        SubstanceGroupKind::Generic(name) => (14, Some(name)),
    };
    w.write_u8(code);
    if let Some(name) = generic_name {
        w.write_string(name);
    }
}

fn read_substance_group_kind(r: &mut PickleReader) -> Result<SubstanceGroupKind, PickleError> {
    match r.read_u8()? {
        0 => Ok(SubstanceGroupKind::Data),
        1 => Ok(SubstanceGroupKind::Superatom),
        2 => Ok(SubstanceGroupKind::MultipleGroup),
        3 => Ok(SubstanceGroupKind::StructuralRepeatUnit),
        4 => Ok(SubstanceGroupKind::Monomer),
        5 => Ok(SubstanceGroupKind::Copolymer),
        6 => Ok(SubstanceGroupKind::Crosslink),
        7 => Ok(SubstanceGroupKind::Graft),
        8 => Ok(SubstanceGroupKind::Modification),
        9 => Ok(SubstanceGroupKind::Mer),
        10 => Ok(SubstanceGroupKind::AnyPolymer),
        11 => Ok(SubstanceGroupKind::MixtureComponent),
        12 => Ok(SubstanceGroupKind::Mixture),
        13 => Ok(SubstanceGroupKind::Formulation),
        14 => {
            let name = r.read_string()?;
            Ok(SubstanceGroupKind::Generic(name.into()))
        }
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "SubstanceGroupKind",
        }),
    }
}

fn write_sgroup_connection(w: &mut PickleWriter, conn: Option<&SGroupConnection>) {
    match conn {
        None => w.write_u8(0),
        Some(SGroupConnection::HeadToHead) => w.write_u8(1),
        Some(SGroupConnection::HeadToTail) => w.write_u8(2),
        Some(SGroupConnection::Either) => w.write_u8(3),
        Some(SGroupConnection::Unknown(s)) => {
            w.write_u8(4);
            w.write_string(s);
        }
    }
}

fn read_sgroup_connection(r: &mut PickleReader) -> Result<Option<SGroupConnection>, PickleError> {
    match r.read_u8()? {
        0 => Ok(None),
        1 => Ok(Some(SGroupConnection::HeadToHead)),
        2 => Ok(Some(SGroupConnection::HeadToTail)),
        3 => Ok(Some(SGroupConnection::Either)),
        4 => Ok(Some(SGroupConnection::Unknown(r.read_string()?.into()))),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "SGroupConnection",
        }),
    }
}

fn write_sgroup_bracket_style(w: &mut PickleWriter, style: Option<&SGroupBracketStyle>) {
    match style {
        None => w.write_u8(0),
        Some(SGroupBracketStyle::Bracket) => w.write_u8(1),
        Some(SGroupBracketStyle::Parenthesis) => w.write_u8(2),
        Some(SGroupBracketStyle::None) => w.write_u8(3),
        Some(SGroupBracketStyle::Unknown(s)) => {
            w.write_u8(4);
            w.write_string(s);
        }
    }
}

fn read_sgroup_bracket_style(
    r: &mut PickleReader,
) -> Result<Option<SGroupBracketStyle>, PickleError> {
    match r.read_u8()? {
        0 => Ok(None),
        1 => Ok(Some(SGroupBracketStyle::Bracket)),
        2 => Ok(Some(SGroupBracketStyle::Parenthesis)),
        3 => Ok(Some(SGroupBracketStyle::None)),
        4 => Ok(Some(SGroupBracketStyle::Unknown(r.read_string()?.into()))),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "SGroupBracketStyle",
        }),
    }
}

fn write_sgroup_display(w: &mut PickleWriter, display: Option<&SGroupDisplay>) {
    match display {
        None => w.write_bool(false),
        Some(d) => {
            w.write_bool(true);
            // Brackets
            w.write_count(d.brackets.len());
            for bracket in &d.brackets {
                w.write_f64(bracket.points[0][0]);
                w.write_f64(bracket.points[0][1]);
                w.write_f64(bracket.points[1][0]);
                w.write_f64(bracket.points[1][1]);
            }
            // Field position
            match d.field_position {
                Some(pos) => {
                    w.write_bool(true);
                    w.write_f64(pos[0]);
                    w.write_f64(pos[1]);
                }
                None => w.write_bool(false),
            }
            // Display tag
            w.write_option_string(d.display_tag.as_ref());
        }
    }
}

fn read_sgroup_display(r: &mut PickleReader) -> Result<Option<SGroupDisplay>, PickleError> {
    if !r.read_bool()? {
        return Ok(None);
    }
    let mut display = SGroupDisplay::default();
    let bracket_count = r.read_count(1)?;
    for _ in 0..bracket_count {
        let p1x = r.read_f64()?;
        let p1y = r.read_f64()?;
        let p2x = r.read_f64()?;
        let p2y = r.read_f64()?;
        display.brackets.push(SGroupBracket::new([
            [p1x, p1y, 0.0],
            [p2x, p2y, 0.0],
            [0.0; 3],
        ]));
    }
    if r.read_bool()? {
        let fx = r.read_f64()?;
        let fy = r.read_f64()?;
        display.field_position = Some([fx, fy]);
    }
    display.display_tag = r.read_option_string()?.map(Into::into);
    Ok(Some(display))
}

fn write_sgroup_data(w: &mut PickleWriter, data: Option<&SGroupData>) {
    match data {
        None => w.write_bool(false),
        Some(d) => {
            w.write_bool(true);
            w.write_option_string(d.field_name.as_ref());
            w.write_option_string(d.field_type.as_ref());
            w.write_option_string(d.field_info.as_ref());
            w.write_option_string(d.field_display.as_ref());
            w.write_option_string(d.units.as_ref());
            w.write_option_string(d.query_type.as_ref());
            w.write_option_string(d.query_op.as_ref());
            w.write_count(d.values.len());
            for v in &d.values {
                w.write_string(v);
            }
        }
    }
}

fn read_sgroup_data(r: &mut PickleReader) -> Result<Option<SGroupData>, PickleError> {
    if !r.read_bool()? {
        return Ok(None);
    }
    let mut data = SGroupData::default();
    data.field_name = r.read_option_string()?.map(Into::into);
    data.field_type = r.read_option_string()?.map(Into::into);
    data.field_info = r.read_option_string()?.map(Into::into);
    data.field_display = r.read_option_string()?.map(Into::into);
    data.units = r.read_option_string()?.map(Into::into);
    data.query_type = r.read_option_string()?.map(Into::into);
    data.query_op = r.read_option_string()?.map(Into::into);
    let val_count = r.read_count(1)?;
    for _ in 0..val_count {
        data.values.push(r.read_string()?.into());
    }
    Ok(Some(data))
}

fn write_substance_group(w: &mut PickleWriter, sg: &SubstanceGroup) {
    // id index
    w.write_count(sg.id().index());

    // rdkit sequence id
    if let Some(seq_id) = sg.rdkit_sequence_id() {
        w.write_bool(true);
        w.write_u32(seq_id);
    } else {
        w.write_bool(false);
    }

    // external id
    if let Some(ext_id) = sg.external_id() {
        w.write_bool(true);
        w.write_u32(ext_id);
    } else {
        w.write_bool(false);
    }

    // kind
    write_substance_group_kind(w, sg.kind());

    // atoms
    let atoms = sg.atoms();
    w.write_count(atoms.len());
    for &a in atoms {
        w.write_count(a.index());
    }

    // bonds
    let bonds = sg.bonds();
    w.write_count(bonds.len());
    for &b in bonds {
        w.write_count(b.index());
    }

    // bond roles — only write non-default roles
    let mut roles_written = 0u32;
    for &b in bonds {
        if sg.bond_role(b) == SGroupBondRole::Contained {
            roles_written += 1;
        }
    }
    w.write_u32(roles_written);
    for &b in bonds {
        if sg.bond_role(b) == SGroupBondRole::Contained {
            w.write_count(b.index());
        }
    }

    // parent atoms
    let parent_atoms = sg.parent_atoms();
    w.write_count(parent_atoms.len());
    for &a in parent_atoms {
        w.write_count(a.index());
    }

    // parent
    if let Some(parent) = sg.parent() {
        w.write_bool(true);
        w.write_count(parent.index());
    } else {
        w.write_bool(false);
    }

    // label
    w.write_option_string(sg.label());

    // connection
    write_sgroup_connection(w, sg.connection());

    // subtype
    w.write_option_string(sg.subtype());

    // bracket style
    write_sgroup_bracket_style(w, sg.bracket_style());

    // expansion state
    w.write_option_string(sg.expansion_state());

    // class
    w.write_option_string(sg.class());

    // component number
    if let Some(cn) = sg.component_number() {
        w.write_bool(true);
        w.write_u32(cn);
    } else {
        w.write_bool(false);
    }

    // display
    write_sgroup_display(w, sg.display());

    // data
    write_sgroup_data(w, sg.data());

    // attach points
    let attach_pts = sg.attach_points();
    w.write_count(attach_pts.len());
    for ap in attach_pts {
        w.write_count(ap.atom.index());
        if let Some(la) = ap.leaving_atom {
            w.write_bool(true);
            w.write_count(la.index());
        } else {
            w.write_bool(false);
        }
        w.write_option_string(ap.label.as_ref());
        if let Some(order) = ap.order {
            w.write_bool(true);
            w.write_u32(order);
        } else {
            w.write_bool(false);
        }
    }

    // cstates
    let cstates = sg.cstates();
    w.write_count(cstates.len());
    for cs in cstates {
        w.write_count(cs.bond.index());
        w.write_f64(cs.vector[0]);
        w.write_f64(cs.vector[1]);
    }

    // props
    w.write_props(sg.props());

    // data fields
    let data_fields = sg.data_fields();
    w.write_count(data_fields.len());
    for df in data_fields {
        w.write_string(df);
    }
}

// ──────────────────────────────────────────────
// Stereo group serialization
// ──────────────────────────────────────────────

fn write_stereo_group_kind(w: &mut PickleWriter, kind: StereoGroupKind) {
    let code: u8 = match kind {
        StereoGroupKind::Absolute => 0,
        StereoGroupKind::Or => 1,
        StereoGroupKind::And => 2,
    };
    w.write_u8(code);
}

fn read_stereo_group_kind(r: &mut PickleReader) -> Result<StereoGroupKind, PickleError> {
    match r.read_u8()? {
        0 => Ok(StereoGroupKind::Absolute),
        1 => Ok(StereoGroupKind::Or),
        2 => Ok(StereoGroupKind::And),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "StereoGroupKind",
        }),
    }
}

fn write_legacy_trust(w: &mut PickleWriter) {
    // COS-native raw-v2/v3 compatibility byte. Current canonical topology has
    // no trust classification; structurally validated graph writes historical1.
    w.write_u8(1);
}

fn read_legacy_trust(r: &mut PickleReader) -> Result<(), PickleError> {
    match r.read_u8()? {
        0..=2 => Ok(()),
        value => Err(PickleError::InvalidEnumValue {
            value,
            type_name: "legacy TopologyTrust",
        }),
    }
}

fn write_stereo_group(w: &mut PickleWriter, sg: &StereoGroup) {
    // id
    if let Some(id) = sg.id() {
        w.write_bool(true);
        w.write_u32(id);
    } else {
        w.write_bool(false);
    }

    // kind
    write_stereo_group_kind(w, sg.kind());

    // atoms
    let atoms = sg.atoms();
    w.write_count(atoms.len());
    for &a in atoms {
        w.write_count(a.index());
    }

    // bonds
    let bonds = sg.bonds();
    w.write_count(bonds.len());
    for &b in bonds {
        w.write_count(b.index());
    }
}

// ──────────────────────────────────────────────
// Main serialization / deserialization
// ──────────────────────────────────────────────

/// Serialize a `Molecule` to a compact binary format.
///
/// Format structure:
/// - Versioned archive header and manifest
/// - Molecule state: atoms, bonds, coordinates, groups, and properties
/// - Derived chemistry state: rings, valence, aromaticity, and stereo validity
///
/// Archive 2.0 stores the complete molecule once, using Müsli storage, with a
/// separate derived-state block. Legacy raw 1..3 and archives 1.0..1.2 remain
/// readable through `decode_molecule_binary`.
///
/// # Errors
///
/// Returns `PickleError` if serialization encounters an internal issue.
pub fn encode_molecule_binary(mol: &BinaryInput<'_>) -> Result<Vec<u8>, PickleError> {
    mol.validate()?;
    archive_v2::encode(mol)
}

#[cfg(test)]
fn mol_to_legacy_binary(mol: &BinaryInput<'_>) -> Result<Vec<u8>, PickleError> {
    mol_to_legacy_binary_version(mol, 3)
}

fn mol_to_legacy_binary_version(
    mol: &BinaryInput<'_>,
    version: u8,
) -> Result<Vec<u8>, PickleError> {
    mol_to_legacy_binary_version_with_provenance(mol, version, None)
}
fn mol_to_legacy_binary_version_with_provenance(
    mol: &BinaryInput<'_>,
    version: u8,
    provenance: Option<&[LegacyStoreProvenance]>,
) -> Result<Vec<u8>, PickleError> {
    let mut w = PickleWriter::new();

    // Version
    w.write_u8(version);

    // ── Atoms ──
    let atoms = mol.atoms();
    if atoms.len() > u32::MAX as usize {
        return Err(PickleError::TooManyAtoms(atoms.len()));
    }
    w.write_count(atoms.len());
    for (index, atom) in atoms.iter().enumerate() {
        w.legacy_current = provenance.map(|p| p[index].clone());
        write_atom(&mut w, atom, version);
    }

    // ── Bonds ──
    let bonds = mol.bonds();
    if bonds.len() > u32::MAX as usize {
        return Err(PickleError::TooManyBonds(bonds.len()));
    }
    w.write_count(bonds.len());
    for (index, bond) in bonds.iter().enumerate() {
        w.legacy_current = provenance.map(|p| p[atoms.len() + index].clone());
        write_bond(&mut w, bond, version);
    }

    w.legacy_current = None;
    // ── 2D Coordinates ──
    if let Some(coords_2d) = mol.coordinates_2d() {
        w.write_bool(true);
        if coords_2d.len() > u32::MAX as usize {
            return Err(PickleError::TooManyAtoms(coords_2d.len()));
        }
        w.write_count(coords_2d.len());
        for &[x, y] in coords_2d {
            w.write_f64(x);
            w.write_f64(y);
        }
    } else {
        w.write_bool(false);
    }

    // ── 3D Conformers ──
    let conformers = mol.conformers_3d();
    w.write_count(conformers.len());
    for (ordinal, conf) in conformers.iter().enumerate() {
        w.write_count(ordinal);
        let coords = conf.coordinates();
        w.write_count(coords.len());
        for &[x, y, z] in coords {
            w.write_f64(x);
            w.write_f64(y);
            w.write_f64(z);
        }
        w.write_bool(conf.is_3d());
        w.write_props(conf.props());
    }

    // ── Source coordinate dimension ──
    match mol.source_coordinate_dim() {
        None => w.write_u8(0),
        Some(CoordinateDimension::TwoD) => w.write_u8(1),
        Some(CoordinateDimension::ThreeD) => w.write_u8(2),
    }

    // ── COSMolKit semantic capabilities ──
    if version >= 2 {
        write_legacy_trust(&mut w);
    }

    // ── Substance Groups ──
    let sgroups = mol.substance_groups();
    w.write_count(sgroups.len());
    for sg in sgroups {
        write_substance_group(&mut w, sg);
    }

    // ── Stereo Groups ──
    let stereo_groups = mol.stereo_groups();
    w.write_count(stereo_groups.len());
    for sg in stereo_groups {
        write_stereo_group(&mut w, sg);
    }

    // ── Molecule Properties ──
    let props = mol.properties();
    w.legacy_current = provenance.map(|p| p[atoms.len() + bonds.len()].clone());
    w.write_option_string(props.name());
    w.write_props(props.props());
    if version >= 3 {
        w.write_computed_props(props.computed_prop_names());
    }

    w.legacy_current = None;
    // SDF data fields
    let sdf_fields = props.sdf_data_fields();
    w.write_count(sdf_fields.len());
    for (key, value) in sdf_fields {
        w.write_string(key);
        w.write_string(value);
    }

    // SDF property lists
    let sdf_prop_lists = props.sdf_property_lists();
    w.write_count(sdf_prop_lists.len());
    for plist in sdf_prop_lists {
        match plist.target() {
            SdfPropertyListTarget::Atom => w.write_u8(0),
            SdfPropertyListTarget::Bond => w.write_u8(1),
        }
        w.write_string(plist.name());
        let values = plist.values();
        w.write_count(values.len());
        for v in values {
            match v {
                Some(value) => {
                    w.write_bool(true);
                    w.write_string(&property_value_to_string(value).map_err(invalid)?);
                }
                None => w.write_bool(false),
            }
        }
    }

    // ── Finish ──
    w.into_inner()
}

/// Deserialize a `Molecule` from binary data produced by `mol_to_binary`.
///
/// # Errors
///
/// Returns `PickleError` if the data is corrupt, has an unsupported version,
/// or produces an invalid molecule state.
pub fn decode_molecule_binary(data: &[u8]) -> Result<BinaryRecord, PickleError> {
    if data.first() == Some(&4) {
        let record = native_state_v2::decode(&data[1..])?;
        if native_state_v2::encode_record(&record)?.as_slice() != &data[1..] {
            return Err(PickleError::InvalidArchive(
                "raw4 disagrees with validated native state".into(),
            ));
        }
        return Ok(record);
    }
    if data.starts_with(archive_v2::MAGIC) {
        archive_v2::decode(data)
    } else if data.starts_with(ARCHIVE_MAGIC) {
        decode_sectioned_archive(data)
    } else {
        mol_from_legacy_binary(data)
    }
}

fn mol_from_legacy_binary(data: &[u8]) -> Result<BinaryRecord, PickleError> {
    Ok(mol_from_legacy_binary_with_provenance(data)?.0)
}
fn mol_from_legacy_binary_with_provenance(
    data: &[u8],
) -> Result<(BinaryRecord, Vec<LegacyStoreProvenance>), PickleError> {
    let mut provenance = Vec::new();
    let mut r = PickleReader::new(data);

    // Version check
    let version = r.read_u8()?;
    if !(1..=3).contains(&version) {
        return Err(PickleError::UnsupportedVersion(version));
    }

    // ── Atoms ──
    let atom_count = r.read_count(1)?;
    if atom_count > 1_000_000 {
        return Err(PickleError::TooManyAtoms(atom_count));
    }

    let mut atom_specs = Vec::with_capacity(atom_count);
    for i in 0..atom_count {
        let atomic_number = r.read_u8()?;
        let formal_charge = r.read_i8()?;
        let isotope = if r.read_bool()? {
            Some(u16::try_from(r.read_u32()?).map_err(invalid)?)
        } else {
            None
        };
        let chiral_tag = read_chiral_tag(&mut r)?;
        let chiral_perm = if r.read_bool()? {
            Some(r.read_u32()?)
        } else {
            None
        };
        let unknown_stereo = r.read_bool()?;
        let mol_parity = if r.read_bool()? {
            Some(r.read_i32()?)
        } else {
            None
        };
        let mol_inv_flag = if r.read_bool()? {
            Some(r.read_i32()?)
        } else {
            None
        };
        let radical_electrons = r.read_u8()?;
        let is_aromatic = r.read_bool()?;
        let hybridization = read_hybridization(&mut r)?;
        let atom_map = if r.read_bool()? {
            Some(r.read_u32()?)
        } else {
            None
        };
        let no_implicit = r.read_bool()?;
        let implicit_hydrogen = r.read_bool()?;
        let explicit_hydrogens = r.read_u8()?;
        let tracked_isotope_count = r.read_count(1)?;
        let mut tracked_isotopic_hydrogens = Vec::with_capacity(tracked_isotope_count);
        for _ in 0..tracked_isotope_count {
            tracked_isotopic_hydrogens.push(u16::try_from(r.read_u32()?).map_err(invalid)?);
        }
        if r.read_bool()? {
            return Err(PickleError::InvalidMolecule(
                "query state in concrete molecule payload".into(),
            ));
        }
        let props = r.read_props()?;
        let computed_props = if version >= 3 {
            r.read_computed_props()?
        } else {
            BTreeSet::new()
        };
        let _has_pdb_info = r.read_bool()?;

        let element =
            Element::from_atomic_number(atomic_number).ok_or(PickleError::InvalidEnumValue {
                value: atomic_number,
                type_name: "Element",
            })?;

        let mut spec = AtomSpec::new(element)
            .with_formal_charge(formal_charge)
            .with_chiral_tag(chiral_tag)
            .with_radical_electrons(radical_electrons)
            .with_aromatic(is_aromatic)
            .with_hybridization(hybridization)
            .with_no_implicit(no_implicit)
            .with_implicit_hydrogen(implicit_hydrogen)
            .with_explicit_hydrogens(explicit_hydrogens)
            .with_unknown_stereo(unknown_stereo);

        if let Some(iso) = isotope {
            spec = spec.with_isotope(iso);
        }
        if let Some(perm) = chiral_perm {
            spec = spec.with_chiral_permutation(perm);
        }
        if let Some(map) = atom_map {
            spec = spec.with_atom_map(map);
        }
        if let Some(parity) = mol_parity {
            spec = spec.with_mol_parity(parity);
        }
        if let Some(inv) = mol_inv_flag {
            spec = spec.with_mol_inversion_flag(inv);
        }
        if !tracked_isotopic_hydrogens.is_empty() {
            spec = spec.with_tracked_isotopic_hydrogens(tracked_isotopic_hydrogens);
        }

        let (rows, original) =
            legacy_map_rows(&props, &computed_props, format!("raw{version} atom {i}"))?;
        for (key, value) in rows {
            spec = spec.with_prop(key, value).map_err(invalid)?;
        }
        provenance.push(original);
        atom_specs.push((i, spec));
    }

    // ── Bonds ──
    let bond_count = r.read_count(1)?;
    if bond_count > 1_000_000 {
        return Err(PickleError::TooManyBonds(bond_count));
    }

    let mut bonds = Vec::with_capacity(bond_count);
    for index in 0..bond_count {
        bonds.push(read_bond(
            &mut r,
            version,
            BondId::new(index),
            &mut provenance,
        )?);
    }

    // ── 2D Coordinates ──
    let mut coords_2d: Option<Vec<[f64; 2]>> = None;
    if r.read_bool()? {
        let coord_count = r.read_count(1)?;
        let mut coords = Vec::with_capacity(coord_count);
        for _ in 0..coord_count {
            let x = r.read_f64()?;
            let y = r.read_f64()?;
            coords.push([x, y]);
        }
        coords_2d = Some(coords);
    }

    // ── 3D Conformers ──
    let conformer_count = r.read_count(1)?;
    let mut conformers_3d = Vec::with_capacity(conformer_count);
    for _ in 0..conformer_count {
        let conf_id = r.read_u32()? as usize;
        let coord_count = r.read_count(1)?;
        let mut coords = Vec::with_capacity(coord_count);
        for _ in 0..coord_count {
            let x = r.read_f64()?;
            let y = r.read_f64()?;
            let z = r.read_f64()?;
            coords.push([x, y, z]);
        }
        let is_3d = r.read_bool()?;
        let props = r.read_props()?;

        let mut conformer = Conformer3D::new(conf_id, coords, is_3d);
        // Replay props using the builder pattern — Conformer3D doesn't expose set_prop directly
        // but has with_prop (consumes self). We must reconstruct.
        // Actually Conformer3D::new creates with empty props, so we use the with_prop pattern
        // but since with_prop takes self, we need to re-architect for the loop.
        // Use crate::Conformer3D constructor then... actually let's just make a new one each time.
        for (k, v) in &props {
            conformer = conformer.with_prop(k.clone(), v.clone());
        }
        conformers_3d.push(conformer);
    }

    // ── Source coordinate dimension ──
    let source_coordinate_dim = match r.read_u8()? {
        0 => None,
        1 => Some(CoordinateDimension::TwoD),
        2 => Some(CoordinateDimension::ThreeD),
        value => {
            return Err(PickleError::InvalidEnumValue {
                value,
                type_name: "CoordinateDimension",
            });
        }
    };

    if version >= 2 {
        read_legacy_trust(&mut r)?;
    }

    // ── Substance Groups ──
    let sgroup_count = r.read_count(1)?;
    let mut sgroups = Vec::with_capacity(sgroup_count);
    for _ in 0..sgroup_count {
        let id = SubstanceGroupId::new(r.read_u32()? as usize);
        let has_rdkit_seq = r.read_bool()?;
        let rdkit_seq = if has_rdkit_seq {
            Some(r.read_u32()?)
        } else {
            None
        };
        let has_ext_id = r.read_bool()?;
        let ext_id = if has_ext_id {
            Some(r.read_u32()?)
        } else {
            None
        };
        let kind = read_substance_group_kind(&mut r)?;
        let atom_count_sg = r.read_count(1)?;
        let mut atoms_sg = Vec::with_capacity(atom_count_sg);
        for _ in 0..atom_count_sg {
            atoms_sg.push(AtomId::new(r.read_u32()? as usize));
        }
        let bond_count_sg = r.read_count(1)?;
        let mut bonds_sg = Vec::with_capacity(bond_count_sg);
        for _ in 0..bond_count_sg {
            bonds_sg.push(BondId::new(r.read_u32()? as usize));
        }
        let role_count = r.read_count(1)?;
        let mut bond_roles = BTreeMap::new();
        for _ in 0..role_count {
            let b = BondId::new(r.read_u32()? as usize);
            bond_roles.insert(b, SGroupBondRole::Contained);
        }
        let parent_atom_count = r.read_count(1)?;
        let mut parent_atoms = Vec::with_capacity(parent_atom_count);
        for _ in 0..parent_atom_count {
            parent_atoms.push(AtomId::new(r.read_u32()? as usize));
        }
        let has_parent = r.read_bool()?;
        let parent = if has_parent {
            Some(SubstanceGroupId::new(r.read_u32()? as usize))
        } else {
            None
        };
        let label: Option<PropertyText> = r.read_option_string()?.map(Into::into);
        let connection = read_sgroup_connection(&mut r)?;
        let subtype = r.read_option_string()?;
        let bracket_style = read_sgroup_bracket_style(&mut r)?;
        let expansion_state = r.read_option_string()?;
        let class = r.read_option_string()?;
        let has_component_number = r.read_bool()?;
        let component_number = if has_component_number {
            Some(r.read_u32()?)
        } else {
            None
        };
        let display = read_sgroup_display(&mut r)?;
        let data = read_sgroup_data(&mut r)?;

        // Attach points
        let ap_count = r.read_count(1)?;
        let mut attach_points = Vec::with_capacity(ap_count);
        for _ in 0..ap_count {
            let ap_atom = AtomId::new(r.read_u32()? as usize);
            let leaving = if r.read_bool()? {
                Some(AtomId::new(r.read_u32()? as usize))
            } else {
                None
            };
            let ap_label = r.read_option_string()?;
            let has_ap_order = r.read_bool()?;
            let ap_order = if has_ap_order {
                Some(r.read_u32()?)
            } else {
                None
            };
            attach_points.push(SGroupAttachPoint {
                atom: ap_atom,
                leaving_atom: leaving,
                label: ap_label.map(Into::into),
                order: ap_order,
            });
        }

        // CStates
        let cs_count = r.read_count(1)?;
        let mut cstates = Vec::with_capacity(cs_count);
        for _ in 0..cs_count {
            let cs_bond = BondId::new(r.read_u32()? as usize);
            let cs_x = r.read_f64()?;
            let cs_y = r.read_f64()?;
            cstates.push(SGroupCState {
                bond: cs_bond,
                vector: [cs_x, cs_y, 0.0],
            });
        }

        let props = r.read_props()?;
        let data_field_count = r.read_count(1)?;
        let mut data_fields = Vec::with_capacity(data_field_count);
        for _ in 0..data_field_count {
            data_fields.push(r.read_string()?);
        }

        let mut sg = SubstanceGroup::new(id, kind)
            .with_atoms(atoms_sg)
            .with_bonds(bonds_sg)
            .with_parent_atoms(parent_atoms)
            .with_attach_points(attach_points)
            .with_cstates(cstates);

        if let Some(seq) = rdkit_seq {
            sg = sg.with_rdkit_sequence_id(seq);
        }
        if let Some(eid) = ext_id {
            sg = sg.with_external_id(eid);
        }
        if let Some(p) = parent {
            sg = sg.with_parent(p);
        }
        if let Some(l) = label {
            sg = sg.with_label(l);
        }
        if let Some(conn) = connection {
            sg = sg.with_connection(conn);
        }
        if let Some(st) = subtype {
            sg = sg.with_subtype(st);
        }
        if let Some(bs) = bracket_style {
            sg = sg.with_bracket_style(bs);
        }
        if let Some(disp) = display {
            sg = sg.with_display(disp);
        }
        if let Some(es) = expansion_state {
            sg = sg.with_expansion_state(es);
        }
        if let Some(c) = class {
            sg = sg.with_class(c);
        }
        if let Some(cn) = component_number {
            sg = sg.with_component_number(cn);
        }
        if let Some(d) = data {
            sg = sg.with_data(d);
        }
        for (key, value) in &props {
            sg = sg.with_prop(key.clone(), value.clone()).map_err(invalid)?;
        }
        for df in &data_fields {
            sg = sg.with_data_field(df.clone());
        }
        // Apply bond roles (need to call after bonds are set)
        for (bond, role) in &bond_roles {
            if *role == SGroupBondRole::Contained {
                sg = sg.with_bond_role(*bond, SGroupBondRole::Contained);
            }
        }

        sgroups.push(sg);
    }

    // ── Stereo Groups ──
    let stereo_group_count = r.read_count(1)?;
    let mut stereo_groups = Vec::with_capacity(stereo_group_count);
    for _ in 0..stereo_group_count {
        let has_id = r.read_bool()?;
        let sg_id = if has_id { Some(r.read_u32()?) } else { None };
        let kind = read_stereo_group_kind(&mut r)?;
        let atom_count_sg = r.read_count(1)?;
        let mut atoms_sg = Vec::with_capacity(atom_count_sg);
        for _ in 0..atom_count_sg {
            atoms_sg.push(AtomId::new(r.read_u32()? as usize));
        }
        let bond_count_sg = r.read_count(1)?;
        let mut bonds_sg = Vec::with_capacity(bond_count_sg);
        for _ in 0..bond_count_sg {
            bonds_sg.push(BondId::new(r.read_u32()? as usize));
        }
        let mut sg = StereoGroup::new(kind, atoms_sg, bonds_sg)?;
        if let Some(id) = sg_id {
            sg = sg.with_id(id);
        }
        stereo_groups.push(sg);
    }

    // ── Molecule Properties ──
    let prop_name = r.read_option_string()?;
    let props = r.read_props()?;
    let computed_props = if version >= 3 {
        r.read_computed_props()?
    } else {
        BTreeSet::new()
    };
    let sdf_field_count = r.read_count(1)?;
    let mut sdf_data_fields = Vec::with_capacity(sdf_field_count);
    for _ in 0..sdf_field_count {
        let key = r.read_string()?;
        let value = r.read_string()?;
        sdf_data_fields.push((key, value));
    }
    let sdf_plist_count = r.read_count(1)?;
    let mut sdf_property_lists = Vec::with_capacity(sdf_plist_count);
    for _ in 0..sdf_plist_count {
        let target = match r.read_u8()? {
            0 => SdfPropertyListTarget::Atom,
            1 => SdfPropertyListTarget::Bond,
            value => {
                return Err(PickleError::InvalidEnumValue {
                    value,
                    type_name: "SdfPropertyListTarget",
                });
            }
        };
        let name = r.read_string()?;
        let val_count = r.read_count(1)?;
        let mut values = Vec::with_capacity(val_count);
        for _ in 0..val_count {
            values.push(
                r.read_option_string()?
                    .map(|value| PropertyValue::String(value.into())),
            );
        }
        sdf_property_lists.push(SdfPropertyList::new(target, name, values));
    }

    let atoms = atom_specs
        .into_iter()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect();
    let topology =
        TopologyBlock::try_from_parts(atoms, bonds, sgroups, stereo_groups).map_err(invalid)?;
    let coordinates = CoordinateBlock {
        conformers_2d: coords_2d
            .into_iter()
            .map(|rows| Conformer2D::new(0, rows))
            .collect(),
        conformers_3d,
        source_coordinate_dim,
        source_conformer_order: None,
    };

    validate_computed(&props, &computed_props)?;
    // Molecule properties — construct complete MoleculeProperties and set it
    let mut mol_props = MoleculeProperties::default();
    if let Some(name) = &prop_name {
        mol_props = mol_props.with_name(name.clone());
    }
    let (rows, original) =
        legacy_map_rows(&props, &computed_props, format!("raw{version} molecule"))?;
    for (key, value) in rows {
        mol_props = mol_props.with_prop(key, value).map_err(invalid)?;
    }
    provenance.push(original);
    for (key, value) in &sdf_data_fields {
        mol_props = mol_props.with_sdf_data_field(key.clone(), value.clone());
    }
    for plist in &sdf_property_lists {
        mol_props = mol_props.with_sdf_property_list(plist.clone());
    }
    if r.remaining() != 0 {
        return Err(PickleError::InvalidArchive(
            "trailing raw molecule bytes".into(),
        ));
    }
    let result = BinaryRecord {
        topology,
        coordinates,
        properties: mol_props,
        derived: BinaryDerivedState::default(),
    };
    result.validate()?;
    Ok((result, provenance))
}

// ──────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────

/// Explicit borrowed input to the single detached codec. No runtime or commit
/// authority crosses the IO boundary.
pub struct BinaryInput<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub properties: &'a MoleculeProperties,
    pub derived: BinaryDerivedView<'a>,
}

#[derive(Default)]
pub struct BinaryDerivedView<'a> {
    pub rings: Option<&'a RingInfo>,
    pub ring_families: Option<&'a RingInfo>,
    pub valence: Option<&'a ValenceAssignment>,
    pub aromaticity_valid: bool,
    pub stereo_valid: bool,
    pub valid_bits: u16,
}

#[derive(Debug, Default)]
pub struct BinaryDerivedState {
    pub rings: Option<RingInfo>,
    pub ring_families: Option<RingInfo>,
    pub valence: Option<ValenceAssignment>,
    pub aromaticity_valid: bool,
    pub stereo_valid: bool,
    /// Historical archives encode only presence and two validity flags.
    /// None keeps that distinction explicit for the private runtime.
    pub valid_bits: Option<u16>,
}

#[derive(Debug)]
pub struct BinaryRecord {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    pub derived: BinaryDerivedState,
}

impl BinaryInput<'_> {
    fn atoms(&self) -> &[Atom] {
        &self.topology.atoms
    }
    fn bonds(&self) -> &[Bond] {
        &self.topology.bonds
    }
    fn num_atoms(&self) -> usize {
        self.topology.atoms.len()
    }
    fn coordinates_2d(&self) -> Option<&[[f64; 2]]> {
        self.coordinates
            .conformers_2d
            .first()
            .map(Conformer2D::coordinates)
    }
    fn conformers_3d(&self) -> &[Conformer3D] {
        &self.coordinates.conformers_3d
    }
    fn source_coordinate_dim(&self) -> Option<CoordinateDimension> {
        self.coordinates.source_coordinate_dim
    }
    fn substance_groups(&self) -> &[SubstanceGroup] {
        &self.topology.substance_groups
    }
    fn stereo_groups(&self) -> &[StereoGroup] {
        &self.topology.stereo_groups
    }
    fn properties(&self) -> &MoleculeProperties {
        self.properties
    }
    fn validate(&self) -> Result<(), PickleError> {
        self.topology.validate().map_err(invalid)?;
        self.coordinates
            .validate_for_atom_count(self.num_atoms())
            .map_err(invalid)?;
        if self.derived.valid_bits & !0xff != 0 {
            return Err(PickleError::InvalidArchive(
                "unknown native derived-validity bits".into(),
            ));
        }
        if self.derived.aromaticity_valid != (self.derived.valid_bits & 8 != 0)
            || self.derived.stereo_valid != (self.derived.valid_bits & 16 != 0)
        {
            return Err(PickleError::InvalidArchive(
                "native and legacy validity flags disagree".into(),
            ));
        }
        Ok(())
    }
}

impl BinaryRecord {
    fn num_atoms(&self) -> usize {
        self.topology.atoms.len()
    }
    fn num_bonds(&self) -> usize {
        self.topology.bonds.len()
    }
    fn validate(&self) -> Result<(), PickleError> {
        self.topology.validate().map_err(invalid)?;
        self.coordinates
            .validate_for_atom_count(self.num_atoms())
            .map_err(invalid)?;
        if let Some(bits) = self.derived.valid_bits {
            if bits & !0xff != 0
                || self.derived.aromaticity_valid != (bits & 8 != 0)
                || self.derived.stereo_valid != (bits & 16 != 0)
            {
                return Err(PickleError::InvalidArchive(
                    "invalid native derived-validity metadata".into(),
                ));
            }
        }
        Ok(())
    }
}

fn invalid(error: impl std::fmt::Display) -> PickleError {
    PickleError::InvalidMolecule(error.to_string())
}

fn validate_computed<T>(
    props: &BTreeMap<String, T>,
    names: &BTreeSet<String>,
) -> Result<(), PickleError> {
    if names.iter().any(|key| !props.contains_key(key)) {
        return Err(PickleError::InvalidArchive(
            "computed property name without value".into(),
        ));
    }
    Ok(())
}

fn postcard_exact<'a, T: Deserialize<'a>>(data: &'a [u8]) -> Result<T, postcard::Error> {
    let (value, remainder) = postcard::take_from_bytes(data)?;
    if !remainder.is_empty() {
        return Err(postcard::Error::DeserializeBadEncoding);
    }
    Ok(value)
}

// Section4/v1 is an explicitly COS-native field extension. Its primitive
// tags and fixed-width float representation are part of this new format,
// not an inference about RDKit or historical string properties.
fn write_value(w: &mut PickleWriter, value: &PropertyValue) {
    match value {
        PropertyValue::String(value) => {
            w.write_u8(0);
            w.write_string(value);
        }
        PropertyValue::Int(value) => {
            w.write_u8(1);
            w.write_i32(*value);
        }
        PropertyValue::UInt(value) => {
            w.write_u8(2);
            w.write_u32(*value);
        }
        PropertyValue::IntVector(value) => {
            w.write_u8(3);
            w.write_count(value.len());
            for v in value {
                w.write_i32(*v);
            }
        }
        PropertyValue::Double(value) => {
            w.write_u8(4);
            w.write_f64(*value);
        }
        PropertyValue::Bool(value) => {
            w.write_u8(5);
            w.write_bool(*value);
        }
        PropertyValue::StringVector(_) => {
            w.error.get_or_insert(PickleError::InvalidArchive(
                "StringVector has no legacy canonical value tag".into(),
            ));
        }
    }
}

fn read_value(r: &mut PickleReader<'_>) -> Result<PropertyValue, PickleError> {
    match r.read_u8()? {
        0 => Ok(PropertyValue::String(r.read_string()?.into())),
        1 => Ok(PropertyValue::Int(r.read_i32()?)),
        2 => Ok(PropertyValue::UInt(r.read_u32()?)),
        3 => {
            let count = r.read_count(4)?;
            let mut values = Vec::with_capacity(count);
            for _ in 0..count {
                values.push(r.read_i32()?);
            }
            Ok(PropertyValue::IntVector(values))
        }
        4 => Ok(PropertyValue::Double(r.read_f64()?)),
        5 => Ok(PropertyValue::Bool(r.read_bool()?)),
        value => Err(PickleError::InvalidEnumValue {
            value,
            type_name: "PropertyValue",
        }),
    }
}

fn write_ordered_props<'a>(
    w: &mut PickleWriter,
    props: impl ExactSizeIterator<Item = (&'a PropertyText, &'a PropertyValue)>,
    computed: Result<Option<&[PropertyText]>, cosmolkit_model::PropertyValueError>,
) {
    let computed = match computed {
        Ok(v) => v.unwrap_or_default(),
        Err(e) => {
            w.error.get_or_insert(invalid(e));
            return;
        }
    };
    w.write_count(props.len());
    for (key, value) in props {
        w.write_string(key);
        write_value(w, value);
        w.write_bool(computed.contains(key));
    }
}

fn read_ordered_props(
    r: &mut PickleReader<'_>,
) -> Result<Vec<(String, PropertyValue, bool)>, PickleError> {
    let count = r.read_count(6)?;
    let mut seen = BTreeSet::new();
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let key = r.read_string()?;
        if key.is_empty() || !seen.insert(key.clone()) {
            return Err(PickleError::InvalidArchive(
                "empty or duplicate ordered property name".into(),
            ));
        }
        rows.push((key, read_value(r)?, r.read_bool()?));
    }
    Ok(rows)
}

fn write_pdb_info(w: &mut PickleWriter, info: Option<&AtomPdbResidueInfo>) {
    w.write_bool(info.is_some());
    if let Some(i) = info {
        w.write_string(i.atom_name());
        w.write_i32(i.serial_number());
        w.write_string(i.alt_loc());
        w.write_string(i.residue_name());
        w.write_i32(i.residue_number());
        w.write_string(i.chain_id());
        w.write_string(i.insertion_code());
        w.write_f64(i.occupancy());
        w.write_f64(i.temp_factor());
        w.write_bool(i.is_hetero_atom());
        w.write_u32(i.secondary_structure());
        w.write_u32(i.segment_number());
        w.write_string(i.monomer_class());
    }
}

fn read_pdb_info(r: &mut PickleReader<'_>) -> Result<Option<AtomPdbResidueInfo>, PickleError> {
    if !r.read_bool()? {
        return Ok(None);
    }
    let name = r.read_string()?;
    let serial = r.read_i32()?;
    let alt = r.read_string()?;
    let residue = r.read_string()?;
    let number = r.read_i32()?;
    let chain = r.read_string()?;
    let insertion = r.read_string()?;
    let occupancy = r.read_f64()?;
    let temperature = r.read_f64()?;
    let hetero = r.read_bool()?;
    let secondary = r.read_u32()?;
    let segment = r.read_u32()?;
    let class = r.read_string()?;
    Ok(Some(
        AtomPdbResidueInfo::new(name, serial, residue, number, chain, hetero)
            .with_alt_loc(alt)
            .with_insertion_code(insertion)
            .with_occupancy(occupancy)
            .with_temp_factor(temperature)
            .with_secondary_structure(secondary)
            .with_segment_number(segment)
            .with_monomer_class(class),
    ))
}

fn write_template_order(w: &mut PickleWriter, order: Option<&TemplateAttachmentOrder>) {
    w.write_bool(order.is_some());
    if let Some(order) = order {
        w.write_count(order.entries().len());
        for e in order.entries() {
            w.write_count(e.target().index());
            w.write_string(e.label());
        }
    }
}

fn read_template_order(
    r: &mut PickleReader<'_>,
) -> Result<Option<TemplateAttachmentOrder>, PickleError> {
    if !r.read_bool()? {
        return Ok(None);
    }
    let count = r.read_count(8)?;
    let mut entries = Vec::with_capacity(count);
    for _ in 0..count {
        entries.push(TemplateAttachment::new(
            AtomId::new(r.read_u32()? as usize),
            r.read_string()?,
        ));
    }
    Ok(Some(
        TemplateAttachmentOrder::new(entries).map_err(invalid)?,
    ))
}

fn read_native_id(r: &mut PickleReader<'_>) -> Result<usize, PickleError> {
    usize::try_from(r.read_u64()?).map_err(invalid)
}

#[cfg(test)]
fn encode_canonical_state(mol: &BinaryInput<'_>) -> Result<Vec<u8>, PickleError> {
    let mut w = PickleWriter::new();
    w.buf
        .extend_from_slice(&mol.derived.valid_bits.to_le_bytes());
    // Native canonical state preserves cache extents independently of topology:
    // AddHs can retain an authoritative cache containing only the original rows.
    // The historical derived-state section stays byte-for-byte unchanged.
    for rings in [mol.derived.rings, mol.derived.ring_families] {
        if let Some(rings) = rings {
            w.write_u32(checked_u32(rings.atom_row_count(), "ring atom extent")?);
            w.write_u32(checked_u32(rings.bond_row_count(), "ring bond extent")?);
        }
    }
    w.write_count(mol.atoms().len());
    for atom in mol.atoms() {
        write_ordered_props(
            &mut w,
            ordered_atom_properties(atom),
            atom.computed_prop_names(),
        );
        w.write_u64(atom.temporary_flags());
        write_pdb_info(&mut w, atom.pdb_residue_info());
        write_template_order(&mut w, atom.template_attachment_order());
    }
    w.write_count(mol.bonds().len());
    for bond in mol.bonds() {
        write_ordered_props(
            &mut w,
            ordered_bond_properties(bond),
            bond.computed_prop_names(),
        );
        w.write_u64(bond.temporary_flags());
    }
    w.write_count(mol.coordinates.conformers_2d.len());
    for conf in &mol.coordinates.conformers_2d {
        w.write_u64(u64::try_from(conf.id()).map_err(invalid)?);
        w.write_count(conf.coordinates().len());
        for row in conf.coordinates() {
            for x in row {
                w.write_f64(*x);
            }
        }
        w.write_props(conf.props());
    }
    w.write_count(mol.conformers_3d().len());
    for conf in mol.conformers_3d() {
        w.write_u64(u64::try_from(conf.id()).map_err(invalid)?);
        w.write_count(conf.coordinates().len());
        for row in conf.coordinates() {
            for x in row {
                w.write_f64(*x);
            }
        }
        w.write_bool(conf.is_3d());
        w.write_props(conf.props());
    }
    w.write_count(mol.substance_groups().len());
    for sg in mol.substance_groups() {
        w.write_bool(sg.display().is_some());
        if let Some(display) = sg.display() {
            w.write_count(display.brackets.len());
            for bracket in &display.brackets {
                for row in bracket.points {
                    for x in row {
                        w.write_f64(x);
                    }
                }
            }
        }
        w.write_count(sg.cstates().len());
        for cs in sg.cstates() {
            w.write_count(cs.bond.index());
            for x in cs.vector {
                w.write_f64(x);
            }
        }
        for list in [sg.head_crossing_bonds(), sg.crossing_bond_correspondence()] {
            w.write_count(list.len());
            for id in list {
                w.write_count(id.index());
            }
        }
    }
    w.write_count(mol.stereo_groups().len());
    for sg in mol.stereo_groups() {
        w.write_u32(sg.write_id());
    }
    w.write_count(mol.properties.sdf_property_lists().len());
    for list in mol.properties.sdf_property_lists() {
        w.write_u8(match list.target() {
            SdfPropertyListTarget::Atom => 0,
            SdfPropertyListTarget::Bond => 1,
        });
        w.write_string(list.name());
        w.write_count(list.values().len());
        for value in list.values() {
            w.write_bool(value.is_some());
            if let Some(value) = value {
                write_value(&mut w, value);
            }
        }
    }
    w.write_bool(mol.coordinates.source_conformer_order.is_some());
    if let Some(order) = &mol.coordinates.source_conformer_order {
        w.write_count(order.len());
        for dimension in order {
            w.write_u8(match dimension {
                CoordinateDimension::TwoD => 2,
                CoordinateDimension::ThreeD => 3,
            });
        }
    }
    w.into_inner()
}

fn exact_count(r: &mut PickleReader<'_>, expected: usize) -> Result<(), PickleError> {
    let actual = r.read_count(1)?;
    if actual != expected {
        return Err(PickleError::DataLengthMismatch { expected, actual });
    }
    Ok(())
}

fn decode_canonical_state(
    data: &[u8],
    version: u16,
    record: &mut BinaryRecord,
) -> Result<(), PickleError> {
    let mut r = PickleReader::new(data);
    let bits = read_u16_le(&mut r)?;
    record.derived.valid_bits = Some(bits);
    let atom_count = record.num_atoms();
    let bond_count = record.num_bonds();
    for slot in [&mut record.derived.rings, &mut record.derived.ring_families] {
        if let Some(rings) = slot.take() {
            let atom_rows = r.read_u32()? as usize;
            let bond_rows = r.read_u32()? as usize;
            if atom_rows > atom_count
                || bond_rows > bond_count
                || (!rings.is_initialized() && (atom_rows != 0 || bond_rows != 0))
            {
                return Err(PickleError::InvalidArchive(
                    "ring cache extents exceed topology or uninitialized state".into(),
                ));
            }
            *slot = Some(
                RingInfo::from_persisted_components(
                    rings.is_initialized(),
                    rings.persisted_find_type(),
                    atom_rows,
                    bond_rows,
                    rings.atom_rings().to_vec(),
                    rings.bond_rings().to_vec(),
                    rings.atom_ring_families().to_vec(),
                    rings.bond_ring_families().to_vec(),
                    rings.persisted_relevant_cycle_count(),
                    rings.persisted_fused_rings().to_vec(),
                    rings.persisted_num_fused_bonds().to_vec(),
                )
                .map_err(invalid)?,
            );
        }
    }
    exact_count(&mut r, record.num_atoms())?;
    for atom in &mut record.topology.atoms {
        let properties = legacy_ordered_rows(
            read_ordered_props(&mut r)?,
            format!("canonical{version} atom {}", atom.id().index()),
        )?;
        atom.replace_property_records(properties).map_err(invalid)?;
        atom.set_temporary_flags(r.read_u64()?);
        atom.set_pdb_residue_info(read_pdb_info(&mut r)?);
        replace_atom_template_attachment_order(atom, read_template_order(&mut r)?);
    }
    exact_count(&mut r, record.num_bonds())?;
    for bond in &mut record.topology.bonds {
        let properties = legacy_ordered_rows(
            read_ordered_props(&mut r)?,
            format!("canonical{version} bond {}", bond.id().index()),
        )?;
        bond.replace_property_records(properties).map_err(invalid)?;
        bond.set_temporary_flags(r.read_u64()?);
    }
    let count = r.read_count(12)?;
    let mut two_d = Vec::with_capacity(count);
    for _ in 0..count {
        let id = read_native_id(&mut r)?;
        let count = r.read_count(16)?;
        if count != record.num_atoms() {
            return Err(PickleError::DataLengthMismatch {
                expected: record.num_atoms(),
                actual: count,
            });
        }
        let mut rows = Vec::with_capacity(count);
        for _ in 0..count {
            rows.push([r.read_f64()?, r.read_f64()?]);
        }
        let mut conf = Conformer2D::new(id, rows);
        for (key, value) in r.read_props()? {
            conf = conf.with_prop(key, value);
        }
        two_d.push(conf);
    }
    let count = r.read_count(13)?;
    let mut three_d = Vec::with_capacity(count);
    for _ in 0..count {
        let id = read_native_id(&mut r)?;
        let count = r.read_count(24)?;
        if count != record.num_atoms() {
            return Err(PickleError::DataLengthMismatch {
                expected: record.num_atoms(),
                actual: count,
            });
        }
        let mut rows = Vec::with_capacity(count);
        for _ in 0..count {
            rows.push([r.read_f64()?, r.read_f64()?, r.read_f64()?]);
        }
        let mut conf = Conformer3D::new(id, rows, r.read_bool()?);
        for (key, value) in r.read_props()? {
            conf = conf.with_prop(key, value);
        }
        three_d.push(conf);
    }
    record.coordinates.conformers_2d = two_d;
    record.coordinates.conformers_3d = three_d;
    exact_count(&mut r, record.topology.substance_groups.len())?;
    let groups = std::mem::take(&mut record.topology.substance_groups);
    for mut sg in groups {
        let display = r.read_bool()?;
        if display != sg.display().is_some() {
            return Err(PickleError::InvalidArchive(
                "display presence mismatch".into(),
            ));
        }
        if display {
            let count = r.read_count(72)?;
            let mut brackets = Vec::with_capacity(count);
            for _ in 0..count {
                let mut points = [[0.0; 3]; 3];
                for row in &mut points {
                    for x in row {
                        *x = r.read_f64()?;
                    }
                }
                brackets.push(SGroupBracket::new(points));
            }
            sg.display_mut().brackets = brackets;
        }
        let count = r.read_count(28)?;
        let mut cstates = Vec::with_capacity(count);
        for _ in 0..count {
            cstates.push(SGroupCState::new(
                BondId::new(r.read_u32()? as usize),
                [r.read_f64()?, r.read_f64()?, r.read_f64()?],
            ));
        }
        sg = sg.with_cstates(cstates);
        let mut lists = Vec::with_capacity(2);
        for _ in 0..2 {
            let count = r.read_count(4)?;
            let mut values = Vec::with_capacity(count);
            for _ in 0..count {
                values.push(BondId::new(r.read_u32()? as usize));
            }
            lists.push(values);
        }
        let mut lists = lists.into_iter();
        sg = sg
            .with_head_crossing_bonds(lists.next().expect("first fixed list"))
            .with_crossing_bond_correspondence(lists.next().expect("second fixed list"));
        record.topology.substance_groups.push(sg);
    }
    exact_count(&mut r, record.topology.stereo_groups.len())?;
    let groups = std::mem::take(&mut record.topology.stereo_groups);
    for sg in groups {
        record
            .topology
            .stereo_groups
            .push(sg.with_write_id(r.read_u32()?));
    }
    exact_count(&mut r, record.properties.sdf_property_lists().len())?;
    let mut props = MoleculeProperties::default();
    if let Some(name) = record.properties.name() {
        props = props.with_name(name);
    }
    for (key, value) in record.properties.ordered_props() {
        props = props.with_prop(key, value).map_err(invalid)?;
    }
    for (key, value) in record.properties.sdf_data_fields() {
        props = props.with_sdf_data_field(key, value);
    }
    for old in record.properties.sdf_property_lists() {
        let target = match r.read_u8()? {
            0 => SdfPropertyListTarget::Atom,
            1 => SdfPropertyListTarget::Bond,
            value => {
                return Err(PickleError::InvalidEnumValue {
                    value,
                    type_name: "SdfPropertyListTarget",
                });
            }
        };
        let name = r.read_string()?;
        let count = r.read_count(1)?;
        if target != old.target()
            || name.as_bytes() != old.name().as_bytes()
            || count != old.values().len()
        {
            return Err(PickleError::InvalidArchive(
                "SDF property list carrier mismatch".into(),
            ));
        }
        let mut values = Vec::with_capacity(count);
        for _ in 0..count {
            values.push(if r.read_bool()? {
                Some(read_value(&mut r)?)
            } else {
                None
            });
        }
        props = props.with_sdf_property_list(SdfPropertyList::new(target, name, values));
    }
    record.properties = props;
    record.coordinates.source_conformer_order = if version >= 2 && r.read_bool()? {
        let count = r.read_count(1)?;
        let mut order = Vec::with_capacity(count);
        for _ in 0..count {
            order.push(match r.read_u8()? {
                2 => CoordinateDimension::TwoD,
                3 => CoordinateDimension::ThreeD,
                value => {
                    return Err(PickleError::InvalidEnumValue {
                        value,
                        type_name: "CoordinateDimension",
                    });
                }
            });
        }
        Some(order)
    } else {
        // Earlier CK canonical archives did not encode interleaving. Do not
        // invent a first conformer when importing that unrecoverable state.
        None
    };
    if r.remaining() != 0 {
        return Err(PickleError::InvalidArchive(
            "trailing canonical state bytes".into(),
        ));
    }
    record.validate()
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{
        AtomSpec, BondOrder, BondSpec, BondStereo, ChiralTag, Element, Hybridization,
        SdfPropertyList, SdfPropertyListTarget, StereoGroup, StereoGroupKind,
    };

    #[test]
    fn source_order_binary_roundtrip_keeps_cross_dimension_appends() {
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        let mut record = builder.build().unwrap();
        record
            .coordinates
            .record_source_conformer_append(CoordinateDimension::ThreeD)
            .unwrap();
        record
            .coordinates
            .conformers_3d
            .push(Conformer3D::new(9, vec![[1.0, 2.0, -0.0]], false));
        record
            .coordinates
            .record_source_conformer_append(CoordinateDimension::TwoD)
            .unwrap();
        record
            .coordinates
            .conformers_2d
            .push(Conformer2D::new(9, vec![[3.0, 4.0]]));
        record
            .coordinates
            .record_source_conformer_append(CoordinateDimension::ThreeD)
            .unwrap();
        record
            .coordinates
            .conformers_3d
            .push(Conformer3D::new(0, vec![[5.0, 6.0, 7.0]], true));
        let bytes = fixture_encode(&record).unwrap();
        let restored = decode_molecule_binary(&bytes).unwrap();
        assert_record_equal(&record, &restored, "source order transport");
        match restored
            .coordinates
            .first_source_conformer()
            .unwrap()
            .unwrap()
        {
            cosmolkit_model::CoordinateSourceConformer::ThreeD(conformer) => {
                assert_eq!(conformer.id(), 9);
                assert!(!conformer.is_3d());
                assert_eq!(
                    conformer.coordinates()[0][2].to_bits(),
                    (-0.0_f64).to_bits()
                );
            }
            _ => panic!("first actual append remains first after decoding"),
        }
        assert_eq!(fixture_encode(&restored).unwrap(), bytes);
    }

    #[test]
    fn native12_preserves_authoritative_sparse_ring_cache_extents() {
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        for atom in 1..5 {
            builder.add_atom(AtomSpec::new(Element::H));
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(atom),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        let mut original = builder.build().unwrap();
        original.derived.rings = Some(RingInfo::new(RingFindType::SymmSssr, 1, 0));
        original.derived.ring_families = Some(RingInfo::new(RingFindType::Fast, 0, 0));
        original.derived.valid_bits = Some(3);
        let bytes = fixture_encode(&original).unwrap();
        let restored = decode_molecule_binary(&bytes).unwrap();
        assert_record_equal(&original, &restored, "authoritative sparse ring extent");
        assert_eq!(restored.derived.rings.as_ref().unwrap().atom_row_count(), 1);
        assert_eq!(restored.derived.rings.as_ref().unwrap().bond_row_count(), 0);
        assert_eq!(fixture_encode(&restored).unwrap(), bytes);
    }

    fn fixture_input(record: &BinaryRecord) -> BinaryInput<'_> {
        let d = &record.derived;
        let bits = d.valid_bits.unwrap_or(
            (if d.rings.is_some() { 1 } else { 0 })
                | (if d.ring_families.is_some() { 2 } else { 0 })
                | (if d.valence.is_some() { 4 } else { 0 })
                | (if d.aromaticity_valid { 8 } else { 0 })
                | (if d.stereo_valid { 16 } else { 0 }),
        );
        BinaryInput {
            topology: &record.topology,
            coordinates: &record.coordinates,
            properties: &record.properties,
            derived: BinaryDerivedView {
                rings: d.rings.as_ref(),
                ring_families: d.ring_families.as_ref(),
                valence: d.valence.as_ref(),
                aromaticity_valid: d.aromaticity_valid,
                stereo_valid: d.stereo_valid,
                valid_bits: bits,
            },
        }
    }
    fn fixture_encode(record: &BinaryRecord) -> Result<Vec<u8>, PickleError> {
        encode_molecule_binary(&fixture_input(record))
    }
    fn fixture_encode_legacy12(record: &BinaryRecord) -> Result<Vec<u8>, PickleError> {
        let input = fixture_input(record);
        let mut archive = ARCHIVE_MAGIC.to_vec();
        write_u16_le(&mut archive, 1);
        write_u16_le(&mut archive, 2);
        write_u16_le(&mut archive, 4);
        for (id, version, flags, codec, payload) in [
            (1, 1, 0, 1, encode_manifest()?),
            (
                2,
                1,
                1,
                1,
                encode_molecule_state(mol_to_legacy_binary(&input)?)?,
            ),
            (3, 1, 1, 0, encode_derived_state(&input)?),
            (4, 2, 1, 0, encode_canonical_state(&input)?),
        ] {
            write_archive_section(&mut archive, id, version, flags, codec, &payload)?;
        }
        Ok(archive)
    }

    #[test]
    fn frozen_legacy_raw_and_archive_layouts_remain_readable_without_writers() {
        // Fixed empty-molecule bytes transcribed from the legacy wire layout:
        // raw1/2/3; archive1.0/1.1; archive1.2 canonical1 and canonical2.
        // Producer is the fixed text "legacy", not the current package version.
        let fixtures = [
            "010000000000000000000000000000000000000000000000000000000000000000000000",
            "02000000000000000000000000000001000000000000000000000000000000000000000000",
            "0300000000000000000000000000000100000000000000000000000000000000000000000000000000",
            "43534d4f4c504b4c01000000020001000100000109000000066c656761637901010200010001012c0000000003290300000000000000000000000000000100000000000000000000000000000000000000000000000000",
            "43534d4f4c504b4c01000100030001000100000109000000066c656761637901010200010001012c0000000003290300000000000000000000000000000100000000000000000000000000000000000000000000000000030001000100050000000000000000",
            "43534d4f4c504b4c01000200040001000100000109000000066c656761637901010200010001012c00000000032903000000000000000000000000000001000000000000000000000000000000000000000000000000000300010001000500000000000000000400010001001e000000000000000000000000000000000000000000000000000000000000000000",
            "43534d4f4c504b4c01000200040001000100000109000000066c656761637901010200010001012c00000000032903000000000000000000000000000001000000000000000000000000000000000000000000000000000300010001000500000000000000000400020001001f00000000000000000000000000000000000000000000000000000000000000000000",
        ];
        for (i, hex) in fixtures.into_iter().enumerate() {
            let bytes = (0..hex.len())
                .step_by(2)
                .map(|j| u8::from_str_radix(&hex[j..j + 2], 16).unwrap())
                .collect::<Vec<_>>();
            let record = decode_molecule_binary(&bytes)
                .unwrap_or_else(|e| panic!("legacy fixture {i}: {e}"));
            assert_record_equal(&BinaryRecord::new(), &record, "fixed legacy wire");
            let upgraded = fixture_encode(&record).unwrap();
            assert_eq!(&upgraded[8..12], &[2, 0, 0, 0]);
            assert_record_equal(
                &record,
                &decode_molecule_binary(&upgraded).unwrap(),
                "legacy to archive 2",
            );
        }
    }

    #[test]
    fn legacy12_companions_use_their_actual_raw_version() {
        let record = BinaryRecord::new();
        let input = fixture_input(&record);
        for version in 1..=3 {
            let raw = mol_to_legacy_binary_version(&input, version).unwrap();
            let state = postcard::to_allocvec(&MoleculeStateV1 {
                encoding: 0,
                encoding_version: version,
                payload: raw,
            })
            .unwrap();
            let manifest = encode_manifest().unwrap();
            let derived = encode_derived_state(&input).unwrap();
            let canonical = encode_canonical_state(&input).unwrap();
            let mut archive = ARCHIVE_MAGIC.to_vec();
            write_u16_le(&mut archive, 1);
            write_u16_le(&mut archive, 2);
            write_u16_le(&mut archive, 4);
            for (id, v, flags, codec, payload) in [
                (1, 1, 0, 1, manifest),
                (2, 1, 1, 1, state),
                (3, 1, 1, 0, derived),
                (4, 2, 1, 0, canonical),
            ] {
                write_archive_section(&mut archive, id, v, flags, codec, &payload).unwrap();
            }
            assert_record_equal(
                &record,
                &decode_molecule_binary(&archive).unwrap(),
                "version-specific legacy companion",
            );
        }
    }

    #[test]
    fn legacy_computed_collision_imports_empty_and_nonconflicting_state_upgrades() {
        let mut original = build_simple_methane();
        original.properties = original.properties.with_prop("mass", "16.043").unwrap();
        original.topology.atoms[0].set_prop("rank", 7).unwrap();
        original.topology.atoms[0].set_pdb_residue_info(Some(AtomPdbResidueInfo::new(
            "CA", 12, "ALA", 4, "A", false,
        )));
        original
            .coordinates
            .conformers_2d
            .push(Conformer2D::new(9, vec![[-0.0, 2.5]; 5]));
        let bytes = fixture_encode_legacy12(&original).unwrap();
        let imported = decode_molecule_binary(&bytes).unwrap();
        assert_record_equal(&original, &imported, "legacy collision import");
        let current = fixture_encode(&imported).unwrap();
        let restored = decode_molecule_binary(&current).unwrap();
        assert_record_equal(&original, &restored, "legacy collision upgrade");
        assert_eq!(
            restored.properties.prop("mass"),
            Some(&PropertyValue::String("16.043".into()))
        );
        assert!(!restored.properties.is_prop_computed("mass").unwrap());

        // Explicit user policy: normalize only the conflicting reserved slot
        // to an empty list; retain the actual values of other properties.
        let old = BinaryRecord::new()
            .with_prop("__computedProps", "opaque")
            .with_prop("mass", "16.043");
        let provenance = [LegacyStoreProvenance {
            reserved: true,
            collision_reserved: None,
            computed: vec!["mass".into()],
            context: "old molecule".into(),
        }];
        let raw = mol_to_legacy_binary_version_with_provenance(
            &fixture_input(&old),
            3,
            Some(&provenance),
        )
        .unwrap();
        let mut expected = BinaryRecord::new();
        expected.properties = expected
            .properties
            .with_prop("__computedProps", PropertyValue::StringVector(vec![]))
            .unwrap()
            .with_prop("mass", "16.043")
            .unwrap();
        let imported = decode_molecule_binary(&raw).unwrap();
        assert_record_equal(&expected, &imported, "raw3 conflict imports empty");
        assert_eq!(
            imported.properties.computed_prop_names().unwrap(),
            Some(&[][..])
        );
        assert!(!imported.properties.is_prop_computed("mass").unwrap());
        let upgraded = fixture_encode(&imported).unwrap();
        assert_record_equal(
            &expected,
            &decode_molecule_binary(&upgraded).unwrap(),
            "empty conflict slot survives archive 2 upgrade",
        );

        // Archives 1.0/1.1/1.2 must retain actual legacy integrity checks;
        // normalization must not create a false companion mismatch.
        let input = fixture_input(&old);
        for minor in 0..=2 {
            let mut archive = ARCHIVE_MAGIC.to_vec();
            write_u16_le(&mut archive, 1);
            write_u16_le(&mut archive, minor);
            write_u16_le(
                &mut archive,
                2 + u16::from(minor >= 1) + u16::from(minor >= 2),
            );
            write_archive_section(&mut archive, 1, 1, 0, 1, &encode_manifest().unwrap()).unwrap();
            write_archive_section(
                &mut archive,
                2,
                1,
                1,
                1,
                &encode_molecule_state(raw.clone()).unwrap(),
            )
            .unwrap();
            if minor >= 1 {
                write_archive_section(
                    &mut archive,
                    3,
                    1,
                    1,
                    0,
                    &encode_derived_state(&input).unwrap(),
                )
                .unwrap();
            }
            if minor >= 2 {
                write_archive_section(
                    &mut archive,
                    4,
                    2,
                    1,
                    0,
                    &encode_canonical_state(&input).unwrap(),
                )
                .unwrap();
            }
            let imported = decode_molecule_binary(&archive).unwrap();
            assert_record_equal(
                &expected,
                &imported,
                "legacy archive conflict imports empty",
            );
        }
    }
    fn assert_record_equal(a: &BinaryRecord, b: &BinaryRecord, message: &str) {
        assert_eq!(a.topology, b.topology, "{message}: whole topology");
        assert_eq!(
            a.coordinates, b.coordinates,
            "{message}: whole coordinate state"
        );
        assert_eq!(
            a.properties, b.properties,
            "{message}: whole molecule properties"
        );
        assert_eq!(a.derived.rings, b.derived.rings, "{message}: rings");
        assert_eq!(
            a.derived.ring_families, b.derived.ring_families,
            "{message}: ring families"
        );
        assert_eq!(a.derived.valence, b.derived.valence, "{message}: valence");
        assert_eq!(
            a.derived.aromaticity_valid, b.derived.aromaticity_valid,
            "{message}: aromaticity validity"
        );
        assert_eq!(
            a.derived.stereo_valid, b.derived.stereo_valid,
            "{message}: stereo validity"
        );
        if let Some(bits) = b.derived.valid_bits {
            assert_eq!(
                fixture_input(a).derived.valid_bits,
                bits,
                "{message}: all native validity"
            );
        }
    }
    impl BinaryRecord {
        fn new() -> Self {
            Self {
                topology: TopologyBlock::default(),
                coordinates: CoordinateBlock::default(),
                properties: MoleculeProperties::default(),
                derived: BinaryDerivedState::default(),
            }
        }
        fn from_blocks(
            topology: TopologyBlock,
            coordinates: CoordinateBlock,
            properties: MoleculeProperties,
        ) -> Result<Self, PickleError> {
            let record = Self {
                topology,
                coordinates,
                properties,
                derived: BinaryDerivedState::default(),
            };
            record.validate()?;
            Ok(record)
        }
        fn atoms(&self) -> &[Atom] {
            &self.topology.atoms
        }
        fn bonds(&self) -> &[Bond] {
            &self.topology.bonds
        }
        fn coordinates_2d(&self) -> Option<&[[f64; 2]]> {
            self.coordinates
                .conformers_2d
                .first()
                .map(Conformer2D::coordinates)
        }
        fn conformers_3d(&self) -> &[Conformer3D] {
            &self.coordinates.conformers_3d
        }
        fn stereo_groups(&self) -> &[StereoGroup] {
            &self.topology.stereo_groups
        }
        fn properties(&self) -> &MoleculeProperties {
            &self.properties
        }
        fn prop(&self, key: &str) -> Option<&str> {
            self.properties.prop(key).map(|value| {
                super::fixture_text(value.as_string().expect("binary fixture String tag"))
            })
        }
        fn with_name(mut self, name: &str) -> Self {
            self.properties = self.properties.with_name(name);
            self
        }
        fn with_prop(mut self, key: &str, value: &str) -> Self {
            self.properties = self
                .properties
                .with_prop(key, value)
                .expect("valid original property");
            self
        }
        fn with_deserialized_derived_cache(
            mut self,
            derived: BinaryDerivedState,
        ) -> Result<Self, PickleError> {
            self.derived = derived;
            self.validate()?;
            Ok(self)
        }
    }
    /// Test fixture construction only; this is not a production builder or a
    /// live molecule. All values are canonical detached model records.
    struct DetachedFixtureBuilder(BinaryRecord);
    impl DetachedFixtureBuilder {
        fn new() -> Self {
            Self(BinaryRecord::new())
        }
        fn add_atom(&mut self, spec: AtomSpec) -> AtomId {
            let id = AtomId::new(self.0.topology.atoms.len());
            self.0.topology.atoms.push(Atom::from_spec(id, spec));
            id
        }
        fn add_bond(&mut self, spec: BondSpec) -> Result<BondId, PickleError> {
            let id = BondId::new(self.0.topology.bonds.len());
            self.0.topology.bonds.push(Bond::from_spec(id, spec));
            Ok(id)
        }
        fn set_2d_coordinates(&mut self, rows: Vec<[f64; 2]>) -> Result<(), PickleError> {
            self.0.coordinates.conformers_2d = vec![Conformer2D::new(0, rows)];
            Ok(())
        }
        fn add_3d_conformer(&mut self, rows: Vec<[f64; 3]>) -> Result<(), PickleError> {
            self.0.coordinates.conformers_3d.push(Conformer3D::new(
                self.0.coordinates.conformers_3d.len(),
                rows,
                true,
            ));
            Ok(())
        }
        fn add_stereo_group(&mut self, group: StereoGroup) -> Result<(), PickleError> {
            self.0.topology.stereo_groups.push(group);
            Ok(())
        }
        fn with_sdf_data_field(mut self, key: &str, value: &str) -> Self {
            self.0.properties = self.0.properties.with_sdf_data_field(key, value);
            self
        }
        fn build(mut self) -> Result<BinaryRecord, PickleError> {
            self.0.topology = TopologyBlock::try_from_parts(
                self.0.topology.atoms,
                self.0.topology.bonds,
                self.0.topology.substance_groups,
                self.0.topology.stereo_groups,
            )
            .map_err(invalid)?;
            self.0.validate()?;
            Ok(self.0)
        }
    }
    fn build_simple_methane() -> BinaryRecord {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        for _ in 0..4 {
            builder.add_atom(AtomSpec::new(h));
        }
        for i in 0..4 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(i + 1),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        builder.build().expect("build methane")
    }

    fn build_simple_ethanol() -> BinaryRecord {
        let c = Element::from_atomic_number(6).unwrap();
        let o = Element::from_atomic_number(8).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        // C-C-O-H backbone: indices 0,1,2,3
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(
            AtomSpec::new(c)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_chiral_permutation(0),
        );
        builder.add_atom(AtomSpec::new(o));
        builder.add_atom(AtomSpec::new(h));
        // 3 H on C0
        for _ in 0..3 {
            builder.add_atom(AtomSpec::new(h));
        }
        // 2 H on C1
        for _ in 0..2 {
            builder.add_atom(AtomSpec::new(h));
        }
        // 1 H on O
        builder.add_atom(AtomSpec::new(h));

        // C-C single
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();
        // C-O single
        builder
            .add_bond(BondSpec::new(
                AtomId::new(1),
                AtomId::new(2),
                BondOrder::Single,
            ))
            .unwrap();
        // O-H single
        builder
            .add_bond(BondSpec::new(
                AtomId::new(2),
                AtomId::new(3),
                BondOrder::Single,
            ))
            .unwrap();
        // 3 C-H on C0
        for i in 0..3 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(4 + i),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        // 2 C-H on C1
        for i in 0..2 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(1),
                    AtomId::new(7 + i),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        // O-H already on index 3
        builder
            .add_bond(BondSpec::new(
                AtomId::new(2),
                AtomId::new(9),
                BondOrder::Single,
            ))
            .unwrap();

        builder.build().expect("build ethanol")
    }

    fn encode_v1_0_sectioned_archive(mol: &BinaryRecord) -> Vec<u8> {
        let manifest = encode_manifest().expect("encode v1.0 manifest");
        let molecule_state = encode_molecule_state(
            mol_to_legacy_binary(&fixture_input(mol)).expect("encode v1.0 legacy molecule state"),
        )
        .expect("encode v1.0 molecule-state envelope");
        let mut data = Vec::new();
        data.extend_from_slice(ARCHIVE_MAGIC);
        write_u16_le(&mut data, ARCHIVE_MAJOR);
        write_u16_le(&mut data, 0);
        write_u16_le(&mut data, 2);
        write_archive_section(
            &mut data,
            SECTION_MANIFEST,
            MANIFEST_VERSION,
            0,
            SECTION_CODEC_POSTCARD,
            &manifest,
        )
        .expect("write v1.0 manifest section");
        write_archive_section(
            &mut data,
            SECTION_MOLECULE_STATE,
            MOLECULE_STATE_VERSION,
            SECTION_FLAG_REQUIRED,
            SECTION_CODEC_POSTCARD,
            &molecule_state,
        )
        .expect("write v1.0 molecule-state section");
        data
    }

    #[test]
    fn test_empty_molecule_roundtrip() {
        let mol = BinaryRecord::new();
        let data = fixture_encode(&mol).unwrap();
        assert!(data.starts_with(archive_v2::MAGIC));
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "empty molecule roundtrip failed");
    }

    #[test]
    fn test_legacy_molecule_state_roundtrip_remains_readable() {
        let mol = build_simple_methane();
        let data = mol_to_legacy_binary(&fixture_input(&mol)).unwrap();
        assert!(!data.starts_with(ARCHIVE_MAGIC));
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "legacy molecule-state roundtrip failed");
    }

    #[test]
    fn test_v1_0_sectioned_archive_remains_readable() {
        let mol = build_simple_methane();
        let data = encode_v1_0_sectioned_archive(&mol);
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "v1.0 sectioned archive decode failed");
    }

    #[test]
    fn test_v1_1_archive_requires_derived_state_section() {
        let mol = build_simple_methane();
        let mut data = encode_v1_0_sectioned_archive(&mol);
        let minor_offset = ARCHIVE_MAGIC.len() + 2;
        data[minor_offset..minor_offset + 2].copy_from_slice(&1u16.to_le_bytes());

        assert_eq!(
            decode_molecule_binary(&data).unwrap_err(),
            PickleError::MissingRequiredSection(SECTION_DERIVED_STATE)
        );
    }

    #[test]
    fn test_derived_state_rejects_valence_rows_for_a_different_graph() {
        let mol = build_simple_methane()
            .with_deserialized_derived_cache(BinaryDerivedState {
                valence: Some(ValenceAssignment {
                    explicit_valence: vec![4, 1, 1, 1, 1],
                    implicit_hydrogens: vec![0; 5],
                }),
                ..BinaryDerivedState::default()
            })
            .expect("build molecule with validated valence state");
        let mut derived_state =
            encode_derived_state(&fixture_input(&mol)).expect("encode derived state");
        assert_eq!(&derived_state[..3], &[0, 0, 1]);
        derived_state[3..7].copy_from_slice(&4_u32.to_le_bytes());

        let error = decode_derived_state(&derived_state, &mol).unwrap_err();
        assert!(
            matches!(&error, PickleError::InvalidArchive(message) if message.contains("explicit valence rows")),
            "unexpected malformed-cache error: {error:?}"
        );
    }

    #[test]
    fn test_sectioned_archive_rejects_unknown_required_section() {
        let mol = build_simple_methane();
        let molecule_state =
            encode_molecule_state(mol_to_legacy_binary(&fixture_input(&mol)).unwrap()).unwrap();
        let manifest = encode_manifest().unwrap();
        let mut data = Vec::new();
        data.extend_from_slice(ARCHIVE_MAGIC);
        write_u16_le(&mut data, ARCHIVE_MAJOR);
        write_u16_le(&mut data, ARCHIVE_MINOR);
        write_u16_le(&mut data, 3);
        write_archive_section(
            &mut data,
            SECTION_MANIFEST,
            MANIFEST_VERSION,
            0,
            SECTION_CODEC_POSTCARD,
            &manifest,
        )
        .unwrap();
        write_archive_section(
            &mut data,
            SECTION_MOLECULE_STATE,
            MOLECULE_STATE_VERSION,
            SECTION_FLAG_REQUIRED,
            SECTION_CODEC_POSTCARD,
            &molecule_state,
        )
        .unwrap();
        write_archive_section(
            &mut data,
            999,
            1,
            SECTION_FLAG_REQUIRED,
            SECTION_CODEC_RAW,
            b"future",
        )
        .unwrap();

        let err = decode_molecule_binary(&data).unwrap_err();
        assert!(matches!(err, PickleError::UnknownRequiredSection(999)));
    }

    #[test]
    fn test_sectioned_archive_rejects_trailing_bytes() {
        let mol = build_simple_methane();
        let mut data = fixture_encode(&mol).unwrap();
        data.push(0);
        let err = decode_molecule_binary(&data).unwrap_err();
        assert!(matches!(err, PickleError::InvalidArchive(_)));
    }

    #[test]
    fn test_methane_roundtrip() {
        let mol = build_simple_methane();
        assert_eq!(mol.num_atoms(), 5);
        assert_eq!(mol.num_bonds(), 4);

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "methane roundtrip failed");
    }

    #[test]
    fn test_methane_with_props_roundtrip() {
        let mut mol = build_simple_methane();
        mol = mol.with_name("methane_test");
        mol = mol.with_prop("key1", "value1");
        mol = mol.with_prop("key2", "value2");

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();

        assert_eq!(
            mol.properties().name().map(super::fixture_text),
            Some("methane_test")
        );
        assert_eq!(
            mol2.properties().name().map(super::fixture_text),
            Some("methane_test")
        );
        assert_eq!(mol2.prop("key1"), Some("value1"));
        assert_eq!(mol2.prop("key2"), Some("value2"));
        assert_record_equal(&mol, &mol2, "methane with properties roundtrip failed");
    }

    #[test]
    fn test_methane_with_2d_coords() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        for _ in 0..4 {
            builder.add_atom(AtomSpec::new(h));
        }
        for i in 0..4 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(i + 1),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        builder
            .set_2d_coordinates(vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [-0.5, 0.866],
                [-0.5, -0.866],
                [0.0, 1.0],
            ])
            .unwrap();
        let mol = builder.build().expect("build methane with coords");
        assert!(mol.coordinates_2d().is_some());
        assert_eq!(mol.coordinates_2d().unwrap().len(), 5);

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "methane with 2D coords roundtrip failed");
        assert!(mol2.coordinates_2d().is_some());
    }

    #[test]
    fn test_methane_with_3d_conformer() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        for _ in 0..4 {
            builder.add_atom(AtomSpec::new(h));
        }
        for i in 0..4 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(i + 1),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        builder
            .add_3d_conformer(vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [-0.5, 0.866, 0.0],
                [-0.5, -0.866, 0.0],
                [0.0, 1.0, 0.0],
            ])
            .unwrap();
        let mol = builder.build().expect("build methane with 3D");
        assert_eq!(mol.conformers_3d().len(), 1);

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "methane with 3D conformer roundtrip failed");
    }

    #[test]
    fn test_ethanol_roundtrip() {
        let mol = build_simple_ethanol();
        assert_eq!(mol.num_atoms(), 10);

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "ethanol roundtrip failed");
    }

    #[test]
    fn test_roundtrip_with_bond_props() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                    .with_prop("wiberg", "0.85")
                    .expect("valid key")
                    .with_stereo(BondStereo::Z)
                    .with_aromatic(false)
                    .with_conjugated(true),
            )
            .unwrap();
        let mol = builder.build().expect("build molecule");

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();

        let bond = &mol2.bonds()[0];
        assert_eq!(
            bond.prop("wiberg"),
            Some(&PropertyValue::String("0.85".into()))
        );
        assert_eq!(bond.stereo(), BondStereo::Z);
        assert!(!bond.is_aromatic());
        assert!(bond.is_conjugated());
        assert_record_equal(&mol, &mol2, "bond props roundtrip failed");
    }

    #[test]
    fn test_roundtrip_with_atom_props() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(
            AtomSpec::new(c)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_chiral_permutation(42)
                .with_isotope(13)
                .with_formal_charge(1)
                .with_radical_electrons(0)
                .with_hybridization(Hybridization::Sp3)
                .with_atom_map(5)
                .with_aromatic(false)
                .with_prop("test_key", "test_val")
                .expect("valid key"),
        );
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();
        let mol = builder.build().expect("build molecule");

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();

        let atom = &mol2.atoms()[0];
        assert_eq!(atom.atomic_number(), 6);
        assert_eq!(atom.chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(atom.chiral_permutation(), Some(42));
        assert_eq!(atom.isotope(), Some(13));
        assert_eq!(atom.formal_charge(), 1);
        assert_eq!(atom.radical_electrons(), 0);
        assert_eq!(atom.hybridization(), Hybridization::Sp3);
        assert_eq!(atom.atom_map(), Some(5));
        assert!(!atom.is_aromatic());
        assert_eq!(
            atom.prop("test_key"),
            Some(&PropertyValue::String("test_val".into()))
        );
        assert_record_equal(&mol, &mol2, "atom props roundtrip failed");
    }

    #[test]
    fn test_roundtrip_with_sdf_property_lists() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();

        let plist = SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "test_list",
            vec![Some("val1".into()), Some("val2".into())],
        );

        // Build MoleculeProperties separately since builder doesn't expose with_sdf_property_list
        let base = builder.build().expect("build molecule");

        // Construct the molecule with properties directly via internal API
        let mut mol_props = MoleculeProperties::default();
        mol_props = mol_props.with_sdf_property_list(plist);
        let topology = TopologyBlock {
            atoms: base.atoms().to_vec(),
            bonds: base.bonds().to_vec(),
            adjacency: cosmolkit_model::AdjacencyList::from_topology(
                base.num_atoms(),
                base.bonds(),
            ),
            substance_groups: vec![],
            stereo_groups: vec![],
        };
        let coord_block = CoordinateBlock {
            conformers_2d: vec![],
            conformers_3d: vec![],
            source_coordinate_dim: None,
            source_conformer_order: None,
        };
        let mol = BinaryRecord::from_blocks(topology, coord_block, mol_props)
            .expect("build molecule with property lists");
        assert_eq!(mol.properties().sdf_property_lists().len(), 1);

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "SDF property list roundtrip failed");
    }

    #[test]
    fn test_roundtrip_with_stereo_groups() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();

        let sg = StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(0)], vec![])
            .expect("valid distinct stereo members")
            .with_id(1);
        builder.add_stereo_group(sg).unwrap();

        let mol = builder.build().expect("build molecule");
        assert_eq!(mol.stereo_groups().len(), 1);

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "stereo group roundtrip failed");
    }

    #[test]
    fn test_invalid_version() {
        let data = vec![0xFF, 0x00, 0x00, 0x00, 0x00];
        let result = decode_molecule_binary(&data);
        assert!(result.is_err(), "expected error for unsupported version");
        match result {
            Err(PickleError::UnsupportedVersion(v)) => assert_eq!(v, 0xFF),
            _ => panic!("expected UnsupportedVersion error"),
        }
    }

    #[test]
    fn test_truncated_data() {
        let data = vec![0x01];
        let result = decode_molecule_binary(&data);
        assert!(result.is_err(), "expected error for truncated data");
        match result {
            Err(PickleError::UnexpectedEof) => {}
            _ => panic!("expected UnexpectedEof error"),
        }
    }

    #[test]
    fn test_methane_with_sdf_data_fields() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        for _ in 0..4 {
            builder.add_atom(AtomSpec::new(h));
        }
        for i in 0..4 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(i + 1),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        builder = builder.with_sdf_data_field("PUBCHEM_IUPAC_NAME", "methane");
        builder = builder.with_sdf_data_field("PUBCHEM_MOLECULAR_FORMULA", "CH4");
        let mol = builder.build().expect("build methane");

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "methane with SDF data fields roundtrip failed");
    }

    #[test]
    fn test_ethanol_with_stereo_atoms() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();

        // Build ethylene-like with stereo atoms on double bond
        // C=C with stereo
        builder.add_atom(AtomSpec::new(c)); // 0
        builder.add_atom(AtomSpec::new(c)); // 1
        builder.add_atom(AtomSpec::new(h)); // 2
        builder.add_atom(AtomSpec::new(h)); // 3
        builder.add_atom(AtomSpec::new(h)); // 4
        builder.add_atom(AtomSpec::new(h)); // 5

        // Bonds for substituents
        builder
            .add_bond(
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double)
                    .with_stereo(BondStereo::Z)
                    .with_stereo_atoms(AtomId::new(2), AtomId::new(4)),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(2),
                BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(3),
                BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                AtomId::new(1),
                AtomId::new(4),
                BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                AtomId::new(1),
                AtomId::new(5),
                BondOrder::Single,
            ))
            .unwrap();

        let mol = builder.build().expect("build ethylene");
        assert_eq!(
            mol.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(4)])
        );

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "ethylene roundtrip failed");
        assert_eq!(
            mol2.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(4)])
        );
    }

    #[test]
    fn test_multiple_conformers() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = DetachedFixtureBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();

        builder
            .add_3d_conformer(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
            .unwrap();
        builder
            .add_3d_conformer(vec![[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]])
            .unwrap();

        let mol = builder.build().expect("build with conformers");
        assert_eq!(mol.conformers_3d().len(), 2);

        let data = fixture_encode(&mol).unwrap();
        let mol2 = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &mol2, "multiple conformers roundtrip failed");
    }

    #[test]
    fn canonical_native_all_property_kinds_order_flags_pdb_and_coordinates() {
        let mut mol = build_simple_methane();
        let values = [
            PropertyValue::String("001 true".into()),
            PropertyValue::Int(i32::MIN),
            PropertyValue::UInt(u32::MAX),
            PropertyValue::IntVector(vec![i32::MIN, 0, i32::MAX]),
            PropertyValue::Double(f64::from_bits(0x7ff8_1234_5678_9abc)),
            PropertyValue::Bool(true),
        ];
        for (i, value) in values.into_iter().enumerate() {
            let key = format!("z{}", 5 - i);
            if i % 2 == 0 {
                mol.topology.atoms[0]
                    .set_computed_prop(&key, value.clone())
                    .unwrap();
                mol.topology.bonds[0]
                    .set_computed_prop(&key, value)
                    .unwrap();
            } else {
                mol.topology.atoms[0].set_prop(&key, value.clone()).unwrap();
                mol.topology.bonds[0].set_prop(&key, value).unwrap();
            }
        }
        mol.topology.atoms[0]
            .set_prop("negative_zero", PropertyValue::Double(-0.0))
            .unwrap();
        mol.topology.atoms[0].set_temporary_flags(u64::MAX);
        mol.topology.bonds[0].set_temporary_flags(0x8123_4567_89ab_cdef);
        let info = AtomPdbResidueInfo::new(" CA ", i32::MIN, "GLY", i32::MAX, "long-chain", true)
            .with_alt_loc("B")
            .with_insertion_code("A")
            .with_occupancy(-0.0)
            .with_temp_factor(f64::from_bits(0x7ff8_1111_2222_3333))
            .with_secondary_structure(u32::MAX)
            .with_segment_number(42)
            .with_monomer_class("amino");
        mol.topology.atoms[0].set_pdb_residue_info(Some(info.clone()));
        replace_atom_template_attachment_order(
            &mut mol.topology.atoms[0],
            Some(
                TemplateAttachmentOrder::new(vec![
                    TemplateAttachment::new(AtomId::new(2), "tail"),
                    TemplateAttachment::new(AtomId::new(1), "head"),
                ])
                .unwrap(),
            ),
        );
        for id in [7, usize::MAX] {
            mol.coordinates
                .conformers_2d
                .push(Conformer2D::new(id, vec![[-0.0, 1.25]; 5]).with_prop("conformer", "two"));
            mol.coordinates.conformers_3d.push(
                Conformer3D::new(id, vec![[1.0, -0.0, 2.5]; 5], id == 7)
                    .with_prop("conformer", "three"),
            );
        }
        mol.coordinates.source_coordinate_dim = Some(CoordinateDimension::ThreeD);
        mol.properties = mol.properties.with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "typed",
            vec![
                Some(PropertyValue::UInt(u32::MAX)),
                None,
                Some(PropertyValue::Double(-0.0)),
                Some(PropertyValue::Bool(false)),
                Some(PropertyValue::IntVector(vec![])),
            ],
        ));
        mol.derived.valid_bits = Some(0xe0);
        let data = fixture_encode(&mol).unwrap();
        let restored = decode_molecule_binary(&data).unwrap();
        assert_record_equal(&mol, &restored, "complete canonical native state");
        assert_eq!(restored.derived.valid_bits, Some(0xe0));
        assert_eq!(
            ordered_atom_properties(&restored.topology.atoms[0])
                .map(|(k, _)| super::fixture_text(k))
                .collect::<Vec<_>>(),
            [
                "__computedProps",
                "z5",
                "z4",
                "z3",
                "z2",
                "z1",
                "z0",
                "negative_zero"
            ]
        );
        assert_eq!(
            ordered_bond_properties(&restored.topology.bonds[0])
                .map(|(k, _)| super::fixture_text(k))
                .collect::<Vec<_>>(),
            ["__computedProps", "z5", "z4", "z3", "z2", "z1", "z0"]
        );
        // RDProps::setProp inserts the reserved vector before the first
        // computed ordinary property; subsequent membership keeps its order.
        let computed = PropertyValue::StringVector(vec!["z5".into(), "z3".into(), "z1".into()]);
        assert_eq!(
            restored.topology.atoms[0].prop("__computedProps"),
            Some(&computed)
        );
        assert_eq!(
            restored.topology.bonds[0].prop("__computedProps"),
            Some(&computed)
        );
        assert_eq!(restored.topology.atoms[0].pdb_residue_info(), Some(&info));
        assert_eq!(
            restored.coordinates.conformers_2d[0].coordinates()[0][0].to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(restored.coordinates.conformers_3d[1].id(), usize::MAX);
        assert_eq!(data, fixture_encode(&restored).unwrap());
    }

    #[test]
    fn canonical_native_sgroup_full_three_points_vectors_crossing_and_stereo_write_id() {
        let mut mol = build_simple_methane();
        let display = SGroupDisplay {
            brackets: vec![SGroupBracket::new([
                [1., 2., 3.],
                [4., 5., 6.],
                [7., 8., 9.],
            ])],
            field_position: Some([10., 11.]),
            display_tag: Some("DISPLAY".into()),
        };
        mol.topology.substance_groups.push(
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                .with_atoms(vec![AtomId::new(0), AtomId::new(1)])
                .with_bonds(vec![BondId::new(0)])
                .with_bond_role(BondId::new(0), SGroupBondRole::Contained)
                .with_display(display)
                .with_cstates(vec![SGroupCState::new(BondId::new(0), [12., 13., 14.])])
                .with_head_crossing_bonds(vec![BondId::new(0)])
                .with_crossing_bond_correspondence(vec![BondId::new(0)]),
        );
        mol.topology.stereo_groups.push(
            StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(0)], vec![])
                .expect("valid distinct stereo members")
                .with_id(4)
                .with_write_id(17),
        );
        let data = fixture_encode(&mol).unwrap();
        let restored = decode_molecule_binary(&data).unwrap();
        assert_record_equal(
            &mol,
            &restored,
            "full current SGroup and stereo group values",
        );
        assert_eq!(data, fixture_encode(&restored).unwrap());
    }

    #[test]
    fn historical_raw_versions_and_archive11_preserve_strings_without_type_inference() {
        let mut mol = build_simple_methane();
        mol.topology.atoms[0]
            .set_prop("numeric-looking", "4294967295")
            .unwrap();
        mol.topology.bonds[0]
            .set_prop("bool-looking", "true")
            .unwrap();
        for version in 1..=3 {
            let data = mol_to_legacy_binary_version(&fixture_input(&mol), version).unwrap();
            let restored = decode_molecule_binary(&data).unwrap();
            assert_record_equal(&mol, &restored, "legacy raw complete state");
            assert_eq!(
                restored.topology.atoms[0].prop("numeric-looking"),
                Some(&PropertyValue::String("4294967295".into()))
            );
            assert_eq!(restored.derived.valid_bits, None);
        }
        let bytes = fixture_encode_legacy12(&mol).unwrap();
        let (_, sections) = read_archive_sections(&bytes).unwrap();
        let mut historical = ARCHIVE_MAGIC.to_vec();
        write_u16_le(&mut historical, 1);
        write_u16_le(&mut historical, 1);
        write_u16_le(&mut historical, 3);
        for section in sections
            .into_iter()
            .filter(|s| s.id != SECTION_CANONICAL_STATE)
        {
            write_archive_section(
                &mut historical,
                section.id,
                section.version,
                section.flags,
                section.codec,
                section.payload,
            )
            .unwrap();
        }
        // Build the raw4 incompatibility explicitly rather than assuming the
        // legacy1.2 fixture writer emits the newest raw version.
        let (_, sections) = read_archive_sections(&bytes).unwrap();
        historical.truncate(14);
        // Only the raw version tag is needed: archive 1.1 rejects raw4 before
        // decoding its body. Do not use the evolving current raw writer here.
        let raw4 = encode_molecule_state(vec![4]).unwrap();
        for section in sections
            .into_iter()
            .filter(|s| s.id != SECTION_CANONICAL_STATE)
        {
            write_archive_section(
                &mut historical,
                section.id,
                section.version,
                section.flags,
                section.codec,
                if section.id == SECTION_MOLECULE_STATE {
                    &raw4
                } else {
                    section.payload
                },
            )
            .unwrap();
        }
        // Preserve the original invalid 1.1/raw4 fixture as a rejection.
        assert_eq!(
            decode_molecule_binary(&historical).unwrap_err(),
            PickleError::InvalidArchive("old archive cannot contain raw4".into())
        );
        // Historical archive 1.1 carries a legacy raw1..3 molecule section.
        // Retain its manifest, derived companion and every string value.
        let legacy_state =
            encode_molecule_state(mol_to_legacy_binary_version(&fixture_input(&mol), 3).unwrap())
                .unwrap();
        let (_, sections) = read_archive_sections(&bytes).unwrap();
        let mut historical = ARCHIVE_MAGIC.to_vec();
        write_u16_le(&mut historical, 1);
        write_u16_le(&mut historical, 1);
        write_u16_le(&mut historical, 3);
        for section in sections
            .into_iter()
            .filter(|section| section.id != SECTION_CANONICAL_STATE)
        {
            write_archive_section(
                &mut historical,
                section.id,
                section.version,
                section.flags,
                section.codec,
                if section.id == SECTION_MOLECULE_STATE {
                    &legacy_state
                } else {
                    section.payload
                },
            )
            .unwrap();
        }
        let restored = decode_molecule_binary(&historical).unwrap();
        assert_record_equal(&mol, &restored, "historical1.1 all retained fields");
        assert_eq!(
            restored.topology.atoms[0].prop("numeric-looking"),
            Some(&PropertyValue::String("4294967295".into()))
        );
        assert_eq!(
            restored.topology.bonds[0].prop("bool-looking"),
            Some(&PropertyValue::String("true".into()))
        );
        assert_eq!(restored.derived.valid_bits, None);
    }

    #[test]
    fn canonical_native_required_section_missing_duplicate_invalid_flags_and_versions() {
        let data = fixture_encode_legacy12(&build_simple_methane()).unwrap();
        let (_, sections) = read_archive_sections(&data).unwrap();
        for (id, version, flags, codec, expected) in [
            (
                SECTION_CANONICAL_STATE,
                CANONICAL_STATE_VERSION + 1,
                1,
                0,
                "version",
            ),
            (SECTION_CANONICAL_STATE, 1, 0, 0, "required"),
            (SECTION_CANONICAL_STATE, 1, 2, 0, "flags"),
            (SECTION_CANONICAL_STATE, 1, 1, 1, "codec"),
        ] {
            let mut mutated = ARCHIVE_MAGIC.to_vec();
            write_u16_le(&mut mutated, 1);
            write_u16_le(&mut mutated, 2);
            write_u16_le(&mut mutated, 4);
            for section in &sections {
                if section.id == id {
                    write_archive_section(&mut mutated, id, version, flags, codec, section.payload)
                        .unwrap();
                } else {
                    write_archive_section(
                        &mut mutated,
                        section.id,
                        section.version,
                        section.flags,
                        section.codec,
                        section.payload,
                    )
                    .unwrap();
                }
            }
            assert!(decode_molecule_binary(&mutated).is_err(), "{expected}");
        }
        let mut missing = data.clone();
        missing[12..14].copy_from_slice(&3u16.to_le_bytes());
        missing.truncate(data.len() - 10 - sections.last().unwrap().payload.len());
        assert_eq!(
            decode_molecule_binary(&missing).unwrap_err(),
            PickleError::MissingRequiredSection(SECTION_CANONICAL_STATE)
        );
        let mut duplicate = data.clone();
        duplicate[12..14].copy_from_slice(&5u16.to_le_bytes());
        let section = &sections[0];
        write_archive_section(
            &mut duplicate,
            section.id,
            section.version,
            section.flags,
            section.codec,
            section.payload,
        )
        .unwrap();
        assert_eq!(
            decode_molecule_binary(&duplicate).unwrap_err(),
            PickleError::DuplicateSection(SECTION_MANIFEST)
        );
    }

    #[test]
    fn canonical_native_truncation_and_strict_boolean_and_value_tags_do_not_panic() {
        let data = fixture_encode(&build_simple_methane()).unwrap();
        for len in 0..data.len() {
            assert!(
                decode_molecule_binary(&data[..len]).is_err(),
                "truncation {len}"
            );
        }
        let mut raw = mol_to_legacy_binary(&fixture_input(&build_simple_methane())).unwrap();
        raw[7] = 2; // isotope presence follows version/count/element/charge.
        assert!(matches!(
            decode_molecule_binary(&raw),
            Err(PickleError::InvalidEnumValue {
                type_name: "bool",
                value: 2
            })
        ));
        assert!(matches!(
            read_value(&mut PickleReader::new(&[6])),
            Err(PickleError::InvalidEnumValue {
                type_name: "PropertyValue",
                value: 6
            })
        ));
        assert!(read_value(&mut PickleReader::new(&[3, 255, 255, 255, 255])).is_err());
        // Raw4 is the supported NativeStateV2 tag; its missing atom count
        // is truncation, as for the supported historical raw1..3 tags.
        for value in 1..=4 {
            assert_eq!(
                decode_molecule_binary(&[value]).unwrap_err(),
                PickleError::UnexpectedEof
            );
        }
        for value in [0, 255] {
            assert!(
                matches!(decode_molecule_binary(&[value]),Err(PickleError::UnsupportedVersion(v)) if v==value)
            );
        }
    }
}

#[cfg(test)]
fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes())
        .expect("original text fixture must retain exact UTF-8 bytes")
}
#[cfg(test)]
fn fixture_writer_text(value: cosmolkit_model::PropertyText) -> String {
    String::from_utf8(value.into_bytes())
        .expect("original writer fixture must retain exact UTF-8 bytes")
}
