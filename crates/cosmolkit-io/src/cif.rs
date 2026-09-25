//! Detached CIF lexical, document, and scalar-formatting primitives.
//!
//! This module is the single CIF representation used by the Gemmi-primary
//! structural readers. It intentionally does not contain BioStructure
//! conversion or a serializer policy.

use std::collections::HashSet;
use std::fmt;

/// Validation performed after the CIF grammar has been parsed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CifCheckLevel {
    /// Grammar and document-building actions only (`check_level == 0`).
    Syntax,
    /// Missing values and duplicate names (`check_level == 1`).
    Default,
    /// Default checks plus bare block names and empty loops (`check_level > 1`).
    Strict,
}

/// Stable class for a CIF read failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CifReadErrorKind {
    Syntax,
    MissingValue,
    DuplicateName,
    InvalidLoop,
    InvalidValue,
    OutOfRange,
}

/// Structured failure from the detached CIF layer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CifReadError {
    kind: CifReadErrorKind,
    source: String,
    line: usize,
    column: usize,
    message: String,
}

impl CifReadError {
    fn new(
        kind: CifReadErrorKind,
        source: &str,
        line: usize,
        column: usize,
        message: impl Into<String>,
    ) -> Self {
        Self {
            kind,
            source: source.to_owned(),
            line,
            column,
            message: message.into(),
        }
    }

    pub fn kind(&self) -> CifReadErrorKind {
        self.kind
    }

    pub fn source(&self) -> &str {
        &self.source
    }

    pub fn line(&self) -> usize {
        self.line
    }

    pub fn column(&self) -> usize {
        self.column
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for CifReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}:{}: {}",
            self.source, self.line, self.column, self.message
        )
    }
}

impl std::error::Error for CifReadError {}

/// One raw CIF value and its source start position.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CifValue {
    raw: String,
    line: usize,
    column: usize,
}

impl CifValue {
    fn new(raw: String, line: usize, column: usize) -> Self {
        Self { raw, line, column }
    }

    pub fn raw(&self) -> &str {
        &self.raw
    }

    pub fn line(&self) -> usize {
        self.line
    }

    pub fn column(&self) -> usize {
        self.column
    }

    pub fn is_null(&self) -> bool {
        cif_is_null(&self.raw)
    }

    pub fn decoded(&self) -> String {
        cif_as_string(&self.raw)
    }
}

/// A CIF tag/value pair. `value == None` is retained at syntax-only level.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CifPair {
    tag: String,
    value: Option<CifValue>,
    line: usize,
}

impl CifPair {
    pub fn tag(&self) -> &str {
        &self.tag
    }

    pub fn value(&self) -> Option<&CifValue> {
        self.value.as_ref()
    }

    pub fn line(&self) -> usize {
        self.line
    }
}

/// A row-major CIF loop.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CifLoop {
    tags: Vec<String>,
    values: Vec<CifValue>,
    line: usize,
}

impl CifLoop {
    pub fn tags(&self) -> &[String] {
        &self.tags
    }

    pub fn values(&self) -> &[CifValue] {
        &self.values
    }

    pub fn line(&self) -> usize {
        self.line
    }

    pub fn width(&self) -> usize {
        self.tags.len()
    }

    pub fn len(&self) -> usize {
        if self.tags.is_empty() {
            0
        } else {
            self.values.len() / self.tags.len()
        }
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    pub fn find_tag(&self, tag: &str) -> Option<usize> {
        self.tags
            .iter()
            .position(|candidate| candidate.eq_ignore_ascii_case(tag))
    }

    pub fn value(&self, row: usize, column: usize) -> Option<&CifValue> {
        if column >= self.width() || row >= self.len() {
            return None;
        }
        self.values.get(row * self.width() + column)
    }

    fn common_prefix(&self) -> String {
        let Some(first) = self.tags.first() else {
            return String::new();
        };
        let mut len = first.len();
        for tag in self.tags.iter().skip(1) {
            len = first
                .as_bytes()
                .iter()
                .zip(tag.as_bytes())
                .take(len)
                .position(|(a, b)| !a.eq_ignore_ascii_case(b))
                .unwrap_or(len.min(tag.len()));
        }
        first[..len].to_owned()
    }
}

/// A stored CIF item. Comments are intentionally not stored by pinned Gemmi.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CifItem {
    Pair(CifPair),
    Loop(CifLoop),
    Frame(CifBlock),
}

impl CifItem {
    pub fn line(&self) -> usize {
        match self {
            Self::Pair(pair) => pair.line,
            Self::Loop(loop_) => loop_.line,
            Self::Frame(frame) => frame.line,
        }
    }
}

/// One CIF data/global/save block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CifBlock {
    name: String,
    items: Vec<CifItem>,
    line: usize,
}

impl CifBlock {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn items(&self) -> &[CifItem] {
        &self.items
    }

    pub fn line(&self) -> usize {
        self.line
    }

    pub fn find_pair(&self, tag: &str) -> Option<&CifPair> {
        self.items.iter().find_map(|item| match item {
            CifItem::Pair(pair) if pair.tag.eq_ignore_ascii_case(tag) => Some(pair),
            _ => None,
        })
    }

    pub fn find_value(&self, tag: &str) -> Option<&CifValue> {
        if let Some(pair) = self.find_pair(tag) {
            return pair.value.as_ref();
        }
        for item in &self.items {
            if let CifItem::Loop(loop_) = item
                && let Some(column) = loop_.find_tag(tag)
                && loop_.len() == 1
            {
                return loop_.value(0, column);
            }
        }
        None
    }

    pub fn find_loop(&self, tag: &str) -> Option<&CifLoop> {
        self.items.iter().find_map(|item| match item {
            CifItem::Loop(loop_) if loop_.find_tag(tag).is_some() => Some(loop_),
            _ => None,
        })
    }

    pub fn find_values(&self, tag: &str) -> Option<CifColumn<'_>> {
        for (item_index, item) in self.items.iter().enumerate() {
            match item {
                CifItem::Loop(loop_) => {
                    if let Some(column) = loop_.find_tag(tag) {
                        return Some(CifColumn {
                            block: self,
                            item_index,
                            column,
                        });
                    }
                }
                CifItem::Pair(pair) if pair.tag.eq_ignore_ascii_case(tag) => {
                    return Some(CifColumn {
                        block: self,
                        item_index,
                        column: 0,
                    });
                }
                _ => {}
            }
        }
        None
    }

    pub fn has_tag(&self, tag: &str) -> bool {
        self.find_values(tag).is_some()
    }

    pub fn has_any_value(&self, tag: &str) -> bool {
        self.find_values(tag)
            .is_some_and(|column| column.iter().any(|value| !value.is_null()))
    }

    pub fn find_frame(&self, name: &str) -> Option<&CifBlock> {
        self.items.iter().find_map(|item| match item {
            CifItem::Frame(frame) if frame.name.eq_ignore_ascii_case(name) => Some(frame),
            _ => None,
        })
    }

    pub fn find<'a>(&'a self, prefix: &str, tags: &[&str]) -> Result<CifTable<'a>, CifReadError> {
        // Gemmi✔️✔️: inline Table Block::find(const std::string& prefix,
        // Gemmi✔️✔️:                          const std::vector<std::string>& tags) {
        // Gemmi✔️✔️:   Item* loop_item = nullptr;
        // Gemmi✔️✔️:   if (!tags.empty()) {
        // Gemmi✔️✔️:     if (tags[0][0] == '?')
        // Gemmi✔️✔️:       fail("The first tag in find() cannot be ?optional.");
        // Gemmi✔️✔️:     loop_item = find_loop(prefix + tags[0]).item();
        // Gemmi✔️✔️:   }
        // Gemmi✔️✔️:   std::vector<int> indices;
        // Gemmi✔️✔️:   indices.reserve(tags.size());
        // Gemmi✔️✔️:   if (loop_item) {
        // Gemmi✔️✔️:     for (const std::string& tag : tags) {
        // Gemmi✔️✔️:       std::string full_tag = prefix + (tag[0] != '?' ? tag : tag.substr(1));
        // Gemmi✔️✔️:       int idx = loop_item->loop.find_tag(full_tag);
        // Gemmi✔️✔️:       if (idx == -1 && tag[0] != '?') {
        // Gemmi✔️✔️:         loop_item = nullptr;
        // Gemmi✔️✔️:         indices.clear();
        // Gemmi✔️✔️:         break;
        // Gemmi✔️✔️:       }
        // Gemmi✔️✔️:       indices.push_back(idx);
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:   } else {
        // Gemmi✔️✔️:     for (const std::string& tag : tags) {
        // Gemmi✔️✔️:       std::string full_tag = prefix + (tag[0] != '?' ? tag : tag.substr(1));
        // Gemmi✔️✔️:       if (const Item* p = find_pair_item(full_tag)) {
        // Gemmi✔️✔️:         indices.push_back(p - items.data());
        // Gemmi✔️✔️:       } else if (tag[0] == '?') {
        // Gemmi✔️✔️:         indices.push_back(-1);
        // Gemmi✔️✔️:       } else {
        // Gemmi✔️✔️:         indices.clear();
        // Gemmi✔️✔️:         break;
        // Gemmi✔️✔️:       }
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:   }
        // Gemmi✔️✔️:   return Table{loop_item, *this, indices, prefix.length()};
        // Gemmi✔️✔️: }
        // Behavior review: optional columns, one-loop identity and pair-table
        // fallback match the source. Complexity review: the same linear item/tag
        // scans and O(number-of-requested-columns) position vector are used.
        if tags.first().is_some_and(|tag| tag.starts_with('?')) {
            return Err(self.error(
                CifReadErrorKind::InvalidValue,
                self.line,
                1,
                "The first tag in find() cannot be ?optional.",
            ));
        }
        let mut positions = Vec::with_capacity(tags.len());
        let loop_index = tags.first().and_then(|tag| {
            let full = format!("{prefix}{tag}");
            self.items.iter().position(
                |item| matches!(item, CifItem::Loop(loop_) if loop_.find_tag(&full).is_some()),
            )
        });
        if let Some(item_index) = loop_index {
            let CifItem::Loop(loop_) = &self.items[item_index] else {
                unreachable!();
            };
            for tag in tags {
                let optional = tag.starts_with('?');
                let suffix = tag.strip_prefix('?').unwrap_or(tag);
                let full = format!("{prefix}{suffix}");
                let position = loop_.find_tag(&full);
                if position.is_none() && !optional {
                    return Ok(CifTable::empty(self, prefix.len()));
                }
                positions.push(position);
            }
            return Ok(CifTable {
                block: self,
                loop_index: Some(item_index),
                positions,
                prefix_len: prefix.len(),
            });
        }
        for tag in tags {
            let optional = tag.starts_with('?');
            let suffix = tag.strip_prefix('?').unwrap_or(tag);
            let full = format!("{prefix}{suffix}");
            let position = self.items.iter().position(
                |item| matches!(item, CifItem::Pair(pair) if pair.tag.eq_ignore_ascii_case(&full)),
            );
            if position.is_none() && !optional {
                return Ok(CifTable::empty(self, prefix.len()));
            }
            positions.push(position);
        }
        Ok(CifTable {
            block: self,
            loop_index: None,
            positions,
            prefix_len: prefix.len(),
        })
    }

    pub fn find_mmcif_category(&self, category: &str) -> Result<CifTable<'_>, CifReadError> {
        let mut prefix = category.to_owned();
        if !prefix.starts_with('_') {
            return Err(self.error(
                CifReadErrorKind::InvalidValue,
                self.line,
                1,
                format!("Category should start with '_', got: {category}"),
            ));
        }
        if !prefix.ends_with('.') {
            prefix.push('.');
        }
        let prefix_lc = prefix.to_ascii_lowercase();
        let mut pair_positions = Vec::new();
        for (index, item) in self.items.iter().enumerate() {
            match item {
                CifItem::Loop(loop_)
                    if loop_
                        .tags
                        .first()
                        .is_some_and(|tag| tag.to_ascii_lowercase().starts_with(&prefix_lc)) =>
                {
                    let mut positions = Vec::with_capacity(loop_.tags.len());
                    for (position, tag) in loop_.tags.iter().enumerate() {
                        if !tag.to_ascii_lowercase().starts_with(&prefix_lc) {
                            return Err(self.error(
                                CifReadErrorKind::InvalidLoop,
                                loop_.line,
                                1,
                                format!("Tag {tag} in loop with {prefix_lc}"),
                            ));
                        }
                        positions.push(Some(position));
                    }
                    return Ok(CifTable {
                        block: self,
                        loop_index: Some(index),
                        positions,
                        prefix_len: prefix_lc.len(),
                    });
                }
                CifItem::Pair(pair) if pair.tag.to_ascii_lowercase().starts_with(&prefix_lc) => {
                    pair_positions.push(Some(index));
                }
                _ => {}
            }
        }
        Ok(CifTable {
            block: self,
            loop_index: None,
            positions: pair_positions,
            prefix_len: prefix_lc.len(),
        })
    }

    pub fn mmcif_category_names(&self) -> Vec<String> {
        let mut categories: Vec<String> = Vec::new();
        for item in &self.items {
            let tag = match item {
                CifItem::Pair(pair) => Some(pair.tag.as_str()),
                CifItem::Loop(loop_) => loop_.tags.first().map(String::as_str),
                CifItem::Frame(_) => None,
            };
            if let Some(tag) = tag
                && let Some(dot) = tag.find('.')
            {
                let category = &tag[..=dot];
                if !categories.iter().any(|known| tag.starts_with(known)) {
                    categories.push(category.to_owned());
                }
            }
        }
        categories
    }

    fn error(
        &self,
        kind: CifReadErrorKind,
        line: usize,
        column: usize,
        message: impl Into<String>,
    ) -> CifReadError {
        CifReadError::new(kind, "cif", line, column, message)
    }
}

/// Parsed CIF document.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CifDocument {
    source: String,
    blocks: Vec<CifBlock>,
}

impl CifDocument {
    pub fn source(&self) -> &str {
        &self.source
    }

    pub fn blocks(&self) -> &[CifBlock] {
        &self.blocks
    }

    pub fn sole_block(&self) -> Result<&CifBlock, CifReadError> {
        if self.blocks.len() > 1 {
            return Err(CifReadError::new(
                CifReadErrorKind::InvalidValue,
                &self.source,
                1,
                1,
                format!("single data block expected, got {}", self.blocks.len()),
            ));
        }
        self.blocks.first().ok_or_else(|| {
            CifReadError::new(
                CifReadErrorKind::InvalidValue,
                &self.source,
                1,
                1,
                "single data block expected, got 0",
            )
        })
    }

    pub fn find_block(&self, name: &str) -> Option<&CifBlock> {
        self.blocks.iter().find(|block| block.name == name)
    }
}

/// A pair or loop column.
#[derive(Debug, Clone, Copy)]
pub struct CifColumn<'a> {
    block: &'a CifBlock,
    item_index: usize,
    column: usize,
}

impl<'a> CifColumn<'a> {
    pub fn len(&self) -> usize {
        match &self.block.items[self.item_index] {
            CifItem::Pair(pair) => usize::from(pair.value.is_some()),
            CifItem::Loop(loop_) => loop_.len(),
            CifItem::Frame(_) => 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn tag(&self) -> &str {
        match &self.block.items[self.item_index] {
            CifItem::Pair(pair) => &pair.tag,
            CifItem::Loop(loop_) => &loop_.tags[self.column],
            CifItem::Frame(_) => unreachable!(),
        }
    }

    pub fn get(&self, index: isize) -> Option<&'a CifValue> {
        let len = self.len() as isize;
        let normalized = if index < 0 { len + index } else { index };
        if !(0..len).contains(&normalized) {
            return None;
        }
        match &self.block.items[self.item_index] {
            CifItem::Pair(pair) => pair.value.as_ref(),
            CifItem::Loop(loop_) => loop_.value(normalized as usize, self.column),
            CifItem::Frame(_) => None,
        }
    }

    pub fn iter(&self) -> CifColumnIter<'a> {
        CifColumnIter {
            column: *self,
            next: 0,
        }
    }
}

pub struct CifColumnIter<'a> {
    column: CifColumn<'a>,
    next: usize,
}

impl<'a> Iterator for CifColumnIter<'a> {
    type Item = &'a CifValue;

    fn next(&mut self) -> Option<Self::Item> {
        let value = self.column.get(self.next as isize);
        self.next += usize::from(value.is_some());
        value
    }
}

/// Read-only selection of columns from one loop or from tag/value pairs.
#[derive(Debug, Clone)]
pub struct CifTable<'a> {
    block: &'a CifBlock,
    loop_index: Option<usize>,
    positions: Vec<Option<usize>>,
    prefix_len: usize,
}

impl<'a> CifTable<'a> {
    fn empty(block: &'a CifBlock, prefix_len: usize) -> Self {
        Self {
            block,
            loop_index: None,
            positions: Vec::new(),
            prefix_len,
        }
    }

    pub fn is_present(&self) -> bool {
        !self.positions.is_empty()
    }

    pub fn width(&self) -> usize {
        self.positions.len()
    }

    pub fn len(&self) -> usize {
        if let Some(index) = self.loop_index {
            match &self.block.items[index] {
                CifItem::Loop(loop_) => loop_.len(),
                _ => 0,
            }
        } else {
            usize::from(!self.positions.is_empty())
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn prefix(&self) -> Result<&str, CifReadError> {
        for (column, position) in self.positions.iter().enumerate() {
            if position.is_some() {
                let tag = self.tag(column).expect("present position has tag");
                return Ok(&tag[..self.prefix_len.min(tag.len())]);
            }
        }
        Err(CifReadError::new(
            CifReadErrorKind::InvalidValue,
            "cif",
            self.block.line,
            1,
            "The table has no columns.",
        ))
    }

    pub fn tag(&self, column: usize) -> Option<&'a str> {
        let position = *self.positions.get(column)?;
        let position = position?;
        if let Some(item_index) = self.loop_index {
            match &self.block.items[item_index] {
                CifItem::Loop(loop_) => loop_.tags.get(position).map(String::as_str),
                _ => None,
            }
        } else {
            match self.block.items.get(position) {
                Some(CifItem::Pair(pair)) => Some(&pair.tag),
                _ => None,
            }
        }
    }

    pub fn has_column(&self, column: usize) -> bool {
        self.positions.get(column).is_some_and(Option::is_some)
    }

    pub fn row(&'a self, index: isize) -> Option<CifRow<'a>> {
        let len = self.len() as isize;
        let normalized = if index < 0 { len + index } else { index };
        if !(0..len).contains(&normalized) {
            return None;
        }
        Some(CifRow {
            table: self,
            row: normalized as usize,
        })
    }

    pub fn one(&'a self) -> Result<CifRow<'a>, CifReadError> {
        if self.len() != 1 {
            return Err(CifReadError::new(
                CifReadErrorKind::InvalidValue,
                "cif",
                self.block.line,
                1,
                format!("Expected one value, found {}", self.len()),
            ));
        }
        Ok(self.row(0).expect("one-row table"))
    }

    pub fn iter(&'a self) -> CifTableIter<'a> {
        CifTableIter {
            table: self,
            next: 0,
        }
    }

    pub fn find_row(&'a self, decoded: &str) -> Result<CifRow<'a>, CifReadError> {
        for row in self.iter() {
            if row.get(0).is_some_and(|value| value.decoded() == decoded) {
                return Ok(row);
            }
        }
        let tag = self.tag(0).unwrap_or("<missing>");
        Err(CifReadError::new(
            CifReadErrorKind::InvalidValue,
            "cif",
            self.block.line,
            1,
            format!("Not found in {tag}: {decoded}"),
        ))
    }
}

pub struct CifTableIter<'a> {
    table: &'a CifTable<'a>,
    next: usize,
}

impl<'a> Iterator for CifTableIter<'a> {
    type Item = CifRow<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let row = self.table.row(self.next as isize);
        self.next += usize::from(row.is_some());
        row
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CifRow<'a> {
    table: &'a CifTable<'a>,
    row: usize,
}

impl<'a> CifRow<'a> {
    pub fn len(&self) -> usize {
        self.table.width()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn has(&self, column: usize) -> bool {
        self.table.has_column(column)
    }

    pub fn has_value(&self, column: usize) -> bool {
        self.get(column).is_some_and(|value| !value.is_null())
    }

    pub fn get(&self, column: usize) -> Option<&'a CifValue> {
        let position = *self.table.positions.get(column)?;
        let position = position?;
        if let Some(item_index) = self.table.loop_index {
            match &self.table.block.items[item_index] {
                CifItem::Loop(loop_) => loop_.value(self.row, position),
                _ => None,
            }
        } else {
            match self.table.block.items.get(position) {
                Some(CifItem::Pair(pair)) => pair.value.as_ref(),
                _ => None,
            }
        }
    }

    pub fn decoded(&self, column: usize) -> Option<String> {
        self.get(column).map(CifValue::decoded)
    }

    pub fn one_of(&self, primary: usize, fallback: usize) -> Result<&'a str, CifReadError> {
        // Gemmi❗✔️: bool has(size_t n) const { return tab.positions.at(n) >= 0; }
        // Gemmi❗✔️: bool has2(size_t n) const { return has(n) && !cif::is_null(operator[](n)); }
        // Gemmi❗✔️: const std::string& one_of(size_t n1, size_t n2) const {
        // Gemmi❗✔️:   static const std::string nul(1, '.');
        // Gemmi❗✔️:   if (has2(n1))
        // Gemmi❗✔️:     return operator[](n1);
        // Gemmi❗✔️:   if (has(n2))
        // Gemmi❗✔️:     return operator[](n2);
        // Gemmi❗✔️:   return nul;
        // Gemmi❗✔️: }
        // Gemmi✔️✔️: explicit Item(std::string&& t)
        // Gemmi✔️✔️:   : type{ItemType::Pair}, pair{{std::move(t), std::string()}} {}
        // Gemmi❗✔️: inline std::string& Table::Row::operator[](size_t n) {
        // Gemmi❗✔️:   int pos = tab.positions[n];
        // Gemmi❗✔️:   if (Loop* loop = tab.get_loop()) {
        // Gemmi❗✔️:     if (row_index == -1) // tags
        // Gemmi❗✔️:       return loop->tags[pos];
        // Gemmi❗✔️:     return loop->values[loop->width() * row_index + pos];
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   return tab.bloc.items[pos].pair[row_index == -1 ? 0 : 1];
        // Gemmi❗✔️: }
        // Behavior review: preserves the source's primary/null/fallback order,
        // including an empty source pair value and static raw dot when a
        // requested optional position is absent. C++ vector::at failures are
        // represented as typed OutOfRange errors rather than C++ exceptions.
        // Complexity review: constant-time position and row access with no
        // allocation on successful source-value branches; only an error
        // constructs an owned message.
        let out_of_range = |column: usize| {
            CifReadError::new(
                CifReadErrorKind::OutOfRange,
                "cif",
                self.table.block.line,
                1,
                format!(
                    "CIF row column index {column} is outside width {}",
                    self.table.width()
                ),
            )
        };

        let primary_position = self
            .table
            .positions
            .get(primary)
            .ok_or_else(|| out_of_range(primary))?;
        if primary_position.is_some() {
            // A pair without an item value is initialized upstream with an
            // empty second string; CifPair retains that state as None.
            let raw = self.get(primary).map(CifValue::raw).unwrap_or("");
            if !cif_is_null(raw) {
                return Ok(raw);
            }
        }

        let fallback_position = self
            .table
            .positions
            .get(fallback)
            .ok_or_else(|| out_of_range(fallback))?;
        if fallback_position.is_some() {
            return Ok(self.get(fallback).map(CifValue::raw).unwrap_or(""));
        }

        Ok(".")
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum TokenKind {
    Data(String),
    Global,
    Loop,
    Save(String),
    Stop,
    Tag(String),
    Value(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Token {
    kind: TokenKind,
    line: usize,
    column: usize,
}

struct Lexer<'a> {
    input: &'a str,
    bytes: &'a [u8],
    source: &'a str,
    index: usize,
    line: usize,
    column: usize,
}

impl<'a> Lexer<'a> {
    fn new(input: &'a str, source: &'a str) -> Self {
        Self {
            input,
            bytes: input.as_bytes(),
            source,
            index: 0,
            line: 1,
            column: 1,
        }
    }

    fn tokenize(mut self) -> Result<Vec<Token>, CifReadError> {
        let mut tokens = Vec::new();
        while self.skip_space_and_comments() {
            tokens.push(self.next_token()?);
        }
        Ok(tokens)
    }

    fn skip_space_and_comments(&mut self) -> bool {
        loop {
            while self.index < self.bytes.len() && char_table(self.bytes[self.index]) == 2 {
                self.bump_ascii();
            }
            if self.index < self.bytes.len() && self.bytes[self.index] == b'#' {
                while self.index < self.bytes.len() && self.bytes[self.index] != b'\n' {
                    self.bump_ascii();
                }
                continue;
            }
            return self.index < self.bytes.len();
        }
    }

    fn next_token(&mut self) -> Result<Token, CifReadError> {
        // Gemmi❗❌: template<typename Q>
        // Gemmi❗❌: struct endq : pegtl::seq<Q, pegtl::at<pegtl::sor<
        // Gemmi❗❌:                                           pegtl::one<' ','\n','\r','\t','#'>,
        // Gemmi❗❌:                                           pegtl::eof>>> {};
        // Gemmi❗❌: template<typename Q>
        // Gemmi❗❌: struct quoted_tail : pegtl::until<endq<Q>, pegtl::not_one<'\n'>> {};
        // Gemmi❗❌: template<typename Q>
        // Gemmi❗❌: struct quoted : pegtl::if_must<Q, quoted_tail<Q>> {};
        // Gemmi❗❌: struct singlequoted : quoted<pegtl::one<'\''>> {};
        // Gemmi❗❌: struct doublequoted : quoted<pegtl::one<'"'>> {};
        // Gemmi❗❌: struct field_sep : pegtl::seq<pegtl::bol, pegtl::one<';'>> {};
        // Gemmi❗❌: struct textfield : pegtl::if_must<field_sep, pegtl::until<field_sep>> {};
        // Gemmi❗❌: struct unquoted : pegtl::seq<pegtl::not_at<keyword>,
        // Gemmi❗❌:                                pegtl::not_at<pegtl::one<'_','$','#'>>,
        // Gemmi❗❌:                                pegtl::plus<nonblank_ch>> {};
        // Behavior review: token dispatch follows these source grammar rules;
        // complete lexical and UTF-8 position parity remains unverified.
        // Complexity review: this token-vector design retains an extra token
        // allocation layer versus PEGTL's direct document-building actions.
        let line = self.line;
        let column = self.column;
        let byte = self.bytes[self.index];
        if byte == b';' && column == 1 {
            return self.text_field(line, column);
        }
        if byte == b'\'' || byte == b'"' {
            return self.quoted(byte, line, column);
        }
        if !(b'!'..=b'~').contains(&byte) {
            return Err(self.error("unexpected non-ASCII byte outside a quoted/text value"));
        }
        let start = self.index;
        while self.index < self.bytes.len() && char_table(self.bytes[self.index]) != 2 {
            if !(b'!'..=b'~').contains(&self.bytes[self.index]) {
                return Err(self.error("unexpected non-ASCII byte in unquoted value"));
            }
            self.bump_ascii();
        }
        let raw = &self.input[start..self.index];
        let lower = raw.to_ascii_lowercase();
        let kind = if let Some(name) = lower.strip_prefix("data_") {
            TokenKind::Data(raw[raw.len() - name.len()..].to_owned())
        } else if lower == "global_" {
            TokenKind::Global
        } else if lower == "loop_" {
            TokenKind::Loop
        } else if let Some(name) = lower.strip_prefix("save_") {
            TokenKind::Save(raw[raw.len() - name.len()..].to_owned())
        } else if lower == "stop_" {
            TokenKind::Stop
        } else if raw.starts_with('_') && raw.len() > 1 {
            TokenKind::Tag(raw.to_owned())
        } else if raw.starts_with('_') || raw.starts_with('$') || raw.starts_with('#') {
            return Err(CifReadError::new(
                CifReadErrorKind::Syntax,
                self.source,
                line,
                column.saturating_sub(1),
                format!("invalid unquoted CIF token {raw:?}"),
            ));
        } else {
            TokenKind::Value(raw.to_owned())
        };
        Ok(Token { kind, line, column })
    }

    fn quoted(&mut self, quote: u8, line: usize, column: usize) -> Result<Token, CifReadError> {
        // Gemmi❗✔️: template<typename Q>
        // Gemmi❗✔️: struct endq : pegtl::seq<Q, pegtl::at<pegtl::sor<
        // Gemmi❗✔️:                                           pegtl::one<' ','\n','\r','\t','#'>,
        // Gemmi❗✔️:                                           pegtl::eof>>> {};
        // Gemmi❗✔️: template<typename Q>
        // Gemmi❗✔️: struct quoted_tail : pegtl::until<endq<Q>, pegtl::not_one<'\n'>> {};
        // Gemmi❗✔️: template<typename Q>
        // Gemmi❗✔️: struct quoted : pegtl::if_must<Q, quoted_tail<Q>> {};
        // Behavior review: matching-quote delimiters and the pinned current-
        // cursor failure location are reproduced; broader UTF-8 location
        // parity remains unverified. Complexity review: one forward scan and
        // one owned raw-token allocation, with no rescanning.
        let start = self.index;
        self.bump_ascii();
        loop {
            if self.index >= self.bytes.len() {
                return Err(self.error(if quote == b'\'' {
                    "unterminated 'string'"
                } else {
                    "unterminated \"string\""
                }));
            }
            let current = self.bytes[self.index];
            if current == b'\n' {
                return Err(self.error(if quote == b'\'' {
                    "unterminated 'string'"
                } else {
                    "unterminated \"string\""
                }));
            }
            if current == quote {
                let next = self.bytes.get(self.index + 1).copied();
                if next.is_none()
                    || next.is_some_and(|byte| matches!(byte, b' ' | b'\t' | b'\r' | b'\n' | b'#'))
                {
                    self.bump_ascii();
                    let raw = self.input[start..self.index].to_owned();
                    return Ok(Token {
                        kind: TokenKind::Value(raw),
                        line,
                        column,
                    });
                }
            }
            self.bump_utf8();
        }
    }

    fn text_field(&mut self, line: usize, column: usize) -> Result<Token, CifReadError> {
        // Gemmi❗✔️: struct field_sep : pegtl::seq<pegtl::bol, pegtl::one<';'>> {};
        // Gemmi❗✔️: struct textfield : pegtl::if_must<field_sep, pegtl::until<field_sep>> {};
        // Behavior review: BOL field delimiters and raw delimiter retention
        // are implemented; malformed EOF location remains outside the cases
        // established by this step. Complexity review: one forward scan and
        // one raw-value allocation.
        let start = self.index;
        self.bump_ascii();
        while self.index < self.bytes.len() {
            if self.column == 1 && self.bytes[self.index] == b';' {
                self.bump_ascii();
                return Ok(Token {
                    kind: TokenKind::Value(self.input[start..self.index].to_owned()),
                    line,
                    column,
                });
            }
            self.bump_utf8();
        }
        Err(CifReadError::new(
            CifReadErrorKind::Syntax,
            self.source,
            line,
            column.saturating_sub(1),
            "unterminated text field",
        ))
    }

    fn bump_ascii(&mut self) {
        let byte = self.bytes[self.index];
        self.index += 1;
        if byte == b'\n' {
            self.line += 1;
            self.column = 1;
        } else {
            self.column += 1;
        }
    }

    fn bump_utf8(&mut self) {
        let byte = self.bytes[self.index];
        if byte.is_ascii() {
            self.bump_ascii();
            return;
        }
        let length = self.input[self.index..]
            .chars()
            .next()
            .expect("index is in input")
            .len_utf8();
        self.index += length;
        self.column += 1;
    }

    fn error(&self, message: impl Into<String>) -> CifReadError {
        // Gemmi❗✔️: template<typename Rule> struct Errors : public pegtl::normal<Rule> {
        // Gemmi❗✔️:   template<typename Input, typename ... States>
        // Gemmi❗✔️:   static void raise(const Input& in, States&& ...) {
        // Gemmi❗✔️:     throw pegtl::parse_error(error_message<Rule>()
        // Gemmi❗✔️:                            //+ " matching " + pegtl::internal::demangle<Rule>()
        // Gemmi❗✔️:                              , in);
        // Gemmi❗✔️:   }
        // Gemmi❗✔️: };
        // Behavior review: lexer columns are converted from the internal
        // one-based cursor to PEGTL's zero-based source column. Complexity
        // review: constant-time error construction.
        CifReadError::new(
            CifReadErrorKind::Syntax,
            self.source,
            self.line,
            self.column.saturating_sub(1),
            message,
        )
    }
}

struct Parser<'a> {
    tokens: Vec<Token>,
    source: &'a str,
    index: usize,
    syntax_check_only: bool,
}

impl<'a> Parser<'a> {
    fn new(tokens: Vec<Token>, source: &'a str, syntax_check_only: bool) -> Self {
        Self {
            tokens,
            source,
            index: 0,
            syntax_check_only,
        }
    }

    fn parse(mut self) -> Result<CifDocument, CifReadError> {
        // Gemmi✔️✔️: struct datablock : pegtl::seq<datablockheading, ws_or_eof,
        // Gemmi✔️✔️:                              pegtl::star<pegtl::sor<dataitem, loop, frame>>> {};
        // Gemmi✔️✔️: struct content : pegtl::plus<datablock> {};
        // Gemmi✔️✔️: struct file : pegtl::seq<pegtl::opt<whitespace>,
        // Gemmi✔️✔️:                            pegtl::if_must<pegtl::not_at<pegtl::eof>,
        // Gemmi✔️✔️:                                           content, pegtl::eof>> {};
        // Behavior review: the parser requires one or more data/global blocks
        // and consumes the complete token stream. Complexity review: each token
        // is consumed once; nested frame parsing uses the same vector cursor.
        if self.tokens.is_empty() {
            return Err(self.error_at(1, 1, "expected block header (data_)"));
        }
        let mut blocks = Vec::new();
        while self.index < self.tokens.len() {
            let heading = self.tokens[self.index].clone();
            let name = match heading.kind {
                TokenKind::Data(name) => {
                    if name.is_empty() {
                        " ".to_owned()
                    } else {
                        name
                    }
                }
                TokenKind::Global => String::new(),
                _ => {
                    return Err(self.error_token(&heading, "expected block header (data_)"));
                }
            };
            self.index += 1;
            let items = self.parse_items(false)?;
            blocks.push(CifBlock {
                name,
                items,
                line: heading.line,
            });
        }
        Ok(CifDocument {
            source: self.source.to_owned(),
            blocks,
        })
    }

    fn parse_items(&mut self, in_frame: bool) -> Result<Vec<CifItem>, CifReadError> {
        let mut items = Vec::new();
        while let Some(token) = self.tokens.get(self.index).cloned() {
            match token.kind {
                TokenKind::Data(_) | TokenKind::Global if !in_frame => break,
                TokenKind::Save(ref name) if in_frame && name.is_empty() => {
                    self.index += 1;
                    return Ok(items);
                }
                TokenKind::Tag(tag) => {
                    self.index += 1;
                    let value = match self.tokens.get(self.index) {
                        Some(next) if matches!(next.kind, TokenKind::Value(_)) => {
                            let next = self.tokens[self.index].clone();
                            self.index += 1;
                            let TokenKind::Value(raw) = next.kind else {
                                unreachable!();
                            };
                            Some(CifValue::new(raw, next.line, next.column))
                        }
                        Some(next) if next.column == 1 => {
                            if self.syntax_check_only {
                                return Err(self.error_token(next, "tag without value"));
                            }
                            None
                        }
                        None => {
                            if self.syntax_check_only {
                                return Err(self.error_at(
                                    token.line,
                                    token.column,
                                    "tag without value",
                                ));
                            }
                            None
                        }
                        Some(next) => {
                            return Err(self.error_token(next, "parse error after tag"));
                        }
                    };
                    items.push(CifItem::Pair(CifPair {
                        tag,
                        value,
                        line: token.line,
                    }));
                }
                TokenKind::Loop => {
                    self.index += 1;
                    items.push(CifItem::Loop(self.parse_loop(token.line)?));
                }
                TokenKind::Save(name) if !name.is_empty() && !in_frame => {
                    self.index += 1;
                    let frame_items = self.parse_items(true)?;
                    items.push(CifItem::Frame(CifBlock {
                        name,
                        items: frame_items,
                        line: token.line,
                    }));
                }
                TokenKind::Save(_) if in_frame => {
                    return Err(self.error_token(&token, "nested or unnamed save_ frame"));
                }
                _ => return Err(self.error_token(&token, "parse error")),
            }
        }
        if in_frame {
            return Err(self.error_at(
                self.tokens.last().map_or(1, |token| token.line),
                1,
                "unterminated save_ frame",
            ));
        }
        Ok(items)
    }

    fn parse_loop(&mut self, line: usize) -> Result<CifLoop, CifReadError> {
        // Gemmi✔️✔️: struct loop : pegtl::if_must<str_loop,
        // Gemmi✔️✔️:                   whitespace,
        // Gemmi✔️✔️:                   pegtl::plus<pegtl::seq<loop_tag, whitespace, pegtl::discard>>,
        // Gemmi✔️✔️:                   pegtl::sor<pegtl::plus<pegtl::seq<loop_value, ws_or_eof,
        // Gemmi✔️✔️:                                                     pegtl::discard>>,
        // Gemmi✔️✔️:                              pegtl::at<pegtl::sor<keyword, pegtl::eof>>>,
        // Gemmi✔️✔️:                   loop_end> {};
        // Gemmi✔️✔️: if (loop.values.size() % loop.tags.size() != 0)
        // Gemmi✔️✔️:   throw pegtl::parse_error(
        // Gemmi✔️✔️:       "Wrong number of values in loop " + loop.common_prefix() + "*", in);
        // Behavior review: at least one tag, source-tolerated empty loops,
        // optional stop_, and read-path width validation are reproduced.
        // Complexity review: tags and values are appended once and width is O(1).
        let mut tags = Vec::new();
        while let Some(Token {
            kind: TokenKind::Tag(tag),
            ..
        }) = self.tokens.get(self.index)
        {
            tags.push(tag.clone());
            self.index += 1;
        }
        if tags.is_empty() {
            let token = self.tokens.get(self.index);
            return Err(match token {
                Some(token) => self.error_token(token, "loop_ without tags"),
                None => self.error_at(line, 1, "loop_ without tags"),
            });
        }
        let mut values = Vec::new();
        while let Some(token) = self.tokens.get(self.index).cloned() {
            match token.kind {
                TokenKind::Value(raw) => {
                    values.push(CifValue::new(raw, token.line, token.column));
                    self.index += 1;
                }
                TokenKind::Stop => {
                    self.index += 1;
                    break;
                }
                _ => break,
            }
        }
        let loop_ = CifLoop { tags, values, line };
        if !self.syntax_check_only && loop_.values.len() % loop_.tags.len() != 0 {
            return Err(self.error_at(
                line,
                1,
                format!("Wrong number of values in loop {}*", loop_.common_prefix()),
            ));
        }
        Ok(loop_)
    }

    fn error_token(&self, token: &Token, message: impl Into<String>) -> CifReadError {
        self.error_at(token.line, token.column, message)
    }

    fn error_at(&self, line: usize, column: usize, message: impl Into<String>) -> CifReadError {
        // Gemmi❗✔️: template<typename Rule> struct Errors : public pegtl::normal<Rule> {
        // Gemmi❗✔️:   template<typename Input, typename ... States>
        // Gemmi❗✔️:   static void raise(const Input& in, States&& ...) {
        // Gemmi❗✔️:     throw pegtl::parse_error(error_message<Rule>()
        // Gemmi❗✔️:                            //+ " matching " + pegtl::internal::demangle<Rule>()
        // Gemmi❗✔️:                              , in);
        // Gemmi❗✔️:   }
        // Gemmi❗✔️: };
        // Behavior review: parser errors use the pinned input cursor location;
        // parser-owned cursor columns are normalized from one-based to zero-
        // based at this boundary. Complexity review: constant-time conversion.
        CifReadError::new(
            CifReadErrorKind::Syntax,
            self.source,
            line,
            column.saturating_sub(1),
            message,
        )
    }
}

/// Parse an in-memory CIF document with the selected Gemmi check level.
pub fn read_cif_document(
    text: &str,
    source_name: &str,
    check_level: CifCheckLevel,
) -> Result<CifDocument, CifReadError> {
    // Gemmi✔️✔️: template<typename Input> Document read_input(Input&& in, int check_level=1) {
    // Gemmi✔️✔️:   Document doc;
    // Gemmi✔️✔️:   doc.source = in.source();
    // Gemmi✔️✔️:   parse_input(doc, in);
    // Gemmi✔️✔️:   if (check_level > 0) {
    // Gemmi✔️✔️:     check_for_missing_values(doc);
    // Gemmi✔️✔️:     check_for_duplicates(doc);
    // Gemmi✔️✔️:     if (check_level > 1) {
    // Gemmi✔️✔️:       for (const cif::Block& block : doc.blocks) {
    // Gemmi✔️✔️:         if (block.name == " ")
    // Gemmi✔️✔️:           fail(doc.source + ": missing block name (bare data_)");
    // Gemmi✔️✔️:         check_empty_loops(block, doc.source);
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return doc;
    // Gemmi✔️✔️: }
    // Behavior review: validation layers are applied only after a fully built
    // document. Complexity review: parsing plus each enabled validation pass is
    // linear in input/items/tags and does not clone document payloads.
    let tokens = Lexer::new(text, source_name).tokenize()?;
    let document = Parser::new(tokens, source_name, false).parse()?;
    match check_level {
        CifCheckLevel::Syntax => {}
        CifCheckLevel::Default => {
            check_missing_values(&document)?;
            check_duplicates(&document)?;
        }
        CifCheckLevel::Strict => {
            check_missing_values(&document)?;
            check_duplicates(&document)?;
            for block in &document.blocks {
                if block.name == " " {
                    return Err(CifReadError::new(
                        CifReadErrorKind::InvalidValue,
                        source_name,
                        block.line,
                        1,
                        "missing block name (bare data_)",
                    ));
                }
                check_empty_loops(block, source_name)?;
            }
        }
    }
    Ok(document)
}

/// Check only the pinned CIF grammar plus the source's tag-without-value rule.
pub fn check_cif_syntax(text: &str, source_name: &str) -> Result<(), CifReadError> {
    // Gemmi✔️✔️: template<typename Input> bool try_parse(Input&& in, std::string* msg) {
    // Gemmi✔️✔️:   try {
    // Gemmi✔️✔️:     return pegtl::parse<rules::file, CheckAction, Errors>(in);
    // Gemmi✔️✔️:   } catch (pegtl::parse_error& e) {
    // Gemmi✔️✔️:     if (msg)
    // Gemmi✔️✔️:       *msg = e.what();
    // Gemmi✔️✔️:     return false;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // Behavior review: this entry applies grammar and missing-value action but
    // deliberately omits document duplicate/empty-loop checks and loop-width
    // Action validation, matching the separate CheckAction source path.
    // Complexity review: it shares the linear tokenizer/parser without cloning
    // a second input buffer.
    let tokens = Lexer::new(text, source_name).tokenize()?;
    Parser::new(tokens, source_name, true).parse()?;
    Ok(())
}

fn check_missing_values(document: &CifDocument) -> Result<(), CifReadError> {
    fn block_check(block: &CifBlock, source: &str) -> Result<(), CifReadError> {
        // Gemmi✔️✔️: for (const Item& item : block.items) {
        // Gemmi✔️✔️:   if (item.type == ItemType::Pair) {
        // Gemmi✔️✔️:     if (item.pair[1].empty())
        // Gemmi✔️✔️:       cif_fail(source, block, item, item.pair[0] + " has no value");
        // Gemmi✔️✔️:   } else if (item.type == ItemType::Frame) {
        // Gemmi✔️✔️:     check_for_missing_values_in_block(item.frame, source);
        // Gemmi✔️✔️:   }
        // Gemmi✔️✔️: }
        // Behavior and complexity review: recursive item-order traversal is
        // source-equivalent and linear with stack depth equal to frame depth.
        for item in &block.items {
            match item {
                CifItem::Pair(pair) if pair.value.is_none() => {
                    return Err(CifReadError::new(
                        CifReadErrorKind::MissingValue,
                        source,
                        pair.line,
                        1,
                        format!("{} has no value", pair.tag),
                    ));
                }
                CifItem::Frame(frame) => block_check(frame, source)?,
                _ => {}
            }
        }
        Ok(())
    }
    for block in &document.blocks {
        block_check(block, &document.source)?;
    }
    Ok(())
}

fn check_duplicates(document: &CifDocument) -> Result<(), CifReadError> {
    // Gemmi✔️✔️: std::unordered_set<std::string> names;
    // Gemmi✔️✔️: for (const Block& block : d.blocks) {
    // Gemmi✔️✔️:   bool ok = names.insert(gemmi::to_lower(block.name)).second;
    // Gemmi✔️✔️:   if (!ok && !block.name.empty())
    // Gemmi✔️✔️:     fail(d.source + ": duplicate block name: ", block.name);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: for (const Block& block : d.blocks) {
    // Gemmi✔️✔️:   names.clear();
    // Gemmi✔️✔️:   frame_names.clear();
    // Gemmi✔️✔️:   for (const Item& item : block.items) {
    // Gemmi✔️✔️:     if (item.type == ItemType::Pair) {
    // Gemmi✔️✔️:       bool ok = names.insert(gemmi::to_lower(item.pair[0])).second;
    // Gemmi✔️✔️:       if (!ok)
    // Gemmi✔️✔️:         cif_fail(d.source, block, item, "duplicate tag " + item.pair[0]);
    // Gemmi✔️✔️:     } else if (item.type == ItemType::Loop) {
    // Gemmi✔️✔️:       for (const std::string& t : item.loop.tags) {
    // Gemmi✔️✔️:         bool ok = names.insert(gemmi::to_lower(t)).second;
    // Gemmi✔️✔️:         if (!ok)
    // Gemmi✔️✔️:           cif_fail(d.source, block, item, "duplicate tag " + t);
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:     } else if (item.type == ItemType::Frame) {
    // Gemmi✔️✔️:       bool ok = frame_names.insert(gemmi::to_lower(item.frame.name)).second;
    // Gemmi✔️✔️:       if (!ok)
    // Gemmi✔️✔️:         cif_fail(d.source, block, item, "duplicate save_" + item.frame.name);
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // Behavior review: duplicate scopes and global-block exception match source.
    // Complexity review: HashSet insertion is the same expected O(1) strategy.
    let mut block_names = HashSet::new();
    for block in &document.blocks {
        if !block_names.insert(block.name.to_ascii_lowercase()) && !block.name.is_empty() {
            return Err(CifReadError::new(
                CifReadErrorKind::DuplicateName,
                &document.source,
                block.line,
                1,
                format!("duplicate block name: {}", block.name),
            ));
        }
    }
    for block in &document.blocks {
        let mut names = HashSet::new();
        let mut frame_names = HashSet::new();
        for item in &block.items {
            match item {
                CifItem::Pair(pair) => {
                    if !names.insert(pair.tag.to_ascii_lowercase()) {
                        return Err(CifReadError::new(
                            CifReadErrorKind::DuplicateName,
                            &document.source,
                            pair.line,
                            1,
                            format!("duplicate tag {}", pair.tag),
                        ));
                    }
                }
                CifItem::Loop(loop_) => {
                    for tag in &loop_.tags {
                        if !names.insert(tag.to_ascii_lowercase()) {
                            return Err(CifReadError::new(
                                CifReadErrorKind::DuplicateName,
                                &document.source,
                                loop_.line,
                                1,
                                format!("duplicate tag {tag}"),
                            ));
                        }
                    }
                }
                CifItem::Frame(frame) => {
                    if !frame_names.insert(frame.name.to_ascii_lowercase()) {
                        return Err(CifReadError::new(
                            CifReadErrorKind::DuplicateName,
                            &document.source,
                            frame.line,
                            1,
                            format!("duplicate save_{}", frame.name),
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

fn check_empty_loops(block: &CifBlock, source: &str) -> Result<(), CifReadError> {
    // Gemmi✔️✔️: for (const cif::Item& item : block.items) {
    // Gemmi✔️✔️:   if (item.type == cif::ItemType::Loop) {
    // Gemmi✔️✔️:     if (item.loop.values.empty() && !item.loop.tags.empty())
    // Gemmi✔️✔️:       cif_fail(source, block, item, "empty loop with " + item.loop.tags[0]);
    // Gemmi✔️✔️:   } else if (item.type == cif::ItemType::Frame) {
    // Gemmi✔️✔️:     check_empty_loops(item.frame, source);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // Behavior and complexity review: exact recursive condition, one linear pass.
    for item in &block.items {
        match item {
            CifItem::Loop(loop_) if loop_.values.is_empty() && !loop_.tags.is_empty() => {
                return Err(CifReadError::new(
                    CifReadErrorKind::InvalidLoop,
                    source,
                    loop_.line,
                    1,
                    format!("empty loop with {}", loop_.tags[0]),
                ));
            }
            CifItem::Frame(frame) => check_empty_loops(frame, source)?,
            _ => {}
        }
    }
    Ok(())
}

/// Return true only for raw CIF null tokens `.` and `?`.
pub fn cif_is_null(value: &str) -> bool {
    value.len() == 1 && matches!(value.as_bytes()[0], b'.' | b'?')
}

/// Decode one raw CIF value using Gemmi's `as_string` rules.
pub fn cif_as_string(value: &str) -> String {
    // Gemmi✔️✔️: inline std::string as_string(const std::string& value) {
    // Gemmi✔️✔️:   if (value.empty() || is_null(value))
    // Gemmi✔️✔️:     return "";
    // Gemmi✔️✔️:   if (value[0] == '"' || value[0] == '\'')
    // Gemmi✔️✔️:     return std::string(value.begin() + 1, value.end() - 1);
    // Gemmi✔️✔️:   if (value[0] == ';' && value.size() > 2 && *(value.end() - 2) == '\n') {
    // Gemmi✔️✔️:     bool crlf = *(value.end() - 3) == '\r';
    // Gemmi✔️✔️:     return std::string(value.begin() + 1, value.end() - (crlf ? 3 : 2));
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return value;
    // Gemmi✔️✔️: }
    // Behavior review: byte delimiters are ASCII and therefore safe UTF-8
    // boundaries. Complexity review: the source also returns an owned string.
    if value.is_empty() || cif_is_null(value) {
        return String::new();
    }
    let bytes = value.as_bytes();
    if matches!(bytes[0], b'\'' | b'"') {
        return value[1..value.len() - 1].to_owned();
    }
    if bytes[0] == b';' && value.len() > 2 && bytes[value.len() - 2] == b'\n' {
        let crlf = value.len() >= 3 && bytes[value.len() - 3] == b'\r';
        return value[1..value.len() - if crlf { 3 } else { 2 }].to_owned();
    }
    value.to_owned()
}

pub fn cif_as_char(value: &str, null: char) -> Result<char, CifReadError> {
    // Gemmi✔️✔️: inline char as_char(const std::string& value, char null) {
    // Gemmi✔️✔️:   if (is_null(value))
    // Gemmi✔️✔️:     return null;
    // Gemmi✔️✔️:   if (value.size() < 2)
    // Gemmi✔️✔️:     return value[0];
    // Gemmi✔️✔️:   const std::string s = as_string(value);
    // Gemmi✔️✔️:   if (s.size() < 2)
    // Gemmi✔️✔️:     return s[0];
    // Gemmi✔️✔️:   fail("Not a single character: " + value);
    // Gemmi✔️✔️: }
    // Behavior review: the length checks are byte-sized like std::string;
    // an empty string yields the C++ terminator NUL, and multi-byte UTF-8 is
    // rejected. Invalid byte sequences cannot enter through Rust `&str`.
    // Complexity review: the decode is one owned linear copy, matching the
    // source's temporary `std::string`; all other checks are constant-time.
    if cif_is_null(value) {
        return Ok(null);
    }
    if value.len() < 2 {
        return Ok(char::from(
            value.as_bytes().first().copied().unwrap_or_default(),
        ));
    }
    let decoded = cif_as_string(value);
    if decoded.len() < 2 {
        Ok(char::from(
            decoded.as_bytes().first().copied().unwrap_or_default(),
        ))
    } else {
        Err(value_error(value, "Not a single character"))
    }
}

pub fn cif_as_i32(value: &str) -> Result<i32, CifReadError> {
    // Gemmi✔️✔️: inline int as_int(const std::string& str) {
    // Gemmi✔️✔️:   return string_to_int(str, true);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: while ((length == 0 || i < length) && is_space(p[i])) ++i;
    // Gemmi✔️✔️: if (p[i] == '-') { mult = 1; ++i; } else if (p[i] == '+') { ++i; }
    // Gemmi✔️✔️: for (; (length == 0 || i < length) && is_digit(p[i]); ++i) {
    // Gemmi✔️✔️:   n = n * 10 - (p[i] - '0'); has_digits = true;
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: if (checked) { while (is_space(p[i])) ++i;
    // Gemmi✔️✔️:   if (!has_digits || p[i] != '\0') throw std::invalid_argument(...); }
    // Behavior review: defined C-locale whitespace/sign/full-consumption behavior
    // is exact. Source signed-overflow is undefined; Rust returns OutOfRange.
    // Complexity review: one linear byte pass, no temporary trimmed string.
    let bytes = value.as_bytes();
    let mut index = 0;
    while bytes.get(index).is_some_and(|byte| is_c_space(*byte)) {
        index += 1;
    }
    let negative = match bytes.get(index) {
        Some(b'-') => {
            index += 1;
            true
        }
        Some(b'+') => {
            index += 1;
            false
        }
        _ => false,
    };
    let mut has_digits = false;
    let mut magnitude: i64 = 0;
    while let Some(byte @ b'0'..=b'9') = bytes.get(index) {
        has_digits = true;
        magnitude = magnitude
            .checked_mul(10)
            .and_then(|number| number.checked_add(i64::from(*byte - b'0')))
            .ok_or_else(|| range_error(value, "integer outside supported i32 range"))?;
        index += 1;
    }
    while bytes.get(index).is_some_and(|byte| is_c_space(*byte)) {
        index += 1;
    }
    if !has_digits || index != bytes.len() {
        return Err(value_error(value, "not an integer"));
    }
    let signed = if negative { -magnitude } else { magnitude };
    i32::try_from(signed).map_err(|_| range_error(value, "integer outside supported i32 range"))
}

pub fn cif_as_f64(value: &str, null: f64) -> Result<f64, CifReadError> {
    // Gemmi❗❌: inline double as_number(const std::string& s, double nan=NAN) {
    // Gemmi❗❌:   const char* start = s.data();
    // Gemmi❗❌:   const char* end = s.data() + s.size();
    // Gemmi❗❌:   if (*start == '+')
    // Gemmi❗❌:     ++start;
    // Gemmi❗❌:   char f = start[int(*start == '-')] | 0x20;
    // Gemmi❗❌:   if (f == 'i' || f == 'n')
    // Gemmi❗❌:     return nan;
    // Gemmi❗❌:   double d;
    // Gemmi❗❌:   auto result = fast_float::from_chars(start, end, d);
    // Gemmi❗❌:   if (result.ec != std::errc())
    // Gemmi❗❌:     return nan;
    // Gemmi❗❌:   if (*result.ptr == '(') {
    // Gemmi❗❌:     const char* p = result.ptr + 1;
    // Gemmi❗❌:     while (*p >= '0' && *p <= '9')
    // Gemmi❗❌:       ++p;
    // Gemmi❗❌:     if (*p == ')')
    // Gemmi❗❌:       result.ptr = p + 1;
    // Gemmi❗❌:   }
    // Gemmi❗❌:   return result.ptr == end ? d : nan;
    // Gemmi❗❌: }
    // fast_float❗❌:   // C++17 20.19.3.(7.1) explicitly forbids '+' sign here
    // fast_float❗❌:   if ((*p == UC('-')) || (uint64_t(fmt & chars_format::allow_leading_plus) &&
    // fast_float❗❌:                           !basic_json_fmt && *p == UC('+'))) {
    // fast_float❗❌:     ++p;
    // fast_float❗❌:   }
    // Behavior review: source-defined failures return the caller's null value,
    // and zero-digit uncertainty is accepted. The range-error condition below
    // mirrors fast_float's nonzero-decimal-mantissa-to-rounded-zero branch;
    // exact-zero mantissas and representable subnormals remain successful.
    // Gemmi removes at most one `+`; its default fast_float general format
    // rejects a second leading `+`. Successful nonzero decimal conversion
    // still uses Rust `str::parse`, so exact decimal bit parity remains
    // unproven outside the tested underflow boundary.
    // Complexity review: Rust scans for `(` before parsing the numeric
    // prefix, adding an extra linear pass versus the source's direct
    // fast_float scan, and scans the significand once more to distinguish an
    // exact zero token from nonzero input rounded to zero; the sign checks
    // themselves are constant-time.
    if cif_is_null(value) {
        return Ok(null);
    }
    let numeric = value.strip_prefix('+').unwrap_or(value);
    if numeric.starts_with('+') {
        return Ok(null);
    }
    let signless = numeric.strip_prefix('-').unwrap_or(numeric);
    if signless
        .as_bytes()
        .first()
        .is_some_and(|byte| matches!(byte.to_ascii_lowercase(), b'i' | b'n'))
    {
        return Ok(null);
    }
    let (number_text, suffix) = match numeric.find('(') {
        Some(open) => (&numeric[..open], &numeric[open..]),
        None => (numeric, ""),
    };
    if !suffix.is_empty()
        && !(suffix.starts_with('(')
            && suffix.ends_with(')')
            && suffix.len() >= 2
            && suffix[1..suffix.len() - 1]
                .bytes()
                .all(|byte| byte.is_ascii_digit()))
    {
        return Ok(null);
    }
    let Ok(parsed) = number_text.parse::<f64>() else {
        return Ok(null);
    };
    if parsed.is_finite() {
        // fast_float❗❌:   to_float(pns.negative, am, value);
        // fast_float❗❌:   // Test for over/underflow.
        // fast_float❗❌:   if ((pns.mantissa != 0 && am.mantissa == 0 && am.power2 == 0) ||
        // fast_float❗❌:       am.power2 == binary_format<T>::infinite_power()) {
        // fast_float❗❌:     answer.ec = std::errc::result_out_of_range;
        // fast_float❗❌:   }
        // fast_float❗❌:   return answer;
        // The input grammar above is a successfully parsed decimal. For that
        // grammar, fast_float's accumulated decimal mantissa is zero exactly
        // when every significand digit is zero; its long-significand path
        // skips leading zeroes and retains a nonzero significant prefix.
        // Thus this check detects source range errors that round nonzero
        // values to signed zero without an exponent threshold or rounding
        // approximation. Exact zero (including signed zero) stays successful.
        if parsed == 0.0 {
            let significand_end = number_text
                .find(|character| matches!(character, 'e' | 'E'))
                .unwrap_or(number_text.len());
            let decimal_mantissa_is_nonzero = number_text[..significand_end]
                .bytes()
                .any(|byte| matches!(byte, b'1'..=b'9'));
            if decimal_mantissa_is_nonzero {
                return Ok(null);
            }
        }
        Ok(parsed)
    } else {
        Ok(null)
    }
}

/// Quote one decoded value using Gemmi's CIF quoting priority.
pub fn quote_cif_value(mut value: String) -> String {
    // Gemmi✔️✔️: inline std::string quote(std::string v) {
    // Gemmi✔️✔️:   if (std::all_of(v.begin(), v.end(), [](char c) { return char_table(c) == 1; })
    // Gemmi✔️✔️:       && !v.empty() && !is_null(v))
    // Gemmi✔️✔️:     return v;
    // Gemmi✔️✔️:   char q = ';';
    // Gemmi✔️✔️:   if (std::memchr(v.c_str(), '\n', v.size()) == nullptr) {
    // Gemmi✔️✔️:     if (std::memchr(v.c_str(), '\'', v.size()) == nullptr) q = '\'';
    // Gemmi✔️✔️:     else if (std::memchr(v.c_str(), '"', v.size()) == nullptr) q = '"';
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   v.insert(v.begin(), q);
    // Gemmi✔️✔️:   if (q == ';') v += '\n';
    // Gemmi✔️✔️:   v += q;
    // Gemmi✔️✔️:   return v;
    // Gemmi✔️✔️: }
    // Behavior review: exact ordinary/null and quote-priority branches.
    // Complexity review: linear scans and one owned output, as in source.
    if !value.is_empty() && !cif_is_null(&value) && value.bytes().all(|byte| char_table(byte) == 1)
    {
        return value;
    }
    let quote = if value.contains('\n') {
        ';'
    } else if !value.contains('\'') {
        '\''
    } else if !value.contains('"') {
        '"'
    } else {
        ';'
    };
    value.insert(0, quote);
    if quote == ';' {
        value.push('\n');
    }
    value.push(quote);
    value
}

/// Gemmi `to_str(double)`, backed by the fixed `% .9g` semantics without locale.
pub fn format_cif_f64(value: f64) -> String {
    // Gemmi❗❌: inline std::string to_str(double d) {
    // Gemmi❗❌:   char buf[24];
    // Gemmi❗❌:   int len = sprintf_z(buf, "%.9g", d);
    // Gemmi❗❌:   return std::string(buf, len > 0 ? len : 0);
    // Gemmi❗❌: }
    // Behavior review: the fixed profile uses the pinned bundled-stb
    // significant-digit conversion; B22 establishes the recorded cases, not
    // a blanket all-input claim. Complexity review: Rust allocates digit and
    // output Strings instead of using the source's fixed stack buffer.
    format_general(value, 9)
}

/// Gemmi `to_str(float)`, preserving the source f32 input rounding first.
pub fn format_cif_f32(value: f32) -> String {
    // Gemmi❗❌: inline std::string to_str(float d) {
    // Gemmi❗❌:   char buf[16];
    // Gemmi❗❌:   int len = sprintf_z(buf, "%.6g", d);
    // Gemmi❗❌:   return std::string(buf, len > 0 ? len : 0);
    // Gemmi❗❌: }
    // Behavior review: f32 is promoted exactly to f64 before the pinned
    // six-significant-digit conversion; sampled oracle coverage is not an
    // all-input claim. Complexity review: Rust allocates Strings versus the
    // source's fixed stack buffer.
    format_general(f64::from(value), 6)
}

pub fn format_cif_f64_precision<const PREC: usize>(value: f64) -> String {
    // Gemmi❗❌: static_assert(Prec >= 0 && Prec < 7, "unsupported precision");
    // Gemmi❗❌: char buf[16];
    // Gemmi❗❌: int len = d > -1e8 && d < 1e8 ? sprintf_z(buf, "%.*f", Prec, d)
    // Gemmi❗❌:                               : sprintf_z(buf, "%g", d);
    // Gemmi❗❌: return std::string(buf, len > 0 ? len : 0);
    // Behavior review: strict interval dispatch is source-shaped; selected
    // fixed/general cases are covered, but not every f64. Complexity review:
    // fixed arithmetic plus allocated Rust output differs from stack format.
    assert!(PREC < 7, "unsupported precision");
    if value > -1e8 && value < 1e8 {
        format_fixed_cif_precision(value, PREC)
    } else {
        format_general(value, 6)
    }
}

pub fn format_cif_i32(value: i32) -> String {
    value.to_string()
}

pub fn format_cif_usize(value: usize) -> String {
    value.to_string()
}

fn format_fixed_cif_precision(value: f64, precision: usize) -> String {
    // Gemmi❗✔️:   if ((frac_digits < 24)) {
    // Gemmi❗✔️:      stbsp__uint32 dg = 1;
    // Gemmi❗✔️:      if ((stbsp__uint64)bits >= stbsp__powten[9])
    // Gemmi❗✔️:         dg = 10;
    // Gemmi❗✔️:      while ((stbsp__uint64)bits >= stbsp__powten[dg]) {
    // Gemmi❗✔️:         ++dg;
    // Gemmi❗✔️:         if (dg == 20)
    // Gemmi❗✔️:            goto noround;
    // Gemmi❗✔️:      }
    // Gemmi❗✔️:      if (frac_digits < dg) {
    // Gemmi❗✔️:         stbsp__uint64 r;
    // Gemmi❗✔️:         // add 0.5 at the right position and round
    // Gemmi❗✔️:         e = dg - frac_digits;
    // Gemmi❗✔️:         if ((stbsp__uint32)e >= 24)
    // Gemmi❗✔️:            goto noround;
    // Gemmi❗✔️:         r = stbsp__powten[e];
    // Gemmi❗✔️:         bits = bits + (r / 2);
    // Gemmi❗✔️:         if ((stbsp__uint64)bits >= stbsp__powten[dg])
    // Gemmi❗✔️:            ++tens;
    // Gemmi❗✔️:         bits /= r;
    // Gemmi❗✔️:      }
    // Gemmi❗✔️:   noround:;
    // Gemmi❗✔️:   }
    // Behavior review: reproduce the source's positive-magnitude half-unit
    // addition using the exact binary value scaled by 10^precision; a negative
    // sign is applied afterward, so ties round away from zero. This is covered
    // by the fixed Gemmi profile matrix, not claimed as an all-input port of
    // stb's double-double decimal conversion.
    // Complexity review: constant-width u128 arithmetic (precision < 7) plus
    // output-sized decimal formatting; no data-dependent search or table scan.
    let bits = value.abs().to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1_u64 << 52) - 1);
    let (significand, binary_exponent) = if exponent == 0 {
        (fraction, -1074)
    } else {
        ((1_u64 << 52) | fraction, exponent - 1023 - 52)
    };

    let numerator = u128::from(significand) * 5_u128.pow(precision as u32);
    let shift = binary_exponent + precision as i32;
    let rounded = if shift >= 0 {
        // This helper is called only for |value| < 1e8 and precision < 7,
        // which bounds this left shift to a value below 1e14.
        numerator << shift as u32
    } else {
        let denominator_shift = (-shift) as u32;
        if denominator_shift >= 128 {
            0
        } else {
            let denominator = 1_u128 << denominator_shift;
            let quotient = numerator / denominator;
            let remainder = numerator % denominator;
            quotient + u128::from(remainder >= denominator / 2 && denominator_shift > 0)
        }
    };

    let negative = value.is_sign_negative();
    let sign = if negative { "-" } else { "" };
    if precision == 0 {
        return format!("{sign}{rounded}");
    }

    let scale = 10_u128.pow(precision as u32);
    let integer = rounded / scale;
    let fractional = rounded % scale;
    format!("{sign}{integer}.{fractional:0width$}", width = precision)
}

fn stb_double_double_product(x: f64, y: f64) -> (f64, f64) {
    // Gemmi❗✔️: #define stbsp__ddmulthi(oh, ol, xh, yh)                            \
    // Gemmi❗✔️:    {                                                               \
    // Gemmi❗✔️:       double ahi = 0, alo, bhi = 0, blo;                           \
    // Gemmi❗✔️:       stbsp__int64 bt;                                             \
    // Gemmi❗✔️:       oh = xh * yh;                                                \
    // Gemmi❗✔️:       STBSP__COPYFP(bt, xh);                                       \
    // Gemmi❗✔️:       bt &= ((~(stbsp__uint64)0) << 27);                           \
    // Gemmi❗✔️:       STBSP__COPYFP(ahi, bt);                                      \
    // Gemmi❗✔️:       alo = xh - ahi;                                              \
    // Gemmi❗✔️:       STBSP__COPYFP(bt, yh);                                       \
    // Gemmi❗✔️:       bt &= ((~(stbsp__uint64)0) << 27);                           \
    // Gemmi❗✔️:       STBSP__COPYFP(bhi, bt);                                      \
    // Gemmi❗✔️:       blo = yh - bhi;                                              \
    // Gemmi❗✔️:       ol = ((ahi * bhi - oh) + ahi * blo + alo * bhi) + alo * blo; \
    // Gemmi❗✔️:    }
    // Gemmi❗❌: #define STBSP__COPYFP(dest, src)                   \
    // Gemmi❗❌:    {                                               \
    // Gemmi❗❌:       int cn;                                      \
    // Gemmi❗❌:       for (cn = 0; cn < 8; cn++)                   \
    // Gemmi❗❌:          ((char *)&dest)[cn] = ((char *)&src)[cn]; \
    // Gemmi❗❌:    }
    // Behavior review: `to_bits`/`from_bits` performs the same binary64
    // truncation; arithmetic order follows the macro. Complexity review:
    // fixed arithmetic, but tuple and scalar operations replace the macro's
    // in-place temporaries.
    let high = x * y;
    let x_high = f64::from_bits(x.to_bits() & (!0_u64 << 27));
    let x_low = x - x_high;
    let y_high = f64::from_bits(y.to_bits() & (!0_u64 << 27));
    let y_low = y - y_high;
    let low = ((x_high * y_high - high) + x_high * y_low + x_low * y_high) + x_low * y_low;
    (high, low)
}

fn stb_double_double_renormalize(high: f64, low: f64) -> (f64, f64) {
    // Gemmi❗✔️: #define stbsp__ddrenorm(oh, ol) \
    // Gemmi❗✔️:    {                            \
    // Gemmi❗✔️:       double s;                 \
    // Gemmi❗✔️:       s = oh + ol;              \
    // Gemmi❗✔️:       ol = ol - (s - oh);       \
    // Gemmi❗✔️:       oh = s;                   \
    // Gemmi❗✔️:    }
    // Behavior review: preserve the source operation order. Complexity
    // review: constant scalar arithmetic, no allocation.
    let sum = high + low;
    let residual = low - (sum - high);
    (sum, residual)
}

fn stb_double_double_to_i64(high: f64, low: f64) -> i64 {
    // Gemmi❗✔️: #define stbsp__ddtoS64(ob, xh, xl)          \
    // Gemmi❗✔️:    {                                        \
    // Gemmi❗✔️:       double ahi = 0, alo, vh, t;           \
    // Gemmi❗✔️:       ob = (stbsp__int64)xh;                \
    // Gemmi❗✔️:       vh = (double)ob;                      \
    // Gemmi❗✔️:       ahi = (xh - vh);                      \
    // Gemmi❗✔️:       t = (ahi - xh);                       \
    // Gemmi❗✔️:       alo = (xh - (ahi - t)) - (vh + t);    \
    // Gemmi❗✔️:       ob += (stbsp__int64)(ahi + alo + xl); \
    // Gemmi❗✔️:    }
    // Behavior review: the formatter's normalized finite range keeps the
    // conversion in signed-64 range for the supported source-defined path.
    // Complexity review: constant arithmetic and one checked-domain cast.
    let mut result = high as i64;
    let value = result as f64;
    let high_residual = high - value;
    let correction = high_residual - high;
    let low_residual = (high - (high_residual - correction)) - (value + correction);
    result += (high_residual + low_residual + low) as i64;
    result
}

fn stb_raise_to_power10(value: f64, power: i32) -> (f64, f64) {
    // Gemmi❗✔️: static double const stbsp__bot[23] = {
    // Gemmi❗❌:    1e+000, 1e+001, 1e+002, 1e+003, 1e+004, 1e+005, 1e+006, 1e+007, 1e+008, 1e+009, 1e+010, 1e+011,
    // Gemmi❗❌:    1e+012, 1e+013, 1e+014, 1e+015, 1e+016, 1e+017, 1e+018, 1e+019, 1e+020, 1e+021, 1e+022
    // Gemmi❗❌: };
    // Gemmi❗✔️: static double const stbsp__negbot[22] = {
    // Gemmi❗❌:    1e-001, 1e-002, 1e-003, 1e-004, 1e-005, 1e-006, 1e-007, 1e-008, 1e-009, 1e-010, 1e-011,
    // Gemmi❗❌:    1e-012, 1e-013, 1e-014, 1e-015, 1e-016, 1e-017, 1e-018, 1e-019, 1e-020, 1e-021, 1e-022
    // Gemmi❗❌: };
    // Gemmi❗✔️: static double const stbsp__negboterr[22] = {
    // Gemmi❗❌:    -5.551115123125783e-018,  -2.0816681711721684e-019, -2.0816681711721686e-020, -4.7921736023859299e-021, -8.1803053914031305e-022, 4.5251888174113741e-023,
    // Gemmi❗❌:    4.5251888174113739e-024,  -2.0922560830128471e-025, -6.2281591457779853e-026, -3.6432197315497743e-027, 6.0503030718060191e-028,  2.0113352370744385e-029,
    // Gemmi❗❌:    -3.0373745563400371e-030, 1.1806906454401013e-032,  -7.7705399876661076e-032, 2.0902213275965398e-033,  -7.1542424054621921e-034, -7.1542424054621926e-035,
    // Gemmi❗❌:    2.4754073164739869e-036,  5.4846728545790429e-037,  9.2462547772103625e-038,  -4.8596774326570872e-039
    // Gemmi❗❌: };
    // Gemmi❗✔️: static double const stbsp__top[13] = {
    // Gemmi❗❌:    1e+023, 1e+046, 1e+069, 1e+092, 1e+115, 1e+138, 1e+161, 1e+184, 1e+207, 1e+230, 1e+253, 1e+276, 1e+299
    // Gemmi❗❌: };
    // Gemmi❗✔️: static double const stbsp__negtop[13] = {
    // Gemmi❗❌:    1e-023, 1e-046, 1e-069, 1e-092, 1e-115, 1e-138, 1e-161, 1e-184, 1e-207, 1e-230, 1e-253, 1e-276, 1e-299
    // Gemmi❗❌: };
    // Gemmi❗✔️: static double const stbsp__toperr[13] = {
    // Gemmi❗❌:    8388608,
    // Gemmi❗❌:    6.8601809640529717e+028,
    // Gemmi❗❌:    -7.253143638152921e+052,
    // Gemmi❗❌:    -4.3377296974619174e+075,
    // Gemmi❗❌:    -1.5559416129466825e+098,
    // Gemmi❗❌:    -3.2841562489204913e+121,
    // Gemmi❗❌:    -3.7745893248228135e+144,
    // Gemmi❗❌:    -1.7356668416969134e+167,
    // Gemmi❗❌:    -3.8893577551088374e+190,
    // Gemmi❗❌:    -9.9566444326005119e+213,
    // Gemmi❗❌:    6.3641293062232429e+236,
    // Gemmi❗❌:    -5.2069140800249813e+259,
    // Gemmi❗❌:    -5.2504760255204387e+282
    // Gemmi❗❌: };
    // Gemmi❗✔️: static double const stbsp__negtoperr[13] = {
    // Gemmi❗❌:    3.9565301985100693e-040,  -2.299904345391321e-063,  3.6506201437945798e-086,  1.1875228833981544e-109,
    // Gemmi❗❌:    -5.0644902316928607e-132, -6.7156837247865426e-155, -2.812077463003139e-178,  -5.7778912386589953e-201,
    // Gemmi❗❌:    7.4997100559334532e-224,  -4.6439668915134491e-247, -6.3691100762962136e-270, -9.436808465446358e-293,
    // Gemmi❗❌:    8.0970921678014997e-317
    // Gemmi❗❌: };
    // Gemmi❗✔️: static double const stbsp__powten[20] = {
    // Gemmi❗❌:    1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000,
    // Gemmi❗❌:    10000000000ULL, 100000000000ULL, 1000000000000ULL, 10000000000000ULL, 100000000000000ULL,
    // Gemmi❗❌:    1000000000000000ULL, 10000000000000000ULL, 100000000000000000ULL, 1000000000000000000ULL,
    // Gemmi❗❌:    10000000000000000000ULL
    // Gemmi❗❌: };
    // Gemmi❗✔️: #define stbsp__tento19th (1000000000000000000ULL)
    // Gemmi❗❗: #define stbsp__ddmultlo(oh, ol, xh, xl, yh, yl) ol = ol + (xh * yl + xl * yh);
    // Gemmi❗❗: #define stbsp__ddmultlos(oh, ol, xh, yl) ol = ol + (xh * yl);
    // Gemmi❗✔️: static void stbsp__raise_to_power10(double *ohi, double *olo, double d, stbsp__int32 power)
    // Gemmi❗❌: {
    // Gemmi❗❌:    double ph, pl;
    // Gemmi❗❌:    if ((power >= 0) && (power <= 22)) {
    // Gemmi❗❌:       stbsp__ddmulthi(ph, pl, d, stbsp__bot[power]);
    // Gemmi❗❗:    } else {
    // Gemmi❗❗:       stbsp__int32 e, et, eb;
    // Gemmi❗❗:       double p2h, p2l;
    // Gemmi❗❗:       e = power;
    // Gemmi❗❗:       if (power < 0)
    // Gemmi❗❗:          e = -e;
    // Gemmi❗❗:       et = (e * 0x2c9) >> 14; /* %23 */
    // Gemmi❗❗:       if (et > 13)
    // Gemmi❗❗:          et = 13;
    // Gemmi❗❗:       eb = e - (et * 23);
    // Gemmi❗❗:       ph = d;
    // Gemmi❗❗:       pl = 0.0;
    // Gemmi❗❗:       if (power < 0) {
    // Gemmi❗❗:          if (eb) {
    // Gemmi❗❗:             --eb;
    // Gemmi❗❗:             stbsp__ddmulthi(ph, pl, d, stbsp__negbot[eb]);
    // Gemmi❗❗:             stbsp__ddmultlos(ph, pl, d, stbsp__negboterr[eb]);
    // Gemmi❗❗:          }
    // Gemmi❗❗:          if (et) {
    // Gemmi❗❗:             stbsp__ddrenorm(ph, pl);
    // Gemmi❗❗:             --et;
    // Gemmi❗❗:             stbsp__ddmulthi(p2h, p2l, ph, stbsp__negtop[et]);
    // Gemmi❗❗:             stbsp__ddmultlo(p2h, p2l, ph, pl, stbsp__negtop[et], stbsp__negtoperr[et]);
    // Gemmi❗❗:             ph = p2h;
    // Gemmi❗❗:             pl = p2l;
    // Gemmi❗❗:          }
    // Gemmi❗❗:       } else {
    // Gemmi❗❗:          if (eb) {
    // Gemmi❗❗:             e = eb;
    // Gemmi❗❗:             if (eb > 22)
    // Gemmi❗❗:                eb = 22;
    // Gemmi❗❗:             e -= eb;
    // Gemmi❗❗:             stbsp__ddmulthi(ph, pl, d, stbsp__bot[eb]);
    // Gemmi❗❗:             if (e) {
    // Gemmi❗❗:                stbsp__ddrenorm(ph, pl);
    // Gemmi❗❗:                stbsp__ddmulthi(p2h, p2l, ph, stbsp__bot[e]);
    // Gemmi❗❗:                stbsp__ddmultlos(p2h, p2l, stbsp__bot[e], pl);
    // Gemmi❗❗:                ph = p2h;
    // Gemmi❗❗:                pl = p2l;
    // Gemmi❗❗:             }
    // Gemmi❗❗:          }
    // Gemmi❗❗:          if (et) {
    // Gemmi❗❗:             stbsp__ddrenorm(ph, pl);
    // Gemmi❗❗:             --et;
    // Gemmi❗❗:             stbsp__ddmulthi(p2h, p2l, ph, stbsp__top[et]);
    // Gemmi❗❗:             stbsp__ddmultlo(p2h, p2l, ph, pl, stbsp__top[et], stbsp__toperr[et]);
    // Gemmi❗❗:             ph = p2h;
    // Gemmi❗❗:             pl = p2l;
    // Gemmi❗❗:          }
    // Gemmi❗❗:       }
    // Gemmi❗❗:    }
    // Gemmi❗❗:    stbsp__ddrenorm(ph, pl);
    // Gemmi❗❗:    *ohi = ph;
    // Gemmi❗❗:    *olo = pl;
    // Gemmi❗❗: }
    // Behavior review: this is the selected signed power-of-ten closure for
    // the fixed Gemmi profile; actual binary64 inputs and outputs are checked
    // separately below. Complexity review: bounded table indexing and a
    // constant number of double-double products, with no heap allocation.
    const BOT: [f64; 23] = [
        1.0e0, 1.0e1, 1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8, 1.0e9, 1.0e10, 1.0e11,
        1.0e12, 1.0e13, 1.0e14, 1.0e15, 1.0e16, 1.0e17, 1.0e18, 1.0e19, 1.0e20, 1.0e21, 1.0e22,
    ];
    const NEGBOT: [f64; 22] = [
        1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10, 1.0e-11,
        1.0e-12, 1.0e-13, 1.0e-14, 1.0e-15, 1.0e-16, 1.0e-17, 1.0e-18, 1.0e-19, 1.0e-20, 1.0e-21,
        1.0e-22,
    ];
    const NEGBOTERR: [f64; 22] = [
        -5.551115123125783e-18,
        -2.0816681711721684e-19,
        -2.0816681711721686e-20,
        -4.7921736023859299e-21,
        -8.1803053914031305e-22,
        4.5251888174113741e-23,
        4.5251888174113739e-24,
        -2.0922560830128471e-25,
        -6.2281591457779853e-26,
        -3.6432197315497743e-27,
        6.0503030718060191e-28,
        2.0113352370744385e-29,
        -3.0373745563400371e-30,
        1.1806906454401013e-32,
        -7.7705399876661076e-32,
        2.0902213275965398e-33,
        -7.1542424054621921e-34,
        -7.1542424054621926e-35,
        2.4754073164739869e-36,
        5.4846728545790429e-37,
        9.2462547772103625e-38,
        -4.8596774326570872e-39,
    ];
    const TOP: [f64; 13] = [
        1.0e23, 1.0e46, 1.0e69, 1.0e92, 1.0e115, 1.0e138, 1.0e161, 1.0e184, 1.0e207, 1.0e230,
        1.0e253, 1.0e276, 1.0e299,
    ];
    const NEGTOP: [f64; 13] = [
        1.0e-23, 1.0e-46, 1.0e-69, 1.0e-92, 1.0e-115, 1.0e-138, 1.0e-161, 1.0e-184, 1.0e-207,
        1.0e-230, 1.0e-253, 1.0e-276, 1.0e-299,
    ];
    const TOPERR: [f64; 13] = [
        8_388_608.0,
        6.8601809640529717e28,
        -7.253143638152921e52,
        -4.3377296974619174e75,
        -1.5559416129466825e98,
        -3.2841562489204913e121,
        -3.7745893248228135e144,
        -1.7356668416969134e167,
        -3.8893577551088374e190,
        -9.9566444326005119e213,
        6.3641293062232429e236,
        -5.2069140800249813e259,
        -5.2504760255204387e282,
    ];
    const NEGTOPERR: [f64; 13] = [
        3.9565301985100693e-40,
        -2.299904345391321e-63,
        3.6506201437945798e-86,
        1.1875228833981544e-109,
        -5.0644902316928607e-132,
        -6.7156837247865426e-155,
        -2.812077463003139e-178,
        -5.7778912386589953e-201,
        7.4997100559334532e-224,
        -4.6439668915134491e-247,
        -6.3691100762962136e-270,
        -9.436808465446358e-293,
        8.0970921678014997e-317,
    ];

    let (mut high, mut low);
    if (0..=22).contains(&power) {
        (high, low) = stb_double_double_product(value, BOT[power as usize]);
    } else {
        let mut e = power.abs();
        let mut et = (e * 0x2c9) >> 14;
        if et > 13 {
            et = 13;
        }
        let mut eb = e - et * 23;
        high = value;
        low = 0.0;
        if power < 0 {
            if eb != 0 {
                eb -= 1;
                (high, low) = stb_double_double_product(value, NEGBOT[eb as usize]);
                low += value * NEGBOTERR[eb as usize];
            }
            if et != 0 {
                (high, low) = stb_double_double_renormalize(high, low);
                et -= 1;
                let (mut product_high, mut product_low) =
                    stb_double_double_product(high, NEGTOP[et as usize]);
                product_low += high * NEGTOPERR[et as usize] + low * NEGTOP[et as usize];
                high = product_high;
                low = product_low;
            }
        } else {
            if eb != 0 {
                e = eb;
                if eb > 22 {
                    eb = 22;
                }
                e -= eb;
                (high, low) = stb_double_double_product(value, BOT[eb as usize]);
                if e != 0 {
                    (high, low) = stb_double_double_renormalize(high, low);
                    let (mut product_high, mut product_low) =
                        stb_double_double_product(high, BOT[e as usize]);
                    product_low += BOT[e as usize] * low;
                    high = product_high;
                    low = product_low;
                }
            }
            if et != 0 {
                (high, low) = stb_double_double_renormalize(high, low);
                et -= 1;
                let (mut product_high, mut product_low) =
                    stb_double_double_product(high, TOP[et as usize]);
                product_low += high * TOPERR[et as usize] + low * TOP[et as usize];
                high = product_high;
                low = product_low;
            }
        }
    }
    stb_double_double_renormalize(high, low)
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum StbGeneralValue {
    Special {
        text: &'static str,
        negative: bool,
    },
    Finite {
        digits: String,
        decimal_position: i32,
        negative: bool,
    },
}

fn general_sign(negative: bool) -> &'static str {
    // Gemmi❗✔️: static void stbsp__lead_sign(stbsp__uint32 fl, char *sign)
    // Gemmi❗✔️: {
    // Gemmi❗✔️:    sign[0] = 0;
    // Gemmi❗✔️:    if (fl & STBSP__NEGATIVE) {
    // Gemmi❗✔️:       sign[0] = 1;
    // Gemmi❗✔️:       sign[1] = '-';
    // Gemmi❗✔️:    } else if (fl & STBSP__LEADINGSPACE) {
    // Gemmi❗✔️:       sign[0] = 1;
    // Gemmi❗✔️:       sign[1] = ' ';
    // Gemmi❗✔️:    } else if (fl & STBSP__LEADINGPLUS) {
    // Gemmi❗✔️:       sign[0] = 1;
    // Gemmi❗✔️:       sign[1] = '+';
    // Gemmi❗✔️:    }
    // Gemmi❗✔️: }
    // Behavior review: the fixed sprintf_z %.ng callers provide no
    // leading-space or leading-plus flags, so the sign bit is the only
    // reachable non-empty prefix. Complexity review: static borrowed output
    // avoids allocation and is constant-time, matching the source helper.
    if negative { "-" } else { "" }
}

fn stb_real_to_str(value: f64, mut fraction_digits: u32) -> StbGeneralValue {
    // Gemmi❗❌: static stbsp__int32 stbsp__real_to_str(char const **start, stbsp__uint32 *len, char *out, stbsp__int32 *decimal_pos, double value, stbsp__uint32 frac_digits)
    // Gemmi❗❌: {
    // Gemmi❗❌:    double d;
    // Gemmi❗❌:    stbsp__int64 bits = 0;
    // Gemmi❗❌:    stbsp__int32 expo, e, ng, tens;
    // Gemmi❗❌:
    // Gemmi❗❌:    d = value;
    // Gemmi❗❌:    STBSP__COPYFP(bits, d);
    // Gemmi❗❌:    expo = (stbsp__int32)((bits >> 52) & 2047);
    // Gemmi❗❌:    ng = (stbsp__int32)((stbsp__uint64) bits >> 63);
    // Gemmi❗❌:    if (ng)
    // Gemmi❗❌:       d = -d;
    // Gemmi❗❌:
    // Gemmi❗❌:    if (expo == 2047) // is nan or inf?
    // Gemmi❗❌:    {
    // Gemmi❗❌:       *start = (bits & ((((stbsp__uint64)1) << 52) - 1)) ? "NaN" : "Inf";
    // Gemmi❗❌:       *decimal_pos = STBSP__SPECIAL;
    // Gemmi❗❌:       *len = 3;
    // Gemmi❗❌:       return ng;
    // Gemmi❗❌:    }
    // Gemmi❗❌:
    // Gemmi❗❌:    if (expo == 0) // is zero or denormal
    // Gemmi❗❌:    {
    // Gemmi❗❌:       if (((stbsp__uint64) bits << 1) == 0) // do zero
    // Gemmi❗❌:       {
    // Gemmi❗❌:          *decimal_pos = 1;
    // Gemmi❗❌:          *start = out;
    // Gemmi❗❌:          out[0] = '0';
    // Gemmi❗❌:          *len = 1;
    // Gemmi❗❌:          return ng;
    // Gemmi❗❌:       }
    // Gemmi❗❌:       // find the right expo for denormals
    // Gemmi❗❌:       {
    // Gemmi❗❌:          stbsp__int64 v = ((stbsp__uint64)1) << 51;
    // Gemmi❗❌:          while ((bits & v) == 0) {
    // Gemmi❗❌:             --expo;
    // Gemmi❗❌:             v >>= 1;
    // Gemmi❗❌:          }
    // Gemmi❗❌:       }
    // Gemmi❗❌:    }
    // Gemmi❗❌:
    // Gemmi❗❌:    // find the decimal exponent as well as the decimal bits of the value
    // Gemmi❗❌:    {
    // Gemmi❗❌:       double ph, pl;
    // Gemmi❗❌:
    // Gemmi❗❌:       // log10 estimate - very specifically tweaked to hit or undershoot by no more than 1 of log10 of all expos 1..2046
    // Gemmi❗❌:       tens = expo - 1023;
    // Gemmi❗❌:       tens = (tens < 0) ? ((tens * 617) / 2048) : (((tens * 1233) / 4096) + 1);
    // Gemmi❗❌:
    // Gemmi❗❌:       // move the significant bits into position and stick them into an int
    // Gemmi❗❌:       stbsp__raise_to_power10(&ph, &pl, d, 18 - tens);
    // Gemmi❗❌:
    // Gemmi❗❌:       // get full as much precision from double-double as possible
    // Gemmi❗❌:       stbsp__ddtoS64(bits, ph, pl);
    // Gemmi❗❌:
    // Gemmi❗❌:       // check if we undershot
    // Gemmi❗❌:       if (((stbsp__uint64)bits) >= stbsp__tento19th)
    // Gemmi❗❌:          ++tens;
    // Gemmi❗❌:    }
    // Gemmi❗❌:
    // Gemmi❗❌:    // now do the rounding in integer land
    // Gemmi❗❌:    frac_digits = (frac_digits & 0x80000000) ? ((frac_digits & 0x7ffffff) + 1) : (tens + frac_digits);
    // Gemmi❗❌:    if ((frac_digits < 24)) {
    // Gemmi❗❌:       stbsp__uint32 dg = 1;
    // Gemmi❗❌:       if ((stbsp__uint64)bits >= stbsp__powten[9])
    // Gemmi❗❌:          dg = 10;
    // Gemmi❗❌:       while ((stbsp__uint64)bits >= stbsp__powten[dg]) {
    // Gemmi❗❌:          ++dg;
    // Gemmi❗❌:          if (dg == 20)
    // Gemmi❗❌:             goto noround;
    // Gemmi❗❌:       }
    // Gemmi❗❌:       if (frac_digits < dg) {
    // Gemmi❗❌:          stbsp__uint64 r;
    // Gemmi❗❌:          // add 0.5 at the right position and round
    // Gemmi❗❌:          e = dg - frac_digits;
    // Gemmi❗❌:          if ((stbsp__uint32)e >= 24)
    // Gemmi❗❌:             goto noround;
    // Gemmi❗❌:          r = stbsp__powten[e];
    // Gemmi❗❌:          bits = bits + (r / 2);
    // Gemmi❗❌:          if ((stbsp__uint64)bits >= stbsp__powten[dg])
    // Gemmi❗❌:             ++tens;
    // Gemmi❗❌:          bits /= r;
    // Gemmi❗❌:       }
    // Gemmi❗❌:    noround:;
    // Gemmi❗❌:    }
    // Gemmi❗❌:
    // Gemmi❗❌:    // kill long trailing runs of zeros
    // Gemmi❗❌:    if (bits) {
    // Gemmi❗❌:       stbsp__uint32 n;
    // Gemmi❗❌:       for (;;) {
    // Gemmi❗❌:          if (bits <= 0xffffffff)
    // Gemmi❗❌:             break;
    // Gemmi❗❌:          if (bits % 1000)
    // Gemmi❗❌:             goto donez;
    // Gemmi❗❌:          bits /= 1000;
    // Gemmi❗❌:       }
    // Gemmi❗❌:       n = (stbsp__uint32)bits;
    // Gemmi❗❌:       while ((n % 1000) == 0)
    // Gemmi❗❌:          n /= 1000;
    // Gemmi❗❌:       bits = n;
    // Gemmi❗❌:    donez:;
    // Gemmi❗❌:    }
    // Gemmi❗❌:
    // Gemmi❗❌:    // convert to string
    // Gemmi❗❌:    out += 64;
    // Gemmi❗❌:    e = 0;
    // Gemmi❗❌:    for (;;) {
    // Gemmi❗❌:       stbsp__uint32 n;
    // Gemmi❗❌:       char *o = out - 8;
    // Gemmi❗❌:       // do the conversion in chunks of U32s (avoid most 64-bit divides, worth it, constant denomiators be damned)
    // Gemmi❗❌:       if (bits >= 100000000) {
    // Gemmi❗❌:          n = (stbsp__uint32)(bits % 100000000);
    // Gemmi❗❌:          bits /= 100000000;
    // Gemmi❗❌:       } else {
    // Gemmi❗❌:          n = (stbsp__uint32)bits;
    // Gemmi❗❌:          bits = 0;
    // Gemmi❗❌:       }
    // Gemmi❗❌:       while (n) {
    // Gemmi❗❌:          out -= 2;
    // Gemmi❗❌:          *(stbsp__uint16 *)out = *(stbsp__uint16 *)&stbsp__digitpair.pair[(n % 100) * 2];
    // Gemmi❗❌:          n /= 100;
    // Gemmi❗❌:          e += 2;
    // Gemmi❗❌:       }
    // Gemmi❗❌:       if (bits == 0) {
    // Gemmi❗❌:          if ((e) && (out[0] == '0')) {
    // Gemmi❗❌:             ++out;
    // Gemmi❗❌:             --e;
    // Gemmi❗❌:          }
    // Gemmi❗❌:          break;
    // Gemmi❗❌:       }
    // Gemmi❗❌:       while (out != o) {
    // Gemmi❗❌:          *--out = '0';
    // Gemmi❗❌:          ++e;
    // Gemmi❗❌:       }
    // Gemmi❗❌:    }
    // Gemmi❗❌:
    // Gemmi❗❌:    *decimal_pos = tens;
    // Gemmi❗❌:    *start = out;
    // Gemmi❗❌:    *len = e;
    // Gemmi❗❌:    return ng;
    // Gemmi❗❌: }
    // Behavior review: the finite decimal-significand and half-up steps
    // follow the pinned stb body; error/carry branches are covered by B22's
    // fixed oracle matrix. Complexity review: Rust materializes a `String`
    // for the digits instead of writing into stb's caller buffer, increasing
    // allocation; integer scaling and significant-digit work remain bounded.
    const POWTEN: [u64; 20] = [
        1,
        10,
        100,
        1_000,
        10_000,
        100_000,
        1_000_000,
        10_000_000,
        100_000_000,
        1_000_000_000,
        10_000_000_000,
        100_000_000_000,
        1_000_000_000_000,
        10_000_000_000_000,
        100_000_000_000_000,
        1_000_000_000_000_000,
        10_000_000_000_000_000,
        100_000_000_000_000_000,
        1_000_000_000_000_000_000,
        10_000_000_000_000_000_000,
    ];
    const TENTO19TH: u64 = 1_000_000_000_000_000_000;
    let raw_bits = value.to_bits();
    let mut exponent = ((raw_bits >> 52) & 2047) as i32;
    let negative = (raw_bits >> 63) != 0;
    let magnitude = if negative { -value } else { value };

    if exponent == 2047 {
        let text = if raw_bits & ((1_u64 << 52) - 1) != 0 {
            "NaN"
        } else {
            "Inf"
        };
        return StbGeneralValue::Special { text, negative };
    }

    if exponent == 0 {
        if raw_bits << 1 == 0 {
            return StbGeneralValue::Finite {
                digits: "0".to_owned(),
                decimal_position: 1,
                negative,
            };
        }
        let mut leading_bit = 1_u64 << 51;
        while raw_bits & leading_bit == 0 {
            exponent -= 1;
            leading_bit >>= 1;
        }
    }

    let mut tens = exponent - 1023;
    tens = if tens < 0 {
        (tens * 617) / 2048
    } else {
        (tens * 1233) / 4096 + 1
    };
    let (high, low) = stb_raise_to_power10(magnitude, 18 - tens);
    let mut bits = stb_double_double_to_i64(high, low) as u64;
    if bits >= TENTO19TH {
        tens += 1;
    }

    fraction_digits = if fraction_digits & 0x8000_0000 != 0 {
        (fraction_digits & 0x07ff_ffff) + 1
    } else {
        tens as u32 + fraction_digits
    };
    if fraction_digits < 24 {
        let mut digit_count = if bits >= POWTEN[9] { 10 } else { 1 };
        while bits >= POWTEN[digit_count as usize] {
            digit_count += 1;
            if digit_count == 20 {
                break;
            }
        }
        if fraction_digits < digit_count {
            let exponent = digit_count - fraction_digits;
            if exponent < 24 {
                let divisor = POWTEN[exponent as usize];
                bits += divisor / 2;
                if bits >= POWTEN[digit_count as usize] {
                    tens += 1;
                }
                bits /= divisor;
            }
        }
    }

    if bits != 0 {
        while bits > 0xffff_ffff && bits % 1000 == 0 {
            bits /= 1000;
        }
        while bits % 1000 == 0 {
            bits /= 1000;
        }
    }
    let digits = bits.to_string();
    StbGeneralValue::Finite {
        digits,
        decimal_position: tens,
        negative,
    }
}

fn format_general(value: f64, precision: usize) -> String {
    // Gemmi❗❌:       case 'G': // float
    // Gemmi❗❌:       case 'g': // float
    // Gemmi❗❌:          h = (f[0] == 'G') ? hexu : hex;
    // Gemmi❗❌:          fv = va_arg(va, double);
    // Gemmi❗❌:          if (pr == -1)
    // Gemmi❗❌:             pr = 6;
    // Gemmi❗❌:          else if (pr == 0)
    // Gemmi❗❌:             pr = 1; // default is 6
    // Gemmi❗❌:          // read the double into a string
    // Gemmi❗❌:          if (stbsp__real_to_str(&sn, &l, num, &dp, fv, (pr - 1) | 0x80000000))
    // Gemmi❗❌:             fl |= STBSP__NEGATIVE;
    // Gemmi❗❌:
    // Gemmi❗❌:          n = pr;
    // Gemmi❗❌:          // clamp the precision and delete extra zeros after clamp
    // Gemmi❗❌:          if (l > (stbsp__uint32)pr)
    // Gemmi❗❌:             l = pr;
    // Gemmi❗❌:          while ((l > 1) && (pr) && (sn[l - 1] == '0')) {
    // Gemmi❗❌:             --pr;
    // Gemmi❗❌:             --l;
    // Gemmi❗❌:          }
    // Gemmi❗❌:
    // Gemmi❗❌:          // should we use %e
    // Gemmi❗❌:          if ((dp <= -4) || (dp > (stbsp__int32)n)) {
    // Gemmi❗❌:             if (pr > (stbsp__int32)l)
    // Gemmi❗❌:                pr = l - 1;
    // Gemmi❗❌:             else if (pr)
    // Gemmi❗❌:                --pr; // when using %e, there is one digit before the decimal
    // Gemmi❗❌:             goto doexpfromg;
    // Gemmi❗❌:          }
    // Gemmi❗❌:          // this is the insane action to get the pr to match %g semantics for %f
    // Gemmi❗❌:          if (dp > 0) {
    // Gemmi❗❌:             pr = (dp < (stbsp__int32)l) ? l - dp : 0;
    // Gemmi❗❌:          } else {
    // Gemmi❗❌:             pr = -dp + ((pr > (stbsp__int32)l) ? (stbsp__int32) l : pr);
    // Gemmi❗❌:          }
    // Gemmi❗❌:          goto dofloatfromg;
    // Behavior review: this is the fixed Gemmi `%g` profile only (no width,
    // grouping, alternate-form, or locale flags); B22 is the frozen oracle
    // matrix. Complexity review: conversion uses bounded arithmetic but
    // materializes digit and output strings rather than stb's stack buffer.
    assert!(precision > 0);
    let significant_digits = precision as u32;
    let parts = stb_real_to_str(value, ((significant_digits - 1) | 0x8000_0000) as u32);
    let (mut digits, decimal_position, negative) = match parts {
        StbGeneralValue::Special { text, negative } => {
            return format!("{}{text}", general_sign(negative));
        }
        StbGeneralValue::Finite {
            digits,
            decimal_position,
            negative,
        } => (digits, decimal_position, negative),
    };
    let mut length = digits.len().min(precision);
    digits.truncate(length);
    let mut significant_precision = precision;
    while length > 1 && significant_precision > 0 && digits.as_bytes()[length - 1] == b'0' {
        significant_precision -= 1;
        length -= 1;
        digits.truncate(length);
    }

    if decimal_position <= -4 || decimal_position > precision as i32 {
        let fractional_precision = if significant_precision > length {
            length - 1
        } else if significant_precision > 0 {
            significant_precision - 1
        } else {
            0
        };
        format_general_exponent(
            &digits,
            decimal_position - 1,
            fractional_precision,
            negative,
        )
    } else {
        let fractional_precision = if decimal_position > 0 {
            if decimal_position < length as i32 {
                length - decimal_position as usize
            } else {
                0
            }
        } else {
            (-decimal_position) as usize + significant_precision.min(length)
        };
        format_general_fixed(&digits, decimal_position, fractional_precision, negative)
    }
}

fn format_general_exponent(
    digits: &str,
    exponent: i32,
    fractional_precision: usize,
    negative: bool,
) -> String {
    // Gemmi❗❌:       doexpfromg:
    // Gemmi❗❌:          tail[0] = 0;
    // Gemmi❗❌:          stbsp__lead_sign(fl, lead);
    // Gemmi❗❌:          if (dp == STBSP__SPECIAL) {
    // Gemmi❗❌:             s = (char *)sn;
    // Gemmi❗❌:             cs = 0;
    // Gemmi❗❌:             pr = 0;
    // Gemmi❗❌:             goto scopy;
    // Gemmi❗❌:          }
    // Gemmi❗❌:          s = num + 64;
    // Gemmi❗❌:          // handle leading chars
    // Gemmi❗❌:          *s++ = sn[0];
    // Gemmi❗❌:
    // Gemmi❗❌:          if (pr)
    // Gemmi❗❌:             *s++ = stbsp__period;
    // Gemmi❗❌:
    // Gemmi❗❌:          // handle after decimal
    // Gemmi❗❌:          if ((l - 1) > (stbsp__uint32)pr)
    // Gemmi❗❌:             l = pr + 1;
    // Gemmi❗❌:          for (n = 1; n < l; n++)
    // Gemmi❗❌:             *s++ = sn[n];
    // Gemmi❗❌:          // trailing zeros
    // Gemmi❗❌:          tz = pr - (l - 1);
    // Gemmi❗❌:          pr = 0;
    // Gemmi❗❌:          // dump expo
    // Gemmi❗❌:          tail[1] = h[0xe];
    // Gemmi❗❌:          dp -= 1;
    // Gemmi❗❌:          if (dp < 0) {
    // Gemmi❗❌:             tail[2] = '-';
    // Gemmi❗❌:             dp = -dp;
    // Gemmi❗❌:          } else
    // Gemmi❗❌:             tail[2] = '+';
    // Gemmi❗❌: #ifdef STB_SPRINTF_MSVC_MODE
    // Gemmi❗❌:          n = 5;
    // Gemmi❗❌: #else
    // Gemmi❗❌:          n = (dp >= 100) ? 5 : 4;
    // Gemmi❗❌: #endif
    // Gemmi❗❌:          tail[0] = (char)n;
    // Gemmi❗❌:          for (;;) {
    // Gemmi❗❌:             tail[n] = '0' + dp % 10;
    // Gemmi❗❌:             if (n <= 3)
    // Gemmi❗❌:                break;
    // Gemmi❗❌:             --n;
    // Gemmi❗❌:             dp /= 10;
    // Gemmi❗❌:          }
    // Gemmi❗❌:          cs = 1 + (3 << 24); // how many tens
    // Gemmi❗❌:          goto flt_lead;
    // Gemmi❗❌:       flt_lead:
    // Gemmi❗❌:          l = (stbsp__uint32)(s - (num + 64));
    // Gemmi❗❌:          s = num + 64;
    // Gemmi❗❌:          goto scopy;
    // Behavior review: fixed Gemmi uses lowercase `e`, a signed exponent, at
    // least two exponent digits, and a source sign bit for the mantissa.
    // Complexity review: output-sized formatting in one allocation; no scan
    // beyond the digit string.
    let mut output = String::with_capacity(digits.len() + 8);
    output.push_str(general_sign(negative));
    output.push(digits.as_bytes()[0] as char);
    let emitted_fraction = fractional_precision.min(digits.len().saturating_sub(1));
    if fractional_precision > 0 {
        output.push('.');
        output.push_str(&digits[1..1 + emitted_fraction]);
        output.extend(std::iter::repeat('0').take(fractional_precision - emitted_fraction));
    }
    output.push('e');
    if exponent < 0 {
        output.push('-');
    } else {
        output.push('+');
    }
    let exponent_digits = exponent.unsigned_abs().to_string();
    if exponent_digits.len() < 2 {
        output.push('0');
    }
    output.push_str(&exponent_digits);
    output
}

fn format_general_fixed(
    digits: &str,
    decimal_position: i32,
    fractional_precision: usize,
    negative: bool,
) -> String {
    // Gemmi❗❌:       dofloatfromg:
    // Gemmi❗❌:          tail[0] = 0;
    // Gemmi❗❌:          stbsp__lead_sign(fl, lead);
    // Gemmi❗❌:          if (dp == STBSP__SPECIAL) {
    // Gemmi❗❌:             s = (char *)sn;
    // Gemmi❗❌:             cs = 0;
    // Gemmi❗❌:             pr = 0;
    // Gemmi❗❌:             goto scopy;
    // Gemmi❗❌:          }
    // Gemmi❗❌:          s = num + 64;
    // Gemmi❗❌:
    // Gemmi❗❌:          // handle the three decimal varieties
    // Gemmi❗❌:          if (dp <= 0) {
    // Gemmi❗❌:             stbsp__int32 i;
    // Gemmi❗❌:             // handle 0.000*000xxxx
    // Gemmi❗❌:             *s++ = '0';
    // Gemmi❗❌:             if (pr)
    // Gemmi❗❌:                *s++ = stbsp__period;
    // Gemmi❗❌:             n = -dp;
    // Gemmi❗❌:             if ((stbsp__int32)n > pr)
    // Gemmi❗❌:                n = pr;
    // Gemmi❗❌:             i = n;
    // Gemmi❗❌:             while (i) {
    // Gemmi❗❌:                if ((((stbsp__uintptr)s) & 3) == 0)
    // Gemmi❗❌:                   break;
    // Gemmi❗❌:                *s++ = '0';
    // Gemmi❗❌:                --i;
    // Gemmi❗❌:             }
    // Gemmi❗❌:             while (i >= 4) {
    // Gemmi❗❌:                *(stbsp__uint32 *)s = 0x30303030;
    // Gemmi❗❌:                s += 4;
    // Gemmi❗❌:                i -= 4;
    // Gemmi❗❌:             }
    // Gemmi❗❌:             while (i) {
    // Gemmi❗❌:                *s++ = '0';
    // Gemmi❗❌:                --i;
    // Gemmi❗❌:             }
    // Gemmi❗❌:             if ((stbsp__int32)(l + n) > pr)
    // Gemmi❗❌:                l = pr - n;
    // Gemmi❗❌:             i = l;
    // Gemmi❗❌:             while (i) {
    // Gemmi❗❌:                *s++ = *sn++;
    // Gemmi❗❌:                --i;
    // Gemmi❗❌:             }
    // Gemmi❗❌:             tz = pr - (n + l);
    // Gemmi❗❌:             cs = 1 + (3 << 24); // how many tens did we write (for commas below)
    // Gemmi❗❌:          } else {
    // Gemmi❗❌:             cs = (fl & STBSP__TRIPLET_COMMA) ? ((600 - (stbsp__uint32)dp) % 3) : 0;
    // Gemmi❗❌:             if ((stbsp__uint32)dp >= l) {
    // Gemmi❗❌:                // handle xxxx000*000.0
    // Gemmi❗❌:                n = 0;
    // Gemmi❗❌:                for (;;) {
    // Gemmi❗❌:                   if ((fl & STBSP__TRIPLET_COMMA) && (++cs == 4)) {
    // Gemmi❗❌:                      cs = 0;
    // Gemmi❗❌:                      *s++ = stbsp__comma;
    // Gemmi❗❌:                   } else {
    // Gemmi❗❌:                      *s++ = sn[n];
    // Gemmi❗❌:                      ++n;
    // Gemmi❗❌:                      if (n >= l)
    // Gemmi❗❌:                         break;
    // Gemmi❗❌:                   }
    // Gemmi❗❌:                }
    // Gemmi❗❌:                if (n < (stbsp__uint32)dp) {
    // Gemmi❗❌:                   n = dp - n;
    // Gemmi❗❌:                   if ((fl & STBSP__TRIPLET_COMMA) == 0) {
    // Gemmi❗❌:                      while (n) {
    // Gemmi❗❌:                         if ((((stbsp__uintptr)s) & 3) == 0)
    // Gemmi❗❌:                            break;
    // Gemmi❗❌:                         *s++ = '0';
    // Gemmi❗❌:                         --n;
    // Gemmi❗❌:                      }
    // Gemmi❗❌:                      while (n >= 4) {
    // Gemmi❗❌:                         *(stbsp__uint32 *)s = 0x30303030;
    // Gemmi❗❌:                         s += 4;
    // Gemmi❗❌:                         n -= 4;
    // Gemmi❗❌:                      }
    // Gemmi❗❌:                   }
    // Gemmi❗❌:                   while (n) {
    // Gemmi❗❌:                      if ((fl & STBSP__TRIPLET_COMMA) && (++cs == 4)) {
    // Gemmi❗❌:                         cs = 0;
    // Gemmi❗❌:                         *s++ = stbsp__comma;
    // Gemmi❗❌:                      } else {
    // Gemmi❗❌:                         *s++ = '0';
    // Gemmi❗❌:                         --n;
    // Gemmi❗❌:                      }
    // Gemmi❗❌:                   }
    // Gemmi❗❌:                }
    // Gemmi❗❌:                cs = (int)(s - (num + 64)) + (3 << 24); // cs is how many tens
    // Gemmi❗❌:                if (pr) {
    // Gemmi❗❌:                   *s++ = stbsp__period;
    // Gemmi❗❌:                   tz = pr;
    // Gemmi❗❌:                }
    // Gemmi❗❌:             } else {
    // Gemmi❗❌:                // handle xxxxx.xxxx000*000
    // Gemmi❗❌:                n = 0;
    // Gemmi❗❌:                for (;;) {
    // Gemmi❗❌:                   if ((fl & STBSP__TRIPLET_COMMA) && (++cs == 4)) {
    // Gemmi❗❌:                      cs = 0;
    // Gemmi❗❌:                      *s++ = stbsp__comma;
    // Gemmi❗❌:                   } else {
    // Gemmi❗❌:                      *s++ = sn[n];
    // Gemmi❗❌:                      ++n;
    // Gemmi❗❌:                      if (n >= (stbsp__uint32)dp)
    // Gemmi❗❌:                         break;
    // Gemmi❗❌:                   }
    // Gemmi❗❌:                }
    // Gemmi❗❌:                cs = (int)(s - (num + 64)) + (3 << 24); // cs is how many tens
    // Gemmi❗❌:                if (pr)
    // Gemmi❗❌:                   *s++ = stbsp__period;
    // Gemmi❗❌:                if ((l - dp) > (stbsp__uint32)pr)
    // Gemmi❗❌:                   l = pr + dp;
    // Gemmi❗❌:                while (n < l) {
    // Gemmi❗❌:                   *s++ = sn[n];
    // Gemmi❗❌:                   ++n;
    // Gemmi❗❌:                }
    // Gemmi❗❌:                tz = pr - (l - dp);
    // Gemmi❗❌:             }
    // Gemmi❗❌:          }
    // Gemmi❗❌:          pr = 0;
    // Gemmi❗❌:
    // Gemmi❗❌:          // handle k,m,g,t
    // Gemmi❗❌:          if (fl & STBSP__METRIC_SUFFIX) {
    // Gemmi❗❌:             char idx;
    // Gemmi❗❌:             idx = 1;
    // Gemmi❗❌:             if (fl & STBSP__METRIC_NOSPACE)
    // Gemmi❗❌:                idx = 0;
    // Gemmi❗❌:             tail[0] = idx;
    // Gemmi❗❌:             tail[1] = ' ';
    // Gemmi❗❌:             {
    // Gemmi❗❌:                if (fl >> 24) {
    // Gemmi❗❌:                   if (fl & STBSP__METRIC_1024)
    // Gemmi❗❌:                      tail[idx + 1] = "_KMGT"[fl >> 24];
    // Gemmi❗❌:                   else
    // Gemmi❗❌:                      tail[idx + 1] = "_kMGT"[fl >> 24];
    // Gemmi❗❌:                   idx++;
    // Gemmi❗❌:                   if (fl & STBSP__METRIC_1024 && !(fl & STBSP__METRIC_JEDEC)) {
    // Gemmi❗❌:                      tail[idx + 1] = 'i';
    // Gemmi❗❌:                      idx++;
    // Gemmi❗❌:                   }
    // Gemmi❗❌:                   tail[0] = idx;
    // Gemmi❗❌:                }
    // Gemmi❗❌:             }
    // Gemmi❗❌:          };
    // Gemmi❗❌:
    // Gemmi❗❌:       flt_lead:
    // Gemmi❗❌:          // get the length that we copied
    // Gemmi❗❌:          l = (stbsp__uint32)(s - (num + 64));
    // Gemmi❗❌:          s = num + 64;
    // Gemmi❗❌:          goto scopy;
    // Behavior review: this is the fixed-point branch selected by `%g` for
    // the frozen profile; the source's alignment/vector writes are reduced
    // to equivalent string appends. Complexity review: linear in output size
    // with one output allocation, versus stb's stack buffer.
    let mut output = String::with_capacity(digits.len() + fractional_precision + 2);
    output.push_str(general_sign(negative));
    if decimal_position <= 0 {
        output.push('0');
        if fractional_precision > 0 {
            output.push('.');
        }
        let leading_zeros = (-decimal_position).max(0) as usize;
        let leading_zeros = leading_zeros.min(fractional_precision);
        output.extend(std::iter::repeat('0').take(leading_zeros));
        let digit_count = digits.len().min(fractional_precision - leading_zeros);
        output.push_str(&digits[..digit_count]);
        output.extend(
            std::iter::repeat('0').take(fractional_precision - leading_zeros - digit_count),
        );
    } else if decimal_position as usize >= digits.len() {
        output.push_str(digits);
        output.extend(std::iter::repeat('0').take(decimal_position as usize - digits.len()));
        if fractional_precision > 0 {
            output.push('.');
            output.extend(std::iter::repeat('0').take(fractional_precision));
        }
    } else {
        let integer_digits = decimal_position as usize;
        output.push_str(&digits[..integer_digits]);
        if fractional_precision > 0 {
            output.push('.');
            let available_fraction = digits.len() - integer_digits;
            let emitted_fraction = available_fraction.min(fractional_precision);
            output.push_str(&digits[integer_digits..integer_digits + emitted_fraction]);
            output.extend(std::iter::repeat('0').take(fractional_precision - emitted_fraction));
        }
    }
    output
}

fn value_error(value: &str, message: &str) -> CifReadError {
    CifReadError::new(
        CifReadErrorKind::InvalidValue,
        "cif-value",
        1,
        1,
        format!("{message}: {value}"),
    )
}

fn range_error(value: &str, message: &str) -> CifReadError {
    CifReadError::new(
        CifReadErrorKind::OutOfRange,
        "cif-value",
        1,
        1,
        format!("{message}: {value}"),
    )
}

fn is_c_space(byte: u8) -> bool {
    matches!(byte, b'\t'..=b'\r' | b' ')
}

fn char_table(byte: u8) -> u8 {
    match byte {
        b'\t' | b'\n' | b'\r' | b' ' => 2,
        b'!'
        | b'%'
        | b'&'
        | b'('..=b':'
        | b'<'..=b'Z'
        | b']'
        | b'\\'
        | b'^'
        | b'`'..=b'z'
        | b'{'
        | b'|'
        | b'}'
        | b'~' => 1,
        _ => 0,
    }
}
