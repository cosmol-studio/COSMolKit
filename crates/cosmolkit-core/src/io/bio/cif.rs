//! Internal CIF document primitives shared by the Gemmi-aligned reader and writer.

use std::io::{self, Write};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CifToken {
    pub(super) value: String,
    pub(super) line_number: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CifItem {
    pub(super) tag: String,
    pub(super) value: CifToken,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CifLoop {
    pub(super) tags: Vec<String>,
    pub(super) values: Vec<CifToken>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum CifEntry {
    Pair(usize),
    Loop(usize),
    Erased,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct CifBlock {
    pub(super) name: String,
    pub(super) items: Vec<CifItem>,
    pub(super) loops: Vec<CifLoop>,
    pub(super) entries: Vec<CifEntry>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct CifDocument {
    pub(super) blocks: Vec<CifBlock>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(super) struct CifWriteOptions {
    pub(super) prefer_pairs: bool,
    pub(super) compact: bool,
    pub(super) misuse_hash: bool,
    pub(super) align_pairs: u16,
    pub(super) align_loops: u16,
}

impl CifLoop {
    pub(super) fn add_row<I, S>(&mut self, values: I) -> Result<(), &'static str>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let values = values.into_iter().map(Into::into).collect::<Vec<_>>();
        if values.len() != self.tags.len() {
            return Err("CIF loop row width differs from its tag count");
        }
        self.values.extend(values.into_iter().map(|value| CifToken {
            value,
            line_number: 0,
        }));
        Ok(())
    }
}

impl CifBlock {
    pub(super) fn push_pair(&mut self, tag: String, value: CifToken) {
        let index = self.items.len();
        self.items.push(CifItem { tag, value });
        self.entries.push(CifEntry::Pair(index));
    }

    pub(super) fn push_loop(&mut self, loop_: CifLoop) {
        let index = self.loops.len();
        self.loops.push(loop_);
        self.entries.push(CifEntry::Loop(index));
    }

    fn entry_has_prefix(&self, entry: CifEntry, prefix: &str) -> bool {
        match entry {
            CifEntry::Pair(index) => self.items[index]
                .tag
                .get(..prefix.len())
                .is_some_and(|head| head.eq_ignore_ascii_case(prefix)),
            CifEntry::Loop(index) => self.loops[index].tags.first().is_some_and(|tag| {
                tag.get(..prefix.len())
                    .is_some_and(|head| head.eq_ignore_ascii_case(prefix))
            }),
            CifEntry::Erased => false,
        }
    }

    fn normalized_category(category: &str) -> String {
        assert!(
            category.starts_with('_'),
            "CIF category must start with '_'"
        );
        if category.ends_with('.') {
            category.to_string()
        } else {
            format!("{category}.")
        }
    }

    // BEGIN GEMMI CPP FUNCTION gemmi::cif::ItemSpan::set_pair
    // Gemmi✔️✔️: assert_tag(tag);
    // Gemmi✔️✔️: std::string lctag = gemmi::to_lower(tag);
    // Gemmi✔️✔️: auto end = items_.begin() + end_;
    // Gemmi✔️✔️: for (auto i = items_.begin() + begin_; i != end; ++i) {
    // Gemmi✔️✔️:   if (i->type == ItemType::Pair && gemmi::iequal(i->pair[0], lctag)) {
    // Gemmi✔️✔️:     i->pair[0] = tag;
    // Gemmi✔️✔️:     i->pair[1] = value;
    // Gemmi✔️✔️:     return;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   if (i->type == ItemType::Loop && i->loop.find_tag_lc(lctag) != -1) {
    // Gemmi✔️✔️:     i->set_value(Item(tag, value));
    // Gemmi✔️✔️:     return;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: items_.emplace(end, tag, value);
    // END GEMMI CPP FUNCTION
    pub(super) fn set_pair(&mut self, tag: &str, value: String) {
        self.set_pair_in_span(None, tag, value);
    }

    pub(super) fn set_pair_in_category(&mut self, category: &str, tag: &str, value: String) {
        let category = Self::normalized_category(category);
        self.set_pair_in_span(Some(&category), tag, value);
    }

    fn set_pair_in_span(&mut self, category: Option<&str>, tag: &str, value: String) {
        assert!(tag.starts_with('_'), "CIF tag must start with '_'");
        let span = category
            .map(|prefix| {
                let begin = self
                    .entries
                    .iter()
                    .position(|entry| self.entry_has_prefix(*entry, prefix))
                    .unwrap_or(self.entries.len());
                let end = self
                    .entries
                    .iter()
                    .rposition(|entry| self.entry_has_prefix(*entry, prefix))
                    .map_or(begin, |index| index + 1);
                begin..end
            })
            .unwrap_or(0..self.entries.len());

        for entry_index in span.clone() {
            match self.entries[entry_index] {
                CifEntry::Pair(item_index)
                    if self.items[item_index].tag.eq_ignore_ascii_case(tag) =>
                {
                    self.items[item_index].tag = tag.to_string();
                    self.items[item_index].value = CifToken {
                        value,
                        line_number: 0,
                    };
                    return;
                }
                CifEntry::Loop(loop_index)
                    if self.loops[loop_index]
                        .tags
                        .iter()
                        .any(|candidate| candidate.eq_ignore_ascii_case(tag)) =>
                {
                    let item_index = self.items.len();
                    self.items.push(CifItem {
                        tag: tag.to_string(),
                        value: CifToken {
                            value,
                            line_number: 0,
                        },
                    });
                    self.entries[entry_index] = CifEntry::Pair(item_index);
                    return;
                }
                _ => {}
            }
        }

        let item_index = self.items.len();
        self.items.push(CifItem {
            tag: tag.to_string(),
            value: CifToken {
                value,
                line_number: 0,
            },
        });
        self.entries.insert(span.end, CifEntry::Pair(item_index));
    }

    // BEGIN GEMMI CPP FUNCTION gemmi::cif::Block::setup_loop_item
    // Gemmi✔️✔️: if (tab.loop_item) {
    // Gemmi✔️✔️:   item = tab.loop_item;
    // Gemmi✔️✔️:   item->loop.clear();
    // Gemmi✔️✔️: } else if (tab.ok()) {
    // Gemmi✔️✔️:   item = &tab.bloc.items.at(tab.positions[0]);
    // Gemmi✔️✔️:   tab.erase();
    // Gemmi✔️✔️:   item->set_value(Item(LoopArg{}));
    // Gemmi✔️✔️: } else {
    // Gemmi✔️✔️:   items.emplace_back(LoopArg{});
    // Gemmi✔️✔️:   item = &items.back();
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: for (std::string& tag : tags) {
    // Gemmi✔️✔️:   tag.insert(0, prefix);
    // Gemmi✔️✔️:   assert_tag(tag);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: item->loop.tags = std::move(tags);
    // END GEMMI CPP FUNCTION
    pub(super) fn init_mmcif_loop(&mut self, category: &str, suffixes: &[&str]) -> &mut CifLoop {
        let category = Self::normalized_category(category);
        let tags = suffixes
            .iter()
            .map(|suffix| format!("{category}{suffix}"))
            .collect::<Vec<_>>();
        let matching = self
            .entries
            .iter()
            .enumerate()
            .filter_map(|(index, entry)| self.entry_has_prefix(*entry, &category).then_some(index))
            .collect::<Vec<_>>();

        if let Some(&entry_index) = matching.first() {
            if let CifEntry::Loop(loop_index) = self.entries[entry_index] {
                self.loops[loop_index].tags = tags;
                self.loops[loop_index].values.clear();
                for &extra in matching.iter().skip(1) {
                    self.entries[extra] = CifEntry::Erased;
                }
                return &mut self.loops[loop_index];
            }

            let loop_index = self.loops.len();
            self.loops.push(CifLoop {
                tags,
                values: Vec::new(),
            });
            self.entries[entry_index] = CifEntry::Loop(loop_index);
            for &extra in matching.iter().skip(1) {
                self.entries[extra] = CifEntry::Erased;
            }
            return &mut self.loops[loop_index];
        }

        let loop_index = self.loops.len();
        self.loops.push(CifLoop {
            tags,
            values: Vec::new(),
        });
        self.entries.push(CifEntry::Loop(loop_index));
        &mut self.loops[loop_index]
    }

    pub(super) fn erase_mmcif_category(&mut self, category: &str) {
        let category = Self::normalized_category(category);
        let matches = self
            .entries
            .iter()
            .map(|entry| self.entry_has_prefix(*entry, &category))
            .collect::<Vec<_>>();
        for (entry, matches) in self.entries.iter_mut().zip(matches) {
            if matches {
                *entry = CifEntry::Erased;
            }
        }
    }
}

fn is_unquoted_char(byte: u8) -> bool {
    matches!(
        byte,
        b'!'
            | b'%'
            | b'&'
            | b'('
            | b')'
            | b'*'
            | b'+'
            | b','
            | b'-'
            | b'.'
            | b'/'
            | b'0'..=b'9'
            | b':'
            | b'<'
            | b'='
            | b'>'
            | b'?'
            | b'@'
            | b'A'..=b'Z'
            | b'\\'
            | b']'
            | b'^'
            | b'`'
            | b'a'..=b'z'
            | b'|'
            | b'~'
    )
}

pub(super) fn is_cif_null(value: &str) -> bool {
    matches!(value.as_bytes(), [b'?'] | [b'.'])
}

pub(super) fn is_cif_text_field(value: &str) -> bool {
    let bytes = value.as_bytes();
    bytes.len() > 2 && bytes[0] == b';' && matches!(bytes.get(bytes.len() - 2), Some(b'\n' | b'\r'))
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::quote
// Gemmi✔️✔️: if (std::all_of(v.begin(), v.end(), [](char c) { return char_table(c) == 1; })
// Gemmi✔️✔️:     && !v.empty() && !is_null(v))
// Gemmi✔️✔️:   return v;
// Gemmi✔️✔️: char q = ';';
// Gemmi✔️✔️: if (std::memchr(v.c_str(), '\n', v.size()) == nullptr) {
// Gemmi✔️✔️:   if (std::memchr(v.c_str(), '\'', v.size()) == nullptr)
// Gemmi✔️✔️:     q = '\'';
// Gemmi✔️✔️:   else if (std::memchr(v.c_str(), '"', v.size()) == nullptr)
// Gemmi✔️✔️:     q = '"';
// Gemmi✔️✔️: }
// Gemmi✔️✔️: v.insert(v.begin(), q);
// Gemmi✔️✔️: if (q == ';')
// Gemmi✔️✔️:   v += '\n';
// Gemmi✔️✔️: v += q;
// Gemmi✔️✔️: return v;
// END GEMMI CPP FUNCTION
pub(super) fn quote_cif_value(value: &str) -> String {
    if !value.is_empty()
        && !is_cif_null(value)
        && value.as_bytes().iter().copied().all(is_unquoted_char)
    {
        return value.to_string();
    }
    if !value.contains('\n') && !value.contains('\'') {
        return format!("'{value}'");
    }
    if !value.contains('\n') && !value.contains('"') {
        return format!("\"{value}\"");
    }
    format!(";{value}\n;")
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::write_text_field
// Gemmi✔️✔️: for (size_t pos = 0, end = 0; end != std::string::npos; pos = end + 1) {
// Gemmi✔️✔️:   end = value.find("\r\n", pos);
// Gemmi✔️✔️:   size_t len = (end == std::string::npos ? value.size() : end) - pos;
// Gemmi✔️✔️:   os.write(value.c_str() + pos, len);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn write_text_field<W: Write>(writer: &mut W, value: &str) -> io::Result<()> {
    let mut position = 0;
    while let Some(relative_end) = value[position..].find("\r\n") {
        let end = position + relative_end;
        writer.write_all(value[position..end].as_bytes())?;
        position = end + 1;
    }
    writer.write_all(value[position..].as_bytes())
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::write_out_pair
// Gemmi✔️✔️: os << name;
// Gemmi✔️✔️: if (is_text_field(value)) {
// Gemmi✔️✔️:   os.put('\n');
// Gemmi✔️✔️:   write_text_field(os, value);
// Gemmi✔️✔️: } else {
// Gemmi✔️✔️:   if (name.size() + value.size() > 120) {
// Gemmi✔️✔️:     os.put('\n');
// Gemmi✔️✔️:   } else {
// Gemmi✔️✔️:     os.put(' ');
// Gemmi✔️✔️:     if (name.size() < options.align_pairs)
// Gemmi✔️✔️:       os.pad(options.align_pairs - name.size());
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   os << value;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: os.put('\n');
// END GEMMI CPP FUNCTION
fn write_pair<W: Write>(
    writer: &mut W,
    item: &CifItem,
    options: CifWriteOptions,
) -> io::Result<()> {
    writer.write_all(item.tag.as_bytes())?;
    if is_cif_text_field(&item.value.value) {
        writer.write_all(b"\n")?;
        write_text_field(writer, &item.value.value)?;
    } else {
        if item.tag.len() + item.value.value.len() > 120 {
            writer.write_all(b"\n")?;
        } else {
            writer.write_all(b" ")?;
            if item.tag.len() < options.align_pairs as usize {
                write_spaces(writer, options.align_pairs as usize - item.tag.len())?;
            }
        }
        writer.write_all(item.value.value.as_bytes())?;
    }
    writer.write_all(b"\n")
}

fn write_spaces<W: Write>(writer: &mut W, count: usize) -> io::Result<()> {
    const SPACES: &[u8; 64] = b"                                                                ";
    let mut remaining = count;
    while remaining != 0 {
        let length = remaining.min(SPACES.len());
        writer.write_all(&SPACES[..length])?;
        remaining -= length;
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::write_out_loop
// Gemmi✔️✔️: if (loop.values.empty())
// Gemmi✔️✔️:   return;
// Gemmi✔️✔️: if (options.prefer_pairs && loop.length() == 1) {
// Gemmi✔️✔️:   for (size_t i = 0; i != loop.tags.size(); ++i)
// Gemmi✔️✔️:     write_out_pair(os, loop.tags[i], loop.values[i], options);
// Gemmi✔️✔️:   return;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: os.write("loop_", 5);
// Gemmi✔️✔️: for (const std::string& tag : loop.tags) {
// Gemmi✔️✔️:   os.put('\n');
// Gemmi✔️✔️:   os << tag;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: size_t ncol = loop.tags.size();
// Gemmi✔️✔️: std::vector<size_t> col_width(ncol, 0);
// Gemmi✔️✔️: if (options.align_loops > 0) {
// Gemmi✔️✔️:   size_t col = 0;
// Gemmi✔️✔️:   for (const std::string& val : loop.values) {
// Gemmi✔️✔️:     if (!is_text_field(val))
// Gemmi✔️✔️:       col_width[col] = std::max(col_width[col], val.size());
// Gemmi✔️✔️:     if (++col == ncol)
// Gemmi✔️✔️:       col = 0;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   for (size_t& w : col_width)
// Gemmi✔️✔️:     w = std::min(w, (size_t)options.align_loops);
// Gemmi✔️✔️: }
// Gemmi✔️✔️: size_t col = 0;
// Gemmi✔️✔️: bool need_new_line = true;
// Gemmi✔️✔️: for (const std::string& val : loop.values) {
// Gemmi✔️✔️:   bool text_field = is_text_field(val);
// Gemmi✔️✔️:   os.put(need_new_line || text_field ? '\n' : ' ');
// Gemmi✔️✔️:   need_new_line = text_field;
// Gemmi✔️✔️:   if (text_field)
// Gemmi✔️✔️:     write_text_field(os, val);
// Gemmi✔️✔️:   else
// Gemmi✔️✔️:     os << val;
// Gemmi✔️✔️:   if (col != ncol - 1) {
// Gemmi✔️✔️:     if (val.size() < col_width[col])
// Gemmi✔️✔️:       os.pad(col_width[col] - val.size());
// Gemmi✔️✔️:     ++col;
// Gemmi✔️✔️:   } else {
// Gemmi✔️✔️:     col = 0;
// Gemmi✔️✔️:     need_new_line = true;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️: }
// Gemmi✔️✔️: os.put('\n');
// END GEMMI CPP FUNCTION
fn write_loop<W: Write>(
    writer: &mut W,
    loop_: &CifLoop,
    options: CifWriteOptions,
) -> io::Result<()> {
    if loop_.values.is_empty() {
        return Ok(());
    }
    let column_count = loop_.tags.len();
    debug_assert_ne!(column_count, 0);
    debug_assert_eq!(loop_.values.len() % column_count, 0);
    if options.prefer_pairs && loop_.values.len() == column_count {
        for (tag, value) in loop_.tags.iter().zip(&loop_.values) {
            write_pair(
                writer,
                &CifItem {
                    tag: tag.clone(),
                    value: value.clone(),
                },
                options,
            )?;
        }
        return Ok(());
    }

    writer.write_all(b"loop_")?;
    for tag in &loop_.tags {
        writer.write_all(b"\n")?;
        writer.write_all(tag.as_bytes())?;
    }

    let mut column_widths = vec![0; column_count];
    if options.align_loops != 0 {
        for (column, value) in loop_.values.iter().enumerate() {
            if !is_cif_text_field(&value.value) {
                let width = &mut column_widths[column % column_count];
                *width = (*width).max(value.value.len());
            }
        }
        for width in &mut column_widths {
            *width = (*width).min(options.align_loops as usize);
        }
    }

    let mut column = 0;
    let mut need_new_line = true;
    for value in &loop_.values {
        let text_field = is_cif_text_field(&value.value);
        writer.write_all(if need_new_line || text_field {
            b"\n"
        } else {
            b" "
        })?;
        need_new_line = text_field;
        if text_field {
            write_text_field(writer, &value.value)?;
        } else {
            writer.write_all(value.value.as_bytes())?;
        }
        if column != column_count - 1 {
            if value.value.len() < column_widths[column] {
                write_spaces(writer, column_widths[column] - value.value.len())?;
            }
            column += 1;
        } else {
            column = 0;
            need_new_line = true;
        }
    }
    writer.write_all(b"\n")
}

fn entry_category_tag<'a>(block: &'a CifBlock, entry: CifEntry) -> Option<&'a str> {
    match entry {
        CifEntry::Pair(index) => Some(block.items[index].tag.as_str()),
        CifEntry::Loop(index) => block.loops[index].tags.first().map(String::as_str),
        CifEntry::Erased => None,
    }
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::should_be_separated_
// Gemmi✔️✔️: if (a.type == ItemType::Comment || b.type == ItemType::Comment)
// Gemmi✔️✔️:   return false;
// Gemmi✔️✔️: if (a.type != ItemType::Pair || b.type != ItemType::Pair)
// Gemmi✔️✔️:   return true;
// Gemmi✔️✔️: auto adot = a.pair[0].find('.');
// Gemmi✔️✔️: if (adot == std::string::npos)
// Gemmi✔️✔️:   return false;
// Gemmi✔️✔️: auto bdot = b.pair[0].find('.');
// Gemmi✔️✔️: return adot != bdot || a.pair[0].compare(0, adot, b.pair[0], 0, adot) != 0;
// END GEMMI CPP FUNCTION
fn should_be_separated(block: &CifBlock, first: CifEntry, second: CifEntry) -> bool {
    let (CifEntry::Pair(first_index), CifEntry::Pair(second_index)) = (first, second) else {
        return true;
    };
    let first_tag = &block.items[first_index].tag;
    let Some(first_dot) = first_tag.find('.') else {
        return false;
    };
    let second_tag = &block.items[second_index].tag;
    let second_dot = second_tag.find('.');
    second_dot != Some(first_dot) || first_tag.get(..first_dot) != second_tag.get(..first_dot)
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::write_cif_block_to_stream
// Gemmi✔️✔️: os.write("data_", 5);
// Gemmi✔️✔️: os << block.name;
// Gemmi✔️✔️: os.put('\n');
// Gemmi✔️✔️: if (options.misuse_hash)
// Gemmi✔️✔️:   os.write("#\n", 2);
// Gemmi✔️✔️: const Item* prev = nullptr;
// Gemmi✔️✔️: for (const Item& item : block.items) {
// Gemmi✔️✔️:   if (item.type == ItemType::Erased)
// Gemmi✔️✔️:     continue;
// Gemmi✔️✔️:   if (prev && !options.compact && should_be_separated_(*prev, item)) {
// Gemmi✔️✔️:     if (options.misuse_hash)
// Gemmi✔️✔️:       os.put('#');
// Gemmi✔️✔️:     os.put('\n');
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   write_out_item(os, item, options);
// Gemmi✔️✔️:   prev = &item;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: if (options.misuse_hash)
// Gemmi✔️✔️:   os.write("#\n", 2);
// END GEMMI CPP FUNCTION
pub(super) fn write_cif_block<W: Write>(
    writer: &mut W,
    block: &CifBlock,
    options: CifWriteOptions,
) -> io::Result<()> {
    writer.write_all(b"data_")?;
    writer.write_all(block.name.as_bytes())?;
    writer.write_all(b"\n")?;
    if options.misuse_hash {
        writer.write_all(b"#\n")?;
    }
    let mut previous = None;
    for &entry in &block.entries {
        if entry == CifEntry::Erased {
            continue;
        }
        if let Some(previous) = previous
            && !options.compact
            && should_be_separated(block, previous, entry)
        {
            if options.misuse_hash {
                writer.write_all(b"#")?;
            }
            writer.write_all(b"\n")?;
        }
        match entry {
            CifEntry::Pair(index) => write_pair(writer, &block.items[index], options)?,
            CifEntry::Loop(index) => write_loop(writer, &block.loops[index], options)?,
            CifEntry::Erased => unreachable!(),
        }
        previous = Some(entry);
    }
    if options.misuse_hash {
        writer.write_all(b"#\n")?;
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::write_cif_to_stream
// Gemmi✔️✔️: bool first = true;
// Gemmi✔️✔️: for (const Block& block : doc.blocks) {
// Gemmi✔️✔️:   if (!first)
// Gemmi✔️✔️:     os.put('\n');
// Gemmi✔️✔️:   write_cif_block_to_stream(os, block, options);
// Gemmi✔️✔️:   first = false;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn write_cif_document<W: Write>(
    writer: W,
    document: &CifDocument,
    options: CifWriteOptions,
) -> io::Result<()> {
    let mut writer = io::BufWriter::with_capacity(4096, writer);
    for (index, block) in document.blocks.iter().enumerate() {
        if index != 0 {
            writer.write_all(b"\n")?;
        }
        write_cif_block(&mut writer, block, options)?;
    }
    writer.flush()
}

pub(super) fn cif_document_to_string(
    document: &CifDocument,
    options: CifWriteOptions,
) -> io::Result<String> {
    let mut bytes = Vec::new();
    write_cif_document(&mut bytes, document, options)?;
    String::from_utf8(bytes).map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn token(value: &str) -> CifToken {
        CifToken {
            value: value.to_string(),
            line_number: 1,
        }
    }

    #[test]
    fn cif_writer_quote_matches_gemmi_delimiter_selection() {
        assert_eq!(quote_cif_value("ALA"), "ALA");
        assert_eq!(quote_cif_value("."), "'.'");
        assert_eq!(quote_cif_value(""), "''");
        assert_eq!(quote_cif_value("two words"), "'two words'");
        assert_eq!(quote_cif_value("O'Brien"), "\"O'Brien\"");
        assert_eq!(
            quote_cif_value("a 'single' and \"double\" quote"),
            ";a 'single' and \"double\" quote\n;"
        );
        assert_eq!(quote_cif_value("line 1\nline 2"), ";line 1\nline 2\n;");
    }

    #[test]
    fn cif_writer_set_pair_updates_existing_pair_without_moving_it() {
        let mut block = CifBlock::default();
        block.push_pair("_entry.id".to_string(), token("old"));
        block.push_pair("_cell.length_a".to_string(), token("1"));

        block.set_pair("_ENTRY.ID", "new".to_string());

        assert_eq!(block.entries, [CifEntry::Pair(0), CifEntry::Pair(1)]);
        assert_eq!(block.items[0].tag, "_ENTRY.ID");
        assert_eq!(block.items[0].value.value, "new");
    }

    #[test]
    fn cif_writer_category_span_keeps_new_pairs_adjacent() {
        let mut block = CifBlock::default();
        block.push_pair("_entry.id".to_string(), token("x"));
        block.push_pair("_cell.length_a".to_string(), token("10"));
        block.push_pair("_struct.title".to_string(), token("title"));

        block.set_pair_in_category("_cell.", "_cell.length_b", "11".to_string());

        assert_eq!(
            block.entries,
            [
                CifEntry::Pair(0),
                CifEntry::Pair(1),
                CifEntry::Pair(3),
                CifEntry::Pair(2),
            ]
        );
    }

    #[test]
    fn cif_writer_init_loop_replaces_category_pairs_at_their_first_position() {
        let mut block = CifBlock::default();
        block.push_pair("_entry.id".to_string(), token("x"));
        block.push_pair("_atom_site.id".to_string(), token("1"));
        block.push_pair("_atom_site.type_symbol".to_string(), token("C"));
        block.push_pair("_cell.length_a".to_string(), token("10"));

        let loop_ = block.init_mmcif_loop("_atom_site", &["id", "type_symbol"]);
        loop_.add_row(["2", "N"]).unwrap();

        assert_eq!(
            block.entries,
            [
                CifEntry::Pair(0),
                CifEntry::Loop(0),
                CifEntry::Erased,
                CifEntry::Pair(3),
            ]
        );
        assert_eq!(
            block.loops[0].tags,
            ["_atom_site.id", "_atom_site.type_symbol"]
        );
        assert_eq!(block.loops[0].values[0].value, "2");
    }

    #[test]
    fn cif_writer_set_pair_replaces_a_loop_containing_the_same_tag() {
        let mut block = CifBlock::default();
        block
            .init_mmcif_loop("_entity", &["id", "type"])
            .add_row(["1", "polymer"])
            .unwrap();

        block.set_pair("_entity.id", "2".to_string());

        assert_eq!(block.entries, [CifEntry::Pair(0)]);
        assert_eq!(block.items[0].tag, "_entity.id");
        assert_eq!(block.items[0].value.value, "2");
    }

    #[test]
    fn cif_writer_loop_rejects_rows_with_the_wrong_width() {
        let mut loop_ = CifLoop {
            tags: vec!["_atom.id".to_string(), "_atom.type".to_string()],
            values: Vec::new(),
        };
        assert_eq!(
            loop_.add_row(["1"]),
            Err("CIF loop row width differs from its tag count")
        );
        assert!(loop_.values.is_empty());
    }

    #[test]
    fn cif_writer_erase_category_preserves_unrelated_order_entries() {
        let mut block = CifBlock::default();
        block.push_pair("_entry.id".to_string(), token("x"));
        block
            .init_mmcif_loop("_atom_site", &["id"])
            .add_row(["1"])
            .unwrap();
        block.push_pair("_cell.length_a".to_string(), token("10"));

        block.erase_mmcif_category("_atom_site.");

        assert_eq!(
            block.entries,
            [CifEntry::Pair(0), CifEntry::Erased, CifEntry::Pair(1)]
        );
    }

    #[test]
    fn cif_serializer_default_separates_categories() {
        let mut block = CifBlock {
            name: "model".to_string(),
            ..CifBlock::default()
        };
        block.push_pair("_entry.id".to_string(), token("X"));
        block.push_pair("_cell.length_a".to_string(), token("10"));
        let document = CifDocument {
            blocks: vec![block],
        };

        assert_eq!(
            cif_document_to_string(&document, CifWriteOptions::default()).unwrap(),
            "data_model\n_entry.id X\n\n_cell.length_a 10\n"
        );
    }

    #[test]
    fn cif_serializer_compact_removes_category_blank_lines() {
        let mut block = CifBlock {
            name: "model".to_string(),
            ..CifBlock::default()
        };
        block.push_pair("_entry.id".to_string(), token("X"));
        block.push_pair("_cell.length_a".to_string(), token("10"));
        let document = CifDocument {
            blocks: vec![block],
        };

        assert_eq!(
            cif_document_to_string(
                &document,
                CifWriteOptions {
                    compact: true,
                    ..CifWriteOptions::default()
                }
            )
            .unwrap(),
            "data_model\n_entry.id X\n_cell.length_a 10\n"
        );
    }

    #[test]
    fn cif_serializer_misuse_hash_matches_gemmi_pdbx_style() {
        let mut block = CifBlock {
            name: "model".to_string(),
            ..CifBlock::default()
        };
        block.push_pair("_entry.id".to_string(), token("X"));
        block.push_pair("_cell.length_a".to_string(), token("10"));
        let document = CifDocument {
            blocks: vec![block],
        };

        assert_eq!(
            cif_document_to_string(
                &document,
                CifWriteOptions {
                    misuse_hash: true,
                    ..CifWriteOptions::default()
                }
            )
            .unwrap(),
            "data_model\n#\n_entry.id X\n#\n_cell.length_a 10\n#\n"
        );
    }

    #[test]
    fn cif_serializer_prefer_pairs_writes_single_row_loop_as_pairs() {
        let mut block = CifBlock {
            name: "model".to_string(),
            ..CifBlock::default()
        };
        block
            .init_mmcif_loop("_atom", &["id", "type"])
            .add_row(["1", "C"])
            .unwrap();
        let document = CifDocument {
            blocks: vec![block],
        };

        assert_eq!(
            cif_document_to_string(
                &document,
                CifWriteOptions {
                    prefer_pairs: true,
                    ..CifWriteOptions::default()
                }
            )
            .unwrap(),
            "data_model\n_atom.id 1\n_atom.type C\n"
        );
    }

    #[test]
    fn cif_serializer_aligns_pairs_and_loop_columns() {
        let mut block = CifBlock {
            name: "x".to_string(),
            ..CifBlock::default()
        };
        block.push_pair("_x.a".to_string(), token("v"));
        let loop_ = block.init_mmcif_loop("_atom", &["id", "type"]);
        loop_.add_row(["1", "C"]).unwrap();
        loop_.add_row(["22", "N"]).unwrap();
        let document = CifDocument {
            blocks: vec![block],
        };

        assert_eq!(
            cif_document_to_string(
                &document,
                CifWriteOptions {
                    align_pairs: 10,
                    align_loops: 30,
                    ..CifWriteOptions::default()
                }
            )
            .unwrap(),
            "data_x\n_x.a       v\n\nloop_\n_atom.id\n_atom.type\n1  C\n22 N\n"
        );
    }

    #[test]
    fn cif_serializer_text_fields_normalize_crlf_without_changing_delimiters() {
        let mut block = CifBlock {
            name: "x".to_string(),
            ..CifBlock::default()
        };
        block.push_pair("_struct.title".to_string(), token(";line 1\r\nline 2\n;"));
        let document = CifDocument {
            blocks: vec![block],
        };

        assert_eq!(
            cif_document_to_string(&document, CifWriteOptions::default()).unwrap(),
            "data_x\n_struct.title\n;line 1\nline 2\n;\n"
        );
    }

    #[test]
    fn cif_serializer_writes_long_pair_values_on_the_next_line() {
        let mut block = CifBlock {
            name: "x".to_string(),
            ..CifBlock::default()
        };
        let value = "A".repeat(121);
        block.push_pair("_entry.id".to_string(), token(&value));
        let document = CifDocument {
            blocks: vec![block],
        };

        assert_eq!(
            cif_document_to_string(&document, CifWriteOptions::default()).unwrap(),
            format!("data_x\n_entry.id\n{value}\n")
        );
    }

    #[test]
    fn cif_serializer_separates_document_blocks_with_one_blank_line() {
        let document = CifDocument {
            blocks: vec![
                CifBlock {
                    name: "a".to_string(),
                    ..CifBlock::default()
                },
                CifBlock {
                    name: "b".to_string(),
                    ..CifBlock::default()
                },
            ],
        };

        assert_eq!(
            cif_document_to_string(&document, CifWriteOptions::default()).unwrap(),
            "data_a\n\ndata_b\n"
        );
    }
}
