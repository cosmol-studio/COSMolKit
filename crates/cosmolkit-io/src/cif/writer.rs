//! Private CIF serializer for the writer pipeline (single `CifBlock` model).

use std::io::{self, Write};

use super::{CifBlock, CifDocument, CifItem, CifLoop};
use crate::bio_write::BioMmcifWriteParams;

/// Gemmi `cif::WriteOptions` layout controls used by this serializer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct CifWriteLayout {
    pub(super) prefer_pairs: bool,
    pub(super) compact: bool,
    pub(super) misuse_hash: bool,
    pub(super) align_pairs: u16,
    pub(super) align_loops: u16,
}

impl From<&BioMmcifWriteParams> for CifWriteLayout {
    fn from(params: &BioMmcifWriteParams) -> Self {
        // Gemmi✔️✔️: struct WriteOptions {
        // Gemmi✔️✔️:   bool prefer_pairs = false;
        // Gemmi✔️✔️:   bool compact = false;
        // Gemmi✔️✔️:   bool misuse_hash = false;
        // Gemmi✔️✔️:   std::uint16_t align_pairs = 0;
        // Gemmi✔️✔️:   std::uint16_t align_loops = 0;
        // Gemmi✔️✔️: };
        // Behavior: the five layout controls pass through one-to-one; the
        // coordinate profile bits (group_pdb/auth_all) are not layout state.
        // Complexity: constant-size copy, no allocation.
        Self {
            prefer_pairs: params.prefer_pairs,
            compact: params.compact,
            misuse_hash: params.misuse_hash,
            align_pairs: params.align_pairs,
            align_loops: params.align_loops,
        }
    }
}

impl Default for CifWriteLayout {
    fn default() -> Self {
        Self {
            prefer_pairs: false,
            compact: false,
            misuse_hash: false,
            align_pairs: 0,
            align_loops: 0,
        }
    }
}

fn write_spaces<W: Write + ?Sized>(writer: &mut W, count: usize) -> io::Result<()> {
    // Gemmi✔️✔️: void pad(size_t n) {
    // Gemmi✔️✔️:   std::memset(ptr, ' ', n);
    // Gemmi✔️✔️:   ptr += n;
    // Gemmi✔️✔️: }
    // Behavior: exactly `count` spaces reach the sink in bounded chunks.
    // Complexity: O(count) bytes, no allocation.
    const SPACES: &[u8; 64] = b"                                                                ";
    let mut remaining = count;
    while remaining != 0 {
        let length = remaining.min(SPACES.len());
        writer.write_all(&SPACES[..length])?;
        remaining -= length;
    }
    Ok(())
}

/// Gemmi `cif::write_out_pair`: emit one tag/value pair; `value` must
/// already be CIF-quoted by the caller.
pub(super) fn write_out_pair<W: Write + ?Sized>(
    writer: &mut W,
    tag: &str,
    value: &str,
    layout: &CifWriteLayout,
) -> io::Result<()> {
    // Gemmi✔️✔️: inline void write_out_pair(BufOstream& os, const std::string& name,
    // Gemmi✔️✔️:                            const std::string& value, WriteOptions options) {
    // Gemmi✔️✔️:   os << name;
    // Gemmi✔️✔️:   if (is_text_field(value)) {
    // Gemmi✔️✔️:     os.put('\n');
    // Gemmi✔️✔️:     write_text_field(os, value);
    // Gemmi✔️✔️:   } else {
    // Gemmi✔️✔️:     if (name.size() + value.size() > 120) {
    // Gemmi✔️✔️:       os.put('\n');
    // Gemmi✔️✔️:     } else {
    // Gemmi✔️✔️:       os.put(' ');
    // Gemmi✔️✔️:       if (name.size() < options.align_pairs)
    // Gemmi✔️✔️:         os.pad(options.align_pairs - name.size());
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     os << value;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   os.put('\n');
    // Gemmi✔️✔️: }
    // Behavior: text fields break to their own line and reuse c04 emission;
    // pairs whose combined tag+value length exceeds 120 bytes break before
    // the value; otherwise one space plus optional align_pairs padding.
    // Complexity: O(tag+value) bytes with constant bookkeeping, one pass.
    writer.write_all(tag.as_bytes())?;
    if is_text_field(value) {
        writer.write_all(b"\n")?;
        write_text_field(writer, value)?;
    } else if tag.len() + value.len() > 120 {
        writer.write_all(b"\n")?;
        writer.write_all(value.as_bytes())?;
    } else {
        writer.write_all(b" ")?;
        if tag.len() < layout.align_pairs as usize {
            write_spaces(writer, layout.align_pairs as usize - tag.len())?;
        }
        writer.write_all(value.as_bytes())?;
    }
    writer.write_all(b"\n")
}

/// Gemmi `cif::write_out_loop`: emit one loop, honoring single-row
/// `prefer_pairs`, `align_loops` column widths and empty-loop suppression.
pub(super) fn write_out_loop<W: Write + ?Sized>(
    writer: &mut W,
    loop_: &CifLoop,
    layout: &CifWriteLayout,
) -> io::Result<()> {
    // Gemmi✔️✔️: inline void write_out_loop(BufOstream& os, const Loop& loop, WriteOptions options) {
    // Gemmi✔️✔️:   if (loop.values.empty())
    // Gemmi✔️✔️:     return;
    // Gemmi✔️✔️:   if (options.prefer_pairs && loop.length() == 1) {
    // Gemmi✔️✔️:     for (size_t i = 0; i != loop.tags.size(); ++i)
    // Gemmi✔️✔️:       write_out_pair(os, loop.tags[i], loop.values[i], options);
    // Gemmi✔️✔️:     return;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   // tags
    // Gemmi✔️✔️:   os.write("loop_", 5);
    // Gemmi✔️✔️:   for (const std::string& tag : loop.tags) {
    // Gemmi✔️✔️:     os.put('\n');
    // Gemmi✔️✔️:     os << tag;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   // values
    // Gemmi✔️✔️:   size_t ncol = loop.tags.size();
    // Gemmi✔️✔️:   std::vector<size_t> col_width(ncol, 0);
    // Gemmi✔️✔️:   if (options.align_loops > 0) {
    // Gemmi✔️✔️:     size_t col = 0;
    // Gemmi✔️✔️:     for (const std::string& val : loop.values) {
    // Gemmi✔️✔️:       if (!is_text_field(val))
    // Gemmi✔️✔️:         col_width[col] = std::max(col_width[col], val.size());
    // Gemmi✔️✔️:       if (++col == ncol)
    // Gemmi✔️✔️:         col = 0;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     for (size_t& w : col_width)
    // Gemmi✔️✔️:       w = std::min(w, (size_t)options.align_loops);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   size_t col = 0;
    // Gemmi✔️✔️:   bool need_new_line = true;
    // Gemmi✔️✔️:   for (const std::string& val : loop.values) {
    // Gemmi✔️✔️:     bool text_field = is_text_field(val);
    // Gemmi✔️✔️:     os.put(need_new_line || text_field ? '\n' : ' ');
    // Gemmi✔️✔️:     need_new_line = text_field;
    // Gemmi✔️✔️:     if (text_field)
    // Gemmi✔️✔️:       write_text_field(os, val);
    // Gemmi✔️✔️:     else
    // Gemmi✔️✔️:       os << val;
    // Gemmi✔️✔️:     if (col != ncol - 1) {
    // Gemmi✔️✔️:       if (val.size() < col_width[col])
    // Gemmi✔️✔️:         os.pad(col_width[col] - val.size());
    // Gemmi✔️✔️:       ++col;
    // Gemmi✔️✔️:     } else {
    // Gemmi✔️✔️:       col = 0;
    // Gemmi✔️✔️:       need_new_line = true;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   os.put('\n');
    // Gemmi✔️✔️: }
    // Behavior: empty loops emit nothing; a single-row loop under
    // prefer_pairs reuses c05 pair emission; column widths consider only
    // non-text-field values, are capped at align_loops, and padding applies
    // only between columns of one row. Rows start on new lines; a text-field
    // cell forces the following separator back to a newline.
    // Complexity: two linear passes over values plus one O(ncol) width
    // vector; no per-cell allocation.
    let values: Vec<&str> = loop_.values.iter().map(|value| value.raw()).collect();
    if values.is_empty() {
        return Ok(());
    }
    let column_count = loop_.tags.len();
    debug_assert_ne!(column_count, 0);
    debug_assert_eq!(values.len() % column_count, 0);
    if layout.prefer_pairs && values.len() == column_count {
        for (tag, value) in loop_.tags.iter().zip(&values) {
            write_out_pair(writer, tag, value, layout)?;
        }
        return Ok(());
    }
    writer.write_all(b"loop_")?;
    for tag in &loop_.tags {
        writer.write_all(b"\n")?;
        writer.write_all(tag.as_bytes())?;
    }

    let mut column_widths = vec![0usize; column_count];
    if layout.align_loops != 0 {
        for (index, value) in values.iter().enumerate() {
            if !is_text_field(value) {
                let width = &mut column_widths[index % column_count];
                *width = (*width).max(value.len());
            }
        }
        for width in &mut column_widths {
            *width = (*width).min(layout.align_loops as usize);
        }
    }

    let mut column = 0;
    let mut need_new_line = true;
    for value in &values {
        let text_field = is_text_field(value);
        writer.write_all(if need_new_line || text_field {
            b"\n"
        } else {
            b" "
        })?;
        need_new_line = text_field;
        if text_field {
            write_text_field(writer, value)?;
        } else {
            writer.write_all(value.as_bytes())?;
        }
        if column != column_count - 1 {
            if value.len() < column_widths[column] {
                write_spaces(writer, column_widths[column] - value.len())?;
            }
            column += 1;
        } else {
            column = 0;
            need_new_line = true;
        }
    }
    writer.write_all(b"\n")
}

/// Gemmi `cif::write_out_item`: dispatch one stored item to its emitter.
pub(super) fn write_out_item<W: Write + ?Sized>(
    writer: &mut W,
    item: &CifItem,
    layout: &CifWriteLayout,
) -> io::Result<()> {
    // Gemmi✔️✔️: inline void write_out_item(BufOstream& os, const Item& item, WriteOptions options) {
    // Gemmi✔️✔️:   switch (item.type) {
    // Gemmi✔️✔️:     case ItemType::Pair:
    // Gemmi✔️✔️:       write_out_pair(os, item.pair[0], item.pair[1], options);
    // Gemmi✔️✔️:       break;
    // Gemmi✔️✔️:     case ItemType::Loop:
    // Gemmi✔️✔️:       write_out_loop(os, item.loop, options);
    // Gemmi✔️✔️:       break;
    // Gemmi✔️✔️:     case ItemType::Frame:
    // Gemmi✔️✔️:       os.write("save_", 5);
    // Gemmi✔️✔️:       os << item.frame.name;
    // Gemmi✔️✔️:       os.put('\n');
    // Gemmi✔️✔️:       for (const Item& inner_item : item.frame.items)
    // Gemmi✔️✔️:         write_out_item(os, inner_item, options);
    // Gemmi✔️✔️:       os.write("save_\n", 6);
    // Gemmi✔️✔️:       break;
    // Gemmi✔️✔️:     case ItemType::Comment:
    // Gemmi✔️✔️:       os << item.pair[1];
    // Gemmi✔️✔️:       os.put('\n');
    // Gemmi✔️✔️:       break;
    // Gemmi✔️✔️:     case ItemType::Erased:
    // Gemmi✔️✔️:       break;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // Behavior: pairs forward to c05 with their stored raw value (a
    // syntax-only `None` value has no Gemmi counterpart and maps to the
    // empty stored string, which c05 emits exactly as Gemmi emits an empty
    // value; the coordinate pipeline never constructs `None`); loops forward
    // to c06; frames recurse between save_ markers. This model intentionally
    // stores no comment or erased items, so those source branches have no
    // corresponding input state.
    // Complexity: one dispatch plus the delegated emitter cost; recursion
    // depth equals frame nesting depth.
    match item {
        CifItem::Pair(pair) => {
            let raw = pair.value.as_ref().map_or("", |value| value.raw());
            write_out_pair(writer, pair.tag(), raw, layout)
        }
        CifItem::Loop(loop_) => write_out_loop(writer, loop_, layout),
        CifItem::Frame(frame) => {
            writer.write_all(b"save_")?;
            writer.write_all(frame.name().as_bytes())?;
            writer.write_all(b"\n")?;
            for inner in frame.items() {
                write_out_item(writer, inner, layout)?;
            }
            writer.write_all(b"save_\n")
        }
    }
}

/// Gemmi `cif::should_be_separated_`: two pairs stay unseparated only when
/// they carry the same case-sensitive category prefix before the first dot.
pub(super) fn should_be_separated(first: &CifItem, second: &CifItem) -> bool {
    // Gemmi✔️✔️: inline bool should_be_separated_(const Item& a, const Item& b) {
    // Gemmi✔️✔️:   if (a.type == ItemType::Comment || b.type == ItemType::Comment)
    // Gemmi✔️✔️:     return false;
    // Gemmi✔️✔️:   if (a.type != ItemType::Pair || b.type != ItemType::Pair)
    // Gemmi✔️✔️:     return true;
    // Gemmi✔️✔️:   // check if we have mmcif-like tags from different categories
    // Gemmi✔️✔️:   auto adot = a.pair[0].find('.');
    // Gemmi✔️✔️:   if (adot == std::string::npos)
    // Gemmi✔️✔️:     return false;
    // Gemmi✔️✔️:   auto bdot = b.pair[0].find('.');
    // Gemmi✔️✔️:   return adot != bdot || a.pair[0].compare(0, adot, b.pair[0], 0, adot) != 0;
    // Gemmi✔️✔️: }
    // Behavior: comments do not exist in this model, so that early return
    // has no corresponding input; any non-pair item separates from its
    // neighbor; dotless tags never separate; otherwise the byte-exact
    // prefixes (up to the first dot) decide.
    // Complexity: O(prefix length) byte comparison, no allocation.
    let (CifItem::Pair(first), CifItem::Pair(second)) = (first, second) else {
        return true;
    };
    let first_tag = first.tag();
    let Some(first_dot) = first_tag.find('.') else {
        return false;
    };
    let second_tag = second.tag();
    let second_dot = second_tag.find('.');
    second_dot != Some(first_dot) || first_tag.get(..first_dot) != second_tag.get(..first_dot)
}

/// Gemmi `cif::write_cif_block_to_stream`: emit one complete block with
/// its header, optional hash fences and category blank-line separation.
pub(super) fn write_cif_block_to_stream<W: Write + ?Sized>(
    writer: &mut W,
    block: &CifBlock,
    layout: &CifWriteLayout,
) -> io::Result<()> {
    // Gemmi✔️✔️: inline void write_cif_block_to_stream(std::ostream& os_, const Block& block,
    // Gemmi✔️✔️:                                       WriteOptions options=WriteOptions()) {
    // Gemmi✔️✔️:   BufOstream os(os_);
    // Gemmi✔️✔️:   os.write("data_", 5);
    // Gemmi✔️✔️:   os << block.name;
    // Gemmi✔️✔️:   os.put('\n');
    // Gemmi✔️✔️:   if (options.misuse_hash)
    // Gemmi✔️✔️:     os.write("#\n", 2);
    // Gemmi✔️✔️:   const Item* prev = nullptr;
    // Gemmi✔️✔️:   for (const Item& item : block.items) {
    // Gemmi✔️✔️:     if (item.type == ItemType::Erased)
    // Gemmi✔️✔️:       continue;
    // Gemmi✔️✔️:     if (prev && !options.compact && should_be_separated_(*prev, item)) {
    // Gemmi✔️✔️:       if (options.misuse_hash)
    // Gemmi✔️✔️:         os.put('#');
    // Gemmi✔️✔️:       os.put('\n');
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     write_out_item(os, item, options);
    // Gemmi✔️✔️:     prev = &item;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   if (options.misuse_hash)
    // Gemmi✔️✔️:     os.write("#\n", 2);
    // Gemmi✔️✔️: }
    // Behavior: data_ header, optional leading/trailing "#\n" fences, and
    // one blank line ("#\n" when misuse_hash) between items that c07 says
    // belong to different categories, suppressed entirely by compact.
    // prefer_pairs reaches items through c06. This model has no erased-item
    // state (category erasure compacts the item vector), so the Erased skip
    // has no corresponding input.
    // Complexity: one pass over items plus delegated emitter costs.
    writer.write_all(b"data_")?;
    writer.write_all(block.name().as_bytes())?;
    writer.write_all(b"\n")?;
    if layout.misuse_hash {
        writer.write_all(b"#\n")?;
    }
    let mut previous: Option<&CifItem> = None;
    for item in block.items() {
        if let Some(prev) = previous
            && !layout.compact
            && should_be_separated(prev, item)
        {
            if layout.misuse_hash {
                writer.write_all(b"#")?;
            }
            writer.write_all(b"\n")?;
        }
        write_out_item(writer, item, layout)?;
        previous = Some(item);
    }
    if layout.misuse_hash {
        writer.write_all(b"#\n")?;
    }
    Ok(())
}

/// Gemmi `cif::write_cif_to_stream`: serialize every block in document
/// order, one blank line between blocks.
pub(super) fn write_cif_document<W: Write + ?Sized>(
    writer: &mut W,
    document: &CifDocument,
    layout: &CifWriteLayout,
) -> io::Result<()> {
    // Gemmi✔️✔️: inline void write_cif_to_stream(std::ostream& os, const Document& doc,
    // Gemmi✔️✔️:                                 WriteOptions options=WriteOptions()) {
    // Gemmi✔️✔️:   bool first = true;
    // Gemmi✔️✔️:   for (const Block& block : doc.blocks) {
    // Gemmi✔️✔️:     if (!first)
    // Gemmi✔️✔️:       os.put('\n'); // extra blank line for readability
    // Gemmi✔️✔️:     write_cif_block_to_stream(os, block, options);
    // Gemmi✔️✔️:     first = false;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // Behavior: blocks serialize in stored order; exactly one blank line
    // separates consecutive blocks; a single-block document has no leading
    // or trailing separation.
    // Complexity: one pass over blocks plus c08 costs; no reparse.
    let mut first = true;
    for block in document.blocks() {
        if !first {
            writer.write_all(b"\n")?;
        }
        write_cif_block_to_stream(writer, block, layout)?;
        first = false;
    }
    Ok(())
}

/// Serialize a document straight to an owned UTF-8 string; no reparse.
pub(super) fn cif_document_to_string(
    document: &CifDocument,
    layout: &CifWriteLayout,
) -> io::Result<String> {
    // Behavior: identical bytes reach the returned String as the stream
    // writer emits; invalid UTF-8 becomes io::ErrorKind::InvalidData because
    // every stored value is a Rust String, this is unreachable for documents
    // built by the mutation API and parsed documents already validated UTF-8.
    // Complexity: one serialization pass plus the String move.
    let mut bytes = Vec::new();
    write_cif_document(&mut bytes, document, layout)?;
    String::from_utf8(bytes).map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
}

/// Gemmi `cif::is_text_field` over one raw stored value.
pub(super) fn is_text_field(value: &str) -> bool {
    // Gemmi✔️✔️: inline bool is_text_field(const std::string& val) {
    // Gemmi✔️✔️:   size_t len = val.size();
    // Gemmi✔️✔️:   return len > 2 && val[0] == ';' && (val[len-2] == '\n' || val[len-2] == '\r');
    // Gemmi✔️✔️: }
    // Behavior: a stored value opening with ';' and closing with a newline (or
    // carriage return) before its final byte is emitted as a multiline field.
    // Complexity: O(1) length/byte checks, no allocation.
    let bytes = value.as_bytes();
    bytes.len() > 2 && bytes[0] == b';' && matches!(bytes[bytes.len() - 2], b'\n' | b'\r')
}

/// Gemmi `cif::write_text_field`: copy a multiline field, dropping the `\r`
/// of every `\r\n` pair; every sink error propagates unchanged.
pub(super) fn write_text_field<W: Write + ?Sized>(writer: &mut W, value: &str) -> io::Result<()> {
    // Gemmi✔️✔️: inline void write_text_field(BufOstream& os, const std::string& value) {
    // Gemmi✔️✔️:   for (size_t pos = 0, end = 0; end != std::string::npos; pos = end + 1) {
    // Gemmi✔️✔️:     end = value.find("\r\n", pos);
    // Gemmi✔️✔️:     size_t len = (end == std::string::npos ? value.size() : end) - pos;
    // Gemmi✔️✔️:     os.write(value.c_str() + pos, len);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // Behavior: each "\r\n" collapses to "\n"; lone "\r" bytes are preserved.
    // Complexity: single pass over the value, no intermediate buffers.
    let mut position = 0;
    while let Some(offset) = value[position..].find("\r\n") {
        let end = position + offset;
        writer.write_all(value[position..end].as_bytes())?;
        position = end + 1;
    }
    writer.write_all(value[position..].as_bytes())
}

#[cfg(test)]
mod tests {
    use super::{
        CifBlock, CifDocument, CifItem, CifWriteLayout, cif_document_to_string, is_text_field,
        should_be_separated, write_cif_block_to_stream, write_cif_document, write_out_item,
        write_out_loop, write_out_pair, write_text_field,
    };
    use crate::bio_write::BioMmcifWriteParams;
    use crate::cif::{CifCheckLevel, read_cif_document};
    use std::io::{self, Write};

    struct FailingSink;
    impl Write for FailingSink {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "sink closed"))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn emit(value: &str) -> String {
        let mut buffer = Vec::new();
        write_text_field(&mut buffer, value).unwrap();
        String::from_utf8(buffer).unwrap()
    }

    fn parsed_loop(text: &str) -> super::CifLoop {
        let mut document = read_cif_document(text, "c06", CifCheckLevel::Syntax).unwrap();
        let block = &mut document.blocks[0];
        let index = block
            .items()
            .iter()
            .position(|item| matches!(item, CifItem::Loop(_)))
            .unwrap();
        let CifItem::Loop(loop_) = &block.items()[index] else {
            unreachable!()
        };
        loop_.clone()
    }

    fn emit_loop(loop_: &super::CifLoop, layout: &CifWriteLayout) -> Vec<u8> {
        let mut buffer = Vec::new();
        write_out_loop(&mut buffer, loop_, layout).unwrap();
        buffer
    }

    fn item_at(text: &str, index: usize) -> CifItem {
        let document = read_cif_document(text, "c07", CifCheckLevel::Syntax).unwrap();
        document.blocks[0].items()[index].clone()
    }

    fn first_block(text: &str) -> CifBlock {
        let document = read_cif_document(text, "c08", CifCheckLevel::Syntax).unwrap();
        document.blocks[0].clone()
    }

    fn emit_block(block: &CifBlock, layout: &CifWriteLayout) -> Vec<u8> {
        let mut buffer = Vec::new();
        write_cif_block_to_stream(&mut buffer, block, layout).unwrap();
        buffer
    }

    #[test]
    fn bio_pdbscope_c09_document_block_counts_and_boundaries() {
        let plain = CifWriteLayout::default();
        // The reader requires at least one data_ block, so the zero-block
        // boundary document is constructed directly (child-module field access).
        let empty = CifDocument {
            source: String::new(),
            blocks: Vec::new(),
        };
        assert_eq!(empty.blocks().len(), 0);
        assert_eq!(cif_document_to_string(&empty, &plain).unwrap(), "");
        let one = read_cif_document("data_first\n_a.b 1\n", "c09", CifCheckLevel::Syntax).unwrap();
        assert_eq!(
            cif_document_to_string(&one, &plain).unwrap(),
            "data_first\n_a.b 1\n"
        );
        let two = read_cif_document(
            "data_first\n_a.b 1\ndata_second\n_d.e 2\n",
            "c09",
            CifCheckLevel::Syntax,
        )
        .unwrap();
        assert_eq!(
            cif_document_to_string(&two, &plain).unwrap(),
            "data_first\n_a.b 1\n\ndata_second\n_d.e 2\n"
        );
        let mut buffer = Vec::new();
        write_cif_document(&mut buffer, &two, &plain).unwrap();
        assert_eq!(
            buffer,
            cif_document_to_string(&two, &plain).unwrap().as_bytes()
        );
        let three = read_cif_document(
            "data_one\ndata_two\ndata_three\n",
            "c09",
            CifCheckLevel::Syntax,
        )
        .unwrap();
        assert_eq!(
            cif_document_to_string(&three, &plain).unwrap(),
            "data_one\n\ndata_two\n\ndata_three\n"
        );
        let mut sink = FailingSink;
        let error = write_cif_document(&mut sink, &two, &plain).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
    }

    #[test]
    fn bio_pdbscope_c08_all_boolean_combinations_exact_bytes() {
        let mixed = first_block("data_demo\n_entry.id X\n_entry.name Y\nloop_\n_atom.id\n1\n2\n");
        for prefer_pairs in [false, true] {
            for compact in [false, true] {
                for misuse_hash in [false, true] {
                    let layout = CifWriteLayout {
                        prefer_pairs,
                        compact,
                        misuse_hash,
                        align_pairs: 0,
                        align_loops: 0,
                    };
                    let fence = if misuse_hash { "#\n" } else { "" };
                    let separation = if compact {
                        ""
                    } else if misuse_hash {
                        "#\n"
                    } else {
                        "\n"
                    };
                    let expected = format!(
                        "data_demo\n{fence}_entry.id X\n_entry.name Y\n{separation}loop_\n_atom.id\n1\n2\n{fence}"
                    );
                    assert_eq!(
                        emit_block(&mixed, &layout),
                        expected.as_bytes(),
                        "prefer_pairs={prefer_pairs} compact={compact} misuse_hash={misuse_hash}"
                    );
                }
            }
        }
        let plain = CifWriteLayout::default();
        let single_row = first_block("data_demo\n_a.b 1\nloop_\n_x.y\n7\n");
        assert_eq!(
            emit_block(&single_row, &plain),
            b"data_demo\n_a.b 1\n\nloop_\n_x.y\n7\n"
        );
        let mut pair_form = plain;
        pair_form.prefer_pairs = true;
        assert_eq!(
            emit_block(&single_row, &pair_form),
            b"data_demo\n_a.b 1\n\n_x.y 7\n"
        );
        let empty = first_block("data_demo\n");
        assert_eq!(emit_block(&empty, &plain), b"data_demo\n");
        let mut fenced = plain;
        fenced.misuse_hash = true;
        assert_eq!(emit_block(&empty, &fenced), b"data_demo\n#\n#\n");
        let mut aligned = plain;
        aligned.align_pairs = 33;
        let expected = format!(
            "data_demo\n_entry.id{}X\n_entry.name{}Y\n\nloop_\n_atom.id\n1\n2\n",
            " ".repeat(25),
            " ".repeat(23)
        );
        assert_eq!(emit_block(&mixed, &aligned), expected.as_bytes());
    }

    #[test]
    fn bio_pdbscope_c07_item_dispatch_pair_loop_and_frame() {
        let plain = CifWriteLayout::default();
        let pair = item_at("data_t\n_a.b 1\n", 0);
        let mut buffer = Vec::new();
        write_out_item(&mut buffer, &pair, &plain).unwrap();
        assert_eq!(buffer, b"_a.b 1\n");
        let loop_item = item_at("data_t\nloop_\n_a.id\n1\n", 0);
        let mut buffer = Vec::new();
        write_out_item(&mut buffer, &loop_item, &plain).unwrap();
        assert_eq!(buffer, b"loop_\n_a.id\n1\n");
        let frame = item_at("data_t\nsave_frameOne\n_a.b 1\nsave_\n", 0);
        let mut buffer = Vec::new();
        write_out_item(&mut buffer, &frame, &plain).unwrap();
        assert_eq!(buffer, b"save_frameOne\n_a.b 1\nsave_\n");
    }

    #[test]
    fn bio_pdbscope_c07_category_separation_and_sink_failure() {
        let same = (
            item_at("data_t\n_a.b 1\n_a.c 2\n", 0),
            item_at("data_t\n_a.b 1\n_a.c 2\n", 1),
        );
        assert!(!should_be_separated(&same.0, &same.1));
        let different = (
            item_at("data_t\n_a.b 1\n_d.e 2\n", 0),
            item_at("data_t\n_a.b 1\n_d.e 2\n", 1),
        );
        assert!(should_be_separated(&different.0, &different.1));
        let dotless = (
            item_at("data_t\n_tag1 1\n_tag2 2\n", 0),
            item_at("data_t\n_tag1 1\n_tag2 2\n", 1),
        );
        assert!(!should_be_separated(&dotless.0, &dotless.1));
        let first_dotless = (
            item_at("data_t\n_plain 1\n_a.c 2\n", 0),
            item_at("data_t\n_plain 1\n_a.c 2\n", 1),
        );
        assert!(!should_be_separated(&first_dotless.0, &first_dotless.1));
        let dot_positions = (
            item_at("data_t\n_a.b.c 1\n_ab.c 2\n", 0),
            item_at("data_t\n_a.b.c 1\n_ab.c 2\n", 1),
        );
        assert!(should_be_separated(&dot_positions.0, &dot_positions.1));
        let case = (
            item_at("data_t\n_A.b 1\n_a.c 2\n", 0),
            item_at("data_t\n_A.b 1\n_a.c 2\n", 1),
        );
        assert!(should_be_separated(&case.0, &case.1));
        let mixed = (
            item_at("data_t\n_a.b 1\nloop_\n_x.y\n7\n", 0),
            item_at("data_t\n_a.b 1\nloop_\n_x.y\n7\n", 1),
        );
        assert!(should_be_separated(&mixed.0, &mixed.1));
        let mut sink = FailingSink;
        let error = write_out_item(&mut sink, &same.0, &CifWriteLayout::default()).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
    }

    #[test]
    fn bio_pdbscope_c06_row_counts_prefer_pairs_and_empty() {
        let plain = CifWriteLayout::default();
        let two_rows = parsed_loop("data_t\nloop_\n_a.id\n_a.type\n1 C\n2 NN\n");
        assert_eq!(
            emit_loop(&two_rows, &plain),
            b"loop_\n_a.id\n_a.type\n1 C\n2 NN\n"
        );
        let mut prefer = plain;
        prefer.prefer_pairs = true;
        assert_eq!(
            emit_loop(&two_rows, &prefer),
            b"loop_\n_a.id\n_a.type\n1 C\n2 NN\n"
        );
        let one_row = parsed_loop("data_t\nloop_\n_a.id\n_a.type\n1 C\n");
        assert_eq!(emit_loop(&one_row, &plain), b"loop_\n_a.id\n_a.type\n1 C\n");
        assert_eq!(emit_loop(&one_row, &prefer), b"_a.id 1\n_a.type C\n");
        let mut document =
            read_cif_document("data_t\n_other.id x\n", "c06", CifCheckLevel::Syntax).unwrap();
        let empty = document.blocks[0]
            .init_mmcif_loop("_a", &["id"])
            .unwrap()
            .clone();
        assert_eq!(emit_loop(&empty, &plain), b"");
        assert_eq!(emit_loop(&empty, &prefer), b"");
    }

    #[test]
    fn bio_pdbscope_c06_alignment_widths_and_multiline_cell() {
        let plain = CifWriteLayout::default();
        let ragged = parsed_loop("data_t\nloop_\n_a.id\n_a.type\n1 C\n22 NN\n");
        assert_eq!(
            emit_loop(&ragged, &plain),
            b"loop_\n_a.id\n_a.type\n1 C\n22 NN\n"
        );
        let mut aligned1 = plain;
        aligned1.align_loops = 1;
        assert_eq!(
            emit_loop(&ragged, &aligned1),
            b"loop_\n_a.id\n_a.type\n1 C\n22 NN\n"
        );
        let mut aligned30 = plain;
        aligned30.align_loops = 30;
        assert_eq!(
            emit_loop(&ragged, &aligned30),
            b"loop_\n_a.id\n_a.type\n1  C\n22 NN\n"
        );
        let mut document = read_cif_document("data_t\n", "c06", CifCheckLevel::Syntax).unwrap();
        let multiline = {
            let loop_ = document.blocks[0]
                .init_mmcif_loop("_a", &["id", "type"])
                .unwrap();
            loop_
                .add_row(vec![";A\r\nB\n;".to_string(), "C".to_string()])
                .unwrap();
            loop_.clone()
        };
        assert_eq!(
            emit_loop(&multiline, &plain),
            b"loop_\n_a.id\n_a.type\n;A\nB\n;\nC\n"
        );
    }

    #[test]
    fn bio_pdbscope_c05_value_kinds_alignment_and_boundaries() {
        let plain = CifWriteLayout::default();
        let pair = |tag: &str, value: &str, layout: &CifWriteLayout| {
            let mut buffer = Vec::new();
            write_out_pair(&mut buffer, tag, value, layout).unwrap();
            buffer
        };
        assert_eq!(pair("_entry.id", "1JKZ", &plain), b"_entry.id 1JKZ\n");
        assert_eq!(
            pair("_struct.title", "'two words'", &plain),
            b"_struct.title 'two words'\n"
        );
        assert_eq!(pair("_atom.label", ".", &plain), b"_atom.label .\n");
        assert_eq!(pair("_atom.alt", "?", &plain), b"_atom.alt ?\n");
        let mut aligned = plain;
        aligned.align_pairs = 33;
        let mut expected = b"_entry.id".to_vec();
        expected.push(b' ');
        expected.extend(std::iter::repeat(b' ').take(24));
        expected.extend_from_slice(b"1JKZ\n");
        assert_eq!(pair("_entry.id", "1JKZ", &aligned), expected);
        let long_tag = "_very_long_category_identifier_here.x";
        assert_eq!(long_tag.len(), 37);
        assert_eq!(
            pair(long_tag, "v", &aligned),
            format!("{long_tag} v\n").into_bytes()
        );
        let value_118 = "x".repeat(118);
        let value_119 = "x".repeat(119);
        assert_eq!(
            pair("_t", &value_118, &plain),
            format!("_t {value_118}\n").into_bytes()
        );
        assert_eq!(
            pair("_t", &value_119, &plain),
            format!("_t\n{value_119}\n").into_bytes()
        );
        assert_eq!(
            pair("_struct.details", ";line 1\r\nline 2\n;", &plain),
            b"_struct.details\n;line 1\nline 2\n;\n"
        );
    }

    #[test]
    fn bio_pdbscope_c05_sink_failure_and_layout_passthrough() {
        let mut sink = FailingSink;
        let error =
            write_out_pair(&mut sink, "_entry.id", "1JKZ", &CifWriteLayout::default()).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
        let params = BioMmcifWriteParams {
            group_pdb: true,
            auth_all: true,
            prefer_pairs: true,
            compact: true,
            misuse_hash: true,
            align_pairs: 33,
            align_loops: 30,
            ..BioMmcifWriteParams::default()
        };
        let layout = CifWriteLayout::from(&params);
        assert_eq!(
            layout,
            CifWriteLayout {
                prefer_pairs: true,
                compact: true,
                misuse_hash: true,
                align_pairs: 33,
                align_loops: 30,
            }
        );
        assert_eq!(
            CifWriteLayout::from(&BioMmcifWriteParams::default()),
            CifWriteLayout::default()
        );
    }

    #[test]
    fn bio_pdbscope_c04_text_field_shape_and_newline_collapsing() {
        assert!(!is_text_field(""));
        assert!(!is_text_field(";x"));
        assert!(!is_text_field(";\n"));
        assert!(is_text_field(";\n;"));
        assert!(is_text_field(";\r;"));
        assert!(!is_text_field("a\n;"));
        assert_eq!(emit(""), "");
        assert_eq!(emit("plain"), "plain");
        assert_eq!(emit("line 1\nline 2"), "line 1\nline 2");
        assert_eq!(emit("line 1\r\nline 2"), "line 1\nline 2");
        assert_eq!(emit("a\r\nb\r\nc"), "a\nb\nc");
        assert_eq!(emit("only\rleft"), "only\rleft");
        assert_eq!(emit("text\n;\n;"), "text\n;\n;");
        let multiline = ";line 1\nO'Brien says \"hi\"\nline 3\n;";
        assert!(is_text_field(multiline));
        assert_eq!(emit(multiline), multiline);
        assert_eq!(emit(";\r\n;"), ";\n;");
        assert_eq!(emit("ends with newline\n"), "ends with newline\n");
    }

    #[test]
    fn bio_pdbscope_c04_sink_failure_propagates_before_any_byte_is_lost() {
        let mut buffer = Vec::new();
        buffer.extend_from_slice(b"head\n");
        let error = write_text_field(&mut FailingSink, ";data\r\nmore\n;").unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
        assert_eq!(error.to_string(), "sink closed");
        let mut partial = Vec::new();
        let _ = write_text_field(&mut partial, ";keep\r\nfail");
        assert_eq!(partial, b";keep\nfail".as_slice());
        let _ = buffer;
    }
}
