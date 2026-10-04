mod parser;

use parser::{JsonArena, JsonArenaValue, JsonCursor, JsonError};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum MmjsonReadError {
    Parse {
        source_name: String,
        error: JsonError,
    },
    Structure {
        source_name: String,
        error: MmjsonStructureError,
    },
}

impl std::fmt::Display for MmjsonReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Parse { source_name, error } => {
                write!(f, "{source_name}:{} error: {error}", error.line())
            }
            Self::Structure { source_name, error } => write!(f, "{source_name}: {error}"),
        }
    }
}

impl std::error::Error for MmjsonReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Parse { error, .. } => Some(error),
            Self::Structure { error, .. } => Some(error),
        }
    }
}

pub(crate) fn read_mmjson_insitu(
    input: &[u8],
    name: &str,
) -> Result<super::CifDocument, MmjsonReadError> {
    // Gemmi❗✔️: Document read_mmjson_insitu(char* buffer, size_t size, const std::string& name) {
    // Gemmi❗✔️:   Document doc;
    // Gemmi❗✔️:   sajson::document json = sajson::parse(sajson::dynamic_allocation(),
    // Gemmi❗✔️:                                     sajson::mutable_string_view(size, buffer));
    // Gemmi❗✔️:   if (!json.is_valid())
    // Gemmi❗✔️:     fail(name + ":", std::to_string(json.get_error_line()), " error: ",
    // Gemmi❗✔️:          json.get_error_message_as_string());
    // Gemmi❗✔️:   fill_document_from_sajson(doc, json);
    // Gemmi❗✔️:   doc.source = name;
    // Gemmi❗✔️:   return doc;
    // Gemmi❗✔️: }
    // Behavior: parser error precedes document lowering, and source name is
    // installed only after successful lowering; errors retain typed causes.
    // Complexity: one owned mutable input buffer and one arena, then a linear
    // block/category traversal; no JSON-to-text CIF round trip.
    let mut cursor = JsonCursor::new(input);
    let (arena, root) = cursor
        .parse_document()
        .map_err(|error| MmjsonReadError::Parse {
            source_name: name.to_owned(),
            error,
        })?;
    let mut document = fill_document_structure(
        &arena,
        root,
        cursor.bytes(),
        |block, category_name, category_id| {
            fill_category_rows(block, category_name, category_id, &arena, cursor.bytes())
        },
    )
    .map_err(|error| MmjsonReadError::Structure {
        source_name: name.to_owned(),
        error,
    })?;
    document.source = name.to_owned();
    Ok(document)
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmjsonStructureError {
    RootNotObject,
    BlockKeyPrefix {
        block: usize,
    },
    BlockNotObject {
        block: usize,
    },
    InvalidCategory {
        block: usize,
        category: usize,
    },
    ExpectedArray {
        column: usize,
        kind: &'static str,
    },
    ArrayLength {
        column: usize,
        expected: usize,
        actual: usize,
    },
    RowCountOverflow,
    MissingArenaValue,
    Text(MmjsonValueError),
}

impl std::fmt::Display for MmjsonStructureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RootNotObject => f.write_str("not mmJSON - the root is not of type object"),
            Self::BlockKeyPrefix { block } => {
                write!(f, "top-level key {block} should start with data_")
            }
            Self::BlockNotObject { block } => write!(f, "block {block} is not an object"),
            Self::InvalidCategory { block, category } => {
                write!(f, "invalid category {category} in block {block}")
            }
            Self::ExpectedArray { column, kind } => {
                write!(f, "column {column}: Expected array, got {kind}")
            }
            Self::ArrayLength {
                column,
                expected,
                actual,
            } => write!(
                f,
                "column {column}: Expected array of length {expected} not {actual}"
            ),
            Self::RowCountOverflow => f.write_str("mmJSON row count overflow"),
            Self::MissingArenaValue => f.write_str("missing private mmJSON arena value"),
            Self::Text(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for MmjsonStructureError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Text(error) => Some(error),
            _ => None,
        }
    }
}

fn fill_category_rows(
    block: &mut super::CifBlock,
    category_name: &str,
    category_id: usize,
    arena: &JsonArena,
    input: &[u8],
) -> Result<(), MmjsonStructureError> {
    // Gemmi❗✔️:       size_t cif_cols = category.get_length();
    // Gemmi❗✔️:       size_t cif_rows = category.get_object_value(0).get_length();
    // Gemmi❗✔️:       if (cif_rows > 1) {
    // Gemmi❗✔️:         items.emplace_back(LoopArg{});
    // Gemmi❗✔️:         Loop& loop = items.back().loop;
    // Gemmi❗✔️:         loop.tags.reserve(cif_cols);
    // Gemmi❗✔️:         loop.values.resize(cif_cols * cif_rows);
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       for (size_t j = 0; j != cif_cols; ++j) {
    // Gemmi❗✔️:         std::string tag = category_name + category.get_object_key(j).as_string();
    // Gemmi❗✔️:         sajson::value arr = category.get_object_value(j);
    // Gemmi❗✔️:         if (arr.get_type() != sajson::TYPE_ARRAY)
    // Gemmi❗✔️:           fail("Expected array, got " + json_type_as_string(arr.get_type()));
    // Gemmi❗✔️:         if (arr.get_length() != cif_rows)
    // Gemmi❗✔️:           fail("Expected array of length ", std::to_string(cif_rows), " not ",
    // Gemmi❗✔️:                std::to_string(arr.get_length()));
    // Gemmi❗✔️:         if (cif_rows == 1) {
    // Gemmi❗✔️:           items.emplace_back(tag, as_cif_value(arr.get_array_element(0)));
    // Gemmi❗✔️:         } else if (cif_rows != 0) {
    // Gemmi❗✔️:           Loop& loop = items.back().loop;
    // Gemmi❗✔️:           loop.tags.emplace_back(std::move(tag));
    // Gemmi❗✔️:           for (size_t k = 0; k != cif_rows; ++k)
    // Gemmi❗✔️:             loop.values[j + k*cif_cols] = as_cif_value(arr.get_array_element(k));
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:       }
    // Behavior: complete source 0/1/many branch, order and first-error checks;
    // malformed private spans and arithmetic overflow are structured errors.
    // Complexity: O(cols * rows) conversion and one row-major values buffer.
    let JsonArenaValue::Object(columns) = arena
        .value(category_id)
        .ok_or(MmjsonStructureError::MissingArenaValue)?
    else {
        return Err(MmjsonStructureError::MissingArenaValue);
    };
    let cif_cols = columns.len();
    let JsonArenaValue::Array(first) = arena
        .object_entry(category_id, 0)
        .ok_or(MmjsonStructureError::MissingArenaValue)?
        .1
    else {
        return Err(MmjsonStructureError::MissingArenaValue);
    };
    let cif_rows = first.len();
    let mut loop_item = if cif_rows > 1 {
        let count = cif_cols
            .checked_mul(cif_rows)
            .ok_or(MmjsonStructureError::RowCountOverflow)?;
        Some(super::CifLoop {
            tags: Vec::with_capacity(cif_cols),
            values: vec![super::CifValue::new(String::new(), 0, 0); count],
            line: 0,
        })
    } else {
        None
    };
    for column in 0..cif_cols {
        let (tag_span, array_id) = arena
            .object_entry_id(category_id, column)
            .ok_or(MmjsonStructureError::MissingArenaValue)?;
        let tag_bytes = input
            .get(tag_span.clone())
            .ok_or(MmjsonStructureError::Text(MmjsonValueError::InvalidSpan {
                start: tag_span.start,
                end: tag_span.end,
            }))?;
        let tag_text = std::str::from_utf8(tag_bytes).map_err(|error| {
            MmjsonStructureError::Text(MmjsonValueError::InvalidUtf8 {
                valid_up_to: tag_span.start + error.valid_up_to(),
            })
        })?;
        let tag = format!("{category_name}{tag_text}");
        let value = arena
            .value(array_id)
            .ok_or(MmjsonStructureError::MissingArenaValue)?;
        let JsonArenaValue::Array(range) = value else {
            return Err(MmjsonStructureError::ExpectedArray {
                column,
                kind: json_kind_name(value),
            });
        };
        if range.len() != cif_rows {
            return Err(MmjsonStructureError::ArrayLength {
                column,
                expected: cif_rows,
                actual: range.len(),
            });
        }
        if cif_rows == 1 {
            let cell = arena
                .array_element(array_id, 0)
                .ok_or(MmjsonStructureError::MissingArenaValue)?;
            block.items.push(super::CifItem::Pair(super::CifPair {
                tag,
                value: Some(super::CifValue::new(as_cif_cell(cell, arena, input)?, 0, 0)),
                line: 0,
            }));
        } else if let Some(loop_item) = loop_item.as_mut() {
            loop_item.tags.push(tag);
            for row in 0..cif_rows {
                let cell = arena
                    .array_element(array_id, row)
                    .ok_or(MmjsonStructureError::MissingArenaValue)?;
                loop_item.values[column + row * cif_cols] =
                    super::CifValue::new(as_cif_cell(cell, arena, input)?, 0, 0);
            }
        }
    }
    if let Some(loop_item) = loop_item {
        block.items.push(super::CifItem::Loop(loop_item));
    }
    Ok(())
}

fn json_kind_name(value: &JsonArenaValue) -> &'static str {
    match value {
        JsonArenaValue::Null => "<null>",
        JsonArenaValue::Boolean(false) => "<false>",
        JsonArenaValue::Boolean(true) => "<true>",
        JsonArenaValue::Number(_) => "<double>",
        JsonArenaValue::String(_) => "<string>",
        JsonArenaValue::Array(_) => "<array>",
        JsonArenaValue::Object(_) => "<object>",
    }
}

fn as_cif_cell(
    value: &JsonArenaValue,
    arena: &JsonArena,
    input: &[u8],
) -> Result<String, MmjsonStructureError> {
    if let JsonArenaValue::Array(range) = value {
        let first = if range.is_empty() {
            None
        } else {
            arena.array_element_by_range(range, 0)
        };
        array_to_cif_value(range.len(), first, input).map_err(MmjsonStructureError::Text)
    } else {
        scalar_to_cif_value(value, input).map_err(MmjsonStructureError::Text)
    }
}

fn fill_document_structure<F>(
    arena: &JsonArena,
    root: usize,
    input: &[u8],
    mut on_category: F,
) -> Result<super::CifDocument, MmjsonStructureError>
where
    F: FnMut(&mut super::CifBlock, &str, usize) -> Result<(), MmjsonStructureError>,
{
    // Gemmi❗✔️: static void fill_document_from_sajson(Document& d, const sajson::document& s) {
    // Gemmi❗✔️:   // assuming mmJSON here, we'll add handling of CIF-JSON later on
    // Gemmi❗✔️:   sajson::value root = s.get_root();
    // Gemmi❗✔️:   if (root.get_type() != sajson::TYPE_OBJECT)
    // Gemmi❗✔️:     fail("not mmJSON - the root is not of type object");
    // Gemmi❗✔️:   for (size_t block_index = 0; block_index < root.get_length(); ++block_index) {
    // Gemmi❗✔️:     std::string block_name = root.get_object_key(block_index).as_string();
    // Gemmi❗✔️:     if (!starts_with(block_name, "data_"))
    // Gemmi❗✔️:       fail("not mmJSON - top level key should start with data_\n"
    // Gemmi❗✔️:            "(if you use gemmi-cif2json to write JSON, use -m for mmJSON)");
    // Gemmi❗✔️:     d.blocks.emplace_back(block_name.substr(5));
    // Gemmi❗✔️:     std::vector<Item>& items = d.blocks[block_index].items;
    // Gemmi❗✔️:     sajson::value top = root.get_object_value(block_index);
    // Gemmi❗✔️:     if (top.get_type() != sajson::TYPE_OBJECT)
    // Gemmi❗✔️:       fail("");
    // Gemmi❗✔️:     for (size_t i = 0; i != top.get_length(); ++i) {
    // Gemmi❗✔️:       std::string category_name = "_" + top.get_object_key(i).as_string() + ".";
    // Gemmi❗✔️:       sajson::value category = top.get_object_value(i);
    // Gemmi❗✔️:       if (category.get_type() != sajson::TYPE_OBJECT ||
    // Gemmi❗✔️:           category.get_length() == 0 ||
    // Gemmi❗✔️:           category.get_object_value(0).get_type() != sajson::TYPE_ARRAY)
    // Gemmi❗✔️:         fail("");
    // The source's column/row body remains J17 and is not certified here.
    // Behavior: each ordered object entry is visited once; exact data_ casing,
    // category shape and first-column array precede downstream lowering.
    // Complexity: O(blocks + categories), no map, sort or deep copy; owned
    // block names and category names mirror source string construction.
    let Some(JsonArenaValue::Object(root_range)) = arena.value(root) else {
        return Err(MmjsonStructureError::RootNotObject);
    };
    let mut document = super::CifDocument {
        source: String::new(),
        blocks: Vec::new(),
    };
    for block_index in 0..root_range.len() {
        let (block_key, block_id) = arena
            .object_entry_id(root, block_index)
            .ok_or(MmjsonStructureError::MissingArenaValue)?;
        let block_bytes = input
            .get(block_key.clone())
            .ok_or(MmjsonStructureError::Text(MmjsonValueError::InvalidSpan {
                start: block_key.start,
                end: block_key.end,
            }))?;
        let block_name = std::str::from_utf8(block_bytes).map_err(|error| {
            MmjsonStructureError::Text(MmjsonValueError::InvalidUtf8 {
                valid_up_to: block_key.start + error.valid_up_to(),
            })
        })?;
        let Some(name) = block_name.strip_prefix("data_") else {
            return Err(MmjsonStructureError::BlockKeyPrefix { block: block_index });
        };
        let Some(JsonArenaValue::Object(categories)) = arena.value(block_id) else {
            return Err(MmjsonStructureError::BlockNotObject { block: block_index });
        };
        let mut block = super::CifBlock {
            name: name.to_owned(),
            items: Vec::new(),
            line: 0,
        };
        for category_index in 0..categories.len() {
            let (category_key, category_id) = arena
                .object_entry_id(block_id, category_index)
                .ok_or(MmjsonStructureError::MissingArenaValue)?;
            let category_bytes =
                input
                    .get(category_key.clone())
                    .ok_or(MmjsonStructureError::Text(MmjsonValueError::InvalidSpan {
                        start: category_key.start,
                        end: category_key.end,
                    }))?;
            let category_text = std::str::from_utf8(category_bytes).map_err(|error| {
                MmjsonStructureError::Text(MmjsonValueError::InvalidUtf8 {
                    valid_up_to: category_key.start + error.valid_up_to(),
                })
            })?;
            let Some(JsonArenaValue::Object(columns)) = arena.value(category_id) else {
                return Err(MmjsonStructureError::InvalidCategory {
                    block: block_index,
                    category: category_index,
                });
            };
            if columns.is_empty()
                || !matches!(
                    arena.object_entry(category_id, 0),
                    Some((_, JsonArenaValue::Array(_)))
                )
            {
                return Err(MmjsonStructureError::InvalidCategory {
                    block: block_index,
                    category: category_index,
                });
            }
            let category_name = format!("_{category_text}.");
            on_category(&mut block, &category_name, category_id)?;
        }
        document.blocks.push(block);
    }
    Ok(document)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MmjsonUnexpectedType {
    Null,
    False,
    True,
    Array,
    Object,
}

impl MmjsonUnexpectedType {
    fn source_name(self) -> &'static str {
        // Gemmi❗✔️: static std::string json_type_as_string(sajson::type t) {
        // Gemmi❗✔️:   switch (t) {
        // Gemmi❗✔️:     case sajson::TYPE_INTEGER: return "<integer>";
        // Gemmi❗✔️:     case sajson::TYPE_DOUBLE:  return "<double>";
        // Gemmi❗✔️:     case sajson::TYPE_NULL:    return "<null>";
        // Gemmi❗✔️:     case sajson::TYPE_FALSE:   return "<false>";
        // Gemmi❗✔️:     case sajson::TYPE_TRUE:    return "<true>";
        // Gemmi❗✔️:     case sajson::TYPE_STRING:  return "<string>";
        // Gemmi❗✔️:     case sajson::TYPE_ARRAY:   return "<array>";
        // Gemmi❗✔️:     case sajson::TYPE_OBJECT:  return "<object>";
        // Gemmi❗✔️:     default:           return "<unknown type>";
        // Gemmi❗✔️:   }
        // Gemmi❗✔️: }
        // Behavior: the scalar converter only routes array/object failures
        // here; all other types are handled directly in the owning match.
        // Complexity: O(1) dispatch, as in the source switch.
        match self {
            Self::Null => "<null>",
            Self::False => "<false>",
            Self::True => "<true>",
            Self::Array => "<array>",
            Self::Object => "<object>",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmjsonValueError {
    UnexpectedType(MmjsonUnexpectedType),
    InvalidArrayElement(MmjsonUnexpectedType),
    MissingArrayFirstElement,
    InvalidSpan { start: usize, end: usize },
    InvalidUtf8 { valid_up_to: usize },
}

impl std::fmt::Display for MmjsonValueError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnexpectedType(kind) => {
                write!(f, "Unexpected {} as value in JSON.", kind.source_name())
            }
            Self::InvalidArrayElement(kind) => {
                write!(
                    f,
                    "array element {} has no source-defined as_string value",
                    kind.source_name()
                )
            }
            Self::MissingArrayFirstElement => f.write_str("missing first mmJSON array element"),
            Self::InvalidSpan { start, end } => {
                write!(f, "invalid mmJSON byte span {start}..{end}")
            }
            Self::InvalidUtf8 { valid_up_to } => {
                write!(f, "invalid UTF-8 in mmJSON value at byte {valid_up_to}")
            }
        }
    }
}

fn array_to_cif_value(
    length: usize,
    first: Option<&JsonArenaValue>,
    input: &[u8],
) -> Result<String, MmjsonValueError> {
    // Gemmi❗✔️:     case sajson::TYPE_ARRAY: {
    // Gemmi❗✔️:       std::string s;
    // Gemmi❗✔️:       for (size_t i = 0; i < val.get_length(); ++i) {
    // Gemmi❗✔️:         if (i != 0)
    // Gemmi❗✔️:           s += ' ';
    // Gemmi❗✔️:         s += val.get_array_element(0).as_string();
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       return quote(s);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:         std::string as_string() const {
    // Gemmi❗✔️: #ifndef SAJSON_NUMBERS_AS_STRINGS
    // Gemmi❗✔️:             assert_type(TYPE_STRING);
    // Gemmi❗✔️: #else
    // Gemmi❗✔️:             assert_type_2(TYPE_STRING, TYPE_DOUBLE);
    // Gemmi❗✔️: #endif
    // Gemmi❗✔️:             return std::string(text + payload[0], text + payload[1]);
    // Gemmi❗✔️:         }
    // Behavior: only string/raw-number first elements have source-defined
    // as_string; invalid source assertions become typed errors, not UB/panic.
    // Empty arrays skip the first-element read and quote an empty string.
    // Complexity: O(length * first-byte-length), matching repeated append;
    // capacity growth is amortized, with one final canonical CIF quote.
    let first_text = if length == 0 {
        ""
    } else {
        let value = first.ok_or(MmjsonValueError::MissingArrayFirstElement)?;
        let span = match value {
            JsonArenaValue::String(span) | JsonArenaValue::Number(span) => span,
            JsonArenaValue::Null => {
                return Err(MmjsonValueError::InvalidArrayElement(
                    MmjsonUnexpectedType::Null,
                ));
            }
            JsonArenaValue::Boolean(false) => {
                return Err(MmjsonValueError::InvalidArrayElement(
                    MmjsonUnexpectedType::False,
                ));
            }
            JsonArenaValue::Boolean(true) => {
                return Err(MmjsonValueError::InvalidArrayElement(
                    MmjsonUnexpectedType::True,
                ));
            }
            JsonArenaValue::Array(_) => {
                return Err(MmjsonValueError::InvalidArrayElement(
                    MmjsonUnexpectedType::Array,
                ));
            }
            JsonArenaValue::Object(_) => {
                return Err(MmjsonValueError::InvalidArrayElement(
                    MmjsonUnexpectedType::Object,
                ));
            }
        };
        let bytes = input
            .get(span.clone())
            .ok_or(MmjsonValueError::InvalidSpan {
                start: span.start,
                end: span.end,
            })?;
        std::str::from_utf8(bytes).map_err(|error| MmjsonValueError::InvalidUtf8 {
            valid_up_to: span.start + error.valid_up_to(),
        })?
    };
    let mut combined = String::new();
    for index in 0..length {
        if index != 0 {
            combined.push(' ');
        }
        combined.push_str(first_text);
    }
    Ok(super::quote_cif_value(combined))
}

impl std::error::Error for MmjsonValueError {}

fn scalar_to_cif_value(value: &JsonArenaValue, input: &[u8]) -> Result<String, MmjsonValueError> {
    // Gemmi❗✔️: static std::string as_cif_value(const sajson::value& val) {
    // Gemmi❗✔️:   switch (val.get_type()) {
    // Gemmi❗✔️:     case sajson::TYPE_DOUBLE:
    // Gemmi❗✔️:       return val.as_string();
    // Gemmi❗✔️:     case sajson::TYPE_NULL:
    // Gemmi❗✔️:       return "?";
    // Gemmi❗✔️:     // mmJSON files from PDBj (this format has no spec) have special support
    // Gemmi❗✔️:     // for boolean YES|NO, which is used only in category _em_specimen.
    // Gemmi❗✔️:     // IMO it's a bad idea, but we must handle it if we want to read mmJSON.
    // Gemmi❗✔️:     case sajson::TYPE_FALSE:
    // Gemmi❗✔️:       return "NO";  // "." in CIF-JSON
    // Gemmi❗✔️:     case sajson::TYPE_TRUE:
    // Gemmi❗✔️:       return "YES";
    // Gemmi❗✔️:     case sajson::TYPE_STRING:
    // Gemmi❗✔️:       return quote(val.as_string());
    // Gemmi❗✔️:     // Another undocumented feature of mmJSON: arrays as values.
    // Gemmi❗✔️:     // It seems that obscure types int-range and float-range are converted to
    // Gemmi❗✔️:     // 2-element arrays. But not only. link_entity_pdbjplus.db_accession has
    // Gemmi❗✔️:     // arrays with strings.
    // Gemmi❗✔️:     case sajson::TYPE_ARRAY: {
    // Gemmi❗✔️:       std::string s;
    // Gemmi❗✔️:       for (size_t i = 0; i < val.get_length(); ++i) {
    // Gemmi❗✔️:         if (i != 0)
    // Gemmi❗✔️:           s += ' ';
    // Gemmi❗✔️:         s += val.get_array_element(0).as_string();
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       return quote(s);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     default:
    // Gemmi❗✔️:       fail("Unexpected ", json_type_as_string(val.get_type()), " as value in JSON.");
    // Gemmi❗✔️:       return "";
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Behavior: this helper implements only the scalar arms; the source array
    // arm is a separate J15 task. Rust text requires explicit UTF-8 rejection
    // for decoded bytes that sajson can retain but CIF String cannot represent.
    // Complexity: one owned copy for text, then the existing linear CIF quote.
    match value {
        JsonArenaValue::Number(span) | JsonArenaValue::String(span) => {
            let bytes = input
                .get(span.clone())
                .ok_or(MmjsonValueError::InvalidSpan {
                    start: span.start,
                    end: span.end,
                })?;
            let text =
                std::str::from_utf8(bytes).map_err(|error| MmjsonValueError::InvalidUtf8 {
                    valid_up_to: span.start + error.valid_up_to(),
                })?;
            if matches!(value, JsonArenaValue::String(_)) {
                Ok(super::quote_cif_value(text.to_owned()))
            } else {
                Ok(text.to_owned())
            }
        }
        JsonArenaValue::Null => Ok("?".to_owned()),
        JsonArenaValue::Boolean(false) => Ok("NO".to_owned()),
        JsonArenaValue::Boolean(true) => Ok("YES".to_owned()),
        JsonArenaValue::Array(_) => Err(MmjsonValueError::UnexpectedType(
            MmjsonUnexpectedType::Array,
        )),
        JsonArenaValue::Object(_) => Err(MmjsonValueError::UnexpectedType(
            MmjsonUnexpectedType::Object,
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::parser::{JsonArenaValue, JsonCursor};
    use super::{
        MmjsonReadError, MmjsonStructureError, MmjsonUnexpectedType, MmjsonValueError,
        array_to_cif_value, fill_category_rows, fill_document_structure, read_mmjson_insitu,
        scalar_to_cif_value,
    };

    fn j16_structure(
        input: &[u8],
    ) -> Result<(super::super::CifDocument, Vec<String>), MmjsonStructureError> {
        let mut cursor = JsonCursor::new(input);
        let (arena, root) = cursor.parse_document().unwrap();
        let mut categories = Vec::new();
        let document = fill_document_structure(&arena, root, cursor.bytes(), |_, name, _| {
            categories.push(name.to_owned());
            Ok(())
        })?;
        Ok((document, categories))
    }

    fn j17_document(input: &[u8]) -> Result<super::super::CifDocument, MmjsonStructureError> {
        let mut cursor = JsonCursor::new(input);
        let (arena, root) = cursor.parse_document().unwrap();
        fill_document_structure(&arena, root, cursor.bytes(), |block, name, category| {
            fill_category_rows(block, name, category, &arena, cursor.bytes())
        })
    }

    #[test]
    fn bio_read_j18_complete_entry_preserves_source_order_and_caller_bytes() {
        use super::super::CifItem;
        let input = br#"{"data_x":{"cat":{"a":[1.2300E+04,2],"b":[["first","other"],null]},"one":{"value":["x"]}}}"#.to_vec();
        let before = input.clone();
        let document = read_mmjson_insitu(&input, "sample.mmjson").unwrap();
        assert_eq!(input, before);
        assert_eq!(document.source(), "sample.mmjson");
        assert_eq!(document.blocks()[0].name(), "x");
        let CifItem::Loop(rows) = &document.blocks()[0].items()[0] else {
            panic!("source requires loop");
        };
        assert_eq!(rows.tags(), ["_cat.a", "_cat.b"]);
        assert_eq!(
            rows.values()
                .iter()
                .map(|value| value.raw())
                .collect::<Vec<_>>(),
            ["1.2300E+04", "'first first'", "2", "?"]
        );
        assert!(
            matches!(&document.blocks()[0].items()[1], CifItem::Pair(pair) if pair.tag() == "_one.value" && pair.value().unwrap().raw() == "x")
        );
    }

    #[test]
    fn bio_read_j18_errors_keep_source_name_and_typed_cause() {
        let parse = read_mmjson_insitu(b"{", "broken.mmjson").unwrap_err();
        assert!(matches!(parse, MmjsonReadError::Parse { .. }));
        assert!(parse.to_string().starts_with("broken.mmjson:1 error:"));
        let structure = read_mmjson_insitu(br#"{"Data_x":{}}"#, "shape.mmjson").unwrap_err();
        assert!(matches!(
            structure,
            MmjsonReadError::Structure {
                error: MmjsonStructureError::BlockKeyPrefix { block: 0 },
                ..
            }
        ));
        assert!(structure.to_string().contains("shape.mmjson"));
    }

    #[test]
    fn bio_read_j17_zero_one_many_rows_and_row_major_values() {
        use super::super::CifItem;
        let zero = j17_document(br#"{"data_x":{"cat":{"a":[],"b":[]}}}"#).unwrap();
        assert!(zero.blocks()[0].items().is_empty());
        let one = j17_document(br#"{"data_x":{"cat":{"a":["word"],"b":[1.20]}}}"#).unwrap();
        assert!(
            matches!(&one.blocks()[0].items()[0], CifItem::Pair(pair) if pair.tag() == "_cat.a" && pair.value().unwrap().raw() == "word")
        );
        assert!(
            matches!(&one.blocks()[0].items()[1], CifItem::Pair(pair) if pair.tag() == "_cat.b" && pair.value().unwrap().raw() == "1.20")
        );
        let many = j17_document(
            br#"{"data_x":{"cat":{"a":["x","y"],"b":[0,true],"b":[null,["first","other"]]}}}"#,
        )
        .unwrap();
        let CifItem::Loop(rows) = &many.blocks()[0].items()[0] else {
            panic!("source requires one loop");
        };
        assert_eq!(rows.tags(), ["_cat.a", "_cat.b", "_cat.b"]);
        assert_eq!(
            rows.values()
                .iter()
                .map(|value| value.raw())
                .collect::<Vec<_>>(),
            ["x", "0", "?", "y", "YES", "'first first'"]
        );
        assert_eq!(many.blocks()[0].items().len(), 1);
    }

    #[test]
    fn bio_read_j17_later_column_shape_and_length_fail_in_source_order() {
        assert_eq!(
            j17_document(br#"{"data_x":{"cat":{"a":[1,2],"b":false}}}"#).unwrap_err(),
            MmjsonStructureError::ExpectedArray {
                column: 1,
                kind: "<false>"
            }
        );
        assert_eq!(
            j17_document(br#"{"data_x":{"cat":{"a":[1,2],"b":[3]}}}"#).unwrap_err(),
            MmjsonStructureError::ArrayLength {
                column: 1,
                expected: 2,
                actual: 1
            }
        );
        assert_eq!(
            j17_document(br#"{"data_x":{"cat":{"a":[],"b":[3]}}}"#).unwrap_err(),
            MmjsonStructureError::ArrayLength {
                column: 1,
                expected: 0,
                actual: 1
            }
        );
    }

    #[test]
    fn bio_read_j16_root_prefix_and_block_shape_follow_source() {
        assert_eq!(
            j16_structure(b"[]").unwrap_err(),
            MmjsonStructureError::RootNotObject
        );
        for input in [b"{\"Data_a\":{}}".as_slice(), b"{\"x\":{}}".as_slice()] {
            assert_eq!(
                j16_structure(input).unwrap_err(),
                MmjsonStructureError::BlockKeyPrefix { block: 0 }
            );
        }
        assert_eq!(
            j16_structure(b"{\"data_a\":[]}").unwrap_err(),
            MmjsonStructureError::BlockNotObject { block: 0 }
        );
        assert_eq!(
            j16_structure(b"{\"data_a\":{}}").unwrap().0.blocks()[0].name(),
            "a"
        );
    }

    #[test]
    fn bio_read_j16_category_shape_and_first_column_follow_source() {
        for input in [
            b"{\"data_a\":{\"cat\":[]}}".as_slice(),
            b"{\"data_a\":{\"cat\":{}}}".as_slice(),
            b"{\"data_a\":{\"cat\":{\"col\":null}}}".as_slice(),
        ] {
            assert_eq!(
                j16_structure(input).unwrap_err(),
                MmjsonStructureError::InvalidCategory {
                    block: 0,
                    category: 0
                }
            );
        }
        assert_eq!(
            j16_structure(b"{\"data_a\":{\"cat\":{\"col\":[]}}}")
                .unwrap()
                .1,
            ["_cat."]
        );
    }

    #[test]
    fn bio_read_j16_duplicate_blocks_and_categories_keep_source_order() {
        let input =
            br#"{"data_b":{"z":{"x":[]},"z":{"x":[]}},"data_a":{"q":{"x":[]}},"data_b":{}}"#;
        let (document, categories) = j16_structure(input).unwrap();
        assert_eq!(
            document
                .blocks()
                .iter()
                .map(|block| block.name())
                .collect::<Vec<_>>(),
            ["b", "a", "b"]
        );
        assert_eq!(categories, ["_z.", "_z.", "_q."]);
    }

    #[test]
    fn bio_read_j15_empty_singleton_and_repeated_first_element() {
        assert_eq!(array_to_cif_value(0, None, b"").unwrap(), "''");
        assert_eq!(
            array_to_cif_value(1, Some(&JsonArenaValue::String(0..5)), b"alpha").unwrap(),
            "alpha"
        );
        assert_eq!(
            array_to_cif_value(2, Some(&JsonArenaValue::String(0..5)), b"alpha beta").unwrap(),
            "'alpha alpha'"
        );
        assert_eq!(
            array_to_cif_value(3, Some(&JsonArenaValue::String(0..1)), b"x y z").unwrap(),
            "'x x x'"
        );
        assert_eq!(
            array_to_cif_value(2, Some(&JsonArenaValue::Number(0..7)), b"-0.0e+2 99").unwrap(),
            "'-0.0e+2 -0.0e+2'"
        );
    }

    #[test]
    fn bio_read_j15_invalid_first_element_has_typed_safe_boundary() {
        for (first, kind) in [
            (JsonArenaValue::Null, MmjsonUnexpectedType::Null),
            (JsonArenaValue::Boolean(false), MmjsonUnexpectedType::False),
            (JsonArenaValue::Boolean(true), MmjsonUnexpectedType::True),
            (JsonArenaValue::Array(0..0), MmjsonUnexpectedType::Array),
            (JsonArenaValue::Object(0..0), MmjsonUnexpectedType::Object),
        ] {
            assert_eq!(
                array_to_cif_value(2, Some(&first), b"").unwrap_err(),
                MmjsonValueError::InvalidArrayElement(kind)
            );
        }
        assert_eq!(
            array_to_cif_value(1, None, b"").unwrap_err(),
            MmjsonValueError::MissingArrayFirstElement
        );
        assert_eq!(
            array_to_cif_value(1, Some(&JsonArenaValue::String(0..1)), b"\xff").unwrap_err(),
            MmjsonValueError::InvalidUtf8 { valid_up_to: 0 }
        );
    }

    #[test]
    fn bio_read_j14_null_boolean_number_and_null_like_string_values() {
        assert_eq!(
            scalar_to_cif_value(&JsonArenaValue::Null, b"").unwrap(),
            "?"
        );
        assert_eq!(
            scalar_to_cif_value(&JsonArenaValue::Boolean(false), b"").unwrap(),
            "NO"
        );
        assert_eq!(
            scalar_to_cif_value(&JsonArenaValue::Boolean(true), b"").unwrap(),
            "YES"
        );
        for (bytes, expected) in [
            (b"".as_slice(), "''"),
            (b".".as_slice(), "'.'"),
            (b"?".as_slice(), "'?'"),
            (b"word".as_slice(), "word"),
            (b"a b".as_slice(), "'a b'"),
        ] {
            assert_eq!(
                scalar_to_cif_value(&JsonArenaValue::String(0..bytes.len()), bytes).unwrap(),
                expected
            );
        }
        for bytes in [
            b"-0".as_slice(),
            b"1.2300E+04".as_slice(),
            b"0.000e-02".as_slice(),
        ] {
            assert_eq!(
                scalar_to_cif_value(&JsonArenaValue::Number(0..bytes.len()), bytes)
                    .unwrap()
                    .as_bytes(),
                bytes
            );
        }
    }

    #[test]
    fn bio_read_j14_quote_priority_errors_and_exact_spans() {
        for (bytes, expected) in [
            (b"a'b".as_slice(), "\"a'b\""),
            (b"a\"b".as_slice(), "'a\"b'"),
            (b"a'b\"c".as_slice(), ";a'b\"c\n;"),
            (b"a\nb".as_slice(), ";a\nb\n;"),
        ] {
            assert_eq!(
                scalar_to_cif_value(&JsonArenaValue::String(0..bytes.len()), bytes).unwrap(),
                expected
            );
        }
        assert_eq!(
            scalar_to_cif_value(&JsonArenaValue::Number(2..9), b"xx-0.0E+3yy").unwrap(),
            "-0.0E+3"
        );
        assert_eq!(
            scalar_to_cif_value(&JsonArenaValue::Object(0..0), b"").unwrap_err(),
            MmjsonValueError::UnexpectedType(MmjsonUnexpectedType::Object),
        );
        assert_eq!(
            scalar_to_cif_value(&JsonArenaValue::String(0..1), b"\xff").unwrap_err(),
            MmjsonValueError::InvalidUtf8 { valid_up_to: 0 },
        );
    }
}
