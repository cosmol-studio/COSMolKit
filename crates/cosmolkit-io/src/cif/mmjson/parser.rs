//! Private source-shaped mmJSON lexical error state.

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum JsonErrorKind {
    NoError,
    OutOfMemory,
    UnexpectedEnd,
    MissingRootElement,
    BadRoot,
    ExpectedComma,
    MissingObjectKey,
    ExpectedColon,
    ExpectedEndOfInput,
    UnexpectedComma,
    ExpectedValue,
    ExpectedNull,
    ExpectedFalse,
    ExpectedTrue,
    InvalidNumber,
    MissingExponent,
    IllegalCodepoint,
    InvalidUnicodeEscape,
    UnexpectedEndOfUtf16,
    ExpectedU,
    InvalidUtf16TrailSurrogate,
    UnknownEscape,
    InvalidUtf8,
}

impl JsonErrorKind {
    pub(super) const fn source_text(self) -> &'static str {
        // sajson.h::internal::get_error_text, original switch arm order.
        // Gemmi✔️✔️: case ERROR_NO_ERROR: return "no error";
        // Gemmi✔️✔️: case ERROR_OUT_OF_MEMORY: return  "out of memory";
        // Gemmi✔️✔️: case ERROR_UNEXPECTED_END: return  "unexpected end of input";
        // Gemmi✔️✔️: case ERROR_MISSING_ROOT_ELEMENT: return  "missing root element";
        // Gemmi✔️✔️: case ERROR_BAD_ROOT: return  "document root must be object or array";
        // Gemmi✔️✔️: case ERROR_EXPECTED_COMMA: return  "expected ,";
        // Gemmi✔️✔️: case ERROR_MISSING_OBJECT_KEY: return  "missing object key";
        // Gemmi✔️✔️: case ERROR_EXPECTED_COLON: return  "expected :";
        // Gemmi✔️✔️: case ERROR_EXPECTED_END_OF_INPUT: return  "expected end of input";
        // Gemmi✔️✔️: case ERROR_UNEXPECTED_COMMA: return  "unexpected comma";
        // Gemmi✔️✔️: case ERROR_EXPECTED_VALUE: return  "expected value";
        // Gemmi✔️✔️: case ERROR_EXPECTED_NULL: return  "expected 'null'";
        // Gemmi✔️✔️: case ERROR_EXPECTED_FALSE: return  "expected 'false'";
        // Gemmi✔️✔️: case ERROR_EXPECTED_TRUE: return  "expected 'true'";
        // Gemmi✔️✔️: case ERROR_INVALID_NUMBER: return "invalid number";
        // Gemmi✔️✔️: case ERROR_MISSING_EXPONENT: return  "missing exponent";
        // Gemmi✔️✔️: case ERROR_ILLEGAL_CODEPOINT: return  "illegal unprintable codepoint in string";
        // Gemmi✔️✔️: case ERROR_INVALID_UNICODE_ESCAPE: return  "invalid character in unicode escape";
        // Gemmi✔️✔️: case ERROR_UNEXPECTED_END_OF_UTF16: return  "unexpected end of input during UTF-16 surrogate pair";
        // Gemmi✔️✔️: case ERROR_EXPECTED_U: return  "expected \\u";
        // Gemmi✔️✔️: case ERROR_INVALID_UTF16_TRAIL_SURROGATE: return  "invalid UTF-16 trail surrogate";
        // Gemmi✔️✔️: case ERROR_UNKNOWN_ESCAPE: return  "unknown escape";
        // Gemmi✔️✔️: case ERROR_INVALID_UTF8: return  "invalid UTF-8";
        // Behavior: exhaustive source category/message mapping, with the
        // illegal-codepoint argument rendered by JsonError below.
        // Complexity: O(1) match, identical to the source switch.
        match self {
            Self::NoError => "no error",
            Self::OutOfMemory => "out of memory",
            Self::UnexpectedEnd => "unexpected end of input",
            Self::MissingRootElement => "missing root element",
            Self::BadRoot => "document root must be object or array",
            Self::ExpectedComma => "expected ,",
            Self::MissingObjectKey => "missing object key",
            Self::ExpectedColon => "expected :",
            Self::ExpectedEndOfInput => "expected end of input",
            Self::UnexpectedComma => "unexpected comma",
            Self::ExpectedValue => "expected value",
            Self::ExpectedNull => "expected 'null'",
            Self::ExpectedFalse => "expected 'false'",
            Self::ExpectedTrue => "expected 'true'",
            Self::InvalidNumber => "invalid number",
            Self::MissingExponent => "missing exponent",
            Self::IllegalCodepoint => "illegal unprintable codepoint in string",
            Self::InvalidUnicodeEscape => "invalid character in unicode escape",
            Self::UnexpectedEndOfUtf16 => "unexpected end of input during UTF-16 surrogate pair",
            Self::ExpectedU => "expected \\u",
            Self::InvalidUtf16TrailSurrogate => "invalid UTF-16 trail surrogate",
            Self::UnknownEscape => "unknown escape",
            Self::InvalidUtf8 => "invalid UTF-8",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct JsonError {
    kind: JsonErrorKind,
    byte_offset: usize,
    line: usize,
    column: usize,
    argument: Option<i32>,
}

impl JsonError {
    pub(super) fn at(
        input: &[u8],
        offset: usize,
        kind: JsonErrorKind,
        argument: Option<i32>,
    ) -> Self {
        // sajson.h::parser::make_error:
        // Gemmi✔️✔️: if (!p) {
        // Gemmi✔️✔️:     p = input_end;
        // Gemmi✔️✔️: }
        // Gemmi✔️✔️: error_line = 1;
        // Gemmi✔️✔️: error_column = 1;
        // Gemmi✔️✔️: char* c = input.get_data();
        // Gemmi✔️✔️: while (c < p) {
        // Gemmi✔️✔️:     if (*c == '\r') {
        // Gemmi✔️✔️:         if (c + 1 < p && c[1] == '\n') {
        // Gemmi✔️✔️:             ++error_line;
        // Gemmi✔️✔️:             error_column = 1;
        // Gemmi✔️✔️:             ++c;
        // Gemmi✔️✔️:         } else {
        // Gemmi✔️✔️:             ++error_line;
        // Gemmi✔️✔️:             error_column = 1;
        // Gemmi✔️✔️:         }
        // Gemmi✔️✔️:     } else if (*c == '\n') {
        // Gemmi✔️✔️:         ++error_line;
        // Gemmi✔️✔️:         error_column = 1;
        // Gemmi✔️✔️:     } else {
        // Gemmi✔️✔️:         // TODO: count UTF-8 characters
        // Gemmi✔️✔️:         ++error_column;
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:     ++c;
        // Gemmi✔️✔️: }
        // Gemmi✔️✔️: error_code = code;
        // Gemmi✔️✔️: error_arg = arg;
        // Behavior: clamp the error pointer to the owned buffer; byte-wise
        // columns and CRLF folding reproduce the source for valid offsets.
        // Complexity: one O(offset) scan with no allocation, as in source.
        let offset = offset.min(input.len());
        let mut line = 1;
        let mut column = 1;
        let mut cursor = 0;
        while cursor < offset {
            match input[cursor] {
                b'\r' => {
                    line += 1;
                    column = 1;
                    if cursor + 1 < offset && input[cursor + 1] == b'\n' {
                        cursor += 1;
                    }
                }
                b'\n' => {
                    line += 1;
                    column = 1;
                }
                _ => column += 1,
            }
            cursor += 1;
        }
        Self {
            kind,
            byte_offset: offset,
            line,
            column,
            argument,
        }
    }

    pub(super) const fn kind(&self) -> JsonErrorKind {
        self.kind
    }
    pub(super) const fn byte_offset(&self) -> usize {
        self.byte_offset
    }
    pub(super) const fn line(&self) -> usize {
        self.line
    }
    pub(super) const fn column(&self) -> usize {
        self.column
    }
}

impl fmt::Display for JsonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Gemmi✔️✔️: int written = has_significant_error_arg()
        // Gemmi✔️✔️:     ? SAJSON_snprintf(formatted_error_message, ERROR_BUFFER_LENGTH - 1, "%s: %d", _internal_get_error_text(), error_arg)
        // Gemmi✔️✔️:     : SAJSON_snprintf(formatted_error_message, ERROR_BUFFER_LENGTH - 1, "%s", _internal_get_error_text());
        // Behavior: only illegal-codepoint carries the decimal byte argument.
        // Complexity: bounded formatting and no parser-side retry.
        if self.kind == JsonErrorKind::IllegalCodepoint {
            write!(
                f,
                "{}: {}",
                self.kind.source_text(),
                self.argument.unwrap_or_default()
            )
        } else {
            f.write_str(self.kind.source_text())
        }
    }
}

impl std::error::Error for JsonError {}

pub(super) struct JsonCursor {
    input: Vec<u8>,
    position: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum JsonLiteralKind {
    Null,
    False,
    True,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum JsonNumberKind {
    Double,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum JsonStringFast {
    Plain(std::ops::Range<usize>),
    Slow { start: usize, at: usize },
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum JsonStringProgress {
    Complete(std::ops::Range<usize>),
    Unicode { start: usize, at: usize, end: usize },
    RawUtf8 { start: usize, at: usize, end: usize },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum JsonArenaValue {
    Null,
    Boolean(bool),
    Number(std::ops::Range<usize>),
    String(std::ops::Range<usize>),
    Array(std::ops::Range<usize>),
    Object(std::ops::Range<usize>),
}

#[derive(Default)]
pub(super) struct JsonArena {
    values: Vec<JsonArenaValue>,
    array_items: Vec<usize>,
    object_items: Vec<(std::ops::Range<usize>, usize)>,
}

impl JsonArena {
    fn push_value(&mut self, value: JsonArenaValue) -> usize {
        let index = self.values.len();
        self.values.push(value);
        index
    }

    fn install_array(&mut self, elements: &[usize]) -> Option<usize> {
        // Gemmi❗✔️:         bool install_array(size_t* array_base, size_t* array_end) {
        // Gemmi❗✔️:             using namespace sajson::internal;
        // Gemmi❗✔️:
        // Gemmi❗✔️:             const size_t length = array_end - array_base;
        // Gemmi❗✔️:             bool success;
        // Gemmi❗✔️:             size_t* const new_base = allocator.reserve(length + 1, &success);
        // Gemmi❗✔️:             if (SAJSON_UNLIKELY(!success)) {
        // Gemmi❗✔️:                 return false;
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:             size_t* out = new_base + length + 1;
        // Gemmi❗✔️:             size_t* const structure_end = allocator.get_write_pointer_of(0);
        // Gemmi❗✔️:
        // Gemmi❗✔️:             while (array_end > array_base) {
        // Gemmi❗✔️:                 size_t element = *--array_end;
        // Gemmi❗✔️:                 type element_type = get_element_type(element);
        // Gemmi❗✔️:                 size_t element_value = get_element_value(element);
        // Gemmi❗✔️:                 size_t* element_ptr = structure_end - element_value;
        // Gemmi❗✔️:                 *--out = make_element(element_type, element_ptr - new_base);
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:             *--out = length;
        // Gemmi❗✔️:             return true;
        // Gemmi❗✔️:         }
        // Behavior: child IDs are already installed; append in source order
        // and keep a bounded range. Invalid private IDs return None safely.
        // Complexity: O(length) append, one contiguous range, no sort or
        // recursive copy; equivalent asymptotic container installation.
        if elements.iter().any(|&index| index >= self.values.len()) {
            return None;
        }
        self.array_items.try_reserve(elements.len()).ok()?;
        self.values.try_reserve(1).ok()?;
        let start = self.array_items.len();
        self.array_items.extend_from_slice(elements);
        Some(self.push_value(JsonArenaValue::Array(start..self.array_items.len())))
    }

    fn install_object(
        &mut self,
        input_len: usize,
        entries: &[(std::ops::Range<usize>, usize)],
    ) -> Option<usize> {
        // Gemmi❗✔️:         bool install_object(size_t* object_base, size_t* object_end) {
        // Gemmi❗✔️:             using namespace internal;
        // Gemmi❗✔️:
        // Gemmi❗✔️:             assert((object_end - object_base) % 3 == 0);
        // Gemmi❗✔️:             const size_t length_times_3 = object_end - object_base;
        // Gemmi❗✔️: #ifndef SAJSON_UNSORTED_OBJECT_KEYS
        // Gemmi❗✔️:             std::sort(
        // Gemmi❗✔️:                 reinterpret_cast<object_key_record*>(object_base),
        // Gemmi❗✔️:                 reinterpret_cast<object_key_record*>(object_end),
        // Gemmi❗✔️:                 object_key_comparator(input.get_data()));
        // Gemmi❗✔️: #endif
        // Gemmi❗✔️:
        // Gemmi❗✔️:             bool success;
        // Gemmi❗✔️:             size_t* const new_base = allocator.reserve(length_times_3 + 1, &success);
        // Gemmi❗✔️:             if (SAJSON_UNLIKELY(!success)) {
        // Gemmi❗✔️:                 return false;
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:             size_t* out = new_base + length_times_3 + 1;
        // Gemmi❗✔️:             size_t* const structure_end = allocator.get_write_pointer_of(0);
        // Gemmi❗✔️:
        // Gemmi❗✔️:             while (object_end > object_base) {
        // Gemmi❗✔️:                 size_t element = *--object_end;
        // Gemmi❗✔️:                 type element_type = get_element_type(element);
        // Gemmi❗✔️:                 size_t element_value = get_element_value(element);
        // Gemmi❗✔️:                 size_t* element_ptr = structure_end - element_value;
        // Gemmi❗✔️:
        // Gemmi❗✔️:                 *--out = make_element(element_type, element_ptr - new_base);
        // Gemmi❗✔️:                 *--out = *--object_end;
        // Gemmi❗✔️:                 *--out = *--object_end;
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:             *--out = length_times_3 / 3;
        // Gemmi❗✔️:             return true;
        // Gemmi❗✔️:         }
        // Behavior: UNSORTED_OBJECT_KEYS disables source sort; duplicate keys
        // remain separate ordered entries. Invalid private spans return None.
        // Complexity: O(length) contiguous append; no map or key sorting.
        if entries.iter().any(|(key, value)| {
            key.start > key.end || key.end > input_len || *value >= self.values.len()
        }) {
            return None;
        }
        self.object_items.try_reserve(entries.len()).ok()?;
        self.values.try_reserve(1).ok()?;
        let start = self.object_items.len();
        self.object_items.extend_from_slice(entries);
        Some(self.push_value(JsonArenaValue::Object(start..self.object_items.len())))
    }

    pub(super) fn value(&self, index: usize) -> Option<&JsonArenaValue> {
        self.values.get(index)
    }
    pub(super) fn array_element(&self, index: usize, position: usize) -> Option<&JsonArenaValue> {
        let JsonArenaValue::Array(range) = self.value(index)? else {
            return None;
        };
        self.value(*self.array_items.get(range.start.checked_add(position)?)?)
            .filter(|_| position < range.len())
    }
    pub(super) fn array_element_by_range(
        &self,
        range: &std::ops::Range<usize>,
        position: usize,
    ) -> Option<&JsonArenaValue> {
        if position >= range.len() {
            return None;
        }
        self.value(*self.array_items.get(range.start.checked_add(position)?)?)
    }
    pub(super) fn object_entry(
        &self,
        index: usize,
        position: usize,
    ) -> Option<(&std::ops::Range<usize>, &JsonArenaValue)> {
        let JsonArenaValue::Object(range) = self.value(index)? else {
            return None;
        };
        if position >= range.len() {
            return None;
        }
        let (key, value) = self.object_items.get(range.start.checked_add(position)?)?;
        Some((key, self.value(*value)?))
    }

    pub(super) fn object_entry_id(
        &self,
        index: usize,
        position: usize,
    ) -> Option<(&std::ops::Range<usize>, usize)> {
        let JsonArenaValue::Object(range) = self.value(index)? else {
            return None;
        };
        if position >= range.len() {
            return None;
        }
        let (key, value) = self.object_items.get(range.start.checked_add(position)?)?;
        Some((key, *value))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum JsonStackKind {
    Array,
    Object,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum JsonContainerTransition {
    Value,
    ObjectKey,
    Closed(usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum JsonStackEntry {
    Key(std::ops::Range<usize>),
    Value(usize),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct JsonStackFrame {
    kind: JsonStackKind,
    start: usize,
}

#[derive(Default)]
struct JsonParseStack {
    entries: Vec<JsonStackEntry>,
    frames: Vec<JsonStackFrame>,
}

impl JsonParseStack {
    fn open(&mut self, kind: JsonStackKind) {
        self.frames.push(JsonStackFrame {
            kind,
            start: self.entries.len(),
        });
    }

    fn push_key(&mut self, key: std::ops::Range<usize>) -> Option<()> {
        if self.frames.last()?.kind != JsonStackKind::Object {
            return None;
        }
        self.entries.push(JsonStackEntry::Key(key));
        Some(())
    }

    fn push_value(&mut self, value: usize) -> Option<()> {
        self.frames.last()?;
        self.entries.push(JsonStackEntry::Value(value));
        Some(())
    }

    fn finish(&mut self, arena: &mut JsonArena, input_len: usize) -> Option<usize> {
        // Gemmi❗❌:             pop_object: {
        // Gemmi❗❌:                 ++p;
        // Gemmi❗❌:                 size_t* base_ptr = stack.get_pointer_from_offset(current_base);
        // Gemmi❗❌:                 pop_element = *base_ptr;
        // Gemmi❗❌:                 if (SAJSON_UNLIKELY(!install_object(base_ptr + 1, stack.get_top()))) {
        // Gemmi❗❌:                     return oom(p);
        // Gemmi❗❌:                 }
        // Gemmi❗❌:                 goto pop;
        // Gemmi❗❌:             }
        // Gemmi❗❌:             pop_array: {
        // Gemmi❗❌:                 ++p;
        // Gemmi❗❌:                 size_t* base_ptr = stack.get_pointer_from_offset(current_base);
        // Gemmi❗❌:                 pop_element = *base_ptr;
        // Gemmi❗❌:                 if (SAJSON_UNLIKELY(!install_array(base_ptr + 1, stack.get_top()))) {
        // Gemmi❗❌:                     return oom(p);
        // Gemmi❗❌:                 }
        // Gemmi❗❌:                 goto pop;
        // Gemmi❗❌:             }
        // Behavior: one explicit stack carries nested container boundaries.
        // An invalid private sequence leaves the frame/entries intact.
        // Complexity: O(children), but a temporary Vec is allocated per
        // container in addition to the arena append; source uses one reserve.
        let frame = *self.frames.last()?;
        let pending = self.entries.get(frame.start..)?;
        let installed = match frame.kind {
            JsonStackKind::Array => {
                let values: Option<Vec<_>> = pending
                    .iter()
                    .map(|entry| match entry {
                        JsonStackEntry::Value(value) => Some(*value),
                        JsonStackEntry::Key(_) => None,
                    })
                    .collect();
                arena.install_array(&values?)?
            }
            JsonStackKind::Object => {
                let mut chunks = pending.chunks_exact(2);
                let pairs: Option<Vec<_>> = chunks
                    .by_ref()
                    .map(|chunk| match chunk {
                        [JsonStackEntry::Key(key), JsonStackEntry::Value(value)] => {
                            Some((key.clone(), *value))
                        }
                        _ => None,
                    })
                    .collect();
                if !chunks.remainder().is_empty() {
                    return None;
                }
                arena.install_object(input_len, &pairs?)?
            }
        };
        self.frames.pop();
        self.entries.truncate(frame.start);
        Some(installed)
    }
}

impl JsonCursor {
    pub(super) fn parse_document(&mut self) -> Result<(JsonArena, usize), JsonError> {
        // Gemmi❗❌:         bool parse() {
        // Gemmi❗❌:             using namespace internal;
        // Gemmi❗❌:
        // Gemmi❗❌:             // p points to the character currently being parsed
        // Gemmi❗❌:             char* p = input.get_data();
        // Gemmi❗❌:
        // Gemmi❗❌:             bool success;
        // Gemmi❗❌:             auto stack = allocator.get_stack_head(&success);
        // Gemmi❗❌:             if (SAJSON_UNLIKELY(!success)) {
        // Gemmi❗❌:                 return oom(p);
        // Gemmi❗❌:             }
        // Gemmi❗❌:
        // Gemmi❗❌:             p = skip_whitespace(p);
        // Gemmi❗❌:             if (SAJSON_UNLIKELY(!p)) {
        // Gemmi❗❌:                 return make_error(p, ERROR_MISSING_ROOT_ELEMENT);
        // Gemmi❗❌:             }
        // Gemmi❗❌:
        // Gemmi❗❌:             // current_base is an offset to the first element of the current structure (object or array)
        // Gemmi❗❌:             size_t current_base = stack.get_size();
        // Gemmi❗❌:             type current_structure_type;
        // Gemmi❗❌:             if (*p == '[') {
        // Gemmi❗❌:                 current_structure_type = TYPE_ARRAY;
        // Gemmi❗❌:                 bool s = stack.push(make_element(current_structure_type, ROOT_MARKER));
        // Gemmi❗❌:                 if (SAJSON_UNLIKELY(!s)) {
        // Gemmi❗❌:                     return oom(p);
        // Gemmi❗❌:                 }
        // Gemmi❗❌:                 goto array_close_or_element;
        // Gemmi❗❌:             } else if (*p == '{') {
        // Gemmi❗❌:                 current_structure_type = TYPE_OBJECT;
        // Gemmi❗❌:                 bool s = stack.push(make_element(current_structure_type, ROOT_MARKER));
        // Gemmi❗❌:                 if (SAJSON_UNLIKELY(!s)) {
        // Gemmi❗❌:                     return oom(p);
        // Gemmi❗❌:                 }
        // Gemmi❗❌:                 goto object_close_or_element;
        // Gemmi❗❌:             } else {
        // Gemmi❗❌:                 return make_error(p, ERROR_BAD_ROOT);
        // Gemmi❗❌:             }
        // Gemmi❗❌:             // BEGIN STATE MACHINE
        // Gemmi❗❌:             size_t pop_element; // used as an argument into the `pop` routine
        // Gemmi❗❌:             if (0) { // purely for structure
        // The first-member, close/comma and key labels are copied inside the
        // actual implementing helpers above; this owner drives their order.
        // Gemmi❗❌:             next_element:
        // Gemmi❗❌:                 p = skip_whitespace(p);
        // Gemmi❗❌:                 if (SAJSON_UNLIKELY(!p)) {
        // Gemmi❗❌:                     return unexpected_end();
        // Gemmi❗❌:                 }
        // Gemmi❗❌:                 type value_type_result;
        // Gemmi❗❌:                 switch (*p) {
        // Gemmi❗❌:                     case 0:
        // Gemmi❗❌:                         return unexpected_end(p);
        // Gemmi❗❌:                     case 'n':
        // Gemmi❗❌:                         p = parse_null(p);
        // Gemmi❗❌:                         if (!p) {
        // Gemmi❗❌:                             return false;
        // Gemmi❗❌:                         }
        // Gemmi❗❌:                         value_type_result = TYPE_NULL;
        // Gemmi❗❌:                         break;
        // Gemmi❗❌:                     case 'f':
        // Gemmi❗❌:                         p = parse_false(p);
        // Gemmi❗❌:                         if (!p) {
        // Gemmi❗❌:                             return false;
        // Gemmi❗❌:                         }
        // Gemmi❗❌:                         value_type_result = TYPE_FALSE;
        // Gemmi❗❌:                         break;
        // Gemmi❗❌:                     case 't':
        // Gemmi❗❌:                         p = parse_true(p);
        // Gemmi❗❌:                         if (!p) {
        // Gemmi❗❌:                             return false;
        // Gemmi❗❌:                         }
        // Gemmi❗❌:                         value_type_result = TYPE_TRUE;
        // Gemmi❗❌:                         break;
        // Gemmi❗❌:                     case '0':
        // Gemmi❗❌:                     case '1':
        // Gemmi❗❌:                     case '2':
        // Gemmi❗❌:                     case '3':
        // Gemmi❗❌:                     case '4':
        // Gemmi❗❌:                     case '5':
        // Gemmi❗❌:                     case '6':
        // Gemmi❗❌:                     case '7':
        // Gemmi❗❌:                     case '8':
        // Gemmi❗❌:                     case '9':
        // Gemmi❗❌:                     case '-': {
        // Gemmi❗❌:                         auto result = parse_number(p);
        // Gemmi❗❌:                         p = result.first;
        // Gemmi❗❌:                         if (!p) {
        // Gemmi❗❌:                             return false;
        // Gemmi❗❌:                         }
        // Gemmi❗❌:                         value_type_result = result.second;
        // Gemmi❗❌:                         break;
        // Gemmi❗❌:                     }
        // Gemmi❗❌:                     case '"': {
        // Gemmi❗❌:                         bool success_;
        // Gemmi❗❌:                         size_t* string_tag = allocator.reserve(2, &success_);
        // Gemmi❗❌:                         if (SAJSON_UNLIKELY(!success_)) {
        // Gemmi❗❌:                             return oom(p);
        // Gemmi❗❌:                         }
        // Gemmi❗❌:                         p = parse_string(p, string_tag);
        // Gemmi❗❌:                         if (!p) {
        // Gemmi❗❌:                             return false;
        // Gemmi❗❌:                         }
        // Gemmi❗❌:                         value_type_result = TYPE_STRING;
        // Gemmi❗❌:                         break;
        // Gemmi❗❌:                     }
        // Gemmi❗❌:                     case '[': {
        // Gemmi❗❌:                         size_t previous_base = current_base;
        // Gemmi❗❌:                         current_base = stack.get_size();
        // Gemmi❗❌:                         bool s = stack.push(make_element(current_structure_type, previous_base));
        // Gemmi❗❌:                         if (SAJSON_UNLIKELY(!s)) {
        // Gemmi❗❌:                             return oom(p);
        // Gemmi❗❌:                         }
        // Gemmi❗❌:                         current_structure_type = TYPE_ARRAY;
        // Gemmi❗❌:                         goto array_close_or_element;
        // Gemmi❗❌:                     }
        // Gemmi❗❌:                     case '{': {
        // Gemmi❗❌:                         size_t previous_base = current_base;
        // Gemmi❗❌:                         current_base = stack.get_size();
        // Gemmi❗❌:                         bool s = stack.push(make_element(current_structure_type, previous_base));
        // Gemmi❗❌:                         if (SAJSON_UNLIKELY(!s)) {
        // Gemmi❗❌:                             return oom(p);
        // Gemmi❗❌:                         }
        // Gemmi❗❌:                         current_structure_type = TYPE_OBJECT;
        // Gemmi❗❌:                         goto object_close_or_element;
        // Gemmi❗❌:                     }
        // Gemmi❗❌:                     pop: {
        // Gemmi❗❌:                         size_t parent = get_element_value(pop_element);
        // Gemmi❗❌:                         if (parent == ROOT_MARKER) {
        // Gemmi❗❌:                             root_type = current_structure_type;
        // Gemmi❗❌:                             p = skip_whitespace(p);
        // Gemmi❗❌:                             if (SAJSON_UNLIKELY(p)) {
        // Gemmi❗❌:                                 return make_error(p, ERROR_EXPECTED_END_OF_INPUT);
        // Gemmi❗❌:                             }
        // Gemmi❗❌:                             return true;
        // Gemmi❗❌:                         }
        // Gemmi❗❌:                         stack.reset(current_base);
        // Gemmi❗❌:                         current_base = parent;
        // Gemmi❗❌:                         value_type_result = current_structure_type;
        // Gemmi❗❌:                         current_structure_type = get_element_type(pop_element);
        // Gemmi❗❌:                         break;
        // Gemmi❗❌:                     }
        // Gemmi❗❌:                     case ',':
        // Gemmi❗❌:                         return make_error(p, ERROR_UNEXPECTED_COMMA);
        // Gemmi❗❌:                     default:
        // Gemmi❗❌:                         return make_error(p, ERROR_EXPECTED_VALUE);
        // Gemmi❗❌:                 }
        // Gemmi❗❌:                 bool s = stack.push(make_element(
        // Gemmi❗❌:                     value_type_result,
        // Gemmi❗❌:                     allocator.get_write_offset()));
        // Gemmi❗❌:                 if (SAJSON_UNLIKELY(!s)) {
        // Gemmi❗❌:                     return oom(p);
        // Gemmi❗❌:                 }
        // Gemmi❗❌:                 goto structure_close_or_comma;
        // Gemmi❗❌:             }
        // Gemmi❗❌:             SAJSON_UNREACHABLE();
        // Gemmi❗❌:         }
        // Behavior: only array/object roots; explicit stack preserves source
        // order and duplicate keys. Raw number spans are never converted.
        // Complexity: O(tokens + bytes), but the J10 stack currently makes
        // an extra temporary Vec per container; this is an acknowledged gap.
        self.position = 0;
        let Some(position) = self.skip_whitespace() else {
            return Err(JsonError::at(
                &self.input,
                self.position,
                JsonErrorKind::MissingRootElement,
                None,
            ));
        };
        let root_kind = match self.input[position] {
            b'[' => JsonStackKind::Array,
            b'{' => JsonStackKind::Object,
            _ => {
                return Err(JsonError::at(
                    &self.input,
                    position,
                    JsonErrorKind::BadRoot,
                    None,
                ));
            }
        };
        let mut arena = JsonArena::default();
        let mut stack = JsonParseStack::default();
        stack.open(root_kind);
        let mut next = self.container_transition(&mut stack, &mut arena, true)?;
        loop {
            match next {
                JsonContainerTransition::Closed(value) => {
                    if stack.frames.is_empty() {
                        if let Some(position) = self.skip_whitespace() {
                            return Err(JsonError::at(
                                &self.input,
                                position,
                                JsonErrorKind::ExpectedEndOfInput,
                                None,
                            ));
                        }
                        return Ok((arena, value));
                    }
                    stack.push_value(value).ok_or_else(|| {
                        JsonError::at(&self.input, self.position, JsonErrorKind::OutOfMemory, None)
                    })?;
                    next = self.container_transition(&mut stack, &mut arena, false)?;
                }
                JsonContainerTransition::ObjectKey => {
                    self.parse_object_key_colon(&mut stack)?;
                    next = JsonContainerTransition::Value;
                }
                JsonContainerTransition::Value => {
                    let Some(position) = self.skip_whitespace() else {
                        return Err(JsonError::at(
                            &self.input,
                            self.position,
                            JsonErrorKind::UnexpectedEnd,
                            None,
                        ));
                    };
                    let scalar = match self.input[position] {
                        0 => {
                            return Err(JsonError::at(
                                &self.input,
                                position,
                                JsonErrorKind::UnexpectedEnd,
                                None,
                            ));
                        }
                        b'n' => {
                            self.parse_null()?;
                            JsonArenaValue::Null
                        }
                        b'f' => {
                            self.parse_false()?;
                            JsonArenaValue::Boolean(false)
                        }
                        b't' => {
                            self.parse_true()?;
                            JsonArenaValue::Boolean(true)
                        }
                        b'0'..=b'9' | b'-' => {
                            let (span, _) = self.parse_number_raw()?;
                            JsonArenaValue::Number(span)
                        }
                        b'"' => {
                            let span = match self.parse_string_fast()? {
                                JsonStringFast::Plain(span) => span,
                                JsonStringFast::Slow { start, at } => {
                                    match self.parse_string_simple(start, at)? {
                                        JsonStringProgress::Complete(span) => span,
                                        JsonStringProgress::Unicode { .. }
                                        | JsonStringProgress::RawUtf8 { .. } => {
                                            return Err(JsonError::at(
                                                &self.input,
                                                self.position,
                                                JsonErrorKind::UnexpectedEnd,
                                                None,
                                            ));
                                        }
                                    }
                                }
                            };
                            JsonArenaValue::String(span)
                        }
                        b'[' | b'{' => {
                            let kind = if self.input[position] == b'[' {
                                JsonStackKind::Array
                            } else {
                                JsonStackKind::Object
                            };
                            stack.open(kind);
                            next = self.container_transition(&mut stack, &mut arena, true)?;
                            continue;
                        }
                        b',' => {
                            return Err(JsonError::at(
                                &self.input,
                                position,
                                JsonErrorKind::UnexpectedComma,
                                None,
                            ));
                        }
                        _ => {
                            return Err(JsonError::at(
                                &self.input,
                                position,
                                JsonErrorKind::ExpectedValue,
                                None,
                            ));
                        }
                    };
                    arena.values.try_reserve(1).map_err(|_| {
                        JsonError::at(&self.input, self.position, JsonErrorKind::OutOfMemory, None)
                    })?;
                    let value = arena.push_value(scalar);
                    stack.push_value(value).ok_or_else(|| {
                        JsonError::at(&self.input, self.position, JsonErrorKind::OutOfMemory, None)
                    })?;
                    next = self.container_transition(&mut stack, &mut arena, false)?;
                }
            }
        }
    }

    fn container_transition(
        &mut self,
        stack: &mut JsonParseStack,
        arena: &mut JsonArena,
        opened: bool,
    ) -> Result<JsonContainerTransition, JsonError> {
        // Gemmi❗✔️:             array_close_or_element:
        // Gemmi❗✔️:                 p = skip_whitespace(p + 1);
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(!p)) {
        // Gemmi❗✔️:                     return unexpected_end();
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 if (*p == ']') {
        // Gemmi❗✔️:                     goto pop_array;
        // Gemmi❗✔️:                 } else {
        // Gemmi❗✔️:                     goto next_element;
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:             object_close_or_element:
        // Gemmi❗✔️:                 p = skip_whitespace(p + 1);
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(!p)) {
        // Gemmi❗✔️:                     return unexpected_end();
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 if (*p == '}') {
        // Gemmi❗✔️:                     goto pop_object;
        // Gemmi❗✔️:                 } else {
        // Gemmi❗✔️:                     goto object_key;
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:             structure_close_or_comma:
        // Gemmi❗✔️:                 p = skip_whitespace(p);
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(!p)) {
        // Gemmi❗✔️:                     return unexpected_end();
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 if (current_structure_type == TYPE_ARRAY) {
        // Gemmi❗✔️:                     if (*p == ']') {
        // Gemmi❗✔️:                         goto pop_array;
        // Gemmi❗✔️:                     } else {
        // Gemmi❗✔️:                         if (SAJSON_UNLIKELY(*p != ',')) {
        // Gemmi❗✔️:                             return make_error(p, ERROR_EXPECTED_COMMA);
        // Gemmi❗✔️:                         }
        // Gemmi❗✔️:                         ++p;
        // Gemmi❗✔️:                         goto next_element;
        // Gemmi❗✔️:                     }
        // Gemmi❗✔️:                 } else {
        // Gemmi❗✔️:                     assert(current_structure_type == TYPE_OBJECT);
        // Gemmi❗✔️:                     if (*p == '}') {
        // Gemmi❗✔️:                         goto pop_object;
        // Gemmi❗✔️:                     } else {
        // Gemmi❗✔️:                         if (SAJSON_UNLIKELY(*p != ',')) {
        // Gemmi❗✔️:                             return make_error(p, ERROR_EXPECTED_COMMA);
        // Gemmi❗✔️:                         }
        // Gemmi❗✔️:                         ++p;
        // Gemmi❗✔️:                         goto object_key;
        // Gemmi❗✔️:                     }
        // Gemmi❗✔️:                 }
        // Behavior: opened distinguishes the two source first-member labels
        // from close-or-comma after a value. The caller next dispatches the
        // returned label, retaining trailing-comma error precedence.
        // Complexity: O(whitespace) scan, O(children) only when closing via
        // the shared stack/arena installer; no duplicate parser or recursion.
        let kind = stack.frames.last().map(|frame| frame.kind).ok_or_else(|| {
            JsonError::at(
                &self.input,
                self.position,
                JsonErrorKind::ExpectedValue,
                None,
            )
        })?;
        if opened {
            self.position += 1;
        }
        let Some(position) = self.skip_whitespace() else {
            return Err(JsonError::at(
                &self.input,
                self.position,
                JsonErrorKind::UnexpectedEnd,
                None,
            ));
        };
        let closing = match kind {
            JsonStackKind::Array => b']',
            JsonStackKind::Object => b'}',
        };
        if self.input[position] == closing {
            self.position += 1;
            let value = stack.finish(arena, self.input.len()).ok_or_else(|| {
                JsonError::at(&self.input, self.position, JsonErrorKind::OutOfMemory, None)
            })?;
            return Ok(JsonContainerTransition::Closed(value));
        }
        if !opened {
            if self.input[position] != b',' {
                return Err(JsonError::at(
                    &self.input,
                    position,
                    JsonErrorKind::ExpectedComma,
                    None,
                ));
            }
            self.position += 1;
            if kind == JsonStackKind::Array {
                // Gemmi❗✔️:             next_element:
                // Gemmi❗✔️:                 p = skip_whitespace(p);
                // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(!p)) {
                // Gemmi❗✔️:                     return unexpected_end();
                // Gemmi❗✔️:                 }
                // Gemmi❗✔️:                     case ',':
                // Gemmi❗✔️:                         return make_error(p, ERROR_UNEXPECTED_COMMA);
                // Gemmi❗✔️:                     default:
                // Gemmi❗✔️:                         return make_error(p, ERROR_EXPECTED_VALUE);
                // Behavior: precheck only the delimiter/error cases reached
                // after a comma; scalar dispatch remains J13's owner.
                let Some(next) = self.skip_whitespace() else {
                    return Err(JsonError::at(
                        &self.input,
                        self.position,
                        JsonErrorKind::UnexpectedEnd,
                        None,
                    ));
                };
                match self.input[next] {
                    0 => {
                        return Err(JsonError::at(
                            &self.input,
                            next,
                            JsonErrorKind::UnexpectedEnd,
                            None,
                        ));
                    }
                    b',' => {
                        return Err(JsonError::at(
                            &self.input,
                            next,
                            JsonErrorKind::UnexpectedComma,
                            None,
                        ));
                    }
                    b']' | b'}' => {
                        return Err(JsonError::at(
                            &self.input,
                            next,
                            JsonErrorKind::ExpectedValue,
                            None,
                        ));
                    }
                    _ => {}
                }
            }
        }
        Ok(match kind {
            JsonStackKind::Array => JsonContainerTransition::Value,
            JsonStackKind::Object => JsonContainerTransition::ObjectKey,
        })
    }

    fn parse_object_key_colon(&mut self, stack: &mut JsonParseStack) -> Result<(), JsonError> {
        // Gemmi❗✔️:             object_key: {
        // Gemmi❗✔️:                 p = skip_whitespace(p);
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(!p)) {
        // Gemmi❗✔️:                     return unexpected_end();
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(*p != '"')) {
        // Gemmi❗✔️:                     return make_error(p, ERROR_MISSING_OBJECT_KEY);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 bool success_;
        // Gemmi❗✔️:                 size_t* out = stack.reserve(2, &success_);
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(!success_)) {
        // Gemmi❗✔️:                     return oom(p);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 p = parse_string(p, out);
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(!p)) {
        // Gemmi❗✔️:                     return false;
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 p = skip_whitespace(p);
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(!p || *p != ':')) {
        // Gemmi❗✔️:                     return make_error(p, ERROR_EXPECTED_COLON);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 ++p;
        // Gemmi❗✔️:                 goto next_element;
        // Gemmi❗✔️:             }
        // Behavior: existing J06–J09 string owner decodes the key. The
        // explicit stack receives its byte span only after colon validation;
        // on parse error the private stack is not externally observable.
        // Complexity: whitespace/key scan is linear, one stack entry, no
        // key copy or sorted-map lookup.
        let Some(position) = self.skip_whitespace() else {
            return Err(JsonError::at(
                &self.input,
                self.position,
                JsonErrorKind::UnexpectedEnd,
                None,
            ));
        };
        if self.input[position] != b'"' {
            return Err(JsonError::at(
                &self.input,
                position,
                JsonErrorKind::MissingObjectKey,
                None,
            ));
        }
        let span = match self.parse_string_fast()? {
            JsonStringFast::Plain(span) => span,
            JsonStringFast::Slow { start, at } => match self.parse_string_simple(start, at)? {
                JsonStringProgress::Complete(span) => span,
                JsonStringProgress::Unicode { .. } | JsonStringProgress::RawUtf8 { .. } => {
                    return Err(JsonError::at(
                        &self.input,
                        self.position,
                        JsonErrorKind::UnexpectedEnd,
                        None,
                    ));
                }
            },
        };
        let Some(position) = self.skip_whitespace() else {
            return Err(JsonError::at(
                &self.input,
                self.position,
                JsonErrorKind::ExpectedColon,
                None,
            ));
        };
        if self.input[position] != b':' {
            return Err(JsonError::at(
                &self.input,
                position,
                JsonErrorKind::ExpectedColon,
                None,
            ));
        }
        stack.push_key(span).ok_or_else(|| {
            JsonError::at(&self.input, position, JsonErrorKind::ExpectedValue, None)
        })?;
        self.position += 1;
        Ok(())
    }

    pub(super) fn new(input: &[u8]) -> Self {
        Self {
            input: input.to_vec(),
            position: 0,
        }
    }

    pub(super) fn bytes(&self) -> &[u8] {
        &self.input
    }

    fn at_eof(&self) -> bool {
        // Gemmi✔️✔️: bool at_eof(const char* p) {
        // Gemmi✔️✔️:     return p == input_end;
        // Gemmi✔️✔️: }
        // Behavior: position is always bounded by input.len(); equality is
        // therefore the source end-pointer check. Complexity: O(1).
        self.position == self.input.len()
    }

    fn has_remaining_characters(&self, remaining: usize) -> bool {
        // Gemmi✔️✔️: bool has_remaining_characters(char* p, ptrdiff_t remaining) {
        // Gemmi✔️✔️:     return input_end - p >= remaining;
        // Gemmi✔️✔️: }
        // Behavior: all active callers pass nonnegative literal widths;
        // subtraction is safe because the cursor invariant keeps position
        // within the input. Complexity: O(1).
        self.input.len() - self.position >= remaining
    }

    fn skip_whitespace(&mut self) -> Option<usize> {
        // Gemmi✔️✔️: char* skip_whitespace(char* p) {
        // Gemmi✔️✔️:     // There is an opportunity to make better use of superscalar
        // Gemmi✔️✔️:     // hardware here* but if someone cares about JSON parsing
        // Gemmi✔️✔️:     // performance the first thing they do is minify, so prefer
        // Gemmi✔️✔️:     // to optimize for code size here.
        // Gemmi✔️✔️:     // * https://github.com/chadaustin/Web-Benchmarks/blob/master/json/third-party/pjson/pjson.h#L1873
        // Gemmi✔️✔️:     for (;;) {
        // Gemmi✔️✔️:         if (SAJSON_UNLIKELY(p == input_end)) {
        // Gemmi✔️✔️:             return 0;
        // Gemmi✔️✔️:         } else if (internal::is_whitespace(*p)) {
        // Gemmi✔️✔️:             ++p;
        // Gemmi✔️✔️:         } else {
        // Gemmi✔️✔️:             return p;
        // Gemmi✔️✔️:         }
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️: }
        // Gemmi✔️✔️: inline bool is_whitespace(char c) {
        // Gemmi✔️✔️:     //return c == '\r' || c == '\n' || c == '\t' || c == ' ';
        // Gemmi✔️✔️:     return (globals::parse_flags[static_cast<unsigned char>(c)] & 2) != 0;
        // Gemmi✔️✔️: }
        // Behavior: exactly the four source JSON whitespace bytes; the
        // source null pointer on end maps to None and leaves position at end.
        // Complexity: one linear scan of skipped bytes, no allocation.
        loop {
            let Some(&byte) = self.input.get(self.position) else {
                return None;
            };
            if matches!(byte, b' ' | b'\t' | b'\r' | b'\n') {
                self.position += 1;
            } else {
                return Some(self.position);
            }
        }
    }

    fn parse_null(&mut self) -> Result<(usize, JsonLiteralKind), JsonError> {
        // Gemmi✔️✔️: char* parse_null(char* p) {
        // Gemmi✔️✔️:     if (SAJSON_UNLIKELY(!has_remaining_characters(p, 4))) {
        // Gemmi✔️✔️:         make_error(p, ERROR_UNEXPECTED_END);
        // Gemmi✔️✔️:         return 0;
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:     char p1 = p[1];
        // Gemmi✔️✔️:     char p2 = p[2];
        // Gemmi✔️✔️:     char p3 = p[3];
        // Gemmi✔️✔️:     if (SAJSON_UNLIKELY(p1 != 'u' || p2 != 'l' || p3 != 'l')) {
        // Gemmi✔️✔️:         make_error(p, ERROR_EXPECTED_NULL);
        // Gemmi✔️✔️:         return 0;
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:     return p + 4;
        // Gemmi✔️✔️: }
        // Behavior: first byte belongs to caller dispatch, as in source;
        // truncated input takes UnexpectedEnd before spelling check.
        // Complexity: O(1), no allocation on success.
        let start = self.position;
        if !self.has_remaining_characters(4) {
            return Err(JsonError::at(
                &self.input,
                start,
                JsonErrorKind::UnexpectedEnd,
                None,
            ));
        }
        if self.input[start + 1..start + 4] != *b"ull" {
            return Err(JsonError::at(
                &self.input,
                start,
                JsonErrorKind::ExpectedNull,
                None,
            ));
        }
        self.position += 4;
        Ok((self.position, JsonLiteralKind::Null))
    }

    fn parse_false(&mut self) -> Result<(usize, JsonLiteralKind), JsonError> {
        // Gemmi✔️✔️: char* parse_false(char* p) {
        // Gemmi✔️✔️:     if (SAJSON_UNLIKELY(!has_remaining_characters(p, 5))) {
        // Gemmi✔️✔️:         return make_error(p, ERROR_UNEXPECTED_END);
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:     char p1 = p[1];
        // Gemmi✔️✔️:     char p2 = p[2];
        // Gemmi✔️✔️:     char p3 = p[3];
        // Gemmi✔️✔️:     char p4 = p[4];
        // Gemmi✔️✔️:     if (SAJSON_UNLIKELY(p1 != 'a' || p2 != 'l' || p3 != 's' || p4 != 'e')) {
        // Gemmi✔️✔️:         return make_error(p, ERROR_EXPECTED_FALSE);
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:     return p + 5;
        // Gemmi✔️✔️: }
        // Behavior: caller checks the initial f; error offset remains p.
        // Complexity: O(1), no allocation on success.
        let start = self.position;
        if !self.has_remaining_characters(5) {
            return Err(JsonError::at(
                &self.input,
                start,
                JsonErrorKind::UnexpectedEnd,
                None,
            ));
        }
        if self.input[start + 1..start + 5] != *b"alse" {
            return Err(JsonError::at(
                &self.input,
                start,
                JsonErrorKind::ExpectedFalse,
                None,
            ));
        }
        self.position += 5;
        Ok((self.position, JsonLiteralKind::False))
    }

    fn parse_true(&mut self) -> Result<(usize, JsonLiteralKind), JsonError> {
        // Gemmi✔️✔️: char* parse_true(char* p) {
        // Gemmi✔️✔️:     if (SAJSON_UNLIKELY(!has_remaining_characters(p, 4))) {
        // Gemmi✔️✔️:         return make_error(p, ERROR_UNEXPECTED_END);
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:     char p1 = p[1];
        // Gemmi✔️✔️:     char p2 = p[2];
        // Gemmi✔️✔️:     char p3 = p[3];
        // Gemmi✔️✔️:     if (SAJSON_UNLIKELY(p1 != 'r' || p2 != 'u' || p3 != 'e')) {
        // Gemmi✔️✔️:         return make_error(p, ERROR_EXPECTED_TRUE);
        // Gemmi✔️✔️:     }
        // Gemmi✔️✔️:     return p + 4;
        // Gemmi✔️✔️: }
        // Behavior: caller checks the initial t; length precedes spelling.
        // Complexity: O(1), no allocation on success.
        let start = self.position;
        if !self.has_remaining_characters(4) {
            return Err(JsonError::at(
                &self.input,
                start,
                JsonErrorKind::UnexpectedEnd,
                None,
            ));
        }
        if self.input[start + 1..start + 4] != *b"rue" {
            return Err(JsonError::at(
                &self.input,
                start,
                JsonErrorKind::ExpectedTrue,
                None,
            ));
        }
        self.position += 4;
        Ok((self.position, JsonLiteralKind::True))
    }

    fn parse_number_raw(&mut self) -> Result<(std::ops::Range<usize>, JsonNumberKind), JsonError> {
        // Gemmi❗✔️:         std::pair<char*, type> parse_number(char* p) {
        // Gemmi❗✔️:             size_t start = p - input.get_data();
        // Gemmi❗✔️:             if ('-' == *p) {
        // Gemmi❗✔️:                 ++p;
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(at_eof(p))) {
        // Gemmi❗✔️:                     return std::make_pair(make_error(p, ERROR_UNEXPECTED_END), TYPE_NULL);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:             if (*p == '0') {
        // Gemmi❗✔️:                 ++p;
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(at_eof(p))) {
        // Gemmi❗✔️:                     return std::make_pair(make_error(p, ERROR_UNEXPECTED_END), TYPE_NULL);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:             } else {
        // Gemmi❗✔️:                 unsigned char c = *p;
        // Gemmi❗✔️:                 if (c < '0' || c > '9') {
        // Gemmi❗✔️:                     return std::make_pair(make_error(p, ERROR_INVALID_NUMBER), TYPE_NULL);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                 do {
        // Gemmi❗✔️:                     ++p;
        // Gemmi❗✔️:                     if (SAJSON_UNLIKELY(at_eof(p))) {
        // Gemmi❗✔️:                         return std::make_pair(make_error(p, ERROR_UNEXPECTED_END), TYPE_NULL);
        // Gemmi❗✔️:                     }
        // Gemmi❗✔️:                     c = *p;
        // Gemmi❗✔️:                 } while (c >= '0' && c <= '9');
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:
        // Gemmi❗✔️:             if ('.' == *p) {
        // Gemmi❗✔️:                 ++p;
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(at_eof(p))) {
        // Gemmi❗✔️:                     return std::make_pair(make_error(p, ERROR_UNEXPECTED_END), TYPE_NULL);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 char c = *p;
        // Gemmi❗✔️:                 if (c < '0' || c > '9') {
        // Gemmi❗✔️:                     return std::make_pair(make_error(p, ERROR_INVALID_NUMBER), TYPE_NULL);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                 do {
        // Gemmi❗✔️:                     ++p;
        // Gemmi❗✔️:                     if (SAJSON_UNLIKELY(at_eof(p))) {
        // Gemmi❗✔️:                         return std::make_pair(make_error(p, ERROR_UNEXPECTED_END), TYPE_NULL);
        // Gemmi❗✔️:                     }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                     c = *p;
        // Gemmi❗✔️:                 } while (c >= '0' && c <= '9');
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:
        // Gemmi❗✔️:             char e = *p;
        // Gemmi❗✔️:             if ('e' == e || 'E' == e) {
        // Gemmi❗✔️:                 ++p;
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(at_eof(p))) {
        // Gemmi❗✔️:                     return std::make_pair(make_error(p, ERROR_UNEXPECTED_END), TYPE_NULL);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                 if ('-' == *p) {
        // Gemmi❗✔️:                     ++p;
        // Gemmi❗✔️:                     if (SAJSON_UNLIKELY(at_eof(p))) {
        // Gemmi❗✔️:                         return std::make_pair(make_error(p, ERROR_UNEXPECTED_END), TYPE_NULL);
        // Gemmi❗✔️:                     }
        // Gemmi❗✔️:                 } else if ('+' == *p) {
        // Gemmi❗✔️:                     ++p;
        // Gemmi❗✔️:                     if (SAJSON_UNLIKELY(at_eof(p))) {
        // Gemmi❗✔️:                         return std::make_pair(make_error(p, ERROR_UNEXPECTED_END), TYPE_NULL);
        // Gemmi❗✔️:                     }
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                 char c = *p;
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(c < '0' || c > '9')) {
        // Gemmi❗✔️:                     return std::make_pair(make_error(p, ERROR_MISSING_EXPONENT), TYPE_NULL);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 for (;;) {
        // Gemmi❗✔️:                     ++p;
        // Gemmi❗✔️:                     if (SAJSON_UNLIKELY(at_eof(p))) {
        // Gemmi❗✔️:                         return std::make_pair(make_error(p, ERROR_UNEXPECTED_END), TYPE_NULL);
        // Gemmi❗✔️:                     }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                     c = *p;
        // Gemmi❗✔️:                     if (c < '0' || c > '9') {
        // Gemmi❗✔️:                         break;
        // Gemmi❗✔️:                     }
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:
        // Gemmi❗✔️:             bool success;
        // Gemmi❗✔️:             size_t* out = allocator.reserve(2, &success);
        // Gemmi❗✔️:             if (SAJSON_UNLIKELY(!success)) {
        // Gemmi❗✔️:                 return std::make_pair(oom(p), TYPE_NULL);
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:             out[0] = start;
        // Gemmi❗✔️:             out[1] = p - input.get_data();
        // Gemmi❗✔️:             return std::make_pair(p, TYPE_DOUBLE);
        // Gemmi❗✔️:         }
        // Behavior review: cursor and original digit bytes are retained; the
        // source requires one lookahead byte after every complete number.
        // Complexity review: one O(number length) scan, no conversion or
        // allocation for the numeric lexeme; source reserves two AST words.
        let start = self.position;
        let mut p = start;
        if self.input.get(p) == Some(&b'-') {
            p += 1;
            if p == self.input.len() {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::UnexpectedEnd,
                    None,
                ));
            }
        }
        if self.input.get(p) == Some(&b'0') {
            p += 1;
            if p == self.input.len() {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::UnexpectedEnd,
                    None,
                ));
            }
        } else {
            let Some(&first) = self.input.get(p) else {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::UnexpectedEnd,
                    None,
                ));
            };
            if !first.is_ascii_digit() {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::InvalidNumber,
                    None,
                ));
            }
            loop {
                p += 1;
                let Some(&next) = self.input.get(p) else {
                    return Err(JsonError::at(
                        &self.input,
                        p,
                        JsonErrorKind::UnexpectedEnd,
                        None,
                    ));
                };
                if !next.is_ascii_digit() {
                    break;
                }
            }
        }
        if self.input[p] == b'.' {
            p += 1;
            let Some(&first) = self.input.get(p) else {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::UnexpectedEnd,
                    None,
                ));
            };
            if !first.is_ascii_digit() {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::InvalidNumber,
                    None,
                ));
            }
            loop {
                p += 1;
                let Some(&next) = self.input.get(p) else {
                    return Err(JsonError::at(
                        &self.input,
                        p,
                        JsonErrorKind::UnexpectedEnd,
                        None,
                    ));
                };
                if !next.is_ascii_digit() {
                    break;
                }
            }
        }
        if matches!(self.input[p], b'e' | b'E') {
            p += 1;
            if p == self.input.len() {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::UnexpectedEnd,
                    None,
                ));
            }
            if matches!(self.input[p], b'+' | b'-') {
                p += 1;
                if p == self.input.len() {
                    return Err(JsonError::at(
                        &self.input,
                        p,
                        JsonErrorKind::UnexpectedEnd,
                        None,
                    ));
                }
            }
            if !self.input[p].is_ascii_digit() {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::MissingExponent,
                    None,
                ));
            }
            loop {
                p += 1;
                let Some(&next) = self.input.get(p) else {
                    return Err(JsonError::at(
                        &self.input,
                        p,
                        JsonErrorKind::UnexpectedEnd,
                        None,
                    ));
                };
                if !next.is_ascii_digit() {
                    break;
                }
            }
        }
        self.position = p;
        Ok((start..p, JsonNumberKind::Double))
    }

    fn read_hex(&mut self) -> Result<u32, JsonError> {
        // Gemmi❗✔️:         char* read_hex(char* p, unsigned& u) {
        // Gemmi❗✔️:             unsigned v = 0;
        // Gemmi❗✔️:             int i = 4;
        // Gemmi❗✔️:             while (i--) {
        // Gemmi❗✔️:                 unsigned char c = *p++;
        // Gemmi❗✔️:                 if (c >= '0' && c <= '9') {
        // Gemmi❗✔️:                     c -= '0';
        // Gemmi❗✔️:                 } else if (c >= 'a' && c <= 'f') {
        // Gemmi❗✔️:                     c = c - 'a' + 10;
        // Gemmi❗✔️:                 } else if (c >= 'A' && c <= 'F') {
        // Gemmi❗✔️:                     c = c - 'A' + 10;
        // Gemmi❗✔️:                 } else {
        // Gemmi❗✔️:                     return make_error(p, ERROR_INVALID_UNICODE_ESCAPE);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:                 v = (v << 4) + c;
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:
        // Gemmi❗✔️:             u = v;
        // Gemmi❗✔️:             return p;
        // Gemmi❗✔️:         }
        // Behavior: caller normally checks four bytes; the private helper
        // also returns UnexpectedEnd safely for a truncated direct call.
        // Complexity: four O(1) reads and no allocation.
        let mut value = 0u32;
        for _ in 0..4 {
            let Some(&byte) = self.input.get(self.position) else {
                return Err(JsonError::at(
                    &self.input,
                    self.position,
                    JsonErrorKind::UnexpectedEnd,
                    None,
                ));
            };
            self.position += 1;
            let digit = match byte {
                b'0'..=b'9' => byte - b'0',
                b'a'..=b'f' => byte - b'a' + 10,
                b'A'..=b'F' => byte - b'A' + 10,
                _ => {
                    return Err(JsonError::at(
                        &self.input,
                        self.position,
                        JsonErrorKind::InvalidUnicodeEscape,
                        None,
                    ));
                }
            };
            value = (value << 4) + u32::from(digit);
        }
        Ok(value)
    }

    fn write_utf8(&self, codepoint: u32, output: &mut [u8]) -> Result<usize, JsonError> {
        // Gemmi❗✔️:         void write_utf8(unsigned codepoint, char*& end) {
        // Gemmi❗✔️:             if (codepoint < 0x80) {
        // Gemmi❗✔️:                 *end++ = codepoint;
        // Gemmi❗✔️:             } else if (codepoint < 0x800) {
        // Gemmi❗✔️:                 *end++ = 0xC0 | (codepoint >> 6);
        // Gemmi❗✔️:                 *end++ = 0x80 | (codepoint & 0x3F);
        // Gemmi❗✔️:             } else if (codepoint < 0x10000) {
        // Gemmi❗✔️:                 *end++ = 0xE0 | (codepoint >> 12);
        // Gemmi❗✔️:                 *end++ = 0x80 | ((codepoint >> 6) & 0x3F);
        // Gemmi❗✔️:                 *end++ = 0x80 | (codepoint & 0x3F);
        // Gemmi❗✔️:             } else {
        // Gemmi❗✔️:                 assert(codepoint < 0x200000);
        // Gemmi❗✔️:                 *end++ = 0xF0 | (codepoint >> 18);
        // Gemmi❗✔️:                 *end++ = 0x80 | ((codepoint >> 12) & 0x3F);
        // Gemmi❗✔️:                 *end++ = 0x80 | ((codepoint >> 6) & 0x3F);
        // Gemmi❗✔️:                 *end++ = 0x80 | (codepoint & 0x3F);
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:         }
        // Behavior: byte emission follows the source even for direct helper
        // calls with surrogate values; invalid source-assert ranges and output
        // capacity are explicit safe errors, not claims of upstream parity.
        // Complexity: O(1), at most four output bytes and no prefix copy.
        let bytes: ([u8; 4], usize) = if codepoint < 0x80 {
            ([codepoint as u8, 0, 0, 0], 1)
        } else if codepoint < 0x800 {
            (
                [
                    0xC0 | (codepoint >> 6) as u8,
                    0x80 | (codepoint & 0x3F) as u8,
                    0,
                    0,
                ],
                2,
            )
        } else if codepoint < 0x10000 {
            (
                [
                    0xE0 | (codepoint >> 12) as u8,
                    0x80 | ((codepoint >> 6) & 0x3F) as u8,
                    0x80 | (codepoint & 0x3F) as u8,
                    0,
                ],
                3,
            )
        } else if codepoint < 0x200000 {
            (
                [
                    0xF0 | (codepoint >> 18) as u8,
                    0x80 | ((codepoint >> 12) & 0x3F) as u8,
                    0x80 | ((codepoint >> 6) & 0x3F) as u8,
                    0x80 | (codepoint & 0x3F) as u8,
                ],
                4,
            )
        } else {
            return Err(JsonError::at(
                &self.input,
                self.position,
                JsonErrorKind::InvalidUnicodeEscape,
                None,
            ));
        };
        if bytes.1 > output.len() {
            return Err(JsonError::at(
                &self.input,
                self.position,
                JsonErrorKind::OutOfMemory,
                None,
            ));
        }
        output[..bytes.1].copy_from_slice(&bytes.0[..bytes.1]);
        Ok(bytes.1)
    }

    fn parse_string_fast(&mut self) -> Result<JsonStringFast, JsonError> {
        // Gemmi❗✔️:         char* parse_string(char* p, size_t* tag) {
        // Gemmi❗✔️:             using namespace internal;
        // Gemmi❗✔️:
        // Gemmi❗✔️:             ++p; // "
        // Gemmi❗✔️:             size_t start = p - input.get_data();
        // Gemmi❗✔️:             char* input_end_local = input_end;
        // Gemmi❗✔️:             while (input_end_local - p >= 4) {
        // Gemmi❗✔️:                 if (!is_plain_string_character(p[0])) { goto found; }
        // Gemmi❗✔️:                 if (!is_plain_string_character(p[1])) { p += 1; goto found; }
        // Gemmi❗✔️:                 if (!is_plain_string_character(p[2])) { p += 2; goto found; }
        // Gemmi❗✔️:                 if (!is_plain_string_character(p[3])) { p += 3; goto found; }
        // Gemmi❗✔️:                 p += 4;
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:             for (;;) {
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(p >= input_end_local)) {
        // Gemmi❗✔️:                     return make_error(p, ERROR_UNEXPECTED_END);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                 if (!is_plain_string_character(*p)) {
        // Gemmi❗✔️:                     break;
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                 ++p;
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:         found:
        // Gemmi❗✔️:             if (SAJSON_LIKELY(*p == '"')) {
        // Gemmi❗✔️:                 tag[0] = start;
        // Gemmi❗✔️:                 tag[1] = p - input.get_data();
        // Gemmi❗✔️:                 *p = '\0';
        // Gemmi❗✔️:                 return p + 1;
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:
        // Gemmi❗✔️:             if (*p >= 0 && *p < 0x20) {
        // Gemmi❗✔️:                 return make_error(p, ERROR_ILLEGAL_CODEPOINT, static_cast<int>(*p));
        // Gemmi❗✔️:             } else {
        // Gemmi❗✔️:                 // backslash or >0x7f
        // Gemmi❗✔️:                 return parse_string_slow(p, tag, start);
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:         inline bool is_plain_string_character(char c) {
        // Gemmi❗✔️:             //return c >= 0x20 && c <= 0x7f && c != 0x22 && c != 0x5c;
        // Gemmi❗✔️:             return (globals::parse_flags[static_cast<unsigned char>(c)] & 1) != 0;
        // Gemmi❗✔️:         }
        // Behavior: input owns its bytes; source's in-place NUL terminator is
        // represented by a span. Slow-path handoff is retained for J07–J09.
        // Complexity: O(length), four-byte unrolled scan, no prefix clone.
        let start = self.position + 1;
        let mut p = start;
        let plain = |byte: u8| byte >= 0x20 && byte < 0x80 && !matches!(byte, b'"' | b'\\');
        while self.input.len().saturating_sub(p) >= 4 {
            let window = &self.input[p..p + 4];
            if let Some(index) = window.iter().position(|&byte| !plain(byte)) {
                p += index;
                break;
            }
            p += 4;
        }
        if p == self.input.len() || self.input.get(p).is_some_and(|&byte| plain(byte)) {
            loop {
                let Some(&byte) = self.input.get(p) else {
                    return Err(JsonError::at(
                        &self.input,
                        p,
                        JsonErrorKind::UnexpectedEnd,
                        None,
                    ));
                };
                if !plain(byte) {
                    break;
                }
                p += 1;
            }
        }
        let byte = self.input[p];
        if byte == b'"' {
            self.input[p] = 0;
            self.position = p + 1;
            return Ok(JsonStringFast::Plain(start..p));
        }
        if byte < 0x20 {
            return Err(JsonError::at(
                &self.input,
                p,
                JsonErrorKind::IllegalCodepoint,
                Some(i32::from(byte)),
            ));
        }
        self.position = p;
        Ok(JsonStringFast::Slow { start, at: p })
    }

    fn parse_string_simple(
        &mut self,
        start: usize,
        at: usize,
    ) -> Result<JsonStringProgress, JsonError> {
        // Gemmi❗✔️:         char* parse_string_slow(char* p, size_t* tag, size_t start) {
        // Gemmi❗✔️:             char* end = p;
        // Gemmi❗✔️:             char* input_end_local = input_end;
        // Gemmi❗✔️:
        // Gemmi❗✔️:             for (;;) {
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(p >= input_end_local)) {
        // Gemmi❗✔️:                     return make_error(p, ERROR_UNEXPECTED_END);
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                 if (SAJSON_UNLIKELY(*p >= 0 && *p < 0x20)) {
        // Gemmi❗✔️:                     return make_error(p, ERROR_ILLEGAL_CODEPOINT, static_cast<int>(*p));
        // Gemmi❗✔️:                 }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                 switch (*p) {
        // Gemmi❗✔️:                     case '"':
        // Gemmi❗✔️:                         tag[0] = start;
        // Gemmi❗✔️:                         tag[1] = end - input.get_data();
        // Gemmi❗✔️:                         *end = '\0';
        // Gemmi❗✔️:                         return p + 1;
        // Gemmi❗✔️:
        // Gemmi❗✔️:                     case '\\':
        // Gemmi❗✔️:                         ++p;
        // Gemmi❗✔️:                         if (SAJSON_UNLIKELY(p >= input_end_local)) {
        // Gemmi❗✔️:                             return make_error(p, ERROR_UNEXPECTED_END);
        // Gemmi❗✔️:                         }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                         char replacement;
        // Gemmi❗✔️:                         switch (*p) {
        // Gemmi❗✔️:                             case '"': replacement = '"'; goto replace;
        // Gemmi❗✔️:                             case '\\': replacement = '\\'; goto replace;
        // Gemmi❗✔️:                             case '/': replacement = '/'; goto replace;
        // Gemmi❗✔️:                             case 'b': replacement = '\b'; goto replace;
        // Gemmi❗✔️:                             case 'f': replacement = '\f'; goto replace;
        // Gemmi❗✔️:                             case 'n': replacement = '\n'; goto replace;
        // Gemmi❗✔️:                             case 'r': replacement = '\r'; goto replace;
        // Gemmi❗✔️:                             case 't': replacement = '\t'; goto replace;
        // Gemmi❗✔️:                             replace:
        // Gemmi❗✔️:                                 *end++ = replacement;
        // Gemmi❗✔️:                                 ++p;
        // Gemmi❗✔️:                                 break;
        // Gemmi❗✔️:                             case 'u': {
        // Gemmi❗✔️:                                 ++p;
        // Gemmi❗✔️:                                 if (SAJSON_UNLIKELY(!has_remaining_characters(p, 4))) {
        // Gemmi❗✔️:                                     return make_error(p, ERROR_UNEXPECTED_END);
        // Gemmi❗✔️:                                 }
        // Gemmi❗✔️:                                 unsigned u = 0; // gcc's complaining that this could be used uninitialized. wrong.
        // Gemmi❗✔️:                                 p = read_hex(p, u);
        // Gemmi❗✔️:                                 if (!p) {
        // Gemmi❗✔️:                                     return 0;
        // Gemmi❗✔️:                                 }
        // Gemmi❗✔️:                                 if (u >= 0xD800 && u <= 0xDBFF) {
        // Gemmi❗✔️:                                     if (SAJSON_UNLIKELY(!has_remaining_characters(p, 6))) {
        // Gemmi❗✔️:                                         return make_error(p, ERROR_UNEXPECTED_END_OF_UTF16);
        // Gemmi❗✔️:                                     }
        // Gemmi❗✔️:                                     char p0 = p[0];
        // Gemmi❗✔️:                                     char p1 = p[1];
        // Gemmi❗✔️:                                     if (p0 != '\\' || p1 != 'u') {
        // Gemmi❗✔️:                                         return make_error(p, ERROR_EXPECTED_U);
        // Gemmi❗✔️:                                     }
        // Gemmi❗✔️:                                     p += 2;
        // Gemmi❗✔️:                                     unsigned v = 0; // gcc's complaining that this could be used uninitialized. wrong.
        // Gemmi❗✔️:                                     p = read_hex(p, v);
        // Gemmi❗✔️:                                     if (!p) {
        // Gemmi❗✔️:                                         return p;
        // Gemmi❗✔️:                                     }
        // Gemmi❗✔️:
        // Gemmi❗✔️:                                     if (v < 0xDC00 || v > 0xDFFF) {
        // Gemmi❗✔️:                                         return make_error(p, ERROR_INVALID_UTF16_TRAIL_SURROGATE);
        // Gemmi❗✔️:                                     }
        // Gemmi❗✔️:                                     u = 0x10000 + (((u - 0xD800) << 10) | (v - 0xDC00));
        // Gemmi❗✔️:                                 }
        // Gemmi❗✔️:                                 write_utf8(u, end);
        // Gemmi❗✔️:                                 break;
        // Gemmi❗✔️:                             }
        // Gemmi❗✔️:                             default:
        // Gemmi❗✔️:                                 return make_error(p, ERROR_UNKNOWN_ESCAPE);
        // Gemmi❗✔️:                         }
        // Gemmi❗✔️:                         break;
        // Behavior: source simple escapes mutate the single owned byte buffer;
        // Unicode and raw multibyte branches retain exact cursor/end state for
        // their later complete packets, with no invented fallback.
        // Complexity: one O(length) scan and O(1) work per simple escape;
        // no prefix copy or secondary string allocation.
        let mut p = at;
        let mut end = at;
        loop {
            let Some(&byte) = self.input.get(p) else {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::UnexpectedEnd,
                    None,
                ));
            };
            if byte < 0x20 {
                return Err(JsonError::at(
                    &self.input,
                    p,
                    JsonErrorKind::IllegalCodepoint,
                    Some(i32::from(byte)),
                ));
            }
            match byte {
                b'"' => {
                    self.input[end] = 0;
                    self.position = p + 1;
                    return Ok(JsonStringProgress::Complete(start..end));
                }
                b'\\' => {
                    p += 1;
                    let Some(&escaped) = self.input.get(p) else {
                        return Err(JsonError::at(
                            &self.input,
                            p,
                            JsonErrorKind::UnexpectedEnd,
                            None,
                        ));
                    };
                    let replacement = match escaped {
                        b'"' => b'"',
                        b'\\' => b'\\',
                        b'/' => b'/',
                        b'b' => 8,
                        b'f' => 12,
                        b'n' => b'\n',
                        b'r' => b'\r',
                        b't' => b'\t',
                        b'u' => {
                            // Gemmi❗❌:                             case 'u': {
                            // Gemmi❗❌:                                 ++p;
                            // Gemmi❗❌:                                 if (SAJSON_UNLIKELY(!has_remaining_characters(p, 4))) {
                            // Gemmi❗❌:                                     return make_error(p, ERROR_UNEXPECTED_END);
                            // Gemmi❗❌:                                 }
                            // Gemmi❗❌:                                 unsigned u = 0; // gcc's complaining that this could be used uninitialized. wrong.
                            // Gemmi❗❌:                                 p = read_hex(p, u);
                            // Gemmi❗❌:                                 if (!p) {
                            // Gemmi❗❌:                                     return 0;
                            // Gemmi❗❌:                                 }
                            // Gemmi❗❌:                                 if (u >= 0xD800 && u <= 0xDBFF) {
                            // Gemmi❗❌:                                     if (SAJSON_UNLIKELY(!has_remaining_characters(p, 6))) {
                            // Gemmi❗❌:                                         return make_error(p, ERROR_UNEXPECTED_END_OF_UTF16);
                            // Gemmi❗❌:                                     }
                            // Gemmi❗❌:                                     char p0 = p[0];
                            // Gemmi❗❌:                                     char p1 = p[1];
                            // Gemmi❗❌:                                     if (p0 != '\\' || p1 != 'u') {
                            // Gemmi❗❌:                                         return make_error(p, ERROR_EXPECTED_U);
                            // Gemmi❗❌:                                     }
                            // Gemmi❗❌:                                     p += 2;
                            // Gemmi❗❌:                                     unsigned v = 0; // gcc's complaining that this could be used uninitialized. wrong.
                            // Gemmi❗❌:                                     p = read_hex(p, v);
                            // Gemmi❗❌:                                     if (!p) {
                            // Gemmi❗❌:                                         return p;
                            // Gemmi❗❌:                                     }
                            // Gemmi❗❌:
                            // Gemmi❗❌:                                     if (v < 0xDC00 || v > 0xDFFF) {
                            // Gemmi❗❌:                                         return make_error(p, ERROR_INVALID_UTF16_TRAIL_SURROGATE);
                            // Gemmi❗❌:                                     }
                            // Gemmi❗❌:                                     u = 0x10000 + (((u - 0xD800) << 10) | (v - 0xDC00));
                            // Gemmi❗❌:                                 }
                            // Gemmi❗✔️:                                 write_utf8(u, end);
                            // Gemmi❗❌:                                 break;
                            // Gemmi❗❌:                             }
                            // Behavior: source checks the first four hex bytes, then
                            // requires a second complete escape only for a high surrogate.
                            // A lone low surrogate follows source UTF-8 emission.
                            // Complexity: the shared writer emits at most four bytes
                            // through a stack slice and copies only the written prefix;
                            // no per-escape heap allocation or prefix rescan.
                            p += 1;
                            if self.input.len().saturating_sub(p) < 4 {
                                return Err(JsonError::at(
                                    &self.input,
                                    p,
                                    JsonErrorKind::UnexpectedEnd,
                                    None,
                                ));
                            }
                            self.position = p;
                            let mut codepoint = self.read_hex()?;
                            p = self.position;
                            if (0xD800..=0xDBFF).contains(&codepoint) {
                                if self.input.len().saturating_sub(p) < 6 {
                                    return Err(JsonError::at(
                                        &self.input,
                                        p,
                                        JsonErrorKind::UnexpectedEndOfUtf16,
                                        None,
                                    ));
                                }
                                if self.input[p] != b'\\' || self.input[p + 1] != b'u' {
                                    return Err(JsonError::at(
                                        &self.input,
                                        p,
                                        JsonErrorKind::ExpectedU,
                                        None,
                                    ));
                                }
                                self.position = p + 2;
                                let trail = self.read_hex()?;
                                p = self.position;
                                if !(0xDC00..=0xDFFF).contains(&trail) {
                                    return Err(JsonError::at(
                                        &self.input,
                                        p,
                                        JsonErrorKind::InvalidUtf16TrailSurrogate,
                                        None,
                                    ));
                                }
                                codepoint =
                                    0x10000 + (((codepoint - 0xD800) << 10) | (trail - 0xDC00));
                            }
                            let mut encoded = [0u8; 4];
                            let written = self.write_utf8(codepoint, &mut encoded)?;
                            let decoded_end = end + written;
                            self.input[end..decoded_end].copy_from_slice(&encoded[..written]);
                            end = decoded_end;
                            continue;
                        }
                        _ => {
                            return Err(JsonError::at(
                                &self.input,
                                p,
                                JsonErrorKind::UnknownEscape,
                                None,
                            ));
                        }
                    };
                    self.input[end] = replacement;
                    end += 1;
                    p += 1;
                }
                0x00..=0x7f => {
                    self.input[end] = byte;
                    end += 1;
                    p += 1;
                }
                _ => {
                    // Gemmi❗✔️:                     default:
                    // Gemmi❗✔️:                         // validate UTF-8
                    // Gemmi❗✔️:                         unsigned char c0 = p[0];
                    // Gemmi❗✔️:                         if (c0 < 128) {
                    // Gemmi❗✔️:                             *end++ = *p++;
                    // Gemmi❗✔️:                         } else if (c0 < 224) {
                    // Gemmi❗✔️:                             if (SAJSON_UNLIKELY(!has_remaining_characters(p, 2))) {
                    // Gemmi❗✔️:                                 return unexpected_end(p);
                    // Gemmi❗✔️:                             }
                    // Gemmi❗✔️:                             unsigned char c1 = p[1];
                    // Gemmi❗✔️:                             if (c1 < 128 || c1 >= 192) {
                    // Gemmi❗✔️:                                 return make_error(p + 1, ERROR_INVALID_UTF8);
                    // Gemmi❗✔️:                             }
                    // Gemmi❗✔️:                             end[0] = c0;
                    // Gemmi❗✔️:                             end[1] = c1;
                    // Gemmi❗✔️:                             end += 2;
                    // Gemmi❗✔️:                             p += 2;
                    // Gemmi❗✔️:                         } else if (c0 < 240) {
                    // Gemmi❗✔️:                             if (SAJSON_UNLIKELY(!has_remaining_characters(p, 3))) {
                    // Gemmi❗✔️:                                 return unexpected_end(p);
                    // Gemmi❗✔️:                             }
                    // Gemmi❗✔️:                             unsigned char c1 = p[1];
                    // Gemmi❗✔️:                             if (c1 < 128 || c1 >= 192) {
                    // Gemmi❗✔️:                                 return make_error(p + 1, ERROR_INVALID_UTF8);
                    // Gemmi❗✔️:                             }
                    // Gemmi❗✔️:                             unsigned char c2 = p[2];
                    // Gemmi❗✔️:                             if (c2 < 128 || c2 >= 192) {
                    // Gemmi❗✔️:                                 return make_error(p + 2, ERROR_INVALID_UTF8);
                    // Gemmi❗✔️:                             }
                    // Gemmi❗✔️:                             end[0] = c0;
                    // Gemmi❗✔️:                             end[1] = c1;
                    // Gemmi❗✔️:                             end[2] = c2;
                    // Gemmi❗✔️:                             end += 3;
                    // Gemmi❗✔️:                             p += 3;
                    // Gemmi❗✔️:                         } else if (c0 < 248) {
                    // Gemmi❗✔️:                             if (SAJSON_UNLIKELY(!has_remaining_characters(p, 4))) {
                    // Gemmi❗✔️:                                 return unexpected_end(p);
                    // Gemmi❗✔️:                             }
                    // Gemmi❗✔️:                             unsigned char c1 = p[1];
                    // Gemmi❗✔️:                             if (c1 < 128 || c1 >= 192) {
                    // Gemmi❗✔️:                                 return make_error(p + 1, ERROR_INVALID_UTF8);
                    // Gemmi❗✔️:                             }
                    // Gemmi❗✔️:                             unsigned char c2 = p[2];
                    // Gemmi❗✔️:                             if (c2 < 128 || c2 >= 192) {
                    // Gemmi❗✔️:                                 return make_error(p + 2, ERROR_INVALID_UTF8);
                    // Gemmi❗✔️:                             }
                    // Gemmi❗✔️:                             unsigned char c3 = p[3];
                    // Gemmi❗✔️:                             if (c3 < 128 || c3 >= 192) {
                    // Gemmi❗✔️:                                 return make_error(p + 3, ERROR_INVALID_UTF8);
                    // Gemmi❗✔️:                             }
                    // Gemmi❗✔️:                             end[0] = c0;
                    // Gemmi❗✔️:                             end[1] = c1;
                    // Gemmi❗✔️:                             end[2] = c2;
                    // Gemmi❗✔️:                             end[3] = c3;
                    // Gemmi❗✔️:                             end += 4;
                    // Gemmi❗✔️:                             p += 4;
                    // Gemmi❗✔️:                         } else {
                    // Gemmi❗✔️:                             return make_error(p, ERROR_INVALID_UTF8);
                    // Gemmi❗✔️:                         }
                    // Gemmi❗✔️:                         break;
                    // Behavior: source validates only continuation-byte ranges,
                    // including its acceptance of overlong and high code-point
                    // byte forms; the later text consumer handles representability.
                    // Complexity: O(width) per sequence, no allocation.
                    let width = if byte < 224 {
                        2
                    } else if byte < 240 {
                        3
                    } else if byte < 248 {
                        4
                    } else {
                        return Err(JsonError::at(
                            &self.input,
                            p,
                            JsonErrorKind::InvalidUtf8,
                            None,
                        ));
                    };
                    if self.input.len().saturating_sub(p) < width {
                        return Err(JsonError::at(
                            &self.input,
                            p,
                            JsonErrorKind::UnexpectedEnd,
                            None,
                        ));
                    }
                    for offset in 1..width {
                        if !(128..192).contains(&self.input[p + offset]) {
                            return Err(JsonError::at(
                                &self.input,
                                p + offset,
                                JsonErrorKind::InvalidUtf8,
                                None,
                            ));
                        }
                    }
                    self.input.copy_within(p..p + width, end);
                    p += width;
                    end += width;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        JsonArena, JsonArenaValue, JsonContainerTransition, JsonCursor, JsonError, JsonErrorKind,
        JsonLiteralKind, JsonNumberKind, JsonParseStack, JsonStackKind, JsonStringFast,
        JsonStringProgress,
    };

    #[test]
    fn bio_read_j13_nested_duplicate_keys_and_raw_number_lexeme() {
        let mut cursor = JsonCursor::new(br#"{"a":-0.0e+2,"a":[true,null,"x"]}"#);
        let (arena, root) = cursor.parse_document().unwrap();
        let (first_key, first_value) = arena.object_entry(root, 0).unwrap();
        assert_eq!(&cursor.input[first_key.clone()], b"a");
        let JsonArenaValue::Number(number) = first_value else {
            panic!("raw number")
        };
        assert_eq!(&cursor.input[number.clone()], b"-0.0e+2");
        let (second_key, second_value) = arena.object_entry(root, 1).unwrap();
        assert_eq!(&cursor.input[second_key.clone()], b"a");
        let JsonArenaValue::Array(_) = second_value else {
            panic!("array")
        };
        let JsonArenaValue::Object(object_range) = arena.value(root).unwrap() else {
            panic!("object")
        };
        assert_eq!(object_range.len(), 2);
        let second_id = arena.object_items[object_range.start + 1].1;
        assert_eq!(
            arena.array_element(second_id, 0),
            Some(&JsonArenaValue::Boolean(true))
        );
        assert_eq!(
            arena.array_element(second_id, 1),
            Some(&JsonArenaValue::Null)
        );
        let JsonArenaValue::String(text) = arena.array_element(second_id, 2).unwrap() else {
            panic!("string")
        };
        assert_eq!(&cursor.input[text.clone()], b"x");
        assert_eq!(arena.object_entry(root, 2), None);

        let mut array = JsonCursor::new(b" [ { } , [ ] ] \r\n");
        let (arena, root) = array.parse_document().unwrap();
        assert!(matches!(arena.value(root), Some(JsonArenaValue::Array(_))));
        assert!(matches!(
            arena.array_element(root, 0),
            Some(JsonArenaValue::Object(_))
        ));
        assert!(matches!(
            arena.array_element(root, 1),
            Some(JsonArenaValue::Array(_))
        ));
    }

    #[test]
    fn bio_read_j13_root_trailing_and_line_error_contract() {
        for (input, kind, offset) in [
            (b"".as_slice(), JsonErrorKind::MissingRootElement, 0),
            (b" \t\n".as_slice(), JsonErrorKind::MissingRootElement, 3),
            (b"null".as_slice(), JsonErrorKind::BadRoot, 0),
            (b"1".as_slice(), JsonErrorKind::BadRoot, 0),
            (b"\"x\"".as_slice(), JsonErrorKind::BadRoot, 0),
            (b"{}x".as_slice(), JsonErrorKind::ExpectedEndOfInput, 2),
            (b"[1,]".as_slice(), JsonErrorKind::ExpectedValue, 3),
            (b"{\"a\":\n}".as_slice(), JsonErrorKind::ExpectedValue, 6),
        ] {
            let error = match JsonCursor::new(input).parse_document() {
                Err(error) => error,
                Ok(_) => panic!("accepted invalid root: {input:?}"),
            };
            assert_eq!(
                (error.kind(), error.byte_offset()),
                (kind, offset),
                "{input:?}"
            );
        }
        let error = match JsonCursor::new(b"{\"a\":\n}").parse_document() {
            Err(error) => error,
            Ok(_) => panic!("accepted invalid value"),
        };
        assert_eq!(
            (error.line(), error.column(), error.to_string()),
            (2, 1, "expected value".to_owned())
        );
    }

    #[test]
    fn bio_read_j12_empty_and_nested_container_closure() {
        for (bytes, kind, expected) in [
            (
                b"[]".as_slice(),
                JsonStackKind::Array,
                JsonArenaValue::Array(0..0),
            ),
            (
                b"{}".as_slice(),
                JsonStackKind::Object,
                JsonArenaValue::Object(0..0),
            ),
        ] {
            let mut cursor = JsonCursor::new(bytes);
            let mut stack = JsonParseStack::default();
            let mut arena = JsonArena::default();
            stack.open(kind);
            let JsonContainerTransition::Closed(root) = cursor
                .container_transition(&mut stack, &mut arena, true)
                .unwrap()
            else {
                panic!("not closed")
            };
            assert_eq!(arena.value(root), Some(&expected));
            assert_eq!(cursor.position, bytes.len());
            assert!(stack.frames.is_empty());
        }

        let mut cursor = JsonCursor::new(b"[{}]");
        let mut stack = JsonParseStack::default();
        let mut arena = JsonArena::default();
        stack.open(JsonStackKind::Array);
        assert_eq!(
            cursor
                .container_transition(&mut stack, &mut arena, true)
                .unwrap(),
            JsonContainerTransition::Value
        );
        stack.open(JsonStackKind::Object);
        let JsonContainerTransition::Closed(child) = cursor
            .container_transition(&mut stack, &mut arena, true)
            .unwrap()
        else {
            panic!("child")
        };
        stack.push_value(child).unwrap();
        let JsonContainerTransition::Closed(root) = cursor
            .container_transition(&mut stack, &mut arena, false)
            .unwrap()
        else {
            panic!("root")
        };
        assert_eq!(
            arena.array_element(root, 0),
            Some(&JsonArenaValue::Object(0..0))
        );
        assert_eq!(cursor.position, 4);
    }

    #[test]
    fn bio_read_j12_comma_closing_and_eof_error_precedence() {
        for (input, kind, position, opened, expected, offset) in [
            (
                b"[".as_slice(),
                JsonStackKind::Array,
                0,
                true,
                JsonErrorKind::UnexpectedEnd,
                1,
            ),
            (
                b"{".as_slice(),
                JsonStackKind::Object,
                0,
                true,
                JsonErrorKind::UnexpectedEnd,
                1,
            ),
            (
                b"[1".as_slice(),
                JsonStackKind::Array,
                2,
                false,
                JsonErrorKind::UnexpectedEnd,
                2,
            ),
            (
                b"{\"a\":1".as_slice(),
                JsonStackKind::Object,
                6,
                false,
                JsonErrorKind::UnexpectedEnd,
                6,
            ),
            (
                b"[1 2]".as_slice(),
                JsonStackKind::Array,
                2,
                false,
                JsonErrorKind::ExpectedComma,
                3,
            ),
            (
                b"[1}".as_slice(),
                JsonStackKind::Array,
                2,
                false,
                JsonErrorKind::ExpectedComma,
                2,
            ),
            (
                b"[1,]".as_slice(),
                JsonStackKind::Array,
                2,
                false,
                JsonErrorKind::ExpectedValue,
                3,
            ),
            (
                b"[1,,2]".as_slice(),
                JsonStackKind::Array,
                2,
                false,
                JsonErrorKind::UnexpectedComma,
                3,
            ),
        ] {
            let mut cursor = JsonCursor::new(input);
            let mut stack = JsonParseStack::default();
            let mut arena = JsonArena::default();
            stack.open(kind);
            cursor.position = position;
            let error = cursor
                .container_transition(&mut stack, &mut arena, opened)
                .unwrap_err();
            assert_eq!(
                (error.kind(), error.byte_offset()),
                (expected, offset),
                "{input:?}"
            );
        }

        let mut cursor = JsonCursor::new(b"{\"a\":1,}");
        let mut stack = JsonParseStack::default();
        let mut arena = JsonArena::default();
        stack.open(JsonStackKind::Object);
        cursor.position = 6;
        assert_eq!(
            cursor
                .container_transition(&mut stack, &mut arena, false)
                .unwrap(),
            JsonContainerTransition::ObjectKey
        );
        let error = cursor.parse_object_key_colon(&mut stack).unwrap_err();
        assert_eq!(
            (error.kind(), error.byte_offset()),
            (JsonErrorKind::MissingObjectKey, 7)
        );
    }

    #[test]
    fn bio_read_j11_quoted_keys_whitespace_and_duplicates() {
        let mut cursor = JsonCursor::new(b"{ \t\"a\" \r\n: 1, \"a\":2 }");
        let mut stack = JsonParseStack::default();
        stack.open(JsonStackKind::Object);
        cursor.position = 1;
        cursor.parse_object_key_colon(&mut stack).unwrap();
        assert_eq!(stack.entries, vec![super::JsonStackEntry::Key(4..5)]);
        assert_eq!(cursor.position, 10);
        cursor.position = 14;
        cursor.parse_object_key_colon(&mut stack).unwrap();
        assert_eq!(
            stack.entries,
            vec![
                super::JsonStackEntry::Key(4..5),
                super::JsonStackEntry::Key(15..16)
            ]
        );
        assert_eq!(cursor.position, 18);

        let mut escaped = JsonCursor::new(br#"{"a\n":0}"#);
        let mut escaped_stack = JsonParseStack::default();
        escaped_stack.open(JsonStackKind::Object);
        escaped.position = 1;
        escaped.parse_object_key_colon(&mut escaped_stack).unwrap();
        let super::JsonStackEntry::Key(span) = &escaped_stack.entries[0] else {
            panic!("missing key")
        };
        assert_eq!(&escaped.input[span.clone()], b"a\n");
    }

    #[test]
    fn bio_read_j11_missing_colon_nonstring_key_and_trailing_comma() {
        for (input, position, kind, offset) in [
            (b"{1:2}".as_slice(), 1, JsonErrorKind::MissingObjectKey, 1),
            (b"{\"a\" 1}".as_slice(), 1, JsonErrorKind::ExpectedColon, 5),
            (b"{\"a\"".as_slice(), 1, JsonErrorKind::ExpectedColon, 4),
            (
                b"{\"a\":1,}".as_slice(),
                7,
                JsonErrorKind::MissingObjectKey,
                7,
            ),
            (b"{\"a\":1,".as_slice(), 7, JsonErrorKind::UnexpectedEnd, 7),
        ] {
            let mut cursor = JsonCursor::new(input);
            let mut stack = JsonParseStack::default();
            stack.open(JsonStackKind::Object);
            cursor.position = position;
            let error = cursor.parse_object_key_colon(&mut stack).unwrap_err();
            assert_eq!(
                (error.kind(), error.byte_offset()),
                (kind, offset),
                "{input:?}"
            );
        }
    }

    #[test]
    fn bio_read_j10_empty_ordered_and_duplicate_containers() {
        let mut arena = JsonArena::default();
        let first = arena.push_value(JsonArenaValue::Boolean(true));
        let second = arena.push_value(JsonArenaValue::Null);
        let empty_array = arena.install_array(&[]).unwrap();
        assert_eq!(arena.value(empty_array), Some(&JsonArenaValue::Array(0..0)));
        assert_eq!(arena.array_element(empty_array, 0), None);
        let array = arena.install_array(&[second, first, second]).unwrap();
        assert_eq!(arena.array_element(array, 0), Some(&JsonArenaValue::Null));
        assert_eq!(
            arena.array_element(array, 1),
            Some(&JsonArenaValue::Boolean(true))
        );
        assert_eq!(arena.array_element(array, 2), Some(&JsonArenaValue::Null));
        assert_eq!(arena.array_element(array, 3), None);
        let empty_object = arena.install_object(2, &[]).unwrap();
        assert_eq!(
            arena.value(empty_object),
            Some(&JsonArenaValue::Object(0..0))
        );
        assert_eq!(arena.object_entry(empty_object, 0), None);
        let object = arena
            .install_object(2, &[(0..1, first), (1..2, second), (0..1, array)])
            .unwrap();
        assert_eq!(
            arena.object_entry(object, 0),
            Some((&(0..1), &JsonArenaValue::Boolean(true)))
        );
        assert_eq!(
            arena.object_entry(object, 1),
            Some((&(1..2), &JsonArenaValue::Null))
        );
        assert_eq!(
            arena.object_entry(object, 2),
            Some((&(0..1), &JsonArenaValue::Array(0..3)))
        );
        assert_eq!(arena.object_entry(object, 3), None);
        assert_eq!(arena.value(usize::MAX), None);
        assert_eq!(arena.install_array(&[usize::MAX]), None);
        assert_eq!(arena.install_object(2, &[(2..3, first)]), None);
        assert_eq!(arena.install_object(2, &[(1..0, first)]), None);
    }

    #[test]
    fn bio_read_j10_explicit_nested_stack_and_bounded_access() {
        let mut arena = JsonArena::default();
        let scalar = arena.push_value(JsonArenaValue::Number(4..5));
        let mut stack = JsonParseStack::default();
        stack.open(JsonStackKind::Object);
        stack.push_key(0..1).unwrap();
        stack.open(JsonStackKind::Array);
        stack.push_value(scalar).unwrap();
        let child = stack.finish(&mut arena, 5).unwrap();
        stack.push_value(child).unwrap();
        let root = stack.finish(&mut arena, 5).unwrap();
        assert_eq!(
            arena.object_entry(root, 0),
            Some((&(0..1), &JsonArenaValue::Array(0..1)))
        );
        assert_eq!(
            arena.array_element(child, 0),
            Some(&JsonArenaValue::Number(4..5))
        );
        assert_eq!(arena.array_element(child, 1), None);
        assert_eq!(arena.object_entry(root, 1), None);
        assert_eq!(stack.finish(&mut arena, 5), None);

        stack.open(JsonStackKind::Object);
        stack.push_key(0..1).unwrap();
        assert_eq!(stack.finish(&mut arena, 5), None);
        assert_eq!(stack.frames.len(), 1);
        assert_eq!(stack.entries.len(), 1);
    }

    #[test]
    fn bio_read_j09_source_widths_and_unvalidated_scalar_boundaries() {
        for bytes in [
            b"\xc2\xa2".as_slice(),
            b"\xe2\x82\xac",
            b"\xf0\x9f\x92\xa9",
            b"\xc0\x80",
            b"\xed\xa0\x80",
            b"\xf7\xbf\xbf\xbf",
            b"\x80\x80",
        ] {
            let mut input = vec![b'"', b'X', b'\\', b'n'];
            input.extend_from_slice(bytes);
            input.extend_from_slice(b"Y\"");
            let mut cursor = JsonCursor::new(&input);
            let JsonStringFast::Slow { start, at } = cursor.parse_string_fast().unwrap() else {
                panic!("escaped prefix did not select slow path");
            };
            let JsonStringProgress::Complete(span) = cursor.parse_string_simple(start, at).unwrap()
            else {
                panic!("raw multibyte input did not complete");
            };
            let mut expected = vec![b'X', b'\n'];
            expected.extend_from_slice(bytes);
            expected.push(b'Y');
            assert_eq!(&cursor.input[span.clone()], expected, "{input:?}");
            assert_eq!(cursor.input[span.end], 0);
            assert_eq!(cursor.position, input.len());
        }
    }

    #[test]
    fn bio_read_j09_invalid_continuation_truncation_and_lead_boundaries() {
        for (bytes, kind, offset) in [
            (b"\xc2".as_slice(), JsonErrorKind::UnexpectedEnd, 1),
            (b"\xe2\x82", JsonErrorKind::UnexpectedEnd, 1),
            (b"\xf0\x9f\x92", JsonErrorKind::UnexpectedEnd, 1),
            (b"\xc2A", JsonErrorKind::InvalidUtf8, 2),
            (b"\xe2\x82A", JsonErrorKind::InvalidUtf8, 3),
            (b"\xf0\x9f\x92A", JsonErrorKind::InvalidUtf8, 4),
            (b"\xf8\x80\x80\x80", JsonErrorKind::InvalidUtf8, 1),
        ] {
            let mut input = vec![b'"'];
            input.extend_from_slice(bytes);
            let mut cursor = JsonCursor::new(&input);
            let JsonStringFast::Slow { start, at } = cursor.parse_string_fast().unwrap() else {
                panic!("raw multibyte input did not select slow path");
            };
            let error = cursor.parse_string_simple(start, at).unwrap_err();
            assert_eq!(
                (error.kind(), error.byte_offset()),
                (kind, offset),
                "{input:?}"
            );
        }
    }

    #[test]
    fn bio_read_j08_bmp_nul_and_surrogate_decoded_spans() {
        for (input, expected) in [
            (b"\"\\u0000\"".as_slice(), b"\0".as_slice()),
            (b"\"A\\u0000B\"", b"A\0B"),
            (b"\"\\u007f\"", b"\x7f"),
            (b"\"\\u0080\"", b"\xc2\x80"),
            (b"\"\\u07ff\"", b"\xdf\xbf"),
            (b"\"\\u0800\"", b"\xe0\xa0\x80"),
            (b"\"\\ud7ff\"", b"\xed\x9f\xbf"),
            (b"\"\\ue000\"", b"\xee\x80\x80"),
            (b"\"\\uffff\"", b"\xef\xbf\xbf"),
            (b"\"\\ud800\\udc00\"", b"\xf0\x90\x80\x80"),
            (b"\"\\udbff\\udfff\"", b"\xf4\x8f\xbf\xbf"),
            (b"\"\\udc00\"", b"\xed\xb0\x80"),
            (b"\"A\\n\\u0000B\\u007f\\u0080\\u07ff\\u0800\\ud7ff\\ue000\\uffff\\ud800\\udc00\\udbff\\udfff\\udc00Z\"", b"A\n\0B\x7f\xc2\x80\xdf\xbf\xe0\xa0\x80\xed\x9f\xbf\xee\x80\x80\xef\xbf\xbf\xf0\x90\x80\x80\xf4\x8f\xbf\xbf\xed\xb0\x80Z"),
            (b"\"\\u0080\\u0080\\u0080\"", b"\xc2\x80\xc2\x80\xc2\x80"),
        ] {
            let mut cursor = JsonCursor::new(input);
            let JsonStringFast::Slow { start, at } = cursor.parse_string_fast().unwrap() else {
                panic!("Unicode input did not select slow path");
            };
            let JsonStringProgress::Complete(span) = cursor.parse_string_simple(start, at).unwrap()
            else {
                panic!("Unicode input did not complete");
            };
            assert_eq!(span.len(), expected.len(), "{input:?}");
            assert_eq!(&cursor.input[span.clone()], expected, "{input:?}");
            assert_eq!(cursor.input[span.end], 0, "{input:?}");
            assert_eq!(cursor.position, input.len(), "{input:?}");
        }
    }

    #[test]
    fn bio_read_j08_malformed_surrogates_and_hex_offsets() {
        for (input, kind, offset) in [
            (
                b"\"\\ud800\"".as_slice(),
                JsonErrorKind::UnexpectedEndOfUtf16,
                7,
            ),
            (b"\"\\ud800x0000\"", JsonErrorKind::ExpectedU, 7),
            (
                b"\"\\ud800\\u0041\"",
                JsonErrorKind::InvalidUtf16TrailSurrogate,
                13,
            ),
            (
                b"\"\\ud800\\u12xz\"",
                JsonErrorKind::InvalidUnicodeEscape,
                12,
            ),
            (b"\"\\u12", JsonErrorKind::UnexpectedEnd, 3),
            (b"\"\\u12xz\"", JsonErrorKind::InvalidUnicodeEscape, 6),
        ] {
            let mut cursor = JsonCursor::new(input);
            let JsonStringFast::Slow { start, at } = cursor.parse_string_fast().unwrap() else {
                panic!("malformed Unicode input did not select slow path");
            };
            let error = cursor.parse_string_simple(start, at).unwrap_err();
            assert_eq!(
                (error.kind(), error.byte_offset()),
                (kind, offset),
                "{input:?}"
            );
        }
    }

    #[test]
    fn bio_read_j07_eight_simple_escapes_and_exact_decoded_spans() {
        for (escaped, decoded) in [
            (b'"', b'"'),
            (b'\\', b'\\'),
            (b'/', b'/'),
            (b'b', 8),
            (b'f', 12),
            (b'n', b'\n'),
            (b'r', b'\r'),
            (b't', b'\t'),
        ] {
            let input = [b'"', b'A', b'\\', escaped, b'B', b'"'];
            let mut cursor = JsonCursor::new(&input);
            assert_eq!(
                cursor.parse_string_fast().unwrap(),
                JsonStringFast::Slow { start: 1, at: 2 }
            );
            let JsonStringProgress::Complete(span) = cursor.parse_string_simple(1, 2).unwrap()
            else {
                panic!("simple escape did not complete");
            };
            assert_eq!(span, 1..4);
            assert_eq!(&cursor.input[span.clone()], &[b'A', decoded, b'B']);
            assert_eq!(cursor.input[span.end], 0);
            assert_eq!(cursor.position, input.len());
        }
        let input = b"\"a\\n\\t\\\"z\"";
        let mut cursor = JsonCursor::new(input);
        assert_eq!(
            cursor.parse_string_fast().unwrap(),
            JsonStringFast::Slow { start: 1, at: 2 }
        );
        let JsonStringProgress::Complete(span) = cursor.parse_string_simple(1, 2).unwrap() else {
            panic!("mixed simple escapes did not complete");
        };
        assert_eq!(&cursor.input[span.clone()], b"a\n\t\"z");
        assert_eq!(span, 1..6);
        assert_eq!(cursor.input[span.end], 0);
    }

    #[test]
    fn bio_read_j07_invalid_escape_final_backslash_and_raw_nul() {
        for escaped in [b'0', b'q'] {
            let input = [b'"', b'\\', escaped, b'"'];
            let mut cursor = JsonCursor::new(&input);
            assert_eq!(
                cursor.parse_string_fast().unwrap(),
                JsonStringFast::Slow { start: 1, at: 1 }
            );
            let error = cursor.parse_string_simple(1, 1).unwrap_err();
            assert_eq!(
                (error.kind(), error.byte_offset()),
                (JsonErrorKind::UnknownEscape, 2)
            );
        }
        let input = b"\"abc\\";
        let mut cursor = JsonCursor::new(input);
        assert_eq!(
            cursor.parse_string_fast().unwrap(),
            JsonStringFast::Slow { start: 1, at: 4 }
        );
        let error = cursor.parse_string_simple(1, 4).unwrap_err();
        assert_eq!(
            (error.kind(), error.byte_offset()),
            (JsonErrorKind::UnexpectedEnd, input.len())
        );

        let input = b"\"a\0b\"";
        let mut cursor = JsonCursor::new(input);
        let error = cursor.parse_string_fast().unwrap_err();
        assert_eq!(
            (error.kind(), error.byte_offset()),
            (JsonErrorKind::IllegalCodepoint, 2)
        );
    }

    #[test]
    fn bio_read_j06_plain_string_spans_and_terminator_positions() {
        for input in [
            b"\"\"".as_slice(),
            b"\"abc\"",
            b"\"a b~\x7f\"",
            b"\"abcdefghijklmnopqrstuvwxyz0123456789\"",
        ] {
            let mut cursor = JsonCursor::new(input);
            let result = cursor.parse_string_fast().unwrap();
            assert_eq!(result, JsonStringFast::Plain(1..input.len() - 1));
            assert_eq!(cursor.position, input.len());
        }
        for input in [b"\"".as_slice(), b"\"abc"] {
            let mut cursor = JsonCursor::new(input);
            let error = cursor.parse_string_fast().unwrap_err();
            assert_eq!(
                (error.kind(), error.byte_offset()),
                (JsonErrorKind::UnexpectedEnd, input.len())
            );
        }
    }

    #[test]
    fn bio_read_j06_controls_and_slow_path_handoff() {
        for control in 0..0x20u8 {
            let input = [b'"', b'A', control, b'"'];
            let mut cursor = JsonCursor::new(&input);
            let error = cursor.parse_string_fast().unwrap_err();
            assert_eq!(
                (error.kind(), error.byte_offset()),
                (JsonErrorKind::IllegalCodepoint, 2)
            );
            assert_eq!(
                error.to_string(),
                format!("illegal unprintable codepoint in string: {control}")
            );
        }
        for (input, expected) in [
            (
                b"\"a\\\"b\"".as_slice(),
                JsonStringFast::Slow { start: 1, at: 2 },
            ),
            (b"\"abc\xc3\xa9\"", JsonStringFast::Slow { start: 1, at: 4 }),
        ] {
            let mut cursor = JsonCursor::new(input);
            assert_eq!(cursor.parse_string_fast().unwrap(), expected);
            assert_eq!(
                cursor.position,
                match expected {
                    JsonStringFast::Slow { at, .. } => at,
                    _ => unreachable!(),
                }
            );
        }
    }

    #[test]
    fn bio_read_j05_hex_digit_classes_and_error_position() {
        for (input, expected) in [
            (b"0123".as_slice(), 0x0123),
            (b"89ab", 0x89ab),
            (b"CDEF", 0xcdef),
            (b"aBcD", 0xabcd),
        ] {
            let mut cursor = JsonCursor::new(input);
            assert_eq!(cursor.read_hex().unwrap(), expected);
            assert_eq!(cursor.position, 4);
        }
        for index in 0..4 {
            let mut input = *b"1234";
            input[index] = b'g';
            let mut cursor = JsonCursor::new(&input);
            let error = cursor.read_hex().unwrap_err();
            assert_eq!(error.kind(), JsonErrorKind::InvalidUnicodeEscape);
            assert_eq!(error.byte_offset(), index + 1);
            assert_eq!(cursor.position, index + 1);
        }
        for length in 0..4 {
            let mut cursor = JsonCursor::new(&b"1234"[..length]);
            let error = cursor.read_hex().unwrap_err();
            assert_eq!(
                (error.kind(), error.byte_offset()),
                (JsonErrorKind::UnexpectedEnd, length)
            );
        }
    }

    #[test]
    fn bio_read_j05_utf8_width_boundaries_and_checked_capacity() {
        let cursor = JsonCursor::new(b"\\u");
        for (codepoint, expected) in [
            (0x00, &[0x00][..]),
            (0x7f, &[0x7f]),
            (0x80, &[0xc2, 0x80]),
            (0x7ff, &[0xdf, 0xbf]),
            (0x800, &[0xe0, 0xa0, 0x80]),
            (0xffff, &[0xef, 0xbf, 0xbf]),
            (0x10000, &[0xf0, 0x90, 0x80, 0x80]),
            (0x1fffff, &[0xf7, 0xbf, 0xbf, 0xbf]),
        ] {
            let mut output = [b'x'; 5];
            let written = cursor.write_utf8(codepoint, &mut output[1..]).unwrap();
            assert_eq!(written, expected.len());
            assert_eq!(&output[1..][..written], expected);
            assert_eq!(output[0], b'x');
            let mut bounded = [b'x'; 5];
            let error = cursor
                .write_utf8(codepoint, &mut bounded[1..expected.len()])
                .unwrap_err();
            assert_eq!(error.kind(), JsonErrorKind::OutOfMemory);
            assert_eq!(bounded, [b'x'; 5]);
        }
        let mut output = [b'x'; 4];
        let error = cursor.write_utf8(0x200000, &mut output).unwrap_err();
        assert_eq!(error.kind(), JsonErrorKind::InvalidUnicodeEscape);
        assert_eq!(output, [b'x'; 4]);
    }

    #[test]
    fn bio_read_j04_raw_number_lexemes_and_lookahead() {
        for lexeme in [
            "-0",
            "0",
            "1",
            "1234567890123456789012345678901234567890",
            "1.0",
            "-0.000",
            "9e0",
            "9E+003",
            "9e-003",
            "12.34E-56",
        ] {
            let input = format!("{lexeme},");
            let mut cursor = JsonCursor::new(input.as_bytes());
            let (span, kind) = cursor.parse_number_raw().unwrap();
            assert_eq!(kind, JsonNumberKind::Double);
            assert_eq!(&input.as_bytes()[span.clone()], lexeme.as_bytes());
            assert_eq!(span, 0..lexeme.len());
            assert_eq!(cursor.position, lexeme.len());
        }
        // Source consumes a zero as the entire integer part; outer state
        // subsequently rejects the extra digit, rather than this scanner.
        let mut cursor = JsonCursor::new(b"01,");
        assert_eq!(cursor.parse_number_raw().unwrap().0, 0..1);
        assert_eq!(cursor.position, 1);
    }

    #[test]
    fn bio_read_j04_invalid_number_and_eof_boundaries() {
        for (input, offset, kind) in [
            ("-", 1, JsonErrorKind::UnexpectedEnd),
            ("-x", 1, JsonErrorKind::InvalidNumber),
            ("1.", 2, JsonErrorKind::UnexpectedEnd),
            ("1.x", 2, JsonErrorKind::InvalidNumber),
            ("1e", 2, JsonErrorKind::UnexpectedEnd),
            ("1e+", 3, JsonErrorKind::UnexpectedEnd),
            ("1e-", 3, JsonErrorKind::UnexpectedEnd),
            ("1ex", 2, JsonErrorKind::MissingExponent),
            ("1e+x", 3, JsonErrorKind::MissingExponent),
            ("0", 1, JsonErrorKind::UnexpectedEnd),
            ("-0", 2, JsonErrorKind::UnexpectedEnd),
            ("123", 3, JsonErrorKind::UnexpectedEnd),
            ("1.2", 3, JsonErrorKind::UnexpectedEnd),
            ("1e2", 3, JsonErrorKind::UnexpectedEnd),
        ] {
            let mut cursor = JsonCursor::new(input.as_bytes());
            let error = cursor.parse_number_raw().unwrap_err();
            assert_eq!(
                (error.byte_offset(), error.kind()),
                (offset, kind),
                "{input}"
            );
            assert_eq!(cursor.position, 0, "{input}");
        }
    }

    #[test]
    fn bio_read_j03_literal_bounds_kinds_and_trailing_bytes() {
        for (input, kind, end) in [
            (b"null".as_slice(), JsonLiteralKind::Null, 4),
            (b"false", JsonLiteralKind::False, 5),
            (b"true", JsonLiteralKind::True, 4),
        ] {
            let mut cursor = JsonCursor::new(input);
            let result = match kind {
                JsonLiteralKind::Null => cursor.parse_null(),
                JsonLiteralKind::False => cursor.parse_false(),
                JsonLiteralKind::True => cursor.parse_true(),
            };
            assert_eq!(result.unwrap(), (end, kind));
            assert!(cursor.at_eof());
        }
        for (input, kind, end) in [
            (b"null,".as_slice(), JsonLiteralKind::Null, 4),
            (b"false]", JsonLiteralKind::False, 5),
            (b"true ", JsonLiteralKind::True, 4),
        ] {
            let mut cursor = JsonCursor::new(input);
            let result = match kind {
                JsonLiteralKind::Null => cursor.parse_null(),
                JsonLiteralKind::False => cursor.parse_false(),
                JsonLiteralKind::True => cursor.parse_true(),
            };
            assert_eq!(result.unwrap(), (end, kind));
            assert!(!cursor.at_eof());
            assert_eq!(cursor.position, end);
        }
    }

    #[test]
    fn bio_read_j03_truncated_and_misspelled_literals_report_start() {
        for (full, expected_kind) in [
            (b"null".as_slice(), JsonErrorKind::ExpectedNull),
            (b"false", JsonErrorKind::ExpectedFalse),
            (b"true", JsonErrorKind::ExpectedTrue),
        ] {
            for len in 0..full.len() {
                let mut cursor = JsonCursor::new(&full[..len]);
                let error = match expected_kind {
                    JsonErrorKind::ExpectedNull => cursor.parse_null().unwrap_err(),
                    JsonErrorKind::ExpectedFalse => cursor.parse_false().unwrap_err(),
                    JsonErrorKind::ExpectedTrue => cursor.parse_true().unwrap_err(),
                    _ => unreachable!(),
                };
                assert_eq!(error.kind(), JsonErrorKind::UnexpectedEnd);
                assert_eq!(error.byte_offset(), 0);
                assert_eq!(cursor.position, 0);
            }
        }
        for (input, expected_kind) in [
            (b"nulL".as_slice(), JsonErrorKind::ExpectedNull),
            (b"falsE", JsonErrorKind::ExpectedFalse),
            (b"truE", JsonErrorKind::ExpectedTrue),
        ] {
            let mut cursor = JsonCursor::new(input);
            let error = match expected_kind {
                JsonErrorKind::ExpectedNull => cursor.parse_null().unwrap_err(),
                JsonErrorKind::ExpectedFalse => cursor.parse_false().unwrap_err(),
                JsonErrorKind::ExpectedTrue => cursor.parse_true().unwrap_err(),
                _ => unreachable!(),
            };
            assert_eq!(error.kind(), expected_kind);
            assert_eq!(error.byte_offset(), 0);
            assert_eq!(cursor.position, 0);
        }
        // The pinned scanners check suffix bytes; their caller dispatches on
        // the first byte, so these direct private calls retain that contract.
        assert_eq!(
            JsonCursor::new(b"Null").parse_null().unwrap().1,
            JsonLiteralKind::Null
        );
    }

    #[test]
    fn bio_read_j02_four_json_whitespace_bytes_and_other_controls() {
        for whitespace in [b' ', b'\t', b'\r', b'\n'] {
            let input = [whitespace, b'x'];
            let mut cursor = JsonCursor::new(&input);
            assert_eq!(cursor.skip_whitespace(), Some(1));
            assert_eq!(cursor.position, 1);
            assert!(!cursor.at_eof());
        }
        for non_whitespace in [b'\x00', b'\x0b', b'\x0c', b'\x1f', b'x', 0x80] {
            let input = [non_whitespace, b'x'];
            let mut cursor = JsonCursor::new(&input);
            assert_eq!(cursor.skip_whitespace(), Some(0));
            assert_eq!(cursor.position, 0);
        }
        let mut cursor = JsonCursor::new(b" \t\r\nx");
        assert_eq!(cursor.skip_whitespace(), Some(4));
        assert_eq!(cursor.position, 4);
    }

    #[test]
    fn bio_read_j02_empty_and_every_eof_boundary() {
        let mut empty = JsonCursor::new(b"");
        assert!(empty.at_eof());
        assert!(empty.has_remaining_characters(0));
        assert!(!empty.has_remaining_characters(1));
        assert_eq!(empty.skip_whitespace(), None);

        for input in [b"".as_slice(), b" ", b"\n\r", b" \t\r\n"] {
            let mut cursor = JsonCursor::new(input);
            assert_eq!(cursor.skip_whitespace(), None);
            assert_eq!(cursor.position, input.len());
            assert!(cursor.at_eof());
            assert!(cursor.has_remaining_characters(0));
            assert!(!cursor.has_remaining_characters(1));
        }
        for position in 0..=4 {
            let mut cursor = JsonCursor::new(b"abcd");
            cursor.position = position;
            assert_eq!(cursor.at_eof(), position == 4);
            for width in 0..=5 {
                assert_eq!(
                    cursor.has_remaining_characters(width),
                    4 - position >= width
                );
            }
        }
    }

    #[test]
    fn bio_read_j01_source_error_categories_and_argument_rendering() {
        let cases = [
            (JsonErrorKind::NoError, "no error"),
            (JsonErrorKind::OutOfMemory, "out of memory"),
            (JsonErrorKind::UnexpectedEnd, "unexpected end of input"),
            (JsonErrorKind::MissingRootElement, "missing root element"),
            (
                JsonErrorKind::BadRoot,
                "document root must be object or array",
            ),
            (JsonErrorKind::ExpectedComma, "expected ,"),
            (JsonErrorKind::MissingObjectKey, "missing object key"),
            (JsonErrorKind::ExpectedColon, "expected :"),
            (JsonErrorKind::ExpectedEndOfInput, "expected end of input"),
            (JsonErrorKind::UnexpectedComma, "unexpected comma"),
            (JsonErrorKind::ExpectedValue, "expected value"),
            (JsonErrorKind::ExpectedNull, "expected 'null'"),
            (JsonErrorKind::ExpectedFalse, "expected 'false'"),
            (JsonErrorKind::ExpectedTrue, "expected 'true'"),
            (JsonErrorKind::InvalidNumber, "invalid number"),
            (JsonErrorKind::MissingExponent, "missing exponent"),
            (
                JsonErrorKind::InvalidUnicodeEscape,
                "invalid character in unicode escape",
            ),
            (
                JsonErrorKind::UnexpectedEndOfUtf16,
                "unexpected end of input during UTF-16 surrogate pair",
            ),
            (JsonErrorKind::ExpectedU, "expected \\u"),
            (
                JsonErrorKind::InvalidUtf16TrailSurrogate,
                "invalid UTF-16 trail surrogate",
            ),
            (JsonErrorKind::UnknownEscape, "unknown escape"),
            (JsonErrorKind::InvalidUtf8, "invalid UTF-8"),
        ];
        for (kind, expected) in cases {
            let error = JsonError::at(b"x", 0, kind, None);
            assert_eq!(error.kind(), kind);
            assert_eq!(error.to_string(), expected);
        }
        let error = JsonError::at(b"\x01", 0, JsonErrorKind::IllegalCodepoint, Some(1));
        assert_eq!(
            error.to_string(),
            "illegal unprintable codepoint in string: 1"
        );
    }

    #[test]
    fn bio_read_j01_byte_offsets_and_crlf_line_columns_follow_source() {
        let input = b"a\r\nb\nc\rd";
        let expected = [
            (0, 1, 1),
            (1, 1, 2),
            (2, 2, 1),
            (3, 2, 1),
            (4, 2, 2),
            (5, 3, 1),
            (6, 3, 2),
            (7, 4, 1),
            (8, 4, 2),
        ];
        for (offset, line, column) in expected {
            let error = JsonError::at(input, offset, JsonErrorKind::ExpectedValue, None);
            assert_eq!(
                (error.byte_offset(), error.line(), error.column()),
                (offset, line, column)
            );
        }
        let unicode = JsonError::at("éx".as_bytes(), 2, JsonErrorKind::InvalidUtf8, None);
        assert_eq!((unicode.line(), unicode.column()), (1, 3));
    }
}
