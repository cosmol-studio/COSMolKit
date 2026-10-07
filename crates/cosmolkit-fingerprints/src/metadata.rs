//! Shared source property-tree fingerprint metadata transport.
use std::fmt;
#[derive(Debug)]
pub enum FingerprintJsonError {
    Parse(serde_json::Error),
    Invalid(String),
    UnsupportedComponent {
        component: &'static str,
        source_type: String,
    },
}
impl fmt::Display for FingerprintJsonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(e) => e.fmt(f),
            Self::Invalid(s) => f.write_str(s),
            Self::UnsupportedComponent {
                component,
                source_type,
            } => write!(
                f,
                "unsupported source generator component {component}: {source_type}"
            ),
        }
    }
}
impl std::error::Error for FingerprintJsonError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Parse(e) => Some(e),
            Self::Invalid(_) | Self::UnsupportedComponent { .. } => None,
        }
    }
}
/// Source property_tree nodes retain ordered, duplicate children and text leaves.
#[derive(Debug)]
pub(crate) struct SourceNode {
    data: String,
    children: Vec<(String, SourceNode)>,
    object: bool,
    array: bool,
    first: std::collections::BTreeMap<String, usize>,
}
impl SourceNode {
    pub(crate) fn get(&self, key: &str) -> Option<&Self> {
        // Boost1.88❗✔️: const_assoc_iterator el = find(fragment);
        // Boost1.88❗✔️: if(el == not_found()) {
        // Boost1.88❗✔️:     return 0;
        // Behavior: find the first equivalent key, never serde's last duplicate.
        // Complexity: source indexed lookup, O(log children).
        self.first.get(key).map(|&index| &self.children[index].1)
    }
    pub(crate) fn is_object(&self) -> bool {
        self.object
    }
    pub(crate) fn as_str(&self) -> Option<&str> {
        Some(&self.data)
    }
    pub(crate) fn children(&self) -> impl Iterator<Item = &Self> {
        self.children.iter().map(|(_, child)| child)
    }
}
pub(crate) fn parse_object(json: &str) -> Result<SourceNode, FingerprintJsonError> {
    // Boost1.88❗❌: void on_number(Range code_units) {
    // Boost1.88❗❌:     new_value().assign(code_units.begin(), code_units.end());
    // Boost1.88❗❌: }
    // Boost1.88❗❌: void on_null() {
    // Boost1.88❗❌:     new_value() = constants::null_value<char_type>();
    // Boost1.88❗❌: }
    // Boost1.88❗❌: void on_boolean(bool b) {
    // Boost1.88❗❌:     new_value() = b ? constants::true_value<char_type>()
    // Boost1.88❗❌:                     : constants::false_value<char_type>();
    // Boost1.88❗❌: }
    // Behavior: source callbacks retain the original numeric spelling, decoded
    // strings, object duplicates and insertion order; root need not be object.
    // Complexity: O(JSON) plus indexed children. The independent RawValue
    // validation pass and duplicated index keys add costs versus source's one
    // parser/indexed nodes; do not claim complexity/performance equivalence.
    let _: &serde_json::value::RawValue =
        serde_json::from_str(json).map_err(FingerprintJsonError::Parse)?;
    // RawValue validates all syntax without materializing or normalizing numbers.
    // This cursor only reads the validated token stream, once, monotonically.
    fn ws(bytes: &[u8], pos: &mut usize) {
        while bytes.get(*pos).is_some_and(u8::is_ascii_whitespace) {
            *pos += 1;
        }
    }
    fn text(json: &str, pos: &mut usize) -> String {
        let start = *pos;
        *pos += 1;
        while json.as_bytes()[*pos] != b'"' {
            if json.as_bytes()[*pos] == b'\\' {
                *pos += 1;
            }
            *pos += 1;
        }
        *pos += 1;
        serde_json::from_str(&json[start..*pos]).expect("validated JSON string")
    }
    fn node(json: &str, pos: &mut usize) -> SourceNode {
        let bytes = json.as_bytes();
        ws(bytes, pos);
        let token = bytes[*pos];
        let mut value = SourceNode {
            data: String::new(),
            children: Vec::new(),
            object: token == b'{',
            array: token == b'[',
            first: Default::default(),
        };
        if value.object || value.array {
            let end = if value.object { b'}' } else { b']' };
            *pos += 1;
            ws(bytes, pos);
            while bytes[*pos] != end {
                let key = if value.object {
                    let key = text(json, pos);
                    ws(bytes, pos);
                    *pos += 1;
                    key
                } else {
                    String::new()
                };
                let child = node(json, pos);
                if value.object {
                    value
                        .first
                        .entry(key.clone())
                        .or_insert(value.children.len());
                }
                value.children.push((key, child));
                ws(bytes, pos);
                if bytes[*pos] == b',' {
                    *pos += 1;
                    ws(bytes, pos);
                } else {
                    break;
                }
            }
            *pos += 1;
        } else if token == b'"' {
            value.data = text(json, pos);
        } else {
            let start = *pos;
            while bytes
                .get(*pos)
                .is_some_and(|b| !b.is_ascii_whitespace() && !b",]}".contains(b))
            {
                *pos += 1;
            }
            value.data = json[start..*pos].into();
        }
        value
    }
    Ok(node(json, &mut 0))
}
pub(crate) fn common_arguments_string(
    count_simulation: bool,
    fp_size: u32,
    bits_per_feature: u32,
    include_chirality: bool,
) -> String {
    // RDKit source: FingerprintGenerator.cpp lines 50-54
    // RDKit❗✔️: std::string FingerprintArguments::commonArgumentsString() const {
    // RDKit❗✔️:   return "Common arguments : countSimulation=" +
    // RDKit❗✔️:          std::to_string(df_countSimulation) +
    // RDKit❗✔️:          " fpSize=" + std::to_string(d_fpSize) +
    // RDKit❗✔️:          " bitsPerFeature=" + std::to_string(d_numBitsPerFeature) +
    // RDKit❗✔️:          " includeChirality=" + std::to_string(df_includeChirality);
    // RDKit❗✔️: }
    format!(
        "Common arguments : countSimulation={} fpSize={} bitsPerFeature={} includeChirality={}",
        count_simulation as u8, fp_size, bits_per_feature, include_chirality as u8
    )
}

pub(crate) fn common_arguments_json(
    count_simulation: bool,
    fp_size: u32,
    bits_per_feature: u32,
    include_chirality: bool,
    count_bounds: &[u32],
) -> String {
    // RDKit source: FingerprintGenerator.cpp lines 58-71
    // RDKit❗✔️: void FingerprintArguments::toJSON(boost::property_tree::ptree &pt) const {
    // RDKit❗✔️:   pt.put("countSimulation", df_countSimulation);
    // RDKit❗✔️:   pt.put("fpSize", d_fpSize);
    // RDKit❗✔️:   pt.put("numBitsPerFeature", d_numBitsPerFeature);
    // RDKit❗✔️:   pt.put("includeChirality", df_includeChirality);
    // RDKit❗✔️:
    // RDKit❗✔️:   boost::property_tree::ptree countBoundsNode;
    // RDKit❗✔️:   for (const auto &bound : d_countBounds) {
    // RDKit❗✔️:     boost::property_tree::ptree boundNode;
    // RDKit❗✔️:     boundNode.put("", bound);
    // RDKit❗✔️:     countBoundsNode.push_back(std::make_pair("", boundNode));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   pt.add_child("countBounds", countBoundsNode);
    // RDKit❗✔️: }
    let count_bounds = count_bounds
        .iter()
        .map(|bound| format!("\"{bound}\""))
        .collect::<Vec<_>>()
        .join(",");
    // Boost property_tree stores every leaf as text, so write_json quotes
    // numeric and boolean values even though fromJSON reads typed values.
    format!(
        "{{\"countSimulation\":\"{}\",\"fpSize\":\"{}\",\"numBitsPerFeature\":\"{}\",\"includeChirality\":\"{}\",\"countBounds\":{}}}",
        count_simulation,
        fp_size,
        bits_per_feature,
        include_chirality,
        if count_bounds.is_empty() {
            "\"\"".to_string()
        } else {
            format!("[{count_bounds}]")
        }
    )
}

pub(crate) fn common_arguments_from_json(
    value: &SourceNode,
    count_simulation: &mut bool,
    fp_size: &mut u32,
    bits_per_feature: &mut u32,
    include_chirality: &mut bool,
    count_bounds: &mut Vec<u32>,
) -> Result<(), FingerprintJsonError> {
    // RDKit source: FingerprintGenerator.cpp lines 73-90
    // RDKit❗✔️: void FingerprintArguments::fromJSON(const boost::property_tree::ptree &pt) {
    // RDKit❗✔️:   df_countSimulation = pt.get<bool>("countSimulation", df_countSimulation);
    // RDKit❗✔️:   d_fpSize = pt.get<std::uint32_t>("fpSize", d_fpSize);
    // RDKit❗✔️:   d_numBitsPerFeature =
    // RDKit❗✔️:       pt.get<std::uint32_t>("numBitsPerFeature", d_numBitsPerFeature);
    // RDKit❗✔️:   df_includeChirality = pt.get<bool>("includeChirality", df_includeChirality);
    // RDKit❗✔️:
    // RDKit❗✔️:   d_countBounds.clear();
    // RDKit❗✔️:   auto countBoundsNode = pt.get_child_optional("countBounds");
    // RDKit❗✔️:   if (countBoundsNode) {
    // RDKit❗✔️:     for (const auto &boundNode : *countBoundsNode) {
    // RDKit❗✔️:       d_countBounds.push_back(boundNode.second.get_value<std::uint32_t>());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    *count_simulation = bool_or(value, "countSimulation", *count_simulation);
    *fp_size = u32_or(value, "fpSize", *fp_size);
    *bits_per_feature = u32_or(value, "numBitsPerFeature", *bits_per_feature);
    *include_chirality = bool_or(value, "includeChirality", *include_chirality);
    count_bounds.clear();
    if let Some(field) = value.get("countBounds") {
        for bound in field.children() {
            // No default overload: a bad child remains a typed conversion error.
            count_bounds.push(json_value_as_u32("countBounds entry", bound)?);
        }
    }
    Ok(())
}
pub(crate) fn bool_or(node: &SourceNode, key: &str, current: bool) -> bool {
    // Boost1.88❗✔️: return get_optional<Type>(path, tr).get_value_or(default_value);
    // Source-defined fallback applies to absent AND untranslatable scalar nodes.
    node.get(key)
        .and_then(|v| json_value_as_bool(key, v).ok())
        .unwrap_or(current)
}
pub(crate) fn u32_or(node: &SourceNode, key: &str, current: u32) -> u32 {
    // Boost1.88❗✔️: return get_optional<Type>(path, tr).get_value_or(default_value);
    // Source-defined fallback; countBounds get_value without default stays strict.
    node.get(key)
        .and_then(|v| json_value_as_u32(key, v).ok())
        .unwrap_or(current)
}

pub(crate) fn json_value_as_bool(
    name: &str,
    value: &SourceNode,
) -> Result<bool, FingerprintJsonError> {
    // Boost1.88❗✔️: s >> e;
    // Boost1.88❗✔️: if(s.fail()) {
    // Boost1.88❗✔️:     // Try again in word form.
    // Boost1.88❗✔️:     s.clear();
    // Boost1.88❗✔️:     s.setf(std::ios_base::boolalpha);
    // Boost1.88❗✔️:     s >> e;
    // Boost1.88❗✔️: }
    // Boost1.88❗✔️: if(iss.fail() || iss.bad() || iss.get() != Traits::eof()) {
    // Boost1.88❗✔️:     return boost::optional<E>();
    // Boost1.88❗✔️: }
    // Behavior: numeric bool extraction accepts signed 0/1, then lower-case
    // words; complete consumption with classic locale whitespace required.
    // Complexity: linear in one scalar text, no graph/old configuration copies.
    let text = value.data.trim_matches(|c: char| c.is_ascii_whitespace());
    // Numeric extraction consumes a sign/digit prefix before setting failbit.
    // boolalpha retries from that CURRENT stream position, without seekg(0).
    // Thus "2true" and "+true" succeed, while "1true" fails the final EOF check.
    let bytes = text.as_bytes();
    let mut consumed = usize::from(bytes.first().is_some_and(|b| *b == b'+' || *b == b'-'));
    let digits_start = consumed;
    while bytes.get(consumed).is_some_and(u8::is_ascii_digit) {
        consumed += 1;
    }
    if consumed > digits_start {
        if let Some(number) = signed_decimal(&text[..consumed]) {
            if number == 0 || number == 1 {
                if consumed == text.len() {
                    return Ok(number == 1);
                }
                return Err(FingerprintJsonError::Invalid(format!(
                    "{name} must be a boolean"
                )));
            }
        }
    }
    let retry = text[consumed..].trim_matches(|c: char| c.is_ascii_whitespace());
    if retry == "true" {
        return Ok(true);
    }
    if retry == "false" {
        return Ok(false);
    }
    Err(FingerprintJsonError::Invalid(format!(
        "{name} must be a boolean"
    )))
}
fn signed_decimal(text: &str) -> Option<i64> {
    // Full decimal extraction, optional sign; stream decimal has no 0x shortcut.
    let digits = text
        .strip_prefix('+')
        .or_else(|| text.strip_prefix('-'))
        .unwrap_or(text);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    text.parse().ok()
}
pub(crate) fn json_value_as_u32(
    name: &str,
    value: &SourceNode,
) -> Result<u32, FingerprintJsonError> {
    // Boost1.88❗✔️: s >> e;
    // Boost1.88❗✔️: if(!s.eof()) {
    // Boost1.88❗✔️:     s >> std::ws;
    // Boost1.88❗✔️: }
    // Boost1.88❗✔️: if(iss.fail() || iss.bad() || iss.get() != Traits::eof()) {
    // Boost1.88❗✔️:     return boost::optional<E>();
    // Boost1.88❗✔️: }
    // Behavior: unsigned C++ decimal extraction accepts sign and wraps negative
    // magnitudes within UINT32_MAX; larger magnitudes fail, including negative.
    // Complexity: one scalar scan, fixed-size arithmetic; source-defined cast.
    let text = value.data.trim_matches(|c: char| c.is_ascii_whitespace());
    let digits = text
        .strip_prefix('+')
        .or_else(|| text.strip_prefix('-'))
        .unwrap_or(text);
    if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
        if let Ok(magnitude) = digits.parse::<u32>() {
            return Ok(if text.starts_with('-') {
                magnitude.wrapping_neg()
            } else {
                magnitude
            });
        }
    }
    Err(FingerprintJsonError::Invalid(format!(
        "{name} must be a 32-bit integer"
    )))
}

/// Append one counted property-tree string leaf, with source byte escaping.
pub(crate) fn append_json_byte_string(
    output: &mut cosmolkit_model::PropertyText,
    value: &cosmolkit_model::PropertyText,
) {
    // Boost1.85❗✔️:     template<class Ch>
    // Boost1.85❗✔️:     std::basic_string<Ch> create_escapes(const std::basic_string<Ch> &s)
    // Boost1.85❗✔️:     {
    // Boost1.85❗✔️:         std::basic_string<Ch> result;
    // Boost1.85❗✔️:         typename std::basic_string<Ch>::const_iterator b = s.begin();
    // Boost1.85❗✔️:         typename std::basic_string<Ch>::const_iterator e = s.end();
    // Boost1.85❗✔️:         while (b != e)
    // Boost1.85❗✔️:         {
    // Boost1.85❗✔️:             typedef typename make_unsigned<Ch>::type UCh;
    // Boost1.85❗✔️:             UCh c(*b);
    // Boost1.85❗✔️:             // This assumes an ASCII superset. But so does everything in PTree.
    // Boost1.85❗✔️:             // We escape everything outside ASCII, because this code can't
    // Boost1.85❗✔️:             // handle high unicode characters.
    // Boost1.85❗✔️:             if (c == 0x20 || c == 0x21 || (c >= 0x23 && c <= 0x2E) ||
    // Boost1.85❗✔️:                 (c >= 0x30 && c <= 0x5B) || (c >= 0x5D && c <= 0xFF))
    // Boost1.85❗✔️:                 result += *b;
    // Boost1.85❗✔️:             else if (*b == Ch('\b')) result += Ch('\\'), result += Ch('b');
    // Boost1.85❗✔️:             else if (*b == Ch('\f')) result += Ch('\\'), result += Ch('f');
    // Boost1.85❗✔️:             else if (*b == Ch('\n')) result += Ch('\\'), result += Ch('n');
    // Boost1.85❗✔️:             else if (*b == Ch('\r')) result += Ch('\\'), result += Ch('r');
    // Boost1.85❗✔️:             else if (*b == Ch('\t')) result += Ch('\\'), result += Ch('t');
    // Boost1.85❗✔️:             else if (*b == Ch('/')) result += Ch('\\'), result += Ch('/');
    // Boost1.85❗✔️:             else if (*b == Ch('"'))  result += Ch('\\'), result += Ch('"');
    // Boost1.85❗✔️:             else if (*b == Ch('\\')) result += Ch('\\'), result += Ch('\\');
    // Boost1.85❗✔️:             else
    // Boost1.85❗✔️:             {
    // Boost1.85❗✔️:                 const char *hexdigits = "0123456789ABCDEF";
    // Boost1.85❗✔️:                 unsigned long u = (std::min)(static_cast<unsigned long>(
    // Boost1.85❗✔️:                                                  static_cast<UCh>(*b)),
    // Boost1.85❗✔️:                                              0xFFFFul);
    // Boost1.85❗✔️:                 unsigned long d1 = u / 4096; u -= d1 * 4096;
    // Boost1.85❗✔️:                 unsigned long d2 = u / 256; u -= d2 * 256;
    // Boost1.85❗✔️:                 unsigned long d3 = u / 16; u -= d3 * 16;
    // Boost1.85❗✔️:                 unsigned long d4 = u;
    // Boost1.85❗✔️:                 result += Ch('\\'); result += Ch('u');
    // Boost1.85❗✔️:                 result += Ch(hexdigits[d1]); result += Ch(hexdigits[d2]);
    // Boost1.85❗✔️:                 result += Ch(hexdigits[d3]); result += Ch(hexdigits[d4]);
    // Boost1.85❗✔️:             }
    // Boost1.85❗✔️:             ++b;
    // Boost1.85❗✔️:         }
    // Boost1.85❗✔️:         return result;
    // Boost1.85❗✔️:     }
    // Behavior: Boost's Ch=char specialization converts each char to unsigned
    // before its ASCII-superset range check. Bytes 0x80..0xff remain bytes;
    // control characters and slash use the exact source escape spelling.
    // Complexity: one linear pass into the existing output, no decoding,
    // secondary Unicode representation or intermediate escaped allocation.
    // Output quotes are the source write_json_helper leaf delimiters.
    // Boost1.85❗✔️:             stream << Ch('"') << data << Ch('"');
    output.push_byte(b'"');
    for &byte in value.as_bytes() {
        match byte {
            0x20 | 0x21 | 0x23..=0x2e | 0x30..=0x5b | 0x5d..=0xff => output.push_byte(byte),
            b'\x08' => output.extend_bytes(b"\\b"),
            b'\x0c' => output.extend_bytes(b"\\f"),
            b'\n' => output.extend_bytes(b"\\n"),
            b'\r' => output.extend_bytes(b"\\r"),
            b'\t' => output.extend_bytes(b"\\t"),
            b'/' => output.extend_bytes(b"\\/"),
            b'"' => output.extend_bytes(b"\\\""),
            b'\\' => output.extend_bytes(b"\\\\"),
            _ => {
                output.extend_bytes(b"\\u00");
                output.push_byte(b"0123456789ABCDEF"[(byte >> 4) as usize]);
                output.push_byte(b"0123456789ABCDEF"[(byte & 15) as usize]);
            }
        }
    }
    output.push_byte(b'"');
}

#[cfg(test)]
mod counted_json_leaf_regressions {
    use super::append_json_byte_string;
    use cosmolkit_model::PropertyText;

    #[test]
    fn property_tree_char_leaf_retains_opaque_bytes_and_source_escape_spelling() {
        for (input, expected) in [
            (b"".as_slice(), b"\"\"".as_slice()),
            (
                b"\0\x01\x1f".as_slice(),
                b"\"\\u0000\\u0001\\u001F\"".as_slice(),
            ),
            (
                b"\x08\x0c\n\r\t".as_slice(),
                b"\"\\b\\f\\n\\r\\t\"".as_slice(),
            ),
            (b"/\"\\".as_slice(), b"\"\\/\\\"\\\\\"".as_slice()),
            (b"\x7f\x80\xff".as_slice(), b"\"\x7f\x80\xff\"".as_slice()),
            (b"R\0S\xff".as_slice(), b"\"R\\u0000S\xff\"".as_slice()),
        ] {
            let value = PropertyText::from_bytes(input);
            let before = value.clone();
            let mut output = PropertyText::from("prefix");
            append_json_byte_string(&mut output, &value);
            assert_eq!(&output.as_bytes()[6..], expected);
            assert_eq!(value, before);
        }
    }
}
