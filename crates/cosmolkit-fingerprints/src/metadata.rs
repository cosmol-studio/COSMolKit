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
